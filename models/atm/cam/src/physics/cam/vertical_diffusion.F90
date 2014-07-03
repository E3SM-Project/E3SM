module vertical_diffusion

  !----------------------------------------------------------------------------------------------------- !
  ! Module to compute vertical diffusion of momentum,  moisture, trace constituents                      !
  ! and static energy. Separate modules compute                                                          !
  !   1. stresses associated with turbulent flow over orography                                          !
  !      ( turbulent mountain stress )                                                                   !
  !   2. eddy diffusivities, including nonlocal tranport terms                                           !
  !   3. molecular diffusivities                                                                         !
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      !
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !  vertical_diffusion_init      Initializes vertical diffustion constants and modules                  !
  !        init_molec_diff        Initializes molecular diffusivity module                               !
  !        init_eddy_diff         Initializes eddy diffusivity module (includes PBL)                     !
  !        init_tms               Initializes turbulent mountain stress module                           !
  !        init_vdiff             Initializes diffusion solver module                                    !
  !  vertical_diffusion_ts_init   Time step initialization (only used for upper boundary condition)      !
  !  vertical_diffusion_tend      Computes vertical diffusion tendencies                                 !
  !        compute_tms            Computes turbulent mountain stresses                                   !
  !        compute_eddy_diff      Computes eddy diffusivities and countergradient terms                  !
  !        compute_vdiff          Solves vertical diffusion equations, including molecular diffusivities !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! J. Rosinski : Jun. 1992                                                                              !
  ! J. McCaa    : Sep. 2004                                                                              !
  ! S. Park     : Aug. 2006, Dec. 2008. Jan. 2010                                                        ! 
  !----------------------------------------------------------------------------------------------------- !

  use shr_kind_mod,     only : r8 => shr_kind_r8, i4=> shr_kind_i4
  use ppgrid,           only : pcols, pver, pverp
  use constituents,     only : pcnst, qmin, cnst_get_ind
  use diffusion_solver, only : vdiff_selector
  use abortutils,       only : endrun
  use error_messages,   only : handle_errmsg
  use physconst,        only :          &
                               cpair  , &     ! Specific heat of dry air
                               gravit , &     ! Acceleration due to gravity
                               rair   , &     ! Gas constant for dry air
                               zvir   , &     ! rh2o/rair - 1
                               latvap , &     ! Latent heat of vaporization
                               latice , &     ! Latent heat of fusion
                               karman , &     ! von Karman constant
                               mwdry  , &     ! Molecular weight of dry air
                               avogad , &     ! Avogadro's number
                               boltz  , &     ! Boltzman's constant
                               tms_orocnst,&  ! turbulent mountain stress parameter
                               tms_z0fac      ! Factor determining z_0 from orographic standard deviation [no unit]
  use cam_history,      only : fieldname_len
  use perf_mod
  use cam_logfile,      only : iulog
  use ref_pres,         only : do_molec_diff
  use phys_control,     only : phys_getopts, waccmx_is
  use time_manager,     only : is_first_step

  implicit none
  private      
  save
  
  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public vd_readnl
  public vd_register                                   ! Register multi-time-level variables with physics buffer
  public vertical_diffusion_init                       ! Initialization
  public vertical_diffusion_ts_init                    ! Time step initialization (only used for upper boundary condition)
  public vertical_diffusion_tend                       ! Full vertical diffusion routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  character(len=16)    :: eddy_scheme                  ! Default set in phys_control.F90, use namelist to change
                                                       !     'HB'       = Holtslag and Boville (default)
                                                       !     'HBR'      = Holtslag and Boville and Rash 
                                                       !     'diag_TKE' = Bretherton and Park ( UW Moist Turbulence Scheme )
  integer, parameter   :: nturb = 5                    ! Number of iterations for solution ( when 'diag_TKE' scheme is selected )
  logical, parameter   :: wstarent = .true.            ! Use wstar (.true.) or TKE (.false.) entrainment closure
                                                       ! ( when 'diag_TKE' scheme is selected )
  logical              :: do_pseudocon_diff = .false.  ! If .true., do pseudo-conservative variables diffusion

  character(len=16)    :: shallow_scheme               ! For checking compatibility between eddy diffusion
                                                       ! and shallow convection schemes
                                                       !     'Hack'     = Hack Shallow Convection Scheme
                                                       !     'UW'       = Park and Bretherton ( UW Shallow Convection Scheme )
  character(len=16)    :: microp_scheme                ! Microphysics scheme

  type(vdiff_selector) :: fieldlist_wet                ! Logical switches for moist mixing ratio diffusion
  type(vdiff_selector) :: fieldlist_dry                ! Logical switches for dry mixing ratio diffusion
  type(vdiff_selector) :: fieldlist_molec              ! Logical switches for molecular diffusion
  integer              :: ntop                         ! Top interface level to which vertical diffusion is applied ( = 1 ).
  integer              :: nbot                         ! Bottom interface level to which vertical diffusion is applied ( = pver ).
  integer              :: tke_idx, kvh_idx, kvm_idx    ! TKE and eddy diffusivity indices for fields in the physics buffer
  integer              :: kvt_idx                      ! Index for kinematic molecular conductivity
  integer              :: turbtype_idx, smaw_idx       ! Turbulence type and instability functions
  integer              :: tauresx_idx, tauresy_idx     ! Redisual stress for implicit surface stress

  character(len=fieldname_len) :: vdiffnam(pcnst)      ! Names of vertical diffusion tendencies
  integer              :: ixcldice, ixcldliq           ! Constituent indices for cloud liquid and ice water
  integer              :: ixnumice, ixnumliq


  logical              :: history_amwg                 ! output the variables used by the AMWG diag package
  logical              :: history_eddy                 ! output the eddy variables
  logical              :: history_budget               ! Output tendencies and state variables for CAM4 T, qv, ql, qi
  integer              :: history_budget_histfile_num  ! output history file number for budget fields
  logical              :: history_waccm                ! output variables of interest for WACCM runs

  integer              :: qrl_idx    = 0               ! pbuf index 
  integer              :: wsedl_idx  = 0               ! pbuf index

  integer              :: pblh_idx, tpert_idx, qpert_idx

  real(r8)             :: kv_top_pressure              ! Pressure defining the bottom of the upper atmosphere for kvh scaling (Pa)
  real(r8)             :: kv_top_scale                 ! Eddy diffusivity scale factor for upper atmosphere
  real(r8)             :: kv_freetrop_scale            ! Eddy diffusivity scale factor for the free troposphere
  real(r8)             :: eddy_lbulk_max               ! Maximum master length for diag_TKE
  real(r8)             :: eddy_leng_max                ! Maximum dissipation length for diag_TKE
  real(r8)             :: eddy_max_bot_pressure        ! Bottom pressure level (hPa) for eddy_leng_max
  logical              :: diff_cnsrv_mass_check        ! do mass conservation check
  logical              :: do_tms                       ! switch for turbulent mountain stress
  logical              :: do_iss                       ! switch for implicit turbulent surface stress
  logical              :: prog_modal_aero = .false.    ! set true if prognostic modal aerosols are present
  integer              :: pmam_ncnst = 0               ! number of prognostic modal aerosol constituents
  integer, allocatable :: pmam_cnst_idx(:)             ! constituent indices of prognostic modal aerosols

contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  subroutine vd_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use spmd_utils,      only: masterproc
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'vd_readnl'
  
    namelist /vert_diff_nl/ kv_top_pressure, kv_top_scale, kv_freetrop_scale, eddy_lbulk_max, eddy_leng_max, &
         eddy_max_bot_pressure, diff_cnsrv_mass_check, do_iss
    !-----------------------------------------------------------------------------
  
    if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'vert_diff_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, vert_diff_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname // ':: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if
  
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(kv_top_pressure,                 1 , mpir8,   0, mpicom)
    call mpibcast(kv_top_scale,                    1 , mpir8,   0, mpicom)
    call mpibcast(kv_freetrop_scale,               1 , mpir8,   0, mpicom)
    call mpibcast(eddy_lbulk_max,                  1 , mpir8,   0, mpicom)
    call mpibcast(eddy_leng_max,                   1 , mpir8,   0, mpicom)
    call mpibcast(eddy_max_bot_pressure,           1 , mpir8,   0, mpicom)
    call mpibcast(diff_cnsrv_mass_check,           1 , mpilog,  0, mpicom)
    call mpibcast(do_iss,                          1 , mpilog,  0, mpicom)
#endif

  end subroutine vd_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vd_register()

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    use physics_buffer,      only : pbuf_add_field, dtype_r8, dtype_i4

    ! Get eddy_scheme setting from phys_control.F90

    call phys_getopts( eddy_scheme_out          =          eddy_scheme, & 
                       shallow_scheme_out       =       shallow_scheme, &
                       microp_scheme_out        =        microp_scheme, &
                       do_tms_out               =               do_tms)

    ! Add fields to physics buffer

    ! kvt is used by gw_drag.  only needs physpkg scope.
    call pbuf_add_field('kvt', 'physpkg', dtype_r8, (/pcols,pverp/), kvt_idx) 


    call pbuf_add_field('pblh',     'global', dtype_r8, (/pcols/),        pblh_idx)
    call pbuf_add_field('tke',      'global', dtype_r8, (/pcols, pverp/), tke_idx) 
    call pbuf_add_field('kvh',      'global', dtype_r8, (/pcols, pverp/), kvh_idx) 
    call pbuf_add_field('kvm',      'global', dtype_r8, (/pcols, pverp/), kvm_idx ) 
    call pbuf_add_field('turbtype', 'global', dtype_i4, (/pcols, pverp/), turbtype_idx)
    call pbuf_add_field('smaw',     'global', dtype_r8, (/pcols, pverp/), smaw_idx) 

    call pbuf_add_field('tauresx',  'global', dtype_r8, (/pcols/),        tauresx_idx)
    call pbuf_add_field('tauresy',  'global', dtype_r8, (/pcols/),        tauresy_idx)

    call pbuf_add_field('tpert', 'global', dtype_r8, (/pcols/),                       tpert_idx)
    call pbuf_add_field('qpert', 'global', dtype_r8, (/pcols,pcnst/),                 qpert_idx)

  end subroutine vd_register

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_init(pbuf2d)

    !------------------------------------------------------------------!
    ! Initialization of time independent fields for vertical diffusion !
    ! Calls initialization routines for subsidiary modules             !
    !----------------------------------------------------------------- !

    use cam_history,       only : addfld, add_default, phys_decomp
    use eddy_diff,         only : init_eddy_diff
    use hb_diff,           only : init_hb_diff
    use molec_diff,        only : init_molec_diff
    use trb_mtn_stress,    only : init_tms
    use diffusion_solver,  only : init_vdiff, new_fieldlist_vdiff, vdiff_select
    use constituents,      only : cnst_get_ind, cnst_get_type_byind, cnst_name, cnst_get_molec_byind
    use spmd_utils,        only : masterproc
    use ref_pres,          only : ntop_molec, nbot_molec, press_lim_idx, pref_mid
    use physics_buffer,    only : pbuf_set_field, pbuf_get_index, physics_buffer_desc
    use rad_constituents,  only : rad_cnst_get_info, rad_cnst_get_mode_num_idx, &
                                  rad_cnst_get_mam_mmr_idx

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    character(128) :: errstring   ! Error status for init_vdiff
    integer        :: ntop_eddy   ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
    integer        :: nbot_eddy   ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
    integer        :: k           ! Vertical loop index

    real(r8), parameter :: ntop_eddy_pres = 1.e-5_r8 ! Pressure below which eddy diffusion is not done in WACCM-X. (Pa)

    integer :: im, l, m, nmodes, nspec

    ! ----------------------------------------------------------------- !

    if (masterproc) then
       write(iulog,*)'Initializing vertical diffusion (vertical_diffusion_init)'
    end if

    ! ----------------------------------------------------------------- !
    ! Get indices of cloud liquid and ice within the constituents array !
    ! ----------------------------------------------------------------- !

    call cnst_get_ind( 'CLDLIQ', ixcldliq )
    call cnst_get_ind( 'CLDICE', ixcldice )
    if( microp_scheme == 'MG' ) then
        call cnst_get_ind( 'NUMLIQ', ixnumliq )
        call cnst_get_ind( 'NUMICE', ixnumice )
    endif

    ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)
    if (prog_modal_aero) then

       ! Get the constituent indices of the number and mass mixing ratios of the modal 
       ! aerosols.
       !
       ! N.B. - This implementation assumes that the prognostic modal aerosols are 
       !        impacting the climate calculation (i.e., can get info from list 0).
       ! 

       ! First need total number of mam constituents
       call rad_cnst_get_info(0, nmodes=nmodes)
       do m = 1, nmodes
          call rad_cnst_get_info(0, m, nspec=nspec)
          pmam_ncnst = pmam_ncnst + 1 + nspec
       end do

       allocate(pmam_cnst_idx(pmam_ncnst))

       ! Get the constituent indicies
       im = 1
       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, pmam_cnst_idx(im))
          im = im + 1
          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, pmam_cnst_idx(im))
             im = im + 1
          end do
       end do
    end if

    ! ---------------------------------------------------------------------------------------- !
    ! Initialize molecular diffusivity module                                                  !
    ! Note that computing molecular diffusivities is a trivial expense, but constituent        !
    ! diffusivities depend on their molecular weights. Decomposing the diffusion matric        !
    ! for each constituent is a needless expense unless the diffusivity is significant.        !
    ! ---------------------------------------------------------------------------------------- !

    !----------------------------------------------------------------------------------------
    ! Initialize molecular diffusion and get top and bottom molecular diffusion limits
    !----------------------------------------------------------------------------------------

    if( do_molec_diff ) then
       call init_molec_diff( r8, pcnst, rair, mwdry, avogad, gravit, &
            cpair, boltz, errstring)

       call handle_errmsg(errstring, subname="init_molec_diff")

       call addfld( 'TTPXMLC', 'K/S', 1, 'A', 'Top interf. temp. flux: molec. viscosity', phys_decomp )
       call add_default ( 'TTPXMLC', 1, ' ' )
       if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =', ntop_molec, 'NBOT_MOLEC =', nbot_molec
    end if

    ! ---------------------------------- !    
    ! Initialize eddy diffusivity module !
    ! ---------------------------------- !
 
    ! ntop_eddy must be 1 or <= nbot_molec
    ! Currently, it is always 1 except for WACCM-X.
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       ntop_eddy  = press_lim_idx(ntop_eddy_pres, top=.true.)
    else
       ntop_eddy = 1
    end if
    nbot_eddy  = pver

    if (masterproc) write(iulog, fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy

    select case ( eddy_scheme )
    case ( 'diag_TKE' ) 
        if( masterproc ) write(iulog,*) &
             'vertical_diffusion_init: eddy_diffusivity scheme: UW Moist Turbulence Scheme by Bretherton and Park'
        ! Check compatibility of eddy and shallow scheme
        if( shallow_scheme .ne. 'UW' ) then
            write(iulog,*) 'ERROR: shallow convection scheme ', shallow_scheme,' is incompatible with eddy scheme ', eddy_scheme
            call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
        endif
        call init_eddy_diff( r8, pver, gravit, cpair, rair, zvir, latvap, latice, &
                             ntop_eddy, nbot_eddy, karman, eddy_lbulk_max, eddy_leng_max, &
                             eddy_max_bot_pressure )
        if( masterproc ) write(iulog,*) 'vertical_diffusion: nturb, ntop_eddy, nbot_eddy ', nturb, ntop_eddy, nbot_eddy
    case ( 'HB', 'HBR')
        if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'
        call init_hb_diff(gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, &
                          karman, eddy_scheme)
        call addfld('HB_ri', 'no',      pver,  'A',  'Richardson Number (HB Scheme), I',  phys_decomp )
    end select
    
    ! The vertical diffusion solver must operate 
    ! over the full range of molecular and eddy diffusion

    ntop = min(ntop_molec,ntop_eddy)
    nbot = max(nbot_molec,nbot_eddy)
    
    ! ------------------------------------------- !
    ! Initialize turbulent mountain stress module !
    ! ------------------------------------------- !

    if( do_tms ) then
       call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring )

       call handle_errmsg(errstring, subname="init_tms")

       call addfld( 'TAUTMSX' ,'N/m2  ',  1,  'A',  'Zonal      turbulent mountain surface stress',  phys_decomp )
       call addfld( 'TAUTMSY' ,'N/m2  ',  1,  'A',  'Meridional turbulent mountain surface stress',  phys_decomp )

       if (history_amwg) then
          call add_default('TAUTMSX ', 1, ' ')
          call add_default('TAUTMSY ', 1, ' ')
       end if

       if (masterproc) then
          write(iulog,*)'Using turbulent mountain stress module'
          write(iulog,*)'  tms_orocnst = ',tms_orocnst
          write(iulog,*)'  tms_z0fac = ',tms_z0fac
       end if
    endif
    
    ! ---------------------------------- !
    ! Initialize diffusion solver module !
    ! ---------------------------------- !

    call init_vdiff( r8, iulog, rair, gravit, do_iss, errstring )
    call handle_errmsg(errstring, subname="init_vdiff")

    ! Use fieldlist_wet to select the fields which will be diffused using moist mixing ratios ( all by default )
    ! Use fieldlist_dry to select the fields which will be diffused using dry   mixing ratios.

    fieldlist_wet = new_fieldlist_vdiff( pcnst)
    fieldlist_dry = new_fieldlist_vdiff( pcnst)
    fieldlist_molec = new_fieldlist_vdiff( pcnst)

    if( vdiff_select( fieldlist_wet, 'u' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'u' ) )
    if( vdiff_select( fieldlist_wet, 'v' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'v' ) )
    if( vdiff_select( fieldlist_wet, 's' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 's' ) )

    constit_loop: do k = 1, pcnst

       if (prog_modal_aero) then
     ! Do not diffuse droplet number - treated in dropmixnuc
          if (k == ixnumliq) cycle constit_loop
          ! Don't diffuse modal aerosol - treated in dropmixnuc
          do m = 1, pmam_ncnst
             if (k == pmam_cnst_idx(m)) cycle constit_loop
       enddo
       end if

       if( cnst_get_type_byind(k) .eq. 'wet' ) then
          if( vdiff_select( fieldlist_wet, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'q', k ) )
       else
          if( vdiff_select( fieldlist_dry, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_dry, 'q', k ) )
       endif
    
       ! ----------------------------------------------- !
       ! Select constituents for molecular diffusion     !
       ! ----------------------------------------------- !
       if ( cnst_get_molec_byind(k) .eq. 'minor' ) then
         if( vdiff_select(fieldlist_molec,'q',k) .ne. '' ) call endrun( vdiff_select( fieldlist_molec,'q',k ) )
       endif

    end do constit_loop
    
    ! ------------------------ !
    ! Diagnostic output fields !
    ! ------------------------ !

    do k = 1, pcnst
       vdiffnam(k) = 'VD'//cnst_name(k)
       if( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
       call addfld( vdiffnam(k), 'kg/kg/s ', pver, 'A', 'Vertical diffusion of '//cnst_name(k), phys_decomp )
    end do

    call addfld( 'TKE'         , 'm2/s2'  , pverp  , 'A', 'Turbulent Kinetic Energy'                          , phys_decomp )
    call addfld( 'PBLH'        , 'm'      , 1      , 'A', 'PBL height'                                        , phys_decomp )
    call addfld( 'TPERT'       , 'K'      , 1      , 'A', 'Perturbation temperature (eddies in PBL)'          , phys_decomp )
    call addfld( 'QPERT'       , 'kg/kg'  , 1      , 'A', 'Perturbation specific humidity (eddies in PBL)'    , phys_decomp )
    call addfld( 'USTAR'       , 'm/s'    , 1      , 'A', 'Surface friction velocity'                         , phys_decomp )
    call addfld( 'KVH'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion diffusivities (heat/moisture)'  , phys_decomp )
    call addfld( 'KVM'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion diffusivities (momentum)'       , phys_decomp )
    call addfld( 'KVT'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion kinematic molecular conductivity', phys_decomp )
    call addfld( 'CGS'         , 's/m2'   , pverp  , 'A', 'Counter-gradient coeff on surface kinematic fluxes', phys_decomp )
    call addfld( 'DTVKE'       , 'K/s'    , pver   , 'A', 'dT/dt vertical diffusion KE dissipation'           , phys_decomp )
    call addfld( 'DTV'         , 'K/s'    , pver   , 'A', 'T vertical diffusion'                              , phys_decomp )
    call addfld( 'DUV'         , 'm/s2'   , pver   , 'A', 'U vertical diffusion'                              , phys_decomp )
    call addfld( 'DVV'         , 'm/s2'   , pver   , 'A', 'V vertical diffusion'                              , phys_decomp )
    call addfld( 'QT'          , 'kg/kg'  , pver   , 'A', 'Total water mixing ratio'                          , phys_decomp )
    call addfld( 'SL'          , 'J/kg'   , pver   , 'A', 'Liquid water static energy'                        , phys_decomp )
    call addfld( 'SLV'         , 'J/kg'   , pver   , 'A', 'Liq wat virtual static energy'                     , phys_decomp )
    call addfld( 'SLFLX'       , 'W/m2'   , pverp  , 'A', 'Liquid static energy flux'                         , phys_decomp ) 
    call addfld( 'QTFLX'       , 'W/m2'   , pverp  , 'A', 'Total water flux'                                  , phys_decomp ) 
    call addfld( 'UFLX'        , 'W/m2'   , pverp  , 'A', 'Zonal momentum flux'                               , phys_decomp ) 
    call addfld( 'VFLX'        , 'W/m2'   , pverp  , 'A', 'Meridional momentm flux'                           , phys_decomp ) 
    call addfld( 'WGUSTD'      , 'm/s'    , 1      , 'A', 'wind gusts from turbulence'                        , phys_decomp )

    ! ---------------------------------------------------------------------------- !
    ! Below ( with '_PBL') are for detailed analysis of UW Moist Turbulence Scheme !
    ! ---------------------------------------------------------------------------- !

    call addfld( 'qt_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qt_prePBL'                                         , phys_decomp )
    call addfld( 'sl_pre_PBL  ', 'J/kg'   , pver   , 'A', 'sl_prePBL'                                         , phys_decomp )
    call addfld( 'slv_pre_PBL ', 'J/kg'   , pver   , 'A', 'slv_prePBL'                                        , phys_decomp )
    call addfld( 'u_pre_PBL   ', 'm/s'    , pver   , 'A', 'u_prePBL'                                          , phys_decomp )
    call addfld( 'v_pre_PBL   ', 'm/s'    , pver   , 'A', 'v_prePBL'                                          , phys_decomp )
    call addfld( 'qv_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qv_prePBL'                                         , phys_decomp )
    call addfld( 'ql_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'ql_prePBL'                                         , phys_decomp )
    call addfld( 'qi_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qi_prePBL'                                         , phys_decomp )
    call addfld( 't_pre_PBL   ', 'K'      , pver   , 'A', 't_prePBL'                                          , phys_decomp )
    call addfld( 'rh_pre_PBL  ', '%'      , pver   , 'A', 'rh_prePBL'                                         , phys_decomp )

    call addfld( 'qt_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qt_afterPBL'                                       , phys_decomp )
    call addfld( 'sl_aft_PBL  ', 'J/kg'   , pver   , 'A', 'sl_afterPBL'                                       , phys_decomp )
    call addfld( 'slv_aft_PBL ', 'J/kg'   , pver   , 'A', 'slv_afterPBL'                                      , phys_decomp )
    call addfld( 'u_aft_PBL   ', 'm/s'    , pver   , 'A', 'u_afterPBL'                                        , phys_decomp )
    call addfld( 'v_aft_PBL   ', 'm/s'    , pver   , 'A', 'v_afterPBL'                                        , phys_decomp )
    call addfld( 'qv_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qv_afterPBL'                                       , phys_decomp )
    call addfld( 'ql_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'ql_afterPBL'                                       , phys_decomp )
    call addfld( 'qi_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qi_afterPBL'                                       , phys_decomp )
    call addfld( 't_aft_PBL   ', 'K'      , pver   , 'A', 't_afterPBL'                                        , phys_decomp )
    call addfld( 'rh_aft_PBL  ', '%'      , pver   , 'A', 'rh_afterPBL'                                       , phys_decomp )

    call addfld( 'slflx_PBL   ', 'J/m2/s' , pverp  , 'A', 'sl flux by PBL'                                    , phys_decomp ) 
    call addfld( 'qtflx_PBL   ', 'kg/m2/s', pverp  , 'A', 'qt flux by PBL'                                    , phys_decomp ) 
    call addfld( 'uflx_PBL    ', 'kg/m/s2', pverp  , 'A', 'u flux by PBL'                                     , phys_decomp ) 
    call addfld( 'vflx_PBL    ', 'kg/m/s2', pverp  , 'A', 'v flux by PBL'                                     , phys_decomp ) 

    call addfld( 'slflx_cg_PBL', 'J/m2/s' , pverp  , 'A', 'sl_cg flux by PBL'                                 , phys_decomp ) 
    call addfld( 'qtflx_cg_PBL', 'kg/m2/s', pverp  , 'A', 'qt_cg flux by PBL'                                 , phys_decomp ) 
    call addfld( 'uflx_cg_PBL ', 'kg/m/s2', pverp  , 'A', 'u_cg flux by PBL'                                  , phys_decomp ) 
    call addfld( 'vflx_cg_PBL ', 'kg/m/s2', pverp  , 'A', 'v_cg flux by PBL'                                  , phys_decomp ) 

    call addfld( 'qtten_PBL   ', 'kg/kg/s', pver   , 'A', 'qt tendency by PBL'                                , phys_decomp )
    call addfld( 'slten_PBL   ', 'J/kg/s' , pver   , 'A', 'sl tendency by PBL'                                , phys_decomp )
    call addfld( 'uten_PBL    ', 'm/s2'   , pver   , 'A', 'u tendency by PBL'                                 , phys_decomp )
    call addfld( 'vten_PBL    ', 'm/s2'   , pver   , 'A', 'v tendency by PBL'                                 , phys_decomp )
    call addfld( 'qvten_PBL   ', 'kg/kg/s', pver   , 'A', 'qv tendency by PBL'                                , phys_decomp )
    call addfld( 'qlten_PBL   ', 'kg/kg/s', pver   , 'A', 'ql tendency by PBL'                                , phys_decomp )
    call addfld( 'qiten_PBL   ', 'kg/kg/s', pver   , 'A', 'qi tendency by PBL'                                , phys_decomp )
    call addfld( 'tten_PBL    ', 'K/s'    , pver   , 'A', 'T tendency by PBL'                                 , phys_decomp )
    call addfld( 'rhten_PBL   ', '%/s'    , pver   , 'A', 'RH tendency by PBL'                                , phys_decomp )

    call addfld ('ustar',     ' ',1, 'A',' ',phys_decomp)
    call addfld ('obklen',    ' ',1, 'A',' ',phys_decomp)

    if( eddy_scheme .eq. 'diag_TKE' ) then    
       call addfld( 'BPROD   ',  'M2/S3   ',pverp,   'A', 'Buoyancy Production'                               ,phys_decomp)
       call addfld( 'SFI     ',  'FRACTION',pverp,   'A', 'Interface-layer sat frac'                          ,phys_decomp)    
       call addfld( 'SPROD   ',  'M2/S3   ',pverp,   'A', 'Shear Production'                                  ,phys_decomp)   
    endif
 
    ! ----------------------------
    ! determine default variables
    ! ----------------------------
 
    call phys_getopts( history_amwg_out = history_amwg, &
                       history_eddy_out = history_eddy, &
                       history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_waccm_out = history_waccm)

    if (history_amwg) then
       call add_default(  vdiffnam(1), 1, ' ' )
       call add_default( 'DTV'       , 1, ' ' )  
       call add_default( 'PBLH'      , 1, ' ' )
       if( eddy_scheme .eq. 'diag_TKE' ) then    
          call add_default( 'WGUSTD  ', 1, ' ' )
       endif
    endif
 
    if (history_eddy) then
       if( eddy_scheme .eq. 'diag_TKE' ) then    
          call add_default( 'UFLX    ', 1, ' ' )
          call add_default( 'VFLX    ', 1, ' ' )
       endif
    endif

    if( history_budget ) then
        call add_default( vdiffnam(ixcldliq), history_budget_histfile_num, ' ' )
        call add_default( vdiffnam(ixcldice), history_budget_histfile_num, ' ' )
        if( history_budget_histfile_num > 1 ) then
           call add_default(  vdiffnam(1), history_budget_histfile_num, ' ' )
           call add_default( 'DTV'       , history_budget_histfile_num, ' ' )
        end if
    end if

    if ( history_waccm ) then
       call add_default( 'DUV'     , 1, ' ' )
       call add_default( 'DVV'     , 1, ' ' )
    end if
     ! ----------------------------
   

     qrl_idx   = pbuf_get_index('QRL')
     wsedl_idx = pbuf_get_index('WSEDL')


     ! Initialization of some pbuf fields
     if (is_first_step()) then
        ! Initialization of pbuf fields tke, kvh, kvm are done in phys_inidat
        call pbuf_set_field(pbuf2d, turbtype_idx, 0    )
        call pbuf_set_field(pbuf2d, smaw_idx,     0._r8)
        call pbuf_set_field(pbuf2d, tauresx_idx,  0._r8)
        call pbuf_set_field(pbuf2d, tauresy_idx,  0._r8)
     end if

  end subroutine vertical_diffusion_init

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ts_init( pbuf2d, state )

    !-------------------------------------------------------------- !
    ! Timestep dependent setting,                                   !
    ! At present only invokes upper bc code for molecular diffusion !
    !-------------------------------------------------------------- !
    use molec_diff    , only : init_timestep_molec_diff
    use physics_types , only : physics_state
    use ppgrid        , only : begchunk, endchunk
    
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

 
    if (do_molec_diff) call init_timestep_molec_diff(pbuf2d, state )

  end subroutine vertical_diffusion_ts_init

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_tend( &
                                      ztodt    , state    ,                  &
                                      taux     , tauy     , shflx    , cflx, &
                                      ustar    , obklen   , ptend    , &
                                      cldn     , ocnfrac  , landfrac , sgh      , pbuf) 
    !---------------------------------------------------- !
    ! This is an interface routine for vertical diffusion !
    !---------------------------------------------------- !
    use physics_buffer,     only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
    use physics_types,      only : physics_state, physics_ptend, physics_ptend_init
    use cam_history,        only : outfld
    
    use trb_mtn_stress,     only : compute_tms
    use eddy_diff,          only : compute_eddy_diff
    use hb_diff,            only : compute_hb_diff
    use wv_saturation,      only : qsat
    use molec_diff,         only : compute_molec_diff, vd_lu_qdecomp
    use constituents,       only : qmincg, qmin
    use diffusion_solver,   only : compute_vdiff, any, operator(.not.)
    use physconst,          only : cpairv, rairv !Needed for calculation of upward H flux
    use time_manager,       only : get_nstep
    use constituents,       only : cnst_get_type_byind, cnst_name, cnst_fixed_ubc, cnst_fixed_ubflx
    use physconst,          only : pi
    use pbl_utils,          only : virtem, calc_obklen

    ! --------------- !
    ! Input Arguments !
    ! --------------- !

    type(physics_state), intent(in)    :: state                     ! Physics state variables

    real(r8),            intent(in)    :: taux(pcols)               ! x surface stress  [ N/m2 ]
    real(r8),            intent(in)    :: tauy(pcols)               ! y surface stress  [ N/m2 ]
    real(r8),            intent(in)    :: shflx(pcols)              ! Surface sensible heat flux  [ w/m2 ]
    real(r8),            intent(in)    :: cflx(pcols,pcnst)         ! Surface constituent flux [ kg/m2/s ]
    real(r8),            intent(in)    :: ztodt                     ! 2 delta-t [ s ]
    real(r8),            intent(in)    :: cldn(pcols,pver)          ! New stratus fraction [ fraction ]
    real(r8),            intent(in)    :: ocnfrac(pcols)            ! Ocean fraction
    real(r8),            intent(in)    :: landfrac(pcols)           ! Land fraction
    real(r8),            intent(in)    :: sgh(pcols)                ! Standard deviation of orography [ unit ? ]

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !
    
    type(physics_ptend), intent(out) :: ptend                       ! Individual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    real(r8),            intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    real(r8),            intent(out)   :: obklen(pcols)             ! Obukhov length [ m ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    logical  :: kvinit                                              ! Tell compute_eddy_diff/ caleddy
                                                                    ! to initialize kvh, kvm (uses kvf)
    character(128) :: errstring                                     ! Error status for compute_vdiff
    real(r8), pointer, dimension(:,:) :: qrl                        ! LW radiative cooling rate
    real(r8), pointer, dimension(:,:) :: wsedl                      ! Sedimentation velocity
                                                                    ! of stratiform liquid cloud droplet [ m/s ]

    integer  :: lchnk                                               ! Chunk identifier
    integer  :: ncol                                                ! Number of atmospheric columns
    integer  :: i, k, l, m                                          ! column, level, constituent indices
    integer  :: ierr                                                ! status for allocate/deallocate

    real(r8) :: dtk(pcols,pver)                                     ! T tendency from KE dissipation
    real(r8), pointer   :: tke(:,:)                                 ! Turbulent kinetic energy [ m2/s2 ]
    integer(i4),pointer :: turbtype(:,:)                            ! Turbulent interface types [ no unit ]
    real(r8), pointer   :: smaw(:,:)                                ! Normalized Galperin instability function
                                                                    ! ( 0<= <=4.964 and 1 at neutral )
    real(r8) :: cgs(pcols,pverp)                                    ! Counter-gradient star  [ cg/flux ]
    real(r8) :: cgh(pcols,pverp)                                    ! Counter-gradient term for heat
    real(r8) :: rztodt                                              ! 1./ztodt [ 1/s ]
    real(r8) :: ksrftms(pcols)                                      ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
    real(r8) :: tautmsx(pcols)                                      ! U component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautmsy(pcols)                                      ! V component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautotx(pcols)                                      ! U component of total surface stress [ N/m2 ]
    real(r8) :: tautoty(pcols)                                      ! V component of total surface stress [ N/m2 ]

    real(r8), pointer :: kvh_in(:,:)                                ! kvh from previous timestep [ m2/s ]
    real(r8), pointer :: kvm_in(:,:)                                ! kvm from previous timestep [ m2/s ]
    real(r8), pointer :: kvt(:,:)                                   ! Molecular kinematic conductivity for temperature [  ]
    real(r8) :: kvq(pcols,pverp)                                    ! Eddy diffusivity for constituents [ m2/s ]
    real(r8) :: kvh(pcols,pverp)                                    ! Eddy diffusivity for heat [ m2/s ]
    real(r8) :: kvm(pcols,pverp)                                    ! Eddy diffusivity for momentum [ m2/s ]
    real(r8) :: bprod(pcols,pverp)                                  ! Buoyancy production of tke [ m2/s3 ]
    real(r8) :: sprod(pcols,pverp)                                  ! Shear production of tke [ m2/s3 ]
    real(r8) :: sfi(pcols,pverp)                                    ! Saturation fraction at interfaces [ fraction ]
    real(r8) :: sl(pcols,pver)
    real(r8) :: qt(pcols,pver)
    real(r8) :: slv(pcols,pver)
    real(r8) :: sl_prePBL(pcols,pver)
    real(r8) :: qt_prePBL(pcols,pver)
    real(r8) :: slv_prePBL(pcols,pver)
    real(r8) :: slten(pcols,pver)
    real(r8) :: qtten(pcols,pver)
    real(r8) :: slvten(pcols,pver)
    real(r8) :: slflx(pcols,pverp)
    real(r8) :: qtflx(pcols,pverp)
    real(r8) :: uflx(pcols,pverp)
    real(r8) :: vflx(pcols,pverp)
    real(r8) :: slflx_cg(pcols,pverp)
    real(r8) :: qtflx_cg(pcols,pverp)
    real(r8) :: uflx_cg(pcols,pverp)
    real(r8) :: vflx_cg(pcols,pverp)
    real(r8) :: th(pcols,pver)                                      ! Potential temperature
    real(r8) :: topflx(pcols)                                       ! Molecular heat flux at top interface
    real(r8) :: wpert(pcols)                                        ! Turbulent wind gusts
    real(r8) :: rhoair

    real(r8) :: ri(pcols,pver)                                      ! richardson number (HB output)
    
    ! for obklen calculation outside HB
    real(r8) :: thvs(pcols)                                         ! Virtual potential temperature at surface
    real(r8) :: rrho(pcols)                                         ! Reciprocal of density at surface
    real(r8) :: khfs(pcols)                                         ! sfc kinematic heat flux [mK/s]
    real(r8) :: kqfs(pcols)                                         ! sfc kinematic water vapor flux [m/s]
    real(r8) :: kbfs(pcols)                                         ! sfc kinematic buoyancy flux [m^2/s^3]

    real(r8) :: ftem(pcols,pver)                                    ! Saturation vapor pressure before PBL
    real(r8) :: ftem_prePBL(pcols,pver)                             ! Saturation vapor pressure before PBL
    real(r8) :: ftem_aftPBL(pcols,pver)                             ! Saturation vapor pressure after PBL
    real(r8) :: tem2(pcols,pver)                                    ! Saturation specific humidity and RH
    real(r8) :: t_aftPBL(pcols,pver)                                ! Temperature after PBL diffusion
    real(r8) :: tten(pcols,pver)                                    ! Temperature tendency by PBL diffusion
    real(r8) :: rhten(pcols,pver)                                   ! RH tendency by PBL diffusion
    real(r8) :: qv_aft_PBL(pcols,pver)                              ! qv after PBL diffusion
    real(r8) :: ql_aft_PBL(pcols,pver)                              ! ql after PBL diffusion
    real(r8) :: qi_aft_PBL(pcols,pver)                              ! qi after PBL diffusion
    real(r8) :: s_aft_PBL(pcols,pver)                               ! s after PBL diffusion
    real(r8) :: u_aft_PBL(pcols,pver)                               ! u after PBL diffusion
    real(r8) :: v_aft_PBL(pcols,pver)                               ! v after PBL diffusion
    real(r8) :: qv_pro(pcols,pver) 
    real(r8) :: ql_pro(pcols,pver)
    real(r8) :: qi_pro(pcols,pver)
    real(r8) :: s_pro(pcols,pver)
    real(r8) :: t_pro(pcols,pver)
    real(r8), pointer :: tauresx(:)                                      ! Residual stress to be added in vdiff to correct
    real(r8), pointer :: tauresy(:)                                      ! for turb stress mismatch between sfc and atm accumulated.
    real(r8) :: ipbl(pcols)
    real(r8) :: kpblh(pcols)
    real(r8) :: wstarPBL(pcols)
    real(r8) :: tpertPBL(pcols)
    real(r8) :: qpertPBL(pcols)
    real(r8) :: rairi(pcols,pver+1)                                 ! interface gas constant needed for compute_vdiff

    real(r8), pointer :: tpert(:)
    real(r8), pointer :: qpert(:)
    real(r8), pointer :: pblh(:)

    real(r8) :: tmp1(pcols)                                         ! Temporary storage

    integer  :: nstep
    real(r8) :: sum1, sum2, sum3, pdelx 
    real(r8) :: sflx

    ! Copy state so we can pass to intent(inout) routines that return
    ! new state instead of a tendency.
    real(r8) :: s_tmp(pcols,pver)
    real(r8) :: u_tmp(pcols,pver)
    real(r8) :: v_tmp(pcols,pver)
    real(r8) :: q_tmp(pcols,pver,pcnst)

    logical  :: lq(pcnst)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    rztodt = 1._r8 / ztodt
    lchnk  = state%lchnk
    ncol   = state%ncol

    call pbuf_get_field(pbuf, tauresx_idx,  tauresx)
    call pbuf_get_field(pbuf, tauresy_idx,  tauresy)
    call pbuf_get_field(pbuf, tpert_idx,    tpert)
    call pbuf_get_field(pbuf, qpert_idx,    qpert)
    call pbuf_get_field(pbuf, pblh_idx,     pblh)
    call pbuf_get_field(pbuf, turbtype_idx, turbtype)

    ! ---------------------------------------- !
    ! Computation of turbulent mountain stress !
    ! ---------------------------------------- !
   
    ! Consistent with the computation of 'normal' drag coefficient, we are using 
    ! the raw input (u,v) to compute 'ksrftms', not the provisionally-marched 'u,v' 
    ! within the iteration loop of the PBL scheme. 

    if( do_tms ) then
        call compute_tms( pcols      , pver     , ncol    ,              &
                          state%u    , state%v  , state%t , state%pmid , & 
                          state%exner, state%zm , sgh     , ksrftms    , & 
                          tautmsx    , tautmsy  , landfrac )
      ! Here, both 'taux, tautmsx' are explicit surface stresses.        
      ! Note that this 'tautotx, tautoty' are different from the total stress
      ! that has been actually added into the atmosphere. This is because both
      ! taux and tautmsx are fully implicitly treated within compute_vdiff.
      ! However, 'tautotx, tautoty' are not used in the actual numerical
      ! computation in this module.   
        tautotx(:ncol) = taux(:ncol) + tautmsx(:ncol)
        tautoty(:ncol) = tauy(:ncol) + tautmsy(:ncol)
    else
        ksrftms(:ncol) = 0._r8
        tautotx(:ncol) = taux(:ncol)
        tautoty(:ncol) = tauy(:ncol)
    endif

    !----------------------------------------------------------------------- !
    !   Computation of eddy diffusivities - Select appropriate PBL scheme    !
    !----------------------------------------------------------------------- !
    call pbuf_get_field(pbuf, kvm_idx,  kvm_in)
    call pbuf_get_field(pbuf, kvh_idx,  kvh_in)
    call pbuf_get_field(pbuf, smaw_idx, smaw)
    call pbuf_get_field(pbuf, tke_idx,  tke)


    select case (eddy_scheme)
    case ( 'diag_TKE' ) 

       ! ---------------------------------------------------------------- !
       ! At first time step, have eddy_diff.F90:caleddy() use kvh=kvm=kvf !
       ! This has to be done in compute_eddy_diff after kvf is calculated !
       ! ---------------------------------------------------------------- !

       if( is_first_step() ) then
           kvinit = .true.
       else
           kvinit = .false.
       endif

       ! ---------------------------------------------- !
       ! Get LW radiative heating out of physics buffer !
       ! ---------------------------------------------- !
       call pbuf_get_field(pbuf, qrl_idx,   qrl   )

       call pbuf_get_field(pbuf, wsedl_idx, wsedl )

       ! Retrieve eddy diffusivities for heat and momentum from physics buffer
       ! from last timestep ( if first timestep, has been initialized by inidat.F90 )

       call compute_eddy_diff( lchnk    ,                                                                    &
                               pcols    , pver        , ncol       , state%t    , state%q(:,:,1) , ztodt   , &
                               state%q(:,:,ixcldliq)  , state%q(:,:,ixcldice)   ,                            &
                               state%s  , state%rpdel , cldn       , qrl        , wsedl          ,           &
                               state%zm , state%zi    , state%pmid , state%pint , state%u        , state%v , &
                               taux     , tauy        , shflx      , cflx(:,1)  , wstarent       , nturb   , &
                               rrho     , ustar       , pblh       , kvm_in     , kvh_in         , kvm     , &
                               kvh      , kvq         , cgh        ,                                         &
                               cgs      , tpert       , qpert      , wpert      , tke            , bprod   , &
                               sprod    , sfi         , kvinit     ,                                         &
                               tauresx  , tauresy     , ksrftms    ,                                         &
                               ipbl(:)  , kpblh(:)    , wstarPBL(:), turbtype   , smaw )

       ! The diag_TKE scheme does not calculate the Monin-Obukhov length, which is used in dry deposition calculations.
       ! Use the routines from pbl_utils to accomplish this. Assumes ustar and rrho have been set.
       th(:ncol,pver) = state%t(:ncol,pver) * state%exner(:ncol,pver)
       thvs(:ncol) = virtem(th(:ncol,pver),state%q(:ncol,pver,1))
       call calc_obklen(th(:ncol,pver), thvs(:ncol), cflx(:ncol,1), shflx(:ncol), rrho(:ncol), ustar(:ncol), &
                        khfs(:ncol),    kqfs(:ncol), kbfs(:ncol),   obklen(:ncol))

       ! ----------------------------------------------- !       
       ! Store TKE in pbuf for use by shallow convection !
       ! ----------------------------------------------- !   

       tpertPBL(:ncol) = tpert(:ncol)
       qpertPBL(:ncol) = qpert(:ncol)


       ! The diffusivities from diag_TKE can be much larger than from HB in the free
       ! troposphere and upper atmosphere. These seem to be larger than observations,
       ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
       ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
       ! or the upper atmosphere.
       !
       ! NOTE: Further investigation should be done as to why the diffusivities are
       ! larger in diag_TKE.
       if ((kv_freetrop_scale /= 1._r8) .or. ((kv_top_scale /= 1._r8) .and. (kv_top_pressure > 0._r8))) then
         do i = 1, ncol
           do k = 1, pverp
           
             ! Outside of the boundary layer?
             if (state%zi(i,k) > pblh(i)) then

               ! In the upper atmosphere?
               if (state%pint(i,k) <= kv_top_pressure) then
                 kvh(i,k) = kvh(i,k) * kv_top_scale
                 kvm(i,k) = kvm(i,k) * kv_top_scale
                 kvq(i,k) = kvq(i,k) * kv_top_scale
               else
                 kvh(i,k) = kvh(i,k) * kv_freetrop_scale
                 kvm(i,k) = kvm(i,k) * kv_freetrop_scale
                 kvq(i,k) = kvq(i,k) * kv_freetrop_scale
               end if
             else
               exit
             end if
           end do
         end do
       end if

       ! Write out fields that are only used by this scheme

       call outfld( 'BPROD   ', bprod(1,1), pcols, lchnk )
       call outfld( 'SPROD   ', sprod(1,1), pcols, lchnk )
       call outfld( 'SFI     ', sfi,        pcols, lchnk )

    case ( 'HB', 'HBR' )

       ! Modification : We may need to use 'taux' instead of 'tautotx' here, for
       !                consistency with the previous HB scheme.

       th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

       call compute_hb_diff( lchnk     , ncol     ,                                &
                             th        , state%t  , state%q , state%zm , state%zi, &
                             state%pmid, state%u  , state%v , tautotx  , tautoty , &
                             shflx     , cflx(:,1), obklen  , ustar    , pblh    , &
                             kvm       , kvh      , kvq     , cgh      , cgs     , &
                             tpert     , qpert    , cldn    , ocnfrac  , tke     , &
                             ri        , &
                             eddy_scheme )
    
       call outfld( 'HB_ri',          ri,         pcols,   lchnk )

       wpert = 0._r8  

    end select

    call outfld( 'ustar',   ustar(:), pcols, lchnk )
    call outfld( 'obklen', obklen(:), pcols, lchnk )

    ! kvh (in pbuf) is used by other physics parameterizations, and as an initial guess in compute_eddy_diff
    ! on the next timestep.  It is not updated by the compute_vdiff call below.
    call pbuf_set_field(pbuf, kvh_idx, kvh)

    ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
    ! The contributions for molecular diffusion made to kvm by the call to compute_vdiff below 
    ! are not included in the pbuf as these are not needed in the initial guess by compute_eddy_diff.
    call pbuf_set_field(pbuf, kvm_idx, kvm)

    call outfld( 'WGUSTD' , wpert, pcols, lchnk )

    !------------------------------------ ! 
    !    Application of diffusivities     !
    !------------------------------------ !

    ! Set arrays from input state.
    q_tmp(:ncol,:,:) = state%q(:ncol,:,:)
    s_tmp(:ncol,:) = state%s(:ncol,:)
    u_tmp(:ncol,:) = state%u(:ncol,:)
    v_tmp(:ncol,:) = state%v(:ncol,:)

    !------------------------------------------------------ !
    ! Write profile output before applying diffusion scheme !
    !------------------------------------------------------ !

    sl_prePBL(:ncol,:pver)  = s_tmp(:ncol,:) -   latvap * q_tmp(:ncol,:,ixcldliq) &
                                             - ( latvap + latice) * q_tmp(:ncol,:,ixcldice)
    qt_prePBL(:ncol,:pver)  = q_tmp(:ncol,:,1) + q_tmp(:ncol,:,ixcldliq) &
                                               + q_tmp(:ncol,:,ixcldice)
    slv_prePBL(:ncol,:pver) = sl_prePBL(:ncol,:pver) * ( 1._r8 + zvir*qt_prePBL(:ncol,:pver) ) 

    call qsat(state%t(:ncol,:), state%pmid(:ncol,:), &
         tem2(:ncol,:), ftem(:ncol,:))
    ftem_prePBL(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8

    call outfld( 'qt_pre_PBL   ', qt_prePBL,                 pcols, lchnk )
    call outfld( 'sl_pre_PBL   ', sl_prePBL,                 pcols, lchnk )
    call outfld( 'slv_pre_PBL  ', slv_prePBL,                pcols, lchnk )
    call outfld( 'u_pre_PBL    ', state%u,                   pcols, lchnk )
    call outfld( 'v_pre_PBL    ', state%v,                   pcols, lchnk )
    call outfld( 'qv_pre_PBL   ', state%q(:ncol,:,1),        pcols, lchnk )
    call outfld( 'ql_pre_PBL   ', state%q(:ncol,:,ixcldliq), pcols, lchnk )
    call outfld( 'qi_pre_PBL   ', state%q(:ncol,:,ixcldice), pcols, lchnk )
    call outfld( 't_pre_PBL    ', state%t,                   pcols, lchnk )
    call outfld( 'rh_pre_PBL   ', ftem_prePBL,               pcols, lchnk )

    ! --------------------------------------------------------------------------------- !
    ! Call the diffusivity solver and solve diffusion equation                          !
    ! The final two arguments are optional function references to                       !
    ! constituent-independent and constituent-dependent moleculuar diffusivity routines !
    ! --------------------------------------------------------------------------------- !

    !------------------------------------------------------------------------ 
    !  Check to see if constituent dependent gas constant needed (WACCM-X) 
    !------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
      rairi(:ncol,1) = rairv(:ncol,1,lchnk)
      do k = 2, pver
        do i = 1, ncol
          rairi(i,k) = 0.5_r8 * (rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
        end do
      end do
    else
      rairi(:ncol,:pver+1) = rair 
    endif

    ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
    !                separately print out as diagnostic output, because these are different from
    !                the explicit 'tautotx, tautoty' computed above. 
    ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

    call pbuf_get_field(pbuf, kvt_idx, kvt)

    if( any(fieldlist_wet) ) then

        call compute_vdiff( state%lchnk   ,                                                                     &
                            pcols         , pver               , pcnst        , ncol          , state%pmid    , &
                            state%pint    , state%rpdel        , state%t      , ztodt         , taux          , &
                            tauy          , shflx              , cflx         , ntop          , nbot          , &
                            kvh           , kvm                , kvq          , cgs           , cgh           , &
                            state%zi      , ksrftms            , qmincg       , fieldlist_wet , fieldlist_molec,&
                            u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
                            tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
                            tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), rairi, &
                            do_molec_diff , compute_molec_diff , vd_lu_qdecomp, kvt )

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_wet call from vertical_diffusion.")

    end if
 
    if( any( fieldlist_dry ) ) then

        if( do_molec_diff ) then
            errstring = "Design flaw: dry vdiff not currently supported with molecular diffusion"
            call endrun(errstring)
        end if

        call compute_vdiff( state%lchnk   ,                                                                     &
                            pcols         , pver               , pcnst        , ncol          , state%pmiddry , &
                            state%pintdry , state%rpdeldry     , state%t      , ztodt         , taux          , &       
                            tauy          , shflx              , cflx         , ntop          , nbot          , &       
                            kvh           , kvm                , kvq          , cgs           , cgh           , &   
                            state%zi      , ksrftms            , qmincg       , fieldlist_dry , fieldlist_molec,&
                            u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
                            tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
                            tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), rairi, &
                            do_molec_diff , compute_molec_diff , vd_lu_qdecomp )

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_dry call from vertical_diffusion.")

    end if

    if (prog_modal_aero) then

       ! Modal aerosol species not diffused, so just add the explicit surface fluxes to the
       ! lowest layer

    tmp1(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
       do m = 1, pmam_ncnst
          l = pmam_cnst_idx(m)
          q_tmp(:ncol,pver,l) = q_tmp(:ncol,pver,l) + tmp1(:ncol) * cflx(:ncol,l)
       enddo
    end if

    ! -------------------------------------------------------- !
    ! Diagnostics and output writing after applying PBL scheme !
    ! -------------------------------------------------------- !

    sl(:ncol,:pver)  = s_tmp(:ncol,:) -   latvap           * q_tmp(:ncol,:,ixcldliq) &
                                      - ( latvap + latice) * q_tmp(:ncol,:,ixcldice)
    qt(:ncol,:pver)  = q_tmp(:ncol,:,1) + q_tmp(:ncol,:,ixcldliq) &
                                    + q_tmp(:ncol,:,ixcldice)
    slv(:ncol,:pver) = sl(:ncol,:pver) * ( 1._r8 + zvir*qt(:ncol,:pver) ) 

    slflx(:ncol,1) = 0._r8
    qtflx(:ncol,1) = 0._r8
    uflx(:ncol,1)  = 0._r8
    vflx(:ncol,1)  = 0._r8

    slflx_cg(:ncol,1) = 0._r8
    qtflx_cg(:ncol,1) = 0._r8
    uflx_cg(:ncol,1)  = 0._r8
    vflx_cg(:ncol,1)  = 0._r8

    do k = 2, pver
       do i = 1, ncol
          rhoair     = state%pint(i,k) / ( rair * ( ( 0.5_r8*(slv(i,k)+slv(i,k-1)) - gravit*state%zi(i,k))/cpair ) )
          slflx(i,k) = kvh(i,k) * &
                          ( - rhoair*(sl(i,k-1)-sl(i,k))/(state%zm(i,k-1)-state%zm(i,k)) &
                            + cgh(i,k) ) 
          qtflx(i,k) = kvh(i,k) * &
                          ( - rhoair*(qt(i,k-1)-qt(i,k))/(state%zm(i,k-1)-state%zm(i,k)) &
                            + rhoair*(cflx(i,1)+cflx(i,ixcldliq)+cflx(i,ixcldice))*cgs(i,k) )
          uflx(i,k)  = kvm(i,k) * &
                          ( - rhoair*(u_tmp(i,k-1)-u_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
          vflx(i,k)  = kvm(i,k) * &
                          ( - rhoair*(v_tmp(i,k-1)-v_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
          slflx_cg(i,k) = kvh(i,k) * cgh(i,k)
          qtflx_cg(i,k) = kvh(i,k) * rhoair * ( cflx(i,1) + cflx(i,ixcldliq) + cflx(i,ixcldice) ) * cgs(i,k)
          uflx_cg(i,k)  = 0._r8
          vflx_cg(i,k)  = 0._r8
       end do
    end do

  ! Modification : I should check whether slflx(:ncol,pverp) is correctly computed.
  !                Note also that 'tautotx' is explicit total stress, different from
  !                the ones that have been actually added into the atmosphere.

    slflx(:ncol,pverp) = shflx(:ncol)
    qtflx(:ncol,pverp) = cflx(:ncol,1)
    uflx(:ncol,pverp)  = tautotx(:ncol)
    vflx(:ncol,pverp)  = tautoty(:ncol)

    slflx_cg(:ncol,pverp) = 0._r8
    qtflx_cg(:ncol,pverp) = 0._r8
    uflx_cg(:ncol,pverp)  = 0._r8
    vflx_cg(:ncol,pverp)  = 0._r8

    ! --------------------------------------------------------------- !
    ! Convert the new profiles into vertical diffusion tendencies.    !
    ! Convert KE dissipative heat change into "temperature" tendency. !
    ! --------------------------------------------------------------- !

    ! All variables are modified by vertical diffusion

    lq(:) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, "vertical diffusion", &
         ls=.true., lu=.true., lv=.true., lq=lq)

    ptend%s(:ncol,:)       = ( s_tmp(:ncol,:) - state%s(:ncol,:) ) * rztodt
    ptend%u(:ncol,:)       = ( u_tmp(:ncol,:) - state%u(:ncol,:) ) * rztodt
    ptend%v(:ncol,:)       = ( v_tmp(:ncol,:) - state%v(:ncol,:) ) * rztodt
    ptend%q(:ncol,:pver,:) = ( q_tmp(:ncol,:pver,:) - state%q(:ncol,:pver,:) ) * rztodt
    slten(:ncol,:)         = ( sl(:ncol,:) - sl_prePBL(:ncol,:) ) * rztodt
    qtten(:ncol,:)         = ( qt(:ncol,:) - qt_prePBL(:ncol,:) ) * rztodt

    ! ----------------------------------------------------------- !
    ! In order to perform 'pseudo-conservative varible diffusion' !
    ! perform the following two stages:                           !
    !                                                             !
    ! I.  Re-set (1) 'qvten' by 'qtten', and 'qlten = qiten = 0'  !
    !            (2) 'sten'  by 'slten', and                      !
    !            (3) 'qlten = qiten = 0'                          !
    !                                                             !
    ! II. Apply 'positive_moisture'                               !
    !                                                             !
    ! ----------------------------------------------------------- !

    if( eddy_scheme .eq. 'diag_TKE' .and. do_pseudocon_diff ) then

         ptend%q(:ncol,:pver,1) = qtten(:ncol,:pver)
         ptend%s(:ncol,:pver)   = slten(:ncol,:pver)
         ptend%q(:ncol,:pver,ixcldliq) = 0._r8         
         ptend%q(:ncol,:pver,ixcldice) = 0._r8         
         ptend%q(:ncol,:pver,ixnumliq) = 0._r8         
         ptend%q(:ncol,:pver,ixnumice) = 0._r8         

         do i = 1, ncol
            do k = 1, pver
               qv_pro(i,k) = state%q(i,k,1)        + ptend%q(i,k,1)             * ztodt       
               ql_pro(i,k) = state%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)      * ztodt
               qi_pro(i,k) = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)      * ztodt              
               s_pro(i,k)  = state%s(i,k)          + ptend%s(i,k)               * ztodt
               t_pro(i,k)  = state%t(i,k)          + (1._r8/cpair)*ptend%s(i,k) * ztodt
            end do 
         end do 

         call positive_moisture( cpair, latvap, latvap+latice, ncol, pver, ztodt, qmin(1), qmin(ixcldliq), qmin(ixcldice),    &
                                 state%pdel(:ncol,pver:1:-1), qv_pro(:ncol,pver:1:-1), ql_pro(:ncol,pver:1:-1), &
                                 qi_pro(:ncol,pver:1:-1), t_pro(:ncol,pver:1:-1), s_pro(:ncol,pver:1:-1),       &
                                 ptend%q(:ncol,pver:1:-1,1), ptend%q(:ncol,pver:1:-1,ixcldliq),                 &
                                 ptend%q(:ncol,pver:1:-1,ixcldice), ptend%s(:ncol,pver:1:-1) )

    end if

    ! ----------------------------------------------------------------- !
    ! Re-calculate diagnostic output variables after vertical diffusion !
    ! ----------------------------------------------------------------- !
 
    qv_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,1)         + ptend%q(:ncol,:pver,1)        * ztodt
    ql_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,ixcldliq)  + ptend%q(:ncol,:pver,ixcldliq) * ztodt
    qi_aft_PBL(:ncol,:pver)  =   state%q(:ncol,:pver,ixcldice)  + ptend%q(:ncol,:pver,ixcldice) * ztodt
    s_aft_PBL(:ncol,:pver)   =   state%s(:ncol,:pver)           + ptend%s(:ncol,:pver)          * ztodt
    t_aftPBL(:ncol,:pver)    = ( s_aft_PBL(:ncol,:pver) - gravit*state%zm(:ncol,:pver) ) / cpair 

    u_aft_PBL(:ncol,:pver)   =  state%u(:ncol,:pver)          + ptend%u(:ncol,:pver)            * ztodt
    v_aft_PBL(:ncol,:pver)   =  state%v(:ncol,:pver)          + ptend%v(:ncol,:pver)            * ztodt

    call qsat(t_aftPBL(:ncol,:pver), state%pmid(:ncol,:pver), &
         tem2(:ncol,:pver), ftem(:ncol,:pver))
    ftem_aftPBL(:ncol,:pver) = qv_aft_PBL(:ncol,:pver) / ftem(:ncol,:pver) * 100._r8

    tten(:ncol,:pver)        = ( t_aftPBL(:ncol,:pver)    - state%t(:ncol,:pver) )              * rztodt     
    rhten(:ncol,:pver)       = ( ftem_aftPBL(:ncol,:pver) - ftem_prePBL(:ncol,:pver) )          * rztodt 

    ! -------------------------------------------------------------- !
    ! mass conservation check.........
    ! -------------------------------------------------------------- !
    if (diff_cnsrv_mass_check) then

       ! Conservation check
       nstep = get_nstep()
       
       do m = 1, pcnst
          fixed_ubc: if ((.not.cnst_fixed_ubc(m)).and.(.not.cnst_fixed_ubflx(m))) then
             col_loop: do i = 1, ncol
                sum1 = 0._r8
                sum2 = 0._r8
                sum3 = 0._r8
                do k = 1, pver
                   if(cnst_get_type_byind(m).eq.'wet') then
                      pdelx = state%pdel(i,k)
                   else
                      pdelx = state%pdeldry(i,k)
                   endif
                   sum1 = sum1 + state%q(i,k,m)*pdelx/gravit                          ! total column
                   sum2 = sum2 +(state%q(i,k,m)+ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! total column after tendancy is applied
                   sum3 = sum3 +(               ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! rate of change in column
                enddo
                sum1 = sum1 + (cflx(i,m) * ztodt) ! add in surface flux (kg/m2)
                sflx = (cflx(i,m) * ztodt) 
                if (sum1>1.e-36_r8) then
                   if( abs((sum2-sum1)/sum1) .gt. 1.e-12_r8  ) then
                      write(iulog,'(a,a8,a,I4,2f8.3,5e25.16)') &
                           'MASSCHECK vert diff : nstep,lon,lat,mass1,mass2,sum3,sflx,rel-diff : ', &
                           trim(cnst_name(m)), ' : ', nstep, state%lon(i)*180._r8/pi, state%lat(i)*180._r8/pi, &
                           sum1, sum2, sum3, sflx, abs(sum2-sum1)/sum1
                      call endrun('vertical_diffusion_tend : mass not conserved' )
                   endif
                endif
             enddo col_loop
          endif fixed_ubc
       enddo
    endif

    ! -------------------------------------------------------------- !
    ! Writing state variables after PBL scheme for detailed analysis !
    ! -------------------------------------------------------------- !

    call outfld( 'sl_aft_PBL'   , sl,                        pcols, lchnk )
    call outfld( 'qt_aft_PBL'   , qt,                        pcols, lchnk )
    call outfld( 'slv_aft_PBL'  , slv,                       pcols, lchnk )
    call outfld( 'u_aft_PBL'    , u_aft_PBL,                 pcols, lchnk )
    call outfld( 'v_aft_PBL'    , v_aft_PBL,                 pcols, lchnk )
    call outfld( 'qv_aft_PBL'   , qv_aft_PBL,                pcols, lchnk )
    call outfld( 'ql_aft_PBL'   , ql_aft_PBL,                pcols, lchnk )
    call outfld( 'qi_aft_PBL'   , qi_aft_PBL,                pcols, lchnk )
    call outfld( 't_aft_PBL '   , t_aftPBL,                  pcols, lchnk )
    call outfld( 'rh_aft_PBL'   , ftem_aftPBL,               pcols, lchnk )
    call outfld( 'slflx_PBL'    , slflx,                     pcols, lchnk )
    call outfld( 'qtflx_PBL'    , qtflx,                     pcols, lchnk )
    call outfld( 'uflx_PBL'     , uflx,                      pcols, lchnk )
    call outfld( 'vflx_PBL'     , vflx,                      pcols, lchnk )
    call outfld( 'slflx_cg_PBL' , slflx_cg,                  pcols, lchnk )
    call outfld( 'qtflx_cg_PBL' , qtflx_cg,                  pcols, lchnk )
    call outfld( 'uflx_cg_PBL'  , uflx_cg,                   pcols, lchnk )
    call outfld( 'vflx_cg_PBL'  , vflx_cg,                   pcols, lchnk )
    call outfld( 'slten_PBL'    , slten,                     pcols, lchnk )
    call outfld( 'qtten_PBL'    , qtten,                     pcols, lchnk )
    call outfld( 'uten_PBL'     , ptend%u(:ncol,:),          pcols, lchnk )
    call outfld( 'vten_PBL'     , ptend%v(:ncol,:),          pcols, lchnk )
    call outfld( 'qvten_PBL'    , ptend%q(:ncol,:,1),        pcols, lchnk )
    call outfld( 'qlten_PBL'    , ptend%q(:ncol,:,ixcldliq), pcols, lchnk )
    call outfld( 'qiten_PBL'    , ptend%q(:ncol,:,ixcldice), pcols, lchnk )
    call outfld( 'tten_PBL'     , tten,                      pcols, lchnk )
    call outfld( 'rhten_PBL'    , rhten,                     pcols, lchnk )

    ! ------------------------------------------- !
    ! Writing the other standard output variables !
    ! ------------------------------------------- !

    call outfld( 'QT'           , qt,                        pcols, lchnk )
    call outfld( 'SL'           , sl,                        pcols, lchnk )
    call outfld( 'SLV'          , slv,                       pcols, lchnk )
    call outfld( 'SLFLX'        , slflx,                     pcols, lchnk )
    call outfld( 'QTFLX'        , qtflx,                     pcols, lchnk )
    call outfld( 'UFLX'         , uflx,                      pcols, lchnk )
    call outfld( 'VFLX'         , vflx,                      pcols, lchnk )
    call outfld( 'TKE'          , tke,                       pcols, lchnk )

    call outfld( 'PBLH'         , pblh,                      pcols, lchnk )
    call outfld( 'TPERT'        , tpert,                     pcols, lchnk )
    call outfld( 'QPERT'        , qpert,                     pcols, lchnk )
    call outfld( 'USTAR'        , ustar,                     pcols, lchnk )
    call outfld( 'KVH'          , kvh,                       pcols, lchnk )
    call outfld( 'KVT'          , kvt,                       pcols, lchnk )
    call outfld( 'KVM'          , kvm,                       pcols, lchnk )
    call outfld( 'CGS'          , cgs,                       pcols, lchnk )
    dtk(:ncol,:) = dtk(:ncol,:) / cpair              ! Normalize heating for history
    call outfld( 'DTVKE'        , dtk,                       pcols, lchnk )
    dtk(:ncol,:) = ptend%s(:ncol,:) / cpair          ! Normalize heating for history using dtk
    call outfld( 'DTV'          , dtk,                       pcols, lchnk ) 
    call outfld( 'DUV'          , ptend%u,                   pcols, lchnk )
    call outfld( 'DVV'          , ptend%v,                   pcols, lchnk )
    do m = 1, pcnst
       call outfld( vdiffnam(m) , ptend%q(1,1,m),            pcols, lchnk )
    end do
    if( do_tms ) then
      ! Here, 'tautmsx,tautmsy' are implicit 'tms' that have been actually
      ! added into the atmosphere.
        call outfld( 'TAUTMSX'  , tautmsx,                   pcols, lchnk )
        call outfld( 'TAUTMSY'  , tautmsy,                   pcols, lchnk )
    end if
    if( do_molec_diff ) then
        call outfld( 'TTPXMLC'  , topflx,                    pcols, lchnk )
    end if

    return
  end subroutine vertical_diffusion_tend

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine positive_moisture( cp, xlv, xls, ncol, mkx, dt, qvmin, qlmin, qimin, & 
                                dp, qv, ql, qi, t, s, qvten, qlten, qiten, sten )
  ! ------------------------------------------------------------------------------- !
  ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
  ! force them to be larger than minimum value by (1) condensating water vapor      !
  ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
  ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
  ! Update final state variables and tendencies associated with this correction.    !
  ! If any condensation happens, update (s,t) too.                                  !
  ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
  ! input tendencies.                                                               !
  ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          ! 
  ! ------------------------------------------------------------------------------- !
    implicit none
    integer,  intent(in)     :: ncol, mkx
    real(r8), intent(in)     :: cp, xlv, xls
    real(r8), intent(in)     :: dt, qvmin, qlmin, qimin
    real(r8), intent(in)     :: dp(ncol,mkx)
    real(r8), intent(inout)  :: qv(ncol,mkx), ql(ncol,mkx), qi(ncol,mkx), t(ncol,mkx), s(ncol,mkx)
    real(r8), intent(inout)  :: qvten(ncol,mkx), qlten(ncol,mkx), qiten(ncol,mkx), sten(ncol,mkx)
    integer   i, k
    real(r8)  dql, dqi, dqv, sum, aa, dum 

  ! Modification : I should check whether this is exactly same as the one used in
  !                shallow convection and cloud macrophysics.

    do i = 1, ncol
       do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
          dql        = max(0._r8,1._r8*qlmin-ql(i,k))
          dqi        = max(0._r8,1._r8*qimin-qi(i,k))
          qlten(i,k) = qlten(i,k) +  dql/dt
          qiten(i,k) = qiten(i,k) +  dqi/dt
          qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
          sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
          ql(i,k)    = ql(i,k) +  dql
          qi(i,k)    = qi(i,k) +  dqi
          qv(i,k)    = qv(i,k) -  dql - dqi
          s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
          t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
          dqv        = max(0._r8,1._r8*qvmin-qv(i,k))
          qvten(i,k) = qvten(i,k) + dqv/dt
          qv(i,k)    = qv(i,k)    + dqv
          if( k .ne. 1 ) then 
              qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
              qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
          endif
          qv(i,k) = max(qv(i,k),qvmin)
          ql(i,k) = max(ql(i,k),qlmin)
          qi(i,k) = max(qi(i,k),qimin)
       end do
       ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
       ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
       ! preserves column moisture. 
       if( dqv .gt. 1.e-20_r8 ) then
           sum = 0._r8
           do k = 1, mkx
              if( qv(i,k) .gt. 2._r8*qvmin ) sum = sum + qv(i,k)*dp(i,k)
           enddo
           aa = dqv*dp(i,1)/max(1.e-20_r8,sum)
           if( aa .lt. 0.5_r8 ) then
               do k = 1, mkx
                  if( qv(i,k) .gt. 2._r8*qvmin ) then
                      dum        = aa*qv(i,k)
                      qv(i,k)    = qv(i,k) - dum
                      qvten(i,k) = qvten(i,k) - dum/dt
                  endif
               enddo 
           else 
               write(iulog,*) 'Full positive_moisture is impossible in vertical_diffusion'
           endif
       endif 
    end do
    return

  end subroutine positive_moisture

end module vertical_diffusion
