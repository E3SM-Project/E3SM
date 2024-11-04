
module vertical_diffusion
  !-------------------------------------------------------------------------------------------------
  ! Module to compute vertical diffusion of momentum,  moisture, trace constituents 
  ! and static energy. Separate modules compute                                     
  !   1. stresses associated with turbulent flow over orography (turbulent mountain stress)
  !   2. eddy diffusivities, including nonlocal tranport terms                    
  !   3. molecular diffusivities                                                  
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by    
  ! differencing the diffused and initial states.                                 
  !                                                                               
  ! Calling sequence:                                                          
  !  vertical_diffusion_init      Initializes vertical diffustion constants and modules 
  !        init_molec_diff        Initializes molecular diffusivity module              
  !        init_tms               Initializes turbulent mountain stress module          
  !        init_vdiff             Initializes diffusion solver module                   
  !  vertical_diffusion_ts_init   Time step init (only used for upper boundary condition)
  !  vertical_diffusion_tend      Computes vertical diffusion tendencies                           
  !        compute_tms            Computes turbulent mountain stresses                             
  !        compute_vdiff          Solves vertical diffusion, including molecular diffusivities
  !-------------------------------------------------------------------------------------------------
  use perf_mod
  use shr_kind_mod,     only : r8 => shr_kind_r8, i4=> shr_kind_i4
  use ppgrid,           only : pcols, pver, pverp
  use constituents,     only : pcnst, qmin, cnst_get_ind
  use diffusion_solver, only : vdiff_selector
  use cam_abortutils,   only : endrun
  use error_messages,   only : handle_errmsg
  use physconst,        only : cpair  , &     ! Specific heat of dry air
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
  use cam_logfile,      only : iulog
  use ref_pres,         only : do_molec_diff
  use phys_control,     only : phys_getopts, waccmx_is
  use time_manager,     only : is_first_step

  implicit none
  private      
  save
  
  !-----------------------------------------------------------------------------
  ! Public interfaces
  !-----------------------------------------------------------------------------
  public vd_readnl
  public vertical_diffusion_register                     ! Register multi-time-level variables with physics buffer
  public vertical_diffusion_init     ! Initialization
  public vertical_diffusion_ts_init  ! Time step initialization (only used for upper boundary condition)
  public vertical_diffusion_tend     ! Full vertical diffusion routine

  !-----------------------------------------------------------------------------
  ! Private data
  !-----------------------------------------------------------------------------
  character(len=16) :: MMF_microphysics_scheme
  logical           :: use_ECPP

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

  real(r8), parameter  :: unset_r8 = huge(1._r8)
  real(r8)             :: kv_top_pressure              ! Pressure defining the bottom of the upper atmosphere for kvh scaling (Pa)
  real(r8)             :: kv_top_scale                 ! Eddy diffusivity scale factor for upper atmosphere
  real(r8)             :: kv_freetrop_scale            ! Eddy diffusivity scale factor for the free troposphere
  real(r8)             :: eddy_lbulk_max               ! Maximum master length for diag_TKE
  real(r8)             :: eddy_leng_max                ! Maximum dissipation length for diag_TKE
  real(r8)             :: eddy_max_bot_pressure        ! Bottom pressure level (hPa) for eddy_leng_max
  real(r8)             :: eddy_moist_entrain_a2l = unset_r8 ! Moist entrainment enhancement param
  logical              :: diff_cnsrv_mass_check        ! do mass conservation check
  logical              :: do_tms                       ! switch for turbulent mountain stress
  logical              :: do_iss                       ! switch for implicit turbulent surface stress
  logical              :: prog_modal_aero = .false.    ! set true if prognostic modal aerosols are present
  integer              :: pmam_ncnst = 0               ! number of prognostic modal aerosol constituents
  integer, allocatable :: pmam_cnst_idx(:)             ! constituent indices of prognostic modal aerosols

contains

!===============================================================================
!===============================================================================

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
       eddy_max_bot_pressure, eddy_moist_entrain_a2l, diff_cnsrv_mass_check, do_iss
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
  call mpibcast(eddy_moist_entrain_a2l,          1 , mpir8,   0, mpicom)
  call mpibcast(diff_cnsrv_mass_check,           1 , mpilog,  0, mpicom)
  call mpibcast(do_iss,                          1 , mpilog,  0, mpicom)
#endif

end subroutine vd_readnl

!===============================================================================
!===============================================================================

subroutine vertical_diffusion_register()

  !------------------------------------------------ !
  ! Register physics buffer fields and constituents !
  !------------------------------------------------ !

  use physics_buffer,      only : pbuf_add_field, dtype_r8, dtype_i4

  ! Get settings from phys_control.F90
  call phys_getopts( MMF_microphysics_scheme_out = MMF_microphysics_scheme)
  call phys_getopts( use_ECPP_out = use_ECPP)
  call phys_getopts( do_tms_out = do_tms)

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

end subroutine vertical_diffusion_register

!===============================================================================
!===============================================================================

subroutine vertical_diffusion_init(pbuf2d)
  !-----------------------------------------------------------------------------
  ! Initialization of time independent fields for vertical diffusion
  ! Calls initialization routines for subsidiary modules            
  !-----------------------------------------------------------------------------
  use cam_history,       only : addfld, horiz_only, add_default
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
  !-----------------------------------------------------------------------------
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  character(128) :: errstring   ! Error status for init_vdiff
  integer        :: ntop_eddy   ! Top    interface to which eddy diffusion is applied ( = 1 )
  integer        :: nbot_eddy   ! Bottom interface to which eddy diffusion is applied ( = pver )
  integer        :: k           ! Vertical loop index
  integer        :: im, l, m, nmodes, nspec
  real(r8), parameter :: ntop_eddy_pres = 1.e-5_r8 ! Pressure below which eddy diffusion is not done in WACCM-X. (Pa)
  !-----------------------------------------------------------------------------
  if (masterproc) write(iulog,*)'Initializing vertical diffusion (vertical_diffusion_init)'

  ! Get indices of cloud liquid and ice within the constituents array
  call cnst_get_ind( 'CLDLIQ', ixcldliq )
  call cnst_get_ind( 'CLDICE', ixcldice )
  if( MMF_microphysics_scheme == 'p3' ) then
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

  !-----------------------------------------------------------------------------
  ! Initialize molecular diffusivity module                                    
  ! Note that computing molecular diffusivities is a trivial expense, but constituent
  ! diffusivities depend on their molecular weights. Decomposing the diffusion matric
  ! for each constituent is a needless expense unless the diffusivity is significant.
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Initialize molecular diffusion and get top and bottom limits
  !-----------------------------------------------------------------------------
  if( do_molec_diff ) then
     call init_molec_diff( r8, pcnst, rair, mwdry, avogad, gravit, cpair, boltz, errstring)

     call handle_errmsg(errstring, subname="init_molec_diff")

     call addfld( 'TTPXMLC', horiz_only, 'A', 'K/S', 'Top interf. temp. flux: molec. viscosity' )
     call add_default ( 'TTPXMLC', 1, ' ' )
     if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =', ntop_molec, 'NBOT_MOLEC =', nbot_molec
  end if

  !-----------------------------------------------------------------------------
  ! Initialize eddy diffusivity module
  !-----------------------------------------------------------------------------
  ! ntop_eddy must be 1 or <= nbot_molec (it is always 1 except for WACCM-X)
  ntop_eddy = 1
  nbot_eddy = pver

  if (masterproc) write(iulog, fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy
  if (masterproc) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'
  
  call init_hb_diff(gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, karman)

  call addfld('HB_ri',      (/ 'lev' /),  'A', 'no',  'Richardson Number (HB Scheme), I' )
  
  ! The vertical diffusion solver must operate over the full range of molecular and eddy diffusion
  ntop = min(ntop_molec,ntop_eddy)
  nbot = max(nbot_molec,nbot_eddy)
  
  !-----------------------------------------------------------------------------
  ! Initialize turbulent mountain stress module
  !-----------------------------------------------------------------------------
  if( do_tms ) then

     call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring )

     call handle_errmsg(errstring, subname="init_tms")

     call addfld( 'TAUTMSX' ,  horiz_only,  'A','N/m2',  'Zonal      turbulent mountain surface stress' )
     call addfld( 'TAUTMSY' ,  horiz_only,  'A','N/m2',  'Meridional turbulent mountain surface stress' )

     if (history_amwg) then
        call add_default('TAUTMSX ', 1, ' ')
        call add_default('TAUTMSY ', 1, ' ')
     end if

     if (masterproc) then
        write(iulog,*)'Using turbulent mountain stress module'
        write(iulog,*)'  tms_orocnst = ',tms_orocnst
        write(iulog,*)'  tms_z0fac = ',tms_z0fac
     end if

  endif ! do_tms
  
  !-----------------------------------------------------------------------------
  ! Initialize diffusion solver module
  !-----------------------------------------------------------------------------
  call init_vdiff( r8, iulog, rair, gravit, do_iss, errstring )
  call handle_errmsg(errstring, subname="init_vdiff")

  ! Use fieldlist_wet to select the fields which will be diffused using moist mixing ratios ( all by default )
  ! Use fieldlist_dry to select the fields which will be diffused using dry   mixing ratios.

  fieldlist_wet = new_fieldlist_vdiff( pcnst)
  fieldlist_dry = new_fieldlist_vdiff( pcnst)
  fieldlist_molec = new_fieldlist_vdiff( pcnst)

  ! Add all fields to fieldlist_wet
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

    ! Convert all constituents to wet before doing diffusion.
    if( vdiff_select( fieldlist_wet, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'q', k ) )
  
    ! Select constituents for molecular diffusion  
    if ( cnst_get_molec_byind(k) .eq. 'minor' ) then
      if( vdiff_select(fieldlist_molec,'q',k) .ne. '' ) call endrun( vdiff_select( fieldlist_molec,'q',k ) )
    endif

  end do constit_loop
  
  !-----------------------------------------------------------------------------
  ! Diagnostic output fields 
  !-----------------------------------------------------------------------------
  do k = 1, pcnst
     vdiffnam(k) = 'VD'//cnst_name(k)
     if( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
     call addfld( vdiffnam(k), (/ 'lev' /), 'A', 'kg/kg/s', 'Vertical diffusion of '//cnst_name(k) )
  end do

  call addfld( 'TKE'         , (/ 'ilev' /)  , 'A', 'm2/s2'  , 'Turbulent Kinetic Energy'                           )
  call addfld( 'PBLH'        , horiz_only    , 'A', 'm'      , 'PBL height'                                         )
  call addfld( 'TPERT'       , horiz_only    , 'A', 'K'      , 'Perturbation temperature (eddies in PBL)'           )
  call addfld( 'QPERT'       , horiz_only    , 'A', 'kg/kg'  , 'Perturbation specific humidity (eddies in PBL)'     )
  call addfld( 'USTAR'       , horiz_only    , 'A', 'm/s'    , 'Surface friction velocity'                          )
  call addfld( 'KVH'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion diffusivities (heat/moisture)'   )
  call addfld( 'KVM'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion diffusivities (momentum)'        )
  call addfld( 'KVT'         , (/ 'ilev' /)  , 'A', 'm2/s'   , 'Vertical diffusion kinematic molecular conductivity')
  call addfld( 'CGS'         , (/ 'ilev' /)  , 'A', 's/m2'   , 'Counter-gradient coeff on surface kinematic fluxes' )
  call addfld( 'DTVKE'       , (/ 'lev' /)   , 'A', 'K/s'    , 'dT/dt vertical diffusion KE dissipation'            )
  call addfld( 'DTV'         , (/ 'lev' /)   , 'A', 'K/s'    , 'T vertical diffusion'                               )
  call addfld( 'DUV'         , (/ 'lev' /)   , 'A', 'm/s2'   , 'U vertical diffusion'                               )
  call addfld( 'DVV'         , (/ 'lev' /)   , 'A', 'm/s2'   , 'V vertical diffusion'                               )
  call addfld( 'QT'          , (/ 'lev' /)   , 'A', 'kg/kg'  , 'Total water mixing ratio'                           )
  call addfld( 'SL'          , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liquid water static energy'                         )
  call addfld( 'SLV'         , (/ 'lev' /)   , 'A', 'J/kg'   , 'Liq wat virtual static energy'                      )
  call addfld( 'SLFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Liquid static energy flux'                          ) 
  call addfld( 'QTFLX'       , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Total water flux'                                   ) 
  call addfld( 'UFLX'        , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Zonal momentum flux'                                ) 
  call addfld( 'VFLX'        , (/ 'ilev' /)  , 'A', 'W/m2'   , 'Meridional momentm flux'                            ) 

  call addfld ('ustar',horiz_only, 'A',     ' ',' ')
  call addfld ('obklen',horiz_only, 'A',    ' ',' ')

  !-----------------------------------------------------------------------------
  ! determine default variables
  !-----------------------------------------------------------------------------
  call phys_getopts( history_amwg_out = history_amwg, &
                     history_eddy_out = history_eddy, &
                     history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num, &
                     history_waccm_out = history_waccm)

  if (history_amwg) then
    call add_default(  vdiffnam(1), 1, ' ' )
    call add_default( 'DTV'       , 1, ' ' )  
    call add_default( 'PBLH'      , 1, ' ' )
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
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

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

!===============================================================================
!===============================================================================

subroutine vertical_diffusion_ts_init( pbuf2d, state )
  !-----------------------------------------------------------------------------
  ! Timestep dependent setting,                                   
  ! At present only invokes upper bc code for molecular diffusion 
  !-----------------------------------------------------------------------------
  use molec_diff    , only : init_timestep_molec_diff
  use physics_types , only : physics_state
  use ppgrid        , only : begchunk, endchunk
  use physics_buffer, only : physics_buffer_desc
  !-----------------------------------------------------------------------------
  type(physics_state), intent(in) :: state(begchunk:endchunk)                 
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !-----------------------------------------------------------------------------
  if (do_molec_diff) call init_timestep_molec_diff(pbuf2d, state )

end subroutine vertical_diffusion_ts_init

!===============================================================================
!===============================================================================

subroutine vertical_diffusion_tend( ztodt    , state    ,                  &
                                    taux     , tauy     , shflx    , cflx, &
                                    ustar    , obklen   , ptend    , &
                                    cldn     , landfrac , sgh , pbuf) 
  !-----------------------------------------------------------------------------
  ! This is an interface routine for vertical diffusion
  !-----------------------------------------------------------------------------
  use physics_buffer,     only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
  use physics_types,      only : physics_state, physics_ptend, physics_ptend_init
  use physics_types,      only : set_dry_to_wet, set_wet_to_dry
  use cam_history,        only : outfld
  use trb_mtn_stress,     only : compute_tms
  use hb_diff,            only : compute_hb_diff
  use wv_saturation,      only : qsat
  use molec_diff,         only : compute_molec_diff, vd_lu_qdecomp
  use constituents,       only : qmincg, qmin, cnst_type
  use co2_cycle,          only : co2_cycle_set_cnst_type
  use diffusion_solver,   only : compute_vdiff, any, operator(.not.)
  use physconst,          only : cpairv, rairv ! Needed for calculation of upward H flux
  use time_manager,       only : get_nstep
  use constituents,       only : cnst_get_type_byind, cnst_name, cnst_fixed_ubc, cnst_fixed_ubflx
  use pbl_utils,          only : virtem, calc_obklen
  use physconst,          only : pi
  !-----------------------------------------------------------------------------
  ! Interface Arguments 
  type(physics_state), intent(inout) :: state                     ! Physics state variables
  real(r8),            intent(in)    :: taux(pcols)               ! x surface stress  [ N/m2 ]
  real(r8),            intent(in)    :: tauy(pcols)               ! y surface stress  [ N/m2 ]
  real(r8),            intent(in)    :: shflx(pcols)              ! Surface sensible heat flux  [ w/m2 ]
  real(r8),            intent(in)    :: cflx(pcols,pcnst)         ! Surface constituent flux [ kg/m2/s ]
  real(r8),            intent(in)    :: ztodt                     ! 2 delta-t [ s ]
  real(r8),            intent(in)    :: cldn(pcols,pver)          ! New stratus fraction [ fraction ]
  real(r8),            intent(in)    :: landfrac(pcols)           ! Land fraction
  real(r8),            intent(in)    :: sgh(pcols)                ! Standard deviation of orography [ unit ? ]
  type(physics_ptend), intent(out)   :: ptend                     ! Individual parameterization tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)
  real(r8),            intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
  real(r8),            intent(out)   :: obklen(pcols)             ! Obukhov length [ m ]
  !-----------------------------------------------------------------------------
  ! Local Variables
  character(128) :: errstring                      ! Error status for compute_vdiff
  real(r8), pointer, dimension(:,:) :: qrl         ! LW radiative cooling rate
  real(r8), pointer, dimension(:,:) :: wsedl       ! Sedimentation velocity of stratiform liquid cloud droplet [ m/s ]
  integer  :: lchnk                                ! Chunk identifier
  integer  :: ncol                                 ! Number of atmospheric columns
  integer  :: i, k, l, m                           ! loop iterators
  integer  :: ierr                                 ! status for allocate/deallocate

  real(r8) :: dtk(pcols,pver)                      ! T tendency from KE dissipation
  real(r8) :: dtvke(pcols,pver)                    ! T tendency from KE dissipation
  real(r8) :: dtv(pcols,pver)                      ! T tendency from KE dissipation
  real(r8), pointer   :: tke(:,:)                  ! Turbulent kinetic energy [ m2/s2 ]
  integer(i4),pointer :: turbtype(:,:)             ! Turbulent interface types [ no unit ]
  real(r8), pointer   :: smaw(:,:)                 ! Normalized Galperin instability function ( 0<= <=4.964 and 1 at neutral )

  real(r8) :: cgs(pcols,pverp)                     ! Counter-gradient star  [ cg/flux ]
  real(r8) :: cgh(pcols,pverp)                     ! Counter-gradient term for heat
  real(r8) :: rztodt                               ! 1./ztodt [ 1/s ]
  real(r8) :: ksrftms(pcols)                       ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
  real(r8) :: tautmsx(pcols)                       ! U component of turbulent mountain stress [ N/m2 ]
  real(r8) :: tautmsy(pcols)                       ! V component of turbulent mountain stress [ N/m2 ]
  real(r8) :: tautotx(pcols)                       ! U component of total surface stress [ N/m2 ]
  real(r8) :: tautoty(pcols)                       ! V component of total surface stress [ N/m2 ]

  real(r8):: shflx_tmp(pcols)                      ! Temporary surface sensible heat flux  [ w/m2 ]
  real(r8):: cflx_tmp(pcols,pcnst)                 ! Temporary surface constituent flux [ kg/m2/s ]

  real(r8), pointer :: kvt(:,:)                    ! Molecular kinematic conductivity for temperature [  ]
  real(r8) :: kvq(pcols,pverp)                     ! Eddy diffusivity for constituents [ m2/s ]
  real(r8) :: kvh(pcols,pverp)                     ! Eddy diffusivity for heat [ m2/s ]
  real(r8) :: kvm(pcols,pverp)                     ! Eddy diffusivity for momentum [ m2/s ]
  real(r8), pointer :: bprod(:,:)                  ! Buoyancy production of tke [ m2/s3 ]
  real(r8) :: sprod(pcols,pverp)                   ! Shear production of tke [ m2/s3 ]
  real(r8) :: sfi(pcols,pverp)                     ! Saturation fraction at interfaces [ fraction ]
  real(r8) :: sl(pcols,pver)
  real(r8) :: qt(pcols,pver)
  real(r8) :: slv(pcols,pver)
  real(r8) :: slvten(pcols,pver)
  real(r8) :: slflx(pcols,pverp)
  real(r8) :: qtflx(pcols,pverp)
  real(r8) :: uflx(pcols,pverp)
  real(r8) :: vflx(pcols,pverp)
  real(r8) :: th(pcols,pver)                        ! Potential temperature
  real(r8) :: topflx(pcols)                         ! Molecular heat flux at top interface
  real(r8) :: rhoair
  real(r8) :: ri(pcols,pver)                        ! richardson number (HB output)
  
  ! for obklen calculation outside HB
  real(r8) :: ftem(pcols,pver)                      ! Saturation vapor pressure before PBL
  real(r8) :: ftem_prePBL(pcols,pver)               ! Saturation vapor pressure before PBL
  real(r8) :: ftem_aftPBL(pcols,pver)               ! Saturation vapor pressure after PBL
  real(r8) :: tem2(pcols,pver)                      ! Saturation specific humidity and RH
  real(r8) :: t_aftPBL(pcols,pver)                  ! Temperature after PBL diffusion
  real(r8) :: tten(pcols,pver)                      ! Temperature tendency by PBL diffusion
  real(r8) :: rhten(pcols,pver)                     ! RH tendency by PBL diffusion
  real(r8) :: rairi(pcols,pver+1)                   ! interface gas constant needed for compute_vdiff
  real(r8),pointer :: tauresx(:)                    ! Residual stress to be added in vdiff to correct
  real(r8),pointer :: tauresy(:)                    ! for turb stress mismatch between sfc and atm accumulated.
  
  real(r8), pointer :: tpert(:)
  real(r8), pointer :: qpert(:)
  real(r8), pointer :: pblh(:)

  real(r8) :: dt_g_rdp(pcols)                       ! dt * g / dp

  integer  :: nstep
  real(r8) :: sum1, sum2, sum3, pdelx 
  real(r8) :: sflx

  ! Copy state so we can pass to intent(inout) routines that return new state instead of a tendency.
  real(r8) :: s_tmp(pcols,pver)
  real(r8) :: u_tmp(pcols,pver)
  real(r8) :: v_tmp(pcols,pver)
  real(r8) :: q_tmp(pcols,pver,pcnst)

  logical  :: lq(pcnst)

  character(len=3), dimension(pcnst) :: cnst_type_loc  ! local override option for constituents cnst_type

  !-----------------------------------------------------------------------------
  ! Main Computation Begins 
  !-----------------------------------------------------------------------------

  ! Assume 'wet' mixing ratios in diffusion code.
  ! don't convert co2 tracers to wet mixing ratios
  cnst_type_loc(:) = cnst_type(:)
  call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
  call set_dry_to_wet(state, cnst_type_loc)

  rztodt = 1._r8 / ztodt
  lchnk  = state%lchnk
  ncol   = state%ncol

  call pbuf_get_field(pbuf, tauresx_idx,  tauresx)
  call pbuf_get_field(pbuf, tauresy_idx,  tauresy)
  call pbuf_get_field(pbuf, tpert_idx,    tpert)
  call pbuf_get_field(pbuf, qpert_idx,    qpert)
  call pbuf_get_field(pbuf, pblh_idx,     pblh)
  call pbuf_get_field(pbuf, turbtype_idx, turbtype)

  !-----------------------------------------------------------------------------
  ! Computation of turbulent mountain stress 
  !-----------------------------------------------------------------------------
 
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

  !-----------------------------------------------------------------------------
  ! Computation of eddy diffusivities
  !-----------------------------------------------------------------------------
  call pbuf_get_field(pbuf, smaw_idx, smaw)
  call pbuf_get_field(pbuf, tke_idx,  tke)

  th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

  call compute_hb_diff( lchnk     , ncol     ,                            &
                        th        , state%t  , state%q , state%zm , state%zi, &
                        state%pmid, state%u  , state%v , tautotx  , tautoty , &
                        shflx     , cflx(:,1), obklen  , ustar    , pblh    , &
                        kvm       , kvh      , kvq     , cgh      , cgs     , &
                        tpert     , qpert    , cldn    , tke      , ri        )

  call outfld( 'HB_ri',         ri, pcols, lchnk )
  call outfld( 'ustar',   ustar(:), pcols, lchnk )
  call outfld( 'obklen', obklen(:), pcols, lchnk )

  !-----------------------------------------------------------------------------
  ! Application of diffusivities    
  !-----------------------------------------------------------------------------

  ! Set arrays from input state.
  q_tmp(:ncol,:,:) = state%q(:ncol,:,:)
  s_tmp(:ncol,:)   = state%s(:ncol,:)
  u_tmp(:ncol,:)   = state%u(:ncol,:)
  v_tmp(:ncol,:)   = state%v(:ncol,:)

  !-----------------------------------------------------------------------------
  ! Call the diffusivity solver and solve diffusion equation    
  ! The final two arguments are optional function references to 
  ! constituent-independent and -dependent moleculuar diffusivity routines
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !  Check to see if constituent dependent gas constant needed (WACCM-X) 
  !-----------------------------------------------------------------------------
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

  ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

  shflx_tmp = shflx
  cflx_tmp = cflx

#if defined( MMF_FLUX_BYPASS )
  shflx_tmp(:)  = 0.
  cflx_tmp(:,1) = 0.
#endif

  call pbuf_get_field(pbuf, kvt_idx, kvt)

  if( any(fieldlist_wet) ) then
      call compute_vdiff( state%lchnk   ,                                                                     &
                          pcols         , pver               , pcnst        , ncol          , state%pmid    , &
                          state%pint    , state%rpdel        , state%t      , ztodt         , taux          , &
                          tauy          , shflx_tmp          , cflx_tmp     , ntop          , nbot          , &
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
      if( do_molec_diff ) call endrun("Design flaw: dry vdiff not currently supported with molecular diffusion")
      call compute_vdiff( state%lchnk   ,                                                                     &
                          pcols         , pver               , pcnst        , ncol          , state%pmiddry , &
                          state%pintdry , state%rpdeldry     , state%t      , ztodt         , taux          , &       
                          tauy          , shflx_tmp          , cflx_tmp     , ntop          , nbot          , &       
                          kvh           , kvm                , kvq          , cgs           , cgh           , &   
                          state%zi      , ksrftms            , qmincg       , fieldlist_dry , fieldlist_molec,&
                          u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
                          tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
                          tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), rairi, &
                          do_molec_diff ,compute_molec_diff  , vd_lu_qdecomp )
      call handle_errmsg(errstring, subname="compute_vdiff", &
           extra_msg="Error in fieldlist_dry call from vertical_diffusion.")
  end if

  if (prog_modal_aero) then
    ! Modal aerosol species not diffused, so just add the explicit surface fluxes to the
    ! lowest layer.  **NOTE** This code assumes wet mmr.
    dt_g_rdp(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
    do m = 1, pmam_ncnst
      l = pmam_cnst_idx(m)
      q_tmp(:ncol,pver,l) = q_tmp(:ncol,pver,l) + dt_g_rdp(:ncol) * cflx(:ncol,l)
    enddo
  end if

  !-----------------------------------------------------------------------------
  ! Diagnostics and output writing after applying PBL scheme 
  !-----------------------------------------------------------------------------

  sl(:ncol,:pver)  = s_tmp(:ncol,:) -   latvap           * q_tmp(:ncol,:,ixcldliq) &
                                    - ( latvap + latice) * q_tmp(:ncol,:,ixcldice)
  qt(:ncol,:pver)  = q_tmp(:ncol,:,1) + q_tmp(:ncol,:,ixcldliq) + q_tmp(:ncol,:,ixcldice)
  slv(:ncol,:pver) = sl(:ncol,:pver) * ( 1._r8 + zvir*qt(:ncol,:pver) ) 

  slflx(:ncol,1)    = 0._r8
  qtflx(:ncol,1)    = 0._r8
  uflx(:ncol,1)     = 0._r8
  vflx(:ncol,1)     = 0._r8

  do k = 2, pver
     do i = 1, ncol
        rhoair     = state%pint(i,k) / ( rair * ( ( 0.5_r8*(slv(i,k)+slv(i,k-1)) - gravit*state%zi(i,k))/cpair ) )
        slflx(i,k) = kvh(i,k) * ( - rhoair*(sl(i,k-1)-sl(i,k))/(state%zm(i,k-1)-state%zm(i,k)) + cgh(i,k) ) 
        qtflx(i,k) = kvh(i,k) * ( - rhoair*(qt(i,k-1)-qt(i,k))/(state%zm(i,k-1)-state%zm(i,k)) &
                                  + rhoair*(cflx(i,1)+cflx(i,ixcldliq)+cflx(i,ixcldice))*cgs(i,k) )
        uflx(i,k)  = kvm(i,k) * ( - rhoair*(u_tmp(i,k-1)-u_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
        vflx(i,k)  = kvm(i,k) * ( - rhoair*(v_tmp(i,k-1)-v_tmp(i,k))/(state%zm(i,k-1)-state%zm(i,k)))
     end do
  end do

  slflx(:ncol,pverp) = shflx(:ncol)
  qtflx(:ncol,pverp) = cflx(:ncol,1)
  uflx(:ncol,pverp)  = tautotx(:ncol)
  vflx(:ncol,pverp)  = tautoty(:ncol)

  !-----------------------------------------------------------------------------
  ! Convert the new profiles into vertical diffusion tendencies.    
  ! Convert KE dissipative heat change into "temperature" tendency. 
  !-----------------------------------------------------------------------------

  lq(:) = .TRUE.
  call physics_ptend_init(ptend,state%psetcols, "vertical diffusion", ls=.true., lu=.true., lv=.true., lq=lq)

  ptend%s(:ncol,:)       = ( s_tmp(:ncol,:) - state%s(:ncol,:) ) * rztodt
  ptend%u(:ncol,:)       = ( u_tmp(:ncol,:) - state%u(:ncol,:) ) * rztodt
  ptend%v(:ncol,:)       = ( v_tmp(:ncol,:) - state%v(:ncol,:) ) * rztodt
  ptend%q(:ncol,:pver,:) = ( q_tmp(:ncol,:pver,:) - state%q(:ncol,:pver,:) ) * rztodt

  ! Convert tendencies of dry constituents to dry basis.
  do m = 1,pcnst
     if (cnst_type(m).eq.'dry') then
        ptend%q(:ncol,:pver,m) = ptend%q(:ncol,:pver,m)*state%pdel(:ncol,:pver)/state%pdeldry(:ncol,:pver)
     endif
  end do

  ! convert wet mmr back to dry before conservation check avoid converting co2 tracers again
  cnst_type_loc(:) = cnst_type(:)
  call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
  call set_wet_to_dry(state, cnst_type_loc)

  !-----------------------------------------------------------------------------
  ! mass conservation check
  !-----------------------------------------------------------------------------
  if (diff_cnsrv_mass_check) then

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
                 endif
              endif
           enddo col_loop
        endif fixed_ubc
     enddo
  endif

  !-----------------------------------------------------------------------------
  ! Writing the other standard output variables 
  !-----------------------------------------------------------------------------
  dtvke(:ncol,:) = dtk(:ncol,:) / cpair      ! Normalize heating for history
  dtv(:ncol,:)   = ptend%s(:ncol,:) / cpair  ! Normalize heating for history using dtk

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
  call outfld( 'DTVKE'        , dtvke,                     pcols, lchnk )
  call outfld( 'DTV'          , dtv,                       pcols, lchnk ) 
  call outfld( 'DUV'          , ptend%u,                   pcols, lchnk )
  call outfld( 'DVV'          , ptend%v,                   pcols, lchnk )
  do m = 1, pcnst
    call outfld( vdiffnam(m) , ptend%q(1,1,m),             pcols, lchnk )
  end do
  ! Here, 'tautmsx,tautmsy' are implicit 'tms' that have been actually added
  if( do_tms )        call outfld( 'TAUTMSX'  , tautmsx,   pcols, lchnk )
  if( do_tms )        call outfld( 'TAUTMSY'  , tautmsy,   pcols, lchnk )
  if( do_molec_diff ) call outfld( 'TTPXMLC'  , topflx,    pcols, lchnk )
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  return
end subroutine vertical_diffusion_tend

!===============================================================================
!===============================================================================

end module vertical_diffusion
