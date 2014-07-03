module clubb_intr

  !----------------------------------------------------------------------------------------------------- !
  ! Module to interface CAM with Cloud Layers Unified by Bi-normals (CLUBB), developed                   !
  !    by the University of Wisconsin Milwaukee Group (UWM).                                             !
  !                                                                                                      !
  ! CLUBB replaces the exisiting turbulence, shallow convection, and macrophysics in CAM5                !  
  !                                                                                                      !  
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      ! 
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! Authors:  P. Bogenschutz, C. Craig, A. Gettelmen                                                     ! 
  !                                                                                                      ! 
  !----------------------------------------------------------------------------------------------------- !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pver, pverp
  use phys_control,  only: phys_getopts
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, rh2o, karman, tms_orocnst, tms_z0fac
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc 
  use constituents,  only: pcnst
  use pbl_utils,     only: calc_ustar, calc_obklen
  use mpishorthand

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, &
#ifdef CLUBB_SGS
            ! This utilizes CLUBB specific variables in its interface
            stats_init_clubb, &
#endif
            stats_end_timestep_clubb, & 
	    clubb_surface

#ifdef CLUBB_SGS	  
  ! Both of these utilize CLUBB specific variables in their interface
  private :: stats_zero, stats_avg
#endif


  ! ------------ !
  ! Private data !
  ! ------------ !

  integer, parameter :: &
      grid_type    = 3, &  		! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0     		! The hydromet array in SAM-CLUBB is currently 0 elements
   
  real(r8), dimension(0) :: &
      sclr_tol = 1.e-8_r8               ! Total water in kg/kg

  character(len=6), parameter :: &
      saturation_equation = "flatau" 	! Flatau polynomial approximation for SVP

  real(r8), parameter :: &
      theta0   = 300._r8, &		! Reference temperature                     [K]
      ts_nudge = 86400._r8, & 		! Time scale for u/v nudging (not used)     [s]
      p0_clubb = 100000._r8
      
  real(r8), parameter :: &
      host_dx = 100000._r8, &		! Host model deltax [m]			     
      host_dy = 100000._r8		! Host model deltay [m]
      
  integer, parameter :: & 
    sclr_dim = 0   			! Higher-order scalars, set to zero

!  Constant parameters
  logical, parameter, private :: &
    l_uv_nudge       = .false.,       &  ! Use u/v nudging (not used)
    l_implemented    = .true.,        &  ! Implemented in a host model (always true)
    l_host_applies_sfc_fluxes = .false.  ! Whether the host model applies the surface fluxes
    
  logical            :: do_tms
  logical            :: lq(pcnst)
  logical            :: prog_modal_aero

  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB
 
!  define physics buffer indicies here       
  integer :: &
    wp2_idx, &         			! vertical velocity variances
    wp3_idx, &         			! third moment of vertical velocity
    wpthlp_idx, &      			! turbulent flux of thetal
    wprtp_idx, &       			! turbulent flux of total water
    rtpthlp_idx, &     			! covariance of thetal and rt
    rtp2_idx, &        			! variance of total water
    thlp2_idx, &       			! variance of thetal
    up2_idx, &         			! variance of east-west wind
    vp2_idx, &         			! variance of north-south wind
    upwp_idx, &        			! east-west momentum flux
    vpwp_idx, &        			! north-south momentum flux
    thlm_idx, &        			! mean thetal
    rtm_idx, &         			! mean total water mixing ratio
    um_idx, &         			! mean of east-west wind
    vm_idx, &           		! mean of north-south wind
    cld_idx, &         			! Cloud fraction
    concld_idx, &       		! Convective cloud fraction
    ast_idx, &          		! Stratiform cloud fraction
    alst_idx, &         		! Liquid stratiform cloud fraction
    aist_idx, & 			! Ice stratiform cloud fraction
    qlst_idx, &         		! Physical in-cloud LWC
    qist_idx, &         		! Physical in-cloud IWC
    dp_frac_idx, &      		! deep convection cloud fraction
    sh_frac_idx, &      		! shallow convection cloud fraction
    rel_idx, &             		! Rel
    kvh_idx, &			        ! CLUBB eddy diffusivity on thermo levels
    kvm_idx, &				! CLUBB eddy diffusivity on mom levels
    pblh_idx, &                         ! PBL pbuf
    icwmrdp_idx, &			! In cloud mixing ratio for deep convection
    tke_idx, &                          ! turbulent kinetic energy
    tpert_idx, &                        ! temperature perturbation from PBL
    fice_idx, &                         ! fice_idx index in physics buffer
    relvar_idx, &                       ! relative cloud water variance
    accre_enhan_idx                     ! optional accretion enhancement factor for MG

  !  Output arrays for CLUBB statistics    
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90

  contains
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_register_cam( )
!-------------------------------------------------------------------------------
! Description:
!   Register the constituents and fields in the physics buffer
! Author: P. Bogenschutz, C. Craig, A. Gettelman
!
!-------------------------------------------------------------------------------
#ifdef CLUBB_SGS

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    !  Add CLUBB fields to pbuf 
    use physics_buffer,  only: pbuf_add_field, dtype_r8, pbuf_times
    use ppgrid,          only: pver, pverp, pcols
    
    call phys_getopts( eddy_scheme_out = eddy_scheme, &
                       do_tms_out      = do_tms)	      

    !  put pbuf_add calls here (see macrop_driver.F90 for sample) use indicies defined at top
    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/),                    pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/),             tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/),             kvh_idx)
    call pbuf_add_field('kvm',        'global', dtype_r8, (/pcols, pverp/),             kvm_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/),                    tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,pbuf_times/),    ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,pbuf_times/),    aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,pbuf_times/),    alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,pbuf_times/),    qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,pbuf_times/),    qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,pbuf_times/),    concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,pbuf_times/),    cld_idx)
    call pbuf_add_field('FICE',       'physpkg',dtype_r8, (/pcols,pver/),               fice_idx)
    
    call pbuf_add_field('WP2',        'global', dtype_r8, (/pcols,pverp,pbuf_times/), wp2_idx)
    call pbuf_add_field('WP3',        'global', dtype_r8, (/pcols,pverp,pbuf_times/), wp3_idx)
    call pbuf_add_field('WPTHLP',     'global', dtype_r8, (/pcols,pverp,pbuf_times/), wpthlp_idx)
    call pbuf_add_field('WPRTP',      'global', dtype_r8, (/pcols,pverp,pbuf_times/), wprtp_idx)
    call pbuf_add_field('RTPTHLP',    'global', dtype_r8, (/pcols,pverp,pbuf_times/), rtpthlp_idx)
    call pbuf_add_field('RTP2',       'global', dtype_r8, (/pcols,pverp,pbuf_times/), rtp2_idx)
    call pbuf_add_field('THLP2',      'global', dtype_r8, (/pcols,pverp,pbuf_times/), thlp2_idx)
    call pbuf_add_field('UP2',        'global', dtype_r8, (/pcols,pverp,pbuf_times/), up2_idx)
    call pbuf_add_field('VP2',        'global', dtype_r8, (/pcols,pverp,pbuf_times/), vp2_idx)
    call pbuf_add_field('UPWP',       'global', dtype_r8, (/pcols,pverp,pbuf_times/), upwp_idx)
    call pbuf_add_field('VPWP',       'global', dtype_r8, (/pcols,pverp,pbuf_times/), vpwp_idx)
    call pbuf_add_field('THLM',       'global', dtype_r8, (/pcols,pverp,pbuf_times/), thlm_idx)
    call pbuf_add_field('RTM',        'global', dtype_r8, (/pcols,pverp,pbuf_times/), rtm_idx)
    call pbuf_add_field('UM',         'global', dtype_r8, (/pcols,pverp,pbuf_times/), um_idx)
    call pbuf_add_field('VM',         'global', dtype_r8, (/pcols,pverp,pbuf_times/), vm_idx)

#endif 

  end subroutine clubb_register_cam
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf2d)
!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------



#ifdef CLUBB_SGS

    !  From CAM libraries
    use physics_types, 		only: physics_state, physics_ptend
    use cam_history,     	only: addfld, add_default, phys_decomp
    use ppgrid, 		only: pver, pverp, pcols
    use ref_pres,               only: pref_mid
    use hb_diff,                only: init_hb_diff
    use trb_mtn_stress,         only: init_tms
    use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx

    !  From the CLUBB libraries
    use clubb_core, 	        only: setup_clubb_core
    use clubb_precision, 	only: time_precision
    use error_code, 		only: set_clubb_debug_level ! Subroutines
    use parameter_indices, 	only: nparams ! Constant
    use parameters_tunable, 	only: read_parameters ! Subroutine
    use stats_variables,        only: l_stats, l_stats_samp, l_grads, l_output_rad_files, &
                                      zt, zm, sfc, rad_zt, rad_zm
    use namelist_utils,         only: find_group_name
    use units,                  only: getunit, freeunit
    use abortutils,             only: endrun
    use error_messages,         only: handle_errmsg
    use time_manager,         	only: is_first_step
    use constants_clubb,        only: em_min, w_tol_sqd, rt_tol, thl_tol  


    !  These are only needed if we're using a passive scalar
    use array_index, 		only: iisclr_rt, iisclr_thl, iisclr_CO2, &    ! [kg/kg]/[K]/[1e6 mol/mol]
      				      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
    use constituents,           only: cnst_get_ind
    use phys_control,           only: phys_getopts

#endif

    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

#ifdef CLUBB_SGS

    real(kind=time_precision) :: dum1, dum2, dum3
    
    real(8), dimension(nparams)  :: clubb_params    ! These adjustable CLUBB parameters (C1, C2 ...)

    logical :: clubb_history, clubb_rad_history ! Stats enabled (T/F)

    character(len=128) :: errstring             ! error status for CLUBB init

    integer :: err_code, iunit 			! Code for when CLUBB fails
    integer :: i, j, k, l                       ! Indices
    integer :: read_status      		! Length of a string
    integer :: ntop_eddy    			! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
    integer :: nbot_eddy   			! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
    integer :: nmodes, nspec, pmam_ncnst, m
    integer :: ixnumliq
    integer :: lptr

    
    real(r8)  :: zt_g(pverp)                        ! Height dummy array
    real(r8)  :: zi_g(pverp)                        ! Height dummy array
      

    namelist /clubb_his_nl/ clubb_history, clubb_rad_history

    !----- Begin Code -----

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that  
    ! CLUBB does not transport aerosol consituents.  Therefore, need to 
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents) 
    ! ----------------------------------------------------------------- !

    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

   !  Select variables to apply tendencies back to CAM
 
   ! Initialize all consituents to true to start
   lq(1:pcnst) = .true.
   edsclr_dim  = pcnst

   if (prog_modal_aero) then
      ! Turn off modal aerosols and decrement edsclr_dim accordingly
      call rad_cnst_get_info(0, nmodes=nmodes)

      do m = 1, nmodes
         call rad_cnst_get_mode_num_idx(m, lptr)
         lq(lptr)=.false.
         edsclr_dim = edsclr_dim-1

         call rad_cnst_get_info(0, m, nspec=nspec)
         do l = 1, nspec
            call rad_cnst_get_mam_mmr_idx(m, l, lptr)
            lq(lptr)=.false.
            edsclr_dim = edsclr_dim-1
         end do
      end do


      !  In addition, if running with MAM, droplet number is transported
      !  in dropmixnuc, therefore we do NOT want CLUBB to apply transport
      !  tendencies to avoid double counted.  Else, we apply tendencies.
      call cnst_get_ind('NUMLIQ',ixnumliq)
      lq(ixnumliq)	= .false.
      edsclr_dim = edsclr_dim-1
   endif

    ! ----------------------------------------------------------------- !
    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    ! ----------------------------------------------------------------- !
    call set_clubb_debug_level( 0 )

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !

    !  Determine if we want clubb_history to be output  
    clubb_history      = .false.   ! Initialize to false
    l_stats            = .false.   ! Initialize to false
    l_output_rad_files = .false.   ! Initialize to false

    !  Read namelist to determine if CLUBB history should be called
    if (masterproc) then
      iunit= getunit()
      open(unit=iunit,file="atm_in",status='old')
      call find_group_name(iunit, 'clubb_his_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_his_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_tend_cam:  error reading namelist')
         end if
      end if
      close(unit=iunit) 
      call freeunit(iunit)
    end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(clubb_history,      1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rad_history,  1,   mpilog,   0, mpicom)
#endif

    !  Overwrite defaults if they are true
    if (clubb_history) l_stats = .true.
    if (clubb_rad_history) l_output_rad_files = .true. 

    !  Defaults
    l_stats_samp = .false.
    l_grads = .false.

    !  Overwrite defaults if needbe     
    if (l_stats) l_stats_samp = .true.

    !  Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC 
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG

    iisclr_rt  = -1
    iisclr_thl = -1
    iisclr_CO2 = -1

    iiedsclr_rt  = -1
    iiedsclr_thl = -1
    iiedsclr_CO2 = -1
    
    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !
    
    !  Read in parameters for CLUBB.  Just read in default values 
      call read_parameters( -99, "", clubb_params )
      
    !  Fill in dummy arrays for height.  Note that these are overwrote
    !  at every CLUBB step to physical values.    
    do k=1,pverp
      zt_g(k) = ((k-1)*1000._r8)-500._r8  !  this is dummy garbage
      zi_g(k) = (k-1)*1000._r8            !  this is dummy garbage
    enddo
   
    !  Set up CLUBB core.  Note that some of these inputs are overwrote
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.     
    call setup_clubb_core     &
         ( pverp, theta0, ts_nudge, &                                 ! In
           hydromet_dim,  sclr_dim, &                                 ! In
           sclr_tol, edsclr_dim, clubb_params, &                      ! In
           l_host_applies_sfc_fluxes, &                               ! In
           l_uv_nudge, saturation_equation,  &                        ! In
           l_implemented, grid_type, zi_g(2), zi_g(1), zi_g(pverp), & ! In
           zi_g(1:pverp), zt_g(1:pverp), &                            ! In
           host_dx, host_dy, zi_g(1), &                               ! In
           err_code )
	   
    ! ----------------------------------------------------------------- !
    ! Set-up HB diffusion.  Only intialized to diagnose PBL depth       !
    ! ----------------------------------------------------------------- !

    ! Initialize eddy diffusivity module
    
    ntop_eddy = 1    ! if >1, must be <= nbot_molec
    nbot_eddy = pver ! currently always pver
    
    call init_hb_diff( gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, karman, eddy_scheme )
    
    ! ----------------------------------------------------------------- !
    ! Initialize turbulent mountain stress module                       !
    ! ------------------------------------------------------------------!

    if ( do_tms) then    
      call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring)
      call handle_errmsg(errstring, subname="init_tms")

      call addfld( 'TAUTMSX' ,'N/m2  ',  1,  'A',  'Zonal      turbulent mountain surface stress',  phys_decomp )
      call addfld( 'TAUTMSY' ,'N/m2  ',  1,  'A',  'Meridional turbulent mountain surface stress',  phys_decomp )
      call add_default( 'TAUTMSX ', 1, ' ' )
      call add_default( 'TAUTMSY ', 1, ' ' )
      if (masterproc) then
         write(iulog,*)'Using turbulent mountain stress module'
         write(iulog,*)'  tms_orocnst = ',tms_orocnst
         write(iulog,*)'  tms_z0fac = ',tms_z0fac
      end if
    endif   

    ! ----------------------------------------------------------------- !
    ! Add output fields for the history files
    ! ----------------------------------------------------------------- !

     !  These are default CLUBB output.  Not the higher order history budgets
     call addfld ('RHO_CLUBB',        'kg/m3',    pverp, 'A', 'Air Density', phys_decomp)
     call addfld ('UP2_CLUBB',        'm2/s2',    pverp, 'A', 'Zonal Velocity Variance', phys_decomp)
     call addfld ('VP2_CLUBB',        'm2/s2',    pverp, 'A', 'Meridional Velocity Variance', phys_decomp)
     call addfld ('WP2_CLUBB',        'm2/s2',    pverp, 'A', 'Vertical Velocity Variance', phys_decomp)
     call addfld ('UPWP_CLUBB',       'm2/s2',    pverp, 'A', 'Zonal Momentum Flux', phys_decomp)
     call addfld ('VPWP_CLUBB',       'm2/s2',    pverp, 'A', 'Meridional Momentum Flux', phys_decomp)
     call addfld ('WP3_CLUBB',        'm3/s3',    pverp, 'A', 'Third Moment Vertical Velocity', phys_decomp)
     call addfld ('WPTHLP_CLUBB',     'W/m2',     pverp, 'A', 'Heat Flux', phys_decomp)
     call addfld ('WPRTP_CLUBB',      'W/m2',     pverp, 'A', 'Moisture Flux', phys_decomp)
     call addfld ('RTP2_CLUBB',       'g^2/kg^2', pverp, 'A', 'Moisture Variance', phys_decomp)
     call addfld ('THLP2_CLUBB',      'K^2',      pverp, 'A', 'Temperature Variance', phys_decomp)
     call addfld ('RTPTHLP_CLUBB',    'K g/kg',   pverp, 'A', 'Temp. Moist. Covariance', phys_decomp)
     call addfld ('RCM_CLUBB',        'g/kg',     pverp, 'A', 'Cloud Water Mixing Ratio', phys_decomp)
     call addfld ('WPRCP_CLUBB',      'W/m2',     pverp, 'A', 'Liquid Water Flux', phys_decomp)
     call addfld ('CLOUDFRAC_CLUBB',  'fraction', pver,  'A', 'Cloud Fraction', phys_decomp)
     call addfld ('RCMINLAYER_CLUBB', 'g/kg',     pverp, 'A', 'Cloud Water in Layer', phys_decomp)
     call addfld ('CLOUDCOVER_CLUBB', 'fraction', pverp, 'A', 'Cloud Cover', phys_decomp) 
     call addfld ('WPTHVP_CLUBB',     'W/m2',     pver,  'A', 'Buoyancy Flux',phys_decomp)
     call addfld ('RVMTEND_CLUBB',    'g/kg /s',  pver,  'A', 'Water vapor tendency',phys_decomp)
     call addfld ('STEND_CLUBB',      'k/s',      pver,  'A', 'Temperature tendency',phys_decomp)
     call addfld ('RCMTEND_CLUBB',    'g/kg /s',  pver,  'A', 'Cloud Liquid Water Tendency',phys_decomp)
     call addfld ('RIMTEND_CLUBB',    'g/kg /s',  pver,  'A', 'Cloud Ice Tendency',phys_decomp)
     call addfld ('UTEND_CLUBB',      'm/s /s',   pver,  'A', 'U-wind Tendency',phys_decomp)
     call addfld ('VTEND_CLUBB',      'm/s /s',   pver,  'A', 'V-wind Tendency',phys_decomp)
     call addfld ('ZT_CLUBB',         'm',        pverp, 'A', 'Thermodynamic Heights',phys_decomp)
     call addfld ('ZM_CLUBB',         'm',        pverp, 'A', 'Momentum Heights',phys_decomp)     
     call addfld ('UM_CLUBB',         'm/s',      pverp, 'A', 'Zonal Wind',phys_decomp)
     call addfld ('VM_CLUBB',         'm/s',      pverp, 'A', 'Meridional Wind',phys_decomp)
     call addfld ('THETAL',           'K',        pver,  'A', 'Liquid Water Potential Temperature',phys_decomp)
     call addfld ('PBLH',             'm',        1,     'A', 'PBL height',phys_decomp)
     call addfld ('QT',               'kg/kg',    pver,  'A', 'Total water mixing ratio',phys_decomp)
     call addfld ('SL',               'J/kg',     pver,  'A', 'Liquid water static energy',phys_decomp)
     call addfld ('CLDST',            'fraction', pver,  'A', 'Stratus cloud fraction',phys_decomp)
     call addfld ('ZMDLF',            'kg/kg/s',  pver,  'A', 'Detrained liquid water from ZM convection',phys_decomp)

     call addfld ('CONCLD   ',        'fraction', pver,  'A', 'Convective cloud cover',phys_decomp)
     
     !  Initialize statistics, below are dummy variables
      dum1 = 300._r8
      dum2 = 1200._r8
      dum3 = 300._r8
			 
      if (l_stats) then 
      
        call stats_init_clubb( .true., dum1, dum2, &
                         pverp, pverp, pverp, dum3 )
			 
	allocate(out_zt(pcols,pverp,zt%nn))
        allocate(out_zm(pcols,pverp,zm%nn))
        allocate(out_sfc(pcols,1,sfc%nn))
	
        allocate(out_radzt(pcols,pverp,rad_zt%nn))
        allocate(out_radzm(pcols,pverp,rad_zm%nn))			 
	
       endif
  
    ! ----------------------------------------------------------------- !
    ! Make all of this output default, this is not CLUBB history
    ! ----------------------------------------------------------------- !
 
     call add_default('RHO_CLUBB',        1, ' ')
     call add_default('UP2_CLUBB',        1, ' ')
     call add_default('VP2_CLUBB',        1, ' ')
     call add_default('WP2_CLUBB',        1, ' ')
     call add_default('WP3_CLUBB',        1, ' ')
     call add_default('UPWP_CLUBB',       1, ' ')
     call add_default('VPWP_CLUBB',       1, ' ')
     call add_default('WPTHLP_CLUBB',     1, ' ')
     call add_default('WPRTP_CLUBB',      1, ' ')
     call add_default('RTP2_CLUBB',       1, ' ')
     call add_default('THLP2_CLUBB',      1, ' ')
     call add_default('RTPTHLP_CLUBB',    1, ' ')
     call add_default('RCM_CLUBB',        1, ' ')
     call add_default('WPRCP_CLUBB',      1, ' ')
     call add_default('CLOUDFRAC_CLUBB',  1, ' ')
     call add_default('RCMINLAYER_CLUBB', 1, ' ')
     call add_default('CLOUDCOVER_CLUBB', 1, ' ')
     call add_default('WPTHVP_CLUBB',     1, ' ')
     call add_default('RVMTEND_CLUBB',    1, ' ')
     call add_default('STEND_CLUBB',      1, ' ')
     call add_default('RCMTEND_CLUBB',    1, ' ')
     call add_default('RIMTEND_CLUBB',    1, ' ')
     call add_default('UTEND_CLUBB',      1, ' ')
     call add_default('VTEND_CLUBB',      1, ' ')
     call add_default('ZT_CLUBB',         1, ' ')
     call add_default('ZM_CLUBB',         1, ' ')   
     call add_default('UM_CLUBB',         1, ' ')
     call add_default('VM_CLUBB',         1, ' ')
     call add_default('PBLH',             1, ' ')
     call add_default('SL',               1, ' ')
     call add_default('QT',               1, ' ')
     call add_default('CONCLD',           1, ' ')
       

    ! --------------- !
    ! First step?     !
    ! Initialization  !
    ! --------------- !

   !  Is this the first time step?  If so then initialize CLUBB variables as follows
   if (is_first_step()) then
 
      call pbuf_set_field(pbuf2d, wp2_idx,     w_tol_sqd)
      call pbuf_set_field(pbuf2d, wp3_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, wpthlp_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d, wprtp_idx,   0.0_r8)
      call pbuf_set_field(pbuf2d, rtpthlp_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, rtp2_idx,    rt_tol**2)
      call pbuf_set_field(pbuf2d, thlp2_idx,   thl_tol**2)
      call pbuf_set_field(pbuf2d, up2_idx,     w_tol_sqd)
      call pbuf_set_field(pbuf2d, vp2_idx,     w_tol_sqd)
      call pbuf_set_field(pbuf2d, upwp_idx,    0.0_r8)
      call pbuf_set_field(pbuf2d, vpwp_idx,    0.0_r8)
      call pbuf_set_field(pbuf2d, tke_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, kvh_idx,     0.0_r8)
      call pbuf_set_field(pbuf2d, fice_idx,    0.0_r8)

   endif
   
    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

#endif
    end subroutine clubb_ini_cam
    
    
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    subroutine clubb_tend_cam( &
                              state,   ptend_all,   pbuf,     hdtime, &
                              cmfmc,   cmfmc2,      cam_in,   sgh30,  dlf,   &
			      det_s, det_ice)
			      
!-------------------------------------------------------------------------------
! Description: Provide tendencies of shallow convection, turbulence, and 
!              macrophysics from CLUBB to CAM
!   
! Author: Cheryl Craig, March 2011
! Modifications: Pete Bogenschutz, March 2011 and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types, 	only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init, &
			      physics_ptend_sum, physics_update

    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              pbuf_set_field, physics_buffer_desc, pbuf_times
			      
    use ppgrid, 	only: pver, pverp, pcols
    use constituents, 	only: cnst_get_ind
    use camsrfexch,     only: cam_in_t
      
#ifdef CLUBB_SGS
    use hb_diff,        only: pblintd
    use scamMOD,	only: single_column,scm_clubb_iop_name
    use phys_grid, 	only: get_lat_p
    use parameter_indices, only: nparams     
    use parameters_tunable, only: read_parameters, setup_parameters ! Subroutine
    use cldwat2m_macro, only: aist_vector
    use clubb_precision,only: time_precision
    use cam_history, 	only: outfld
    use clubb_core, 	only: advance_clubb_core
    use grid_class, 	only: zt2zm, gr, setup_grid, cleanup_grid 
    use constants_clubb,only: em_min, w_tol_sqd, rt_tol, thl_tol  
    use model_flags, 	only: l_use_boussinesq 
    use stats_variables,        only: l_stats, stats_tsamp, stats_tout, zt, &
                                      sfc, zm, rad_zt, rad_zm, l_output_rad_files
    use pdf_parameter_module, only: pdf_parameter ! Type
    use saturation, 	only: sat_mixrat_liq
    use trb_mtn_stress, only: compute_tms
    use stats_subs,     only: stats_begin_timestep
    
#endif

    implicit none
   
    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state), intent(in)    :: state                    ! Physics state variables  			[vary]
    type(cam_in_t),      intent(in)    :: cam_in
    real(r8), 		 intent(in)    :: hdtime                   ! Host model timestep			[s]
    real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection	[kg/ks/s]
    real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c		[kg/m2/s]
    real(r8),            intent(in)    :: cmfmc2(pcols,pverp)      ! shallow convective mass flux--m subc 	[kg/m2/s]
    real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography			[m]
    
    ! ---------------------- !
    ! Input-Output Auguments !
    ! ---------------------- !
    
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! ---------------------- !
    ! Output Auguments !
    ! ---------------------- !

    type(physics_ptend), intent(out)   :: ptend_all 		           ! package tendencies

    ! These two variables are needed for energy check    
    real(r8),            intent(out)   :: det_s(pcols)               ! Integral of detrained static energy from ice
    real(r8),            intent(out)   :: det_ice(pcols)             ! Integral of detrained ice for energy check
        
    ! --------------- !
    ! Local Variables !
    ! --------------- !

#if CLUBB_SGS

   type(physics_state) :: state1         	! Local copy of state variable
   type(physics_ptend) :: ptend_loc      	! Local tendency from processes, added up to return as ptend_all
   
   integer :: i, j, k, t, ixind, nadv
   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
   integer :: itim
   integer :: ncol, lchnk                   	! # of columns, and chunk identifier
   integer :: err_code				! Diagnostic, for if some calculation goes amiss.
   integer :: begin_height, end_height
   integer :: icnt
   
   real(r8) :: frac_limit, ic_limit

   real(r8) :: dtime				! CLUBB time step                               [s]   
   real(r8) :: edsclr_in(pverp,edsclr_dim)      ! Scalars to be diffused through CLUBB 		[units vary]
   real(r8) :: wp2_in(pverp)			! vertical velocity variance (CLUBB)		[m^2/s^2]
   real(r8) :: wp3_in(pverp)			! third moment vertical velocity		[m^3/s^3]
   real(r8) :: wpthlp_in(pverp)			! turbulent flux of thetal			[K m/s]
   real(r8) :: wprtp_in(pverp)			! turbulent flux of total water			[kg/kg m/s]
   real(r8) :: rtpthlp_in(pverp)		! covariance of thetal and qt			[kg/kg K]
   real(r8) :: rtp2_in(pverp)			! total water variance				[kg^2/k^2]
   real(r8) :: thlp2_in(pverp)			! thetal variance				[K^2]
   real(r8) :: up2_in(pverp)			! meridional wind variance			[m^2/s^2]
   real(r8) :: vp2_in(pverp)			! zonal wind variance				[m^2/s^2]
   real(r8) :: upwp_in(pverp)			! meridional wind flux 				[m^2/s^2]
   real(r8) :: vpwp_in(pverp)			! zonal wind flux				[m^2/s^2]
   real(r8) :: thlm_in(pverp)			! liquid water potential temperature (thetal)	[K]
   real(r8) :: rtm_in(pverp)			! total water mixing ratio			[kg/kg]
   real(r8) :: um_in(pverp)			! meridional wind				[m/s]
   real(r8) :: vm_in(pverp)			! zonal wind					[m/s]
   real(r8) :: rho_in(pverp)			! mid-point density				[kg/m^3]
   real(r8) :: rcm_out(pverp)			! CLUBB output of liquid water mixing ratio	[kg/kg]
   real(r8) :: wprcp_out(pverp)			! CLUBB output of flux of liquid water		[kg/kg m/s]
   real(r8) :: cloud_frac_out(pverp)		! CLUBB output of cloud fraction		[fraction]
   real(r8) :: rcm_in_layer_out(pverp)		! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg]
   real(r8) :: cloud_cover_out(pverp)		! CLUBB output of in-cloud cloud fraction	[fraction]
   real(r8) :: rho_ds_zm(pverp)			! Dry, static density on momentum levels      	[kg/m^3]
   real(r8) :: rho_ds_zt(pverp)			! Dry, static density on thermodynamic levels 	[kg/m^3]
   real(r8) :: invrs_rho_ds_zm(pverp)		! Inv. dry, static density on momentum levels 	[m^3/kg]
   real(r8) :: invrs_rho_ds_zt(pverp)		! Inv. dry, static density on thermo. levels  	[m^3/kg]
   real(r8) :: thv_ds_zm(pverp)			! Dry, base-state theta_v on momentum levels  	[K]
   real(r8) :: thv_ds_zt(pverp)			! Dry, base-state theta_v on thermo. levels   	[K]
   real(r8) :: zt_g(pverp)			! Thermodynamic grid of CLUBB		      	[m]
   real(r8) :: zi_g(pverp)			! Momentum grid of CLUBB		      	[m]
   real(r8) :: zt_out(pcols,pverp)              ! output for the thermo CLUBB grid            	[m] 
   real(r8) :: zi_out(pcols,pverp)              ! output for momentum CLUBB grid              	[m]
   real(r8) :: fcor				! Coriolis forcing 			      	[s^-1]
   real(r8) :: sfc_elevation    		! Elevation of ground			      	[m AMSL]			      	[m]
   real(r8) :: ubar				! surface wind                                  [m/s]
   real(r8) :: ustar				! surface stress				[m/s]								
   real(r8) :: z0				! roughness height				[m]
   real(r8) :: thlm_forcing(pverp)		! theta_l forcing (thermodynamic levels)        [K/s]
   real(r8) :: rtm_forcing(pverp)		! r_t forcing (thermodynamic levels)            [(kg/kg)/s]	
   real(r8) :: um_forcing(pverp)		! u wind forcing (thermodynamic levels)     	[m/s/s]
   real(r8) :: vm_forcing(pverp)		! v wind forcing (thermodynamic levels)     	[m/s/s]
   real(r8) :: wm_zm(pverp)			! w mean wind component on momentum levels  	[m/s]
   real(r8) :: wm_zt(pverp)			! w mean wind component on thermo. levels   	[m/s]
   real(r8) :: p_in_Pa(pverp)			! Air pressure (thermodynamic levels)       	[Pa]
   real(r8) :: rho_zt(pverp)                    ! Air density on thermo levels                  [kt/m^3]
   real(r8) :: rho_zm(pverp)			! Air density on momentum levels            	[kg/m^3]
   real(r8) :: exner(pverp)			! Exner function (thermodynamic levels)     	[-]
   real(r8) :: wpthlp_sfc			! w' theta_l' at surface   			[(m K)/s]
   real(r8) :: wprtp_sfc			! w' r_t' at surface       			[(kg m)/( kg s)]
   real(r8) :: upwp_sfc				! u'w' at surface          			[m^2/s^2]
   real(r8) :: vpwp_sfc				! v'w' at surface          			[m^2/s^2]   
   real(r8) :: sclrm_forcing(pverp,sclr_dim)	! Passive scalar forcing			[{units vary}/s]
   real(r8) :: wpsclrp_sfc(sclr_dim)		! Scalar flux at surface			[{units vary} m/s]
   real(r8) :: edsclrm_forcing(pverp,edsclr_dim)! Eddy passive scalar forcing			[{units vary}/s]
   real(r8) :: wpedsclrp_sfc(edsclr_dim)	! Eddy-scalar flux at surface			[{units vary} m/s]
   real(r8) :: sclrm(pverp,sclr_dim)		! Passive scalar mean (thermo. levels) 		[units vary]
   real(r8) :: wpsclrp(pverp,sclr_dim)		! w'sclr' (momentum levels)            		[{units vary} m/s]
   real(r8) :: sclrp2(pverp,sclr_dim)		! sclr'^2 (momentum levels)            		[{units vary}^2]
   real(r8) :: sclrprtp(pverp,sclr_dim)		! sclr'rt' (momentum levels)           		[{units vary} (kg/kg)]
   real(r8) :: sclrpthlp(pverp,sclr_dim)	! sclr'thlp' (momentum levels)         		[{units vary} (K)]  
   real(r8) :: bflx22                           ! Variable for buoyancy flux for pbl            [K m/s]
   real(r8) :: C_10                             ! transfer coefficient                          [-]
   real(r8) :: khzm_out(pverp)                  ! eddy diffusivity on momentum grids            [m^2/s]
   real(r8) :: khzt_out(pverp)                  ! eddy diffusivity on thermo grids              [m^2/s]
   real(r8) :: qclvar_out(pverp)                ! cloud water variance                          [kg^2/kg^2]
   real(r8) :: qclvar(pcols,pverp)              ! cloud water variance                          [kg^2/kg^2]
   real(r8) :: zo                               ! roughness height                              [m]
   real(r8) :: dz_g(pver)                       ! thickness of layer                            [m]
   real(r8) :: newfice(pcols,pver)              ! fraction of ice in cloud at CLUBB start       [-]
   real(r8) :: minqn                            ! minimum total cloud liquid + ice threshold    [kg/kg]
   real(r8) :: tempqn                           ! temporary total cloud liquid + ice            [kg/kg]
   real(r8) :: cldthresh                        ! threshold to determin cloud fraction          [kg/kg]
   
   ! Variables below are needed to compute energy integrals for conservation
   real(r8) :: ke_a(pcols), ke_b(pcols), te_a(pcols), te_b(pcols)
   real(r8) :: wv_a(pcols), wv_b(pcols), wl_b(pcols), wl_a(pcols)
   real(r8) :: se_dis, se_a(pcols), se_b(pcols), clubb_s(pver)

   real(r8) :: exner_clubb(pcols,pverp)         ! Exner function consistent with CLUBB          [-]
   real(r8) :: wpthlp_output(pcols,pverp)       ! Heat flux output variable                     [W/m2]
   real(r8) :: wprtp_output(pcols,pverp)        ! Total water flux output variable              [W/m2]
   real(r8) :: qt_output(pcols,pver)            ! Total water mixing ratio for output           [kg/kg]
   real(r8) :: thetal_output(pcols,pver)        ! Liquid water potential temperature output     [K]
   real(r8) :: sl_output(pcols,pver)            ! Liquid water static energy                    [J/kg]
   real(r8) :: ustar2(pcols)                    ! Surface stress for PBL height                 [m2/s2]
   real(r8) :: rho(pcols,pverp)     		! Midpoint density in CAM      			[kg/m^3]
   real(r8) :: thv(pcols,pver)   		! virtual potential temperature			[K]
   real(r8) :: edsclr_out(pverp,edsclr_dim)     ! Scalars to be diffused through CLUBB		[units vary]
   real(r8) :: rcm(pcols,pverp)			! CLUBB cloud water mixing ratio		[kg/kg]
   real(r8) :: cloud_frac(pcols,pverp)		! CLUBB cloud fraction				[fraction]
   real(r8) :: rcm_in_layer(pcols,pverp)	! CLUBB in-cloud liquid water mixing ratio	[kg/kg]
   real(r8) :: cloud_cover(pcols,pverp)		! CLUBB in-cloud cloud fraction			[fraction]
   real(r8) :: wprcp(pcols,pverp)		! CLUBB liquid water flux			[m/s kg/kg]
   real(r8) :: wpthvp(pcols,pverp)		! CLUBB buoyancy flux				[W/m^2]
   real(r8) :: dlf2(pcols,pver)                 ! Detraining cld H20 from shallow convection    [kg/kg/day]
   real(r8) :: eps				! Rv/Rd                                         [-]
   real(r8) :: dum1				! dummy variable                                [units vary]
   real(r8) :: obklen(pcols)			! Obukov length					[m]
   real(r8) :: kbfs(pcols)			! Kinematic Surface heat flux			[K m/s]
   real(r8) :: th(pcols,pver)			! potential temperature				[K]
   real(r8) :: dummy2(pcols)		        ! dummy variable				[units vary]
   real(r8) :: dummy3(pcols)                    ! dummy variable				[units vary]
   real(r8) :: kinheat(pcols)                   ! Kinematic Surface heat flux			[K m/s]
   real(r8) :: ksrftms(pcols)			! Turbulent mountain stress surface drag 	[kg/s/m2]
   real(r8) :: tautmsx(pcols)			! U component of turbulent mountain stress	[N/m2]
   real(r8) :: tautmsy(pcols)			! V component of turbulent mountain stress      [N/m2]
   real(r8) :: rrho 				! Inverse of air density                        [1/kg/m^3]
   real(r8) :: kinwat(pcols)                    ! Kinematic water vapor flux                    [m/s]
   
   real(kind=time_precision) :: time_elapsed    ! time keep track of stats                      [s]
   real(r8), dimension(nparams)  :: clubb_params    ! These adjustable CLUBB parameters (C1, C2 ...)
   real(r8), dimension(sclr_dim) :: sclr_tol 	! Tolerance on passive scalar 			[units vary]
   type(pdf_parameter), dimension(pverp) :: pdf_params      ! PDF parameters   			[units vary]
   character(len=200) :: temp1, sub             ! Strings needed for CLUBB output

    ! --------------- !
    ! Pointers        !
    ! --------------- !
    
   real(r8), pointer, dimension(:,:) :: wp2      ! vertical velocity variance			[m^2/s^2]
   real(r8), pointer, dimension(:,:) :: wp3      ! third moment of vertical velocity		[m^3/s^3]
   real(r8), pointer, dimension(:,:) :: wpthlp   ! turbulent flux of thetal			[m/s K]
   real(r8), pointer, dimension(:,:) :: wprtp    ! turbulent flux of moisture			[m/s kg/kg]
   real(r8), pointer, dimension(:,:) :: rtpthlp  ! covariance of thetal and qt			[kg/kg K]
   real(r8), pointer, dimension(:,:) :: rtp2     ! moisture variance				[kg^2/kg^2]
   real(r8), pointer, dimension(:,:) :: thlp2    ! temperature variance				[K^2]
   real(r8), pointer, dimension(:,:) :: up2      ! east-west wind variance			[m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vp2      ! north-south wind variance			[m^2/s^2]
   real(r8), pointer, dimension(:,:) :: upwp     ! east-west momentum flux			[m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp     ! north-south momentum flux			[m^2/s^2]
   real(r8), pointer, dimension(:,:) :: thlm     ! mean temperature				[K]
   real(r8), pointer, dimension(:,:) :: rtm      ! mean moisture mixing ratio			[kg/kg]
   real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind				[m/s]
   real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind			[m/s]
   real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction 				[fraction]
   real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction			[fraction]
   real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction			[fraction]
   real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction		[fraction]
   real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction		[fraction]
   real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC			[kg/kg]
   real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC			[kg/kg]
   real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction		[fraction]
   real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction 		[fraction]    
   real(r8), pointer, dimension(:,:) :: khzt     ! eddy diffusivity on thermo levels            [m^2/s]
   real(r8), pointer, dimension(:,:) :: khzm     ! eddy diffusivity on momentum levels          [m^2/s]
   real(r8), pointer, dimension(:,:) :: pblh     ! planetary boundary layer height              [m]
   real(r8), pointer, dimension(:,:) :: tke      ! turbulent kinetic energy                     [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: dp_icwmr ! deep convection in cloud mixing ratio        [kg/kg]
   real(r8), pointer, dimension(:,:) :: relvar   ! relative cloud water variance                [-]
   real(r8), pointer, dimension(:,:) :: accre_enhan ! accretion enhancement factor              [-]

   logical            :: lqice(pcnst)

   intrinsic :: selected_real_kind, max

#endif
   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8
#if CLUBB_SGS

   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !       MAIN COMPUTATION BEGINS HERE               						   !
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!
   !-----------------------------------------------------------------------------------------------!

   frac_limit = 0.01_r8
   ic_limit   = 1.e-12_r8

 !  Get indicees for cloud and ice mass and cloud and ice number

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

 !  Initialize physics tendency arrays, copy the state to state1 array to use in this routine
    
   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq)

   call physics_state_copy(state,state1)

   !  Determine number of columns and which chunk computation is to be performed on

   ncol = state%ncol
   lchnk = state%lchnk

   !  Determine time step of physics buffer
   
   itim = pbuf_old_tim_idx() 

   !  Establish associations between pointers and physics buffer fields   

   call pbuf_get_field(pbuf, wp2_idx,     wp2,     start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wp3_idx,     wp3,     start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wpthlp_idx,  wpthlp,  start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wprtp_idx,   wprtp,   start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtpthlp_idx, rtpthlp, start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtp2_idx,    rtp2,    start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlp2_idx,   thlp2,   start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, up2_idx,     up2,     start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vp2_idx,     vp2,     start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, upwp_idx,    upwp,    start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vpwp_idx,    vpwp,    start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlm_idx,    thlm,    start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtm_idx,     rtm,     start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, um_idx,      um,      start=(/1,1,itim/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vm_idx,      vm,      start=(/1,1,itim/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, tke_idx,     tke)

   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
   call pbuf_get_field(pbuf, relvar_idx,  relvar)
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
   call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
   call pbuf_get_field(pbuf, kvm_idx,     khzt)
   call pbuf_get_field(pbuf, kvh_idx,     khzm)
   call pbuf_get_field(pbuf, pblh_idx,    pblh)
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)

   !  Determine CLUBB time step based on host model time step 
   !  Current algorithm is to always allow at least 4 CLUBB timesteps per host model 
   !  timestep.  However, a maximum timestep of 300 s and a minimum timestep of 60 s 
   !  is imposed.  
   dtime=max(min((1.0_r8*hdtime)/4.0_r8,300.0_r8),60.0_r8)   

   !  Initialize forcings for transported scalars to zero
   
   sclrm_forcing(:,:)   = 0._r8
   edsclrm_forcing(:,:) = 0._r8
   sclrm(:,:)           = 0._r8
   
   minqn = 0._r8
   newfice(:,:) = 0._r8
   where(state1%q(:ncol,:pver,3) .gt. minqn) &
       newfice(:ncol,:pver) = state1%q(:ncol,:pver,3)/(state1%q(:ncol,:pver,2)+state1%q(:ncol,:pver,3))

   !  determine number of timesteps CLUBB core should be advanced, 
   !  host time step divided by CLUBB time step  
   nadv = max(hdtime/dtime,1._r8)

   !  Compute exner function consistent with CLUBB's definition, which uses a constant
   !  surface pressure.  CAM's exner (in state does not).  Therefore, for consistent 
   !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables 
   !  (such as thlm), use "exner_clubb" other wise use the exner in state

   do k=1,pver
     do i=1,ncol
       exner_clubb(i,k) = 1._r8/((state1%pmid(i,k)/p0_clubb)**(rair/cpair))
     enddo
   enddo
   
   !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state 
   !  from the CAM state 

   do k=1,pver   ! loop over levels
     do i=1,ncol ! loop over columns

       rtm(i,k)  = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)
       um(i,k)   = state1%u(i,k)
       vm(i,k)   = state1%v(i,k)
       thlm(i,k) = state1%t(i,k)*exner_clubb(i,k)-(latvap/cpair)*state1%q(i,k,ixcldliq)
	 
     enddo
   enddo

   rtm(1:ncol,pverp)  = rtm(1:ncol,pver)
   um(1:ncol,pverp)   = state1%u(1:ncol,pver)
   vm(1:ncol,pverp)   = state1%v(1:ncol,pver)
   thlm(1:ncol,pverp) = thlm(1:ncol,pver)

   ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
   ! for the computation of total energy before CLUBB is called.  This is for an 
   ! effort to conserve energy since liquid water potential temperature (which CLUBB 
   ! conserves) and static energy (which CAM conserves) are not exactly equal.   
   se_b = 0._r8
   ke_b = 0._r8
   wv_b = 0._r8
   wl_b = 0._r8
   do k=1,pver
     do i=1,ncol
       se_b(i) = se_b(i) + state1%s(i,k)*state1%pdel(i,k)/gravit
       ke_b(i) = ke_b(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)/gravit
       wv_b(i) = wv_b(i) + state1%q(i,k,ixq)*state1%pdel(i,k)/gravit
       wl_b(i) = wl_b(i) + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     enddo
   enddo

   !  Compute virtual potential temperature, which is needed for CLUBB  
   do k=1,pver
     do i=1,ncol	 
       thv(i,k) = state1%t(i,k)*exner_clubb(i,k)*(1._r8+zvir*state1%q(i,k,ixq)&
                  -state1%q(i,k,ixcldliq))
     enddo
   enddo
   
   ! ------------------------------------------------- !
   ! Begin module to compute turbulent mountain stress	!
   ! ------------------------------------------------- !
   
   call compute_tms( pcols,        pver,      ncol,                   &
                     state1%u,     state1%v,  state1%t,  state1%pmid, &
	             state1%exner, state1%zm, sgh30,     ksrftms,     &
		     tautmsx,      tautmsy,   cam_in%landfrac ) 	
		     
   ! ------------------------------------------------- !
   ! End module to compute turbulent mountain stress	!
   ! ------------------------------------------------- !				       

   !  Loop over all columns in lchnk to advance CLUBB core
   do i=1,ncol   ! loop over columns

      !  Set time_elapsed to host model time step, this is for 
      !  CLUBB's budget stats
      time_elapsed = hdtime

      !  Determine Coriolis force at given latitude.  This is never used
      !  when CLUBB is implemented in a host model, therefore just set
      !  to zero.
      fcor = 0._r8 

      !  Define the CLUBB momentum grid (in height, units of m)
      do k=1,pverp
         zi_g(k) = state1%zi(i,pverp-k+1)-state1%zi(i,pver+1)
      enddo 

      !  Define the CLUBB thermodynamic grid (in units of m)
      do k=1,pver
         zt_g(k+1) = state1%zm(i,pver-k+1)-state1%zi(i,pver+1)
	 dz_g(k) = state1%zi(i,k)-state1%zi(i,k+1)  ! compute thickness
      enddo
 
      !  Thermodynamic ghost point is below surface 
      zt_g(1) = -1._r8*zt_g(2)

      !  Set the elevation of the surface
      sfc_elevation = state1%zi(i,pver+1)

      !  Compute thermodynamic stuff needed for CLUBB on thermo levels.  
      !  Inputs for the momentum levels are set below setup_clubb core
      do k=1,pver
         p_in_Pa(k+1)         = state1%pmid(i,pver-k+1)                              ! Pressure profile
         exner(k+1)           = 1._r8/exner_clubb(i,pver-k+1)
         rho_ds_zt(k+1)       = (1._r8/gravit)*(state1%pdel(i,pver-k+1)/dz_g(pver-k+1))
         invrs_rho_ds_zt(k+1) = 1._r8/(rho_ds_zt(k+1))                               ! Inverse ds rho at thermo
         rho(i,k+1)           = rho_ds_zt(k+1)                                       ! rho on thermo 
         thv_ds_zt(k+1)       = thv(i,pver-k+1)                                      ! thetav on thermo
      enddo

      !  Below computes the same stuff for the ghost point.  May or may
      !  not be needed, just to be safe to avoid NaN's
      rho_ds_zt(1)       = rho_ds_zt(2)
      invrs_rho_ds_zt(1) = invrs_rho_ds_zt(2)
      rho(i,1)           = rho_ds_zt(2)
      thv_ds_zt(1)       = thv_ds_zt(2)
      rho_zt(:)          = rho(i,:)
      p_in_Pa(1)         = p_in_Pa(2)
      exner(1)           = exner(2)

      !  Compute mean w wind on thermo grid, convert from omega to w 
      wm_zt(1) = 0._r8
      do k=1,pver
         wm_zt(k+1) = -1._r8*state1%omega(i,pver-k+1)/(rho(i,k+1)*gravit)
      enddo
    
      !  Surface fluxes provided by host model
      wpthlp_sfc = cam_in%shf(i)/(cpair*rho(i,1))       ! Sensible heat flux
      wprtp_sfc  = cam_in%lhf(i)/(latvap*rho(i,1)) 	! Latent heat flux
      upwp_sfc   = cam_in%wsx(i)/rho(i,1)               ! Surface meridional momentum flux
      vpwp_sfc   = cam_in%wsy(i)/rho(i,1)               ! Surface zonal momentum flux

     ! ------------------------------------------------- !
     ! Begin case specific code for SCAM cases.	        !
     ! This section of code block NOT called in          !
     ! global simulations				!
     ! ------------------------------------------------- !

     if (single_column) then

        !  Initialize zo if variable ustar is used

        if (cam_in%landfrac(i) .ge. 0.5_r8) then
           zo = 0.035_r8
        else
           zo = 0.0001_r8
        endif

        !  Compute surface wind (ubar)
        ubar = sqrt(um(i,pver)**2+vm(i,pver)**2)
        if (ubar .lt. 0.25_r8) ubar = 0.25_r8
    
        !  Below denotes case specifics for surface momentum
        !  and thermodynamic fluxes, depending on the case

        !  Define ustar (based on case, if not variable)     
        ustar = 0.25_r8   ! Initialize ustar in case no case
    
        if(trim(scm_clubb_iop_name) .eq. 'DYCOMSrf01_4day') then
           ustar = 0.25_r8
        endif
    
        if(trim(scm_clubb_iop_name) .eq. 'DYCOMSrf02_06hr') then
           ustar = 0.25_r8
        endif
    
        if(trim(scm_clubb_iop_name) .eq. 'BOMEX_5day') then
           ustar = 0.28_r8
        endif
    
        if(trim(scm_clubb_iop_name) .eq. 'ATEX_48hr') then
           ustar      = 0.30_r8
           wpthlp_sfc = 3.0_r8/(cpair*rho(i,1)) 
           wprtp_sfc  = 110.0_r8/(latvap*rho(i,1))
        endif
    
        if(trim(scm_clubb_iop_name) .eq. 'RICO_3day') then
           ustar      = 0.28_r8
           wpthlp_sfc = 9.5_r8/(cpair*rho(i,1))  
           wprtp_sfc  = 138.0_r8/(latvap*rho(i,1))
        endif

        if(trim(scm_clubb_iop_name) .eq. 'arm97' .or. trim(scm_clubb_iop_name) .eq. 'gate' .or. &
           trim(scm_clubb_iop_name) .eq. 'toga' .or. trim(scm_clubb_iop_name) .eq. 'mpace' .or. &
           trim(scm_clubb_iop_name) .eq. 'ARM_CC') then
       
             bflx22 = (gravit/theta0)*wpthlp_sfc
             ustar  = diag_ustar(zt_g(2),bflx22,ubar,zo)      
        endif
    
        if(trim(scm_clubb_iop_name) .eq. 'gate') then
           C_10       = 0.0013_r8
           wpthlp_sfc = -C_10*ubar*(thlm(1,pver)-300.5_r8*(1000._r8/1015._r8)**(rair/cpair))
           wprtp_sfc  = -C_10*ubar*(rtm(1,pver)-0.0198293_r8)
        endif
    
        !  Compute the surface momentum fluxes, if this is a SCAM simulation       
        upwp_sfc = -um(i,pver)*ustar**2/ubar
        vpwp_sfc = -vm(i,pver)*ustar**2/ubar
    
     endif

     ! ------------------------------------------------- !
     ! End case specific code for SCAM cases	         !
     ! ------------------------------------------------- !
    
     ! ------------------------------------------------- !
     ! Apply TMS	                                 !
     ! ------------------------------------------------- !    
      
     upwp_sfc = upwp_sfc-((ksrftms(i)*state1%u(i,pver))/rho(i,1))
     vpwp_sfc = vpwp_sfc-((ksrftms(i)*state1%v(i,pver))/rho(i,1))    

     !  Define surface sources for transported variables for diffusion, will 
     !  be zero as these tendencies are done in clubb_surface
     do ixind=1,edsclr_dim
        wpedsclrp_sfc(ixind) = 0._r8
     enddo 

     !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
     !  forcings already applied through CAM
     thlm_forcing(1:pverp) = 0._r8
     rtm_forcing(1:pverp)  = 0._r8
     um_forcing(1:pverp)   = 0._r8
     vm_forcing(1:pverp)   = 0._r8

     !  Set stats output and increment equal to CLUBB and host dt
     stats_tsamp = dtime
     stats_tout  = hdtime

     !  Heights need to be set at each timestep.  Therefore, recall 
     !  setup_grid and setup_parameters for this.  
    
     !  Read in parameters for CLUBB.  Just read in default values 
     call read_parameters( -99, "", clubb_params )

     !  Set-up CLUBB core at each CLUBB call because heights can change 	   
     call setup_grid(pverp, sfc_elevation, l_implemented, grid_type, &
       zi_g(2), zi_g(1), zi_g(pverp), zi_g(1:pverp), zt_g(1:pverp), &
       begin_height, end_height)
	
     call setup_parameters(zi_g(2), clubb_params, pverp, grid_type, &
       zi_g(begin_height:end_height), zt_g(begin_height:end_height), err_code)
	
     !  Compute some inputs from the thermodynamic grid
     !  to the momentum grid	   
     rho_ds_zm       = zt2zm(rho_ds_zt)	   
     rho_zm          = zt2zm(rho_zt)
     invrs_rho_ds_zm = zt2zm(invrs_rho_ds_zt)
     thv_ds_zm       = zt2zm(thv_ds_zt)
     wm_zm           = zt2zm(wm_zt)
 
     !  Need to flip arrays around for CLUBB core	   
     do k=1,pverp
        um_in(k)      = um(i,pverp-k+1)
        vm_in(k)      = vm(i,pverp-k+1)
        upwp_in(k)    = upwp(i,pverp-k+1)
        vpwp_in(k)    = vpwp(i,pverp-k+1)
        up2_in(k)     = up2(i,pverp-k+1)
        vp2_in(k)     = vp2(i,pverp-k+1)
        wp2_in(k)     = wp2(i,pverp-k+1)
        wp3_in(k)     = wp3(i,pverp-k+1)
        rtp2_in(k)    = rtp2(i,pverp-k+1)
        thlp2_in(k)   = thlp2(i,pverp-k+1)
        thlm_in(k)    = thlm(i,pverp-k+1)
        rtm_in(k)     = rtm(i,pverp-k+1)
        wprtp_in(k)   = wprtp(i,pverp-k+1)
        wpthlp_in(k)  = wpthlp(i,pverp-k+1)
        rtpthlp_in(k) = rtpthlp(i,pverp-k+1)
	     
        !  Initialize these to prevent crashing behavior
        rcm_out(k)          = 0._r8
        wprcp_out(k)        = 0._r8
        cloud_frac_out(k)   = 0._r8
        rcm_in_layer_out(k) = 0._r8
        cloud_cover_out(k)  = 0._r8
        edsclr_in(k,:)      = 0._r8
        edsclr_out(k,:)     = 0._r8
        khzm_out(k)         = 0._r8
        khzt_out(k)         = 0._r8
	     
        !  higher order scalar stuff, put to zero
        sclrm(k,:)          = 0._r8
        wpsclrp(k,:)        = 0._r8
        sclrp2(k,:)         = 0._r8
        sclrprtp(k,:)       = 0._r8
        sclrpthlp(k,:)      = 0._r8
        wpsclrp_sfc(:)      = 0._r8
	     	     
     enddo

     !  Do the same for tracers 
     icnt=0
     do ixind=1,pcnst
        if (lq(ixind))  then 
           icnt=icnt+1
           do k=1,pver
              edsclr_in(k+1,icnt) = state1%q(i,pver-k+1,ixind)
           enddo
           edsclr_in(1,icnt) = edsclr_in(2,icnt)
        end if
     enddo
	    
     rho_in(:) = rho(i,:)	   

     do t=1,nadv    ! do needed number of "sub" timesteps for each CAM step
    
       !  Increment the statistics then being stats timestep
       if (l_stats) then
          time_elapsed = time_elapsed+dtime
          call stats_begin_timestep(time_elapsed)
       endif 

       !  Advance CLUBB CORE one timestep in the future
       call advance_clubb_core &
          ( l_implemented, dtime, fcor, sfc_elevation, &
          thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
          sclrm_forcing, edsclrm_forcing,&      
          wm_zm, wm_zt, &
          wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
          wpsclrp_sfc, wpedsclrp_sfc, &       
          p_in_Pa, rho_zm, rho_in, exner, &
          rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
          invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
          um_in, vm_in, upwp_in, &
	  vpwp_in, up2_in, vp2_in, &
          thlm_in, rtm_in, wprtp_in, wpthlp_in, &
          wp2_in, wp3_in, rtp2_in, &
	  thlp2_in, rtpthlp_in, &
          sclrm, sclrp2, sclrprtp, sclrpthlp, &        
          wpsclrp, edsclr_in, err_code, &		   
          rcm_out, wprcp_out, cloud_frac_out, &
          rcm_in_layer_out, cloud_cover_out, &
          khzm_out, khzt_out, qclvar_out, &
          pdf_params)	
	  
       !  Check to see if stats should be output, here stats are read into
       !  output arrays to make them conformable to CAM output
       if (l_stats) call stats_end_timestep_clubb(lchnk,i,out_zt,out_zm,&
                                out_radzt,out_radzm,out_sfc)      

     enddo  ! end time loop
    
     call cleanup_grid()
 
     !  Arrays need to be "flipped" to CAM grid 
     do k=1,pverp
     
        um(i,k)           = um_in(pverp-k+1)
        vm(i,k)           = vm_in(pverp-k+1)
        upwp(i,k)         = upwp_in(pverp-k+1)
        vpwp(i,k)         = vpwp_in(pverp-k+1)
        up2(i,k)          = up2_in(pverp-k+1)
        vp2(i,k)          = vp2_in(pverp-k+1)
        thlm(i,k)         = thlm_in(pverp-k+1)
        rtm(i,k)          = rtm_in(pverp-k+1)
        wprtp(i,k)        = wprtp_in(pverp-k+1)
        wpthlp(i,k)       = wpthlp_in(pverp-k+1)
        wp2(i,k)          = wp2_in(pverp-k+1)
        wp3(i,k)          = wp3_in(pverp-k+1)
        rtp2(i,k)         = rtp2_in(pverp-k+1)
        thlp2(i,k)        = thlp2_in(pverp-k+1)
        rtpthlp(i,k)      = rtpthlp_in(pverp-k+1)
        rcm(i,k)          = rcm_out(pverp-k+1)
        wprcp(i,k)        = wprcp_out(pverp-k+1)
        cloud_frac(i,k)   = min(cloud_frac_out(pverp-k+1),1._r8)
        rcm_in_layer(i,k) = rcm_in_layer_out(pverp-k+1)
        cloud_cover(i,k)  = min(cloud_cover_out(pverp-k+1),1._r8)
        zt_out(i,k)       = zt_g(pverp-k+1)
        zi_out(i,k)       = zi_g(pverp-k+1)
        khzm(i,k)         = khzm_out(pverp-k+1)
        khzt(i,k)         = khzt_out(pverp-k+1)
        qclvar(i,k)       = min(1._r8,qclvar_out(pverp-k+1))
     
        do ixind=1,edsclr_dim
            edsclr_out(k,ixind) = edsclr_in(pverp-k+1,ixind)
        enddo

     enddo 

     zi_out(i,1) = 0._r8
     
     ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
     ! after CLUBB is called.  This is for energy conservation purposes.
     se_a = 0._r8
     ke_a = 0._r8
     wv_a = 0._r8
     wl_a = 0._r8
     do k=1,pver
       clubb_s(k) = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner_clubb(i,k))+ &
                    gravit*state1%zm(i,k)+state1%phis(i)
       se_a(i) = se_a(i) + clubb_s(k)*state1%pdel(i,k)/gravit
       ke_a(i) = ke_a(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)/gravit
       wv_a(i) = wv_a(i) + (rtm(i,k)-rcm(i,k))*state1%pdel(i,k)/gravit
       wl_a(i) = wl_a(i) + (rcm(i,k))*state1%pdel(i,k)/gravit
     enddo
     
     ! Based on these integrals, compute the total energy before and after CLUBB call
     do k=1,pver
       te_a(i) = se_a(i) + ke_a(i) + (latvap+latice)*wv_a(i)+latice*wl_a(i)
       te_b(i) = se_b(i) + ke_b(i) + (latvap+latice)*wv_b(i)+latice*wl_b(i)
     enddo
     
     ! Take into account the surface fluxes of heat and moisture
     te_b(i) = te_b(i)+(cam_in%shf(i)+(cam_in%lhf(i)/latvap)*(latvap+latice))*hdtime
     
     ! Compute the disbalance of total energy
     se_dis = (te_a(i) - te_b(i))/(state1%ps(i)-state1%pint(i,1))
     
     ! Fix the total energy coming out of CLUBB so it achieves enery conservation.
     ! Apply this fixer throughout the column evenly.
     do k=1,pver
       clubb_s(k) = clubb_s(k) - se_dis*gravit
     enddo    

     !  Now compute the tendencies of CLUBB to CAM, note that pverp is the ghost point
     !  for all variables and therefore is never called in this loop
     do k=1,pver
  
        ptend_loc%u(i,k)   = (um(i,k)-state1%u(i,k))/hdtime             ! east-west wind
        ptend_loc%v(i,k)   = (vm(i,k)-state1%v(i,k))/hdtime             ! north-south wind
        ptend_loc%q(i,k,ixq) = (rtm(i,k)-rcm(i,k)-state1%q(i,k,ixq))/hdtime ! water vapor
        ptend_loc%q(i,k,ixcldliq) = (rcm(i,k)-state1%q(i,k,ixcldliq))/hdtime   ! Tendency of liquid water
        ptend_loc%s(i,k)   = (clubb_s(k)-state1%s(i,k))/hdtime          ! Tendency of static energy

       !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
       !  Loading up this array doesn't mean the tendencies are applied.  
       ! edsclr_out is compressed with just the constituents being used, ptend and state are not compressed

        icnt=0
        do ixind=1,pcnst
           if (lq(ixind)) then
              icnt=icnt+1
              if ((ixind /= ixq) .and. (ixind /= ixcldliq)) then
                 ptend_loc%q(i,k,ixind) = (edsclr_out(k,icnt)-state1%q(i,k,ixind))/hdtime ! transported constituents 
              end if
           end if
        enddo

     enddo

   enddo  ! end column loop
 
    ! ------------------------------------------------- !
    ! End column computation of CLUBB, begin to apply   !
    ! and compute output, etc		 	        !
    ! ------------------------------------------------- !

   !  Output CLUBB tendencies 
   call outfld( 'RVMTEND_CLUBB', ptend_loc%q(:,:,ixq)*1000._r8, pcols, lchnk)
   call outfld( 'RCMTEND_CLUBB', ptend_loc%q(:,:,ixcldliq)*1000._r8, pcols, lchnk)
   call outfld( 'RIMTEND_CLUBB', ptend_loc%q(:,:,ixcldice)*1000._r8, pcols, lchnk)
   call outfld( 'STEND_CLUBB',   ptend_loc%s,pcols, lchnk)
   call outfld( 'UTEND_CLUBB',   ptend_loc%u,pcols, lchnk)
   call outfld( 'VTEND_CLUBB',   ptend_loc%v,pcols, lchnk)

   !  Update physics tendencies     
   call physics_ptend_init(ptend_all, state%psetcols, 'clubb')
   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)
   
   ! --------------------------------------------------------------------------------- !  
   !  COMPUTE THE ICE CLOUD DETRAINMENT				                       !
   !  Detrainment of convective condensate into the environment or stratiform cloud    !
   ! --------------------------------------------------------------------------------- ! 

   !  Initialize the shallow convective detrainment rate, will always be zero
   dlf2(:,:) = 0.0_r8

   lqice(:)        = .false.
   lqice(ixcldliq) = .true.
   lqice(ixcldice) = .true.
   lqice(ixnumliq) = .true.
   lqice(ixnumice) = .true.
    
   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb', ls=.true., lq=lqice)
   
   do k=1,pver
      do i=1,ncol
        if( state1%t(i,k) > 268.15_r8 ) then
	  dum1 = 0.0_r8
	elseif ( state1%t(i,k) < 238.15_r8 ) then
	  dum1 = 1.0_r8
	else
	  dum1 = ( 268.15_r8 - state1%t(i,k) ) / 30._r8 
	endif
        
	ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
        ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
        ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
	                          / (4._r8*3.14_r8* 8.e-6_r8**3*997._r8) + & ! Deep    Convection
                                  3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
				  / (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection 
        ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
	                          / (4._r8*3.14_r8*25.e-6_r8**3*500._r8) + & ! Deep    Convection
                                  3._r8 * (                         dlf2(i,k)    *  dum1 ) &
				  / (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
        ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

        ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
	!   track of the integrals of ice and static energy that is effected from conversion to ice
	!   so that the energy checker doesn't complain.	
	det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)/gravit
	det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)/gravit
	
      enddo
   enddo
   
   det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water
  
   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)

   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ ! 
   ! ------------------------------------------------------------ !
   ! The rest of the code deals with diagnosing variables         !
   ! for microphysics/radiation computation and macrophysics      !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !

   ! ------------------------------------------------- !
   ! Diagnose relative cloud water variance	       !
   ! ------------------------------------------------- !

   
   relvar(:,:) = 1.0_r8  ! default
   
   where (rcm(:ncol,:pver) /= 0 .and. qclvar(:ncol,:pver) /= 0) &
         relvar(:ncol,:pver) = min(1.0_r8,max(0.001_r8,rcm(:ncol,:pver)**2/qclvar(:ncol,:pver)))
   
   ! ------------------------------------------------- !
   ! Optional Accretion enhancement factor	       !
   ! ------------------------------------------------- !   
   
   accre_enhan(:ncol,:pver) = 1._r8+0.65_r8*(1.0_r8/relvar(:ncol,:pver))
   
   ! ------------------------------------------------- !
   ! Diagnose some output variables		       !
   ! ------------------------------------------------- !

   !  density
   rho(:ncol,1:pver) = state1%pmid(:ncol,1:pver)/(rair*state1%t(:ncol,1:pver))
   rho(:ncol,pverp)  = state1%ps(:ncol)/(rair*state1%t(:ncol,pver))

   eps = rair/rh2o
   wpthvp(:,:) = 0.0_r8
   do k=1,pver
     do i=1,ncol
        !  buoyancy flux
        wpthvp(i,k) = wpthlp(i,k)+((1._r8-eps)/eps)*theta0*wprtp(i,k)+((latvap/cpair)* &
          state1%exner(i,k)-(1._r8/eps)*theta0)*wprcp(i,k)

        !  total water mixing ratio
        qt_output(i,k) = state1%q(i,k,ixq)+state1%q(i,k,ixcldliq)+state1%q(i,k,ixcldice)
        !  liquid water potential temperature
        thetal_output(i,k) = (state1%t(i,k)*state1%exner(i,k))-(latvap/cpair)*state1%q(i,k,ixcldliq)
        !  liquid water static energy
        sl_output(i,k) = cpair*state1%t(i,k)+gravit*state1%zm(i,k)-latvap*state1%q(i,k,ixcldliq)
     enddo
   enddo
   
   do k=1,pverp
     do i=1,ncol
        !  liquid water potential temperature flux
        wpthlp_output(i,k) = wpthlp(i,k)*rho(i,k)*cpair
        !  total water mixig ratio flux
        wprtp_output(i,k)  = wprtp(i,k)*rho(i,k)*latvap
	!  turbulent kinetic energy
	tke(i,k) = 0.5_r8*(up2(i,k)+vp2(i,k)+wp2(i,k))
     enddo
   enddo
   
   ! --------------------------------------------------------------------------------- ! 
   !  Diagnose some quantities that are computed in macrop_tend here.                  !
   !  These are inputs required for the microphysics calculation.                      !
   !                              						       !	
   !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
   ! --------------------------------------------------------------------------------- ! 
 
   !  initialize variables 
   alst(:,:) = 0.0_r8
   qlst(:,:) = 0.0_r8 
 
   do k=1,pver
     do i=1,ncol
       alst(i,k) = cloud_frac(i,k)   
       qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
     enddo
   enddo
 
   ! --------------------------------------------------------------------------------- !  
   !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION 		       !
   ! --------------------------------------------------------------------------------- ! 
 
   deepcu(:,pver) = 0.0_r8
   shalcu(:,pver) = 0.0_r8
 
   do k=1,pver-1
     do i=1,ncol
       !  diagnose the deep convective cloud fraction, as done in macrophysics based on the 
       !  deep convective mass flux, read in from pbuf.  Since shallow convection is never 
       !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
       !  fraction is purely from deep convection scheme.  
       deepcu(i,k) = max(0.0_r8,min(0.1_r8*log(1.0_r8+500.0_r8*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.6_r8))
       shalcu(i,k) = 0._r8
       
       if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
         deepcu(i,k) = 0._r8
       endif
             
       !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable 
       !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
       !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud 
       !  from CLUBB plus the deep convective cloud fraction
       concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+deepcu(i,k),0.80_r8)
     enddo
   enddo
   
   if (single_column) then
     if (trim(scm_clubb_iop_name) .eq. 'ATEX_48hr' .or. trim(scm_clubb_iop_name) .eq. 'BOMEX_5day' .or. &
         trim(scm_clubb_iop_name) .eq. 'DYCOMSrf01_4day' .or. &
	 trim(scm_clubb_iop_name) .eq. 'DYCOMSrf02_06hr' .or. &
	 trim(scm_clubb_iop_name) .eq. 'RICO_3day' .or. &
	 trim(scm_clubb_iop_name) .eq. 'ARM_CC') then
       
           deepcu(:,:) = 0.0_r8
           concld(:,:) = 0.0_r8
       
     endif       
   endif
   
   ! --------------------------------------------------------------------------------- !  
   !  COMPUTE THE ICE CLOUD FRACTION PORTION					       !
   !  use the aist_vector function to compute the ice cloud fraction                   !
   ! --------------------------------------------------------------------------------- !  

    do k=1,pver
      call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
           cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol)
    enddo
  
   ! --------------------------------------------------------------------------------- !  
   !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION				       !
   !										       !
   !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
   !  because CLUBB does nothing with ice.  Here I simply overwrite the liquid stratus ! 
   !  fraction that was coded in macrop_driver 					       !	
   ! --------------------------------------------------------------------------------- !  
 
   !  Recompute net stratus fraction using maximum over-lapping assumption, as done
   !  in macrophysics code, using alst computed above and aist read in from physics buffer
   
   cldthresh=1.e-18_r8

   do k=1,pver
     do i=1,ncol
       ast(i,k) = 0._r8  ! init AST

       if (newfice(i,k) .le. 0.5_r8 .and. state1%q(i,k,2) .gt. cldthresh) then
          ast(i,k) = alst(i,k)
       else if ((newfice(i,k) .gt. 0.5_r8 .and. state1%q(i,k,3) .gt. cldthresh) .or. &
                (newfice(i,k) .le. 0.5_r8 .and. state1%q(i,k,2) .lt. cldthresh &
                                          .and. state1%q(i,k,3) .gt. cldthresh)) then
          ast(i,k) = aist(i,k)
       end if
	 
       qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k)) 
     enddo
   enddo
   
   !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just 
   !  be outputting the shallow convective cloud fraction 

   do k=1,pver
     do i=1,ncol
       cloud_frac(i,k) = min(ast(i,k)+deepcu(i,k),1.0_r8)
     enddo
   enddo
   
   ! --------------------------------------------------------------------------------- !  
   !  DIAGNOSE THE PBL DEPTH					                       !
   !  this is needed for aerosol code                                                  !
   ! --------------------------------------------------------------------------------- ! 
    
   do i=1,ncol
     do k=1,pver
       th(i,k) = state1%t(i,k)*state1%exner(i,k)
       thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq))
     enddo
   enddo
 
   ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
   do i=1,ncol
       call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                        rrho, ustar2(i) )
       call calc_obklen( th(i,pver), thv(i,pver), cam_in%lhf(i)/latvap, cam_in%shf(i), rrho, ustar2(i), &
                        kinheat(i), kinwat(i), kbfs(i), obklen(i) )  		
   enddo
   
   dummy2(:) = 0._r8
   dummy3(:) = 0._r8

   !  Compute PBL depth according to Holtslag-Boville Scheme
   call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kinheat, pblh, dummy2, &
        	state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)
		 
   !  Output the PBL depth
   call outfld('PBLH', pblh, pcols, lchnk)
 
   ! Assign the first pver levels of cloud_frac back to cld
   cld(:,1:pver) = cloud_frac(:,1:pver)

   ! --------------------------------------------------------------------------------- !   
   !  END CLOUD FRACTION DIAGNOSIS, begin to store variables back into buffer          !
   ! --------------------------------------------------------------------------------- !  
 
   !  Output calls of variables goes here
   call outfld( 'RHO_CLUBB',        rho,                   pcols, lchnk )
   call outfld( 'WP2_CLUBB',        wp2,                   pcols, lchnk )
   call outfld( 'UP2_CLUBB',        up2,                   pcols, lchnk )
   call outfld( 'VP2_CLUBB',        vp2,                   pcols, lchnk )
   call outfld( 'WP3_CLUBB',        wp3,                   pcols, lchnk )
   call outfld( 'UPWP_CLUBB',       upwp,                  pcols, lchnk )
   call outfld( 'VPWP_CLUBB',       vpwp,                  pcols, lchnk )
   call outfld( 'WPTHLP_CLUBB',     wpthlp_output,         pcols, lchnk )
   call outfld( 'WPRTP_CLUBB',      wprtp_output,          pcols, lchnk )
   call outfld( 'RTP2_CLUBB',       rtp2*1000._r8,         pcols, lchnk )
   call outfld( 'THLP2_CLUBB',      thlp2,                 pcols, lchnk )
   call outfld( 'RTPTHLP_CLUBB',    rtpthlp*1000._r8,      pcols, lchnk )
   call outfld( 'RCM_CLUBB',        rcm*1000._r8,          pcols, lchnk )
   call outfld( 'WPRCP_CLUBB',      wprcp*latvap,          pcols, lchnk )
   call outfld( 'CLOUDFRAC_CLUBB',  alst,                  pcols, lchnk )
   call outfld( 'RCMINLAYER_CLUBB', rcm_in_layer*1000._r8, pcols, lchnk )
   call outfld( 'CLOUDCOVER_CLUBB', cloud_frac,            pcols, lchnk )
   call outfld( 'WPTHVP_CLUBB',     wpthvp*cpair,          pcols, lchnk )
   call outfld( 'ZT_CLUBB',         1._r8*zt_out,          pcols, lchnk )
   call outfld( 'ZM_CLUBB',         1._r8*zi_out,          pcols, lchnk )
   call outfld( 'UM_CLUBB',         um,                    pcols, lchnk )
   call outfld( 'VM_CLUBB',         vm,                    pcols, lchnk )
   call outfld( 'THETAL',           thetal_output,         pcols, lchnk )
   call outfld( 'QT',               qt_output,             pcols, lchnk )
   call outfld( 'SL',               sl_output,             pcols, lchnk )

   !  Output CLUBB history here
   if (l_stats) then 
      
   do i=1,zt%nn
   
     temp1 = trim(zt%f%var(i)%name)
     sub   = temp1
     if (len(temp1) .gt. 16) sub = temp1(1:16)
   
     call outfld(trim(sub), out_zt(:,:,i), pcols, lchnk )
   enddo
   
   do i=1,zm%nn
   
     temp1 = trim(zm%f%var(i)%name)
     sub   = temp1
     if (len(temp1) .gt. 16) sub = temp1(1:16)
   
     call outfld(trim(sub),out_zm(:,:,i), pcols, lchnk)
   enddo

   if (l_output_rad_files) then  
     do i=1,rad_zt%nn
       call outfld(trim(rad_zt%f%var(i)%name), out_radzt(:,:,i), pcols, lchnk)
     enddo
   
     do i=1,rad_zm%nn
       call outfld(trim(rad_zm%f%var(i)%name), out_radzm(:,:,i), pcols, lchnk)
     enddo
   endif
   
   do i=1,sfc%nn
     call outfld(trim(sfc%f%var(i)%name), out_sfc(:,:,i), pcols, lchnk)
   enddo
   
   endif
   
   return
#endif
  end subroutine clubb_tend_cam
    
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
    subroutine clubb_surface ( &
    			     state, ptend, ztodt, cam_in, ustar, obklen)
    
!-------------------------------------------------------------------------------
! Description: Provide the obukov length and the surface friction velocity 
!              for the dry deposition code in routine tphysac.  Since University
!              of Washington Moist Turbulence (UWMT) scheme is not called when 
!              CLUBB is turned on the obukov length and ustar are never initialized
!              nor computed (sometimes never updated from NaN).  In addition, surface
!              fluxes are applied to the constituents.  
!   
! Author: Peter Bogenschutz, August 2011
! Origin: Based heavily on UWMT code (eddy_diff.F90)
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types,		only: physics_state, physics_ptend, physics_ptend_init
    use physconst, 		only: gravit, zvir, latvap
    use ppgrid, 	        only: pver, pcols
    use constituents, 	        only: pcnst, cnst_get_ind
    use camsrfexch,             only: cam_in_t
#ifdef MODAL_AERO
    use modal_aero_data
#endif
    
    implicit none
    
    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state), intent(in)    	:: state		! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in
    
    real(r8),  		 intent(in)     :: ztodt		! 2 delta-t        [ s ] 

    ! ---------------- !
    ! Output Auguments !
    ! ---------------- !
    
    type(physics_ptend), intent(out)    :: ptend		! Individual parameterization tendencies
    real(r8), 		 intent(out)  	:: obklen(pcols)    	! Obukhov length [ m ]
    real(r8), 		 intent(out)	:: ustar(pcols)		! Surface friction velocity [ m/s ]
    
#ifdef CLUBB_SGS 

    ! --------------- !
    ! Local Variables !
    ! --------------- !
    
    integer :: i   						! indicees					
    integer :: ncol						! # of atmospheric columns
    
    real(r8) :: th(pcols)                                       ! surface potential temperature
    real(r8) :: thv(pcols)                                      ! surface virtual potential temperature
    real(r8) :: kinheat						! kinematic surface heat flux
    real(r8) :: kinwat						! kinematic surface vapor flux
    real(r8) :: kbfs						! kinematic surface buoyancy flux
    real(r8) :: tmp1(pcols)
    real(r8) :: rztodt                                          ! 1./ztodt
    integer  :: m
    integer  :: ixq
    real(r8) :: rrho						! Inverse air density
    
    logical  :: lq(pcnst)

#endif
    obklen(pcols) = 0.0_r8
    ustar(pcols)  = 0.0_r8
#ifdef CLUBB_SGS

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !


    call cnst_get_ind('Q',ixq)
    
    lq(:) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'clubb', lq=lq)
    
    ncol = state%ncol
    
    ! Compute the surface friction velocity and obukov length
    
    do i = 1, ncol
      th(i) = state%t(i,pver)*state%exner(i,pver)         ! diagnose potential temperature
      thv(i) = th(i)*(1._r8+zvir*state%q(i,pver,ixq))  ! diagnose virtual potential temperature
    enddo
    
    do i = 1, ncol
      call calc_ustar( state%t(i,pver), state%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                       rrho, ustar(i) )
      call calc_obklen( th(i), thv(i), cam_in%lhf(i)/latvap, cam_in%shf(i), rrho, ustar(i), &
                       kinheat, kinwat, kbfs, obklen(i) )
    enddo

    rztodt = 1._r8/ztodt
    ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
    tmp1(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
    do m = 2, pcnst
      ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) * cam_in%cflx(:ncol,m)
    enddo
    
    ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) * rztodt
    
    return

#endif

    end subroutine clubb_surface

#ifdef CLUBB_SGS
! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!

real(r8) function diag_ustar( z, bflx, wnd, z0 ) 

use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

implicit none

real(r8), parameter      :: am   =  4.8_r8   !   "          "         "
real(r8), parameter      :: bm   = 19.3_r8  !   "          "         "

real(r8), parameter      :: grav = shr_const_g
real(r8), parameter      :: vonk = shr_const_karman
real(r8), parameter      :: pi   = shr_const_pi

real(r8), intent (in)    :: z             ! height where u locates
real(r8), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real(r8), intent (in)    :: wnd           ! wind speed at z
real(r8), intent (in)    :: z0            ! momentum roughness height


integer :: iterate
real(r8)    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0_r8 - 3.0_r8*log( 2.0_r8 )

ustar =  wnd*klnz
if (abs(bflx) > 1.e-6_r8) then
  do iterate=1,4

    if (ustar > 1.e-6_r8) then
        lmo   = -ustar**3 / ( vonk * bflx )
        zeta  = z/lmo
        if (zeta > 0._r8) then
          ustar =  vonk*wnd  /(lnz + am*zeta)
        else
          x     = sqrt( sqrt( 1.0_r8 - bm*zeta ) )
          psi1  = 2._r8*log( 1.0_r8+x ) + log( 1.0_r8+x*x ) - 2._r8*atan( x ) + c1
          ustar = wnd*vonk/(lnz - psi1)
        end if

    endif

  end do
end if


diag_ustar = ustar

return


end function diag_ustar
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS
			 
  subroutine stats_init_clubb( l_stats_in, stats_tsamp_in, stats_tout_in, &
                         nnzp, nnrad_zt,nnrad_zm, delt )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.   
    
    !-----------------------------------------------------------------------


    use stats_variables, only: & 
      zt,      & ! Variables
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      ztscr07, & 
      ztscr08, & 
      ztscr09, & 
      ztscr10, & 
      ztscr11, & 
      ztscr12, & 
      ztscr13, & 
      ztscr14, & 
      ztscr15, & 
      ztscr16, & 
      ztscr17, & 
      ztscr18, & 
      ztscr19, & 
      ztscr20, & 
      ztscr21

    use stats_variables, only: & 
      zm,      & 
      zmscr01, & 
      zmscr02, & 
      zmscr03, & 
      zmscr04, & 
      zmscr05, & 
      zmscr06, & 
      zmscr07, & 
      zmscr08, & 
      zmscr09, & 
      zmscr10, & 
      zmscr11, & 
      zmscr12, & 
      zmscr13, & 
      zmscr14, & 
      zmscr15, &
      zmscr16, &
      zmscr17, &
      rad_zt,  &
      rad_zm,  &
      sfc,     & 
      l_stats, &
      l_output_rad_files, & 
      stats_tsamp,   & 
      stats_tout,    & 
      l_stats_samp,  & 
      l_stats_last, & 
      fname_rad_zt, &
      fname_rad_zm, & 
      fname_sfc, & 
      l_netcdf, & 
      l_grads

    use clubb_precision, only:         time_precision   ! 
    use stats_zm, only:                nvarmax_zm, stats_init_zm ! 
    use stats_zt, only:                nvarmax_zt, stats_init_zt ! 
    use stats_rad_zt, only:            nvarmax_rad_zt, stats_init_rad_zt ! 
    use stats_rad_zm, only:            nvarmax_rad_zm, stats_init_rad_zm !        
    use stats_sfc, only:               nvarmax_sfc, stats_init_sfc ! 
    use error_code, only:              clubb_at_least_debug_level ! 
    use constants_clubb, only:         fstderr, var_length !     
    use cam_history, only:             addfld, phys_decomp
    use namelist_utils,only:           find_group_name
    use units,only:                    getunit, freeunit
    use abortutils,only:               endrun

    implicit none

    ! Input Variables

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count] 
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Timestep (dtmain in CLUBB)         [s]


    !  Local Variables

    !  Namelist Variables

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      clubb_vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_zm) ::  & 
      clubb_vars_zm  ! Variables on the momentum levels

    character(len=var_length), dimension(nvarmax_rad_zt) ::  & 
      clubb_vars_rad_zt  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_rad_zm) ::  & 
      clubb_vars_rad_zm  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_sfc) ::  &
      clubb_vars_sfc ! Variables at the model surface

    namelist /clubb_stats_nl/ & 
      clubb_vars_zt, & 
      clubb_vars_zm, &
      clubb_vars_rad_zt, &
      clubb_vars_rad_zm, & 
      clubb_vars_sfc

    !  Local Variables

    logical :: l_error

    character(len=200) :: fname, temp1, sub

    integer :: i, ntot, read_status
    integer :: iunit

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in
    
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_last  = .false.
      return
    end if

    !  Initialize namelist variables

    clubb_vars_zt     = ''
    clubb_vars_zm     = ''
    clubb_vars_rad_zt = ''
    clubb_vars_rad_zm = ''
    clubb_vars_sfc    = ''

    !  Read variables to compute from the namelist    
    if (masterproc) then
      iunit= getunit()
      open(unit=iunit,file="atm_in",status='old')
      call find_group_name(iunit, 'clubb_stats_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_stats_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_tend_cam:  error reading namelist')
         end if
      end if
      close(unit=iunit)
      call freeunit(iunit)
    end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(clubb_vars_zt,      var_length*nvarmax_zt,       mpichar,   0, mpicom)
      call mpibcast(clubb_vars_zm,      var_length*nvarmax_zm,       mpichar,   0, mpicom)
      call mpibcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpichar,   0, mpicom)
      call mpibcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpichar,   0, mpicom)
      call mpibcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpichar,   0, mpicom)
#endif

    !  Hardcode these for use in CAM-CLUBB, don't want either
    l_netcdf = .false.
    l_grads  = .false.

    !  Check sampling and output frequencies

    !  The model time step length, delt (which is dtmain), should multiply
    !  evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) )  & 
           > 1.e-8_r8 ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                       'delt (which is dtmain).  Check the appropriate ',  &
                       'model.in file.'
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
      write(fstderr,*) 'delt = ', delt
    endif

    !  Initialize zt (mass points)

    i = 1
    do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_zt(i)) /= 0 & 
               .and. i <= nvarmax_zt )
       i = i + 1
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      stop "stats_init_clubb:  number of zt statistical variables exceeds limit"
    endif

    zt%nn = ntot
    zt%kk = nnzp

    allocate( zt%z( zt%kk ) )

    allocate( zt%x( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%n( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%l_in_update( 1, 1, zt%kk, zt%nn ) )
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )

    allocate( zt%f%var( zt%nn ) )
    allocate( zt%f%z( zt%kk ) )

    !  Allocate scratch space

    allocate( ztscr01(zt%kk) )
    allocate( ztscr02(zt%kk) )
    allocate( ztscr03(zt%kk) )
    allocate( ztscr04(zt%kk) )
    allocate( ztscr05(zt%kk) )
    allocate( ztscr06(zt%kk) )
    allocate( ztscr07(zt%kk) )
    allocate( ztscr08(zt%kk) )
    allocate( ztscr09(zt%kk) )
    allocate( ztscr10(zt%kk) )
    allocate( ztscr11(zt%kk) )
    allocate( ztscr12(zt%kk) )
    allocate( ztscr13(zt%kk) )
    allocate( ztscr14(zt%kk) )
    allocate( ztscr15(zt%kk) )
    allocate( ztscr16(zt%kk) )
    allocate( ztscr17(zt%kk) )
    allocate( ztscr18(zt%kk) )
    allocate( ztscr19(zt%kk) )
    allocate( ztscr20(zt%kk) )
    allocate( ztscr21(zt%kk) )

    ztscr01 = 0.0_r8
    ztscr02 = 0.0_r8
    ztscr03 = 0.0_r8
    ztscr04 = 0.0_r8
    ztscr05 = 0.0_r8
    ztscr06 = 0.0_r8
    ztscr07 = 0.0_r8
    ztscr08 = 0.0_r8
    ztscr09 = 0.0_r8
    ztscr10 = 0.0_r8
    ztscr11 = 0.0_r8
    ztscr12 = 0.0_r8
    ztscr13 = 0.0_r8
    ztscr14 = 0.0_r8
    ztscr15 = 0.0_r8
    ztscr16 = 0.0_r8
    ztscr17 = 0.0_r8
    ztscr18 = 0.0_r8
    ztscr19 = 0.0_r8
    ztscr20 = 0.0_r8
    ztscr21 = 0.0_r8

    !  Default initialization for array indices for zt

    call stats_init_zt( clubb_vars_zt, l_error )

    !  Initialize zm (momentum points)

    i = 1
    do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_zm(i)) /= 0 & 
               .and. i <= nvarmax_zm )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      stop "stats_init_clubb:  number of zm statistical variables exceeds limit"
    endif

    zm%nn = ntot
    zm%kk = nnzp

    allocate( zm%z( zm%kk ) )

    allocate( zm%x( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%n( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%l_in_update( 1, 1, zm%kk, zm%nn ) )

    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )

    allocate( zm%f%var( zm%nn ) )
    allocate( zm%f%z( zm%kk ) )

    !  Allocate scratch space

    allocate( zmscr01(zm%kk) )
    allocate( zmscr02(zm%kk) )
    allocate( zmscr03(zm%kk) )
    allocate( zmscr04(zm%kk) )
    allocate( zmscr05(zm%kk) )
    allocate( zmscr06(zm%kk) )
    allocate( zmscr07(zm%kk) )
    allocate( zmscr08(zm%kk) )
    allocate( zmscr09(zm%kk) )
    allocate( zmscr10(zm%kk) )
    allocate( zmscr11(zm%kk) )
    allocate( zmscr12(zm%kk) )
    allocate( zmscr13(zm%kk) )
    allocate( zmscr14(zm%kk) )
    allocate( zmscr15(zm%kk) )
    allocate( zmscr16(zm%kk) )
    allocate( zmscr17(zm%kk) )

    zmscr01 = 0.0_r8
    zmscr02 = 0.0_r8
    zmscr03 = 0.0_r8
    zmscr04 = 0.0_r8
    zmscr05 = 0.0_r8
    zmscr06 = 0.0_r8
    zmscr07 = 0.0_r8
    zmscr08 = 0.0_r8
    zmscr09 = 0.0_r8
    zmscr10 = 0.0_r8
    zmscr11 = 0.0_r8
    zmscr12 = 0.0_r8
    zmscr13 = 0.0_r8
    zmscr14 = 0.0_r8
    zmscr15 = 0.0_r8
    zmscr16 = 0.0_r8
    zmscr17 = 0.0_r8

    call stats_init_zm( clubb_vars_zm, l_error )

    !  Initialize rad_zt (radiation points)

    if (l_output_rad_files) then
    
      i = 1
      do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_rad_zt(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        stop "stats_init_clubb:  number of rad_zt statistical variables exceeds limit"
      endif

      rad_zt%nn = ntot
      rad_zt%kk = nnrad_zt

      allocate( rad_zt%z( rad_zt%kk ) )

      allocate( rad_zt%x( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%n( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%l_in_update( 1, 1, rad_zt%kk, rad_zt%nn ) )

      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )

      allocate( rad_zt%f%var( rad_zt%nn ) )
      allocate( rad_zt%f%z( rad_zt%kk ) )

      fname = trim( fname_rad_zt )

      call stats_init_rad_zt( clubb_vars_rad_zt, l_error )

      !  Initialize rad_zm (radiation points)

      i = 1
      do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_rad_zm(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zm )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        stop "stats_init_clubb:  number of rad_zm statistical variables exceeds limit"
      endif

      rad_zm%nn = ntot
      rad_zm%kk = nnrad_zm

      allocate( rad_zm%z( rad_zm%kk ) )

      allocate( rad_zm%x( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%n( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%l_in_update( 1, 1, rad_zm%kk, rad_zm%nn ) )

      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )

      allocate( rad_zm%f%var( rad_zm%nn ) )
      allocate( rad_zm%f%z( rad_zm%kk ) )


      fname = trim( fname_rad_zm )

      call stats_init_rad_zm( clubb_vars_rad_zm, l_error )
    end if ! l_output_rad_files


    !  Initialize sfc (surface point)

    i = 1
    do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_sfc(i)) /= 0 & 
               .and. i <= nvarmax_sfc )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      stop "stats_init_clubb:  number of sfc statistical variables exceeds limit"
    endif

    sfc%nn = ntot
    sfc%kk = 1

    allocate( sfc%z( sfc%kk ) )

    allocate( sfc%x( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%n( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%l_in_update( 1, 1, sfc%kk, sfc%nn ) )

    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    allocate( sfc%f%var( sfc%nn ) )
    allocate( sfc%f%z( sfc%kk ) )

    fname = trim( fname_sfc )

    call stats_init_sfc( clubb_vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      stop
    endif

!   Now call add fields
    do i = 1, zt%nn
    
      temp1 = trim(zt%f%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)
     
       call addfld(trim(sub),trim(zt%f%var(i)%units),nnzp,&
            'A',trim(zt%f%var(i)%description),phys_decomp)
    enddo
    
    do i = 1, zm%nn
    
      temp1 = trim(zm%f%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)
    
      call addfld(trim(sub),trim(zm%f%var(i)%units),nnzp,&
           'A',trim(zm%f%var(i)%description),phys_decomp)
    enddo

    if (l_output_rad_files) then     
      do i = 1, rad_zt%nn
        call addfld(trim(rad_zt%f%var(i)%name),trim(rad_zt%f%var(i)%units),nnzp,&
           'A',trim(rad_zt%f%var(i)%description),phys_decomp)
      enddo
    
      do i = 1, rad_zm%nn
        call addfld(trim(rad_zm%f%var(i)%name),trim(rad_zm%f%var(i)%units),nnzp,&
           'A',trim(rad_zm%f%var(i)%description),phys_decomp)
      enddo
    endif 
    
    do i = 1, sfc%nn
      call addfld(trim(sfc%f%var(i)%name),trim(sfc%f%var(i)%units),1,&
           'A',trim(sfc%f%var(i)%description),phys_decomp)
    enddo

    return


  end subroutine stats_init_clubb 
  
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  
    !-----------------------------------------------------------------------
  subroutine stats_end_timestep_clubb(lchnk,thecol,out_zt,out_zm,out_radzt,out_radzm,out_sfc)

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

#ifdef CLUBB_SGS

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        zm, & 
        rad_zt, &
        rad_zm, &
        sfc, & 
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout, &
        l_output_rad_files

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)
	
    use cam_history, 	only: outfld
    
    use ppgrid, only: pcols, pverp

    implicit none


#endif

    integer :: lchnk
    integer :: thecol
    
    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pverp,rad_zt%nn)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%nn)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%nn)

#ifdef CLUBB_SGS
    ! Local Variables

    integer :: i, k
    logical :: l_error

    !  Check if it is time to write to file

    if ( .not. l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zt statistics at each vertical level.
    do i = 1, zt%nn
      do k = 1, zt%kk

        if ( zt%n(1,1,k,i) /= 0 .and.  &
             zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zt%f%var(i)%name), ' in zt ',  &
                             'at k = ', k,  &
                             '; zt%n(',k,',',i,') = ', zt%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zm statistics at each vertical level.
    do i = 1, zm%nn
      do k = 1, zm%kk

        if ( zm%n(1,1,k,i) /= 0 .and.  &
             zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zm%f%var(i)%name), ' in zm ',  &
                             'at k = ', k,  &
                             '; zm%n(',k,',',i,') = ', zm%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    if (l_output_rad_files) then
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zt statistics at each vertical level.
      do i = 1, rad_zt%nn
        do k = 1, rad_zt%kk

          if ( rad_zt%n(1,1,k,i) /= 0 .and.  &
               rad_zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zt%f%var(i)%name), ' in rad_zt ',  &
                               'at k = ', k,  &
                               '; rad_zt%n(',k,',',i,') = ', rad_zt%n(1,1,k,i)
            endif

          endif

        enddo
      enddo
    
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zm statistics at each vertical level.
      do i = 1, rad_zm%nn
        do k = 1, rad_zm%kk

          if ( rad_zm%n(1,1,k,i) /= 0 .and.  &
               rad_zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zm%f%var(i)%name), ' in rad_zm ',  &
                               'at k = ', k,  &
                               '; rad_zm%n(',k,',',i,') = ', rad_zm%n(1,1,k,i)
            endif

          endif

        enddo
      enddo
    end if ! l_output_rad_files

    !  Look for errors by checking the number of sampling points
    !  for each variable in the sfc statistics at each vertical level.
    do i = 1, sfc%nn
      do k = 1, sfc%kk

        if ( sfc%n(1,1,k,i) /= 0 .and.  &
             sfc%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(sfc%f%var(i)%name), ' in sfc ',  &
                             'at k = ', k,  &
                             '; sfc%n(',k,',',i,') = ', sfc%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Stop the run if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      stop 'stats_end_timestep:  error(s) found'
    endif

    !  Compute averages
    call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
    call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
    if (l_output_rad_files) then
      call stats_avg( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n )
      call stats_avg( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n )
    end if
    call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

   !  Here we are not outputting the data, rather reading the stats into 
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.    
    do i = 1, zt%nn
      do k = 1, zt%kk 
         out_zt(thecol,k,i) = zt%x(1,1,zt%kk-k+1,i)
	 if(out_zt(thecol,k,i) .ne. out_zt(thecol,k,i)) out_zt(thecol,k,i) = 0.0_r8
      enddo   
    enddo
    
    do i = 1, zm%nn
      do k = 1, zt%kk 
         out_zm(thecol,k,i) = zm%x(1,1,zt%kk-k+1,i)
	 if(out_zm(thecol,k,i) .ne. out_zm(thecol,k,i)) out_zm(thecol,k,i) = 0.0_r8
      enddo   
    enddo

    if (l_output_rad_files) then 
      do i = 1, rad_zt%nn
        do k = 1, rad_zt%kk 
          out_radzt(thecol,k,i) = rad_zt%x(1,1,zt%kk-k+1,i)
	  if(out_radzt(thecol,k,i) .ne. out_radzt(thecol,k,i)) out_radzt(thecol,k,i) = 0.0_r8
        enddo   
      enddo
    
      do i = 1, rad_zm%nn
        do k = 1, rad_zm%kk 
          out_radzm(thecol,k,i) = rad_zm%x(1,1,zt%kk-k+1,i)
	  if(out_radzm(thecol,k,i) .ne. out_radzm(thecol,k,i)) out_radzm(thecol,k,i) = 0.0_r8
        enddo   
      enddo
    endif
    
    do i = 1, sfc%nn
      out_sfc(thecol,1,i) = sfc%x(1,1,1,i)   
      if(out_sfc(thecol,1,i) .ne. out_sfc(thecol,1,i)) out_sfc(thecol,1,i) = 0.0_r8
    enddo

    !  Reset sample fields
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )
    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )
      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )
    end if
    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )


    return

#endif

  end subroutine stats_end_timestep_clubb
  
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS
  
    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, nn

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical, dimension(1,1,kk,nn), intent(out)                 :: l_in_update

    !  Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0_r8
      n(:,:,:,:) = 0
      l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero
  
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  
#ifdef CLUBB_SGS
    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: nn, kk
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m=1,nn
      do k=1,kk

        if ( n(1,1,k,m) > 0 ) then
          x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
        end if

      end do
    end do

    return

  end subroutine stats_avg

#endif
  
end module clubb_intr
