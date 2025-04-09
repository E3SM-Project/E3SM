
module zm_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective 
! precip from the ZM scheme
!
! Apr 2006: RBN: Code added to perform a dilute ascent for closure of the CM mass flux
!                based on an entraining plume a la Raymond and Blythe (1992)
!
! Author: Byron Boville, from code in tphysbc
!
! April 2021: X. Song added code for convective microphysics
! April 2022: X. Song added code for mass flux adjustment
!---------------------------------------------------------------------------------
  use shr_kind_mod,           only: r8 => shr_kind_r8
  use spmd_utils,             only: masterproc
  use ppgrid,                 only: pcols, pver, pverp
  use cloud_fraction,         only: cldfrc_fice
  use physconst,              only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                                    cpwv, cpliq, rh2o
  use cam_abortutils,         only: endrun
  use cam_logfile,            only: iulog
  use zm_conv_cape,           only: compute_dilute_cape
  use zm_conv_types,          only: zm_const_t, zm_param_t
  use zm_conv_util,           only: qsat_hpa ! remove after moving cldprp to new module
  use zm_aero,                only: zm_aero_t
  use zm_microphysics,        only: zm_mphy
  use zm_microphysics_state,  only: zm_microp_st, zm_microp_st_alloc, zm_microp_st_dealloc, zm_microp_st_ini, zm_microp_st_gb

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zmconv_readnl            ! read zmconv_nl namelist
  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! evaporation of precip from ZM schemea
  public trigdcape_ull            ! true if to use dcape-ULL trigger
  public trig_dcape_only          ! true if to use dcape only trigger
  public trig_ull_only            ! true if to ULL along with default CAPE-based trigger
  public zm_microp                ! true for convective microphysics
  public MCSP                     ! true if running MCSP
  public MCSP_heat_coeff          ! MCSP coefficient setting degree of dry static energy transport
  public MCSP_moisture_coeff      ! MCSP coefficient setting degree of moisture transport
  public MCSP_uwind_coeff         ! MCSP coefficient setting degree of zonal wind transport
  public MCSP_vwind_coeff         ! MCSP coefficient setting degree of meridional wind transport

!
! PUBLIC: data
!

   type(zm_const_t), public :: zm_const ! derived type to hold ZM constants
   type(zm_param_t), public :: zm_param ! derived type to hold ZM tunable parameters

!
! Private data
!
   real(r8), parameter :: unset_r8   = huge(1.0_r8)
   integer , parameter :: unset_int  = huge(1)
   real(r8) :: zmconv_c0_lnd         = unset_r8    
   real(r8) :: zmconv_c0_ocn         = unset_r8    
   real(r8) :: zmconv_ke             = unset_r8 
   real(r8) :: zmconv_tau            = unset_r8   
   real(r8) :: zmconv_dmpdz          = unset_r8   
   real(r8) :: zmconv_alfa           = unset_r8   
   real(r8) :: zmconv_tiedke_add     = unset_r8   
   logical  :: zmconv_trigdcape_ull  = .false.
   logical  :: zmconv_trig_dcape_only  = .false.
   logical  :: zmconv_trig_ull_only  = .false.

   logical  :: zmconv_microp         = .false.   ! switch for convective microphysics 
   logical  :: zmconv_clos_dyn_adj   = .false.   ! switch for mass flux adjustment    
   logical  :: zmconv_tpert_fix      = .false.   ! switch for tpert impacts           
   integer  :: zmconv_cape_cin       = unset_int
   integer  :: zmconv_mx_bot_lyr_adj = unset_int
   real(r8) :: zmconv_tp_fac         = unset_r8
   real(r8) :: zmconv_auto_fac       = unset_r8
   real(r8) :: zmconv_accr_fac       = unset_r8
   real(r8) :: zmconv_micro_dcs      = unset_r8
   real(r8) :: zmconv_MCSP_heat_coeff = 0._r8
   real(r8) :: zmconv_MCSP_moisture_coeff = 0._r8
   real(r8) :: zmconv_MCSP_uwind_coeff = 0._r8
   real(r8) :: zmconv_MCSP_vwind_coeff = 0._r8

   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.

!DCAPE-ULL, including options for DCAPE_only and ull_only
   real(r8), parameter :: trigdcapelmt = 0._r8  ! threshold value of dcape for deep convection
   logical :: trigdcape_ull    = .false. !true to use DCAPE trigger and ULL
   logical :: trig_dcape_only  = .false. !true to use DCAPE trigger, ULL not used
   logical :: trig_ull_only    = .false. !true to use ULL along with default CAPE-based trigger

   real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   real(r8) :: alfa_scalar  ! maximum downdraft mass flux fraction  
   logical  :: zm_microp    = .false.   ! switch for convective microphysics                   
   logical  :: clos_dyn_adj = .false.   ! true if apply mass flux adjustment to CAPE closure   

   real(r8) tau   ! convective time scale
   real(r8),parameter :: c1 = 6.112_r8
   real(r8),parameter :: c2 = 17.67_r8
   real(r8),parameter :: c3 = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1
      

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 
   

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair
   
   ! integer,protected :: limcnv ! upper pressure interface level to limit deep convection

   real(r8) :: auto_fac = unset_r8
   real(r8) :: accr_fac = unset_r8
   real(r8) :: micro_dcs= unset_r8

   logical  :: MCSP
   real(r8) :: MCSP_heat_coeff = unset_r8
   real(r8) :: MCSP_moisture_coeff = unset_r8
   real(r8) :: MCSP_uwind_coeff = unset_r8
   real(r8) :: MCSP_vwind_coeff = unset_r8

!===================================================================================================
contains
!===================================================================================================

subroutine zmconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use zm_conv_types,   only: zm_param_mpi_broadcast
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zmconv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_tau, & 
           zmconv_dmpdz, zmconv_alfa, zmconv_tiedke_add,     &
           zmconv_cape_cin, zmconv_mx_bot_lyr_adj, zmconv_tp_fac, zmconv_trigdcape_ull, &
           zmconv_trig_dcape_only, zmconv_trig_ull_only, zmconv_microp, zmconv_auto_fac,&
           zmconv_accr_fac, zmconv_micro_dcs, zmconv_clos_dyn_adj, zmconv_tpert_fix,    &
           zmconv_MCSP_heat_coeff, zmconv_MCSP_moisture_coeff, &
           zmconv_MCSP_uwind_coeff, zmconv_MCSP_vwind_coeff   
   !-----------------------------------------------------------------------------

   zmconv_tau = 3600._r8
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      c0_lnd           = zmconv_c0_lnd
      c0_ocn           = zmconv_c0_ocn
      ke               = zmconv_ke
      tau              = zmconv_tau
      trigdcape_ull    = zmconv_trigdcape_ull
      trig_dcape_only  = zmconv_trig_dcape_only
      trig_ull_only    = zmconv_trig_ull_only
      zm_microp        = zmconv_microp
      clos_dyn_adj     = zmconv_clos_dyn_adj
      auto_fac         = zmconv_auto_fac
      accr_fac         = zmconv_accr_fac
      micro_dcs        = zmconv_micro_dcs
      MCSP_heat_coeff  = zmconv_MCSP_heat_coeff
      MCSP_moisture_coeff = zmconv_MCSP_moisture_coeff
      MCSP_uwind_coeff = zmconv_MCSP_uwind_coeff
      MCSP_vwind_coeff = zmconv_MCSP_vwind_coeff

      ! set zm_param values
      zm_param%trig_dcape      = trigdcape_ull .or. trig_dcape_only
      zm_param%trig_ull        = trigdcape_ull .or. trig_ull_only
      zm_param%tiedke_add      = zmconv_tiedke_add
      zm_param%dmpdz           = zmconv_dmpdz
      zm_param%num_cin         = zmconv_cape_cin
      zm_param%tpert_fix       = zmconv_tpert_fix
      zm_param%tpert_fac       = zmconv_tp_fac
      zm_param%mx_bot_lyr_adj  = zmconv_mx_bot_lyr_adj

      if( abs(MCSP_heat_coeff)+abs(MCSP_moisture_coeff)+abs(MCSP_uwind_coeff)+abs(MCSP_vwind_coeff) > 0._r8 ) then
           MCSP = .true.
      else
           MCSP = .false.
      end if 

      if ( zmconv_alfa /= unset_r8 ) then
           alfa_scalar = zmconv_alfa
      else
           alfa_scalar = 0.1_r8
      end if

      if (trigdcape_ull) then
         write(iulog,*)'**** ZMCONV-DCAPE trigger with unrestricted launch level:', trigdcape_ull
      endif

      if (trig_dcape_only) then
         write(iulog,*)'**** ZM scheme uses DCAPE-only trigger:', trig_dcape_only
      endif

      if (trig_ull_only) then
         write(iulog,*)'**** ZM scheme uses unrestricted launch level along with default CAPE-based trigger:', trig_ull_only
      endif

      if(MCSP) then
         write(iulog,*)'**** ZM scheme uses multiscale coherent structure parameterization (MCSP):',MCSP
      end if

      if(zm_microp) then
         write(iulog,*)'**** ZM scheme uses convective microphysics scheme:',zm_microp
      end if

      if(clos_dyn_adj) then
         write(iulog,*)'**** ZM scheme uses cloud-base mass flux adjustment:',clos_dyn_adj
      end if

      if(zm_param%tpert_fix) then
         write(iulog,*)'**** ZM scheme uses tpert_fix:',zm_param%tpert_fix
      end if
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(ke,                1, mpir8,  0, mpicom)
   call mpibcast(tau,               1, mpir8,  0, mpicom)
   call mpibcast(alfa_scalar,       1, mpir8,  0, mpicom)
   call mpibcast(trigdcape_ull,     1, mpilog, 0, mpicom)
   call mpibcast(trig_dcape_only,   1, mpilog, 0, mpicom)
   call mpibcast(trig_ull_only,     1, mpilog, 0, mpicom)
   call mpibcast(zm_microp,         1, mpilog, 0, mpicom)
   call mpibcast(clos_dyn_adj,      1, mpilog, 0, mpicom)
   call mpibcast(auto_fac,          1, mpir8,  0, mpicom)
   call mpibcast(accr_fac,          1, mpir8,  0, mpicom)
   call mpibcast(micro_dcs,         1, mpir8,  0, mpicom)  
   call mpibcast(MCSP,              1, mpilog, 0, mpicom)
   call mpibcast(MCSP_heat_coeff,   1, mpir8,  0, mpicom)
   call mpibcast(MCSP_moisture_coeff,1, mpir8,  0, mpicom)
   call mpibcast(MCSP_uwind_coeff,  1, mpir8,  0, mpicom)
   call mpibcast(MCSP_vwind_coeff,  1, mpir8,  0, mpicom)

   call zm_param_mpi_broadcast(zm_param)

#endif

end subroutine zmconv_readnl


subroutine zm_convi(limcnv_in, no_deep_pbl_in)

   use dycore,       only: dycore_is, get_resolution
   use zm_conv_types, only: zm_const_set_to_global

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL 

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of ZM constants
   zm_param%limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   if ( present(no_deep_pbl_in) )  then
      no_deep_pbl = no_deep_pbl_in
   else
      no_deep_pbl = .false.
   endif

   ! tau=4800. were used in canadian climate center. however, in echam3 t42, 
   ! convection is too weak, thus adjusted to 2400.

   hgrid = get_resolution()

   ! set zm_const using global values
   call zm_const_set_to_global(zm_const)

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convi: tau',tau
      write(iulog,*) 'tuning parameters zm_convi: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn 
      write(iulog,*) 'tuning parameters zm_convi: ke',ke
      write(iulog,*) 'tuning parameters zm_convi: dmpdz',zm_param%dmpdz
      write(iulog,*) 'tuning parameters zm_convi: alfa',alfa_scalar
      write(iulog,*) 'tuning parameters zm_convi: no_deep_pbl',no_deep_pbl
   endif

   if (masterproc) write(iulog,*)'**** ZM: DILUTE Buoyancy Calculation ****'

end subroutine zm_convi



subroutine zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     ,omega   , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac, &
                    t_star  ,q_star, dcape,   &
                    aero    ,qi      ,dif     ,dnlf    ,dnif    , & 
                    dsf     ,dnsf    ,sprd    ,rice    ,frz     , &
                    mudpcu  ,lambdadpcu, microp_st, wuc)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
! 
!-----------------------------------------------------------------------
   use time_manager, only: is_first_step 
!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qh(pcols,pver)   ! grid slice of specific humidity.
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: omega(pcols,pver)      ! Vertical velocity Pa/s
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
   type(zm_aero_t), intent(inout) :: aero         ! aerosol object. intent(inout) because the
                                                  ! gathered arrays are set here
                                                  ! before passing object
                                                  ! to microphysics
   type(zm_microp_st), intent(inout) :: microp_st ! state and tendency of convective microphysics

!DCAPE-ULL
   real(r8), intent(in), pointer, dimension(:,:) :: t_star ! intermediate T between n and n-1 time step
   real(r8), intent(in), pointer, dimension(:,:) :: q_star ! intermediate q between n and n-1 time step


!
! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
   real(r8), intent(out) :: sprd(pcols,pver)       ! snow production rate
   real(r8), intent(out) :: dif(pcols,pver)        ! detrained convective cloud ice mixing ratio.
   real(r8), intent(out) :: dnlf(pcols,pver)       ! detrained convective cloud water num concen.
   real(r8), intent(out) :: dnif(pcols,pver)       ! detrained convective cloud ice num concen.
   real(r8), intent(out) :: dsf(pcols,pver)        ! detrained convective snow mixing ratio.
   real(r8), intent(out) :: dnsf(pcols,pver)       ! detrained convective snow num concen.
   real(r8), intent(out) :: lambdadpcu(pcols,pver) ! slope of cloud liquid size distr
   real(r8), intent(out) :: mudpcu(pcols,pver)     ! width parameter of droplet size distr
   real(r8), intent(out) :: frz(pcols,pver)        ! freezing heating
   real(r8), intent(out) :: rice(pcols)            ! reserved ice (not yet in cldce) for energy integrals
   real(r8), intent(out) :: qi(pcols,pver)         ! cloud ice mixing ratio. 
   real(r8), intent(inout),optional :: wuc(pcols,pver) ! vertical velocity from ZMmp
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8), intent(out) :: eu(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8), intent(out) :: ed(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   integer,  intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   integer,  intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: rliq(pcols)   ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: dcape(pcols)           ! output dynamical CAPE


   real(r8) zs(pcols)
   real(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) pflxg(pcols,pverp) ! gather precip flux at each level
   real(r8) cug(pcols,pver)    ! gathered condensation rate
   real(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) mumax(pcols)
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.

   integer pblt(pcols)           ! i row of pbl top indices.
   integer pbltg(pcols)          ! i row of pbl top indices.



!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
   real(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
   real(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.
   real(r8) tpm1(pcols,pver)           ! w parcel temperatures at n-1 time step. 
   real(r8) qstpm1(pcols,pver)         ! w parcel temp. saturation mixing ratio at n-1 time step
   real(r8) tlm1(pcols)                ! w LCL parcel Temperature at n-1 time step
   real(r8) capem1(pcols)              ! w CAPE at n-1 time step 
   integer lclm1(pcols)                ! w base level index of deep cumulus convection.
   integer lelm1(pcols)                ! w index of highest theoretical convective plume.
   integer lonm1(pcols)                ! w index of onset level for deep convection.
   integer maxim1(pcols)               ! w index of level with largest moist static energy.

   logical iclosure                    ! switch on sequence of call to buoyan_dilute to derive DCAPE
   real(r8) capelmt_wk                 ! work capelmt to allow diff values passed to closure with trigdcape

   integer lcl(pcols)                  ! w  base level index of deep cumulus convection.
   integer lel(pcols)                  ! w  index of highest theoretical convective plume.

   integer lon(pcols)                  ! w  index of onset level for deep convection.
   integer maxi(pcols)                 ! w  index of level with largest moist static energy.
   integer index(pcols)
   real(r8) precip
!
! gathered work fields:
!
   real(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
   real(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
   real(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
   real(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
   real(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
   real(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
   real(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
   real(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
   real(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
   real(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
   real(r8) cmeg(pcols,pver)
   real(r8) omegag(pcols,pver)         ! wg grid slice of gathered values of omega.
   
   real(r8) rprdg(pcols,pver)          ! wg gathered rain production rate
   real(r8) capeg(pcols)               ! wg gathered convective available potential energy.
   real(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
   real(r8) landfracg(pcols)           ! wg grid slice of landfrac  
   real(r8) tpertg(pcols)              ! wg grid slice of gathered values of tpert (temperature perturbation from PBL)

   integer lclg(pcols)       ! wg gathered values of lcl.
   integer lelg(pcols)
!
! work fields arising from gathered calculations.
!
   real(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
   real(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.
!      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
   real(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
   real(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
   real(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
   real(r8) hmn(pcols,pver)            ! wg moist static energy.
   real(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
   real(r8) qlg(pcols,pver)
   real(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
   real(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.
!      real(r8) ud(pcols,pver)
!      real(r8) vd(pcols,pver)
   real(r8) :: lambdadpcug(pcols,pver) ! slope of cloud liquid size distr
   real(r8) :: mudpcug(pcols,pver)     ! width parameter of droplet size distr

   real(r8) sprdg(pcols,pver)          ! wg gathered snow production rate
   real(r8) dig(pcols,pver)
   real(r8) dsg(pcols,pver)
   real(r8) dnlg(pcols,pver)
   real(r8) dnig(pcols,pver)
   real(r8) dnsg(pcols,pver)

   real(r8) qldeg(pcols,pver)          ! cloud water mixing ratio for detrainment (kg/kg)
   real(r8) qideg(pcols,pver)          ! cloud ice mixing ratio for detrainment (kg/kg)
   real(r8) qsdeg(pcols,pver)          ! snow mixing ratio for detrainment (kg/kg)
   real(r8) ncdeg(pcols,pver)          ! cloud water number concentration for detrainment (1/kg)
   real(r8) nideg(pcols,pver)          ! cloud ice number concentration for detrainment (1/kg)
   real(r8) nsdeg(pcols,pver)          ! snow concentration for detrainment (1/kg)

   real(r8) dsfmg(pcols,pver)          ! mass tendency due to detrainment of snow
   real(r8) dsfng(pcols,pver)          ! num tendency due to detrainment of snow
   real(r8) frzg(pcols,pver)           ! gathered heating rate due to freezing

   type(zm_microp_st)  :: loc_microp_st ! state and tendency of convective microphysics


   real(r8) mb(pcols)                  ! wg cloud base mass flux.

   integer jlcl(pcols)
   integer j0(pcols)                 ! wg detrainment initiation level index.
   integer jd(pcols)                 ! wg downdraft initiation level index.

   real(r8) delt                     ! length of model time-step in seconds.

   integer i
   integer ii
   integer k, kk
   integer msg                      !  ic number of missing moisture levels at the top of model.
   integer ierror

   real(r8) qdifr
   real(r8) sdifr

   integer l, m
   real(r8), parameter :: dcon  = 25.e-6_r8
   real(r8), parameter :: mucon = 5.3_r8
   real(r8) negadq
   logical doliq

   integer dcapemx(pcols)  ! launching level index saved from 1st call for CAPE calculation;  used in 2nd call when DCAPE-ULL active

!
!--------------------------Data statements------------------------------
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = zm_param%limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in cam
!
   qtnd(:,:) = 0._r8
   heat(:,:) = 0._r8
   mcon(:,:) = 0._r8
   rliq(:ncol)   = 0._r8
   rice(:ncol)   = 0._r8    
!
! Allocate microphysics arrays 
   if (zm_microp) call zm_microp_st_alloc(loc_microp_st)
!
! initialize convective tendencies
!
   prec(:ncol) = 0._r8
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k)  = 0._r8
         dsdt(i,k)  = 0._r8
         dudt(i,k)  = 0._r8
         dvdt(i,k)  = 0._r8
         pflx(i,k)  = 0._r8
         pflxg(i,k) = 0._r8
         cme(i,k)   = 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         qlg(i,k)   = 0._r8
         dlf(i,k)   = 0._r8
         dlg(i,k)   = 0._r8
         ! Convective microphysics
         dif(i,k)   = 0._r8
         dsf(i,k)   = 0._r8
         dnlf(i,k)  = 0._r8
         dnif(i,k)  = 0._r8
         dnsf(i,k)  = 0._r8

         dig(i,k)   = 0._r8
         dsg(i,k)   = 0._r8
         dnlg(i,k)  = 0._r8
         dnig(i,k)  = 0._r8
         dnsg(i,k)  = 0._r8

         qi(i,k)    = 0._r8
         sprd(i,k)  = 0._r8
         frz(i,k)   = 0._r8
         
         sprdg(i,k) = 0._r8
         frzg(i,k)   = 0._r8 

         qldeg(i,k)  = 0._r8
         qideg(i,k)  = 0._r8
         qsdeg(i,k)  = 0._r8
         ncdeg(i,k) = 0._r8
         nideg(i,k) = 0._r8
         nsdeg(i,k) = 0._r8
          
         dsfmg(i,k) = 0._r8
         dsfng(i,k) = 0._r8
      end do
   end do

! Initialize microphysics arrays
   if (zm_microp) then
      call zm_microp_st_ini(microp_st,    ncol)
      call zm_microp_st_ini(loc_microp_st,ncol)
   end if

   lambdadpcu  = (mucon + 1._r8)/dcon
   mudpcu      = mucon
   lambdadpcug = lambdadpcu
   mudpcug     = mudpcu

   do i = 1,ncol
      pflx(i,pverp) = 0
      pflxg(i,pverp) = 0
   end do
!
   do i = 1,ncol
      pblt(i) = pver
      dsubcld(i) = 0._r8

      jctop(i) = pver
      jcbot(i) = 1
   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01_r8
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01_r8
         pf(i,k) = paph(i,k)*0.01_r8
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
   do k = pver - 1,msg + 1,-1
      do i = 1,ncol
         if (abs(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5_r8) pblt(i) = k
      end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qh(i,k)
         s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
         tp(i,k)=0.0_r8
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
      end do
   end do

   do i = 1,ncol
      capeg(i) = 0._r8
      lclg(i) = 1
      lelg(i) = pver
      maxg(i) = 1
      tlg(i) = 400._r8
      dsubcld(i) = 0._r8
   end do

   ! Evaluate Tparcel, qs(Tparcel), buoyancy, CAPE, lcl, lel, parcel launch level at index maxi()=hmax
   ! - call #1, iclosure=.true.   standard calculation using state of current step
   ! - call #2, iclosure=.false.  use state from previous step and launch level from call #1
   ! DCAPE is the difference in CAPE between the two calls using the same launch level

   iclosure = .true.
   call compute_dilute_cape( pcols, ncol, pver, pverp, &
                             zm_param%num_cin, msg, &
                             q, t, z, p, pf, &
                             pblt, tpert, &
                             tp, qstp, maxi, tl, &
                             lcl, lel, cape, &
                             zm_const, zm_param, &
                             iclosure )
   if (trigdcape_ull .or. trig_dcape_only) dcapemx(:ncol) = maxi(:ncol)

   ! Calculate dcape trigger condition
   if ( .not.is_first_step() .and. zm_param%trig_dcape ) then
      iclosure = .false.
      call compute_dilute_cape( pcols, ncol, pver, pverp, &
                                zm_param%num_cin, msg, &
                                q_star, t_star, z, p, pf, &
                                pblt, tpert, &
                                tpm1, qstpm1, maxim1, tlm1, &
                                lclm1, lelm1, capem1, &
                                zm_const, zm_param, &
                                iclosure, dcapemx )
      dcape(:ncol) = (cape(:ncol)-capem1(:ncol))/(delt*2._r8)
   endif

!
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).
!
   capelmt_wk = capelmt   ! capelmt_wk default to capelmt for default trigger

   if ((trigdcape_ull .or. trig_dcape_only) .and. (.not. is_first_step()) )  &
      capelmt_wk = 0.0_r8

   lengath = 0
   do i=1,ncol
     if (trigdcape_ull .or. trig_dcape_only) then
     ! DCAPE-ULL
      if (is_first_step()) then
         !Will this cause restart to be non-BFB
           if (cape(i) > capelmt) then
              lengath = lengath + 1
              index(lengath) = i
           end if
       else if (cape(i) > 0.0_r8 .and. dcape(i) > trigdcapelmt) then
           ! use constant 0 or a separate threshold for capt because capelmt is for default trigger
           lengath = lengath + 1
           index(lengath) = i
       endif
     else
      if (cape(i) > capelmt) then
         lengath = lengath + 1
         index(lengath) = i
      end if
     end if
   end do

   if (lengath.eq.0) return
   do ii=1,lengath
      i=index(ii)
      ideep(ii)=i
   end do
!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,pver
      do i = 1,lengath
         dp(i,k) = 0.01_r8*dpp(ideep(i),k)
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         pg(i,k) = p(ideep(i),k)
         zg(i,k) = z(ideep(i),k)
         sg(i,k) = s(ideep(i),k)
         tpg(i,k) = tp(ideep(i),k)
         zfg(i,k) = zf(ideep(i),k)
         qstpg(i,k) = qstp(ideep(i),k)
         omegag(i,k) = omega(ideep(i),k)
         ug(i,k) = 0._r8
         vg(i,k) = 0._r8
      end do
   end do

   if (zm_microp) then

      if (aero%scheme == 'modal') then

         do m = 1, aero%nmodes

            do k = 1,pver
               do i = 1,lengath
                  aero%numg_a(i,k,m) = aero%num_a(m)%val(ideep(i),k)
                  aero%dgnumg(i,k,m) = aero%dgnum(m)%val(ideep(i),k)
               end do
            end do

            do l = 1, aero%nspec(m)
               do k = 1,pver
                  do i = 1,lengath
                     aero%mmrg_a(i,k,l,m) = aero%mmr_a(l,m)%val(ideep(i),k)
                  end do
               end do
            end do

         end do

      else if (aero%scheme == 'bulk') then

         do m = 1, aero%nbulk
            do k = 1,pver
               do i = 1,lengath
                  aero%mmrg_bulk(i,k,m) = aero%mmr_bulk(m)%val(ideep(i),k)
               end do
            end do
         end do

      end if

   end if

!
   do i = 1,lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
   end do
   do i = 1,lengath
      capeg(i) = cape(ideep(i))
      lclg(i) = lcl(ideep(i))
      lelg(i) = lel(ideep(i))
      maxg(i) = maxi(ideep(i))
      tlg(i) = tl(ideep(i))
      landfracg(i) = landfrac(ideep(i))
      pbltg(i) = pblt(ideep(i))
      tpertg(i) = tpert(ideep(i))
   end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + dp(i,k)
         end if
      end do
   end do
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
!            alpha(i,k) = 0.5
         sdifr = 0._r8
         qdifr = 0._r8
         if (sg(i,k) > 0._r8 .or. sg(i,k-1) > 0._r8) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._r8 .or. qg(i,k-1) > 0._r8) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > 1.E-6_r8) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > 1.E-6_r8) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do
!
! obtain cloud properties.
!

   call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg     , &
               cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
               maxg    ,j0      ,jd      ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     , &
               pflxg   ,evpg    ,cug     ,rprdg   ,zm_param%limcnv  , &
               landfracg, tpertg, &  
               aero    ,qhat ,lambdadpcug,mudpcug ,sprdg   ,frzg ,  &
               qldeg   ,qideg   ,qsdeg   ,ncdeg   ,nideg   ,nsdeg,  &
               dsfmg   ,dsfng   ,loc_microp_st )   


!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du   (i,k) = du   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu   (i,k) = eu   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed   (i,k) = ed   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug  (i,k) = cug  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg (i,k) = cmeg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg(i,k) = rprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         sprdg(i,k) = sprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         frzg (i,k) = frzg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         evpg (i,k) = evpg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   call closure(lchnk   , &
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt_wk )
!
! limit cloud base mass flux to theoretical upper bound.
!
   do i=1,lengath
      mumax(i) = 0
   end do
   do k=msg + 2,pver
      do i=1,lengath
        mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
   end do

   do i=1,lengath
      if (mumax(i) > 0._r8) then
         mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
      else
         mb(i) = 0._r8
      endif
   end do

   if (clos_dyn_adj) then 
      do i = 1,lengath
         mb(i) = max(mb(i) - omegag(i,pbltg(i))*0.01_r8, 0._r8)
      end do
   end if

   ! If no_deep_pbl = .true., don't allow convection entirely 
   ! within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (no_deep_pbl) then
      do i=1,lengath
         if (zm(ideep(i),jt(i)) < pblh(ideep(i))) mb(i) = 0
      end do
   end if


   do k=msg+1,pver
      do i=1,lengath
         mu   (i,k)  = mu   (i,k)*mb(i)
         md   (i,k)  = md   (i,k)*mb(i)
         mc   (i,k)  = mc   (i,k)*mb(i)
         du   (i,k)  = du   (i,k)*mb(i)
         eu   (i,k)  = eu   (i,k)*mb(i)
         ed   (i,k)  = ed   (i,k)*mb(i)
         cmeg (i,k)  = cmeg (i,k)*mb(i)
         rprdg(i,k)  = rprdg(i,k)*mb(i)
         cug  (i,k)  = cug  (i,k)*mb(i)
         evpg (i,k)  = evpg (i,k)*mb(i)
         pflxg(i,k+1)= pflxg(i,k+1)*mb(i)*100._r8/grav
         sprdg(i,k)  = sprdg(i,k)*mb(i)
         frzg(i,k)   = frzg(i,k)*mb(i)


         if ( zm_microp .and. mb(i).eq.0._r8) then
            qlg (i,k) = 0._r8
            dsfmg(i,k) = 0._r8
            dsfng(i,k) = 0._r8 
            frzg (i,k) = 0._r8
            loc_microp_st%wu(i,k) = 0._r8
            loc_microp_st%qliq (i,k) = 0._r8
            loc_microp_st%qice (i,k) = 0._r8
            loc_microp_st%qrain(i,k) = 0._r8
            loc_microp_st%qsnow(i,k) = 0._r8
            loc_microp_st%qgraupel(i,k)= 0._r8
            loc_microp_st%qnl (i,k) = 0._r8
            loc_microp_st%qni (i,k) = 0._r8
            loc_microp_st%qnr (i,k) = 0._r8
            loc_microp_st%qns (i,k) = 0._r8
            loc_microp_st%qng (i,k)  = 0._r8

            loc_microp_st%autolm(i,k) = 0._r8
            loc_microp_st%accrlm(i,k) = 0._r8
            loc_microp_st%bergnm(i,k) = 0._r8
            loc_microp_st%fhtimm(i,k) = 0._r8
            loc_microp_st%fhtctm(i,k) = 0._r8
            loc_microp_st%fhmlm (i,k) = 0._r8
            loc_microp_st%hmpim (i,k) = 0._r8
            loc_microp_st%accslm(i,k) = 0._r8
            loc_microp_st%dlfm  (i,k) = 0._r8

            loc_microp_st%autoln(i,k) = 0._r8
            loc_microp_st%accrln(i,k) = 0._r8
            loc_microp_st%bergnn(i,k) = 0._r8
            loc_microp_st%fhtimn(i,k) = 0._r8
            loc_microp_st%fhtctn(i,k) = 0._r8
            loc_microp_st%fhmln (i,k) = 0._r8
            loc_microp_st%accsln(i,k) = 0._r8
            loc_microp_st%activn(i,k) = 0._r8
            loc_microp_st%dlfn  (i,k) = 0._r8
            loc_microp_st%cmel  (i,k) = 0._r8


            loc_microp_st%autoim(i,k) = 0._r8
            loc_microp_st%accsim(i,k) = 0._r8
            loc_microp_st%difm  (i,k) = 0._r8
            loc_microp_st%cmei  (i,k) = 0._r8

            loc_microp_st%nuclin(i,k) = 0._r8
            loc_microp_st%autoin(i,k) = 0._r8
            loc_microp_st%accsin(i,k) = 0._r8
            loc_microp_st%hmpin (i,k) = 0._r8
            loc_microp_st%difn  (i,k) = 0._r8

            loc_microp_st%trspcm(i,k) = 0._r8
            loc_microp_st%trspcn(i,k) = 0._r8
            loc_microp_st%trspim(i,k) = 0._r8
            loc_microp_st%trspin(i,k) = 0._r8

            loc_microp_st%accgrm(i,k) = 0._r8
            loc_microp_st%accglm(i,k) = 0._r8
            loc_microp_st%accgslm(i,k)= 0._r8
            loc_microp_st%accgsrm(i,k)= 0._r8
            loc_microp_st%accgirm(i,k)= 0._r8
            loc_microp_st%accgrim(i,k)= 0._r8
            loc_microp_st%accgrsm(i,k)= 0._r8

            loc_microp_st%accgsln(i,k)= 0._r8
            loc_microp_st%accgsrn(i,k)= 0._r8
            loc_microp_st%accgirn(i,k)= 0._r8

            loc_microp_st%accsrim(i,k)= 0._r8
            loc_microp_st%acciglm(i,k)= 0._r8
            loc_microp_st%accigrm(i,k)= 0._r8
            loc_microp_st%accsirm(i,k)= 0._r8

            loc_microp_st%accigln(i,k)= 0._r8
            loc_microp_st%accigrn(i,k)= 0._r8
            loc_microp_st%accsirn(i,k)= 0._r8
            loc_microp_st%accgln(i,k) = 0._r8
            loc_microp_st%accgrn(i,k) = 0._r8

            loc_microp_st%accilm(i,k) = 0._r8
            loc_microp_st%acciln(i,k) = 0._r8

            loc_microp_st%fallrm(i,k) = 0._r8
            loc_microp_st%fallsm(i,k) = 0._r8
            loc_microp_st%fallgm(i,k) = 0._r8
            loc_microp_st%fallrn(i,k) = 0._r8
            loc_microp_st%fallsn(i,k) = 0._r8
            loc_microp_st%fallgn(i,k) = 0._r8
            loc_microp_st%fhmrm (i,k) = 0._r8
         end if
      end do
   end do
!
! compute temperature and moisture changes due to convection.
!
   call q1q2_pjr(lchnk   , &
                 dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                 su      ,du      ,qhat    ,shat    ,dp      , &
                 mu      ,md      ,sd      ,qd      ,qldeg   , &
                 dsubcld ,jt      ,maxg    ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg     ,evpg    ,cug     ,qideg   ,dig     , &
                 ncdeg   ,nideg   ,dnlg    ,dnig    ,frzg    , &
                 qsdeg   ,nsdeg   ,dsg     ,dnsg    )

!
! Conservation check
!
  if (zm_microp) then
    do k = msg + 1,pver
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
      do i = 1,lengath
         if (dqdt(i,k)*2._r8*delt+qg(i,k)<0._r8) then
            negadq = dqdt(i,k)+0.5_r8*qg(i,k)/delt
            dqdt(i,k) = dqdt(i,k)-negadq

            ! First evaporate precipitation from k layer to cloud top assuming that the preciptation 
            ! above will fall down and evaporate at k layer. So dsdt will be applied at k layer.
            do kk=k,jt(i),-1
              if (negadq<0._r8) then
                 if (rprdg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                 ! precipitation is enough
                    dsdt(i,k) = dsdt(i,k) + negadq*rl/cpres
                    if (rprdg(i,kk)>sprdg(i,kk)) then
                    ! if there is rain, evaporate it first
                       if(rprdg(i,kk)-sprdg(i,kk)<-negadq*dp(i,k)/dp(i,kk)) then
                       ! if rain is not enough, evaporate snow and graupel
                         dsdt(i,k) = dsdt(i,k) + (negadq+ (rprdg(i,kk)-sprdg(i,kk))*dp(i,kk)/dp(i,k))*latice/cpres
                         sprdg(i,kk) = negadq*dp(i,k)/dp(i,kk)+rprdg(i,kk)
                       end if
                    else
                    ! if there is not rain, evaporate snow and graupel
                       sprdg(i,kk) = sprdg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                       dsdt(i,k) = dsdt(i,k) + negadq*latice/cpres
                    end if
                    rprdg(i,kk) = rprdg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                    negadq = 0._r8
                 else
                 ! precipitation is not enough. calculate the residue and evaporate next layer
                    negadq = rprdg(i,kk)*dp(i,kk)/dp(i,k)+negadq
                    dsdt(i,k) = dsdt(i,k) - rprdg(i,kk)*rl/cpres*dp(i,kk)/dp(i,k)
                    dsdt(i,k) = dsdt(i,k) - sprdg(i,kk)*latice/cpres*dp(i,kk)/dp(i,k)
                    sprdg(i,kk) = 0._r8
                    rprdg(i,kk) = 0._r8
                 end if


                 if (negadq<0._r8) then
                     if (dlg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                    ! first evaporate (detrained) cloud liquid water
                        dsdt(i,k) = dsdt(i,k) + negadq*rl/cpres
                        dnlg(i,kk) = dnlg(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dlg(i,kk))
                        dlg(i,kk)  = dlg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                        negadq = 0._r8
                     else
                    ! if cloud liquid water is not enough then calculate the residual and evaporate the detrained cloud ice
                        negadq = negadq + dlg(i,kk)*dp(i,kk)/dp(i,k)
                        dsdt(i,k) = dsdt(i,k) - dlg(i,kk)*dp(i,kk)/dp(i,k)*rl/cpres
                        dlg(i,kk) = 0._r8
                        dnlg(i,kk) = 0._r8
                       if (dig(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                         dsdt(i,k) = dsdt(i,k) + negadq*(rl+latice)/cpres
                         dnig(i,kk) = dnig(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dig(i,kk))
                         dig(i,kk)  = dig(i,kk)+negadq*dp(i,k)/dp(i,kk)
                         negadq = 0._r8
                       else
                       ! if cloud ice is not enough, then calculate the residual and evaporate the detrained snow
                         negadq = negadq + dig(i,kk)*dp(i,kk)/dp(i,k)
                         dsdt(i,k) = dsdt(i,k) - dig(i,kk)*dp(i,kk)/dp(i,k)*(rl+latice)/cpres
                         dig(i,kk) = 0._r8
                         dnig(i,kk) = 0._r8
                         if (dsg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                           dsdt(i,k) = dsdt(i,k) + negadq*(rl+latice)/cpres
                           dnsg(i,kk) = dnsg(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dsg(i,kk))
                           dsg(i,kk)  = dsg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                           negadq = 0._r8
                         else
                         ! if cloud ice is not enough, then calculate the residual and evaporate next layer
                           negadq = negadq + dsg(i,kk)*dp(i,kk)/dp(i,k)
                           dsdt(i,k) = dsdt(i,k) - dsg(i,kk)*dp(i,kk)/dp(i,k)*(rl+latice)/cpres
                           dsg(i,kk) = 0._r8
                           dnsg(i,kk) = 0._r8
                         end if
                        end if
                     end if
                 end if

              end if
            end do

            if (negadq<0._r8) then
               dqdt(i,k) = dqdt(i,k) - negadq
            end if

          end if
      end do
   end do
  end if



! gather back temperature and mixing ratio.
!
   do k = msg + 1,pver
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
      do i = 1,lengath
!
! q is updated to compute net precip.
!
         q(ideep(i),k) = qh(ideep(i),k) + 2._r8*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         cme (ideep(i),k) = cmeg (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*cpres
         dlf (ideep(i),k) = dlg  (i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql  (ideep(i),k) = qlg  (i,k)

         sprd(ideep(i),k) = sprdg(i,k)
         dif (ideep(i),k) = dig  (i,k)
         dsf (ideep(i),k) = dsg  (i,k)
         dnlf(ideep(i),k) = dnlg  (i,k)
         dnif(ideep(i),k) = dnig  (i,k)
         dnsf(ideep(i),k) = dnsg  (i,k)
         lambdadpcu(ideep(i),k) = lambdadpcug(i,k)
         mudpcu(ideep(i),k)     = mudpcug(i,k)
         frz(ideep(i),k)  = frzg(i,k)*latice/cpres
         if (zm_microp) qi(ideep(i),k) = loc_microp_st%qice(i,k)
      end do
   end do

! Gather back microphysics arrays.
   if (zm_microp)  call zm_microp_st_gb(microp_st,loc_microp_st,ideep,lengath)
   
   if (zm_microp) then
     do k = msg + 1,pver
       do i = 1,ncol

         !Interpolate variable from interface to mid-layer.

         if(k.lt.pver) then
            microp_st%qice (i,k) = 0.5_r8*(microp_st%qice(i,k)+microp_st%qice(i,k+1))
            microp_st%qliq (i,k) = 0.5_r8*(microp_st%qliq(i,k)+microp_st%qliq(i,k+1))
            microp_st%qrain (i,k) = 0.5_r8*(microp_st%qrain(i,k)+microp_st%qrain(i,k+1))
            microp_st%qsnow (i,k) = 0.5_r8*(microp_st%qsnow(i,k)+microp_st%qsnow(i,k+1))
            microp_st%qgraupel(i,k) = 0.5_r8*(microp_st%qgraupel(i,k)+microp_st%qgraupel(i,k+1))
            microp_st%qni (i,k) = 0.5_r8*(microp_st%qni(i,k)+microp_st%qni(i,k+1))
            microp_st%qnl (i,k) = 0.5_r8*(microp_st%qnl(i,k)+microp_st%qnl(i,k+1))
            microp_st%qnr (i,k) = 0.5_r8*(microp_st%qnr(i,k)+microp_st%qnr(i,k+1))
            microp_st%qns (i,k) = 0.5_r8*(microp_st%qns(i,k)+microp_st%qns(i,k+1))
            microp_st%qng (i,k) = 0.5_r8*(microp_st%qng(i,k)+microp_st%qng(i,k+1))
            microp_st%wu(i,k)   = 0.5_r8*(microp_st%wu(i,k)+microp_st%wu(i,k+1))
            wuc(i,k) = microp_st%wu(i,k)
         end if

         if (t(i,k).gt.tmelt .and. t(i,k-1).le.tmelt) then
             microp_st%qice (i,k-1) = microp_st%qice (i,k-1) + microp_st%qice (i,k)
             microp_st%qice (i,k) = 0._r8
             microp_st%qni (i,k-1) = microp_st%qni (i,k-1) + microp_st%qni (i,k)
             microp_st%qni (i,k) = 0._r8
             microp_st%qsnow (i,k-1) = microp_st%qsnow (i,k-1) + microp_st%qsnow (i,k)
             microp_st%qsnow (i,k) = 0._r8
             microp_st%qns (i,k-1) = microp_st%qns (i,k-1) + microp_st%qns (i,k)
             microp_st%qns (i,k) = 0._r8
             microp_st%qgraupel (i,k-1) = microp_st%qgraupel (i,k-1) + microp_st%qgraupel (i,k)
             microp_st%qgraupel (i,k) = 0._r8
             microp_st%qng (i,k-1) = microp_st%qng (i,k-1) + microp_st%qng (i,k)
             microp_st%qng (i,k) = 0._r8
         end if
         !Convert it from units of "kg/kg" to "g/m3"
         microp_st%qice (i,k) = microp_st%qice(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         microp_st%qliq (i,k) = microp_st%qliq(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         microp_st%qrain (i,k) = microp_st%qrain(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         microp_st%qsnow (i,k) = microp_st%qsnow(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         microp_st%qgraupel (i,k) = microp_st%qgraupel(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         microp_st%qni (i,k) = microp_st%qni(i,k) * pap(i,k)/t(i,k)/rgas
         microp_st%qnl (i,k) = microp_st%qnl(i,k) * pap(i,k)/t(i,k)/rgas
         microp_st%qnr (i,k) = microp_st%qnr(i,k) * pap(i,k)/t(i,k)/rgas
         microp_st%qns (i,k) = microp_st%qns(i,k) * pap(i,k)/t(i,k)/rgas
         microp_st%qng (i,k) = microp_st%qng(i,k) * pap(i,k)/t(i,k)/rgas 
       end do
     end do
   end if
!
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
   do i = 1,lengath
      jctop(ideep(i)) = jt(i)
!++bee
      jcbot(ideep(i)) = maxg(i)
!--bee
      pflx(ideep(i),pverp) = pflxg(i,pverp)
   end do

! Compute precip by integrating change in water vapor minus detrained cloud water
   do k = pver,msg + 1,-1
      do i = 1,ncol
          prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k)) - dpp(i,k)*(dlf(i,k)+dif(i,k)+dsf(i,k))*2._r8*delt
      end do
   end do

! obtain final precipitation rate in m/s.
   do i = 1,ncol
      prec(i) = rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
   end do

! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
   do k = 1, pver
      do i = 1, ncol
          rliq(i) = rliq(i) + (dlf(i,k)+dif(i,k)+dsf(i,k))*dpp(i,k)/gravit
          rice(i) = rice(i) + (dif(i,k)+dsf(i,k))*dpp(i,k)/gravit
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8
   rice(:ncol) = rice(:ncol) /1000._r8    

! Deallocate microphysics arrays.
   if (zm_microp) call zm_microp_st_dealloc(loc_microp_st)

   return
end subroutine zm_convr

!===============================================================================
subroutine zm_conv_evap(ncol,lchnk, &
     t,pmid,pdel,q, &
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat,  &
     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, prdsnow, old_snow )
!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

    use wv_saturation,  only: qsat
    use phys_grid, only: get_rlat_all_p

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol, lchnk             ! number of columns and chunk index
    real(r8),intent(in), dimension(pcols,pver) :: t          ! temperature (K)
    real(r8),intent(in), dimension(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    real(r8),intent(in), dimension(pcols,pver) :: pdel       ! layer thickness (Pa)
    real(r8),intent(in), dimension(pcols,pver) :: q          ! water vapor (kg/kg)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwprd ! Heating rate of snow production
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    


    real(r8), intent(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step

    real(r8), intent(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    real(r8), intent(out)   :: snow(pcols)        ! Convective-scale snowfall rate

    ! Convective microphysics
    real(r8), intent(in   ) :: prdsnow(pcols,pver)! snow production (kg/ks/s)
    logical,  intent(in)    :: old_snow           ! true for old estimate of snow production
!
!---------------------------Local storage-------------------------------

    real(r8) :: es    (pcols,pver)    ! Saturation vapor pressure
    real(r8) :: fice   (pcols,pver)    ! ice fraction in precip production
    real(r8) :: fsnow_conv(pcols,pver) ! snow fraction in precip production
    real(r8) :: qs   (pcols,pver)    ! saturation specific humidity
    real(r8),intent(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),intent(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8),intent(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    real(r8),intent(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
    real(r8) :: work1                  ! temp variable (pjr)
    real(r8) :: work2                  ! temp variable (pjr)

    real(r8) :: evpvint(pcols)         ! vertical integral of evaporation
    real(r8) :: evpprec(pcols)         ! evaporation of precipitation (kg/kg/s)
    real(r8) :: evpsnow(pcols)         ! evaporation of snowfall (kg/kg/s)
    real(r8) :: snowmlt(pcols)         ! snow melt tendency in layer
    real(r8) :: flxsntm(pcols)         ! flux of snow into layer, after melting

    real(r8) :: evplimit               ! temp variable for evaporation limits
    real(r8) :: rlat(pcols)

    real(r8) :: dum
    real(r8) :: omsm
    integer :: i,k                     ! longitude,level indices


!-----------------------------------------------------------------------

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._r8

! determine saturation vapor pressure
    call qsat(t(1:ncol, 1:pver), pmid(1:ncol, 1:pver), &
         es(1:ncol, 1:pver), qs(1:ncol, 1:pver))

! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldfrc_fice(ncol, t, fice, fsnow_conv)

! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0._r8
    flxsnow(:ncol,1) = 0._r8
    evpvint(:ncol)   = 0._r8
    omsm=0.99999_r8          ! to prevent problems due to round off error

    do k = 1, pver
       do i = 1, ncol

! Melt snow falling into layer, if necessary. 
        if( old_snow ) then
          if (t(i,k) > tmelt) then
             flxsntm(i) = 0._r8
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
          end if
        else
        ! make sure melting snow doesn't reduce temperature below threshold
          if (t(i,k) > tmelt) then
              dum = -latice/cpres*flxsnow(i,k)*gravit/pdel(i,k)*deltat
              if (t(i,k) + dum .le. tmelt) then
                dum = (t(i,k)-tmelt)*cpres/latice/deltat
                dum = dum/(flxsnow(i,k)*gravit/pdel(i,k))
                dum = max(0._r8,dum)
                dum = min(1._r8,dum)
              else
                dum = 1._r8
              end if
              dum = dum*omsm
              flxsntm(i) = flxsnow(i,k)*(1.0_r8-dum)
              snowmlt(i) = dum*flxsnow(i,k)*gravit/ pdel(i,k)
           else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
           end if
         end if


! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._r8 - q(i,k)/qs(i,k), 0._r8)

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))
!**********************************************************
!!          evpprec(i) = 0.    ! turn off evaporation for now
!**********************************************************

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qs
          evplimit   = max(0._r8, (qs(i,k)-q(i,k)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
          evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = min(evplimit, evpprec(i))

          if( .not.old_snow ) then
            evpprec(i) = max(0._r8, evpprec(i))
            evpprec(i) = evpprec(i)*omsm
          end if

! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(i,k) > 0._r8) then
!            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
!            prevent roundoff problems
             work1 = min(max(0._r8,flxsntm(i)/flxprec(i,k)),1._r8)
             if (.not.old_snow .and. prdsnow(i,k)>prdprec(i,k)) work1 = 1._r8
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._r8
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
!pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
!pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
! the small amount added to flxprec in the work1 expression has been increased from 
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.

      if( old_snow ) then
#ifdef PERGRO
          work1 = min(max(0._r8,flxsnow(i,k)/(flxprec(i,k)+8.64e-11_r8)),1._r8)
#else
          if (flxprec(i,k).gt.0._r8) then
             work1 = min(max(0._r8,flxsnow(i,k)/flxprec(i,k)),1._r8)
          else
             work1 = 0._r8
          endif
#endif
          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i).gt.0._r8) work2 = 0._r8
!         work2 = fsnow_conv(i,k)
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k)*work2*latice
          tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*latice
       else
          ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i))
          tend_s_snwprd  (i,k) = prdsnow(i,k)*latice
          tend_s_snwevmlt(i,k) = -min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i) )*latice
       end if

! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

! protect against rounding error
          flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
          flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)
! more protection (pjr)
!         flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

! heating (cooling) and moistening due to evaporation 
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
          if( old_snow ) then
             tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          else
             tend_s(i,k)   =-evpprec(i)*latvap + tend_s_snwevmlt(i,k)
          end if

          tend_q(i,k) = evpprec(i)
       end do
    end do

! protect against rounding error
    if( .not.old_snow ) then
      do i = 1, ncol
        if(flxsnow(i,pverp).gt.flxprec(i,pverp)) then
           dum = (flxsnow(i,pverp)-flxprec(i,pverp))*gravit
           do k = pver, 1, -1
             if (ntsnprd(i,k)>ntprprd(i,k).and. dum > 0._r8) then
                ntsnprd(i,k) = ntsnprd(i,k) - dum/pdel(i,k)
                tend_s_snwevmlt(i,k) = tend_s_snwevmlt(i,k) - dum/pdel(i,k)*latice
                tend_s(i,k)  = tend_s(i,k) - dum/pdel(i,k)*latice
                dum = 0._r8
             end if
           end do
           flxsnow(i,pverp) = flxprec(i,pverp)
        end if
      end do
    end if


! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000._r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000._r8

!**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
!**********************************************************

  end subroutine zm_conv_evap

!=========================================================================================

subroutine cldprp(lchnk   , &
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      , &
                  cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
                  mx      ,j0      ,jd      ,rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     , &
                  pflx    ,evp     ,cu      ,rprd    ,limcnv  , &
                  landfrac,tpertg  , &
                  aero    ,qhat ,lambdadpcu ,mudpcu  ,sprd   ,frz1 , &
                  qcde    ,qide   ,qsde     ,ncde    ,nide   ,nsde , &
                  dsfm    ,dsfn   ,loc_microp_st )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
! 
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------

   implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier

   real(r8), intent(in) :: q(pcols,pver)         ! spec. humidity of env
   real(r8), intent(in) :: t(pcols,pver)         ! temp of env
   real(r8), intent(in) :: p(pcols,pver)         ! pressure of env
   real(r8), intent(in) :: z(pcols,pver)         ! height of env
   real(r8), intent(in) :: s(pcols,pver)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
   real(r8), intent(in) :: u(pcols,pver)         ! zonal velocity of env
   real(r8), intent(in) :: v(pcols,pver)         ! merid. velocity of env

   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac

   integer, intent(in) :: jb(pcols)              ! updraft base level
   integer, intent(in) :: lel(pcols)             ! updraft launch level
   integer, intent(out) :: jt(pcols)              ! updraft plume top
   integer, intent(out) :: jlcl(pcols)            ! updraft lifting cond level
   integer, intent(in) :: mx(pcols)              ! updraft base level (same is jb)
   integer, intent(out) :: j0(pcols)              ! level where updraft begins detraining
   integer, intent(out) :: jd(pcols)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pcols,pver)      ! interface values of dry stat energy
   real(r8), intent(in) :: tpertg(pcols)

   real(r8), intent(in) :: qhat(pcols,pver)      ! wg grid slice of upper interface mixing ratio.
   type(zm_aero_t), intent(in) :: aero           ! aerosol object

!
! output
!
   real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
   real(r8), intent(out) :: qst(pcols,pver)      ! saturation mixing ratio of env.
   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft

   ! Convective microphysics
   type(zm_microp_st)  :: loc_microp_st ! state and tendency of convective microphysics

   real(r8), intent(out) :: qcde(pcols,pver)     ! cloud water mixing ratio for detrainment (kg/kg)
   real(r8), intent(out) :: qide(pcols,pver)     ! cloud ice mixing ratio for detrainment (kg/kg)
   real(r8), intent(out) :: qsde(pcols,pver)     ! snow mixing ratio for detrainment (kg/kg)
   real(r8), intent(out) :: ncde(pcols,pver)     ! cloud water number concentration for detrainment (1/kg)
   real(r8), intent(out) :: nide(pcols,pver)     ! cloud ice number concentration for detrainment (1/kg)
   real(r8), intent(out) :: nsde(pcols,pver)     ! snow number concentration for detrainment (1/kg)
   real(r8), intent(out) :: sprd(pcols,pver)     ! rate of production of snow at that layer

   ! tendency for output

   real(r8), intent(out) :: dsfm  (pcols,pver)   !mass tendency due to detrainment of snow
   real(r8), intent(out) :: dsfn  (pcols,pver)   !num tendency due to detrainment of snow

   real(r8), intent(inout) :: lambdadpcu(pcols,pver) ! slope of cloud liquid size distr
   real(r8), intent(inout) :: mudpcu(pcols,pver)     ! width parameter of droplet size distr

   real(r8) rd                   ! gas constant for dry air
   real(r8) grav                 ! gravity
   real(r8) cp                   ! heat capacity of dry air

!
! Local workspace
!
   real(r8) gamma(pcols,pver)
   real(r8) dz(pcols,pver)
   real(r8) iprm(pcols,pver)
   real(r8) hu(pcols,pver)
   real(r8) hd(pcols,pver)
   real(r8) eps(pcols,pver)
   real(r8) f(pcols,pver)
   real(r8) k1(pcols,pver)
   real(r8) i2(pcols,pver)
   real(r8) ihat(pcols,pver)
   real(r8) i3(pcols,pver)
   real(r8) idag(pcols,pver)
   real(r8) i4(pcols,pver)
   real(r8) qsthat(pcols,pver)
   real(r8) hsthat(pcols,pver)
   real(r8) gamhat(pcols,pver)
   real(r8) cu(pcols,pver)
   real(r8) evp(pcols,pver)
   real(r8) cmeg(pcols,pver)
   real(r8) qds(pcols,pver)
! RBN For c0mask
   real(r8) c0mask(pcols)

   real(r8) hmin(pcols)
   real(r8) expdif(pcols)
   real(r8) expnum(pcols)
   real(r8) ftemp(pcols)
   real(r8) eps0(pcols)
   real(r8) rmue(pcols)
   real(r8) zuef(pcols)
   real(r8) zdef(pcols)
   real(r8) epsm(pcols)
   real(r8) ratmjb(pcols)
   real(r8) est(pcols)
   real(r8) totpcp(pcols)
   real(r8) totevp(pcols)
   real(r8) alfa(pcols)
   real(r8) ql1
   real(r8) tu
   real(r8) estu
   real(r8) qstu

   real(r8) small
   real(r8) mdt

   ! Convective microphysics
   real(r8) fice(pcols,pver)        ! ice fraction in precip production
   real(r8) tug(pcols,pver)

   real(r8) totfrz(pcols)
   real(r8) frz1(pcols,pver)        ! rate of freezing
   real(r8) frz (pcols,pver)        ! rate of freezing
   real(r8) pflxs(pcols,pverp)      ! frozen precipitation flux thru layer
   real(r8) dum, sdum 

   real(r8), parameter :: omsm=0.99999_r8       ! to prevent problems due to round off error
   real(r8), parameter :: mu_min = 0.02_r8      ! minimum value of mu
   real(r8), parameter :: t_homofrz = 233.15_r8 ! homogeneous freezing temperature
   real(r8), parameter :: t_mphase  = 40._r8    ! mixed phase temperature = tmelt-t_homofrz = 273.15K - 233.15K

   integer  jto(pcols)              ! updraft plume old top
   integer  tmplel(pcols)

   integer  iter, itnum
   integer  m

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical doit(pcols)
   logical done(pcols)

!
!------------------------------------------------------------------------------
   dsfm  (:il2g,:) = 0._r8
   dsfn  (:il2g,:) = 0._r8
   if (zm_microp) then
      loc_microp_st%autolm(:il2g,:) = 0._r8
      loc_microp_st%accrlm(:il2g,:) = 0._r8
      loc_microp_st%bergnm(:il2g,:) = 0._r8
      loc_microp_st%fhtimm(:il2g,:) = 0._r8
      loc_microp_st%fhtctm(:il2g,:) = 0._r8
      loc_microp_st%fhmlm (:il2g,:) = 0._r8
      loc_microp_st%hmpim (:il2g,:) = 0._r8
      loc_microp_st%accslm(:il2g,:) = 0._r8
      loc_microp_st%dlfm  (:il2g,:) = 0._r8

      loc_microp_st%autoln(:il2g,:) = 0._r8
      loc_microp_st%accrln(:il2g,:) = 0._r8
      loc_microp_st%bergnn(:il2g,:) = 0._r8
      loc_microp_st%fhtimn(:il2g,:) = 0._r8
      loc_microp_st%fhtctn(:il2g,:) = 0._r8
      loc_microp_st%fhmln (:il2g,:) = 0._r8
      loc_microp_st%accsln(:il2g,:) = 0._r8
      loc_microp_st%activn(:il2g,:) = 0._r8
      loc_microp_st%dlfn  (:il2g,:) = 0._r8

      loc_microp_st%autoim(:il2g,:) = 0._r8
      loc_microp_st%accsim(:il2g,:) = 0._r8
      loc_microp_st%difm  (:il2g,:) = 0._r8


      loc_microp_st%nuclin(:il2g,:) = 0._r8
      loc_microp_st%autoin(:il2g,:) = 0._r8
      loc_microp_st%accsin(:il2g,:) = 0._r8
      loc_microp_st%hmpin (:il2g,:) = 0._r8
      loc_microp_st%difn  (:il2g,:) = 0._r8

      loc_microp_st%trspcm(:il2g,:) = 0._r8
      loc_microp_st%trspcn(:il2g,:) = 0._r8
      loc_microp_st%trspim(:il2g,:) = 0._r8
      loc_microp_st%trspin(:il2g,:) = 0._r8

      do k = 1,pver
        do i = 1,il2g
         loc_microp_st%accgrm(i,k) = 0._r8
         loc_microp_st%accglm(i,k) = 0._r8
         loc_microp_st%accgslm(i,k)= 0._r8
         loc_microp_st%accgsrm(i,k)= 0._r8
         loc_microp_st%accgirm(i,k)= 0._r8
         loc_microp_st%accgrim(i,k)= 0._r8
         loc_microp_st%accgrsm(i,k)= 0._r8

         loc_microp_st%accgsln(i,k)= 0._r8
         loc_microp_st%accgsrn(i,k)= 0._r8
         loc_microp_st%accgirn(i,k)= 0._r8

         loc_microp_st%accsrim(i,k)= 0._r8
         loc_microp_st%acciglm(i,k)= 0._r8
         loc_microp_st%accigrm(i,k)= 0._r8
         loc_microp_st%accsirm(i,k)= 0._r8

         loc_microp_st%accigln(i,k)= 0._r8
         loc_microp_st%accigrn(i,k)= 0._r8
         loc_microp_st%accsirn(i,k)= 0._r8
         loc_microp_st%accgln(i,k) = 0._r8
         loc_microp_st%accgrn(i,k) = 0._r8

         loc_microp_st%accilm(i,k) = 0._r8
         loc_microp_st%acciln(i,k) = 0._r8

         loc_microp_st%fallrm(i,k) = 0._r8
         loc_microp_st%fallsm(i,k) = 0._r8
         loc_microp_st%fallgm(i,k) = 0._r8
         loc_microp_st%fallrn(i,k) = 0._r8
         loc_microp_st%fallsn(i,k) = 0._r8
         loc_microp_st%fallgn(i,k) = 0._r8
         loc_microp_st%fhmrm (i,k) = 0._r8      
       end do
      end do
   end if
!
   do i = 1,il2g
      ftemp(i) = 0._r8
      expnum(i) = 0._r8
      expdif(i) = 0._r8
      c0mask(i)  = c0_ocn * (1._r8-landfrac(i)) +   c0_lnd * landfrac(i) 
   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
   do k = 1,pver
      do i = 1,il2g
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!
! initialize many output and work variables to zero
!
   pflx(:il2g,1) = 0

   do k = 1,pver
      do i = 1,il2g
         k1(i,k) = 0._r8
         i2(i,k) = 0._r8
         i3(i,k) = 0._r8
         i4(i,k) = 0._r8
         mu(i,k) = 0._r8
         f(i,k) = 0._r8
         eps(i,k) = 0._r8
         eu(i,k) = 0._r8
         du(i,k) = 0._r8
         ql(i,k) = 0._r8
         cu(i,k) = 0._r8
         evp(i,k) = 0._r8
         cmeg(i,k) = 0._r8
         qds(i,k) = q(i,k)
         md(i,k) = 0._r8
         ed(i,k) = 0._r8
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0._r8
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
         call qsat_hPa(t(i,k), p(i,k), est(i), qst(i,k))
!++bee
         if ( p(i,k)-est(i) <= 0._r8 ) then
            qst(i,k) = 1.0_r8
         end if
!--bee
         gamma(i,k) = qst(i,k)*(1._r8 + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         rprd(i,k) = 0._r8
         ! Convective microphysics
         sprd(i,k) = 0._r8
         fice(i,k) = 0._r8
         tug(i,k)  = 0._r8
         qcde(i,k) = 0._r8
         qide(i,k) = 0._r8
         qsde(i,k) = 0._r8
         ncde(i,k) = 0._r8
         nide(i,k) = 0._r8
         nsde(i,k) = 0._r8
         frz(i,k)  = 0._r8
         frz1(i,k) = 0._r8
         if (zm_microp) then
           loc_microp_st%cmel(i,k) = 0._r8
           loc_microp_st%cmei(i,k) = 0._r8
           loc_microp_st%wu(i,k)   = 0._r8
           loc_microp_st%qliq(i,k) = 0._r8
           loc_microp_st%qice(i,k) = 0._r8
           loc_microp_st%qrain(i,k)= 0._r8
           loc_microp_st%qsnow(i,k)= 0._r8
           loc_microp_st%qgraupel(i,k) = 0._r8
           loc_microp_st%qnl(i,k)  = 0._r8
           loc_microp_st%qni(i,k)  = 0._r8
           loc_microp_st%qnr(i,k)  = 0._r8
           loc_microp_st%qns(i,k)  = 0._r8
           loc_microp_st%qng(i,k)  = 0._r8
         end if
      end do
   end do
!
!jr Set to zero things which make this routine blow up
!
   do k=1,msg
      do i=1,il2g
         rprd(i,k) = 0._r8
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do k = 1, msg+1
      do i = 1,il2g
         hsthat(i,k) = hsat(i,k)
         qsthat(i,k) = qst(i,k)
         gamhat(i,k) = gamma(i,k)
      end do
   end do
   do i = 1,il2g
      totpcp(i) = 0._r8
      totevp(i) = 0._r8
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6_r8) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_r8) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                                (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
      end do
   end do
!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   jt(:) = pver
   jto(:)= pver              
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.E6_r8
   end do
!
! find the level of minimum hsat, where detrainment starts
!

   do k = msg + 1,pver
      do i = 1,il2g
         if (hsat(i,k) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
            hmin(i) = hsat(i,k)
            j0(i) = k
         end if
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
      j0(i) = min(j0(i),pver)
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
           ! Tunable temperature perturbation (tiedke_add) was already added to parcel hu/su to 
           ! represent subgrid temperature perturbation. If PBL temperature perturbation (tpert)
           ! is used to represent subgrid temperature perturbation, tiedke_add may need to be 
           ! removed. In addition, current calculation of PBL temperature perturbation is not 
           ! accurate enough so that a new tunable parameter tpert_fac was introduced. This introduced
           ! new uncertainties into the ZM scheme. The original code of ZM scheme will be used 
           ! when tpert_fix=.true.  
           if(zm_param%tpert_fix) then
             hu(i,k) = hmn(i,mx(i)) + cp*zm_param%tiedke_add 
             su(i,k) = s(i,mx(i))   +    zm_param%tiedke_add
           else
             hu(i,k) = hmn(i,mx(i)) + cp*(zm_param%tiedke_add+zm_param%tpert_fac*tpertg(i)) 
             su(i,k) = s(i,mx(i))   +     zm_param%tiedke_add+zm_param%tpert_fac*tpertg(i)
           end if
         end if
      end do
   end do
!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5_r8* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5_r8* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5_r8* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
         end if
      end do
   end do
!
! re-initialize hmin array for ensuing calculation.
!
   do i = 1,il2g
      hmin(i) = 1.E6_r8
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
         end if
      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0._r8
         ftemp(i) = 0._r8
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k) = 0._r8
            expnum(i) = 0._r8
         else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) + &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
         end if
         if ((expdif(i) > 100._r8 .and. expnum(i) > 0._r8) .and. &
	     k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                     (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                     ftemp(i)**3 + (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ &
                     5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                     k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._r8)
            f(i,k) = min(f(i,k),0.0002_r8)
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(i,j0(i)) < 1.E-6_r8 .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
            eps(i,k) = f(i,j0(i))
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
   end do


   itnum = 1
   if (zm_microp) itnum = 2

   do iter=1, itnum

      do k = pver,msg + 1,-1
        do i = 1,il2g
           cu(i,k) = 0._r8
           if (zm_microp)  loc_microp_st%qliq(i,k) = 0._r8
           if (zm_microp)  loc_microp_st%qice(i,k) = 0._r8
           ql(i,k) = 0._r8
           frz1(i,k) = 0._r8
        end do
      end do
      do i = 1,il2g
          totpcp(i) = 0._r8
          if (zm_microp)  hu(i,jb(i)) = hmn(i,jb(i)) + cp*zm_param%tiedke_add
      end do

!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
   do i = 1,il2g
      if (eps0(i) > 0._r8) then
         mu(i,jb(i)) = 1._r8
         eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
      end if
      if (zm_microp) then
        tmplel(i) = lel(i)
      else
        tmplel(i) = jt(i)
      end if
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
          if (eps0(i) > 0._r8 .and. (k >= tmplel(i) .and. k < jb(i))) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1._r8/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._r8)/zuef(i)
            mu(i,k) = (1._r8/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._r8)/zuef(i)
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
         end if
      end do
   end do
!
   khighest = pverp
   klowest = 1
   do i=1,il2g
      khighest = min(khighest,lel(i))
      klowest = max(klowest,jb(i))
   end do

   do k = klowest-1,khighest,-1
      do i = 1,il2g
         if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0._r8) then
            if (mu(i,k) < 0.02_r8) then
               hu(i,k) = hmn(i,k)
               mu(i,k) = 0._r8
               eu(i,k) = 0._r8
               du(i,k) = mu(i,k+1)/dz(i,k)
            else
              if (zm_microp) then
                hu(i,k) = (mu(i,k+1)*hu(i,k+1) + dz(i,k)*(eu(i,k)*hmn(i,k) +   &
                            latice*frz(i,k)))/(mu(i,k)+ dz(i,k)*du(i,k))
              else
                hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)- du(i,k)*hsat(i,k))
              end if
            end if
         end if
      end do
   end do
!
!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
   do i=1,il2g
      doit(i) = .true.
      totfrz(i)= 0._r8
      do k = pver,msg + 1,-1
         totfrz(i)= totfrz(i)+ frz(i,k)*dz(i,k)
      end do
   end do
   do k=klowest-2,khighest-1,-1
      do i=1,il2g
         if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
            if (hu(i,k) <= hsthat(i,k) .and. hu(i,k+1) > hsthat(i,k+1) .and. mu(i,k) >= mu_min) then
               if (hu(i,k)-hsthat(i,k) < -2000._r8) then
                  jt(i) = k + 1
                  doit(i) = .false.
               else
                  jt(i) = k
                  doit(i) = .false.
               end if
             else if ( (hu(i,k) > hu(i,jb(i)) .and. totfrz(i)<=0._r8) .or. mu(i,k) < mu_min) then
               jt(i) = k + 1
               doit(i) = .false.
            end if
         end if
      end do
   end do

   do i = 1,il2g
      if (iter == 1)  jto(i) = jt(i)
   end do

   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0._r8) then
            mu(i,k) = 0._r8
            eu(i,k) = 0._r8
            du(i,k) = 0._r8
            hu(i,k) = hmn(i,k)
         end if
         if (k == jt(i) .and. eps0(i) > 0._r8) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            eu(i,k) = 0._r8
            mu(i,k) = 0._r8
         end if
      end do
   end do
!
   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k == jb(i) .and. eps0(i) > 0._r8) then
            qu(i,k) = q(i,mx(i))
            su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
         end if
         if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0._r8) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                      dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                            du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cp*zf(i,k)
            call qsat_hPa(tu, (p(i,k)+p(i,k-1))/2._r8, estu, qstu)
            if (qu(i,k) >= qstu) then
               jlcl(i) = k
               kount = kount + 1
               done(i) = .true.
            end if
         end if
      end do
      if (kount >= il2g) goto 690
   end do
690 continue
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1._r8+gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                     (rl* (1._r8+gamhat(i,k)))
         end if
      end do
   end do

   do i = 1,il2g
     if (zm_microp) then
        tmplel(i) = jlcl(i)+1
     else
        tmplel(i) = jb(i)
     end if
   end do

! compute condensation in updraft
   do k = pver,msg + 2,-1
      do i = 1,il2g
        if (k >= jt(i) .and. k < tmplel(i) .and. eps0(i) > 0._r8) then
          if (zm_microp) then
             cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                     dz(i,k)- eu(i,k)*s(i,k)+du(i,k)*su(i,k))/(rl/cp)  &
                       - latice*frz(i,k)/rl
          else

            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                      dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/(rl/cp)
          end if
            if (k == jt(i)) cu(i,k) = 0._r8
            cu(i,k) = max(0._r8,cu(i,k))
         end if
      end do
   end do

   if (zm_microp) then

      tug(:il2g,:) = t(:il2g,:)
      fice(:,:)    = 0._r8

      do k = pver, msg+2, -1
         do i = 1, il2g
            tug(i,k) = su(i,k) - grav/cp*zf(i,k)
         end do
      end do
      do k = 1, pver-1
         do i = 1, il2g

            if (tug(i,k+1) > tmelt) then
               ! If warmer than tmelt then water phase
               fice(i,k) = 0._r8

            else if (tug(i,k+1) < t_homofrz) then
               ! If colder than t_homofrz then ice phase
               fice(i,k) = 1._r8

            else
               ! Otherwise mixed phase, with ice fraction decreasing linearly
               ! from t_homofrz to tmelt
               fice(i,k) =(tmelt - tug(i,k+1)) / t_mphase
            end if
         end do
      end do

      do k = 1, pver
         do i = 1,il2g
            loc_microp_st%cmei(i,k) = cu(i,k)* fice(i,k)
            loc_microp_st%cmel(i,k) = cu(i,k) * (1._r8-fice(i,k))
         end do
      end do

      call  zm_mphy(su,   qu,    mu,    du,   eu,  loc_microp_st%cmel,  loc_microp_st%cmei,   zf,    p,   & 
                     t,    q,  eps0,  jb,  jt,  jlcl,  msg,   il2g,  grav,   cp,    rd,   aero, gamhat,   &
                    loc_microp_st%qliq,   loc_microp_st%qice,    loc_microp_st%qnl,   loc_microp_st%qni,  &
                    qcde,  qide,  ncde,  nide,  rprd, sprd, frz, loc_microp_st%wu,  loc_microp_st%qrain,  &                 
                    loc_microp_st%qsnow,  loc_microp_st%qnr,  loc_microp_st%qns, loc_microp_st%qgraupel,  &
                    loc_microp_st%qng,   qsde,   nsde,    loc_microp_st%autolm,    loc_microp_st%accrlm,  &
                    loc_microp_st%bergnm, loc_microp_st%fhtimm, loc_microp_st%fhtctm, loc_microp_st%fhmlm,   &
                    loc_microp_st%hmpim,  loc_microp_st%accslm, loc_microp_st%dlfm,  loc_microp_st%autoln,   &
                    loc_microp_st%accrln, loc_microp_st%bergnn, loc_microp_st%fhtimn, loc_microp_st%fhtctn,  &
                    loc_microp_st%fhmln,  loc_microp_st%accsln, loc_microp_st%activn, loc_microp_st%dlfn,    &
                    loc_microp_st%autoim, loc_microp_st%accsim, loc_microp_st%difm,   loc_microp_st%nuclin,  &
                    loc_microp_st%autoin, loc_microp_st%accsin, loc_microp_st%hmpin,  loc_microp_st%difn,    &
                    loc_microp_st%trspcm, loc_microp_st%trspcn, loc_microp_st%trspim, loc_microp_st%trspin,  &
                    lambdadpcu, mudpcu,   &
                    loc_microp_st%accgrm, loc_microp_st%accglm, loc_microp_st%accgslm,loc_microp_st%accgsrm, &
                    loc_microp_st%accgirm,loc_microp_st%accgrim,loc_microp_st%accgrsm,loc_microp_st%accgsln, &
                    loc_microp_st%accgsrn,loc_microp_st%accgirn,loc_microp_st%accsrim,loc_microp_st%acciglm, &
                    loc_microp_st%accigrm,loc_microp_st%accsirm,loc_microp_st%accigln,loc_microp_st%accigrn, &
                    loc_microp_st%accsirn,loc_microp_st%accgln ,loc_microp_st%accgrn ,loc_microp_st%accilm , &
                    loc_microp_st%acciln ,loc_microp_st%fallrm ,loc_microp_st%fallsm ,loc_microp_st%fallgm , &
                    loc_microp_st%fallrn ,loc_microp_st%fallsn ,loc_microp_st%fallgn ,loc_microp_st%fhmrm  , &
                    dsfm,   dsfn, auto_fac, accr_fac, micro_dcs) 


      do k = pver,msg + 2,-1
         do i = 1,il2g
            ! In the original ZM scheme, which does not consider ice phase, ql actually represents total cloud 
            ! water. With convective microphysics, loc_microp_st%qliq and loc_microp_st%qice represent cloud
            ! liquid water and cloud ice, respectively. Since ql is still used in other subroutines as total
            ! cloud water, here ql is calculated as total cloud water for consistency.   
            ql(i,k) = loc_microp_st%qliq(i,k)+ loc_microp_st%qice(i,k)
            frz1(i,k) = frz(i,k)         
         end do
      end do

      do i = 1,il2g
        if (iter == 2 .and. jt(i)> jto(i)) then
          do k = jt(i), jto(i), -1
             frz1(i,k) = 0.0_r8
             cu(i,k)=0.0_r8
          end do
        end if
      end do

      do k = pver,msg + 2,-1
         do i = 1,il2g
            if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
               totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*(qcde(i,k+1)+qide(i,k+1)+qsde(i,k+1) ))
            end if
         end do
      end do

   else  ! no microphysics

! compute condensed liquid, rain production rate
! accumulate total precipitation (condensation - detrainment of liquid)
! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
! consistently applied.
!    mu, ql are interface quantities
!    cu, du, eu, rprd are midpoint quantites
     do k = pver,msg + 2,-1
       do i = 1,il2g
         rprd(i,k) = 0._r8
         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
            if (mu(i,k) > 0._r8) then
               ql1 = 1._r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                     dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
               ql(i,k) = ql1/ (1._r8+dz(i,k)*c0mask(i))
            else
               ql(i,k) = 0._r8
            end if
            totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))
            rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
            ! reset convective microphysics variables
            qcde(i,k) = ql(i,k)
            qide(i,k) = 0._r8
            qsde(i,k) = 0._r8
            ncde(i,k) = 0._r8
            nide(i,k) = 0._r8
            nsde(i,k) = 0._r8
            sprd(i,k) = 0._r8
            frz1(i,k) = 0._r8
         end if
       end do
     end do
   end if  ! zm_microp

 end do   !iter

! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
      alfa(i) = alfa_scalar
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0._r8) then
         epsm(i) = eps0(i)
         md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
      end if
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i)
         end if
      end do
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)
         end if
      end do
   end do

   small = 1.e-20_r8
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0._r8) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
         end if
      end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ &
               (rl*(1._r8 + gamhat(i,k)))
         end if
      end do
   end do
!
   do i = 1,il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            mdt = min(md(i,k+1),-small)
            if (zm_microp)   evp(i,k) = min(evp(i,k),rprd(i,k))
            sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do
   do i = 1,il2g
!*guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
   end do
!!$   if (.true.) then
   if (.false.) then
      do i = 1,il2g
         k = jb(i)
         if (eps0(i) > 0._r8) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   endif

   do i = 1,il2g
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0._r8 .and. totpcp(i) > 0._r8) then
            md(i,k)  = md (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(i,k)  = ed (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(i,k) = 0._r8
            ed(i,k) = 0._r8
            evp(i,k) = 0._r8
         end if
! cmeg is the cloud water condensed - rain water evaporated
! rprd is the cloud water converted to rain - (rain evaporated)
         cmeg(i,k) = cu(i,k) - evp(i,k)
         if (zm_microp) then
           if (rprd(i,k)> 0._r8)  then
              frz1(i,k) = frz1(i,k)- evp(i,k)*min(1._r8,sprd(i,k)/rprd(i,k))
              sprd(i,k) = sprd(i,k)- evp(i,k)*min(1._r8,sprd(i,k)/rprd(i,k))
           end if
         end if
         rprd(i,k) = rprd(i,k)-evp(i,k)
      end do
   end do

! compute the net precipitation flux across interfaces
   pflx(:il2g,1) = 0._r8
   if (zm_microp)   pflxs(:,:) = 0._r8
   do k = 2,pverp
      do i = 1,il2g
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
         if (zm_microp) pflxs(i,k) = pflxs(i,k-1) + sprd(i,k-1)*dz(i,k-1)
      end do
   end do
! protect against rounding error
   if (zm_microp) then
     do i = 1,il2g
       if(pflxs(i,pverp).gt.pflx(i,pverp)) then
         dum = (pflxs(i,pverp)-pflx(i,pverp))/omsm
         do k = pver, msg+2, -1
           if (sprd(i,k) > 0._r8 .and. dum > 0._r8) then
             sdum = min(sprd(i,k),dum/dz(i,k))
             sprd(i,k) = sprd(i,k)- sdum
             frz1(i,k) = frz1(i,k)- sdum
             dum = dum - sdum*dz(i,k)
           end if
         end do
       end if
     end do   
   end if

   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do
!
   do i = 1,il2g
     if ( zm_microp .and. jt(i)>=jlcl(i)) then
       do k = msg + 1,pver  
          mu(i,k)   = 0._r8
          eu(i,k)   = 0._r8
          du(i,k)   = 0._r8
          ql(i,k)   = 0._r8
          cu(i,k)   = 0._r8
          evp(i,k)  = 0._r8
          cmeg(i,k) = 0._r8
          md(i,k)   = 0._r8
          ed(i,k)   = 0._r8
          mc(i,k)   = 0._r8
          rprd(i,k) = 0._r8
          sprd(i,k) = 0._r8
          fice(i,k) = 0._r8
          qcde(i,k) = 0._r8
          qide(i,k) = 0._r8
          qsde(i,k) = 0._r8
          ncde(i,k) = 0._r8
          nide(i,k) = 0._r8
          nsde(i,k) = 0._r8
          frz(i,k)  = 0._r8
          frz1(i,k) = 0._r8
          loc_microp_st%wu(i,k)   = 0._r8      
          loc_microp_st%cmel(i,k) = 0._r8
          loc_microp_st%cmei(i,k) = 0._r8
          loc_microp_st%qliq(i,k) = 0._r8
          loc_microp_st%qice(i,k) = 0._r8
          loc_microp_st%qrain(i,k)= 0._r8
          loc_microp_st%qsnow(i,k)= 0._r8
          loc_microp_st%qgraupel(i,k) = 0._r8
          loc_microp_st%qnl(i,k)  = 0._r8
          loc_microp_st%qni(i,k)  = 0._r8
          loc_microp_st%qnr(i,k)  = 0._r8
          loc_microp_st%qns(i,k)  = 0._r8
          loc_microp_st%qng(i,k)  = 0._r8
       end do
     end if
   end do       
   return
end subroutine cldprp

subroutine closure(lchnk   , &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,capelmt )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
   use dycore,    only: dycore_is, get_resolution

   implicit none

!
!-----------------------------Arguments---------------------------------
!
   integer, intent(in) :: lchnk                 ! chunk identifier

   real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
   real(r8), intent(inout) :: t(pcols,pver)        ! temperature
   real(r8), intent(inout) :: p(pcols,pver)        ! pressure (mb)
   real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
   real(r8), intent(in) :: z(pcols,pver)        ! height (m)
   real(r8), intent(in) :: s(pcols,pver)        ! normalized dry static energy
   real(r8), intent(in) :: tp(pcols,pver)       ! parcel temp
   real(r8), intent(in) :: qs(pcols,pver)       ! sat spec humidity
   real(r8), intent(in) :: qu(pcols,pver)       ! updraft spec. humidity
   real(r8), intent(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   real(r8), intent(in) :: mc(pcols,pver)       ! net convective mass flux
   real(r8), intent(in) :: du(pcols,pver)       ! detrainment from updraft
   real(r8), intent(in) :: mu(pcols,pver)       ! mass flux of updraft
   real(r8), intent(in) :: md(pcols,pver)       ! mass flux of downdraft
   real(r8), intent(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
   real(r8), intent(in) :: sd(pcols,pver)       ! dry static energy of downdraft
   real(r8), intent(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
   real(r8), intent(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
   real(r8), intent(in) :: dp(pcols,pver)       ! pressure thickness of layers
   real(r8), intent(in) :: qstp(pcols,pver)     ! spec humidity of parcel
   real(r8), intent(in) :: zf(pcols,pver+1)     ! height of interface levels
   real(r8), intent(in) :: ql(pcols,pver)       ! liquid water mixing ratio

   real(r8), intent(in) :: cape(pcols)          ! available pot. energy of column
   real(r8), intent(in) :: tl(pcols)
   real(r8), intent(in) :: dsubcld(pcols)       ! thickness of subcloud layer

   integer, intent(in) :: lcl(pcols)        ! index of lcl
   integer, intent(in) :: lel(pcols)        ! index of launch leve
   integer, intent(in) :: jt(pcols)         ! top of updraft
   integer, intent(in) :: mx(pcols)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real(r8) dtpdt(pcols,pver)
   real(r8) dqsdtp(pcols,pver)
   real(r8) dtmdt(pcols,pver)
   real(r8) dqmdt(pcols,pver)
   real(r8) dboydt(pcols,pver)
   real(r8) thetavp(pcols,pver)
   real(r8) thetavm(pcols,pver)

   real(r8) dtbdt(pcols),dqbdt(pcols),dtldt(pcols)
   real(r8) beta
   real(r8) capelmt
   real(r8) cp
   real(r8) dadt(pcols)
   real(r8) debdt
   real(r8) dltaa
   real(r8) eb
   real(r8) grav

   integer i
   integer il1g
   integer il2g
   integer k, kmin, kmax
   integer msg

   real(r8) rd
   real(r8) rl
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
      mb(i) = 0._r8
      eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      dtbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
      dqbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
                 md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
      debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
      dtldt(i) = -2840._r8* (3.5_r8/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t(i,mx(i)))-log(eb)-4.805_r8)**2
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(i,k) = 0._r8
         dqmdt(i,k) = 0._r8
      end do
   end do
!
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
            dqmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                         qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
         end if
      end do
   end do
!
   beta = 0._r8
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                         dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
            dqsdtp(i,k) = qstp(i,k)* (1._r8+qstp(i,k)/eps1)*eps1*rl/(rd*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(i,k) = tp(i,k)/ (1._r8+rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                        (dtbdt(i)/t(i,mx(i))+rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                         tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._r8/(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))* &
                          (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_r8/ &
                          (1._r8+0.608_r8*q(i,k))*dqmdt(i,k)))*grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+0.608_r8/ (1._r8+0.608_r8*q(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t(i,k)-0.608_r8/ (1._r8+0.608_r8*q(i,k))*dqmdt(i,k))* &
                          grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   dadt(il1g:il2g)  = 0._r8
   kmin = minval(lel(il1g:il2g))
   kmax = maxval(mx(il1g:il2g)) - 1
   do k = kmin, kmax
      do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         endif
      end do
   end do
   do i = il1g,il2g
      dltaa = -1._r8* (cape(i)-capelmt)
      if (dadt(i) /= 0._r8) mb(i) = max(dltaa/tau/dadt(i),0._r8)
      if (zm_microp .and. mx(i)-jt(i) < 2._r8) mb(i) =0.0_r8
   end do
!
   return
end subroutine closure

subroutine q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
                    dl      ,evp     ,cu      ,qice    ,dice    , &
                    qnl     ,qni     ,dnl     ,dni     ,frz     , &
                    qsde    ,nsde    ,dsnow   ,dns    )

   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
! 
!-----------------------------------------------------------------------


   real(r8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   real(r8), intent(in) :: q(pcols,pver)
   real(r8), intent(in) :: qs(pcols,pver)
   real(r8), intent(in) :: qu(pcols,pver)
   real(r8), intent(in) :: su(pcols,pver)
   real(r8), intent(in) :: du(pcols,pver)
   real(r8), intent(in) :: qhat(pcols,pver)
   real(r8), intent(in) :: shat(pcols,pver)
   real(r8), intent(in) :: dp(pcols,pver)
   real(r8), intent(in) :: mu(pcols,pver)
   real(r8), intent(in) :: md(pcols,pver)
   real(r8), intent(in) :: sd(pcols,pver)
   real(r8), intent(in) :: qd(pcols,pver)
   real(r8), intent(in) :: ql(pcols,pver)
   real(r8), intent(in) :: evp(pcols,pver)
   real(r8), intent(in) :: cu(pcols,pver)
   real(r8), intent(in) :: dsubcld(pcols)

   ! Convective microphysics
   real(r8), intent(in) :: frz(pcols,pver)
   real(r8), intent(in) :: qice(pcols,pver)
   real(r8), intent(in) :: qnl(pcols,pver)
   real(r8), intent(in) :: qni(pcols,pver)
   real(r8), intent(in) :: qsde(pcols,pver)
   real(r8), intent(in) :: nsde(pcols,pver)

   real(r8),intent(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
   real(r8),intent(out) :: dl(pcols,pver)

   ! Convective microphysics
   real(r8),intent(out) :: dice(pcols,pver)
   real(r8),intent(out) :: dnl(pcols,pver)
   real(r8),intent(out) :: dni(pcols,pver)
   real(r8),intent(out) :: dsnow(pcols,pver)
   real(r8),intent(out) :: dns(pcols,pver)


   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
!
! work fields:
!
   integer i
   integer k

   real(r8) emc
   real(r8) rl
!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0._r8
         dqdt(i,k) = 0._r8
         dl(i,k) = 0._r8
         ! Convective microphysics
         dice(i,k) = 0._r8
         dnl(i,k)  = 0._r8
         dni(i,k)  = 0._r8
         dsnow(i,k) = 0._r8
         dns(i,k)  = 0._r8
      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         if (zm_microp) dsdt(i,k) = dsdt(i,k) + latice/cp*frz(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = du(i,k)*ql(i,k+1)
         if (zm_microp) then
           dice(i,k) = du(i,k)*qice(i,k+1)
           dnl(i,k)  = du(i,k)*qnl(i,k+1)
           dni(i,k)  = du(i,k)*qni(i,k+1)
           dsnow(i,k) = du(i,k)*qsde(i,k+1)
           dns(i,k)  = du(i,k)*nsde(i,k+1)
         end if

      end do
   end do

!
#ifdef CPRCRAY
!DIR$ NOINTERCHANGE!
#endif
   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        )
            dqdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        )
         else if (k > mx(i)) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do
!
   return
end subroutine q1q2_pjr

end module zm_conv
