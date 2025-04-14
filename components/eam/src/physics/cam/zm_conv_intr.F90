module zm_conv_intr
   !----------------------------------------------------------------------------
   ! 
   ! Interface to the Zhang-McFarlane deep convection scheme
   !
   ! Author: D.B. Coleman
   ! January 2010 modified by J. Kay to add COSP simulator fields to pbuf
   ! July 2015 B. Singh Added code for unified convective trasport
   !----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8=>shr_kind_r8
   use spmd_utils,            only: masterproc
   use perf_mod,              only: t_startf, t_stopf
   use cam_abortutils,        only: endrun
   use cam_history,           only: outfld, addfld, horiz_only, add_default
   use cam_logfile,           only: iulog
   use physconst,             only: pi, cpair
   use ppgrid,                only: pver, pcols, pverp, begchunk, endchunk
   use rad_constituents,      only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr
   use rad_constituents,      only: rad_cnst_get_aer_props, rad_cnst_get_mode_props
   use ndrop_bam,             only: ndrop_bam_init
   use zm_conv,               only: zm_conv_evap, zm_convr, trigdcape_ull, trig_dcape_only
   use zm_conv,               only: MCSP, MCSP_heat_coeff, MCSP_moisture_coeff, MCSP_uwind_coeff, MCSP_vwind_coeff
   use zm_conv,               only: zm_microp
   use zm_transport,          only: zm_transport_tracer, zm_transport_momentum
   use zm_aero,               only: zm_aero_t
   use zm_microphysics_state, only: zm_microp_st

   implicit none
   private
   save

   ! public methods
   public :: zm_conv_register  ! register fields in physics buffer
   public :: zm_conv_init      ! initialize donner_deep module
   public :: zm_conv_tend      ! return tendencies
   public :: zm_conv_tend_2    ! return tendencies

   ! physics buffer field indices
   integer :: dp_flxprc_idx    ! deep conv flux of precipitation from deep convection (kg/m2/s)
   integer :: dp_flxsnw_idx    ! deep conv flux of snow from deep convection (kg/m2/s)
   integer :: dp_cldliq_idx    ! deep conv cloud liq water (kg/kg)
   integer :: dp_cldice_idx    ! deep conv cloud ice water (kg/kg)
   integer :: dlfzm_idx        ! detrained convective cloud water mixing ratio
   integer :: difzm_idx        ! detrained convective cloud ice mixing ratio
   integer :: prec_dp_idx      ! total surface precipitation rate from deep conv
   integer :: snow_dp_idx      ! frozen surface precipitation rate from deep conv
   integer :: t_star_idx       ! DCAPE temperature from previous time step
   integer :: q_star_idx       ! DCAPE water vapor from previous time step
   integer :: cld_idx          ! cloud fraction
   integer :: icwmrdp_idx      ! in-cloud water mixing ratio
   integer :: rprddp_idx       ! rain production
   integer :: fracis_idx       ! fraction of transported species that are insoluble
   integer :: nevapr_dpcu_idx  ! evaporation of deep conv precipitation
   integer :: dgnum_idx        ! dry aerosol mode diameter
   integer :: lambdadpcu_idx   ! slope of cloud liquid size distribution
   integer :: mudpcu_idx       ! width parameter of droplet size distr
   integer :: icimrdp_idx      ! in-cloud ice mixing ratio

   ! other private module data
   logical :: convproc_do_aer 
   logical :: convproc_do_gas 
   logical :: clim_modal_aero
   logical :: old_snow  = .true.   ! flag to use old estimate of snow production in zm_conv_evap
                                   ! set false to use snow production from zm microphysics
   integer :: nmodes
   integer :: nbulk
   type(zm_aero_t), allocatable :: aero(:) ! object contains aerosol information
   
   real(r8), parameter :: ZM_upper_limit_pref   = 40e2_r8  ! pressure limit above which deep convection is not allowed [Pa]
   real(r8), parameter :: MCSP_storm_speed_pref = 600e2_r8 ! pressure level for winds in MCSP calculation [Pa]
   real(r8), parameter :: MCSP_conv_depth_min   = 700e2_r8 ! pressure thickness of convective heating [Pa]
   real(r8), parameter :: MCSP_shear_min        = 3.0_r8   ! min shear value for MCSP to be active
   real(r8), parameter :: MCSP_shear_max        = 200.0_r8 ! max shear value for MCSP to be active

!===================================================================================================
contains
!===================================================================================================

subroutine zm_conv_register
   !----------------------------------------------------------------------------
   ! Purpose: register fields with the physics buffer
   !----------------------------------------------------------------------------
   use physics_buffer,  only: pbuf_add_field, dtype_r8
   use misc_diagnostics,only: dcape_diags_register
   use zm_microphysics, only: zm_microphysics_register
   implicit none

   integer idx

   ! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx)

   ! Flux of snow from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx)

   ! deep conv cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)

   ! deep conv cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)

   ! previous time step data for DCAPE calculation
   if (trigdcape_ull .or. trig_dcape_only) then
      call pbuf_add_field('T_STAR','global',dtype_r8,(/pcols,pver/), t_star_idx)
      call pbuf_add_field('Q_STAR','global',dtype_r8,(/pcols,pver/), q_star_idx)
   end if

   ! detrained convective cloud water mixing ratio.
   call pbuf_add_field('DLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dlfzm_idx)
   ! detrained convective cloud ice mixing ratio.
   call pbuf_add_field('DIFZM', 'physpkg', dtype_r8, (/pcols,pver/), difzm_idx)

   ! Only add the number conc fields if the microphysics is active.
   if (zm_microp) call zm_microphysics_register()

   ! Register variables for dCAPE diagnosis and decomposition
   call dcape_diags_register( pcols )

end subroutine zm_conv_register

!===================================================================================================

subroutine zm_conv_init(pref_edge)
   !----------------------------------------------------------------------------
   ! Purpose: declare output fields, initialize variables needed by convection
   !----------------------------------------------------------------------------
   use zm_conv,                 only: zm_convi
   use pmgrid,                  only: plev,plevp
   use spmd_utils,              only: masterproc
   use error_messages,          only: alloc_err
   use phys_control,            only: phys_deepconv_pbl, phys_getopts
   use physics_buffer,          only: pbuf_get_index
   use rad_constituents,        only: rad_cnst_get_info
   use zm_microphysics,         only: zm_mphyi
   use zm_aero,                 only: zm_aero_init
   use zm_microphysics_history, only: zm_microphysics_history_init
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces
   !----------------------------------------------------------------------------
   ! Local variables
   logical :: no_deep_pbl                 ! if true, no deep convection in PBL
   integer :: limcnv                      ! top interface level limit for convection
   logical :: history_budget              ! output tendencies and state variables for 
                                          ! temperature, water vapor, cloud ice/liq budgets
   integer :: history_budget_histfile_num ! output history file number for budget fields
   integer i, k, istat
   character(len=*), parameter :: routine = 'zm_conv_init'
   !----------------------------------------------------------------------------
   ! Allocate the basic aero structure outside the zmconv_microp logical
   ! This allows the aero structure to be passed
   ! Note that all of the arrays inside this structure are conditionally allocated
   allocate(aero(begchunk:endchunk))

   ! Register fields with the output buffer

   call addfld('PRECZ',        horiz_only, 'A', 'm/s',      'ZM total precipitation rate')
   call addfld('ZMDT',         (/ 'lev'/), 'A', 'K/s',      'ZM T tendency')
   call addfld('ZMDQ',         (/ 'lev'/), 'A', 'kg/kg/s',  'ZM Q tendency')
   call addfld('ZMDICE',       (/ 'lev'/), 'A', 'kg/kg/s',  'ZM cloud ice tendency')
   call addfld('ZMDLIQ',       (/ 'lev'/), 'A', 'kg/kg/s',  'ZM cloud liq tendency')
   call addfld('EVAPTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from evaporation/snow prod')
   call addfld('FZSNTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from rain to snow conversion')
   call addfld('EVSNTZM',      (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from snow to rain prod')
   call addfld('EVAPQZM',      (/ 'lev'/), 'A', 'kg/kg/s',  'ZM Q tendency from evaporation')
   call addfld('ZMFLXPRC',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM flux of precipitation')
   call addfld('ZMFLXSNW',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM flux of snow')
   call addfld('ZMNTPRPD',     (/ 'lev'/), 'A', 'kg/kg/s',  'ZM net precipitation production')
   call addfld('ZMNTSNPD',     (/ 'lev'/), 'A', 'kg/kg/s',  'ZM net snow production')
   call addfld('ZMEIHEAT',     (/ 'lev'/), 'A', 'W/kg',     'ZM heating by ice and evaporation')
   call addfld('CMFMCDZM',     (/'ilev'/), 'A', 'kg/m2/s',  'ZM convection mass flux')
   call addfld('PRECCDZM',     horiz_only, 'A', 'm/s',      'ZM convective precipitation rate')
   call addfld('PCONVB',       horiz_only, 'A', 'Pa',       'ZM convection base pressure')
   call addfld('PCONVT',       horiz_only, 'A', 'Pa',       'ZM convection top  pressure')
   call addfld('MAXI',         horiz_only, 'A', 'level',    'ZM model level of launching parcel')
   call addfld('CAPE_ZM',      horiz_only, 'A', 'J/kg',     'ZM convectively available potential energy')
   call addfld('DCAPE',        horiz_only, 'A', 'J/kg',     'ZM rate of change of CAPE')
   call addfld('FREQZM',       horiz_only, 'A', 'fraction', 'ZM fractional occurrence of convection')
   call addfld('ZMMTT',        (/ 'lev'/), 'A', 'K/s',      'ZM T tendency from convective momentum transport')
   call addfld('ZMMTU',        (/ 'lev'/), 'A', 'm/s2',     'ZM U tendency from convective momentum transport')
   call addfld('ZMMTV',        (/ 'lev'/), 'A', 'm/s2',     'ZM V tendency from convective momentum transport')
   call addfld('ZMMU',         (/ 'lev'/), 'A', 'kg/m2/s',  'ZM convection updraft mass flux')
   call addfld('ZMMD',         (/ 'lev'/), 'A', 'kg/m2/s',  'ZM convection downdraft mass flux')
   call addfld('ZMUPGU',       (/ 'lev'/), 'A', 'm/s2',     'ZM zonal force from updraft pressure gradient term')
   call addfld('ZMUPGD',       (/ 'lev'/), 'A', 'm/s2',     'ZM zonal force from downdraft pressure gradient term')
   call addfld('ZMVPGU',       (/ 'lev'/), 'A', 'm/s2',     'ZM meridional force from updraft pressure gradient term')
   call addfld('ZMVPGD',       (/ 'lev'/), 'A', 'm/s2',     'ZM merdional force from downdraft pressure gradient term')
   call addfld('ZMICUU',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud U updrafts')
   call addfld('ZMICUD',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud U downdrafts')
   call addfld('ZMICVU',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud V updrafts')
   call addfld('ZMICVD',       (/ 'lev'/), 'A', 'm/s',      'ZM in-cloud V downdrafts')

   if (MCSP) then 
      call addfld('MCSP_DT',   (/ 'lev'/), 'A', 'K/s',      'MCSP T tendency')
      call addfld('MCSP_freq', horiz_only, 'A', '1',        'MCSP frequency of activation')
      call addfld('MCSP_DU',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP U tendency')
      call addfld('MCSP_DV',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP V tendency')
      call addfld('MCSP_shear',horiz_only, 'A', 'm/s',      'MCSP vertical zonal wind shear')
      call addfld('ZM_depth',  horiz_only, 'A', 'Pa',       'ZM convection depth')
   end if

   call phys_getopts( history_budget_out = history_budget, &
                      history_budget_histfile_num_out = history_budget_histfile_num, &
                      convproc_do_aer_out = convproc_do_aer, &
                      convproc_do_gas_out = convproc_do_gas)
   
   ! Determine whether modal aerosols are used
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   if ( history_budget ) then
      call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
      call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
      call add_default('ZMDT     ', history_budget_histfile_num, ' ')
      call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
      call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
      call add_default('ZMDICE   ', history_budget_histfile_num, ' ')
      call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
   end if

   ! Limit deep convection to regions below ZM_upper_limit_pref
   limcnv = 0 ! initialize to null value to check against below
   if (pref_edge(1) >= ZM_upper_limit_pref) then
      limcnv = 1
   else
      do k = 1,plev
         if (pref_edge(k)   <  ZM_upper_limit_pref .and. &
             pref_edge(k+1) >= ZM_upper_limit_pref) then
            limcnv = k
            exit
         end if
      end do
      if ( limcnv == 0 ) limcnv = plevp
   end if
    
   if (masterproc) write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at ', &
                                 'intfc ',limcnv,' which is ',pref_edge(limcnv),' pascals'

   no_deep_pbl = phys_deepconv_pbl()
   call zm_convi( limcnv, no_deep_pbl_in=no_deep_pbl )

   dp_flxprc_idx   = pbuf_get_index('DP_FLXPRC')
   dp_flxsnw_idx   = pbuf_get_index('DP_FLXSNW')
   dp_cldliq_idx   = pbuf_get_index('DP_CLDLIQ')
   dp_cldice_idx   = pbuf_get_index('DP_CLDICE')
   cld_idx         = pbuf_get_index('CLD')
   icwmrdp_idx     = pbuf_get_index('ICWMRDP')
   rprddp_idx      = pbuf_get_index('RPRDDP')
   fracis_idx      = pbuf_get_index('FRACIS')
   nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
   prec_dp_idx     = pbuf_get_index('PREC_DP')
   snow_dp_idx     = pbuf_get_index('SNOW_DP')
   lambdadpcu_idx  = pbuf_get_index('LAMBDADPCU')
   mudpcu_idx      = pbuf_get_index('MUDPCU')
   icimrdp_idx     = pbuf_get_index('ICIMRDP')

   ! Initialization for convective microphysics
   if (zm_microp) then

      call zm_microphysics_history_init()

      call zm_mphyi()

      ! use old estimate of snow production in zm_conv_evap
      old_snow = .false. 

      ! Initialize the aerosol object with data from the modes/species
      ! affecting climate, i.e., the list index is hardcoded to 0
      call rad_cnst_get_info( 0, nmodes=nmodes, naero=nbulk )

      do i = begchunk, endchunk
         call zm_aero_init(nmodes, nbulk, aero(i))
      end do

      if (nmodes > 0) then
         dgnum_idx = pbuf_get_index('DGNUM')
      else if (nbulk > 0) then
         call ndrop_bam_init()
      end if

   end if ! zmconv_microp

end subroutine zm_conv_init

!===================================================================================================

subroutine zm_conv_tend(pblh, mcon, cme, tpert, dlftot, pflx, zdu, &
                        rliq, rice, ztodt, jctop, jcbot, &
                        state, ptend_all, landfrac, pbuf, mu, eu, &
                        du, md, ed, dp, dsubcld, jt, maxg, ideep, lengath )
   !----------------------------------------------------------------------------
   ! Purpose: Primary interface with ZM parameterization
   !----------------------------------------------------------------------------
   use physics_types,           only: physics_state, physics_ptend
   use physics_types,           only: physics_ptend_init
   use physics_update_mod,      only: physics_update
   use physics_types,           only: physics_state_copy, physics_state_dealloc
   use physics_types,           only: physics_ptend_sum, physics_ptend_dealloc
   use phys_grid,               only: get_lat_p, get_lon_p
   use time_manager,            only: get_nstep, is_first_step
   use physics_buffer,          only: pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,            only: pcnst, cnst_get_ind, cnst_is_convtran1
   use physconst,               only: gravit
   use time_manager,            only: get_curr_date
   use interpolate_data,        only: vertinterp
   use zm_microphysics,         only: dnlfzm_idx, dnifzm_idx, dsfzm_idx, dnsfzm_idx, wuc_idx
   use zm_microphysics_state,   only: zm_microp_st_alloc, zm_microp_st_dealloc
   use zm_microphysics_history, only: zm_microphysics_history_out
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),target,       intent(in)  :: state      ! Physics state variables
   type(physics_ptend),              intent(out) :: ptend_all  ! individual parameterization tendencies
   type(physics_buffer_desc),        pointer     :: pbuf(:)    ! physics buffer
   real(r8),                         intent(in)  :: ztodt      ! 2 delta t (model time increment)
   real(r8), dimension(pcols),       intent(in)  :: pblh       ! Planetary boundary layer height
   real(r8), dimension(pcols),       intent(in)  :: tpert      ! Thermal temperature excess
   real(r8), dimension(pcols),       intent(in)  :: landfrac   ! Land fraction
   real(r8), dimension(pcols,pverp), intent(out) :: mcon       ! Convective mass flux--m sub c
   real(r8), dimension(pcols,pver ), intent(out) :: dlftot     ! scattrd version of the detraining cld h2o tend
   real(r8), dimension(pcols,pverp), intent(out) :: pflx       ! scattered precip flux at each level
   real(r8), dimension(pcols,pver ), intent(out) :: cme        ! cmf condensation - evaporation
   real(r8), dimension(pcols,pver ), intent(out) :: zdu        ! detraining mass flux
   real(r8), dimension(pcols),       intent(out) :: rliq       ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), dimension(pcols),       intent(out) :: rice       ! reserved ice (not yet in cldice) for energy integrals
   real(r8), dimension(pcols,pver ), intent(out) :: mu         ! upward cloud mass flux
   real(r8), dimension(pcols,pver ), intent(out) :: eu         ! entrainment in updraft
   real(r8), dimension(pcols,pver ), intent(out) :: du         ! detrainment in updraft
   real(r8), dimension(pcols,pver ), intent(out) :: md         ! downward cloud mass flux
   real(r8), dimension(pcols,pver ), intent(out) :: ed         ! entrainment in downdraft
   real(r8), dimension(pcols,pver ), intent(out) :: dp         ! layer thickness
   real(r8), dimension(pcols),       intent(out) :: dsubcld    ! wg layer thickness in mbs (between upper/lower interface)
   integer,  dimension(pcols),       intent(out) :: jt         ! wg layer thickness in mbs between lcl and maxi
   integer,  dimension(pcols),       intent(out) :: maxg       ! wg top  level index of deep cumulus convection
   integer,  dimension(pcols),       intent(out) :: ideep      ! wg gathered values of maxi
   integer,                          intent(out) :: lengath    ! gathered points vs longitude index
   !----------------------------------------------------------------------------
   ! Local variables
   type(zm_microp_st) :: microp_st     ! ZM microphysics data structure

   integer :: i,k,l,m                  ! loop iterators
   integer :: nstep                    ! model time step number
   integer :: ixcldice, ixcldliq       ! constituent indices for cloud liquid and ice water
   integer :: lchnk                    ! chunk identifier
   integer :: ncol                     ! number of atmospheric columns
   integer :: itim_old                 ! for physics buffer fields
   logical :: lq(pcnst)                ! flags for ptend initialization

   real(r8), dimension(pcols,pver) :: ftem            ! Temporary workspace for outfld variables
   real(r8), dimension(pcols,pver) :: ntprprd         ! evap outfld: net precip production in layer
   real(r8), dimension(pcols,pver) :: ntsnprd         ! evap outfld: net snow production in layer
   real(r8), dimension(pcols,pver) :: tend_s_snwprd   ! Heating rate of snow production
   real(r8), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
   real(r8), dimension(pcols,pver) :: fake_dpdry      ! used in convtran call

   ! physics types
   type(physics_state)        :: state1               ! copy of state for evaporation
   type(physics_ptend),target :: ptend_loc            ! output tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)     :: prec        ! total precipitation
   real(r8), pointer, dimension(:)     :: snow        ! snow from ZM convection 
   real(r8), pointer, dimension(:,:)   :: cld         ! cloud fraction
   real(r8), pointer, dimension(:,:)   :: ql          ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:)   :: rprd        ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis      ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:)   :: evapcdp     ! evaporation of precipitation
   real(r8), pointer, dimension(:,:)   :: flxprec     ! convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:)   :: flxsnow     ! convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:)   :: dp_cldliq   ! cloud liq water
   real(r8), pointer, dimension(:,:)   :: dp_cldice   ! cloud ice water
   real(r8), pointer, dimension(:,:)   :: wuc         ! vertical velocity

   ! DCAPE-ULL
   real(r8), pointer, dimension(:,:) :: t_star        ! DCAPE T from time step n-1
   real(r8), pointer, dimension(:,:) :: q_star        ! DCAPE q from time step n-1
   real(r8), dimension(pcols)        :: dcape         ! DCAPE cape change
   real(r8), dimension(pcols)        :: maxgsav       ! DCAPE tmp array for recording and outfld to MAXI

   real(r8), pointer, dimension(:,:) :: dlf           ! detrained convective cloud water mixing ratio
   real(r8), pointer, dimension(:,:) :: dif           ! detrained convective cloud ice mixing ratio
   real(r8), pointer, dimension(:,:) :: dsf           ! detrained convective snow mixing ratio
   real(r8), pointer, dimension(:,:) :: dnlf          ! detrained convective cloud water num concen
   real(r8), pointer, dimension(:,:) :: dnif          ! detrained convective cloud ice num concen
   real(r8), pointer, dimension(:,:) :: dnsf          ! detrained convective snow num concen
   real(r8), pointer, dimension(:,:) :: lambdadpcu    ! slope of cloud liquid size distr
   real(r8), pointer, dimension(:,:) :: mudpcu        ! width parameter of droplet size distr
   real(r8), pointer, dimension(:,:) :: qi            ! grid slice of cloud ice

   integer , dimension(pcols) :: jctop                ! row of top-of-deep-convection indices passed out
   integer , dimension(pcols) :: jcbot                ! row of base of cloud indices passed out
   real(r8), dimension(pcols) :: pcont                ! convection top pressure
   real(r8), dimension(pcols) :: pconb                ! convection base pressure
   real(r8), dimension(pcols) :: freqzm               ! fractional occurrence of ZM convection

   ! history output fields
   real(r8), dimension(pcols)      :: cape            ! convective available potential energy
   real(r8), dimension(pcols,pver) :: mu_out          ! updraft mass flux for output
   real(r8), dimension(pcols,pver) :: md_out          ! downdraft mass flux for output

   ! used in momentum transport calculations
   real(r8), dimension(pcols,pver,2) :: winds
   real(r8), dimension(pcols,pver,2) :: wind_tends
   real(r8), dimension(pcols,pver,2) :: pguall
   real(r8), dimension(pcols,pver,2) :: pgdall
   real(r8), dimension(pcols,pver,2) :: icwu
   real(r8), dimension(pcols,pver,2) :: icwd
   real(r8), dimension(pcols,pver)   :: seten

   real(r8), dimension(pcols,pver) :: sprd
   real(r8), dimension(pcols,pver) :: frz

   ! MCSP
   logical  :: doslop
   logical  :: doslop_heat
   logical  :: doslop_moisture
   logical  :: doslop_uwind
   logical  :: doslop_vwind
   real(r8) :: alpha2, alpha_moisture, alphau, alphav
   real(r8) :: mcsp_top, mcsp_bot
   real(r8) :: dpg
   real(r8) :: Qcq_adjust
   real(r8), dimension(pcols)      :: Q_dis
   real(r8), dimension(pcols)      :: Qq_dis
   real(r8), dimension(pcols,pver) :: Qm
   real(r8), dimension(pcols,pver) :: Qmq
   real(r8), dimension(pcols,pver) :: Qmu
   real(r8), dimension(pcols,pver) :: Qmv
   real(r8), dimension(pcols)      :: Qm_int_end
   real(r8), dimension(pcols)      :: Qmq_int_end
   real(r8), dimension(pcols)      :: Pa_int_end
   real(r8), dimension(pcols)      :: Qs_zmconv
   real(r8), dimension(pcols)      :: Qv_zmconv
   real(r8), dimension(pcols)      :: MCSP_freq
   real(r8), dimension(pcols,pver) :: MCSP_DT
   real(r8), dimension(pcols)      :: ZM_depth
   real(r8), dimension(pcols)      :: MCSP_shear
   real(r8), dimension(pcols)      :: du600
   real(r8), dimension(pcols)      :: dv600

   !----------------------------------------------------------------------------

   if (zm_microp) call zm_microp_st_alloc(microp_st)

   doslop          = .false.
   doslop_heat     = .false.
   doslop_moisture = .false.
   doslop_uwind    = .false.
   doslop_vwind    = .false.
   alphau          = 0
   alphav          = 0
   if ( MCSP ) then
      if ( MCSP_heat_coeff > 0 ) then
         doslop_heat = .true.
         alpha2 = MCSP_heat_coeff
      end if
      if ( MCSP_moisture_coeff > 0 ) then
         doslop_moisture = .true.
         alpha_moisture = MCSP_moisture_coeff
      end if
      if ( MCSP_uwind_coeff > 0 ) then
         doslop_uwind = .true.
         alphau = MCSP_uwind_coeff
      end if
      if ( MCSP_vwind_coeff > 0 ) then
         doslop_vwind = .true.
         alphav = MCSP_vwind_coeff
      end if
   end if

   if (doslop_heat .or. doslop_moisture .or. doslop_uwind .or. doslop_vwind) then
      doslop = .true.
   end if

   !----------------------------------------------------------------------------
   ! Initialize stuff

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem      (1:ncol,1:pver)   = 0
   mu_out    (1:ncol,1:pver)   = 0
   md_out    (1:ncol,1:pver)   = 0
   dlftot    (1:ncol,1:pver)   = 0
   wind_tends(1:ncol,1:pver,:) = 0

   call physics_state_copy(state,state1) ! make local copy of state

   lq(:) = .FALSE.
   lq(1) = .TRUE.
 
   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq, lu=doslop_uwind, lv=doslop_vwind)

   !----------------------------------------------------------------------------
   ! Associate pointers with physics buffer fields

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,           cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, fracis_idx,        fracis, start=(/1,1,1/),        kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, icwmrdp_idx,       ql)
   call pbuf_get_field(pbuf, rprddp_idx,        rprd)
   call pbuf_get_field(pbuf, nevapr_dpcu_idx,   evapcdp)
   call pbuf_get_field(pbuf, prec_dp_idx,       prec)
   call pbuf_get_field(pbuf, snow_dp_idx,       snow)
   call pbuf_get_field(pbuf, icimrdp_idx,       qi)
   call pbuf_get_field(pbuf, dlfzm_idx,         dlf)
   call pbuf_get_field(pbuf, difzm_idx,         dif)

   ! DCAPE-ULL
   if (trigdcape_ull .or. trig_dcape_only) then
      call pbuf_get_field(pbuf, t_star_idx, t_star)
      call pbuf_get_field(pbuf, q_star_idx, q_star)
      if ( is_first_step()) then
         q_star(1:ncol,1:pver) = state%q(1:ncol,1:pver,1)
         t_star(1:ncol,1:pver) = state%t(1:ncol,1:pver)
      end if
   end if

   if (zm_microp) then
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlf)
      call pbuf_get_field(pbuf, dnifzm_idx, dnif)
      call pbuf_get_field(pbuf, dsfzm_idx,  dsf)
      call pbuf_get_field(pbuf, dnsfzm_idx, dnsf)
      call pbuf_get_field(pbuf, wuc_idx,    wuc)
   else
      allocate(dnlf(pcols,pver), &
               dnif(pcols,pver), &
               dsf(pcols,pver),  &
               dnsf(pcols,pver), &
               wuc(pcols,pver)   )
   end if
   wuc(1:pcols,1:pver) = 0

   call pbuf_get_field(pbuf, lambdadpcu_idx, lambdadpcu)
   call pbuf_get_field(pbuf, mudpcu_idx,     mudpcu)

   if (zm_microp) then
      if (nmodes > 0) then

         ! Associate pointers with the modes and species that affect the climate (list 0)
         do m = 1, nmodes
            call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, aero(lchnk)%num_a(m)%val)
            call pbuf_get_field(pbuf, dgnum_idx, aero(lchnk)%dgnum(m)%val, start=(/1,1,m/), kount=(/pcols,pver,1/))
            do l = 1, aero(lchnk)%nspec(m)
               call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, aero(lchnk)%mmr_a(l,m)%val)
            end do
         end do

      else if (nbulk > 0) then

         ! Associate pointers with the bulk aerosols that affect the climate (list 0)
         do m = 1, nbulk
            call rad_cnst_get_aer_mmr(0, m, state, pbuf, aero(lchnk)%mmr_bulk(m)%val)
         end do

      end if
   end if

   ! Call the primary Zhang-McFarlane convection parameterization
   call t_startf ('zm_convr')
   call zm_convr( lchnk, ncol, state%t, state%q(:,:,1), prec, jctop, jcbot, &
                  pblh, state%zm, state%phis, state%zi, ptend_loc%q(:,:,1), &
                  ptend_loc%s, state%pmid, state%pint, state%pdel, state%omega, &
                  0.5*ztodt, mcon, cme, cape, tpert, dlf, pflx, zdu, rprd, mu, md, du, eu, ed, &
                  dp, dsubcld, jt, maxg, ideep, lengath, ql, rliq, landfrac, &
                  t_star, q_star, dcape, &  
                  aero(lchnk), qi, dif, dnlf, dnif, dsf, dnsf, sprd, rice, frz, mudpcu, &
                  lambdadpcu, microp_st, wuc )
   call t_stopf ('zm_convr')

   if (zm_microp) then
      dlftot(1:ncol,1:pver) = dlf(1:ncol,1:pver) + dif(1:ncol,1:pver) + dsf(1:ncol,1:pver)
   else
      dlftot(1:ncol,1:pver) = dlf(1:ncol,1:pver)
   end if

   !----------------------------------------------------------------------------
   ! Begin MCSP parameterization here
   !----------------------------------------------------------------------------
   if (doslop) then
      du600              = 0
      dv600              = 0
      ZM_depth           = 0
      MCSP_shear         = 0
      Qs_zmconv (1:ncol) = 0
      Qv_zmconv (1:ncol) = 0
      Pa_int_end(1:ncol) = 0
      Qm (1:ncol,1:pver) = 0
      Qmq(1:ncol,1:pver) = 0
      Qmu(1:ncol,1:pver) = 0
      Qmv(1:ncol,1:pver) = 0

      call vertinterp(ncol, pcols, pver, state%pmid, MCSP_storm_speed_pref, state%u,du600)
      call vertinterp(ncol, pcols, pver, state%pmid, MCSP_storm_speed_pref, state%v,dv600)

      do i = 1,ncol
         if (state%pmid(i,pver).gt.MCSP_storm_speed_pref) then
            du600(i) = du600(i)-state%u(i,pver)
            dv600(i) = dv600(i)-state%v(i,pver)
         else
            du600(i) = -999
            dv600(i) = -999
         end if
      end do

      MCSP_shear = du600
      do i = 1,ncol
         if ( jctop(i).ne.pver ) ZM_depth(i) = state%pint(i,pver+1) - state%pmid(i,jctop(i))
      end do

      ! Define parameters
      do i = 1,ncol
         do k = jctop(i),pver
            Qs_zmconv(i) = Qs_zmconv(i) + ptend_loc%s(i,k) * state%pdel(i,k)
            Qv_zmconv(i) = Qv_zmconv(i) + ptend_loc%q(i,k,1) * state%pdel(i,k)
            Pa_int_end(i) = Pa_int_end(i) + state%pdel(i,k)
         end do
      end do

      do i = 1,ncol
         if (jctop(i) .ne. pver) then
            Qs_zmconv(i) = Qs_zmconv(i) / Pa_int_end(i)
            Qv_zmconv(i) = Qv_zmconv(i) / Pa_int_end(i)
         else
            Qs_zmconv(i) = 0
            Qv_zmconv(i) = 0
         end if
      end do

      do i = 1,ncol
         Qm_int_end(i)  = 0
         Qmq_int_end(i) = 0
         Pa_int_end(i)  = 0
         Q_dis(i)       = 0
         Qq_dis(i)      = 0

         if ( (state%pint(i,pver+1)-state%pmid(i,jctop(i))) >= MCSP_conv_depth_min ) then
            if ( abs(du600(i)).ge.MCSP_shear_min .and. &
                 abs(du600(i)).lt.MCSP_shear_max .and. &
                 Qs_zmconv(i).gt.0 ) then
               do k = jctop(i),pver
                  mcsp_top = state%pint(i,pver+1) - state%pmid(i,k)
                  mcsp_bot = state%pint(i,pver+1) - state%pmid(i,jctop(i))

                  Qm(i,k)  = -1 * Qs_zmconv(i) * alpha2 * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
                  Qmq(i,k) = -1 * Qv_zmconv(i) * alpha2 * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
                  Qmq(i,k) = Qm(i,k)/2500000.0_r8/4.0_r8

                  Qmu(i,k) = alphau * (cos(pi*(mcsp_top/mcsp_bot)))
                  Qmv(i,k) = alphav * (cos(pi*(mcsp_top/mcsp_bot)))

                  dpg = state%pdel(i,k)/gravit

                  Qm_int_end(i)  = Qm_int_end(i) + Qm(i,k) * dpg
                  Qm_int_end(i)  = Qm_int_end(i) + (2.0_r8*Qmu(i,k)*ztodt*state%u(i,k)+ &
                                                   Qmu(i,k)*Qmu(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt
                  Qm_int_end(i)  = Qm_int_end(i) + (2.0_r8*Qmv(i,k)*ztodt*state%v(i,k)+ &
                                                   Qmv(i,k)*Qmv(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt
                  Qmq_int_end(i) = Qmq_int_end(i) + Qmq(i,k) * dpg
                  Pa_int_end(i)  = Pa_int_end(i) + state%pdel(i,k)
               end do
               Q_dis(i)  = Qm_int_end(i)  / Pa_int_end(i)
               Qq_dis(i) = Qmq_int_end(i) / Pa_int_end(i)
            end if
         end if
      end do

      MCSP_DT   = 0
      MCSP_freq = 0

      do i = 1,ncol
         do k = jctop(i),pver
            Qcq_adjust = ptend_loc%q(i,k,1) - Qq_dis(i) * gravit

            ptend_loc%s(i,k) = ptend_loc%s(i,k) - Q_dis(i)*gravit ! energy fixer

            MCSP_DT(i,k) = -Q_dis(i)*gravit+Qm(i,k)
            if (abs(Qm(i,k)).gt.0 .and. abs(Qmu(i,k)).gt.0) MCSP_freq(i) = 1

            if (doslop_heat)     ptend_loc%s(i,k)   = ptend_loc%s(i,k) + Qm(i,k)
            if (doslop_moisture) ptend_loc%q(i,k,1) = Qcq_adjust + Qmq(i,k)
            if (doslop_uwind)    ptend_loc%u(i,k)   = Qmu(i,k)
            if (doslop_vwind)    ptend_loc%v(i,k)   = Qmv(i,k)
         end do
      end do

      MCSP_DT(1:ncol,1:pver) = MCSP_DT(1:ncol,1:pver)/cpair
      call outfld('MCSP_DT    ',MCSP_DT,   pcols, lchnk )
      call outfld('MCSP_freq  ',MCSP_freq, pcols, lchnk )
      if (doslop_uwind) call outfld('MCSP_DU    ',ptend_loc%u*86400.0_r8, pcols, lchnk )
      if (doslop_vwind) call outfld('MCSP_DV    ',ptend_loc%v*86400.0_r8, pcols, lchnk )
      call outfld('ZM_depth   ',ZM_depth,   pcols, lchnk )
      call outfld('MCSP_shear ',MCSP_shear, pcols, lchnk )

   end if
   !----------------------------------------------------------------------------
   ! End of MCSP parameterization calculations
   !----------------------------------------------------------------------------

   call outfld('DCAPE',  dcape, pcols, lchnk )
   call outfld('CAPE_ZM', cape, pcols, lchnk )

   ! Output fractional occurrence of ZM convection
   freqzm(:) = 0
   do i = 1,lengath
      freqzm(ideep(i)) = 1.0_r8
   end do
   call outfld('FREQZM  ', freqzm, pcols, lchnk )

   ! Convert mass flux from reported mb/s to kg/m^2/s
   mcon(1:ncol,1:pver) = mcon(1:ncol,1:pver) * 100.0_r8/gravit

   call outfld('CMFMCDZM', mcon, pcols, lchnk )

   ! Store upward and downward mass fluxes in un-gathered arrays + convert from mb/s to kg/m^2/s
   do i = 1,lengath
      do k = 1,pver
         mu_out(ideep(i),k) = mu(i,k) * 100.0_r8/gravit
         md_out(ideep(i),k) = md(i,k) * 100.0_r8/gravit
      end do
   end do

   if (convproc_do_aer .or. convproc_do_gas) then 
      call outfld('ZMMU', mu_out, pcols, lchnk )
      call outfld('ZMMD', md_out, pcols, lchnk )
   else
      call outfld('ZMMU', mu_out(1,1), pcols, lchnk )
      call outfld('ZMMD', md_out(1,1), pcols, lchnk )
   end if

   ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,1:pver)/cpair
   call outfld('ZMDT    ', ftem, pcols, lchnk )
   call outfld('ZMDQ    ', ptend_loc%q(1,1,1), pcols, lchnk )

   maxgsav(1:ncol) = 0 ! zero if no convection. true mean to be MAXI/FREQZM
   pcont(1:ncol) = state%ps(1:ncol)
   pconb(1:ncol) = state%ps(1:ncol)
   do i = 1,lengath
      if (maxg(i).gt.jt(i)) then
         pcont(ideep(i)) = state%pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
         pconb(ideep(i)) = state%pmid(ideep(i),maxg(i))! gathered array
         maxgsav(ideep(i)) = dble(maxg(i))             ! gathered array for launching level
      end if
   end do
   call outfld('PCONVT  ', pcont, pcols, lchnk )
   call outfld('PCONVB  ', pconb, pcols, lchnk )
   call outfld('MAXI  ', maxgsav, pcols, lchnk )

   call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')

   ! add tendency from this process to tendencies from other processes
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc 
   call physics_update(state1, ptend_loc, ztodt)

   ! initialize ptend for next process
   lq(:) = .FALSE.
   lq(1) = .TRUE.
   call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

   ! Determine the phase of the precipitation produced and add latent heat of fusion
   ! Evaporate some of the precip directly into the environment
   ! Allow this to use the updated copy of state (state1) and the fresh ptend_loc type
   ! heating and specific humidity tendencies produced

   call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
   call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
   call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
   call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
   dp_cldliq(1:ncol,1:pver) = 0
   dp_cldice(1:ncol,1:pver) = 0

   call t_startf ('zm_conv_evap')
   call zm_conv_evap(state1%ncol, state1%lchnk, &
                     state1%t, state1%pmid, state1%pdel, state1%q(1:pcols,1:pver,1), &
                     ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
                     rprd, cld, ztodt, prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, sprd, old_snow)
   call t_stopf ('zm_conv_evap')

   evapcdp(1:ncol,1:pver) = ptend_loc%q(1:ncol,1:pver,1)

   ! Write out variables from zm_conv_evap
   ftem(1:ncol,1:pver) = ptend_loc%s(1:ncol,1:pver)/cpair
   call outfld('EVAPTZM ', ftem, pcols, lchnk )
   ftem(1:ncol,1:pver) = tend_s_snwprd  (1:ncol,1:pver)/cpair
   call outfld('FZSNTZM ', ftem, pcols, lchnk )
   ftem(1:ncol,1:pver) = tend_s_snwevmlt(1:ncol,1:pver)/cpair
   call outfld('EVSNTZM ', ftem, pcols, lchnk )

   call outfld('EVAPQZM ', ptend_loc%q(1,1,1), pcols, lchnk )
   call outfld('ZMFLXPRC', flxprec,            pcols, lchnk )
   call outfld('ZMFLXSNW', flxsnow,            pcols, lchnk )
   call outfld('ZMNTPRPD', ntprprd,            pcols, lchnk )
   call outfld('ZMNTSNPD', ntsnprd,            pcols, lchnk )
   call outfld('ZMEIHEAT', ptend_loc%s,        pcols, lchnk )
   call outfld('PRECCDZM', prec,               pcols, lchnk )
   call outfld('PRECZ   ', prec,               pcols, lchnk )

   if (zm_microp) call zm_microphysics_history_out( lchnk, ncol, microp_st, prec, dlf, dif, dnlf, dnif, frz )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc
   call physics_update(state1, ptend_loc, ztodt)

   !----------------------------------------------------------------------------
   ! convective momentum transport

   call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_transport_momentum', ls=.true., lu=.true., lv=.true.)

   winds(1:ncol,1:pver,1) = state1%u(1:ncol,1:pver)
   winds(1:ncol,1:pver,2) = state1%v(1:ncol,1:pver)

   call t_startf ('zm_transport_momentum')
   call zm_transport_momentum( ncol, winds, 2, &
                               mu, md, du, eu, ed, dp, &
                               jt, maxg, ideep, 1, lengath, &
                               wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )
   call t_stopf ('zm_transport_momentum')

   ptend_loc%u(1:ncol,1:pver) = wind_tends(1:ncol,1:pver,1)
   ptend_loc%v(1:ncol,1:pver) = wind_tends(1:ncol,1:pver,2)
   ptend_loc%s(1:ncol,1:pver) = seten     (1:ncol,1:pver)

   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   ! update physics state type state1 with ptend_loc
   call physics_update(state1, ptend_loc, ztodt)

   ftem(1:ncol,1:pver) = seten(1:ncol,1:pver)/cpair
   call outfld('ZMMTT', ftem             , pcols, lchnk )
   call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk )
   call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk )

   ! Output apparent force from  pressure gradient
   call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk )
   call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk )
   call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk )
   call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk )

   ! Output in-cloud winds
   call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk )
   call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk )
   call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk )
   call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk )

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   !----------------------------------------------------------------------------
   ! convective tracer transport

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_transport_tracer_1', lq=lq)

   ! dpdry is not used in next convtran call since cloud liq/ice mixing ratios are moist
   fake_dpdry(1:ncol,1:pver) = 0

   call t_startf ('zm_transport_tracer_1')
   call zm_transport_tracer( ptend_loc%lq, state1%q, pcnst, &
                             mu, md, du, eu, ed, dp, &
                             jt, maxg, ideep, 1, lengath, &
                             fracis, ptend_loc%q, fake_dpdry, ztodt)  
   call t_stopf ('zm_transport_tracer_1')

   call outfld('ZMDICE ', ptend_loc%q(1,1,ixcldice), pcols, lchnk )
   call outfld('ZMDLIQ ', ptend_loc%q(1,1,ixcldliq), pcols, lchnk )

   ! add tendency from this process to tendency from other processes
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   !----------------------------------------------------------------------------
   ! deallocate local copies

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

   if (zm_microp) then
      call zm_microp_st_dealloc(microp_st)
   else
      deallocate(dnlf, dnif, dsf, dnsf)
   end if

end subroutine zm_conv_tend

!===================================================================================================

subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf, mu, eu, du, md, ed, dp, &
                           dsubcld, jt, maxg, ideep, lengath, species_class)
   !----------------------------------------------------------------------------
   ! Purpose: Secondary interface with ZM for additional convective transport
   !----------------------------------------------------------------------------
   use physics_types,      only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,       only: get_nstep
   use physics_buffer,     only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,       only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages,     only: alloc_err
   use physconst,          only: spec_class_aerosol, spec_class_gas
   !----------------------------------------------------------------------------
   ! Arguments
   type(physics_state),             intent(in ):: state        ! Physics state variables
   type(physics_ptend),             intent(out):: ptend        ! indivdual parameterization tendencies
   type(physics_buffer_desc),          pointer :: pbuf(:)      ! physics buffer
   real(r8),                        intent(in) :: ztodt        ! 2 delta t (model time increment)
   real(r8), dimension(pcols,pver), intent(in) :: mu           ! upward cloud mass flux
   real(r8), dimension(pcols,pver), intent(in) :: eu           ! entrainment in updraft
   real(r8), dimension(pcols,pver), intent(in) :: du           ! detrainment in updraft
   real(r8), dimension(pcols,pver), intent(in) :: md           ! downward cloud mass flux
   real(r8), dimension(pcols,pver), intent(in) :: ed           ! entrainment in downdraft
   real(r8), dimension(pcols,pver), intent(in) :: dp           ! layer thickness
   real(r8), dimension(pcols),      intent(in) :: dsubcld      ! wg layer thickness in mbs (between interfaces)
   integer,  dimension(pcols),      intent(in) :: jt           ! wg layer thickness in mbs between lcl and maxi
   integer,  dimension(pcols),      intent(in) :: maxg         ! wg top  level index of deep cumulus convection
   integer,  dimension(pcols),      intent(in) :: ideep        ! wg gathered values of maxi
   integer,                         intent(in) :: lengath      ! number of gathered columns
   integer,  dimension(:),          intent(in) :: species_class! constituent tracer type
   !----------------------------------------------------------------------------
   ! Local variables
   integer :: i, lchnk, ncol, istat, m
   integer :: nstep
   logical :: lq(pcnst)
   integer :: ifld
   real(r8), dimension(pcols,pver) :: dpdry
   real(r8), pointer, dimension(:,:,:) :: fracis ! fraction of transported species that are insoluble
   !----------------------------------------------------------------------------
   ! Initialize

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   !----------------------------------------------------------------------------
   ! convective tracer transport

   lq(:) = .FALSE.
   lq(:) = .not. cnst_is_convtran1(:)
   call physics_ptend_init(ptend, state%psetcols, 'zm_transport_tracer_2', lq=lq )

   ! Associate pointers with physics buffer fields
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

   ! Transport all constituents except cloud water and ice
   if ( (convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero ) then
      do m = 1, pcnst
         if ( (species_class(m)==spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m)==spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
         end if
      end do
   end if

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran for tracers of dry mixing ratio type
      dpdry(1:ncol,1:pver) = 0
      do i = 1,lengath
         dpdry(i,1:pver) = state%pdeldry(ideep(i),1:pver)/100_r8
      end do
      call t_startf ('zm_transport_tracer_2')
      call zm_transport_tracer( ptend%lq, state%q, pcnst,  &
                                mu, md, du, eu, ed, dp,  &
                                jt, maxg, ideep, 1, lengath, &
                                fracis, ptend%q, dpdry, ztodt)
      call t_stopf ('zm_transport_tracer_2')
   end if

   if ( (convproc_do_aer .or. convproc_do_gas) .and. clim_modal_aero ) then
      do m = 1, pcnst
         if ( (species_class(m)==spec_class_aerosol .and. convproc_do_aer) .or. &
              (species_class(m)==spec_class_gas     .and. convproc_do_gas) ) then
            ptend%lq(m) = .false.
            ptend%q(1:ncol,1:pver,m) = 0
         end if
      end do
   end if

end subroutine zm_conv_tend_2

!===================================================================================================

end module zm_conv_intr
