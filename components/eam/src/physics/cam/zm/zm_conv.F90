module zm_conv
   !----------------------------------------------------------------------------
   ! Purpose: primary methods for the Zhang-McFarlane convection scheme
   !----------------------------------------------------------------------------
   ! Contributors: Guang Zhang, Norman McFarlane, Michael Lazare, Phil Rasch,
   !               Rich Neale, Byron Boville, Xiaoliang Song, Walter Hannah
   !----------------------------------------------------------------------------
   ! Relevant literature:
   !
   ! ZM95 => Zhang, G. J., & N. A. McFarlane (1995): Sensitivity of climate
   !   simulations to the parameterization of cumulus convection in the Canadian
   !   climate centre general circulation model. Atmosphere-Ocean, 33(3),
   !   407–446. https://doi.org/10.1080/07055900.1995.9649539
   !
   ! Z02 => Zhang, G. J. (2002): Convective quasi-equilibrium in midlatitude
   !   continental environment and its effect on convective parameterization,
   !   J. Geophys. Res., 107(D14), doi:10.1029/2001JD001005
   !
   ! AS74 => Arakawa, A., and W. H. Schubert (1974): Interaction of a Cumulus
   !   Cloud Ensemble with the Large-Scale Environment, Part I. J. Atmos. Sci.,
   !   31, 674–701, https://doi.org/10.1175/1520-0469(1974)031<0674:IOACCE>2.0.CO;2
   !----------------------------------------------------------------------------
   ! for equations and details the best reference is the CAM5 tehcnical description:
   !   Neale, R., Gettelman, A., Park, S., Chen, C.-C., Lauritzen, P. H., et al.
   !     (2012). Description of the NCAR Community Atmosphere Model (CAM 5.0).
   !     https://doi.org/10.5065/wgtk-4g06
   ! which can also be found here:
   !   https://opensky.ucar.edu/islandora/object/technotes%3A594?search_api_fulltext=CAM5
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_params, only: r8
   use zm_eamxx_bridge_methods,only: cldfrc_fice
#else
   use shr_kind_mod,           only: r8 => shr_kind_r8
   use cloud_fraction,         only: cldfrc_fice
   use zm_microphysics,        only: zm_mphy, zm_microphysics_adjust
#endif
   use zm_conv_cape,           only: compute_dilute_cape
   use zm_conv_types,          only: zm_const_t, zm_param_t
   use zm_conv_util,           only: qsat_hpa
   use zm_aero_type,           only: zm_aero_t
   use zm_microphysics_state,  only: zm_microp_st, zm_microp_st_alloc, zm_microp_st_dealloc
   use zm_microphysics_state,  only: zm_microp_st_ini, zm_microp_st_zero, zm_microp_st_scatter
   !----------------------------------------------------------------------------
   implicit none
   save
   private
   !----------------------------------------------------------------------------
   ! public methods
   public zm_conv_main_init        ! ZM scheme initialization
   public zm_conv_main             ! ZM scheme calculations
   public zm_conv_evap             ! ZM scheme evaporation of precip
   !----------------------------------------------------------------------------
   ! public variables
   type(zm_const_t), public :: zm_const ! derived type to hold ZM constants
   type(zm_param_t), public :: zm_param ! derived type to hold ZM tunable parameters
   !----------------------------------------------------------------------------
   ! private variables
   real(r8), parameter :: cape_threshold_old = 70._r8       ! threshold value of cape for deep convection (old value before DCAPE)
   real(r8), parameter :: cape_threshold_new = 0._r8        ! threshold value of cape for deep convection
   real(r8), parameter :: dcape_threshold    = 0._r8        ! threshold value of dcape for deep convection
   real(r8), parameter :: interp_diff_min    = 1.E-6_r8     ! minimum threshold for interpolation method - see eq (4.109), (4.118), (4.119)
   real(r8), parameter :: omsm               = 0.99999_r8   ! to prevent problems due to round off error
   real(r8), parameter :: small              = 1.e-20_r8    ! small number to limit blowup when normalizing by mass flux
!===================================================================================================
contains
!===================================================================================================

subroutine zm_conv_main_init(limcnv_in, no_deep_pbl_in)
   !----------------------------------------------------------------------------
   ! Purpose: initialize quantities for ZM convection scheme
   !----------------------------------------------------------------------------
   use zm_conv_types, only: zm_const_set_to_global, zm_param_print
   !----------------------------------------------------------------------------
   ! Arguments
   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   logical, intent(in), optional :: no_deep_pbl_in  ! flag to eliminate deep convection within PBL
   !----------------------------------------------------------------------------
   zm_param%limcnv = limcnv_in

   if ( present(no_deep_pbl_in) )  then
      zm_param%no_deep_pbl = no_deep_pbl_in
   else
      zm_param%no_deep_pbl = .false.
   end if

   ! set zm_const using global values
   call zm_const_set_to_global(zm_const)

   ! print parameter values to the log file
   call zm_param_print(zm_param)

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_main_init

!===================================================================================================

subroutine zm_get_gather_index(pcols, ncol, pver, pverp, is_first_step, cape, dcape, &
                               cape_threshold_loc, gather_index, lengath)
   !----------------------------------------------------------------------------
   ! Purpose: determine length of gathered arrays
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols              ! maximum number of columns
   integer,                         intent(in   ) :: ncol               ! actual number of columns
   integer,                         intent(in   ) :: pver               ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp              ! number of interface vertical levels
   logical,                         intent(in   ) :: is_first_step      ! flag for first step of run
   real(r8), dimension(pcols),      intent(in   ) :: cape               ! conv. avail. potential energy     [J]
   real(r8), dimension(pcols),      intent(in   ) :: dcape              ! CAPE generated by dycor (dCAPE)   [J]
   real(r8),                        intent(  out) :: cape_threshold_loc ! cape threshold
   integer,  dimension(pcols),      intent(  out) :: gather_index       ! flag for active columns
   integer,                         intent(  out) :: lengath            ! number of columns for gathered arrays
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, ii
   !----------------------------------------------------------------------------
   ! set local threshold to be used for zm_closure()
   if ( zm_param%trig_dcape .and. (.not.is_first_step) ) then
      cape_threshold_loc = cape_threshold_new
   else
      cape_threshold_loc = cape_threshold_old
   end if

   !----------------------------------------------------------------------------
   ! determine number of active columns
   lengath = 0
   do i = 1,ncol
      if ( zm_param%trig_dcape .and. (.not.is_first_step) ) then
         if ( cape(i)>cape_threshold_loc .and. dcape(i)>dcape_threshold ) then
            lengath = lengath + 1
            gather_index(lengath) = i
         end if
      else
         if (cape(i) > cape_threshold_loc) then
            lengath = lengath + 1
            gather_index(lengath) = i
         end if
      end if
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_get_gather_index

!===================================================================================================

subroutine zm_conv_main(pcols, ncol, pver, pverp, is_first_step, time_step, &
                        t_mid, q_mid_in, omega, p_mid_in, p_int_in, p_del_in, &
                        geos, z_mid_in, z_int_in, pbl_hgt, &
                        tpert, landfrac, t_star, q_star, &
                        lengath, gather_index, maxg, jctop, jcbot, jt, &
                        prec, heat, qtnd, cape, dcape, &
                        mcon, pflx, zdu, mu, eu, du, md, ed, p_del, dsubcld, &
                        ql, rliq, rprd, dlf, aero, microp_st )
   !----------------------------------------------------------------------------
   ! Purpose: Main driver for Zhang-Mcfarlane convection scheme
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols           ! maximum number of columns
   integer,                         intent(in   ) :: ncol            ! actual number of columns
   integer,                         intent(in   ) :: pver            ! number of mid-point levels
   integer,                         intent(in   ) :: pverp           ! number of interface levels
   logical,                         intent(in   ) :: is_first_step   ! flag for first step of run
   real(r8),                        intent(in   ) :: time_step       ! model time-step                   [s]
   real(r8), dimension(pcols,pver), intent(in   ) :: t_mid           ! temperature                       [K]
   real(r8), dimension(pcols,pver), intent(in   ) :: q_mid_in        ! specific humidity                 [kg/kg]
   real(r8), dimension(pcols,pver), intent(in   ) :: omega           ! vertical pressure velocity        [Pa/s]
   real(r8), dimension(pcols,pver), intent(in   ) :: p_mid_in        ! mid-point pressure                [Pa]
   real(r8), dimension(pcols,pverp),intent(in   ) :: p_int_in        ! interface pressure                [Pa]
   real(r8), dimension(pcols,pver), intent(in   ) :: p_del_in        ! pressure thickness                [Pa]
   real(r8), dimension(pcols),      intent(in   ) :: geos            ! surface geopotential              [m2/s2]
   real(r8), dimension(pcols,pver), intent(in   ) :: z_mid_in        ! mid-point geopotential            [m2/s2]
   real(r8), dimension(pcols,pverp),intent(in   ) :: z_int_in        ! interface geopotential            [m2/s2]
   real(r8), dimension(pcols),      intent(in   ) :: pbl_hgt         ! boundary layer height             [m]
   real(r8), dimension(pcols),      intent(in   ) :: tpert           ! parcel temperature perturbation   [K]
   real(r8), dimension(pcols),      intent(in   ) :: landfrac        ! land fraction
   real(r8),pointer,dimension(:,:), intent(in   ) :: t_star          ! for DCAPE - prev temperature      [K]
   real(r8),pointer,dimension(:,:), intent(in   ) :: q_star          ! for DCAPE - prev sp. humidity     [kg/kg]
   integer,                         intent(  out) :: lengath         ! number of active columns in chunk for gathering
   integer,  dimension(pcols),      intent(  out) :: gather_index    ! flag for active columns
   integer,  dimension(pcols),      intent(  out) :: maxg            ! gathered level indices of max MSE (maxi)
   integer,  dimension(pcols),      intent(  out) :: jctop           ! top-of-deep-convection indices
   integer,  dimension(pcols),      intent(  out) :: jcbot           ! base of cloud indices
   integer,  dimension(pcols),      intent(  out) :: jt              ! gathered top level index of deep cumulus convection
   real(r8), dimension(pcols),      intent(  out) :: prec            ! output precipitation
   real(r8), dimension(pcols,pver), intent(  out) :: heat            ! dry static energy tendency        [W/kg]
   real(r8), dimension(pcols,pver), intent(  out) :: qtnd            ! specific humidity tendency        [kg/kg/s]
   real(r8), dimension(pcols),      intent(  out) :: cape            ! conv. avail. potential energy     [J]
   real(r8), dimension(pcols),      intent(  out) :: dcape           ! CAPE generated by dycor (dCAPE)   [J]
   real(r8), dimension(pcols,pverp),intent(  out) :: mcon            ! convective mass flux              [mb/s]
   real(r8), dimension(pcols,pverp),intent(  out) :: pflx            ! precip flux at each level         [?]
   real(r8), dimension(pcols,pver), intent(  out) :: zdu             ! detraining mass flux              [?]
   real(r8), dimension(pcols,pver), intent(  out) :: mu              ! updraft mass flux                 [?]
   real(r8), dimension(pcols,pver), intent(  out) :: eu              ! updraft entrainment               [?]
   real(r8), dimension(pcols,pver), intent(  out) :: du              ! updraft detrainment               [?]
   real(r8), dimension(pcols,pver), intent(  out) :: md              ! downdraft mass flux               [?]
   real(r8), dimension(pcols,pver), intent(  out) :: ed              ! downdraft entrainment             [?]
   real(r8), dimension(pcols,pver), intent(  out) :: p_del           ! layer thickness                   [mb]
   real(r8), dimension(pcols),      intent(  out) :: dsubcld         ! thickness between lcl and maxi    [mb]
   real(r8), dimension(pcols,pver), intent(  out) :: ql              ! cloud liquid water for chem/wetdep
   real(r8), dimension(pcols),      intent(  out) :: rliq            ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), dimension(pcols,pver), intent(  out) :: rprd            ! rain production rate
   real(r8), dimension(pcols,pver), intent(  out) :: dlf             ! detrained cloud liq mixing ratio
   type(zm_aero_t),                 intent(inout) :: aero            ! aerosol object
   type(zm_microp_st),              intent(inout) :: microp_st       ! convective microphysics state and tendencies
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8), dimension(pcols,pver) :: s_mid        ! scaled dry static energy (t+gz/cp)      [K]
   real(r8), dimension(pcols,pver) :: q_mid        ! local copy of specific humidity         [kg/kg]
   real(r8), dimension(pcols,pver) :: p_mid        ! local copy of mid-point pressure        [mb]
   real(r8), dimension(pcols,pverp):: p_int        ! local copy of interface pressure        [mb]
   real(r8), dimension(pcols,pver) :: z_mid        ! local copy of mid-point altitude        [m]
   real(r8), dimension(pcols,pverp):: z_int        ! local copy of interface altitude        [m]
   real(r8), dimension(pcols)      :: z_srf        ! surface altitude                        [m]
   real(r8), dimension(pcols)      :: mumax        ! max value of mu/dp

   integer,  dimension(pcols)      :: pbl_top      ! pbl top indices
   integer,  dimension(pcols)      :: pbl_top_g    ! gathered pbl top indices

   real(r8), dimension(pcols,pver) :: t_pcl        ! parcel temperature                      [K]
   real(r8), dimension(pcols,pver) :: q_pcl_sat    ! parcel saturation specific humidity     [kg/kg]
   real(r8), dimension(pcols)      :: tl           ! parcel temperature at lcl               [K]
   integer,  dimension(pcols)      :: lcl          ! base level index of deep cumulus convection
   integer,  dimension(pcols)      :: lel          ! index of highest theoretical convective plume
   integer,  dimension(pcols)      :: lon          ! index of onset level for deep convection
   integer,  dimension(pcols)      :: maxi         ! index of level with largest moist static energy

   real(r8), dimension(pcols,pver) :: t_pcl_m1     ! time n-1 parcel temperatures
   real(r8), dimension(pcols,pver) :: q_pcl_sat_m1 ! time n-1 parcel saturation specific humidity
   real(r8), dimension(pcols)      :: tl_m1        ! time n-1 parcel Temperature at LCL
   integer,  dimension(pcols)      :: lcl_m1       ! time n-1 base level index of deep cumulus convection
   integer,  dimension(pcols)      :: lel_m1       ! time n-1 index of highest theoretical convective plume
   integer,  dimension(pcols)      :: lon_m1       ! time n-1 index of onset level for deep convection
   integer,  dimension(pcols)      :: maxi_m1      ! time n-1 index of level with largest moist static energy
   real(r8), dimension(pcols)      :: cape_m1      ! time n-1 CAPE

   logical  :: cape_calc_msemax_klev               ! flag for compute_dilute_cape()
   real(r8) :: cape_threshold_loc                  ! local CAPE threshold

   real(r8), dimension(pcols,pver) :: q_mid_g      ! gathered mid-point env specific humidity
   real(r8), dimension(pcols,pver) :: q_int_g      ! gathered interface env specific humidity
   real(r8), dimension(pcols,pver) :: t_mid_g      ! gathered temperature at interface
   real(r8), dimension(pcols,pver) :: p_mid_g      ! gathered values of p_mid
   real(r8), dimension(pcols,pver) :: z_mid_g      ! gathered values of z_mid
   real(r8), dimension(pcols,pver) :: s_mid_g      ! gathered values of s
   real(r8), dimension(pcols,pver) :: s_int_g      ! gathered upper interface dry static energy
   real(r8), dimension(pcols,pverp):: z_int_g      ! gathered values of z_int
   real(r8), dimension(pcols,pver) :: t_pcl_g      ! gathered values of t_pcl
   real(r8), dimension(pcols,pver) :: q_pcl_sat_g  ! gathered values of q_pcl_sat
   real(r8), dimension(pcols,pver) :: omega_g      ! gathered values of omega
   real(r8), dimension(pcols,pver) :: rprd_g       ! gathered rain production rate
   real(r8), dimension(pcols,pver) :: ql_g         ! gathered cloud liquid water
   real(r8), dimension(pcols)      :: cape_g       ! gathered convective available potential energy
   real(r8), dimension(pcols)      :: tl_g         ! gathered values of tl
   real(r8), dimension(pcols)      :: landfrac_g   ! gathered landfrac
   real(r8), dimension(pcols)      :: tpert_g      ! gathered values of tpert (temperature perturbation from PBL)
   integer,  dimension(pcols)      :: lcl_g        ! gathered values of lcl level index
   integer,  dimension(pcols)      :: lel_g        ! gathered values of equilibrium level index

   real(r8), dimension(pcols,pver) :: dqdt         ! gathered specific humidity tendency
   real(r8), dimension(pcols,pver) :: dsdt         ! gathered dry static energy ("temp") tendency at gathered points
   real(r8), dimension(pcols,pver) :: sd           ! gathered downdraft dry static energy
   real(r8), dimension(pcols,pver) :: qd           ! gathered downdraft specific humidity
   real(r8), dimension(pcols,pver) :: mc           ! gathered net upward (scaled by mb) cloud mass flux
   real(r8), dimension(pcols,pver) :: qu           ! gathered updraft specific humidity
   real(r8), dimension(pcols,pver) :: su           ! gathered updraft dry static energy
   real(r8), dimension(pcols,pver) :: q_mid_sat_g  ! gathered mid-point saturation specific humidity
   real(r8), dimension(pcols,pver) :: dl_g         ! gathered detraining cld h2o tend
   real(r8), dimension(pcols,pverp):: pflx_g       ! gathered precip flux at each level
   real(r8), dimension(pcols,pver) :: cu_g         ! gathered condensation rate
   real(r8), dimension(pcols,pver) :: evp_g        ! gathered evap rate of rain in downdraft

   integer,  dimension(pcols)      :: jlcl         ! updraft lifting cond level
   integer,  dimension(pcols)      :: j0           ! detrainment initiation level index
   integer,  dimension(pcols)      :: jd           ! downdraft initiation level index

   real(r8), dimension(pcols)      :: cld_bass_mass_flux ! cloud base mass flux determined from zm_closure()

   type(zm_microp_st) :: loc_microp_st ! local (gathered) convective microphysics state and tendencies

   integer i, ii, k, kk             ! loop iterators
   integer msg                      ! number of missing moisture levels at the top of model

   real(r8) qdifr
   real(r8) sdifr

   integer l, m
   real(r8), parameter :: dcon  = 25.e-6_r8
   real(r8), parameter :: mucon = 5.3_r8

   integer prev_msemax_klev(pcols)  ! launching level index saved from 1st call for CAPE calculation;  used in 2nd call when DCAPE-ULL active

   !----------------------------------------------------------------------------
   ! Set upper limit of convection to "limcnv-1"
   msg = zm_param%limcnv - 1

   !----------------------------------------------------------------------------
   ! initialize various arrays
   do i = 1,ncol
      do k = 1,pver
         qtnd(i,k)  = 0._r8
         heat(i,k)  = 0._r8
         mcon(i,k)  = 0._r8
         dqdt(i,k)  = 0._r8
         dsdt(i,k)  = 0._r8
         pflx(i,k)  = 0._r8
         pflx_g(i,k)= 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         ql_g(i,k)  = 0._r8
         dlf(i,k)   = 0._r8
         dl_g(i,k)  = 0._r8
      end do
      prec(i)        = 0._r8
      rliq(i)        = 0._r8
      pflx(i,pverp)  = 0
      pflx_g(i,pverp)= 0
      pbl_top(i)     = pver
      dsubcld(i)     = 0._r8
      jctop(i)       = pver
      jcbot(i)       = 1
   end do

   !----------------------------------------------------------------------------
   ! Allocate and/or Initialize microphysics state/tend derived types
   if (zm_param%zm_microp) then
      call zm_microp_st_alloc(loc_microp_st, pcols, pver)
      call zm_microp_st_ini(loc_microp_st, pcols, pver)
      call zm_microp_st_ini(microp_st, pcols, pver)
      loc_microp_st%lambdadpcu  = (mucon + 1._r8)/dcon
      loc_microp_st%mudpcu      = mucon
   end if

   !----------------------------------------------------------------------------
   ! calculate local pressure [mb] and height [m] for both interface and mid-point
   do i = 1,ncol
      z_srf(i)        = geos(i)*zm_const%rgrav
      p_int(i,pver+1) = p_int_in(i,pver+1)*0.01_r8
      z_int(i,pver+1) = z_int_in(i,pver+1) + z_srf(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p_mid(i,k) = p_mid_in(i,k)*0.01_r8
         p_int(i,k) = p_int_in(i,k)*0.01_r8
         z_mid(i,k) = z_mid_in(i,k) + z_srf(i)
         z_int(i,k) = z_int_in(i,k) + z_srf(i)
      end do
   end do

   !----------------------------------------------------------------------------
   ! locate PBL top index
   do k = pver-1, msg+1, -1
      do i = 1,ncol
         if (abs(z_mid(i,k)-z_srf(i)-pbl_hgt(i)) < (z_int(i,k)-z_int(i,k+1))*0.5_r8) then
            pbl_top(i) = k
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! store input sp. humidity for calculation of precip (through change in storage)
   ! define dry static energy normalized by cp

   do k = 1,pver
      do i = 1,ncol
         q_mid(i,k)     = q_mid_in(i,k)
         s_mid(i,k)     = t_mid(i,k) + (zm_const%grav/zm_const%cpair)*z_mid(i,k)
         t_pcl(i,k)     = 0.0_r8
         s_int_g(i,k)   = s_mid(i,k)
         q_int_g(i,k)   = q_mid(i,k)
      end do
   end do

   do i = 1,ncol
      cape_g(i)  = 0._r8
      lcl_g(i)   = 1
      lel_g(i)   = pver
      maxg(i)    = 1
      tl_g(i)    = 400._r8
      dsubcld(i) = 0._r8
   end do

   !----------------------------------------------------------------------------
   ! Evaluate Tparcel, qs(Tparcel), buoyancy, CAPE, lcl, lel, parcel launch level at index maxi()=hmax
   ! - call #1, cape_calc_msemax_klev=.true.   standard calculation using state of current step
   ! - call #2, cape_calc_msemax_klev=.false.  use launch level (msemax_klev) from previous call
   ! DCAPE is the difference in CAPE between the two calls using the same launch level
   cape_calc_msemax_klev = .true.
   call compute_dilute_cape( pcols, ncol, pver, pverp, &
                             zm_param%num_cin, msg, &
                             q_mid, t_mid, z_mid, p_mid, p_int, &
                             pbl_top, tpert, &
                             t_pcl, q_pcl_sat, maxi, tl, &
                             lcl, lel, cape, &
                             zm_const, zm_param, &
                             cape_calc_msemax_klev )

   ! Calculate dcape trigger condition
   if ( .not.is_first_step .and. zm_param%trig_dcape ) then
      cape_calc_msemax_klev = .false.
      prev_msemax_klev(:ncol) = maxi(:ncol)
      call compute_dilute_cape( pcols, ncol, pver, pverp, &
                                zm_param%num_cin, msg, &
                                q_star, t_star, z_mid, p_mid, p_int, &
                                pbl_top, tpert, &
                                t_pcl_m1, q_pcl_sat_m1, maxi_m1, tl_m1, &
                                lcl_m1, lel_m1, cape_m1, &
                                zm_const, zm_param, &
                                cape_calc_msemax_klev, prev_msemax_klev )
      dcape(:ncol) = (cape(:ncol)-cape_m1(:ncol))/time_step
   end if

   !----------------------------------------------------------------------------
   ! determine whether active columns for gathering
   call zm_get_gather_index(pcols, ncol, pver, pverp, is_first_step, cape, dcape, &
                            cape_threshold_loc, gather_index, lengath)

   if (lengath.eq.0) then
      ! Deallocate local microphysics arrays before returning
      if (zm_param%zm_microp) call zm_microp_st_dealloc(loc_microp_st)
      return
   end if

   !----------------------------------------------------------------------------
   ! copy data to gathered arrays
   do i = 1,lengath
      do k = 1,pver
         p_del(i,k)        = 0.01_r8*p_del_in(gather_index(i),k)
         q_mid_g(i,k)      = q_mid(gather_index(i),k)
         t_mid_g(i,k)      = t_mid(gather_index(i),k)
         p_mid_g(i,k)      = p_mid(gather_index(i),k)
         z_mid_g(i,k)      = z_mid(gather_index(i),k)
         s_mid_g(i,k)      = s_mid(gather_index(i),k)
         t_pcl_g(i,k)      = t_pcl(gather_index(i),k)
         z_int_g(i,k)      = z_int(gather_index(i),k)
         q_pcl_sat_g(i,k)  = q_pcl_sat(gather_index(i),k)
         omega_g(i,k)      = omega(gather_index(i),k)
      end do
      z_int_g(i,pverp)     = z_int(gather_index(i),pver+1)
      cape_g(i)            = cape(gather_index(i))
      lcl_g(i)             = lcl(gather_index(i))
      lel_g(i)             = lel(gather_index(i))
      maxg(i)              = maxi(gather_index(i))
      tl_g(i)              = tl(gather_index(i))
      landfrac_g(i)        = landfrac(gather_index(i))
      pbl_top_g(i)         = pbl_top(gather_index(i))
      tpert_g(i)           = tpert(gather_index(i))
   end do

   !----------------------------------------------------------------------------
   ! copy aerosol data to gathered arrays
   if (zm_param%zm_microp) then
      if (aero%scheme == 'modal') then
         do m = 1, aero%nmodes
            do i = 1,lengath
               do k = 1,pver
                  aero%numg_a(i,k,m) = aero%num_a(m)%val(gather_index(i),k)
                  aero%dgnumg(i,k,m) = aero%dgnum(m)%val(gather_index(i),k)
                  do l = 1, aero%nspec(m)
                     aero%mmrg_a(i,k,l,m) = aero%mmr_a(l,m)%val(gather_index(i),k)
                  end do
               end do
            end do
         end do
      else if (aero%scheme == 'bulk') then
         do m = 1, aero%nbulk
            do k = 1,pver
               do i = 1,lengath
                  aero%mmrg_bulk(i,k,m) = aero%mmr_bulk(m)%val(gather_index(i),k)
               end do
            end do
         end do
      end if
   end if

   !----------------------------------------------------------------------------
   ! calculate sub-cloud layer pressure "thickness" for closure and tendency calculations
   do k = msg+1, pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + p_del(i,k)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! define interfacial values for (q,s) used in subsequent routines
   do k = msg+2, pver
      do i = 1,lengath
         sdifr = 0._r8
         qdifr = 0._r8
         if (s_mid_g(i,k) > 0._r8 .or. s_mid_g(i,k-1) > 0._r8) &
            sdifr = abs((s_mid_g(i,k)-s_mid_g(i,k-1))/max(s_mid_g(i,k-1),s_mid_g(i,k)))
         if (q_mid_g(i,k) > 0._r8 .or. q_mid_g(i,k-1) > 0._r8) &
            qdifr = abs((q_mid_g(i,k)-q_mid_g(i,k-1))/max(q_mid_g(i,k-1),q_mid_g(i,k)))
         if (sdifr > interp_diff_min) then
            s_int_g(i,k) = log(s_mid_g(i,k-1)/s_mid_g(i,k))*s_mid_g(i,k-1)*s_mid_g(i,k)/(s_mid_g(i,k-1)-s_mid_g(i,k))
         else
            s_int_g(i,k) = 0.5_r8*(s_mid_g(i,k)+s_mid_g(i,k-1))
         end if
         if (qdifr > interp_diff_min) then
            q_int_g(i,k) = log(q_mid_g(i,k-1)/q_mid_g(i,k))*q_mid_g(i,k-1)*q_mid_g(i,k)/(q_mid_g(i,k-1)-q_mid_g(i,k))
         else
            q_int_g(i,k) = 0.5_r8*(q_mid_g(i,k)+q_mid_g(i,k-1))
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! calculate updraft and downdraft properties
   call zm_cloud_properties(pcols, lengath, pver, pverp, msg, zm_param%limcnv, &
                            p_mid_g, z_mid_g, z_int_g, t_mid_g, s_mid_g, s_int_g, q_mid_g, &
                            landfrac_g, tpert_g, &
                            maxg, lel_g, jt, jlcl, j0, jd, &
                            mu, eu, du, md, ed, mc, &
                            su, qu, ql_g, sd, qd,  &
                            q_mid_sat_g, cu_g, evp_g, pflx_g, rprd_g, &
                            aero, loc_microp_st )

   !---------------------------------------------------------------------------
   ! convert from units of "per length" [1/m] to "per pressure" [1/mb].
   do k = msg+1, pver
      do i = 1,lengath
         du    (i,k) = du    (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         eu    (i,k) = eu    (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         ed    (i,k) = ed    (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         cu_g  (i,k) = cu_g  (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         rprd_g(i,k) = rprd_g(i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         evp_g  (i,k)= evp_g (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         if (zm_param%zm_microp) then
            loc_microp_st%frz (i,k) = loc_microp_st%frz (i,k) * (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
            loc_microp_st%sprd(i,k) = loc_microp_st%sprd(i,k) * (z_int_g(i,k)-z_int_g(i,k+1))/p_del(i,k)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! apply closure assumption to calculate cloud base mass flux
   call zm_closure(pcols, lengath, pver, pverp, msg, cape_threshold_loc, &
                   lcl_g, lel_g, jt, maxg, dsubcld, &
                   z_mid_g, z_int_g, p_mid_g, p_del, t_mid_g, &
                   s_mid_g, q_mid_g, q_mid_sat_g, ql_g, s_int_g, q_int_g, &
                   tl_g, t_pcl_g, q_pcl_sat_g, su, qu, &
                   mc, du, mu, md, qd, sd, cape_g, cld_bass_mass_flux )

   !----------------------------------------------------------------------------
   ! limit cloud base mass flux to theoretical upper bound.
   do i = 1,lengath
      mumax(i) = 0
      do k = msg+2, pver
        mumax(i) = max(mumax(i), mu(i,k)/p_del(i,k))
      end do
      if (mumax(i) > 0._r8) then
         cld_bass_mass_flux(i) = min(cld_bass_mass_flux(i),1/(time_step*mumax(i)))
      else
         cld_bass_mass_flux(i) = 0._r8
      end if
      if (zm_param%clos_dyn_adj) then
         cld_bass_mass_flux(i) = max(cld_bass_mass_flux(i) - omega_g(i,pbl_top_g(i))*0.01_r8, 0._r8)
      end if
   end do

   !----------------------------------------------------------------------------
   ! don't allow convection within PBL (suggestion of Bjorn Stevens, 8-2000)
   if (zm_param%no_deep_pbl) then
      do i = 1,lengath
         if (z_mid_in(gather_index(i),jt(i)) < pbl_hgt(gather_index(i))) then
            cld_bass_mass_flux(i) = 0
         end if
      end do
   end if

   !----------------------------------------------------------------------------
   ! apply cloud base mass flux scaling
   do i = 1,lengath
      ! zero out micro data for inactive columns
      if ( zm_param%zm_microp .and. cld_bass_mass_flux(i).eq.0._r8) call zm_microp_st_zero(loc_microp_st,i,pver)
      ! scale variables
      do k = msg+1,pver
         mu    (i,k)   = mu    (i,k)  *cld_bass_mass_flux(i)
         md    (i,k)   = md    (i,k)  *cld_bass_mass_flux(i)
         mc    (i,k)   = mc    (i,k)  *cld_bass_mass_flux(i)
         du    (i,k)   = du    (i,k)  *cld_bass_mass_flux(i)
         eu    (i,k)   = eu    (i,k)  *cld_bass_mass_flux(i)
         ed    (i,k)   = ed    (i,k)  *cld_bass_mass_flux(i)
         rprd_g(i,k)   = rprd_g(i,k)  *cld_bass_mass_flux(i)
         cu_g  (i,k)   = cu_g  (i,k)  *cld_bass_mass_flux(i)
         evp_g  (i,k)  = evp_g (i,k)  *cld_bass_mass_flux(i)
         pflx_g(i,k+1) = pflx_g(i,k+1)*cld_bass_mass_flux(i)*100._r8/zm_const%grav
         ! scale microphysics variables
         if (zm_param%zm_microp) then
            loc_microp_st%sprd(i,k)  = loc_microp_st%sprd(i,k)*cld_bass_mass_flux(i)
            loc_microp_st%frz (i,k)  = loc_microp_st%frz (i,k)*cld_bass_mass_flux(i)
            if (cld_bass_mass_flux(i).eq.0._r8) then
               ql_g(i,k) = 0._r8
            end if
         end if
      end do ! k
   end do ! i

   !----------------------------------------------------------------------------
   ! compute temperature and moisture changes due to convection.
   call zm_calc_output_tend(pcols, lengath, pver, pverp, msg, &
                            jt, maxg, dsubcld, p_del, s_int_g, q_int_g, su, qu, &
                            mu, du, md, sd, qd, ql_g, evp_g, cu_g, &
                            dsdt, dqdt, dl_g, &
                            loc_microp_st)

   !----------------------------------------------------------------------------
   ! conservation check and adjusment
#ifndef SCREAM_CONFIG_IS_CMAKE
   if (zm_param%zm_microp) then
      call zm_microphysics_adjust(pcols, lengath, pver, jt, msg, time_step, zm_const, &
                                  p_del, q_mid_g, dl_g, dsdt, dqdt, rprd_g, loc_microp_st)
   end if
#endif

   !----------------------------------------------------------------------------
   ! scatter data (i.e. undo the gathering)
   do k = msg+1, pver
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
      do i = 1,lengath
         ! q is updated to compute net precip.
         q_mid(gather_index(i),k) = q_mid_in(gather_index(i),k) + time_step*dqdt(i,k)
         qtnd(gather_index(i),k)  = dqdt   (i,k)
         rprd(gather_index(i),k)  = rprd_g (i,k)
         zdu (gather_index(i),k)  = du     (i,k)
         mcon(gather_index(i),k)  = mc     (i,k)
         heat(gather_index(i),k)  = dsdt   (i,k)*zm_const%cpair
         dlf (gather_index(i),k)  = dl_g   (i,k)
         pflx(gather_index(i),k)  = pflx_g (i,k)
         ql  (gather_index(i),k)  = ql_g   (i,k)
      end do
   end do
   ! scatter 1d variables
   do i = 1,lengath
      jctop(gather_index(i)) = jt(i)
      jcbot(gather_index(i)) = maxg(i)
      pflx(gather_index(i),pverp) = pflx_g(i,pverp)
   end do

   !----------------------------------------------------------------------------
   ! scatter microphysics data (i.e. undo the gathering)
   if (zm_param%zm_microp) then
      call zm_microp_st_scatter(loc_microp_st,microp_st,pcols,lengath,pver,gather_index)
   end if

#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif

   !----------------------------------------------------------------------------
   ! Compute precip by integrating change in water vapor minus detrained cloud water
   do i = 1,ncol
      do k = pver, msg+1, -1
         if (zm_param%zm_microp) then
            prec(i) = prec(i) - p_del_in(i,k)*(q_mid(i,k)-q_mid_in(i,k)) - p_del_in(i,k)*(dlf(i,k)+microp_st%dif(i,k)+microp_st%dsf(i,k))*time_step
         else
            prec(i) = prec(i) - p_del_in(i,k)*(q_mid(i,k)-q_mid_in(i,k)) - p_del_in(i,k)*(dlf(i,k))*time_step
         end if
      end do
      ! obtain final precipitation rate in m/s
      prec(i) = zm_const%rgrav*max(prec(i),0._r8)/ time_step/1000._r8
   end do

   !----------------------------------------------------------------------------
   ! Compute reserved liquid (and ice) that is not yet in cldliq for energy integrals
   ! Treat rliq as flux out bottom, to be added back later
   do k = 1,pver
      do i = 1,ncol
         if (zm_param%zm_microp) then
            rliq(i) = rliq(i) + (dlf(i,k)+microp_st%dif(i,k)+microp_st%dsf(i,k))*p_del_in(i,k)/zm_const%grav
            microp_st%rice(i) = microp_st%rice(i) &
                              + (microp_st%dif(i,k)+microp_st%dsf(i,k))*p_del_in(i,k)/zm_const%grav
         else
            rliq(i) = rliq(i) + dlf(i,k)*p_del_in(i,k)/zm_const%grav
         end if
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8
   if (zm_param%zm_microp) microp_st%rice(:ncol) = microp_st%rice(:ncol) /1000._r8

   !----------------------------------------------------------------------------
   ! Deallocate microphysics arrays
   if (zm_param%zm_microp) call zm_microp_st_dealloc(loc_microp_st)

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_main

!===================================================================================================

subroutine zm_conv_evap(pcols, ncol, pver, pverp, time_step, &
                        p_mid, p_del, t_mid, q_mid, prdprec, cldfrc, &
                        tend_s, tend_q, tend_s_snwprd, tend_s_snwevmlt, &
                        prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, microp_st )
   !----------------------------------------------------------------------------
   ! Purpose: - compute tendencies due to evaporation of rain from ZM scheme,
   !          - compute total precip and snow fluxes at the surface
   !          - add the latent heat of fusion for snow formation and melt
   !          - evaporate some precip directly into the environment using a Sundqvist type algorithm
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_wv_saturation, only: qsat
#else
   use wv_saturation,  only: qsat
#endif
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols              ! maximum number of columns
   integer,                         intent(in   ) :: ncol               ! actual number of columns
   integer,                         intent(in   ) :: pver               ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp              ! number of interface vertical levels
   real(r8),                        intent(in   ) :: time_step          ! model time step                         [s]
   real(r8), dimension(pcols,pver), intent(in   ) :: p_mid              ! midpoint pressure                       [Pa]
   real(r8), dimension(pcols,pver), intent(in   ) :: p_del              ! layer thickness                         [Pa]
   real(r8), dimension(pcols,pver), intent(in   ) :: t_mid              ! temperature                             [K]
   real(r8), dimension(pcols,pver), intent(in   ) :: q_mid              ! water vapor                             [kg/kg]
   real(r8), dimension(pcols,pver), intent(in   ) :: prdprec            ! precipitation production                [kg/kg/s]
   real(r8), dimension(pcols,pver), intent(in   ) :: cldfrc             ! cloud fraction
   real(r8), dimension(pcols,pver), intent(inout) :: tend_s             ! heating rate                            [J/kg/s]
   real(r8), dimension(pcols,pver), intent(inout) :: tend_q             ! water vapor tendency                    [kg/kg/s]
   real(r8), dimension(pcols,pver), intent(out  ) :: tend_s_snwprd      ! Heating rate of snow production         [J/kg/s]
   real(r8), dimension(pcols,pver), intent(out  ) :: tend_s_snwevmlt    ! Heating rate of snow evap/melt          [J/kg/s]
   real(r8), dimension(pcols),      intent(inout) :: prec               ! Convective-scale prec rate              [m/s]
   real(r8), dimension(pcols),      intent(out  ) :: snow               ! Convective-scale snow rate              [m/s]
   real(r8), dimension(pcols,pver), intent(out  ) :: ntprprd            ! net precip production in layer          [?]
   real(r8), dimension(pcols,pver), intent(out  ) :: ntsnprd            ! net snow production in layer            [?]
   real(r8), dimension(pcols,pverp),intent(out  ) :: flxprec            ! Convective flux of prec at interfaces   [kg/m2/s]
   real(r8), dimension(pcols,pverp),intent(out  ) :: flxsnow            ! Convective flux of snow at interfaces   [kg/m2/s]
   type(zm_microp_st),              intent(inout) :: microp_st          ! ZM microphysics data structure
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k ! loop iterators
   real(r8), dimension(pcols,pver) :: es           ! Saturation vapor pressure
   real(r8), dimension(pcols,pver) :: fice         ! ice fraction in precip production
   real(r8), dimension(pcols,pver) :: fsnow_conv   ! snow fraction in precip production
   real(r8), dimension(pcols,pver) :: qs           ! saturation specific humidity
   real(r8), dimension(pcols,pver) :: prdsnow      ! snow production                   [kg/ks/s]
   real(r8), dimension(pcols)      :: evpvint      ! vertical integral of evaporation
   real(r8), dimension(pcols)      :: evpprec      ! evaporation of precipitation      [kg/kg/s]
   real(r8), dimension(pcols)      :: evpsnow      ! evaporation of snowfall           [kg/kg/s]
   real(r8), dimension(pcols)      :: snowmlt      ! snow melt tendency in layer
   real(r8), dimension(pcols)      :: flxsntm      ! flux of snow into layer, after melting
   real(r8) :: work1    ! temporary work variable
   real(r8) :: work2    ! temporary work variable
   real(r8) :: evplimit ! temporary work variable for evaporation limits
   real(r8) :: dum      ! temporary work variable

   logical :: pergro_active ! flag for perturbation growth test (pergro)
   real(r8), parameter :: pergro_perturbation = 8.64e-11_r8 ! value used when pergro_active is true
   !----------------------------------------------------------------------------
   ! set flag for perturbation growth test
#ifdef PERGRO
   pergro_active = .true.
#else
   pergro_active = .false.
#endif
   !----------------------------------------------------------------------------
   if (zm_param%zm_microp) then
      prdsnow(1:ncol,1:pver) = microp_st%sprd(1:ncol,1:pver)
   else
      prdsnow(1:ncol,1:pver) = 0._r8
   end if

   ! convert input precip to kg/m2/s
   prec(:ncol) = prec(:ncol)*1000._r8

   ! determine saturation vapor pressure
   call qsat( t_mid(1:ncol,1:pver), p_mid(1:ncol,1:pver), es(1:ncol,1:pver), qs(1:ncol,1:pver))

   ! determine ice fraction in rain production (use cloud water parameterization fraction at present)
   call cldfrc_fice(ncol, t_mid, fice, fsnow_conv)

   ! zero the flux integrals on the top boundary
   flxprec(:ncol,1) = 0._r8
   flxsnow(:ncol,1) = 0._r8
   evpvint(:ncol)   = 0._r8

   do k = 1,pver
      do i = 1,ncol

         ! Melt snow falling into layer, if necessary.
         if (zm_param%old_snow) then
            if (t_mid(i,k) > zm_const%tfreez) then
               flxsntm(i) = 0._r8
               snowmlt(i) = flxsnow(i,k) * zm_const%grav/ p_del(i,k)
            else
               flxsntm(i) = flxsnow(i,k)
               snowmlt(i) = 0._r8
            end if
         else
            ! make sure melting snow doesn't reduce temperature below threshold
            if (t_mid(i,k) > zm_const%tfreez) then
               dum = -zm_const%latice/zm_const%cpair*flxsnow(i,k)*zm_const%grav/p_del(i,k)*time_step
               if (t_mid(i,k) + dum .le. zm_const%tfreez) then
                  dum = (t_mid(i,k)-zm_const%tfreez)*zm_const%cpair/zm_const%latice/time_step
                  dum = dum/(flxsnow(i,k)*zm_const%grav/p_del(i,k))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if
               dum = dum*omsm
               flxsntm(i) = flxsnow(i,k)*(1.0_r8-dum)
               snowmlt(i) = dum*flxsnow(i,k)*zm_const%grav/ p_del(i,k)
            else
               flxsntm(i) = flxsnow(i,k)
               snowmlt(i) = 0._r8
            end if
         end if

         ! relative humidity depression must be > 0 for evaporation
         evplimit = max(1._r8 - q_mid(i,k)/qs(i,k), 0._r8)

         ! total evaporation depends on flux in the top of the layer
         ! flux prec is the net production above layer minus evaporation into environmet
         evpprec(i) = zm_param%ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))

         ! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
         ! Currently does not include heating/cooling change to qs
         evplimit = max(0._r8, (qs(i,k)-q_mid(i,k)) / time_step)

         ! Don't evaporate more than is falling into the layer from above.
         ! Don't evaporate rain formed in this layer, but if precip production
         ! is negative, remove from the available precip. Negative precip
         ! production occurs because of evaporation in downdrafts.
         evplimit = min(evplimit, flxprec(i,k) * zm_const%grav / p_del(i,k))

         ! Total evaporation cannot exceed input precipitation
         evplimit = min(evplimit, (prec(i) - evpvint(i)) * zm_const%grav / p_del(i,k))

         evpprec(i) = min(evplimit, evpprec(i))

         if (.not.zm_param%old_snow) then
            evpprec(i) = max(0._r8, evpprec(i))
            evpprec(i) = evpprec(i)*omsm
         end if

         ! evaporation of snow depends on snow fraction of total precipitation in the top after melting
         if (flxprec(i,k) > 0._r8) then
            ! use limiters to prevent roundoff problems
            work1 = min( max(0._r8, flxsntm(i)/flxprec(i,k) ), 1._r8)
            if (.not.zm_param%old_snow .and. prdsnow(i,k)>prdprec(i,k)) work1 = 1._r8
            evpsnow(i) = evpprec(i) * work1
         else
            evpsnow(i) = 0._r8
         end if

         ! vertically integrated evaporation
         evpvint(i) = evpvint(i) + evpprec(i) * p_del(i,k)/zm_const%grav

         ! net precip production => production - evaporation
         ntprprd(i,k) = prdprec(i,k) - evpprec(i)

         ! net snow production => precip production * ice fraction - evaporation - melting

         ! the small amount added to flxprec in the work1 expression was increased
         ! from 1e-36 to 8.64e-11 (1e-5 mm/day) to address error growth problems.
         ! This causes temperature partitioning to be used for small flxprec amounts.

         if (zm_param%old_snow) then
            if (pergro_active) then
               work1 = min(max(0._r8,flxsnow(i,k)/(flxprec(i,k)+pergro_perturbation)),1._r8)
            else
               if (flxprec(i,k).gt.0._r8) then
                  work1 = min(max(0._r8,flxsnow(i,k)/flxprec(i,k)),1._r8)
               else
                  work1 = 0._r8
               end if
            end if
            work2 = max(fsnow_conv(i,k), work1)
            if (snowmlt(i).gt.0._r8) work2 = 0._r8
            ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
            tend_s_snwprd  (i,k) = prdprec(i,k)*work2*zm_const%latice
            tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*zm_const%latice
         else
            ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k)*zm_const%grav/p_del(i,k), evpsnow(i)+snowmlt(i))
            tend_s_snwprd  (i,k) = prdsnow(i,k)*zm_const%latice
            tend_s_snwevmlt(i,k) = -min(flxsnow(i,k)*zm_const%grav/p_del(i,k), evpsnow(i)+snowmlt(i) )*zm_const%latice
         end if

         ! precipitation fluxes
         flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * p_del(i,k)/zm_const%grav
         flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * p_del(i,k)/zm_const%grav

         ! protect against rounding error
         flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
         flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)

         ! heating (cooling) and moistening due to evaporation
         ! - latent heat of vaporization for precip production has already been accounted for
         ! - snow is contained in prec
         if (zm_param%old_snow) then
            tend_s(i,k) =-evpprec(i)*zm_const%latvap + ntsnprd(i,k)*zm_const%latice
         else
            tend_s(i,k) =-evpprec(i)*zm_const%latvap + tend_s_snwevmlt(i,k)
         end if

         tend_q(i,k) = evpprec(i)
      end do ! i
   end do ! k

   ! protect against rounding error
   if (.not.zm_param%old_snow) then
      do i = 1,ncol
         if (flxsnow(i,pverp).gt.flxprec(i,pverp)) then
            dum = (flxsnow(i,pverp)-flxprec(i,pverp))*zm_const%grav
            do k = pver, 1, -1
               if (ntsnprd(i,k)>ntprprd(i,k).and. dum > 0._r8) then
                  ntsnprd(i,k) = ntsnprd(i,k) - dum/p_del(i,k)
                  tend_s_snwevmlt(i,k) = tend_s_snwevmlt(i,k) - dum/p_del(i,k)*zm_const%latice
                  tend_s(i,k)  = tend_s(i,k) - dum/p_del(i,k)*zm_const%latice
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

end subroutine zm_conv_evap

!===================================================================================================

subroutine zm_calc_fractional_entrainment(pcols, ncol, pver, pverp, msg, &
                                          jb, jt, j0, z_mid, z_int, dz, &
                                          h_env, h_env_sat, h_env_min, &
                                          lambda, lambda_max)
   !----------------------------------------------------------------------------
   ! Purpose: Determine properties of ZM updrafts and downdrafts
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols        ! maximum number of columns
   integer,                         intent(in   ) :: ncol         ! actual number of columns
   integer,                         intent(in   ) :: pver         ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp        ! number of interface vertical levels
   integer,                         intent(in   ) :: msg          ! number of levels to ignore at model top
   integer,  dimension(pcols),      intent(in   ) :: jb           ! updraft base level
   integer,  dimension(pcols),      intent(in   ) :: jt           ! updraft top level
   integer,  dimension(pcols),      intent(inout) :: j0           ! level where updraft begins detraining
   real(r8), dimension(pcols,pver), intent(in   ) :: z_mid        ! env altitude at mid-point
   real(r8), dimension(pcols,pverp),intent(in   ) :: z_int        ! env altitude at interface
   real(r8), dimension(pcols,pver) ,intent(in   ) :: dz           ! layer thickness
   real(r8), dimension(pcols,pver), intent(in   ) :: h_env        ! env moist stat energy
   real(r8), dimension(pcols,pver), intent(in   ) :: h_env_sat    ! env saturated moist stat energy
   real(r8), dimension(pcols)     , intent(inout) :: h_env_min    ! mid-tropospheric MSE minimum
   real(r8), dimension(pcols,pver), intent(  out) :: lambda       ! fractional entrainment
   real(r8), dimension(pcols)     , intent(  out) :: lambda_max   ! fractional entrainment maximum
   !----------------------------------------------------------------------------
   ! Local variables

   integer :: i,k ! loop iterators

   ! variables used for Taylor series expansion when solving eq (4.78) for lamda_D (i.e. fractional entrainment)
   real(r8), dimension(pcols,pver) :: lambda_tmp   ! fractional entrainment work variable
   real(r8), dimension(pcols,pver) :: k1           ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i2           ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i3           ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i4           ! term for Taylor series
   real(r8), dimension(pcols,pver) :: ihat         ! term for Taylor series
   real(r8), dimension(pcols,pver) :: idag         ! term for Taylor series
   real(r8), dimension(pcols,pver) :: iprm         ! term for Taylor series
   real(r8) :: tmp          ! term for Taylor series
   real(r8) :: expnum       ! term for Taylor series

   real(r8), parameter :: lambda_limit_min = 0._r8     ! limiter
   real(r8), parameter :: lambda_limit_max = 0.0002_r8 ! limiter
   real(r8), parameter :: lambda_threshold = 1.E-6_r8  ! threshold for moving detrainment level
   !----------------------------------------------------------------------------
   ! initialize variables
   do k = 1,pver
      do i = 1,ncol
         k1(i,k) = 0._r8
         i2(i,k) = 0._r8
         i3(i,k) = 0._r8
         i4(i,k) = 0._r8
         lambda_tmp(i,k) = 0._r8
         lambda(i,k)     = 0._r8
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! compute taylor series for approximate lambda(z) below
   do k = pver-1, msg+1, -1
      do i = 1,ncol
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (h_env(i,jb(i))-h_env(i,k))*dz(i,k)
            ihat(i,k) = 0.5_r8* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5_r8* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5_r8* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! re-initialize minimum MSE for ensuing calculation
   h_env_min(1:ncol) = 1.E6_r8
   do k = msg+1, pver
      do i = 1,ncol
         if (k >= j0(i) .and. k <= jb(i) .and. h_env(i,k) <= h_env_min(i)) then
            h_env_min(i) = h_env(i,k)
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! compute approximate lambda(z) using above taylor series - see eq (A6) in ZM95
   do k = msg+2, pver
      do i = 1,ncol
         expnum = 0._r8
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k)  = 0._r8
            expnum = 0._r8
         else
            expnum = h_env(i,jb(i)) - (h_env_sat(i,k-1)*(z_int(i,k)-z_mid(i,k)) + &
                     h_env_sat(i,k)*(z_mid(i,k-1)-z_int(i,k)))/(z_mid(i,k-1)-z_mid(i,k))
         end if
         if ( ( (h_env(i,jb(i))-h_env_min(i)) > 100._r8 .and. expnum > 0._r8 ) .and. k1(i,k) > expnum*dz(i,k)) then
            tmp = expnum / k1(i,k)
            lambda_tmp(i,k) = tmp + &
                              i2(i,k)/k1(i,k) * tmp**2 + &
                              (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2 * tmp**3 + &
                              (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/k1(i,k)**3 * tmp**4
            lambda_tmp(i,k) = max(lambda_tmp(i,k),lambda_limit_min)
            lambda_tmp(i,k) = min(lambda_tmp(i,k),lambda_limit_max)
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! move detrainment level downward if fractional entrainment is too low
   do i = 1,ncol
      if (j0(i) < jb(i)) then
         if ( lambda_tmp(i,j0(i))<lambda_threshold .and. lambda_tmp(i,j0(i)+1)>lambda_tmp(i,j0(i)) ) then
            j0(i) = j0(i) + 1
         end if
      end if
   end do ! i

   !----------------------------------------------------------------------------
   ! ensure that entrainment does not increase above the level that detrainment starts
   do k = msg+2, pver
      do i = 1,ncol
         if (k >= jt(i) .and. k <= j0(i)) then
            lambda_tmp(i,k) = max(lambda_tmp(i,k),lambda_tmp(i,k-1))
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! specify maximum fractional entrainment
   do i = 1,ncol
      lambda_max(i) = lambda_tmp(i,j0(i))
      lambda(i,jb(i)) = lambda_max(i)
   end do

   !----------------------------------------------------------------------------
   ! The modification below comes from:
   !   Rasch, P. J., J. E. Kristjánsson, A comparison of the CCM3 model climate
   !   using diagnosed and predicted condensate parameterizations, J. Clim., 1997.
   do k = pver, msg+1, -1
      do i = 1,ncol
         if ( k >=j0(i) .and. k<=jb(i) ) lambda(i,k) = lambda_tmp(i,j0(i))
         if ( k  <j0(i) .and. k>=jt(i) ) lambda(i,k) = lambda_tmp(i,k)
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   return

end subroutine zm_calc_fractional_entrainment

!===================================================================================================

subroutine zm_downdraft_properties(pcols, ncol, pver, pverp, msg, &
                                   jb, jt, j0, jd, z_int, dz, s_mid, q_mid, h_env, &
                                   lambda, lambda_max, qsthat, hsthat, gamhat, rprd, &
                                   mu, md, ed, sd, qd, hd, qds, evp, totevp )
   !----------------------------------------------------------------------------
   ! Purpose: Calculate properties of ZM downdrafts
   ! Notes:
   ! - Downward mass flux is scaled so that net flux (up-down) at cloud base in not negative
   ! - No downdrafts if jd>=jb
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols        ! maximum number of columns
   integer,                         intent(in   ) :: ncol         ! actual number of columns
   integer,                         intent(in   ) :: pver         ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp        ! number of interface vertical levels
   integer,                         intent(in   ) :: msg          ! number of levels to ignore at model top
   integer,  dimension(pcols),      intent(in   ) :: jb           ! updraft base level
   integer,  dimension(pcols),      intent(inout) :: jt           ! updraft top level
   integer,  dimension(pcols),      intent(in   ) :: j0           ! level where updraft begins detraining
   integer,  dimension(pcols),      intent(inout) :: jd           ! level of downdraft
   real(r8), dimension(pcols,pverp),intent(in   ) :: z_int        ! env altitude at interface
   real(r8), dimension(pcols,pver) ,intent(in   ) :: dz           ! layer thickness
   real(r8), dimension(pcols,pver), intent(in   ) :: s_mid        ! env dry static energy of env [K] (normalized)
   real(r8), dimension(pcols,pver), intent(in   ) :: q_mid        ! env specific humidity
   real(r8), dimension(pcols,pver), intent(in   ) :: h_env        ! ambient env moist stat energy
   real(r8), dimension(pcols,pver), intent(in   ) :: lambda       ! fractional entrainment
   real(r8), dimension(pcols),      intent(in   ) :: lambda_max   ! fractional entrainment max
   real(r8), dimension(pcols,pver), intent(in   ) :: qsthat       ! interface interpolated qst
   real(r8), dimension(pcols,pver), intent(in   ) :: hsthat       ! interface interpolated hst
   real(r8), dimension(pcols,pver), intent(in   ) :: gamhat       ! interface interpolated gamma
   real(r8), dimension(pcols,pver), intent(in   ) :: rprd         ! rate of production of precip at that layer
   real(r8), dimension(pcols,pver), intent(in   ) :: mu           ! updraft mass flux
   real(r8), dimension(pcols,pver), intent(inout) :: md           ! downdraft mass flux
   real(r8), dimension(pcols,pver), intent(inout) :: ed           ! downdraft entrainment rate
   real(r8), dimension(pcols,pver), intent(inout) :: sd           ! dndraft dry static energy [K] (normalized)
   real(r8), dimension(pcols,pver), intent(inout) :: qd           ! dndraft specific humidity [kg/kg]
   real(r8), dimension(pcols,pver), intent(inout) :: hd           ! dndraft moist static energy
   real(r8), dimension(pcols,pver), intent(inout) :: qds          ! dndraft saturation specific humdity
   real(r8), dimension(pcols,pver), intent(inout) :: evp          ! evaporation rate
   real(r8), dimension(pcols),      intent(inout) :: totevp       ! total evap   for dndraft proportionality factor - see eq (4.106)
   !----------------------------------------------------------------------------
   ! Local variables
   integer :: i,k ! loop iterators
   real(r8), dimension(pcols) :: ratmjb ! ?
   real(r8) dz_tmp   ! temporary vertical thickness for downdraft mass flux calculation
   real(r8) mdt      ! temporary downdraft mass flux for normalization
   !----------------------------------------------------------------------------
   ! calculate downdraft mass flux
   do i = 1,ncol
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = h_env(i,jd(i)-1)
      if (jd(i) < jb(i) .and. lambda_max(i) > 0._r8) then
         ! NOTE - this nonsensical lambda_max/lambda_max factor
         ! was retained to preserve BFB results during ZM refactoring
         md(i,jd(i)) = -zm_param%alfa * lambda_max(i) / lambda_max(i)
      end if
   end do ! i
   do k = msg+1, pver
      do i = 1,ncol
         if ((k > jd(i) .and. k <= jb(i)) .and. lambda_max(i) > 0._r8) then
            dz_tmp = z_int(i,jd(i)) - z_int(i,k)
            md(i,k) = -zm_param%alfa / (2._r8*lambda_max(i))*(exp(2._r8*lambda_max(i)*dz_tmp)-1._r8)/dz_tmp
         end if
      end do ! i
   end do ! k
   do k = msg+1, pver
      do i = 1,ncol
         if ((k >= jt(i) .and. k <= jb(i)) .and. lambda_max(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! calculate downdraft entrainment and MSE
   do k = msg+1, pver
      do i = 1,ncol
         if ((k >= jt(i) .and. k <= pver) .and. lambda_max(i) > 0._r8) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*h_env(i,k-1))/mdt
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! calculate downdraft specific humidity
   do k = msg+2, pver
      do i = 1,ncol
         if ((k >= jd(i) .and. k <= jb(i)) .and. lambda_max(i) > 0._r8 .and. jd(i) < jb(i)) then
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k)) / (zm_const%latvap*(1._r8+gamhat(i,k)))
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! downdraft quantities at source level
   do i = 1,ncol
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - zm_const%latvap*qd(i,jd(i)))/zm_const%cpair
   end do

   !----------------------------------------------------------------------------
   ! calculate downdraft evaporation
   do k = msg+2, pver
      do i = 1,ncol
         if (k >= jd(i) .and. k < jb(i) .and. lambda_max(i) > 0._r8) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q_mid(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            mdt = min(md(i,k+1),-small)
            if (zm_param%zm_microp)   evp(i,k) = min(evp(i,k),rprd(i,k))
            sd(i,k+1) = ((zm_const%latvap/zm_const%cpair*evp(i,k)-ed(i,k)*s_mid(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q_mid(i,k)
         end if
      end do ! i
   end do ! k

   do i = 1,ncol
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_downdraft_properties

!===================================================================================================

subroutine zm_cloud_properties(pcols, ncol, pver, pverp, msg, limcnv, &
                               p_mid, z_mid, z_int, t_mid, s_mid, s_int, q_mid, &
                               landfrac, tpert_g, &
                               jb, lel, jt, jlcl, j0, jd, &
                               mu, eu, du, md, ed, mc, &
                               su, qu, ql, sd, qd,  &
                               qst, cu, evp, pflx, rprd, &
                               aero, loc_microp_st )
   !----------------------------------------------------------------------------
   ! Purpose: Determine properties of ZM updrafts and downdrafts
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in ) :: pcols          ! declared number of columns
   integer,                         intent(in ) :: ncol           ! actual number of columns for iteration
   integer,                         intent(in ) :: pver           ! number of mid-point vertical levels
   integer,                         intent(in ) :: pverp          ! number of interface vertical levels
   integer,                         intent(in ) :: msg            ! number of levels to ignore at model top
   integer,                         intent(in ) :: limcnv         ! convection limiting level
   real(r8), dimension(pcols,pver), intent(in ) :: p_mid          ! env pressure at mid-point
   real(r8), dimension(pcols,pver), intent(in ) :: z_mid          ! env altitude at mid-point
   real(r8), dimension(pcols,pverp),intent(in ) :: z_int          ! env altitude at interface
   real(r8), dimension(pcols,pver), intent(in ) :: t_mid          ! env temperature
   real(r8), dimension(pcols,pver), intent(in ) :: s_mid          ! env dry static energy of env [K] (normalized)
   real(r8), dimension(pcols,pver), intent(in ) :: s_int          ! interface values of dry stat energy
   real(r8), dimension(pcols,pver), intent(in ) :: q_mid          ! env specific humidity
   real(r8), dimension(pcols),      intent(in ) :: landfrac       ! Land fraction
   real(r8), dimension(pcols),      intent(in ) :: tpert_g        ! PBL temperature perturbation
   integer,  dimension(pcols),      intent(in ) :: jb             ! updraft base level
   integer,  dimension(pcols),      intent(in ) :: lel            ! updraft parcel launch level
   integer,  dimension(pcols),      intent(out) :: jt             ! updraft plume top
   integer,  dimension(pcols),      intent(out) :: jlcl           ! updraft lifting cond level
   integer,  dimension(pcols),      intent(out) :: j0             ! level where detrainment begins (starting at h_env_min)
   integer,  dimension(pcols),      intent(out) :: jd             ! level of downdraft
   real(r8), dimension(pcols,pver), intent(out) :: mu             ! updraft mass flux
   real(r8), dimension(pcols,pver), intent(out) :: eu             ! entrainment rate of updraft
   real(r8), dimension(pcols,pver), intent(out) :: du             ! detrainement rate of updraft
   real(r8), dimension(pcols,pver), intent(out) :: md             ! downdraft mass flux
   real(r8), dimension(pcols,pver), intent(out) :: ed             ! downdraft entrainment rate
   real(r8), dimension(pcols,pver), intent(out) :: mc             ! net mass flux
   real(r8), dimension(pcols,pver), intent(out) :: su             ! updraft dry static energy [K] (normalized)
   real(r8), dimension(pcols,pver), intent(out) :: qu             ! updraft specific humidity [kg/kg]
   real(r8), dimension(pcols,pver), intent(out) :: ql             ! updraft liq water
   real(r8), dimension(pcols,pver), intent(out) :: sd             ! dndraft dry static energy [K] (normalized)
   real(r8), dimension(pcols,pver), intent(out) :: qd             ! dndraft specific humidity [kg/kg]
   real(r8), dimension(pcols,pver), intent(out) :: qst            ! env saturation mixing ratio
   real(r8), dimension(pcols,pver), intent(out) :: cu             ! condensation rate
   real(r8), dimension(pcols,pver), intent(out) :: evp            ! evaporation rate
   real(r8), dimension(pcols,pverp),intent(out) :: pflx           ! precipitation flux thru layer
   real(r8), dimension(pcols,pver), intent(out) :: rprd           ! rate of production of precip at that layer
   type(zm_aero_t),                 intent(in ) :: aero           ! aerosol object
   type(zm_microp_st)                           :: loc_microp_st  ! state and tendency of convective microphysics
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8), dimension(pcols,pver) :: gamma        ! latent-heating correction for pseudo-adiabatic parcel lifting
   real(r8), dimension(pcols,pver) :: dz           ! layer thickness
   real(r8), dimension(pcols,pver) :: h_env        ! ambient env moist stat energy
   real(r8), dimension(pcols,pver) :: h_env_sat    ! ambient env saturated moist stat energy
   real(r8), dimension(pcols)      :: h_env_min    ! mid-tropospheric MSE minimum
   real(r8), dimension(pcols,pver) :: hu           ! updraft moist static energy
   real(r8), dimension(pcols,pver) :: hd           ! dndraft moist static energy
   real(r8), dimension(pcols,pver) :: qsthat       ! interface interpolated qst
   real(r8), dimension(pcols,pver) :: hsthat       ! interface interpolated hst
   real(r8), dimension(pcols,pver) :: gamhat       ! interface interpolated gamma
   real(r8), dimension(pcols,pver) :: qds          ! dndraft saturation specific humdity
   real(r8), dimension(pcols)      :: c0mask       ! land/ocean modification
   real(r8), dimension(pcols)      :: totpcp       ! total precip for dndraft proportionality factor - see eq (4.106)
   real(r8), dimension(pcols)      :: totevp       ! total evap   for dndraft proportionality factor - see eq (4.106)
   real(r8), dimension(pcols,pver) :: lambda       ! fractional entrainment
   real(r8), dimension(pcols)      :: lambda_max   ! fractional entrainment max

   ! Convective microphysics
   real(r8), dimension(pcols,pver) :: fice         ! ice fraction in precip production
   real(r8), dimension(pcols,pver) :: tug          ! temporary updraft temperature
   real(r8), dimension(pcols,pver) :: tmp_frz      ! temporary rate of freezing
   real(r8), dimension(pcols)      :: tot_frz      ! total column freezing rate
   real(r8), dimension(pcols,pverp):: pflxs        ! frozen precipitation flux thru layer
   integer, dimension(pcols)       :: jto          ! updraft plume old top

   real(r8) :: zuef           ! temporary vertical thickness for updraft calculations
   real(r8) :: rmue           ! temporary for updraft entrainment calculations
   real(r8) :: ql1            ! temporary cloud water
   real(r8) :: tu             ! updraft temperature
   real(r8) :: est            ! saturation vapor pressure
   real(r8) :: estu           ! updraft saturation vapor pressure
   real(r8) :: qstu           ! updraft saturation specific humidity
   real(r8) :: dum, sdum      ! dummy variables for round-off error protection

   integer  :: tmp_k_limit    ! temporary limit on k index in various contexts
   integer  :: iter           ! iteration counter
   integer  :: itnum          ! iteration number
   integer  :: khighest       ! k iteration limit
   integer  :: klowest        ! k iteration limit
   integer  :: kount          ! counter for LCL determination
   integer  :: i,k            ! loop iterator

   logical, dimension(pcols) :: doit ! flag to reset cloud
   logical, dimension(pcols) :: done ! flag for LCL determination

   real(r8), parameter :: mu_min       = 0.02_r8      ! minimum value of mu
   real(r8), parameter :: t_homofrz    = 233.15_r8    ! homogeneous freezing temperature
   real(r8), parameter :: t_mphase     = 40._r8       ! mixed phase temperature = tfreez-t_homofrz = 273.15K - 233.15K
   real(r8), parameter :: hu_diff_min  = -2000._r8    ! limit on hu(i,k)-hsthat(i,k)

   !----------------------------------------------------------------------------
   ! initialize 1D variables
   do i = 1,ncol
      totpcp(i) = 0._r8
      totevp(i) = 0._r8
      c0mask(i) = zm_param%c0_ocn*(1._r8-landfrac(i)) + zm_param%c0_lnd*landfrac(i)
   end do

   ! initialize 2D variables
   do k = 1,pver
      do i = 1,ncol
         ! mass fluxes & mixing variables
         mu(i,k)     = 0._r8
         eu(i,k)     = 0._r8
         du(i,k)     = 0._r8
         mc(i,k)     = 0._r8
         md(i,k)     = 0._r8
         ed(i,k)     = 0._r8
         ! cloud process variables
         ql(i,k)     = 0._r8
         evp(i,k)    = 0._r8
         cu(i,k)     = 0._r8
         rprd(i,k)   = 0._r8
         pflx(i,k)   = 0._r8
         ! calculate layer thickness
         dz(i,k)     = z_int(i,k) - z_int(i,k+1)
         ! calculate saturation specific humidity
         call qsat_hPa(t_mid(i,k), p_mid(i,k), est, qst(i,k))
         if ( p_mid(i,k)-est <= 0._r8 ) qst(i,k) = 1.0_r8
         ! compute gamma - see eq. (4.117)
         gamma(i,k)  = qst(i,k)*(1._r8 + qst(i,k)/zm_const%epsilo) &
                       * zm_const%epsilo*zm_const%latvap/(zm_const%rdair*t_mid(i,k)**2) &
                       * zm_const%latvap/zm_const%cpair
         ! cloud thermodynamic variables
         su(i,k)        = s_mid(i,k)
         sd(i,k)        = s_mid(i,k)
         qd(i,k)        = q_mid(i,k)
         qds(i,k)       = q_mid(i,k)
         qu(i,k)        = q_mid(i,k)
         h_env(i,k)     = zm_const%cpair*t_mid(i,k) + zm_const%grav*z_mid(i,k) + zm_const%latvap*q_mid(i,k)
         h_env_sat(i,k) = zm_const%cpair*t_mid(i,k) + zm_const%grav*z_mid(i,k) + zm_const%latvap*qst(i,k)
         hu(i,k)        = h_env(i,k)
         hd(i,k)        = h_env(i,k)
         ! convective microphysics
         pflxs(i,k)  = 0._r8
         fice(i,k)   = 0._r8
         tug(i,k)    = 0._r8
         tmp_frz(i,k)= 0._r8
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! interpolate the mid-point values to interfaces
   do k = msg+1, pver
      do i = 1,ncol
         hsthat(i,k) = h_env_sat(i,k)
         qsthat(i,k) = qst(i,k)
         gamhat(i,k) = gamma(i,k)
         if (k>msg+1) then
            if (abs(qst(i,k-1)-qst(i,k)) > interp_diff_min) then
               qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
            else
               qsthat(i,k) = qst(i,k)
            end if
            hsthat(i,k) = zm_const%cpair*s_int(i,k) + zm_const%latvap*qsthat(i,k)
            if (abs(gamma(i,k-1)-gamma(i,k)) > interp_diff_min) then
               gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ (gamma(i,k-1)-gamma(i,k))
            else
               gamhat(i,k) = gamma(i,k)
            end if
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! initialize cloud top to highest plume top (upper/lower limit => pver / limcnv+1)
   jt(:) = pver
   jto(:)= pver
   do i = 1,ncol
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      h_env_min(i) = 1.E6_r8
   end do ! i

   !----------------------------------------------------------------------------
   ! find the level of minimum saturated MSE (h_env_sat), where detrainment starts
   do i = 1,ncol
      do k = msg+1, pver
         if (h_env_sat(i,k) <= h_env_min(i) .and. k >= jt(i) .and. k <= jb(i)) then
            h_env_min(i) = h_env_sat(i,k)
            j0(i) = k
         end if
      end do ! k
      ! apply limiters to j0
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
      ! don't let j0 be greater than pver
      j0(i) = min(j0(i),pver)
   end do ! i

   !----------------------------------------------------------------------------
   ! Initialize cloud moist and dry static energies (hu=MSE & su=DSE)
   do k = msg+1, pver
      do i = 1,ncol
         if (k >= jt(i) .and. k <= jb(i)) then
            ! Tunable temperature perturbation (tiedke_add) was already added to parcel hu/su to
            ! represent subgrid temperature perturbation. If PBL temperature pert (tpert) also
            ! represents subgrid temperature pert, tiedke_add may need to be removed. Additionally,
            ! current calculation of PBL temperature perturbation is not accurate enough so that a
            ! new tunable parameter tpert_fac was introduced. This introduced new uncertainties into
            ! the ZM scheme. The original code of ZM scheme will be used when tpert_fix=.true.
            if (zm_param%tpert_fix) then
               hu(i,k) = h_env(i,jb(i)) + zm_const%cpair*zm_param%tiedke_add
               su(i,k) = s_mid(i,jb(i)) +                zm_param%tiedke_add
            else
               hu(i,k) = h_env(i,jb(i)) + zm_const%cpair*(zm_param%tiedke_add+zm_param%tpert_fac*tpert_g(i))
               su(i,k) = s_mid(i,jb(i)) +                 zm_param%tiedke_add+zm_param%tpert_fac*tpert_g(i)
            end if
         end if
      end do ! i
   end do ! k

   !----------------------------------------------------------------------------
   ! calculate fractional entrainment (i.e. "lambda D") - see eq (4.78) from Neale et al. (2012)
   call zm_calc_fractional_entrainment(pcols, ncol, pver, pverp, msg, &
                                       jb, jt, j0, z_mid, z_int, dz, &
                                       h_env, h_env_sat, h_env_min, &
                                       lambda, lambda_max)

   !----------------------------------------------------------------------------
   ! iteration to set cloud properties
   itnum = 1
   if (zm_param%zm_microp) itnum = 2

   do iter = 1,itnum

      do k = pver, msg+1, -1
         do i = 1,ncol
            cu(i,k) = 0._r8
            ql(i,k) = 0._r8
            if (zm_param%zm_microp) then
               loc_microp_st%qliq(i,k) = 0._r8
               loc_microp_st%qice(i,k) = 0._r8
               loc_microp_st%frz(i,k)  = 0._r8
            end if
         end do ! i
      end do ! k
      do i = 1,ncol
         totpcp(i) = 0._r8
         if (zm_param%zm_microp) hu(i,jb(i)) = h_env(i,jb(i)) + zm_const%cpair*zm_param%tiedke_add
      end do

      do k = pver, msg+1, -1
         do i = 1,ncol
            ! intialize updraft mass flux variables - here and below all normalized by cloud base mass flux (mb)
            if (lambda_max(i) > 0._r8) then
               mu(i,jb(i)) = 1._r8
               eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
               if (     zm_param%zm_microp) tmp_k_limit = lel(i)
               if (.not.zm_param%zm_microp) tmp_k_limit = jt(i)
               ! compute profiles of updraft mass fluxes - see eq (4.79) - (4.81)
               if ( k>=tmp_k_limit .and. k<jb(i) ) then
                  zuef = z_int(i,k) - z_int(i,jb(i))
                  rmue = (1._r8/lambda_max(i))* (exp(lambda(i,k+1)*zuef)-1._r8)/zuef
                  mu(i,k) = (1._r8/lambda_max(i))* (exp(lambda(i,k  )*zuef)-1._r8)/zuef
                  eu(i,k) = (rmue-mu(i,k+1))/dz(i,k)
                  du(i,k) = (rmue-mu(i,k))/dz(i,k)
               end if
            end if ! lambda_max(i)>0._r8
         end do ! i
      end do ! k

      khighest = pverp
      klowest = 1
      do i = 1,ncol
         khighest = min(khighest,lel(i))
         klowest = max(klowest,jb(i))
      end do

      do k = klowest-1,khighest,-1
         do i = 1,ncol
            if (k <= jb(i)-1 .and. k >= lel(i) .and. lambda_max(i) > 0._r8) then
               if (mu(i,k) < 0.02_r8) then
                  hu(i,k) = h_env(i,k)
                  mu(i,k) = 0._r8
                  eu(i,k) = 0._r8
                  du(i,k) = mu(i,k+1)/dz(i,k)
               else
                 if (zm_param%zm_microp) then
                   hu(i,k) = ( mu(i,k+1)*hu(i,k+1) &
                              +dz(i,k)*( eu(i,k)*h_env(i,k) &
                                        +zm_const%latice*tmp_frz(i,k) ) &
                             ) / ( mu(i,k) + dz(i,k)*du(i,k) )
                 else
                   hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                            dz(i,k)/mu(i,k)* (eu(i,k)*h_env(i,k)- du(i,k)*h_env_sat(i,k))
                 end if
               end if
            end if
         end do ! i
      end do ! k

      ! reset cloud top index beginning from two layers above the
      ! cloud base (i.e. if cloud is only one layer thick, top is not reset
      do i = 1,ncol
         doit(i) = .true.
         tot_frz(i)= 0._r8
         do k = pver, msg+1, -1
            tot_frz(i)= tot_frz(i) + tmp_frz(i,k)*dz(i,k)
         end do ! k
      end do ! i

      do k = klowest-2, khighest-1, -1
         do i = 1,ncol
            if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
               if (hu(i,k)  <= hsthat(i,k) .and. &
                   hu(i,k+1) > hsthat(i,k+1) .and. &
                   mu(i,k)  >= mu_min) then
                  if ( hu(i,k)-hsthat(i,k) < hu_diff_min) then
                     jt(i) = k + 1
                     doit(i) = .false.
                  else
                     jt(i) = k
                     doit(i) = .false.
                  end if
                else if ( (hu(i,k) > hu(i,jb(i)) .and. tot_frz(i)<=0._r8) .or. mu(i,k) < mu_min) then
                  jt(i) = k + 1
                  doit(i) = .false.
               end if
            end if
         end do ! i
      end do ! k

      do i = 1,ncol
         if (iter == 1) jto(i) = jt(i)
      end do

      do k = pver, msg+1, -1
         do i = 1,ncol
            if (k >= lel(i) .and. k <= jt(i) .and. lambda_max(i) > 0._r8) then
               mu(i,k) = 0._r8
               eu(i,k) = 0._r8
               du(i,k) = 0._r8
               hu(i,k) = h_env(i,k)
            end if
            if (k == jt(i) .and. lambda_max(i) > 0._r8) then
               du(i,k) = mu(i,k+1)/dz(i,k)
               eu(i,k) = 0._r8
               mu(i,k) = 0._r8
            end if
         end do ! i
      end do ! k

      ! determine LCL - see eq (4.127)- (4.130) of Neale et al. (2012)
      done(1:ncol) = .false.
      kount = 0
      do k = pver, msg+2, -1
         do i = 1,ncol
            if (k == jb(i) .and. lambda_max(i) > 0._r8) then
               qu(i,k) = q_mid(i,jb(i))
               su(i,k) = (hu(i,k)-zm_const%latvap*qu(i,k))/zm_const%cpair
            end if
            if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. lambda_max(i) > 0._r8) then
               su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s_mid(i,k)
               qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q_mid(i,k) - du(i,k)*qst(i,k))
               tu = su(i,k) - zm_const%grav/zm_const%cpair*z_int(i,k)
               call qsat_hPa(tu, (p_mid(i,k)+p_mid(i,k-1))/2._r8, estu, qstu)
               if (qu(i,k) >= qstu) then
                  jlcl(i) = k
                  kount = kount + 1
                  done(i) = .true.
               end if
            end if
         end do ! i
         if (kount >= ncol) goto 690
      end do ! k

690 continue

      do k = msg+2, pver
         do i = 1,ncol
            if ((k > jt(i) .and. k <= jlcl(i)) .and. lambda_max(i) > 0._r8) then
               su(i,k) = s_int(i,k)  +             (hu(i,k)-hsthat(i,k)) / (zm_const%cpair* (1._r8+gamhat(i,k)))
               qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k)) / (zm_const%latvap* (1._r8+gamhat(i,k)))
            end if
         end do ! i
      end do ! k

      ! compute condensation in updraft
      do k = pver, msg+2, -1
         do i = 1,ncol
            if (lambda_max(i)>0._r8) then
               if (     zm_param%zm_microp) tmp_k_limit = jlcl(i)+1
               if (.not.zm_param%zm_microp) tmp_k_limit = jb(i)
               if ( k>=jt(i) .and. k<tmp_k_limit ) then
                  if (zm_param%zm_microp) then
                     cu(i,k) = ( ( mu(i,k)*su(i,k) - mu(i,k+1)*su(i,k+1) )/dz(i,k) - eu(i,k)*s_mid(i,k) + du(i,k)*su(i,k) &
                               )/(zm_const%latvap/zm_const%cpair) - zm_const%latice*tmp_frz(i,k)/zm_const%latvap
                  else
                     cu(i,k) = ( ( mu(i,k)*su(i,k) - mu(i,k+1)*su(i,k+1) )/dz(i,k) - ( eu(i,k) - du(i,k) )*s_mid(i,k) &
                               )/(zm_const%latvap/zm_const%cpair)
                  end if
                  ! apply limiters
                  if (k == jt(i)) cu(i,k) = 0._r8
                  cu(i,k) = max(0._r8,cu(i,k))
               end if
            end if ! lambda_max(i)>0._r8
         end do ! i
      end do ! k

      !-------------------------------------------------------------------------
      ! microphysical calculations
      if (zm_param%zm_microp) then

         ! calculate updraft temperature
         tug(1:ncol,:) = t_mid(1:ncol,:)
         do k = pver, msg+2, -1
            do i = 1,ncol
               tug(i,k) = su(i,k) - zm_const%grav/zm_const%cpair*z_int(i,k)
            end do ! i
         end do ! k

         ! specify ice fraction
         fice(1:ncol,:) = 0._r8
         do k = 1,pver-1
            do i = 1,ncol
               if (tug(i,k+1) > zm_const%tfreez) then
                  ! If warmer than zm_const%tfreez then water phase
                  fice(i,k) = 0._r8
               else if (tug(i,k+1) < t_homofrz) then
                  ! If colder than t_homofrz then ice phase
                  fice(i,k) = 1._r8
               else
                  ! mixed phase - ice frac decreasing linearly from t_homofrz to zm_const%tfreez
                  fice(i,k) =(zm_const%tfreez - tug(i,k+1)) / t_mphase
               end if
            end do ! i
         end do ! k

         do k = 1,pver
            do i = 1,ncol
               loc_microp_st%cmei(i,k) = cu(i,k) * fice(i,k)
               loc_microp_st%cmel(i,k) = cu(i,k) * (1._r8-fice(i,k))
            end do ! i
         end do ! k

#ifndef SCREAM_CONFIG_IS_CMAKE
         call  zm_mphy( pcols, ncol, msg, &
                        zm_const%grav, zm_const%cpair, zm_const%rdair, &
                        zm_param%auto_fac, zm_param%accr_fac, zm_param%micro_dcs, &
                        jb, jt, jlcl, su, qu, mu, du, eu, z_int, p_mid, t_mid, q_mid, gamhat, lambda_max, &
                        loc_microp_st%cmel,  loc_microp_st%cmei, aero, &
                        loc_microp_st%qliq,     loc_microp_st%qice,     loc_microp_st%qnl,      loc_microp_st%qni,     &
                        loc_microp_st%qcde,     loc_microp_st%qide,     loc_microp_st%ncde,     loc_microp_st%nide,    &
                        rprd,                   loc_microp_st%sprd,     tmp_frz,                loc_microp_st%wu,      &
                        loc_microp_st%qrain,    loc_microp_st%qsnow,    loc_microp_st%qnr,      loc_microp_st%qns,     &
                        loc_microp_st%qgraupel, loc_microp_st%qng,      loc_microp_st%qsde,     loc_microp_st%nsde,    &
                        loc_microp_st%autolm,   loc_microp_st%accrlm,   &
                        loc_microp_st%bergnm,   loc_microp_st%fhtimm,   loc_microp_st%fhtctm,   loc_microp_st%fhmlm,   &
                        loc_microp_st%hmpim,    loc_microp_st%accslm,   loc_microp_st%dlfm,     loc_microp_st%autoln,  &
                        loc_microp_st%accrln,   loc_microp_st%bergnn,   loc_microp_st%fhtimn,   loc_microp_st%fhtctn,  &
                        loc_microp_st%fhmln,    loc_microp_st%accsln,   loc_microp_st%activn,   loc_microp_st%dlfn,    &
                        loc_microp_st%autoim,   loc_microp_st%accsim,   loc_microp_st%difm,     loc_microp_st%nuclin,  &
                        loc_microp_st%autoin,   loc_microp_st%accsin,   loc_microp_st%hmpin,    loc_microp_st%difn,    &
                        loc_microp_st%trspcm,   loc_microp_st%trspcn,   loc_microp_st%trspim,   loc_microp_st%trspin,  &
                        loc_microp_st%lambdadpcu,loc_microp_st%mudpcu,  &
                        loc_microp_st%accgrm,   loc_microp_st%accglm,   loc_microp_st%accgslm,  loc_microp_st%accgsrm, &
                        loc_microp_st%accgirm,  loc_microp_st%accgrim,  loc_microp_st%accgrsm,  loc_microp_st%accgsln, &
                        loc_microp_st%accgsrn,  loc_microp_st%accgirn,  loc_microp_st%accsrim,  loc_microp_st%acciglm, &
                        loc_microp_st%accigrm,  loc_microp_st%accsirm,  loc_microp_st%accigln,  loc_microp_st%accigrn, &
                        loc_microp_st%accsirn,  loc_microp_st%accgln,   loc_microp_st%accgrn,   loc_microp_st%accilm,  &
                        loc_microp_st%acciln,   loc_microp_st%fallrm,   loc_microp_st%fallsm,   loc_microp_st%fallgm,  &
                        loc_microp_st%fallrn,   loc_microp_st%fallsn,   loc_microp_st%fallgn,   loc_microp_st%fhmrm,   &
                        loc_microp_st%dsfm,     loc_microp_st%dsfn )
#endif

         do k = pver, msg+2, -1
            do i = 1,ncol
               ! In the original ZM scheme, which does not consider ice phase, ql actually represents total cloud
               ! water. With convective microphysics, loc_microp_st%qliq and loc_microp_st%qice represent cloud
               ! liquid water and cloud ice, respectively. Since ql is still used in other subroutines as total
               ! cloud water, here ql is calculated as total cloud water for consistency.
               ql(i,k) = loc_microp_st%qliq(i,k)+ loc_microp_st%qice(i,k)
               loc_microp_st%frz(i,k) = tmp_frz(i,k)
            end do ! i
         end do ! k

         do i = 1,ncol
           if (iter == 2 .and. jt(i)> jto(i)) then
             do k = jt(i), jto(i), -1
                loc_microp_st%frz(i,k) = 0.0_r8
                cu(i,k)=0.0_r8
             end do ! k
           end if
         end do ! i

         do k = pver, msg+2, -1
            do i = 1,ncol
               if (k >= jt(i) .and. k < jb(i) .and. lambda_max(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
                  totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*( loc_microp_st%qcde(i,k+1) &
                                                                    +loc_microp_st%qide(i,k+1) &
                                                                    +loc_microp_st%qsde(i,k+1) ))
               end if
            end do ! i
         end do ! k

      else  ! no microphysics

         ! compute condensed liquid, rain production rate
         ! accumulate total precipitation (condensation - detrainment of liquid)
         ! Note: ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
         ! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is consistently applied.
         !   interface quantities => mu, ql are 
         !   mid-point quantities => cu, du, eu, rprd
         do k = pver, msg+2, -1
            do i = 1,ncol
               rprd(i,k) = 0._r8
               if (k >= jt(i) .and. k < jb(i) .and. lambda_max(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
                  if (mu(i,k) > 0._r8) then
                     ql1 = 1._r8/mu(i,k) * ( mu(i,k+1)*ql(i,k+1) - dz(i,k)*du(i,k)*ql(i,k+1) + dz(i,k)*cu(i,k) )
                     ql(i,k) = ql1 / (1._r8+dz(i,k)*c0mask(i))
                  else
                     ql(i,k) = 0._r8
                  end if
                  totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))
                  rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
               end if
            end do ! i
         end do ! k

      end if  ! zm_param%zm_microp

   end do ! iter = 1,itnum

   !----------------------------------------------------------------------------
   ! calculate downdraft properties
   call zm_downdraft_properties(pcols, ncol, pver, pverp, msg, &
                                jb, jt, j0, jd, z_int, dz, s_mid, q_mid, h_env, &
                                lambda, lambda_max, qsthat, hsthat, gamhat, rprd, &
                                mu, md, ed, sd, qd, hd, qds, evp, totevp)

   !----------------------------------------------------------------------------
   ! ensure totpcp and totevp are non-negative
   do i = 1,ncol
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do

   do k = msg+2, pver
      do i = 1,ncol
         ! also ensure that downdraft strength is consistent with precipitation availability - see eq (4.106)
         if (totevp(i) > 0._r8 .and. totpcp(i) > 0._r8) then
            md(i,k)  = md (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(i,k)  = ed (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(i,k)  = 0._r8
            ed(i,k)  = 0._r8
            evp(i,k) = 0._r8
         end if
         ! rprd is the cloud water converted to rain - (rain evaporated)
         if (zm_param%zm_microp) then
           if (rprd(i,k)> 0._r8)  then
              loc_microp_st%frz(i,k)  = loc_microp_st%frz(i,k) - evp(i,k)*min(1._r8,loc_microp_st%sprd(i,k)/rprd(i,k))
              loc_microp_st%sprd(i,k) = loc_microp_st%sprd(i,k)- evp(i,k)*min(1._r8,loc_microp_st%sprd(i,k)/rprd(i,k))
           end if
         end if
         rprd(i,k) = rprd(i,k)-evp(i,k)
      end do ! i
   end do ! k

   ! compute the net precipitation flux across interfaces
   do k = 2,pverp
      do i = 1,ncol
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
         if (zm_param%zm_microp) then
            pflxs(i,k) = pflxs(i,k-1) + loc_microp_st%sprd(i,k-1)*dz(i,k-1)
         end if
      end do ! i
   end do ! k

   ! calculate net mass flux
   do k = msg+1, pver
      do i = 1,ncol
         mc(i,k) = mu(i,k) + md(i,k)
      end do ! i
   end do ! k

   if (zm_param%zm_microp) then
      do i = 1,ncol
         ! protect against rounding error
         if (pflxs(i,pverp).gt.pflx(i,pverp)) then
            dum = (pflxs(i,pverp)-pflx(i,pverp))/omsm
            do k = pver, msg+2, -1
               if (loc_microp_st%sprd(i,k) > 0._r8 .and. dum > 0._r8) then
                  sdum = min(loc_microp_st%sprd(i,k),dum/dz(i,k))
                  loc_microp_st%sprd(i,k) = loc_microp_st%sprd(i,k)- sdum
                  loc_microp_st%frz(i,k)  = loc_microp_st%frz(i,k) - sdum
                  dum = dum - sdum*dz(i,k)
               end if
            end do ! k
         end if
         ! disable columns if top is at or below LCL if using ZM microphysics
         if ( jt(i)>=jlcl(i) ) then
            do k = msg+1, pver
               mu(i,k)   = 0._r8
               eu(i,k)   = 0._r8
               du(i,k)   = 0._r8
               ql(i,k)   = 0._r8
               cu(i,k)   = 0._r8
               evp(i,k)  = 0._r8
               md(i,k)   = 0._r8
               ed(i,k)   = 0._r8
               mc(i,k)   = 0._r8
               rprd(i,k) = 0._r8
               fice(i,k) = 0._r8
            end do ! k
            call zm_microp_st_zero(loc_microp_st,i,pver)
         end if
      end do ! i
   end if ! zm_microp

   !----------------------------------------------------------------------------
   return

end subroutine zm_cloud_properties

!===================================================================================================

subroutine zm_closure(pcols, ncol, pver, pverp, msg, cape_threshold_in, &
                      lcl, lel, jt, mx, dsubcld, &
                      z_mid, z_int, p_mid, p_del, t_mid, &
                      s_mid, q_mid, qs, ql, s_int, q_int,  &
                      tl, t_pcl, q_pcl_sat, su, qu, &
                      mc, du, mu, md, qd, sd, cape, cld_bass_mass_flux )
   !----------------------------------------------------------------------------
   ! Purpose: calculate closure condition for ZM convection scheme using the
   !          revised quasi-equilibrium hypothesis of Z02, in which a
   !          quasi-equilibrium exists between the convective and large-scale
   !          modifications of the free-tropospheric CAPE, such that the net
   !          contribution is negilgible. This differs notably from AS74, where
   !          they assumed that CAPE changes from free-tropospheric and boundary
   !          layer changes are in balance. The Z02 revised closure is based on
   !          the observation that the total CAPE change is comparable to the
   !          CAPE change due to boundary layer thermodynamic changes.
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in ) :: pcols             ! maximum number of columns
   integer,                         intent(in ) :: ncol              ! actual number of columns
   integer,                         intent(in ) :: pver              ! number of mid-point vertical levels
   integer,                         intent(in ) :: pverp             ! number of interface vertical levels
   integer,                         intent(in ) :: msg               ! number of levels to ignore at model top
   real(r8),                        intent(in ) :: cape_threshold_in ! CAPE threshold for "cloud work function" (i.e. A)
   integer,  dimension(pcols),      intent(in ) :: lcl               ! index of lcl
   integer,  dimension(pcols),      intent(in ) :: lel               ! index of launch leve
   integer,  dimension(pcols),      intent(in ) :: jt                ! top of updraft
   integer,  dimension(pcols),      intent(in ) :: mx                ! base of updraft
   real(r8), dimension(pcols),      intent(in ) :: dsubcld           ! thickness of subcloud layer
   real(r8), dimension(pcols,pver), intent(in ) :: z_mid             ! altitude (m)
   real(r8), dimension(pcols,pverp),intent(in ) :: z_int             ! height of interface levels
   real(r8), dimension(pcols,pver), intent(in ) :: p_mid             ! ambient pressure (mb)
   real(r8), dimension(pcols,pver), intent(in ) :: p_del             ! pressure thickness of layers
   real(r8), dimension(pcols,pver), intent(in ) :: t_mid             ! ambient temperature
   real(r8), dimension(pcols,pver), intent(in ) :: s_mid             ! ambient dry static energy (normalized)
   real(r8), dimension(pcols,pver), intent(in ) :: q_mid             ! ambient specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: qs                ! ambient saturation specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: ql                ! ambient liquid water mixing ratio
   real(r8), dimension(pcols,pver), intent(in ) :: s_int             ! env. normalized dry static energy at intrfcs
   real(r8), dimension(pcols,pver), intent(in ) :: q_int             ! environment specific humidity at interfaces
   real(r8), dimension(pcols),      intent(in ) :: tl                ! parcel temperature at lcl
   real(r8), dimension(pcols,pver), intent(in ) :: t_pcl             ! parcel temperature
   real(r8), dimension(pcols,pver), intent(in ) :: q_pcl_sat         ! parcel specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: su                ! updraft dry static energy (normalized)
   real(r8), dimension(pcols,pver), intent(in ) :: qu                ! updraft specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: mc                ! net convective mass flux
   real(r8), dimension(pcols,pver), intent(in ) :: du                ! detrainment from updraft
   real(r8), dimension(pcols,pver), intent(in ) :: mu                ! updraft mass flux
   real(r8), dimension(pcols,pver), intent(in ) :: md                ! dndraft mass flux
   real(r8), dimension(pcols,pver), intent(in ) :: qd                ! dndraft specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: sd                ! dndraft dry static energy
   real(r8), dimension(pcols),      intent(in ) :: cape              ! convective available potential energy
   real(r8), dimension(pcols),      intent(out) :: cld_bass_mass_flux! cloud base mass flux
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8), dimension(pcols,pver) :: dboydt    ! integrand of cape change
   real(r8), dimension(pcols,pver) :: dtmdt     ! free tropospheric tendencies
   real(r8), dimension(pcols,pver) :: dqmdt     ! free tropospheric tendencies
   real(r8), dimension(pcols)      :: dtbdt     ! sub-cloud layer tendencies
   real(r8), dimension(pcols)      :: dqbdt     ! sub-cloud layer tendencies
   real(r8), dimension(pcols)      :: dtldt     ! sub-cloud layer tendencies
   real(r8), dimension(pcols)      :: dadt      ! CAPE consumption rate per unit cloud base mass flux (i.e. "F")
   real(r8) :: dtpdt          ! ?
   real(r8) :: dqsdtp         ! ?
   real(r8) :: thetavp        ! cloud virtual potential temperature
   real(r8) :: thetavm        ! ambient virtual potential temperature
   real(r8) :: dltaa          ! Analogous to the "cloud work function" described by AS74 (i.e. "A")
   real(r8) :: eb             ! ?
   real(r8) :: debdt          ! ?
   integer  :: i, k           ! loop iterators
   integer  :: kmin, kmax     ! vertical iteration limits for vertical integration

   real(r8), parameter :: beta = 0._r8   ! proportion to use liquid water from layer below
   !----------------------------------------------------------------------------
   ! initialization
   dadt(1:ncol) = 0._r8
   dtmdt(1:ncol,(msg+1):pver) = 0._r8
   dqmdt(1:ncol,(msg+1):pver) = 0._r8

   kmin = minval(lel(1:ncol))
   kmax = maxval( mx(1:ncol)) - 1

   !----------------------------------------------------------------------------
   ! Calculate sub-cloud tendencies of virtual temperature and humidity
   do i = 1,ncol
      cld_bass_mass_flux(i) = 0._r8
      eb = p_mid(i,mx(i))*q_mid(i,mx(i))/ (zm_const%epsilo+q_mid(i,mx(i)))
      dtbdt(i) = (1._r8/dsubcld(i)) &
                 *( mu(i,mx(i))*(s_int(i,mx(i))-su(i,mx(i))) &
                   +md(i,mx(i))*(s_int(i,mx(i))-sd(i,mx(i))) )
      dqbdt(i) = (1._r8/dsubcld(i)) &
                 *( mu(i,mx(i))*(q_int(i,mx(i))-qu(i,mx(i))) &
                   +md(i,mx(i))*(q_int(i,mx(i))-qd(i,mx(i))) )
      debdt = zm_const%epsilo*p_mid(i,mx(i)) / (zm_const%epsilo+q_mid(i,mx(i)))**2 * dqbdt(i)
      dtldt(i) = -2840._r8 * (3.5_r8/t_mid(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t_mid(i,mx(i)))-log(eb)-4.805_r8)**2
   end do

   !----------------------------------------------------------------------------
   ! Calculate dtmdt & dqmdt
   do k = msg+1, pver-1
      do i = 1,ncol
         ! cloud top
         if (k==jt(i)) then
            dtmdt(i,k) = (1._r8/p_del(i,k)) &
                         *(mu(i,k+1)*(su(i,k+1)-s_int(i,k+1)-zm_const%latvap/zm_const%cpair*ql(i,k+1)) &
                         + md(i,k+1)*(sd(i,k+1)-s_int(i,k+1)))
            dqmdt(i,k) = (1._r8/p_del(i,k)) &
                         *(mu(i,k+1)*(qu(i,k+1)-q_int(i,k+1)+ql(i,k+1) ) &
                         + md(i,k+1)*(qd(i,k+1)-q_int(i,k+1)))
         end if
         ! below cloud top
         if ( k>jt(i) .and. k<mx(i) ) then
            dtmdt(i,k) = ( mc(i,k  )*(s_int(i,k)-s_mid(i,k)    ) &
                          +mc(i,k+1)*(s_mid(i,k)-s_int(i,k+1)) ) / p_del(i,k) &
                         - zm_const%latvap/zm_const%cpair * du(i,k)*( beta*ql(i,k) + (1-beta)*ql(i,k+1) )
            dqmdt(i,k) = ( mu(i,k+1)*(qu(i,k+1)-q_int(i,k+1)+zm_const%cpair/zm_const%latvap*(su(i,k+1)-s_mid(i,k))) &
                          -mu(i,k  )*(qu(i,k  )-q_int(i,k  )+zm_const%cpair/zm_const%latvap*(su(i,k  )-s_mid(i,k))) &
                          +md(i,k+1)*(qd(i,k+1)-q_int(i,k+1)+zm_const%cpair/zm_const%latvap*(sd(i,k+1)-s_mid(i,k))) &
                          -md(i,k  )*(qd(i,k  )-q_int(i,k  )+zm_const%cpair/zm_const%latvap*(sd(i,k  )-s_mid(i,k))))/p_del(i,k) &
                         + du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! Calculate dboydt (integrand of cape change)
   do k = msg+1, pver
      do i = 1,ncol
         ! levels between parcel launch and LCL
         if ( k>lcl(i) .and. k<mx(i) ) then
            thetavp = t_pcl(i,k)* (1000._r8/p_mid(i,k))**(zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q_mid(i,mx(i)))
            thetavm = t_mid(i,k)* (1000._r8/p_mid(i,k))**(zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q_mid(i,k))
            dboydt(i,k) = (dtbdt(i)/t_mid(i,mx(i))+0.608_r8/ (1._r8+0.608_r8*q_mid(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t_mid(i,k)-0.608_r8/ (1._r8+0.608_r8*q_mid(i,k))*dqmdt(i,k))* &
                          zm_const%grav*thetavp/thetavm
         end if
         ! levels between LCL and cloud top
         if ( k>=lel(i) .and. k<=lcl(i) ) then
            thetavp = t_pcl(i,k)* (1000._r8/p_mid(i,k))**(zm_const%rdair/zm_const%cpair)*(1._r8+1.608_r8*q_pcl_sat(i,k)-q_mid(i,mx(i)))
            thetavm = t_mid(i,k)* (1000._r8/p_mid(i,k))**(zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q_mid(i,k))
            dqsdtp  = q_pcl_sat(i,k) * (1._r8+q_pcl_sat(i,k)/zm_const%epsilo)*zm_const%epsilo*zm_const%latvap/(zm_const%rdair*t_pcl(i,k)**2)
            dtpdt   = t_pcl(i,k)/ (1._r8+zm_const%latvap/zm_const%cpair* (dqsdtp-q_pcl_sat(i,k)/t_pcl(i,k)))* &
                          (dtbdt(i)/t_mid(i,mx(i))+zm_const%latvap/zm_const%cpair* (dqbdt(i)/tl(i)-q_mid(i,mx(i))/tl(i)**2*dtldt(i)))
            dboydt(i,k) = ((dtpdt/t_pcl(i,k)+1._r8/(1._r8+1.608_r8*q_pcl_sat(i,k)-q_mid(i,mx(i)))* &
                          (1.608_r8 * dqsdtp * dtpdt -dqbdt(i))) - (dtmdt(i,k)/t_mid(i,k)+0.608_r8/ &
                          (1._r8+0.608_r8*q_mid(i,k))*dqmdt(i,k)))*zm_const%grav*thetavp/thetavm
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! vertically integrate buoyancy change
   do k = kmin, kmax
      do i = 1,ncol
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)*(z_int(i,k)-z_int(i,k+1))
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! Calculate cloud base mass flux - see eq (8) in Z02
   do i = 1,ncol
      dltaa = -1._r8*( cape(i) - cape_threshold_in )
      if (dadt(i) /= 0._r8) then
         cld_bass_mass_flux(i) = max( dltaa/zm_param%tau/dadt(i), 0._r8)
      end if
      if (zm_param%zm_microp .and. mx(i)-jt(i) < 2._r8) then
         cld_bass_mass_flux(i) = 0.0_r8
      end if
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_closure

!===================================================================================================

subroutine zm_calc_output_tend(pcols, ncol, pver, pverp, msg, &
                               jt, mx, dsubcld, p_del, s_int, q_int, su, qu, &
                               mu, du, md, sd, qd, ql, evp, cu, &
                               dsdt, dqdt, dl, &
                               loc_microp_st)
   !----------------------------------------------------------------------------
   ! Purpose: calculate final output tendencies for the ZM convection scheme
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols        ! maximum number of columns
   integer,                         intent(in   ) :: ncol         ! actual number of columns
   integer,                         intent(in   ) :: pver         ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp        ! number of interface vertical levels
   integer,                         intent(in   ) :: msg          ! number of levels to ignore at model top
   integer,  dimension(pcols),      intent(in   ) :: jt           ! level index of updraft top
   integer,  dimension(pcols),      intent(in   ) :: mx           ! level index of updraft base
   real(r8), dimension(pcols),      intent(in   ) :: dsubcld      ! sub-cloud layer thickness
   real(r8), dimension(pcols,pver), intent(in   ) :: p_del        ! pressure thickness
   real(r8), dimension(pcols,pver), intent(in   ) :: s_int        ! ambient interface dry static energy
   real(r8), dimension(pcols,pver), intent(in   ) :: q_int        ! ambient interface specific humidity
   real(r8), dimension(pcols,pver), intent(in   ) :: su           ! updraft dry static energy
   real(r8), dimension(pcols,pver), intent(in   ) :: qu           ! updraft specific humidity
   real(r8), dimension(pcols,pver), intent(in   ) :: mu           ! updraft mass flux
   real(r8), dimension(pcols,pver), intent(in   ) :: du           ! updraft detrainment
   real(r8), dimension(pcols,pver), intent(in   ) :: md           ! downdraft mass flux
   real(r8), dimension(pcols,pver), intent(in   ) :: sd           ! downdraft dry static energy
   real(r8), dimension(pcols,pver), intent(in   ) :: qd           ! downdraft specific humidity
   real(r8), dimension(pcols,pver), intent(in   ) :: ql           ! cloud liquid water
   real(r8), dimension(pcols,pver), intent(in   ) :: evp          ! evaporation
   real(r8), dimension(pcols,pver), intent(in   ) :: cu           ! updraft condensation
   real(r8), dimension(pcols,pver), intent(  out) :: dqdt         ! output tendency for specific humidity
   real(r8), dimension(pcols,pver), intent(  out) :: dsdt         ! output tendency for dry static energy
   real(r8), dimension(pcols,pver), intent(  out) :: dl           ! output tendency for cloud liquid water
   type(zm_microp_st),              intent(inout) :: loc_microp_st! convective microphysics state and tendencies
   !----------------------------------------------------------------------------
   ! Local variables
   integer i,k
   integer kbm
   integer ktm
   real(r8) emc
   !----------------------------------------------------------------------------
   ! initialize variables
   do k = msg+1, pver
      do i = 1,ncol
         dsdt(i,k) = 0._r8
         dqdt(i,k) = 0._r8
         dl(i,k) = 0._r8
         ! Convective microphysics
         if (zm_param%zm_microp) then
            loc_microp_st%dif(i,k)   = 0._r8
            loc_microp_st%dsf(i,k)   = 0._r8
            loc_microp_st%dnlf(i,k)  = 0._r8
            loc_microp_st%dnif(i,k)  = 0._r8
            loc_microp_st%dnsf(i,k)  = 0._r8
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = 1,ncol
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   !----------------------------------------------------------------------------
   ! calculate large-scale tendencies
   do k = ktm,pver-1
      do i = 1,ncol
         emc = -cu(i,k) + evp(i,k) ! condensation in updraft and evaporating rain in downdraft

         dsdt(i,k) = -zm_const%latvap/zm_const%cpair*emc + &
                     (+mu(i,k+1)*(su(i,k+1)-s_int(i,k+1)) - mu(i,k)*(su(i,k)-s_int(i,k)) &
                      +md(i,k+1)*(sd(i,k+1)-s_int(i,k+1)) - md(i,k)*(sd(i,k)-s_int(i,k)) &
                     )/p_del(i,k)

         dqdt(i,k) = emc + &
                     (+mu(i,k+1)*(qu(i,k+1)-q_int(i,k+1)) - mu(i,k)*(qu(i,k)-q_int(i,k)) &
                      +md(i,k+1)*(qd(i,k+1)-q_int(i,k+1)) - md(i,k)*(qd(i,k)-q_int(i,k)) &
                     )/p_del(i,k)

         if (zm_param%zm_microp) then
            dsdt(i,k) = dsdt(i,k) + zm_const%latice/zm_const%cpair*loc_microp_st%frz(i,k)
            loc_microp_st%dif (i,k) = du(i,k)*loc_microp_st%qide(i,k+1)
            loc_microp_st%dnlf(i,k) = du(i,k)*loc_microp_st%ncde(i,k+1)
            loc_microp_st%dnif(i,k) = du(i,k)*loc_microp_st%nide(i,k+1)
            loc_microp_st%dsf (i,k) = du(i,k)*loc_microp_st%qsde(i,k+1)
            loc_microp_st%dnsf(i,k) = du(i,k)*loc_microp_st%nsde(i,k+1)
            dl(i,k) = du(i,k)*loc_microp_st%qcde(i,k+1)
         else
            dl(i,k) = du(i,k)*ql(i,k+1)
         end if

      end do
   end do

   !----------------------------------------------------------------------------
   ! calculate large-scale tendencies at and below cloud base
#ifdef CPRCRAY
!DIR$ NOINTERCHANGE!
#endif
   do k = kbm,pver
      do i = 1,ncol
         if (k == mx(i)) then
            dsdt(i,k) = (1._r8/dsubcld(i))* (-mu(i,k)*(su(i,k)-s_int(i,k)) &
                                             -md(i,k)*(sd(i,k)-s_int(i,k)) )
            dqdt(i,k) = (1._r8/dsubcld(i))* (-mu(i,k)*(qu(i,k)-q_int(i,k)) &
                                             -md(i,k)*(qd(i,k)-q_int(i,k)) )
         else if (k > mx(i)) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_calc_output_tend

!===================================================================================================

end module zm_conv
