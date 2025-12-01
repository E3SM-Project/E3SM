module zm_conv
   !----------------------------------------------------------------------------
   ! Purpose: primary methods for the Zhang-McFarlane convection scheme
   !----------------------------------------------------------------------------
   ! Contributors: Guang Zhang, Norman McFarlane, Michael Lazare, Phil Rasch,
   !               Rich Neale, Byron Boville, Xiaoliang Song, Walter Hannah
   !----------------------------------------------------------------------------
   ! for equations and technical details the best reference is:
   !   Neale, R., Gettelman, A., Park, S., Chen, C.-C., Lauritzen, P. H., et al.
   !   (2012). Description of the NCAR Community Atmosphere Model (CAM 5.0).
   !    https://doi.org/10.5065/wgtk-4g06
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
   public zm_convi                 ! ZM scheme initialization
   public zm_convr                 ! ZM scheme calculations
   public zm_conv_evap             ! ZM scheme evaporation of precip
   !----------------------------------------------------------------------------
   ! public variables
   type(zm_const_t), public :: zm_const ! derived type to hold ZM constants
   type(zm_param_t), public :: zm_param ! derived type to hold ZM tunable parameters
   !----------------------------------------------------------------------------
   ! private variables
   real(r8), parameter :: capelmt         = 70._r8       ! threshold value for cape for deep convection
   real(r8), parameter :: trigdcapelmt    = 0._r8        ! threshold value of dcape for deep convection
   real(r8), parameter :: omsm            = 0.99999_r8   ! to prevent problems due to round off error
   real(r8), parameter :: interp_diff_min = 1.E-6_r8     ! minimum threshold for interpolation method - see eq (4.109), (4.118), (4.119)
!===================================================================================================
contains
!===================================================================================================

subroutine zm_convi(limcnv_in, no_deep_pbl_in)
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
   endif

   ! set zm_const using global values
   call zm_const_set_to_global(zm_const)

   ! print parameter values to the log file
   call zm_param_print(zm_param)

   !----------------------------------------------------------------------------
   return

end subroutine zm_convi

!===================================================================================================

subroutine zm_convr( pcols, ncol, pver, pverp, is_first_step, delt, &
                     t_mid, qh, omega, p_mid_in, p_int_in, p_del_in, &
                     geos, z_mid_in, z_int_in, pbl_hgt, &
                     tpert, landfrac, t_star, q_star, &
                     lengath, ideep, maxg, jctop, jcbot, jt, &
                     prec, heat, qtnd, cape, dcape, &
                     mcon, pflx, zdu, mu, eu, du, md, ed, dp, dsubcld, &
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
   real(r8),                        intent(in   ) :: delt            ! model time-step                   [s]
   real(r8), dimension(pcols,pver), intent(in   ) :: t_mid           ! temperature                       [K]
   real(r8), dimension(pcols,pver), intent(in   ) :: qh              ! specific humidity                 [kg/kg]
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
   integer,  dimension(pcols),      intent(  out) :: ideep           ! flag for active columns
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
   real(r8), dimension(pcols,pver), intent(  out) :: dp              ! layer thickness                   [mb]
   real(r8), dimension(pcols),      intent(  out) :: dsubcld         ! thickness between lcl and maxi    [mb]
   real(r8), dimension(pcols,pver), intent(  out) :: ql              ! cloud liquid water for chem/wetdep
   real(r8), dimension(pcols),      intent(  out) :: rliq            ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), dimension(pcols,pver), intent(  out) :: rprd            ! rain production rate
   real(r8), dimension(pcols,pver), intent(  out) :: dlf             ! detrained cloud liq mixing ratio
   type(zm_aero_t),                 intent(inout) :: aero            ! aerosol object
   type(zm_microp_st),              intent(inout) :: microp_st       ! convective microphysics state and tendencies
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8), dimension(pcols,pver) :: s            ! scaled dry static energy (t+gz/cp)      [K]
   real(r8), dimension(pcols,pver) :: q            ! local copy of specific humidity         [kg/kg]
   real(r8), dimension(pcols,pver) :: p_mid        ! local copy of mid-point pressure        [mb]
   real(r8), dimension(pcols,pverp):: p_int        ! local copy of interface pressure        [mb]
   real(r8), dimension(pcols,pver) :: z_mid        ! local copy of mid-point altitude        [m]
   real(r8), dimension(pcols,pverp):: z_int        ! local copy of interface altitude        [m]
   real(r8), dimension(pcols)      :: z_srf        ! surface altitude                        [m]
   real(r8), dimension(pcols,pver) :: dlg          ! gathered detraining cld h2o tend
   real(r8), dimension(pcols,pverp):: pflxg        ! gathered precip flux at each level
   real(r8), dimension(pcols,pver) :: cug          ! gathered condensation rate
   real(r8), dimension(pcols,pver) :: evpg         ! gathered evap rate of rain in downdraft
   real(r8), dimension(pcols)      :: mumax        ! max value of mu/dp

   integer,  dimension(pcols)      :: pbl_top      ! pbl top indices
   integer,  dimension(pcols)      :: pbl_top_g    ! gathered pbl top indices

   real(r8), dimension(pcols,pver) :: tp           ! parcel temperature                      [K]
   real(r8), dimension(pcols,pver) :: qstp         ! parcel saturation specific humidity     [kg/kg]
   real(r8), dimension(pcols)      :: tl           ! parcel temperature at lcl               [K]
   integer,  dimension(pcols)      :: lcl          ! base level index of deep cumulus convection
   integer,  dimension(pcols)      :: lel          ! index of highest theoretical convective plume
   integer,  dimension(pcols)      :: lon          ! index of onset level for deep convection
   integer,  dimension(pcols)      :: maxi         ! index of level with largest moist static energy

   real(r8), dimension(pcols,pver) :: tpm1         ! time n-1 parcel temperatures
   real(r8), dimension(pcols,pver) :: qstpm1       ! time n-1 parcel saturation specific humidity
   real(r8), dimension(pcols)      :: tlm1         ! time n-1 parcel Temperature at LCL
   integer,  dimension(pcols)      :: lclm1        ! time n-1 base level index of deep cumulus convection
   integer,  dimension(pcols)      :: lelm1        ! time n-1 index of highest theoretical convective plume
   integer,  dimension(pcols)      :: lonm1        ! time n-1 index of onset level for deep convection
   integer,  dimension(pcols)      :: maxim1       ! time n-1 index of level with largest moist static energy
   real(r8), dimension(pcols)      :: capem1       ! time n-1 CAPE

   logical  iclosure                               ! flag for compute_dilute_cape()
   real(r8) capelmt_wk                             ! local capelmt to allow exceptions when calling closure() with trigdcape

   integer,  dimension(pcols)      :: gather_index ! temporary variable used to set ideep

   real(r8), dimension(pcols,pver) :: qg           ! gathered specific humidity
   real(r8), dimension(pcols,pver) :: tg           ! gathered temperature at interface
   real(r8), dimension(pcols,pver) :: pg           ! gathered values of p
   real(r8), dimension(pcols,pver) :: z_mid_g      ! gathered values of z_mid
   real(r8), dimension(pcols,pver) :: sg           ! gathered values of s
   real(r8), dimension(pcols,pver) :: tpg          ! gathered values of tp
   real(r8), dimension(pcols,pverp):: z_int_g      ! gathered values of z_int
   real(r8), dimension(pcols,pver) :: qstpg        ! gathered values of qstp
   real(r8), dimension(pcols,pver) :: ug           ! gathered values of u
   real(r8), dimension(pcols,pver) :: vg           ! gathered values of v
   real(r8), dimension(pcols,pver) :: omegag       ! gathered values of omega
   real(r8), dimension(pcols,pver) :: rprdg        ! gathered rain production rate
   real(r8), dimension(pcols)      :: capeg        ! gathered convective available potential energy
   real(r8), dimension(pcols)      :: tlg          ! gathered values of tl
   real(r8), dimension(pcols)      :: landfracg    ! gathered landfrac
   real(r8), dimension(pcols)      :: tpertg       ! gathered values of tpert (temperature perturbation from PBL)
   integer,  dimension(pcols)      :: lclg         ! gathered values of lcl level index
   integer,  dimension(pcols)      :: lelg         ! gathered values of equilibrium level index

   ! work fields arising from gathered calculations
   real(r8), dimension(pcols,pver) :: dqdt         ! gathered specific humidity tendency
   real(r8), dimension(pcols,pver) :: dsdt         ! gathered dry static energy ("temp") tendency at gathered points
   real(r8), dimension(pcols,pver) :: sd           ! gathered downdraft dry static energy
   real(r8), dimension(pcols,pver) :: qd           ! gathered downdraft specific humidity
   real(r8), dimension(pcols,pver) :: mc           ! gathered net upward (scaled by mb) cloud mass flux
   real(r8), dimension(pcols,pver) :: qhat         ! gathered upper interface specific humidity
   real(r8), dimension(pcols,pver) :: qu           ! gathered updraft specific humidity
   real(r8), dimension(pcols,pver) :: su           ! gathered updraft dry static energy
   real(r8), dimension(pcols,pver) :: qs           ! gathered saturation specific humidity
   real(r8), dimension(pcols,pver) :: shat         ! gathered upper interface dry static energy
   real(r8), dimension(pcols,pver) :: qlg          ! gathered cloud liquid water
   real(r8), dimension(pcols,pver) :: dudt         ! gathered u-wind tendency at gathered points
   real(r8), dimension(pcols,pver) :: dvdt         ! gathered v-wind tendency at gathered points

   real(r8), dimension(pcols)      :: mb           ! cloud base mass flux
   integer,  dimension(pcols)      :: jlcl         ! updraft lifting cond level
   integer,  dimension(pcols)      :: j0           ! detrainment initiation level index
   integer,  dimension(pcols)      :: jd           ! downdraft initiation level index

   type(zm_microp_st) :: loc_microp_st ! local (gathered) convective microphysics state and tendencies

   integer i, ii, k, kk             ! loop iterators
   integer msg                      ! number of missing moisture levels at the top of model

   real(r8) qdifr
   real(r8) sdifr

   integer l, m
   real(r8), parameter :: dcon  = 25.e-6_r8
   real(r8), parameter :: mucon = 5.3_r8

   integer dcapemx(pcols)  ! launching level index saved from 1st call for CAPE calculation;  used in 2nd call when DCAPE-ULL active

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
         dudt(i,k)  = 0._r8
         dvdt(i,k)  = 0._r8
         pflx(i,k)  = 0._r8
         pflxg(i,k) = 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         qlg(i,k)   = 0._r8
         dlf(i,k)   = 0._r8
         dlg(i,k)   = 0._r8
      end do
      prec(i)        = 0._r8
      rliq(i)        = 0._r8
      pflx(i,pverp)  = 0
      pflxg(i,pverp) = 0
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
   ! calculate local pressure (mbs) and height (m) for both interface and mid-point

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
   do k = pver - 1,msg + 1,-1
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
         q(i,k)    = qh(i,k)
         s(i,k)    = t_mid(i,k) + (zm_const%grav/zm_const%cpair)*z_mid(i,k)
         tp(i,k)   = 0.0_r8
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
      end do
   end do

   do i = 1,ncol
      capeg(i)   = 0._r8
      lclg(i)    = 1
      lelg(i)    = pver
      maxg(i)    = 1
      tlg(i)     = 400._r8
      dsubcld(i) = 0._r8
   end do

   !----------------------------------------------------------------------------
   ! Evaluate Tparcel, qs(Tparcel), buoyancy, CAPE, lcl, lel, parcel launch level at index maxi()=hmax
   ! - call #1, iclosure=.true.   standard calculation using state of current step
   ! - call #2, iclosure=.false.  use state from previous step and launch level from call #1
   ! DCAPE is the difference in CAPE between the two calls using the same launch level

   iclosure = .true.
   call compute_dilute_cape( pcols, ncol, pver, pverp, &
                             zm_param%num_cin, msg, &
                             q, t_mid, z_mid, p_mid, p_int, &
                             pbl_top, tpert, &
                             tp, qstp, maxi, tl, &
                             lcl, lel, cape, &
                             zm_const, zm_param, &
                             iclosure )

   ! Calculate dcape trigger condition
   if ( .not.is_first_step .and. zm_param%trig_dcape ) then
      iclosure = .false.
      dcapemx(:ncol) = maxi(:ncol)
      call compute_dilute_cape( pcols, ncol, pver, pverp, &
                                zm_param%num_cin, msg, &
                                q_star, t_star, z_mid, p_mid, p_int, &
                                pbl_top, tpert, &
                                tpm1, qstpm1, maxim1, tlm1, &
                                lclm1, lelm1, capem1, &
                                zm_const, zm_param, &
                                iclosure, dcapemx )
      dcape(:ncol) = (cape(:ncol)-capem1(:ncol))/(delt*2._r8)
   endif

   !----------------------------------------------------------------------------
   ! determine whether ZM is active in each column (ideep=1) or not (ideep=0),
   ! based on requirement that cape>0 and lel<lcl

   capelmt_wk = capelmt   ! capelmt_wk default to capelmt for default trigger

   if ( zm_param%trig_dcape .and. (.not.is_first_step) )  capelmt_wk = 0.0_r8

   lengath = 0
   do i=1,ncol
     if (zm_param%trig_dcape) then
     ! DCAPE-ULL
      if (is_first_step) then
         !Will this cause restart to be non-BFB
           if (cape(i) > capelmt) then
              lengath = lengath + 1
              gather_index(lengath) = i
           end if
       else if (cape(i) > 0.0_r8 .and. dcape(i) > trigdcapelmt) then
           ! use constant 0 or a separate threshold for capt because capelmt is for default trigger
           lengath = lengath + 1
           gather_index(lengath) = i
       endif
     else
      if (cape(i) > capelmt) then
         lengath = lengath + 1
         gather_index(lengath) = i
      end if
     end if
   end do

   if (lengath.eq.0) then
      ! Deallocate local microphysics arrays before returning
      if (zm_param%zm_microp) call zm_microp_st_dealloc(loc_microp_st)
      return
   end if

   do ii=1,lengath
      ideep(ii)=gather_index(ii)
   end do

   !----------------------------------------------------------------------------
   ! copy data to gathered arrays

   do i = 1,lengath
      do k = 1,pver
         dp(i,k)        = 0.01_r8*p_del_in(ideep(i),k)
         qg(i,k)        = q(ideep(i),k)
         tg(i,k)        = t_mid(ideep(i),k)
         pg(i,k)        = p_mid(ideep(i),k)
         z_mid_g(i,k)   = z_mid(ideep(i),k)
         sg(i,k)        = s(ideep(i),k)
         tpg(i,k)       = tp(ideep(i),k)
         z_int_g(i,k)   = z_int(ideep(i),k)
         qstpg(i,k)     = qstp(ideep(i),k)
         omegag(i,k)    = omega(ideep(i),k)
         ug(i,k)        = 0._r8
         vg(i,k)        = 0._r8
      end do
      z_int_g(i,pverp)  = z_int(ideep(i),pver+1)
      capeg(i)          = cape(ideep(i))
      lclg(i)           = lcl(ideep(i))
      lelg(i)           = lel(ideep(i))
      maxg(i)           = maxi(ideep(i))
      tlg(i)            = tl(ideep(i))
      landfracg(i)      = landfrac(ideep(i))
      pbl_top_g(i)      = pbl_top(ideep(i))
      tpertg(i)         = tpert(ideep(i))
   end do

   !----------------------------------------------------------------------------
   ! copy aerosol data to gathered arrays

   if (zm_param%zm_microp) then
      if (aero%scheme == 'modal') then
         do m = 1, aero%nmodes
            do i = 1,lengath
               do k = 1,pver
                  aero%numg_a(i,k,m) = aero%num_a(m)%val(ideep(i),k)
                  aero%dgnumg(i,k,m) = aero%dgnum(m)%val(ideep(i),k)
                  do l = 1, aero%nspec(m)
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

   !----------------------------------------------------------------------------
   ! calculate sub-cloud layer pressure "thickness" for use in
   ! closure and tendency routines.

   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + dp(i,k)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! define interfacial values for (q,s) used in subsequent routines.
   do k = msg + 2,pver
      do i = 1,lengath
         sdifr = 0._r8
         qdifr = 0._r8
         if (sg(i,k) > 0._r8 .or. sg(i,k-1) > 0._r8) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._r8 .or. qg(i,k-1) > 0._r8) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > interp_diff_min) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > interp_diff_min) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! calculate updraft and downdraft properties

   call cldprp(pcols, lengath, pver, pverp, msg, zm_param%limcnv, &
               pg, z_mid_g, z_int_g, tg, sg, shat, qg, ug, vg, landfracg, tpertg, &
               maxg, lelg, jt, jlcl, j0, jd, &
               mu, eu, du, md, ed, mc, &
               su, qu, qlg, sd, qd,  &
               qs, cug, evpg, pflxg, rprdg, &
               aero, loc_microp_st )

   !----------------------------------------------------------------------------
   ! convert detrainment from units of "1/m" to "1/mb".

   do k = msg + 1,pver
      do i = 1,lengath
         du   (i,k) = du   (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         eu   (i,k) = eu   (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         ed   (i,k) = ed   (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         cug  (i,k) = cug  (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         rprdg(i,k) = rprdg(i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         evpg (i,k) = evpg (i,k)* (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         if (zm_param%zm_microp) then
            loc_microp_st%frz (i,k) = loc_microp_st%frz (i,k) * (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
            loc_microp_st%sprd(i,k) = loc_microp_st%sprd(i,k) * (z_int_g(i,k)-z_int_g(i,k+1))/dp(i,k)
         end if
      end do
   end do

   !----------------------------------------------------------------------------

   call closure(pcols, ncol, pver, pverp, &
                qg      ,tg      ,pg      ,z_mid_g      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,z_int_g     , &
                qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,msg     ,capelmt_wk )

   !----------------------------------------------------------------------------
   ! limit cloud base mass flux to theoretical upper bound.

   do i=1,lengath
      mumax(i) = 0
      do k=msg + 2,pver
        mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
      if (mumax(i) > 0._r8) then
         mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
      else
         mb(i) = 0._r8
      endif
      if (zm_param%clos_dyn_adj) mb(i) = max(mb(i) - omegag(i,pbl_top_g(i))*0.01_r8, 0._r8)
   end do

   !----------------------------------------------------------------------------
   ! don't allow convection within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (zm_param%no_deep_pbl) then
      do i=1,lengath
         if (z_mid_in(ideep(i),jt(i)) < pbl_hgt(ideep(i))) mb(i) = 0
      end do
   end if

   !----------------------------------------------------------------------------
   ! apply cloud base mass flux scaling

   do i=1,lengath

      ! zero out micro data for inactive columns
      if ( zm_param%zm_microp .and. mb(i).eq.0._r8) call zm_microp_st_zero(loc_microp_st,i,pver)

      do k=msg+1,pver
         mu   (i,k)  = mu   (i,k)*mb(i)
         md   (i,k)  = md   (i,k)*mb(i)
         mc   (i,k)  = mc   (i,k)*mb(i)
         du   (i,k)  = du   (i,k)*mb(i)
         eu   (i,k)  = eu   (i,k)*mb(i)
         ed   (i,k)  = ed   (i,k)*mb(i)
         rprdg(i,k)  = rprdg(i,k)*mb(i)
         cug  (i,k)  = cug  (i,k)*mb(i)
         evpg (i,k)  = evpg (i,k)*mb(i)
         pflxg(i,k+1)= pflxg(i,k+1)*mb(i)*100._r8/zm_const%grav

         if (zm_param%zm_microp) then
            loc_microp_st%sprd(i,k)  = loc_microp_st%sprd(i,k)*mb(i)
            loc_microp_st%frz (i,k)  = loc_microp_st%frz (i,k)*mb(i)
            if (mb(i).eq.0._r8) then
               qlg (i,k)  = 0._r8
            end if
         end if

      end do
   end do

   !----------------------------------------------------------------------------
   ! compute temperature and moisture changes due to convection.

   call zm_calculate_output_tend(pcols, ncol, pver, pverp, &
                                 dqdt, dsdt, qg, qs, qu, &
                                 su, du, qhat, shat, dp, &
                                 mu, md, sd, qd, qlg, &
                                 dsubcld, jt, maxg, 1, lengath, msg, &
                                 dlg, evpg, cug, &
                                 loc_microp_st)

   !----------------------------------------------------------------------------
   ! conservation check and adjusment
#ifndef SCREAM_CONFIG_IS_CMAKE
   if (zm_param%zm_microp) call zm_microphysics_adjust(pcols, lengath, pver, jt, msg, delt, zm_const, &
                                                       dp, qg, dlg, dsdt, dqdt, rprdg, loc_microp_st)
#endif

   !----------------------------------------------------------------------------
   ! scatter data (i.e. undo the gathering)

   do k = msg + 1,pver
#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif
      do i = 1,lengath
         ! q is updated to compute net precip.
         q(ideep(i),k) = qh(ideep(i),k) + 2._r8*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*zm_const%cpair
         dlf (ideep(i),k) = dlg  (i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql  (ideep(i),k) = qlg  (i,k)
      end do
   end do

   do i = 1,lengath
      jctop(ideep(i)) = jt(i)
      jcbot(ideep(i)) = maxg(i)
      pflx(ideep(i),pverp) = pflxg(i,pverp)
   end do

   !----------------------------------------------------------------------------
   ! scatter microphysics data (i.e. undo the gathering)

   if (zm_param%zm_microp) call zm_microp_st_scatter(loc_microp_st,microp_st,pcols,lengath,pver,ideep)

#ifdef CPRCRAY
!DIR$ CONCURRENT
#endif

   !----------------------------------------------------------------------------
   ! Compute precip by integrating change in water vapor minus detrained cloud water
   do i = 1,ncol
      do k = pver,msg + 1,-1
         if (zm_param%zm_microp) then
            prec(i) = prec(i) - p_del_in(i,k)*(q(i,k)-qh(i,k)) - p_del_in(i,k)*(dlf(i,k)+microp_st%dif(i,k)+microp_st%dsf(i,k))*2._r8*delt
         else
            prec(i) = prec(i) - p_del_in(i,k)*(q(i,k)-qh(i,k)) - p_del_in(i,k)*(dlf(i,k))*2._r8*delt
         end if
      end do
      ! obtain final precipitation rate in m/s
      prec(i) = zm_const%rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
   end do

   !----------------------------------------------------------------------------
   ! Compute reserved liquid (and ice) that is not yet in cldliq for energy integrals
   ! Treat rliq as flux out bottom, to be added back later
   do k = 1, pver
      do i = 1, ncol
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
end subroutine zm_convr

!===================================================================================================

subroutine zm_conv_evap(pcols, ncol, pver, pverp, deltat, &
                        pmid, pdel, t_mid, q, prdprec, cldfrc, &
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
   real(r8),                        intent(in   ) :: deltat             ! time step                               [s]
   real(r8), dimension(pcols,pver), intent(in   ) :: pmid               ! midpoint pressure                       [Pa]
   real(r8), dimension(pcols,pver), intent(in   ) :: pdel               ! layer thickness                         [Pa]
   real(r8), dimension(pcols,pver), intent(in   ) :: t_mid              ! temperature                             [K]
   real(r8), dimension(pcols,pver), intent(in   ) :: q                  ! water vapor                             [kg/kg]
   real(r8), dimension(pcols,pver), intent(in   ) :: prdprec            ! precipitation production                [kg/kg/s]
   real(r8), dimension(pcols,pver), intent(in   ) :: cldfrc             ! cloud fraction
   real(r8), dimension(pcols,pver), intent(inout) :: tend_s             ! heating rate                            [J/kg/s]
   real(r8), dimension(pcols,pver), intent(inout) :: tend_q             ! water vapor tendency                    [kg/kg/s]
   real(r8), dimension(pcols,pver), intent(out  ) :: tend_s_snwprd      ! Heating rate of snow production         [J/kg/s]
   real(r8), dimension(pcols,pver), intent(out  ) :: tend_s_snwevmlt    ! Heating rate of snow evap/melt          [J/kg/s]
   real(r8), dimension(pcols),      intent(inout) :: prec(pcols)        ! Convective-scale prec rate              [m/s]
   real(r8), dimension(pcols),      intent(out  ) :: snow(pcols)        ! Convective-scale snow rate              [m/s]
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
   !----------------------------------------------------------------------------
   if (zm_param%zm_microp) then
      prdsnow(1:ncol,1:pver) = microp_st%sprd(1:ncol,1:pver)
   else
      prdsnow(1:ncol,1:pver) = 0._r8
   end if

   ! convert input precip to kg/m2/s
   prec(:ncol) = prec(:ncol)*1000._r8

   ! determine saturation vapor pressure
   call qsat(t_mid(1:ncol, 1:pver), pmid(1:ncol, 1:pver), es(1:ncol, 1:pver), qs(1:ncol, 1:pver))

   ! determine ice fraction in rain production (use cloud water parameterization fraction at present)
   call cldfrc_fice(ncol, t_mid, fice, fsnow_conv)

   ! zero the flux integrals on the top boundary
   flxprec(:ncol,1) = 0._r8
   flxsnow(:ncol,1) = 0._r8
   evpvint(:ncol)   = 0._r8

   do k = 1, pver
      do i = 1, ncol

         ! Melt snow falling into layer, if necessary.
         if( zm_param%old_snow ) then
            if (t_mid(i,k) > zm_const%tfreez) then
               flxsntm(i) = 0._r8
               snowmlt(i) = flxsnow(i,k) * zm_const%grav/ pdel(i,k)
            else
               flxsntm(i) = flxsnow(i,k)
               snowmlt(i) = 0._r8
            end if
         else
            ! make sure melting snow doesn't reduce temperature below threshold
            if (t_mid(i,k) > zm_const%tfreez) then
               dum = -zm_const%latice/zm_const%cpair*flxsnow(i,k)*zm_const%grav/pdel(i,k)*deltat
               if (t_mid(i,k) + dum .le. zm_const%tfreez) then
                  dum = (t_mid(i,k)-zm_const%tfreez)*zm_const%cpair/zm_const%latice/deltat
                  dum = dum/(flxsnow(i,k)*zm_const%grav/pdel(i,k))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if
               dum = dum*omsm
               flxsntm(i) = flxsnow(i,k)*(1.0_r8-dum)
               snowmlt(i) = dum*flxsnow(i,k)*zm_const%grav/ pdel(i,k)
            else
               flxsntm(i) = flxsnow(i,k)
               snowmlt(i) = 0._r8
            end if
         end if

         ! relative humidity depression must be > 0 for evaporation
         evplimit = max(1._r8 - q(i,k)/qs(i,k), 0._r8)

         ! total evaporation depends on flux in the top of the layer
         ! flux prec is the net production above layer minus evaporation into environmet
         evpprec(i) = zm_param%ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))

         ! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
         ! Currently does not include heating/cooling change to qs
         evplimit = max(0._r8, (qs(i,k)-q(i,k)) / deltat)

         ! Don't evaporate more than is falling into the layer from above.
         ! Don't evaporate rain formed in this layer, but if precip production
         ! is negative, remove from the available precip. Negative precip
         ! production occurs because of evaporation in downdrafts.
         evplimit = min(evplimit, flxprec(i,k) * zm_const%grav / pdel(i,k))

         ! Total evaporation cannot exceed input precipitation
         evplimit = min(evplimit, (prec(i) - evpvint(i)) * zm_const%grav / pdel(i,k))

         evpprec(i) = min(evplimit, evpprec(i))

         if( .not.zm_param%old_snow ) then
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
         evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/zm_const%grav

         ! net precip production => production - evaporation
         ntprprd(i,k) = prdprec(i,k) - evpprec(i)

         ! net snow production => precip production * ice fraction - evaporation - melting

         ! the small amount added to flxprec in the work1 expression was increased
         ! from 1e-36 to 8.64e-11 (1e-5 mm/day) to address error growth problems.
         ! This causes temperature partitioning to be used for small flxprec amounts.

         if( zm_param%old_snow ) then
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
            ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
            tend_s_snwprd  (i,k) = prdprec(i,k)*work2*zm_const%latice
            tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*zm_const%latice
         else
            ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k)*zm_const%grav/pdel(i,k), evpsnow(i)+snowmlt(i))
            tend_s_snwprd  (i,k) = prdsnow(i,k)*zm_const%latice
            tend_s_snwevmlt(i,k) = -min(flxsnow(i,k)*zm_const%grav/pdel(i,k), evpsnow(i)+snowmlt(i) )*zm_const%latice
         end if

         ! precipitation fluxes
         flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/zm_const%grav
         flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/zm_const%grav

         ! protect against rounding error
         flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
         flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)

         ! heating (cooling) and moistening due to evaporation
         ! - latent heat of vaporization for precip production has already been accounted for
         ! - snow is contained in prec
         if( zm_param%old_snow ) then
            tend_s(i,k) =-evpprec(i)*zm_const%latvap + ntsnprd(i,k)*zm_const%latice
         else
            tend_s(i,k) =-evpprec(i)*zm_const%latvap + tend_s_snwevmlt(i,k)
         end if

         tend_q(i,k) = evpprec(i)
      end do ! i
   end do ! k

   ! protect against rounding error
   if( .not.zm_param%old_snow ) then
      do i = 1, ncol
         if(flxsnow(i,pverp).gt.flxprec(i,pverp)) then
            dum = (flxsnow(i,pverp)-flxprec(i,pverp))*zm_const%grav
            do k = pver, 1, -1
               if (ntsnprd(i,k)>ntprprd(i,k).and. dum > 0._r8) then
                  ntsnprd(i,k) = ntsnprd(i,k) - dum/pdel(i,k)
                  tend_s_snwevmlt(i,k) = tend_s_snwevmlt(i,k) - dum/pdel(i,k)*zm_const%latice
                  tend_s(i,k)  = tend_s(i,k) - dum/pdel(i,k)*zm_const%latice
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

subroutine calculate_fractional_entrainment(pcols, ncol, pver, pverp, msg, &
                                            jb, jt, j0, z_mid, z_int, dz, &
                                            h_env, h_env_sat, h_env_min, eps, eps0)
   !----------------------------------------------------------------------------
   ! Purpose: Determine properties of ZM updrafts and downdrafts
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols     ! maximum number of columns
   integer,                         intent(in   ) :: ncol      ! actual number of columns
   integer,                         intent(in   ) :: pver      ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp     ! number of interface vertical levels
   integer,                         intent(in   ) :: msg       ! missing moisture levels
   integer,  dimension(pcols),      intent(in   ) :: jb        ! updraft base level
   integer,  dimension(pcols),      intent(in   ) :: jt        ! updraft top level
   integer,  dimension(pcols),      intent(inout) :: j0        ! level where updraft begins detraining
   real(r8), dimension(pcols,pver), intent(in   ) :: z_mid     ! env altitude at mid-point
   real(r8), dimension(pcols,pverp),intent(in   ) :: z_int     ! env altitude at interface
   real(r8), dimension(pcols,pver) ,intent(in   ) :: dz        ! layer thickness
   real(r8), dimension(pcols,pver), intent(in   ) :: h_env     ! env moist stat energy
   real(r8), dimension(pcols,pver), intent(in   ) :: h_env_sat ! env saturated moist stat energy
   real(r8), dimension(pcols)     , intent(inout) :: h_env_min ! mid-tropospheric MSE minimum
   real(r8), dimension(pcols,pver), intent(  out) :: eps       ! fractional entrainment
   real(r8), dimension(pcols)     , intent(  out) :: eps0      ! fractional entrainment maximum
   !----------------------------------------------------------------------------
   ! Local variables

   integer :: i,k ! loop iterators

   ! variables used for Taylor series expansion when solving eq (4.78) for lamda_D (i.e. fractional entrainment)
   real(r8), dimension(pcols,pver) :: f         ! fractional entrainment work variable
   real(r8), dimension(pcols)      :: ftemp     ! temporary value of f
   real(r8), dimension(pcols)      :: expdif    ! diff between env MSE low-level max and mid-level min
   real(r8), dimension(pcols)      :: expnum    ! ?
   real(r8), dimension(pcols,pver) :: k1        ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i2        ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i3        ! term for Taylor series
   real(r8), dimension(pcols,pver) :: i4        ! term for Taylor series
   real(r8), dimension(pcols,pver) :: ihat      ! term for Taylor series
   real(r8), dimension(pcols,pver) :: idag      ! term for Taylor series
   real(r8), dimension(pcols,pver) :: iprm      ! term for Taylor series

   real(r8), parameter :: lambda_min = 0._r8     ! limiter
   real(r8), parameter :: lambda_max = 0.0002_r8 ! limiter
   real(r8), parameter :: lambda_thrs = 1.E-6_r8 ! threshold
   !----------------------------------------------------------------------------

   ! initialize 1D variables
   do i = 1,ncol
      ftemp(i)  = 0._r8
      expnum(i) = 0._r8
      expdif(i) = 0._r8
   end do

   ! initialize 2D variables
   do k = 1,pver
      do i = 1,ncol
         ! misc work variables
         k1(i,k)     = 0._r8
         i2(i,k)     = 0._r8
         i3(i,k)     = 0._r8
         i4(i,k)     = 0._r8
         f(i,k)      = 0._r8
         eps(i,k)    = 0._r8
      end do
   end do

   !----------------------------------------------------------------------------
   ! compute taylor series for approximate eps(z) below
   do k = pver - 1,msg + 1,-1
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
      end do
   end do

   !----------------------------------------------------------------------------
   ! re-initialize minimum MSE for ensuing calculation
   h_env_min(1:ncol) = 1.E6_r8
   do k = msg + 1,pver
      do i = 1,ncol
         if (k >= j0(i) .and. k <= jb(i) .and. h_env(i,k) <= h_env_min(i)) then
            h_env_min(i) = h_env(i,k)
            expdif(i) = h_env(i,jb(i)) - h_env_min(i)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! compute approximate eps(z) using above taylor series
   do k = msg + 2,pver
      do i = 1,ncol
         expnum(i)   = 0._r8
         ftemp(i)    = 0._r8
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k)  = 0._r8
            expnum(i)= 0._r8
         else
            expnum(i) = h_env(i,jb(i)) - (h_env_sat(i,k-1)*(z_int(i,k)-z_mid(i,k)) + &
                        h_env_sat(i,k)*(z_mid(i,k-1)-z_int(i,k)))/(z_mid(i,k-1)-z_mid(i,k))
         end if
         if ((expdif(i) > 100._r8 .and. expnum(i) > 0._r8) .and. k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i) / k1(i,k)
            f(i,k) = ftemp(i) + &
                     i2(i,k)/k1(i,k) * ftemp(i)**2 + &
                     (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2 * ftemp(i)**3 + &
                     (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ 5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/k1(i,k)**3 * ftemp(i)**4
            f(i,k) = max(f(i,k),lambda_min)
            f(i,k) = min(f(i,k),lambda_max)
         end if
      end do
   end do

   ! move detrainment level downward if fractional entrainment is too low
   do i = 1,ncol
      if (j0(i) < jb(i)) then
         if ( f(i,j0(i))<1.E-6_r8 .and. f(i,j0(i)+1)>f(i,j0(i)) ) then
            j0(i) = j0(i) + 1
         end if
      end if
   end do

   ! ensure that entrainment does not increase above the level that detrainment starts
   do k = msg + 2,pver
      do i = 1,ncol
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do

   ! specify maximum fractional entrainment
   do i = 1,ncol
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do

   ! The modification below comes from:
   !   Rasch, P. J., J. E. Kristjánsson, A comparison of the CCM3 model climate
   !   using diagnosed and predicted condensate parameterizations, J. Clim., 1997.
   do k = pver,msg + 1,-1
      do i = 1,ncol
         if ( k >=j0(i) .and. k<=jb(i) ) eps(i,k) = f(i,j0(i))
         if ( k  <j0(i) .and. k>=jt(i) ) eps(i,k) = f(i,k)
      end do
   end do

   !----------------------------------------------------------------------------
   return

end subroutine calculate_fractional_entrainment

!===================================================================================================

subroutine cldprp(pcols, ncol, pver, pverp, msg, limcnv, &
                  p_mid, z_mid, z_int, t_mid, s, shat, q, u, v, landfrac, tpertg, &
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
   integer,                         intent(in ) :: msg            ! missing moisture levels
   integer,                         intent(in ) :: limcnv         ! convection limiting level
   real(r8), dimension(pcols,pver), intent(in ) :: p_mid          ! env pressure at mid-point
   real(r8), dimension(pcols,pver), intent(in ) :: z_mid          ! env altitude at mid-point
   real(r8), dimension(pcols,pverp),intent(in ) :: z_int          ! env altitude at interface
   real(r8), dimension(pcols,pver), intent(in ) :: t_mid          ! env temperature
   real(r8), dimension(pcols,pver), intent(in ) :: s              ! env dry static energy of env [K] (normalized)
   real(r8), dimension(pcols,pver), intent(in ) :: shat           ! interface values of dry stat energy
   real(r8), dimension(pcols,pver), intent(in ) :: q              ! env specific humidity
   real(r8), dimension(pcols,pver), intent(in ) :: u              ! env zonal wind
   real(r8), dimension(pcols,pver), intent(in ) :: v              ! env meridional wind
   real(r8), dimension(pcols),      intent(in ) :: landfrac       ! Land fraction
   real(r8), dimension(pcols),      intent(in ) :: tpertg         ! PBL temperature perturbation
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
   real(r8), dimension(pcols,pver), intent(out) :: ed             ! entrainment rate of downdraft
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
   real(r8), dimension(pcols,pver) :: gamma     ! latent-heating correction for pseudo-adiabatic parcel lifting
   real(r8), dimension(pcols,pver) :: dz        ! layer thickness
   real(r8), dimension(pcols,pver) :: h_env     ! ambient env moist stat energy
   real(r8), dimension(pcols,pver) :: h_env_sat ! ambient env saturated moist stat energy
   real(r8), dimension(pcols)      :: h_env_min ! mid-tropospheric MSE minimum
   real(r8), dimension(pcols,pver) :: hu        ! updraft moist static energy
   real(r8), dimension(pcols,pver) :: hd        ! dndraft moist static energy
   real(r8), dimension(pcols,pver) :: qsthat    ! interface interpolated qst
   real(r8), dimension(pcols,pver) :: hsthat    ! interface interpolated hst
   real(r8), dimension(pcols,pver) :: gamhat    ! interface interpolated gamma
   real(r8), dimension(pcols,pver) :: qds       ! dndraft saturation specific humdity
   real(r8), dimension(pcols)      :: c0mask    ! land/ocean modification

   real(r8), dimension(pcols)      :: rmue      !
   real(r8), dimension(pcols)      :: zuef      !
   real(r8), dimension(pcols)      :: zdef      !
   real(r8), dimension(pcols)      :: epsm      !
   real(r8), dimension(pcols)      :: ratmjb    !
   real(r8), dimension(pcols)      :: est       !
   real(r8), dimension(pcols)      :: totpcp    ! total precip for dndraft proportionality factor - see eq (4.106)
   real(r8), dimension(pcols)      :: totevp    ! total evap   for dndraft proportionality factor - see eq (4.106)

   real(r8), dimension(pcols,pver) :: eps       ! fractional entrainment
   real(r8), dimension(pcols)      :: eps0      ! fractional entrainment max

   real(r8) ql1   ! ?
   real(r8) tu    ! ?
   real(r8) estu  ! ?
   real(r8) qstu  ! ?
   real(r8) mdt   ! ?

   ! Convective microphysics
   real(r8), dimension(pcols,pver) :: fice       ! ice fraction in precip production
   real(r8), dimension(pcols,pver) :: tug        ! temporary updraft temperature
   real(r8), dimension(pcols,pver) :: tmp_frz    ! temporary rate of freezing
   real(r8), dimension(pcols)      :: tot_frz    ! total column freezing rate
   real(r8), dimension(pcols,pverp):: pflxs      ! frozen precipitation flux thru layer
   real(r8) dum, sdum

   integer, dimension(pcols) :: jto              ! updraft plume old top

   integer :: tmp_k_limit  ! temporary limit on k index in various contexts
   integer :: iter         ! ?
   integer :: itnum        ! ?
   integer :: m

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical, dimension(pcols) :: doit ! flag to reset cloud
   logical, dimension(pcols) :: done ! ?

   real(r8), parameter :: small        = 1.e-20_r8    ! small number to limit blowup when normalizing by mass flux
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
         call qsat_hPa(t_mid(i,k), p_mid(i,k), est(i), qst(i,k))
         if ( p_mid(i,k)-est(i) <= 0._r8 ) qst(i,k) = 1.0_r8
         ! compute gamma - see eq. (4.117)
         gamma(i,k)  = qst(i,k)*(1._r8 + qst(i,k)/zm_const%epsilo) &
                       * zm_const%epsilo*zm_const%latvap/(zm_const%rdair*t_mid(i,k)**2) &
                       * zm_const%latvap/zm_const%cpair
         ! cloud thermodynamic variables
         su(i,k)        = s(i,k)
         sd(i,k)        = s(i,k)
         qd(i,k)        = q(i,k)
         qds(i,k)       = q(i,k)
         qu(i,k)        = q(i,k)
         h_env(i,k)     = zm_const%cpair*t_mid(i,k) + zm_const%grav*z_mid(i,k) + zm_const%latvap*q(i,k)
         h_env_sat(i,k) = zm_const%cpair*t_mid(i,k) + zm_const%grav*z_mid(i,k) + zm_const%latvap*qst(i,k)
         hu(i,k)        = h_env(i,k)
         hd(i,k)        = h_env(i,k)
         ! convective microphysics
         pflxs(i,k)  = 0._r8
         fice(i,k)   = 0._r8
         tug(i,k)    = 0._r8
         tmp_frz(i,k)= 0._r8
      end do
   end do

   !----------------------------------------------------------------------------
   ! interpolate the mid-point values to interfaces
   do k = msg+1,pver
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
            hsthat(i,k) = zm_const%cpair*shat(i,k) + zm_const%latvap*qsthat(i,k)
            if (abs(gamma(i,k-1)-gamma(i,k)) > interp_diff_min) then
               gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ (gamma(i,k-1)-gamma(i,k))
            else
               gamhat(i,k) = gamma(i,k)
            end if
         end if
      end do
   end do

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
   end do

   !----------------------------------------------------------------------------
   ! find the level of minimum saturated MSE (h_env_sat), where detrainment starts
   do i = 1,ncol
      do k = msg + 1,pver
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
   do k = msg + 1,pver
      do i = 1,ncol
         if (k >= jt(i) .and. k <= jb(i)) then
            ! Tunable temperature perturbation (tiedke_add) was already added to parcel hu/su to
            ! represent subgrid temperature perturbation. If PBL temperature pert (tpert) also
            ! represents subgrid temperature pert, tiedke_add may need to be removed. Additionally,
            ! current calculation of PBL temperature perturbation is not accurate enough so that a
            ! new tunable parameter tpert_fac was introduced. This introduced new uncertainties into
            ! the ZM scheme. The original code of ZM scheme will be used when tpert_fix=.true.
            if(zm_param%tpert_fix) then
               hu(i,k) = h_env(i,jb(i)) + zm_const%cpair*zm_param%tiedke_add
               su(i,k) = s(i,jb(i))     +                zm_param%tiedke_add
            else
               hu(i,k) = h_env(i,jb(i)) + zm_const%cpair*(zm_param%tiedke_add+zm_param%tpert_fac*tpertg(i))
               su(i,k) = s(i,jb(i))     +                 zm_param%tiedke_add+zm_param%tpert_fac*tpertg(i)
            end if
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! calculate fractional entrainment (i.e. "lambda D") - see eq (4.78) from Neale et al. (2012)
   call calculate_fractional_entrainment(pcols, ncol, pver, pverp, msg, &
                                         jb, jt, j0, z_mid, z_int, dz, &
                                         h_env, h_env_sat, h_env_min, eps, eps0)

   !----------------------------------------------------------------------------
   ! iteration to set cloud properties

   itnum = 1
   if (zm_param%zm_microp) itnum = 2

   do iter = 1,itnum

      do k = pver,msg + 1,-1
         do i = 1,ncol
            cu(i,k) = 0._r8
            ql(i,k) = 0._r8
            if (zm_param%zm_microp) then
               loc_microp_st%qliq(i,k) = 0._r8
               loc_microp_st%qice(i,k) = 0._r8
               loc_microp_st%frz(i,k)  = 0._r8
            end if
         end do
      end do
      do i = 1,ncol
         totpcp(i) = 0._r8
         if (zm_param%zm_microp) hu(i,jb(i)) = h_env(i,jb(i)) + zm_const%cpair*zm_param%tiedke_add
      end do

      do k = pver,msg + 1,-1
         do i = 1,ncol
            ! intialize updraft mass flux variables - here and below all normalized by cloud base mass flux (mb)
            if (eps0(i) > 0._r8) then
               mu(i,jb(i)) = 1._r8
               eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
               if (     zm_param%zm_microp) tmp_k_limit = lel(i)
               if (.not.zm_param%zm_microp) tmp_k_limit = jt(i)
               ! compute profiles of updraft mass fluxes - see eq (4.79) - (4.81)
               if ( k>=tmp_k_limit .and. k<jb(i) ) then
                  zuef(i) = z_int(i,k) - z_int(i,jb(i))
                  rmue(i) = (1._r8/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._r8)/zuef(i)
                  mu(i,k) = (1._r8/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._r8)/zuef(i)
                  eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
                  du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
               end if
            end if ! eps0(i)>0._r8
         end do
      end do

      khighest = pverp
      klowest = 1
      do i=1,ncol
         khighest = min(khighest,lel(i))
         klowest = max(klowest,jb(i))
      end do

      do k = klowest-1,khighest,-1
         do i = 1,ncol
            if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0._r8) then
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
         end do
      end do

      ! reset cloud top index beginning from two layers above the
      ! cloud base (i.e. if cloud is only one layer thick, top is not reset
      do i=1,ncol
         doit(i) = .true.
         tot_frz(i)= 0._r8
         do k = pver,msg + 1,-1
            tot_frz(i)= tot_frz(i) + tmp_frz(i,k)*dz(i,k)
         end do
      end do

      do k=klowest-2,khighest-1,-1
         do i=1,ncol
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
         end do
      end do

      do i = 1,ncol
         if (iter == 1) jto(i) = jt(i)
      end do

      do k = pver,msg + 1,-1
         do i = 1,ncol
            if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0._r8) then
               mu(i,k) = 0._r8
               eu(i,k) = 0._r8
               du(i,k) = 0._r8
               hu(i,k) = h_env(i,k)
            end if
            if (k == jt(i) .and. eps0(i) > 0._r8) then
               du(i,k) = mu(i,k+1)/dz(i,k)
               eu(i,k) = 0._r8
               mu(i,k) = 0._r8
            end if
         end do
      end do

      ! determine LCL - see eq (4.127)- (4.130) of Neale et al. (2012)
      done(1:ncol) = .false.
      kount = 0
      do k = pver,msg + 2,-1
         do i = 1,ncol
            if (k == jb(i) .and. eps0(i) > 0._r8) then
               qu(i,k) = q(i,jb(i))
               su(i,k) = (hu(i,k)-zm_const%latvap*qu(i,k))/zm_const%cpair
            end if
            if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0._r8) then
               su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
               qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k) - du(i,k)*qst(i,k))
               tu = su(i,k) - zm_const%grav/zm_const%cpair*z_int(i,k)
               call qsat_hPa(tu, (p_mid(i,k)+p_mid(i,k-1))/2._r8, estu, qstu)
               if (qu(i,k) >= qstu) then
                  jlcl(i) = k
                  kount = kount + 1
                  done(i) = .true.
               end if
            end if
         end do
         if (kount >= ncol) goto 690
      end do

690 continue

      do k = msg + 2,pver
         do i = 1,ncol
            if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
               su(i,k) = shat(i,k)   +             (hu(i,k)-hsthat(i,k)) / (zm_const%cpair* (1._r8+gamhat(i,k)))
               qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k)) / (zm_const%latvap* (1._r8+gamhat(i,k)))
            end if
         end do
      end do

      ! compute condensation in updraft
      do k = pver,msg + 2,-1
         do i = 1,ncol
            if (eps0(i)>0._r8) then
               if (     zm_param%zm_microp) tmp_k_limit = jlcl(i)+1
               if (.not.zm_param%zm_microp) tmp_k_limit = jb(i)
               if ( k>=jt(i) .and. k<tmp_k_limit ) then
                  if (zm_param%zm_microp) then
                     cu(i,k) = ( ( mu(i,k)*su(i,k) - mu(i,k+1)*su(i,k+1) )/dz(i,k) - eu(i,k)*s(i,k) + du(i,k)*su(i,k) &
                               )/(zm_const%latvap/zm_const%cpair) - zm_const%latice*tmp_frz(i,k)/zm_const%latvap
                  else
                     cu(i,k) = ( ( mu(i,k)*su(i,k) - mu(i,k+1)*su(i,k+1) )/dz(i,k) - ( eu(i,k) - du(i,k) )*s(i,k) &
                               )/(zm_const%latvap/zm_const%cpair)
                  end if
                  ! apply limiters
                  if (k == jt(i)) cu(i,k) = 0._r8
                  cu(i,k) = max(0._r8,cu(i,k))
               end if
            end if ! eps0(i)>0._r8
         end do ! i
      end do ! k

      if (zm_param%zm_microp) then

         tug(1:ncol,:) = t_mid(1:ncol,:)
         fice(1:ncol,:) = 0._r8

         do k = pver, msg+2, -1
            do i = 1, ncol
               tug(i,k) = su(i,k) - zm_const%grav/zm_const%cpair*z_int(i,k)
            end do
         end do

         do k = 1, pver-1
            do i = 1, ncol
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
            end do
         end do

         do k = 1, pver
            do i = 1,ncol
               loc_microp_st%cmei(i,k) = cu(i,k)* fice(i,k)
               loc_microp_st%cmel(i,k) = cu(i,k) * (1._r8-fice(i,k))
            end do
         end do

#ifndef SCREAM_CONFIG_IS_CMAKE
         call  zm_mphy( pcols, ncol, msg, &
                        zm_const%grav, zm_const%cpair, zm_const%rdair, &
                        zm_param%auto_fac, zm_param%accr_fac, zm_param%micro_dcs, &
                        jb, jt, jlcl, su, qu, mu, du, eu, z_int, p_mid, t_mid, q, gamhat, eps0, &
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

         do k = pver,msg + 2,-1
            do i = 1,ncol
               ! In the original ZM scheme, which does not consider ice phase, ql actually represents total cloud
               ! water. With convective microphysics, loc_microp_st%qliq and loc_microp_st%qice represent cloud
               ! liquid water and cloud ice, respectively. Since ql is still used in other subroutines as total
               ! cloud water, here ql is calculated as total cloud water for consistency.
               ql(i,k) = loc_microp_st%qliq(i,k)+ loc_microp_st%qice(i,k)
               loc_microp_st%frz(i,k) = tmp_frz(i,k)
            end do
         end do

         do i = 1,ncol
           if (iter == 2 .and. jt(i)> jto(i)) then
             do k = jt(i), jto(i), -1
                loc_microp_st%frz(i,k) = 0.0_r8
                cu(i,k)=0.0_r8
             end do
           end if
         end do

         do k = pver,msg + 2,-1
            do i = 1,ncol
               if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
                  totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*( loc_microp_st%qcde(i,k+1) &
                                                                    +loc_microp_st%qide(i,k+1) &
                                                                    +loc_microp_st%qsde(i,k+1) ))
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
            do i = 1,ncol
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
               end if
            end do
         end do

      end if  ! zm_param%zm_microp

   end do ! iter = 1,itnum

   !----------------------------------------------------------------------------
   ! specify downdraft properties (note - no downdrafts if jd>=jb)
   ! downward mass flux profile is scaled so that net flux (up-down) at cloud base in not negative
   do i = 1,ncol
      ! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = h_env(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0._r8) then
         epsm(i) = eps0(i)
         md(i,jd(i)) = -zm_param%alfa * epsm(i) / eps0(i)
      end if
   end do
   do k = msg + 1,pver
      do i = 1,ncol
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8) then
            zdef(i) = z_int(i,jd(i)) - z_int(i,k)
            md(i,k) = -zm_param%alfa / (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i)
         end if
      end do
   end do
   do k = msg + 1,pver
      do i = 1,ncol
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)
         end if
      end do
   end do

   do k = msg + 1,pver
      do i = 1,ncol
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0._r8) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*h_env(i,k-1))/mdt
         end if
      end do
   end do

   ! calculate updraft and downdraft properties
   do k = msg + 2,pver
      do i = 1,ncol
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k)) / (zm_const%latvap*(1._r8+gamhat(i,k)))
         end if
      end do
   end do

   do i = 1,ncol
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - zm_const%latvap*qd(i,jd(i)))/zm_const%cpair
   end do

   do k = msg + 2,pver
      do i = 1,ncol
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            mdt = min(md(i,k+1),-small)
            if (zm_param%zm_microp)   evp(i,k) = min(evp(i,k),rprd(i,k))
            sd(i,k+1) = ((zm_const%latvap/zm_const%cpair*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do

   ! calculate total precip for downdraft scaling
   do i = 1,ncol
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do

   do k = msg + 2,pver
      do i = 1,ncol
         ! ensure that downdraft strength is consistent with precipitation availability - see eq (4.106)
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
      end do
   end do

   ! compute the net precipitation flux across interfaces
   do k = 2,pverp
      do i = 1,ncol
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
         if (zm_param%zm_microp) pflxs(i,k) = pflxs(i,k-1) + loc_microp_st%sprd(i,k-1)*dz(i,k-1)
      end do
   end do

   ! calculate net mass flux
   do k = msg + 1,pver
      do i = 1,ncol
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do

   if (zm_param%zm_microp) then
      do i = 1,ncol
         ! protect against rounding error
         if(pflxs(i,pverp).gt.pflx(i,pverp)) then
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
            do k = msg + 1,pver
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
   end if

   return
end subroutine cldprp

!===================================================================================================

subroutine closure(pcols, ncol, pver, pverp, &
                   q       ,t_mid   ,p_mid       ,z_mid   ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,z_int   , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,msg     ,capelmt )
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

   implicit none

!
!-----------------------------Arguments---------------------------------
!
   integer,          intent(in) :: pcols        ! maximum number of columns
   integer,          intent(in) :: ncol         ! actual number of columns
   integer,          intent(in) :: pver         ! number of mid-point vertical levels
   integer,          intent(in) :: pverp        ! number of interface vertical levels
   real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
   real(r8), intent(inout) :: t_mid(pcols,pver)        ! temperature
   real(r8), intent(inout) :: p_mid(pcols,pver)        ! pressure (mb)
   real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
   real(r8), intent(in) :: z_mid(pcols,pver)        ! height (m)
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
   real(r8), intent(in) :: z_int(pcols,pver+1)     ! height of interface levels
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
   real(r8) dadt(pcols)
   real(r8) debdt
   real(r8) dltaa
   real(r8) eb

   integer i
   integer il1g
   integer il2g
   integer k, kmin, kmax
   integer msg

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
      eb = p_mid(i,mx(i))*q(i,mx(i))/ (zm_const%epsilo+q(i,mx(i)))
      dtbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
      dqbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
                 md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
      debdt = zm_const%epsilo*p_mid(i,mx(i))/ (zm_const%epsilo+q(i,mx(i)))**2*dqbdt(i)
      dtldt(i) = -2840._r8* (3.5_r8/t_mid(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t_mid(i,mx(i)))-log(eb)-4.805_r8)**2
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
                          zm_const%latvap/zm_const%cpair*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
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
                         dp(i,k) - zm_const%latvap/zm_const%cpair*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+zm_const%cpair/zm_const%latvap* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+zm_const%cpair/zm_const%latvap*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+zm_const%cpair/zm_const%latvap*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+zm_const%cpair/zm_const%latvap*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p_mid(i,k))** (zm_const%rdair/zm_const%cpair)*(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t_mid(i,k)* (1000._r8/p_mid(i,k))** (zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q(i,k))
            dqsdtp(i,k) = qstp(i,k)* (1._r8+qstp(i,k)/zm_const%epsilo)*zm_const%epsilo*zm_const%latvap/(zm_const%rdair*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(i,k) = tp(i,k)/ (1._r8+zm_const%latvap/zm_const%cpair* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                        (dtbdt(i)/t_mid(i,mx(i))+zm_const%latvap/zm_const%cpair* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                         tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._r8/(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))* &
                          (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t_mid(i,k)+0.608_r8/ &
                          (1._r8+0.608_r8*q(i,k))*dqmdt(i,k)))*zm_const%grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)   * (1000._r8/p_mid(i,k))** (zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q(i,mx(i)))
            thetavm(i,k) = t_mid(i,k)* (1000._r8/p_mid(i,k))** (zm_const%rdair/zm_const%cpair)*(1._r8+0.608_r8*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (dtbdt(i)/t_mid(i,mx(i))+0.608_r8/ (1._r8+0.608_r8*q(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t_mid(i,k)-0.608_r8/ (1._r8+0.608_r8*q(i,k))*dqmdt(i,k))* &
                          zm_const%grav*thetavp(i,k)/thetavm(i,k)
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
            dadt(i) = dadt(i) + dboydt(i,k)* (z_int(i,k)-z_int(i,k+1))
         endif
      end do
   end do
   do i = il1g,il2g
      dltaa = -1._r8* (cape(i)-capelmt)
      if (dadt(i) /= 0._r8) mb(i) = max(dltaa/zm_param%tau/dadt(i),0._r8)
      if (zm_param%zm_microp .and. mx(i)-jt(i) < 2._r8) mb(i) =0.0_r8
   end do
!
   return
end subroutine closure

!===================================================================================================

subroutine zm_calculate_output_tend(pcols, ncol, pver, pverp, &
                                    dqdt, dsdt, q, qs, qu, &
                                    su, du, qhat, shat, dp, &
                                    mu, md, sd, qd, ql, &
                                    dsubcld, jt, mx, il1g, il2g, msg, &
                                    dl, evp, cu, &
                                    loc_microp_st)
   !----------------------------------------------------------------------------
   ! Purpose: initialize quantities for ZM convection scheme
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols        ! maximum number of columns
   integer,                         intent(in   ) :: ncol         ! actual number of columns
   integer,                         intent(in   ) :: pver         ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp        ! number of interface vertical levels
   integer,                         intent(in   ) :: il1g         !
   integer,                         intent(in   ) :: il2g         !
   integer,                         intent(in   ) :: msg          !
   real(r8), dimension(pcols,pver), intent(in   ) :: q            !
   real(r8), dimension(pcols,pver), intent(in   ) :: qs           !
   real(r8), dimension(pcols,pver), intent(in   ) :: qu           !
   real(r8), dimension(pcols,pver), intent(in   ) :: su           !
   real(r8), dimension(pcols,pver), intent(in   ) :: du           !
   real(r8), dimension(pcols,pver), intent(in   ) :: qhat         !
   real(r8), dimension(pcols,pver), intent(in   ) :: shat         !
   real(r8), dimension(pcols,pver), intent(in   ) :: dp           !
   real(r8), dimension(pcols,pver), intent(in   ) :: mu           !
   real(r8), dimension(pcols,pver), intent(in   ) :: md           !
   real(r8), dimension(pcols,pver), intent(in   ) :: sd           !
   real(r8), dimension(pcols,pver), intent(in   ) :: qd           !
   real(r8), dimension(pcols,pver), intent(in   ) :: ql           !
   real(r8), dimension(pcols,pver), intent(in   ) :: evp          !
   real(r8), dimension(pcols,pver), intent(in   ) :: cu           !
   real(r8), dimension(pcols),      intent(in   ) :: dsubcld      !
   real(r8), dimension(pcols,pver), intent(  out) :: dqdt         !
   real(r8), dimension(pcols,pver), intent(  out) :: dsdt         !
   real(r8), dimension(pcols,pver), intent(  out) :: dl           !
   type(zm_microp_st),              intent(inout) :: loc_microp_st! convective microphysics state and tendencies
   !----------------------------------------------------------------------------
   ! Local variables
   integer i,k
   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
   real(r8) emc
   !----------------------------------------------------------------------------
   ! initialize variables
   do k = msg + 1,pver
      do i = il1g,il2g
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
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   !----------------------------------------------------------------------------
   ! calculate large-scale tendencies
   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu(i,k) + evp(i,k) ! condensation in updraft and evaporating rain in downdraft

         dsdt(i,k) = -zm_const%latvap/zm_const%cpair*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         if (zm_param%zm_microp) then
            dsdt(i,k) = dsdt(i,k) + zm_const%latice/zm_const%cpair*loc_microp_st%frz(i,k)
            dl                (i,k) = du(i,k)*loc_microp_st%qcde(i,k+1)
            loc_microp_st%dif (i,k) = du(i,k)*loc_microp_st%qide(i,k+1)
            loc_microp_st%dnlf(i,k) = du(i,k)*loc_microp_st%ncde(i,k+1)
            loc_microp_st%dnif(i,k) = du(i,k)*loc_microp_st%nide(i,k+1)
            loc_microp_st%dsf (i,k) = du(i,k)*loc_microp_st%qsde(i,k+1)
            loc_microp_st%dnsf(i,k) = du(i,k)*loc_microp_st%nsde(i,k+1)
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
   !----------------------------------------------------------------------------
   return
end subroutine zm_calculate_output_tend

!===================================================================================================

end module zm_conv
