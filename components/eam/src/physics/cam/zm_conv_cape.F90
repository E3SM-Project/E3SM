module zm_conv_cape
   !----------------------------------------------------------------------------
   ! Purpose: CAPE calculation methods for ZM deep convection scheme
   !----------------------------------------------------------------------------
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use cam_abortutils,   only: endrun
   use zm_conv_util,     only: entropy, ientropy, qsat_hPa
   use zm_conv_types,    only: zm_const_t, zm_param_t

   implicit none
   private

   public :: compute_dilute_cape ! calculate convective available potential energy (CAPE) with dilute parcel ascent

   real(r8), parameter :: lcl_pressure_threshold     = 600._r8   ! if LCL pressure is lower => no convection and cape is zero
   real(r8), parameter :: ull_upper_launch_pressure  = 600._r8   ! upper search limit for unrestricted launch level (ULL)
   real(r8), parameter :: pergro_rhd_threshold       = -1.e-4_r8 ! MSE difference threshold for perturbation growth test
!===================================================================================================
contains
!===================================================================================================

subroutine compute_dilute_cape( pcols, ncol, pver, pverp, &
                                num_cin, num_msg, &
                                sp_humidity_in, temperature_in, &
                                zmid, pmid, pint, pblt, tpert, &
                                parcel_temp, parcel_qsat, msemax_klev, &
                                lcl_temperature, lcl_klev, &
                                eql_klev, cape, &
                                zm_const, zm_param, &
                                iclosure, dcapemx, &
                                use_input_tq_mx, q_mx, t_mx )
   !----------------------------------------------------------------------------
   ! Purpose: calculate convective available potential energy (CAPE), lifting
   !          condensation level (LCL), and convective top with dilute parcel ascent
   ! Method: parcel temperature based on a plume model with constant entrainment
   ! Original Author: Richard Neale - September 2004
   ! References:
   !   Raymond, D. J., and A. M. Blyth, 1986: A Stochastic Mixing Model for
   !     Nonprecipitating Cumulus Clouds. J. Atmos. Sci., 43, 2708–2718
   !   Raymond, D. J., and A. M. Blyth, 1992: Extension of the Stochastic Mixing
   !     Model to Cumulonimbus Clouds. J. Atmos. Sci., 49, 1968–1983
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                              intent(in   ) :: pcols               ! number of atmospheric columns (max)
   integer,                              intent(in   ) :: ncol                ! number of atmospheric columns (actual)
   integer,                              intent(in   ) :: pver                ! number of mid-point vertical levels
   integer,                              intent(in   ) :: pverp               ! number of interface vertical levels
   integer,                              intent(in   ) :: num_cin             ! num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
   integer,                              intent(in   ) :: num_msg             ! index of highest level convection is allowed
   real(r8), dimension(pcols,pver),      intent(in   ) :: sp_humidity_in      ! specific humidity
   real(r8), dimension(pcols,pver),      intent(in   ) :: temperature_in      ! temperature
   real(r8), dimension(pcols,pver),      intent(in   ) :: zmid                ! altitude/height at mid-levels
   real(r8), dimension(pcols,pver),      intent(in   ) :: pmid                ! pressure at mid-levels
   real(r8), dimension(pcols,pverp),     intent(in   ) :: pint                ! pressure at interfaces
   integer,  dimension(pcols),           intent(in   ) :: pblt                ! index of pbl top used as upper limit index of max MSE search
   real(r8), dimension(pcols),           intent(in   ) :: tpert               ! perturbation temperature by pbl processes
   real(r8), dimension(pcols,pver),      intent(  out) :: parcel_temp         ! parcel temperature
   real(r8), dimension(pcols,pver),      intent(inout) :: parcel_qsat         ! parcel saturation mixing ratio
   integer,  dimension(pcols),           intent(inout) :: msemax_klev         ! index of max MSE at parcel launch level
   real(r8), dimension(pcols),           intent(  out) :: lcl_temperature     ! lifting condensation level (LCL) temperature
   integer,  dimension(pcols),           intent(inout) :: lcl_klev            ! index of lifting condensation level (i.e. cloud bottom)
   integer,  dimension(pcols),           intent(inout) :: eql_klev            ! index of equilibrium level (i.e. cloud top)
   real(r8), dimension(pcols),           intent(inout) :: cape                ! convective available potential energy
   type(zm_const_t),                     intent(in   ) :: zm_const            ! derived type to hold ZM constants
   type(zm_param_t),                     intent(in   ) :: zm_param            ! derived type to hold ZM tunable parameters
   logical,                              intent(in   ) :: iclosure            ! true for normal procedure, otherwise use dcapemx from 1st call
   integer,  dimension(pcols), optional, intent(in   ) :: dcapemx             ! values of msemax_klev from previous call for dcape closure
   logical,                    optional, intent(in   ) :: use_input_tq_mx     ! if .true., use input values of dcapemx, q_mx, t_mx in the CAPE calculation
   real(r8), dimension(pcols), optional, intent(inout) :: q_mx                ! specified sp humidity to apply at level of max MSE if use_input_tq_mx=.true.
   real(r8), dimension(pcols), optional, intent(inout) :: t_mx                ! specified temperature to apply at level of max MSE if use_input_tq_mx=.true.
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8), dimension(pcols,pver)    :: sp_humidity     ! local version of specific humidity
   real(r8), dimension(pcols,pver)    :: temperature     ! local version of temperature
   real(r8), dimension(pcols,pver)    :: tv              ! virtual temperature
   real(r8), dimension(pcols,pver)    :: parcel_vtemp    ! parcel virtual temperature
   real(r8), dimension(pcols)         :: lcl_pmid        ! lifting condensation level (LCL) pressure
   real(r8), dimension(pcols)         :: mse_max_val     ! value of max MSE at parcel launch level
   integer,  dimension(pcols)         :: pblt_ull        ! upper limit index of max MSE search for ULL
   integer,  dimension(pcols)         :: msemax_top_k    ! upper limit index of max MSE search

   integer  :: i, k                       ! loop iterators
   logical  :: pergro_active              ! flag for perturbation growth test (pergro)
   logical  :: use_input_tq_mx_loc    ! flag to use input parcel temperature and specific humidity
   !----------------------------------------------------------------------------
   ! set flag for perturbation growth test
#ifdef PERGRO
   pergro_active = .true.
#else
   pergro_active = .false.
#endif
   !----------------------------------------------------------------------------
   use_input_tq_mx_loc = .false.
   if (present(use_input_tq_mx)) use_input_tq_mx_loc = use_input_tq_mx

   if ( use_input_tq_mx_loc  .and. &
       ((.not.present(t_mx)) .or.  &
        (.not.present(q_mx)) .or.  &
        (.not.present(dcapemx)) )  ) then
      call endrun('compute_dilute_cape: use_input_tq_mx = .true. but dcapemx, t_mx or q_mx is not provided')
   end if

   !----------------------------------------------------------------------------
   ! Copy the incoming temperature and specific humidity values to local arrays 
   temperature(:ncol,:) = temperature_in(:ncol,:)
   sp_humidity(:ncol,:) = sp_humidity_in(:ncol,:)

   !----------------------------------------------------------------------------
   ! initialize msemax_klev and potentially modify T/q
   if (use_input_tq_mx_loc) then
      ! note - in this case we expect:
      ! (1) the incoming array dcapemx contains prev identified launching level index, and 
      ! (2) the arrays q_mx and t_mx contain q and T values at the old launching level 
      !     at the time when the old launching level was identified. 
      ! Copy the old values to work arrays for calculations in the rest of this subroutine
      msemax_klev(:ncol) = dcapemx(:ncol)
      do i=1,ncol
         sp_humidity(i,msemax_klev(i)) = q_mx(i)
         temperature(i,msemax_klev(i)) = t_mx(i)
      end do
   else
      msemax_klev(:) = pver
   end if

   !----------------------------------------------------------------------------
   ! Initialize parcel properties
   parcel_temp(1:ncol,1:pver) = temperature(1:ncol,1:pver)
   parcel_qsat(1:ncol,1:pver) = sp_humidity(1:ncol,1:pver)
   tv  (1:ncol,1:pver) = temperature(1:ncol,1:pver) &
                         * ( 1._r8+zm_const%zvir*sp_humidity(1:ncol,1:pver) ) &
                         / ( 1._r8+sp_humidity(1:ncol,1:pver) )
   parcel_vtemp (1:ncol,1:pver) = tv(1:ncol,1:pver)

   !----------------------------------------------------------------------------
   ! Find new upper bound for parcel starting level - unrestricted launch level (ULL)
   if (zm_param%trig_ull) then
      pblt_ull(:) = 1
      do k = pver-1, num_msg+1, -1
         do i = 1,ncol
            if ( (pmid(i,k)  .le.ull_upper_launch_pressure) .and. &
                 (pmid(i,k+1).gt.ull_upper_launch_pressure) ) then
               pblt_ull(i) = k
            end if
         end do
      end do
   endif

   !----------------------------------------------------------------------------
   ! Set level of max moist static energy for parcel initialization
   if ( zm_param%trig_dcape .and. (.not.iclosure) ) then
      ! Use max moist static energy level that is passed in
      if (.not.present(dcapemx)) call endrun('** ZM CONV compute_dilute_cape: dcapemx not present **')
      msemax_klev(1:ncol) = dcapemx(1:ncol)
   elseif (.not.use_input_tq_mx_loc) then
      if (     zm_param%trig_ull) msemax_top_k(:ncol) = pblt_ull(:ncol)
      if (.not.zm_param%trig_ull) msemax_top_k(:ncol) = pblt(:ncol)
      call find_mse_max( pcols, ncol, pver, num_msg, msemax_top_k, pergro_active, &
                         temperature, zmid, sp_humidity, zm_const, zm_param, &
                         msemax_klev, mse_max_val )
   end if

   !----------------------------------------------------------------------------
   do i=1,ncol
      ! Save launching level T, q for output
      if ( .not.use_input_tq_mx_loc .and. present(q_mx) .and. present(t_mx) ) then
         q_mx(i) = sp_humidity(i,msemax_klev(i))
         t_mx(i) = temperature(i,msemax_klev(i))
      end if
      ! save LCL values for compute_dilute_parcel()
      lcl_klev(i) = msemax_klev(i)
      lcl_pmid(i) = pmid(i,msemax_klev(i))
      lcl_temperature(i) = temperature(i,msemax_klev(i))
   end do

   !----------------------------------------------------------------------------
   ! entraining parcel calculation
   call compute_dilute_parcel( pcols, ncol, pver, num_msg, msemax_klev, &
                               pmid, temperature, sp_humidity, tpert, pblt, &
                               zm_const, zm_param, &
                               parcel_temp, parcel_vtemp, parcel_qsat, &
                               lcl_pmid, lcl_temperature, lcl_klev )

   !----------------------------------------------------------------------------
   ! calculate CAPE
   call compute_cape_from_parcel( pcols, ncol, pver, pverp, num_cin, num_msg, &
                                  temperature, tv, zmid, sp_humidity, pint, &
                                  msemax_klev, lcl_pmid, lcl_klev, &
                                  zm_const, zm_param, &
                                  parcel_qsat, parcel_temp, parcel_vtemp, &
                                  eql_klev, cape )

   !----------------------------------------------------------------------------
   return

end subroutine compute_dilute_cape

!===================================================================================================

subroutine find_mse_max( pcols, ncol, pver, num_msg, msemax_top_k, pergro_active, &
                         temperature, zmid, sp_humidity, zm_const, zm_param, &
                         msemax_klev, mse_max_val)
   !----------------------------------------------------------------------------
   ! Purpose: find level of max moist static energy for parcel initialization
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols           ! number of atmospheric columns (max)
   integer,                         intent(in   ) :: ncol            ! number of atmospheric columns (actual)
   integer,                         intent(in   ) :: pver            ! number of mid-point vertical levels
   integer,                         intent(in   ) :: num_msg         ! number of missing moisture levels at the top of model
   integer,  dimension(pcols),      intent(in   ) :: msemax_top_k    ! upper limit index of max MSE search
   logical,                         intent(in   ) :: pergro_active   ! flag for perturbation growth test (pergro)
   real(r8), dimension(pcols,pver), intent(in   ) :: temperature     ! environement temperature
   real(r8), dimension(pcols,pver), intent(in   ) :: zmid            ! height/altitude at mid-levels
   real(r8), dimension(pcols,pver), intent(in   ) :: sp_humidity     ! specific humidity
   type(zm_const_t),                intent(in   ) :: zm_const        ! derived type to hold ZM constants
   type(zm_param_t),                intent(in   ) :: zm_param        ! derived type to hold ZM tunable parameters
   integer,  dimension(pcols),      intent(inout) :: msemax_klev     ! index of max MSE at parcel launch level
   real(r8), dimension(pcols),      intent(inout) :: mse_max_val     ! value of max MSE at parcel launch level
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i,k                       ! loop iterators
   integer  :: bot_layer                 ! lower limit to search for parcel launch level
   real(r8) :: pergro_rhd                ! relative MSE (h) difference for perturbation growth test (pergro)
   real(r8), dimension(pcols) :: mse_env ! env moist static energy
   !----------------------------------------------------------------------------
   ! initialize values
   mse_max_val(1:ncol) = 0._r8
   bot_layer = pver - zm_param%mx_bot_lyr_adj ! set lower limit to search for launch level with max MSE
   !----------------------------------------------------------------------------
   do k = bot_layer, num_msg+1, -1
      do i = 1,ncol
         ! calculate moist static energy
         mse_env(i) = zm_const%cpair*temperature(i,k) + zm_const%grav*zmid(i,k) + zm_const%latvap*sp_humidity(i,k)
         if (pergro_active) then
            ! Reset max moist static energy level when relative difference exceeds 1.e-4
            pergro_rhd = (mse_env(i) - mse_max_val(i))/(mse_env(i) + mse_max_val(i))
            if (k >= msemax_top_k(i) .and. pergro_rhd > pergro_rhd_threshold) then
               mse_max_val(i) = mse_env(i)
               msemax_klev(i) = k
            end if
         else
            ! find level and value of max moist static energy
            if (k >= msemax_top_k(i) .and. mse_env(i) > mse_max_val(i)) then
               mse_max_val(i) = mse_env(i)
               msemax_klev(i) = k
            end if
         end if
      end do
   end do
   !----------------------------------------------------------------------------
end subroutine find_mse_max

!===================================================================================================

subroutine compute_dilute_parcel( pcols, ncol, pver, num_msg, klaunch, &
                                  pmid, temperature, sp_humidity, tpert, pblt, &
                                  zm_const, zm_param, &
                                  parcel_temp, parcel_vtemp, parcel_qsat, &
                                  lcl_pmid, lcl_temperature, lcl_klev )
   !----------------------------------------------------------------------------
   ! Purpose: Calculate thermodynamic properties of an entraining air parcel 
   !          lifted from the PBL using fractional mass entrainment rate 
   !          specified by zm_param%dmpdz
   !----------------------------------------------------------------------------
   implicit none
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                         intent(in   ) :: pcols           ! number of atmospheric columns (max)
   integer,                         intent(in   ) :: ncol            ! number of atmospheric columns (actual)
   integer,                         intent(in   ) :: pver            ! number of mid-point vertical levels
   integer,                         intent(in   ) :: num_msg         ! number of missing moisture levels at the top of model
   integer,  dimension(pcols),      intent(in   ) :: klaunch         ! index of parcel launch level based on max MSE
   real(r8), dimension(pcols,pver), intent(in   ) :: pmid            ! ambient env pressure at cell center
   real(r8), dimension(pcols,pver), intent(in   ) :: temperature     ! ambient env temperature at cell center
   real(r8), dimension(pcols,pver), intent(in   ) :: sp_humidity     ! ambient env specific humidity at cell center
   real(r8), dimension(pcols),      intent(in   ) :: tpert           ! PBL temperature perturbation
   integer,  dimension(pcols),      intent(in   ) :: pblt            ! index of pbl depth 
   type(zm_const_t),                intent(in   ) :: zm_const        ! derived type to hold ZM constants
   type(zm_param_t),                intent(in   ) :: zm_param        ! derived type to hold ZM tunable parameters
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_temp     ! Parcel temperature
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_vtemp    ! Parcel virtual temperature
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_qsat     ! Parcel water vapour (sat value above lcl)
   real(r8), dimension(pcols)     , intent(inout) :: lcl_pmid        ! lifting condensation level (LCL) pressure
   real(r8), dimension(pcols)     , intent(inout) :: lcl_temperature ! lifting condensation level (LCL) temperature
   integer,  dimension(pcols)     , intent(inout) :: lcl_klev        ! lifting condensation level (LCL) vertical index
   !----------------------------------------------------------------------------
   ! Local variables
   integer i,k,ii   ! loop iterators

   real(r8), dimension(pcols,pver) :: tmix        ! tempertaure of the entraining parcel.
   real(r8), dimension(pcols,pver) :: qtmix       ! total water of the entraining parcel.
   real(r8), dimension(pcols,pver) :: qsmix       ! saturated mixing ratio at the tmix.
   real(r8), dimension(pcols,pver) :: smix        ! entropy of the entraining parcel.
   real(r8), dimension(pcols,pver) :: xsh2o       ! precipitate lost from parcel.
   real(r8), dimension(pcols,pver) :: ds_xsh2o    ! entropy change due to loss of condensate.
   real(r8), dimension(pcols,pver) :: ds_freeze   ! entropy change sue to freezing of precip.

   real(r8), dimension(pcols) :: mp               ! parcel mass flux
   real(r8), dimension(pcols) :: qtp              ! parcel total water
   real(r8), dimension(pcols) :: sp               ! parcel entropy
   real(r8), dimension(pcols) :: sp0              ! parcel launch entropy
   real(r8), dimension(pcols) :: qtp0             ! parcel launch total water
   real(r8), dimension(pcols) :: mp0              ! parcel launch relative mass [0-1]

   real(r8), dimension(pcols) :: tpert_loc           ! gather parcel temperature perturbation

   real(r8) dmpdp    ! parcel fractional mass entrainment rate [1/mb]
   real(r8) dpdz     ! hydrstatic relation
   real(r8) dzdp     ! inverse hydrstatic relation
   real(r8) senv     ! environmental entropy
   real(r8) qtenv    ! environmental total water
   real(r8) penv     ! environmental total pressure
   real(r8) tenv     ! environmental total temperature
   real(r8) new_s    ! hold value for entropy after condensation/freezing adjustments
   real(r8) new_q    ! hold value for total water after condensation/freezing adjustments
   real(r8) dp       ! layer thickness (center to center)
   real(r8) tfguess  ! first guess for entropy inversion
   real(r8) tscool   ! super cooled temperature offset from freezing temperature when cloud water loading freezes

   real(r8) qxsk     ! LCL excess water @ k
   real(r8) qxskp1   ! LCL excess water @ k+1
   real(r8) dsdp     ! LCL entropy gradient @ k
   real(r8) dqtdp    ! LCL total water gradient @ k
   real(r8) dqxsdp   ! LCL excess water gradient @ k+1
   real(r8) slcl     ! LCL entropy
   real(r8) qtlcl    ! LCL total water
   real(r8) qslcl    ! LCL saturated vapor mixing ratio

   integer rcall     ! ientropy call id for error message
   
   integer,  parameter :: nit_lheat = 2      ! Number of iterations for condensation/freezing loop
   real(r8), parameter :: lwmax = 1.e-3_r8   ! maximum condesate that can be held in cloud before rainout

   !----------------------------------------------------------------------------
   ! initialize values
   tscool      = 0._r8
   qtmix       = 0._r8
   smix        = 0._r8
   qtenv       = 0._r8
   senv        = 0._r8
   tenv        = 0._r8
   penv        = 0._r8
   qtp0        = 0._r8
   sp0         = 0._r8
   mp0         = 0._r8
   qtp         = 0._r8
   sp          = 0._r8
   mp          = 0._r8
   new_q       = 0._r8
   new_s       = 0._r8
   !----------------------------------------------------------------------------
   ! The original ZM scheme only treated PBL-rooted convection. A PBL temperature 
   ! perturbation (tpert) was then used to increase the parcel temperatue at launch
   ! level, which is in PBL. The dcape_ull or ull triggr enables ZM scheme to treat
   ! elevated convection with launch level above PBL. If parcel launch level is
   ! above PBL top, tempeature perturbation in PBL should not be able to influence 
   ! it. In this situation, the temporary variable tpert_loc is reset to zero.  
   do i=1,ncol
      tpert_loc(i) = tpert(i)
      if ( zm_param%tpert_fix .and. klaunch(i)<pblt(i) ) then
         tpert_loc(i) = 0._r8
      end if
   end do

   !----------------------------------------------------------------------------
   ! entrainment loop
   do k = pver, num_msg+1, -1
      do i = 1,ncol

         if ( k == klaunch(i) ) then 

            ! initialize values at launch level
            mp0(i)      = 1._r8            ! initial relative mass - value of 1.0 does not change for undilute (dmpdp=0)
            qtp0(i)     = sp_humidity(i,k) ! initial total water - assuming subsaturated
            sp0(i)      = entropy( temperature(i,k), pmid(i,k), qtp0(i), zm_const )
            smix(i,k)   = sp0(i)
            qtmix(i,k)  = qtp0(i)
            tfguess     = temperature(i,k)
            rcall = 1
            call ientropy( rcall, smix(i,k), pmid(i,k), qtmix(i,k), tmix(i,k), qsmix(i,k), tfguess, zm_const )
         
         elseif ( k < klaunch(i) ) then 

            ! set environmental values for this level
            dp    = pmid(i,k) - pmid(i,k+1)
            qtenv = 0.5_r8*(sp_humidity(i,k)+sp_humidity(i,k+1))
            tenv  = 0.5_r8*(temperature(i,k)+temperature(i,k+1)) 
            penv  = 0.5_r8*(pmid(i,k)+pmid(i,k+1))
            senv  = entropy( tenv, penv, qtenv, zm_const )

            ! determine fractional entrainment rate 1/pa given value 1/m
            dpdz = -(penv*zm_const%grav)/(zm_const%rdair*tenv) ! [mb/m] since pressure should be [mb]
            dzdp = 1._r8/dpdz                  ! [m/mb]
            dmpdp = zm_param%dmpdz*dzdp   ! Fractional entrainment [1/mb]

            ! sum entrainment to current level - entrain q,s out of intervening dp layers, 
            ! assuming linear variation (i.e. entrain the mean of the 2 stored values)
            sp(i)  = sp(i)  - dmpdp*dp*senv 
            qtp(i) = qtp(i) - dmpdp*dp*qtenv 
            mp(i)  = mp(i)  - dmpdp*dp
            
            ! entrain s and qt to next level
            smix(i,k)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
            qtmix(i,k) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

            ! Invert entropy from s and q to determine T and saturation-capped q of mixture
            ! temperature(i,k) used as a first guess so that it converges faster

            tfguess = tmix(i,k+1)
            rcall = 2
            call ientropy( rcall, smix(i,k), pmid(i,k), qtmix(i,k), tmix(i,k), qsmix(i,k), tfguess, zm_const )

            ! determine if we are at the LCL if this is first level where qsmix<=qtmix on ascending
            if ( qsmix(i,k)<=qtmix(i,k) .and. qsmix(i,k+1)>qtmix(i,k+1) ) then
               lcl_klev(i) = k
               qxsk        = qtmix(i,k) - qsmix(i,k)
               qxskp1      = qtmix(i,k+1) - qsmix(i,k+1)
               dqxsdp      = (qxsk - qxskp1)/dp
               lcl_pmid(i) = pmid(i,k+1) - qxskp1/dqxsdp
               dsdp        = (smix(i,k)  - smix(i,k+1))/dp
               dqtdp       = (qtmix(i,k) - qtmix(i,k+1))/dp
               slcl        = smix(i,k+1)  + dsdp* (lcl_pmid(i)-pmid(i,k+1))  
               qtlcl       = qtmix(i,k+1) + dqtdp*(lcl_pmid(i)-pmid(i,k+1))
               tfguess     = tmix(i,k)
               rcall = 3
               call ientropy( rcall, slcl, lcl_pmid(i), qtlcl, lcl_temperature(i), qslcl, tfguess, zm_const )
            endif

         end if !  k < klaunch

      end do ! i = 1,ncol
   end do ! k = pver, num_msg+1, -1
   !----------------------------------------------------------------------------
   ! end of entrainment loop
   
   !----------------------------------------------------------------------------
   ! We now have a profile of entropy and total water of the entraining parcel
   ! Varying with height from the launch level klaunch parcel=environment. To the
   ! top allowed level for the existence of convection. If we stop now it will 
   ! provide some estimate of buoyancy without the effects of freezing/condensation.
   ! 
   ! Instead, we will adjust these values such that the water held in vapor is
   ! <=qsmix. We assume that the cloud holds a certain amount of condensate (lwmax)
   ! and the rest is rained out (xsh2o). This provides latent heating to the
   ! mixed parcel and so this has to be added back to it.

   !----------------------------------------------------------------------------
   ! precipitation/freezing loop - iterate twice for accuracy
   xsh2o = 0._r8
   ds_xsh2o = 0._r8
   ds_freeze = 0._r8
   do k = pver, num_msg+1, -1
      do i = 1,ncol      

         if ( k == klaunch(i) ) then

            ! initialize values at launch level - assume no liquid water
            parcel_temp(i,k)  = tmix(i,k)
            parcel_qsat(i,k)  = sp_humidity(i,k) 
            parcel_vtemp(i,k) = ( parcel_temp(i,k) +zm_param%tpert_fac*tpert_loc(i) ) & 
                                * (1._r8+zm_const%zvir*parcel_qsat(i,k)) / (1._r8+parcel_qsat(i,k))

         elseif ( k < klaunch(i) ) then

            ! iterate nit_lheat times for s,qt changes
            do ii = 0,nit_lheat-1

               ! rain (xsh2o) is excess condensate, bar lwmax (accumulated loss from qtmix)
               xsh2o(i,k) = max (0._r8, qtmix(i,k) - qsmix(i,k) - lwmax)

               ! contribution to ds from precip loss of condensate (accumulated change from smix)
               ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - zm_const%cpliq * log (tmix(i,k)/zm_const%tfreez) * max(0._r8,(xsh2o(i,k)-xsh2o(i,k+1)))

               ! calculate entropy of freezing => ( latice x amount of water involved ) / T

               ! one off freezing of condensate
               if (tmix(i,k) <= (zm_const%tfreez+tscool) .and. ds_freeze(i,k+1) == 0._r8) then 
                  ! entropy change from latent heat
                  ds_freeze(i,k) = (zm_const%latice/tmix(i,k)) * max(0._r8,qtmix(i,k)-qsmix(i,k)-xsh2o(i,k)) 
               end if
            
               if (tmix(i,k) <= zm_const%tfreez+tscool .and. ds_freeze(i,k+1) /= 0._r8) then 
                  ! continual freezing of additional condensate
                  ds_freeze(i,k) = ds_freeze(i,k+1)+(zm_const%latice/tmix(i,k)) * max(0._r8,(qsmix(i,k+1)-qsmix(i,k)))
               end if
            
               ! adjust entropy and accordingly to sum of ds (be careful of signs)
               new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k) 

               ! adjust liquid water and accordingly to xsh2o
               new_q = qtmix(i,k) - xsh2o(i,k)

               ! invert entropy to get updated Tmix and qsmix of parcel
               tfguess = tmix(i,k)
               rcall =4
               call ientropy( rcall, new_s, pmid(i,k), new_q, tmix(i,k), qsmix(i,k), tfguess, zm_const )
               
            end do  ! iteration loop for freezing processes

            ! parcel temp is temp of mixture
            ! parcel virtual temp should be density temp with new_q total water
            parcel_temp(i,k) = tmix(i,k)

            ! parcel_vtemp=tprho in the presence of condensate (i.e. when new_q > qsmix)
            if (new_q > qsmix(i,k)) then  ! super-saturated so condensate present - reduces buoyancy
               parcel_qsat(i,k) = qsmix(i,k)
            else                          ! just saturated/sub-saturated - no condensate virtual effects
               parcel_qsat(i,k) = new_q
            end if

            parcel_vtemp(i,k) = ( parcel_temp(i,k) +zm_param%tpert_fac*tpert_loc(i) ) &
                                * (1._r8+zm_const%zvir*parcel_qsat(i,k)) / (1._r8+ new_q)

         end if ! k < klaunch
      
      end do ! i = 1,ncol
   end do ! k = pver, num_msg+1, -1

   !----------------------------------------------------------------------------
   return

end subroutine compute_dilute_parcel

!===================================================================================================

subroutine compute_cape_from_parcel( pcols, ncol, pver, pverp, num_cin, num_msg, &
                                     temperature, tv, zmid, sp_humidity, pint, &
                                     msemax_klev, lcl_pmid, lcl_klev, &
                                     zm_const, zm_param, &
                                     parcel_qsat, parcel_temp, parcel_vtemp, &
                                     eql_klev, cape )
   !----------------------------------------------------------------------------
   ! Purpose: calculate convective available potential energy (CAPE)
   !          from parcel thermodynamic properties from compute_dilute_parcel()
   !----------------------------------------------------------------------------
   integer,                         intent(in   ) :: pcols        ! number of atmospheric columns (max)
   integer,                         intent(in   ) :: ncol         ! number of atmospheric columns (actual)
   integer,                         intent(in   ) :: pver         ! number of mid-point vertical levels
   integer,                         intent(in   ) :: pverp        ! number of interface vertical levels
   integer,                         intent(in   ) :: num_cin      ! num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
   integer,                         intent(in   ) :: num_msg      ! number of missing moisture levels at the top of model
   real(r8), dimension(pcols,pver), intent(in   ) :: temperature  ! temperature
   real(r8), dimension(pcols,pver), intent(in   ) :: tv           ! virtual temperature
   real(r8), dimension(pcols,pver), intent(in   ) :: zmid         ! height/altitude at mid-levels
   real(r8), dimension(pcols,pver), intent(in   ) :: sp_humidity  ! specific humidity
   real(r8), dimension(pcols,pverp),intent(in   ) :: pint         ! pressure at interfaces
   integer,  dimension(pcols),      intent(in   ) :: msemax_klev  ! index of max MSE at parcel launch level
   real(r8), dimension(pcols),      intent(in   ) :: lcl_pmid     ! lifting condensation level (LCL) pressure
   integer,  dimension(pcols),      intent(in   ) :: lcl_klev     ! lifting condensation level (LCL) index
   type(zm_const_t),                intent(in   ) :: zm_const     ! derived type to hold ZM constants
   type(zm_param_t),                intent(in   ) :: zm_param     ! derived type to hold ZM tunable parameters
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_qsat  ! parcel saturation mixing ratio
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_temp  ! parcel temperature
   real(r8), dimension(pcols,pver), intent(inout) :: parcel_vtemp ! parcel virtual temperature
   integer,  dimension(pcols),      intent(inout) :: eql_klev     ! index of equilibrium level (i.e. cloud top)
   real(r8), dimension(pcols),      intent(inout) :: cape         ! convective available potential energy
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, k, n                                      ! loop iterators
   real(r8), dimension(pcols,pver)    :: buoyancy           ! parcel buoyancy
   real(r8), dimension(pcols,num_cin) :: cape_tmp           ! provisional value of cape
   integer,  dimension(pcols,num_cin) :: eql_klev_tmp       ! provisional value of equilibrium level index
   integer,  dimension(pcols)         :: neg_buoyancy_cnt   ! counter for levels with negative bounancy

   logical plge600(pcols) ! for testing - remove!
   !----------------------------------------------------------------------------
   ! Initialize variables
   eql_klev     (1:ncol)               = pver
   eql_klev_tmp (1:ncol,1:num_cin)     = pver
   cape            (1:ncol)            = 0._r8
   cape_tmp        (1:ncol,1:num_cin)  = 0._r8
   buoyancy        (1:ncol,1:pver)     = 0._r8
   neg_buoyancy_cnt(1:ncol)            = 0
   !----------------------------------------------------------------------------
   ! Calculate buoyancy
   do k = pver, num_msg+1, -1
      do i=1,ncol
         ! Define buoyancy from launch level to equilibrium level
         if ( k <= msemax_klev(i) .and. lcl_pmid(i).ge.lcl_pressure_threshold ) then
            buoyancy(i,k) = parcel_vtemp(i,k) - tv(i,k) + zm_param%tiedke_add
         else
            parcel_qsat(i,k)  = sp_humidity(i,k)
            parcel_temp(i,k)  = temperature(i,k)            
            parcel_vtemp(i,k) = tv(i,k)
         endif
      end do
   end do

   !----------------------------------------------------------------------------
   ! find convective equilibrium level accounting for negative buoyancy levels
   do k = num_msg+2, pver
      do i = 1,ncol
         if ( k < lcl_klev(i) .and. lcl_pmid(i).ge.lcl_pressure_threshold ) then
            if ( buoyancy(i,k+1) > 0._r8 .and. &
                 buoyancy(i,k)   <=0._r8 ) then
               neg_buoyancy_cnt(i) = min( num_cin, neg_buoyancy_cnt(i)+1 )
               eql_klev_tmp(i,neg_buoyancy_cnt(i)) = k
            end if
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! integrate buoyancy to obtain possible CAPE values
   do n = 1,num_cin
      do k = num_msg+1, pver
         do i = 1,ncol
            if ( lcl_pmid(i).ge.lcl_pressure_threshold .and. &
                 k <= msemax_klev(i) .and. k > eql_klev_tmp(i,n)) then
               cape_tmp(i,n) = cape_tmp(i,n) + zm_const%rdair*buoyancy(i,k)*log(pint(i,k+1)/pint(i,k))
            end if
         end do
      end do
   end do

   !----------------------------------------------------------------------------
   ! find maximum cape from all possible tentative CAPE values
   ! and use it as the final cape (April 26, 1995)
   do n = 1,num_cin
      do i = 1,ncol
         if (cape_tmp(i,n) > cape(i)) then
            cape(i) = cape_tmp(i,n)
            eql_klev(i) = eql_klev_tmp(i,n)
         end if
      end do
   end do

   !----------------------------------------------------------------------------
   ! apply limiter to ensure CAPE is positive
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
   
   !----------------------------------------------------------------------------
   return

end subroutine compute_cape_from_parcel

!===================================================================================================

end module zm_conv_cape
