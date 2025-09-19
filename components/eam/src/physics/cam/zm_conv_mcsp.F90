module zm_conv_mcsp
   !----------------------------------------------------------------------------
   ! Purpose: methods for mesoscale coherent structure parameterization (MCSP)
   !          This scheme essentially redistributes the ZM heating and drying
   !          tendencies vertically in order to mimic the effects of mesoscale
   !          organization. It can also add momentum tendencies, although this
   !          capability has not been extensively tested
   !----------------------------------------------------------------------------
   ! References:
   ! 
   !   Moncrieff, M. W., & Liu, C. (2006). Representing convective organization in
   !     prediction models by a hybrid strategy. J. Atmos. Sci., 63, 3404–3420.
   !     https://doi.org/10.1175/JAS3812.1
   !
   !   Chen, C.-C., Richter, J. H., Liu, C., Moncrieff, M. W., Tang, Q., Lin, W.,
   !     et al. (2021). Effects of organized convection parameterization on the MJO
   !     and precipitation in E3SMv1. Part I: Mesoscale heating. J. Adv.
   !     Mod. Earth Sys., 13, e2020MS002401, https://doi.org/10.1029/2020MS002401
   !
   !   Moncrieff, M. W., C. Liu, and P. Bogenschutz, 2017: Simulation, Modeling, 
   !     and Dynamically Based Parameterization of Organized Tropical Convection 
   !     for Global Climate Models. J. Atmos. Sci., 74, 1363–1380, 
   !     https://doi.org/10.1175/JAS-D-16-0166.1.
   !
   !   Moncrieff, M. W. (2019). Toward a Dynamical Foundation for Organized Convection
   !     Parameterization in GCMs. Geophys. Res. Lett., 46, 14103–14108.
   !     https://doi.org/10.1029/2019GL085316
   !
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_params, only: r8, pcols, pver, pverp
#else
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use cam_abortutils,   only: endrun
   use cam_logfile,      only: iulog
#endif
   use zm_conv_types,    only: zm_const_t, zm_param_t

   implicit none
   private

#ifndef SCREAM_CONFIG_IS_CMAKE
   public :: zm_conv_mcsp_init ! Initialize MCSP output fields
#endif
   public :: zm_conv_mcsp_tend ! Perform MCSP tendency calculations

   real(r8), parameter :: MCSP_storm_speed_pref = 600e2_r8 ! pressure level for winds in MCSP calculation [Pa]
   real(r8), parameter :: MCSP_conv_depth_min   = 700e2_r8 ! pressure thickness of convective heating [Pa]
   real(r8), parameter :: mcsp_shear_min        = 3.0_r8   ! min shear value for MCSP to be active
   real(r8), parameter :: mcsp_shear_max        = 200.0_r8 ! max shear value for MCSP to be active

!===================================================================================================
contains
!===================================================================================================

! We need to avoid building this for now when bridging from EAMxx
#ifndef SCREAM_CONFIG_IS_CMAKE

subroutine zm_conv_mcsp_init()
   !----------------------------------------------------------------------------
   ! Purpose: initialize MCSP output fields
   !----------------------------------------------------------------------------
   use cam_history,     only: addfld, horiz_only
   use mpishorthand
   !----------------------------------------------------------------------------
   call addfld('MCSP_DT', (/'lev'/), 'A', 'K/s',      'MCSP T tendency')
   call addfld('MCSP_DQ', (/'lev'/), 'A', 'kg/kg/s',  'MCSP qv tendency')
   call addfld('MCSP_DU', (/'lev'/), 'A', 'm/s/day',  'MCSP U wind tendency')
   call addfld('MCSP_DV', (/'lev'/), 'A', 'm/s/day',  'MCSP V wind tendency')

   call addfld('MCSP_freq',     horiz_only, 'A', '1',        'MCSP frequency of activation')
   call addfld('MCSP_shear',    horiz_only, 'A', 'm/s',      'MCSP vertical shear of zonal wind')
   call addfld('MCSP_zm_depth', horiz_only, 'A', 'Pa',  'ZM convection depth for MCSP')

end subroutine zm_conv_mcsp_init

#endif /* SCREAM_CONFIG_IS_CMAKE */

!===================================================================================================

subroutine zm_conv_mcsp_calculate_shear( pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear)
   !----------------------------------------------------------------------------
   ! Purpose: calculate shear for MCSP
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_methods, only: vertinterp
#else
   use interpolate_data, only: vertinterp
#endif
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pmid ! physics state mid-point pressure
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_u    ! physics state u momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_v    ! physics state v momentum
   real(r8), dimension(pcols),            intent(  out) :: mcsp_shear
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i
   real(r8), dimension(pcols) :: storm_u         ! u wind at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_v         ! v wind at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_u_shear   ! u shear at storm reference level set by MCSP_storm_speed_pref
   real(r8), dimension(pcols) :: storm_v_shear   ! v shear at storm reference level set by MCSP_storm_speed_pref
   !----------------------------------------------------------------------------
   ! Interpolate wind to pressure level specified by MCSP_storm_speed_pref
   call vertinterp( ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_u, storm_u )
   call vertinterp( ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_v, storm_v )

   !----------------------------------------------------------------------------
   ! calculate low-level shear
   do i = 1,ncol
      if (state_pmid(i,pver).gt.MCSP_storm_speed_pref) then
         storm_u_shear(i) = storm_u(i)-state_u(i,pver)
         storm_v_shear(i) = storm_v(i)-state_v(i,pver)
      else
         storm_u_shear(i) = -999
         storm_v_shear(i) = -999
      end if
      mcsp_shear(i) = storm_u_shear(i)
   end do

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_calculate_shear

!===================================================================================================

subroutine zm_conv_mcsp_tend( lchnk, pcols, ncol, pver, pverp, &
                              ztodt, jctop, zm_const, zm_param, &
                              state_pmid, state_pint, state_pdel, &
                              state_s, state_q, state_u, state_v, &
                              ptend_zm_s, ptend_zm_q, &
                              ptend_s, ptend_q, ptend_u, ptend_v )
   !----------------------------------------------------------------------------
   ! Purpose: perform MCSP tendency calculations
   !----------------------------------------------------------------------------
#ifndef SCREAM_CONFIG_IS_CMAKE
   use cam_history,      only: outfld
#endif
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                               intent(in   ) :: lchnk      ! chunk identifier
   integer,                               intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                               intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                               intent(in   ) :: pver       ! number of mid-point vertical levels
   integer,                               intent(in   ) :: pverp      ! number of interface vertical levels
   real(r8),                              intent(in   ) :: ztodt      ! 2x physics time step
   integer,  dimension(pcols),            intent(in   ) :: jctop      ! cloud top level indices
   type(zm_const_t),                      intent(in   ) :: zm_const   ! derived type to hold ZM constants
   type(zm_param_t),                      intent(in   ) :: zm_param   ! derived type to hold ZM parameters
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pmid ! physics state mid-point pressure
   real(r8), dimension(pcols,pverp),      intent(in   ) :: state_pint ! physics state interface pressure
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_pdel ! physics state pressure thickness
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_s    ! physics state dry static energy
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_q    ! physics state specific humidity
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_u    ! physics state u momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: state_v    ! physics state v momentum
   real(r8), dimension(pcols,pver),       intent(in   ) :: ptend_zm_s ! input ZM tendency for dry static energy (DSE)
   real(r8), dimension(pcols,pver),       intent(in   ) :: ptend_zm_q ! input ZM tendency for specific humidity (qv)
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_s    ! output tendency of DSE
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_q    ! output tendency of qv
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_u    ! output tendency of u-wind
   real(r8), dimension(pcols,pver),       intent(inout) :: ptend_v    ! output tendency of v-wind
   !----------------------------------------------------------------------------
   ! Local variables
   integer  :: i, k

   real(r8) :: tend_k       ! temporary variable for kinetic energy tendency
   real(r8) :: pdepth_mid_k ! temporary pressure depth used for vertical structure
   real(r8) :: pdepth_total ! temporary pressure depth used for vertical structure

   real(r8), dimension(pcols)      :: zm_avg_tend_s   ! mass weighted column average DSE tendency from ZM
   real(r8), dimension(pcols)      :: zm_avg_tend_q   ! mass weighted column average qv tendency from ZM
   real(r8), dimension(pcols)      :: zm_depth        ! pressure depth of ZM heating
   real(r8), dimension(pcols)      :: mcsp_shear      ! shear used to check against threshold
   real(r8), dimension(pcols)      :: pdel_sum        ! column integrated pressure thickness

   real(r8), dimension(pcols,pver) :: mcsp_tend_s     ! MCSP tendency before energy fixer for DSE
   real(r8), dimension(pcols,pver) :: mcsp_tend_q     ! MCSP tendency before energy fixer for qv
   real(r8), dimension(pcols,pver) :: mcsp_tend_u     ! MCSP tendency before energy fixer for u wind
   real(r8), dimension(pcols,pver) :: mcsp_tend_v     ! MCSP tendency before energy fixer for v wind

   real(r8), dimension(pcols)      :: mcsp_avg_tend_s ! mass weighted column average MCSP tendency of DSE
   real(r8), dimension(pcols)      :: mcsp_avg_tend_q ! mass weighted column average MCSP tendency of qv 
   real(r8), dimension(pcols)      :: mcsp_avg_tend_k ! mass weighted column average MCSP tendency of kinetic energy

   real(r8), dimension(pcols)      :: mcsp_freq       ! MSCP frequency for output
   real(r8), dimension(pcols,pver) :: mcsp_dt_out     ! final MCSP tendency for DSE
   real(r8), dimension(pcols,pver) :: mcsp_dq_out     ! final MCSP tendency for qv
   real(r8), dimension(pcols,pver) :: mcsp_du_out     ! final MCSP tendency for u wind
   real(r8), dimension(pcols,pver) :: mcsp_dv_out     ! final MCSP tendency for v wind

   logical :: do_mcsp_t = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_q = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_u = .false.   ! internal flag to enable tendency calculations
   logical :: do_mcsp_v = .false.   ! internal flag to enable tendency calculations

   !----------------------------------------------------------------------------

   if (.not.zm_param%mcsp_enabled) return

   !----------------------------------------------------------------------------
   ! initialize variables

   if (zm_param%mcsp_t_coeff>0) do_mcsp_t = .true.
   if (zm_param%mcsp_q_coeff>0) do_mcsp_q = .true.
   if (zm_param%mcsp_u_coeff>0) do_mcsp_u = .true.
   if (zm_param%mcsp_v_coeff>0) do_mcsp_v = .true.

   zm_avg_tend_s(1:ncol) = 0
   zm_avg_tend_q(1:ncol) = 0

   pdel_sum(1:ncol) = 0

   mcsp_avg_tend_s(1:ncol) = 0
   mcsp_avg_tend_q(1:ncol) = 0
   mcsp_avg_tend_k(1:ncol) = 0
   
   mcsp_tend_s(1:ncol,1:pver) = 0
   mcsp_tend_q(1:ncol,1:pver) = 0
   mcsp_tend_u(1:ncol,1:pver) = 0
   mcsp_tend_v(1:ncol,1:pver) = 0

   !----------------------------------------------------------------------------
   ! calculate shear

   call zm_conv_mcsp_calculate_shear( pcols, ncol, pver, state_pmid, state_u, state_v, mcsp_shear )

   !----------------------------------------------------------------------------
   ! calculate mass weighted column average tendencies from ZM

   zm_depth(1:ncol) = 0
   do i = 1,ncol
      if (jctop(i) .ne. pver) then
         ! integrate pressure and ZM tendencies over column
         do k = jctop(i),pver
            zm_avg_tend_s(i) = zm_avg_tend_s(i) + ptend_zm_s(i,k) * state_pdel(i,k)
            zm_avg_tend_q(i) = zm_avg_tend_q(i) + ptend_zm_q(i,k) * state_pdel(i,k)
            pdel_sum(i) = pdel_sum(i) + state_pdel(i,k)
         end do
         ! normalize integrated ZM tendencies by total mass
         zm_avg_tend_s(i) = zm_avg_tend_s(i) / pdel_sum(i)
         zm_avg_tend_q(i) = zm_avg_tend_q(i) / pdel_sum(i)
         ! calculate diagnostic zm_depth
         zm_depth(i) = state_pint(i,pver+1) - state_pmid(i,jctop(i))
      else
         zm_avg_tend_s(i) = 0
         zm_avg_tend_q(i) = 0
         zm_depth(i)  = 0
      end if
   end do

   !----------------------------------------------------------------------------
   ! Note: To conserve total energy we need to account for the kinteic energy tendency
   ! which we can obtain from the velocity tendencies based on the following:
   !   KE_new = (u_new^2 + v_new^2)/2 
   !          = [ (u_old+du)^2 + (v_old+dv)^2 ]/2
   !          = [ ( u_old^2 + 2*u_old*du + du^2 ) + ( v_old^2 + 2*v_old*dv + dv^2 ) ]/2
   !          = ( u_old^2 + v_old^2 )/2 + ( 2*u_old*du + du^2 + 2*v_old*dv + dv^2 )/2
   !          = KE_old + [ 2*u_old*du + du^2 + 2*v_old*dv + dv^2 ] /2

   !----------------------------------------------------------------------------
   ! calculate MCSP tendencies

   do i = 1,ncol

      ! check that ZM produced tendencies over a depth that exceeds the threshold
      if ( zm_depth(i) >= MCSP_conv_depth_min ) then
         ! check that ZM provided a non-zero column total heating
         if ( zm_avg_tend_s(i) > 0 ) then
            ! check that there is sufficient wind shear to justify coherent organization
            if ( abs(mcsp_shear(i)).ge.mcsp_shear_min .and. &
                 abs(mcsp_shear(i)).lt.mcsp_shear_max ) then
               do k = jctop(i),pver

                  ! See eq 7-8 of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
                  pdepth_mid_k = state_pint(i,pver+1) - state_pmid(i,k)
                  pdepth_total = state_pint(i,pver+1) - state_pmid(i,jctop(i))

                  ! specify the assumed vertical structure
                  if (do_mcsp_t) mcsp_tend_s(i,k) = -1*zm_param%mcsp_t_coeff * sin(2.0_r8*zm_const%pi*(pdepth_mid_k/pdepth_total))
                  if (do_mcsp_q) mcsp_tend_q(i,k) = -1*zm_param%mcsp_q_coeff * sin(2.0_r8*zm_const%pi*(pdepth_mid_k/pdepth_total))
                  if (do_mcsp_u) mcsp_tend_u(i,k) =    zm_param%mcsp_u_coeff * (cos(zm_const%pi*(pdepth_mid_k/pdepth_total)))
                  if (do_mcsp_v) mcsp_tend_v(i,k) =    zm_param%mcsp_v_coeff * (cos(zm_const%pi*(pdepth_mid_k/pdepth_total)))

                  ! scale the vertical structure by the ZM heating/drying tendencies
                  if (do_mcsp_t) mcsp_tend_s(i,k) = zm_avg_tend_s(i) * mcsp_tend_s(i,k)
                  if (do_mcsp_q) mcsp_tend_q(i,k) = zm_avg_tend_q(i) * mcsp_tend_q(i,k)

                  ! integrate the DSE/qv tendencies for energy/mass fixer
                  if (do_mcsp_t) mcsp_avg_tend_s(i) = mcsp_avg_tend_s(i) + mcsp_tend_s(i,k) * state_pdel(i,k) / pdel_sum(i)
                  if (do_mcsp_q) mcsp_avg_tend_q(i) = mcsp_avg_tend_q(i) + mcsp_tend_q(i,k) * state_pdel(i,k) / pdel_sum(i)

                  ! integrate the change in kinetic energy (KE) for energy fixer
                  if (do_mcsp_u.or.do_mcsp_v) then
                     tend_k = ( 2.0_r8*mcsp_tend_u(i,k)*ztodt*state_u(i,k) + mcsp_tend_u(i,k)*mcsp_tend_u(i,k)*ztodt*ztodt &
                               +2.0_r8*mcsp_tend_v(i,k)*ztodt*state_v(i,k) + mcsp_tend_v(i,k)*mcsp_tend_v(i,k)*ztodt*ztodt )/2.0_r8/ztodt
                     mcsp_avg_tend_k(i) = mcsp_avg_tend_k(i) + tend_k*state_pdel(i,k) / pdel_sum(i)
                  end if

               end do ! k = jctop(i),pver
            end if ! shear threshold
         end if ! zm_avg_tend_s(i) > 0
      end if ! zm_depth(i) >= MCSP_conv_depth_min
   end do

   !----------------------------------------------------------------------------
   ! calculate final output tendencies

   mcsp_dt_out(1:ncol,1:pver) = 0
   mcsp_dq_out(1:ncol,1:pver) = 0
   mcsp_du_out(1:ncol,1:pver) = 0
   mcsp_dv_out(1:ncol,1:pver) = 0

   mcsp_freq(1:ncol) = 0

   do i = 1,ncol
      do k = jctop(i),pver

         ! update frequency if MCSP contributes any tendency in the column
         if ( abs(mcsp_tend_s(i,k)).gt.0 .or. abs(mcsp_tend_q(i,k)).gt.0 .or.&
              abs(mcsp_tend_u(i,k)).gt.0 .or. abs(mcsp_tend_v(i,k)).gt.0 ) then
            mcsp_freq(i) = 1
         end if

         ! subtract mass weighted average tendencies for energy/mass conservation
         mcsp_dt_out(i,k) = mcsp_tend_s(i,k) - mcsp_avg_tend_s(i) 
         mcsp_dq_out(i,k) = mcsp_tend_q(i,k) - mcsp_avg_tend_q(i)
         mcsp_du_out(i,k) = mcsp_tend_u(i,k)
         mcsp_dv_out(i,k) = mcsp_tend_v(i,k)

         ! make sure kinetic energy correction is added to DSE tendency
         ! to conserve total energy whenever momentum tendencies are calculated
         if (do_mcsp_u.or.do_mcsp_v) then
            mcsp_dt_out(i,k) = mcsp_dt_out(i,k) - mcsp_avg_tend_k(i)
         end if

         ! update output tendencies
         if (do_mcsp_t) ptend_s(i,k) = ptend_s(i,k) + mcsp_dt_out(i,k)
         if (do_mcsp_q) ptend_q(i,k) = ptend_q(i,k) + mcsp_dq_out(i,k)
         if (do_mcsp_u) ptend_u(i,k) = ptend_u(i,k) + mcsp_du_out(i,k)
         if (do_mcsp_v) ptend_v(i,k) = ptend_v(i,k) + mcsp_dv_out(i,k)

         ! adjust units for diagnostic outputs
         if (do_mcsp_t) mcsp_dt_out(i,k) = mcsp_dt_out(i,k)/zm_const%cpair

      end do
   end do

   !----------------------------------------------------------------------------
   ! write out MCSP diagnostic history fields
#ifndef SCREAM_CONFIG_IS_CMAKE
   call outfld('MCSP_DT',       mcsp_dt_out, pcols, lchnk )
   call outfld('MCSP_DQ',       mcsp_dq_out, pcols, lchnk )
   call outfld('MCSP_DU',       mcsp_du_out, pcols, lchnk )
   call outfld('MCSP_DV',       mcsp_dv_out, pcols, lchnk )

   call outfld('MCSP_freq',     mcsp_freq,   pcols, lchnk )
   call outfld('MCSP_shear',    mcsp_shear,  pcols, lchnk )
   call outfld('MCSP_zm_depth', zm_depth,    pcols, lchnk )
#endif
   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_tend

!===================================================================================================

end module zm_conv_mcsp
