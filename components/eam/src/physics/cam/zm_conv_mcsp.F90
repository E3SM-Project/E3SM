module zm_conv_mcsp
   !----------------------------------------------------------------------------
   ! Purpose: 
   ! 
   ! References:
   ! 
   ! Moncrieff, M. W., & Liu, C. (2006). Representing convective organization in
   ! prediction models by a hybrid strategy. Journal of the Atmospheric Sciences, 
   ! 63, 3404–3420. https://doi.org/10.1175/JAS3812.1
   ! 
   ! Chen, C.-C., Richter, J. H., Liu, C., Moncrieff, M. W., Tang, Q., Lin, W.,
   !   et al. (2021). Effects of organized convection parameterization on the MJO
   !   and precipitation in E3SMv1. Part I: Mesoscale heating. J. Adv.
   !   Mod. Earth Sys., 13, e2020MS002401, https://doi.org/10.1029/2020MS002401
   ! 
   ! Moncrieff, M. W., C. Liu, and P. Bogenschutz, 2017: Simulation, Modeling, 
   !   and Dynamically Based Parameterization of Organized Tropical Convection 
   !   for Global Climate Models. J. Atmos. Sci., 74, 1363–1380, https://doi.org/10.1175/JAS-D-16-0166.1.
   ! 
   !----------------------------------------------------------------------------
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use cam_abortutils,   only: endrun
   use zm_conv_types,    only: zm_const_t

   implicit none
   private

   public :: zm_conv_mcsp_init ! ???
   public :: zm_conv_mcsp_tend ! ???
   ! public :: zm_conv_mcsp_outfld ! ???

   real(r8), parameter :: unset_r8   = huge(1.0_r8)
   integer , parameter :: unset_int  = huge(1)

   ! public MCSP                     ! true if running MCSP
   ! public MCSP_heat_coeff          ! MCSP coefficient setting degree of dry static energy transport
   ! public MCSP_moisture_coeff      ! MCSP coefficient setting degree of moisture transport
   ! public MCSP_uwind_coeff         ! MCSP coefficient setting degree of zonal wind transport
   ! public MCSP_vwind_coeff         ! MCSP coefficient setting degree of meridional wind transport

   ! logical,  public :: MCSP
   ! real(r8), public :: MCSP_heat_coeff = unset_r8
   ! real(r8), public :: MCSP_moisture_coeff = unset_r8
   ! real(r8), public :: MCSP_uwind_coeff = unset_r8
   ! real(r8), public :: MCSP_vwind_coeff = unset_r8


   logical,  public :: MCSP_enabled = .false.
   real(r8), public :: MCSP_t_coeff = unset_r8
   real(r8), public :: MCSP_q_coeff = unset_r8
   real(r8), public :: MCSP_u_coeff = unset_r8
   real(r8), public :: MCSP_v_coeff = unset_r8

   logical :: do_MCSP_t = .false.
   logical :: do_MCSP_q = .false.
   logical :: do_MCSP_u = .false.
   logical :: do_MCSP_v = .false.

   real(r8), parameter :: MCSP_storm_speed_pref = 600e2_r8 ! pressure level for winds in MCSP calculation [Pa]
   real(r8), parameter :: MCSP_conv_depth_min   = 700e2_r8 ! pressure thickness of convective heating [Pa]
   real(r8), parameter :: MCSP_shear_min        = 3.0_r8   ! min shear value for MCSP to be active
   real(r8), parameter :: MCSP_shear_max        = 200.0_r8 ! max shear value for MCSP to be active


   ! MCSP
   ! logical  :: doslop
   ! logical  :: doslop_heat
   ! logical  :: doslop_moisture
   ! logical  :: doslop_uwind
   ! logical  :: doslop_vwind
   ! real(r8) :: alpha2, alpha_moisture, alphau, alphav

   ! Notes:
   ! - add pi to zm_const
   ! - combine ptend updates into single line (i.e. no seperate line for energy fix)
   ! - replace *_int_end with *_dis and replace with more sensible name
   ! - replace MCSP_DT with ptend_s - they seem to be equivalent except for cpair normalization
   ! - create MCSP_DQ

!===================================================================================================
contains
!===================================================================================================

subroutine zm_conv_mcsp_init()
   !----------------------------------------------------------------------------
   ! Purpose:
   !----------------------------------------------------------------------------
   use cam_history, only: addfld, horiz_only
   !----------------------------------------------------------------------------
   ! Arguments

   !----------------------------------------------------------------------------
   ! Local variables

   !----------------------------------------------------------------------------

   call addfld('MCSP_freq', horiz_only, 'A', '1',        'MCSP frequency of activation')
   call addfld('MCSP_DT',   (/ 'lev'/), 'A', 'K/s',      'MCSP T tendency')
   call addfld('MCSP_DQ',   (/ 'lev'/), 'A', '???',      'MCSP Q tendency')
   call addfld('MCSP_DU',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP U tendency')
   call addfld('MCSP_DV',   (/ 'lev'/), 'A', 'm/s/day',  'MCSP V tendency')
   call addfld('MCSP_shear',horiz_only, 'A', 'm/s',      'MCSP vertical shear of zonal wind')
   call addfld('ZM_depth',  horiz_only, 'A', 'Pa',       'ZM convection depth for MCSP')

   if ( MCSP_t_coeff > 0 ) do_MCSP_t = .true.
   if ( MCSP_q_coeff > 0 ) do_MCSP_q = .true.
   if ( MCSP_u_coeff > 0 ) do_MCSP_u = .true.
   if ( MCSP_v_coeff > 0 ) do_MCSP_v = .true.

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_init

!===================================================================================================

subroutine zm_conv_mcsp_tend( pcols, ncol, pver, ztodt, jctop, zm_const, &
                              state_pmid, state_pint, state_pdel, &
                              state_s, state_q, state_u, state_v, &
                              ptend_s, ptend_q, ptend_u, ptend_v )
   !----------------------------------------------------------------------------
   ! Purpose:
   !----------------------------------------------------------------------------
   use interpolate_data, only: vertinterp
   !----------------------------------------------------------------------------
   ! Arguments
   integer,                              intent(in   ) :: pcols      ! number of atmospheric columns (max)
   integer,                              intent(in   ) :: ncol       ! number of atmospheric columns (actual)
   integer,                              intent(in   ) :: pver       ! number of mid-point vertical levels
   ! type(zm_param_t),                     intent(in   ) :: zm_param   ! derived type to hold ZM tunable parameters
   real(r8),                             intent(in   ) :: ztodt
   integer,  dimension(pcols)            intent(in   ) :: jctop      ! cloud top level indices
   type(zm_const_t),                     intent(in   ) :: zm_const   ! derived type to hold ZM constants
   real(r8), dimension(pcols,pver),      intent(inout) :: state_pmid ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_pint ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_pdel ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_s    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_q    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_u    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: state_v    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: ptend_s    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: ptend_q    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: ptend_u    ! ?
   real(r8), dimension(pcols,pver),      intent(inout) :: ptend_v    ! ?
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8) :: mcsp_top, mcsp_bot
   real(r8) :: dpg
   real(r8) :: Qcq_adjust
   real(r8), dimension(pcols)      :: Qs_dis       ! tendency?
   real(r8), dimension(pcols)      :: Qq_dis       ! tendency?
   real(r8), dimension(pcols,pver) :: Qms          ! tendency?
   real(r8), dimension(pcols,pver) :: Qmq          ! tendency?
   real(r8), dimension(pcols,pver) :: Qmu          ! tendency?
   real(r8), dimension(pcols,pver) :: Qmv          ! tendency?
   real(r8), dimension(pcols)      :: Qms_int_end
   real(r8), dimension(pcols)      :: Qmq_int_end
   real(r8), dimension(pcols)      :: Pa_int_end   ! column integrated pressure (i.e. total mass)
   real(r8), dimension(pcols)      :: Qs_zmconv    ! column integrated dse tendency from ZM
   real(r8), dimension(pcols)      :: Qv_zmconv    ! column integrated qv tendency from ZM
   real(r8), dimension(pcols)      :: MCSP_freq
   real(r8), dimension(pcols,pver) :: MCSP_DT      ! tendency for output
   real(r8), dimension(pcols,pver) :: MCSP_DQ      ! tendency for output
   real(r8), dimension(pcols,pver) :: MCSP_DU      ! tendency for output
   real(r8), dimension(pcols,pver) :: MCSP_DV      ! tendency for output
   real(r8), dimension(pcols)      :: ZM_depth
   real(r8), dimension(pcols)      :: MCSP_shear
   real(r8), dimension(pcols)      :: storm_u
   real(r8), dimension(pcols)      :: storm_v
   real(r8), dimension(pcols)      :: storm_u_shear
   real(r8), dimension(pcols)      :: storm_v_shear
   !----------------------------------------------------------------------------
   ! initialize variables
   storm_u(1:ncol)    = 0
   storm_v(1:ncol)    = 0
   MCSP_shear         = 0
   ZM_depth           = 0

   Qs_zmconv (1:ncol) = 0
   Qv_zmconv (1:ncol) = 0
   Pa_int_end(1:ncol) = 0
   
   Qms(1:ncol,1:pver) = 0
   Qmq(1:ncol,1:pver) = 0
   Qmu(1:ncol,1:pver) = 0
   Qmv(1:ncol,1:pver) = 0

   !----------------------------------------------------------------------------
   ! Interpolate wind
   call vertinterp(ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_u,storm_u)
   call vertinterp(ncol, pcols, pver, state_pmid, MCSP_storm_speed_pref, state_v,storm_v)

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
   end do
   MCSP_shear = storm_u_shear

   !----------------------------------------------------------------------------
   ! calculate shear
   do i = 1,ncol
      if ( jctop(i).ne.pver ) ZM_depth(i) = state_pint(i,pver+1) - state_pmid(i,jctop(i))
   end do

   !----------------------------------------------------------------------------
   ! calculate column integrated DSE/qv tendencies from ZM
   do i = 1,ncol
      do k = jctop(i),pver
         Qs_zmconv(i) = Qs_zmconv(i) + ptend_s(i,k) * state_pdel(i,k)
         Qv_zmconv(i) = Qv_zmconv(i) + ptend_q(i,k) * state_pdel(i,k)
         Pa_int_end(i) = Pa_int_end(i) + state_pdel(i,k)
      end do
   end do
   
   ! normalize by total mass
   do i = 1,ncol
      if (jctop(i) .ne. pver) then
         Qs_zmconv(i) = Qs_zmconv(i) / Pa_int_end(i)
         Qv_zmconv(i) = Qv_zmconv(i) / Pa_int_end(i)
      else
         Qs_zmconv(i) = 0
         Qv_zmconv(i) = 0
      end if
   end do

   !----------------------------------------------------------------------------
   ! ???
   do i = 1,ncol
      Qms_int_end(i) = 0
      Qmq_int_end(i) = 0
      ! Pa_int_end(i)  = 0
      Qs_dis(i)      = 0
      Qq_dis(i)      = 0

      if ( (state_pint(i,pver+1)-state_pmid(i,jctop(i))) >= MCSP_conv_depth_min ) then
         if ( abs(MCSP_shear(i)).ge.MCSP_shear_min .and. &
              abs(MCSP_shear(i)).lt.MCSP_shear_max .and. &
              Qs_zmconv(i).gt.0 ) then
            do k = jctop(i),pver
               
               ! See eq (7) of Moncrieff et al. (2017) - also eq (5) of Moncrieff & Liu (2006)
               mcsp_top = state_pint(i,pver+1) - state_pmid(i,k)
               mcsp_bot = state_pint(i,pver+1) - state_pmid(i,jctop(i))

               Qms(i,k) = -1 * Qs_zmconv(i) * MCSP_t_coeff * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
               Qmq(i,k) = -1 * Qv_zmconv(i) * MCSP_q_coeff * sin(2.0_r8*pi*(mcsp_top/mcsp_bot))
                  
               ! convert units from kg/kg/s to J/s ????
               Qmq(i,k) = Qmq(i,k)/2500000.0_r8/4.0_r8

               ! See eq (8) of Moncrieff et al. (2017)
               Qmu(i,k) = MCSP_u_coeff * (cos(pi*(mcsp_top/mcsp_bot)))
               Qmv(i,k) = MCSP_v_coeff * (cos(pi*(mcsp_top/mcsp_bot)))

               dpg = state_pdel(i,k)/zm_const%gravit

               Qms_int_end(i) = Qms_int_end(i) + Qms(i,k) * dpg
               Qms_int_end(i) = Qms_int_end(i) + (2.0_r8*Qmu(i,k)*ztodt*state_u(i,k)+ Qmu(i,k)*Qmu(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt
               Qms_int_end(i) = Qms_int_end(i) + (2.0_r8*Qmv(i,k)*ztodt*state_v(i,k)+ Qmv(i,k)*Qmv(i,k)*ztodt*ztodt)/2.0_r8 * dpg/ztodt

               Qmq_int_end(i) = Qmq_int_end(i) + Qmq(i,k) * dpg
               ! Pa_int_end(i)  = Pa_int_end(i) + state_pdel(i,k)
            end do
            Qs_dis(i) = Qms_int_end(i) / Pa_int_end(i)
            Qq_dis(i) = Qmq_int_end(i) / Pa_int_end(i)
         end if
      end if
   end do

   !----------------------------------------------------------------------------
   ! ???
   MCSP_DT(1:ncol,1:pver) = 0
   MCSP_DQ(1:ncol,1:pver) = 0
   MCSP_DU(1:ncol,1:pver) = 0
   MCSP_DV(1:ncol,1:pver) = 0

   MCSP_freq(1:ncol) = 0

   do i = 1,ncol
      do k = jctop(i),pver
         
         Qcq_adjust   = ptend_q(i,k) - Qq_dis(i)*zm_const%gravit
         ptend_s(i,k) = ptend_s(i,k) - Qs_dis(i)*zm_const%gravit ! energy fixer

         MCSP_DT(i,k) = -Qs_dis(i)*zm_const%gravit + Qms(i,k)
         
         if (abs(Qms(i,k)).gt.0 .and. abs(Qmu(i,k)).gt.0) MCSP_freq(i) = 1

         ! if (    abs(Qms(i,k)).gt.0 &
         !    .or. abs(Qmq(i,k)).gt.0 &
         !    .or. abs(Qmu(i,k)).gt.0 &
         !    .or. abs(Qmv(i,k)).gt.0 ) then
         !    MCSP_freq(i) = 1
         ! end if

         if (do_MCSP_t) ptend_s(i,k) = ptend_s(i,k) + Qms(i,k)
         if (do_MCSP_q) ptend_q(i,k) = Qcq_adjust + Qmq(i,k)
         if (do_MCSP_u) ptend_u(i,k) = Qmu(i,k)
         if (do_MCSP_v) ptend_v(i,k) = Qmv(i,k)
      end do
   end do

   MCSP_DT(1:ncol,1:pver) = MCSP_DT(1:ncol,1:pver)/cpair
   MCSP_DQ(1:ncol,1:pver) = ptend_q(1:ncol,1:pver)
   MCSP_DU(1:ncol,1:pver) = ptend_u(1:ncol,1:pver)*86400.0_r8
   MCSP_DV(1:ncol,1:pver) = ptend_v(1:ncol,1:pver)*86400.0_r8
   !----------------------------------------------------------------------------
   ! call zm_conv_mcsp_outfld()

   if (do_MCSP_t) call outfld('MCSP_DT    ',MCSP_DT, pcols, lchnk )
   if (do_MCSP_q) call outfld('MCSP_DQ    ',MCSP_DQ, pcols, lchnk )
   if (do_MCSP_u) call outfld('MCSP_DU    ',MCSP_DU, pcols, lchnk )
   if (do_MCSP_v) call outfld('MCSP_DV    ',MCSP_DV, pcols, lchnk )
   call outfld('MCSP_freq  ',MCSP_freq,  pcols, lchnk )
   call outfld('MCSP_shear ',MCSP_shear, pcols, lchnk )
   call outfld('ZM_depth   ',ZM_depth,   pcols, lchnk )

   !----------------------------------------------------------------------------
   return

end subroutine zm_conv_mcsp_tend

!===================================================================================================

! subroutine zm_conv_mcsp_outfld( pcols, ncol, pver, ztodt, jctop, &
!                                 )
!    !----------------------------------------------------------------------------
!    ! Purpose:
!    !----------------------------------------------------------------------------
!    use cam_history, only: outfld
!    !----------------------------------------------------------------------------
!    ! Arguments
!    integer,                              intent(in   ) :: pcols      ! number of atmospheric columns (max)
!    integer,                              intent(in   ) :: ncol       ! number of atmospheric columns (actual)
!    integer,                              intent(in   ) :: pver       ! number of mid-point vertical levels
!    integer,  dimension(pcols)            intent(in   ) :: jctop      ! cloud top level indices

!    !----------------------------------------------------------------------------
!    ! Local variables

!    !----------------------------------------------------------------------------
!    if (do_MCSP_t) call outfld('MCSP_DT    ',MCSP_DT,   pcols, lchnk )
!    if (do_MCSP_q) call outfld('MCSP_DQ    ',MCSP_DQ,   pcols, lchnk )
!    if (do_MCSP_u) call outfld('MCSP_DU    ',ptend_u*86400.0_r8, pcols, lchnk )
!    if (do_MCSP_v) call outfld('MCSP_DV    ',ptend_v*86400.0_r8, pcols, lchnk )
!    call outfld('MCSP_freq  ',MCSP_freq, pcols, lchnk )
!    call outfld('MCSP_shear ',MCSP_shear, pcols, lchnk )
!    call outfld('ZM_depth   ',ZM_depth,   pcols, lchnk )

!    !----------------------------------------------------------------------------
!    return

! end subroutine zm_conv_mcsp_outfld

!===================================================================================================

! subroutine xxx()
!    !----------------------------------------------------------------------------
!    ! Purpose:
!    !----------------------------------------------------------------------------
!    ! use ???
!    !----------------------------------------------------------------------------
!    ! Arguments

!    !----------------------------------------------------------------------------
!    ! Local variables

!    !----------------------------------------------------------------------------
   
!    !----------------------------------------------------------------------------
!    return

! end subroutine xxx

!===================================================================================================

end module zm_conv_mcsp
