module wetdep

!----------------------------------------------------------------------- 
!
! Wet deposition routines for both aerosols and gas phase constituents.
! 
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid,       only: pcols, pver
use physconst,    only: gravit, rair, tmelt
use phys_control, only: cam_physpkg_is
use cam_logfile,  only: iulog

implicit none
save
private

public :: wetdepa_v1  ! scavenging codes for very soluble aerosols -- CAM4 version
public :: wetdepa_v2  ! scavenging codes for very soluble aerosols -- CAM5 version
public :: wetdepg     ! scavenging of gas phase constituents by henry's law
public :: clddiag     ! calc of cloudy volume and rain mixing ratio
public :: faer_resusp_vs_fprec_evap_mpln

public :: wetdep_inputs_t
public :: wetdep_init
public :: wetdep_inputs_set

real(r8), parameter :: cmftau = 3600._r8
real(r8), parameter :: rhoh2o = 1000._r8            ! density of water
real(r8), parameter :: molwta = 28.97_r8            ! molecular weight dry air gm/mole

type wetdep_inputs_t
   real(r8), pointer :: cldt(:,:) => null()  ! cloud fraction
   real(r8), pointer :: qme(:,:) => null()
   real(r8), pointer :: prain(:,:) => null()
   real(r8), pointer :: evapr(:,:) => null()
   real(r8) :: cldcu(pcols,pver)     ! convective cloud fraction, currently empty
   real(r8) :: evapc(pcols,pver)     ! Evaporation rate of convective precipitation
   real(r8) :: cmfdqr(pcols,pver)    ! convective production of rain
   real(r8) :: conicw(pcols,pver)    ! convective in-cloud water
   real(r8) :: totcond(pcols, pver)  ! total condensate
   real(r8) :: cldv(pcols,pver)      ! cloudy volume undergoing wet chem and scavenging
   real(r8) :: cldvcu(pcols,pver)    ! Convective precipitation area at the top interface of current layer
   real(r8) :: cldvst(pcols,pver)    ! Stratiform precipitation area at the top interface of current layer 
end type wetdep_inputs_t

integer :: cld_idx             = 0
integer :: qme_idx             = 0 
integer :: prain_idx           = 0 
integer :: nevapr_idx          = 0 

integer :: icwmrdp_idx         = 0 
integer :: icwmrsh_idx         = 0 
integer :: rprddp_idx          = 0 
integer :: rprdsh_idx          = 0 
integer :: sh_frac_idx         = 0 
integer :: dp_frac_idx         = 0 
integer :: nevapr_shcu_idx     = 0 
integer :: nevapr_dpcu_idx     = 0 
integer :: ixcldice, ixcldliq

logical :: pergro_mods         = .false.

!==============================================================================
contains
!==============================================================================

!==============================================================================
!==============================================================================
subroutine wetdep_init()
  use physics_buffer, only: pbuf_get_index
  use constituents,   only: cnst_get_ind
  use phys_control,   only: phys_getopts

  cld_idx             = pbuf_get_index('CLD')    
  qme_idx             = pbuf_get_index('QME')    
  prain_idx           = pbuf_get_index('PRAIN')  
  nevapr_idx          = pbuf_get_index('NEVAPR') 

  icwmrdp_idx         = pbuf_get_index('ICWMRDP') 
  rprddp_idx          = pbuf_get_index('RPRDDP')  
  icwmrsh_idx         = pbuf_get_index('ICWMRSH') 
  rprdsh_idx          = pbuf_get_index('RPRDSH')  
  sh_frac_idx         = pbuf_get_index('SH_FRAC' )
  dp_frac_idx         = pbuf_get_index('DP_FRAC') 
  nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU') 
  nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU') 

  call cnst_get_ind('CLDICE', ixcldice)
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call phys_getopts(pergro_mods_out = pergro_mods)

endsubroutine wetdep_init

!==============================================================================
! gathers up the inputs needed for the wetdepa routines
!==============================================================================
subroutine wetdep_inputs_set( state, pbuf, inputs )
  use phys_control,   only: cam_physpkg_is
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

  ! args

  type(physics_state),  intent(in )  :: state           !! physics state
  type(physics_buffer_desc), pointer :: pbuf(:)         !! physics buffer
  type(wetdep_inputs_t), intent(out) :: inputs          !! collection of wetdepa inputs

  ! local vars

  real(r8), pointer :: icwmrdp(:,:)    ! in cloud water mixing ratio, deep convection
  real(r8), pointer :: rprddp(:,:)     ! rain production, deep convection
  real(r8), pointer :: icwmrsh(:,:)    ! in cloud water mixing ratio, deep convection
  real(r8), pointer :: rprdsh(:,:)     ! rain production, deep convection
  real(r8), pointer :: sh_frac(:,:)    ! Shallow convective cloud fraction
  real(r8), pointer :: dp_frac(:,:)    ! Deep convective cloud fraction
  real(r8), pointer :: evapcsh(:,:)    ! Evaporation rate of shallow convective precipitation >=0.
  real(r8), pointer :: evapcdp(:,:)    ! Evaporation rate of deep    convective precipitation >=0.

  real(r8) :: rainmr(pcols,pver)       ! mixing ratio of rain within cloud volume
  real(r8) :: cldst(pcols,pver)        ! Stratiform cloud fraction

  integer :: itim, ncol

  ncol = state%ncol
  itim = pbuf_old_tim_idx()

  call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qme_idx,         inputs%qme     )
  call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
  call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
  call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
  call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
  call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
  call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
  call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
  call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
  call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
  call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

  inputs%cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
  cldst(:ncol,:)          = inputs%cldt(:ncol,:) - inputs%cldcu(:ncol,:)       ! Stratiform cloud fraction
  inputs%evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)
  inputs%cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

  ! sum deep and shallow convection contributions
  if (cam_physpkg_is('cam5')) then
     ! Dec.29.2009. Sungsu
     inputs%conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                              max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
  else
     inputs%conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
  end if

  inputs%totcond(:ncol,:) = state%q(:ncol,:,ixcldliq) + state%q(:ncol,:,ixcldice)

  call clddiag( state%t,     state%pmid,   state%pdel,   inputs%cmfdqr, inputs%evapc, &
               inputs%cldt,  inputs%cldcu,       cldst,  inputs%qme,    inputs%evapr, &
               inputs%prain, inputs%cldv, inputs%cldvcu, inputs%cldvst,       rainmr, &
                state%ncol )

end subroutine wetdep_inputs_set

subroutine clddiag(t, pmid, pdel, cmfdqr, evapc, &
                   cldt, cldcu, cldst, cme, evapr, &
                   prain, cldv, cldvcu, cldvst, rain, &
                   ncol)

   ! ------------------------------------------------------------------------------------ 
   ! Estimate the cloudy volume which is occupied by rain or cloud water as
   ! the max between the local cloud amount or the
   ! sum above of (cloud*positive precip production)      sum total precip from above
   !              ----------------------------------   x ------------------------
   ! sum above of     (positive precip           )        sum positive precip from above
   ! Author: P. Rasch
   !         Sungsu Park. Mar.2010 
   ! ------------------------------------------------------------------------------------

   ! Input arguments:
   real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
   real(r8), intent(in) :: pmid(pcols,pver)     ! pressure at layer midpoints
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across layers
   real(r8), intent(in) :: cmfdqr(pcols,pver)   ! dq/dt due to convective rainout 
   real(r8), intent(in) :: evapc(pcols,pver)    ! Evaporation rate of convective precipitation ( >= 0 ) 
   real(r8), intent(in) :: cldt(pcols,pver)    ! total cloud fraction
   real(r8), intent(in) :: cldcu(pcols,pver)    ! Cumulus cloud fraction
   real(r8), intent(in) :: cldst(pcols,pver)    ! Stratus cloud fraction
   real(r8), intent(in) :: cme(pcols,pver)      ! rate of cond-evap within the cloud
   real(r8), intent(in) :: evapr(pcols,pver)    ! rate of evaporation of falling precipitation (kg/kg/s)
   real(r8), intent(in) :: prain(pcols,pver)    ! rate of conversion of condensate to precipitation (kg/kg/s)
   integer, intent(in) :: ncol

   ! Output arguments:
   real(r8), intent(out) :: cldv(pcols,pver)     ! fraction occupied by rain or cloud water 
   real(r8), intent(out) :: cldvcu(pcols,pver)   ! Convective precipitation volume
   real(r8), intent(out) :: cldvst(pcols,pver)   ! Stratiform precipitation volume
   real(r8), intent(out) :: rain(pcols,pver)     ! mixing ratio of rain (kg/kg)

   ! Local variables:
   integer  i, k
   real(r8) convfw         ! used in fallspeed calculation; taken from findmcnew
   real(r8) sumppr(pcols)        ! precipitation rate (kg/m2-s)
   real(r8) sumpppr(pcols)       ! sum of positive precips from above
   real(r8) cldv1(pcols)         ! precip weighted cloud fraction from above
   real(r8) lprec                ! local production rate of precip (kg/m2/s)
   real(r8) lprecp               ! local production rate of precip (kg/m2/s) if positive
   real(r8) rho                  ! air density
   real(r8) vfall
   real(r8) sumppr_cu(pcols)     ! Convective precipitation rate (kg/m2-s)
   real(r8) sumpppr_cu(pcols)    ! Sum of positive convective precips from above
   real(r8) cldv1_cu(pcols)      ! Convective precip weighted convective cloud fraction from above
   real(r8) lprec_cu             ! Local production rate of convective precip (kg/m2/s)
   real(r8) lprecp_cu            ! Local production rate of convective precip (kg/m2/s) if positive
   real(r8) sumppr_st(pcols)     ! Stratiform precipitation rate (kg/m2-s)
   real(r8) sumpppr_st(pcols)    ! Sum of positive stratiform precips from above
   real(r8) cldv1_st(pcols)      ! Stratiform precip weighted stratiform cloud fraction from above
   real(r8) lprec_st             ! Local production rate of stratiform precip (kg/m2/s)
   real(r8) lprecp_st            ! Local production rate of stratiform precip (kg/m2/s) if positive
   ! -----------------------------------------------------------------------

   convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
   do i=1,ncol
      sumppr(i) = 0._r8
      cldv1(i) = 0._r8
      sumpppr(i) = 1.e-36_r8
      sumppr_cu(i)  = 0._r8
      cldv1_cu(i)   = 0._r8
      sumpppr_cu(i) = 1.e-36_r8
      sumppr_st(i)  = 0._r8
      cldv1_st(i)   = 0._r8
      sumpppr_st(i) = 1.e-36_r8
   end do

   do k = 1,pver
      do i = 1,ncol
         cldv(i,k) = &
            max(min(1._r8, &
            cldv1(i)/sumpppr(i) &
            )*sumppr(i)/sumpppr(i), &
            cldt(i,k) &
            )
         lprec = pdel(i,k)/gravit &
            *(prain(i,k)+cmfdqr(i,k)-evapr(i,k))
         lprecp = max(lprec,1.e-30_r8)
         cldv1(i) = cldv1(i)  + cldt(i,k)*lprecp
         sumppr(i) = sumppr(i) + lprec
         sumpppr(i) = sumpppr(i) + lprecp

         ! For convective precipitation volume at the top interface of each layer. Neglect the current layer.
         cldvcu(i,k)   = max(min(1._r8,cldv1_cu(i)/sumpppr_cu(i))*(sumppr_cu(i)/sumpppr_cu(i)),0._r8)
         lprec_cu      = (pdel(i,k)/gravit)*(cmfdqr(i,k)-evapc(i,k))
         lprecp_cu     = max(lprec_cu,1.e-30_r8)
         cldv1_cu(i)   = cldv1_cu(i) + cldcu(i,k)*lprecp_cu
         sumppr_cu(i)  = sumppr_cu(i) + lprec_cu
         sumpppr_cu(i) = sumpppr_cu(i) + lprecp_cu

         ! For stratiform precipitation volume at the top interface of each layer. Neglect the current layer.
         cldvst(i,k)   = max(min(1._r8,cldv1_st(i)/sumpppr_st(i))*(sumppr_st(i)/sumpppr_st(i)),0._r8)
         lprec_st      = (pdel(i,k)/gravit)*(prain(i,k)-evapr(i,k))
         lprecp_st     = max(lprec_st,1.e-30_r8)
         cldv1_st(i)   = cldv1_st(i) + cldst(i,k)*lprecp_st
         sumppr_st(i)  = sumppr_st(i) + lprec_st
         sumpppr_st(i) = sumpppr_st(i) + lprecp_st

         rain(i,k) = 0._r8
         if(t(i,k) .gt. tmelt) then
            rho = pmid(i,k)/(rair*t(i,k))
            vfall = convfw/sqrt(rho)
            rain(i,k) = sumppr(i)/(rho*vfall)
            if (rain(i,k).lt.1.e-14_r8) rain(i,k) = 0._r8
         endif
      end do
   end do

end subroutine clddiag

!==============================================================================

! REASTER 08/05/2015
! changed arguments
!    put them in a more logical order
!    optional arguments are now mandatory, and commented out:
!       all "if ( present(xx) )" tests
!       any code for ".not. present(xx)" cases
!    eliminated the sol_fact**_in - now just use sol_fact**

! old argument order
! ubroutine wetdepa_v2( t, p, q, pdel, &
!                       cldt, cldc, cmfdqr, evapc, conicw, precs, conds, &
!                       evaps, cwat, tracer, deltat, &
!                       scavt, iscavt, cldv, cldvcu, cldvst, dlf, fracis, sol_fact, ncol, &
!                       scavcoef, is_strat_cloudborne, rate1ord_cw2pr_st, qqcw, f_act_conv, &
!                       icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt, &  
!                       sol_facti_in, sol_factbi_in, sol_factii_in, &
!                       sol_factic_in, sol_factiic_in, resus_fix ) 
! note - p, q, cldv not needed

! new argument order
subroutine wetdepa_v2( ncol, deltat, &
                       t, p, q, pdel, &
                       cmfdqr, evapc, dlf, conicw, &
                       precs, conds, evaps, cwat, &
                       cldt, cldc, cldv, cldvcu, cldvst, &
                       sol_factb, sol_factbi, sol_facti, sol_factii, sol_factic, sol_factiic, &
                       mam_prevap_resusp_optcc, is_strat_cloudborne, scavcoef, rate1ord_cw2pr_st, f_act_conv, &
                       tracer, qqcw, &
                       fracis, scavt, iscavt, &
                       icscavt, isscavt, bcscavt, bsscavt, rcscavt, rsscavt )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging code for very soluble aerosols
      ! 
      ! Author: P. Rasch
      ! Modified by T. Bond 3/2003 to track different removals
      ! Sungsu Park. Mar.2010 : Impose consistencies with a few changes in physics.
      !-----------------------------------------------------------------------

      use phys_control, only: phys_getopts

      implicit none

      integer, intent(in) :: ncol

      real(r8), intent(in) ::&
         deltat,               &! time step
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cmfdqr(pcols,pver),   &! rate of production of convective precip
! Sungsu
         evapc(pcols,pver),    &! Evaporation rate of convective precipitation
         dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
! Sungsu
         conicw(pcols,pver),   &! convective cloud water
         precs(pcols,pver),    &! rate of production of stratiform precip
         conds(pcols,pver),    &! rate of production of condensate
         evaps(pcols,pver),    &! rate of evaporation of precip
         cwat(pcols,pver),     &! cloud water amount 
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cldv(pcols,pver),     &! total cloud fraction
! Sungsu
         cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer
         cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer
! Sungsu
         tracer(pcols,pver)     ! trace species
      ! If subroutine is called with just sol_fact:
            ! sol_fact is used for both in- and below-cloud scavenging
      ! If subroutine is called with optional argument sol_facti_in:
            ! sol_fact  is used for below cloud scavenging
            ! sol_facti is used for in cloud scavenging

      integer, intent(in) :: mam_prevap_resusp_optcc

!     logical, intent(in) :: resus_fix
      ! rce 2010/05/01
      ! is_strat_cloudborne = .true. if tracer is stratiform-cloudborne aerosol; else .false. 
      logical, intent(in) :: is_strat_cloudborne   
      real(r8), intent(in) :: scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
      ! rate1ord_cw2pr_st = 1st order rate for strat cw to precip (1/s) 
      real(r8), intent(in) :: rate1ord_cw2pr_st(pcols,pver)
      ! qqcw = strat-cloudborne aerosol corresponding to tracer when is_strat_cloudborne==.false.; else 0.0 
      ! f_act_conv = conv-cloud activation fraction when is_strat_cloudborne==.false.; else 0.0 
      real(r8), intent(in) :: f_act_conv(pcols,pver)

      real(r8), intent(in) :: qqcw(pcols,pver)
      ! end rce 2010/05/01

!     real(r8), intent(in) :: sol_fact ! solubility factor (fraction of aer scavenged below & in, or just below or sol_facti is provided)
      real(r8), intent(in) :: sol_factb   ! solubility factor (frac of aerosol scavenged below cloud)
      real(r8), intent(in) :: sol_factbi  ! solubility factor (frac of aerosol scavenged below cloud by ice)
      real(r8), intent(in) :: sol_facti   ! solubility factor (frac of aerosol scavenged in cloud)
      real(r8), intent(in) :: sol_factii  ! solubility factor (frac of aerosol scavenged in cloud by ice)
      real(r8), intent(in) :: sol_factic(pcols,pver)  ! sol_facti for convective clouds
      real(r8), intent(in) :: sol_factiic ! sol_factii for convective clouds
         
      real(r8), intent(out) :: fracis(pcols,pver)  ! fraction of species not scavenged
      real(r8), intent(out) :: scavt(pcols,pver)   ! scavenging tend 
      real(r8), intent(out) :: iscavt(pcols,pver)  ! incloud scavenging tends

      real(r8), intent(out) :: icscavt(pcols,pver)  ! incloud, convective
      real(r8), intent(out) :: isscavt(pcols,pver)  ! incloud, stratiform
      real(r8), intent(out) :: bcscavt(pcols,pver)  ! below cloud, convective
      real(r8), intent(out) :: bsscavt(pcols,pver)  ! below cloud, stratiform
      real(r8), intent(out) :: rcscavt(pcols,pver)  ! resuspension, convective 
      real(r8), intent(out) :: rsscavt(pcols,pver)  ! resuspension, stratiform 

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) fracev(pcols)        ! fraction of precip from above that is evaporating
! Sungsu
      real(r8) fracev_cu(pcols)     ! Fraction of convective precip from above that is evaporating
! Sungsu
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precabc(pcols)       ! conv precip from above (work array)
      real(r8) precabs(pcols)       ! strat precip from above (work array)
      real(r8) precbl               ! precip falling out of level (work array)
      real(r8) precmin              ! minimum convective precip causing scavenging
      real(r8) rat(pcols)           ! ratio of amount available to amount removed
      real(r8) scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(pcols)       ! scavenged tracer flux from above (work array)
      real(r8) srcc                 ! tend for convective rain
      real(r8) srcs                 ! tend for stratiform rain
      real(r8) srct(pcols)          ! work variable
      real(r8) tracab(pcols)        ! column integrated tracer amount
!      real(r8) vfall                ! fall speed of precip
      real(r8) fins                 ! fraction of rem. rate by strat rain
      real(r8) finc                 ! fraction of rem. rate by conv. rain
      real(r8) srcs1                ! work variable
      real(r8) srcs2                ! work variable
      real(r8) tc                   ! temp in celcius
      real(r8) weight               ! fraction of condensate which is ice
      real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
      real(r8) odds                 ! limit on removal rate (proportional to prec)
      real(r8) dblchek(pcols)
      logical :: found

    ! Jan.16.2009. Sungsu for wet scavenging below clouds.
    ! real(r8) cldovr_cu(pcols)     ! Convective precipitation area at the base of each layer
    ! real(r8) cldovr_st(pcols)     ! Stratiform precipitation area at the base of each layer

      real(r8) tracer_incu
      real(r8) tracer_mean

    ! End by Sungsu

!     real(r8) sol_facti,  sol_factb  ! in cloud and below cloud fraction of aerosol scavenged
!     real(r8) sol_factii, sol_factbi ! in cloud and below cloud fraction of aerosol scavenged by ice
!     real(r8) sol_factic(pcols,pver)             ! sol_facti for convective clouds
!     real(r8) sol_factiic            ! sol_factii for convective clouds
      ! sol_factic & solfact_iic added for MODAL_AERO.  
      ! For stratiform cloud, cloudborne aerosol is treated explicitly,
      !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
      ! For convective cloud, cloudborne aerosol is not treated explicitly,
      !    and sol_factic is 1.0 for both cloudborne and interstitial.

      integer  jstrcnv

      real(r8), parameter :: prec_smallaa = 1.0e-30_r8  ! 1e-30 kg/m2/s (or mm/s) = 3.2e-23 mm/yr
      real(r8), parameter :: x_smallaa = 1.0e-30_r8

      real(r8) arainx
      real(r8) evapx
      real(r8) pprdx
      real(r8) precabc_base(pcols)  ! conv precip at an effective cloud base for calculations in a particular layer
      real(r8) precabs_base(pcols)  ! strat precip at an effective cloud base for calculations in a particular layer
      real(r8) precabx_old, precabx_tmp, precabx_new
      real(r8) precabx_base_old, precabx_base_tmp, precabx_base_new
      real(r8) precnums_base(pcols)  ! stratiform precip number flux at the bottom of a particular layer
      real(r8) precnumc_base(pcols)  ! convective precip number flux at the bottom of a particular layer
      real(r8) precnumx_base_old, precnumx_base_tmp, precnumx_base_new
      real(r8) resusp_c          ! aerosol mass re-suspension in a particular layer from convective rain
      real(r8) resusp_s          ! aerosol mass re-suspension in a particular layer from stratiform rain

      real(r8) resusp_x
      real(r8) resusp_c_sv(pcols)
      real(r8) resusp_s_sv(pcols)
      real(r8) scavabx_old, scavabx_tmp, scavabx_new
      real(r8) srcx
      real(r8) tmpa, tmpb
      real(r8) u_old, u_tmp
      real(r8) x_old, x_tmp, x_ratio

      
#ifdef CRM_NZ
      ! crm_nz is used to disable warnings above the CRM with MMF
      integer, parameter :: crm_nz = CRM_NZ   
#else
      ! if not MMF, CRM_NZ is not defined, so set to zero to avoid build error
      integer, parameter :: crm_nz = 0
#endif
      logical use_SPCAM
      call phys_getopts( use_SPCAM_out = use_SPCAM)

! ------------------------------------------------------------------------
!      omsm = 1.-1.e-10          ! used to prevent roundoff errors below zero
      omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero
      precmin =  0.1_r8/8.64e4_r8      ! set critical value to 0.1 mm/day in kg/m2/s

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
!      vfall = 4.
	
      ! default (if other sol_facts aren't in call, set all to required sol_fact
!     sol_facti = sol_fact
!     sol_factb = sol_fact
!     sol_factii = sol_fact
!     sol_factbi = sol_fact

!     if ( present(sol_facti_in) )  sol_facti = sol_facti_in
!     if ( present(sol_factii_in) )  sol_factii = sol_factii_in
!     if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in
!     sol_facti = sol_facti_in
!     sol_factii = sol_factii_in
!     sol_factbi = sol_factbi_in

!     sol_factic  = sol_facti
!     sol_factiic = sol_factii
!     if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
!     if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in
!     sol_factic  = sol_factic_in
!     sol_factiic = sol_factiic_in

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds, 
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.

      do i = 1,pcols
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0

         precabs_base(i) = 0.0_r8
         precabc_base(i) = 0.0_r8
         precnums_base(i) = 0.0_r8
         precnumc_base(i) = 0.0_r8
         resusp_c_sv(i) = 0.0_r8
         resusp_s_sv(i) = 0.0_r8

       ! Jan.16. Sungsu 
       ! I added below to compute vertically projected cumulus and stratus fractions from the top to the
       ! current model layer by assuming a simple independent maximum overlapping assumption for 
       ! each cloud.
       ! cldovr_cu(i) = 0._r8
       ! cldovr_st(i) = 0._r8
       ! End by Sungsu

      end do

main_k_loop: &
      do k = 1,pver
main_i_loop: &
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
            weight = 0._r8                                 ! assume no ice

            pdog = pdel(i,k)/gravit

            ! ****************** Evaporation **************************
            ! calculate the fraction of strat precip from above 
            !                 which evaporates within this layer
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            ! trap to ensure reasonable ratio bounds
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))

! Sungsu : Same as above but convective precipitation part
            fracev_cu(i) = evapc(i,k)*pdel(i,k)/gravit/max(1.e-12_r8,precabc(i))
            fracev_cu(i) = max(0._r8,min(1._r8,fracev_cu(i)))
! Sungsu

            if (mam_prevap_resusp_optcc <= 0) then
               fracev(i) = 0.0_r8
               fracev_cu(i) = 0.0_r8
            endif

            ! ****************** Convection ***************************
            ! now do the convective scavenging

            ! set odds proportional to fraction of the grid box that is swept by the 
            ! precipitation =precabc/rhoh20*(area of sphere projected on plane
            !                                /volume of sphere)*deltat
            ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
            ! unless the fraction of the area that is cloud is less than odds, in which
            ! case use the cloud fraction (assumes precabs is in kg/m2/s)
            ! is really: precabs*3/4/1000./1e-3*deltat
            ! here I use .1 from Balkanski
            !
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/gravit*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1.-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1.-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1.-weight) &
            !         /deltat 

            ! fraction of convective cloud water converted to rain
            ! Dec.29.2009 : Sungsu multiplied cldc(i,k) to conicw(i,k) below
            ! fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            ! fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,cldc(i,k)*conicw(i,k))
            ! Sungsu: Below new formula of 'fracp' is necessary since 'conicw' is a LWC/IWC
            !         that has already precipitated out, that is, 'conicw' does not contain
            !         precipitation at all ! 
              fracp = cmfdqr(i,k)*deltat/max(1.e-12_r8,cldc(i,k)*conicw(i,k)+(cmfdqr(i,k)+dlf(i,k))*deltat) ! Sungsu.Mar.19.2010.
            ! Dec.29.2009
            ! Note cmfdrq can be negative from evap of rain, so constrain it <-- This is wrong. cmfdqr does not
            ! contain evaporation of precipitation.
            fracp = max(min(1._r8,fracp),0._r8)
            ! remove that amount from within the convective area
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat ! liquid only
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)/deltat             ! any condensation
!           srcs1 = 0.
!           Jan.02.2010. Sungsu : cldt --> cldc below.
            ! rce 2010/05/01
!           if (present(is_strat_cloudborne)) then  ! Tianyi, 2011/03/29
               if ( is_strat_cloudborne ) then
                  ! only strat in-cloud removal affects strat-cloudborne aerosol
                  srcs1 = 0._r8
               else
                  tracer_incu = f_act_conv(i,k)*(tracer(i,k)+& 
                       min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k)))))))              
                  srcs1 = sol_factic(i,k)*cldc(i,k)*fracp*tracer_incu*(1._r8-weight)/deltat &  ! Liquid
                       + sol_factiic    *cldc(i,k)*fracp*tracer_incu*(weight)/deltat          ! Ice
               end if
!           else
!              srcs1 = sol_factic(i,k)*cldc(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat &  ! liquid
!                   +  sol_factiic*cldc(i,k)*fracp*tracer(i,k)*(weight)/deltat      ! ice
!           end if


            !--mcb

            ! scavenge below cloud

            !            cldmabc(i) = max(cldc(i,k),cldmabc(i))
            !            cldmabc(i) = max(cldt(i,k),cldmabc(i))
            ! cldmabc(i) = max(cldv(i,k),cldmabc(i))
            ! cldmabc(i) = cldv(i,k)
            cldmabc(i) = cldvcu(i,k)

            ! Jan. 16. 2010. Sungsu
            ! cldmabc(i) = cldmabc(i) * cldovr_cu(i) / max( 0.01_r8, cldovr_cu(i) + cldovr_st(i) )
            ! End by Sungsu

            ! rce 2010/05/01
!           if (present(is_strat_cloudborne)) then  ! Tianyi, 2011/03/29
               if ( is_strat_cloudborne ) then
                  ! only strat in-cloud removal affects strat-cloudborne aerosol
                  srcs2 = 0._r8
               else
                  tracer_mean = tracer(i,k)*(1._r8-cldc(i,k)*f_act_conv(i,k))-cldc(i,k)*f_act_conv(i,k)*&
                       min(qqcw(i,k),tracer(i,k)*((cldt(i,k)-cldc(i,k))/max(0.01_r8,(1._r8-(cldt(i,k)-cldc(i,k))))))
                  tracer_mean = max(0._r8,tracer_mean) 
                  odds  = max(min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8)*scavcoef(i,k)*deltat),0._r8) ! Dana and Hales coefficient (/mm)
                  srcs2 = sol_factb *cldmabc(i)*odds*tracer_mean*(1._r8-weight)/deltat & ! Liquid
                       + sol_factbi*cldmabc(i)*odds*tracer_mean*(weight)/deltat         ! Ice
               end if
!           else
!              odds=max( &
!                   min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8) &
!                   *scavcoef(i,k)*deltat),0._r8) ! Dana and Hales coefficient (/mm)
!              srcs2 = sol_factb*cldmabc(i)*odds*tracer(i,k)*(1._r8-weight)/deltat & ! liquid
!                   +  sol_factbi*cldmabc(i)*odds*tracer(i,k)*(weight)/deltat    !ice
!           end if


            !Note that using the temperature-determined weight doesn't make much sense here


            srcc = srcs1 + srcs2  ! convective tend by both processes
            finc = srcs1/(srcc + 1.e-36_r8) ! fraction in-cloud

            ! ****************** Stratiform ***********************
            ! now do the stratiform scavenging

            ! incloud scavenging

            ! rce 2010/05/01
!           if(present(is_strat_cloudborne)) then  ! Tianyi 2011/03/29
               if ( is_strat_cloudborne ) then
                  ! new code for stratiform incloud scav of cloudborne (modal) aerosol 
                  ! >> use the 1st order cw to precip rate calculated in microphysics routine
                  ! >> cloudborne aerosol resides in cloudy portion of grid cell, so do not apply "cldt" factor
                  ! fracp = rate1ord_cw2pr_st(i,k)*deltat
                  ! fracp = max(0._r8,min(1._r8,fracp))
                  fracp = precs(i,k)*deltat/max(cwat(i,k)+precs(i,k)*deltat,1.e-12_r8) ! Sungsu. Mar.19.2010.
                  fracp = max(0._r8,min(1._r8,fracp))
                  srcs1 = sol_facti *fracp*tracer(i,k)/deltat*(1._r8-weight) &  ! Liquid
                       + sol_factii*fracp*tracer(i,k)/deltat*(weight)          ! Ice
               else
                  ! strat in-cloud removal only affects strat-cloudborne aerosol
                  srcs1 = 0._r8
               end if
!           else 
!              ! fracp is the fraction of cloud water converted to precip
!              ! Sungsu modified fracp as the convectiv case.
!              !        Below new formula by Sungsu of 'fracp' is necessary since 'cwat' is a LWC/IWC
!              !        that has already precipitated out, that is, 'cwat' does not contain
!              !        precipitation at all ! 
!              !            fracp =  precs(i,k)*deltat/max(cwat(i,k),1.e-12_r8)
!              fracp =  precs(i,k)*deltat/max(cwat(i,k)+precs(i,k)*deltat,1.e-12_r8) ! Sungsu. Mar.19.2010.
!              fracp = max(0._r8,min(1._r8,fracp))
!              !            fracp = 0.     ! for debug
!              
!              ! assume the corresponding amnt of tracer is removed
!              !++mcb -- remove cldc; change cldt to cldv 
!              !            srcs1 = (cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat
!              !            srcs1 = cldv(i,k)*fracp*tracer(i,k)/deltat &
!              !            srcs1 = cldt(i,k)*fracp*tracer(i,k)/deltat            ! all condensate
!              !            Jan.02.2010. Sungsu : cldt --> cldt - cldc below.
!              srcs1 = sol_facti*(cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat*(1._r8-weight) &  ! liquid
!                   + sol_factii*(cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat*(weight)       ! ice
!           end if
            ! end rce 2010/05/01


            ! below cloud scavenging

!           volume undergoing below cloud scavenging
!           cldmabs(i) = cldv(i,k)   ! precipitating volume
!           cldmabs(i) = cldt(i,k)   ! local cloud volume
            cldmabs(i) = cldvst(i,k) ! Stratiform precipitation area at the top interface of current layer

            ! Jan. 16. 2010. Sungsu
            ! cldmabs(i) = cldmabs(i) * cldovr_st(i) / max( 0.01_r8, cldovr_cu(i) + cldovr_st(i) )
            ! End by Sungsu

            ! rce 2010/05/01
!           if (present(is_strat_cloudborne)) then  ! Tianyi 2011/03/29
               if ( is_strat_cloudborne ) then
                  ! only strat in-cloud removal affects strat-cloudborne aerosol
                  srcs2 = 0._r8
               else
                  odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
                  odds = max(min(1._r8,odds),0._r8)
                  srcs2 = sol_factb *cldmabs(i)*odds*tracer_mean*(1._r8-weight)/deltat & ! Liquid
                       + sol_factbi*cldmabs(i)*odds*tracer_mean*(weight)/deltat         ! Ice
               end if
!           else
!              odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
!              odds = max(min(1._r8,odds),0._r8)
!              srcs2 =sol_factb*(cldmabs(i)*odds) *tracer(i,k)*(1._r8-weight)/deltat & ! liquid
!                   + sol_factbi*(cldmabs(i)*odds) *tracer(i,k)*(weight)/deltat       ! ice
!           end if
                        
            !Note that using the temperature-determined weight doesn't make much sense here

            srcs = srcs1 + srcs2             ! total stratiform scavenging
            fins=srcs1/(srcs + 1.e-36_r8)    ! fraction taken by incloud processes

            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i).lt.1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm

            
            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            if(pergro_mods) then
               fracp = deltat*srct(i)/(max(cldmabs(i),1.e-4_r8)*max(tracer(i,k),1.e-36_r8))  ! amount removed !BSINGH - phil suggested 2nd approach
            else
               fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  ! amount removed ! original code
            endif
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1._r8 - fracp

            ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            ! Sungsu added cumulus contribution in the below 3 blocks
         

! mam_prevap_resusp_optcc values:
!     0 = no resuspension
!     1 = linear resuspension of aerosol mass or number following original mam coding
!     2 = same as 1 but resuspension tendencies are in rc/sscavt rather than combined with bc/sscavt
!     3 = same as 2 but with some added "xxx = max( 0, xxx)" lines
!   130 = non-linear resuspension of aerosol mass   based on scavenged aerosol mass
!   230 = non-linear resuspension of aerosol number based on raindrop number
resusp_block_aa: &
            if ( mam_prevap_resusp_optcc >= 100) then

jstrcnv_loop_aa: &
            do jstrcnv = 1, 2

! step 1 - load working ("x") variables from stratiform ("s") or convective ("c") variables
            if (jstrcnv == 1) then
               pprdx = precs(i,k)
               evapx = evaps(i,k)
               precabx_old = precabs(i)
               precabx_base_old = precabs_base(i)
               if ( mam_prevap_resusp_optcc <= 130) then
                  scavabx_old = scavab(i)
                  srcx = srcs
               else
                  precnumx_base_old = precnums_base(i)
                  arainx = cldvst(i,min(k+1,pver))
               end if
            else
               pprdx = cmfdqr(i,k)
               evapx = evapc(i,k)
               precabx_old = precabc(i)
               precabx_base_old = precabc_base(i)
               if ( mam_prevap_resusp_optcc <= 130) then
                  scavabx_old = scavabc(i)
                  srcx = srcc
               else
                  precnumx_base_old = precnumc_base(i)
                  arainx = cldvcu(i,min(k+1,pver))
               end if
            end if

! force these to be non-negative
            precabx_base_old = max( 0.0_r8, precabx_base_old )
            precabx_old  = max( 0.0_r8, precabx_old )
            precabx_old  = min( precabx_base_old, precabx_old )
            if ( mam_prevap_resusp_optcc <= 130) then
               scavabx_old = max( 0.0_r8, scavabx_old )
            else
               precnumx_base_old = max( 0.0_r8, precnumx_base_old )
               precnumx_base_tmp = precnumx_base_old
            end if

! step 2 - do evaporation and resuspension
            precabx_base_tmp = precabx_base_old
            tmpa = max( 0.0_r8, evapx*pdel(i,k)/gravit )
            precabx_tmp = max( 0.0_r8, precabx_old - tmpa )
            precabx_tmp  = min( precabx_base_tmp, precabx_tmp )

            if (precabx_tmp < prec_smallaa) then
               ! precip rate is essentially zero so do complete resuspension
               if ( mam_prevap_resusp_optcc <= 130) then
                  ! linear resuspension based on scavenged aerosol mass or number
                  scavabx_tmp = 0.0_r8
                  resusp_x = scavabx_old
               else
                  ! non-linear resuspension of aerosol number based on raindrop number
                  if (precabx_base_old < prec_smallaa) then
                     resusp_x = 0.0_r8
                  else
                     u_old = precabx_old/precabx_base_old
                     u_old = max( 0.0_r8, min( 1.0_r8, u_old ) )
                     x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln( 1.0_r8-u_old, jstrcnv )
                     x_old = max( 0.0_r8, min( 1.0_r8, x_old ) )
                     x_tmp = 0.0_r8
                     resusp_x = max( 0.0_r8, precnumx_base_tmp*(x_old - x_tmp) )
                  end if
               end if
               ! setting both these precip rates to zero causes the resuspension 
               ! calculations to start fresh if there is any more precip production
               precabx_tmp = 0.0_r8
               precabx_base_tmp = 0.0_r8

            else if (evapx <= 0.0_r8) then
               ! no evap so no resuspension
               if ( mam_prevap_resusp_optcc <= 130) then
                  scavabx_tmp = scavabx_old
               end if
               resusp_x = 0.0_r8

            else
               u_old = precabx_old/precabx_base_old
               u_old = max( 0.0_r8, min( 1.0_r8, u_old ) )
               if ( mam_prevap_resusp_optcc <= 130) then
                  ! non-linear resuspension of aerosol mass
                  x_old = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln( 1.0_r8-u_old, jstrcnv )
               else 
                  ! non-linear resuspension of aerosol number based on raindrop number
                  x_old = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln( 1.0_r8-u_old, jstrcnv )
               end if
               x_old = max( 0.0_r8, min( 1.0_r8, x_old ) )

               if (x_old < x_smallaa) then
                  x_tmp = 0.0_r8
                  x_ratio = 0.0_r8
               else
                  u_tmp = precabx_tmp/precabx_base_tmp
                  u_tmp = max( 0.0_r8, min( 1.0_r8, u_tmp ) )
                  u_tmp = min( u_tmp, u_old )
                  if ( mam_prevap_resusp_optcc <= 130) then
                     ! non-linear resuspension of aerosol mass
                     x_tmp = 1.0_r8 - faer_resusp_vs_fprec_evap_mpln( 1.0_r8-u_tmp, jstrcnv )
                  else
                     ! non-linear resuspension of aerosol number based on raindrop number
                     x_tmp = 1.0_r8 - fprecn_resusp_vs_fprec_evap_mpln( 1.0_r8-u_tmp, jstrcnv )
                  end if
                  x_tmp = max( 0.0_r8, min( 1.0_r8, x_tmp ) )
                  x_tmp = min( x_tmp, x_old )
                  x_ratio = x_tmp/x_old
                  x_ratio = max( 0.0_r8, min( 1.0_r8, x_ratio ) )
               end if

               if ( mam_prevap_resusp_optcc <= 130) then
                  ! aerosol mass resuspension
                  scavabx_tmp = max( 0.0_r8, scavabx_old * x_ratio )
                  resusp_x = max( 0.0_r8, scavabx_old - scavabx_tmp )
               else
                  ! number resuspension
                  resusp_x = max( 0.0_r8, precnumx_base_tmp*(x_old - x_tmp) )
               end if 
            end if

! step 3 - do precip production and scavenging
            tmpa = max( 0.0_r8, pprdx*pdel(i,k)/gravit )
            precabx_base_new = max( 0.0_r8, precabx_base_tmp + tmpa )
            precabx_new = max( 0.0_r8, precabx_tmp + tmpa )
            precabx_new = min( precabx_base_new, precabx_new )

            if ( mam_prevap_resusp_optcc <= 130) then
               ! aerosol mass scavenging
               tmpa = max( 0.0_r8, srcx*pdel(i,k)/gravit )
               scavabx_new = max( 0.0_r8, scavabx_tmp + tmpa )
            else
               ! raindrop number increase
               if (precabx_base_new < prec_smallaa) then
                  precnumx_base_new = 0.0_r8
               else if (precabx_base_new > precabx_base_tmp) then
                  ! note - calc rainshaft number flux from rainshaft water flux, 
                  !    then multiply by rainshaft area to get grid-average number flux
                  arainx = max( arainx, 0.01_r8 )
                  tmpa = arainx * flux_precnum_vs_flux_prec_mpln( (precabx_base_new/arainx), jstrcnv )
                  precnumx_base_new = max( 0.0_r8, tmpa )
               else
                  precnumx_base_new = precnumx_base_tmp
               end if
            end if

! step 4 - update stratiform ("s") or convective ("c") variables from working ("x") variables
            if (jstrcnv == 1) then
               resusp_s = resusp_x
               precabs(i) = precabx_new
               precabs_base(i) = precabx_base_new
               if ( mam_prevap_resusp_optcc <= 130) then
                  scavab(i) = scavabx_new
               else
                  precnums_base(i) = precnumx_base_new
               end if
            else
               resusp_c = resusp_x
               precabc(i) = precabx_new
               precabc_base(i) = precabx_base_new
               if ( mam_prevap_resusp_optcc <= 130) then
                  scavabc(i) = scavabx_new
               else
                  precnumc_base(i) = precnumx_base_new
               end if
            end if

            end do jstrcnv_loop_aa


            else resusp_block_aa
               resusp_c = fracev_cu(i)*scavabc(i)
               resusp_s = fracev(i)*scavab(i)

            end if resusp_block_aa

            resusp_s_sv(i) = resusp_s
            resusp_c_sv(i) = resusp_c


            if ( mam_prevap_resusp_optcc <= 3) then
               scavt(i,k) = -srct(i) + (fracev(i)*scavab(i)+fracev_cu(i)*scavabc(i))*gravit/pdel(i,k)
            else
               scavt(i,k) = -srct(i) + (resusp_s+resusp_c)*gravit/pdel(i,k)
            endif

            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

!           if ( present(icscavt) ) icscavt(i,k) = -(srcc*finc) * omsm
!           if ( present(isscavt) ) isscavt(i,k) = -(srcs*fins) * omsm
            icscavt(i,k) = -(srcc*finc) * omsm
            isscavt(i,k) = -(srcs*fins) * omsm

!           if(.not.present(resus_fix)) then
!              if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm +  &
!                   fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)
!              if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
!                   fracev(i)*scavab(i)*gravit/pdel(i,k)
!           endif

!           if(present(resus_fix)) then
!              if ( .not. resus_fix ) then
               if (mam_prevap_resusp_optcc <= 1) then
!                 if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm +  &
!                      fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)
!                 if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
!                      fracev(i)*scavab(i)*gravit/pdel(i,k)
                  bcscavt(i,k) = -(srcc * (1-finc)) * omsm +  &
                       fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)
                  bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
                       fracev(i)*scavab(i)*gravit/pdel(i,k)
                  rcscavt(i,k) = 0.0
                  rsscavt(i,k) = 0.0
               else if (mam_prevap_resusp_optcc == 2 .or. mam_prevap_resusp_optcc == 3) then
!                 if ( present(bcscavt) ) then
!                    if ( present(rcscavt) ) then
                        bcscavt(i,k) = -(srcc * (1-finc)) * omsm                 !RCE
                        rcscavt(i,k) = fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)  !RCE
!                    else
!                       bcscavt(i,k) = -(srcc * (1-finc)) * omsm  &
!                            + fracev_cu(i)*scavabc(i)*gravit/pdel(i,k)
!                    end if
!                 end if
!                 if ( present(bsscavt) ) then
!                    if ( present(rsscavt) ) then
                        bsscavt(i,k) = -(srcs * (1-fins)) * omsm               !RCE
                        rsscavt(i,k) = + fracev(i)*scavab(i)*gravit/pdel(i,k)  !RCE
!                    else
!                       bsscavt(i,k) = -(srcs * (1-fins)) * omsm  &
!                            + fracev(i)*scavab(i)*gravit/pdel(i,k)
!                    end if
!                 end if
               else ! here mam_prevap_resusp_optcc == 130, 210, 230
                  bcscavt(i,k) = -(srcc * (1-finc)) * omsm
                  rcscavt(i,k) = resusp_c*gravit/pdel(i,k)
                  bsscavt(i,k) = -(srcs * (1-fins)) * omsm
                  rsscavt(i,k) = resusp_s*gravit/pdel(i,k)
               endif
!           endif

            dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

            ! now keep track of scavenged mass and precip
            if (mam_prevap_resusp_optcc <= 3) then
               scavab(i) = scavab(i)*(1-fracev(i)) + srcs*pdel(i,k)/gravit
               precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
               scavabc(i) = scavabc(i)*(1-fracev_cu(i)) + srcc*pdel(i,k)/gravit
               precabc(i) = precabc(i) + (cmfdqr(i,k) - evapc(i,k))*pdel(i,k)/gravit
               if (mam_prevap_resusp_optcc == 3) then
                  scavab(i)  = max( 0.0_r8, scavab(i)  )
                  scavabc(i) = max( 0.0_r8, scavabc(i) )
                  precabs(i) = max( 0.0_r8, precabs(i) )
                  precabc(i) = max( 0.0_r8, precabc(i) )
               endif
            endif

            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit

       ! Jan.16.2010. Sungsu
       ! Compute convective and stratiform precipitation areas at the base interface
       ! of current layer. These are for computing 'below cloud scavenging' in the 
       ! next layer below.

       ! cldovr_cu(i) = max( cldovr_cu(i), cldc(i,k) )
       ! cldovr_st(i) = max( cldovr_st(i), max( 0._r8, cldt(i,k) - cldc(i,k) ) )

       ! cldovr_cu(i) = max( 0._r8, min ( 1._r8, cldovr_cu(i) ) )
       ! cldovr_st(i) = max( 0._r8, min ( 1._r8, cldovr_st(i) ) )

       ! End by Sungsu

         end do main_i_loop ! End of i = 1, ncol

         found = .false.
         do i = 1,ncol
            if ( dblchek(i) < 0._r8 ) then
               found = .true.
               exit
            end if
         end do

         ! the log files from MMF runs were getting really large with these warnings due 
         ! to tiny negative values (~ -1e-300) at the top of the model, well above the CRM,
         ! so this block avoids this problem when running the MMF
         if ( use_SPCAM ) then
            if ( found ) then
               if ( k < (pver-crm_nz) ) found = .false.
            end if
         end if

         if ( found ) then
            do i = 1,ncol
               if (dblchek(i) .lt. 0._r8) then
                  write(iulog,*) ' wetdepa_v2: negative value ', i, k, tracer(i,k), &
                     dblchek(i), scavt(i,k), srct(i), rat(i), fracev(i)
                  write(iulog,*) ' wetdepa_v2: negative value ', i, k, &
                       mam_prevap_resusp_optcc, is_strat_cloudborne, tracer(i,k), &
                       dblchek(i), deltat*scavt(i,k), deltat*srct(i), &
                       deltat*resusp_s_sv(i)*gravit/pdel(i,k), deltat*resusp_c_sv(i)*gravit/pdel(i,k)
               endif
            end do
         endif

      end do main_k_loop ! End of k = 1, pver

   end subroutine wetdepa_v2


!==============================================================================
      function flux_precnum_vs_flux_prec_mpln( flux_prec, jstrcnv )
      real(r8) :: flux_precnum_vs_flux_prec_mpln
      real(r8), intent(in) :: flux_prec
      integer,  intent(in) :: jstrcnv

      if (jstrcnv <= 1) then
         flux_precnum_vs_flux_prec_mpln = flux_precnum_vs_flux_prec_mp( flux_prec )
      else
         flux_precnum_vs_flux_prec_mpln = flux_precnum_vs_flux_prec_ln( flux_prec )
      end if

      return
      end function flux_precnum_vs_flux_prec_mpln


!==============================================================================
      function faer_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
      real(r8) :: faer_resusp_vs_fprec_evap_mpln
      real(r8), intent(in) :: fprec_evap
      integer,  intent(in) :: jstrcnv

      if (jstrcnv <= 1) then
         faer_resusp_vs_fprec_evap_mpln = faer_resusp_vs_fprec_evap_mp( fprec_evap )
      else
         faer_resusp_vs_fprec_evap_mpln = faer_resusp_vs_fprec_evap_ln( fprec_evap )
      end if

      return
      end function faer_resusp_vs_fprec_evap_mpln


!==============================================================================
      function fprecn_resusp_vs_fprec_evap_mpln( fprec_evap, jstrcnv )
      real(r8) :: fprecn_resusp_vs_fprec_evap_mpln
      real(r8), intent(in) :: fprec_evap
      integer,  intent(in) :: jstrcnv

      if (jstrcnv <= 1) then
         fprecn_resusp_vs_fprec_evap_mpln = fprecn_resusp_vs_fprec_evap_mp( fprec_evap )
      else
         fprecn_resusp_vs_fprec_evap_mpln = fprecn_resusp_vs_fprec_evap_ln( fprec_evap )
      end if

      return
      end function fprecn_resusp_vs_fprec_evap_mpln


!==============================================================================
      function flux_precnum_vs_flux_prec_mp( flux_prec )
!
!  flux_prec = precipitation mass flux at the cloud base (kg/m^2/s)
!  flux_precnum_vs_flux_prec_mp = precipitation number flux
!     at the cloud base (drops/m^2/s), assuming marshall-palmer raindrop size distribution
!
!
      real(r8) :: flux_precnum_vs_flux_prec_mp
      real(r8), intent(in) :: flux_prec

      real(r8), parameter :: a0 =  1.0885896550304022E+01_r8
      real(r8), parameter :: a1 =  4.3660645528167907E-01_r8

      real(r8) :: x, y
   
      if (flux_prec >= 1.0e-36_r8) then
         x = log( flux_prec )
         y = exp( a0 + a1*x )    
      else
         y = 0.0_r8
      end if
      flux_precnum_vs_flux_prec_mp = y

      return
      end function flux_precnum_vs_flux_prec_mp


!==============================================================================
      function flux_precnum_vs_flux_prec_ln( flux_prec )
!
!  flux_prec = precipitation mass flux at the cloud base (kg/m^2/s)
!  flux_precnum_vs_flux_prec_ln = precipitation number flux
!     at the cloud base (drops/m^2/s), assuming log-normal raindrop size distribution
!
!
      real(r8) :: flux_precnum_vs_flux_prec_ln
      real(r8), intent(in) :: flux_prec

      real(r8), parameter :: a0 =  9.9067806476181524E+00_r8
      real(r8), parameter :: a1 =  4.2690709912134056E-01_r8

      real(r8) :: x, y
   
      if (flux_prec >= 1.0e-36_r8) then
         x = log( flux_prec )
         y = exp( a0 + a1*x )    
      else
         y = 0.0_r8
      end if
      flux_precnum_vs_flux_prec_ln = y

      return
      end function flux_precnum_vs_flux_prec_ln


!==============================================================================
      function faer_resusp_vs_fprec_evap_mp( fprec_evap )
!
!  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
!  faer_resusp_vs_fprec_evap_mp = corresponding fraction of precipitation-borne aerosol
!     flux that is resuspended, assuming marshall-palmer raindrop size distribution
!
!  note that these fractions are relative to the cloud-base fluxes,
!      and not to the layer immediately above fluxes
!
      real(r8) :: faer_resusp_vs_fprec_evap_mp
      real(r8), intent(in) :: fprec_evap

      real(r8), parameter :: a01 =  8.6591133737322856E-02_r8
      real(r8), parameter :: a02 = -1.7389168499601941E+00_r8
      real(r8), parameter :: a03 =  2.7401882373663732E+01_r8
      real(r8), parameter :: a04 = -1.5861714653209464E+02_r8
      real(r8), parameter :: a05 =  5.1338179363011193E+02_r8
      real(r8), parameter :: a06 = -9.6835933124501412E+02_r8
      real(r8), parameter :: a07 =  1.0588489932213311E+03_r8
      real(r8), parameter :: a08 = -6.2184513459217271E+02_r8
      real(r8), parameter :: a09 =  1.5184126886039758E+02_r8
      real(r8), parameter :: x_lox_lin =  5.0000000000000003E-02_r8
      real(r8), parameter :: y_lox_lin =  2.5622471203221014E-03_r8

      real(r8) :: x, y

      x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
      if (x < x_lox_lin) then
         y = y_lox_lin * (x/x_lox_lin)
      else
         y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
           + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
      end if
      faer_resusp_vs_fprec_evap_mp = y

      return
      end function faer_resusp_vs_fprec_evap_mp


!==============================================================================
      function faer_resusp_vs_fprec_evap_ln( fprec_evap )
!
!  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
!  faer_resusp_vs_fprec_evap_ln = corresponding fraction of precipitation-borne aerosol
!     flux that is resuspended, assuming log-normal raindrop size distribution
!
!  note that these fractions are relative to the cloud-base fluxes,
!      and not to the layer immediately above fluxes
!
      real(r8) :: faer_resusp_vs_fprec_evap_ln
      real(r8), intent(in) :: fprec_evap

      real(r8), parameter :: a01 =  6.1944215103685640E-02_r8
      real(r8), parameter :: a02 = -2.0095166685965378E+00_r8
      real(r8), parameter :: a03 =  2.3882460251821236E+01_r8
      real(r8), parameter :: a04 = -1.2695611774753374E+02_r8
      real(r8), parameter :: a05 =  4.0086943562320101E+02_r8
      real(r8), parameter :: a06 = -7.4954272875943707E+02_r8
      real(r8), parameter :: a07 =  8.1701055892023624E+02_r8
      real(r8), parameter :: a08 = -4.7941894659538502E+02_r8
      real(r8), parameter :: a09 =  1.1710291076059025E+02_r8
      real(r8), parameter :: x_lox_lin =  1.0000000000000001E-01_r8
      real(r8), parameter :: y_lox_lin =  6.2227889828044350E-04_r8

      real(r8) :: x, y

      x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
      if (x < x_lox_lin) then
         y = y_lox_lin * (x/x_lox_lin)
      else
         y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
           + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
      end if
      faer_resusp_vs_fprec_evap_ln = y

      return
      end function faer_resusp_vs_fprec_evap_ln


!==============================================================================
      function fprecn_resusp_vs_fprec_evap_mp( fprec_evap )
!
!  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
!  fprecn_resusp_vs_fprec_evap_mp = Rain number evaporation fraction, 
!                                  assuming marshall-palmer raindrop size distribution
!
!  note that these fractions are relative to the cloud-base fluxes,
!      and not to the layer immediately above fluxes
!
      real(r8) :: fprecn_resusp_vs_fprec_evap_mp
      real(r8), intent(in) :: fprec_evap

      real(r8), parameter :: a01 =  4.5461070198414655E+00_r8
      real(r8), parameter :: a02 = -3.0381753620077529E+01_r8
      real(r8), parameter :: a03 =  1.7959619926085665E+02_r8
      real(r8), parameter :: a04 = -6.7152282193785618E+02_r8
      real(r8), parameter :: a05 =  1.5651931323557126E+03_r8
      real(r8), parameter :: a06 = -2.2743927701175126E+03_r8
      real(r8), parameter :: a07 =  2.0004645897056735E+03_r8
      real(r8), parameter :: a08 = -9.7351466279626209E+02_r8
      real(r8), parameter :: a09 =  2.0101198012962413E+02_r8
      real(r8), parameter :: x_lox_lin =  5.0000000000000003E-02_r8
      real(r8), parameter :: y_lox_lin =  1.7005858490684875E-01_r8

      real(r8) :: x, y

      x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
      if (x < x_lox_lin) then
         y = y_lox_lin * (x/x_lox_lin)
      else
         y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
           + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
      end if
      fprecn_resusp_vs_fprec_evap_mp = y

      return
      end function fprecn_resusp_vs_fprec_evap_mp


!==============================================================================
      function fprecn_resusp_vs_fprec_evap_ln( fprec_evap )
!
!  fprec_evap = fraction of precipitation flux that has evaporated (below cloud base)
!  fprecn_resusp_vs_fprec_evap_ln = Rain number evaporation fraction, 
!                                  assuming log-normal raindrop size distribution
!
!  note that these fractions are relative to the cloud-base fluxes,
!      and not to the layer immediately above fluxes
!
      real(r8) :: fprecn_resusp_vs_fprec_evap_ln
      real(r8), intent(in) :: fprec_evap

      real(r8), parameter :: a01 = -5.2335291116884175E-02_r8
      real(r8), parameter :: a02 =  2.7203158069178226E+00_r8
      real(r8), parameter :: a03 =  9.4730878152409375E+00_r8
      real(r8), parameter :: a04 = -5.0573187592544798E+01_r8
      real(r8), parameter :: a05 =  9.4732631441282862E+01_r8
      real(r8), parameter :: a06 = -8.8265926556465814E+01_r8
      real(r8), parameter :: a07 =  3.5247835268269142E+01_r8
      real(r8), parameter :: a08 =  1.5404586576716444E+00_r8
      real(r8), parameter :: a09 = -3.8228795492549068E+00_r8
      real(r8), parameter :: x_lox_lin =  1.0000000000000001E-01_r8
      real(r8), parameter :: y_lox_lin =  2.7247994766566485E-02_r8

      real(r8) :: x, y

      x = max( 0.0_r8, min( 1.0_r8, fprec_evap ) )
      if (x < x_lox_lin) then
         y = y_lox_lin * (x/x_lox_lin)
      else
         y = x*( a01 + x*( a02 + x*( a03 + x*( a04 + x*( a05 &
           + x*( a06 + x*( a07 + x*( a08 + x*a09 ))))))))
      end if
      fprecn_resusp_vs_fprec_evap_ln = y

      return
      end function fprecn_resusp_vs_fprec_evap_ln


!==============================================================================

! This is the frozen CAM4 version of wetdepa.


   subroutine wetdepa_v1( t, p, q, pdel, &
                       cldt, cldc, cmfdqr, conicw, precs, conds, &
                       evaps, cwat, tracer, deltat, &
                       scavt, iscavt, cldv, fracis, sol_fact, ncol, &
                       scavcoef,icscavt, isscavt, bcscavt, bsscavt, &
                       sol_facti_in, sol_factbi_in, sol_factii_in, &
                       sol_factic_in, sol_factiic_in )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging code for very soluble aerosols
      ! 
      ! Author: P. Rasch
      ! Modified by T. Bond 3/2003 to track different removals
      !-----------------------------------------------------------------------

      implicit none

      real(r8), intent(in) ::&
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cldt(pcols,pver),    &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cmfdqr(pcols,pver),   &! rate of production of convective precip
         conicw(pcols,pver),   &! convective cloud water
         cwat(pcols,pver),     &! cloud water amount 
         precs(pcols,pver),    &! rate of production of stratiform precip
         conds(pcols,pver),    &! rate of production of condensate
         evaps(pcols,pver),    &! rate of evaporation of precip
         cldv(pcols,pver),     &! total cloud fraction
         deltat,               &! time step
         tracer(pcols,pver)     ! trace species
      ! If subroutine is called with just sol_fact:
            ! sol_fact is used for both in- and below-cloud scavenging
      ! If subroutine is called with optional argument sol_facti_in:
            ! sol_fact  is used for below cloud scavenging
            ! sol_facti is used for in cloud scavenging
         real(r8), intent(in) :: sol_fact ! solubility factor (fraction of aer scavenged below & in, or just below or sol_facti_in is provided)
         real(r8), intent(in), optional :: sol_facti_in   ! solubility factor (frac of aerosol scavenged in cloud)
         real(r8), intent(in), optional :: sol_factbi_in  ! solubility factor (frac of aerosol scavenged below cloud by ice)
         real(r8), intent(in), optional :: sol_factii_in  ! solubility factor (frac of aerosol scavenged in cloud by ice)
         real(r8), intent(in), optional :: sol_factic_in(pcols,pver)  ! sol_facti_in for convective clouds
         real(r8), intent(in), optional :: sol_factiic_in ! sol_factii_in for convective clouds
         real(r8), intent(in) :: scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
         
      integer, intent(in) :: ncol

      real(r8), intent(out) ::&
         scavt(pcols,pver),    &! scavenging tend 
         iscavt(pcols,pver),   &! incloud scavenging tends
         fracis(pcols,pver)     ! fraction of species not scavenged

      real(r8), intent(out), optional ::    icscavt(pcols,pver)     ! incloud, convective
      real(r8), intent(out), optional ::    isscavt(pcols,pver)     ! incloud, stratiform
      real(r8), intent(out), optional ::    bcscavt(pcols,pver)     ! below cloud, convective
      real(r8), intent(out), optional ::    bsscavt(pcols,pver)     ! below cloud, stratiform

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) fracev(pcols)        ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precabc(pcols)       ! conv precip from above (work array)
      real(r8) precabs(pcols)       ! strat precip from above (work array)
      real(r8) precbl               ! precip falling out of level (work array)
      real(r8) precmin              ! minimum convective precip causing scavenging
      real(r8) rat(pcols)           ! ratio of amount available to amount removed
      real(r8) scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(pcols)       ! scavenged tracer flux from above (work array)
      real(r8) srcc                 ! tend for convective rain
      real(r8) srcs                 ! tend for stratiform rain
      real(r8) srct(pcols)          ! work variable
      real(r8) tracab(pcols)        ! column integrated tracer amount
!      real(r8) vfall                ! fall speed of precip
      real(r8) fins                 ! fraction of rem. rate by strat rain
      real(r8) finc                 ! fraction of rem. rate by conv. rain
      real(r8) srcs1                ! work variable
      real(r8) srcs2                ! work variable
      real(r8) tc                   ! temp in celcius
      real(r8) weight               ! fraction of condensate which is ice
      real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
      real(r8) odds                 ! limit on removal rate (proportional to prec)
      real(r8) dblchek(pcols)
      logical :: found

      real(r8) sol_facti,  sol_factb  ! in cloud and below cloud fraction of aerosol scavenged
      real(r8) sol_factii, sol_factbi ! in cloud and below cloud fraction of aerosol scavenged by ice
      real(r8) sol_factic(pcols,pver)             ! sol_facti for convective clouds
      real(r8) sol_factiic            ! sol_factii for convective clouds
      ! sol_factic & solfact_iic added for MODAL_AERO.  
      ! For stratiform cloud, cloudborne aerosol is treated explicitly,
      !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
      ! For convective cloud, cloudborne aerosol is not treated explicitly,
      !    and sol_factic is 1.0 for both cloudborne and interstitial.

      ! ------------------------------------------------------------------------
!      omsm = 1.-1.e-10          ! used to prevent roundoff errors below zero
      omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero
      precmin =  0.1_r8/8.64e4_r8      ! set critical value to 0.1 mm/day in kg/m2/s

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
!      vfall = 4.
	
      ! default (if other sol_facts aren't in call, set all to required sol_fact
      sol_facti = sol_fact
      sol_factb = sol_fact
      sol_factii = sol_fact
      sol_factbi = sol_fact

      if ( present(sol_facti_in) )  sol_facti = sol_facti_in
      if ( present(sol_factii_in) )  sol_factii = sol_factii_in
      if ( present(sol_factbi_in) )  sol_factbi = sol_factbi_in

      sol_factic  = sol_facti
      sol_factiic = sol_factii
      if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in
      if ( present(sol_factiic_in) )  sol_factiic = sol_factiic_in

      ! this section of code is for highly soluble aerosols,
      ! the assumption is that within the cloud that
      ! all the tracer is in the cloud water
      !
      ! for both convective and stratiform clouds, 
      ! the fraction of cloud water converted to precip defines
      ! the amount of tracer which is pulled out.
      !

      do i = 1,pcols
         precabs(i) = 0
         precabc(i) = 0
         scavab(i) = 0
         scavabc(i) = 0
         tracab(i) = 0
         cldmabs(i) = 0
         cldmabc(i) = 0
      end do

      do k = 1,pver
         do i = 1,ncol
            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
            weight = 0._r8                                 ! assume no ice

            pdog = pdel(i,k)/gravit

            ! ****************** Evaporation **************************
            ! calculate the fraction of strat precip from above 
            !                 which evaporates within this layer
            fracev(i) = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))

            ! trap to ensure reasonable ratio bounds
            fracev(i) = max(0._r8,min(1._r8,fracev(i)))

            ! ****************** Convection ***************************
            ! now do the convective scavenging

            ! set odds proportional to fraction of the grid box that is swept by the 
            ! precipitation =precabc/rhoh20*(area of sphere projected on plane
            !                                /volume of sphere)*deltat
            ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
            ! unless the fraction of the area that is cloud is less than odds, in which
            ! case use the cloud fraction (assumes precabs is in kg/m2/s)
            ! is really: precabs*3/4/1000./1e-3*deltat
            ! here I use .1 from Balkanski
            !
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/gravit*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1.-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1.-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1.-weight) &
            !         /deltat 

            ! fraction of convective cloud water converted to rain
            fracp = cmfdqr(i,k)*deltat/max(1.e-8_r8,conicw(i,k))
            ! note cmfdrq can be negative from evap of rain, so constrain it
            fracp = max(min(1._r8,fracp),0._r8)
            ! remove that amount from within the convective area
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat ! liquid only
!           srcs1 = cldc(i,k)*fracp*tracer(i,k)/deltat             ! any condensation
!           srcs1 = 0.
            srcs1 = sol_factic(i,k)*cldt(i,k)*fracp*tracer(i,k)*(1._r8-weight)/deltat &  ! liquid
                 +  sol_factiic*cldt(i,k)*fracp*tracer(i,k)*(weight)/deltat      ! ice


            !--mcb

            ! scavenge below cloud

            !            cldmabc(i) = max(cldc(i,k),cldmabc(i))
            !            cldmabc(i) = max(cldt(i,k),cldmabc(i))
            cldmabc(i) = max(cldv(i,k),cldmabc(i))
            cldmabc(i) = cldv(i,k)

            odds=max( &
                 min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8) &
                 *scavcoef(i,k)*deltat),0._r8) ! Dana and Hales coefficient (/mm)
            srcs2 = sol_factb*cldmabc(i)*odds*tracer(i,k)*(1._r8-weight)/deltat & ! liquid
                 +  sol_factbi*cldmabc(i)*odds*tracer(i,k)*(weight)/deltat    !ice
            !Note that using the temperature-determined weight doesn't make much sense here


            srcc = srcs1 + srcs2  ! convective tend by both processes
            finc = srcs1/(srcc + 1.e-36_r8) ! fraction in-cloud

            ! ****************** Stratiform ***********************
            ! now do the stratiform scavenging

            ! incloud scavenging

            ! fracp is the fraction of cloud water converted to precip
            fracp =  precs(i,k)*deltat/max(cwat(i,k),1.e-12_r8)
            fracp = max(0._r8,min(1._r8,fracp))
!           fracp = 0.     ! for debug

            ! assume the corresponding amnt of tracer is removed
            !++mcb -- remove cldc; change cldt to cldv 
            !            srcs1 = (cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat
            !            srcs1 = cldv(i,k)*fracp*tracer(i,k)/deltat &
!           srcs1 = cldt(i,k)*fracp*tracer(i,k)/deltat            ! all condensate
            srcs1 = sol_facti*cldt(i,k)*fracp*tracer(i,k)/deltat*(1._r8-weight) &  ! liquid
                 + sol_factii*cldt(i,k)*fracp*tracer(i,k)/deltat*(weight)       ! ice


            ! below cloud scavenging

!           volume undergoing below cloud scavenging
            cldmabs(i) = cldv(i,k)   ! precipitating volume
!           cldmabs(i) = cldt(i,k)   ! local cloud volume

            odds = precabs(i)/max(cldmabs(i),1.e-5_r8)*scavcoef(i,k)*deltat
            odds = max(min(1._r8,odds),0._r8)
            srcs2 =sol_factb*(cldmabs(i)*odds) *tracer(i,k)*(1._r8-weight)/deltat & ! liquid
                 + sol_factbi*(cldmabs(i)*odds) *tracer(i,k)*(weight)/deltat       ! ice
            !Note that using the temperature-determined weight doesn't make much sense here


            srcs = srcs1 + srcs2             ! total stratiform scavenging
            fins=srcs1/(srcs + 1.e-36_r8)    ! fraction taken by incloud processes

            ! make sure we dont take out more than is there
            ! ratio of amount available to amount removed
            rat(i) = tracer(i,k)/max(deltat*(srcc+srcs),1.e-36_r8)
            if (rat(i).lt.1._r8) then
               srcs = srcs*rat(i)
               srcc = srcc*rat(i)
            endif
            srct(i) = (srcc+srcs)*omsm

            
            ! fraction that is not removed within the cloud
            ! (assumed to be interstitial, and subject to convective transport)
            fracp = deltat*srct(i)/max(cldmabs(i)*tracer(i,k),1.e-36_r8)  ! amount removed
            fracp = max(0._r8,min(1._r8,fracp))
            fracis(i,k) = 1._r8 - fracp

            ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            scavt(i,k) = -srct(i) + fracev(i)*scavab(i)*gravit/pdel(i,k)
            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

            if ( present(icscavt) ) icscavt(i,k) = -(srcc*finc) * omsm
            if ( present(isscavt) ) isscavt(i,k) = -(srcs*fins) * omsm
            if ( present(bcscavt) ) bcscavt(i,k) = -(srcc * (1-finc)) * omsm
            if ( present(bsscavt) ) bsscavt(i,k) = -(srcs * (1-fins)) * omsm +  &
                 fracev(i)*scavab(i)*gravit/pdel(i,k)

            dblchek(i) = tracer(i,k) + deltat*scavt(i,k)

            ! now keep track of scavenged mass and precip
            scavab(i) = scavab(i)*(1-fracev(i)) + srcs*pdel(i,k)/gravit
            precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
            scavabc(i) = scavabc(i) + srcc*pdel(i,k)/gravit
            precabc(i) = precabc(i) + (cmfdqr(i,k))*pdel(i,k)/gravit
            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit

         end do

         found = .false.
         do i = 1,ncol
            if ( dblchek(i) < 0._r8 ) then
               found = .true.
               exit
            end if
         end do

         if ( found ) then
            do i = 1,ncol
               if (dblchek(i) .lt. 0._r8) then
                  write(iulog,*) ' wetdapa: negative value ', i, k, tracer(i,k), &
                       dblchek(i), scavt(i,k), srct(i), rat(i), fracev(i)
               endif
            end do
         endif

      end do

   end subroutine wetdepa_v1

!==============================================================================

! wetdepg is currently being used for both CAM4 and CAM5 by making use of the
! cam_physpkg_is method.

   subroutine wetdepg( t, p, q, pdel, &
                       cldt, cldc, cmfdqr, evapc, precs, evaps, &
                       rain, cwat, tracer, deltat, molwt, &
                       solconst, scavt, iscavt, cldv, icwmr1, &
                       icwmr2, fracis, ncol )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! scavenging of gas phase constituents by henry's law
      ! 
      ! Author: P. Rasch
      !-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cmfdqr(pcols,pver),   &! rate of production of convective precip
         rain (pcols,pver),    &! total rainwater mixing ratio
         cwat(pcols,pver),     &! cloud water amount 
         precs(pcols,pver),    &! rate of production of stratiform precip
         evaps(pcols,pver),    &! rate of evaporation of precip
! Sungsu
         evapc(pcols,pver),    &! Rate of evaporation of convective precipitation
! Sungsu 
         cldv(pcols,pver),     &! estimate of local volume occupied by clouds
         icwmr1 (pcols,pver),  &! in cloud water mixing ration for zhang scheme
         icwmr2 (pcols,pver),  &! in cloud water mixing ration for hack  scheme
         deltat,               &! time step
         tracer(pcols,pver),   &! trace species
         molwt                  ! molecular weights

      integer, intent(in) :: ncol

      real(r8) &
         solconst(pcols,pver)   ! Henry's law coefficient

      real(r8), intent(out) ::&
         scavt(pcols,pver),    &! scavenging tend 
         iscavt(pcols,pver),   &! incloud scavenging tends
         fracis(pcols, pver)    ! fraction of constituent that is insoluble

      ! local variables

      integer i                 ! x index
      integer k                 ! z index

      real(r8) adjfac               ! factor stolen from cmfmca
      real(r8) aqfrac               ! fraction of tracer in aqueous phase
      real(r8) cwatc                ! local convective total water amount 
      real(r8) cwats                ! local stratiform total water amount 
      real(r8) cwatl                ! local cloud liq water amount 
      real(r8) cwatp                ! local water amount falling from above precip
      real(r8) cwatpl               ! local water amount falling from above precip (liq)
      real(r8) cwatt                ! local sum of strat + conv total water amount 
      real(r8) cwatti               ! cwatt/cldv = cloudy grid volume mixing ratio
      real(r8) fracev               ! fraction of precip from above that is evaporating
      real(r8) fracp                ! fraction of cloud water converted to precip
      real(r8) gafrac               ! fraction of tracer in gas phasea
      real(r8) hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) mpla                 ! moles / liter H2O entering the layer from above
      real(r8) mplb                 ! moles / liter H2O leaving the layer below
      real(r8) omsm                 ! 1 - (a small number)
      real(r8) part                 !  partial pressure of tracer in atmospheres
      real(r8) patm                 ! total pressure in atmospheres
      real(r8) pdog                 ! work variable (pdel/gravit)
      real(r8) precab(pcols)        ! precip from above (work array)
      real(r8) precbl               ! precip work variable
      real(r8) precxx               ! precip work variable
      real(r8) precxx2               !
      real(r8) precic               ! precip work variable
      real(r8) rat                  ! ratio of amount available to amount removed
      real(r8) scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) scavabc(pcols)       ! scavenged tracer flux from above (work array)
      !      real(r8) vfall                ! fall speed of precip
      real(r8) scavmax              ! an estimate of the max tracer avail for removal
      real(r8) scavbl               ! flux removed at bottom of layer
      real(r8) fins                 ! in cloud fraction removed by strat rain
      real(r8) finc                 ! in cloud fraction removed by conv rain
      real(r8) rate                 ! max removal rate estimate
      real(r8) scavlimt             ! limiting value 1
      real(r8) scavt1               ! limiting value 2
      real(r8) scavin               ! scavenging by incloud processes
      real(r8) scavbc               ! scavenging by below cloud processes
      real(r8) tc
      real(r8) weight               ! ice fraction
      real(r8) wtpl                 ! work variable
      real(r8) cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) cldmabc(pcols)       ! maximum cloud at or above this level
      !-----------------------------------------------------------

      omsm = 1._r8-2*epsilon(1._r8)   ! used to prevent roundoff errors below zero

      adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

      ! assume 4 m/s fall speed currently (should be improved)
      !      vfall = 4.

      ! zero accumulators
      do i = 1,pcols
         precab(i) = 1.e-36_r8
         scavab(i) = 0._r8
         cldmabs(i) = 0._r8
      end do

      do k = 1,pver
         do i = 1,ncol

            tc     = t(i,k) - tmelt
            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

            cldmabs(i) = max(cldmabs(i),cldt(i,k))

            ! partitioning coefs for gas and aqueous phase
            !              take as a cloud water amount, the sum of the stratiform amount
            !              plus the convective rain water amount 

            ! convective amnt is just the local precip rate from the hack scheme
            !              since there is no storage of water, this ignores that falling from above
            !            cwatc = cmfdqr(i,k)*deltat/adjfac
            !++mcb -- test cwatc
            cwatc = (icwmr1(i,k) + icwmr2(i,k)) * (1._r8-weight)
            !--mcb 

            ! strat cloud water amount and also ignore the part falling from above
            cwats = cwat(i,k) 

            ! cloud water as liq
            !++mcb -- add cwatc later (in cwatti)
            !            cwatl = (1.-weight)*(cwatc+cwats)
            cwatl = (1._r8-weight)*cwats
            ! cloud water as ice
            !*not used        cwati = weight*(cwatc+cwats)

            ! total suspended condensate as liquid
            cwatt = cwatl + rain(i,k)

            ! incloud version 
            !++mcb -- add cwatc here
            cwatti = cwatt/max(cldv(i,k), 0.00001_r8) + cwatc

            ! partitioning terms
            patm = p(i,k)/1.013e5_r8 ! pressure in atmospheres
            hconst = molwta*patm*solconst(i,k)*cwatti/rhoh2o
            aqfrac = hconst/(1._r8+hconst)
            gafrac = 1/(1._r8+hconst)
            fracis(i,k) = gafrac


            ! partial pressure of the tracer in the gridbox in atmospheres
            part = patm*gafrac*tracer(i,k)*molwta/molwt

            ! use henrys law to give moles tracer /liter of water
            ! in this volume 
            ! then convert to kg tracer /liter of water (kg tracer / kg water)
            mplb = solconst(i,k)*part*molwt/1000._r8


            pdog = pdel(i,k)/gravit

            ! this part of precip will be carried downward but at a new molarity of mpl 
            precic = pdog*(precs(i,k) + cmfdqr(i,k))

            ! we cant take out more than entered, plus that available in the cloud
            !                  scavmax = scavab(i)+tracer(i,k)*cldt(i,k)/deltat*pdog
            scavmax = scavab(i)+tracer(i,k)*cldv(i,k)/deltat*pdog

            ! flux of tracer by incloud processes
            scavin = precic*(1._r8-weight)*mplb

            ! fraction of precip which entered above that leaves below
            if (cam_physpkg_is('cam5')) then
               ! Sungsu added evaporation of convective precipitation below.
               precxx = precab(i)-pdog*(evaps(i,k)+evapc(i,k))
            else
               precxx = precab(i)-pdog*evaps(i,k)
            end if
            precxx = max (precxx,0.0_r8)

            ! flux of tracer by below cloud processes
            !++mcb -- removed wtpl because it is now not assigned and previously
            !          when it was assigned it was unnecessary:  if(tc.gt.0)wtpl=1
            if (tc.gt.0) then
               !               scavbc = precxx*wtpl*mplb ! if liquid
               scavbc = precxx*mplb ! if liquid
            else
               precxx2=max(precxx,1.e-36_r8)
               scavbc = scavab(i)*precxx2/(precab(i)) ! if ice
            endif

            scavbl = min(scavbc + scavin, scavmax)

            ! first guess assuming that henries law works
            scavt1 = (scavab(i)-scavbl)/pdog*omsm

            ! pjr this should not be required, but we put it in to make sure we cant remove too much
            ! remember, scavt1 is generally negative (indicating removal)
            scavt1 = max(scavt1,-tracer(i,k)*cldv(i,k)/deltat)

            !++mcb -- remove this limitation for gas species
            !c use the dana and hales or balkanski limit on scavenging
            !c            rate = precab(i)*0.1
            !            rate = (precic + precxx)*0.1
            !            scavlimt = -tracer(i,k)*cldv(i,k)
            !     $           *rate/(1.+rate*deltat)

            !            scavt(i,k) = max(scavt1, scavlimt)

            ! instead just set scavt to scavt1
            scavt(i,k) = scavt1
            !--mcb

            ! now update the amount leaving the layer
            scavbl = scavab(i) - scavt(i,k)*pdog 

            ! in cloud amount is that formed locally over the total flux out bottom
            fins = scavin/(scavin + scavbc + 1.e-36_r8)
            iscavt(i,k) = scavt(i,k)*fins

            scavab(i) = scavbl
            precab(i) = max(precxx + precic,1.e-36_r8)

        
            
         end do
      end do
      
   end subroutine wetdepg

!##############################################################################

end module wetdep
