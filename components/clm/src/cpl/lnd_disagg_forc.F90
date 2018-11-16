module lnd_disagg_forc

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Disaggregates gridcell quantities to topounits
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_megan_mod  , only : shr_megan_mechcomps_n
  use clm_varpar     , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon     , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl     , only : iulog, use_c13, use_cn, use_lch4, iulog
  use clm_cpl_indices
  use seq_drydep_mod , only : n_drydep, drydep_method, DD_XLND
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use atm2lndType    , only : atm2lnd_type
  use lnd2atmType    , only : lnd2atm_type
  use glc2lndMod     , only : glc2lnd_type
  use GridcellType   , only : grc_pp
  use TopounitType   , only : top_pp, top_as, top_es        
  use LandunitType   , only : lun_pp                
  use ColumnType     , only : col_pp                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: downscale_atmo_state_to_topounit    ! Downscale atm state forcing fields from gridcell to column
  !
  ! !PRIVATE MEMBER FUNCTIONS:
!  private :: downscale_longwave_top      ! Downscale longwave radiation from gridcell to column
!  private :: build_normalization         ! Compute normalization factors so that downscaled fields are conservative
!  private :: check_downscale_consistency ! Check consistency of downscaling
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine downscale_atmo_state_to_topounit(g, i, x2l)
    !
    ! !DESCRIPTION:
    ! Downscale atmospheric forcing fields from gridcell to topounit
    !
    ! Downscaling is done over topounits if # topounits > 1.
    !
    ! !USES:
    use clm_time_manager, only : get_nstep
    use clm_varcon      , only : rair, cpair, grav, lapse_glcmec
    use clm_varcon      , only : glcmec_rain_snow_threshold
    use landunit_varcon , only : istice_mec 
    use clm_varctl      , only : glcmec_downscale_rain_snow_convert
    use domainMod       , only : ldomain
    use QsatMod         , only : Qsat
    !
    ! !ARGUMENTS:
    type(integer)                    , intent(in)    :: g  
    type(integer)                    , intent(in)    :: i
    real(r8)                         , intent(in)    :: x2l(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: t, l, c, fc         ! indices
    integer :: clo, cc

    ! temporaries for topo downscaling
    real(r8) :: hsurf_g,hsurf_t,Hbot
    real(r8) :: zbot_g, tsfc_g, tbot_g, pbot_g, thbot_g, qbot_g, qs_g, es_g
    real(r8) :: zbot_t, tsfc_t, tbot_t, pbot_t, thbot_t, qbot_t, qs_t, es_t
    real(r8) :: egcm_t, rhos_t
    real(r8) :: dum1,   dum2

!    real(r8), dimension(bounds%begg : bounds%endg) :: sum_qbot_g    ! weighted sum of column-level lwrad
!    real(r8), dimension(bounds%begg : bounds%endg) :: sum_wts_g      ! sum of weights that contribute to sum_lwrad_g
!    real(r8), dimension(bounds%begg : bounds%endg) :: qbot_norm_g   ! normalization factors

    logical  :: do_lapse_downscaling = .true.

    character(len=*), parameter :: subname = 'downscale_forcings'
    !-----------------------------------------------------------------------

!    nstep = get_nstep()
    
    ! Downscale forc_t, forc_th, forc_q, forc_pbot, and forc_rho to columns.
    ! For glacier_mec columns the downscaling is based on surface elevation.
    ! For other columns the downscaling is a simple copy (above).

    do t = grc_pp%topi(g), grc_pp%topf(g)

      top_as%tbot(t)    = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
      top_as%thbot(t)   = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
      top_as%pbot(t)    = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy    Atm state Pa
      top_as%qbot(t)    = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
      top_as%ubot(t)    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
      top_as%vbot(t)    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
      top_as%zbot(t)    = x2l(index_x2l_Sa_z,i)         ! zgcmxy    Atm state m

         ! This is a simple downscaling procedure 
         ! Note that forc_hgt, forc_u, and forc_v are not downscaled.

      hsurf_g = ldomain%topo(g)                       ! gridcell sfc elevation
!         hsurf_t = top_pp%elevation(t)                  ! topounit sfc elevation
      hsurf_t = ldomain%topo(g)                       ! topounit sfc elevation
      tbot_g  = x2l(index_x2l_Sa_tbot,i)              ! atm sfc temp
      thbot_g = x2l(index_x2l_Sa_ptem,i)              ! atm sfc pot temp
      qbot_g  = x2l(index_x2l_Sa_shum,i)              ! atm sfc spec humid
      pbot_g  = x2l(index_x2l_Sa_pbot,i)              ! atm sfc pressure
      zbot_g  = x2l(index_x2l_Sa_z,i)                 ! atm ref height

         zbot_t  = zbot_g
         tbot_t  = tbot_g-lapse_glcmec*(hsurf_t-hsurf_g) ! sfc temp for column
 
         Hbot    = rair*0.5_r8*(tbot_g+tbot_t)/grav      ! scale ht at avg temp
         pbot_t  = pbot_g*exp(-(hsurf_t-hsurf_g)/Hbot)   ! column sfc press

         ! Derivation of potential temperature calculation:
         ! 
         ! The textbook definition would be:
         ! thbot_c = tbot_c * (p0/pbot_c)^(rair/cpair)
         ! 
         ! Note that pressure is related to scale height as:
         ! pbot_c = p0 * exp(-zbot_c/H)
         !
         ! Using Hbot in place of H, we get:
         ! pbot_c = p0 * exp(-zbot_c/Hbot)
         !
         ! Plugging this in to the textbook definition, then manipulating, we get:
         ! thbot_c = tbot_c * (p0/(p0*exp(-zbot_c/Hbot)))^(rair/cpair)
         !         = tbot_c * (1/exp(-zbot_c/Hbot))^(rair/cpair)
         !         = tbot_c * (exp(zbot_c/Hbot))^(rair/cpair)
         !         = tbot_c * exp((zbot_c/Hbot) * (rair/cpair))

         thbot_t= tbot_t*exp((zbot_t/Hbot)*(rair/cpair))  ! pot temp calc

         call Qsat(tbot_g,pbot_g,es_g,dum1,qs_g,dum2)
         call Qsat(tbot_t,pbot_t,es_t,dum1,qs_t,dum2)

         qbot_t = qbot_g*(qs_t/qs_g)
         egcm_t = qbot_t*pbot_t/(0.622+0.378*qbot_t)
         rhos_t = (pbot_t-0.378*egcm_t) / (rair*tbot_t)

       top_as%tbot(t) = tbot_t
       top_as%thbot(t) = thbot_t
!	 top_as%vp_atm(t) = es_t
       top_as%qbot(t) = qbot_t
       top_as%pbot(t) = pbot_t

     end do

!      call downscale_longwave_top(bounds, atm2lnd_vars, top_pp, top_es)

!      call check_downscale_consistency(bounds, atm2lnd_vars)

  end subroutine downscale_atmo_state_to_topounit

end module lnd_disagg_forc
