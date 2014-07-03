module BalanceCheckMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BalanceCheckMod
!
! !DESCRIPTION:
! Water and energy balance check.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
  use clm_varctl,   only: iulog

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginWaterBalance  ! Initialize water balance check
  public :: BalanceCheck       ! Water and energy balance check
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginWaterBalance
!
! !INTERFACE:
  subroutine BeginWaterBalance(lbc, ubc, lbp, ubp, &
             num_nolakec, filter_nolakec, num_lakec, filter_lakec, &
             num_hydrologyc, filter_hydrologyc)
!
! !DESCRIPTION:
! Initialize column-level water balance at beginning of time step
!
! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use clm_varpar   , only : nlevgrnd, nlevsoi, nlevurb
    use subgridAveMod, only : p2c
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, &
                              icol_road_imperv
    use clm_varcon   , only : denh2o, denice
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_lakec                   ! number of column non-lake points in column filter
    integer, intent(in) :: filter_lakec(ubc-lbc+1)     ! column filter for non-lake points
    integer , intent(in)  :: num_hydrologyc               ! number of column soil points in column filter
    integer , intent(in)  :: filter_hydrologyc(ubc-lbc+1) ! column filter for soil points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! Created by Peter Thornton
!
!EOP
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in variables
!
    real(r8), pointer :: h2osfc(:)             ! surface water (mm)
    real(r8), pointer :: londeg(:)             ! longitude
    real(r8), pointer :: latdeg(:)             ! latitude
    integer , pointer :: cgridcell(:)          ! column's gridcell index
    integer , pointer :: clandunit(:)          ! column's landunit
    integer , pointer :: ltype(:)              ! landunit type
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    real(r8), pointer :: h2osoi_ice(:,:)       ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:)       ! liquid water (kg/m2)
    real(r8), pointer :: h2ocan_pft(:)         ! canopy water (mm H2O) (pft-level) 
    real(r8), pointer :: wa(:)                 ! water in the unconfined aquifer (mm)
    integer , pointer :: ctype(:)              ! column type 
    real(r8), pointer :: zwt(:)                ! water table depth (m)
    real(r8), pointer :: zi(:,:)               ! interface level below a "z" level (m)
!
! local pointers to original implicit out variables
!
    real(r8), pointer :: h2ocan_col(:)         ! canopy water (mm H2O) (column level)
    real(r8), pointer :: begwb(:)              ! water mass begining of the time step
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, p, f, j, fc                  ! indices
    real(r8):: h2osoi_vol
    real(r8), pointer :: dz(:,:), watsat(:,:)
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (column-level)

    h2osfc             => clm3%g%l%c%cws%h2osfc
    londeg             => clm3%g%londeg
    latdeg             => clm3%g%latdeg
    cgridcell          => clm3%g%l%c%gridcell
    clandunit          => clm3%g%l%c%landunit
    ltype              => clm3%g%l%itype
    dz                 => clm3%g%l%c%cps%dz
    watsat             => clm3%g%l%c%cps%watsat
    h2osno             => clm3%g%l%c%cws%h2osno
    h2osoi_ice         => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq         => clm3%g%l%c%cws%h2osoi_liq
    begwb              => clm3%g%l%c%cwbal%begwb
    h2ocan_col         => clm3%g%l%c%cws%pws_a%h2ocan
    wa                 => clm3%g%l%c%cws%wa
    ctype              => clm3%g%l%c%itype
    zwt                => clm3%g%l%c%cws%zwt
    zi                 => clm3%g%l%c%cps%zi

    ! Assign local pointers to derived type members (pft-level)

    h2ocan_pft         => clm3%g%l%c%p%pws%h2ocan

    ! Determine beginning water balance for time step
    ! pft-level canopy water averaged to column
    call p2c(num_nolakec, filter_nolakec, h2ocan_pft, h2ocan_col)

    do f = 1, num_hydrologyc
       c = filter_hydrologyc(f)
       if(zwt(c) <= zi(c,nlevsoi)) then
          wa(c) = 5000._r8
       end if
    end do

    do f = 1, num_nolakec
       c = filter_nolakec(f)
       if (ctype(c) == icol_roof .or. ctype(c) == icol_sunwall &
          .or. ctype(c) == icol_shadewall .or. ctype(c) == icol_road_imperv) then
         begwb(c) = h2ocan_col(c) + h2osno(c)
       else
         begwb(c) = h2ocan_col(c) + h2osno(c) + h2osfc(c) + wa(c)
       end if

    end do
    do j = 1, nlevgrnd
      do f = 1, num_nolakec
         c = filter_nolakec(f)
         if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
          .or. ctype(c) == icol_roof) .and. j > nlevurb) then
         else
            begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
         end if
      end do
    end do

    do f = 1, num_lakec
       c = filter_lakec(f)
       begwb(c) = h2osno(c)
    end do

  end subroutine BeginWaterBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BalanceCheck
!
! !INTERFACE:
  subroutine BalanceCheck(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg)
!
! !DESCRIPTION:
! This subroutine accumulates the numerical truncation errors of the water
! and energy balance calculation. It is helpful to see the performance of
! the process of integration.
!
! The error for energy balance:
!
! error = abs(Net radiation - change of internal energy - Sensible heat
!             - Latent heat)
!
! The error for water balance:
!
! error = abs(precipitation - change of water storage - evaporation - runoff)
!
! !USES:
    use clmtype
    use clm_atmlnd   , only : clm_a2l
    use subgridAveMod
    use clm_time_manager , only : get_step_size, get_nstep
    use clm_varcon   , only : isturb, icol_roof, icol_sunwall, icol_shadewall, &
                              spval, icol_road_perv, icol_road_imperv, istice_mec, &
                              istdlak, istslak,istsoil,istcrop,istwet
    use clm_varctl   , only : glc_dyntopo, create_glacier_mec_landunit
!
! !ARGUMENTS:
    implicit none
    integer :: lbp, ubp ! pft-index bounds
    integer :: lbc, ubc ! column-index bounds
    integer :: lbl, ubl ! landunit-index bounds
    integer :: lbg, ubg ! grid-index bounds
!
! !CALLED FROM:
! subroutine clm_driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 10 November 2000: Mariana Vertenstein
! Migrated to new data structures by Mariana Vertenstein and
! Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arguments
!
    real(r8), pointer :: tws(:)                !total water storage (mm H2O)
    real(r8), pointer :: volr(:)               !river water storage (m3)
    real(r8), pointer :: area(:)               !gridcell area (km2)
    logical , pointer :: do_capsnow(:)         ! true => do snow capping
    real(r8), pointer :: qflx_rain_grnd_col(:) ! rain on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_grnd_col(:) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), pointer :: qflx_snow_h2osfc(:)   ! snow falling on surface water (mm/s)
    real(r8), pointer :: frac_sno_eff(:)       ! effective snow fraction
    real(r8), pointer :: qflx_h2osfc_to_ice(:) ! conversion of h2osfc to ice
    real(r8), pointer :: qflx_snow_melt(:)     ! snow melt (net)
    real(r8), pointer :: frac_sno(:)           ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: qflx_drain_perched(:) ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_floodc(:)        ! total runoff due to flooding
    real(r8), pointer :: qflx_evap_soi(:)      ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_h2osfc_surf(:)!surface water runoff (mm/s)
    real(r8), pointer :: sabg_soil(:)       ! solar radiation absorbed by soil (W/m**2)
    real(r8), pointer :: sabg_snow(:)       ! solar radiation absorbed by snow (W/m**2)
    real(r8), pointer :: sabg_chk(:)        ! sum of soil/snow using current fsno, for balance check
    integer , pointer :: pcolumn(:)         ! pft's column index
    logical , pointer :: pactive(:)         ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: cactive(:)         ! true=>do computations on this column (see reweightMod for details)
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: plandunit(:)       ! pft's landunit index
    integer , pointer :: cgridcell(:)       ! column's gridcell index
    integer , pointer :: clandunit(:)       ! column's landunit index
    integer , pointer :: ltype(:)           ! landunit type 
    integer , pointer :: ctype(:)           ! column type 
    real(r8), pointer :: forc_rain(:)       ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)       ! snow rate [mm/s]
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: endwb(:)           ! water mass end of the time step
    real(r8), pointer :: begwb(:)           ! water mass begining of the time step
    real(r8), pointer :: fsa(:)             ! solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsr(:)             ! solar radiation reflected (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: sabv(:)            ! solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_totg(:)    ! total sensible heat flux at grid level (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_dynbal(:)     ! energy conversion flux due to dynamic land cover change(W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: qflx_irrig(:)      ! irrigation flux (mm H2O /s)
    real(r8), pointer :: qflx_surf(:)       ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)      ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_drain(:)      ! sub-surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_runoff(:)     ! total runoff (mm H2O /s)
    real(r8), pointer :: qflx_runoffg(:)    ! total runoff at gridcell level inc land cover change flux (mm H2O /s)
    real(r8), pointer :: qflx_liq_dynbal(:) ! liq runoff due to dynamic land cover change (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice(:)  ! excess snowfall due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: qflx_glcice(:)     ! flux of new glacier ice (mm H2O /s) [+ if ice grows]
    real(r8), pointer :: qflx_glcice_frz(:) ! ice growth (mm H2O/s) [+]
    real(r8), pointer :: qflx_snwcp_iceg(:) ! excess snowfall due to snow cap inc land cover change flux (mm H20/s)
    real(r8), pointer :: qflx_ice_dynbal(:) ! ice runoff due to dynamic land cover change (mm H2O /s)
    real(r8), pointer :: forc_solad(:,:)    ! direct beam radiation (vis=forc_sols , nir=forc_soll )
    real(r8), pointer :: forc_solai(:,:)    ! diffuse radiation     (vis=forc_solsd, nir=forc_solld)
    real(r8), pointer :: eflx_traffic_pft(:)    ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat_pft(:)  ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: canyon_hwr(:)      ! ratio of building height to street width
    real(r8), pointer :: eflx_heat_from_ac_pft(:) !sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: h2osno(:)             ! snow water (mm H2O)
    real(r8), pointer :: h2osno_old(:)         ! snow water (mm H2O) at previous time step
    real(r8), pointer :: qflx_dew_snow(:)      ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_sub_snow(:)      ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_top_soil(:)      ! net water input into soil from top (mm/s)
    real(r8), pointer :: qflx_dew_grnd(:)      ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: qflx_evap_grnd(:)     ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_prec_grnd(:)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), pointer :: qflx_snwcp_liq(:)     ! excess liquid water due to snow capping (mm H2O /s) [+]`
    real(r8), pointer :: qflx_sl_top_soil(:)   ! liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
    integer , pointer :: snl(:)                ! number of snow layers
!
! local pointers to original implicit out arguments
!
    real(r8), pointer :: errh2o(:)          ! water conservation error (mm H2O)
    real(r8), pointer :: errsol(:)          ! solar radiation conservation error (W/m**2)
    real(r8), pointer :: errlon(:)          ! longwave radiation conservation error (W/m**2)
    real(r8), pointer :: errseb(:)          ! surface energy conservation error (W/m**2)
    real(r8), pointer :: netrad(:)          ! net radiation (positive downward) (W/m**2)
    real(r8), pointer :: errsoi_col(:)      ! column-level soil/lake energy conservation error (W/m**2)
    real(r8), pointer :: snow_sources(:)    ! snow sources (mm H2O /s)
    real(r8), pointer :: snow_sinks(:)      ! snow sinks (mm H2O /s)
    real(r8), pointer :: errh2osno(:)       ! error in h2osno (kg m-2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: p,c,l,g                     ! indices
    real(r8) :: dtime                       ! land model time step (sec)
    integer  :: nstep                       ! time step number
    logical  :: found                       ! flag in search loop
    integer  :: indexp,indexc,indexl,indexg ! index of first found in search loop
    real(r8) :: forc_rain_col(lbc:ubc)      ! column level rain rate [mm/s]
    real(r8) :: forc_snow_col(lbc:ubc)      ! column level snow rate [mm/s]
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (gridcell-level)

    tws                 => clm3%g%tws
    area                => clm3%g%area
    volr                => clm_a2l%volr
    do_capsnow          => clm3%g%l%c%cps%do_capsnow
    qflx_rain_grnd_col  => clm3%g%l%c%cwf%pwf_a%qflx_rain_grnd
    qflx_snow_grnd_col  => clm3%g%l%c%cwf%pwf_a%qflx_snow_grnd
    qflx_snow_h2osfc    => clm3%g%l%c%cwf%qflx_snow_h2osfc
    frac_sno_eff        => clm3%g%l%c%cps%frac_sno_eff
    qflx_h2osfc_to_ice  => clm3%g%l%c%cwf%qflx_h2osfc_to_ice
    frac_sno            => clm3%g%l%c%cps%frac_sno 
    qflx_drain_perched  => clm3%g%l%c%cwf%qflx_drain_perched
    qflx_floodc         => clm3%g%l%c%cwf%qflx_floodc
    qflx_evap_soi       => clm3%g%l%c%cwf%pwf_a%qflx_evap_soi
    qflx_h2osfc_surf    => clm3%g%l%c%cwf%qflx_h2osfc_surf
    qflx_snow_melt      => clm3%g%l%c%cwf%qflx_snow_melt
    sabg_soil           => clm3%g%l%c%p%pef%sabg_soil
    sabg_snow           => clm3%g%l%c%p%pef%sabg_snow
    sabg_chk            => clm3%g%l%c%p%pef%sabg_chk
    pcolumn             => clm3%g%l%c%p%column
    forc_rain           => clm_a2l%forc_rain
    forc_snow           => clm_a2l%forc_snow
    forc_lwrad          => clm_a2l%forc_lwrad
    forc_solad          => clm_a2l%forc_solad
    forc_solai          => clm_a2l%forc_solai

    ! Assign local pointers to derived type scalar members (landunit-level)

    ltype             => clm3%g%l%itype
    canyon_hwr        => clm3%g%l%canyon_hwr

    ! Assign local pointers to derived type scalar members (column-level)

    cactive           => clm3%g%l%c%active
    ctype             => clm3%g%l%c%itype
    cgridcell         => clm3%g%l%c%gridcell
    clandunit         => clm3%g%l%c%landunit
    endwb             => clm3%g%l%c%cwbal%endwb
    begwb             => clm3%g%l%c%cwbal%begwb
    qflx_irrig        => clm3%g%l%c%cwf%qflx_irrig
    qflx_surf         => clm3%g%l%c%cwf%qflx_surf
    qflx_qrgwl        => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_drain        => clm3%g%l%c%cwf%qflx_drain
    qflx_runoff       => clm3%g%l%c%cwf%qflx_runoff
    qflx_snwcp_ice    => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_ice
    qflx_evap_tot     => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    qflx_glcice       => clm3%g%l%c%cwf%qflx_glcice
    qflx_glcice_frz   => clm3%g%l%c%cwf%qflx_glcice_frz
    errh2o            => clm3%g%l%c%cwbal%errh2o
    errsoi_col        => clm3%g%l%c%cebal%errsoi
    h2osno             => clm3%g%l%c%cws%h2osno
    h2osno_old         => clm3%g%l%c%cws%h2osno_old
    qflx_dew_snow      => clm3%g%l%c%cwf%pwf_a%qflx_dew_snow
    qflx_sub_snow      => clm3%g%l%c%cwf%pwf_a%qflx_sub_snow
    qflx_top_soil      => clm3%g%l%c%cwf%qflx_top_soil
    qflx_evap_grnd     => clm3%g%l%c%cwf%pwf_a%qflx_evap_grnd
    qflx_dew_grnd      => clm3%g%l%c%cwf%pwf_a%qflx_dew_grnd
    qflx_prec_grnd     => clm3%g%l%c%cwf%pwf_a%qflx_prec_grnd
    qflx_snwcp_liq     => clm3%g%l%c%cwf%pwf_a%qflx_snwcp_liq
    qflx_sl_top_soil   => clm3%g%l%c%cwf%qflx_sl_top_soil
    snow_sources       => clm3%g%l%c%cws%snow_sources
    snow_sinks         => clm3%g%l%c%cws%snow_sinks
    errh2osno          => clm3%g%l%c%cws%errh2osno
    snl                => clm3%g%l%c%cps%snl

    ! Assign local pointers to derived type scalar members (pft-level)

    pactive           => clm3%g%l%c%p%active
    pgridcell         => clm3%g%l%c%p%gridcell
    plandunit         => clm3%g%l%c%p%landunit
    fsa               => clm3%g%l%c%p%pef%fsa
    fsr               => clm3%g%l%c%p%pef%fsr
    eflx_lwrad_out    => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net    => clm3%g%l%c%p%pef%eflx_lwrad_net
    sabv              => clm3%g%l%c%p%pef%sabv
    sabg              => clm3%g%l%c%p%pef%sabg
    eflx_sh_tot       => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lh_tot       => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_soil_grnd    => clm3%g%l%c%p%pef%eflx_soil_grnd
    errsol            => clm3%g%l%c%p%pebal%errsol
    errseb            => clm3%g%l%c%p%pebal%errseb
    errlon            => clm3%g%l%c%p%pebal%errlon
    netrad            => clm3%g%l%c%p%pef%netrad
    eflx_wasteheat_pft => clm3%g%l%c%p%pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => clm3%g%l%c%p%pef%eflx_heat_from_ac_pft
    eflx_traffic_pft  => clm3%g%l%c%p%pef%eflx_traffic_pft

    ! Assign local pointers to derived type scalar members (gridcell-level)

    qflx_runoffg       => clm3%g%gwf%qflx_runoffg
    qflx_liq_dynbal    => clm3%g%gwf%qflx_liq_dynbal
    qflx_snwcp_iceg    => clm3%g%gwf%qflx_snwcp_iceg
    qflx_ice_dynbal    => clm3%g%gwf%qflx_ice_dynbal
    eflx_sh_totg       => clm3%g%gef%eflx_sh_totg
    eflx_dynbal        => clm3%g%gef%eflx_dynbal

    ! Get step size and time step

    nstep = get_nstep()
    dtime = get_step_size()

    ! Determine column level incoming snow and rain
    ! Assume no incident precipitation on urban wall columns (as in Hydrology1Mod.F90).

    do c = lbc,ubc
       g = cgridcell(c)
       if (ctype(c) == icol_sunwall .or.  ctype(c) == icol_shadewall) then
          forc_rain_col(c) = 0.
          forc_snow_col(c) = 0.
       else
          forc_rain_col(c) = forc_rain(g)
          forc_snow_col(c) = forc_snow(g)
       end if
    end do

    ! Water balance check

    do c = lbc, ubc
       g = cgridcell(c)
       l = clandunit(c)
      
       ! add qflx_drain_perched and qflx_flood
       if (cactive(c))then
          errh2o(c) = endwb(c) - begwb(c) &
               - (forc_rain_col(c) + forc_snow_col(c)  + qflx_floodc(c) + qflx_irrig(c) &
                 - qflx_evap_tot(c) - qflx_surf(c)  - qflx_h2osfc_surf(c) &
                 - qflx_qrgwl(c) - qflx_drain(c) - qflx_drain_perched(c) - qflx_snwcp_ice(c)) * dtime

          ! Suppose glc_dyntopo = T:   
          ! (1) We have qflx_snwcp_ice = 0, and excess snow has been incorporated in qflx_glcice.  
          !     This flux must be included here to complete the water balance.
          ! (2) Meltwater from ice is allowed to run off and is included in qflx_qrgwl,
          !     but the water content of the ice column has not changed (at least for now) because
          !     an equivalent ice mass has been "borrowed" from the base of the column.  That
          !     meltwater is included in qflx_glcice.
          !
          ! Note that qflx_glcice is only valid over ice_mec landunits; elsewhere it is spval

          if (glc_dyntopo .and. ltype(l)==istice_mec) then
             errh2o(c) = errh2o(c) + qflx_glcice(c)*dtime
          end if

       else

          errh2o(c) = 0.0_r8

       end if

    end do

    found = .false.
    do c = lbc, ubc
       if (abs(errh2o(c)) > 1e-7_r8) then
          found = .true.
          indexc = c
       end if
    end do

    if ( found ) then
       write(iulog,*)'WARNING:  water balance error ',&
            ' nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
       if ((ctype(indexc) .eq. icol_roof .or. ctype(indexc) .eq. icol_road_imperv .or. &
            ctype(indexc) .eq. icol_road_perv) .and. abs(errh2o(indexc)) > 1.e-1 .and. (nstep > 2) ) then
          write(iulog,*)'clm urban model is stopping - error is greater than 1.e-1'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'forc_rain    = ',forc_rain_col(indexc)
          write(iulog,*)'forc_snow    = ',forc_snow_col(indexc)
          write(iulog,*)'endwb        = ',endwb(indexc)
          write(iulog,*)'begwb        = ',begwb(indexc)
          write(iulog,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
          write(iulog,*)'qflx_irrig   = ',qflx_irrig(indexc)
          write(iulog,*)'qflx_surf    = ',qflx_surf(indexc)
          write(iulog,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
          write(iulog,*)'qflx_drain   = ',qflx_drain(indexc)
          write(iulog,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun()
       else if (abs(errh2o(indexc)) > .10_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping - error is greater than .10'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2o= ',errh2o(indexc)
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'forc_rain    = ',forc_rain_col(indexc)
          write(iulog,*)'forc_snow    = ',forc_snow_col(indexc)
          write(iulog,*)'endwb        = ',endwb(indexc)
          write(iulog,*)'begwb        = ',begwb(indexc)
          write(iulog,*)'qflx_evap_tot= ',qflx_evap_tot(indexc)
          write(iulog,*)'qflx_irrig   = ',qflx_irrig(indexc)
          write(iulog,*)'qflx_surf    = ',qflx_surf(indexc)
          write(iulog,*)'qflx_h2osfc_surf    = ',qflx_h2osfc_surf(indexc)
          write(iulog,*)'qflx_qrgwl   = ',qflx_qrgwl(indexc)
          write(iulog,*)'qflx_drain   = ',qflx_drain(indexc)
          write(iulog,*)'qflx_drain_perched   = ',qflx_drain_perched(indexc)
          write(iulog,*)'qflx_flood   = ',qflx_floodc(indexc)
          write(iulog,*)'qflx_snwcp_ice   = ',qflx_snwcp_ice(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Snow balance check

    do c = lbc, ubc
       g = cgridcell(c)
       l = clandunit(c)
          ! As defined here, snow_sources - snow_sinks will equal the change in h2osno at 
          ! any given time step but only if there is at least one snow layer.  h2osno 
          ! also includes snow that is part of the soil column (an initial snow layer is 
          ! only created if h2osno > 10mm).
          if (snl(c) .lt. 0) then
             snow_sources(c) = qflx_prec_grnd(c) + qflx_dew_snow(c) + qflx_dew_grnd(c)
             snow_sinks(c)   = qflx_sub_snow(c) + qflx_evap_grnd(c) + qflx_snow_melt(c) &
                               + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) + qflx_sl_top_soil(c)

             if (ltype(l) == istdlak) then 
                if ( do_capsnow(c) ) then
                   snow_sources(c) = qflx_snow_grnd_col(c) &
                        + frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) ) 
                   
                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                        + (qflx_snwcp_ice(c) + qflx_snwcp_liq(c) - qflx_prec_grnd(c))  &
                        + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
                else
                   snow_sources(c) = qflx_snow_grnd_col(c) &
                        + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                        +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) 
                   
                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c) ) &
                        + qflx_snow_melt(c)  + qflx_sl_top_soil(c)
                endif
             endif

             if (ltype(l) == istsoil .or. ltype(l) == istcrop .or. ltype(l) == istwet ) then
                 if ( do_capsnow(c) ) then
                    snow_sources(c) = frac_sno_eff(c) * (qflx_dew_snow(c) + qflx_dew_grnd(c) ) &
                         + qflx_h2osfc_to_ice(c) + qflx_prec_grnd(c)

                    snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                         + qflx_snwcp_ice(c) + qflx_snwcp_liq(c) &
                         + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                 else
                    snow_sources(c) = (qflx_snow_grnd_col(c) - qflx_snow_h2osfc(c) ) &
                         + frac_sno_eff(c) * (qflx_rain_grnd_col(c) &
                         +  qflx_dew_snow(c) + qflx_dew_grnd(c) ) + qflx_h2osfc_to_ice(c)

                   snow_sinks(c)   = frac_sno_eff(c) * (qflx_sub_snow(c) + qflx_evap_grnd(c)) &
                        + qflx_snow_melt(c) + qflx_sl_top_soil(c)
                endif
             endif

             ! For ice_mec landunits, if glc_dyntopo is true, then qflx_snwcp_ice = 0,
             ! and qflx_glcice_frz instead stores this flux
             if (ltype(l) == istice_mec .and. glc_dyntopo) then
                snow_sinks(c) = snow_sinks(c) + qflx_glcice_frz(c)
             end if

             errh2osno(c) = (h2osno(c) - h2osno_old(c)) - (snow_sources(c) - snow_sinks(c)) * dtime
          else
             snow_sources(c) = 0._r8
             snow_sinks(c) = 0._r8
             errh2osno(c) = 0._r8
          end if
    end do

    found = .false.
    do c = lbc, ubc
       if (cactive(c) .and. abs(errh2osno(c)) > 1.0e-7_r8) then
          found = .true.
          indexc = c
       end if
    end do
    if ( found ) then
       write(iulog,*)'WARNING:  snow balance error ',&
            ' nstep = ',nstep,' indexc= ',indexc,'ltype: ', ltype(clandunit(indexc)),' errh2osno= ',errh2osno(indexc)
       if (abs(errh2osno(indexc)) > 0.1_r8 .and. (nstep > 2) ) then
          write(iulog,*)'clm model is stopping - error is greater than .10'
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errh2osno= ',errh2osno(indexc)
          write(iulog,*)'ltype: ', ltype(clandunit(indexc))
          write(iulog,*)'ctype(indexc): ',ctype(indexc)
          write(iulog,*)'snl: ',snl(indexc)
          write(iulog,*)'h2osno: ',h2osno(indexc)
          write(iulog,*)'h2osno_old: ',h2osno_old(indexc)
          write(iulog,*)'snow_sources: ', snow_sources(indexc)
          write(iulog,*)'snow_sinks: ', snow_sinks(indexc)
          write(iulog,*)'qflx_prec_grnd: ',qflx_prec_grnd(indexc)*dtime
          write(iulog,*)'qflx_sub_snow: ',qflx_sub_snow(indexc)*dtime
          write(iulog,*)'qflx_evap_grnd: ',qflx_evap_grnd(indexc)*dtime
          write(iulog,*)'qflx_top_soil: ',qflx_top_soil(indexc)*dtime
          write(iulog,*)'qflx_dew_snow: ',qflx_dew_snow(indexc)*dtime
          write(iulog,*)'qflx_dew_grnd: ',qflx_dew_grnd(indexc)*dtime
          write(iulog,*)'qflx_snwcp_ice: ',qflx_snwcp_ice(indexc)*dtime
          write(iulog,*)'qflx_snwcp_liq: ',qflx_snwcp_liq(indexc)*dtime
          write(iulog,*)'qflx_sl_top_soil: ',qflx_sl_top_soil(indexc)*dtime
          if (create_glacier_mec_landunit) &
          write(iulog,*)'qflx_glcice_frz: ',qflx_glcice_frz(indexc)*dtime
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Energy balance checks

    do p = lbp, ubp
       if (pactive(p)) then
          l = plandunit(p)
          g = pgridcell(p)

          ! Solar radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (ltype(l) /= isturb) then
             errsol(p) = fsa(p) + fsr(p) &
                  - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2))
          else
             errsol(p) = spval
          end if
          
          ! Longwave radiation energy balance
          ! Do not do this check for an urban pft since it will not balance on a per-column
          ! level because of interactions between columns and since a separate check is done
          ! in the urban radiation module
          if (ltype(l) /= isturb) then
             errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(g)
          else
             errlon(p) = spval
          end if
          
          ! Surface energy balance
          ! Changed to using (eflx_lwrad_net) here instead of (forc_lwrad - eflx_lwrad_out) because
          ! there are longwave interactions between urban columns (and therefore pfts). 
          ! For surfaces other than urban, (eflx_lwrad_net) equals (forc_lwrad - eflx_lwrad_out),
          ! and a separate check is done above for these terms.
          
          if (ltype(l) /= isturb) then
             c=pcolumn(p)
             errseb(p) = sabv(p) + sabg_chk(p) + forc_lwrad(g) - eflx_lwrad_out(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
          else
             errseb(p) = sabv(p) + sabg(p) &
                         - eflx_lwrad_net(p) &
                         - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p) &
                         + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
          end if
          netrad(p) = fsa(p) - eflx_lwrad_net(p)
       end if
    end do

    ! Solar radiation energy balance check

    found = .false.
    do p = lbp, ubp
       if (pactive(p)) then
          if ( (errsol(p) /= spval) .and. (abs(errsol(p)) > .10_r8) ) then
             found = .true.
             indexp = p
             indexg = pgridcell(p)
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: solar radiation balance error', nstep, indexp, errsol(indexp)
       write(iulog,*)'fsa          = ',fsa(indexp)
       write(iulog,*)'fsr          = ',fsr(indexp)
       write(iulog,*)'forc_solad(1)= ',forc_solad(indexg,1)
       write(iulog,*)'forc_solad(2)= ',forc_solad(indexg,2)
       write(iulog,*)'forc_solai(1)= ',forc_solai(indexg,1)
       write(iulog,*)'forc_solai(2)= ',forc_solai(indexg,2)
       write(iulog,*)'forc_tot     = ',forc_solad(indexg,1)+forc_solad(indexg,2)&
                                  +forc_solai(indexg,1)+forc_solai(indexg,2)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Longwave radiation energy balance check

    found = .false.
    do p = lbp, ubp
       if (pactive(p)) then
          if ( (errlon(p) /= spval) .and. (abs(errlon(p)) > .10_r8) ) then
             found = .true.
             indexp = p
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: longwave enery balance error',nstep,indexp,errlon(indexp)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Surface energy balance check

    found = .false.
    do p = lbp, ubp
       if (pactive(p)) then
          if (abs(errseb(p)) > .10_r8 ) then
             found = .true.
             indexp = p
          end if
       end if
    end do
    if ( found  .and. (nstep > 2) ) then
       write(iulog,100)'BalanceCheck: surface flux energy balance error',nstep,indexp,errseb(indexp)
       write(iulog,*)' sabv           = ',sabv(indexp)
       c=pcolumn(indexp)
       write(iulog,*)' column      = ',c
       write(iulog,*)' sabg           = ',sabg(indexp), ((1._r8- frac_sno(c))*sabg_soil(indexp) + frac_sno(c)*sabg_snow(indexp)),sabg_chk(indexp)
       write(iulog,*)' eflx_lwrad_net = ',eflx_lwrad_net(indexp)
       write(iulog,*)' eflx_sh_tot    = ',eflx_sh_tot(indexp)
       write(iulog,*)' eflx_lh_tot    = ',eflx_lh_tot(indexp)
       write(iulog,*)' eflx_soil_grnd = ',eflx_soil_grnd(indexp)
       write(iulog,*)'clm model is stopping'
       call endrun()
    end if

    ! Soil energy balance check

    found = .false.
    do c = lbc, ubc
       if (abs(errsoi_col(c)) > 1.0e-7_r8 ) then
          found = .true.
          indexc = c
       end if
    end do
    if ( found ) then
       if (abs(errsoi_col(indexc)) > .10_r8 .and. (nstep > 2) ) then
          write(iulog,100)'BalanceCheck: soil balance error',nstep,indexc,errsoi_col(indexc)
          write(iulog,*)'nstep = ',nstep,' indexc= ',indexc,' errsoi_col= ',errsoi_col(indexc)
          write(iulog,*)'clm model is stopping'
          call endrun()
       end if
    end if

    ! Update SH and RUNOFF for dynamic land cover change energy and water fluxes
    call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                &
              qflx_runoff(lbc:ubc), qflx_runoffg(lbg:ubg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = lbg, ubg
       qflx_runoffg(g) = qflx_runoffg(g) - qflx_liq_dynbal(g)
    enddo

    call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                      &
              qflx_snwcp_ice(lbc:ubc), qflx_snwcp_iceg(lbg:ubg), &
              c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
    do g = lbg, ubg
       qflx_snwcp_iceg(g) = qflx_snwcp_iceg(g) - qflx_ice_dynbal(g)
    enddo

    call p2g( lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg,      &
              eflx_sh_tot(lbp:ubp), eflx_sh_totg(lbg:ubg), &
              p2c_scale_type='unity',c2l_scale_type='urbanf',l2g_scale_type='unity')
    do g = lbg, ubg
       eflx_sh_totg(g) =  eflx_sh_totg(g) - eflx_dynbal(g)
    enddo

! calculate total water storage for history files
! first set tws to gridcell total endwb
    call c2g( lbc, ubc, lbl, ubl, lbg, ubg,                &
         endwb(lbc:ubc), tws(lbg:ubg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity' )
! second add river storage as gridcell average depth
! 1.e-3 converts [m3/km2] to [mm]
    do g = lbg, ubg
       tws(g) = tws(g) + volr(g) / area(g) * 1.e-3_r8
    enddo
100 format (1x,a,' nstep =',i10,' point =',i6,' imbalance =',f12.6,' W/m2')
200 format (1x,a,' nstep =',i10,' point =',i6,' imbalance =',f12.6,' mm')

  end subroutine BalanceCheck

end module BalanceCheckMod
