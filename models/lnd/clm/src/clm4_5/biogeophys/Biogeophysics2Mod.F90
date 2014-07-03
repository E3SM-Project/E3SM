module Biogeophysics2Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Biogeophysics2Mod
!
! !DESCRIPTION:
! Performs the calculation of soil/snow and ground temperatures
! and updates surface fluxes based on the new ground temperature.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
  use clm_varctl,   only: iulog
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics2   ! Calculate soil/snow and ground temperatures
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
! !IROUTINE: Biogeophysics2
!
! !INTERFACE:
  subroutine Biogeophysics2 (lbl, ubl, lbc, ubc, lbp, ubp, &
             num_urbanl, filter_urbanl, num_nolakec, filter_nolakec, &
             num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of soil/snow and
! ground temperatures and update surface fluxes based on the new ground
! temperature
!
! Calling sequence is:
! Biogeophysics2:             surface biogeophysics driver
!    -> SoilTemperature:      soil/snow and ground temperatures
!          -> SoilTermProp    thermal conductivities and heat capacities
!          -> Tridiagonal     tridiagonal matrix solution
!          -> PhaseChange     phase change of liquid/ice contents
!
! (1) Snow and soil temperatures
!     o The volumetric heat capacity is calculated as a linear combination
!       in terms of the volumetric fraction of the constituent phases.
!     o The thermal conductivity of soil is computed from
!       the algorithm of Johansen (as reported by Farouki 1981), and the
!       conductivity of snow is from the formulation used in
!       SNTHERM (Jordan 1991).
!     o Boundary conditions:
!       F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
!     o Soil / snow temperature is predicted from heat conduction
!       in 10 soil layers and up to 5 snow layers.
!       The thermal conductivities at the interfaces between two
!       neighboring layers (j, j+1) are derived from an assumption that
!       the flux across the interface is equal to that from the node j
!       to the interface and the flux from the interface to the node j+1.
!       The equation is solved using the Crank-Nicholson method and
!       results in a tridiagonal system equation.
!
! (2) Phase change (see PhaseChange.F90)
!
! !USES:
    use clmtype
    use clm_atmlnd        , only : clm_a2l
    use clm_time_manager  , only : get_step_size
    use clm_varcon        , only : hvap, cpair, grav, vkc, tfrz, sb, icol_road_perv, &
                                   isturb, icol_roof, icol_sunwall, icol_shadewall, istsoil, &
                                   istice_mec
    use clm_varcon        , only : istcrop
    use clm_varpar        , only : nlevsno, nlevgrnd, nlevurb, max_pft_per_col
    use SoilTemperatureMod, only : SoilTemperature
    use subgridAveMod     , only : p2c
    use perf_mod          , only : t_startf, t_stopf
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                    ! pft bounds
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: lbl, ubl                    ! landunit bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in) :: num_urbanl                  ! number of urban landunits in clump
    integer, intent(in) :: filter_urbanl(ubl-lbl+1)    ! urban landunit filter
    integer, intent(in) :: num_nolakep                 ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)   ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! Migrated to clm2.0 by Keith Oleson and Mariana Vertenstein
! Migrated to clm2.1 new data structures by Peter Thornton and M. Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: frac_sno_eff(:)    !eff. fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno(:)        ! fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: h2osfc(:)          ! surface water (mm)
    real(r8), pointer :: t_h2osfc(:) 	    ! surface water temperature
    real(r8), pointer :: t_h2osfc_bef(:)    ! saved surface water temperature
    real(r8), pointer :: frac_h2osfc(:)     ! fraction of ground covered by surface water (0 to 1)
    real(r8), pointer :: qflx_ev_snow(:)    ! evaporation flux from snow (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_ev_soil(:)    ! evaporation flux from soil (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_ev_h2osfc(:)  ! evaporation flux from soil (W/m**2) [+ to atm]
    real(r8), pointer :: sabg_soil(:)       ! solar radiation absorbed by soil (W/m**2)
    real(r8), pointer :: sabg_snow(:)       ! solar radiation absorbed by snow (W/m**2)
    integer , pointer :: ctype(:)           ! column type
    integer , pointer :: ltype(:)           ! landunit type
    logical , pointer :: pactive(:)         ! true=>do computations on this pft (see reweightMod for details)
    integer , pointer :: pcolumn(:)         ! pft's column index
    integer , pointer :: plandunit(:)       ! pft's landunit index
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: npfts(:)           ! column's number of pfts 
    integer , pointer :: pfti(:)            ! column's beginning pft index 
    integer , pointer :: snl(:)             ! number of snow layers
    logical , pointer :: do_capsnow(:)      ! true => do snow capping
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: emg(:)             ! ground emissivity
    real(r8), pointer :: htvp(:)            ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: t_grnd(:)          ! ground temperature (Kelvin)
    integer , pointer :: frac_veg_nosno(:)  ! fraction of vegetation not covered by snow (0 OR 1 now) [-]
    real(r8), pointer :: cgrnds(:)          ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)          ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: dlrad(:)           ! downward longwave radiation below the canopy [W/m2]
    real(r8), pointer :: ulrad(:)           ! upward longwave radiation above the canopy [W/m2]
    real(r8), pointer :: eflx_sh_veg(:)     ! sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_veg(:)   ! vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:)   ! vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_evap_can(:)   ! evaporation from leaves and stems (mm H2O/s) (+ = to atm)
    real(r8), pointer :: wtcol(:)           ! pft weight relative to column
    real(r8), pointer :: tssbef(:,:)        ! soil/snow temperature before update
    real(r8), pointer :: t_soisno(:,:)      ! soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_ice(:,:)    ! ice lens (kg/m2) (new)
    real(r8), pointer :: h2osoi_liq(:,:)    ! liquid water (kg/m2) (new)
    real(r8), pointer :: eflx_building_heat(:)   ! heat flux from urban building interior to walls, roof
    real(r8), pointer :: eflx_traffic_pft(:)    ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat_pft(:)  ! sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac_pft(:) ! sensible heat flux put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: canyon_hwr(:)      ! ratio of building height to street width (-)
 
! local pointers to implicit inout arguments
!
    real(r8), pointer :: eflx_sh_grnd(:)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_soi(:)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_snwcp_liq(:)  ! excess rainfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_ice(:)  ! excess snowfall due to snow capping (mm H2O /s)
! 
! local pointers to implicit out arguments
! 
    real(r8), pointer :: dt_grnd(:)         ! change in t_grnd, last iteration (Kelvin)
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_soil_grnd_u(:)! urban soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_soil_grnd_r(:)! rural soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_u(:)   ! urban total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_r(:)   ! rural total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_tot_u(:)   ! urban total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_tot_r(:)   ! rural total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: qflx_evap_grnd(:)  ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), pointer :: qflx_sub_snow(:)   ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_snow(:)   ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), pointer :: qflx_dew_grnd(:)   ! ground surface dew formation (mm H2O /s) [+]
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_u(:)  ! urban net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lwrad_net_r(:)  ! rural net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: eflx_lh_vege(:)    ! veg evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_vegt(:)    ! veg transpiration heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_grnd(:)    ! ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: errsoi_pft(:)      ! pft-level soil/lake energy conservation error (W/m**2)
    real(r8), pointer :: errsoi_col(:)      ! column-level soil/lake energy conservation error (W/m**2)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: p,c,g,j,pi,l         ! indices
    integer  :: fc,fp                ! lake filtered column and pft indices
    real(r8) :: dtime                ! land model time step (sec)
    real(r8) :: egsmax(lbc:ubc)      ! max. evaporation which soil can provide at one time step
    real(r8) :: egirat(lbc:ubc)      ! ratio of topsoil_evap_tot : egsmax
    real(r8) :: tinc(lbc:ubc)        ! temperature difference of two time step
    real(r8) :: xmf(lbc:ubc)         ! total latent heat of phase change of ground water
    real(r8) :: sumwt(lbc:ubc)       ! temporary
    real(r8) :: evaprat(lbp:ubp)     ! ratio of qflx_evap_soi/topsoil_evap_tot
    real(r8) :: save_qflx_evap_soi   ! temporary storage for qflx_evap_soi
    real(r8) :: topsoil_evap_tot(lbc:ubc)          ! column-level total evaporation from top soil layer
    real(r8) :: fact(lbc:ubc, -nlevsno+1:nlevgrnd) ! used in computing tridiagonal matrix
    real(r8) :: eflx_lwrad_del(lbp:ubp)            ! update due to eflx_lwrad
    real(r8) :: t_grnd0(lbc:ubc)    !t_grnd of previous time step
    real(r8) :: c_h2osfc(lbc:ubc)   !heat capacity of surface water
    real(r8) :: xmf_h2osfc(lbc:ubc) !latent heat of phase change of surface water
    real(r8) :: eflx_temp(lbc:ubc)
    real(r8) :: xmf_temp(lbc:ubc)
    real(r8) :: dq_temp(lbc:ubc)
    real(r8) :: lw_grnd
    real(r8) :: fsno_eff
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (gridcell-level)

    forc_lwrad => clm_a2l%forc_lwrad

    ! Assign local pointers to derived subtypes components (landunit-level)

    ltype      => clm3%g%l%itype
    canyon_hwr     => clm3%g%l%canyon_hwr

    ! Assign local pointers to derived subtypes components (column-level)

    frac_sno_eff  => clm3%g%l%c%cps%frac_sno_eff 
    frac_sno      => clm3%g%l%c%cps%frac_sno 
    h2osfc        => clm3%g%l%c%cws%h2osfc
    frac_h2osfc   => clm3%g%l%c%cps%frac_h2osfc
    t_h2osfc      => clm3%g%l%c%ces%t_h2osfc
    t_h2osfc_bef  => clm3%g%l%c%ces%t_h2osfc_bef
    qflx_ev_snow  => clm3%g%l%c%p%pwf%qflx_ev_snow
    qflx_ev_soil  => clm3%g%l%c%p%pwf%qflx_ev_soil
    qflx_ev_h2osfc=> clm3%g%l%c%p%pwf%qflx_ev_h2osfc
    sabg_soil     => clm3%g%l%c%p%pef%sabg_soil
    sabg_snow     => clm3%g%l%c%p%pef%sabg_snow
    ctype      => clm3%g%l%c%itype
    npfts      => clm3%g%l%c%npfts
    pfti       => clm3%g%l%c%pfti
    snl        => clm3%g%l%c%cps%snl
    do_capsnow => clm3%g%l%c%cps%do_capsnow
    htvp       => clm3%g%l%c%cps%htvp
    emg        => clm3%g%l%c%cps%emg
    t_grnd     => clm3%g%l%c%ces%t_grnd
    dt_grnd    => clm3%g%l%c%ces%dt_grnd
    t_soisno   => clm3%g%l%c%ces%t_soisno
    tssbef     => clm3%g%l%c%ces%tssbef
    h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
    errsoi_col => clm3%g%l%c%cebal%errsoi
    eflx_building_heat => clm3%g%l%c%cef%eflx_building_heat

    ! Assign local pointers to derived subtypes components (pft-level)

    pactive        => clm3%g%l%c%p%active
    pcolumn        => clm3%g%l%c%p%column
    plandunit      => clm3%g%l%c%p%landunit
    pgridcell      => clm3%g%l%c%p%gridcell
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    sabg           => clm3%g%l%c%p%pef%sabg
    dlrad          => clm3%g%l%c%p%pef%dlrad
    ulrad          => clm3%g%l%c%p%pef%ulrad
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_veg    => clm3%g%l%c%p%pef%eflx_sh_veg
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_veg  => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
    qflx_evap_can  => clm3%g%l%c%p%pwf%qflx_evap_can
    qflx_snwcp_liq => clm3%g%l%c%p%pwf%qflx_snwcp_liq
    qflx_snwcp_ice => clm3%g%l%c%p%pwf%qflx_snwcp_ice
    qflx_evap_tot  => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_grnd => clm3%g%l%c%p%pwf%qflx_evap_grnd
    qflx_sub_snow  => clm3%g%l%c%p%pwf%qflx_sub_snow
    qflx_dew_snow  => clm3%g%l%c%p%pwf%qflx_dew_snow
    qflx_dew_grnd  => clm3%g%l%c%p%pwf%qflx_dew_grnd
    eflx_soil_grnd => clm3%g%l%c%p%pef%eflx_soil_grnd
    eflx_soil_grnd_u => clm3%g%l%c%p%pef%eflx_soil_grnd_u
    eflx_soil_grnd_r => clm3%g%l%c%p%pef%eflx_soil_grnd_r
    eflx_sh_tot    => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_sh_tot_u  => clm3%g%l%c%p%pef%eflx_sh_tot_u
    eflx_sh_tot_r  => clm3%g%l%c%p%pef%eflx_sh_tot_r
    eflx_lh_tot    => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_lh_tot_u    => clm3%g%l%c%p%pef%eflx_lh_tot_u
    eflx_lh_tot_r    => clm3%g%l%c%p%pef%eflx_lh_tot_r
    eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net => clm3%g%l%c%p%pef%eflx_lwrad_net
    eflx_lwrad_net_u => clm3%g%l%c%p%pef%eflx_lwrad_net_u
    eflx_lwrad_net_r => clm3%g%l%c%p%pef%eflx_lwrad_net_r
    eflx_lh_vege   => clm3%g%l%c%p%pef%eflx_lh_vege
    eflx_lh_vegt   => clm3%g%l%c%p%pef%eflx_lh_vegt
    eflx_lh_grnd   => clm3%g%l%c%p%pef%eflx_lh_grnd
    cgrnds         => clm3%g%l%c%p%pef%cgrnds
    cgrndl         => clm3%g%l%c%p%pef%cgrndl
    eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
    qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
    errsoi_pft     => clm3%g%l%c%p%pebal%errsoi
    wtcol          => clm3%g%l%c%p%wtcol
    eflx_wasteheat_pft => clm3%g%l%c%p%pef%eflx_wasteheat_pft
    eflx_heat_from_ac_pft => clm3%g%l%c%p%pef%eflx_heat_from_ac_pft
    eflx_traffic_pft => clm3%g%l%c%p%pef%eflx_traffic_pft

    ! Get step size

    dtime = get_step_size()

    ! Determine soil temperatures including surface soil temperature

    call t_startf('soiltemperature')
    call SoilTemperature(lbl, ubl, lbc, ubc, num_urbanl, filter_urbanl, &
                         num_nolakec, filter_nolakec, xmf , fact, c_h2osfc, xmf_h2osfc)
    call t_stopf('soiltemperature')

    call t_startf('bgp2_loop_1')
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       j = snl(c)+1

       ! Calculate difference in soil temperature from last time step, for
       ! flux corrections

       if (snl(c) < 0) then
          t_grnd0(c) = frac_sno_eff(c) * tssbef(c,snl(c)+1) &
               + (1 - frac_sno_eff(c) - frac_h2osfc(c)) * tssbef(c,1) &
               + frac_h2osfc(c) * t_h2osfc_bef(c)
       else
          t_grnd0(c) = (1 - frac_h2osfc(c)) * tssbef(c,1) + frac_h2osfc(c) * t_h2osfc_bef(c)
       endif

       tinc(c) = t_grnd(c) - t_grnd0(c)

       ! Determine ratio of topsoil_evap_tot

       egsmax(c) = (h2osoi_ice(c,j)+h2osoi_liq(c,j)) / dtime

       ! added to trap very small negative soil water,ice

       if (egsmax(c) < 0._r8) then
          egsmax(c) = 0._r8
       end if
    end do

    ! A preliminary pft loop to determine if corrections are required for
    ! excess evaporation from the top soil layer... Includes new logic
    ! to distribute the corrections between pfts on the basis of their
    ! evaporative demands.
    ! egirat holds the ratio of demand to availability if demand is
    ! greater than availability, or 1.0 otherwise.
    ! Correct fluxes to present soil temperature

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       eflx_sh_grnd(p) = eflx_sh_grnd(p) + tinc(c)*cgrnds(p)
       qflx_evap_soi(p) = qflx_evap_soi(p) + tinc(c)*cgrndl(p)

       ! set ev_snow, ev_soil for urban landunits here
       l = plandunit(p)
       if (ltype(l) == isturb) then
          qflx_ev_snow(p) = qflx_evap_soi(p)
          qflx_ev_soil(p) = 0._r8
          qflx_ev_h2osfc(p) = 0._r8
       else
          qflx_ev_snow(p) = qflx_ev_snow(p) + tinc(c)*cgrndl(p)
          qflx_ev_soil(p) = qflx_ev_soil(p) + tinc(c)*cgrndl(p)
          qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) + tinc(c)*cgrndl(p)
       endif
    end do

    ! Set the column-average qflx_evap_soi as the weighted average over all pfts
    ! but only count the pfts that are evaporating

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       topsoil_evap_tot(c) = 0._r8
       sumwt(c) = 0._r8
    end do

    do pi = 1,max_pft_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                topsoil_evap_tot(c) = topsoil_evap_tot(c) + qflx_evap_soi(p) * wtcol(p)
             end if
          end if
       end do
    end do
    call t_stopf('bgp2_loop_1')
    call t_startf('bgp2_loop_2')

    ! Calculate ratio for rescaling pft-level fluxes to meet availability

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       if (topsoil_evap_tot(c) > egsmax(c)) then
          egirat(c) = (egsmax(c)/topsoil_evap_tot(c))
       else
          egirat(c) = 1.0_r8
       end if
    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       ! Correct soil fluxes for possible evaporation in excess of top layer water
       ! excess energy is added to the sensible heat flux from soil

       if (egirat(c) < 1.0_r8) then
          save_qflx_evap_soi = qflx_evap_soi(p)
          qflx_evap_soi(p) = qflx_evap_soi(p) * egirat(c)
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + (save_qflx_evap_soi - qflx_evap_soi(p))*htvp(c)
          qflx_ev_snow(p) = qflx_ev_snow(p) * egirat(c)
          qflx_ev_soil(p) = qflx_ev_soil(p) * egirat(c)
          qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) * egirat(c)
       end if

       ! Ground heat flux

       if (ltype(l) /= isturb) then
          lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
               +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
               +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

          eflx_soil_grnd(p) = ((1._r8- frac_sno_eff(c))*sabg_soil(p) + frac_sno_eff(c)*sabg_snow(p)) + dlrad(p) &
               + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(g) &
            - emg(c)*sb*lw_grnd - emg(c)*sb*t_grnd0(c)**3*(4._r8*tinc(c)) &
               - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
          
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
            eflx_soil_grnd_r(p) = eflx_soil_grnd(p)
          end if
       else
          ! For all urban columns we use the net longwave radiation (eflx_lwrad_net) since
          ! the term (emg*sb*tssbef(snl+1)**4) is not the upward longwave flux because of 
          ! interactions between urban columns.

          eflx_lwrad_del(p) = 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)

          ! Include transpiration term because needed for pervious road
          ! and wasteheat and traffic flux
          eflx_soil_grnd(p) = sabg(p) + dlrad(p) &
                              - eflx_lwrad_net(p) - eflx_lwrad_del(p) &
                              - (eflx_sh_grnd(p) + qflx_evap_soi(p)*htvp(c) + qflx_tran_veg(p)*hvap) &
                              + eflx_wasteheat_pft(p) + eflx_heat_from_ac_pft(p) + eflx_traffic_pft(p)
          eflx_soil_grnd_u(p) = eflx_soil_grnd(p)
       end if

       ! Total fluxes (vegetation + ground)

       eflx_sh_tot(p) = eflx_sh_veg(p) + eflx_sh_grnd(p)
       qflx_evap_tot(p) = qflx_evap_veg(p) + qflx_evap_soi(p)
       eflx_lh_tot(p)= hvap*qflx_evap_veg(p) + htvp(c)*qflx_evap_soi(p)
       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         eflx_lh_tot_r(p)= eflx_lh_tot(p)
         eflx_sh_tot_r(p)= eflx_sh_tot(p)
       else if (ltype(l) == isturb) then
         eflx_lh_tot_u(p)= eflx_lh_tot(p)
         eflx_sh_tot_u(p)= eflx_sh_tot(p)
       end if

       ! Assign ground evaporation to sublimation from soil ice or to dew
       ! on snow or ground

       qflx_evap_grnd(p) = 0._r8
       qflx_sub_snow(p) = 0._r8
       qflx_dew_snow(p) = 0._r8
       qflx_dew_grnd(p) = 0._r8

       if (qflx_ev_snow(p) >= 0._r8) then
          ! for evaporation partitioning between liquid evap and ice sublimation, 
	  ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
	  if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0.) then
             qflx_evap_grnd(p) = max(qflx_ev_snow(p)*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
	  else
	     qflx_evap_grnd(p) = 0.
	  end if
          qflx_sub_snow(p) = qflx_ev_snow(p) - qflx_evap_grnd(p)
       else
          if (t_grnd(c) < tfrz) then
             qflx_dew_snow(p) = abs(qflx_ev_snow(p))
          else
             qflx_dew_grnd(p) = abs(qflx_ev_snow(p))
          end if
       end if

       ! Update the pft-level qflx_snwcp
       ! This was moved in from Hydrology2 to keep all pft-level
       ! calculations out of Hydrology2

       if (snl(c) < 0 .and. do_capsnow(c)) then
          qflx_snwcp_liq(p) = qflx_snwcp_liq(p)+frac_sno_eff(c)*qflx_dew_grnd(p)
          qflx_snwcp_ice(p) = qflx_snwcp_ice(p)+frac_sno_eff(c)*qflx_dew_snow(p)
       end if

       ! Variables needed by history tape

       qflx_evap_can(p)  = qflx_evap_veg(p) - qflx_tran_veg(p)
       eflx_lh_vege(p)   = (qflx_evap_veg(p) - qflx_tran_veg(p)) * hvap
       eflx_lh_vegt(p)   = qflx_tran_veg(p) * hvap
       eflx_lh_grnd(p)   = qflx_evap_soi(p) * htvp(c)

    end do
    call t_stopf('bgp2_loop_2')
    call t_startf('bgp2_loop_3')

    ! Soil Energy balance check

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       errsoi_pft(p) = eflx_soil_grnd(p) - xmf(c) - xmf_h2osfc(c) &
            - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c)) &
            *(c_h2osfc(c)/dtime)

       ! For urban sunwall, shadewall, and roof columns, the "soil" energy balance check
       ! must include the heat flux from the interior of the building.
       if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall .or. ctype(c)==icol_roof) then
          errsoi_pft(p) = errsoi_pft(p) + eflx_building_heat(c) 
       end if
    end do
    do j = -nlevsno+1,nlevgrnd
       do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          if ((ctype(c) /= icol_sunwall .and. ctype(c) /= icol_shadewall &
               .and. ctype(c) /= icol_roof) .or. ( j <= nlevurb)) then
            ! area weight heat absorbed by snow layers
             if (j >= snl(c)+1 .and. j < 1) errsoi_pft(p) = errsoi_pft(p) &
                  - frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
             if (j >= 1) errsoi_pft(p) = errsoi_pft(p) &
                  - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
          end if
       end do
    end do
    call t_stopf('bgp2_loop_3')
    call t_startf('bgp2_loop_4')

    ! Outgoing long-wave radiation from vegetation + ground
    ! For conservation we put the increase of ground longwave to outgoing
    ! For urban pfts, ulrad=0 and (1-fracveg_nosno)=1, and eflx_lwrad_out and eflx_lwrad_net 
    ! are calculated in UrbanRadiation. The increase of ground longwave is added directly 
    ! to the outgoing longwave and the net longwave.

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       if (ltype(l) /= isturb) then
          lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
               +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
               +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

          eflx_lwrad_out(p) = ulrad(p) &
               + (1-frac_veg_nosno(p))*(1.-emg(c))*forc_lwrad(g) &
               + (1-frac_veg_nosno(p))*emg(c)*sb*lw_grnd &
               + 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)

          eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(g)
          if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
            eflx_lwrad_net_r(p) = eflx_lwrad_out(p) - forc_lwrad(g)
          end if
       else
          eflx_lwrad_out(p) = eflx_lwrad_out(p) + eflx_lwrad_del(p)
          eflx_lwrad_net(p) = eflx_lwrad_net(p) + eflx_lwrad_del(p)
          eflx_lwrad_net_u(p) = eflx_lwrad_net_u(p) + eflx_lwrad_del(p)
       end if
    end do

    ! lake balance for errsoi is not over pft
    ! therefore obtain column-level radiative temperature

    call p2c(num_nolakec, filter_nolakec, errsoi_pft, errsoi_col)
    call t_stopf('bgp2_loop_4')

  end subroutine Biogeophysics2

end module Biogeophysics2Mod
