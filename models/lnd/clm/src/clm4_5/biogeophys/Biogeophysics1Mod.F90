module Biogeophysics1Mod

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Biogeophysics1Mod
!
! !DESCRIPTION:
! Performs calculation of leaf temperature and surface fluxes.
! Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varctl  , only : iulog
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: Biogeophysics1   ! Calculate leaf temperature and surface fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Biogeophysics1
!
! !INTERFACE:
  subroutine Biogeophysics1(lbg, ubg, lbc, ubc, lbp, ubp, &
       num_nolakec, filter_nolakec, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! This is the main subroutine to execute the calculation of leaf temperature
! and surface fluxes. Biogeophysics2.F90 then determines soil/snow and ground
! temperatures and updates the surface fluxes for the new ground
! temperature.
!
! Calling sequence is:
! Biogeophysics1:           surface biogeophysics driver
!  -> QSat:                 saturated vapor pressure, specific humidity, and
!                           derivatives at ground surface and derivatives at
!                           leaf surface using updated leaf temperature
! Leaf temperature
! Foliage energy conservation is given by the foliage energy budget
! equation:
!                Rnet - Hf - LEf = 0
! The equation is solved by Newton-Raphson iteration, in which this
! iteration includes the calculation of the photosynthesis and
! stomatal resistance, and the integration of turbulent flux profiles.
! The sensible and latent heat transfer between foliage and atmosphere
! and ground is linked by the equations:
!                Ha = Hf + Hg and Ea = Ef + Eg
!
! !USES:
    use clmtype
    use clm_atmlnd         , only : clm_a2l
    use clm_varcon         , only : denh2o, denice, roverg, hvap, hsub, &
                                    istice, istice_mec, istwet, istsoil, isturb, istdlak, &
                                    zlnd, zsno, tfrz, &
                                    icol_roof, icol_sunwall, icol_shadewall,     &
                                    icol_road_imperv, icol_road_perv, tfrz, spval, istdlak
    use clm_varcon         , only : istcrop
    use clm_varpar         , only : nlevgrnd, nlevurb, nlevsno, max_pft_per_gcell, nlevsoi
    use QSatMod            , only : QSat
    use shr_const_mod      , only : SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg                    ! gridcell-index bounds
    integer, intent(in) :: lbc, ubc                    ! column-index bounds
    integer, intent(in) :: lbp, ubp                    ! pft-index bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
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
! 27 February 2008: Keith Oleson; weighted soil/snow emissivity
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(r8), pointer :: frac_sno_eff(:)  ! eff. fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_h2osfc(:)   ! fraction of ground covered by surface water (0 to 1)
    real(r8), pointer :: h2osfc(:)        ! surface water (mm)
    real(r8), pointer :: t_h2osfc(:) 	  ! surface water temperature
    real(r8), pointer :: t_h2osfc_bef(:)  ! saved surface water temperature
    real(r8), pointer :: qg_snow(:)       ! specific humidity at snow surface [kg/kg]
    real(r8), pointer :: qg_soil(:)       ! specific humidity at soil surface [kg/kg]
    real(r8), pointer :: qg_h2osfc(:)     ! specific humidity at h2osfc surface [kg/kg]
    logical , pointer :: pactive(:)       !true=>do computations on this pft (see reweightMod for details)
    integer , pointer :: ivt(:)           !pft vegetation type
    integer , pointer :: ityplun(:)       !landunit type
    integer , pointer :: clandunit(:)     !column's landunit index
    integer , pointer :: cgridcell(:)     !column's gridcell index
    integer , pointer :: ctype(:)         !column type
    real(r8), pointer :: forc_pbot(:)     !atmospheric pressure (Pa)
    real(r8), pointer :: forc_q(:)        !atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_t(:)        !atmospheric temperature (Kelvin)
    real(r8), pointer :: forc_hgt_t(:)    !observational height of temperature [m]
    real(r8), pointer :: forc_hgt_u(:)    !observational height of wind [m]
    real(r8), pointer :: forc_hgt_q(:)    !observational height of specific humidity [m]
    integer , pointer :: npfts(:)         !number of pfts on gridcell
    integer , pointer :: pfti(:)          !initial pft on gridcell
    integer , pointer :: plandunit(:)     !pft's landunit index
    real(r8), pointer :: forc_hgt_u_pft(:) !observational height of wind at pft level [m]
    real(r8), pointer :: forc_hgt_t_pft(:) !observational height of temperature at pft level [m]
    real(r8), pointer :: forc_hgt_q_pft(:) !observational height of specific humidity at pft level [m]
    integer , pointer :: frac_veg_nosno(:) !fraction of vegetation not covered by snow (0 OR 1) [-]
    integer , pointer :: pgridcell(:)      !pft's gridcell index
    integer , pointer :: pcolumn(:)        !pft's column index
    real(r8), pointer :: z_0_town(:)      !momentum roughness length of urban landunit (m)
    real(r8), pointer :: z_d_town(:)      !displacement height of urban landunit (m)
    real(r8), pointer :: forc_th(:)       !atmospheric potential temperature (Kelvin)
    real(r8), pointer :: forc_u(:)        !atmospheric wind speed in east direction (m/s)
    real(r8), pointer :: forc_v(:)        !atmospheric wind speed in north direction (m/s)
    real(r8), pointer :: smpmin(:)        !restriction for min of soil potential (mm)
    integer , pointer :: snl(:)           !number of snow layers
    real(r8), pointer :: frac_sno(:)      !fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: h2osno(:)        !snow water (mm H2O)
    real(r8), pointer :: elai(:)          !one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)          !one-sided stem area index with burying by snow
    real(r8), pointer :: z0mr(:)          !ratio of momentum roughness length to canopy top height (-)
    real(r8), pointer :: displar(:)       !ratio of displacement height to canopy top height (-)
    real(r8), pointer :: htop(:)          !canopy top (m)
    real(r8), pointer :: dz(:,:)          !layer depth (m)
    real(r8), pointer :: t_soisno(:,:)    !soil temperature (Kelvin)
    real(r8), pointer :: h2osoi_liq(:,:)  !liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  !ice lens (kg/m2)
    real(r8), pointer :: watsat(:,:)      !volumetric soil water at saturation (porosity)
    real(r8), pointer :: sucsat(:,:)      !minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)         !Clapp and Hornberger "b"
    real(r8), pointer :: watfc(:,:)       !volumetric soil water at field capacity
    real(r8), pointer :: watopt(:,:)      !volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
    real(r8), pointer :: watdry(:,:)      !volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
    real(r8), pointer :: rootfr_road_perv(:,:) !fraction of roots in each soil layer for urban pervious road
    real(r8), pointer :: rootr_road_perv(:,:) !effective fraction of roots in each soil layer for urban pervious road
!
! local pointers to implicit out arguments
!
    real(r8), pointer :: t_grnd(:)        !ground temperature (Kelvin)
    real(r8), pointer :: qg(:)            !ground specific humidity [kg/kg]
    real(r8), pointer :: dqgdT(:)         !d(qg)/dT
    real(r8), pointer :: emg(:)           !ground emissivity
    real(r8), pointer :: htvp(:)          !latent heat of vapor of water (or sublimation) [j/kg]
    real(r8), pointer :: beta(:)          !coefficient of convective velocity [-]
    real(r8), pointer :: zii(:)           !convective boundary height [m]
    real(r8), pointer :: thm(:)           !intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(r8), pointer :: thv(:)           !virtual potential temperature (kelvin)
    real(r8), pointer :: z0mg(:)          !roughness length over ground, momentum [m]
    real(r8), pointer :: z0hg(:)          !roughness length over ground, sensible heat [m]
    real(r8), pointer :: z0qg(:)          !roughness length over ground, latent heat [m]
    real(r8), pointer :: emv(:)           !vegetation emissivity
    real(r8), pointer :: z0m(:)           !momentum roughness length (m)
    real(r8), pointer :: displa(:)        !displacement height (m)
    real(r8), pointer :: z0mv(:)          !roughness length over vegetation, momentum [m]
    real(r8), pointer :: z0hv(:)          !roughness length over vegetation, sensible heat [m]
    real(r8), pointer :: z0qv(:)          !roughness length over vegetation, latent heat [m]
    real(r8), pointer :: eflx_sh_tot(:)   !total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_u(:) !urban total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_sh_tot_r(:) !rural total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)   !total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_tot_u(:) !urban total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_lh_tot_r(:) !rural total latent heat flux (W/m**2)  [+ to atm]
    real(r8), pointer :: eflx_sh_veg(:)   !sensible heat flux from leaves (W/m**2) [+ to atm]
    real(r8), pointer :: qflx_evap_tot(:) !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(r8), pointer :: qflx_evap_veg(:) !vegetation evaporation (mm H2O/s) (+ = to atm)
    real(r8), pointer :: qflx_tran_veg(:) !vegetation transpiration (mm H2O/s) (+ = to atm)
    real(r8), pointer :: cgrnd(:)         !deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(r8), pointer :: cgrnds(:)        !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
    real(r8), pointer :: cgrndl(:)        !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
    real(r8) ,pointer :: tssbef(:,:)      !soil/snow temperature before update
    real(r8) ,pointer :: soilalpha(:)     !factor that reduces ground saturated specific humidity (-)
    real(r8) ,pointer :: soilbeta(:)      !factor that reduces ground evaporation
    real(r8) ,pointer :: soilalpha_u(:)   !Urban factor that reduces ground saturated specific humidity (-)

!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: g,l,c,p !indices
    integer  :: j       !soil/snow level index
    integer  :: fp      !lake filter pft index
    integer  :: fc      !lake filter column index
    real(r8) :: qred    !soil surface relative humidity
    real(r8) :: avmuir  !ir inverse optical depth per unit leaf area
    real(r8) :: eg      !water vapor pressure at temperature T [pa]
    real(r8) :: qsatg   !saturated humidity [kg/kg]
    real(r8) :: degdT   !d(eg)/dT
    real(r8) :: qsatgdT !d(qsatg)/dT
    real(r8) :: fac     !soil wetness of surface layer
    real(r8) :: psit    !negative potential of soil
    real(r8) :: hr      !relative humidity
    real(r8) :: hr_road_perv  !relative humidity for urban pervious road
    real(r8) :: wx      !partial volume of ice and water of surface layer
    real(r8) :: fac_fc        !soil wetness of surface layer relative to field capacity
    real(r8) :: eff_porosity  ! effective porosity in layer
    real(r8) :: vol_ice       ! partial volume of ice lens in layer
    real(r8) :: vol_liq       ! partial volume of liquid water in layer
    integer  :: pi            !index
    real(r8) :: fh2o_eff(lbc:ubc) ! effective surface water fraction (i.e. seen by atm)
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

    frac_sno_eff   => clm3%g%l%c%cps%frac_sno_eff
    frac_h2osfc   => clm3%g%l%c%cps%frac_h2osfc
    h2osfc        => clm3%g%l%c%cws%h2osfc
    t_h2osfc      => clm3%g%l%c%ces%t_h2osfc
    t_h2osfc_bef  => clm3%g%l%c%ces%t_h2osfc_bef
    qg_snow       => clm3%g%l%c%cws%qg_snow 
    qg_soil       => clm3%g%l%c%cws%qg_soil
    qg_h2osfc     => clm3%g%l%c%cws%qg_h2osfc
    forc_hgt_t    => clm_a2l%forc_hgt_t
    forc_u        => clm_a2l%forc_u
    forc_v        => clm_a2l%forc_v
    forc_hgt_u    => clm_a2l%forc_hgt_u
    forc_hgt_q    => clm_a2l%forc_hgt_q
    npfts         => clm3%g%npfts
    pfti          => clm3%g%pfti

    ! Assign local pointers to derived type members (landunit-level)

    ityplun       => clm3%g%l%itype
    z_0_town      => clm3%g%l%z_0_town
    z_d_town      => clm3%g%l%z_d_town

    ! Assign local pointers to derived type members (column-level)

    forc_pbot     => clm3%g%l%c%cps%forc_pbot
    forc_q        => clm3%g%l%c%cws%forc_q
    forc_t        => clm3%g%l%c%ces%forc_t
    forc_th       => clm3%g%l%c%ces%forc_th

    cgridcell     => clm3%g%l%c%gridcell
    clandunit     => clm3%g%l%c%landunit
    ctype         => clm3%g%l%c%itype
    beta          => clm3%g%l%c%cps%beta
    dqgdT         => clm3%g%l%c%cws%dqgdT
    emg           => clm3%g%l%c%cps%emg
    frac_sno      => clm3%g%l%c%cps%frac_sno
    h2osno        => clm3%g%l%c%cws%h2osno
    htvp          => clm3%g%l%c%cps%htvp
    qg            => clm3%g%l%c%cws%qg
    smpmin        => clm3%g%l%c%cps%smpmin
    snl           => clm3%g%l%c%cps%snl
    t_grnd        => clm3%g%l%c%ces%t_grnd
    thv           => clm3%g%l%c%ces%thv
    z0hg          => clm3%g%l%c%cps%z0hg
    z0mg          => clm3%g%l%c%cps%z0mg
    z0qg          => clm3%g%l%c%cps%z0qg
    zii           => clm3%g%l%c%cps%zii
    bsw           => clm3%g%l%c%cps%bsw
    dz            => clm3%g%l%c%cps%dz
    h2osoi_ice    => clm3%g%l%c%cws%h2osoi_ice
    h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
    soilalpha     => clm3%g%l%c%cws%soilalpha
    soilbeta      => clm3%g%l%c%cws%soilbeta
    soilalpha_u   => clm3%g%l%c%cws%soilalpha_u
    sucsat        => clm3%g%l%c%cps%sucsat
    t_soisno      => clm3%g%l%c%ces%t_soisno
    tssbef        => clm3%g%l%c%ces%tssbef
    watsat        => clm3%g%l%c%cps%watsat
    watfc         => clm3%g%l%c%cps%watfc
    watdry        => clm3%g%l%c%cps%watdry
    watopt        => clm3%g%l%c%cps%watopt
    rootfr_road_perv => clm3%g%l%c%cps%rootfr_road_perv
    rootr_road_perv  => clm3%g%l%c%cps%rootr_road_perv

    ! Assign local pointers to derived type members (pft-level)

    pactive       => clm3%g%l%c%p%active
    ivt           => clm3%g%l%c%p%itype
    elai          => clm3%g%l%c%p%pps%elai
    esai          => clm3%g%l%c%p%pps%esai
    htop          => clm3%g%l%c%p%pps%htop
    emv           => clm3%g%l%c%p%pps%emv
    z0m           => clm3%g%l%c%p%pps%z0m
    displa        => clm3%g%l%c%p%pps%displa
    z0mv          => clm3%g%l%c%p%pps%z0mv
    z0hv          => clm3%g%l%c%p%pps%z0hv
    z0qv          => clm3%g%l%c%p%pps%z0qv
    eflx_sh_tot   => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_sh_tot_u => clm3%g%l%c%p%pef%eflx_sh_tot_u
    eflx_sh_tot_r => clm3%g%l%c%p%pef%eflx_sh_tot_r
    eflx_lh_tot   => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_lh_tot_u => clm3%g%l%c%p%pef%eflx_lh_tot_u
    eflx_lh_tot_r => clm3%g%l%c%p%pef%eflx_lh_tot_r
    eflx_sh_veg   => clm3%g%l%c%p%pef%eflx_sh_veg
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    qflx_evap_veg => clm3%g%l%c%p%pwf%qflx_evap_veg
    qflx_tran_veg => clm3%g%l%c%p%pwf%qflx_tran_veg
    cgrnd         => clm3%g%l%c%p%pef%cgrnd
    cgrnds        => clm3%g%l%c%p%pef%cgrnds
    cgrndl        => clm3%g%l%c%p%pef%cgrndl
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
    forc_hgt_t_pft => clm3%g%l%c%p%pps%forc_hgt_t_pft
    forc_hgt_q_pft => clm3%g%l%c%p%pps%forc_hgt_q_pft
    plandunit      => clm3%g%l%c%p%landunit
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    thm            => clm3%g%l%c%p%pes%thm
    pgridcell      => clm3%g%l%c%p%gridcell
    pcolumn        => clm3%g%l%c%p%column

    ! Assign local pointers to derived type members (ecophysiological)

    z0mr          => pftcon%z0mr
    displar       => pftcon%displar

    do j = -nlevsno+1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ((ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall &
               .or. ctype(c) == icol_roof) .and. j > nlevurb) then
            tssbef(c,j) = spval 
          else
            tssbef(c,j) = t_soisno(c,j)
          end if
          ! record t_h2osfc prior to updating
          t_h2osfc_bef(c) = t_h2osfc(c)   
       end do
    end do

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)

       if (ctype(c) == icol_road_perv) then
          hr_road_perv = 0._r8
       end if

       ! begin calculations that relate only to the column level
       ! Ground and soil temperatures from previous time step

       ! ground temperature is weighted average of exposed soil, snow, and h2osfc
       if (snl(c) < 0) then
          t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
               + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
               + frac_h2osfc(c) * t_h2osfc(c)
       else
          t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
       endif
       ! Saturated vapor pressure, specific humidity and their derivatives
       ! at ground surface

       qred = 1._r8
       if (ityplun(l)/=istwet .AND. ityplun(l)/=istice  &
                              .AND. ityplun(l)/=istice_mec) then
          if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
             wx   = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)
             fac  = min(1._r8, wx/watsat(c,1))
             fac  = max( fac, 0.01_r8 )
             psit = -sucsat(c,1) * fac ** (-bsw(c,1))
             psit = max(smpmin(c), psit)
             ! modify qred to account for h2osfc
             hr   = exp(psit/roverg/t_soisno(c,1))
             qred = (1._r8 - frac_sno(c) - frac_h2osfc(c))*hr &
                  + frac_sno(c) + frac_h2osfc(c)

             !! Lee and Pielke 1992 beta, added by K.Sakaguchi
             if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
                fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
                fac_fc  = max( fac_fc, 0.01_r8 )
                ! modify soil beta by snow cover. soilbeta for snow surface is one
                soilbeta(c) = (1._r8-frac_sno(c)-frac_h2osfc(c)) &
                     *0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                              + frac_sno(c)+ frac_h2osfc(c)
             else   !when water content of ths top layer is more than that at F.C.
                soilbeta(c) = 1._r8
             end if

             soilalpha(c) = qred
          ! Pervious road depends on water in total soil column
          else if (ctype(c) == icol_road_perv) then
             do j = 1, nlevsoi
                if (t_soisno(c,j) >= tfrz) then
                   vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
                   eff_porosity = watsat(c,j)-vol_ice
                   vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
                   fac = min( max(vol_liq-watdry(c,j),0._r8) / (watopt(c,j)-watdry(c,j)), 1._r8 )
                else
                   fac = 0._r8
                end if
                rootr_road_perv(c,j) = rootfr_road_perv(c,j)*fac
                hr_road_perv = hr_road_perv + rootr_road_perv(c,j)
             end do
             ! Allows for sublimation of snow or dew on snow
             qred = (1.-frac_sno(c))*hr_road_perv + frac_sno(c)

             ! Normalize root resistances to get layer contribution to total ET
             if (hr_road_perv .gt. 0._r8) then
                do j = 1, nlevsoi
                   rootr_road_perv(c,j) = rootr_road_perv(c,j)/hr_road_perv
                end do
             end if
             soilalpha_u(c) = qred
             soilbeta(c) = 0._r8
          else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
             qred = 0._r8
             soilbeta(c) = 0._r8
             soilalpha_u(c) = spval
          else if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
             qred = 1._r8
             soilbeta(c) = 0._r8
             soilalpha_u(c) = spval
          end if
       else
          soilalpha(c) = spval
          soilbeta(c) =   1._r8
       end if

       ! compute humidities individually for snow, soil, h2osfc for vegetated landunits
       if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then

          call QSat(t_soisno(c,snl(c)+1), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
          if (qsatg > forc_q(c) .and. forc_q(c) > qsatg) then
             qsatg = forc_q(c)
             qsatgdT = 0._r8
          end if

          qg_snow(c) = qsatg
          dqgdT(c) = frac_sno(c)*qsatgdT

          call QSat(t_soisno(c,1) , forc_pbot(c), eg, degdT, qsatg, qsatgdT)
          if (qsatg > forc_q(c) .and. forc_q(c) > hr*qsatg) then
             qsatg = forc_q(c)
             qsatgdT = 0._r8
          end if
          qg_soil(c) = hr*qsatg
          dqgdT(c) = dqgdT(c) + (1._r8 - frac_sno(c) - frac_h2osfc(c))*hr*qsatgdT

          ! to be consistent with hs_top values in SoilTemp, set qg_snow to qg_soil for snl = 0 case
          ! this ensures hs_top_snow will equal hs_top_soil
          if (snl(c) >= 0) then
             qg_snow(c) = qg_soil(c)
             dqgdT(c) = (1._r8 - frac_h2osfc(c))*hr*dqgdT(c)
          endif

          call QSat(t_h2osfc(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
          if (qsatg > forc_q(c) .and. forc_q(c) > qsatg) then
             qsatg = forc_q(c)
             qsatgdT = 0._r8
          end if
          qg_h2osfc(c) = qsatg
          dqgdT(c) = dqgdT(c) + frac_h2osfc(c) * qsatgdT

!          qg(c) = frac_sno(c)*qg_snow(c) + (1._r8 - frac_sno(c) - frac_h2osfc(c))*qg_soil(c) &
          qg(c) = frac_sno_eff(c)*qg_snow(c) + (1._r8 - frac_sno_eff(c) - frac_h2osfc(c))*qg_soil(c) &
               + frac_h2osfc(c) * qg_h2osfc(c)

       else
          call QSat(t_grnd(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)
          qg(c) = qred*qsatg
          dqgdT(c) = qred*qsatgdT

          if (qsatg > forc_q(c) .and. forc_q(c) > qred*qsatg) then
             qg(c) = forc_q(c)
             dqgdT(c) = 0._r8
          end if

          qg_snow(c) = qg(c)
          qg_soil(c) = qg(c)
          qg_h2osfc(c) = qg(c)
       endif

       ! Ground emissivity - only calculate for non-urban landunits 
       ! Urban emissivities are currently read in from data file

       if (ityplun(l) /= isturb) then
          if (ityplun(l)==istice .or. ityplun(l)==istice_mec) then
             emg(c) = 0.97_r8
          else
             emg(c) = (1._r8-frac_sno(c))*0.96_r8 + frac_sno(c)*0.97_r8
          end if
       end if

       ! Latent heat. We arbitrarily assume that the sublimation occurs
       ! only as h2osoi_liq = 0

       htvp(c) = hvap
       if (h2osoi_liq(c,snl(c)+1) <= 0._r8 .and. h2osoi_ice(c,snl(c)+1) > 0._r8) htvp(c) = hsub

       ! Switch between vaporization and sublimation causes rapid solution
       ! separation in perturbation growth test

#if (defined PERGRO)
       htvp(c) = hvap
#endif

       ! Ground roughness lengths over non-lake columns (includes bare ground, ground
       ! underneath canopy, wetlands, etc.)

       if (frac_sno(c) > 0._r8) then
          z0mg(c) = zsno
       else
          z0mg(c) = zlnd
       end if
       z0hg(c) = z0mg(c)            ! initial set only
       z0qg(c) = z0mg(c)            ! initial set only

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height

       beta(c) = 1._r8
       zii(c)  = 1000._r8
       thv(c)  = forc_th(c)*(1._r8+0.61_r8*forc_q(c))

    end do ! (end of columns loop)
    
    ! Initialization

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)

       ! Initial set (needed for history tape fields)

       eflx_sh_tot(p) = 0._r8
       l = plandunit(p)
       if (ityplun(l) == isturb) then
         eflx_sh_tot_u(p) = 0._r8
       else if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then 
         eflx_sh_tot_r(p) = 0._r8
       end if
       eflx_lh_tot(p) = 0._r8
       if (ityplun(l) == isturb) then
         eflx_lh_tot_u(p) = 0._r8
       else if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then 
         eflx_lh_tot_r(p) = 0._r8
       end if
       eflx_sh_veg(p) = 0._r8
       qflx_evap_tot(p) = 0._r8
       qflx_evap_veg(p) = 0._r8
       qflx_tran_veg(p) = 0._r8

       ! Initial set for calculation

       cgrnd(p)  = 0._r8
       cgrnds(p) = 0._r8
       cgrndl(p) = 0._r8

       ! Vegetation Emissivity

       avmuir = 1._r8
       emv(p) = 1._r8-exp(-(elai(p)+esai(p))/avmuir)

       ! Roughness lengths over vegetation

       z0m(p)    = z0mr(ivt(p)) * htop(p)
       displa(p) = displar(ivt(p)) * htop(p)

       z0mv(p)   = z0m(p)
       z0hv(p)   = z0mv(p)
       z0qv(p)   = z0mv(p)
    end do

    ! Make forcing height a pft-level quantity that is the atmospheric forcing 
    ! height plus each pft's z0m+displa
    do pi = 1,max_pft_per_gcell
       do g = lbg, ubg
          if (pi <= npfts(g)) then
            p = pfti(g) + pi - 1
            if (pactive(p)) then
              l = plandunit(p)
              c = pcolumn(p)
              if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
                if (frac_veg_nosno(p) == 0) then
                  forc_hgt_u_pft(p) = forc_hgt_u(g) + z0mg(c) + displa(p)
                  forc_hgt_t_pft(p) = forc_hgt_t(g) + z0mg(c) + displa(p)
                  forc_hgt_q_pft(p) = forc_hgt_q(g) + z0mg(c) + displa(p)
                else
                  forc_hgt_u_pft(p) = forc_hgt_u(g) + z0m(p) + displa(p)
                  forc_hgt_t_pft(p) = forc_hgt_t(g) + z0m(p) + displa(p)
                  forc_hgt_q_pft(p) = forc_hgt_q(g) + z0m(p) + displa(p)
                end if
              else if (ityplun(l) == istwet .or. ityplun(l) == istice      &
                                            .or. ityplun(l) == istice_mec) then
                forc_hgt_u_pft(p) = forc_hgt_u(g) + z0mg(c)
                forc_hgt_t_pft(p) = forc_hgt_t(g) + z0mg(c)
                forc_hgt_q_pft(p) = forc_hgt_q(g) + z0mg(c)
              ! Appropriate momentum roughness length will be added in SLakeFLuxesMod.
              else if (ityplun(l) == istdlak) then
                forc_hgt_u_pft(p) = forc_hgt_u(g)
                forc_hgt_t_pft(p) = forc_hgt_t(g)
                forc_hgt_q_pft(p) = forc_hgt_q(g)
              else if (ityplun(l) == isturb) then
                forc_hgt_u_pft(p) = forc_hgt_u(g) + z_0_town(l) + z_d_town(l)
                forc_hgt_t_pft(p) = forc_hgt_t(g) + z_0_town(l) + z_d_town(l)
                forc_hgt_q_pft(p) = forc_hgt_q(g) + z_0_town(l) + z_d_town(l)
              end if
            end if
          end if
       end do
    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       thm(p)  = forc_t(c) + 0.0098_r8*forc_hgt_t_pft(p)
    end do

  end subroutine Biogeophysics1

end module Biogeophysics1Mod
