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
    integer , pointer :: ivt(:)           !pft vegetation type
    integer , pointer :: ityplun(:)       !landunit type
    integer , pointer :: clandunit(:)     !column's landunit index
    integer , pointer :: cgridcell(:)     !column's gridcell index
    real(r8), pointer :: pwtgcell(:)      !weight relative to gridcell for each pft
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
!------------------------------------------------------------------------------

   ! Assign local pointers to derived type members (gridcell-level)

    forc_hgt_t    => clm_a2l%forc_hgt_t
    forc_u        => clm_a2l%forc_u
    forc_v        => clm_a2l%forc_v
    forc_hgt_u    => clm_a2l%forc_hgt_u
    forc_hgt_q    => clm_a2l%forc_hgt_q
    npfts         => grc%npfts
    pfti          => grc%pfti

    ! Assign local pointers to derived type members (landunit-level)

    ityplun       => lun%itype
    z_0_town      => lun%z_0_town
    z_d_town      => lun%z_d_town

    ! Assign local pointers to derived type members (column-level)

    forc_pbot     => cps%forc_pbot
    forc_q        => cws%forc_q
    forc_t        => ces%forc_t
    forc_th       => ces%forc_th

    cgridcell     => col%gridcell
    clandunit     => col%landunit
    ctype         => col%itype
    beta          => cps%beta
    dqgdT         => cws%dqgdT
    emg           => cps%emg
    frac_sno      => cps%frac_sno
    h2osno        => cws%h2osno
    htvp          => cps%htvp
    qg            => cws%qg
    smpmin        => cps%smpmin
    snl           => cps%snl
    t_grnd        => ces%t_grnd
    thv           => ces%thv
    z0hg          => cps%z0hg
    z0mg          => cps%z0mg
    z0qg          => cps%z0qg
    zii           => cps%zii
    bsw           => cps%bsw
    dz            => cps%dz
    h2osoi_ice    => cws%h2osoi_ice
    h2osoi_liq    => cws%h2osoi_liq
    soilalpha     => cws%soilalpha
    soilbeta      => cws%soilbeta
    soilalpha_u   => cws%soilalpha_u
    sucsat        => cps%sucsat
    t_soisno      => ces%t_soisno
    tssbef        => ces%tssbef
    watsat        => cps%watsat
    watfc         => cps%watfc
    watdry        => cps%watdry
    watopt        => cps%watopt
    rootfr_road_perv => cps%rootfr_road_perv
    rootr_road_perv  => cps%rootr_road_perv

    ! Assign local pointers to derived type members (pft-level)

    ivt           => pft%itype
    elai          => pps%elai
    esai          => pps%esai
    htop          => pps%htop
    emv           => pps%emv
    z0m           => pps%z0m
    displa        => pps%displa
    z0mv          => pps%z0mv
    z0hv          => pps%z0hv
    z0qv          => pps%z0qv
    eflx_sh_tot   => pef%eflx_sh_tot
    eflx_sh_tot_u => pef%eflx_sh_tot_u
    eflx_sh_tot_r => pef%eflx_sh_tot_r
    eflx_lh_tot   => pef%eflx_lh_tot
    eflx_lh_tot_u => pef%eflx_lh_tot_u
    eflx_lh_tot_r => pef%eflx_lh_tot_r
    eflx_sh_veg   => pef%eflx_sh_veg
    qflx_evap_tot => pwf%qflx_evap_tot
    qflx_evap_veg => pwf%qflx_evap_veg
    qflx_tran_veg => pwf%qflx_tran_veg
    cgrnd         => pef%cgrnd
    cgrnds        => pef%cgrnds
    cgrndl        => pef%cgrndl
    forc_hgt_u_pft => pps%forc_hgt_u_pft
    forc_hgt_t_pft => pps%forc_hgt_t_pft
    forc_hgt_q_pft => pps%forc_hgt_q_pft
    plandunit      => pft%landunit
    frac_veg_nosno => pps%frac_veg_nosno
    thm            => pes%thm
    pgridcell      => pft%gridcell
    pcolumn        => pft%column
    pwtgcell       => pft%wtgcell

    ! Assign local pointers to derived type members (ecophysiological)

    z0mr          => pftcon%z0mr
    displar       => pftcon%displar

    do j = -nlevsno+1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          tssbef(c,j) = t_soisno(c,j)
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

       t_grnd(c) = t_soisno(c,snl(c)+1)

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
             hr   = exp(psit/roverg/t_grnd(c))
             qred = (1.-frac_sno(c))*hr + frac_sno(c)

             !! Lee and Pielke 1992 beta, added by K.Sakaguchi
             if (wx < watfc(c,1) ) then  !when water content of ths top layer is less than that at F.C.
                fac_fc  = min(1._r8, wx/watfc(c,1))  !eqn5.66 but divided by theta at field capacity
                fac_fc  = max( fac_fc, 0.01_r8 )
                ! modifiy soil beta by snow cover. soilbeta for snow surface is one
                soilbeta(c) = (1._r8-frac_sno(c))*0.25_r8*(1._r8 - cos(SHR_CONST_PI*fac_fc))**2._r8 &
                              + frac_sno(c)
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

       call QSat(t_grnd(c), forc_pbot(c), eg, degdT, qsatg, qsatgdT)

       qg(c) = qred*qsatg
       dqgdT(c) = qred*qsatgdT

       if (qsatg > forc_q(c) .and. forc_q(c) > qred*qsatg) then
          qg(c) = forc_q(c)
          dqgdT(c) = 0._r8
       end if

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
            l = plandunit(p)
            ! Note: Some glacier_mec pfts may have zero weight  
            if (pwtgcell(p) > 0._r8 .or. ityplun(l)==istice_mec) then 
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
              else if (ityplun(l) == istdlak) then
                ! Should change the roughness lengths to shared constants
                if (t_grnd(c) >= tfrz) then
                  forc_hgt_u_pft(p) = forc_hgt_u(g) + 0.01_r8
                  forc_hgt_t_pft(p) = forc_hgt_t(g) + 0.01_r8
                  forc_hgt_q_pft(p) = forc_hgt_q(g) + 0.01_r8
                else
                  forc_hgt_u_pft(p) = forc_hgt_u(g) + 0.04_r8
                  forc_hgt_t_pft(p) = forc_hgt_t(g) + 0.04_r8
                  forc_hgt_q_pft(p) = forc_hgt_q(g) + 0.04_r8
                end if
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
