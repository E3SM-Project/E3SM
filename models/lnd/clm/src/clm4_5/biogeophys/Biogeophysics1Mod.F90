module Biogeophysics1Mod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs calculation of leaf temperature and surface fluxes.
  ! Biogeophysics2.F90 then determines soil/snow and ground
  ! temperatures and updates the surface fluxes for the new ground
  ! temperature.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only : iulog
  use decompMod   , only: bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Biogeophysics1   ! Calculate leaf temperature and surface fluxes
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine Biogeophysics1(bounds, num_nolakec, filter_nolakec, num_nolakep, filter_nolakep)
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
    use clm_atmlnd         , only : clm_a2l, a2l_downscaled_col
    use clm_varcon         , only : denh2o, denice, roverg, hvap, hsub, &
         istice, istice_mec, istwet, istsoil, istdlak, &
         zlnd, zsno, tfrz, icol_roof, icol_sunwall, icol_shadewall,     &
         icol_road_imperv, icol_road_perv, tfrz, spval, istdlak
    use clm_varcon         , only : istcrop
    use clm_varpar         , only : nlevgrnd, nlevurb, nlevsno, nlevsoi
    use QSatMod            , only : QSat
    use shr_const_mod      , only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds    ! bounds
    integer, intent(in) :: num_nolakec         ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:)   ! column filter for non-lake points
    integer, intent(in) :: num_nolakep         ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(:)   ! pft filter for non-lake points
    !
    ! !LOCAL VARIABLES:
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
    real(r8) :: fh2o_eff(bounds%begc:bounds%endc) ! effective surface water fraction (i.e. seen by atm)
    !------------------------------------------------------------------------------

   associate(& 
   frac_sno_eff              =>    cps%frac_sno_eff             , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_h2osfc               =>    cps%frac_h2osfc              , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc                    =>    cws%h2osfc                   , & ! Input:  [real(r8) (:)]  surface water (mm)                      
   t_h2osfc                  =>    ces%t_h2osfc                 , & ! Input:  [real(r8) (:)]  surface water temperature               
   t_h2osfc_bef              =>    ces%t_h2osfc_bef             , & ! Input:  [real(r8) (:)]  saved surface water temperature         
   qg_snow                   =>    cws%qg_snow                  , & ! Input:  [real(r8) (:)]  specific humidity at snow surface [kg/kg]
   qg_soil                   =>    cws%qg_soil                  , & ! Input:  [real(r8) (:)]  specific humidity at soil surface [kg/kg]
   qg_h2osfc                 =>    cws%qg_h2osfc                , & ! Input:  [real(r8) (:)]  specific humidity at h2osfc surface [kg/kg]
   forc_hgt_t                =>    clm_a2l%forc_hgt_t           , & ! Input:  [real(r8) (:)] observational height of temperature [m]  
   forc_u                    =>    clm_a2l%forc_u               , & ! Input:  [real(r8) (:)] atmospheric wind speed in east direction (m/s)
   forc_v                    =>    clm_a2l%forc_v               , & ! Input:  [real(r8) (:)] atmospheric wind speed in north direction (m/s)
   forc_hgt_u                =>    clm_a2l%forc_hgt_u           , & ! Input:  [real(r8) (:)] observational height of wind [m]         
   forc_hgt_q                =>    clm_a2l%forc_hgt_q           , & ! Input:  [real(r8) (:)] observational height of specific humidity [m]
   ityplun                   =>    lun%itype                    , & ! Input:  [integer (:)] landunit type                             
   urbpoi                    =>    lun%urbpoi                   , & ! Input:  [logical (:)]  true => landunit is an urban point       
   z_0_town                  =>   lun%z_0_town                  , & ! Input:  [real(r8) (:)] momentum roughness length of urban landunit (m)
   z_d_town                  =>   lun%z_d_town                  , & ! Input:  [real(r8) (:)] displacement height of urban landunit (m)
   forc_pbot                 =>    a2l_downscaled_col%forc_pbot , & ! Input:  [real(r8) (:)] atmospheric pressure (Pa)                
   forc_q                    =>    a2l_downscaled_col%forc_q    , & ! Input:  [real(r8) (:)] atmospheric specific humidity (kg/kg)    
   forc_t                    =>    a2l_downscaled_col%forc_t    , & ! Input:  [real(r8) (:)] atmospheric temperature (Kelvin)         
   forc_th                   =>    a2l_downscaled_col%forc_th   , & ! Input:  [real(r8) (:)] atmospheric potential temperature (Kelvin)
   cgridcell                 =>   col%gridcell                  , & ! Input:  [integer (:)] column's gridcell index                   
   clandunit                 =>   col%landunit                  , & ! Input:  [integer (:)] column's landunit index                   
   ctype                     =>    col%itype                    , & ! Input:  [integer (:)] column type                               
   beta                      =>    cps%beta                     , & ! Output: [real(r8) (:)] coefficient of convective velocity [-]   
   dqgdT                     =>    cws%dqgdT                    , & ! Output: [real(r8) (:)] d(qg)/dT                                 
   emg                       =>    cps%emg                      , & ! Output: [real(r8) (:)] ground emissivity                        
   frac_sno                  =>    cps%frac_sno                 , & ! Input:  [real(r8) (:)] fraction of ground covered by snow (0 to 1)
   h2osno                    =>    cws%h2osno                   , & ! Input:  [real(r8) (:)] snow water (mm H2O)                      
   htvp                      =>    cps%htvp                     , & ! Output: [real(r8) (:)] latent heat of vapor of water (or sublimation) [j/kg]
   qg                        =>    cws%qg                       , & ! Output: [real(r8) (:)] ground specific humidity [kg/kg]         
   smpmin                    =>    cps%smpmin                   , & ! Input:  [real(r8) (:)] restriction for min of soil potential (mm)
   snl                       =>    cps%snl                      , & ! Input:  [integer (:)] number of snow layers                     
   t_grnd                    =>    ces%t_grnd                   , & ! Output: [real(r8) (:)] ground temperature (Kelvin)              
   thv                       =>    ces%thv                      , & ! Output: [real(r8) (:)] virtual potential temperature (kelvin)   
   z0hg                      =>    cps%z0hg                     , & ! Output: [real(r8) (:)] roughness length over ground, sensible heat [m]
   z0mg                      =>    cps%z0mg                     , & ! Output: [real(r8) (:)] roughness length over ground, momentum [m]
   z0qg                      =>    cps%z0qg                     , & ! Output: [real(r8) (:)] roughness length over ground, latent heat [m]
   zii                       =>    cps%zii                      , & ! Output: [real(r8) (:)] convective boundary height [m]           
   bsw                       =>    cps%bsw                      , & ! Input:  [real(r8) (:,:)] Clapp and Hornberger "b"               
   dz                        =>    cps%dz                       , & ! Input:  [real(r8) (:,:)] layer depth (m)                        
   h2osoi_ice                =>    cws%h2osoi_ice               , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2)                       
   h2osoi_liq                =>    cws%h2osoi_liq               , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2)                   
   soilalpha                 =>    cws%soilalpha                , & ! Output: [real(r8) (:)] factor that reduces ground saturated specific humidity (-)
   soilbeta                  =>    cws%soilbeta                 , & ! Output: [real(r8) (:)] factor that reduces ground evaporation   
   soilalpha_u               =>    cws%soilalpha_u              , & ! Output: [real(r8) (:)] Urban factor that reduces ground saturated specific humidity (-)
   sucsat                    =>    cps%sucsat                   , & ! Input:  [real(r8) (:,:)] minimum soil suction (mm)              
   t_soisno                  =>    ces%t_soisno                 , & ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)              
   tssbef                    =>    ces%tssbef                   , & ! Output: [real(r8) (:,:)] soil/snow temperature before update    
   watsat                    =>    cps%watsat                   , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
   watfc                     =>    cps%watfc                    , & ! Input:  [real(r8) (:,:)] volumetric soil water at field capacity
   watdry                    =>    cps%watdry                   , & ! Input:  [real(r8) (:,:)] volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
   watopt                    =>    cps%watopt                   , & ! Input:  [real(r8) (:,:)] volumetric soil moisture corresponding to no restriction on ET from urban pervious surface
   rootfr_road_perv          =>    cps%rootfr_road_perv         , & ! Input:  [real(r8) (:,:)] fraction of roots in each soil layer for urban pervious road
   rootr_road_perv           =>    cps%rootr_road_perv          , & ! Input:  [real(r8) (:,:)] effective fraction of roots in each soil layer for urban pervious road
   pactive                   =>    pft%active                   , & ! Input:  [logical (:)] true=>do computations on this pft 
   ivt                       =>   pft%itype                     , & ! Input:  [integer (:)] pft vegetation type                       
   elai                      =>    pps%elai                     , & ! Input:  [real(r8) (:)] one-sided leaf area index with burying by snow
   esai                      =>    pps%esai                     , & ! Input:  [real(r8) (:)] one-sided stem area index with burying by snow
   htop                      =>    pps%htop                     , & ! Input:  [real(r8) (:)] canopy top (m)                           
   emv                       =>    pps%emv                      , & ! Output: [real(r8) (:)] vegetation emissivity                    
   z0m                       =>    pps%z0m                      , & ! Output: [real(r8) (:)] momentum roughness length (m)            
   displa                    =>    pps%displa                   , & ! Output: [real(r8) (:)] displacement height (m)                  
   z0mv                      =>    pps%z0mv                     , & ! Output: [real(r8) (:)] roughness length over vegetation, momentum [m]
   z0hv                      =>    pps%z0hv                     , & ! Output: [real(r8) (:)] roughness length over vegetation, sensible heat [m]
   z0qv                      =>    pps%z0qv                     , & ! Output: [real(r8) (:)] roughness length over vegetation, latent heat [m]
   eflx_sh_tot               =>    pef%eflx_sh_tot              , & ! Output: [real(r8) (:)] total sensible heat flux (W/m**2) [+ to atm]
   eflx_sh_tot_u             =>    pef%eflx_sh_tot_u            , & ! Output: [real(r8) (:)] urban total sensible heat flux (W/m**2) [+ to atm]
   eflx_sh_tot_r             =>    pef%eflx_sh_tot_r            , & ! Output: [real(r8) (:)] rural total sensible heat flux (W/m**2) [+ to atm]
   eflx_lh_tot               =>    pef%eflx_lh_tot              , & ! Output: [real(r8) (:)] total latent heat flux (W/m**2)  [+ to atm]
   eflx_lh_tot_u             =>    pef%eflx_lh_tot_u            , & ! Output: [real(r8) (:)] urban total latent heat flux (W/m**2)  [+ to atm]
   eflx_lh_tot_r             =>    pef%eflx_lh_tot_r            , & ! Output: [real(r8) (:)] rural total latent heat flux (W/m**2)  [+ to atm]
   eflx_sh_veg               =>    pef%eflx_sh_veg              , & ! Output: [real(r8) (:)] sensible heat flux from leaves (W/m**2) [+ to atm]
   qflx_evap_tot             =>    pwf%qflx_evap_tot            , & ! Output: [real(r8) (:)] qflx_evap_soi + qflx_evap_can + qflx_tran_veg
   qflx_evap_veg             =>    pwf%qflx_evap_veg            , & ! Output: [real(r8) (:)] vegetation evaporation (mm H2O/s) (+ = to atm)
   qflx_tran_veg             =>    pwf%qflx_tran_veg            , & ! Output: [real(r8) (:)] vegetation transpiration (mm H2O/s) (+ = to atm)
   cgrnd                     =>    pef%cgrnd                    , & ! Output: [real(r8) (:)] deriv. of soil energy flux wrt to soil temp [w/m2/k]
   cgrnds                    =>    pef%cgrnds                   , & ! Output: [real(r8) (:)] deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
   cgrndl                    =>    pef%cgrndl                   , & ! Output: [real(r8) (:)] deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
   forc_hgt_u_pft            =>    pps%forc_hgt_u_pft           , & ! Input:  [real(r8) (:)] observational height of wind at pft level [m]
   forc_hgt_t_pft            =>    pps%forc_hgt_t_pft           , & ! Input:  [real(r8) (:)] observational height of temperature at pft level [m]
   forc_hgt_q_pft            =>    pps%forc_hgt_q_pft           , & ! Input:  [real(r8) (:)] observational height of specific humidity at pft level [m]
   plandunit                 =>   pft%landunit                  , & ! Input:  [integer (:)] pft's landunit index                      
   frac_veg_nosno            =>    pps%frac_veg_nosno           , & ! Input:  [integer (:)] fraction of vegetation not covered by snow (0 OR 1) [-]
   thm                       =>    pes%thm                      , & ! Output: [real(r8) (:)] intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
   pgridcell                 =>   pft%gridcell                  , & ! Input:  [integer (:)] pft's gridcell index                      
   pcolumn                   =>   pft%column                    , & ! Input:  [integer (:)] pft's column index                        
   z0mr                      =>    pftcon%z0mr                  , & ! Input:  [real(r8) (:)] ratio of momentum roughness length to canopy top height (-)
   displar                   =>    pftcon%displar                 & ! Input:  [real(r8) (:)] ratio of displacement height to canopy top height (-)
   )

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

       if (.not. urbpoi(l)) then
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
       if (urbpoi(l)) then
         eflx_sh_tot_u(p) = 0._r8
       else if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then 
         eflx_sh_tot_r(p) = 0._r8
       end if
       eflx_lh_tot(p) = 0._r8
       if (urbpoi(l)) then
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
    do p = bounds%begp,bounds%endp
       if (pactive(p)) then
          g = pgridcell(p)
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
          else if (urbpoi(l)) then
             forc_hgt_u_pft(p) = forc_hgt_u(g) + z_0_town(l) + z_d_town(l)
             forc_hgt_t_pft(p) = forc_hgt_t(g) + z_0_town(l) + z_d_town(l)
             forc_hgt_q_pft(p) = forc_hgt_q(g) + z_0_town(l) + z_d_town(l)
          end if
       end if
    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       thm(p)  = forc_t(c) + 0.0098_r8*forc_hgt_t_pft(p)
    end do

    end associate 
   end subroutine Biogeophysics1

end module Biogeophysics1Mod
