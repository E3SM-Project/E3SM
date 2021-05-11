module LakeTemperatureMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates surface fluxes and temperature for lakes.
  ! Created by Zack Subin, 2009
  !
  ! !USES
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use CH4Mod            , only : ch4_type
  use EnergyFluxType    , only : energyflux_type
  use LakeStateType     , only : lakestate_type
  use SoilStateType     , only : soilstate_type
  use SolarAbsorbedType , only : solarabs_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ef, col_ws, col_wf  
  use VegetationType    , only : veg_pp
  use VegetationDataType, only : veg_ef  
  !    
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LakeTemperature              ! Calculates Lake related temperature
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: SoilThermProp_Lake   ! Set therm conductivities and heat cap of snow/soil layers
  private :: PhaseChange_Lake     ! Calculation of the phase change within snow/soil/lake layers
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LakeTemperature(bounds, num_lakec, filter_lakec, num_lakep, filter_lakep, &
       solarabs_vars, soilstate_vars, waterstate_vars, waterflux_vars, ch4_vars, &
       energyflux_vars, temperature_vars, lakestate_vars)
    !
    ! !DESCRIPTION:
    ! Calculates temperatures in the 25-45 layer column of (possible) snow,
    ! lake water, soil, and bedrock beneath lake.
    ! Snow and soil temperatures are determined as in Laketemperature, except
    ! for appropriate boundary conditions at the top of the snow (the flux is fixed
    ! to be the ground heat flux calculated in LakeFluxes), the bottom of the snow
    ! (adjacent to top lake layer), and the top of the soil (adjacent to the bottom
    ! lake layer). Also, the soil is kept fully saturated.
    ! The whole column is solved simultaneously as one tridiagonal matrix.
    ! Major changes from CLM4:
    !    i) Lake water layers can freeze by any fraction and release latent heat; thermal
    !       and mechanical properties are adjusted for ice fraction.
    !   ii) Convective mixing (though not eddy diffusion) still occurs for frozen lakes.
    !  iii) No sunlight is absorbed in the lake if there are snow layers (except for that allowed through
    !                                                                     to the top "soil" layer by SNICAR)
    !   iv) Light is allowed to reach the top soil layer (where it is assumed to be completely absorbed).
    !    v) Lakes have variable depth, set ultimately in surface data set but now in initLakeMod.
    !   vi) The extinction coefficient is now variable with depth.
    !   vii) The fraction of shortwave absorbed at the surface is now the NIR fraction, rather than a fixed parameter.
    !   viii) Enhanced background diffusion and option for increased mixing for deep lakes is added.
    !   See discussion in Subin et al. 2011
    !   
    !   Lakes are allowed to have variable depth, set in initLakeMod.
    !
    ! Use the Crank-Nicholson method to set up tridiagonal system of equations to
    ! solve for ts at time n+1, where the temperature equation for layer i is
    ! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
    !
    ! The solution conserves energy as:
    !
    !
    ! cv*([ts(top soisno layer)] n+1 - [ts(top soisno layer)] n)*dz(top soisno layer)/dt + ... +
    ! cv*([ts(nlevlak+nlevgrnd)] n+1 - [ts(nlevlak+nlevgrnd)] n)*dz(nlevlak+nlevgrnd)/dt = eflx_soil_grnd
    ! cv is not constant.
    ! ts is with respect to freezing temperature, as there is a discontinuity in cv at this temperature.
    !
    ! where:
    ! [ts] n   = old temperature (kelvin)
    ! [ts] n+1 = new temperature (kelvin)
    ! eflx_soil_grnd      = energy flux into lake (w/m**2)
    !          = betaprime*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot
    !          (This is now the same as the ground heat flux.)
    !            + phi(1) + ... + phi(nlevlak) + phi(top soil level)
    ! betaprime = beta (NIR fraction) for no snow layers, and sabg_lyr(p,jtop)/sabg(p) for snow layers.
    !
    ! WARNING: This subroutine assumes lake columns have one and only one pft.
    !
    ! Outline:
    ! 1  ) Initialization
    ! 2  ) Lake density
    ! 3  ) Diffusivity
    ! 4  ) Heat source term from solar radiation penetrating lake
    ! 5  ) Set thermal props and find initial energy content
    ! 6  ) Set up vectors for tridiagonal matrix solution
    ! 7  ) Solve tridiagonal and back-substitute
    ! 8  ) (Optional) Do first energy check using temperature change at constant heat capacity.
    ! 9  ) Phase change
    ! 9.5) (Optional) Do second energy check using temperature change and latent heat, 
    !      considering changed heat capacity. Also do soil water balance check.
    ! 10 ) Convective mixing 
    ! 11 ) Do final energy check to detect small numerical errors (especially from convection)
    !      and dump small imbalance into sensible heat, or pass large errors to BalanceCheckMod for abort.
    !
    ! !USES:
    use LakeCon           , only : betavis, za_lake, n2min, tdmax, pudz, depthcrit, mixfact
    use LakeCon           , only : lakepuddling, lake_no_ed
    use QSatMod            , only : QSat
    use TridiagonalMod     , only : Tridiagonal
    use elm_varpar         , only : nlevlak, nlevgrnd, nlevsno
    use clm_time_manager   , only : get_step_size
    use elm_varcon         , only : hfus, cpliq, cpice, tkwat, tkice, denice
    use elm_varcon         , only : vkc, grav, denh2o, tfrz, cnfac
    use elm_varctl         , only : iulog, use_lch4
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds	  
    integer                , intent(in)    :: num_lakec       ! number of column non-lake points in column filter
    integer                , intent(in)    :: filter_lakec(:) ! column filter for non-lake points
    integer                , intent(in)    :: num_lakep       ! number of column non-lake points in pft filter
    integer                , intent(in)    :: filter_lakep(:) ! patch filter for non-lake points
    type(solarabs_type)    , intent(in)    :: solarabs_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    type(waterflux_type)   , intent(inout) :: waterflux_vars
    type(ch4_type)         , intent(inout) :: ch4_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(temperature_type) , intent(inout) :: temperature_vars
    type(lakestate_type)   , intent(inout) :: lakestate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: p0 = 1._r8                                      ! neutral value of turbulent prandtl number
    integer  :: i,j,fc,fp,g,c,p                                            ! do loop or array index
    real(r8) :: dtime                                                      ! land model time step (sec)
    real(r8) :: beta(bounds%begc:bounds%endc)                              ! fraction of solar rad absorbed at surface: equal to NIR fraction
                                                                           ! of surface absorbed shortwave
    real(r8) :: eta                                                        ! light extinction coefficient (/m): depends on lake type
    real(r8) :: cwat                                                       ! specific heat capacity of water (j/m**3/kelvin)
    real(r8) :: cice_eff                                                   ! effective heat capacity of ice (using density of
                                                                           ! water because layer depth is not adjusted when freezing
    real(r8) :: cfus                                                       ! effective heat of fusion per unit volume
                                                                           ! using water density as above
    real(r8) :: km                                                         ! molecular diffusion coefficient (m**2/s)
    real(r8) :: tkice_eff                                                  ! effective conductivity since layer depth is constant
    real(r8) :: a(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)     ! "a" vector for tridiagonal matrix
    real(r8) :: b(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)     ! "b" vector for tridiagonal matrix
    real(r8) :: c1(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)    ! "c" vector for tridiagonal matrix
    real(r8) :: r(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)     ! "r" vector for tridiagonal solution
    real(r8) :: rhow(bounds%begc:bounds%endc,nlevlak)                      ! density of water (kg/m**3)
    real(r8) :: phi(bounds%begc:bounds%endc,nlevlak)                       ! solar radiation absorbed by layer (w/m**2)
    real(r8) :: kme(bounds%begc:bounds%endc,nlevlak)                       ! molecular + eddy diffusion coefficient (m**2/s)
    real(r8) :: rsfin                                                      ! relative flux of solar radiation into layer
    real(r8) :: rsfout                                                     ! relative flux of solar radiation out of layer
    real(r8) :: phi_soil(bounds%begc:bounds%endc)                          ! solar radiation into top soil layer (W/m**2)
    real(r8) :: ri                                                         ! richardson number
    real(r8) :: fin(bounds%begc:bounds%endc)                               ! net heat flux into lake at ground interface (w/m**2)
    real(r8) :: ocvts(bounds%begc:bounds%endc)                             ! (cwat*(t_lake[n  ])*dz
    real(r8) :: ncvts(bounds%begc:bounds%endc)                             ! (cwat*(t_lake[n+1])*dz
    real(r8) :: ke                                                         ! eddy diffusion coefficient (m**2/s)
    real(r8) :: zin                                                        ! depth at top of layer (m)
    real(r8) :: zout                                                       ! depth at bottom of layer (m)
    real(r8) :: drhodz                                                     ! d [rhow] /dz (kg/m**4)
    real(r8) :: n2                                                         ! brunt-vaisala frequency (/s**2)
    real(r8) :: num                                                        ! used in calculating ri
    real(r8) :: den                                                        ! used in calculating ri
    real(r8) :: tav_froz(bounds%begc:bounds%endc)                          ! used in aver temp for convectively mixed layers (C)
    real(r8) :: tav_unfr(bounds%begc:bounds%endc)                          ! "
    real(r8) :: nav(bounds%begc:bounds%endc)                               ! used in aver temp for convectively mixed layers
    real(r8) :: phidum                                                     ! temporary value of phi
    real(r8) :: iceav(bounds%begc:bounds%endc)                             ! used in calc aver ice for convectively mixed layers
    real(r8) :: qav(bounds%begc:bounds%endc)                               ! used in calc aver heat content for conv. mixed layers
    integer  :: jtop(bounds%begc:bounds%endc)                              ! top level for each column (no longer all 1)
    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)           ! heat capacity of soil/snow [J/(m2 K)]
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)           ! thermal conductivity of soil/snow [W/(m K)]
                                                                           ! (at interface below, except for j=0)
    real(r8) :: cv_lake (bounds%begc:bounds%endc,1:nlevlak)                ! heat capacity [J/(m2 K)]
    real(r8) :: tk_lake (bounds%begc:bounds%endc,1:nlevlak)                ! thermal conductivity at layer node [W/(m K)]
    real(r8) :: cvx (bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)  ! heat capacity for whole column [J/(m2 K)]
    real(r8) :: tkix(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)  ! thermal conductivity at layer interfaces
                                                                           ! for whole column [W/(m K)]
    real(r8) :: tx(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)    ! temperature of whole column [K]
    real(r8) :: tktopsoillay(bounds%begc:bounds%endc)                      ! thermal conductivity [W/(m K)]
    real(r8) :: fnx(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)   ! heat diffusion through the layer interface below [W/m2]
    real(r8) :: phix(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)  ! solar source term for whole column [W/m**2]
    real(r8) :: zx(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd)    ! interface depth (+ below surface) for whole column [m]
    real(r8) :: dzm                                                        ! used in computing tridiagonal matrix [m]
    real(r8) :: dzp                                                        ! used in computing tridiagonal matrix [m]
    integer  :: jprime                                                     ! j - nlevlak
    real(r8) :: factx(bounds%begc:bounds%endc,-nlevsno+1:nlevlak+nlevgrnd) ! coefficient used in computing tridiagonal matrix
    real(r8) :: t_lake_bef(bounds%begc:bounds%endc,1:nlevlak)              ! beginning lake temp for energy conservation check [K]
    real(r8) :: t_soisno_bef(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)  ! beginning soil temp for E cons. check [K]
    real(r8) :: lhabs(bounds%begc:bounds%endc)                             ! total per-column latent heat abs. from phase change  (J/m^2)
    real(r8) :: esum1(bounds%begc:bounds%endc)                             ! temp for checking energy (J/m^2)
    real(r8) :: esum2(bounds%begc:bounds%endc)                             ! ""
    real(r8) :: zsum(bounds%begc:bounds%endc)                              ! temp for putting ice at the top during convection (m)
    real(r8) :: wsum(bounds%begc:bounds%endc)                              ! temp for checking water (kg/m^2)
    real(r8) :: wsum_end(bounds%begc:bounds%endc)                          ! temp for checking water (kg/m^2)
    real(r8) :: sabg_col(bounds%begc:bounds%endc)                          ! absorbed ground solar for column (W/m^2)
    real(r8) :: sabg_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:1)         ! absorbed ground solar in layer for column (W/m^2)
    real(r8) :: sabg_nir                                                   ! NIR that is absorbed (W/m^2)
    integer  :: jconvect(bounds%begc:bounds%endc)                          ! Lowest level where convection occurs
    integer  :: jconvectbot(bounds%begc:bounds%endc)                       ! Hightest level where bottom-originating convection occurs
    logical  :: bottomconvect(bounds%begc:bounds%endc)                     ! Convection originating in bottom layer of lake triggers special convection loop
    real(r8) :: fangkm                                                     ! (m^2/s) extra diffusivity based on Fang & Stefan 1996, citing Ellis, 1991

    ! They think that mixing energy will generally get into lake to make
    ! diffusivity exceed molecular; the energy is damped out according to the Brunt-Vaisala
    ! frequency, yielding a maximum diffusivity for neutral stability of about 6 times molecular
    ! For puddling
    logical  :: puddle(bounds%begc:bounds%endc)
    real(r8) :: icesum(bounds%begc:bounds%endc) ! m
    logical  :: frzn(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    associate(                                                       & 
         dz_lake         =>   col_pp%dz_lake                          , & ! Input:  [real(r8) (:,:) ]  layer thickness for lake (m)          
         z_lake          =>   col_pp%z_lake                           , & ! Input:  [real(r8) (:,:) ]  layer depth for lake (m)              
         dz              =>   col_pp%dz                               , & ! Input:  [real(r8) (:,:) ]  layer thickness for snow & soil (m)   
         z               =>   col_pp%z                                , & ! Input:  [real(r8) (:,:) ]  layer depth for snow & soil (m)       
         snl             =>   col_pp%snl                              , & ! Input:  [integer  (:)   ]  negative of number of snow layers        
         lakedepth       =>   col_pp%lakedepth                        , & ! Input:  [real(r8) (:)   ]  column lake depth (m)                   

         sabg            =>   solarabs_vars%sabg_patch             , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)
         sabg_lyr        =>   solarabs_vars%sabg_lyr_patch         , & ! Input:  [real(r8) (:,:) ]  absorbed solar radiation (pft,lyr) [W/m2]
         fsds_nir_d      =>   solarabs_vars%fsds_nir_d_patch       , & ! Input:  [real(r8) (:)   ]  incident direct beam nir solar radiation (W/m**2)
         fsds_nir_i      =>   solarabs_vars%fsds_nir_i_patch       , & ! Input:  [real(r8) (:)   ]  incident diffuse nir solar radiation (W/m**2)
         fsr_nir_d       =>   solarabs_vars%fsr_nir_d_patch        , & ! Input:  [real(r8) (:)   ]  reflected direct beam nir solar radiation (W/m**2)
         fsr_nir_i       =>   solarabs_vars%fsr_nir_i_patch        , & ! Input:  [real(r8) (:)   ]  reflected diffuse nir solar radiation (W/m**2)

         etal            =>   lakestate_vars%etal_col              , & ! Input:  [real(r8) (:)   ]  extinction coefficient from surface data (1/m)
         ks              =>   lakestate_vars%ks_col                , & ! Input:  [real(r8) (:)   ]  coefficient passed to LakeTemperature            
         ws              =>   lakestate_vars%ws_col                , & ! Input:  [real(r8) (:)   ]  surface friction velocity (m/s)                   
         lake_raw       =>    lakestate_vars%lake_raw_col          , & ! Input:  [real(r8) (:)   ]  aerodynamic resistance for moisture (s/m)   
         
         h2osno          =>   col_ws%h2osno           , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                     
         h2osoi_liq      =>   col_ws%h2osoi_liq       , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) [for snow & soil layers]
         h2osoi_ice      =>   col_ws%h2osoi_ice       , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2) [for snow & soil layers]
         frac_iceold     =>   col_ws%frac_iceold      , & ! Output: [real(r8) (:,:) ]  fraction of ice relative to the tot water

         qflx_snofrz_col =>   col_wf%qflx_snofrz       , & ! Output: [real(r8) (:)   ]  column-integrated snow freezing rate (kg m-2 s-1) [+]

         t_grnd          =>   col_es%t_grnd          , & ! Input:  [real(r8) (:)   ]  ground temperature (Kelvin)             
         t_soisno        =>   col_es%t_soisno        , & ! Output: [real(r8) (:,:) ]  soil (or snow) temperature (Kelvin)   
         t_lake          =>   col_es%t_lake          , & ! Output: [real(r8) (:,:) ]  col lake temperature (Kelvin)             
         hc_soi          =>   col_es%hc_soi          , & ! Output: [real(r8) (:)   ]  soil heat content (MJ/m2)               
         hc_soisno       =>   col_es%hc_soisno       , & ! Output: [real(r8) (:)   ]  soil plus snow plus lake heat content (MJ/m2)

         beta            =>   lakestate_vars%betaprime_col         , & ! Output: [real(r8) (:)   ]  col effective beta: sabg_lyr(p,jtop) for snow layers, beta otherwise
         lake_icefrac    =>   lakestate_vars%lake_icefrac_col      , & ! Output: [real(r8) (:,:) ]  col mass fraction of lake layer that is frozen
         lake_icethick   =>   lakestate_vars%lake_icethick_col     , & ! Output: [real(r8) (:)   ]  col ice thickness (m) (integrated if lakepuddling)
         savedtke1       =>   lakestate_vars%savedtke1_col         , & ! Output: [real(r8) (:)   ]  col top level eddy conductivity (W/mK)      
         lakeresist      =>   lakestate_vars%lakeresist_col        , & ! Output: [real(r8) (:)   ]  col (Needed for calc. of grnd_ch4_cond) (s/m) 

         grnd_ch4_cond   =>   ch4_vars%grnd_ch4_cond_col           , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s] (only over lake points)

         eflx_soil_grnd  =>   veg_ef%eflx_soil_grnd , & ! Output: [real(r8) (:)   ]  heat flux into snow / lake (W/m**2) [+ = into soil]
         eflx_sh_grnd    =>   veg_ef%eflx_sh_grnd   , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot     =>   veg_ef%eflx_sh_tot    , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_gnet       =>   veg_ef%eflx_gnet      , & ! Output: [real(r8) ( :)  ]  net heat flux into ground (W/m**2) at the surface interface
         errsoi          =>   col_ef%errsoi             & ! Output: [real(r8) (:)   ]  soil/lake energy conservation error (W/m**2)
         )

    ! 1!) Initialization
    ! Determine step size

    dtime = get_step_size()

    ! Initialize constants
    cwat = cpliq*denh2o ! water heat capacity per unit volume
    cice_eff = cpice*denh2o !use water density because layer depth is not adjusted
                              !for freezing
    cfus = hfus*denh2o  ! latent heat per unit volume
    tkice_eff = tkice * denice/denh2o !effective conductivity since layer depth is constant
    km = tkwat/cwat     ! a constant (molecular diffusivity)

    ! Needed for Lahey compiler which doesn't seem to allow shortcircuit logic for undefined variables.
    puddle(bounds%begc:bounds%endc) = .false.
    frzn(bounds%begc:bounds%endc) = .false.

    ! Begin calculations

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       ! Initialize quantities needed below

       ocvts(c) = 0._r8
       ncvts(c) = 0._r8
       esum1(c) = 0._r8
       esum2(c) = 0._r8
       hc_soisno(c) = 0._r8
       hc_soi(c)    = 0._r8
       if (use_lch4) then
          jconvect(c) = 0
          jconvectbot(c) = nlevlak+1
          lakeresist(c) = 0._r8
       end if
       bottomconvect(bounds%begc:bounds%endc) = .false.

       qflx_snofrz_col(c) = 0._r8
    end do

    ! Initialize set of previous time-step variables as in DriverInit,
    ! which is currently not called over lakes. This has to be done
    ! here because phase change will occur in this routine.
    ! Ice fraction of snow at previous time step

    do j = -nlevsno+1,0
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         if (j >= snl(c) + 1) then
            frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
         end if
      end do
    end do

    ! Prepare for lake layer temperature calculations below

    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = veg_pp%column(p)

       ! fin(c) = betaprime * sabg(p) + forc_lwrad(c) - (eflx_lwrad_out(p) + &
       !     eflx_sh_tot(p) + eflx_lh_tot(p)) 
       ! fin(c) now passed from LakeFluxes as eflx_gnet
       fin(c) = eflx_gnet(p)

       ! Calculate the NIR fraction of absorbed solar. This will now be the "beta" parameter.
       ! Total NIR absorbed:
       sabg_nir = fsds_nir_d(p) + fsds_nir_i(p) - fsr_nir_d(p) - fsr_nir_i(p)
       sabg_nir = min(sabg_nir, sabg(p))
       beta(c) = sabg_nir/max(1.e-5_r8,sabg(p))
       beta(c) = beta(c) + (1._r8-beta(c))*betavis

    end do

    ! 2!) Lake density

    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * & 
                      1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-tdmax))**1.68_r8 ) &
                    + lake_icefrac(c,j)*denice
                    ! Allow for ice fraction; assume constant ice density.
                    ! This is not the correct average-weighting but that's OK because the density will only
                    ! be used for convection for lakes with ice, and the ice fraction will dominate the
                    ! density differences between layers.
                    ! Using this average will make sure that surface ice is treated properly during
                    ! convective mixing.
       end do
    end do

    ! 3!) Diffusivity and implied thermal "conductivity" = diffusivity * cwat
    do j = 1, nlevlak-1
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          drhodz = (rhow(c,j+1)-rhow(c,j)) / (z_lake(c,j+1)-z_lake(c,j))
          n2 = grav / rhow(c,j) * drhodz

          ! Fixed sign error here: our z goes up going down into the lake, so no negative
          ! sign is needed to make this positive unlike in Hostetler. --ZS
          num = 40._r8 * n2 * (vkc*z_lake(c,j))**2._r8
          den = max( (ws(c)**2._r8) * exp(-2._r8*ks(c)*z_lake(c,j)), 1.e-10_r8 )
          ri = ( -1._r8 + sqrt( max(1._r8+num/den, 0._r8) ) ) / 20._r8

          if (lakepuddling .and. j == 1) frzn(c) = .false.

          if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0 .and. &
              (.not. lakepuddling .or. (lake_icefrac(c,j) == 0._r8 .and. .not. frzn(c))) ) then
             ke = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
             kme(c,j) = km + ke

             if (.not. lake_no_ed) then
                fangkm = 1.039e-8_r8 * max(n2,n2min)**(-0.43_r8)  ! Fang & Stefan 1996, citing Ellis et al 1991
                kme(c,j) = kme(c,j) + fangkm
             end if
             if (lakedepth(c) >= depthcrit) then
                kme(c,j) = kme(c,j) * mixfact
             end if

             tk_lake(c,j) = kme(c,j)*cwat
          else
             kme(c,j) = km
             if (.not. lake_no_ed) then
                fangkm = 1.039e-8_r8 * max(n2,n2min)**(-0.43_r8)
                kme(c,j) = kme(c,j) + fangkm
                if (lakedepth(c) >= depthcrit) then
                   kme(c,j) = kme(c,j) * mixfact
                end if
                tk_lake(c,j) = kme(c,j)*cwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                               + kme(c,j)*cwat*lake_icefrac(c,j) )
             else
                tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                               + tkwat*lake_icefrac(c,j) )
             ! Assume the resistances add as for the calculation of conductivities at layer interfaces.
             end if
             if (lakepuddling) frzn(c) = .true.
                ! Prevent eddy mixing beneath frozen layers even when surface is unfrozen.
          end if
       end do
    end do

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       j = nlevlak
       kme(c,nlevlak) = kme(c,nlevlak-1)

       if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0 .and. &
           (.not. lakepuddling .or. (lake_icefrac(c,j) == 0._r8 .and. .not. frzn(c)) ) ) then
          tk_lake(c,j) = tk_lake(c,j-1)
       else
          if (.not. lake_no_ed) then
             tk_lake(c,j) = kme(c,j)*cwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                               + kme(c,j)*cwat*lake_icefrac(c,j) )
          else
             tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                               + tkwat*lake_icefrac(c,j) )
          end if
       end if

       ! Use in surface flux calculation for next timestep.
       savedtke1(c) = kme(c,1)*cwat ! Will only be used if unfrozen
       ! set number of column levels for use by Tridiagonal below
       jtop(c) = snl(c) + 1
    end do

    ! 4!) Heat source term
    do j = 1, nlevlak
       do fp = 1, num_lakep
          p = filter_lakep(fp)
          c = veg_pp%column(p)

          ! If no eta from surface data,
          ! Set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
          ! (regression of Secchi Depth with lake depth for small glacial basin lakes), and the
          ! Poole & Atkins expression for extinction coeffient of 1.7 / Secchi Depth (m).
          if (etal(c) > 0._r8) then ! use eta from surface data
             eta = etal(c)
          else
             eta = 1.1925_r8*max(lakedepth(c),1._r8)**(-0.424_r8)
          end if

          zin  = z_lake(c,j) - 0.5_r8*dz_lake(c,j)
          zout = z_lake(c,j) + 0.5_r8*dz_lake(c,j)
          rsfin  = exp( -eta*max(  zin-za_lake,0._r8 ) )
          rsfout = exp( -eta*max( zout-za_lake,0._r8 ) )

          ! Let rsfout for bottom layer go into soil.
          ! This looks like it should be robust even for pathological cases,
            ! like lakes thinner than za_lake.

          if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
             phidum = (rsfin-rsfout) * sabg(p) * (1._r8-beta(c))
             if (j == nlevlak) then
                phi_soil(c) = rsfout * sabg(p) * (1._r8-beta(c))
             end if
          else if (j == 1 .and. snl(c) == 0) then !if frozen but no snow layers
             phidum = sabg(p) * (1._r8-beta(c))
             ! This should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
             ! cause spring convection.
          else if (j == 1) then
             phidum = sabg_lyr(p,j)
          !some radiation absorbed in snow layers, the rest in the top layer of lake
          !radiation absorbed in snow layers will be applied below
          else
             phidum = 0._r8
             if (j == nlevlak) phi_soil(c) = 0._r8
          end if
          phi(c,j) = phidum

       end do
    end do

    ! 5!) Set thermal properties and check initial energy content.

    ! For lake
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
       end do
    end do

    ! For snow / soil
    call SoilThermProp_Lake(bounds, num_lakec, filter_lakec, &
         tk(bounds%begc:bounds%endc, :), &
         cv(bounds%begc:bounds%endc, :), &
         tktopsoillay(bounds%begc:bounds%endc), &
         soilstate_vars, waterstate_vars, temperature_vars)
    
    ! Sum cv*t_lake for energy check
    ! Include latent heat term, and use tfrz as reference temperature
    ! to prevent abrupt change in heat content due to changing heat capacity with phase change.

    ! This will need to be over all soil / lake / snow layers. Lake is below.
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          ocvts(c) = ocvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
                   + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) 
          t_lake_bef(c,j) = t_lake(c,j)
       end do
    end do

    ! Now do for soil / snow layers
    do j = -nlevsno + 1, nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if (j >= jtop(c)) then
             ocvts(c) = ocvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                      + hfus*h2osoi_liq(c,j)
             if (j == 1 .and. h2osno(c) > 0._r8 .and. j == jtop(c)) then
                ocvts(c) = ocvts(c) - h2osno(c)*hfus
             end if
             t_soisno_bef(c,j) = t_soisno(c,j)
          end if
       end do
    end do

    !!!!!!!!!!!!!!!!!!!
    ! 6!) Set up vector r and vectors a, b, c1 that define tridiagonal matrix

    ! Heat capacity and resistance of snow without snow layers (<1cm) is ignored during diffusion,
    ! but its capacity to absorb latent heat may be used during phase change.

    ! Transfer sabg and sabg_lyr to column level
    do j = -nlevsno+1,1
       do fp = 1, num_lakep
          p = filter_lakep(fp)
          c = veg_pp%column(p)

          if (j >= jtop(c)) then
             if (j == jtop(c)) sabg_col(c) = sabg(p)
             sabg_lyr_col(c,j) = sabg_lyr(p,j)
          end if
       end do
    end do

    ! Set up interface depths, zx, heat capacities, cvx, solar source terms, phix, and temperatures, tx.
    do j = -nlevsno+1, nlevlak+nlevgrnd
       do fc = 1,num_lakec
          c = filter_lakec(fc)

          jprime = j - nlevlak

          if (j >= jtop(c)) then
             if (j < 1) then !snow layer
                zx(c,j) = z(c,j)
                cvx(c,j) = cv(c,j)
                if (j == jtop(c)) then ! no absorption because it has already been assigned to the surface
                                       ! interface
                   phix(c,j) = 0._r8
                else
                   phix(c,j) = sabg_lyr_col(c,j) !New for SNICAR
                end if
                tx(c,j) = t_soisno(c,j)
             else if (j <= nlevlak) then !lake layer
                zx(c,j) = z_lake(c,j)
                cvx(c,j) = cv_lake(c,j)
                phix(c,j) = phi(c,j)
                tx(c,j) = t_lake(c,j)
             else !soil layer
                zx(c,j) = zx(c,nlevlak) + dz_lake(c,nlevlak)/2._r8 + z(c,jprime)
                cvx(c,j) = cv(c,jprime)
                if (j == nlevlak + 1) then !top soil layer
                   phix(c,j) = phi_soil(c)
                else !middle or bottom soil layer
                   phix(c,j) = 0._r8
                end if
                tx(c,j) = t_soisno(c,jprime)
             end if
          end if

       end do
    end do

    ! Determine interface thermal conductivities, tkix

    do j = -nlevsno+1, nlevlak+nlevgrnd
       do fc = 1,num_lakec
          c = filter_lakec(fc)

          jprime = j - nlevlak

          if (j >= jtop(c)) then
             if (j < 0) then !non-bottom snow layer
                tkix(c,j) = tk(c,j)
             else if (j == 0) then !bottom snow layer
                dzp = zx(c,j+1) - zx(c,j)
                tkix(c,j) = tk_lake(c,1)*tk(c,j)*dzp / &
                      (tk(c,j)*z_lake(c,1) + tk_lake(c,1)*(-z(c,j)) )
                ! tk(c,0) is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
             else if (j < nlevlak) then !non-bottom lake layer
                tkix(c,j) = ( tk_lake(c,j)*tk_lake(c,j+1) * (dz_lake(c,j+1)+dz_lake(c,j)) ) &
                           / ( tk_lake(c,j)*dz_lake(c,j+1) + tk_lake(c,j+1)*dz_lake(c,j) )
             else if (j == nlevlak) then !bottom lake layer
                dzp = zx(c,j+1) - zx(c,j)
                tkix(c,j) = (tktopsoillay(c)*tk_lake(c,j)*dzp / &
                    (tktopsoillay(c)*dz_lake(c,j)/2._r8 + tk_lake(c,j)*z(c,1) ) )
                    ! tktopsoillay is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
             else !soil layer
                tkix(c,j) = tk(c,jprime)
             end if
         end if

      end do 
   end do


    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c1 that define tridiagonal
    ! matrix and solve system

    do j = -nlevsno+1, nlevlak+nlevgrnd
       do fc = 1,num_lakec
          c = filter_lakec(fc)
          if (j >= jtop(c)) then
             if (j < nlevlak+nlevgrnd) then !top or interior layer
                factx(c,j) = dtime/cvx(c,j)
                fnx(c,j) = tkix(c,j)*(tx(c,j+1)-tx(c,j))/(zx(c,j+1)-zx(c,j))
             else !bottom soil layer
                factx(c,j) = dtime/cvx(c,j)
                fnx(c,j) = 0._r8 !not used
             end if
          end if
       enddo
    end do

    do j = -nlevsno+1,nlevlak+nlevgrnd
       do fc = 1,num_lakec
          c = filter_lakec(fc)
          if (j >= jtop(c)) then
             if (j == jtop(c)) then !top layer
                dzp    = zx(c,j+1)-zx(c,j)
                a(c,j) = 0._r8
                b(c,j) = 1+(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
                c1(c,j) =  -(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
                r(c,j) = tx(c,j) + factx(c,j)*( fin(c) + phix(c,j) + cnfac*fnx(c,j) )
             else if (j < nlevlak+nlevgrnd) then !middle layer
                dzm    = (zx(c,j)-zx(c,j-1))
                dzp    = (zx(c,j+1)-zx(c,j))
                a(c,j) =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j-1)/dzm
                b(c,j) = 1._r8+ (1._r8-cnfac)*factx(c,j)*(tkix(c,j)/dzp + tkix(c,j-1)/dzm)
                c1(c,j) =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j)/dzp
                r(c,j) = tx(c,j) + cnfac*factx(c,j)*( fnx(c,j) - fnx(c,j-1) ) + factx(c,j)*phix(c,j)
             else  !bottom soil layer
                dzm     = (zx(c,j)-zx(c,j-1))
                a(c,j) =   - (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
                b(c,j) = 1._r8+ (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
                c1(c,j) = 0._r8
                r(c,j) = tx(c,j) - cnfac*factx(c,j)*fnx(c,j-1)
             end if
          end if
       enddo
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! 7!) Solve for tdsolution

    call Tridiagonal(bounds, -nlevsno + 1, nlevlak + nlevgrnd, &
         jtop(bounds%begc:bounds%endc), &
         num_lakec, filter_lakec, &
         a(bounds%begc:bounds%endc, :), &
         b(bounds%begc:bounds%endc, :), &
         c1(bounds%begc:bounds%endc, :), &
         r(bounds%begc:bounds%endc, :), &
         tx(bounds%begc:bounds%endc, :))
 
    ! Set t_soisno and t_lake
    do j = -nlevsno+1, nlevlak + nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          jprime = j - nlevlak

          ! Don't do anything with invalid snow layers.
          if (j >= jtop(c)) then
             if (j < 1) then !snow layer
                t_soisno(c,j) = tx(c,j)
             else if (j <= nlevlak) then !lake layer
                t_lake(c,j)   = tx(c,j)
             else !soil layer
                t_soisno(c,jprime) = tx(c,j)
             end if
          end if
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!

    ! 8!) Sum energy content and total energy into lake for energy check. Any errors will be from the
    !     Tridiagonal solution.
    ! This section was for debugging only and has been removed.  See original "ch4" branch code.

    !!!!!!!!!!!!!!!!!!!!!!!

    ! 9!) Phase change
    call PhaseChange_Lake(bounds, num_lakec, filter_lakec, &
         cv(bounds%begc:bounds%endc, :), &
         cv_lake(bounds%begc:bounds%endc, :), &
         lhabs(bounds%begc:bounds%endc), &
         waterstate_vars, waterflux_vars, temperature_vars, &
         energyflux_vars, lakestate_vars)

    !!!!!!!!!!!!!!!!!!!!!!!

    ! 9.5!) Second energy check and water check.  Now check energy balance before and after phase
    !       change, considering the possibility of changed heat capacity during phase change, by
    !       using initial heat capacity in the first step, final heat capacity in the second step,
    !       and differences from tfrz only to avoid enthalpy correction for (cpliq-cpice)*melt*tfrz.
    !       Also check soil water sum.
    ! This section was for debugging only and has been removed.  See original "ch4" branch code.

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 10!) Convective mixing: make sure fracice*dz is conserved, heat content c*dz*T is conserved, and
    ! all ice ends up at the top. Done over all lakes even if frozen.
    ! Either an unstable density profile or ice in a layer below an incompletely frozen layer will trigger.

    !Recalculate density
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * &
                      1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-tdmax))**1.68_r8 ) &
                    + lake_icefrac(c,j)*denice
       end do
    end do

    if (lakepuddling) then
    ! For sensitivity tests
       do j = 1, nlevlak
          do fc = 1, num_lakec
             c = filter_lakec(fc)
   
             if (j == 1) then
                icesum(c) = 0._r8
                puddle(c) = .false.
             end if
   
             icesum(c) = icesum(c) + lake_icefrac(c,j)*dz(c,j)
   
             if (j == nlevlak) then
                if (icesum(c) >= pudz) puddle(c) = .true.
             end if
          end do
       end do
    end if

    ! ZMS: The algorithm here, from CLM2-4 and presumably Hostetler, starts at the top and proceeds down;
    ! each time it finds an unstable density profile the lake is mixed from this point to the top fully.
    ! Occasionally in deep temperate lakes, the unstable profile can originate at the bottom because of small amounts of
    ! heat coming from the sediments. To prevent an unrealistic complete overturning, convection starting in the bottom
    ! layer is treated separately, mixing from the bottom up only one level at a time until
    ! a stable density profile is attained, rather than mixing all the way to the top immediately.

    ! First examine top nlevlak-1 layers.
    do j = 1, nlevlak-2
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          qav(c) = 0._r8
          nav(c) = 0._r8
          iceav(c) = 0._r8
       end do

       do i = 1, j+1
          do fc = 1, num_lakec
             c = filter_lakec(fc)
             if ( (.not. lakepuddling .or. .not. puddle(c) ) .and. (rhow(c,j) > rhow(c,j+1) .or. &
                (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) ) then
                qav(c) = qav(c) + dz_lake(c,i)*(t_lake(c,i)-tfrz) * & 
                        ((1._r8 - lake_icefrac(c,i))*cwat + lake_icefrac(c,i)*cice_eff)
                !tav(c) = tav(c) + t_lake(c,i)*dz_lake(c,i)
                iceav(c) = iceav(c) + lake_icefrac(c,i)*dz_lake(c,i)
                nav(c) = nav(c) + dz_lake(c,i)
                if (use_lch4) then
                   jconvect(c) = j+1
                end if
             end if
          end do
       end do

       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if ( (.not. lakepuddling .or. .not. puddle(c) ) .and. (rhow(c,j) > rhow(c,j+1) .or. &
                (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) ) then
             qav(c) = qav(c)/nav(c)
             iceav(c) = iceav(c)/nav(c)
             !If the average temperature is above freezing, put the extra energy into the water.
             !If it is below freezing, take it away from the ice.
             if (qav(c) > 0._r8) then
                tav_froz(c) = 0._r8 !Celsius
                tav_unfr(c) = qav(c) / ((1._r8 - iceav(c))*cwat)
             else if (qav(c) < 0._r8) then
                tav_froz(c) = qav(c) / (iceav(c)*cice_eff)
                tav_unfr(c) = 0._r8 !Celsius
             else
                tav_froz(c) = 0._r8
                tav_unfr(c) = 0._r8
             end if
          end if
       end do

       do i = 1, j+1
          do fc = 1, num_lakec
             c = filter_lakec(fc)
             if (nav(c) > 0._r8) then

                !Put all the ice at the top.!
                !If the average temperature is above freezing, put the extra energy into the water.
                !If it is below freezing, take it away from the ice.
                !For the layer with both ice & water, be careful to use the average temperature
                !that preserves the correct total heat content given what the heat capacity of that
                !layer will actually be.
                if (i == 1) zsum(c) = 0._r8
                if ((zsum(c)+dz_lake(c,i))/nav(c) <= iceav(c)) then
                   lake_icefrac(c,i) = 1._r8
                   t_lake(c,i) = tav_froz(c) + tfrz
                else if (zsum(c)/nav(c) < iceav(c)) then
                   lake_icefrac(c,i) = (iceav(c)*nav(c) - zsum(c)) / dz_lake(c,i)
                   ! Find average value that preserves correct heat content.
                   t_lake(c,i) = ( lake_icefrac(c,i)*tav_froz(c)*cice_eff &
                               + (1._r8 - lake_icefrac(c,i))*tav_unfr(c)*cwat ) &
                               / ( lake_icefrac(c,i)*cice_eff + (1-lake_icefrac(c,i))*cwat ) + tfrz
                else
                   lake_icefrac(c,i) = 0._r8
                   t_lake(c,i) = tav_unfr(c) + tfrz
                end if
                zsum(c) = zsum(c) + dz_lake(c,i)

                rhow(c,i) = (1._r8 - lake_icefrac(c,i)) * & 
                            1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,i)-tdmax))**1.68_r8 ) &
                          + lake_icefrac(c,i)*denice
             end if
          end do
       end do
    end do

    ! Now check bottom layer
    j = nlevlak-1
    do fc = 1, num_lakec
       c = filter_lakec(fc)

       if ( (.not. lakepuddling .or. .not. puddle(c) ) .and. (rhow(c,j) > rhow(c,j+1) .or. &
            (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) ) then
          ! convection originating in bottom layer. Could be coming from sediments-- be careful not to
          ! unnecessarily mix all the way to the top of the lake
          bottomconvect(c) = .true.
       end if
    end do

    ! Start mixing from bottom up. Only mix as high as the unstable density profile persists.
    do j = nlevlak-1, 1, -1
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          qav(c) = 0._r8
          nav(c) = 0._r8
          iceav(c) = 0._r8
       end do

       do i = j, nlevlak
          do fc = 1, num_lakec
             c = filter_lakec(fc)
             if ( bottomconvect(c) .and. &
                  (.not. lakepuddling .or. .not. puddle(c) ) .and. (rhow(c,j) > rhow(c,j+1) .or. &
                  (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) ) then
                qav(c) = qav(c) + dz_lake(c,i)*(t_lake(c,i)-tfrz) * & 
                        ((1._r8 - lake_icefrac(c,i))*cwat + lake_icefrac(c,i)*cice_eff)
                !tav(c) = tav(c) + t_lake(c,i)*dz_lake(c,i)
                iceav(c) = iceav(c) + lake_icefrac(c,i)*dz_lake(c,i)
                nav(c) = nav(c) + dz_lake(c,i)
                if (use_lch4) then
                   jconvectbot(c) = j
                end if
             end if
          end do
       end do

       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if ( bottomconvect(c) .and. &
               (.not. lakepuddling .or. .not. puddle(c) ) .and. (rhow(c,j) > rhow(c,j+1) .or. &
               (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) ) then
             qav(c) = qav(c)/nav(c)
             iceav(c) = iceav(c)/nav(c)
             !If the average temperature is above freezing, put the extra energy into the water.
             !If it is below freezing, take it away from the ice.
             if (qav(c) > 0._r8) then
                tav_froz(c) = 0._r8 !Celsius
                tav_unfr(c) = qav(c) / ((1._r8 - iceav(c))*cwat)
             else if (qav(c) < 0._r8) then
                tav_froz(c) = qav(c) / (iceav(c)*cice_eff)
                tav_unfr(c) = 0._r8 !Celsius
             else
                tav_froz(c) = 0._r8
                tav_unfr(c) = 0._r8
             end if
          end if
       end do

       do i = j, nlevlak
          do fc = 1, num_lakec
             c = filter_lakec(fc)
             if (bottomconvect(c) .and. nav(c) > 0._r8) then

                !Put all the ice at the top.!
                !If the average temperature is above freezing, put the extra energy into the water.
                !If it is below freezing, take it away from the ice.
                !For the layer with both ice & water, be careful to use the average temperature
                !that preserves the correct total heat content given what the heat capacity of that
                !layer will actually be.

                if (i == j) zsum(c) = 0._r8
                if ((zsum(c)+dz_lake(c,i))/nav(c) <= iceav(c)) then
                   lake_icefrac(c,i) = 1._r8
                   t_lake(c,i) = tav_froz(c) + tfrz
                else if (zsum(c)/nav(c) < iceav(c)) then
                   lake_icefrac(c,i) = (iceav(c)*nav(c) - zsum(c)) / dz_lake(c,i)
                   ! Find average value that preserves correct heat content.
                   t_lake(c,i) = ( lake_icefrac(c,i)*tav_froz(c)*cice_eff &
                               + (1._r8 - lake_icefrac(c,i))*tav_unfr(c)*cwat ) &
                               / ( lake_icefrac(c,i)*cice_eff + (1-lake_icefrac(c,i))*cwat ) + tfrz
                else
                   lake_icefrac(c,i) = 0._r8
                   t_lake(c,i) = tav_unfr(c) + tfrz
                end if
                zsum(c) = zsum(c) + dz_lake(c,i)

                rhow(c,i) = (1._r8 - lake_icefrac(c,i)) * & 
                            1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,i)-tdmax))**1.68_r8 ) &
                          + lake_icefrac(c,i)*denice
             end if
          end do
       end do
    end do

    ! Calculate lakeresist and grnd_ch4_cond for CH4 Module
    ! The CH4 will diffuse directly from the top soil layer to the atmosphere, so 
    ! the whole lake resistance is included.

    if (use_lch4) then
       do j = 1, nlevlak
          do fc = 1, num_lakec
             c = filter_lakec(fc)

             if (j > jconvect(c) .and. j < jconvectbot(c)) then  ! Assume resistance is zero for levels that convect
                lakeresist(c) = lakeresist(c) + dz(c,j)/kme(c,j) ! dz/eddy or molecular diffusivity
             end if

             if (j == nlevlak) then ! Calculate grnd_ch4_cond
                grnd_ch4_cond(c) = 1._r8 / (lakeresist(c) + lake_raw(c))

                ! Lake water R +  aerodynamic R
                ! Snow will be considered in methane routine
                ! No methane conduction through frozen lake
                if (lake_icefrac(c,1) > 0.1_r8) grnd_ch4_cond(c) = 0._r8
             end if

          end do
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    ! 11!) Re-evaluate thermal properties and sum energy content.
    ! For lake
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
       end do
    end do

    ! For snow / soil
    call SoilThermProp_Lake(bounds, num_lakec, filter_lakec, &
         tk(bounds%begc:bounds%endc, :), &
         cv(bounds%begc:bounds%endc, :), &
         tktopsoillay(bounds%begc:bounds%endc), &
         soilstate_vars, waterstate_vars, temperature_vars)


    ! Do as above to sum energy content
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          ncvts(c) = ncvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
                   + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) 
          fin(c) = fin(c) + phi(c,j)
          ! New for CLM 4
          hc_soisno(c) = hc_soisno(c) + cv_lake(c,j)*t_lake(c,j)/1.e6
       end do
    end do

    do j = -nlevsno + 1, nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if (j >= jtop(c)) then
             ncvts(c) = ncvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                      + hfus*h2osoi_liq(c,j) 
             if (j < 1) fin(c) = fin(c) + phix(c,j) !For SNICAR
             if (j == 1 .and. h2osno(c) > 0._r8 .and. j == jtop(c)) then
                ncvts(c) = ncvts(c) - h2osno(c)*hfus
             end if
             hc_soisno(c) = hc_soisno(c) + cv(c,j)*t_soisno(c,j)/1.e6
             if (j >= 1) hc_soi(c) = hc_soi(c) + cv(c,j)*t_soisno(c,j)/1.e6
          end if
          if (j == 1) fin(c) = fin(c) + phi_soil(c)
       end do
    end do


    ! Check energy conservation.
    do fp = 1, num_lakep
       p = filter_lakep(fp)
       c = veg_pp%column(p)
       errsoi(c) = (ncvts(c)-ocvts(c)) / dtime - fin(c)
       if (abs(errsoi(c)) < 0.10_r8) then ! else send to Balance Check and abort
          eflx_sh_tot(p)    = eflx_sh_tot(p)    - errsoi(c)
          eflx_sh_grnd(p)   = eflx_sh_grnd(p)   - errsoi(c)
          eflx_soil_grnd(p) = eflx_soil_grnd(p) + errsoi(c)
          eflx_gnet(p)      = eflx_gnet(p)      + errsoi(c)
          if (abs(errsoi(c)) > 1.e-3_r8) then
             write(iulog,*)'errsoi incorporated into sensible heat in LakeTemperature: c, (W/m^2):', c, errsoi(c)
          end if
          errsoi(c) = 0._r8
       end if

    end do
    ! This loop assumes only one point per column.

    ! lake_icethick diagnostic.
    do j = 1, nlevlak
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          if (j == 1) lake_icethick(c) = 0._r8

          lake_icethick(c) = lake_icethick(c) + lake_icefrac(c,j)*dz_lake(c,j)*denh2o/denice
                                                           ! Convert from nominal to physical thickness
       end do
    end do

    end associate 
   end subroutine LakeTemperature

   !-----------------------------------------------------------------------
   subroutine SoilThermProp_Lake (bounds,  num_lakec, filter_lakec, tk, cv, tktopsoillay, &
        soilstate_vars, waterstate_vars, temperature_vars)
     !
     ! !DESCRIPTION:
     ! Calculation of thermal conductivities and heat capacities of
     ! snow/soil layers
     ! (1) The volumetric heat capacity is calculated as a linear combination
     !     in terms of the volumetric fraction of the constituent phases.
     !
     ! (2) The thermal conductivity of soil is computed from the algorithm of
     !     Johansen (as reported by Farouki 1981), and of snow is from the
     !     formulation used in SNTHERM (Jordan 1991).
     ! The thermal conductivities at the interfaces between two neighboring
     ! layers (j, j+1) are derived from an assumption that the flux across
     ! the interface is equal to that from the node j to the interface and the
     ! flux from the interface to the node j+1.
     !
     ! For lakes, the proper soil layers (not snow) should always be saturated.
     !
     ! !USES:
     use elm_varcon  , only : denh2o, denice, tfrz, tkwat, tkice, tkair
     use elm_varcon  , only : cpice,  cpliq, thk_bedrock
     use elm_varpar  , only : nlevsno, nlevsoi, nlevgrnd
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)  :: bounds  
     integer                , intent(in)  :: num_lakec                        ! number of column lake points in column filter
     integer                , intent(in)  :: filter_lakec(:)                  ! column filter for lake points
     real(r8)               , intent(out) :: cv( bounds%begc: , -nlevsno+1: ) ! heat capacity [J/(m2 K)] [col, lev]
     real(r8)               , intent(out) :: tk( bounds%begc: , -nlevsno+1: ) ! thermal conductivity [W/(m K)] [col, lev]
     real(r8)               , intent(out) :: tktopsoillay( bounds%begc: )     ! thermal conductivity [W/(m K)] [col]
     type(soilstate_type)   , intent(in)  :: soilstate_vars
     type(waterstate_type)  , intent(in)  :: waterstate_vars
     type(temperature_type) , intent(in)  :: temperature_vars
    
     !
     ! !LOCAL VARIABLES:
     integer  :: l,c,j                     ! indices
     integer  :: fc                        ! lake filtered column indices
     real(r8) :: bw                        ! partial density of water (ice + liquid)
     real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
     real(r8) :: dke                       ! kersten number
     real(r8) :: fl                        ! fraction of liquid or unfrozen water to total water
     real(r8) :: satw                      ! relative total water content of soil.
     real(r8) :: thk(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd) ! thermal conductivity of layer
     real(r8) :: xicevol                   ! (virtual excess ice volume per nominal soil volume)
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(cv)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(tk)           == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(tktopsoillay) == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

     associate(                                           & 
          snl         => col_pp%snl                        , & ! Input:  [integer (:)]  number of snow layers                    
          dz          => col_pp%dz                         , & ! Input:  [real(r8) (:,:)]  layer thickness (m)                   
          zi          => col_pp%zi                         , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m) 
          z           => col_pp%z                          , & ! Input:  [real(r8) (:,:)]  layer depth (m)                       

          watsat      => soilstate_vars%watsat_col      , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
          tksatu      => soilstate_vars%tksatu_col      , & ! Input:  [real(r8) (:,:)]  thermal conductivity, saturated soil [W/m-K]
          tkmg        => soilstate_vars%tkmg_col        , & ! Input:  [real(r8) (:,:)]  thermal conductivity, soil minerals  [W/m-K]
          tkdry       => soilstate_vars%tkdry_col       , & ! Input:  [real(r8) (:,:)]  thermal conductivity, dry soil (W/m/Kelvin)
          csol        => soilstate_vars%csol_col        , & ! Input:  [real(r8) (:,:)]  heat capacity, soil solids (J/m**3/Kelvin)

          h2osoi_liq  => col_ws%h2osoi_liq , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                  
          h2osoi_ice  => col_ws%h2osoi_ice , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                      

          t_soisno    => col_es%t_soisno    & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)             
          )

       ! Thermal conductivity of soil from Farouki (1981)

       do j = -nlevsno+1,nlevgrnd
          do fc = 1, num_lakec
             c = filter_lakec(fc)

             ! Only examine levels from 1->nlevsoi
             if (j >= 1 .and. j <= nlevsoi) then
                !             l = clandunit(c)
                !             if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
                ! This could be altered later for allowing this to be over glaciers.

                ! Soil should be saturated in LakeHydrology
                satw = 1._r8
                ! ZMS: Note the following needs to be updated for the corrections to be merged into SoilTemperature.
                ! They are especially important here because of no supercooled water.
                fl = h2osoi_liq(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j))
                if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                   dke = max(0._r8, log10(satw) + 1.0_r8)
                   dksat = tksatu(c,j)
                else                               ! Frozen soil
                   dke = satw
                   dksat = tkmg(c,j)*0.249_r8**(fl*watsat(c,j))*2.29_r8**watsat(c,j)
                endif
                thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                ! Hydrology routine won't let the excess be liquid.
                if (satw > 1._r8) then
                   xicevol = (satw-1._r8)*watsat(c,j)
                   thk(c,j) = (thk(c,j) + xicevol*tkice)/(1._r8 + xicevol)/(1._r8 + xicevol)
                   ! Use simple mean because excess ice is likely to be in chunks, thus conductivities add
                   ! rather than the usual addition of resistances.
                   ! Conductivity is reduced by the extra virtual volume fraction, as dz is not changing.
                end if
                ! See discussion in LakeHydrology.  This is the simplest way to treat thermal cycling near
                ! freezing even if not modeling excess ice, to assume up to 10% excess ice over watsat after refreezing.
                ! Assume some ground heaving rather than drainage of the water while it's freezing.
                ! This allows for energy conservation and for the final allowed ice volume to be independent
                ! of the timestep or precise freezing trajectory.
                ! Does real sediment under lakes heave or drain? Could revisit later.
             elseif (j > nlevsoi) then
                thk(c,j) = thk_bedrock
             endif

             ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
             ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
             if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
                bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dz(c,j)
                thk(c,j) = tkair + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkair)
             end if

          end do
       end do

       ! Thermal conductivity at the layer interface

       ! Have to correct for the fact that bottom snow layer and top soil layer border lake.
       ! For the first case, the snow layer conductivity for the middle of the layer will be returned.
       ! Because the interfaces are below the soil layers, the conductivity for the top soil layer
       ! will have to be returned separately.
       do j = -nlevsno+1,nlevgrnd
          do fc = 1,num_lakec
             c = filter_lakec(fc)
             if (j >= snl(c)+1 .AND. j <= nlevgrnd-1 .AND. j /= 0) then
                tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                     /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
             else if (j == 0 .and. j >= snl(c)+1) then
                tk(c,j) = thk(c,j)
             else if (j == nlevgrnd) then
                tk(c,j) = 0._r8
             end if
             ! For top soil layer.
             if (j == 1) tktopsoillay(c) = thk(c,j)
          end do
       end do

       ! Soil heat capacity, from de Vires (1963)

       do j = 1, nlevgrnd
          do fc = 1,num_lakec
             c = filter_lakec(fc)
             cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) +   &
                  (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
             !   if (j == 1) then
             !      if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
             !         cv(c,j) = cv(c,j) + cpice*h2osno(c)
             !      end if
             !   end if
             ! Won't worry about heat capacity for thin snow on lake with no snow layers.
             ! Its temperature will be assigned based on air temperature anyway if a new node is formed.
          enddo
       end do

       ! Snow heat capacity

       do j = -nlevsno+1,0
          do fc = 1,num_lakec
             c = filter_lakec(fc)
             if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
                cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
             end if
          end do
       end do

     end associate

   end subroutine SoilThermProp_Lake

   !-----------------------------------------------------------------------
   subroutine PhaseChange_Lake (bounds, num_lakec, filter_lakec, cv, cv_lake, lhabs, &
        waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars, lakestate_vars)
     !
     ! !DESCRIPTION:
     ! Calculation of the phase change within snow, soil, & lake layers:
     ! (1) Check the conditions for which the phase change may take place,
     !     i.e., the layer temperature is great than the freezing point
     !     and the ice mass is not equal to zero (i.e. melting),
     !     or the layer temperature is less than the freezing point
     !     and the liquid water mass is greater than zero (i.e. melting). No supercooling.
     ! (2) Assess the amount of phase change from the energy excess (or deficit)
     !     after setting the layer temperature to freezing point, depending on
     !     how much water or ice is available.
     ! (3) Re-adjust the ice and liquid mass, and the layer temperature: either to
     !     the freezing point if enough water or ice is available to fully compensate,
     !     or to a remaining temperature.
     ! The specific heats are adjusted during phase change for precise energy conservation.
     ! Errors will be trapped at the end of LakeTemperature.
     !
     ! !USES:
     use clm_time_manager , only : get_step_size
     use elm_varcon       , only : tfrz, hfus, denh2o, denice, cpliq, cpice
     use elm_varpar       , only : nlevsno, nlevgrnd, nlevlak
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds                     
     integer                , intent(in)    :: num_lakec                        ! number of lake columns
     integer                , intent(in)    :: filter_lakec(:)                  ! column filter for lake points
     real(r8)               , intent(inout) :: cv( bounds%begc: , -nlevsno+1: ) ! heat capacity [J/(m2 K)] [col, lev]
     real(r8)               , intent(inout) :: cv_lake( bounds%begc: , 1: )     ! heat capacity [J/(m2 K)] [col, levlak]
     real(r8)               , intent(out)   :: lhabs( bounds%begc: )            ! total per-column latent heat abs. (J/m^2) [col]
     type(waterstate_type)  , intent(inout) :: waterstate_vars
     type(waterflux_type)   , intent(inout) :: waterflux_vars
     type(temperature_type) , intent(inout) :: temperature_vars
     type(energyflux_type)  , intent(inout) :: energyflux_vars
     type(lakestate_type)   , intent(inout) :: lakestate_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: j,c,g                              ! do loop index
     integer  :: fc                                 ! lake filtered column indices
     real(r8) :: dtime                              ! land model time step (sec)
     real(r8) :: heatavail                          ! available energy for melting or freezing (J/m^2)
     real(r8) :: heatrem                            ! energy residual or loss after melting or freezing
     real(r8) :: melt                               ! actual melting (+) or freezing (-) [kg/m2]
     !real(r8), parameter :: smallnumber = 1.e-7_r8 ! to prevent tiny residuals from rounding error
     real(r8), parameter :: smallnumber = 1.e-12_r8 ! The above actually was enough to cause a 0.1 W/m^2 energy imbalance
     ! when the bottom lake layer started freezing in a 50m Arctic lake
     logical  :: dophasechangeflag
     !-----------------------------------------------------------------------

     ! Enforce expected array sizes
     SHR_ASSERT_ALL((ubound(cv)      == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(cv_lake) == (/bounds%endc, nlevlak/)),  errMsg(__FILE__, __LINE__))
     SHR_ASSERT_ALL((ubound(lhabs)   == (/bounds%endc/)),           errMsg(__FILE__, __LINE__))

     associate(                                                   & 
          dz_lake         => col_pp%dz_lake                        , & ! Input:  [real(r8)  (:,:) ] lake layer thickness (m)              
          dz              => col_pp%dz                             , & ! Input:  [real(r8)  (:,:) ] layer thickness (m)                   
          snl             => col_pp%snl                            , & ! Input:  [integer   (:)   ] number of snow layers                    

          snow_depth      => col_ws%snow_depth     , & ! Output: [real(r8)  (:)   ] snow height (m)                         
          h2osno          => col_ws%h2osno         , & ! Output: [real(r8)  (:)   ] snow water (mm H2O)                     
          h2osoi_liq      => col_ws%h2osoi_liq     , & ! Output: [real(r8)  (:,:) ] liquid water (kg/m2)                  
          h2osoi_ice      => col_ws%h2osoi_ice     , & ! Output: [real(r8)  (:,:) ] ice lens (kg/m2)                      

          lake_icefrac    => lakestate_vars%lake_icefrac_col    , & ! Input:  [real(r8)  (:,:) ] mass fraction of lake layer that is frozen
          
          qflx_snofrz_lyr => col_wf%qflx_snofrz_lyr , & ! Input:  [real(r8)  (:,:) ] snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
          qflx_snow_melt  => col_wf%qflx_snow_melt  , & ! Output: [real(r8)  (:)   ] net snow melt                           
          qflx_snomelt    => col_wf%qflx_snomelt    , & ! Output: [real(r8)  (:)   ] snow melt (mm H2O /s)                   
          qflx_snofrz_col => col_wf%qflx_snofrz     , & ! Output: [real(r8)  (:)   ] column-integrated snow freezing rate (kg m-2 s-1) [+]

          t_soisno        => col_es%t_soisno      , & ! Input:  [real(r8)  (:,:) ] soil temperature (Kelvin)             
          t_lake          => col_es%t_lake        , & ! Input:  [real(r8)  (:,:) ] lake temperature (Kelvin)             
          imelt           => col_ef%imelt         , & ! Output: [integer   (:,:) ] flag for melting (=1), freezing (=2), Not=0 (new)

          eflx_snomelt    => col_ef%eflx_snomelt     & ! Output: [real(r8)  (:)   ] snow melt heat flux (W/m**2)            
          )

       ! Get step size

       dtime = get_step_size()

       ! Initialization

       do fc = 1,num_lakec
          c = filter_lakec(fc)

          qflx_snomelt(c)   = 0._r8
          eflx_snomelt(c)   = 0._r8
          lhabs(c)          = 0._r8
          qflx_snow_melt(c) = 0._r8
       end do

       do j = -nlevsno+1,0
          do fc = 1,num_lakec
             c = filter_lakec(fc)

             qflx_snofrz_lyr(c,j) = 0._r8
             ! Do for all possible snow layers in case snl changes over timestep. ! Bug corrected ZMS 10/14/11
             imelt(c,j) = 0
          end do
       end do

       ! Check for case of snow without snow layers and top lake layer temp above freezing.

       do fc = 1,num_lakec
          c = filter_lakec(fc)

          if (snl(c) == 0 .and. h2osno(c) > 0._r8 .and. t_lake(c,1) > tfrz) then
             heatavail = (t_lake(c,1) - tfrz) * cv_lake(c,1)
             melt = min(h2osno(c), heatavail/hfus)
             heatrem = max(heatavail - melt*hfus, 0._r8)
             !catch small negative value to keep t at tfrz
             t_lake(c,1) = tfrz + heatrem/(cv_lake(c,1))
             snow_depth(c) = snow_depth(c)*(1._r8 - melt/h2osno(c))
             h2osno(c) = h2osno(c) - melt
             lhabs(c) = lhabs(c) + melt*hfus
             qflx_snomelt(c)   = qflx_snomelt(c)   + melt/dtime
             qflx_snow_melt(c) = qflx_snow_melt(c) + melt/dtime
             ! Prevent tiny residuals
             if (h2osno(c) < smallnumber) h2osno(c) = 0._r8
             if (snow_depth(c) < smallnumber) snow_depth(c) = 0._r8
          end if
       end do

       ! Lake phase change

       do j = 1,nlevlak
          do fc = 1,num_lakec
             c = filter_lakec(fc)

             dophasechangeflag = .false.
             if (t_lake(c,j) > tfrz .and. lake_icefrac(c,j) > 0._r8) then ! melting
                dophasechangeflag = .true.
                heatavail = (t_lake(c,j) - tfrz) * cv_lake(c,j)
                melt = min(lake_icefrac(c,j)*denh2o*dz_lake(c,j), heatavail/hfus)
                !denh2o is used because layer thickness is not adjusted for freezing
                heatrem = max(heatavail - melt*hfus, 0._r8)
                !catch small negative value to keep t at tfrz
             else if (t_lake(c,j) < tfrz .and. lake_icefrac(c,j) < 1._r8) then !freezing
                dophasechangeflag = .true.
                heatavail = (t_lake(c,j) - tfrz) * cv_lake(c,j)
                melt = max(-(1._r8-lake_icefrac(c,j))*denh2o*dz_lake(c,j), heatavail/hfus)
                !denh2o is used because layer thickness is not adjusted for freezing
                heatrem = min(heatavail - melt*hfus, 0._r8)
                !catch small positive value to keep t at tfrz
             end if
             ! Update temperature and ice fraction.
             if (dophasechangeflag) then
                lake_icefrac(c,j) = lake_icefrac(c,j) - melt/(denh2o*dz_lake(c,j))
                lhabs(c) = lhabs(c) + melt*hfus
                ! Update heat capacity
                cv_lake(c,j) = cv_lake(c,j) + melt*(cpliq-cpice)
                t_lake(c,j) = tfrz + heatrem/cv_lake(c,j)
                ! Prevent tiny residuals
                if (lake_icefrac(c,j) > 1._r8 - smallnumber) lake_icefrac(c,j) = 1._r8
                if (lake_icefrac(c,j) < smallnumber)         lake_icefrac(c,j) = 0._r8
             end if
          end do
       end do

       ! Snow & soil phase change
       ! Currently, does not do freezing point depression.

       do j = -nlevsno+1,nlevgrnd
          do fc = 1,num_lakec
             c = filter_lakec(fc)
             dophasechangeflag = .false.

             if (j >= snl(c) + 1) then

                if (t_soisno(c,j) > tfrz .and. h2osoi_ice(c,j) > 0._r8) then ! melting
                   dophasechangeflag = .true.
                   heatavail = (t_soisno(c,j) - tfrz) * cv(c,j)
                   melt = min(h2osoi_ice(c,j), heatavail/hfus)
                   heatrem = max(heatavail - melt*hfus, 0._r8)
                   !catch small negative value to keep t at tfrz
                   if (j <= 0) then !snow
                      imelt(c,j) = 1
                      qflx_snomelt(c) = qflx_snomelt(c) + melt/dtime
                   end if
                else if (t_soisno(c,j) < tfrz .and. h2osoi_liq(c,j) > 0._r8) then !freezing
                   dophasechangeflag = .true.
                   heatavail = (t_soisno(c,j) - tfrz) * cv(c,j)
                   melt = max(-h2osoi_liq(c,j), heatavail/hfus)
                   heatrem = min(heatavail - melt*hfus, 0._r8)
                   !catch small positive value to keep t at tfrz
                   if (j <= 0) then !snow
                      imelt(c,j) = 2
                      !qflx_snomelt(c) = qflx_snomelt(c) + melt
                      ! Does this works for both signs of melt in SnowHydrology? I think
                      ! qflx_snomelt(c) is just output.
                      ! It looks like qflx_snomelt is just supposed to be positive.
                      ! New variable for CLM 4
                      qflx_snofrz_lyr(c,j) = -melt/dtime
                   end if
                end if

                ! Update temperature and soil components.
                if (dophasechangeflag) then
                   h2osoi_ice(c,j) = h2osoi_ice(c,j) - melt
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + melt
                   lhabs(c) = lhabs(c) + melt*hfus
                   ! Update heat capacity
                   cv(c,j) = cv(c,j) + melt*(cpliq-cpice)
                   t_soisno(c,j) = tfrz + heatrem/cv(c,j)
                   ! Prevent tiny residuals
                   if (h2osoi_ice(c,j) < smallnumber) h2osoi_ice(c,j) = 0._r8
                   if (h2osoi_liq(c,j) < smallnumber) h2osoi_liq(c,j) = 0._r8
                end if

             end if
          end do
       end do

       ! Update eflx_snomelt(c)
       do fc = 1,num_lakec
          c = filter_lakec(fc)
          eflx_snomelt(c) = qflx_snomelt(c)*hfus
       end do
       do j = -nlevsno+1,0
          do fc = 1,num_lakec
             c = filter_lakec(fc)
             qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
          end do
       end do

     end associate

   end subroutine PhaseChange_Lake

 end module LakeTemperatureMod

