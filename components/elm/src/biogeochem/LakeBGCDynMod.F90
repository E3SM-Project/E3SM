module LakeBGCDynMod

   ! ==============================================================================================
   ! This module creates subroutines to represent methane-centered
   ! biogeochemical processes in the lake system (Tan et al., 2015, 2017).
   ! ==============================================================================================
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use decompMod         , only : bounds_type
   use abortutils        , only : endrun
   use LakeStateType     , only : lakestate_type
   use SoilStateType     , only : soilstate_type
   use GridcellType      , only : grc_pp
   use TopounitDataType  , only : top_as, top_af
   use LakeBGCType       , only : lakebgc_type
   use lnd2atmType       , only : lnd2atm_type
   use ColumnType        , only : col_pp
   use ColumnDataType    , only : col_es, col_ws

   implicit none
   save
   private

   ! !PUBLIC MEMBER FUNCTIONS:
   public  :: readLakeBGCParams
   public  :: LakeBGCDynamics

   ! !PRIVATE MEMBER FUNCTIONS:
   private :: Photosynthesis
   private :: AutotrophicR
   private :: HeterotrophicR
   private :: Methanogenesis 
   private :: Methanotrophy
   private :: BubbleDynamics
   private :: SedimentEbullition 

   type, private :: LakeBGCParamsType
      ! Chla-specific light saturated growth rate (gC (g Chl)-1 d-1)
      real(r8) :: Vchs
      real(r8) :: Vchl
      ! half-saturation for phosphorus limitation (mol/m3)
      real(r8) :: Ksrps
      real(r8) :: Ksrpl
      ! metabolic loss rate (d-1)
      real(r8) :: Klrs
      real(r8) :: Klrl
      ! phytoplankton resuspension fraction
      real(r8) :: frcResusps
      real(r8) :: frcResuspl
      ! ebullition rate (s-1) 
      real(r8) :: Re 
      ! relative concentration at which ebullition begins
      real(r8) :: Ae
      ! sediment OC dampening rate (m-1)
      real(r8) :: csedDMP
      ! ice bubble flux and dissolution rate (s-1)
      real(r8) :: icebflux
      real(r8) :: icebloss
      ! sediment OC decomposition rate through aerobic respiration (s-1) 
      real(r8) :: Rcapas
      real(r8) :: Rcaact
      ! sediment OC decomposition rate through methanogenesis (s-1)
      real(r8) :: Rcpas
      real(r8) :: Rcact
      real(r8) :: Rcold
      ! sediment OC decomposition Q10 through methanogenesis
      real(r8) :: PQ10pas
      real(r8) :: PQ10act
      ! O2 suppression coefficient for mathanogenesis (m3 water mol^-1)
      real(r8) :: etaO2
      ! CH4 oxidation potential (mol/m3/s)
      real(r8) :: Qch4
      ! CH4 oxidation Q10
      real(r8) :: OQ10
      ! CH4 concentration exponent coefficient of CH4 oxidation (unitless)
      real(r8) :: betaCH4
      ! O2 inhibition coefficient of CH4 oxidation (unitless)
      real(r8) :: lamdaO2
      ! oxic CH4 production rate to photosynthesis (mol CH4/m3/s (mol O2/m3/s)-1)
      real(r8) :: Rcoxic

   end type
   type(LakeBGCParamsType), private ::  LakeBGCParamsInst

contains

   !----------------------------------------------------------------------- 
   subroutine LakeBGCDynamics(bounds, num_lakec, filter_lakec, num_lakep, &
      filter_lakep, soilstate_vars, lakestate_vars, lakebgc_vars, &
      lnd2atm_vars)
      !
      ! !DESCRIPTION:
      ! Calculate lake biogeochemical dynamics in both water and sediment columns 
      !
      ! !USES:
      use subgridAveMod      , only : c2g
      use clm_time_manager   , only : get_step_size
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak, nphytolak, nsoilclak
      use elm_varcon         , only : vkc, grav, denh2o, tfrz, cnfac
      use elm_varcon         , only : catomw, patomw
      use elm_varctl         , only : iulog
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak, wsrplak
      use LakeBGCType        , only : small_phyto, large_phyto
      use LakeBGCType        , only : actC, pasC, nsubstep
      use LakeBGCType        , only : Kbg, Kc
      use TridiagonalMod     , only : Tridiagonal
      !
      ! !ARGUMENTS:
      type(bounds_type)      , intent(in)    :: bounds
      integer                , intent(in)    :: num_lakec       ! number of column lake points in column filter
      integer                , intent(in)    :: filter_lakec(:) ! column filter for non-lake points
      integer                , intent(in)    :: num_lakep       ! number of column non-lake points in pft filter
      integer                , intent(in)    :: filter_lakep(:) ! patch filter for non-lake points
      type(soilstate_type)   , intent(in)    :: soilstate_vars
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      type(lnd2atm_type)     , intent(inout) :: lnd2atm_vars
      !
      ! !LOCAL VARIABLES:
      integer  :: c, t, g, fc                         ! indices
      integer  :: j, k, r, j0, js                     ! indices
      integer  :: iter                                ! time indices
      logical  :: isFullTStep                         ! flag for full time step
      real(r8) :: dtime                               ! timestep size [seconds]
      real(r8) :: temp, por, depth                    ! auxiliary variables
      real(r8) :: solu_avg, biomas_avg                ! auxiliary variables 
      real(r8) :: zsum, dzm, dzp, dzb                 ! auxiliary variables
      real(r8) :: deficit, deficitCum, biomas_tot
      real(r8) :: blamda, bsusp, frcResusp
      real(r8) :: ztot_act, cdep, cdep_tmp
      integer  :: jtop(bounds%begc:bounds%endc)                      ! top level for each column
      integer  :: jwat(bounds%begc:bounds%endc)                      ! top water layer index
      integer  :: jmix(bounds%begc:bounds%endc)                      ! mixing layer index
      real(r8) :: Hcline(bounds%begc:bounds%endc)                    ! thermocline depth (m)
      real(r8) :: zTopblayer(bounds%begc:bounds%endc)                ! lower bound of top boundary layer (m)
      real(r8) :: zBotblayer(bounds%begc:bounds%endc)                ! upper bound of bottom boundary layer (m)
      real(r8) :: eb_sed(bounds%begc:bounds%endc,1:ngaslak)          ! gas ebullition from sediment (mol/m2/s)
      real(r8) :: eb_surf(bounds%begc:bounds%endc,1:ngaslak)         ! gas ebullition from surface (mol/m2/s)
      real(r8) :: df_surf(bounds%begc:bounds%endc,1:nsolulak)        ! diffusion from surface (mol/m2/s)
      real(r8) :: df_sed(bounds%begc:bounds%endc,1:nsolulak)         ! diffusion from sediment (mol/m2/s)
      real(r8) :: errch4(bounds%begc:bounds%endc)                    ! error (Mol CH4 /m^2) [+ = too much CH4]
      real(r8) :: errbiomas(bounds%begc:bounds%endc)                 ! error (gC/m2) [+ = too much biomass]
      real(r8) :: delta_ch4(bounds%begc:bounds%endc)                 ! net change (Mol CH4 /m^2) 
      real(r8) :: delta_iceb(bounds%begc:bounds%endc)                ! net ice bubble change (Mol CH4 /m^2)
      real(r8) :: soilc_act_tot(bounds%begc:bounds%endc)             ! total active OC pool (gC/m2)
      real(r8) :: bphyto_dep_tot(bounds%begc:bounds%endc)            ! total phytoplankton deposition (gC/m2/s)
      real(r8) :: bphyto_dead(bounds%begc:bounds%endc,1:nphytolak)   ! vertically-integrated dead phytoplankton (gC/m2/s)
      real(r8) :: conc_eq(bounds%begc:bounds%endc,1:nsolulak)        ! equilibrium concentration (mol/m^3)
      real(r8) :: conc_iceb_old(bounds%begc:bounds%endc,1:ngaslak)   ! ice-trapped bubble gases at the last time step (mol/m2)
      real(r8) :: ex_iceb(bounds%begc:bounds%endc,1:ngaslak)         ! ice-trapped bubble oxidation rate (mol/m2/s)
      real(r8) :: kg(bounds%begc:bounds%endc,1:nsolulak)             ! transfer velocity (m/s)
      real(r8) :: zx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! interface depth (+ below surface) for whole column (m)
      real(r8) :: dzx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)     ! cell thickness (+ below surface) for whole column (m)
      real(r8) :: qx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! porosity (+ below surface) for whole column
      real(r8) :: dx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! diffusivity (+ below surface) for whole column (m^2/s)
      real(r8) :: va(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "a" vector for tridiagonal matrix
      real(r8) :: vb(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "b" vector for tridiagonal matrix
      real(r8) :: vc(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "c" vector for tridiagonal matrix
      real(r8) :: vr(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "r" vector for tridiagonal solution
      real(r8) :: ex(bounds%begc:bounds%endc,1:nlevlak,1:ngaslak)                ! bubble dissolution rate (mol/m3/s)
      real(r8) :: sx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:nsolulak)       ! source (+ below surface) for whole column (mol/m^3/s)
      real(r8) :: csrc(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:nsolulak)     ! solute production rate (mol/m3/s)
      real(r8) :: csnk(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:nsolulak)     ! solute loss rate (mol/m3/s)
      real(r8) :: gebul(bounds%begc:bounds%endc,1:nlevsoi,1:ngaslak)             ! sediment ebullition rate (mol/m3/s) 
      real(r8) :: conc_old(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:nsolulak) ! solute concentration at the last time step (mol/m3)
      real(r8) :: conc_new(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:nsolulak) ! updated solute concentration (mol/m3)
      real(r8) :: vpx(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)             ! phytoplankton settling velocity (m/s) [+ = sinking]
      real(r8) :: biomas_old(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)      ! phytoplankton biomass at the last time step (gC/m3)
      real(r8) :: biomas_new(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)      ! updated phytoplankton biomass (gC/m3)
      real(r8) :: hr_vr(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)               ! vertically-resolved heterotrophic respiration (gC/m3/s)
      real(r8) :: gpp_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)          ! vertically-resolved primary production (gC/m3/s)
      real(r8) :: bloss_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)        ! vertically-resolved metabolic loss (gC/m3/s)
      real(r8) :: ar_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)           ! vertically-resolved autotrophic respiration (gC/m3/s)
      real(r8) :: pch4_vr(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)             ! vertically-resolved methane production (mol/m3/s)
      real(r8) :: och4_vr(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)             ! vertically-resolved methane oxidation (mol/m3/s)
      real(r8) :: bmotion_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)      ! vertically-resolved time-averaged biomass motion (gC/m2/s)
      real(r8) :: soilc_loss(bounds%begc:bounds%endc,1:nlevsoi,1:nsoilclak)      ! vertically-resolved soil C loss (gC/m3/s)
      !----------------------------------------------------------------------- 
      
      associate(                                                     &
            dz_lake              => col_pp%dz_lake                      , & ! Input: [real(r8) (:,:)] layer thickness for lake [m]       
            z_lake               => col_pp%z_lake                       , & ! Input: [real(r8) (:,:)] layer depth for lake [m]
            dz                   => col_pp%dz                           , & ! Input: [real(r8) (:,:)] layer thickness [m]              
            z                    => col_pp%z                            , & ! Input: [real(r8) (:,:)] layer depth [m]                       
            lakedepth            => col_pp%lakedepth                    , & ! Input: [real(r8) (:)] column lake depth [m]

            t_lake               => col_es%t_lake                       , & ! Input: [real(r8) (:,:)] col lake temperature [Kelvin] 
            t_soisno             => col_es%t_soisno                     , & ! Input: [real(r8) (:,:)] soil (or snow) temperature [Kelvin]

            forc_wind            => top_as%windbot                      , & ! Input: [real(r8) (:)] horizontal component of wind at atmospheric forcing height [m/s]

            watsat               => soilstate_vars%watsat_col           , & ! Input: [real(r8) (:,:)] volumetric soil water at saturation (porosity)

            icethick             => lakestate_vars%lake_icethick_col    , & ! Input: [real(r8) (:)] ice thickness [m] (integrated if lakepuddling)
            icefrac              => lakestate_vars%lake_icefrac_col     , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            fsds_vis             => lakestate_vars%lake_fsds_vis_col    , & ! Input: [real(r8) (:,:)] incident vis radiation (W/m^2)

            lake_csed            => lakebgc_vars%csed_col               , & ! Input: [real(r8) (:)] lake sediment OC density [kgC/m3]
            lake_cdep            => lakebgc_vars%cdep_col               , & ! Input: [real(r8) (:)] allocthonous OC deposition [gC/m2/yr]
            lake_soilc_old       => lakebgc_vars%soilc_old_col          , & ! Input: [real(r8) (:,:,:)] lake sediment Pleistocene carbon [gC/m3]

            nem_lake_grc         => lnd2atm_vars%nem_lake_grc           , & ! Output: [real(r8) (:)] gridcell average net methane correction to CO2 flux [gC/m^2/s]

            soilpor              => lakebgc_vars%soilpor_col            , & ! Output: [real(r8) (:,:)] lake sediment porosity
            fsds_vis_bgc         => lakebgc_vars%fsds_vis_col           , & ! Output: [real(r8) (:,:)] incident vis radiation for BGC (W/m^2)
            conc_wat             => lakebgc_vars%conc_wat_col           , & ! Output: [real(r8) (:,:,:)] water solute conc [mol/m3]
            conc_sed             => lakebgc_vars%conc_sed_col           , & ! Output: [real(r8) (:,:,:)] porewater solute conc [mol/m3]
            conc_bubl            => lakebgc_vars%conc_bubl_col          , & ! Output: [real(r8) (:,:,:)] lake bubble gas [mol/m3]
            conc_iceb            => lakebgc_vars%conc_iceb_col          , & ! Output: [real(r8) (:,:)] lake bubble trapped in ice [mol/m2]
            biomas_phyto         => lakebgc_vars%biomas_phyto_col       , & ! Output: [real(r8) (:,:,:)] lake phytoplankton biomass [gC/m3]
            ctot_dep             => lakebgc_vars%ctot_dep_col           , & ! Output: [real(r8) (:)] lake OC deposition rate [gC/m2/s]
            ch4_sed_diff         => lakebgc_vars%ch4_sed_diff_col       , & ! Output: [real(r8) (:)] lake ch4 diffusion at water-sediment interface [mol/m2/s]
            ch4_sed_ebul         => lakebgc_vars%ch4_sed_ebul_col       , & ! Output: [real(r8) (:)] lake ch4 ebullition at water-sediment interface [mol/m2/s]
            ch4_surf_diff        => lakebgc_vars%ch4_surf_diff_col      , & ! Output: [real(r8) (:)] lake ch4 diffusion at air-water interface [mol/m2/s]
            ch4_surf_ebul        => lakebgc_vars%ch4_surf_ebul_col      , & ! Output: [real(r8) (:)] lake ch4 ebullition at air-water interface [mol/m2/s]
            ch4_surf_flux        => lakebgc_vars%ch4_surf_flux_col      , & ! Output: [real(r8) (:)] lake ch4 total flux at air-water interface [kgC/m2/s]
            ch4_prod_wat         => lakebgc_vars%ch4_prod_wat_col       , & ! Output: [real(r8) (:,:)] lake ch4 production rate [mol/m3/s]
            ch4_prod_sed         => lakebgc_vars%ch4_prod_sed_col       , & ! Output: [real(r8) (:,:)] sediment ch4 production rate [mol/m3/s]
            ch4_oxid_wat         => lakebgc_vars%ch4_oxid_wat_col       , & ! Output: [real(r8) (:,:)] lake ch4 oxidation rate [mol/m3/s]
            ch4_oxid_sed         => lakebgc_vars%ch4_oxid_sed_col       , & ! Output: [real(r8) (:,:)] sediment ch4 oxidation rate [mol/m3/s]
            ch4_prod_tot         => lakebgc_vars%ch4_prod_tot_col       , & ! Output: [real(r8) (:)] lake total CH4 production rate [mol/m2/s]
            ch4_oxid_tot         => lakebgc_vars%ch4_oxid_tot_col       , & ! Output: [real(r8) (:)] lake total CH4 oxidation rate [mol/m2/s]
            hr_wat               => lakebgc_vars%hr_wat_vr_col          , & ! Output: [real(r8) (:,:)] lake heterotrophic respiration [gC/m3/s]
            hr_sed               => lakebgc_vars%hr_sed_vr_col          , & ! Output: [real(r8) (:,:)] sediment heterotrophic respiration [gC/m3/s]
            nem_col              => lakebgc_vars%nem_col                , & ! Output: [real(r8) (:)] net adjustment to atm. C flux from methane production [g C/m**2/s]
            lake_gpp             => lakebgc_vars%gpp_vr_col             , & ! Output: [real(r8) (:,:)] depth-resolved gross primary production [gC/m3/s]
            lake_npp             => lakebgc_vars%npp_vr_col             , & ! Output: [real(r8) (:,:)] depth-resolved net primary production [gC/m3/s]
            lake_gpp_tot         => lakebgc_vars%gpp_tot_col            , & ! Output: [real(r8) (:)] depth-integrated gross primary production [gC/m2/s]
            lake_npp_tot         => lakebgc_vars%npp_tot_col            , & ! Output: [real(r8) (:)] depth-integrated net primary production [gC/m2/s]
            lake_chla            => lakebgc_vars%chla_col               , & ! Output: [real(r8) (:,:)] chlorophyll-a conc [g/m3]
            lake_soilc           => lakebgc_vars%soilc_col              , & ! Output: [real(r8) (:,:,:)] lake sediment carbon pools [gC/m3]
            totsoilc             => lakebgc_vars%totsoilc_col           , & ! Output: [real(r8) (:)] sediment total OC [gC/m2]
            totphytoc            => lakebgc_vars%totphytoc_col          & ! Output: [real(r8) (:)] lake total phytoplankton biomass [gC/m2]
            )

      ! Get time step
      dtime = get_step_size() / DBLE(nsubstep)

      ! physical and BGC processes
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         t = col_pp%topounit(c)

         ! initialize flux/rate variables
         ctot_dep(c)          = 0._r8
         ch4_sed_diff(c)      = 0._r8
         ch4_sed_ebul(c)      = 0._r8
         ch4_surf_diff(c)     = 0._r8
         ch4_surf_ebul(c)     = 0._r8 
         ch4_surf_flux(c)     = 0._r8
         ch4_prod_wat(c,:)    = 0._r8
         ch4_prod_sed(c,:)    = 0._r8
         ch4_oxid_wat(c,:)    = 0._r8
         ch4_oxid_sed(c,:)    = 0._r8
         ch4_prod_tot(c)      = 0._r8
         ch4_oxid_tot(c)      = 0._r8
         hr_wat(c,:)          = 0._r8
         hr_sed(c,:)          = 0._r8
         lake_gpp(c,:)        = 0._r8
         lake_npp(c,:)        = 0._r8
         lake_gpp_tot(c)      = 0._r8
         lake_npp_tot(c)      = 0._r8

         ! calculate incident radiation for BGC
         do j = 1, nlevlak
            if (j==1 .or. icefrac(c,j)>0._r8) then
               fsds_vis_bgc(c,j) = fsds_vis(c,j)
            else
               biomas_tot = sum(biomas_phyto(c,j,:))
               fsds_vis_bgc(c,j) = fsds_vis_bgc(c,j-1) * &
                     exp(-(Kbg+Kc*biomas_tot)*dz_lake(c,j))
            end if
         end do

         ! define specified water layers
         call DefineWaterLayers(lakestate_vars, c, jtop(c), jwat(c), &
               jmix(c), zTopblayer(c), zBotblayer(c), Hcline(c))

         ! set lake sediment porosity
         do j = 1, nlevgrnd
            ! update sediment porosity
            if (z(c,j)<0.15_r8) then
               soilpor(c,j) = max(watsat(c,j),0.8_r8)
            else
               soilpor(c,j) = watsat(c,j)
            end if
         end do

         ! update Pleistocene-age C
         call SetYedomaLakeOldC(lakebgc_vars, c)

         ! diffusivity
         call CalcSoluteDiffusivity(lakestate_vars, c, dx(c,:))

         ! equilibrium concentration
         call CalcGasEQConc(lakebgc_vars, c, conc_eq(c,:))

         ! flux velocity
         temp = t_lake(c,1)
         do k = 1, nsolulak
            if (icethick(c)<1.e-8_r8) then
               if (k<=ngaslak) then
                  kg(c,k) = GasFluxVelocity(forc_wind(t), temp, k)
               else
                  kg(c,k) = 1.1574e-5_r8  ! 1/24/3600 
               end if
            else
               kg(c,k) = 0._r8
            end if
         end do
      end do

      ! sub-step cycle for solute, bubble, and biomass dynamics
      do iter = 1, nsubstep, 1

         do fc = 1, num_lakec
            c = filter_lakec(fc)
            t = col_pp%topounit(c)

            ! record solute conc at the last time step
            do j = 1, nlevlak+nlevsoi
               if (j<=nlevlak) then
                  conc_old(c,j,:) = conc_wat(c,j,:)
                  ! O2 supply from surface & groundwater
                  if (j==jwat(c) .and. icethick(c)>1.e-8_r8 .and. &
                        lakedepth(c)>=2._r8) then
                     conc_old(c,j,wo2lak) = 0.35_r8
                  end if
               else
                  conc_old(c,j,:) = conc_sed(c,j-nlevlak,:)
               end if
            end do
            conc_iceb_old(c,:) = conc_iceb(c,:) 

            ! record phytoplankton biomass at the last time step
            do j = 1, nlevlak
               biomas_old(c,j,:) = biomas_phyto(c,j,:)
            end do

            ! record soil C at the last time step 
            soilc_act_tot(c) = 0._r8
            do j = 1, nlevsoi
               soilc_act_tot(c) = soilc_act_tot(c) + lake_soilc(c,j,actC)*dz(c,j)
            end do 

            ! surface and sediment diffusive fluxes 
            dzb = 0.5_r8 * (dz_lake(c,nlevlak) + dz(c,1))
            do k = 1, nsolulak
               df_surf(c,k) = kg(c,k) * (conc_old(c,1,k) - conc_eq(c,k))
               df_sed(c,k) = dx(c,nlevlak) / dzb * (conc_old(c,nlevlak+1,k)/ &
                     soilpor(c,1) - conc_old(c,nlevlak,k))
            end do
            ch4_sed_diff(c) = ch4_sed_diff(c) + df_sed(c,wch4lak)*dtime

            ! initialize source and sink terms of solutes 
            csrc(c,:,:)       = 0._r8
            csnk(c,:,:)       = 0._r8
            gebul(c,:,:)      = 0._r8
            eb_surf(c,:)      = 0._r8

            ! initialize source and sink terms of c pools
            soilc_loss(c,:,:) = 0._r8
            bphyto_dead(c,:)  = 0._r8
            bphyto_dep_tot(c) = 0._r8
            gpp_vr(c,:,:)     = 0._r8
            bloss_vr(c,:,:)   = 0._r8
            ar_vr(c,:,:)      = 0._r8
            hr_vr(c,:)        = 0._r8
            pch4_vr(c,:)      = 0._r8
            och4_vr(c,:)      = 0._r8
            bmotion_vr(c,:,:) = 0._r8

            ! phytoplankton settling velocity
            call PhytoSettleVelocity(lakestate_vars, c, jwat(c), zTopblayer(c), &
                     zBotblayer(c), conc_old(c,:,wsrplak), vpx(c,:,:))

            ! sediment bubbling
            call SedimentEbullition(soilstate_vars, lakebgc_vars, c, &
                     conc_old(c,:,1:ngaslak), gebul(c,:,:))
            do k = 1, ngaslak         
               eb_sed(c,k) = sum(gebul(c,:,k)*dz(c,1:nlevsoi))
            end do
            ! CH4 ebullition at the water-sediment interface
            ch4_sed_ebul(c) = ch4_sed_ebul(c) + eb_sed(c,wch4lak)*dtime

            ! bubble transport in water (11 size bins)
            if (iter==1) then
               isFullTStep = .True.
            else
               isFullTStep = .False.
            end if
            call BubbleDynamics(lakestate_vars, lakebgc_vars, c, jwat(c), isFullTStep, &
                     dtime, conc_old(c,:,1:ngaslak), eb_sed(c,:), eb_surf(c,:), &
                     ex(c,:,:), ex_iceb(c,:), conc_bubl(c,:,:), conc_iceb(c,:))
            ch4_surf_ebul(c) = ch4_surf_ebul(c) + eb_surf(c,wch4lak)*dtime

            ! photosynthesis
            call Photosynthesis(lakestate_vars, lakebgc_vars, c, conc_old(c,:,:), &
                     biomas_old(c,:,:), csrc(c,:,:), csnk(c,:,:), gpp_vr(c,:,:), &
                     lake_chla(c,:))

            ! autotrophic respiration
            call AutotrophicR(lakestate_vars, c, conc_old(c,:,:), biomas_old(c,:,:), &
                     csrc(c,:,:), csnk(c,:,:), bloss_vr(c,:,:), ar_vr(c,:,:))
            do k = 1, nphytolak
               bphyto_dead(c,k) = bphyto_dead(c,k) + sum(dz_lake(c,1:nlevlak)* &
                     (bloss_vr(c,1:nlevlak,k)-ar_vr(c,1:nlevlak,k)))
            end do

            ! heterotrophic respiration
            call HeterotrophicR(lakestate_vars, c, conc_old(c,:,:), lake_soilc(c,:,:), &
                     csrc(c,:,:), csnk(c,:,:), soilc_loss(c,:,:), hr_vr(c,:))
            hr_wat(c,1:nlevlak) = hr_wat(c,1:nlevlak) + hr_vr(c,1:nlevlak)*dtime
            hr_sed(c,1:nlevsoi) = hr_sed(c,1:nlevsoi) + hr_vr(c,nlevlak+1:nlevlak+nlevsoi)*dtime

            ! methanogenesis 
            call Methanogenesis(lakestate_vars, lakebgc_vars, c, conc_old(c,:,:), &
                     lake_soilc(c,:,:), lake_soilc_old(c,:), csrc(c,:,:), &
                     csnk(c,:,:), soilc_loss(c,:,:), pch4_vr(c,:)) 
            ch4_prod_wat(c,1:nlevlak) = ch4_prod_wat(c,1:nlevlak) + &
                  pch4_vr(c,1:nlevlak)*dtime
            ch4_prod_sed(c,1:nlevsoi) = ch4_prod_sed(c,1:nlevsoi) + &
                  pch4_vr(c,nlevlak+1:nlevlak+nlevsoi)*dtime 

            ! methanotrophy
            call Methanotrophy(lakestate_vars, c, conc_old(c,:,:), csrc(c,:,:), &
                     csnk(c,:,:), och4_vr(c,:))
            ch4_oxid_wat(c,1:nlevlak) = ch4_oxid_wat(c,1:nlevlak) + &
                  och4_vr(c,1:nlevlak)*dtime
            ch4_oxid_sed(c,1:nlevsoi) = ch4_oxid_sed(c,1:nlevsoi) + &
                  och4_vr(c,nlevlak+1:nlevlak+nlevsoi)*dtime
         end do

         ! Set up interface depths, zx, gas diffusivity, dx, and porosity, qx.
         do j = 1, nlevlak+nlevsoi
            do fc = 1, num_lakec
               c = filter_lakec(fc)
            
               if (j<=nlevlak) then
                  zx(c,j) = z_lake(c,j)
                  dzx(c,j) = dz_lake(c,j)
               else
                  zx(c,j) = zx(c,nlevlak) + 0.5_r8*dz_lake(c,nlevlak) + z(c,j-nlevlak)
                  dzx(c,j) = dz(c,j-nlevlak)
               end if
               if (j<=nlevlak) then
                  qx(c,j) = 1.0_r8
               else
                  qx(c,j) = soilpor(c,j-nlevlak)
               end if

               do k = 1, nsolulak
                  if (j==1) then
                     if (k<=ngaslak) then
                        sx(c,j,k) = ex(c,j,k) + csrc(c,j,k) - csnk(c,j,k) - &
                              df_surf(c,k)/dz_lake(c,j)
                     else
                        sx(c,j,k) = csrc(c,j,k) - csnk(c,j,k) - &
                              df_surf(c,k)/dz_lake(c,j)
                     end if
                  else if (j<=nlevlak) then
                     if (k<=ngaslak) then
                        sx(c,j,k) = ex(c,j,k) + csrc(c,j,k) - csnk(c,j,k)
                     else
                        sx(c,j,k) = csrc(c,j,k) - csnk(c,j,k)
                     end if
                  else
                     if (k<=ngaslak) then
                        sx(c,j,k) = csrc(c,j,k) - csnk(c,j,k) - gebul(c,j-nlevlak,k)
                     else
                        sx(c,j,k) = csrc(c,j,k) - csnk(c,j,k)
                     end if
                  end if
                  if (j==jwat(c) .and. k<=ngaslak) then
                     sx(c,j,k) = sx(c,j,k) + ex_iceb(c,k)/dz_lake(c,j)
                  end if
               end do
            end do
         end do

         ! solve governing equations 
         do k = 1, nsolulak
         
            ! Set up vector r and vectors a, b, c that define tridiagonal matrix 
            do j = 1, nlevlak+nlevsoi
               do fc = 1, num_lakec
                  c = filter_lakec(fc)
              
                  if (j==1) then ! top layer
                     dzp      = zx(c,j+1) - zx(c,j)
                     va(c,j)  = 0._r8          
                     vb(c,j)  = 1._r8 + (1._r8-cnfac)/qx(c,j)*dtime*dx(c,j)/dzp/dzx(c,j)
                     vc(c,j)  = -(1._r8-cnfac)/qx(c,j+1)*dtime*dx(c,j)/dzp/dzx(c,j)
                     vr(c,j)  = conc_old(c,j,k) + dtime*sx(c,j,k) + dtime*cnfac*dx(c,j)* &
                           (conc_old(c,j+1,k)/qx(c,j+1)-conc_old(c,j,k)/qx(c,j))/dzp/dzx(c,j)
                  else if (j<nlevlak+nlevsoi) then
                     dzp      = zx(c,j+1) - zx(c,j)
                     dzm      = zx(c,j) - zx(c,j-1)
                     va(c,j)  = -(1._r8-cnfac)/qx(c,j-1)*dtime*dx(c,j-1)/dzm/dzx(c,j)
                     vb(c,j)  = 1._r8 + (1._r8-cnfac)/qx(c,j)*dtime*(dx(c,j)/dzp+ &
                           dx(c,j-1)/dzm)/dzx(c,j)
                     vc(c,j)  = -(1._r8-cnfac)/qx(c,j+1)*dtime*dx(c,j)/dzp/dzx(c,j)
                     vr(c,j)  = conc_old(c,j,k) + dtime*sx(c,j,k) + dtime*cnfac/dzx(c,j)* &
                           (dx(c,j)*(conc_old(c,j+1,k)/qx(c,j+1)-conc_old(c,j,k)/qx(c,j))/dzp - &
                           dx(c,j-1)*(conc_old(c,j,k)/qx(c,j)-conc_old(c,j-1,k)/qx(c,j-1))/dzm)
                  else  ! bottom layer
                     dzm      = zx(c,j) - zx(c,j-1)
                     va(c,j)  = -(1._r8-cnfac)/qx(c,j-1)*dtime*dx(c,j-1)/dzm/dzx(c,j)
                     vb(c,j)  = 1._r8 + (1._r8-cnfac)/qx(c,j)*dtime*dx(c,j-1)/dzm/dzx(c,j)
                     vc(c,j)  = 0._r8
                     vr(c,j)  = conc_old(c,j,k) + dtime*sx(c,j,k) - dtime*cnfac*dx(c,j-1)* &
                           (conc_old(c,j,k)/qx(c,j)-conc_old(c,j-1,k)/qx(c,j-1))/dzm/dzx(c,j)
                  end if
               end do
            end do

            call Tridiagonal(bounds, 1, nlevlak + nlevsoi, &
                     jtop(bounds%begc:bounds%endc), &
                     num_lakec, filter_lakec, &
                     va(bounds%begc:bounds%endc, :), &
                     vb(bounds%begc:bounds%endc, :), &
                     vc(bounds%begc:bounds%endc, :), &
                     vr(bounds%begc:bounds%endc, :), &
                     conc_new(bounds%begc:bounds%endc, :, k))
         end do

         ! update negative concentrations 
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            do k = 1, nsolulak
               ! mixing
               j0 = jwat(c)
               solu_avg = sum(conc_new(c,j0:jmix(c),k)*dz_lake(c,j0:jmix(c))) / &
                     sum(dz_lake(c,j0:jmix(c)))
               conc_new(c,j0:jmix(c),k) = solu_avg

               deficitCum = 0._r8
               ! correct negative concentration
               do j = 1, nlevlak+nlevsoi
                  if (conc_new(c,j,k)<0._r8) then
                     deficit = -conc_new(c,j,k) * dzx(c,j)  ! mol/m2
                     deficitCum = deficitCum + deficit
                     if (k==wch4lak) then ! for CH4
                        if (deficit > 1.e-3_r8) then
                           if (deficit > 1.e-2_r8) then
                              write(iulog,*)'Note: sink > source in LakeBGCDynamics, sources are changing '// &
                                    ' quickly relative to diffusion timestep, and/or diffusion is rapid.'
                              g = col_pp%gridcell(c)
                              write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
                              write(iulog,*)'This typically occurs when there is a larger than normal '// &
                                    ' diffusive flux.'
                              write(iulog,*)'If this occurs frequently, consider reducing lake bgc model '// &
                                    'timestep, or reducing the max. sink per timestep.'
                           end if
                           !write(iulog,*) 'Negative conc. in LakeBGCDynamics g,c,j,deficit (mol):',g,c,j,deficit
                        end if
                     end if
                     conc_new(c,j,k) = 0._r8
                  end if
               end do
               ! add negative CH4 to bottom layer and negative conc of other
               ! gases to surface layer
               if (k==wch4lak .and. deficitCum>0._r8) then
                  ! add negative CH4 to bottom layer
                  conc_new(c,nlevlak,k) = conc_new(c,nlevlak,k) - deficitCum/dzx(c,nlevlak)
                  if (conc_new(c,nlevlak,k)<0._r8) then
                     deficit = -conc_new(c,nlevlak,k) * dzx(c,nlevlak)
                     if (icethick(c)<1.e-8_r8) then
                        df_surf(c,k) = df_surf(c,k) - deficit/dtime
                     else
                        conc_iceb(c,k) = conc_iceb(c,k) - deficit
                     end if
                     conc_new(c,nlevlak,k) = 0._r8
                  end if
               else if (k<=ngaslak .and. deficitCum>0._r8) then
                  conc_new(c,j0,k) = conc_new(c,j0,k) - deficitCum/dzx(c,j0)
                  if (conc_new(c,j0,k)<0._r8) then
                     deficit = -conc_new(c,j0,k) * dzx(c,j0)
                     if (icethick(c)<1.e-8_r8) then
                        df_surf(c,k) = df_surf(c,k) - deficit/dtime
                     else
                        conc_iceb(c,k) = conc_iceb(c,k) - deficit
                     end if
                     conc_new(c,j0,k) = 0._r8
                  end if
               end if

               ! other solute corrections in water column
               do j = 1, nlevlak
                  ! correct solutes in ice layers
                  if (conc_new(c,j,k)>0._r8 .and. icefrac(c,j)>=1._r8) then
                     if (k<=ngaslak) then
                        conc_iceb(c,k) = conc_iceb(c,k) + conc_new(c,j,k)*dzx(c,j)
                     end if
                     conc_new(c,j,k) = 0._r8
                  end if
               end do
               ! assume dissolved N2 always replete
               conc_new(c,1:nlevlak,wn2lak) = conc_eq(c,wn2lak)

               ! set dissolved gas conc for outputs
               conc_wat(c,1:nlevlak,k) = conc_new(c,1:nlevlak,k)
               conc_sed(c,1:nlevsoi,k) = conc_new(c,nlevlak+1:nlevlak+nlevsoi,k)
            end do
         end do

         ! Do Balance Check
         do j = 1, nlevlak+nlevsoi
            do fc = 1, num_lakec
               c = filter_lakec(fc)
            
               if (j==1) then
                  errch4(c) = conc_iceb(c,wch4lak) - conc_iceb_old(c,wch4lak)
                  delta_iceb(c) = errch4(c) 
                  delta_ch4(c) = 0._r8
               end if
               if (j<=nlevlak) then
                  errch4(c) = errch4(c) + dzx(c,j)*(conc_new(c,j,wch4lak)- &
                        conc_old(c,j,wch4lak))
                  errch4(c) = errch4(c) - pch4_vr(c,j)*dzx(c,j)*dtime
                  errch4(c) = errch4(c) + och4_vr(c,j)*dzx(c,j)*dtime
                  delta_ch4(c) = delta_ch4(c) + dzx(c,j)*(conc_new(c,j,wch4lak)- &
                        conc_old(c,j,wch4lak))
               else
                  errch4(c) = errch4(c) + dzx(c,j)*(conc_new(c,j,wch4lak)- &
                        conc_old(c,j,wch4lak))
                  errch4(c) = errch4(c) - pch4_vr(c,j)*dzx(c,j)*dtime
                  errch4(c) = errch4(c) + och4_vr(c,j)*dzx(c,j)*dtime
                  delta_ch4(c) = delta_ch4(c) + dzx(c,j)*(conc_new(c,j,wch4lak)- &
                        conc_old(c,j,wch4lak))
               end if
            end do
         end do
      
         do fc = 1, num_lakec
            c = filter_lakec(fc)
      
            errch4(c) = errch4(c) + (df_surf(c,wch4lak) + eb_surf(c,wch4lak))*dtime
            if (abs(errch4(c)) < 1.e-8_r8) then
               df_surf(c,wch4lak) = df_surf(c,wch4lak) - errch4(c)/dtime
            else ! errch4 > 1e-8 mol / m^2 / timestep
               write(iulog,*)'CH4 Conservation Error in LakeBGCDynamics during c, errch4 (mol /m^2.timestep)', &
                  c,errch4(c)
               g = col_pp%gridcell(c)
               write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: CH4 Conservation Error in LakeBGCDynamics during diffusion'//&
                     errMsg(__FILE__, __LINE__))
            end if
            
            ch4_surf_diff(c) = ch4_surf_diff(c) + df_surf(c,wch4lak)*dtime
         end do

         ! solve phytoplankton equations
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            do k = 1, nphytolak
               do j = 1, nlevlak
                  if (biomas_old(c,j,k)>1.e-9_r8) then
                     blamda = (bloss_vr(c,j,k) - gpp_vr(c,j,k)) / biomas_old(c,j,k)
                     biomas_new(c,j,k) = biomas_old(c,j,k)*exp(-blamda*dtime)
                     ! time-averaged gain/loss
                     if (abs(blamda)>0._r8) then
                        gpp_vr(c,j,k) = gpp_vr(c,j,k)*(1._r8-exp(-blamda*dtime))/ &
                              blamda/dtime
                        bloss_vr(c,j,k) = bloss_vr(c,j,k)*(1._r8-exp(-blamda*dtime))/ &
                              blamda/dtime 
                        ar_vr(c,j,k) = ar_vr(c,j,k)*(1._r8-exp(-blamda*dtime))/ &
                              blamda/dtime
                     end if
                  else
                     blamda = bloss_vr(c,j,k) - gpp_vr(c,j,k)
                     biomas_new(c,j,k) = biomas_old(c,j,k) - blamda*dtime 
                  end if 
                  lake_gpp(c,j) = lake_gpp(c,j) + gpp_vr(c,j,k)*dtime
                  lake_npp(c,j) = lake_npp(c,j) + (gpp_vr(c,j,k)-ar_vr(c,j,k))*dtime

                  if (vpx(c,j,k)>=0._r8) then
                     bmotion_vr(c,j,k) = min(vpx(c,j,k),dzx(c,j)/dtime)*biomas_new(c,j,k)
                  else
                     bmotion_vr(c,j,k) = max(vpx(c,j,k),-dzx(c,j)/dtime)*biomas_new(c,j,k)
                  end if
               end do

               if (k==small_phyto) then
                  frcResusp = LakeBGCParamsInst%frcResusps
               else if (k==large_phyto) then
                  frcResusp = LakeBGCParamsInst%frcResuspl
               end if

               ! mixing
               j0 = jwat(c)
               biomas_avg = sum(biomas_new(c,j0:jmix(c),k)*dzx(c,j0:jmix(c))) / &
                     sum(dzx(c,j0:jmix(c)))
               biomas_new(c,j0:jmix(c),k) = biomas_avg

               ! vertical move including actively swim
               do j = 1, nlevlak
                  if (j==1) then
                     biomas_new(c,j,k) = biomas_new(c,j,k) - (abs(bmotion_vr(c,j,k)) - &
                           max(-bmotion_vr(c,j+1,k),0._r8)) / dzx(c,j) * dtime
                  else if (j==nlevlak) then
                     biomas_new(c,j,k) = biomas_new(c,j,k) - (abs(bmotion_vr(c,j,k)) - &
                           max(bmotion_vr(c,j-1,k),0._r8)) / dzx(c,j) * dtime
                  else
                     biomas_new(c,j,k) = biomas_new(c,j,k) - (abs(bmotion_vr(c,j,k)) - &
                           max(-bmotion_vr(c,j+1,k),0._r8) - max(bmotion_vr(c,j-1,k),0._r8)) / &
                           dzx(c,j) * dtime
                  end if
               end do

               ! deposition and resuspension
               if (lakedepth(c)>zBotblayer(c)) then
                  bphyto_dep_tot(c) = bphyto_dep_tot(c) + (1._r8 - frcResusp) * &
                        max(bmotion_vr(c,nlevlak,k),0._r8)
                  bsusp = max(bmotion_vr(c,nlevlak,k),0._r8) * frcResusp / &
                        (lakedepth(c) - zBotblayer(c))
                  do j = jwat(c), nlevlak
                     if (z_lake(c,j)>=zBotblayer(c)) then
                        biomas_new(c,j,k) = biomas_new(c,j,k) + bsusp*dtime
                     end if
                  end do
               else
                  zsum = sum(dzx(c,jwat(c):nlevlak)) 
                  bsusp = max(bmotion_vr(c,nlevlak,k),0._r8) / zsum
                  do j = jwat(c), nlevlak
                     biomas_new(c,j,k) = biomas_new(c,j,k) + bsusp*dtime
                  end do
               end if
            end do

         end do

         ! check negative values and balance
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            totphytoc(c) = 0._r8 
            do k = 1, nphytolak
               do j = 1, nlevlak
                  ! correct negative biomas
                  if (biomas_new(c,j,k)<0._r8) then
                     deficit = -biomas_new(c,j,k) * dzx(c,j)  ! gC/m2
                     if (deficit > 1.e-2_r8) then
                        if (deficit > 1.e-1_r8) then
                           write(iulog,*)'Note: sink > source in LakeBGCDynamics, sources are changing '// &
                                 ' quickly relative to deposition timestep, and/or deposition is rapid.'
                           g = col_pp%gridcell(c)
                           write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
                           write(iulog,*)'This typically occurs when there is a larger than normal '// &
                                 ' phytoplankton deposition.'
                           write(iulog,*)'If this occurs frequently, consider reducing lake bgc model '// &
                                 'timestep, or reducing the max. deposition per timestep.'
                        end if
                        !write(iulog,*) 'Negative biomas. in LakeBGCDynamics c,j,deficit (mol):',c,j,deficit
                     end if
                     biomas_new(c,j,k) = 0._r8
                     bphyto_dep_tot(c) = bphyto_dep_tot(c) - deficit/dtime
                  end if
               end do

               ! set phytoplankton biomass for outputs
               do j = 1, nlevlak
                  biomas_phyto(c,j,k) = biomas_new(c,j,k)
                  totphytoc(c) = totphytoc(c) + biomas_phyto(c,j,k)*dzx(c,j)
               end do
            end do
         end do

         ! Do Balance Check
         do j = 1, nlevlak
            do fc = 1, num_lakec
               c = filter_lakec(fc)
            
               if (j==1) errbiomas(c) = 0._r8
               do k = 1, nphytolak
                  errbiomas(c) = errbiomas(c) + dzx(c,j)*(biomas_new(c,j,k)-biomas_old(c,j,k))
                  errbiomas(c) = errbiomas(c) - gpp_vr(c,j,k)*dzx(c,j)*dtime
                  errbiomas(c) = errbiomas(c) + bloss_vr(c,j,k)*dzx(c,j)*dtime
               end do
            end do
         end do
      
         do fc = 1, num_lakec
            c = filter_lakec(fc)
      
            errbiomas(c) = errbiomas(c) + bphyto_dep_tot(c)*dtime
            if (abs(errbiomas(c)) < 1.e-8_r8) then
               bphyto_dep_tot(c) = bphyto_dep_tot(c) - errbiomas(c)/dtime
            else ! errbiomass > 1e-8 mol / m^2 / timestep
               write(iulog,*)'Biomass Conservation Error in LakeBGCDynamics during c, errbiomas (gC/m^2.timestep)', &
                     c,errbiomas(c)
               g = col_pp%gridcell(c)
               write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
               call endrun(msg=' ERROR: Biomass Conservation Error in LakeBGCDynamics during deposition'//&
                     errMsg(__FILE__, __LINE__))
            end if
         end do

         ! update sediment C pools
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            ! OC deposition (terrestrial + aquatic)
            cdep = lake_cdep(c) / 3.1536e7_r8
            cdep_tmp = max(bphyto_dep_tot(c),0._r8) + cdep
            do k = 1, nphytolak
               cdep_tmp = cdep_tmp + bphyto_dead(c,k)
            end do
            ctot_dep(c) = ctot_dep(c) + cdep_tmp*dtime
            soilc_act_tot(c) = soilc_act_tot(c) + (cdep_tmp - 0.95_r8*cdep)*dtime 

            ! decomposition and transformation (only active pool changed)
            do j = 1, nlevsoi
               soilc_act_tot(c) = soilc_act_tot(c) - soilc_loss(c,j,actC)* &
                     dz(c,j)*dtime
            end do
            soilc_act_tot(c) = max(0._r8, soilc_act_tot(c))

            ! active C redistribution
            ztot_act = 0._r8
            do j = 1, nlevsoi
               if (ztot_act>=0.2_r8) then
                  exit
               end if
               ztot_act = ztot_act + dz(c,j)
               js = j
            end do
            do j = 1, nlevsoi
               if (j<=js) then
                  lake_soilc(c,j,actC) = soilc_act_tot(c) / ztot_act
               else
                  lake_soilc(c,j,actC) = 0._r8
               end if
            end do

            totsoilc(c) = 0._r8
            do j = 1, nlevsoi
               totsoilc(c) = totsoilc(c) + sum(lake_soilc(c,j,:))*dz(c,j)
            end do
         end do
      end do ! end substep cycling

      ! average fluxes/rates over full time step 
      do fc = 1, num_lakec
         c = filter_lakec(fc)

         ! average flux/rate variables
         ctot_dep(c)          = ctot_dep(c) / (nsubstep*dtime)
         ch4_sed_diff(c)      = ch4_sed_diff(c) / (nsubstep*dtime)
         ch4_sed_ebul(c)      = ch4_sed_ebul(c) / (nsubstep*dtime)
         ch4_surf_diff(c)     = ch4_surf_diff(c) / (nsubstep*dtime)
         ch4_surf_ebul(c)     = ch4_surf_ebul(c) / (nsubstep*dtime)
         ch4_prod_wat(c,:)    = ch4_prod_wat(c,:) / (nsubstep*dtime)
         ch4_prod_sed(c,:)    = ch4_prod_sed(c,:) / (nsubstep*dtime)
         ch4_oxid_wat(c,:)    = ch4_oxid_wat(c,:) / (nsubstep*dtime)
         ch4_oxid_sed(c,:)    = ch4_oxid_sed(c,:) / (nsubstep*dtime)
         hr_wat(c,:)          = hr_wat(c,:) / (nsubstep*dtime)
         hr_sed(c,:)          = hr_sed(c,:) / (nsubstep*dtime)
         lake_gpp(c,:)        = lake_gpp(c,:) / (nsubstep*dtime)
         lake_npp(c,:)        = lake_npp(c,:) / (nsubstep*dtime)

         ch4_surf_flux(c)     = 1.e-3_r8 * catomw * (ch4_surf_diff(c) + ch4_surf_ebul(c))

         do j = 1, nlevlak+nlevsoi
            if (j<=nlevlak) then
               ch4_prod_tot(c) = ch4_prod_tot(c) + ch4_prod_wat(c,j)*dz_lake(c,j)
               ch4_oxid_tot(c) = ch4_oxid_tot(c) + ch4_oxid_wat(c,j)*dz_lake(c,j)
               lake_gpp_tot(c) = lake_gpp_tot(c) + lake_gpp(c,j)*dz_lake(c,j)
               lake_npp_tot(c) = lake_npp_tot(c) + lake_npp(c,j)*dz_lake(c,j)
            else
               ch4_prod_tot(c) = ch4_prod_tot(c) + ch4_prod_sed(c,j-nlevlak)*dz(c,j-nlevlak)
               ch4_oxid_tot(c) = ch4_oxid_tot(c) + ch4_oxid_sed(c,j-nlevlak)*dz(c,j-nlevlak)
            end if
         end do
         ! Adjustment to NEE flux to atm. for methane production/oxidation
         nem_col(c) = (ch4_oxid_tot(c) - ch4_prod_tot(c)) * catomw
      end do

      ! Now average up to gridcell for fluxes
      call c2g( bounds, &
           nem_col(bounds%begc:bounds%endc), nem_lake_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )

      end associate
   end subroutine LakeBGCDynamics

   !-----------------------------------------------------------------------
   real(r8) function HenrySolubility(temp, pH, gas)
      !
      ! !DESCRIPTION:
      ! Calculate gas solubility in water (mol/m3/Pa) 
      !
      ! !USES:
      use elm_varcon         , only : tfrz
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak 
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      real(r8), intent(in) :: pH       ! pH
      integer,  intent(in) :: gas      ! gas index
      !
      ! !LOCAL VARIABLES:
      real(r8), parameter :: SOLO2(11) = (/14.6_r8, 12.8_r8, 11.3_r8, 10.1_r8, &
                              9.1_r8, 8.3_r8, 7.6_r8, 7.0_r8, 6.5_r8, 6.0_r8, 5.6_r8/)
      real(r8) :: hi, kc1, kc2, par
      integer  :: indx 
      !--------------------------------------------------------------------

      if (gas==wn2lak) then
         HenrySolubility = 6.1e-6_r8*exp(-1300._r8*(1._r8/temp-1._r8/298._r8)) 
      else if (gas==wo2lak) then
         if (temp>=tfrz .and. temp<=tfrz+50._r8) then 
            indx = min( int((temp-tfrz)/5._r8) + 1, 10 )
            par = (temp - tfrz - 5*indx + 5) / 5._r8
            HenrySolubility = ( SOLO2(indx+1) * par + &
               SOLO2(indx) * (1._r8-par) ) / 6.596e5_r8 
         else
            HenrySolubility = 1.3e-5_r8*exp(-1500._r8*(1._r8/temp-1._r8/298._r8))
         end if
      else if (gas==wco2lak) then
         hi = 10._r8**(-pH)    ! Concentration of hydrogen ion
         ! rate constant of dissolved CO2 for first and second dissolution
         kc1 = 4.3e-7_r8*exp(-921.4_r8*(1._r8/temp-1._r8/298._r8))
         kc2 = 4.7e-11_r8*exp(-1787.4_r8*(1._r8/temp-1._r8/298._r8))
         HenrySolubility = 3.4e-4_r8*exp(-2400._r8*(1._r8/temp-1._r8/298._r8))* &
            (1._r8 + kc1/hi + kc1*kc2/hi**2._r8)
      else if (gas==wch4lak) then
         HenrySolubility = 1.3e-5_r8*exp(-1700._r8*(1._r8/temp-1._r8/298._r8))
      else
         call endrun(msg=' ERROR: no gas type index is found.'//&
               errMsg(__FILE__, __LINE__))
      end if

   end function HenrySolubility

   !-----------------------------------------------------------------------
   real(r8) function BubbleGasDiffusivity(temp, gas)
      !
      ! !DESCRIPTION:
      ! Calculate bubble gas diffusivity (m^2/s) 
      !
      ! !USES:
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      integer,  intent(in) :: gas      ! gas index
      !
      ! !LOCAL VARIABLES:
      !--------------------------------------------------------------------

      if (gas==wn2lak) then
         BubbleGasDiffusivity = 2.57e-7_r8 * (temp/273.0_r8)
      else if (gas==wo2lak) then
         BubbleGasDiffusivity = 2.4e-7_r8 * (temp/298.0_r8)
      else if (gas==wco2lak) then
         BubbleGasDiffusivity = 1.81e-3_r8 * exp(-2032.6_r8/temp)
      else if (gas==wch4lak) then
         BubbleGasDiffusivity = 1.5e-7_r8 * (temp/298.0_r8)
      else
         call endrun(msg=' ERROR: no gas type index is found.'//&
               errMsg(__FILE__, __LINE__))
      end if

   end function BubbleGasDiffusivity

   !-----------------------------------------------------------------------
   real(r8) function SchmidtNumber(temp, gas)
      !
      ! !DESCRIPTION:
      ! Calculate gas Schmidt number (m^2/s) 
      !
      ! !USES:
      use elm_varcon         , only : tfrz
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      integer,  intent(in) :: gas      ! gas index
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw
      !--------------------------------------------------------------------

      tw = min(temp-tfrz, 30._r8)
      if (tw<0._r8) then
         SchmidtNumber = 1.e30_r8
         return
      end if
      if (gas==wn2lak) then
         SchmidtNumber = 1970.7_r8 - 131.45_r8*tw + 4.139_r8*(tw**2._r8) - &
            0.052106_r8*(tw**3.0_r8)
      else if (gas==wo2lak) then
         SchmidtNumber = 1800.6_r8 - 120.1_r8*tw + 3.7818_r8*(tw**2._r8) - &
            0.047608_r8*(tw**3.0_r8)
      else if (gas==wco2lak) then
         SchmidtNumber = 1911.0_r8 - 113.7_r8*tw + 2.967_r8*(tw**2._r8) - &
            0.02943_r8*(tw**3.0_r8)
      else if (gas==wch4lak) then
         SchmidtNumber = 1898_r8 - 110.1_r8*tw + 2.834_r8*(tw**2._r8) - &
            0.02791_r8*(tw**3.0_r8)
      else
         call endrun(msg=' ERROR: no gas type index is found.'//&
               errMsg(__FILE__, __LINE__))
      end if

   end function SchmidtNumber

   !-----------------------------------------------------------------------
   real(r8) function SurfaceTension(temp)
      !
      ! !DESCRIPTION:
      ! Calculate bubble surface tension (N/m) 
      !
      ! !USES:
      use elm_varcon         , only : tfrz
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      !
      ! !LOCAL VARIABLES:
      real(r8), parameter :: r0 = 75.64e-3_r8
      real(r8), parameter :: r25 = 71.97e-3_r8
      real(r8), parameter :: r50 = 67.91e-3_r8
      real(r8) :: par
      !--------------------------------------------------------------------

      if (temp<tfrz) then
         SurfaceTension = 0._r8
      else if (temp>=tfrz .and. temp<tfrz+25._r8) then
         par = (temp-tfrz)/25.0_r8
         SurfaceTension = r0*(1-par) + r25*par
      else if (temp>=tfrz+25._r8 .and. temp<tfrz+50._r8) then
         par = (temp-tfrz-25.0_r8)/25.0_r8
         SurfaceTension = r25*(1-par) + r50*par
      else
         SurfaceTension = r50
      end if

   end function SurfaceTension

   !-----------------------------------------------------------------------
   real(r8) function KinematicViscosity(temp)
      !
      ! !DESCRIPTION:
      ! Calculate water kinematic viscosity (m^2/s) 
      !
      ! !USES:
      use elm_varcon         , only : tfrz, denh2o
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      !
      ! !LOCAL VARIABLES:
      !--------------------------------------------------------------------

      if (temp<tfrz) then
         KinematicViscosity = 1.e30_r8 
      else
         KinematicViscosity = 1._r8/denh2o * ( 2.414e-5_r8 * &
            10.0_r8**(247.8_r8/(temp-140.0_r8)) )
      end if

   end function KinematicViscosity 

   !-----------------------------------------------------------------------
   real(r8) function BuoyantVelocity(radius, vsc)
      !
      ! !DESCRIPTION:
      ! Calculate bubble rising velocity (m/s) 
      !
      ! !USES:
      use elm_varcon         , only : grav
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: radius   ! bubble radius (m)
      real(r8), intent(in) :: vsc      ! viscosity (m^2/s)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: xx, yy 
      !--------------------------------------------------------------------

      if (vsc>1.0e10_r8) then
         BuoyantVelocity = 1.e-30_r8
         return
      end if
      
      xx = grav * (radius**3.0_r8) / (vsc**2.0_r8)
      yy = 10.82_r8 / xx
      BuoyantVelocity = (2.0_r8*(radius**2.0_r8)*grav/9.0_r8/vsc) * &
         ((yy**2.0_r8+2.0_r8*yy)**0.5_r8-yy)

   end function BuoyantVelocity 

   !-----------------------------------------------------------------------
   real(r8) function GasFluxVelocity(wind, temp, gas)
      !
      ! !DESCRIPTION:
      ! Calculate gas flux velocity (m/s) 
      !
      ! !USES:
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: wind     ! wind speed (m/s)
      real(r8), intent(in) :: temp     ! temperature (K)
      integer,  intent(in) :: gas      ! gas id
      !
      ! !LOCAL VARIABLES:
      real(r8) :: schmidt, k600
      !--------------------------------------------------------------------

      schmidt = SchmidtNumber(temp, gas)
      k600 = 2.778e-6_r8 * (2.07_r8 + 0.215_r8*wind**1.7_r8)
      GasFluxVelocity = k600 * (600.0_r8/Schmidt)**0.5_r8

   end function GasFluxVelocity

   !-----------------------------------------------------------------------
   subroutine CalcGasEQConc(lakebgc_vars, c, conc_eq)
      !
      ! !DESCRIPTION:
      ! Calculate equilibrium gas concentration (mol/m^3) 
      !
      ! !USES:
      use elm_varcon         , only : grav, patomw
      use elm_varpar         , only : ngaslak, nsolulak
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak, wsrplak
      ! !ARGUMENTS:
      implicit none
      type(lakebgc_type)     , intent(in)  :: lakebgc_vars
      integer                , intent(in)  :: c             ! column index
      real(r8)               , intent(out) :: conc_eq(1:nsolulak)
      ! !CONSTANTS
      real(r8), parameter :: Xn2 = 0.78_r8            ! N2 mixing ratio
      !
      ! !LOCAL VARIABLES:
      real(r8) :: henry, temp, pressure
      integer  :: k, t
      !--------------------------------------------------------------------

      associate(                                            &
            t_lake         => col_es%t_lake                       , & ! Input: [real(r8) (:,:)] col lake temperature [Kelvin]

            forc_pbot      => top_as%pbot                         , & ! Input: [real(r8) (:)] atmospheric pressure [Pa]
            forc_po2       => top_as%po2bot                       , & ! Input: [real(r8) (:)] O2 partial pressure [Pa]
            forc_pco2      => top_as%pco2bot                      , & ! Input: [real(r8) (:)] CO2 partial pressure [Pa]      
            forc_pch4      => top_as%pch4bot                      , & ! Input: [real(r8) (:)] CH4 partial pressure [Pa]

            lake_tp        => lakebgc_vars%tp_col                 , & ! Input: [real(r8) (:)] epilimnion average total phosphorus [gP/m3]
            lake_ph        => lakebgc_vars%ph_col                 & ! Input: [real(r8) (:)] lake mean water pH
            )

      t = col_pp%topounit(c)
      temp = t_lake(c,1)
      do k = 1, nsolulak
         if (k==wn2lak) then
            pressure = Xn2*forc_pbot(t)
            henry = HenrySolubility(temp, lake_ph(c), wn2lak)
            conc_eq(k) = henry * pressure 
         else if (k==wo2lak) then
            pressure = forc_po2(t)
            henry = HenrySolubility(temp, lake_ph(c), wo2lak)
            conc_eq(k) = henry * pressure
         else if (k==wco2lak) then
            pressure = forc_pco2(t)
            henry = HenrySolubility(temp, lake_ph(c), wco2lak)
            conc_eq(k) = henry * pressure
         else if (k==wch4lak) then
            pressure = forc_pch4(t)
            henry = HenrySolubility(temp, lake_ph(c), wch4lak)
            conc_eq(k) = henry * pressure
         else if (k==wsrplak) then
            conc_eq(k) = 0.1_r8*lake_tp(c)/patomw
         else
            call endrun(msg=' ERROR: no solute index is found.'//&
                        errMsg(__FILE__, __LINE__)) 
         end if 
      end do

      end associate
   end subroutine CalcGasEQConc

   subroutine DefineWaterLayers(lakestate_vars, c, jtop, jwat, jmix, &
                                zTopblayer, zBotblayer, Hcline)
      !
      ! !DESCRIPTION:
      ! Define the indices or thickness of different water layers. 
      !
      ! !USES:
      use elm_varpar         , only : nlevlak
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)  :: lakestate_vars
      integer                , intent(in)  :: c             ! column index
      integer                , intent(out) :: jtop          ! top index of water column
      integer                , intent(out) :: jwat          ! top index of unfrozon water layers
      integer                , intent(out) :: jmix          ! bottom index of surface mixing layer
      real(r8)               , intent(out) :: zTopblayer    ! top position of thermocline
      real(r8)               , intent(out) :: zBotblayer    ! bottom position of thermocline
      real(r8)               , intent(out) :: Hcline        ! middle position of thermocline 
      !
      ! !CONSTANTS
      !
      ! !LOCAL VARIABLES:
      real(r8) :: zsum
      integer  :: j
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                      , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                       , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                    , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                       , & ! Input: [real(r8) (:,:)] col lake temperature [Kelvin]

            hmix           => lakestate_vars%lake_hmix_col        , & ! Input: [real(r8) (:)] lake mixing layer depth [m]
            icethick       => lakestate_vars%lake_icethick_col    & ! Input: [real(r8) (:)] ice thickness [m] (integrated if lakepuddling)
            )

      jtop = 1    ! not include snow layers
      jwat = COUNT(z_lake(c,:)+0.5_r8*dz_lake(c,:)<=icethick(c)) + 1

      ! mixing layer bottom index
      zsum = 0._r8
      jmix = 1
      do j = 1, nlevlak
         if (zsum>=hmix(c)) then
            exit
         end if
         zsum = zsum + dz_lake(c,j)
         jmix = j
      end do

      ! top boundary layer lower position 
      zTopblayer = sum(dz_lake(c,1:jwat-1))
      do j = jwat, nlevlak, 1
         if (abs(t_lake(c,j)-t_lake(c,jwat))>0.5_r8) then
            exit
         end if
         zTopblayer = zTopblayer + dz_lake(c,j)
      end do
      ! bottom boundary layer upper position
      zBotblayer = lakedepth(c)
      do j = nlevlak, jwat, -1
         if (abs(t_lake(c,j)-t_lake(c,nlevlak))>0.5_r8) then
            exit
         end if
         zBotblayer = zBotblayer - dz_lake(c,j)
      end do
      ! thermocline depth
      if (zTopblayer<zBotblayer) then
         Hcline = 0.5_r8 * (zTopblayer + zBotblayer)
      else
         zTopblayer = lakedepth(c)
         zBotblayer = lakedepth(c)
         Hcline = lakedepth(c)
      end if       

      end associate
   end subroutine DefineWaterLayers

   subroutine CalcSoluteDiffusivity(lakestate_vars, c, dx)
      !
      ! !DESCRIPTION:
      ! Calculate the diffusivity of solutes in both water and sediment. 
      !
      ! !USES:
      use elm_varpar         , only : nlevlak, nlevsoi
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)  :: lakestate_vars
      integer                , intent(in)  :: c                         ! column index
      real(r8)               , intent(out) :: dx(1:nlevlak+nlevsoi)     ! diffusivity (m^2/s)
      !
      ! !CONSTANTS
      !
      ! !LOCAL VARIABLES:
      real(r8) :: kmg(1:nlevlak+nlevsoi)     ! diffusivity (m2/s)
      real(r8) :: dzp
      integer  :: j
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                      , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            dz             => col_pp%dz                           , & ! Input: [real(r8) (:,:)] layer thickness [m]
            z              => col_pp%z                            , & ! Input: [real(r8) (:,:)] layer depth [m]

            kme            => lakestate_vars%lake_kme_col         , & ! Input: [real(r8) (:,:)] lake water and sediment heat diffusivity [m2/s]
            icefrac        => lakestate_vars%lake_icefrac_col     & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            )

      do j = 1, nlevlak+nlevsoi
         if (j<nlevlak) then
            kmg(j) = 1.25_r8 * kme(c,j)
         else if (j==nlevlak) then
            kmg(j) = max(1.25_r8*kme(c,j), 1.5e-6_r8)
         else
            kmg(j) = 2.57e-7_r8
         end if
      end do

      do j = 1, nlevlak+nlevsoi
         if (j<=nlevlak) then
            if (icefrac(c,j)<1._r8) then
               if (j<nlevlak) then
                  dx(j) = ( kmg(j)*kmg(j+1) * (dz_lake(c,j+1)+dz_lake(c,j)) ) &
                           / ( kmg(j)*dz_lake(c,j+1) + kmg(j+1)*dz_lake(c,j) )
               else
                  dzp = 0.5_r8 * (dz_lake(c,j) + dz(c,1))
                  dx(j) = ( kmg(j+1)*kmg(j)*dzp / &
                          (kmg(j+1)*dz_lake(c,j)/2._r8 + kmg(j)*z(c,1) ) )
               end if
            else
               dx(j) = 0._r8
            end if
         else
            dx(j) = kmg(j)
         end if
      end do

      end associate
   end subroutine CalcSoluteDiffusivity

   subroutine PhytoSettleVelocity(lakestate_vars, c, jwat, zTopblayer, &
                                  zBotblayer, conc_srp, vpx)
      !
      ! !DESCRIPTION:
      ! Calculate phytoplankton settling velocity, including phototaxis 
      ! and chemotaxis driven swim
      !
      ! !USES:
      use elm_varpar         , only : nlevlak, nphytolak 
      use LakeBGCType        , only : small_phyto, large_phyto
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)  :: lakestate_vars
      integer                , intent(in)  :: c             ! column index
      integer                , intent(in)  :: jwat          ! top index of unfrozon water layers
      real(r8)               , intent(in)  :: zTopblayer    ! top position of thermocline
      real(r8)               , intent(in)  :: zBotblayer    ! bottom position of thermocline
      real(r8)               , intent(in)  :: conc_srp(1:nlevlak) ! SRP concentration (mol/m3)
      real(r8)               , intent(out) :: vpx(1:nlevlak,1:nphytolak) ! settling velocity
      !
      ! !CONSTANTS
      real(r8), parameter :: Vset_small = 9.84e-8_r8  ! settling velocity of small phytoplankton (m/s)
      real(r8), parameter :: Vset_large = 8.68e-6_r8  ! settling velocity of large phytoplankton (m/s)
      real(r8), parameter :: Vswim = 2.89e-6_r8       ! swim velocity of small phytoplankton (m/s)
      real(r8), parameter :: fsrp_vmdown = 0.67_r8    ! threshold of nutrient limitation for chemotaxis
      real(r8), parameter :: fsrp_vmup = 0.75_r8      ! threshold of nutrient limition for phototaxis
      real(r8), parameter :: ipar_crit = 0.1_r8       ! threshold of light for phototaxis (umol/m2/s)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: fsrp, ipar
      integer  :: j
      !-------------------------------------------------------------------- 

      associate(                                            &
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            fsds_vis       => lakestate_vars%lake_fsds_vis_col & ! Input: [real(r8) (:,:)] incident vis radiation (W/m^2)
            )

      do j = 1, nlevlak
         ipar = 4.6_r8 * fsds_vis(c,j) ! umol/m2/s
         if (j<=jwat-1) then
            vpx(j,small_phyto) = 0._r8 
            vpx(j,large_phyto) = 0._r8
         else
            if (z_lake(c,j)<=zTopblayer) then
               ! surface mixing layer
               vpx(j,small_phyto) = 0._r8
               vpx(j,large_phyto) = 0._r8
            else if (z_lake(c,j)<zBotblayer) then
               ! thermocline
               fsrp = conc_srp(j) / (LakeBGCParamsInst%Ksrps + conc_srp(j))
               if (fsrp<fsrp_vmdown) then
                  vpx(j,small_phyto) = Vswim
               else if (fsrp>fsrp_vmup .and. ipar>ipar_crit) then
                  vpx(j,small_phyto) = -Vswim
               else
                  vpx(j,small_phyto) = 0._r8
               end if
               fsrp = conc_srp(j) / (LakeBGCParamsInst%Ksrpl + conc_srp(j)) 
               if (fsrp<fsrp_vmdown) then
                  vpx(j,large_phyto) = Vset_large
               else
                  vpx(j,large_phyto) = 0._r8
               end if
            else
               ! bottom mixing layer
               vpx(j,small_phyto) = Vset_small
               vpx(j,large_phyto) = Vset_large
            end if
         end if
         if (j==jwat) then
            vpx(j,small_phyto) = max(vpx(j,small_phyto),0._r8)
            vpx(j,large_phyto) = max(vpx(j,large_phyto),0._r8)
         end if
      end do

      end associate
   end subroutine PhytoSettleVelocity 

   subroutine SedimentEbullition(soilstate_vars, lakebgc_vars, c, conc_gas, gebul) 
      !
      ! !DESCRIPTION:
      ! Simulate bubble generation 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak
      ! !ARGUMENTS:
      implicit none
      type(soilstate_type)   , intent(in)  :: soilstate_vars
      type(lakebgc_type)     , intent(in)  :: lakebgc_vars
      integer                , intent(in)  :: c                            ! column index
      real(r8)               , intent(in)  :: conc_gas(1:nlevlak+nlevsoi,1:ngaslak) ! dissolved gas (mol/m3)
      real(r8)               , intent(out) :: gebul(1:nlevsoi,1:ngaslak)   ! ebullition rate (mol/m3/s)
      !
      ! !CONSTANTS
      !
      ! !LOCAL VARIABLES:
      real(r8) :: temp, por, depth
      real(r8) :: c_n2, c_ch4
      real(r8) :: Hn2, Hch4
      real(r8) :: Prn2, Prch4, Prtot, Prnet, Psat 
      integer  :: j, t
      !-------------------------------------------------------------------- 

      associate(                                            &
            z              => col_pp%z                         , & ! Input: [real(r8) (:,:)] layer depth [m]
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_soisno       => col_es%t_soisno                  , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            forc_pbot      => top_as%pbot                      , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)

            soilpor        => lakebgc_vars%soilpor_col         , & ! Input: [real(r8) (:,:)] sediment porosity 
            lake_ph        => lakebgc_vars%ph_col              & ! Input: [real(r8) (:)] lake column mean water pH      
            )

      t = col_pp%topounit(c)
      
      do j = 1, nlevsoi
         temp = t_soisno(c,j)
         por = soilpor(c,j)
         depth = lakedepth(c) + z(c,j)
         c_n2 = conc_gas(j+nlevlak,wn2lak)
         c_ch4 = conc_gas(j+nlevlak,wch4lak)
         Hn2 = HenrySolubility(temp, lake_ph(c), wn2lak)
         Hch4 = HenrySolubility(temp, lake_ph(c), wch4lak)
         Prn2 = c_n2/Hn2                  ! partial pressure of N2
         Prch4 = c_ch4/Hch4               ! partial pressure of CH4
         Prtot = Prn2 + Prch4
         Psat = LakeBGCParamsInst%Ae * por * (forc_pbot(t) + grav*denh2o*depth)
         Prnet = Prtot - Psat             ! excessive pressure
         if (Prnet>1.e-8_r8) then
            gebul(j,wn2lak) = LakeBGCParamsInst%Re * (Prn2/Prtot) * Prnet * Hn2
            gebul(j,wch4lak) = LakeBGCParamsInst%Re * (Prch4/Prtot) * Prnet * Hch4
            gebul(j,wo2lak) = 0.0_r8
            gebul(j,wco2lak) = 0.0_r8
         else
            gebul(j,1:ngaslak) = 0.0_r8
         end if
      end do

      end associate
   end subroutine SedimentEbullition

   !-----------------------------------------------------------------------
   subroutine BubbleDynamics(lakestate_vars, lakebgc_vars, c, jwat, &
                             isFullTStep, dtime, conc_gas, ebb, ebt, &
                             ex, ex_iceb, conc_bubl, conc_iceb) 
      !
      ! !DESCRIPTION:
      ! Simulate bubble transport in the water column 
      !
      ! !USES:
      use elm_varcon         , only : tfrz, denh2o, grav, rgas, rpi
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak
      use LakeBGCType        , only : wn2lak, wo2lak, wco2lak, wch4lak
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(in)    :: lakebgc_vars
      integer                , intent(in)    :: c
      integer                , intent(in)    :: jwat                    ! unfrozen water top index
      logical                , intent(in)    :: isFullTStep             ! flag for full time step
      real(r8)               , intent(in)    :: dtime                   ! timestep size [seconds]
      real(r8)               , intent(in)    :: conc_gas(1:nlevlak+nlevsoi,1:ngaslak)  ! dissolved gas (mol/m3)
      real(r8)               , intent(in)    :: ebb(1:ngaslak)          ! bottom ebullition (mol/m2/s)
      real(r8)               , intent(out)   :: ebt(1:ngaslak)          ! surface ebullition (mol/m2/s)
      real(r8)               , intent(inout) :: ex(1:nlevlak,1:ngaslak) ! bubble gas dissolution (mol/m3/s)
      real(r8)               , intent(inout) :: ex_iceb(1:ngaslak)      ! bubble gas oxidation (mol/m2/s)
      real(r8)               , intent(inout) :: conc_bubl(1:nlevlak,1:ngaslak)   ! bubble gas concentration (mol/m3)
      real(r8)               , intent(inout) :: conc_iceb(1:ngaslak)    ! gas concentration of ice-trapped bubbles (mol/m2)
      !
      ! !CONSTANTS:
      integer,  parameter :: nrlev = 11   ! bubble size number
      real(r8), parameter :: Rb0(11) = (/2.5_r8, 3.25_r8, 4.0_r8, 4.75_r8, 5.5_r8, &
                              6.25_r8, 7.0_r8, 7.75_r8, 8.5_r8, 9.25_r8, 10.0_r8/)  ! mm
      !
      ! !LOCAL VARIABLES:
      integer  :: rindx, locindx
      integer  :: j, k, t
      real(r8) :: pos, pr, vb, temp
      real(r8) :: rr, Atmp
      real(r8) :: dr_tmp1, dr_tmp2, dt
      real(r8) :: numb, pr_tmp
      real(r8) :: vsc(1:nlevlak)                                  ! viscosity
      real(r8) :: gama(1:nlevlak)                                 ! bubble tension
      real(r8) :: peclet(1:ngaslak)                               ! Peclet number
      real(r8) :: reynold(1:ngaslak)                              ! Reynold number
      real(r8) :: nusselt(1:ngaslak)                              ! Nusselt number
      real(r8) :: k_ndm(1:ngaslak)                                ! gas transfer coefficient
      real(r8) :: ndm(1:ngaslak)                                  ! gas exchange rate of single bubble
      real(r8) :: vbg(1:ngaslak)                                  ! bubble gas conc (mol/m^3)
      real(r8) :: exlayer(1:ngaslak)                              ! gas exchange at layer level (mol/m2/s)
      real(r8) :: vrr(1:nlevlak+1,1:nrlev)                        ! bubble radius at each lake level (m)
      real(r8) :: schmidt(1:nlevlak,1:ngaslak)                    ! Schmidt number
      real(r8) :: gdiff(1:nlevlak,1:ngaslak)                      ! gas diffusivity
      real(r8) :: gsolu(1:nlevlak,1:ngaslak)                      ! gas solubility
      real(r8) :: vb_rr(1:nlevlak+1,1:nrlev)                      ! bubble rising velocity (m/s)
      real(r8) :: bgas_rr(1:nlevlak+1,1:ngaslak,1:nrlev)          ! size-resolved bubble gas conc (mol/m^3)
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 

            forc_pbot      => top_as%pbot                      , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            lake_ph        => lakebgc_vars%ph_col              & ! Input: [real(r8) (:)] lake column mean water pH      
            )

      t = col_pp%topounit(c)

      ! only run the bubble model at the full time step
      if (isFullTStep) then
         ! Calculate bubble release rate at the water-sediment interface
         do j = 1, nlevlak
            temp = max(t_lake(c,j),tfrz)
            gama(j) = SurfaceTension(temp) 
            vsc(j) = KinematicViscosity(temp)
            do k = 1, ngaslak
               gsolu(j,k) = HenrySolubility(temp, lake_ph(c), k)
               gdiff(j,k) = BubbleGasDiffusivity(temp, k)
               schmidt(j,k) = SchmidtNumber(temp, k)
            end do
         end do
         pr = forc_pbot(t) + denh2o*grav*lakedepth(c)
         Atmp = 0._r8
         do rindx = 1, nrlev
            rr = 1.e-3_r8 * Rb0(rindx)
            vrr(nlevlak+1,rindx) = rr
            vb = BuoyantVelocity(rr, vsc(nlevlak)) 
            vb_rr(nlevlak+1,rindx) = vb
            Atmp = Atmp + (pr*rr + 2._r8*gama(nlevlak)) * vb
            bgas_rr(nlevlak+1,:,rindx) = ebb * (pr*rr+2._r8*gama(nlevlak))
         end do
         bgas_rr(nlevlak+1,:,:) = bgas_rr(nlevlak+1,:,:) / Atmp

         ! time step (s)
         dt = 4.0_r8 / (0.1418_r8*(10.0*Rb0(1))**2._r8 + &
            0.05579_r8*(10._r8*Rb0(1)) + 0.7794_r8)

         do rindx = 1, nrlev
            ! no sediment bubbling
            if (ebb(wch4lak)<1.e-10_r8) then
               bgas_rr(1:nlevlak,:,rindx) = 0._r8
               vb_rr(1:nlevlak,rindx) = 0._r8
               vrr(1:nlevlak,rindx) = 1.e-3_r8 * Rb0(rindx)
               cycle
            end if
            locindx = nlevlak
            pos = -lakedepth(c) 
            rr = 1.e-3_r8 * Rb0(rindx)
            vbg = bgas_rr(nlevlak+1,:,rindx)
            ! keep bubble number constant (missing some nonlinear effect)
            temp = t_lake(c,locindx)
            pr = forc_pbot(t) - denh2o*grav*pos + 2._r8*gama(locindx)/rr
            numb = 0.75e-3_r8*rgas*temp/(rpi*pr*rr**3.0_r8)*sum(vbg)
            do while (pos<-icethick(c) .and. rr>1.e-8_r8 .and. &
                     sum(vbg)>1.e-10_r8)
               if (pos>=-z_lake(c,locindx)+0.5_r8*dz_lake(c,locindx)) then
                  bgas_rr(locindx,:,rindx) = vbg 
                  vb_rr(locindx,rindx) = vb
                  vrr(locindx,rindx) = rr 
                  locindx = locindx - 1 
               end if
               temp = t_lake(c,locindx)
               vb = BuoyantVelocity(rr, vsc(locindx))
               pr = forc_pbot(t) - denh2o*grav*pos + 2._r8*gama(locindx)/rr
               !numb = 0.75e-3_r8*rgas*temp/(rpi*pr*rr**3.0)*sum(vbg)
               ! gas exchange rate
               peclet = rr * vb / gdiff(locindx,:)
               reynold = peclet / schmidt(locindx,:) 
               where (reynold<=1._r8)
                  nusselt = sqrt(2._r8*rpi*peclet/3.0_r8)
               elsewhere
                  nusselt = 0.45_r8*(reynold**(1._r8/6._r8))*(peclet**(1._r8/3._r8)) 
               end where
               k_ndm = -4._r8 * rpi * rr * gdiff(locindx,:) * nusselt
               ndm = k_ndm * (gsolu(locindx,:)*pr*vbg/sum(vbg) - conc_gas(locindx,:))
               ! new bubble gas
               vbg = vbg + numb * ndm * dt
               where (vbg<0._r8) vbg = 0._r8    ! remove small negatives
               pr_tmp = 3._r8*forc_pbot(t) - 3._r8*denh2o*grav*pos + &
                     4._r8*gama(locindx)/rr
               dr_tmp1 = 0.75e-3_r8*rgas*temp/(rpi*rr**2._r8)/pr_tmp 
               dr_tmp2 = rr*denh2o*grav*vb/pr_tmp
               ! new bubble radius
               rr = rr + (dr_tmp1*sum(ndm)+dr_tmp2)*dt
               pos = pos + vb * dt     ! new position
            end do 
            ! including ice layers
            do j = 1, locindx, 1
               bgas_rr(j,:,rindx) = vbg
               vb_rr(j,rindx) = vb
               vrr(j,rindx) = max(0._r8,rr)
            end do
         end do

         ! update state and flux variables 
         ebt = ebb
         conc_bubl(1:nlevlak,:) = 0._r8
         ex(1:nlevlak,:) = 0._r8
         do rindx = 1, nrlev
            ! size-integrated bubble gas and gas exchange
            do j = nlevlak, jwat, -1
               conc_bubl(j,:) = conc_bubl(j,:) + 0.5_r8 * &
                     (bgas_rr(j,:,rindx)+bgas_rr(j+1,:,rindx))
               exlayer = (bgas_rr(j+1,:,rindx) - bgas_rr(j,:,rindx)) * &
                     vb_rr(nlevlak+1,rindx)
               ex(j,:) = ex(j,:) + exlayer / dz_lake(c,j)
               ebt = ebt - exlayer
            end do
         end do 
      else
         ! update flux variable (assuming unchanged ex and conc_bubl)
         ebt = ebb
         ! size-integrated bubble gas and gas exchange
         do j = nlevlak, jwat, -1
            ebt = ebt - ex(j,:) * dz_lake(c,j)
         end do
      end if

      if (icethick(c)>1.e-8_r8) then
         ! bubbles trapped in ice layers
         conc_iceb = conc_iceb + ebt * dtime
         ebt = 0._r8
      else
         ! CH4 ebullition or loss from trapped ice bubbles
         ebt = ebt + conc_iceb * LakeBGCParamsInst%icebflux
         do k = 1, ngaslak
            conc_iceb(k) = conc_iceb(k) * max(0._r8, &
                  (1._r8 - LakeBGCParamsInst%icebflux*dtime))
         end do 
      end if

      ! CH4 loss in ice-trapped bubbles
      do k = 1, ngaslak
         if (k==wch4lak) then
            ex_iceb(k) = conc_iceb(k) * LakeBGCParamsInst%icebloss
            conc_iceb(k) = conc_iceb(k) * max(0._r8, &
                  (1._r8 - LakeBGCParamsInst%icebloss*dtime))
         else
            conc_iceb(k) = 0._r8
         end if
      end do

      end associate
   end subroutine BubbleDynamics

   !----------------------------------------------------------------------- 
   subroutine Photosynthesis(lakestate_vars, lakebgc_vars, c, conc_solu, &
                             biomas_phyto, csrc, csnk, gpp_vr, chla_vr)
      !
      ! !DESCRIPTION:
      ! Simulate phytoplankton photosynthesis 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz
      use elm_varcon         , only : catomw, patomw
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak, nphytolak
      use LakeBGCType        , only : small_phyto, large_phyto 
      use LakeBGCType        , only : wo2lak, wco2lak, wsrplak
      use LakeBGCType        , only : YC2P_POM
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(in)    :: lakebgc_vars
      integer                , intent(in)    :: c
      real(r8)               , intent(in)    :: conc_solu(1:nlevlak+nlevsoi,1:nsolulak)   ! solute (mol/m3)
      real(r8)               , intent(in)    :: biomas_phyto(1:nlevlak,1:nphytolak)       ! biomass (gC/m3)
      real(r8)               , intent(inout) :: csrc(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: csnk(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s   
      real(r8)               , intent(out)   :: gpp_vr(1:nlevlak,1:nphytolak)             ! gC/m3/s
      real(r8)               , intent(out)   :: chla_vr(1:nlevlak)                        ! g/m3
      !
      ! !CONSTANTS
      real(r8), parameter :: ThetaG = 1.08_r8      ! temperature multiplier for phytoplankton growth
      real(r8), parameter :: kt(nphytolak) = (/1.93841_r8, 12.76830_r8/)
      real(r8), parameter :: at(nphytolak) = (/29.27777_r8, 21.67022_r8/)
      real(r8), parameter :: bt(nphytolak) = (/0.28991_r8, 0.21632_r8/)
      real(r8), parameter :: Kco2 = 6.163e-2_r8    ! CO2 half-saturation constant (mol/m3)
      real(r8), parameter :: C2Chlmin(nphytolak) = (/1.44e2_r8, 1.44e2_r8/)   ! minimum carbon to chlorophill ratio (gC g chl-1)
      real(r8), parameter :: C2Chlmax(nphytolak) = (/4.8e2_r8, 3.6e2_r8/)     ! maximum carbon to chlorophill ratio (gC g chl-1)
      real(r8), parameter :: mu0(nphytolak) = (/0.4_r8, 1.0_r8/)  ! the maximum growth rate of phytoplankton at 0 celcius (d-1)
      real(r8), parameter :: Kpc2chl(nphytolak) = (/95._r8, 70._r8/)   ! The slope of C:Chl ratio vs. growth rate (gC g chl-1 d)
      real(r8), parameter :: phAlpha(nphytolak) = (/2.0e3_r8, 6.0e3_r8/)   ! initial slope of P-E curve ((mol photons m-2 s-1)-1 d-1)
      real(r8), parameter :: phBeta(nphytolak) = (/4.0e2_r8, 4.0e2_r8/)    ! photoinhibition parameter ((mol photons m-2 s-1)-1 d-1) 
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, c2chl, chla, ipar0, ipar
      real(r8) :: fpar, ftemp, fco2, fsrp
      real(r8) :: Vch, Ksrp
      real(r8) :: C2Chl0
      real(r8) :: gpp_phyto
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            fsds_vis       => lakebgc_vars%fsds_vis_col        & ! Input: [real(r8) (:,:)] incident vis radiation for BGC (W/m^2)
            )

      ipar0 = 4.6e-6_r8 * fsds_vis(c,1)   ! surface incident vis radiation (mol/m2/s) 

      do j = 1, nlevlak
         ! no radiation
         if (ipar0<1.e-8_r8) then
            chla_vr(j) = sum( biomas_phyto(j,:)/C2Chlmax )
            gpp_vr(j,:)  = 0._r8
            cycle
         end if

         ! ice layers
         if (icefrac(c,j)>=1._r8) then
            chla_vr(j) = 0._r8
            gpp_vr(j,:) = 0._r8
            cycle
         end if

         chla_vr(j) = 0._r8

         tw = t_lake(c,j) - tfrz
         ipar = 4.6e-6_r8 * fsds_vis(c,j)

         do k = 1, nphytolak
            if (k==small_phyto) then
               Vch = LakeBGCParamsInst%Vchs
               Ksrp = LakeBGCParamsInst%Ksrps
            else if (k==large_phyto) then
               Vch = LakeBGCParamsInst%Vchl
               Ksrp = LakeBGCParamsInst%Ksrpl
            end if

            ! temperature factor
            ftemp = max( 0.1_r8, ThetaG**(tw-20._r8) - &
                  ThetaG**(kt(k)*(tw-at(k))) + bt(k) )

            ! CO2 factor (not used now)
            fco2 = conc_solu(j,wco2lak) / (Kco2 + conc_solu(j,wco2lak)) 

            ! nutrient factor
            fsrp = conc_solu(j,wsrplak) / (Ksrp + conc_solu(j,wsrplak))

            ! C:Chla ratio
            C2Chl0 = C2Chlmax(k) - Kpc2chl(k) * mu0(k) * ftemp * fsrp
            if (icefrac(c,j)<1._r8 .and. ipar>0.01_r8*ipar0) then
               c2chl = C2Chl0 - (C2Chl0 - C2Chlmin(k)) * log(ipar0/ipar) / 4.605_r8
               c2chl = min(max(c2chl, C2Chlmin(k)), C2Chlmax(k))
            else
               c2chl = C2Chlmax(k)
            end if

            chla = biomas_phyto(j,k)*(1._r8-icefrac(c,j))/c2chl
            chla_vr(j) = chla_vr(j) + chla

            ! radiation factor
            fpar = (1._r8 - exp(-phAlpha(k)*ipar/mu0(k))) * exp(-phBeta(k)*ipar/mu0(k))

            gpp_phyto = Vch/8.64e4_r8 * fpar * ftemp * fsrp * chla / catomw
            gpp_vr(j,k)       = catomw * gpp_phyto

            csrc(j,wo2lak)    = csrc(j,wo2lak)  + gpp_phyto
            csnk(j,wco2lak)   = csnk(j,wco2lak) + gpp_phyto
            csnk(j,wsrplak)   = csnk(j,wsrplak) + gpp_phyto/YC2P_POM
         end do 

      end do
      
      end associate
   end subroutine Photosynthesis

   subroutine AutotrophicR(lakestate_vars, c, conc_solu, biomas_phyto, &
                           csrc, csnk, bloss_vr, ar_vr)
      !
      ! !DESCRIPTION:
      ! Simulate phytoplankton autotrophic respiration 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak, nphytolak
      use LakeBGCType        , only : small_phyto, large_phyto
      use LakeBGCType        , only : wo2lak, wco2lak, wsrplak
      use LakeBGCType        , only : YC2P_POM
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      integer                , intent(in)    :: c
      real(r8)               , intent(in)    :: conc_solu(1:nlevlak+nlevsoi,1:nsolulak)   ! solute (mol/m3)
      real(r8)               , intent(in)    :: biomas_phyto(1:nlevlak,1:nphytolak)       ! biomass (gC/m3)
      real(r8)               , intent(inout) :: csrc(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: csnk(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(out)   :: bloss_vr(1:nlevlak,1:nphytolak)           ! gC/m3/s
      real(r8)               , intent(out)   :: ar_vr(1:nlevlak,1:nphytolak)              ! gC/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: ThetaML = 1.073_r8 ! temperature multiplier for metabolic loss 
      real(r8), parameter :: frres(nphytolak) = (/0.9_r8, 0.75_r8/)  ! fraction of respiration in total metabolic loss
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, bloss 
      real(r8) :: Klr
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            )

      do j = 1, nlevlak
         ! ice layers
         if (icefrac(c,j)>=1._r8) then
            bloss_vr(j,:) = 0._r8
            ar_vr(j,:) = 0._r8
            cycle
         end if

         tw = t_lake(c,j) - tfrz

         do k = 1, nphytolak
            if (k==small_phyto) then
               Klr = LakeBGCParamsInst%Klrs
            else if (k==large_phyto) then
               Klr = LakeBGCParamsInst%Klrl
            end if
            
            bloss = Klr/8.64e4_r8 / catomw * (ThetaML**(tw-20)) * &
                  biomas_phyto(j,k) * (1._r8 - icefrac(c,j))
            bloss_vr(j,k) = bloss * catomw
            ar_vr(j,k) = bloss * catomw * frres(k)
            
            csnk(j,wo2lak)    = csnk(j,wo2lak)  + bloss*frres(k)
            csrc(j,wco2lak)   = csrc(j,wco2lak) + bloss*frres(k)
            csrc(j,wsrplak)   = csrc(j,wsrplak) + bloss*frres(k)/YC2P_POM
         end do 
      
      end do

      end associate
   end subroutine AutotrophicR

   subroutine HeterotrophicR(lakestate_vars, c, conc_solu, soilc, &
                             csrc, csnk, soilc_loss, hr_vr)
      !
      ! !DESCRIPTION:
      ! Simulate microbial heterotrophic respiration 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak, nsoilclak
      use LakeBGCType        , only : wo2lak, wco2lak, wsrplak
      use LakeBGCType        , only : pasC, actC, frac_resp
      use LakeBGCType        , only : YC2P_POM, YC2P_DOM
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      integer                , intent(in)    :: c
      real(r8)               , intent(in)    :: conc_solu(1:nlevlak+nlevsoi,1:nsolulak)   ! solute (mol/m3)
      real(r8)               , intent(in)    :: soilc(1:nlevsoi,1:nsoilclak)              ! soil C (gC/m3)
      real(r8)               , intent(inout) :: csrc(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: csnk(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: soilc_loss(1:nlevsoi,1:nsoilclak)         ! gC/m3/s
      real(r8)               , intent(out)   :: hr_vr(1:nlevlak+nlevsoi)                  ! gC/m3/s 
      !
      ! !CONSTANTS
      real(r8), parameter :: Ko2CM = 4.6875e-2_r8  ! O2 half-saturation constant for HR (mol/m3)
      real(r8), parameter :: ThetaCM = 1.073_r8    ! temperature multiplier for heterotrophic respiration
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, ts, c_o2
      real(r8) :: Rca, pco2_soilc, pco2_BOD
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            dz             => col_pp%dz                        , & ! Input: [real(r8) (:,:)] layer thickness (m)                   
            z              => col_pp%z                         , & ! Input: [real(r8) (:,:)] layer depth (m)
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 
            t_soisno       => col_es%t_soisno                  , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            )

      ! for water column
      do j = 1, nlevlak
         ! ice layers
         if (icefrac(c,j)>=1._r8) then
            hr_vr(j) = 0._r8
            cycle
         end if

         tw = t_lake(c,j) - tfrz
         c_o2 = conc_solu(j,wo2lak)

         pco2_BOD = 3.62e-9_r8 * (ThetaCM**(tw-20._r8)) * &
               (c_o2/(Ko2CM+c_o2)) * (1._r8 - icefrac(c,j)) 
         hr_vr(j) = pco2_BOD * catomw

         csrc(j,wco2lak)   = csrc(j,wco2lak) + pco2_BOD 
         csrc(j,wsrplak)   = csrc(j,wsrplak) + pco2_BOD/YC2P_DOM
         csnk(j,wo2lak)    = csnk(j,wo2lak)  + pco2_BOD
      end do

      ! for sediment
      do j = 1, nlevsoi
         ts = t_soisno(c,j) - tfrz
         c_o2 = conc_solu(j+nlevlak,wo2lak)

         hr_vr(j+nlevlak) = 0._r8
         do k = 1, nsoilclak
            if (k==actC) then
               Rca = LakeBGCParamsInst%Rcaact
            else if (k==pasC) then
               Rca = LakeBGCParamsInst%Rcapas
            end if
            pco2_soilc = Rca * (ThetaCM**(ts-20._r8)) * (c_o2/(Ko2CM+c_o2)) * &
                  soilc(j,k) / catomw
            hr_vr(j+nlevlak) = hr_vr(j+nlevlak) + pco2_soilc*frac_resp(k)*catomw
            soilc_loss(j,k) = soilc_loss(j,k) + pco2_soilc*catomw

            csrc(j+nlevlak,wco2lak) = csrc(j+nlevlak,wco2lak) + &
                  pco2_soilc*frac_resp(k)
            csrc(j+nlevlak,wsrplak) = csrc(j+nlevlak,wsrplak) + &
                  pco2_soilc*frac_resp(k)/YC2P_POM
            csnk(j+nlevlak,wo2lak)  = csnk(j+nlevlak,wo2lak)  + &
                  pco2_soilc*frac_resp(k)
         end do
      end do

      end associate
   end subroutine HeterotrophicR

   subroutine Methanogenesis(lakestate_vars, lakebgc_vars, c, conc_solu, &
                             soilc, soilc_old, csrc, csnk, soilc_loss, pch4_vr)
      !
      ! !DESCRIPTION:
      ! Simulate methane production 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak, nsoilclak
      use LakeBGCType        , only : wo2lak, wco2lak, wch4lak, wsrplak
      use LakeBGCType        , only : pasC, actC
      use LakeBGCType        , only : yedoma_lake, thaw_lake 
      use LakeBGCType        , only : YC2P_POM
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(in)    :: lakebgc_vars
      integer                , intent(in)    :: c
      real(r8)               , intent(in)    :: conc_solu(1:nlevlak+nlevsoi,1:nsolulak)   ! solute (mol/m3)
      real(r8)               , intent(in)    :: soilc(1:nlevsoi,1:nsoilclak)              ! soil C (gC/m3)
      real(r8)               , intent(in)    :: soilc_old(1:nlevsoi)                      ! soil C (gC/m3)
      real(r8)               , intent(inout) :: csrc(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s 
      real(r8)               , intent(inout) :: csnk(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: soilc_loss(1:nlevsoi,1:nsoilclak)         ! gC/m3/s
      real(r8)               , intent(out)   :: pch4_vr(1:nlevlak+nlevsoi)                ! mol/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: Tpr(nsoilclak) = (/273.15_r8,276.65_r8/)   ! CH4 production reference temperature (K)
      real(r8), parameter :: pHmin = 2.2_r8  ! minimum allowable pH for CH4 production
      real(r8), parameter :: pHmax = 9.0_r8  ! maximum allowable pH for CH4 production
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, ts
      real(r8) :: Rc, PQ10, etaO2
      real(r8) :: c_o2, fo2, ftemp, fph
      real(r8) :: pch4_soilc, pch4_soilc_yedoma
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 
            t_soisno       => col_es%t_soisno                  , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            lake_tp        => lakebgc_vars%tp_col              , & ! Input: [real(r8) (:)] lake column mean total phosphorus [gP/m3]
            lake_ph        => lakebgc_vars%ph_col              , & ! Input: [real(r8) (:)] lake column mean water pH
            lake_type      => lakebgc_vars%ltype_col           & ! Input: [integer  (:)] lake type identifier 
            )

      ! CH4 production in oxic water of oligotrophic and mesotrophic lakes
      !if (lake_tp(c)<0.024) then
      !   do j = 1, nlevlak
      !      if (icefrac(c,j)>=1._r8) pch4_vr(j) = 0._r8
      !      csrc(j,wch4lak) = csrc(j,wch4lak) + LakeBGCParamsInst%Rcoxic*csrc(j,wo2lak) 
      !      pch4_vr(j) = LakeBGCParamsInst%Rcoxic*csrc(j,wo2lak)
      !   end do
      !else
         csrc(1:nlevlak,wch4lak) = 0._r8
         pch4_vr(1:nlevlak) = 0._r8
      !end if

      ! pH factor (borrowed from wetland CH4 model)
      if (lake_ph(c)>pHmin .and. lake_ph(c)<pHmax) then
         fph = 10._r8**(-0.2235_r8*lake_ph(c)**2._r8 + 2.7727_r8*lake_ph(c) - 8.6_r8)
         fph = min(1._r8, max(0._r8, fph))
      else
         fph = 0._r8
      end if

      ! CH4 production in sediment
      do j = 1, nlevsoi
         ! sediment temperature
         ts = t_soisno(c,j)
         
         ! O2 suppression (Tang et al., 2010; Biogeosciences)
         etaO2 = LakeBGCParamsInst%etaO2
         c_o2 = conc_solu(j+nlevlak,wo2lak)
         fo2 = 1._r8 / (1._r8 + etaO2*c_o2)
         
         pch4_vr(j+nlevlak) = 0._r8
         ! yedoma thermokarst lakes
         if (lake_type(c)==yedoma_lake) then
            Rc = LakeBGCParamsInst%Rcold

            pch4_soilc_yedoma = 0.5_r8 * Rc * soilc_old(j) / catomw 
            pch4_vr(j+nlevlak) = pch4_vr(j+nlevlak) + pch4_soilc_yedoma

            csrc(j+nlevlak,wch4lak) = csrc(j+nlevlak,wch4lak) + pch4_soilc_yedoma
            csrc(j+nlevlak,wco2lak) = csrc(j+nlevlak,wco2lak) + pch4_soilc_yedoma
            csrc(j+nlevlak,wsrplak) = csrc(j+nlevlak,wsrplak) + &
                  2._r8*pch4_soilc_yedoma/YC2P_POM
         end if
         ! other lakes
         do k = 1, nsoilclak
            if (k==actC) then
               ! active C mainly through acetoclastic methanogenesis
               Rc = LakeBGCParamsInst%Rcact
               PQ10 = LakeBGCParamsInst%PQ10act
               ftemp = PQ10**(0.1_r8*(ts-Tpr(k)))
               
               pch4_soilc = 0.5_r8 * Rc * ftemp * fo2 * fph * soilc(j,k) / catomw
               pch4_vr(j+nlevlak) = pch4_vr(j+nlevlak) + pch4_soilc
               soilc_loss(j,k) = soilc_loss(j,k) + 2._r8*catomw*pch4_soilc

               csrc(j+nlevlak,wch4lak) = csrc(j+nlevlak,wch4lak) + pch4_soilc
               csrc(j+nlevlak,wco2lak) = csrc(j+nlevlak,wco2lak) + pch4_soilc
               csrc(j+nlevlak,wsrplak) = csrc(j+nlevlak,wsrplak) + &
                  2._r8*pch4_soilc/YC2P_POM
            else if (k==pasC) then
               ! passive C mainly through hydrogenotrophic methanogenesis
               Rc = LakeBGCParamsInst%Rcpas
               PQ10 = LakeBGCParamsInst%PQ10pas
               ftemp = PQ10**(0.1_r8*(ts-Tpr(k)))

               pch4_soilc = 0.25_r8 * Rc * ftemp * fo2 * fph * soilc(j,k) / catomw
               pch4_vr(j+nlevlak) = pch4_vr(j+nlevlak) + pch4_soilc
               soilc_loss(j,k) = soilc_loss(j,k) + 4._r8*catomw*pch4_soilc

               csrc(j+nlevlak,wch4lak) = csrc(j+nlevlak,wch4lak) + pch4_soilc
               csrc(j+nlevlak,wco2lak) = csrc(j+nlevlak,wco2lak) + 3._r8*pch4_soilc
               csrc(j+nlevlak,wsrplak) = csrc(j+nlevlak,wsrplak) + &
                  4._r8*pch4_soilc/YC2P_POM
            end if
         end do
         
      end do

      end associate
   end subroutine Methanogenesis 

   subroutine SetYedomaLakeOldC(lakebgc_vars, c)
      !
      ! !DESCRIPTION:
      ! Set Pleistocene-age C pool for yedoma lake 
      !
      ! !USES:
      use elm_varpar         , only : nlevsoi, nlevgrnd
      use LakeBGCType        , only : yedoma_lake, thaw_lake
      ! !ARGUMENTS:
      implicit none
      type(lakebgc_type)     , intent(in)    :: lakebgc_vars
      integer                , intent(in)    :: c
      !
      ! !CONSTANTS
      real(r8), parameter :: oldcarb0 = 29.3e3_r8  ! yedoma permafrost C (gC/m3)
      real(r8), parameter :: Ct = 0.77_r8    ! permafrost thawing rate (m/(yr^0.5))
      real(r8), parameter :: Rco = 4.2983e-10_r8   ! yedoma C reference loss rate (s-1)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: talik, Rc, tthaw 
      integer  :: j
      !-------------------------------------------------------------------- 

      associate(                                            &
            z                    => col_pp%z                   , & ! Input: [real(r8) (:,:)] layer depth (m)
            lakedepth            => col_pp%lakedepth           , & ! Input: [real(r8) (:)] column lake depth [m] 

            lake_type            => lakebgc_vars%ltype_col     , & ! Input: [integer  (:)] lake type identifier

            lake_soilc_old       => lakebgc_vars%soilc_old_col   & ! Output: [real(r8) (:,:)] column lake Pleistocene-age C (gC/m3) 
            )

      talik = (lakedepth(c) - 2._r8) / 0.75_r8  ! talik thickness 
      do j = 1, nlevgrnd
         ! yedoma thermokarst lakes
         if (lake_type(c)==yedoma_lake) then
            if (j<nlevsoi .and. talik>=z(c,j) .and. talik<z(c,j+1)) then
               tthaw = 3.1536e7_r8 * (talik**2._r8 - z(c,j)**2._r8) / Ct**2._r8
               lake_soilc_old(c,j) = oldcarb0 * exp(-Rco*tthaw) / 3._r8
            else if (j==nlevsoi .and. talik>=z(c,j)) then
               tthaw = 3.1536e7_r8 * (talik**2._r8 - z(c,j)**2._r8) / Ct**2._r8
               lake_soilc_old(c,j) = oldcarb0 * exp(-Rco*tthaw) / 3._r8
            else
               lake_soilc_old(c,j) = 0._r8
            end if
         else
            lake_soilc_old(c,j) = 0._r8
         end if
      end do 

      end associate
   end subroutine

   subroutine Methanotrophy(lakestate_vars, c, conc_solu, csrc, csnk, och4_vr)
      !
      ! !DESCRIPTION:
      ! Simulate methane oxidation
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz
      use elm_varpar         , only : nlevlak, nlevsoi
      use elm_varpar         , only : ngaslak, nsolulak
      use LakeBGCType        , only : wo2lak, wch4lak, wco2lak
      ! !ARGUMENTS:
      implicit none
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      integer                , intent(in)    :: c
      real(r8)               , intent(in)    :: conc_solu(1:nlevlak+nlevsoi,1:nsolulak)   ! solute (mol/m3)
      real(r8)               , intent(inout) :: csrc(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(inout) :: csnk(1:nlevlak+nlevsoi,1:nsolulak)        ! mol/m3/s
      real(r8)               , intent(out)   :: och4_vr(1:nlevlak+nlevsoi)                ! mol/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: Tor = 267.65_r8 ! CH4 oxidation reference temperature (K)
      real(r8), parameter :: Och4_cr = 1.7e-3_r8 ! critical CH4 concentration for O2-inhibitation (mol/m3)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, ts
      real(r8) :: Qch4, OQ10, betaCH4, lamdaO2
      real(r8) :: c_ch4, c_o2, och4_oxic
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 
            t_soisno       => col_es%t_soisno                  , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            )

      ! CH4 oxidation in the water column
      do j = 1, nlevlak
         ! ice layers
         if (icefrac(c,j)>=1._r8) then
            och4_vr(j) = 0._r8
            cycle
         end if
         
         tw = t_lake(c,j)
         c_ch4 = conc_solu(j,wch4lak)
         c_o2  = conc_solu(j,wo2lak)

         Qch4 = LakeBGCParamsInst%Qch4
         OQ10 = LakeBGCParamsInst%OQ10
         betaCH4 = LakeBGCParamsInst%betaCH4

         lamdaO2 = LakeBGCParamsInst%lamdaO2 / (1._r8 + c_ch4/Och4_cr)

         och4_oxic = (1._r8-icefrac(c,j)) * Qch4 * (OQ10**(0.1_r8*(tw-Tor))) * &
            (c_ch4**betaCH4) * exp(-lamdaO2*c_o2) * (1._r8 - exp(-180._r8*c_o2))
         och4_vr(j) = och4_oxic

         csnk(j,wch4lak) = csnk(j,wch4lak) + och4_oxic
         csnk(j,wo2lak)  = csnk(j,wo2lak)  + 2._r8*och4_oxic
         csrc(j,wco2lak) = csrc(j,wco2lak) + och4_oxic
      end do      

      ! no CH4 oxidation in sediment
      och4_vr(nlevlak+1:nlevlak+nlevsoi) = 0._r8

      end associate
   end subroutine Methanotrophy

   !----------------------------------------------------------------------- 
   subroutine readLakeBGCParams ( ncid )
      !
      ! !USES:
      use shr_kind_mod , only : r8 => shr_kind_r8
      use ncdio_pio    , only : file_desc_t,ncd_io
      use elm_varpar   , only : nsoilclak
      !
      ! !ARGUMENTS:
      implicit none
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      !
      ! !LOCAL VARIABLES:
      character(len=32)  :: subname = 'CH4ParamsType'
      character(len=100) :: errCode = '-Error reading in parameters file:'
      logical            :: readv ! has variable been read in or not
      real(r8)           :: tempr ! temporary to read in constant
      character(len=100) :: tString ! temp. var for reading
   !-----------------------------------------------------------------------
      tString='Vchs'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Vchs = tempr

      tString='Vchl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Vchl = tempr

      tString='Ksrps'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ksrps = tempr

      tString='Ksrpl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ksrpl = tempr

      tString='Klrs'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Klrs = tempr

      tString='Klrl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Klrl = tempr

      tString = 'frcResusps'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%frcResusps = tempr

      tString = 'frcResuspl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%frcResuspl = tempr

      tString='Re'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Re = tempr

      tString='Ae'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ae = tempr

      tString='csedDMP'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%csedDMP = tempr

      tString='icebflux'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%icebflux = tempr

      tString='icebloss'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%icebloss = tempr

      tString='Rcapas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcapas = tempr

      tString='Rcaact'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcaact = tempr

      tString='Rcpas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcpas = tempr

      tString='Rcact'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcact = tempr

      tString='Rcold'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcold = tempr

      tString='PQ10pas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%PQ10pas = tempr

      tString='PQ10act'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%PQ10act = tempr

      tString='etaO2'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%etaO2 = tempr

      tString='Qch4'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Qch4 = tempr

      tString='OQ10'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%OQ10 = tempr

      tString='betaCH4'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%betaCH4 = tempr

      tString='lamdaO2'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%lamdaO2 = tempr

      tString='Rcoxic'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcoxic = tempr

   end subroutine 

end module LakeBGCDynMod
