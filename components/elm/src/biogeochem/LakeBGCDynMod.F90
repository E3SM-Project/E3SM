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
      ! initial slope of P-E curve ((mol photons m-2 s-1)-1 d-1)
      real(r8) :: phAlphas
      real(r8) :: phAlphal
      ! photoinhibition parameter ((mol photons m-2 s-1)-1 d-1)
      real(r8) :: phBetas
      real(r8) :: phBetal
      ! half-saturation for phosphorus limitation (mol/m3)
      real(r8) :: Ktps
      real(r8) :: Ktpl
      ! metabolic loss rate (d-1)
      real(r8) :: Klrs
      real(r8) :: Klrl
      ! phytoplankton resuspension fraction
      real(r8) :: frcResusp
      ! ebullition rate (s-1) 
      real(r8) :: Re 
      ! ice bubble flux and dissolution rate (s-1)
      real(r8) :: icebflux
      real(r8) :: icebloss
      ! sediment OC decomposition rate through aerobic respiration (s-1) 
      real(r8) :: Rcapas
      real(r8) :: Rcaact
      ! sediment OC decomposition rate through methanogenesis (s-1)
      real(r8) :: Rcpas
      real(r8) :: Rcact
      ! sediment OC decomposition Q10 through methanogenesis
      real(r8) :: PQ10pas
      real(r8) :: PQ10act
      ! CH4 oxidation potential (mol/m3/s)
      real(r8) :: Qch4
      ! CH4 oxidation Q10
      real(r8) :: OQ10
      ! Michaelis-Menten constant (unit: mol/m3)
      real(r8) :: Kch4
      real(r8) :: Ko2
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
      use elm_varpar         , only : ngaslak, nphytolak, nsoilclak
      use elm_varcon         , only : vkc, grav, denh2o, tfrz, cnfac
      use elm_varcon         , only : catomw
      use elm_varctl         , only : iulog
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak
      use LakeBGCType        , only : small_phyto, large_phyto
      use LakeBGCType        , only : actC, pasC, nsubstep
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
      ! !CONSTANTS
      real(r8), parameter :: Xn2 = 0.78_r8            ! N2 mixing ratio
      real(r8), parameter :: Vset_small = 9.84e-8_r8  ! settling velocity of small phytoplankton (m/s)
      real(r8), parameter :: Vset_large = 8.68e-6_r8  ! settling velocity of large phytoplankton (m/s)
      ! !LOCAL VARIABLES:
      integer  :: c, t, g, fc                         ! indices
      integer  :: j, k, r, j0                         ! indices
      integer  :: it                                  ! time indices
      real(r8) :: dtime                               ! timestep size [seconds]
      real(r8) :: temp, por, depth                    ! auxiliary variables
      real(r8) :: cgas_avg, zsum                      ! auxiliary variables
      real(r8) :: cdep, zsoi_inc                      ! auxiliary variables
      real(r8) :: c_o2, c_n2, c_co2, c_ch4            ! concentration (mol/m3)
      real(r8) :: Hn2, Hch4                           ! solubility
      real(r8) :: Prn2, Prch4, Prtot, Prnet           ! partial pressure (Pa)
      real(r8) :: Psat, Pgas                          ! partial pressure (Pa)
      real(r8) :: dzm, dzp, vzm, vzp
      real(r8) :: deficit
      real(r8) :: Hcline, Hsum1, Hsum2
      integer  :: jtop(bounds%begc:bounds%endc)                      ! top level for each column
      integer  :: jwat(bounds%begc:bounds%endc)                      ! top water layer index
      integer  :: jmix(bounds%begc:bounds%endc)                      ! mixing layer index
      real(r8) :: errch4(bounds%begc:bounds%endc)                    ! error (Mol CH4 /m^2) [+ = too much CH4]
      real(r8) :: errbiomas(bounds%begc:bounds%endc)                 ! error (gC/m2) [+ = too much biomass]
      real(r8) :: bphyto_dep_tot(bounds%begc:bounds%endc)            ! total phytoplankton deposition (gC/m2/s)
      real(r8) :: bphyto_dep(bounds%begc:bounds%endc,1:nphytolak)    ! phytoplankton deposition (gC/m2/s)
      real(r8) :: bphyto_ssp(bounds%begc:bounds%endc,1:nphytolak)    ! phytoplankton suspension (gC/m3/s)
      real(r8) :: bphyto_dead(bounds%begc:bounds%endc,1:nphytolak)   ! vertically-integrated dead phytoplankton (gC/m2/s)
      real(r8) :: eb_sed(bounds%begc:bounds%endc,1:ngaslak)          ! gas ebullition from sediment (mol/m2/s)
      real(r8) :: df_top(bounds%begc:bounds%endc,1:ngaslak)          ! gas diffusion from lake surface (mol/m2/s)
      real(r8) :: conc_eq(bounds%begc:bounds%endc,1:ngaslak)         ! equilibrium gas concentration (mol/m^3)
      real(r8) :: conc_iceb_old(bounds%begc:bounds%endc,1:ngaslak)   ! ice-trapped bubbles at the last time step (mol/m2)
      real(r8) :: kg(bounds%begc:bounds%endc,1:ngaslak)              ! gas transfer velocity (m/s)
      real(r8) :: zx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! interface depth (+ below surface) for whole column (m)
      real(r8) :: dzx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)     ! cell thickness (+ below surface) for whole column (m)
      real(r8) :: qx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! porosity (+ below surface) for whole column
      real(r8) :: dx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! diffusivity (+ below surface) for whole column (m^2/s)
      real(r8) :: va(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "a" vector for tridiagonal matrix
      real(r8) :: vb(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "b" vector for tridiagonal matrix
      real(r8) :: vc(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "c" vector for tridiagonal matrix
      real(r8) :: vr(bounds%begc:bounds%endc,1:nlevlak+nlevsoi)      ! "r" vector for tridiagonal solution
      real(r8) :: ex(bounds%begc:bounds%endc,1:nlevlak,1:ngaslak)    ! bubble dissolution rate (mol/m3/s)
      real(r8) :: sx(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:ngaslak)        ! source (+ below surface) for whole column (mol/m^3/s)
      real(r8) :: gsrc(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:ngaslak)      ! dissolved gas production rate (mol/m3/s)
      real(r8) :: gsnk(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:ngaslak)      ! dissolved gas loss rate (mol/m3/s)
      real(r8) :: conc_old(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:ngaslak)  ! gas concentration at the last time step (mol/m3)
      real(r8) :: conc_new(bounds%begc:bounds%endc,1:nlevlak+nlevsoi,1:ngaslak)  ! updated gas concentration (mol/m3)
      real(r8) :: vba(bounds%begc:bounds%endc,1:nlevlak)             ! "a" vector for tridiagonal matrix
      real(r8) :: vbb(bounds%begc:bounds%endc,1:nlevlak)             ! "b" vector for tridiagonal matrix
      real(r8) :: vbc(bounds%begc:bounds%endc,1:nlevlak)             ! "c" vector for tridiagonal matrix
      real(r8) :: vbr(bounds%begc:bounds%endc,1:nlevlak)             ! "r" vector for tridiagonal solution
      real(r8) :: vpx(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)             ! phytoplankton settling velocity (m/s)
      real(r8) :: sbx(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)             ! source for whole column (gC/m^3/s)
      real(r8) :: biomas_old(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)      ! phytoplankton biomass at the last time step (gC/m3)
      real(r8) :: biomas_new(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)      ! updated phytoplankton biomass (gC/m3)
      real(r8) :: soilc_old(bounds%begc:bounds%endc,1:nlevsoi,1:nsoilclak)       ! lake sediment C at the last time step (gC/m3)
      real(r8) :: gpp_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)          ! vertically-resolved primary production (gC/m3/s)
      real(r8) :: bloss_vr(bounds%begc:bounds%endc,1:nlevlak,1:nphytolak)        ! vertically-resolved phytoplankton loss (gC/m3/s)
      real(r8) :: soilc_loss(bounds%begc:bounds%endc,1:nlevsoi,1:nsoilclak)      ! vertically-resolved soil C loss (gC/m3/s)
      !----------------------------------------------------------------------- 
      
      associate(                                                     &
            dz_lake              => col_pp%dz_lake                      , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)          
            z_lake               => col_pp%z_lake                       , & ! Input: [real(r8) (:,:)] layer depth for lake (m)
            dz                   => col_pp%dz                           , & ! Input: [real(r8) (:,:)] layer thickness (m)                   
            z                    => col_pp%z                            , & ! Input: [real(r8) (:,:)] layer depth (m)                       
            lakedepth            => col_pp%lakedepth                    , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake               => col_es%t_lake                       , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 
            t_soisno             => col_es%t_soisno                     , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            forc_pbot            => top_as%pbot                         , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)
            forc_wind            => top_as%windbot                      , & ! Input: [real(r8) (:)] horizontal component of wind at atmospheric forcing height (m/s)
            forc_po2             => top_as%po2bot                       , & ! Input: [real(r8) (:)] O2 partial pressure (Pa)
            forc_pco2            => top_as%pco2bot                      , & ! Input: [real(r8) (:)] CO2 partial pressure (Pa)      
            forc_pch4            => top_as%pch4bot                      , & ! Input: [real(r8) (:)] CH4 partial pressure (Pa)

            watsat               => soilstate_vars%watsat_col           , & ! Input: [real(r8) (:,:)] volumetric soil water at saturation (porosity)

            kme                  => lakestate_vars%lake_kme_col         , & ! Input: [real(r8) (:,:)] lake water and sediment heat diffusivity [m2/s]
            hmix                 => lakestate_vars%lake_hmix_col        , & ! Input: [real(r8) (:)] lake mixing layer depth [m]
            icethick             => lakestate_vars%lake_icethick_col    , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac              => lakestate_vars%lake_icefrac_col     , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen
            fsds_vis             => lakestate_vars%lake_fsds_vis_col    , & ! Input: [real(r8) (:,:)] incident vis radiation (W/m^2)

            lake_tp              => lakebgc_vars%tp_col                 , & ! Input: [real(r8) (:)] lake column mean total phosphorus [gP/m3]
            lake_ph              => lakebgc_vars%ph_col                 , & ! Input: [real(r8) (:)] lake column mean water pH
            lake_sdep            => lakebgc_vars%sdep_col               , & ! Input: [real(r8) (:)] lake column mean sediment deposition [g/m2/yr]

            nem_lake_grc         => lnd2atm_vars%nem_lake_grc           , & ! Output: [real(r8) (:)] gridcell average net methane correction to CO2 flux (g C/m^2/s)
          
            conc_gas_wat         => lakebgc_vars%conc_wat_col           , & ! Output: [real(r8) (:,:,:)] lake dissolved gas [mol/m3]
            conc_gas_sed         => lakebgc_vars%conc_sed_col           , & ! Output: [real(r8) (:,:,:)] sediment dissolved gas [mol/m3]
            conc_bubl            => lakebgc_vars%conc_bubl_col          , & ! Output: [real(r8) (:,:,:)] lake bubble gas [mol/m3]
            conc_iceb            => lakebgc_vars%conc_iceb_col          , & ! Output: [real(r8) (:,:)] lake bubble trapped in ice [mol/m2]
            biomas_phyto         => lakebgc_vars%biomas_phyto_col       , & ! Output: [real(r8) (:,:,:)] lake phytoplankton biomass [gC/m3]
            lake_chla            => lakebgc_vars%chla_col               , & ! Output: [real(r8) (:,:)] lake chlorophyll-a [g/m3]
            ch4_sed_diff         => lakebgc_vars%ch4_sed_diff_col       , & ! Output: [real(r8) (:)] lake ch4 diffusion at water-sediment interface [mol/m2/s]
            ch4_surf_diff        => lakebgc_vars%ch4_surf_diff_col      , & ! Output: [real(r8) (:)] lake ch4 diffusion at air-water interface [mol/m2/s]
            ch4_sed_ebul         => lakebgc_vars%ch4_sed_ebul_col       , & ! Output: [real(r8) (:)] lake ch4 ebullition at water-sediment interface [mol/m2/s]
            ch4_surf_ebul        => lakebgc_vars%ch4_surf_ebul_col      , & ! Output: [real(r8) (:)] lake ch4 ebullition at air-water interface [mol/m2/s]
            ch4_surf_flux_tot    => lakebgc_vars%ch4_surf_flux_tot_col  , & ! Output: [real(r8) (:)] lake ch4 total flux at air-water interface [kgC/m2/s]
            ch4_prod_wat         => lakebgc_vars%ch4_prod_wat_col       , & ! Output: [real(r8) (:,:)] lake ch4 production rate [mol/m3/s]
            ch4_prod_sed         => lakebgc_vars%ch4_prod_sed_col       , & ! Output: [real(r8) (:,:)] sediment ch4 production rate [mol/m3/s]
            ch4_oxid_wat         => lakebgc_vars%ch4_oxid_wat_col       , & ! Output: [real(r8) (:,:)] lake ch4 oxidation rate [mol/m3/s]
            ch4_oxid_sed         => lakebgc_vars%ch4_oxid_sed_col       , & ! Output: [real(r8) (:,:)] sediment ch4 oxidation rate [mol/m3/s]
            ch4_prod_tot         => lakebgc_vars%ch4_prod_tot_col       , & ! Output: [real(r8) (:)] lake total CH4 production rate [gC/m2/s]
            ch4_oxid_tot         => lakebgc_vars%ch4_oxid_tot_col       , & ! Output: [real(r8) (:)] lake total CH4 oxidation rate [gC/m2/s]
            nem_col              => lakebgc_vars%nem_col                , & ! Output: [real(r8) (:)] net adjustment to atm. C flux from methane production (g C/m**2/s)
            lake_soilc           => lakebgc_vars%soilc_col              & ! Output: [real(r8) (:,:,:)] lake sediment carbon pools (gC/m3)
            )

      ! Get time step
      dtime = get_step_size() / DBLE(nsubstep)

      ! substep cycle
      do it = 1, nsubstep

      ! physical and BGC processes
      do fc = 1, num_lakec
         c = filter_lakec(fc)
         t = col_pp%topounit(c)

         ch4_surf_flux_tot(c) = 0._r8
         ch4_prod_tot(c) = 0._r8
         ch4_oxid_tot(c) = 0._r8

         nem_col(c) = 0._r8

         jtop(c) = 1    ! not include snow layers

         jwat(c) = COUNT(z_lake(c,:)+0.5*dz_lake(c,:)<=icethick(c)) + 1

         ! mixing layer bottom index
         zsum = 0._r8
         jmix(c) = 1
         do j = 1, nlevlak
            if (zsum>=hmix(c)) then
               exit
            end if
            zsum = zsum + dz_lake(c,j)
            jmix(c) = j
         end do

         ! thermocline depth
         if (icethick(c)<1e-8_r8) then 
            ! top boundary layer thickness
            Hsum1 = 0._r8
            do j = 1, nlevlak, 1
               if (abs(t_lake(c,j)-t_lake(c,1))>1._r8) then
                  exit
               end if
               Hsum1 = Hsum1 + dz_lake(c,j)
            end do
            ! bottom boundary layer thickness
            Hsum2 = 0._r8
            do j = nlevlak, 1, -1
               if (abs(t_lake(c,j)-t_lake(c,nlevlak))>1._r8) then
                  exit
               end if
               Hsum2 = Hsum2 + dz_lake(c,j)
            end do
            if (Hsum1+Hsum2>lakedepth(c)) then
               Hcline = lakedepth(c)
            else
               Hcline = 0.5_r8 * (Hsum1 + (lakedepth(c)-Hsum2))
            end if
         else
            Hcline = lakedepth(c)
         end if

         ! record gas conc at the last time step
         do j = 1, nlevlak+nlevsoi
            if (j<=nlevlak) then
               conc_old(c,j,:) = conc_gas_wat(c,j,:)
            else
               conc_old(c,j,:) = conc_gas_sed(c,j-nlevlak,:)
            end if
         end do
         conc_iceb_old(c,:) = conc_iceb(c,:) 

         ! record phytoplankton biomass at the last time step
         do j = 1, nlevlak
            biomas_old(c,j,:) = biomas_phyto(c,j,:)
         end do

         do j = 1, nlevsoi
            soilc_old(c,j,:) = lake_soilc(c,j,:)
         end do 

         ! gas tranfer velocity
         temp = t_lake(c,1)
         do k = 1, ngaslak
            if (icethick(c)<1e-8_r8) then
               kg(c,k) = GasFluxVelocity(forc_wind(t), temp, k)
            else
               kg(c,k) = 0._r8
            end if
            if (k==gn2lak) then
               Pgas = Xn2 * forc_pbot(t)
            else if (k==go2lak) then
               Pgas = forc_po2(t)
            else if (k==gco2lak) then
               Pgas = forc_pco2(t)
            else if (k==gch4lak) then
               Pgas = forc_pch4(t)
            else
               call endrun(msg=' ERROR: no gas type index is found.'//&
                  errMsg(__FILE__, __LINE__))
            end if
            conc_eq(c,k) = GasEQConc(Pgas, temp, lake_ph(c), k) 
            df_top(c,k) = kg(c,k) * (conc_gas_wat(c,1,k) - conc_eq(c,k))
         end do

         ! phytoplankton settling velocity
         do j = 1, nlevlak
            if (t_lake(c,j)>=tfrz) then
               if (z_lake(c,j)<Hcline) then
                  vpx(c,j,:) = 0._r8
               else
                  vpx(c,j,small_phyto) = Vset_small
                  vpx(c,j,large_phyto) = Vset_large
               end if 
            else
               vpx(c,j,:) = 0._r8
            end if
         end do 

         ! initialize source and sink terms of dissolved gases
         ex(c,:,:) = 0._r8
         gsrc(c,:,:) = 0._r8
         gsnk(c,:,:) = 0._r8
         soilc_loss(c,:,:) = 0._r8

         ! sediment bubbling
         call SedimentEbullition(c, soilstate_vars, lakebgc_vars, eb_sed(c,:))

         ! bubble transport in water (11 size bins)
         call BubbleDynamics(c, lakestate_vars, lakebgc_vars, eb_sed(c,:), &
               ex(c,:,:), dtime)

         ! photosynthesis
         call Photosynthesis(c, lakestate_vars, lakebgc_vars, gpp_vr(c,:,:), &
               gsrc(c,:,:), gsnk(c,:,:))

         ! autotrophic respiration
         call AutotrophicR(c, lakestate_vars, lakebgc_vars, bloss_vr(c,:,:), &
               bphyto_dead(c,:), gsrc(c,:,:), gsnk(c,:,:))

         ! heterotrophic respiration
         call HeterotrophicR(c, lakestate_vars, lakebgc_vars, soilc_loss(c,:,:), &
               gsrc(c,:,:), gsnk(c,:,:))

         ! methanogenesis 
         call Methanogenesis(c, lakestate_vars, lakebgc_vars, soilc_loss(c,:,:), &
               gsrc(c,:,:), gsnk(c,:,:)) 

         ! methanotrophy
         call Methanotrophy(c, lakestate_vars, lakebgc_vars, gsrc(c,:,:), gsnk(c,:,:))

         ! phytoplankton deposition
         bphyto_dep_tot(c) = 0._r8
         do k = 1, nphytolak
            bphyto_dep(c,k) = vpx(c,nlevlak,k)*biomas_old(c,nlevlak,k)
            bphyto_dep_tot(c) = bphyto_dep_tot(c) + bphyto_dep(c,k)*(1._r8-LakeBGCParamsInst%frcResusp)
            bphyto_ssp(c,k) = bphyto_dep(c,k)*LakeBGCParamsInst%frcResusp/lakedepth(c)
         end do

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
               qx(c,j) = watsat(c,j-nlevlak)
            end if

            do k = 1, ngaslak
               if (j==1) then
                  sx(c,j,k) = ex(c,j,k) + gsrc(c,j,k) - gsnk(c,j,k) - df_top(c,k)/dz_lake(c,j)
               else if (j<=nlevlak) then
                  sx(c,j,k) = ex(c,j,k) + gsrc(c,j,k) - gsnk(c,j,k)
               else
                  sx(c,j,k) = gsrc(c,j,k) - gsnk(c,j,k)
               end if
            end do
         end do
      end do

      do j = 1, nlevlak+nlevsoi
         do fc = 1, num_lakec
            c = filter_lakec(fc)

            if (j<nlevlak) then
               dx(c,j) = ( kme(c,j)*kme(c,j+1) * (dz_lake(c,j+1)+dz_lake(c,j)) ) &
                           / ( kme(c,j)*dz_lake(c,j+1) + kme(c,j+1)*dz_lake(c,j) )
            else if (j==nlevlak) then
               dzp = zx(c,j+1) - zx(c,j)
               dx(c,j) = ( kme(c,j+1)*kme(c,j)*dzp / &
                    (kme(c,j+1)*dz_lake(c,j)/2._r8 + kme(c,j)*z(c,1) ) )
            else
               dx(c,j) = kme(c,j)
            end if
         end do
      end do

      ! solve governing equations 
      do k = 1, ngaslak
         
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
                  vb(c,j)  = 1._r8 + (1._r8-cnfac)/qx(c,j)*dtime*(dx(c,j)/dzp+dx(c,j-1)/dzm)/dzx(c,j)
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

         ! mixing
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            j0 = jwat(c) 
            cgas_avg = sum(conc_new(c,j0:jmix(c),k)*dz_lake(c,j0:jmix(c))) / &
               sum(dz_lake(c,j0:jmix(c)))
            conc_new(c,j0:jmix(c),k) = cgas_avg
         end do
      end do

      ! calculate CH4 flux at the lake surface
      do fc = 1, num_lakec
         c = filter_lakec(fc)

         ch4_surf_diff(c) = df_top(c,gch4lak) 

         do k = 1, ngaslak
            do j = 1, nlevlak+nlevsoi
               ! correct negative concentration
               if (conc_new(c,j,k)<0._r8) then
                  if (k==gch4lak) then ! for CH4
                     deficit = -conc_new(c,j,k) * dzx(c,j)  ! mol/m2
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
                        write(iulog,*) 'Negative conc. in LakeBGCDynamics g,c,j,deficit (mol):',g,c,j,deficit
                     end if
                     conc_new(c,j,k) = 0._r8
                     ch4_surf_diff(c) = ch4_surf_diff(c) - deficit/dtime
                  else  ! for other gases
                     conc_new(c,j,k) = 0._r8
                  end if
               end if
               
               ! set dissolved gas conc for outputs
               if (j<=nlevlak) then
                  conc_gas_wat(c,j,k) = conc_new(c,j,k)
               else
                  conc_gas_sed(c,j-nlevlak,k) = conc_new(c,j,k)
               end if
            end do
         end do
      end do

      do fc = 1, num_lakec
         c = filter_lakec(fc)

         ch4_surf_flux_tot(c) = 1e-3_r8 * catomw * (ch4_surf_diff(c) + ch4_surf_ebul(c))

         do j = 1, nlevlak
            ch4_prod_tot(c) = ch4_prod_tot(c) + ch4_prod_wat(c,j)*dz_lake(c,j)*catomw
            ch4_oxid_tot(c) = ch4_oxid_tot(c) + ch4_oxid_wat(c,j)*dz_lake(c,j)*catomw 
         end do

         do j = 1, nlevsoi
            ch4_prod_tot(c) = ch4_prod_tot(c) + ch4_prod_sed(c,j)*dz(c,j)*catomw
            ch4_oxid_tot(c) = ch4_oxid_tot(c) + ch4_oxid_sed(c,j)*dz(c,j)*catomw 
            if (j==nlevsoi) then
               ! Adjustment to NEE flux to atm. for methane production
               nem_col(c) = nem_col(c) - ch4_prod_tot(c)
               ! Adjustment to NEE flux to atm. for methane oxidation
               nem_col(c) = nem_col(c) + ch4_oxid_tot(c)
            end if
         end do
      end do

      ! Do Balance Check
      do j = 1, nlevlak+nlevsoi
         do fc = 1, num_lakec
            c = filter_lakec(fc)
            
            if (j==1) errch4(c) = conc_iceb(c,gch4lak) - conc_iceb_old(c,gch4lak)
            if (j<=nlevlak) then
               errch4(c) = errch4(c) + dzx(c,j)*(conc_new(c,j,gch4lak)-conc_old(c,j,gch4lak))
               errch4(c) = errch4(c) - ch4_prod_wat(c,j)*dzx(c,j)*dtime
               errch4(c) = errch4(c) + ch4_oxid_wat(c,j)*dzx(c,j)*dtime
            else
               errch4(c) = errch4(c) + dzx(c,j)*(conc_new(c,j,gch4lak)-conc_old(c,j,gch4lak))
               errch4(c) = errch4(c) - ch4_prod_sed(c,j-nlevlak)*dzx(c,j)*dtime
               errch4(c) = errch4(c) + ch4_oxid_sed(c,j-nlevlak)*dzx(c,j)*dtime
            end if
         end do
      end do
      
      do fc = 1, num_lakec
         c = filter_lakec(fc)
      
         errch4(c) = errch4(c) + (ch4_surf_diff(c) + ch4_surf_ebul(c))*dtime
         if (abs(errch4(c)) < 1.e-8_r8) then
            ch4_surf_diff(c) = ch4_surf_diff(c) - errch4(c)/dtime
         else ! errch4 > 1e-8 mol / m^2 / timestep
            write(iulog,*)'CH4 Conservation Error in LakeBGCDynamics during c, errch4 (mol /m^2.timestep)', &
                 c,errch4(c)
            g = col_pp%gridcell(c)
            write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(g),grc_pp%londeg(g)
            call endrun(msg=' ERROR: CH4 Conservation Error in LakeBGCDynamics during diffusion'//&
                 errMsg(__FILE__, __LINE__))
         end if
      end do

      ! Now average up to gridcell for fluxes
      call c2g( bounds, &
           nem_col(bounds%begc:bounds%endc), nem_lake_grc(bounds%begg:bounds%endg), &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )

      ! solve phytoplankton equations
      do k = 1, nphytolak

         ! Set up vector r and vectors a, b, c that define tridiagonal matrix 
         do j = 1, nlevlak
            do fc = 1, num_lakec
               c = filter_lakec(fc)
               
               ! source/sink term
               sbx(c,j,k) = gpp_vr(c,j,k) - bloss_vr(c,j,k) + bphyto_ssp(c,k)
               if (j==nlevlak) then
                  sbx(c,j,k) = sbx(c,j,k) - bphyto_dep(c,k)/dzx(c,j)
               end if

               if (j==1) then ! top layer
                  vzp       = 0.5_r8*(vpx(c,j,k)+vpx(c,j+1,k)) 
                  vba(c,j)  = 0._r8
                  vbb(c,j)  = 1._r8 + 0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzp
                  vbc(c,j)  = 0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzp
                  vbr(c,j)  = biomas_old(c,j,k) + dtime*sbx(c,j,k) - 0.5_r8*dtime/dzx(c,j)* &
                     cnfac*vzp*(biomas_old(c,j,k)+biomas_old(c,j+1,k))
               else if (j<nlevlak) then
                  vzm       = 0.5_r8*(vpx(c,j,k)+vpx(c,j-1,k))
                  vzp       = 0.5_r8*(vpx(c,j,k)+vpx(c,j+1,k))
                  vba(c,j)  = -0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzm
                  vbb(c,j)  = 1._r8 + 0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzp - &
                     0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzm
                  vbc(c,j)  = 0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzp
                  vbr(c,j)  = biomas_old(c,j,k) + dtime*sbx(c,j,k) - 0.5_r8*dtime/dzx(c,j)* &
                     cnfac*vzp*(biomas_old(c,j,k)+biomas_old(c,j+1,k)) + 0.5_r8*dtime/dzx(c,j)* &
                     cnfac*vzm*(biomas_old(c,j,k)+biomas_old(c,j-1,k))
               else  ! bottom layer
                  vzm       = 0.5_r8*(vpx(c,j,k)+vpx(c,j-1,k))
                  vba(c,j)  = -0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzm
                  vbb(c,j)  = 1._r8 - 0.5_r8*dtime/dzx(c,j)*(1._r8-cnfac)*vzm
                  vbc(c,j)  = 0._r8
                  vbr(c,j)  = biomas_old(c,j,k) + dtime*sbx(c,j,k) + 0.5_r8*dtime/dzx(c,j)* &
                     cnfac*vzm*(biomas_old(c,j,k)+biomas_old(c,j-1,k))
               end if
            end do
         end do

         call Tridiagonal(bounds, 1, nlevlak, &
               jtop(bounds%begc:bounds%endc), &
               num_lakec, filter_lakec, &
               vba(bounds%begc:bounds%endc, :), &
               vbb(bounds%begc:bounds%endc, :), &
               vbc(bounds%begc:bounds%endc, :), &
               vbr(bounds%begc:bounds%endc, :), &
               biomas_new(bounds%begc:bounds%endc, :, k))
      end do

      ! check negative values and balance
      do fc = 1, num_lakec
         c = filter_lakec(fc)

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
                     write(iulog,*) 'Negative biomas. in LakeBGCDynamics c,j,deficit (mol):',c,j,deficit
                  end if
                  biomas_new(c,j,k) = 0._r8
                  bphyto_dep_tot(c) = bphyto_dep_tot(c) - deficit/dtime
               end if

               ! set phytoplankton biomass for outputs
               biomas_phyto(c,j,k) = biomas_new(c,j,k) 
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

         ! decomposition
         do k = 1, nsoilclak
            do j = 1, nlevsoi
               lake_soilc(c,j,k) = max(0._r8, soilc_old(c,j,k)-soilc_loss(c,j,k)*dtime)
            end do
         end do

         ! POC deposition (Tan et al., 2017)
         cdep = 0.061_r8 * lake_sdep(c)**0.766_r8 / 3.1536e7_r8

         ! deposition thickness (m)
         zsoi_inc = lake_sdep(c)/3.1536e7_r8 + max(bphyto_dep_tot(c),0._r8)
         do k = 1, nphytolak
            zsoi_inc = zsoi_inc + bphyto_dead(c,k)
         end do
         zsoi_inc = 1e-3_r8 * zsoi_inc / 1.4e3_r8 * dtime 

         ! deposition and compaction
         do j = 1, nlevsoi
            if (j==1) then
               lake_soilc(c,j,actC) = lake_soilc(c,j,actC)*(1._r8-zsoi_inc/dz(c,j)) + &
                  0.05_r8*cdep/dz(c,j)*dtime
               lake_soilc(c,j,pasC) = lake_soilc(c,j,pasC)*(1._r8-zsoi_inc/dz(c,j)) + &
                  0.95_r8*cdep/dz(c,j)*dtime
               do k = 1, nphytolak
                  lake_soilc(c,j,actC) = lake_soilc(c,j,actC) + bphyto_dead(c,k)/dz(c,j)*dtime
               end do 
               lake_soilc(c,j,actC) = lake_soilc(c,j,actC) + bphyto_dep_tot(c)/dz(c,j)*dtime
            else
               lake_soilc(c,j,:) = lake_soilc(c,j,:)*(1._r8-zsoi_inc/dz(c,j)) + &
                  lake_soilc(c,j-1,:)*zsoi_inc/dz(c,j)
            end if
            lake_soilc(c,j,actC) = max(0._r8, lake_soilc(c,j,actC))
            lake_soilc(c,j,pasC) = max(0._r8, lake_soilc(c,j,pasC))
         end do

      end do

      end do ! end substep cycling

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
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak 
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      real(r8), intent(in) :: pH       ! pH
      integer,  intent(in) :: gas      ! gas index
      !
      ! !LOCAL VARIABLES:
      real(r8), parameter :: SOLO2(11) = (/14.6, 12.8, 11.3, 10.1, &
                              9.1, 8.3, 7.6, 7.0, 6.5, 6.0, 5.6/)
      real(r8) :: hi, kc1, kc2, par
      integer  :: indx 
      !--------------------------------------------------------------------

      if (gas==gn2lak) then
         HenrySolubility = 6.1e-6_r8*exp(-1300.0*(1.0/temp-1.0/298.0)) 
      else if (gas==go2lak) then
         if (temp>=tfrz .and. temp<=tfrz+50.0) then 
            indx = min( int((temp-tfrz)/5) + 1, 10 )
            par = (temp - tfrz - 5*indx + 5) / 5.0_r8
            HenrySolubility = ( SOLO2(indx+1) * par + &
               SOLO2(indx) * (1-par) ) / 6.596e5_r8 
         else
            HenrySolubility = 1.3e-5_r8*exp(-1500.0*(1/temp-1/298.0))
         end if
      else if (gas==gco2lak) then
         hi = 10.0_r8**(-pH)    ! Concentration of hydrogen ion
         ! rate constant of dissolved CO2 for first and second dissolution
         kc1 = 4.3e-7_r8*exp(-921.4*(1/temp-1/298.0))
         kc2 = 4.7e-11_r8*exp(-1787.4*(1/temp-1/298.0))
         HenrySolubility = 3.4e-4_r8*exp(-2400.0*(1/temp-1/298.0))* &
            (1.0 + kc1/hi + kc1*kc2/hi**2.0)
      else if (gas==gch4lak) then
         HenrySolubility = 1.3e-5_r8*exp(-1700.0*(1/temp-1/298.0))
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
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: temp     ! temperature (K)
      integer,  intent(in) :: gas      ! gas index
      !
      ! !LOCAL VARIABLES:
      !--------------------------------------------------------------------

      if (gas==gn2lak) then
         BubbleGasDiffusivity = 2.57e-9_r8 * (temp/273.0)
      else if (gas==go2lak) then
         BubbleGasDiffusivity = 2.4e-9_r8 * (temp/298.0)
      else if (gas==gco2lak) then
         BubbleGasDiffusivity = 1.81e-6_r8 * exp(-2032.6/temp)
      else if (gas==gch4lak) then
         BubbleGasDiffusivity = 1.5e-9_r8 * (temp/298.0)
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
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak
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
      if (gas==gn2lak) then
         SchmidtNumber = 1970.7_r8 - 131.45_r8*tw + 4.139_r8*(tw**2.0) - &
            0.052106_r8*(tw**3.0)
      else if (gas==go2lak) then
         SchmidtNumber = 1800.6_r8 - 120.1_r8*tw + 3.7818_r8*(tw**2.0) - &
            0.047608_r8*(tw**3.0)
      else if (gas==gco2lak) then
         SchmidtNumber = 1911.0_r8 - 113.7_r8*tw + 2.967_r8*(tw**2.0) - &
            0.02943_r8*(tw**3.0)
      else if (gas==gch4lak) then
         SchmidtNumber = 1898_r8 - 110.1_r8*tw + 2.834_r8*(tw**2.0) - &
            0.02791_r8*(tw**3.0)
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
         KinematicViscosity = 1e30_r8 
      else
         KinematicViscosity = 1._r8/denh2o * ( 2.414e-5_r8 * &
            10.0**(247.8/(temp-140.0)) )
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
      real(r8) :: rr, xx, tmp1, tmp2, par
      !--------------------------------------------------------------------

      if (vsc>1.0e10_r8) then
         BuoyantVelocity = 1e-30_r8
         return
      end if
      
      rr = 1e6_r8 * radius      ! change meter to micron
      if (rr<80._r8) then
         BuoyantVelocity = grav * (radius**2.0) / (3.0*vsc)
      else if (rr>150._r8) then
         xx = grav * (radius**3.0) / (vsc**2.0)
         BuoyantVelocity = grav * (radius**2.0) / &
            (vsc*18._r8*(1.0-2.0/(1.0+sqrt(1.0+0.091*xx))))
      else
         xx = grav * (radius**3.0) / (vsc**2.0)
         tmp1 = grav * (radius**2.0) / (3.0*vsc)
         tmp2 = grav * (radius**2.0) / &
            (vsc*18.0_r8*(1.0-2.0/(1.0+sqrt(1+0.091*xx))))
         par = (rr-80._r8) / 70._r8
         BuoyantVelocity = tmp1*(1-par) + tmp2*par
      end if

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
   real(r8) function GasEQConc(pressure, temp, pH, gas)
      !
      ! !DESCRIPTION:
      ! Calculate equilibrium gas concentration (mol/m^3) 
      !
      ! !USES:
      use elm_varcon         , only : grav
      ! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: pressure ! partial pressure (Pa)
      real(r8), intent(in) :: temp     ! temperature (K)
      real(r8), intent(in) :: pH       ! pH
      integer,  intent(in) :: gas      ! gas id
      !
      ! !LOCAL VARIABLES:
      real(r8) :: henry
      !--------------------------------------------------------------------

      henry = HenrySolubility(temp, pH, gas)
      GasEQConc = henry * pressure     
 
   end function GasEQConc

   subroutine SedimentEbullition(c, soilstate_vars, lakebgc_vars, eb_sed) 
      !
      ! !DESCRIPTION:
      ! Simulate bubble generation 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)  :: c                         ! column index
      type(soilstate_type)   , intent(in)  :: soilstate_vars
      type(lakebgc_type)     , intent(in)  :: lakebgc_vars
      real(r8)               , intent(out) :: eb_sed(1:ngaslak)         ! ebullition rate (mol/m3/s)
      !
      ! !CONSTANTS
      real(r8), parameter :: Ae = 0.4                 ! relative concentration at which ebullition begins
      !
      ! !LOCAL VARIABLES:
      real(r8) :: temp, por, depth
      real(r8) :: c_n2, c_ch4
      real(r8) :: Hn2, Hch4
      real(r8) :: Pn2, Pch4
      real(r8) :: Prn2, Prch4, Prtot, Prnet, Psat 
      integer  :: j, t, k
      real(r8) :: eb(1:nlevsoi,1:ngaslak)
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz             => col_pp%dz                        , & ! Input: [real(r8) (:,:)] layer thickness (m)                   
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_soisno       => col_es%t_soisno                  , & ! Input: [real(r8) (:,:)] soil (or snow) temperature (Kelvin)

            forc_pbot      => top_as%pbot                      , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)

            watsat         => soilstate_vars%watsat_col        , & ! Input: [real(r8) (:,:)] volumetric soil water at saturation (porosity)

            lake_ph        => lakebgc_vars%ph_col              , & ! Input: [real(r8) (:)] lake column mean water pH      
            conc_gas_sed   => lakebgc_vars%conc_sed_col        , & ! Input: [real(r8) (:,:,:)] sediment dissolved gas [mol/m3]

            ch4_sed_ebul   => lakebgc_vars%ch4_sed_ebul_col    & ! Output:  [real(r8) (:)] sediment CH4 ebulltion (mol/m2/s)
            )

      t = col_pp%topounit(c)
      
      do j = nlevlak+1, nlevlak+nlevsoi
         if (j-nlevlak>nlevsoi) then
            eb(j-nlevlak,1:ngaslak) = 0.0_r8
            cycle
         end if
         temp = t_soisno(c,j-nlevlak)
         por = watsat(c,j-nlevlak)
         depth = lakedepth(c) + 0.5_r8 * dz(c,j-nlevlak)
         c_n2 = conc_gas_sed(c,j-nlevlak,gn2lak)
         c_ch4 = conc_gas_sed(c,j-nlevlak,gch4lak)
         Hn2 = HenrySolubility(temp, lake_ph(c), gn2lak)
         Hch4 = HenrySolubility(temp, lake_ph(c), gch4lak)
         Prn2 = c_n2/Hn2                  ! partial pressure of N2
         Prch4 = c_ch4/Hch4               ! partial pressure of CH4
         Prtot = Pn2 + Pch4
         Psat = Ae * por * (forc_pbot(t) + grav*denh2o*depth)
         Prnet = Prtot - Psat             ! excessive pressure
         if (Prnet>1e-8_r8) then
            eb(j-nlevlak,gn2lak) = LakeBGCParamsInst%Re * &
                  (Prn2/Prtot) * Prnet * Hn2
            eb(j-nlevlak,gch4lak) = LakeBGCParamsInst%Re * &
                  (Prch4/Prtot) * Prnet * Hch4
            eb(j-nlevlak,go2lak) = 0.0_r8
            eb(j-nlevlak,gco2lak) = 0.0_r8
         else
            eb(j-nlevlak,1:ngaslak) = 0.0_r8
         end if
      end do

      ! ebullition at the water-sediment interface
      do k = 1, ngaslak
         eb_sed(k) = sum(eb(:,k)*dz(c,1:nlevsoi))
      end do
      ch4_sed_ebul(c) = eb_sed(gch4lak)

      end associate
   end subroutine SedimentEbullition

   !-----------------------------------------------------------------------
   subroutine BubbleDynamics(c, lakestate_vars, lakebgc_vars, ebb, ex, dtime) 
      !
      ! !DESCRIPTION:
      ! Simulate bubble transport in the water column 
      !
      ! !USES:
      use elm_varcon         , only : tfrz, denh2o, grav, rgas, rpi
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak
      use LakeBGCType        , only : gn2lak, go2lak, gco2lak, gch4lak
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(in)    :: ebb(1:ngaslak)          ! bottom ebullition (mol/m2/s)
      real(r8)               , intent(out)   :: ex(1:nlevlak,1:ngaslak) ! gas dissolution (mol/m3/s)
      real(r8)               , intent(in)    :: dtime                   ! timestep size [seconds]
      !
      ! !CONSTANTS:
      integer,  parameter :: nrlev = 11   ! bubble size number
      real(r8), parameter :: Rb0(11) = (/2.5, 3.25, 4.0, 4.75, 5.5, &
                              6.25, 7.0, 7.75, 8.5, 9.25, 10.0/)  ! mm
      !
      ! !LOCAL VARIABLES:
      integer  :: rindx, locindx
      integer  :: j, k, t, jtop
      real(r8) :: pos, pr, vb, temp
      real(r8) :: rr, rmin, rmax
      real(r8) :: drr, dr_tmp1, dr_tmp2
      real(r8) :: dt, tcur, tini
      real(r8) :: numb, pr_tmp
      real(r8) :: vsc(1:nlevlak)                                  ! viscosity
      real(r8) :: gama(1:nlevlak)                                 ! bubble tension
      real(r8) :: peclet(1:ngaslak)                               ! Peclet number
      real(r8) :: reynold(1:ngaslak)                              ! Reynold number
      real(r8) :: nusselt(1:ngaslak)                              ! Nusselt number
      real(r8) :: ndm_tmp1(1:ngaslak)
      real(r8) :: ndm_tmp2(1:ngaslak)
      real(r8) :: ndm(1:ngaslak)                                  ! gas exchange rate of single bubble
      real(r8) :: vbg(1:ngaslak)                                  ! bubble gas conc (mol/m^3/m)
      real(r8) :: vbg_init(1:ngaslak)
      real(r8) :: vrr(1:nrlev)                                    ! bubble radius at lake surface
      real(r8) :: schmidt(1:nlevlak,1:ngaslak)                    ! Schmidt number
      real(r8) :: gdiff(1:nlevlak,1:ngaslak)                      ! gas diffusivity
      real(r8) :: gsolu(1:nlevlak,1:ngaslak)                      ! gas solubility
      real(r8) :: ex_rr(1:nlevlak,1:ngaslak,1:nrlev)              ! size-resolved gas exchange (mol/m^3/s/m)
      real(r8) :: bgas_rr(1:nlevlak,1:ngaslak,1:nrlev)            ! size-resolved bubble gas conc (mol/m^3/m)
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 

            forc_pbot      => top_as%pbot                      , & ! Input: [real(r8) (:)] atmospheric pressure (Pa)

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            lake_ph        => lakebgc_vars%ph_col              , & ! Input: [real(r8) (:)] lake column mean water pH      
            conc_gas_wat   => lakebgc_vars%conc_wat_col        , & ! Input: [real(r8) (:,:,:)] lake dissolved gas [mol/m3]

            ch4_surf_ebul  => lakebgc_vars%ch4_surf_ebul_col   , & ! Output: [real(r8) (:)] lake ch4 ebullition at air-water interface [mol/m2/s] 
            conc_bubl      => lakebgc_vars%conc_bubl_col       , & ! Output: [real(r8) (:,:,:)] lake bubble gas [mol/m3]
            conc_iceb      => lakebgc_vars%conc_iceb_col       & ! Output: [real(r8) (:,:)] lake bubble trapped in ice [mol/m2]
            )

      t = col_pp%topounit(c)

      ! top water layer index
      jtop = COUNT(z_lake(c,:)+0.5*dz_lake(c,:)<=icethick(c)) + 1  

      ! Calculate bubble release rate at the water-sediment interface
      do j = 1, nlevlak
         temp = t_lake(c,nlevlak)
         gama(j) = SurfaceTension(temp) 
         vsc(j) = KinematicViscosity(temp)
         do k = 1, ngaslak
            gsolu(j,k) = HenrySolubility(temp, lake_ph(c), k)
            gdiff(j,k) = BubbleGasDiffusivity(temp, k)
            schmidt(j,k) = SchmidtNumber(temp, k)
         end do
      end do
      pr = forc_pbot(t) + denh2o*grav*lakedepth(c)
      rmin = 1e-3_r8 * Rb0(1)
      rmax = 1e-3_r8 * Rb0(nrlev)
      do rindx = 1, nrlev
         rr = Rb0(rindx)
         vb = BuoyantVelocity(rr, vsc(nlevlak)) 
         bgas_rr(nlevlak,:,rindx) = ebb * (pr*rr+gama(nlevlak)) / vb / &
            (0.5*(rmax**2.0-rmin**2.0)*pr+gama(nlevlak)*(rmax-rmin))
      end do

      ! time step (s)
      dt = 0.4_r8 / (0.1418_r8*(10.0*Rb0(1))**2.0 + &
         0.05579_r8*(10.0*Rb0(1)) + 0.7794_r8)

      do rindx = 1, nrlev
         ! no sediment bubbling
         if (sum(ebb)<1e-8_r8) then
            bgas_rr(:,:,rindx) = 0._r8
            ex_rr(:,:,rindx) = 0._r8
            cycle
         end if
         locindx = nlevlak
         pos = -lakedepth(c) 
         rr = 1e-3_r8 * Rb0(rindx)
         vbg = bgas_rr(nlevlak,:,rindx)
         vbg_init = vbg
         tini = 0._r8
         tcur = 0._r8
         do while (pos<-icethick(c) .and. rr>1e-8_r8 .and. &
               sum(vbg)>1e-8_r8)
            if (pos>=-z_lake(c,locindx)+0.5_r8*dz_lake(c,locindx)) then
               bgas_rr(locindx,:,rindx) = 0.5_r8 * (vbg + vbg_init) 
               ex_rr(locindx,:,rindx) = (vbg_init - vbg) / (tcur - tini)
               vbg_init = vbg
               tini = tcur
               locindx = locindx - 1 
            end if
            temp = t_lake(c,locindx)
            vb = BuoyantVelocity(rr, vsc(locindx))
            pr = forc_pbot(t) - denh2o*grav*pos + 2*gama(locindx)/rr
            ! bubble number
            numb = 0.75e-3_r8*rgas*temp/(rpi*pr*rr**3.0)*sum(vbg)
            ! gas exchange rate
            peclet = rr * vb / gdiff(locindx,:)
            reynold = peclet / schmidt(locindx,:) 
            where (reynold<=1._r8)
               nusselt = sqrt(2._r8*rpi*peclet/3.0_r8)
            elsewhere
               nusselt = 0.45_r8*(reynold**(1.0/6.0))*(peclet**(1.0/3.0)) 
            end where
            ndm_tmp1 = -4._r8*rpi*rr*gdiff(locindx,:)*nusselt*gsolu(locindx,:)*pr 
            ndm_tmp2 = 4._r8*rpi*rr*gdiff(locindx,:)*nusselt*conc_gas_wat(c,locindx,:)
            ndm = ndm_tmp1*vbg/sum(vbg) + ndm_tmp2
            ! new bubble gas
            vbg = vbg + numb * ndm * dt
            where (vbg<0._r8) vbg = 0._r8    ! remove small negatives
            pr_tmp = 3._r8*forc_pbot(t) - 3._r8*denh2o*grav*pos + &
               4._r8*gama(locindx)/rr
            dr_tmp1 = 0.75e-3_r8*rgas*temp/(rpi*rr**2.0)/pr_tmp 
            dr_tmp2 = rr*denh2o*grav*vb/pr_tmp
            ! new bubble radius
            rr = rr + (dr_tmp1*sum(ndm)+dr_tmp2)*dt
            tcur = tcur + dt
            pos = pos + vb * dt     ! new position
         end do 
         if (locindx>=1) then
            bgas_rr(locindx,:,rindx) = 0.5_r8 * (vbg + vbg_init)
            ex_rr(locindx,:,rindx) = (vbg_init - vbg) / (tcur - tini)
         end if
         vrr(rindx) = max(0._r8,rr) 
         ! for ice layers
         bgas_rr(1:locindx-1,:,rindx) = 0._r8
         ex_rr(1:locindx-1,:,rindx) = 0._r8   
      end do

      ! update state and flux variables 
      ch4_surf_ebul(c) = 0._r8
      conc_bubl(c,:,:) = 0._r8
      ex = 0._r8
      do rindx = 1, nrlev
         if (rindx==1) then
            drr = 0.5e-3_r8 * (Rb0(rindx+1) - Rb0(rindx))
         else if (rindx==nrlev) then
            drr = 0.5e-3_r8 * (Rb0(rindx) - Rb0(rindx-1))
         else
            drr = 0.5e-3_r8 * (Rb0(rindx+1) - Rb0(rindx-1))
         end if 
         
         ! size-integrated bubble gas and gas exchange
         do j = 1, nlevlak
            conc_bubl(c,j,:) = conc_bubl(c,j,:) + &
               bgas_rr(j,:,rindx) * drr 
            ex(j,:) = ex(j,:) + ex_rr(j,:,rindx) * drr
         end do

         ! CH4 ebullition at the lake surface
         if (icethick(c)<1.e-8_r8) then
            vb = BuoyantVelocity(vrr(rindx), vsc(nlevlak)) 
            ch4_surf_ebul(c) = ch4_surf_ebul(c) + vb * &
               bgas_rr(nlevlak,gch4lak,rindx) * drr
         end if

         ! bubbles trapped in ice layers
         if (icethick(c)>=1.e-8_r8) then
            vb = BuoyantVelocity(vrr(rindx), vsc(jtop))
            conc_iceb(c,:) = conc_iceb(c,:) + &
               vb * bgas_rr(jtop,:,rindx) * drr * dtime 
         end if
      end do 

      ! CH4 ebullition or loss from trapped ice bubbles 
      if (icethick(c)<1.e-8_r8) then
         ch4_surf_ebul(c) = ch4_surf_ebul(c) + conc_iceb(c,gch4lak) * &
            LakeBGCParamsInst%icebflux
         do k = 1, ngaslak
            conc_iceb(c,k) = conc_iceb(c,k) * max(0._r8, &
               (1._r8 - LakeBGCParamsInst%icebflux*dtime))
         end do
      else
         do k = 1, ngaslak
            ex(jtop,k) = ex(jtop,k) + conc_iceb(c,k) / &
               dz_lake(c,jtop) * LakeBGCParamsInst%icebloss
            conc_iceb(c,k) = conc_iceb(c,k) * max(0._r8, &
               (1._r8 - LakeBGCParamsInst%icebloss*dtime)) 
         end do
      end if

      end associate
   end subroutine BubbleDynamics

   !----------------------------------------------------------------------- 
   subroutine Photosynthesis(c, lakestate_vars, lakebgc_vars, gpp_vr, &
               gsrc, gsnk)
      !
      ! !DESCRIPTION:
      ! Simulate phytoplankton photosynthesis 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak, nphytolak
      use LakeBGCType        , only : small_phyto, large_phyto 
      use LakeBGCType        , only : go2lak, gco2lak
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(out)   :: gpp_vr(1:nlevlak,1:nphytolak)       ! gC/m3/s
      real(r8)               , intent(inout) :: gsrc(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
      real(r8)               , intent(inout) :: gsnk(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s   
      !
      ! !CONSTANTS
      real(r8), parameter :: ThetaG = 1.08   ! temperature multiplier for phytoplankton growth
      real(r8), parameter :: kt(nphytolak) = (/1.93841_r8, 12.76830_r8/)
      real(r8), parameter :: at(nphytolak) = (/29.27777_r8, 21.67022_r8/)
      real(r8), parameter :: bt(nphytolak) = (/0.28991_r8, 0.21632_r8/)
      real(r8), parameter :: Kco2 = 6.163e-2_r8  ! CO2 half-saturation constant (mol/m3) 
      real(r8), parameter :: C2Chlmin(nphytolak) = (/1.44e2_r8, 1.44e2_r8/)   ! minimum carbon to chlorophill ratio (gC g chl-1)
      real(r8), parameter :: C2Chlmax(nphytolak) = (/4.8e2_r8, 3.6e2_r8/)     ! maximum carbon to chlorophill ratio (gC g chl-1)
      real(r8), parameter :: mu0(nphytolak) = (/0.4_r8, 1.0_r8/)  ! the maximum growth rate of phytoplankton at 0 celcius (d-1)
      real(r8), parameter :: Kpc2chl(nphytolak) = (/95_r8, 70_r8/)   ! The slope of C:Chl ratio vs. growth rate (gC g chl-1 d)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, c2chl, chla, ipar0, ipar
      real(r8) :: fpar, ftemp, fco2, ftp
      real(r8) :: Vch, phAlpha, phBeta, Ktp
      real(r8) :: C2Chl0, Vm0
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
            fsds_vis       => lakestate_vars%lake_fsds_vis_col , & ! Input: [real(r8) (:,:)] incident vis radiation (W/m^2)

            lake_tp        => lakebgc_vars%tp_col              , & ! Input: [real(r8) (:)] epilimnion average total phosphorus [gP/m3]
            conc_gas_wat   => lakebgc_vars%conc_wat_col        , & ! Input: [real(r8) (:,:,:)] lake dissolved gas [mol/m3]
            biomas_phyto   => lakebgc_vars%biomas_phyto_col    , & ! Input: [real(r8) (:,:,:)] phytoplankton biomass (gC/m3)

            lake_gpp    => lakebgc_vars%gpp_col                , & ! Output: [real(r8) (:)] depth-integrated gross primary production (gC/m2/s)
            lake_chla   => lakebgc_vars%chla_col               & ! Output: [real(r8) (:,:)] chlorophyll-a conc (g/m3)
            )

      ipar0 = 1e-6 * 4.6 * fsds_vis(c,1)   ! surface incident vis radiation (mol/m2/s) 

      lake_gpp(c) = 0._r8
      do j = 1, nlevlak

         ! no radiation
         if (ipar0<1e-8_r8) then
            lake_chla(c,j) = sum( biomas_phyto(c,j,:)/C2Chlmax )
            gpp_vr(j,:)  = 0._r8
            cycle
         end if

         lake_chla(c,j) = 0._r8

         tw = t_lake(c,j) - tfrz
         ipar = 1e-6 * 4.6 * fsds_vis(c,j)

         do k = 1, nphytolak
            if (k==small_phyto) then
               Vch = LakeBGCParamsInst%Vchs
               phAlpha = LakeBGCParamsInst%phAlphas
               phBeta = LakeBGCParamsInst%phBetas
               Ktp = LakeBGCParamsInst%Ktps
            else if (k==large_phyto) then
               Vch = LakeBGCParamsInst%Vchl
               phAlpha = LakeBGCParamsInst%phAlphal
               phBeta = LakeBGCParamsInst%phBetal
               Ktp = LakeBGCParamsInst%Ktpl
            end if

            ! temperature factor
            ftemp = max( 0._r8, ThetaG**(tw-20._r8) - &
               ThetaG**(kt(k)*(tw-at(k))) + bt(k) )

            ! CO2 factor
            fco2 = conc_gas_wat(c,j,gco2lak) / (Kco2 + conc_gas_wat(c,j,gco2lak)) 

            ! nutrient factor
            ftp = lake_tp(c) / (Ktp + lake_tp(c))

            ! C:Chla ratio
            C2Chl0 = C2Chlmax(k) - Kpc2chl(k) * mu0(k) * ftemp * ftp
            if (icefrac(c,j)<1._r8 .and. ipar>0.01*ipar0) then
               c2chl = C2Chl0 - (C2Chl0 - C2Chlmin(k)) * log(ipar0/ipar) / 4.605_r8
               c2chl = min(max(c2chl, C2Chlmin(k)), C2Chlmax(k))
            else
               c2chl = C2Chlmax(k)
            end if

            chla = biomas_phyto(c,j,k) * (1._r8 - icefrac(c,j)) / c2chl
            lake_chla(c,j) = lake_chla(c,j) + chla

            ! radiation factor
            Vm0 = Vch / c2chl
            fpar = (1.0 - exp(-phAlpha*ipar/Vm0)) * exp(-phBeta*ipar/Vm0)

            gpp_phyto = Vch/8.64e4_r8 * fpar * ftemp * fco2 * ftp * chla

            gpp_vr(j,k)  = catomw * gpp_phyto
            lake_gpp(c)  = lake_gpp(c) + catomw * gpp_phyto * dz_lake(c,j) 

            gsrc(j,go2lak)    = gsrc(j,go2lak)  + gpp_phyto
            gsnk(j,gco2lak)   = gsnk(j,gco2lak) + gpp_phyto
         end do 

      end do
      
      end associate
   end subroutine Photosynthesis

   subroutine AutotrophicR(c, lakestate_vars, lakebgc_vars, bloss_vr, &
               bphyto_dead, gsrc, gsnk)
      !
      ! !DESCRIPTION:
      ! Simulate phytoplankton autotrophic respiration 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak, nphytolak
      use LakeBGCType        , only : small_phyto, large_phyto
      use LakeBGCType        , only : go2lak, gco2lak
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(out)   :: bloss_vr(1:nlevlak,1:nphytolak)     ! gC/m3/s
      real(r8)               , intent(out)   :: bphyto_dead(1:nphytolak)            ! gC/m2/s
      real(r8)               , intent(inout) :: gsrc(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
      real(r8)               , intent(inout) :: gsnk(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: ThetaML = 1.073 ! temperature multiplier for metabolic loss 
      real(r8), parameter :: frres(nphytolak) = (/0.9_r8, 0.75_r8/)  ! fraction of respiration in total metabolic loss
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, phyto_ar 
      real(r8) :: Klr
      integer  :: j, k
      !-------------------------------------------------------------------- 

      associate(                                            &
            dz_lake        => col_pp%dz_lake                   , & ! Input: [real(r8) (:,:)] layer thickness for lake (m)
            z_lake         => col_pp%z_lake                    , & ! Input: [real(r8) (:,:)] layer depth for lake (m) 
            lakedepth      => col_pp%lakedepth                 , & ! Input: [real(r8) (:)] column lake depth (m)

            t_lake         => col_es%t_lake                    , & ! Input: [real(r8) (:,:)] col lake temperature (Kelvin) 

            icethick       => lakestate_vars%lake_icethick_col , & ! Input: [real(r8) (:)] ice thickness (m) (integrated if lakepuddling)
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            lake_gpp       => lakebgc_vars%gpp_col             , & ! Input: [real(r8) (:)] depth-integrated gross primary production (gC/m2/s)
            biomas_phyto   => lakebgc_vars%biomas_phyto_col    , & ! Input: [real(r8) (:,:,:)] phytoplankton biomass (gC/m3)

            lake_npp       => lakebgc_vars%npp_col             & ! Output: [real(r8) (:)] depth-integrated net primary production (gC/m2/s)
            )

      lake_npp(c) = lake_gpp(c)
      bphyto_dead = 0._r8
      
      do j = 1, nlevlak
         
         bloss_vr(j,:) = 0._r8

         tw = t_lake(c,j) - tfrz

         do k = 1, nphytolak
            if (k==small_phyto) then
               Klr = LakeBGCParamsInst%Klrs
            else if (k==large_phyto) then
               Klr = LakeBGCParamsInst%Klrl
            end if
            
            phyto_ar = Klr/8.64e4_r8 * (ThetaML**(tw-20)) * frres(k) * &
               biomas_phyto(c,j,k) * (1._r8 - icefrac(c,j)) / catomw
            
            gsnk(j,go2lak)    = gsnk(j,go2lak)  + phyto_ar
            gsrc(j,gco2lak)   = gsrc(j,gco2lak) + phyto_ar

            bloss_vr(j,k)     = catomw * phyto_ar / frres(k)
            lake_npp(c)       = lake_npp(c) - catomw * phyto_ar * dz_lake(c,j)
            bphyto_dead(k)    = bphyto_dead(k) + bloss_vr(j,k) * (1._r8-frres(k)) * dz_lake(c,j)
         end do 
      
      end do

      end associate
   end subroutine AutotrophicR

   subroutine HeterotrophicR(c, lakestate_vars, lakebgc_vars, soilc_loss, &
               gsrc, gsnk)
      !
      ! !DESCRIPTION:
      ! Simulate microbial heterotrophic respiration 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak, nsoilclak
      use LakeBGCType        , only : go2lak, gco2lak
      use LakeBGCType        , only : pasC, actC
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(inout) :: soilc_loss(1:nlevsoi,1:nsoilclak)   ! gC/m3/s
      real(r8)               , intent(inout) :: gsrc(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
      real(r8)               , intent(inout) :: gsnk(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
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
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            conc_gas_wat   => lakebgc_vars%conc_wat_col        , & ! Input: [real(r8) (:,:,:)] lake dissolved gas [mol/m3]
            conc_gas_sed   => lakebgc_vars%conc_sed_col        , & ! Input: [real(r8) (:,:,:)] sediment dissolved gas [mol/m3] 
            lake_soilc     => lakebgc_vars%soilc_col           & ! Input: [real(r8) (:,:,:)] lake sediment carbon pools (gC/m3)
            )

      ! for water column
      do j = 1, nlevlak

         tw = t_lake(c,j) - tfrz
         c_o2 = conc_gas_wat(c,j,go2lak)

         pco2_BOD = 3.62e-9_r8 * (ThetaCM**(tw-20._r8)) * &
            (c_o2/(Ko2CM+c_o2)) * (1._r8 - icefrac(c,j)) 
         gsrc(j,gco2lak)   = gsrc(j,gco2lak) + pco2_BOD 
         gsnk(j,go2lak)    = gsnk(j,go2lak)  + pco2_BOD
      end do

      ! for sediment
      do j = 1, nlevsoi

         ts = t_soisno(c,j) - tfrz
         c_o2 = conc_gas_sed(c,j,go2lak)

         do k = 1, nsoilclak
            if (k==pasC) then
               Rca = LakeBGCParamsInst%Rcapas
            else if (k==actC) then
               Rca = LakeBGCParamsInst%Rcaact
            end if
            pco2_soilc = Rca * (ThetaCM**(ts-20._r8)) * (c_o2/(Ko2CM+c_o2)) * &
               lake_soilc(c,j,k) / catomw

            gsrc(j+nlevlak,gco2lak) = gsrc(j+nlevlak,gco2lak)  + pco2_soilc
            gsnk(j+nlevlak,go2lak)  = gsnk(j+nlevlak,go2lak)   + pco2_soilc

            soilc_loss(j,k) = soilc_loss(j,k) + catomw * pco2_soilc
         end do
      end do

      end associate
   end subroutine HeterotrophicR

   subroutine Methanogenesis(c, lakestate_vars, lakebgc_vars, soilc_loss, &
               gsrc, gsnk)
      !
      ! !DESCRIPTION:
      ! Simulate methane production 
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz, catomw
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak, nsoilclak
      use LakeBGCType        , only : go2lak, gco2lak, gch4lak
      use LakeBGCType        , only : pasC, actC 
      use LakeBGCType        , only : yedoma_lake, thaw_lake 
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(inout) :: soilc_loss(1:nlevsoi,1:nsoilclak)   ! gC/m3/s
      real(r8)               , intent(inout) :: gsrc(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s 
      real(r8)               , intent(inout) :: gsnk(1:nlevlak+nlevsoi,1:ngaslak)   ! mol/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: etaO2 = 4e2_r8  ! O2 suppression coefficient (m3 water mol-1) 
      real(r8), parameter :: Tpr(nsoilclak) = (/276.65_r8,273.15_r8/)   ! CH4 production reference temperature (K)
      real(r8), parameter :: pHmin = 2.2_r8  ! minimum allowable pH for CH4 production
      real(r8), parameter :: pHmax = 9.0_r8  ! maximum allowable pH for CH4 production
      real(r8), parameter :: oldcarb0 = 29.3e3_r8  ! yedoma permafrost C (gC/m3)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, ts
      real(r8) :: Rc, PQ10
      real(r8) :: c_o2, fo2, ftemp, fph
      real(r8) :: pch4_soilc, pch4_soilc_yedoma
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
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            lake_tp        => lakebgc_vars%tp_col              , & ! Input: [real(r8) (:)] lake column mean total phosphorus [gP/m3]
            lake_ph        => lakebgc_vars%ph_col              , & ! Input: [real(r8) (:)] lake column mean water pH
            lake_type      => lakebgc_vars%ltype_col           , & ! Input: [integer  (:)] lake type identifier 
            conc_gas_sed   => lakebgc_vars%conc_sed_col        , & ! Input: [real(r8) (:,:,:)] sediment dissolved gas [mol/m3] 
            lake_soilc     => lakebgc_vars%soilc_col           , & ! Input: [real(r8) (:,:,:)] lake sediment carbon pools (gC/m3) 

            ch4_prod_wat   => lakebgc_vars%ch4_prod_wat_col    , & ! Output: [real(r8) (:,:)] lake ch4 production rate [mol/m3/s]
            ch4_prod_sed   => lakebgc_vars%ch4_prod_sed_col    & ! Output: [real(r8) (:,:)] sediment ch4 production rate [mol/m3/s]
            )

      ! CH4 production in oxic water of oligotrophic and mesotrophic lakes
      !if (lake_tp(c)<0.024) then
      !   do j = 1, nlevlak
      !      gsrc(j,gch4lak) = gsrc(j,gch4lak) + LakeBGCParamsInst%Rcoxic*gsrc(j,go2lak) 
      !      ch4_prod_wat(c,j) = LakeBGCParamsInst%Rcoxic*gsrc(j,go2lak)
      !   end do
      !else
         gsrc(1:nlevlak,gch4lak) = 0._r8
         ch4_prod_wat(c,:) = 0._r8
      !end if

      ! pH factor (borrowed from wetland CH4 model)
      if (lake_ph(c)>pHmin .and. lake_ph(c)<pHmax) then
         fph = 10._r8**(-0.2235_r8*lake_ph(c)**2_r8 + 2.7727_r8*lake_ph(c) - 8.6_r8)
         fph = min(1._r8, max(0._r8, fph))
      else
         fph = 0._r8
      end if

      ! CH4 production in sediment
      do j = 1, nlevsoi

         ts = t_soisno(c,j)
         
         ! O2 suppression (Tang et al., 2010; Biogeosciences)
         c_o2 = conc_gas_sed(c,j,go2lak)
         fo2 = 1._r8 / (1._r8 + etaO2*c_o2)
         
         ch4_prod_sed(c,j) = 0._r8
         do k = 1, nsoilclak
            if (k==pasC) then
               ! passive C mainly through hydrogenotrophic methanogenesis
               Rc = LakeBGCParamsInst%Rcpas 
               PQ10 = LakeBGCParamsInst%PQ10pas
               ftemp = PQ10**(0.1_r8*(ts-Tpr(k)))
               
               pch4_soilc = 0.25_r8 * Rc * ftemp * fo2 * fph * lake_soilc(c,j,k) / catomw
               pch4_soilc_yedoma = 0._r8

               gsrc(j+nlevlak,gch4lak) = gsrc(j+nlevlak,gch4lak) + pch4_soilc
               gsrc(j+nlevlak,gco2lak) = gsrc(j+nlevlak,gco2lak) + 3._r8*pch4_soilc

               soilc_loss(j,k) = soilc_loss(j,k) + 4._r8*catomw*pch4_soilc
            else if (k==actC) then
               ! active C mainly through acetoclastic methanogenesis
               Rc = LakeBGCParamsInst%Rcact
               PQ10 = LakeBGCParamsInst%PQ10act
               ftemp = PQ10**(0.1_r8*(ts-Tpr(k)))
               
               pch4_soilc = 0.5_r8 * Rc * ftemp * fo2 * fph * lake_soilc(c,j,k) / catomw
               if (lake_type(c)==yedoma_lake .and. j==nlevsoi) then
                  pch4_soilc_yedoma = 0.5_r8 * Rc * ftemp * fo2 * fph * (oldcarb0/3.0) / catomw
               else
                  pch4_soilc_yedoma = 0._r8 
               end if

               gsrc(j+nlevlak,gch4lak) = gsrc(j+nlevlak,gch4lak) + pch4_soilc + pch4_soilc_yedoma
               gsrc(j+nlevlak,gco2lak) = gsrc(j+nlevlak,gco2lak) + pch4_soilc + pch4_soilc_yedoma

               soilc_loss(j,k) = soilc_loss(j,k) + 2._r8*catomw*pch4_soilc
            end if

            ch4_prod_sed(c,j) = ch4_prod_sed(c,j) + pch4_soilc + pch4_soilc_yedoma
         end do
         
      end do

      end associate
   end subroutine Methanogenesis 

   subroutine Methanotrophy(c, lakestate_vars, lakebgc_vars, gsrc, gsnk)
      !
      ! !DESCRIPTION:
      ! Simulate methane oxidation
      !
      ! !USES:
      use elm_varcon         , only : grav, denh2o, tfrz
      use elm_varpar         , only : nlevlak, nlevgrnd, nlevsoi
      use elm_varpar         , only : ngaslak
      use LakeBGCType        , only : go2lak, gch4lak, gco2lak
      ! !ARGUMENTS:
      implicit none
      integer                , intent(in)    :: c
      type(lakestate_type)   , intent(in)    :: lakestate_vars
      type(lakebgc_type)     , intent(inout) :: lakebgc_vars
      real(r8)               , intent(inout) :: gsrc(1:nlevlak+nlevsoi,1:ngaslak)  ! mol/m3/s
      real(r8)               , intent(inout) :: gsnk(1:nlevlak+nlevsoi,1:ngaslak)  ! mol/m3/s
      !
      ! !CONSTANTS
      real(r8), parameter :: Tor = 267.65 ! CH4 oxidation reference temperature (K)
      !
      ! !LOCAL VARIABLES:
      real(r8) :: tw, ts
      real(r8) :: Qch4, OQ10, Kch4, Ko2
      real(r8) :: c_ch4, c_o2, och4_oxic
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
            icefrac        => lakestate_vars%lake_icefrac_col  , & ! Input: [real(r8) (:,:)] mass fraction of lake layer that is frozen

            conc_gas_wat   => lakebgc_vars%conc_wat_col        , & ! Input: [real(r8) (:,:,:)] lake dissolved gas [mol/m3]
            conc_gas_sed   => lakebgc_vars%conc_sed_col        , & ! Input: [real(r8) (:,:,:)] sediment dissolved gas [mol/m3] 

            ch4_oxid_wat   => lakebgc_vars%ch4_oxid_wat_col    , & ! Output: [real(r8) (:,:)] lake ch4 oxidation rate [mol/m3/s]
            ch4_oxid_sed   => lakebgc_vars%ch4_oxid_sed_col    & ! Output: [real(r8) (:,:)] sediment ch4 oxidation rate [mol/m3/s] 
            )

      ! CH4 oxidation in the water column
      do j = 1, nlevlak
         
         tw = t_lake(c,j)
         c_ch4 = conc_gas_wat(c,j,gch4lak)
         c_o2  = conc_gas_sed(c,j,go2lak)

         if (icefrac(c,j)<1e-8_r8) then
            Qch4 = LakeBGCParamsInst%Qch4
            OQ10 = LakeBGCParamsInst%OQ10
            Kch4 = LakeBGCParamsInst%Kch4
            Ko2  = LakeBGCParamsInst%Ko2

            och4_oxic = Qch4 * (OQ10**(0.1_r8*(tw-Tor))) * &
               (c_ch4/(Kch4+c_ch4)) * (c_o2/(Ko2+c_o2))
         else
            och4_oxic = 0._r8
         end if

         gsnk(j,gch4lak) = gsnk(j,gch4lak) + och4_oxic 
         gsnk(j,go2lak)  = gsnk(j,go2lak)  + 2._r8*och4_oxic
         gsrc(j,gco2lak) = gsrc(j,gco2lak) + och4_oxic
         ch4_oxid_wat(c,j) = och4_oxic
      end do      

      ! no CH4 oxidation in sediment
      ch4_oxid_sed(c,:) = 0._r8

      end associate
   end subroutine Methanotrophy

   !----------------------------------------------------------------------- 
   subroutine readLakeBGCParams ( ncid )
      !
      ! !USES:
      use shr_kind_mod , only : r8 => shr_kind_r8
      use ncdio_pio    , only : file_desc_t,ncd_io
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

      tString='phAlphas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%phAlphas = tempr

      tString='phAlphal'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%phAlphal = tempr

      tString='phBetas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%phBetas = tempr

      tString='phBetal'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%phBetal = tempr

      tString='Ktps'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ktps = tempr

      tString='Ktpl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ktpl = tempr

      tString='Klrs'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Klrs = tempr

      tString='Klrl'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Klrl = tempr

      tString = 'frcResusp'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%frcResusp = tempr

      tString='Re'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Re = tempr

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

      tString='PQ10pas'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%PQ10pas = tempr

      tString='PQ10act'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%PQ10act = tempr

      tString='Qch4'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Qch4 = tempr

      tString='OQ10'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%OQ10 = tempr

      tString='Kch4'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Kch4 = tempr

      tString='Ko2'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Ko2 = tempr

      tString='Rcoxic'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      LakeBGCParamsInst%Rcoxic = tempr

   end subroutine 

end module LakeBGCDynMod
