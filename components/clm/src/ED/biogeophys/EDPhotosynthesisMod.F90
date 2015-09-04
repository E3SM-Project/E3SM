module EDPhotosynthesisMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates the photosynthetic fluxes for the ED model
  ! This code is equivalent to the 'photosynthesis' subroutine in PhotosynthesisMod.F90.
  ! We have split this out to reduce merge conflicts until we can pull out
  ! common code used in both the ED and CLM versions.
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use clm_varctl         , only : iulog
  use clm_varcon         , only : namep 
  use decompMod          , only : bounds_type
  use EcophysConType     , only : ecophyscon
  use atm2lndType        , only : atm2lnd_type
  use CanopyStateType    , only : canopystate_type
  use PhotosynthesisType , only : photosyns_type
  use TemperatureType    , only : temperature_type
  use PatchType          , only : pft
  use EDtypesMod         , only : patch, cohort, gridCellEdState, numpft_ed, dinc_ed
  use EDEcophysContype   , only : EDecophyscon
  use EDVecPatchType     , only : EDpft
  use EDBioType          , only : EDbio_type

  implicit none
  save
  private

  !
  ! PRIVATE MEMBER FUNCTIONS:
  public :: Photosynthesis_ED !ED specific photosynthesis routine
 
  !ED variables
  type (patch) , pointer :: currentPatch
  type (cohort), pointer :: currentCohort

  !------------------------------------------------------------------------------

contains
 
  !---------------------------------------------------------
  ! !INTERFACE:
  subroutine Photosynthesis_ED (bounds, fn, filterp, esat_tv, eair, oair, cair, &
       rb, dayl_factor, &
       atm2lnd_vars, temperature_vars, canopystate_vars, photosyns_vars, EDbio_vars)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    !
    ! !REVISION HISTORY:

    ! !USES:
    use clm_time_manager      , only : get_step_size
    use clm_varcon            , only : rgas, tfrz 
    use clm_varpar            , only : nlevcan_ed, nclmax, nlevsoi, mxpft
    use EDPhysiologyMod       , only : root_fraction
    use EDParamsMod           , only : ED_val_grperc
    use perf_mod              , only : t_startf, t_stopf
    use clm_varctl            , only : use_cn
    use quadraticMod          , only : quadratic
    use CNSharedParamsMod     , only : CNParamsShareInst

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                      ! bounds
    integer                , intent(in)    :: fn                          ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)    :: esat_tv(bounds%begp: )      ! saturation vapor pressure at t_veg (Pa)
    real(r8)               , intent(in)    :: eair( bounds%begp: )        ! vapor pressure of canopy air (Pa)
    real(r8)               , intent(in)    :: oair( bounds%begp: )        ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)    :: cair( bounds%begp: )        ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(inout) :: rb( bounds%begp: )          ! boundary layer resistance (s/m)
    real(r8)               , intent(in)    :: dayl_factor( bounds%begp: ) ! scalar (0-1) for daylength
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    type(canopystate_type) , intent(inout) :: canopystate_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    type(EDbio_type)       , intent(inout) :: EDbio_vars

    ! !CALLED FROM:
    ! subroutine CanopyFluxes 

    ! !LOCAL VARIABLES:
    !
    ! local pointers to implicit in variables
    integer , parameter :: psn_type = 2 !c3 or c4. 

    ! Leaf photosynthesis parameters
    real(r8) :: vcmax_z(nclmax,mxpft,nlevcan_ed)  ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8) :: jmax_z(nclmax,mxpft,nlevcan_ed)   ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: tpu_z(nclmax,mxpft,nlevcan_ed)    ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(r8) :: kp_z(nclmax,mxpft,nlevcan_ed)     ! initial slope of CO2 response curve (C4 plants)
    real(r8) :: lmr_z(nclmax,mxpft,nlevcan_ed)    ! initial slope of CO2 response curve (C4 plants)
    real(r8) :: rs_z(nclmax,mxpft,nlevcan_ed)     ! stomatal resistance s/m
    real(r8) :: gs_z(nclmax,mxpft,nlevcan_ed)     ! stomatal conductance m/s

    real(r8) :: ci(nclmax,mxpft,nlevcan_ed)       ! intracellular leaf CO2 (Pa)
    real(r8) :: lnc(mxpft)                        ! leaf N concentration (gN leaf/m^2)
    real(r8) :: kc( bounds%begp:bounds%endp )     ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: ko( bounds%begp:bounds%endp )     ! Michaelis-Menten constant for O2 (Pa)
    real(r8) :: co2_cp( bounds%begp:bounds%endp ) ! CO2 compensation point (Pa)
    real(r8) :: bbbopt(psn_type)                  ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(r8) :: bbb(mxpft)                        ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8) :: mbbopt(psn_type)                  ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    real(r8) :: mbb(mxpft)                        ! Ball-Berry slope of conductance-photosynthesis relationship

    real(r8) :: kn(mxpft)                         ! leaf nitrogen decay coefficient
    real(r8) :: vcmax25top(mxpft)                 ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top(mxpft)                  ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top(mxpft)                   ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top(mxpft)                   ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top(mxpft)                    ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

    real(r8) :: vcmax25                           ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25                            ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25                             ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25                             ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25                              ! leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: kc25                              ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25                              ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: cp25                              ! CO2 compensation point at 25C (Pa)

    real(r8) :: vcmaxha                           ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha                            ! activation energy for jmax (J/mol)
    real(r8) :: tpuha                             ! activation energy for tpu (J/mol)
    real(r8) :: lmrha                             ! activation energy for lmr (J/mol)
    real(r8) :: kcha                              ! activation energy for kc (J/mol)
    real(r8) :: koha                              ! activation energy for ko (J/mol)
    real(r8) :: cpha                              ! activation energy for cp (J/mol)

    real(r8) :: vcmaxhd                           ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd                            ! deactivation energy for jmax (J/mol)
    real(r8) :: tpuhd                             ! deactivation energy for tpu (J/mol)
    real(r8) :: lmrhd                             ! deactivation energy for lmr (J/mol)

    real(r8) :: vcmaxse                           ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse                            ! entropy term for jmax (J/mol/K)
    real(r8) :: tpuse                             ! entropy term for tpu (J/mol/K)
    real(r8) :: lmrse                             ! entropy term for lmr (J/mol/K)

    real(r8) :: vcmaxc                            ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc                             ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: tpuc                              ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: lmrc                              ! scaling factor for high temperature inhibition (25 C = 1.0)

    real(r8) :: qe(psn_type)                      ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(r8) :: fnps                              ! fraction of light absorbed by non-photosynthetic pigments
    real(r8) :: theta_psii                        ! empirical curvature parameter for electron transport rate

    real(r8) :: theta_cj(psn_type)                ! empirical curvature parameter for ac, aj photosynthesis co-limitation
    real(r8) :: theta_ip                          ! empirical curvature parameter for ap photosynthesis co-limitation

                                                  ! Other
    integer  :: c,CL,f,g,iv,j,p,ps,ft             ! indices
    integer  :: NCL_p                             ! number of canopy layers in patch
    real(r8) :: cf                                ! s m**2/umol -> s/m
    real(r8) :: rsmax0                            ! maximum stomatal resistance [s/m]
    real(r8) :: gb                                ! leaf boundary layer conductance (m/s)
    real(r8) :: gb_mol                            ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8) :: cs                                ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: gs_mol                            ! leaf stomatal conductance (umol H2O/m**2/s)
    real(r8) :: gs                                ! leaf stomatal conductance (m/s)
    real(r8) :: hs                                ! fractional humidity at leaf surface (dimensionless)
    real(r8) :: sco                               ! relative specificity of rubisco
    real(r8) :: tl                                ! leaf temperature in photosynthesis temperature function (K)
    real(r8) :: ha                                ! activation energy in photosynthesis temperature function (J/mol)
    real(r8) :: hd                                ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8) :: se                                ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8) :: cc2                               ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: ciold                             ! previous value of Ci for convergence check
    real(r8) :: gs_mol_err                        ! gs_mol for error check
    real(r8) :: je                                ! electron transport rate (umol electrons/m**2/s)
    real(r8) :: qabs                              ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: aquad,bquad,cquad                 ! terms for quadratic equations
    real(r8) :: r1,r2                             ! roots of quadratic equation
    real(r8) :: ceair                             ! vapor pressure of air, constrained (Pa)
    real(r8) :: act25                             ! (umol/mgRubisco/min) Rubisco activity at 25 C
    integer  :: niter                             ! iteration loop index
    real(r8) :: nscaler                           ! leaf nitrogen scaling coefficient
    real(r8) :: leaf_frac                         ! ratio of to leaf biomass to total alive biomass

    real(r8) :: ac                                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: aj                                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: ap                                ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: ag(nclmax,mxpft,nlevcan_ed)       ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(r8) :: an(nclmax,mxpft,nlevcan_ed)       ! net leaf photosynthesis (umol CO2/m**2/s)
    real(r8) :: an_av(nclmax,mxpft,nlevcan_ed)    ! net leaf photosynthesis (umol CO2/m**2/s) averaged over sun and shade leaves.  
    real(r8) :: ai                                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)

    real(r8) :: laican                            ! canopy sum of lai_z
    real(r8) :: vai                               ! leaf and steam area in ths layer. 
    integer  :: exitloop   
    real(r8) :: laifrac
    real(r8) :: tcsoi                             ! Temperature response function for root respiration. 
    real(r8) :: tc                                ! Temperature response function for wood

    real(r8) :: br                                ! Base rate of root respiration.  (gC/gN/s) 
    real(r8) :: q10                               ! temperature dependence of root respiration  
    integer  :: sunsha                            ! sun (1) or shaded (2)  leaves...
    real(r8) :: dr(2)
    real(r8) :: coarse_wood_frac                  ! amount of woody biomass that is coarse... 
    real(r8) :: tree_area
    real(r8) :: gs_cohort

    ! FIX(SPM, 040714) [I]- these should be proper functions...
    real(r8) :: ft1    ! photosynthesis temperature response (statement function)
    real(r8) :: fth    ! photosynthesis temperature inhibition (statement function)
    real(r8) :: fth25  ! scaling factor for photosynthesis temperature inhibition (statement function)
    ! ... get rid of function statements [I]

    real(r8) dtime ! stepsize in seconds

    !------------------------------------------------------------------------------

    !
    ! FIX(SPM, 040714) [I]- these should be proper functions...Jinyun might be doing this in his refactor...check.
    !
    ! Temperature and soil water response functions
    ft1(tl,ha) = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )
    fth(tl,hd,se,cc2) = cc2 / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )
    fth25(hd,se) = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )
    ! ... get rid of function statements [I]

    associate(                                                &
         ED_patch  => EDpft%ED_patch                        , & ! Input:  does this 'p' have any vegetation associated with it?

         c3psn     => ecophyscon%c3psn                      , & ! photosynthetic pathway: 0. = c4, 1. = c3
         slatop    => ecophyscon%slatop                     , & ! specific leaf area at top of canopy, projected area basis [m^2/gC]
         flnr      => ecophyscon%flnr                       , & ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
         woody     => ecophyscon%woody                      , & ! Is vegetation woody or not? 
         fnitr     => ecophyscon%fnitr                      , & ! foliage nitrogen limitation factor (-)
         leafcn    => ecophyscon%leafcn                     , & ! leaf C:N (gC/gN)

         bb_slope  => EDEcophyscon%BB_slope                 , & ! slope of BB relationship

         forc_pbot => atm2lnd_vars%forc_pbot_downscaled_col , & ! Input:  [real(r8) (:)   ] atmospheric pressure (Pa)

         t_soisno  => temperature_vars%t_soisno_col         , & ! Input: [real(r8)  (:,:) ] soil temperature (Kelvin)                                           
         t_veg     => temperature_vars%t_veg_patch          , & ! Input: [real(r8)  (:)   ] vegetation temperature (Kelvin)                                       
         tgcm      => temperature_vars%thm_patch            , & ! Input: [real(r8)  (:)   ] air temperature at agcm reference height (kelvin) 

         elai      => canopystate_vars%elai_patch           , & ! Input:  [real(r8) (:)  ]  one-sided leaf area index with burying by snow                        
         tlai      => canopystate_vars%tlai_patch           , & ! Input:  [real(r8) (:)   ] one-sided leaf area index
         rscanopy  => canopystate_vars%rscanopy_patch       , & ! Output: [real(r8) (:,:)]  canopy resistance s/m  

         psncanopy => photosyns_vars%psncanopy_patch        , & ! Output: [real(r8) (:,:)]  canopy scale photosynthesis umol CO2 /m**2/ s
         lmrcanopy => photosyns_vars%lmrcanopy_patch        , & ! Output: [real(r8) (:,:)]  canopy scale leaf maintenance respiration umol CO2 /m**2/ s

         GCcanopy  => EDbio_vars%GCcanopy_patch               & ! Output:
         )

      !set timestep
      dtime = get_step_size()

      ! Assign local pointers to derived type members (gridcell-level)
      dr(1) = 0.025_r8; dr(2) = 0.015_r8

      ! Peter Thornton: 3/13/09 
      ! Q10 was originally set to 2.0, an arbitrary choice, but reduced to 1.5 as part of the tuning
      ! to improve seasonal cycle of atmospheric CO2 concentration in global
      ! simulatoins
      q10 = 1.5_r8
      Q10 = CNParamsShareInst%Q10

      !==============================================================================!
      ! Photosynthesis and stomatal conductance parameters, from:
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      !==============================================================================!

      ! vcmax25 parameters, from CN

      act25 = 3.6_r8   !umol/mgRubisco/min
      ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
      act25 = act25 * 1000.0_r8 / 60.0_r8

      ! Activation energy, from:
      ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
      ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

      kcha    = 79430._r8
      koha    = 36380._r8
      cpha    = 37830._r8
      vcmaxha = 65330._r8
      jmaxha  = 43540._r8
      tpuha   = 53100._r8
      lmrha   = 46390._r8

      ! High temperature deactivation, from:
      ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
      ! The factor "c" scales the deactivation to a value of 1.0 at 25C

      vcmaxhd = 149250._r8
      jmaxhd  = 152040._r8
      tpuhd   = 150650._r8
      lmrhd   = 150650._r8

      vcmaxse = 485._r8
      jmaxse  = 495._r8
      tpuse   = 490._r8
      lmrse   = 490._r8

      vcmaxc = fth25(vcmaxhd, vcmaxse)
      jmaxc  = fth25(jmaxhd, jmaxse)
      tpuc   = fth25(tpuhd, tpuse)
      lmrc   = fth25(lmrhd, lmrse)

      ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593

      fnps = 0.15_r8
      theta_psii = 0.7_r8
      theta_ip = 0.95_r8

      qe(1) = 0._r8
      theta_cj(1) = 0.98_r8
      bbbopt(1) = 10000._r8
      mbbopt(1) = 9._r8

      qe(2) = 0.05_r8
      theta_cj(2) = 0.80_r8
      bbbopt(2) = 40000._r8
      mbbopt(2) = 4._r8

      do f = 1,fn
         p = filterp(f)
         call t_startf('edfluxes')
         !these are zeroed even if there's no patch. 
         psncanopy(p) = 0._r8
         lmrcanopy(p) = 0._r8
         rscanopy(p) = 0._r8
         GCcanopy(p) = 0._r8  

         if(ED_patch(p) == 1)then
            g = pft%gridcell(p)

            currentPatch => gridCellEdState(g)%spnt%oldest_patch   
            do while(p /= currentPatch%clm_pno)
               currentPatch => currentPatch%younger
            enddo


            currentPatch%ncan(:,:) = 0
            !redo the canopy structure algorithm to get round a bug that is happening for site 125, FT13. 
            currentCohort => currentPatch%tallest
            do while(associated(currentCohort))

               currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft) = &
                    max(currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft),currentCohort%NV)

               currentCohort => currentCohort%shorter

            enddo !cohort   

            currentPatch%nrad = currentPatch%ncan
            do CL = 1,nclmax
               do ft = 1,numpft_ed
                  currentPatch%present(CL,ft) = 0
                  do iv = 1, currentPatch%nrad(CL,ft);
                     if(currentPatch%canopy_area_profile(CL,ft,iv) > 0._r8)then
                        currentPatch%present(CL,ft) = 1
                     end if
                  end do !iv     
               enddo !ft
            enddo !CL


            g = pft%gridcell(p)
            c = pft%column(p)
            ! Soil water stress applied to Ball-Berry parameters
            do FT = 1,numpft_ed
               if (nint(c3psn(FT)) == 1)then
                  ps = 1
               else
                  ps = 2
               end if
               bbb(FT) = max (bbbopt(ps)*currentPatch%btran_ft(FT), 1._r8)

               mbb(FT) = bb_slope(ft) ! mbbopt(ps)
            end do

            ! kc, ko, currentPatch, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
            !
            !       kc25 = 404.9 umol/mol
            !       ko25 = 278.4 mmol/mol
            !       cp25 = 42.75 umol/mol
            !
            ! Derive sco from currentPatch and O2 using present-day O2 (0.209 mol/mol) and re-calculate
            ! currentPatch to account for variation in O2 using currentPatch = 0.5 O2 / sco
            !

            kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(c)
            ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(c)
            sco  = 0.5_r8 * 0.209_r8 / (42.75_r8 / 1.e06_r8)
            cp25 = 0.5_r8 * oair(p) / sco

            if(t_veg(p).gt.150_r8.and.t_veg(p).lt.350_r8)then
               kc(p) = kc25 * ft1(t_veg(p), kcha)
               ko(p) = ko25 * ft1(t_veg(p), koha)
               co2_cp(p) = cp25 * ft1(t_veg(p), cpha)
            else
               kc(p) = 1
               ko(p) = 1
               co2_cp(p) = 1
               write(iulog,*) 'something wrong with temperature',t_veg(p),p,elai(p),tlai(p)
            end if

         end if
      end do

      ! Multi-layer parameters scaled by leaf nitrogen profile.
      ! Loop through each canopy layer to calculate nitrogen profile using
      ! cumulative lai at the midpoint of the layer


      do f = 1,fn
         p = filterp(f)
         c = pft%column(p)
         if(ED_patch(p) == 1)then
            g = pft%gridcell(p)
            currentPatch => gridCellEdState(g)%spnt%oldest_patch   
            do while(p /= currentPatch%clm_pno)
               currentPatch => currentPatch%younger
            enddo

            do FT = 1,numpft_ed
               if (nint(c3psn(FT)) == 1)then
                  ps = 1
               else
                  ps = 2
               end if
               bbb(FT) = max (bbbopt(ps)*currentPatch%btran_ft(FT), 1._r8)
               mbb(FT) = mbbopt(ps)

               if (nint(c3psn(FT)) == 1)then
                  ci(:,FT,:) = 0.7_r8 * cair(p)
               else
                  ci(:,FT,:) = 0.4_r8 * cair(p)
               end if
            enddo

            NCL_p = currentPatch%NCL_p

            do FT = 1,numpft_ed !calculate patch and pft specific propserties at canopy top. 
               ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
               lnc(FT) = 1._r8 / (slatop(FT) * leafcn(FT))

               !at the moment in ED we assume that there is no active N cycle. This should change, of course. FIX(RF,032414) Sep2011. 
               vcmax25top(FT) = fnitr(FT) !fudge - shortcut using fnitr as a proxy for vcmax... 

               ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
               ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.
               ! Here use a factor "1.67", from Medlyn et al (2002) Plant, Cell and Environment 25:1167-1179

               jmax25top(FT) = 1.67_r8    * vcmax25top(FT)
               !
               ! FIX(RF, 082914) - Changed to 0.167 in ED branch...should be a parameter
               !
               tpu25top(FT)  = 0.167_r8   * jmax25top(FT)
               kp25top(FT)   = 20000._r8  * vcmax25top(FT)

               ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
               ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
               ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25

               if (dayl_factor(p)  ==  0._r8) then
                  kn(FT) =  0._r8
               else
                  kn(FT) = exp(0.00963_r8 * vcmax25top(FT) - 2.43_r8)
               end if

               if (use_cn) then
                  ! Leaf maintenance respiration to match the base rate used in CN
                  ! but with the new temperature functions for C3 and C4 plants.
                  !
                  ! Base rate for maintenance respiration is from:
                  ! M. Ryan, 1991. Effects of climate change on plant respiration.
                  ! Ecological Applications, 1(2), 157-167.
                  ! Original expression is br = 0.0106 molC/(molN h)
                  ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
                  !
                  ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
                  !
                  ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
                  ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
                  !
                  ! Then scale this value at the top of the canopy for canopy depth

                  lmr25top(FT) = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
                  lmr25top(FT) = lmr25top(FT) * lnc(FT) / 12.e-06_r8
               else
                  ! Leaf maintenance respiration in proportion to vcmax25top
                  lmr25top(FT) = vcmax25top(FT) * dr(nint(c3psn(FT)+1))
               end if
            end do !FT 

            !==============================================================================!   
            ! Calculate Nitrogen scaling factors and photosynthetic parameters.         
            !==============================================================================!
            do CL = 1, NCL_p
               do FT = 1,numpft_ed

                  do iv = 1, currentPatch%nrad(CL,FT)
                     if(currentPatch%canopy_area_profile(CL,FT,iv)>0._r8.and.currentPatch%present(CL,FT) /= 1)then
                        write(iulog,*) 'CF: issue with present structure',CL,FT,iv, &
                             currentPatch%canopy_area_profile(CL,FT,iv),currentPatch%present(CL,FT), &
                             currentPatch%nrad(CL,FT),currentPatch%ncl_p,nclmax
                        currentPatch%present(CL,FT) = 1
                     end if
                  enddo

                  if(currentPatch%present(CL,FT) == 1)then ! are there any leaves of this pft in this layer?     

                     if(CL==NCL_p)then !are we in the top canopy layer or a shaded layer?      
                        laican = 0._r8
                     else
                        laican = sum(currentPatch%canopy_layer_lai(CL+1:NCL_p)) 
                     end if

                     ! Loop through canopy layers (above snow). Respiration needs to be
                     ! calculated every timestep. Others are calculated only if daytime    
                     do iv = 1, currentPatch%nrad(CL,FT)
                        vai = (currentPatch%elai_profile(CL,FT,iv)+currentPatch%esai_profile(CL,FT,iv)) !vegetation area index. 
                        if (iv == 1) then
                           laican = laican + 0.5_r8 * vai
                        else
                           laican = laican + 0.5_r8 * (currentPatch%elai_profile(CL,FT,iv-1)+ &
                                currentPatch%esai_profile(CL,FT,iv-1))+vai
                        end if

                        ! Scale for leaf nitrogen profile
                        nscaler = exp(-kn(FT) * laican)


                        ! Maintenance respiration: umol CO2 / m**2 [leaf] / s
                        lmr25 = lmr25top(FT) * nscaler

                        if (nint(c3psn(FT)) == 1)then
                           lmr_z(CL,FT,iv) = lmr25 * ft1(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
                        else
                           lmr_z(CL,FT,iv) = lmr25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                           lmr_z(CL,FT,iv) = lmr_z(CL,FT,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
                        end if


                        if (currentPatch%ed_parsun_z(CL,FT,iv) <= 0._r8) then           ! night time
                           vcmax_z(CL,FT,iv) = 0._r8
                           jmax_z(CL,FT,iv) = 0._r8
                           tpu_z(CL,FT,iv) = 0._r8
                           kp_z(CL,FT,iv) = 0._r8
                        else                                     ! day time
                           vcmax25 = vcmax25top(FT) * nscaler
                           jmax25 = jmax25top(FT) * nscaler
                           tpu25 = tpu25top(FT) * nscaler
                           kp25 = kp25top(FT) * nscaler

                           ! Adjust for temperature
                           vcmax_z(CL,FT,iv) = vcmax25 * ft1(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
                           jmax_z(CL,FT,iv)  = jmax25 * ft1(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
                           tpu_z(CL,FT,iv)   = tpu25 * ft1(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

                           if (nint(c3psn(FT))  /=  1) then
                              vcmax_z(CL,FT,iv) = vcmax25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                              vcmax_z(CL,FT,iv) = vcmax_z(CL,FT,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                              vcmax_z(CL,FT,iv) = vcmax_z(CL,FT,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
                           end if
                           kp_z(CL,FT,iv) = kp25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8) !q10 response of product limited psn. 
                        end if
                        ! Adjust for soil water:(umol co2/m**2/s)

                        vcmax_z(CL,FT,iv) = vcmax_z(CL,FT,iv) * currentPatch%btran_ft(FT)
                        ! completely removed respiration drought response 
                        ! - (lmr_z(CL,FT,iv) * (1.0_r8-currentPatch%btran_ft(FT))  *ecophyscon%resp_drought_response(FT))
                        lmr_z(CL,FT,iv) = lmr_z(CL,FT,iv) 

                     end do ! iv
                  end if !present
               enddo !PFT 
            enddo !CL

            !==============================================================================!
            ! Leaf-level photosynthesis and stomatal conductance
            !==============================================================================!

            rsmax0 = 2.e4_r8

            ! Leaf boundary layer conductance, umol/m**2/s

            cf = forc_pbot(c)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
            gb = 1._r8/rb(p)
            gb_mol = gb * cf
            ! Constrain eair >= 0.05*esat_tv so that solution does not blow up. This ensures
            ! that hs does not go to zero. Also eair <= esat_tv so that hs <= 1

            ceair = min( max(eair(p), 0.05_r8*esat_tv(p)), esat_tv(p) ) 
            ! Loop through canopy layers (above snow). Only do calculations if daytime
            do CL = 1, NCL_p
               do FT = 1,numpft_ed
                  if (nint(c3psn(FT)) == 1)then
                     ps = 1
                  else
                     ps = 2
                  end if
                  if(currentPatch%present(CL,FT) == 1)then ! are there any leaves of this pft in this layer?     
                     do iv = 1, currentPatch%nrad(CL,FT)
                        if (currentPatch%ed_parsun_z(CL,FT,iv) <= 0._r8) then  ! night time
                           ac = 0._r8
                           aj = 0._r8
                           ap = 0._r8
                           ag(CL,FT,iv) = 0._r8
                           an(CL,FT,iv) = ag(CL,FT,iv) - lmr_z(CL,FT,iv)
                           an_av(cl,ft,iv) = 0._r8
                           currentPatch%psn_z(cl,ft,iv) = 0._r8
                           rs_z(CL,FT,iv) = min(rsmax0, 1._r8/bbb(FT) * cf)


                        else ! day time
                           !is there leaf area? - (NV can be larger than 0 with only stem area if deciduous)
                           if(currentPatch%ed_laisun_z(CL,ft,iv)+currentPatch%ed_laisha_z(cl,ft,iv) > 0._r8)then 
                              !Loop aroun shaded and unshaded leaves          
                              currentPatch%psn_z(CL,ft,iv) = 0._r8    ! psn is accumulated across sun and shaded leaves. 
                              rs_z(CL,FT,iv)  = 0._r8                 ! 1/rs is accumulated across sun and shaded leaves. 
                              gs_z(CL,FT,iv) = 0._r8
                              an_av(CL,FT,iv) = 0._r8
                              do  sunsha = 1,2      
                                 ! Electron transport rate for C3 plants. Convert par from W/m2 to umol photons/m**2/s 
                                 ! using the factor 4.6
                                 ! Convert from units of par absorbed per unit ground area to par absorbed per unit leaf area. 

                                 if(sunsha == 1)then !sunlit
                                    if((currentPatch%ed_laisun_z(CL,FT,iv) * currentPatch%canopy_area_profile(CL,FT,iv)) >  &
                                         0.0000000001_r8)then

                                       qabs = currentPatch%ed_parsun_z(CL,FT,iv) / (currentPatch%ed_laisun_z(CL,FT,iv) * &
                                            currentPatch%canopy_area_profile(CL,FT,iv))   
                                       qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 

                                    else
                                       qabs = 0.0_r8
                                    end if
                                 else

                                    qabs = currentPatch%ed_parsha_z(CL,FT,iv) / (currentPatch%ed_laisha_z(CL,FT,iv) * &
                                         currentPatch%canopy_area_profile(CL,FT,iv))  
                                    qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 

                                 end if

                                 !convert the absorbed par into absorbed par per m2 of leaf, 
                                 ! so it is consistant with the vcmax and lmr numbers. 
                                 aquad = theta_psii
                                 bquad = -(qabs + jmax_z(cl,ft,iv))
                                 cquad = qabs * jmax_z(cl,ft,iv)
                                 call quadratic (aquad, bquad, cquad, r1, r2)
                                 je = min(r1,r2)

                                 ! Iterative loop for ci beginning with initial guess
                                 if (nint(c3psn(FT)) == 1)then
                                    ci(cl,ft,iv) = 0.7_r8 * cair(p)
                                 else
                                    ci(cl,ft,iv) = 0.4_r8 * cair(p)
                                 end if

                                 niter = 0
                                 exitloop = 0
                                 do while(exitloop == 0)                 
                                    ! Increment iteration counter. Stop if too many iterations
                                    niter = niter + 1

                                    ! Save old ci
                                    ciold = ci(cl,ft,iv)

                                    ! Photosynthesis limitation rate calculations 
                                    if (nint(c3psn(FT)) == 1)then
                                       ! C3: Rubisco-limited photosynthesis
                                       ac = vcmax_z(cl,ft,iv) * max(ci(cl,ft,iv)-co2_cp(p), 0._r8) / (ci(cl,ft,iv)+kc(p)* &
                                            (1._r8+oair(p)/ko(p)))
                                       ! C3: RuBP-limited photosynthesis
                                       aj = je * max(ci(cl,ft,iv)-co2_cp(p), 0._r8) / (4._r8*ci(cl,ft,iv)+8._r8*co2_cp(p))
                                       ! C3: Product-limited photosynthesis 
                                       ap = 3._r8 * tpu_z(cl,ft,iv)
                                    else
                                       ! C4: Rubisco-limited photosynthesis
                                       ac = vcmax_z(cl,ft,iv)
                                       ! C4: RuBP-limited photosynthesis
                                       if(sunsha == 1)then !sunlit
                                          if((currentPatch%ed_laisun_z(cl,ft,iv) * currentPatch%canopy_area_profile(cl,ft,iv)) > &
                                               0.0000000001_r8)then !guard against /0's in the night.             
                                             aj = qe(ps) * currentPatch%ed_parsun_z(cl,ft,iv) * 4.6_r8
                                             !convert from per cohort to per m2 of leaf)
                                             aj = aj / (currentPatch%ed_laisun_z(cl,ft,iv) * &
                                                        currentPatch%canopy_area_profile(cl,ft,iv))
                                          else
                                             aj = 0._r8
                                          end if
                                       else
                                          aj = qe(ps) * currentPatch%ed_parsha_z(cl,ft,iv) * 4.6_r8
                                          aj = aj / (currentPatch%ed_laisha_z(cl,ft,iv) * &
                                                     currentPatch%canopy_area_profile(cl,ft,iv))         
                                       end if

                                       ! C4: PEP carboxylase-limited (CO2-limited)
                                       ap = kp_z(cl,ft,iv) * max(ci(cl,ft,iv), 0._r8) / forc_pbot(c)
                                    end if
                                    ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap
                                    aquad = theta_cj(ps)
                                    bquad = -(ac + aj)
                                    cquad = ac * aj
                                    call quadratic (aquad, bquad, cquad, r1, r2)
                                    ai = min(r1,r2)

                                    aquad = theta_ip
                                    bquad = -(ai + ap)
                                    cquad = ai * ap
                                    call quadratic (aquad, bquad, cquad, r1, r2)
                                    ag(cl,ft,iv) = min(r1,r2)

                                    ! Net carbon assimilation. Exit iteration if an < 0
                                    an(cl,ft,iv) = ag(cl,ft,iv) - lmr_z(cl,ft,iv)
                                    if (an(cl,ft,iv) < 0._r8) then
                                       exitloop = 1
                                    end if

                                    ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
                                    ! With an <= 0, then gs_mol = bbb

                                    cs = cair(p) - 1.4_r8/gb_mol * an(cl,ft,iv) * forc_pbot(c)
                                    cs = max(cs,1.e-06_r8)
                                    aquad = cs
                                    bquad = cs*(gb_mol - bbb(FT)) - mbb(FT)*an(cl,ft,iv)*forc_pbot(c)
                                    cquad = -gb_mol*(cs*bbb(FT) + mbb(FT)*an(cl,ft,iv)*forc_pbot(c)*ceair/esat_tv(p))
                                    call quadratic (aquad, bquad, cquad, r1, r2)
                                    gs_mol = max(r1,r2)

                                    ! Derive new estimate for ci
                                    ci(cl,ft,iv) = cair(p) - an(cl,ft,iv) * forc_pbot(c) * &
                                         (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)

                                    ! Check for ci convergence. Delta ci/pair = mol/mol. Multiply by 10**6 to
                                    ! convert to umol/mol (ppm). Exit iteration if convergence criteria of +/- 1 x 10**-6 ppm
                                    ! is met OR if at least ten iterations (niter=10) are completed

                                    if ((abs(ci(cl,ft,iv)-ciold)/forc_pbot(c)*1.e06_r8 <=  2.e-06_r8) .or. niter == 5) then
                                       exitloop = 1
                                    end if
                                 end do !iteration loop

                                 ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb
                                 if (an(cl,ft,iv) < 0._r8) then
                                    gs_mol = bbb(FT)
                                 end if

                                 ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)
                                 cs = cair(p) - 1.4_r8/gb_mol * an(cl,ft,iv) * forc_pbot(c)
                                 cs = max(cs,1.e-06_r8)
                                 ci(cl,ft,iv) = cair(p) - an(cl,ft,iv) * forc_pbot(c) * (1.4_r8*gs_mol+1.6_r8*gb_mol) / &
                                      (gb_mol*gs_mol)
                                 ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)
                                 gs = gs_mol / cf

                                 !accumulate total photosynthesis umol/m2 ground/s-1. weight per unit sun and sha leaves.  
                                 if(sunsha == 1)then !sunlit       

                                    currentPatch%psn_z(cl,ft,iv) = currentPatch%psn_z(cl,ft,iv) + ag(cl,ft,iv) * &
                                         currentPatch%f_sun(cl,ft,iv)
                                    an_av(cl,ft,iv) = an_av(cl,ft,iv) + an(cl,ft,iv)                  * &
                                         currentPatch%f_sun(cl,ft,iv) 
                                    gs_z(cl,ft,iv)  = gs_z(cl,ft,iv)  + 1._r8/(min(1._r8/gs, rsmax0)) * &
                                         currentPatch%f_sun(cl,ft,iv) 

                                 else

                                    currentPatch%psn_z(cl,ft,iv) = currentPatch%psn_z(cl,ft,iv) + ag(cl,ft,iv) &
                                         * (1.0_r8-currentPatch%f_sun(cl,ft,iv))                 
                                    an_av(cl,ft,iv)    = an_av(cl,ft,iv)    + an(cl,ft,iv) &
                                         * (1.0_r8-currentPatch%f_sun(cl,ft,iv)) 
                                    gs_z(cl,ft,iv)     = gs_z(cl,ft,iv)     + &
                                         1._r8/(min(1._r8/gs, rsmax0)) * (1.0_r8-currentPatch%f_sun(cl,ft,iv)) 

                                 end if

                                 ! Make sure iterative solution is correct
                                 if (gs_mol < 0._r8) then
                                    write (iulog,*)'Negative stomatal conductance:'
                                    write (iulog,*)'p,iv,gs_mol= ',p,iv,gs_mol
                                    call endrun(decomp_index=p, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
                                 end if

                                 ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b
                                 hs = (gb_mol*ceair + gs_mol*esat_tv(p)) / ((gb_mol+gs_mol)*esat_tv(p))
                                 gs_mol_err = mbb(FT)*max(an(cl,ft,iv), 0._r8)*hs/cs*forc_pbot(c) + bbb(FT)

                                 if (abs(gs_mol-gs_mol_err) > 1.e-01_r8) then
                                    write (iulog,*) 'CF: Ball-Berry error check - stomatal conductance error:'
                                    write (iulog,*) gs_mol, gs_mol_err
                                 end if

                              enddo !sunsha loop
                              !average leaf-level stomatal resistance rate over sun and shade leaves... 
                              rs_z(cl,ft,iv)  = 1._r8/gs_z(cl,ft,iv) 
                           end if !is there leaf area? 
                        end if    ! night or day 
                     end do   ! iv canopy layer 
                  end if    ! present(L,ft) ? rd_array
               end do  ! PFT loop
            end do  !canopy layer

            call t_stopf('edfluxes')
            call t_startf('edunpack')

            !==============================================================================!
            ! Unpack fluxes from arrays into cohorts
            !==============================================================================!
            call root_fraction(currentPatch)

            if(currentPatch%countcohorts > 0.0)then  !avoid errors caused by empty patches 

               currentCohort => currentPatch%tallest  ! Cohort loop

               do while (associated(currentCohort)) ! Cohort loop
                  call t_startf('edfluxunpack1')
                  if(currentCohort%n > 0._r8)then   
                     ! Zero cohort flux accumulators.
                     currentCohort%npp_clm    = 0._r8
                     currentCohort%resp_clm   = 0._r8

                     ! Select canopy layer and PFT.
                     FT = currentCohort%pft  !are we going to have ftindex?
                     CL  = currentCohort%canopy_layer
                     !------------------------------------------------------------------------------
                     ! Accumulate fluxes over the sub-canopy layers of each cohort.
                     !------------------------------------------------------------------------------
                     ! Convert from umolC/m2leaf/s to umolC/indiv/s ( x canopy area x 1m2 leaf area). 
                     tree_area = currentCohort%c_area/currentCohort%n
                     if(currentCohort%nv > 1)then

                        currentCohort%gpp_clm  = sum(currentPatch%psn_z(cl,ft,1:currentCohort%nv-1) * &
                             currentPatch%elai_profile(cl,ft,1:currentCohort%nv-1)) * tree_area
                        currentCohort%rd       = sum(lmr_z(cl,ft,1:currentCohort%nv-1)    * &
                             currentPatch%elai_profile(cl,ft,1:currentCohort%nv-1)) * tree_area 

                        currentCohort%gscan    = sum((1.0_r8/(rs_z(cl,ft,1:currentCohort%nv-1)+rb(p)))) *  tree_area 
                        currentCohort%ts_net_uptake(1:currentCohort%nv) = an_av(cl,ft,1:currentCohort%nv) * 12E-9 * dtime 

                     else

                        currentCohort%gpp_clm = 0.0_r8 
                        currentCohort%rd = 0._r8 
                        currentCohort%gscan = 0._r8 
                        currentCohort%ts_net_uptake(:) = 0._r8

                     end if

                     laifrac = (currentCohort%treelai+currentCohort%treesai)-(currentCohort%nv-1)*dinc_ed

                     gs_cohort = 1.0_r8/(rs_z(cl,ft,currentCohort%nv)+rb(p))*laifrac*tree_area   
                     currentCohort%gscan = currentCohort%gscan+gs_cohort

                     currentCohort%gpp_clm  = currentCohort%gpp_clm + currentPatch%psn_z(cl,ft,currentCohort%nv) * &
                          currentPatch%elai_profile(cl,ft,currentCohort%nv) * laifrac * tree_area
                     currentCohort%rd       = currentCohort%rd      + lmr_z(cl,ft,currentCohort%nv)    * &
                          currentPatch%elai_profile(cl,ft,currentCohort%nv) * laifrac * tree_area 

                     call t_stopf('edfluxunpack1')      
                     call t_startf('edfluxunpack2')

                     !------------------------------------------------------------------------------
                     ! Calculate Whole Plant Respiration (this doesn't really need to be in this iteration at all, surely?)    
                     ! Leaf respn needs to be in the sub-layer loop to account for changing N through canopy. 
                     !
                     ! base rate for maintenance respiration is from:
                     ! M. Ryan, 1991. Effects of climate change on plant respiration.
                     ! Ecological Applications, 1(2), 157-167.
                     ! Original expression is br = 0.0106 molC/(molN h)
                     ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
                     !------------------------------------------------------------------------------

                     br = 2.525e-6_r8

                     leaf_frac = 1.0_r8/(currentCohort%canopy_trim + EDecophyscon%sapwood_ratio(currentCohort%pft) * &
                          currentCohort%hite + ecophyscon%froot_leaf(currentCohort%pft))
                     currentCohort%bsw = EDecophyscon%sapwood_ratio(currentCohort%pft) * currentCohort%hite * &
                          (currentCohort%balive + currentCohort%laimemory)*leaf_frac
                     currentCohort%livestemn  = currentCohort%bsw  / ecophyscon%leafcn(currentCohort%pft)

                     currentCohort%livestem_mr  = 0._r8
                     currentCohort%livecroot_mr = 0._r8
                     if (woody(FT) == 1) then
                        tc = q10**((t_veg(p)-tfrz - 20.0_r8)/10.0_r8) 
                        currentCohort%livestem_mr  = currentCohort%livestemn  * br * tc !*currentPatch%btran_ft(currentCohort%pft)  
                        currentCohort%livecroot_mr = currentCohort%livecrootn * br * tc !*currentPatch%btran_ft(currentCohort%pft) 

                        !convert from gC /indiv/s-1 to kgC/indiv/s-1
                        currentCohort%livestem_mr  = currentCohort%livestem_mr /1000
                        currentCohort%livecroot_mr = currentCohort%livecroot_mr /1000
                     else
                        tc = 1.0_r8
                        currentCohort%livestem_mr  = 0._r8
                        currentCohort%livecroot_mr = 0._r8    
                     end if

                     if (ecophyscon%woody(currentCohort%pft) == 1) then
                        coarse_wood_frac = 0.5_r8
                     else
                        coarse_wood_frac = 0.0_r8
                     end if

                     ! Soil temperature. 
                     currentCohort%froot_mr = 0._r8

                     do j = 1,nlevsoi
                        tcsoi  = q10**((t_soisno(c,j)-tfrz - 20.0_r8)/10.0_r8)
                        !fine root respn. 
                        currentCohort%froot_mr = currentCohort%froot_mr + (1.0_r8 - coarse_wood_frac) * &
                             currentCohort%br*br*tcsoi * currentPatch%rootfr_ft(ft,j)/leafcn(currentCohort%pft) 
                        ! convert from gC/indiv/s-1 to kgC/indiv/s-1
                        currentCohort%froot_mr  =  currentCohort%froot_mr /1000_r8
                     enddo

                     call t_stopf('edfluxunpack2')
                     call t_startf('edfluxunpack3')
                     ! convert gpp and resp from umol/indiv/s-1 to kgC/indiv/s-1  = X * 12 *10-6 * 10-3
                     !currentCohort%resp_m  = currentCohort%rd      * 12.0E-9
                     currentCohort%gpp_clm = currentCohort%gpp_clm * 12.0E-9    
                     ! add on whole plant respiration values in kgC/indiv/s-1  
                     currentCohort%resp_m = currentCohort%livestem_mr + currentCohort%livecroot_mr + currentCohort%froot_mr
                     ! no drought response * (1.0_r8 - currentPatch%btran_ft(currentCohort%pft)*ecophyscon%resp_drought_response(FT))   
                     currentCohort%resp_m = currentCohort%resp_m + currentCohort%rd * 12.0E-9 !this was already corrected fo BTRAN     

                     ! convert from kgC/indiv/s to kgC/indiv/timestep       
                     currentCohort%resp_m   = currentCohort%resp_m  * dtime 
                     currentCohort%gpp_clm  = currentCohort%gpp_clm * dtime          
                     currentCohort%resp_g   = ED_val_grperc * (max(0._r8,currentCohort%gpp_clm - currentCohort%resp_m))
                     currentCohort%resp_clm = currentCohort%resp_m + currentCohort%resp_g ! kgC/indiv/ts    
                     currentCohort%npp_clm  = currentCohort%gpp_clm - currentCohort%resp_clm  ! kgC/indiv/ts

                     !------------------------------------------------------------------------------
                     ! Remove whole plant respiration from net uptake.  (kgC/indiv/ts)
                     if(currentCohort%treelai > 0._r8)then     
                        ! do iv =1,currentCohort%NV
                        ! currentCohort%year_net_uptake(iv) = currentCohort%year_net_uptake(iv) - &
                        ! (timestep_secs*(currentCohort%livestem_mr + currentCohort%livecroot_mr &
                        ! minus contribution to whole plant respn.
                        ! + currentCohort%froot_mr))/(currentCohort%treelai*currentCohort%c_area/currentCohort%n)
                        ! enddo
                     else !lai<0   
                        currentCohort%gpp_clm = 0._r8
                        currentCohort%resp_m = 0._r8
                        currentCohort%gscan = 0._r8
                     end if
                  else !pft<0 n<0
                     write(iulog,*) 'CF: pft 0 or n 0',currentCohort%pft,currentCohort%n,currentCohort%indexnumber
                     currentCohort%gpp_clm = 0._r8
                     currentCohort%resp_m  = 0._r8
                     currentCohort%gscan   = 0._r8   
                     currentCohort%ts_net_uptake(1:currentCohort%nv) = 0._r8    
                  end if !pft<0 n<0

                  psncanopy(p) = psncanopy(p) + currentCohort%gpp_clm
                  lmrcanopy(p) = lmrcanopy(p) + currentCohort%resp_m
                  ! accumulate cohort level canopy conductances over whole area before dividing by total area.  
                  GCcanopy(p)  = GCcanopy(p) + currentCohort%gscan * currentCohort%n /currentPatch%total_canopy_area  

                  currentCohort => currentCohort%shorter

               enddo  ! end cohort loop.   
            end if !count_cohorts is more than zero.

            psncanopy(p) = psncanopy(p) / currentPatch%area
            lmrcanopy(p) = lmrcanopy(p) / currentPatch%area
            if(GCcanopy(p) > 1._r8/rsmax0.and.elai(p) > 0.0_r8)then
               rscanopy(p)  = (1.0_r8/GCcanopy(p))-rb(p)/elai(p) ! this needs to be resistance per unit leaf area.  
            else
               rscanopy(p) = rsmax0
            end if
            GCcanopy(p)  = 1.0_r8/rscanopy(p) *cf /1000 !convert into umol m02 s-1 then mmol m-2 s-1. 
         else !EDpatch
            rscanopy(p) = rsmax0  
         end if !edpatch

         call t_stopf('edfluxunpack3')
         call t_stopf('edunpack')  

      end do !patch loop

    end associate

  end subroutine Photosynthesis_ED

end module EDPhotosynthesisMod
