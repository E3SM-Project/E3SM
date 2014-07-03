module CNNitrifDenitrifMod

  !-----------------------------------------------------------------------
  ! Calculate nitrification and denitrification rates
  !
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_TKFRZ
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  use clm_varcon   , only: secspday
  use clm_varctl   , only: use_lch4
  use abortutils   , only: endrun
  use decompMod    , only: bounds_type
  !
  implicit none
  save
  private
  !
  public :: nitrif_denitrif
  public :: readCNNitrifDenitrifParams
  !
  type, private :: CNNitrifDenitrifParamsType
   real(r8) :: k_nitr_max               !  maximum nitrification rate constant (1/s)
   real(r8) :: surface_tension_water    !  surface tension of water(J/m^2), Arah an and Vinten 1995
   real(r8) :: rij_kro_a                !  Arah and Vinten 1995)
   real(r8) :: rij_kro_alpha            !  parameter to calculate anoxic fraction of soil  (Arah and Vinten 1995)
   real(r8) :: rij_kro_beta             !  (Arah and Vinten 1995)
   real(r8) :: rij_kro_gamma            !  (Arah and Vinten 1995)
   real(r8) :: rij_kro_delta            !  (Arah and Vinten 1995)
  end type CNNitrifDenitrifParamsType

  type(CNNitrifDenitrifParamsType),private ::  CNNitrifDenitrifParamsInst

  logical, public :: no_frozen_nitrif_denitrif = .false.  ! stop nitrification and denitrification in frozen soils
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------  
  subroutine readCNNitrifDenitrifParams ( ncid )
    !
    !
    use ncdio_pio    , only: file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNNitrifDenitrifParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------
    !
    ! read in constants
    !
    tString='k_nitr_max'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%k_nitr_max=tempr

    tString='surface_tension_water'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%surface_tension_water=tempr

    tString='rij_kro_a'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_a=tempr

    tString='rij_kro_alpha'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_alpha=tempr

    tString='rij_kro_beta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_beta=tempr

    tString='rij_kro_gamma'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_gamma=tempr

    tString='rij_kro_delta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNNitrifDenitrifParamsInst%rij_kro_delta=tempr

  end subroutine readCNNitrifDenitrifParams

  !-----------------------------------------------------------------------
  subroutine nitrif_denitrif(bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    !  calculate nitrification and denitrification rates
    !
    ! !USES:
    use clmtype
    use clm_varpar         , only: nlevgrnd,nlevdecomp
    use clm_time_manager   , only: get_curr_date, get_step_size
    use shr_const_mod      , only: SHR_CONST_TKFRZ
    use clm_varctl         , only: iulog
    use clm_varcon         , only: rpi, denh2o, dzsoi, zisoi, grav
    use clm_varcon         , only: d_con_g, d_con_w
    use CNSharedParamsMod  , only: anoxia_wtsat,CNParamsShareInst
    use clm_varcon         , only: spval
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_soilc         ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, reflev, j
    real(r8) :: soil_hr_vr(bounds%begc:bounds%endc,1:nlevdecomp) ! total soil respiration rate (g C / m3 / s)
    real(r8) :: g_per_m3__to__ug_per_gsoil
    real(r8) :: g_per_m3_sec__to__ug_per_gsoil_day
    real(r8) :: k_nitr_max                          ! maximum nitrification rate constant (1/s)
    real(r8) :: mu, sigma
    real(r8) :: t
    real(r8) :: pH(bounds%begc:bounds%endc)
    !debug-- put these in clmtype for outing to hist files
    real(r8) :: co2diff_con(2)                      ! diffusion constants for CO2
    real(r8) :: eps
    real(r8) :: f_a
    real(r8) :: surface_tension_water ! (J/m^2), Arah and Vinten 1995
    real(r8) :: rij_kro_a             !  Arah and Vinten 1995
    real(r8) :: rij_kro_alpha         !  Arah and Vinten 1995
    real(r8) :: rij_kro_beta          !  Arah and Vinten 1995
    real(r8) :: rij_kro_gamma         !  Arah and Vinten 1995
    real(r8) :: rij_kro_delta         !  Arah and Vinten 1995
    real(r8) :: rho_w  = 1.e3_r8                   ! (kg/m3)
    real(r8) :: r_max
    real(r8) :: r_min(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: ratio_diffusivity_water_gas(bounds%begc:bounds%endc,1:nlevdecomp)
    real(r8) :: om_frac
    real(r8) :: anaerobic_frac_sat, r_psi_sat, r_min_sat ! scalar values in sat portion for averaging
    real(r8) :: organic_max              ! organic matter content (kg/m3) where
                                         ! soil is assumed to act like peat
    character(len=32) :: subname='nitrif_denitrif' ! subroutine name
    !-----------------------------------------------------------------------

   associate(& 
   phr_vr                              =>    ccf%phr_vr                                  , & ! Input:  [real(r8) (:,:)]  potential hr (not N-limited)                    
   w_scalar                            =>    ccf%w_scalar                                , & ! Input:  [real(r8) (:,:)]  soil water scalar for decomp                    
   t_scalar                            =>    ccf%t_scalar                                , & ! Input:  [real(r8) (:,:)]  temperature scalar for decomp                   
   h2osoi_vol                          =>    cws%h2osoi_vol                              , & ! Input:  [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
   h2osoi_liq                          =>    cws%h2osoi_liq                              , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
   watsat                              =>    cps%watsat                                  , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
   t_soisno                            =>    ces%t_soisno                                , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   smin_nh4_vr                         =>    cns%smin_nh4_vr                             , & ! Input:  [real(r8) (:,:)]  (gN/m3) soil mineral NH4 pool                   
   smin_no3_vr                         =>    cns%smin_no3_vr                             , & ! Input:  [real(r8) (:,:)]  (gN/m3) soil mineral NO3 pool                   
   bd                                  =>    cps%bd                                      , & ! Input:  [real(r8) (:,:)]  bulk density of dry soil material [kg/m3]       
   dz                                  =>    cps%dz                                      , & ! Input:  [real(r8) (:,:)]  layer thickness (m)  (-nlevsno+1:nlevgrnd)      
   watfc                               =>    cps%watfc                                   , & ! Input:  [real(r8) (:,:)]  volumetric soil water at field capacity (nlevsoi)
   bsw                                 =>    cps%bsw                                     , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)             
   soilpsi                             =>    cps%soilpsi                                 , & ! Input:  [real(r8) (:,:)]  soil water potential in each soil layer (MPa)   
   o2_decomp_depth_unsat               =>    cch4%o2_decomp_depth_unsat                  , & ! Input:  [real(r8) (:,:)]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   conc_o2_unsat                       =>    cch4%conc_o2_unsat                          , & ! Input:  [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   o2_decomp_depth_sat                 =>    cch4%o2_decomp_depth_sat                    , & ! Input:  [real(r8) (:,:)]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
   conc_o2_sat                         =>    cch4%conc_o2_sat                            , & ! Input:  [real(r8) (:,:)]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
   finundated                          =>    cws%finundated                              , & ! Input:  [real(r8) (:)]  fractional inundated area in soil column (excluding dedicated wetland columns)
   sucsat                              =>    cps%sucsat                                  , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   r_psi                               =>    cnf%r_psi                                   , & ! Input:  [real(r8) (:,:)]                                                  
   anaerobic_frac                      =>    cnf%anaerobic_frac                          , & ! Input:  [real(r8) (:,:)]                                                  
   ! ! subsets of the n flux calcs (for diagnostic/debugging purposes)
   smin_no3_massdens_vr                =>    cnf%smin_no3_massdens_vr                    , & ! Input:  [real(r8) (:,:)]  (ugN / g soil) soil nitrate concentration       
   k_nitr_t_vr                         =>    cnf%k_nitr_t_vr                             , & ! Input:  [real(r8) (:,:)]                                                  
   k_nitr_ph_vr                        =>    cnf%k_nitr_ph_vr                            , & ! Input:  [real(r8) (:,:)]                                                  
   k_nitr_h2o_vr                       =>    cnf%k_nitr_h2o_vr                           , & ! Input:  [real(r8) (:,:)]                                                  
   k_nitr_vr                           =>    cnf%k_nitr_vr                               , & ! Input:  [real(r8) (:,:)]                                                  
   wfps_vr                             =>    cnf%wfps_vr                                 , & ! Input:  [real(r8) (:,:)]                                                  
   fmax_denit_carbonsubstrate_vr       =>    cnf%fmax_denit_carbonsubstrate_vr           , & ! Input:  [real(r8) (:,:)]                                                  
   fmax_denit_nitrate_vr               =>    cnf%fmax_denit_nitrate_vr                   , & ! Input:  [real(r8) (:,:)]                                                  
   f_denit_base_vr                     =>    cnf%f_denit_base_vr                         , & ! Input:  [real(r8) (:,:)]                                                  
   diffus                              =>    cnf%diffus                                  , & ! Input:  [real(r8) (:,:)] diffusivity (unitless fraction of total diffusivity)
   ratio_k1                            =>    cnf%ratio_k1                                , & ! Input:  [real(r8) (:,:)]                                                  
   ratio_no3_co2                       =>    cnf%ratio_no3_co2                           , & ! Input:  [real(r8) (:,:)]                                                  
   soil_co2_prod                       =>    cnf%soil_co2_prod                           , & ! Input:  [real(r8) (:,:)]  (ug C / g soil / day)                           
   fr_WFPS                             =>    cnf%fr_WFPS                                 , & ! Input:  [real(r8) (:,:)]                                                  
   soil_bulkdensity                    =>    cnf%soil_bulkdensity                        , & ! Input:  [real(r8) (:,:)]  (kg soil / m3) bulk density of soil (including water)
   cellorg                             =>    cps%cellorg                                 , & ! Input:  [real(r8) (:,:)]  column 3D org (kg/m3 organic matter) (nlevgrnd) 
   pot_f_nit_vr                        =>    cnf%pot_f_nit_vr                            , & ! Input:  [real(r8) (:,:)]  (gN/m3/s) potential soil nitrification flux     
   pot_f_denit_vr                      =>    cnf%pot_f_denit_vr                          , & ! Input:  [real(r8) (:,:)]  (gN/m3/s) potential soil denitrification flux   
   n2_n2o_ratio_denit_vr               =>    cnf%n2_n2o_ratio_denit_vr                     & ! Input:  [real(r8) (:,:)]  ratio of N2 to N2O production by denitrification [gN/gN]
   )

   ! Set maximum nitrification rate constant 
   k_nitr_max =  0.1_r8 / secspday   ! [1/sec] 10%/day  Parton et al., 2001 
   ! Todo:  SPM - the explicit divide gives different results than when that
   ! value is placed in the parameters netcdf file.  To get bfb, keep the 
   ! divide in source.
   !k_nitr_max = CNNitrifDenitrifParamsInst%k_nitr_max

   surface_tension_water = CNNitrifDenitrifParamsInst%surface_tension_water

   ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
   rij_kro_a     = CNNitrifDenitrifParamsInst%rij_kro_a
   rij_kro_alpha = CNNitrifDenitrifParamsInst%rij_kro_alpha
   rij_kro_beta  = CNNitrifDenitrifParamsInst%rij_kro_beta
   rij_kro_gamma = CNNitrifDenitrifParamsInst%rij_kro_gamma
   rij_kro_delta = CNNitrifDenitrifParamsInst%rij_kro_delta

   organic_max = CNParamsShareInst%organic_max

   pH(bounds%begc:bounds%endc) = 6.5  !!! set all soils with the same pH as placeholder here
   co2diff_con(1) =   0.1325_r8
   co2diff_con(2) =   0.0009_r8

   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
                  
         !---------------- calculate soil anoxia state
         ! calculate gas diffusivity of soil at field capacity here
         ! use expression from methane code, but neglect OM for now
         f_a = 1._r8 - watfc(c,j) / watsat(c,j)
         eps =  watsat(c,j)-watfc(c,j) ! Air-filled fraction of total soil volume

         ! use diffusivity calculation including peat
         if (use_lch4) then

            if (organic_max > 0._r8) then
               om_frac = min(cellorg(c,j)/organic_max, 1._r8)
               ! Use first power, not square as in iniTimeConst
            else
               om_frac = 1._r8
            end if
            diffus (c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j)) * 1.e-4_r8 * &
                 (om_frac * f_a**(10._r8/3._r8) / watsat(c,j)**2 + &
                 (1._r8-om_frac) * eps**2 * f_a**(3._r8 / bsw(c,j)) ) 

            ! calculate anoxic fraction of soils
            ! use rijtema and kroess model after Riley et al., 2000
            ! caclulated r_psi as a function of psi
            r_min(c,j) = 2 * surface_tension_water / (rho_w * grav * abs(soilpsi(c,j)))
            r_max = 2 * surface_tension_water / (rho_w * grav * 0.1_r8)
            r_psi(c,j) = sqrt(r_min(c,j) * r_max)
            ratio_diffusivity_water_gas(c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j) ) * 1.e-4_r8 / &
                 ((d_con_w(2,1) + d_con_w(2,2)*t_soisno(c,j) + d_con_w(2,3)*t_soisno(c,j)**2) * 1.e-9_r8)

            if (o2_decomp_depth_unsat(c,j) .ne. spval .and. conc_o2_unsat(c,j) .ne. spval .and.  & 
                 o2_decomp_depth_unsat(c,j) .gt. 0._r8) then
               anaerobic_frac(c,j) = exp(-rij_kro_a * r_psi(c,j)**(-rij_kro_alpha) * &
                    o2_decomp_depth_unsat(c,j)**(-rij_kro_beta) * &
                    conc_o2_unsat(c,j)**rij_kro_gamma * (h2osoi_vol(c,j) + ratio_diffusivity_water_gas(c,j) * &
                    watsat(c,j))**rij_kro_delta)
            else
               anaerobic_frac(c,j) = 0._r8
            endif

            if (anoxia_wtsat) then ! Average saturated fraction values into anaerobic_frac(c,j).
               r_min_sat = 2._r8 * surface_tension_water / (rho_w * grav * abs(grav * 1.e-6_r8 * sucsat(c,j)))
               r_psi_sat = sqrt(r_min_sat * r_max)
               if (o2_decomp_depth_sat(c,j) .ne. spval .and. conc_o2_sat(c,j) .ne. spval .and. &
                    o2_decomp_depth_sat(c,j) .gt. 0._r8) then
                  anaerobic_frac_sat = exp(-rij_kro_a * r_psi_sat**(-rij_kro_alpha) * &
                       o2_decomp_depth_sat(c,j)**(-rij_kro_beta) * &
                       conc_o2_sat(c,j)**rij_kro_gamma * (watsat(c,j) + ratio_diffusivity_water_gas(c,j) * &
                       watsat(c,j))**rij_kro_delta)
               else
                  anaerobic_frac_sat = 0._r8
               endif
               anaerobic_frac(c,j) = (1._r8 - finundated(c))*anaerobic_frac(c,j) + finundated(c)*anaerobic_frac_sat
            end if

         else
            ! NITRIF_DENITRIF requires Methane model to be active, 
            ! otherwise diffusivity will be zeroed out here. EBK CDK 10/18/2011
            anaerobic_frac(c,j) = 0._r8
            diffus (c,j) = 0._r8
            !call endrun(msg=' ERROR: NITRIF_DENITRIF requires Methane model to be active'//errMsg(__FILE__, __LINE__) )
         end if

         
         !---------------- nitrification
         ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

         ! assume nitrification temp function equal to the HR scalar
         k_nitr_t_vr(c,j) = min(t_scalar(c,j), 1._r8)
         
         ! ph function from Parton et al., (2001, 1996)
         k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c)))/rpi
         
         ! moisture function-- assume the same moisture function as limits heterotrophic respiration
         ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
         k_nitr_h2o_vr(c,j) = w_scalar(c,j)
         
         ! nitrification constant is a set scalar * temp, moisture, and ph scalars
         k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

         ! first-order decay of ammonium pool with scalar defined above
         pot_f_nit_vr(c,j) = max(smin_nh4_vr(c,j) * k_nitr_vr(c,j), 0._r8)

         ! limit to oxic fraction of soils
         pot_f_nit_vr(c,j)  = pot_f_nit_vr(c,j) * (1._r8 - anaerobic_frac(c,j))

         ! limit to non-frozen soil layers
         if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
            pot_f_nit_vr(c,j) = 0._r8
         endif


         !---------------- denitrification
         ! first some input variables an unit conversions
         soil_hr_vr(c,j) = phr_vr(c,j)

         ! CENTURY papers give denitrification in units of per gram soil; need to convert from volumetric to mass-based units here
         soil_bulkdensity(c,j) = bd(c,j) + h2osoi_liq(c,j)/dz(c,j)         

         g_per_m3__to__ug_per_gsoil = 1.e3_r8 / soil_bulkdensity(c,j)
         
         g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * secspday
         
         smin_no3_massdens_vr(c,j) = max(smin_no3_vr(c,j), 0._r8) * g_per_m3__to__ug_per_gsoil

         soil_co2_prod(c,j) = (soil_hr_vr(c,j) * (g_per_m3_sec__to__ug_per_gsoil_day))

         !! maximum potential denitrification rates based on heterotrophic respiration rates or nitrate concentrations, 
         !! from (del Grosso et al., 2000)
         fmax_denit_carbonsubstrate_vr(c,j) = (0.1_r8 * (soil_co2_prod(c,j)**1.3_r8)) &
              / g_per_m3_sec__to__ug_per_gsoil_day
         !  
         fmax_denit_nitrate_vr(c,j) = (1.15_r8 * smin_no3_massdens_vr(c,j)**0.57_r8)  &
              / g_per_m3_sec__to__ug_per_gsoil_day
         
         ! find limiting denitrification rate
         f_denit_base_vr(c,j) = max(min(fmax_denit_carbonsubstrate_vr(c,j), fmax_denit_nitrate_vr(c,j)),0._r8) 

         ! limit to non-frozen soil layers
         if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif ) then
            f_denit_base_vr(c,j) = 0._r8
         endif
         
         ! limit to anoxic fraction of soils
         pot_f_denit_vr(c,j) = f_denit_base_vr(c,j) * anaerobic_frac(c,j)

         ! now calculate the ratio of N2O to N2 from denitrifictaion, following Del Grosso et al., 2000
         ! diffusivity constant (figure 6b)
         ratio_k1(c,j) = max(1.7_r8, 38.4_r8 - 350._r8 * diffus(c,j))
         
         ! ratio function (figure 7c)
         if ( soil_co2_prod(c,j) .gt. 0 ) then
            ratio_no3_co2(c,j) = smin_no3_massdens_vr(c,j) / soil_co2_prod(c,j)
         else
            ! fucntion saturates at large no3/co2 ratios, so set as some nominally large number
            ratio_no3_co2(c,j) = 100._r8
         endif

         ! total water limitation function (Del Grosso et al., 2000, figure 7a)
         wfps_vr(c,j) = max(min(h2osoi_vol(c,j)/watsat(c, j), 1._r8), 0._r8) * 100._r8
         fr_WFPS(c,j) = max(0.1_r8, 0.015_r8 * wfps_vr(c,j) - 0.32_r8)
         if (use_lch4) then
            if (anoxia_wtsat) then
               fr_WFPS(c,j) = fr_WFPS(c,j)*(1._r8 - finundated(c)) + finundated(c)*1.18_r8
            end if
         end if

         ! final ratio expression 
         n2_n2o_ratio_denit_vr(c,j) = max(0.16*ratio_k1(c,j), ratio_k1(c,j)*exp(-0.8 * ratio_no3_co2(c,j))) * fr_WFPS(c,j)

      end do
            
   end do

 end associate
 end subroutine nitrif_denitrif

end module CNNitrifDenitrifMod
