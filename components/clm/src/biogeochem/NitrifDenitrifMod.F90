module CNNitrifDenitrifMod

  !-----------------------------------------------------------------------
  ! Calculate nitrification and denitrification rates
  !
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_const_mod       , only : SHR_CONST_TKFRZ
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_const_mod       , only : SHR_CONST_TKFRZ
  use clm_varpar          , only : nlevgrnd,nlevdecomp
  use clm_varcon          , only : rpi, denh2o, dzsoi, zisoi, grav
  use clm_varcon          , only : d_con_g, d_con_w, spval, secspday
  use clm_varctl          , only : use_lch4, iulog
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use SoilStatetype       , only : soilstate_type
  use WaterStateType      , only : waterstate_type
  use TemperatureType     , only : temperature_type
  use CNCarbonfluxType    , only : carbonflux_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CH4Mod              , only : ch4_type
  use ColumnType          , only : col_pp                
  !
  implicit none
  save
  private
  !
  public :: nitrif_denitrif
  public :: readNitrifDenitrifParams
  !
  type, private :: NitrifDenitrifParamsType
   real(r8) :: k_nitr_max               !  maximum nitrification rate constant (1/s)
   real(r8) :: surface_tension_water    !  surface tension of water(J/m^2), Arah an and Vinten 1995
   real(r8) :: rij_kro_a                !  Arah and Vinten 1995)
   real(r8) :: rij_kro_alpha            !  parameter to calculate anoxic fraction of soil  (Arah and Vinten 1995)
   real(r8) :: rij_kro_beta             !  (Arah and Vinten 1995)
   real(r8) :: rij_kro_gamma            !  (Arah and Vinten 1995)
   real(r8) :: rij_kro_delta            !  (Arah and Vinten 1995)
  end type NitrifDenitrifParamsType

  type(NitrifDenitrifParamsType),private ::  NitrifDenitrifParamsInst

  logical, public :: no_frozen_nitrif_denitrif = .false.  ! stop nitrification and denitrification in frozen soils
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------  
  subroutine readNitrifDenitrifParams ( ncid )
    !
    use ncdio_pio    , only: file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'NitrifDenitrifParamsType'
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
    NitrifDenitrifParamsInst%k_nitr_max=tempr

    tString='surface_tension_water'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%surface_tension_water=tempr

    tString='rij_kro_a'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%rij_kro_a=tempr

    tString='rij_kro_alpha'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%rij_kro_alpha=tempr

    tString='rij_kro_beta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%rij_kro_beta=tempr

    tString='rij_kro_gamma'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%rij_kro_gamma=tempr

    tString='rij_kro_delta'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    NitrifDenitrifParamsInst%rij_kro_delta=tempr

  end subroutine readNitrifDenitrifParams

  !-----------------------------------------------------------------------
  subroutine nitrif_denitrif(bounds, num_soilc, filter_soilc, &
       soilstate_vars, waterstate_vars, temperature_vars, ch4_vars, &
       carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    !  calculate nitrification and denitrification rates
    !
    ! !USES:
    use clm_time_manager  , only : get_curr_date, get_step_size
    use SharedParamsMod , only : anoxia_wtsat,ParamsShareInst
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(carbonflux_type)    , intent(in)    :: carbonflux_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
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
    !debug-- put these type structure for outing to hist files
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

    associate(                                                                                     & 
         watsat                        =>    soilstate_vars%watsat_col                           , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity) (nlevgrnd)
         bd                            =>    soilstate_vars%bd_col                               , & ! Input:  [real(r8) (:,:)  ]  bulk density of dry soil material [kg/m3]       
         watfc                         =>    soilstate_vars%watfc_col                            , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at field capacity (nlevsoi)
         bsw                           =>    soilstate_vars%bsw_col                              , & ! Input:  [real(r8) (:,:)  ]  Clapp and Hornberger "b" (nlevgrnd)             
         cellorg                       =>    soilstate_vars%cellorg_col                          , & ! Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd) 
         sucsat                        =>    soilstate_vars%sucsat_col                           , & ! Input:  [real(r8) (:,:)  ]  minimum soil suction (mm)                       
         soilpsi                       =>    soilstate_vars%soilpsi_col                          , & ! Input:  [real(r8) (:,:)  ]  soil water potential in each soil layer (MPa)   
         
         h2osoi_vol                    =>    waterstate_vars%h2osoi_vol_col                      , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
         h2osoi_liq                    =>    waterstate_vars%h2osoi_liq_col                      , & ! Input:  [real(r8) (:,:)  ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
         
         t_soisno                      =>    temperature_vars%t_soisno_col                       , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         
         o2_decomp_depth_unsat         =>    ch4_vars%o2_decomp_depth_unsat_col                  , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
         conc_o2_unsat                 =>    ch4_vars%conc_o2_unsat_col                          , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
         o2_decomp_depth_sat           =>    ch4_vars%o2_decomp_depth_sat_col                    , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
         conc_o2_sat                   =>    ch4_vars%conc_o2_sat_col                            , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)   
         finundated                    =>    ch4_vars%finundated_col                             , & ! Input:  [real(r8) (:)    ]  fractional inundated area in soil column (excluding dedicated wetland columns)
         
         phr_vr                        =>    carbonflux_vars%phr_vr_col                          , & ! Input:  [real(r8) (:,:)  ]  potential hr (not N-limited)                    
         w_scalar                      =>    carbonflux_vars%w_scalar_col                        , & ! Input:  [real(r8) (:,:)  ]  soil water scalar for decomp                    
         t_scalar                      =>    carbonflux_vars%t_scalar_col                        , & ! Input:  [real(r8) (:,:)  ]  temperature scalar for decomp                   

         smin_nh4_vr                   =>    nitrogenstate_vars%smin_nh4_vr_col                  , & ! Input:  [real(r8) (:,:)  ]  (gN/m3) soil mineral NH4 pool                   
         smin_no3_vr                   =>    nitrogenstate_vars%smin_no3_vr_col                  , & ! Input:  [real(r8) (:,:)  ]  (gN/m3) soil mineral NO3 pool                   

         r_psi                         =>    nitrogenflux_vars%r_psi_col                         , & ! Output:  [real(r8) (:,:)  ]                                                  
         anaerobic_frac                =>    nitrogenflux_vars%anaerobic_frac_col                , & ! Output:  [real(r8) (:,:)  ]                                                  
         ! ! subsets of the n flux calcs (for diagnostic/debugging purposes)
         smin_no3_massdens_vr          =>    nitrogenflux_vars%smin_no3_massdens_vr_col          , & ! Output:  [real(r8) (:,:) ]  (ugN / g soil) soil nitrate concentration       
         k_nitr_t_vr                   =>    nitrogenflux_vars%k_nitr_t_vr_col                   , & ! Output:  [real(r8) (:,:) ]                                                  
         k_nitr_ph_vr                  =>    nitrogenflux_vars%k_nitr_ph_vr_col                  , & ! Output:  [real(r8) (:,:) ]                                                  
         k_nitr_h2o_vr                 =>    nitrogenflux_vars%k_nitr_h2o_vr_col                 , & ! Output:  [real(r8) (:,:) ]                                                  
         k_nitr_vr                     =>    nitrogenflux_vars%k_nitr_vr_col                     , & ! Output:  [real(r8) (:,:) ]                                                  
         wfps_vr                       =>    nitrogenflux_vars%wfps_vr_col                       , & ! Output:  [real(r8) (:,:) ]                                                  
         fmax_denit_carbonsubstrate_vr =>    nitrogenflux_vars%fmax_denit_carbonsubstrate_vr_col , & ! Output:  [real(r8) (:,:) ]                                                  
         fmax_denit_nitrate_vr         =>    nitrogenflux_vars%fmax_denit_nitrate_vr_col         , & ! Output:  [real(r8) (:,:) ]                                                  
         f_denit_base_vr               =>    nitrogenflux_vars%f_denit_base_vr_col               , & ! Output:  [real(r8) (:,:) ]                                                  
         diffus                        =>    nitrogenflux_vars%diffus_col                        , & ! Output:  [real(r8) (:,:) ] diffusivity (unitless fraction of total diffusivity)
         ratio_k1                      =>    nitrogenflux_vars%ratio_k1_col                      , & ! Output:  [real(r8) (:,:) ]                                                  
         ratio_no3_co2                 =>    nitrogenflux_vars%ratio_no3_co2_col                 , & ! Output:  [real(r8) (:,:) ]                                                  
         soil_co2_prod                 =>    nitrogenflux_vars%soil_co2_prod_col                 , & ! Output:  [real(r8) (:,:) ]  (ug C / g soil / day)                           
         fr_WFPS                       =>    nitrogenflux_vars%fr_WFPS_col                       , & ! Output:  [real(r8) (:,:) ]                                                  
         soil_bulkdensity              =>    nitrogenflux_vars%soil_bulkdensity_col              , & ! Output:  [real(r8) (:,:) ]  (kg soil / m3) bulk density of soil (including water)
         pot_f_nit_vr                  =>    nitrogenflux_vars%pot_f_nit_vr_col                  , & ! Output:  [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux     
         pot_f_denit_vr                =>    nitrogenflux_vars%pot_f_denit_vr_col                , & ! Output:  [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux   
         n2_n2o_ratio_denit_vr         =>    nitrogenflux_vars%n2_n2o_ratio_denit_vr_col           & ! Output:  [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         )

      ! Set maximum nitrification rate constant 
      k_nitr_max =  0.1_r8 / secspday   ! [1/sec] 10%/day  Parton et al., 2001 

      ! Todo:  FIX(SPM,032414) - the explicit divide gives different results than when that
      ! value is placed in the parameters netcdf file.  To get bfb, keep the 
      ! divide in source.
      !k_nitr_max = NitrifDenitrifParamsInst%k_nitr_max

      surface_tension_water = NitrifDenitrifParamsInst%surface_tension_water

      ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
      rij_kro_a     = NitrifDenitrifParamsInst%rij_kro_a
      rij_kro_alpha = NitrifDenitrifParamsInst%rij_kro_alpha
      rij_kro_beta  = NitrifDenitrifParamsInst%rij_kro_beta
      rij_kro_gamma = NitrifDenitrifParamsInst%rij_kro_gamma
      rij_kro_delta = NitrifDenitrifParamsInst%rij_kro_delta

      organic_max = ParamsShareInst%organic_max

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

               if (o2_decomp_depth_unsat(c,j) /= spval .and. conc_o2_unsat(c,j) /= spval .and.  & 
                    o2_decomp_depth_unsat(c,j) > 0._r8) then
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
                  if (o2_decomp_depth_sat(c,j) /= spval .and. conc_o2_sat(c,j) /= spval .and. &
                       o2_decomp_depth_sat(c,j) > 0._r8) then
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
            soil_bulkdensity(c,j) = bd(c,j) + h2osoi_liq(c,j)/col_pp%dz(c,j)         

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
            if ( soil_co2_prod(c,j) > 0 ) then
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
