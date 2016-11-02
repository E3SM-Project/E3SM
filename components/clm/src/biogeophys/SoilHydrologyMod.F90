module SoilHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil hydrology
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varctl        , only : iulog, use_vichydro
  use clm_varcon        , only : e_ice, denh2o, denice, rpi
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use TemperatureType   , only : temperature_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft     
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRunoff        ! Calculate surface runoff
  public :: Infiltration         ! Calculate infiltration into surface soil layer
  public :: WaterTable           ! Calculate water table before imposing drainage
  public :: Drainage             ! Calculate subsurface drainage
  public :: DrainageVSFM         ! Calculate subsurface drainage for VSFM
  public :: CLMVICMap
  ! !PRIVATE DATA MEMBERS:
  ! Get these from SoilWaterMovementMod eventually
  integer, parameter :: zengdecker2009 = 0
  integer, parameter :: vsfm = 1
  integer :: soilroot_water_method=0     !0: use the Zeng and deck method, this will be readin from namelist in the future
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SurfaceRunoff (bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, waterflux_vars, &
       waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff
    !
    ! !USES:
    use clm_varcon      , only : denice, denh2o, wimp, pondmx_urban
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_imperv, icol_road_perv
    use clm_varpar      , only : nlevsoi, nlevgrnd, maxpatch_pft
    use clm_time_manager, only : get_step_size
    use clm_varpar      , only : nlayer, nlayert
     use clm_varctl     , only : do_varsoil
    use abortutils      , only : endrun
!    use SoilWaterMovementMod, only : soilroot_water_method, zengdecker_2009, vsfm
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds               
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,g,l,i                               !indices
    integer  :: nlevbed                                    !# levels to bedrock
    real(r8) :: dtime                                      !land model time step (sec)
    real(r8) :: xs(bounds%begc:bounds%endc)                !excess soil water above urban ponding limit
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevgrnd) !partial volume of ice lens in layer
    real(r8) :: fff(bounds%begc:bounds%endc)               !decay factor (m-1)
    real(r8) :: s1                                         !variable to calculate qinmax
    real(r8) :: su                                         !variable to calculate qinmax
    real(r8) :: v                                          !variable to calculate qinmax
    real(r8) :: qinmax                                     !maximum infiltration capacity (mm/s)
    real(r8) :: A(bounds%begc:bounds%endc)                 !fraction of the saturated area
    real(r8) :: ex(bounds%begc:bounds%endc)                !temporary variable (exponent)
    real(r8) :: top_moist(bounds%begc:bounds%endc)         !temporary, soil moisture in top VIC layers
    real(r8) :: top_max_moist(bounds%begc:bounds%endc)     !temporary, maximum soil moisture in top VIC layers
    real(r8) :: top_ice(bounds%begc:bounds%endc)           !temporary, ice len in top VIC layers
    character(len=32) :: subname = 'SurfaceRunoff'         !subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        & 
         snl              =>    col%snl                             , & ! Input:  [integer  (:)   ]  minus number of snow layers                        
         dz               =>    col%dz                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         nlev2bed         =>    col%nlev2bed                        , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     

         sucsat           =>    soilstate_vars%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         watsat           =>    soilstate_vars%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         wtfact           =>    soilstate_vars%wtfact_col           , & ! Input:  [real(r8) (:)   ]  maximum saturated fraction for a gridcell         
         hksat            =>    soilstate_vars%hksat_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw              =>    soilstate_vars%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        

         h2osoi_ice       =>    waterstate_vars%h2osoi_ice_col      , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq       =>    waterstate_vars%h2osoi_liq_col      , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            

         qflx_snow_h2osfc =>    waterflux_vars%qflx_snow_h2osfc_col , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)              
         qflx_floodc      =>    waterflux_vars%qflx_floodc_col      , & ! Input:  [real(r8) (:)   ]  column flux of flood water from RTM               
         qflx_evap_grnd   =>    waterflux_vars%qflx_evap_grnd_col   , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]    
         qflx_top_soil    =>    waterflux_vars%qflx_top_soil_col    , & ! Output: [real(r8) (:)   ]  net water input into soil from top (mm/s)         
         qflx_surf        =>    waterflux_vars%qflx_surf_col        , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)                        

         zwt              =>    soilhydrology_vars%zwt_col          , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         max_moist        =>    soilhydrology_vars%max_moist_col    , & ! Input:  [real(r8) (:,:) ]  maximum soil moisture (ice + liq, mm)            
         frost_table      =>    soilhydrology_vars%frost_table_col  , & ! Input:  [real(r8) (:)   ]  frost table depth (m)                             
         zwt_perched      =>    soilhydrology_vars%zwt_perched_col  , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)                     
         b_infil          =>    soilhydrology_vars%b_infil_col      , & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter                       
         moist            =>    soilhydrology_vars%moist_col        , & ! Input:  [real(r8) (:,:) ]  soil moisture in each VIC layers (liq, mm)       
         hkdepth          =>    soilhydrology_vars%hkdepth_col      , & ! Input:  [real(r8) (:)   ]  decay factor (m)                                  
         origflag         =>    soilhydrology_vars%origflag         , & ! Input:  logical
         fcov             =>    soilhydrology_vars%fcov_col         , & ! Output: [real(r8) (:)   ]  fractional impermeable area                       
         fsat             =>    soilhydrology_vars%fsat_col         , & ! Output: [real(r8) (:)   ]  fractional area with water table at surface       
         fracice          =>    soilhydrology_vars%fracice_col      , & ! Output: [real(r8) (:,:) ]  fractional impermeability (-)                    
         icefrac          =>    soilhydrology_vars%icefrac_col      , & ! Output: [real(r8) (:,:) ]                                                  
         ice              =>    soilhydrology_vars%ice_col          , & ! Output: [real(r8) (:,:) ]  ice len in each VIC layers(ice, mm)              
         max_infil        =>    soilhydrology_vars%max_infil_col    , & ! Output: [real(r8) (:)   ]  maximum infiltration capacity in VIC (mm)          
         i_0              =>    soilhydrology_vars%i_0_col            & ! Output: [real(r8) (:)   ]  column average soil moisture in top VIC layers (mm)
         )

      ! Get time step

      dtime = get_step_size()

      do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
        if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
        else
            nlevbed = nlevsoi
        end if
         do j = 1,nlevbed

            ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
            ! fractional impermeability

            vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
            if (origflag == 1) then
               icefrac(c,j) = min(1._r8,h2osoi_ice(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j)))
            else
               icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
            endif

            fracice(c,j) = max(0._r8,exp(-3._r8*(1._r8-icefrac(c,j)))- exp(-3._r8))/(1.0_r8-exp(-3._r8))
         end do
      end do

      ! Saturated fraction

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         fff(c) = 0.5_r8
         if (use_vichydro) then 
            top_moist(c) = 0._r8
            top_ice(c) = 0._r8
            top_max_moist(c) = 0._r8
            do j = 1, nlayer - 1
               top_ice(c) = top_ice(c) + ice(c,j)
               top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
               top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
            end do
            if(top_moist(c)> top_max_moist(c)) top_moist(c)= top_max_moist(c)
            top_ice(c)     = max(0._r8,top_ice(c))
            max_infil(c)   = (1._r8+b_infil(c)) * top_max_moist(c)
            ex(c)          = b_infil(c) / (1._r8 + b_infil(c))
            A(c)           = 1._r8 - (1._r8 - top_moist(c) / top_max_moist(c))**ex(c)
            i_0(c)         = max_infil(c) * (1._r8 - (1._r8 - A(c))**(1._r8/b_infil(c)))
            fsat(c)        = A(c)  !for output
         else
            fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt(c))
         end if

         ! use perched water table to determine fsat (if present)
         if ( frost_table(c) > zwt(c)) then 
            if (use_vichydro) then
               fsat(c) =  A(c)
            else
               fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt(c))
            end if
         else
            if ( frost_table(c) > zwt_perched(c)) then 
               fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt_perched(c))!*( frost_table(c) - zwt_perched(c))/4.0
            endif
         endif
         if (origflag == 1) then
            if (use_vichydro) then
               call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
            else
               fcov(c) = (1._r8 - fracice(c,1)) * fsat(c) + fracice(c,1)
            end if
         else
            fcov(c) = fsat(c)
         endif
      end do

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         ! assume qinmax large relative to qflx_top_soil in control
         if (origflag == 1) then
            qflx_surf(c) =  fcov(c) * qflx_top_soil(c)
         else
            ! only send fast runoff directly to streams
            qflx_surf(c) =   fsat(c) * qflx_top_soil(c)
         endif
      end do

      ! Determine water in excess of ponding limit for urban roof and impervious road.
      ! Excess goes to surface runoff. No surface runoff for sunwall and shadewall.

      do fc = 1, num_urbanc
         c = filter_urbanc(fc)
         if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then

            ! If there are snow layers then all qflx_top_soil goes to surface runoff
            if (snl(c) < 0) then
               qflx_surf(c) = max(0._r8,qflx_top_soil(c))
            else
               xs(c) = max(0._r8, &
                    h2osoi_liq(c,1)/dtime + qflx_top_soil(c) - qflx_evap_grnd(c) - &
                    pondmx_urban/dtime)
               if (xs(c) > 0.) then
                  h2osoi_liq(c,1) = pondmx_urban
               else
                  h2osoi_liq(c,1) = max(0._r8,h2osoi_liq(c,1)+ &
                       (qflx_top_soil(c)-qflx_evap_grnd(c))*dtime)
               end if
               qflx_surf(c) = xs(c)
            end if
         else if (col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall) then
            qflx_surf(c) = 0._r8
         end if
         ! send flood water flux to runoff for all urban columns
         qflx_surf(c) = qflx_surf(c)  + qflx_floodc(c)

      end do

      ! remove stormflow and snow on h2osfc from qflx_top_soil
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         ! add flood water flux to qflx_top_soil
         qflx_top_soil(c) = qflx_top_soil(c) + qflx_snow_h2osfc(c) + qflx_floodc(c)
      end do

    end associate

   end subroutine SurfaceRunoff

   !-----------------------------------------------------------------------
   subroutine Infiltration(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        energyflux_vars, soilhydrology_vars, soilstate_vars, temperature_vars, &
        waterflux_vars, waterstate_vars)
     !
     ! !DESCRIPTION:
     ! Calculate infiltration into surface soil layer (minus the evaporation)
     !
     ! !USES:
     use shr_const_mod    , only : shr_const_pi
     use clm_varpar       , only : nlayer, nlayert
     use clm_varpar       , only : nlevsoi, nlevgrnd
     use clm_varctl       , only : do_varsoil
     use clm_varcon       , only : denh2o, denice, roverg, wimp, pc, mu, tfrz
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_sunwall, icol_shadewall, icol_road_perv
     use landunit_varcon  , only : istsoil, istcrop
     use clm_time_manager , only : get_step_size
!     use SoilWaterMovementMod, only : soilroot_water_method, zengdecker_2009, vsfm
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     type(energyflux_type)    , intent(in)    :: energyflux_vars 
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(soilstate_type)     , intent(inout) :: soilstate_vars
     type(temperature_type)   , intent(in)    :: temperature_vars
     type(waterstate_type)    , intent(inout) :: waterstate_vars
     type(waterflux_type)     , intent(inout) :: waterflux_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,l,fc                                   ! indices
     integer  :: nlevbed                                    !# levels to bedrock
     real(r8) :: dtime                                      ! land model time step (sec)
     real(r8) :: s1,su,v                                    ! variable to calculate qinmax
     real(r8) :: qinmax                                     ! maximum infiltration capacity (mm/s)
     real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevgrnd) ! partial volume of ice lens in layer
     real(r8) :: alpha_evap(bounds%begc:bounds%endc)        ! fraction of total evap from h2osfc
     real(r8) :: qflx_evap(bounds%begc:bounds%endc)         ! local evaporation array
     real(r8) :: qflx_h2osfc_drain(bounds%begc:bounds%endc) ! bottom drainage from h2osfc
     real(r8) :: qflx_in_h2osfc(bounds%begc:bounds%endc)    ! surface input to h2osfc
     real(r8) :: qflx_in_soil(bounds%begc:bounds%endc)      ! surface input to soil
     real(r8) :: qflx_infl_excess(bounds%begc:bounds%endc)  ! infiltration excess runoff -> h2osfc
     real(r8) :: frac_infclust                              ! fraction of submerged area that is connected
     real(r8) :: fsno                                       ! copy of frac_sno
     real(r8) :: k_wet                                      ! linear reservoir coefficient for h2osfc
     real(r8) :: fac                                        ! soil wetness of surface layer
     real(r8) :: psit                                       ! negative potential of soil
     real(r8) :: hr                                         ! relative humidity
     real(r8) :: wx                                         ! partial volume of ice and water of surface layer
     real(r8) :: z_avg
     real(r8) :: rho_avg
     real(r8) :: fmelt
     real(r8) :: f_sno
     real(r8) :: imped
     real(r8) :: d
     real(r8) :: h2osoi_vol                                 
     real(r8) :: basis                                      ! temporary, variable soil moisture holding capacity 
     ! in top VIC layers for runoff calculation
     real(r8) :: rsurf_vic                                  ! temp VIC surface runoff
     real(r8) :: top_moist(bounds%begc:bounds%endc)         ! temporary, soil moisture in top VIC layers
     real(r8) :: top_max_moist(bounds%begc:bounds%endc)     ! temporary, maximum soil moisture in top VIC layers
     real(r8) :: top_ice(bounds%begc:bounds%endc)           ! temporary, ice len in top VIC layers
     real(r8) :: top_icefrac                                ! temporary, ice fraction in top VIC layers
     !-----------------------------------------------------------------------

     associate(                                                                & 
          snl                  =>    col%snl                                 , & ! Input:  [integer  (:)   ]  minus number of snow layers                        
          dz                   =>    col%dz                                  , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         nlev2bed         =>    col%nlev2bed                                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     

          t_soisno             =>    temperature_vars%t_soisno_col           , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          frac_h2osfc          =>    waterstate_vars%frac_h2osfc_col         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          frac_sno             =>    waterstate_vars%frac_sno_eff_col        , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
          h2osoi_ice           =>    waterstate_vars%h2osoi_ice_col          , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_liq           =>    waterstate_vars%h2osoi_liq_col          , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osno               =>    waterstate_vars%h2osno_col              , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                               
          snow_depth           =>    waterstate_vars%snow_depth_col          , & ! Input:  [real(r8) (:)   ]  snow height (m)                                   
          h2osfc               =>    waterstate_vars%h2osfc_col              , & ! Output: [real(r8) (:)   ]  surface water (mm)                                

          qflx_ev_soil         =>    waterflux_vars%qflx_ev_soil_col         , & ! Input:  [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]    
          qflx_evap_soi        =>    waterflux_vars%qflx_evap_soi_col        , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]    
          qflx_evap_grnd       =>    waterflux_vars%qflx_evap_grnd_col       , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]    
          qflx_top_soil        =>    waterflux_vars%qflx_top_soil_col        , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)         
          qflx_ev_h2osfc       =>    waterflux_vars%qflx_ev_h2osfc_col       , & ! Input:  [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]                      
          qflx_surf            =>    waterflux_vars%qflx_surf_col            , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
          qflx_h2osfc_surf     =>    waterflux_vars%qflx_h2osfc_surf_col     , & ! Output: [real(r8) (:)   ]  surface water runoff (mm/s)                       
          qflx_infl            =>    waterflux_vars%qflx_infl_col            , & ! Output: [real(r8) (:)   ] infiltration (mm H2O /s)                           
          qflx_gross_infl_soil =>    waterflux_vars%qflx_gross_infl_soil_col , & ! Output: [real(r8) (:)] gross infiltration (mm H2O/s)
          qflx_gross_evap_soil =>    waterflux_vars%qflx_gross_evap_soil_col , & ! Output: [real(r8) (:)] gross evaporation (mm H2O/s)

          smpmin               =>    soilstate_vars%smpmin_col               , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)        
          sucsat               =>    soilstate_vars%sucsat_col               , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat               =>    soilstate_vars%watsat_col               , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
          bsw                  =>    soilstate_vars%bsw_col                  , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat                =>    soilstate_vars%hksat_col                , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          eff_porosity         =>    soilstate_vars%eff_porosity_col         , & ! Output: [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         

          h2osfc_thresh        =>    soilhydrology_vars%h2osfc_thresh_col    , & ! Input:  [real(r8) (:)   ]  level at which h2osfc "percolates"                
          zwt                  =>    soilhydrology_vars%zwt_col              , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched          =>    soilhydrology_vars%zwt_perched_col      , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)                     
          fcov                 =>    soilhydrology_vars%fcov_col             , & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface       
          b_infil              =>    soilhydrology_vars%b_infil_col          , & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter                       
          frost_table          =>    soilhydrology_vars%frost_table_col      , & ! Input:  [real(r8) (:)   ]  frost table depth (m)                             
          fsat                 =>    soilhydrology_vars%fsat_col             , & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface       
          moist                =>    soilhydrology_vars%moist_col            , & ! Input:  [real(r8) (:,:) ]  soil moisture in each VIC layers (liq, mm)       
          max_moist            =>    soilhydrology_vars%max_moist_col        , & ! Input:  [real(r8) (:,:) ]  maximum soil moisture (ice + liq, mm)            
          max_infil            =>    soilhydrology_vars%max_infil_col        , & ! Input:  [real(r8) (:)   ]  maximum infiltration capacity in VIC (mm)          
          ice                  =>    soilhydrology_vars%ice_col              , & ! Input:  [real(r8) (:,:) ]  ice len in each VIC layers(ice, mm)              
          i_0                  =>    soilhydrology_vars%i_0_col              , & ! Input:  [real(r8) (:)   ]  column average soil moisture in top VIC layers (mm)
          h2osfcflag           =>    soilhydrology_vars%h2osfcflag           , & ! Input:  logical
          icefrac              =>    soilhydrology_vars%icefrac_col            & ! Output: [real(r8) (:,:) ]  fraction of ice                                 
              )

       dtime = get_step_size()

       ! Infiltration into surface soil layer (minus the evaporation)
       do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if
          do j = 1,nlevbed
             ! Porosity of soil, partial volume of ice and liquid
             vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             eff_porosity(c,j) = max(0.01_r8,watsat(c,j)-vol_ice(c,j))
             icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          ! partition moisture fluxes between soil and h2osfc       
          if (lun%itype(col%landunit(c)) == istsoil .or. lun%itype(col%landunit(c))==istcrop) then

             ! explicitly use frac_sno=0 if snl=0
             if (snl(c) >= 0) then
                fsno=0._r8
                ! if no snow layers, sublimation is removed from h2osoi_ice in drainage
                qflx_evap(c)=qflx_evap_grnd(c)
             else
                fsno=frac_sno(c)
                qflx_evap(c)=qflx_ev_soil(c)
             endif

             !1. partition surface inputs between soil and h2osfc
             qflx_in_soil(c) = (1._r8 - frac_h2osfc(c)) * (qflx_top_soil(c)  - qflx_surf(c))
             qflx_in_h2osfc(c) = frac_h2osfc(c) * (qflx_top_soil(c)  - qflx_surf(c))          
             qflx_gross_infl_soil(c) = (1._r8 - frac_h2osfc(c)) * (qflx_top_soil(c)  - qflx_surf(c))
             
             !2. remove evaporation (snow treated in SnowHydrology)
             qflx_in_soil(c) = qflx_in_soil(c) - (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             qflx_in_h2osfc(c) =  qflx_in_h2osfc(c)  - frac_h2osfc(c) * qflx_ev_h2osfc(c)

             if (qflx_evap(c)>0._r8) then
                qflx_gross_evap_soil(c) = (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             else
                qflx_gross_evap_soil(c) = 0._r8
                qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)-(1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             endif

             !3. determine maximum infiltration rate
             if (use_vichydro) then
                top_moist(c)= 0._r8
                top_ice(c)=0._r8
                top_max_moist(c)= 0._r8
                do j = 1, nlayer - 1
                   top_ice(c) = top_ice(c) + ice(c,j)
                   top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
                   top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
                end do
                top_icefrac = min(1._r8,top_ice(c)/top_max_moist(c))
                if(qflx_in_soil(c) <= 0._r8) then
                   rsurf_vic = 0._r8
                else if(max_infil(c) <= 0._r8) then
                   rsurf_vic = qflx_in_soil(c)
                else if((i_0(c) + qflx_in_soil(c)*dtime) > max_infil(c)) then             !(Eq.(3a) Wood et al. 1992)
                   rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c))/dtime
                else                                                                      !(Eq.(3b) Wood et al. 1992)
                   basis = 1._r8 - (i_0(c) + qflx_in_soil(c)*dtime)/max_infil(c)
                   rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c)    &
                        + top_max_moist(c) * basis**(1._r8 + b_infil(c)))/dtime
                end if
                rsurf_vic = min(qflx_in_soil(c), rsurf_vic)
                qinmax = (1._r8 - fsat(c)) * 10._r8**(-e_ice*top_icefrac)*(qflx_in_soil(c) - rsurf_vic)
             else
                qinmax=(1._r8 - fsat(c)) * minval(10._r8**(-e_ice*(icefrac(c,1:3)))*hksat(c,1:3))
             end if
             qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - frac_h2osfc(c))*qinmax)

             !4. soil infiltration and h2osfc "run-on"
             qflx_infl(c) = qflx_in_soil(c) - qflx_infl_excess(c)
             qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) + qflx_infl_excess(c)
             qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)- qflx_infl_excess(c)

             !5. surface runoff from h2osfc
             if (h2osfcflag==1) then
                ! calculate runoff from h2osfc  -------------------------------------
                if (frac_h2osfc(c) <= pc) then 
                   frac_infclust=0.0_r8
                else
                   frac_infclust=(frac_h2osfc(c)-pc)**mu
                endif
             endif

             ! limit runoff to value of storage above S(pc)
             if(h2osfc(c) >= h2osfc_thresh(c) .and. h2osfcflag/=0) then
                ! spatially variable k_wet
                k_wet=1.0_r8 * sin((rpi/180.) * col%topo_slope(c))
                qflx_h2osfc_surf(c) = k_wet * frac_infclust * (h2osfc(c) - h2osfc_thresh(c))

                qflx_h2osfc_surf(c)=min(qflx_h2osfc_surf(c),(h2osfc(c) - h2osfc_thresh(c))/dtime)
             else
                qflx_h2osfc_surf(c)= 0._r8
             endif

             ! cutoff lower limit
             if ( qflx_h2osfc_surf(c) < 1.0e-8) qflx_h2osfc_surf(c) = 0._r8 

             ! use this for non-h2osfc code
             if(h2osfcflag==0) then 
                qflx_h2osfc_surf(c)= 0._r8
                ! shift infiltration excess from h2osfc input to surface runoff
                qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) - qflx_infl_excess(c)
                qflx_surf(c)= qflx_surf(c) + qflx_infl_excess(c) 
                qflx_infl_excess(c) = 0._r8
             endif

             qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) - qflx_h2osfc_surf(c) 

             !6. update h2osfc prior to calculating bottom drainage from h2osfc
             h2osfc(c) = h2osfc(c) + qflx_in_h2osfc(c) * dtime

             !--  if all water evaporates, there will be no bottom drainage
             if (h2osfc(c) < 0.0) then
                qflx_infl(c) = qflx_infl(c) + h2osfc(c)/dtime
                qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c) - h2osfc(c)/dtime                
                h2osfc(c) = 0.0
                qflx_h2osfc_drain(c)= 0._r8
             else
                qflx_h2osfc_drain(c)=min(frac_h2osfc(c)*qinmax,h2osfc(c)/dtime)
             endif

             if(h2osfcflag==0) then 
                qflx_h2osfc_drain(c)= max(0._r8,h2osfc(c)/dtime) !ensure no h2osfc
             endif

             !7. remove drainage from h2osfc and add to qflx_infl
             h2osfc(c) = h2osfc(c) - qflx_h2osfc_drain(c) * dtime
             qflx_infl(c) = qflx_infl(c) + qflx_h2osfc_drain(c)
             qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) + qflx_h2osfc_drain(c)             
          else
             ! non-vegetated landunits (i.e. urban) use original CLM4 code
             if (snl(c) >= 0) then
                ! when no snow present, sublimation is removed in Drainage
                qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) - qflx_evap_grnd(c)
                qflx_gross_infl_soil(c) = qflx_top_soil(c) - qflx_surf(c)
                qflx_gross_evap_soil(c) = qflx_evap_grnd(c)                
             else
                qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) &
                     - (1.0_r8 - frac_sno(c)) * qflx_ev_soil(c)
                qflx_gross_infl_soil(c) = qflx_top_soil(c) - qflx_surf(c)
                qflx_gross_evap_soil(c) = (1.0_r8 - frac_sno(c)) * qflx_ev_soil(c)                     
             end if
             qflx_h2osfc_surf(c) = 0._r8
          endif

       enddo

       ! No infiltration for impervious urban surfaces

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) /= icol_road_perv) then
             qflx_infl(c) = 0._r8
          end if
       end do
    
    end associate

   end subroutine Infiltration

   !-----------------------------------------------------------------------
   subroutine WaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        soilhydrology_vars, soilstate_vars, temperature_vars, waterstate_vars, waterflux_vars) 
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use clm_time_manager , only : get_step_size
     use clm_varcon       , only : pondmx, tfrz, watmin,denice,denh2o
     use clm_varpar       , only : nlevsoi, nlevgrnd
     use column_varcon    , only : icol_roof, icol_road_imperv
     use clm_varctl       , only : use_vsfm, do_varsoil
!     use SoilWaterMovementMod, only : soilroot_water_method, zengdecker_2009, vsfm
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds  
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(temperature_type)   , intent(in)    :: temperature_vars
     type(waterstate_type)    , intent(inout) :: waterstate_vars
     type(waterflux_type)     , intent(inout) :: waterflux_vars
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                                ! indices
     integer  :: nlevbed                                 ! # layers to bedrock
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: dflag=0._r8
     !-----------------------------------------------------------------------

     associate(                                                            & 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers                              
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         nlev2bed            =>    col%nlev2bed                          , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     

          t_soisno           =>    temperature_vars%t_soisno_col         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       

          h2osfc             =>    waterstate_vars%h2osfc_col            , & ! Input:  [real(r8) (:)   ]  surface water (mm)                                
          h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)                                
          h2osoi_vol         =>    waterstate_vars%h2osoi_vol_col        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          frac_h2osfc        =>    waterstate_vars%frac_h2osfc_col       , & ! Input:  [real(r8) (:)   ]                                                    

          qflx_dew_grnd      =>    waterflux_vars%qflx_dew_grnd_col      , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterflux_vars%qflx_dew_snow_col      , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]    

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         

          zwt                =>    soilhydrology_vars%zwt_col            , & ! Output: [real(r8) (:)   ]  water table depth (m)                             
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)                     
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Output: [real(r8) (:)   ]  frost table depth (m)                             
          wa                 =>    soilhydrology_vars%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)              
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical  
          
          qflx_sub_snow      =>    waterflux_vars%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]   
          qflx_drain         =>    waterflux_vars%qflx_drain_col         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
          qflx_drain_perched =>    waterflux_vars%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s)         
          qflx_rsub_sat      =>    waterflux_vars%qflx_rsub_sat_col        & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]                 
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

       do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if
          do j = 1,nlevbed
             dzmm(c,j) = dz(c,j)*1.e3_r8
          end do
       end do

       if (.not.use_vsfm) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             qflx_drain(c)    = 0._r8
             qflx_rsub_sat(c) = 0._r8
             qflx_drain_perched(c)  = 0._r8
          end do
       endif

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if
          jwt(c) = nlevbed
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevbed
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       !============================== QCHARGE =========================================
       ! Water table changes due to qcharge
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if

          !scs: use analytical expression for aquifer specific yield
          rous = watsat(c,nlevbed) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
          rous=max(rous,0.02_r8)

          !--  water table is below the soil column  --------------------------------------
          if(jwt(c) == nlevbed) then             
	    if(.not. (soilroot_water_method .eq. zengdecker2009 .and. do_varsoil)) then
             wa(c)  = wa(c) + qcharge(c)  * dtime
             zwt(c) = zwt(c) - (qcharge(c)  * dtime)/1000._r8/rous
            end if
          else                                
             !-- water table within soil layers 1-9  -------------------------------------
             ! try to raise water table to account for qcharge
             qcharge_tot = qcharge(c) * dtime
             if(qcharge_tot > 0.) then !rising water table
                do j = jwt(c)+1, 1,-1
                   !scs: use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=min(qcharge_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                   qcharge_layer=max(qcharge_layer,0._r8)

                   if(s_y > 0._r8) zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8

                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot <= 0.) exit
                enddo
             else ! deepening water table (negative qcharge)
                do j = jwt(c)+1, nlevbed
                   !scs: use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=max(qcharge_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                   qcharge_layer=min(qcharge_layer,0._r8)
                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot >= 0.) then 
                      zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8
                      exit
                   else
                      zwt(c) = zi(c,j)
                   endif

                enddo
                if (qcharge_tot > 0.) zwt(c) = zwt(c) - qcharge_tot/1000._r8/rous
             endif

             !-- recompute jwt for following calculations  ---------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevbed
             do j = 1,nlevbed
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          endif
       enddo


       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then 
             k_frz=nlevbed
          else
             k_frz=1
          endif

          do k=2, nlevbed
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
          else 
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then 
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

             endif !k_frz > k_perch 
          endif
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Renew the ice and liquid mass due to condensation

          if (snl(c)+1 >= 1) then

             ! make consistent with how evap_grnd removed in infiltration
             if (.not.use_vsfm) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtime
                end if
             end if
          end if
       end do


       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

          if (col%itype(c) == icol_roof .or. col%itype(c) == icol_road_imperv) then
             if (snl(c)+1 >= 1) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtime)
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtime)
                end if
             end if
          end if

       end do

     end associate

   end subroutine WaterTable

   !-----------------------------------------------------------------------
   subroutine Drainage(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc,  &
        temperature_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
     use clm_time_manager , only : get_step_size
     use clm_varpar       , only : nlevsoi, nlevgrnd, nlayer, nlayert
     use clm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     use abortutils       , only : endrun
     use clm_varctl       , only : use_vsfm, do_varsoil
!     use SoilWaterMovementMod, only : soilroot_water_method, zengdecker_2009, vsfm
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds               
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(temperature_type)   , intent(in)    :: temperature_vars
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(waterstate_type)    , intent(inout) :: waterstate_vars
     type(waterflux_type)     , intent(inout) :: waterflux_vars
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'Drainage'           ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     integer  :: nlevbed                                 ! # layers to bedrock
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevsoi)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: wh_zwt                                  ! water head at the water table depth (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: rsub_top_max
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     !-----------------------------------------------------------------------

     associate(                                                            & 
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)           
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)                                 
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers                              
         nlev2bed            =>    col%nlev2bed                          , & ! Input:  [integer  (:)   ]  number of layers to bedrock                     

          t_soisno           =>    temperature_vars%t_soisno_col         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)                       

          h2osfc             =>    waterstate_vars%h2osfc_col            , & ! Input:  [real(r8) (:)   ] surface water (mm)                                

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"                        
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)                       
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)  
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice         
          hk_l               =>    soilstate_vars%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)                    
          smp_l              =>    soilstate_vars%smp_l_col              , & ! Input:  [real(r8) (:,:) ] soil matric potential (mm)                    

          depth              =>    soilhydrology_vars%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth                                   
          c_param            =>    soilhydrology_vars%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)                             
          Dsmax              =>    soilhydrology_vars%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_vars%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_vars%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          Ds                 =>    soilhydrology_vars%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_vars%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_vars%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer                         
          hkdepth            =>    soilhydrology_vars%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)                                  
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)                             
          zwt                =>    soilhydrology_vars%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)                             
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)                     
          wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)              
          ice                =>    soilhydrology_vars%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)                         
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)                      
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_vars%h2osfcflag         , & ! Input:  logical
          
          qflx_snwcp_liq     =>    waterflux_vars%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_snwcp_ice     =>    waterflux_vars%qflx_snwcp_ice_col     , & ! Output: [real(r8) (:)   ] excess snowfall due to snow capping (mm H2O /s) [+]
          qflx_dew_grnd      =>    waterflux_vars%qflx_dew_grnd_col      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]      
          qflx_dew_snow      =>    waterflux_vars%qflx_dew_snow_col      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]    
          qflx_sub_snow      =>    waterflux_vars%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]   
          qflx_drain         =>    waterflux_vars%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)                    
          qflx_qrgwl         =>    waterflux_vars%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterflux_vars%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]                 
          qflx_drain_perched =>    waterflux_vars%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)         

          h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)                            
          h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)                                
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

        do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)
          else
            nlevbed = nlevsoi
          end if
          do j = 1,nlevbed
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)    = 0._r8 
          rsub_bot(c)      = 0._r8
          qflx_rsub_sat(c) = 0._r8
          rsub_top(c)      = 0._r8
          fracice_rsub(c)  = 0._r8
          qflx_qrgwl(c)    = 0._r8
       end do

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)             ! Added by MAB, 5/13/16
          else
            nlevbed = nlevsoi
          end if
          jwt(c) = nlevbed
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevbed
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       rous = 0.2_r8

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)             ! Added by MAB, 5/13/16
          else
            nlevbed = nlevsoi
          end if

          !  specify maximum drainage rate
          q_perch_max = 1.e-5_r8 * sin(col%topo_slope(c) * (rpi/180._r8))

          ! if layer containing water table is frozen, compute the following:
          !     frost table, perched water table, and drainage from perched saturated layer

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then 
             k_frz=nlevbed
          else
             k_frz=1
          endif

          do k=2, nlevbed
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)
          qflx_drain_perched(c) = 0._r8

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation

          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub

             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

             ! remove drainage from perched saturated layers
             rsub_top_tot = -  qflx_drain_perched(c) * dtime
             do k = jwt(c)+1, k_frz
                rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                if (use_vsfm) then
                   rsub_top_layer = 0._r8
                endif
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

                if (rsub_top_tot >= 0.) then 
                   zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                   exit
                else
                   zwt(c) = zi(c,k)
                endif
             enddo

             ! if rsub_top_tot is greater than available water (above frost table), 
             !     then decrease qflx_drain_perched by residual amount for water balance
             qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             !-- recompute jwt  ---------------------------------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevbed
             do j = 1,nlevbed
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          else 
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol/watsat(c,k) <= sat_lev) then 
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

                ! compute drainage from perched saturated region
                wtsub = 0._r8
                q_perch = 0._r8
                do k = k_perch, k_frz
                   imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

                ! no perched water table drainage if using original formulation
                if(origflag == 1) qflx_drain_perched(c) = 0._r8

                ! remove drainage from perched saturated layers
                rsub_top_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   if (use_vsfm) rsub_top_layer = 0._r8
                   rsub_top_tot = rsub_top_tot - rsub_top_layer

                   h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer

                   if (rsub_top_tot >= 0.) then 
                      zwt_perched(c) = zwt_perched(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                      exit
                   else
                      zwt_perched(c) = zi(c,k)
                   endif

                enddo

                ! if rsub_top_tot is greater than available water (above frost table), 
                !     then decrease qflx_drain_perched by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch

             !-- Topographic runoff  ----------------------------------------------------------------------
             fff(c)         = 1._r8/ hkdepth(c)
             dzsum = 0._r8
             icefracsum = 0._r8
             do j = max(jwt(c),1), nlevsoi
                dzsum  = dzsum + dzmm(c,j)
                icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
             end do
             ! add ice impedance factor to baseflow
             if(origflag == 1) then 
                if (use_vichydro) then
                   call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
                else
                   fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum))) &
                        - exp(-3._r8))/(1.0_r8-exp(-3._r8))
                   imped=(1._r8 - fracice_rsub(c))
                   rsub_top_max = 5.5e-3_r8
                end if
             else
                if (use_vichydro) then
                   imped=10._r8**(-e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                   dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                   rsub_top_max = dsmax_tmp(c)
                else
                   imped=10._r8**(-e_ice*(icefracsum/dzsum))
                   rsub_top_max = 10._r8 * sin((rpi/180.) * col%topo_slope(c))
                end if
             endif
             if (use_vichydro) then
                ! ARNO model for the bottom soil layer (based on bottom soil layer 
                ! moisture from previous time step
                ! use watmin instead for resid_moist to be consistent with default hydrology
                rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin) 
                frac = (Ds(c) * rsub_top_max )/Wsvic(c)
                rsub_tmp = (frac * rel_moist)/dtime
                if(rel_moist > Wsvic(c))then
                   frac = (rel_moist - Wsvic(c))/(1.0_r8 - Wsvic(c))
                   rsub_tmp = rsub_tmp + (rsub_top_max * (1.0_r8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtime
                end if
                rsub_top(c) = imped * rsub_tmp
                ! make sure baseflow isn't negative
                rsub_top(c) = max(0._r8, rsub_top(c))
             else
	        if(jwt(c) == nlevbed .and. soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
                  rsub_top(c)    = 0._r8
                else
                  rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*zwt(c))
		end if
             end if

             if (use_vsfm) rsub_top(c) = 0._r8

             ! use analytical expression for aquifer specific yield
             rous = watsat(c,nlevbed) &
                  * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
             rous=max(rous,0.02_r8)

             !--  water table is below the soil column  --------------------------------------
             if(jwt(c) == nlevbed) then             
	      if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
         	 if(-1._r8 * smp_l(c,nlevbed) < 0.5_r8 * dzmm(c,nlevbed)) then
           	   zwt(c) = z(c,nlevbed) - (smp_l(c,nlevbed) / 1000._r8)
		 end if
                 rsub_top(c) = imped * rsub_top_max * exp(-fff(c) * zwt(c))
                 rsub_top_tot = - rsub_top(c) * dtime
                 s_y = watsat(c,nlevbed) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
                 s_y=max(s_y,0.02_r8)     
                 rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,nlevbed) - zwt(c))*1.e3))
                 rsub_top_layer=min(rsub_top_layer,0._r8)
                 h2osoi_liq(c,nlevbed) = h2osoi_liq(c,nlevbed) + rsub_top_layer
                 rsub_top_tot = rsub_top_tot - rsub_top_layer
                 if (rsub_top_tot >= 0.) then 
                   zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8
                 else
                   zwt(c) = zi(c,nlevbed)
                 end if
	         if(rsub_top_tot < 0.) then
	           rsub_top(c) = rsub_top(c) + rsub_top_tot / dtime
                   rsub_top_tot = 0.
                 end if
              else
                wa(c)  = wa(c) - rsub_top(c) * dtime
                zwt(c)     = zwt(c) + (rsub_top(c) * dtime)/1000._r8/rous
                h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) + max(0._r8,(wa(c)-5000._r8))
                wa(c)  = min(wa(c), 5000._r8)
              end if
             else                                
                !-- water table within soil layers 1-9  -------------------------------------
                !============================== RSUB_TOP =========================================
                !--  Now remove water via rsub_top
                rsub_top_tot = - rsub_top(c) * dtime
                !should never be positive... but include for completeness
                if(rsub_top_tot > 0.) then !rising water table

                   call endrun(msg="RSUB_TOP IS POSITIVE in Drainage!"//errmsg(__FILE__, __LINE__))

                else ! deepening water table
                   if (use_vichydro) then
                      wtsub_vic = 0._r8
                      do j = (nlvic(1)+nlvic(2)+1), nlevbed
                         wtsub_vic = wtsub_vic + hk_l(c,j)*dzmm(c,j)
                      end do

                      do j = (nlvic(1)+nlvic(2)+1), nlevbed
                         rsub_top_layer=max(rsub_top_tot, rsub_top_tot*hk_l(c,j)*dzmm(c,j)/wtsub_vic)
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         if (use_vsfm) rsub_top_layer = 0._r8
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                         rsub_top_tot = rsub_top_tot - rsub_top_layer
                      end do
                   else
                      do j = jwt(c)+1, nlevbed
                         ! use analytical expression for specific yield
                         s_y = watsat(c,j) &
                              * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                         s_y=max(s_y,0.02_r8)

                         rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         if (use_vsfm) rsub_top_layer = 0._r8
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer

                         rsub_top_tot = rsub_top_tot - rsub_top_layer

                         if (rsub_top_tot >= 0.) then 
                            zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8

                            exit
                         else
                            zwt(c) = zi(c,j)
                         endif
                      enddo
                   end if

                   !--  remove residual rsub_top  ---------------------------------------------
!                   zwt(c) = zwt(c) - rsub_top_tot/1000._r8/rous
!                   wa(c) = wa(c) + rsub_top_tot
            	    if(rsub_top_tot < 0.) then
              	      rsub_top(c) = rsub_top(c) + rsub_top_tot / dtime
              	      rsub_top_tot = 0._r8
            	    end if
                endif

                !-- recompute jwt  ---------------------------------------------------------
                ! allow jwt to equal zero when zwt is in top layer
                jwt(c) = nlevbed
                do j = 1,nlevbed
                   if(zwt(c) <= zi(c,j)) then
                      jwt(c) = j-1
                      exit
                   end if
                enddo
             end if! end of jwt if construct

             zwt(c) = max(0.0_r8,zwt(c))
             zwt(c) = min(80._r8,zwt(c))

          endif

       end do

       !  excessive water above saturation added to the above unsaturated layer like a bucket
       !  if column fully saturated, excess water goes to runoff

       do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)             ! Added by MAB, 5/13/16
          else
            nlevbed = nlevsoi
          end if
          do j = nlevbed,2,-1
             xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
             if (use_vsfm) then
                xsi(c) = 0._r8
             else
                h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
                h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
             endif
          end do
       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          !scs: watmin addition to fix water balance errors
          xs1(c)          = max(max(h2osoi_liq(c,1)-watmin,0._r8)- &
               max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)-watmin)),0._r8)
          if (use_vsfm) xs1(c) = 0._r8
          h2osoi_liq(c,1) = h2osoi_liq(c,1) - xs1(c)

          if (lun%urbpoi(col%landunit(c))) then
             qflx_rsub_sat(c)     = xs1(c) / dtime
          else
             if(h2osfcflag == 1) then
                ! send this water up to h2osfc rather than sending to drainage
                h2osfc(c) = h2osfc(c) + xs1(c)
                qflx_rsub_sat(c)     = 0._r8
             else
                ! use original code to send water to drainage (non-h2osfc case)
                qflx_rsub_sat(c)     = xs1(c) / dtime
             endif
          endif

          if (use_vsfm .or. do_varsoil) qflx_rsub_sat(c) = 0._r8

          ! add in ice check
          xs1(c)          = max(max(h2osoi_ice(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1))),0._r8)
          h2osoi_ice(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1)), h2osoi_ice(c,1))
          qflx_snwcp_ice(c) = qflx_snwcp_ice(c) + xs1(c) / dtime
       end do

       ! Limit h2osoi_liq to be greater than or equal to watmin.
       ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
       ! If insufficient water in soil layers, get from aquifer water

       do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)             ! Added by MAB, 5/13/16
          else
            nlevbed = nlevsoi
          end if
       	  do j = 1, nlevbed-1
             if (h2osoi_liq(c,j) < watmin) then
                xs(c) = watmin - h2osoi_liq(c,j)
                ! deepen water table if water is passed from below zwt layer
                if(j == jwt(c)) then 
                   zwt(c) = zwt(c) + xs(c)/eff_porosity(c,j)/1000._r8
                endif
             else
                xs(c) = 0._r8
             end if
             h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
          end do
       end do

       ! Get water for bottom layer from layers above if possible
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if(soilroot_water_method .eq. zengdecker2009 .and. do_varsoil) then
       	    nlevbed = nlev2bed(c)             ! Added by MAB, 5/13/16
          else
            nlevbed = nlevsoi
          end if
       	  j = nlevbed
          if (h2osoi_liq(c,j) < watmin) then
             xs(c) = watmin-h2osoi_liq(c,j)
             searchforwater: do i = nlevbed-1, 1, -1
                available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._r8)
                if (available_h2osoi_liq >= xs(c)) then
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - xs(c)
                   xs(c) = 0._r8
                   exit searchforwater
                else
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + available_h2osoi_liq
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - available_h2osoi_liq
                   xs(c) = xs(c) - available_h2osoi_liq
                end if
             end do searchforwater
          else
             xs(c) = 0._r8
          end if
          ! Needed in case there is no water to be found
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
          ! Instead of removing water from aquifer where it eventually
          ! shows up as excess drainage to the ocean, take it back out of 
          ! drainage
          rsub_top(c) = rsub_top(c) - xs(c)/dtime

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          ! Set imbalance for snow capping

          qflx_qrgwl(c) = qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine Drainage

   !-----------------------------------------------------------------------
   subroutine DrainageVSFM(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc,  &
        temperature_vars, soilhydrology_vars, soilstate_vars, waterstate_vars, waterflux_vars)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
     use clm_time_manager , only : get_step_size
     use clm_varpar       , only : nlevsoi, nlevgrnd, nlayer, nlayert
     use clm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     use abortutils       , only : endrun
     use clm_varctl       , only : use_vsfm
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(temperature_type)   , intent(in)    :: temperature_vars
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(waterstate_type)    , intent(inout) :: waterstate_vars
     type(waterflux_type)     , intent(inout) :: waterflux_vars
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'Drainage'           ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: dtime                                   ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevsoi)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: wh_zwt                                  ! water head at the water table depth (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: rsub_top_max
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     integer  :: idx                      ! 1D index for VSFM
     !-----------------------------------------------------------------------

     associate(                                                            &
          z                  =>    col%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          zi                 =>    col%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)
          dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          snl                =>    col%snl                               , & ! Input:  [integer  (:)   ] number of snow layers

          t_soisno           =>    temperature_vars%t_soisno_col         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)

          h2osfc             =>    waterstate_vars%h2osfc_col            , & ! Input:  [real(r8) (:)   ] surface water (mm)

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice
          hk_l               =>    soilstate_vars%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)

          depth              =>    soilhydrology_vars%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth
          c_param            =>    soilhydrology_vars%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)
          Dsmax              =>    soilhydrology_vars%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_vars%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_vars%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          Ds                 =>    soilhydrology_vars%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_vars%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_vars%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer
          hkdepth            =>    soilhydrology_vars%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)
          zwt                =>    soilhydrology_vars%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)
          wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
          ice                =>    soilhydrology_vars%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_vars%h2osfcflag         , & ! Input:  logical

          qflx_snwcp_liq     =>    waterflux_vars%qflx_snwcp_liq_col     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_snwcp_ice     =>    waterflux_vars%qflx_snwcp_ice_col     , & ! Output: [real(r8) (:)   ] excess snowfall due to snow capping (mm H2O /s) [+]
          qflx_dew_grnd      =>    waterflux_vars%qflx_dew_grnd_col      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]
          qflx_dew_snow      =>    waterflux_vars%qflx_dew_snow_col      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow      =>    waterflux_vars%qflx_sub_snow_col      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]
          qflx_drain         =>    waterflux_vars%qflx_drain_col         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)
          qflx_qrgwl         =>    waterflux_vars%qflx_qrgwl_col         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    waterflux_vars%qflx_rsub_sat_col      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]
          qflx_drain_perched =>    waterflux_vars%qflx_drain_perched_col , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)

          mflx_drain_perched_col_1d =>    waterflux_vars%mflx_drain_perched_col_1d   , & ! Input:  [real(r8) (:)   ]  drainage from perched water table (kg H2O /s)

          h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
          )

       ! Get time step

       dtime = get_step_size()

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevgrnd
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)          = 0._r8
          rsub_bot(c)            = 0._r8
          qflx_rsub_sat(c)       = 0._r8
          rsub_top(c)            = 0._r8
          fracice_rsub(c)        = 0._r8
          qflx_qrgwl(c)          = 0._r8
          qflx_drain_perched(c)  = 0._r8
       end do

       mflx_drain_perched_col_1d(:) = 0._r8

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevgrnd
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevgrnd
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       rous = 0.2_r8

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          !  specify maximum drainage rate
          q_perch_max = 1.e-5_r8 * sin(col%topo_slope(c) * (rpi/180._r8))

          ! if layer containing water table is frozen, compute the following:
          !     frost table, perched water table, and drainage from perched saturated layer

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then
             k_frz=nlevsoi
          else
             k_frz=1
          endif

          do k=2, nlevsoi
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)
          qflx_drain_perched(c) = 0._r8

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then

             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub

             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

             ! remove drainage from perched saturated layers
             rsub_top_tot = -  qflx_drain_perched(c) * dtime
             do k = jwt(c)+1, k_frz
                rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                ! save the flux for VSFM
                idx = (c-bounds%begc)*nlevgrnd + k
                mflx_drain_perched_col_1d(idx) = rsub_top_layer/dtime

                if (rsub_top_tot >= 0.) then
                   zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                   exit
                else
                   zwt(c) = zi(c,k)
                endif
             enddo

             ! if rsub_top_tot is greater than available water (above frost table),
             !     then decrease qflx_drain_perched by residual amount for water balance
             qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             !-- recompute jwt  ---------------------------------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevsoi
             do j = 1,nlevsoi
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          else
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol/watsat(c,k) <= sat_lev) then
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then

                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

                ! compute drainage from perched saturated region
                wtsub = 0._r8
                q_perch = 0._r8
                do k = k_perch, k_frz
                   imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

                ! no perched water table drainage if using original formulation
                if(origflag == 1) qflx_drain_perched(c) = 0._r8

                ! remove drainage from perched saturated layers
                rsub_top_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   rsub_top_tot = rsub_top_tot - rsub_top_layer

                   ! save the flux for VSFM
                   idx = (c-bounds%begc)*nlevgrnd + k
                   mflx_drain_perched_col_1d(idx) = rsub_top_layer/dtime

                   if (rsub_top_tot >= 0.) then
                      zwt_perched(c) = zwt_perched(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                      exit
                   else
                      zwt_perched(c) = zi(c,k)
                   endif

                enddo

                ! if rsub_top_tot is greater than available water (above frost table),
                !     then decrease qflx_drain_perched by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch

             !-- Topographic runoff  ----------------------------------------------------------------------
             fff(c)         = 1._r8/ hkdepth(c)
             dzsum = 0._r8
             icefracsum = 0._r8
             do j = max(jwt(c),1), nlevgrnd
                dzsum  = dzsum + dzmm(c,j)
                icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
             end do
             ! add ice impedance factor to baseflow
             if(origflag == 1) then
                if (use_vichydro) then
                   call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
                else
                   fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum))) &
                        - exp(-3._r8))/(1.0_r8-exp(-3._r8))
                   imped=(1._r8 - fracice_rsub(c))
                   rsub_top_max = 5.5e-3_r8
                end if
             else
                if (use_vichydro) then
                   imped=10._r8**(-e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                   dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                   rsub_top_max = dsmax_tmp(c)
                else
                   imped=10._r8**(-e_ice*(icefracsum/dzsum))
                   rsub_top_max = 10._r8 * sin((rpi/180.) * col%topo_slope(c))
                end if
             endif
             if (use_vichydro) then
                ! ARNO model for the bottom soil layer (based on bottom soil layer
                ! moisture from previous time step
                ! use watmin instead for resid_moist to be consistent with default hydrology
                rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin)
                frac = (Ds(c) * rsub_top_max )/Wsvic(c)
                rsub_tmp = (frac * rel_moist)/dtime
                if(rel_moist > Wsvic(c)) then
                   frac = (rel_moist - Wsvic(c))/(1.0_r8 - Wsvic(c))
                   rsub_tmp = rsub_tmp + (rsub_top_max * (1.0_r8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtime
                end if
                rsub_top(c) = imped * rsub_tmp
                ! make sure baseflow isn't negative
                rsub_top(c) = max(0._r8, rsub_top(c))
             else
                rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*zwt(c))
             end if
          endif

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          ! Set imbalance for snow capping

          qflx_qrgwl(c) = qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine DrainageVSFM

  !-----------------------------------------------------------------------
  subroutine CLMVICMap(bounds, numf, filter, &
       soilhydrology_vars, waterstate_vars)
     !
     ! !DESCRIPTION:
     ! Performs  the mapping from CLM layers to VIC layers
     ! Specifically, 10 (or 23 when more_vertlayers == .true.) 
     ! CLM hydrologically active soil layers are mapped to three VIC layers
     ! by assigning the first nlvic(1) layers to VIC layer 1
     !              the next nlvic(2) layers  to VIC alyer 2
     !              and the remaining to VIC layer 3
     ! mapping from VIC to CLM layers, M.Huang
     !
     ! !USES:
     use clm_varcon  , only : denh2o, denice, watmin
     use clm_varpar  , only : nlevsoi, nlayer, nlayert, nlevgrnd 
     use decompMod   , only : bounds_type
     !
     ! !REVISION HISTORY:
     ! Created by Maoyi Huang
     ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC 
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds    
     integer                  , intent(in)    :: numf      ! number of column soil points in column filter
     integer                  , intent(in)    :: filter(:) ! column filter for soil points
     type(waterstate_type)    , intent(in)    :: waterstate_vars 
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     !
     ! !LOCAL VARIABLES
     real(r8) :: ice0(1:nlayer)            ! last step ice lens (mm)  (new)
     real(r8) :: moist0(1:nlayer)          ! last step soil water (mm)  (new)
     integer  :: i, j, c, fc
     ! note: in CLM4 h2osoil_liq unit is kg/m2, in VIC moist is mm
     ! h2osoi_ice is actually water equivalent ice content.
     !-----------------------------------------------------------------------

     associate(                                                   & 
          dz            => col%dz                               , & ! Input:  [real(r8) (:,:)   ] layer depth (m)                        
          zi            => col%zi                               , & ! Input:  [real(r8) (:,:)   ] interface level below a "z" level (m)  
          z             => col%z                                , & ! Input:  [real(r8) (:,:)   ] layer thickness (m)                    

          h2osoi_liq    => waterstate_vars%h2osoi_liq_col       , & ! Input:  [real(r8) (:,:)   ] liquid water (kg/m2)                   
          h2osoi_ice    => waterstate_vars%h2osoi_ice_col       , & ! Input:  [real(r8) (:,:)   ] ice lens (kg/m2)                       
          h2osoi_vol    => waterstate_vars%h2osoi_vol_col       , & ! Input:  [real(r8) (:,:)   ] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)

          depth         => soilhydrology_vars%depth_col         , & ! Input:  [real(r8) (:,:)   ] layer depth of upper layer (m)         
          porosity      => soilhydrology_vars%porosity_col      , & ! Input:  [real(r8) (:,:)   ] soil porisity (1-bulk_density/soil_density)
          max_moist     => soilhydrology_vars%max_moist_col     , & ! Input:  [real(r8) (:,:)   ] max layer moist + ice (mm)             
          vic_clm_fract => soilhydrology_vars%vic_clm_fract_col , & ! Input:  [real(r8) (:,:,:) ] fraction of VIC layers in each CLM layer
          moist         => soilhydrology_vars%moist_col         , & ! Output: [real(r8) (:,:)   ] liquid water (mm)                      
          ice           => soilhydrology_vars%ice_col           , & ! Output: [real(r8) (:,:)   ] ice lens (mm)                          
          moist_vol     => soilhydrology_vars%moist_vol_col       & ! Output: [real(r8) (:,:)   ] volumetric soil moisture for VIC soil layers
          )

       ! map CLM to VIC
       do fc = 1, numf
          c = filter(fc)
          do i = 1, nlayer
             ice0(i) = ice(c,i)
             moist0(i) = moist(c,i)
             ice(c,i) = 0._r8
             moist(c,i) = 0._r8
             do j = 1, nlevsoi
                ice(c,i) = ice(c,i) + h2osoi_ice(c,j) * vic_clm_fract(c,i,j)
                moist(c,i) = moist(c,i) + h2osoi_liq(c,j) * vic_clm_fract(c,i,j)
             end do
             ice(c,i)       = min((moist0(i) + ice0(i)), ice(c,i))
             ice(c,i)       = max(0._r8, ice(c,i))
             moist(c,i)     = max(watmin, moist(c,i))
             moist(c,i)     = min(max_moist(c,i)-ice(c,i), moist(c,i))
             moist_vol(c,i) = moist(c,i)/(depth(c,i)*denice) + ice(c,i)/(depth(c,i)*denh2o)
             moist_vol(c,i) = min(porosity(c,i), moist_vol(c,i))
             moist_vol(c,i) = max(0.01_r8, moist_vol(c,i))
          end do

          ! hydrologic inactive layers
          ice(c, nlayer+1:nlayert)       = h2osoi_ice(c, nlevsoi+1:nlevgrnd)
          moist(c, nlayer+1:nlayert)     = h2osoi_liq(c, nlevsoi+1:nlevgrnd)
          moist_vol(c, nlayer+1:nlayert) = h2osoi_vol(c, nlevsoi+1:nlevgrnd)
       end do

     end associate

   end subroutine CLMVICMap

end module SoilHydrologyMod
