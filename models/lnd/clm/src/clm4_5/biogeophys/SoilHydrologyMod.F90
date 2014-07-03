module SoilHydrologyMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil hydrology
  !
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  use clm_varctl   , only: iulog, use_vichydro
  use clm_varcon   , only: e_ice,denh2o, denice,rpi
  use decompMod    , only: bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilHydrology_readnl ! Initialization for Soil Hydrology
  public :: SurfaceRunoff        ! Calculate surface runoff
  public :: Infiltration         ! Calculate infiltration into surface soil layer
  public :: SoilWater            ! Calculate soil hydrology
  public :: Drainage             ! Calculate subsurface drainage
  public :: WaterTable           ! Calculate water table before imposing drainage
  !
  ! !PUBLIC DATA MEMBERS:
  integer, public :: h2osfcflag=1               !If surface water is active or not
  !
  ! !PRIVATE DATA MEMBERS:
  integer         :: origflag=0                 !use control soil hydraulic properties
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilHydrology_readnl( NLFilename )
    !
    ! !DESCRIPTION:
    ! Read namelist for SoilHydrology
    !
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use clm_nlUtilsMod, only : find_nlgroup_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
    !
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'SoilHydrology_readnl'  ! subroutine name
    !-----------------------------------------------------------------------

    namelist / clm_soilhydrology_inparm / h2osfcflag, origflag

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input. 
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm_soilhydrology_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call find_nlgroup_name(unitn, 'clm_soilhydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_soilhydrology_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading clm_soilhydrology_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(h2osfcflag, mpicom)
    call shr_mpi_bcast(origflag,   mpicom)

  end subroutine SoilHydrology_readnl

  !-----------------------------------------------------------------------
  subroutine SurfaceRunoff (bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff
    !
    ! !USES:
    use clmtype
    use clm_varcon      , only : denice, denh2o, wimp, pondmx_urban
    use clm_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varcon      , only : icol_road_imperv, icol_road_perv
    use clm_varpar      , only : nlevsoi, maxpatch_pft
    use clm_time_manager, only : get_step_size
    use clm_varpar      , only : nlayer, nlayert
    use abortutils      , only : endrun
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in) :: bounds               ! bounds
    integer           , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer           , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    integer           , intent(in) :: num_urbanc           ! number of column urban points in column filter
    integer           , intent(in) :: filter_urbanc(:)     ! column filter for urban points
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,g,l,i                               !indices
    real(r8) :: dtime                                      !land model time step (sec)
    real(r8) :: xs(bounds%begc:bounds%endc)                !excess soil water above urban ponding limit
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi) !partial volume of ice lens in layer
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

   associate(& 
   qflx_snow_h2osfc  =>    cwf%qflx_snow_h2osfc  , & ! Input:  [real(r8) (:)]  snow falling on surface water (mm/s)              
   frac_h2osfc       =>    cps%frac_h2osfc       , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   frost_table       =>    cws%frost_table       , & ! Input:  [real(r8) (:)]  frost table depth (m)                             
   zwt_perched       =>    cws%zwt_perched       , & ! Input:  [real(r8) (:)]  perched water table depth (m)                     
   qflx_floodc       =>    cwf%qflx_floodc       , & ! Input:  [real(r8) (:)]  column flux of flood water from RTM               
   cgridcell         =>    col%gridcell          , & ! Input:  [integer (:)]  gridcell index for each column                     
   clandunit         =>    col%landunit          , & ! Input:  [integer (:)]  column's landunit                                  
   ltype             =>    lun%itype             , & ! Input:  [integer (:)]  landunit type                                      
   ctype             =>    col%itype             , & ! Input:  [integer (:)]  column type index                                  
   qflx_top_soil     =>    cwf%qflx_top_soil     , & ! Input:  [real(r8) (:)]  net water input into soil from top (mm/s)         
   qflx_surf         =>    cwf%qflx_surf         , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)                        
   watsat            =>    cps%watsat            , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   hkdepth           =>    cps%hkdepth           , & ! Input:  [real(r8) (:)]  decay factor (m)                                  
   dz                =>    cps%dz                , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   h2osoi_ice        =>    cws%h2osoi_ice        , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq        =>    cws%h2osoi_liq        , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   fcov              =>    cws%fcov              , & ! Input:  [real(r8) (:)]  fractional impermeable area                       
   fsat              =>    cws%fsat              , & ! Input:  [real(r8) (:)]  fractional area with water table at surface       
   eff_porosity      =>    cps%eff_porosity      , & ! Output: [real(r8) (:,:)]  effective porosity = porosity - vol_ice         
   wtfact            =>    cps%wtfact            , & ! Input:  [real(r8) (:)]  maximum saturated fraction for a gridcell         
   zwt               =>    cws%zwt               , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   fracice           =>    cps%fracice           , & ! Output: [real(r8) (:,:)] fractional impermeability (-)                    
   icefrac           =>    cps%icefrac           , & ! Output: [real(r8) (:,:)]                                                  
   hksat             =>    cps%hksat             , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s)
   bsw               =>    cps%bsw               , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   sucsat            =>    cps%sucsat            , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   snl               =>    cps%snl               , & ! Input:  [integer (:)]  minus number of snow layers                        
   qflx_evap_grnd    =>    pwf_a%qflx_evap_grnd  , & ! Input:  [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]    
   b_infil           =>    cps%b_infil           , & ! Output: [real(r8) (:)] VIC b infiltration parameter                       
   max_moist         =>    cps%max_moist         , & ! Output: [real(r8) (:,:)] maximum soil moisture (ice + liq, mm)            
   moist             =>    cws%moist             , & ! Output: [real(r8) (:,:)] soil moisture in each VIC layers (liq, mm)       
   ice               =>    cws%ice               , & ! Output: [real(r8) (:,:)] ice len in each VIC layers(ice, mm)              
   max_infil         =>    cws%max_infil         , & ! Output: [real(r8) (:)] maximum infiltration capacity in VIC (mm)          
   i_0               =>    cws%i_0               , & ! Output: [real(r8) (:)] column average soil moisture in top VIC layers (mm)
   zi                =>    cps%zi                  & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)           
   )

    ! Get time step

    dtime = get_step_size()

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

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
       if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then

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
       else if (ctype(c) == icol_sunwall .or. ctype(c) == icol_shadewall) then
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
   subroutine Infiltration(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Calculate infiltration into surface soil layer (minus the evaporation)
     !
     ! !USES:
     use clmtype
     use clm_varcon      , only : icol_roof, icol_road_imperv, icol_sunwall
     use clm_varcon      , only : icol_shadewall, icol_road_perv,denh2o, denice, roverg, wimp
     use clm_varcon      , only : istsoil,pc,mu,tfrz, istcrop
     use clm_time_manager, only : get_step_size
     use clm_varpar      , only : nlevsoi
     use H2OSfcMod       , only : FracH2oSfc
     use shr_const_mod   , only : shr_const_pi
     use clm_varpar      , only : nlayer, nlayert
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds               ! bounds
     integer, intent(in) :: num_hydrologyc                 ! number of column soil points in column filter
     integer, intent(in) :: filter_hydrologyc(:)           ! column filter for soil points
     integer, intent(in) :: num_urbanc                     ! number of column urban points in column filter
     integer, intent(in) :: filter_urbanc(:)               ! column filter for urban points
     !
     ! !LOCAL VARIABLES:
     integer :: c,j,l, fc                                   ! indices
     real(r8) :: dtime                                      ! land model time step (sec)
     real(r8) :: s1,su,v                                    ! variable to calculate qinmax
     real(r8) :: qinmax                                     ! maximum infiltration capacity (mm/s)
     real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi) ! partial volume of ice lens in layer
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

   associate(& 
   frost_table        =>    cws%frost_table       , & ! Input:  [real(r8) (:)]  frost table depth (m)                             
   zwt_perched        =>    cws%zwt_perched       , & ! Input:  [real(r8) (:)]  perched water table depth (m)                     
   fsat               =>    cws%fsat              , & ! Input:  [real(r8) (:)]  fractional area with water table at surface       
   h2osfc_thresh      =>    cps%h2osfc_thresh     , & ! Input:  [real(r8) (:)]  level at which h2osfc "percolates"                
   frac_sno           =>    cps%frac_sno_eff      , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)       
   qflx_evap_soi      =>    pwf_a%qflx_evap_soi   , & ! Input:  [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]    
   qflx_h2osfc_surf   =>    cwf%qflx_h2osfc_surf  , & ! Input:  [real(r8) (:)]  surface water runoff (mm/s)                       
   frac_h2osfc        =>    cps%frac_h2osfc       , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc             =>    cws%h2osfc            , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   h2osoi_ice         =>    cws%h2osoi_ice        , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_liq         =>    cws%h2osoi_liq        , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   sucsat             =>    cps%sucsat            , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   watsat             =>    cps%watsat            , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   icefrac            =>    cps%icefrac           , & ! Output: [real(r8) (:,:)]  fraction of ice                                 
   bsw                =>    cps%bsw               , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   t_soisno           =>    ces%t_soisno          , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)                       
   smpmin             =>    cps%smpmin            , & ! Input:  [real(r8) (:)]  restriction for min of soil potential (mm)        
   fcov               =>    cws%fcov              , & ! Input:  [real(r8) (:)]  fractional area with water table at surface       
   eff_porosity       =>    cps%eff_porosity      , & ! Input:  [real(r8) (:,:)]  effective porosity = porosity - vol_ice         
   hksat              =>    cps%hksat             , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s)
   hksat_min          =>    cps%hksat_min         , & ! Input:  [real(r8) (:,:)]  mineral hksat                                   
   dz                 =>    cps%dz                , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   h2osno             =>    cws%h2osno            , & ! Input:  [real(r8) (:)]  snow water (mm H2O)                               
   snow_depth         =>    cps%snow_depth        , & ! Input:  [real(r8) (:)]  snow height (m)                                   
   t_h2osfc           =>    ces%t_h2osfc          , & ! Input:  [real(r8) (:)]  soil temperature (Kelvin)                         
   clandunit          =>   col%landunit           , & ! Input:  [integer (:)]  column's landunit                                  
   ltype              =>    lun%itype             , & ! Input:  [integer (:)]  landunit type                                      
   topo_slope         =>    cps%topo_slope        , & ! Input:  [real(r8) (:)]  topographic slope                                 
   qflx_ev_snow       =>    pwf_a%qflx_ev_snow    , & ! Input:  [real(r8) (:)]  evaporation flux from snow (W/m**2) [+ to atm]    
   qflx_ev_soil       =>    pwf_a%qflx_ev_soil    , & ! Input:  [real(r8) (:)]  evaporation flux from soil (W/m**2) [+ to atm]    
   qflx_ev_h2osfc     =>    pwf_a%qflx_ev_h2osfc  , & ! Input:  [real(r8) (:)]  evaporation flux from h2osfc (W/m**2) [+ to atm]  
   zwt                =>    cws%zwt               , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   ctype              =>    col%itype             , & ! Input:  [integer (:)]  column type index                                  
   snl                =>    cps%snl               , & ! Input:  [integer (:)]  minus number of snow layers                        
   qflx_top_soil      =>    cwf%qflx_top_soil     , & ! Input:  [real(r8) (:)]  net water input into soil from top (mm/s)         
   qflx_surf          =>    cwf%qflx_surf         , & ! Input:  [real(r8) (:)]  surface runoff (mm H2O /s)                        
   qflx_infl          =>    cwf%qflx_infl         , & ! Output: [real(r8) (:)] infiltration (mm H2O /s)                           
   b_infil            =>    cps%b_infil           , & ! Input:  [real(r8) (:)] VIC b infiltration parameter                       
   max_moist          =>    cps%max_moist         , & ! Input:  [real(r8) (:,:)] maximum soil moisture (ice + liq, mm)            
   moist              =>    cws%moist             , & ! Input:  [real(r8) (:,:)] soil moisture in each VIC layers (liq, mm)       
   ice                =>    cws%ice               , & ! Input:  [real(r8) (:,:)] ice len in each VIC layers(ice, mm)              
   max_infil          =>    cws%max_infil         , & ! Input:  [real(r8) (:)] maximum infiltration capacity in VIC (mm)          
   i_0                =>    cws%i_0               , & ! Input:  [real(r8) (:)] column average soil moisture in top VIC layers (mm)
   qflx_evap_grnd     =>    pwf_a%qflx_evap_grnd    & ! Input:  [real(r8) (:)]  ground surface evaporation rate (mm H2O/s) [+]    
   )

    dtime = get_step_size()

    ! Infiltration into surface soil layer (minus the evaporation)
    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          ! Porosity of soil, partial volume of ice and liquid
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
          eff_porosity(c,j) = max(0.01_r8,watsat(c,j)-vol_ice(c,j))
          icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       ! partition moisture fluxes between soil and h2osfc       
       if (ltype(clandunit(c)) == istsoil .or. ltype(clandunit(c))==istcrop) then

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

          !2. remove evaporation (snow treated in SnowHydrology)
          qflx_in_soil(c) = qflx_in_soil(c) - (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
          qflx_in_h2osfc(c) =  qflx_in_h2osfc(c)  - frac_h2osfc(c) * qflx_ev_h2osfc(c)
       
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
             k_wet=1.0_r8 * sin((rpi/180.) * topo_slope(c))
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

       else
          ! non-vegetated landunits (i.e. urban) use original CLM4 code
          if (snl(c) >= 0) then
             ! when no snow present, sublimation is removed in Drainage
             qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) - qflx_evap_grnd(c)
          else
             qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) &
               - (1.0_r8 - frac_sno(c)) * qflx_ev_soil(c)
          end if
          qflx_h2osfc_surf(c) = 0._r8
       endif

    enddo

    ! No infiltration for impervious urban surfaces

    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       if (ctype(c) /= icol_road_perv) then
          qflx_infl(c) = 0._r8
       end if
    end do
    
    end associate 
   end subroutine Infiltration

   !-----------------------------------------------------------------------
   subroutine SoilWater(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Soil hydrology
     ! Soil moisture is predicted from a 10-layer model (as with soil
     ! temperature), in which the vertical soil moisture transport is governed
     ! by infiltration, runoff, gradient diffusion, gravity, and root
     ! extraction through canopy transpiration.  The net water applied to the
     ! surface layer is the snowmelt plus precipitation plus the throughfall
     ! of canopy dew minus surface runoff and evaporation.
     ! CLM3.5 uses a zero-flow bottom boundary condition.
     !
     ! The vertical water flow in an unsaturated porous media is described by
     ! Darcy's law, and the hydraulic conductivity and the soil negative
     ! potential vary with soil water content and soil texture based on the work
     ! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
     ! integrated over the layer thickness, in which the time rate of change in
     ! water mass must equal the net flow across the bounding interface, plus the
     ! rate of internal source or sink. The terms of water flow across the layer
     ! interfaces are linearly expanded by using first-order Taylor expansion.
     ! The equations result in a tridiagonal system equation.
     !
     ! Note: length units here are all millimeter
     ! (in temperature subroutine uses same soil layer
     ! structure required but lengths are m)
     !
     ! Richards equation:
     !
     ! d wat      d     d wat d psi
     ! ----- = - -- [ k(----- ----- - 1) ] + S
     !   dt      dz       dz  d wat
     !
     ! where: wat = volume of water per volume of soil (mm**3/mm**3)
     ! psi = soil matrix potential (mm)
     ! dt  = time step (s)
     ! z   = depth (mm)
     ! dz  = thickness (mm)
     ! qin = inflow at top (mm h2o /s)
     ! qout= outflow at bottom (mm h2o /s)
     ! s   = source/sink flux (mm h2o /s)
     ! k   = hydraulic conductivity (mm h2o /s)
     !
     !                       d qin                  d qin
     ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
     !                       d wat(j-1)             d wat(j)
     !                ==================|=================
     !                                  < qin
     !
     !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
     !
     !                                  > qout
     !                ==================|=================
     !                        d qout               d qout
     ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
     !                        d wat(j)             d wat(j+1)
     !
     !
     ! Solution: linearize k and psi about d wat and use tridiagonal
     ! system of equations to solve for d wat,
     ! where for layer j
     !
     !
     ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
     !
     ! !USES:
     use clmtype
     use clm_varcon    , only : wimp, icol_roof, icol_road_imperv,grav,hfus,tfrz
     use clm_varpar    , only : nlevsoi, max_pft_per_col, nlevgrnd
     use shr_const_mod , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
     use TridiagonalMod, only : Tridiagonal
     use clm_time_manager  , only : get_step_size
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds             ! bounds
     integer , intent(in)  :: num_hydrologyc             ! number of column soil points in column filter
     integer , intent(in)  :: filter_hydrologyc(:)       ! column filter for soil points
     integer , intent(in)  :: num_urbanc                 ! number of column urban points in column filter
     integer , intent(in)  :: filter_urbanc(:)           ! column filter for urban points
     !
     ! !LOCAL VARIABLES:
     integer  :: p,c,fc,j                  ! do loop indices
     integer  :: jtop(bounds%begc:bounds%endc)             ! top level at each column
     real(r8) :: dtime                     ! land model time step (sec)
     real(r8) :: hk(bounds%begc:bounds%endc,1:nlevsoi)     ! hydraulic conductivity [mm h2o/s]
     real(r8) :: dhkdw(bounds%begc:bounds%endc,1:nlevsoi)  ! d(hk)/d(vol_liq)
     real(r8) :: amx(bounds%begc:bounds%endc,1:nlevsoi+1)  ! "a" left off diagonal of tridiagonal matrix
     real(r8) :: bmx(bounds%begc:bounds%endc,1:nlevsoi+1)  ! "b" diagonal column for tridiagonal matrix
     real(r8) :: cmx(bounds%begc:bounds%endc,1:nlevsoi+1)  ! "c" right off diagonal tridiagonal matrix
     real(r8) :: rmx(bounds%begc:bounds%endc,1:nlevsoi+1)  ! "r" forcing term of tridiagonal matrix
     real(r8) :: zmm(bounds%begc:bounds%endc,1:nlevsoi+1)  ! layer depth [mm]
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi+1) ! layer thickness [mm]
     real(r8) :: den                       ! used in calculating qin, qout
     real(r8) :: dqidw0(bounds%begc:bounds%endc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i-1))
     real(r8) :: dqidw1(bounds%begc:bounds%endc,1:nlevsoi+1) ! d(qin)/d(vol_liq(i))
     real(r8) :: dqodw1(bounds%begc:bounds%endc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i))
     real(r8) :: dqodw2(bounds%begc:bounds%endc,1:nlevsoi+1) ! d(qout)/d(vol_liq(i+1))
     real(r8) :: dsmpdw(bounds%begc:bounds%endc,1:nlevsoi+1) ! d(smp)/d(vol_liq)
     real(r8) :: num                         ! used in calculating qin, qout
     real(r8) :: qin(bounds%begc:bounds%endc,1:nlevsoi+1)    ! flux of water into soil layer [mm h2o/s]
     real(r8) :: qout(bounds%begc:bounds%endc,1:nlevsoi+1)   ! flux of water out of soil layer [mm h2o/s]
     real(r8) :: s_node                    ! soil wetness
     real(r8) :: s1                        ! "s" at interface of layer
     real(r8) :: s2                        ! k*s**(2b+2)
     real(r8) :: smp(bounds%begc:bounds%endc,1:nlevsoi)    ! soil matrix potential [mm]
     real(r8) :: sdamp                     ! extrapolates soiwat dependence of evaporation
     integer  :: pi                        ! pft index
     real(r8) :: temp(bounds%begc:bounds%endc)             ! accumulator for rootr weighting
     integer  :: jwt(bounds%begc:bounds%endc)              ! index of the soil layer right above the water table (-)
     real(r8) :: smp1,dsmpdw1,wh,wh_zwt,ka
     real(r8) :: dwat2(bounds%begc:bounds%endc,1:nlevsoi+1)
     real(r8) :: dzq                         ! used in calculating qin, qout (difference in equilbirium matric potential)
     real(r8) :: zimm(bounds%begc:bounds%endc,0:nlevsoi)     ! layer interface depth [mm]
     real(r8) :: zq(bounds%begc:bounds%endc,1:nlevsoi+1)     ! equilibrium matric potential for each layer [mm]
     real(r8) :: vol_eq(bounds%begc:bounds%endc,1:nlevsoi+1) ! equilibrium volumetric water content
     real(r8) :: tempi                       ! temp variable for calculating vol_eq
     real(r8) :: temp0                       ! temp variable for calculating vol_eq
     real(r8) :: voleq1                      ! temp variable for calculating vol_eq
     real(r8) :: zwtmm(bounds%begc:bounds%endc)              ! water table depth [mm]
     real(r8) :: imped(bounds%begc:bounds%endc,1:nlevsoi)             
     real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi)
     real(r8) :: z_mid
     real(r8) :: vwc_zwt(bounds%begc:bounds%endc)
     real(r8) :: vwc_liq(bounds%begc:bounds%endc,1:nlevsoi+1) ! liquid volumetric water content
     real(r8) :: smp_grad(bounds%begc:bounds%endc,1:nlevsoi+1)
     !-----------------------------------------------------------------------

   associate(& 
   h2osoi_ice       =>    cws%h2osoi_ice     , & ! Input:  [real(r8) (:,:)]  ice water (kg/m2)                               
   qcharge          =>    cws%qcharge        , & ! Input:  [real(r8) (:)]  aquifer recharge rate (mm/s)                      
   qflx_deficit     =>    cwf%qflx_deficit   , & ! Input:  [real(r8) (:)]  water deficit to keep non-negative liquid water content
   hkdepth          =>    cps%hkdepth        , & ! Input:  [real(r8) (:)]  decay factor (m)                                  
   zi               =>    cps%zi             , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)           
   zwt              =>    cws%zwt            , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   ctype            =>    col%itype          , & ! Input:  [integer (:)]  column type index                                  
   npfts            =>    col%npfts          , & ! Input:  [integer (:)]  column's number of pfts - ADD                      
   z                =>    cps%z              , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   dz               =>    cps%dz             , & ! Input:  [real(r8) (:,:)]  layer thickness (m)                             
   smpmin           =>    cps%smpmin         , & ! Input:  [real(r8) (:)]  restriction for min of soil potential (mm)        
   watsat           =>    cps%watsat         , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   hksat            =>    cps%hksat          , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s)
   bsw              =>    cps%bsw            , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   icefrac          =>    cps%icefrac        , & ! Input:  [real(r8) (:,:)]  fraction of ice                                 
   sucsat           =>    cps%sucsat         , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   eff_porosity     =>    cps%eff_porosity   , & ! Input:  [real(r8) (:,:)]  effective porosity = porosity - vol_ice         
   rootr_col        =>    cps%rootr_column   , & ! Input:  [real(r8) (:,:)]  effective fraction of roots in each soil layer  
   t_soisno         =>    ces%t_soisno       , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)                       
   h2osoi_liq       =>    cws%h2osoi_liq     , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_vol       =>    cws%h2osoi_vol     , & ! Input:  [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   qflx_infl        =>    cwf%qflx_infl      , & ! Input:  [real(r8) (:)]  infiltration (mm H2O /s)                          
   fracice          =>    cps%fracice        , & ! Input:  [real(r8) (:,:)]  fractional impermeability (-)                   
   qflx_tran_veg_col=>    pwf_a%qflx_tran_veg, & ! Input:  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)  
   pfti             =>    col%pfti           , & ! Input:  [integer (:)]  beginning pft index for each column                
   smp_l            =>    cws%smp_l          , & ! Input:  [real(r8) (:,:)]  soil matrix potential [mm]                      
   hk_l             =>    cws%hk_l           , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity (mm/s)                   
   pactive          =>    pft%active         , & ! Input:  [logical (:)]  true=>do computations on this pft 
   qflx_tran_veg_pft=>    pwf%qflx_tran_veg  , & ! Input:  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)  
   rootr_pft        =>    pps%rootr          , & ! Input:  [real(r8) (:,:)]  effective fraction of roots in each soil layer  
   pwtcol           =>    pft%wtcol            & ! Input:  [real(r8) (:)]  weight relative to column for each pft            
   )

    ! Get time step

    dtime = get_step_size()

    ! Because the depths in this routine are in mm, use local
    ! variable arrays instead of pointers

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          zmm(c,j) = z(c,j)*1.e3_r8
          dzmm(c,j) = dz(c,j)*1.e3_r8
          zimm(c,j) = zi(c,j)*1.e3_r8
          ! calculate icefrac up here
          vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
          icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
          vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_r8)/(dz(c,j)*denh2o)
       end do
    end do

    do fc = 1, num_hydrologyc 
       c = filter_hydrologyc(fc)
       zimm(c,0) = 0.0_r8
       zwtmm(c)  = zwt(c)*1.e3_r8
    end do

    ! First step is to calculate the column-level effective rooting
    ! fraction in each soil layer. This is done outside the usual
    ! PFT-to-column averaging routines because it is not a simple
    ! weighted average of the PFT level rootr arrays. Instead, the
    ! weighting depends on both the per-unit-area transpiration
    ! of the PFT and the PFTs area relative to all PFTs.

    temp(bounds%begc : bounds%endc) = 0._r8

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          rootr_col(c,j) = 0._r8
       end do
    end do

    do pi = 1,max_pft_per_col
       do j = 1,nlevsoi
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             if (pi <= npfts(c)) then
                p = pfti(c) + pi - 1
                if (pactive(p)) then
                   rootr_col(c,j) = rootr_col(c,j) + rootr_pft(p,j) * qflx_tran_veg_pft(p) * pwtcol(p)
                end if
             end if
          end do
       end do
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (pi <= npfts(c)) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                temp(c) = temp(c) + qflx_tran_veg_pft(p) * pwtcol(p)
             end if
          end if
       end do
    end do

    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (temp(c) /= 0._r8) then
             rootr_col(c,j) = rootr_col(c,j)/temp(c)
          end if
       end do
    end do

    !compute jwt index
    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       ! allow jwt to equal zero when zwt is in top layer
       do j = 1,nlevsoi
          if(zwt(c) <= zi(c,j)) then
             jwt(c) = j-1
             exit
          end if
       enddo

      ! compute vwc at water table depth (mainly for case when t < tfrz)
      !     this will only be used when zwt is below the soil column
       vwc_zwt(c) = watsat(c,nlevsoi)
       if(t_soisno(c,jwt(c)+1) < tfrz) then
          vwc_zwt(c) = vwc_liq(c,nlevsoi)
          do j = nlevsoi,nlevgrnd
             if(zwt(c) <= zi(c,j)) then
                smp1 = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                !smp1 = max(0._r8,smp1)
                smp1 = max(sucsat(c,nlevsoi),smp1)
                vwc_zwt(c) = watsat(c,nlevsoi)*(smp1/sucsat(c,nlevsoi))**(-1._r8/bsw(c,nlevsoi))
                ! for temperatures close to tfrz, limit vwc to total water content 
                vwc_zwt(c) = min(vwc_zwt(c), 0.5*(watsat(c,nlevsoi) + h2osoi_vol(c,nlevsoi)) )
                exit
             endif
          enddo
       endif
    end do

    ! calculate the equilibrium water content based on the water table depth
            
    do j=1,nlevsoi 
       do fc=1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if ((zwtmm(c) .le. zimm(c,j-1))) then 
             vol_eq(c,j) = watsat(c,j)
            
          ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
          ! rest of the layer

          else if ((zwtmm(c) .lt. zimm(c,j)) .and. (zwtmm(c) .gt. zimm(c,j-1))) then
             tempi = 1.0_r8
             temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             voleq1 = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
             vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
             vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
             vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
          else
             tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
             vol_eq(c,j) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
             vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
             vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
          endif
          zq(c,j) = -sucsat(c,j)*(max(vol_eq(c,j)/watsat(c,j),0.01_r8))**(-bsw(c,j))
          zq(c,j) = max(smpmin(c), zq(c,j))
       end do
    end do

    ! If water table is below soil column calculate zq for the 11th layer
    j = nlevsoi
    do fc=1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(jwt(c) == nlevsoi) then 
          tempi = 1._r8
          temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
          vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
          vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_r8)
          vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
          zq(c,j+1) = -sucsat(c,j)*(max(vol_eq(c,j+1)/watsat(c,j),0.01_r8))**(-bsw(c,j))
          zq(c,j+1) = max(smpmin(c), zq(c,j+1))
       end if
    end do

    ! Hydraulic conductivity and soil matric potential and their derivatives

    sdamp = 0._r8
    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          ! compute hydraulic conductivity based on liquid water content only

          if(origflag == 1) then
             s1 = 0.5_r8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
                  (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
          else
             s1 = 0.5_r8*(vwc_liq(c,j) + vwc_liq(c,min(nlevsoi, j+1))) / &
                  (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
          endif
          s1 = min(1._r8, s1)
          s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

          ! replace fracice with impedance factor, as in zhao 97,99
          if(origflag == 1) then
             imped(c,j)=(1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))
          else
             imped(c,j)=10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi, j+1)))))
          endif
          hk(c,j) = imped(c,j)*s1*s2
          dhkdw(c,j) = imped(c,j)*(2._r8*bsw(c,j)+3._r8)*s2* &
                       (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))

          ! compute matric potential and derivative based on liquid water content only
          if(origflag == 1) then
             s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          else
             s_node = max(vwc_liq(c,j)/watsat(c,j), 0.01_r8)
          endif
          s_node = min(1.0_r8, s_node)

          smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
          smp(c,j) = max(smpmin(c), smp(c,j))

          if(origflag == 1) then
             dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))
          else
             dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/vwc_liq(c,j)
          endif

          smp_l(c,j) = smp(c,j)
          hk_l(c,j) = hk(c,j)

       end do
    end do

    ! aquifer (11th) layer
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       zmm(c,nlevsoi+1) = 0.5*(1.e3_r8*zwt(c) + zmm(c,nlevsoi))
       if(jwt(c) < nlevsoi) then
         dzmm(c,nlevsoi+1) = dzmm(c,nlevsoi)
       else
         dzmm(c,nlevsoi+1) = (1.e3_r8*zwt(c) - zmm(c,nlevsoi))
       end if
    end do

    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qin(c,j)    = qflx_infl(c)
       den    = (zmm(c,j+1)-zmm(c,j))
       dzq    = (zq(c,j+1)-zq(c,j))
       num    = (smp(c,j+1)-smp(c,j)) - dzq
       qout(c,j)   = -hk(c,j)*num/den
       dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
       dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
       rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c) * rootr_col(c,j)
       amx(c,j) =  0._r8
       bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1(c,j)
       cmx(c,j) =  dqodw2(c,j)
    end do

    ! Nodes j=2 to j=nlevsoi-1

    do j = 2, nlevsoi - 1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          den    = (zmm(c,j) - zmm(c,j-1))
          dzq    = (zq(c,j)-zq(c,j-1))
          num    = (smp(c,j)-smp(c,j-1)) - dzq
          qin(c,j)    = -hk(c,j-1)*num/den
          dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
          dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
          den    = (zmm(c,j+1)-zmm(c,j))
          dzq    = (zq(c,j+1)-zq(c,j))
          num    = (smp(c,j+1)-smp(c,j)) - dzq
          qout(c,j)   = -hk(c,j)*num/den
          dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
          dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
          rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
          amx(c,j)    = -dqidw0(c,j)
          bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
          cmx(c,j)    =  dqodw2(c,j)
       end do
    end do

    ! Node j=nlevsoi (bottom)

    j = nlevsoi
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(j > jwt(c)) then !water table is in soil column
         den    = (zmm(c,j) - zmm(c,j-1))
         dzq    = (zq(c,j)-zq(c,j-1))
         num    = (smp(c,j)-smp(c,j-1)) - dzq
         qin(c,j)    = -hk(c,j-1)*num/den
         dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
         dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
         qout(c,j)   =  0._r8
         dqodw1(c,j) =  0._r8
         rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
         amx(c,j)    = -dqidw0(c,j)
         bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
         cmx(c,j)    =  0._r8

         ! next set up aquifer layer; hydrologically inactive
         rmx(c,j+1) = 0._r8
         amx(c,j+1) = 0._r8
         bmx(c,j+1) = dzmm(c,j+1)/dtime
         cmx(c,j+1) = 0._r8
       else ! water table is below soil column

         ! compute aquifer soil moisture as average of layer 10 and saturation
         if(origflag == 1) then
            s_node = max(0.5*(1.0_r8+h2osoi_vol(c,j)/watsat(c,j)), 0.01_r8)
         else
            s_node = max(0.5*((vwc_zwt(c)+vwc_liq(c,j))/watsat(c,j)), 0.01_r8)
         endif
         s_node = min(1.0_r8, s_node)

         ! compute smp for aquifer layer
         smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
         smp1 = max(smpmin(c), smp1)

         ! compute dsmpdw for aquifer layer
         dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

         ! first set up bottom layer of soil column
         den    = (zmm(c,j) - zmm(c,j-1))
         dzq    = (zq(c,j)-zq(c,j-1))
         num    = (smp(c,j)-smp(c,j-1)) - dzq
         qin(c,j)    = -hk(c,j-1)*num/den
         dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
         dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
         den    = (zmm(c,j+1)-zmm(c,j))
         dzq    = (zq(c,j+1)-zq(c,j))
         num    = (smp1-smp(c,j)) - dzq
         qout(c,j)   = -hk(c,j)*num/den
         dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
         dqodw2(c,j) = -( hk(c,j)*dsmpdw1 + num*dhkdw(c,j))/den

         rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
         amx(c,j) = -dqidw0(c,j)
         bmx(c,j) =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
         cmx(c,j) =  dqodw2(c,j)

         ! next set up aquifer layer; den/num unchanged, qin=qout
         qin(c,j+1)    = qout(c,j)
         dqidw0(c,j+1) = -(-hk(c,j)*dsmpdw(c,j) + num*dhkdw(c,j))/den
         dqidw1(c,j+1) = -( hk(c,j)*dsmpdw1   + num*dhkdw(c,j))/den
         qout(c,j+1)   =  0._r8  ! zero-flow bottom boundary condition
         dqodw1(c,j+1) =  0._r8  ! zero-flow bottom boundary condition
         rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
         amx(c,j+1) = -dqidw0(c,j+1)
         bmx(c,j+1) =  dzmm(c,j+1)/dtime - dqidw1(c,j+1) + dqodw1(c,j+1)
         cmx(c,j+1) =  0._r8
       endif
    end do

    ! Solve for dwat

    jtop(bounds%begc : bounds%endc) = 1
    call Tridiagonal(bounds, 1, nlevsoi+1, &
         jtop(bounds%begc:bounds%endc), &
         num_hydrologyc, filter_hydrologyc, &
         amx(bounds%begc:bounds%endc, :), &
         bmx(bounds%begc:bounds%endc, :), &
         cmx(bounds%begc:bounds%endc, :), &
         rmx(bounds%begc:bounds%endc, :), &
         dwat2(bounds%begc:bounds%endc, :) )

    ! Renew the mass of liquid water
    ! also compute qcharge from dwat in aquifer layer
    ! update in drainage for case jwt < nlevsoi

    do fc = 1,num_hydrologyc
       c = filter_hydrologyc(fc)
       do j = 1, nlevsoi
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat2(c,j)*dzmm(c,j)
       end do

       ! calculate qcharge for case jwt < nlevsoi
       if(jwt(c) < nlevsoi) then
          wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat
          
          ! Recharge rate qcharge to groundwater (positive to aquifer)
          s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_r8)
          s1 = min(1._r8, s_node)
          
          !scs: this is the expression for unsaturated hk
          ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
               *s1**(2._r8*bsw(c,jwt(c)+1)+3._r8)
          
          ! Recharge rate qcharge to groundwater (positive to aquifer)
          smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
          wh      = smp1 - zq(c,max(1,jwt(c)))
          
          !scs: original formulation
          if(jwt(c) == 0) then
             qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._r8)
          else
             !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8)
             !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8*2.0)
          endif
          
          ! To limit qcharge  (for the first several timesteps)
          qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
          qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
       else
          ! if water table is below soil column, compute qcharge from dwat2(11)
          qcharge(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtime
       endif
    end do

! compute the water deficit and reset negative liquid water content
!  Jinyun Tang
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_deficit(c) = 0._r8
       do j = 1, nlevsoi
          if(h2osoi_liq(c,j)<0._r8)then
             qflx_deficit(c) = qflx_deficit(c) - h2osoi_liq(c,j)
          endif
       enddo
    enddo

    end associate 
   end subroutine SoilWater

   !-----------------------------------------------------------------------
   subroutine WaterTable(bounds, num_hydrologyc, filter_hydrologyc, &
        num_urbanc, filter_urbanc)
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
     use clmtype
     use clm_time_manager, only : get_step_size
     use clm_varcon  , only : pondmx, tfrz, icol_roof, icol_road_imperv, &
                              watmin,denice,denh2o
     use clm_varpar  , only : nlevsoi
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds  ! bounds
     integer , intent(in) :: num_hydrologyc               ! number of column soil points in column filter
     integer , intent(in) :: num_urbanc                   ! number of column urban points in column filter
     integer , intent(in) :: filter_urbanc(:)     ! column filter for urban points
     integer , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i                 ! indices
     real(r8) :: dtime                    ! land model time step (sec)
     real(r8) :: xs(bounds%begc:bounds%endc)              ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi)  ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)             ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)        ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)        ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)             ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)             ! excess soil water above saturation at layer i (mm)
     real(r8) :: rous                     ! aquifer yield (-)
     real(r8) :: wh                       ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: ws                       ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                   ! soil wetness (-)
     real(r8) :: dzsum                    ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum               ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)    ! fractional impermeability of soil layers (-)
     real(r8) :: ka                       ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                      ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq     ! available soil liquid water in a layer
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

   associate(& 
   h2osfc                              =>    cws%h2osfc                                  , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   frac_h2osfc                         =>    cps%frac_h2osfc                             , & ! Input:  [real(r8) (:)]                                                    
   frost_table                         =>    cws%frost_table                             , & ! Input:  [real(r8) (:)]  frost table depth (m)                             
   zwt_perched                         =>    cws%zwt_perched                             , & ! Input:  [real(r8) (:)]  perched water table depth (m)                     
   qflx_drain_perched                  =>    cwf%qflx_drain_perched                      , & ! Input:  [real(r8) (:)]  perched wt sub-surface runoff (mm H2O /s)         
   qflx_deficit                        =>    cwf%qflx_deficit                            , & ! Input:  [real(r8) (:)]  water deficit to keep non-negative liquid water content
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type index                                  
   snl                                 =>    cps%snl                                     , & ! Input:  [integer (:)]  number of snow layers                              
   dz                                  =>    cps%dz                                      , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   bsw                                 =>    cps%bsw                                     , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   t_soisno                            =>    ces%t_soisno                                , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)                       
   hksat                               =>    cps%hksat                                   , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s)
   sucsat                              =>    cps%sucsat                                  , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   z                                   =>    cps%z                                       , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   zi                                  =>    cps%zi                                      , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)           
   watsat                              =>    cps%watsat                                  , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   zwt                                 =>    cws%zwt                                     , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   wa                                  =>    cws%wa                                      , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)              
   qcharge                             =>    cws%qcharge                                 , & ! Input:  [real(r8) (:)]  aquifer recharge rate (mm/s)                      
   eff_porosity                        =>    cps%eff_porosity                            , & ! Input:  [real(r8) (:,:)]  effective porosity = porosity - vol_ice         
   qflx_dew_grnd                       =>    pwf_a%qflx_dew_grnd                         , & ! Input:  [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]      
   qflx_dew_snow                       =>    pwf_a%qflx_dew_snow                         , & ! Input:  [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]    
   qflx_sub_snow                       =>    pwf_a%qflx_sub_snow                         , & ! Input:  [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]   
   qflx_drain                          =>    cwf%qflx_drain                              , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)                    
   qflx_rsub_sat                       =>    cwf%qflx_rsub_sat                           , & ! Output: [real(r8) (:)]  soil saturation excess [mm h2o/s]                 
   h2osoi_liq                          =>    cws%h2osoi_liq                              , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_ice                          =>    cws%h2osoi_ice                              , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   h2osoi_vol                          =>    cws%h2osoi_vol                                & ! Input:  [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
   )

    ! Get time step

    dtime = get_step_size()

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          dzmm(c,j) = dz(c,j)*1.e3_r8
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_drain(c)    = 0._r8
       qflx_rsub_sat(c) = 0._r8
       qflx_drain_perched(c)  = 0._r8       
    end do

    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       ! allow jwt to equal zero when zwt is in top layer
       do j = 1,nlevsoi
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

!scs: use analytical expression for aquifer specific yield
       rous = watsat(c,nlevsoi) &
            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
       rous=max(rous,0.02_r8)

!--  water table is below the soil column  --------------------------------------
       if(jwt(c) == nlevsoi) then             
          wa(c)  = wa(c) + qcharge(c)  * dtime
          zwt(c)     = zwt(c) - (qcharge(c)  * dtime)/1000._r8/rous
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
             do j = jwt(c)+1, nlevsoi
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
          jwt(c) = nlevsoi
          do j = 1,nlevsoi
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
          h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
          h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
          if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
             qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
             h2osoi_ice(c,1) = 0._r8
          else
             h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtime
          end if
       end if
    end do


    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

       if (ctype(c) == icol_roof .or. ctype(c) == icol_road_imperv) then
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
    subroutine Drainage(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc)
      !
      ! !DESCRIPTION:
      ! Calculate subsurface drainage
      !
      ! !USES:
      use clmtype
      use clm_time_manager, only : get_step_size
      use clm_varcon  , only : pondmx, tfrz, icol_roof, icol_road_imperv, icol_road_perv, watmin,rpi
      use clm_varpar  , only : nlevsoi,nlevgrnd
      use clm_varcon  , only : secspday,nlvic
      use clm_varpar  , only : nlayer, nlayert
      use CLMVICMapMod , only : CLMVICMap
      use abortutils  , only : endrun
      !
      ! !ARGUMENTS:
      implicit none
      type(bounds_type), intent(in) :: bounds  ! bounds
      integer , intent(in) :: num_hydrologyc               ! number of column soil points in column filter
      integer , intent(in) :: num_urbanc                   ! number of column urban points in column filter
      integer , intent(in) :: filter_urbanc(:)     ! column filter for urban points
      integer , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
      !
      ! !LOCAL VARIABLES:
      character(len=32) :: subname = 'Drainage'  ! subroutine name
      integer  :: c,j,fc,i                 ! indices
      real(r8) :: dtime                    ! land model time step (sec)
      real(r8) :: xs(bounds%begc:bounds%endc)              ! water needed to bring soil moisture to watmin (mm)
      real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi)  ! layer thickness (mm)
      integer  :: jwt(bounds%begc:bounds%endc)             ! index of the soil layer right above the water table (-)
      real(r8) :: rsub_bot(bounds%begc:bounds%endc)        ! subsurface runoff - bottom drainage (mm/s)
      real(r8) :: rsub_top(bounds%begc:bounds%endc)        ! subsurface runoff - topographic control (mm/s)
      real(r8) :: fff(bounds%begc:bounds%endc)             ! decay factor (m-1)
      real(r8) :: xsi(bounds%begc:bounds%endc)             ! excess soil water above saturation at layer i (mm)
      real(r8) :: xsia(bounds%begc:bounds%endc)            ! available pore space at layer i (mm)
      real(r8) :: xs1(bounds%begc:bounds%endc)             ! excess soil water above saturation at layer 1 (mm)
      real(r8) :: smpfz(1:nlevsoi)         ! matric potential of layer right above water table (mm)
      real(r8) :: wtsub                    ! summation of hk*dzmm for layers below water table (mm**2/s)
      real(r8) :: rous                     ! aquifer yield (-)
      real(r8) :: wh                       ! smpfz(jwt)-z(jwt) (mm)
      real(r8) :: wh_zwt                   ! water head at the water table depth (mm)
      real(r8) :: ws                       ! summation of pore space of layers below water table (mm)
      real(r8) :: s_node                   ! soil wetness (-)
      real(r8) :: dzsum                    ! summation of dzmm of layers below water table (mm)
      real(r8) :: icefracsum               ! summation of icefrac*dzmm of layers below water table (-)
      real(r8) :: fracice_rsub(bounds%begc:bounds%endc)    ! fractional impermeability of soil layers (-)
      real(r8) :: ka                       ! hydraulic conductivity of the aquifer (mm/s)
      real(r8) :: dza                      ! fff*(zwt-z(jwt)) (-)
      real(r8) :: available_h2osoi_liq     ! available soil liquid water in a layer
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


   associate(& 
   h2osfc                              =>    cws%h2osfc                                  , & ! Input:  [real(r8) (:)]  surface water (mm)                                
   frac_h2osfc                         =>    cps%frac_h2osfc                             , & ! Input:  [real(r8) (:)]                                                    
   topo_ndx                            =>    cps%topo_ndx                                , & ! Input:  [real(r8) (:)]  topographic index                                 
   topo_slope                          =>    cps%topo_slope                              , & ! Input:  [real(r8) (:)]  topographic slope                                 
   frost_table                         =>    cws%frost_table                             , & ! Input:  [real(r8) (:)]  frost table depth (m)                             
   zwt_perched                         =>    cws%zwt_perched                             , & ! Input:  [real(r8) (:)]  perched water table depth (m)                     
   qflx_drain_perched                  =>    cwf%qflx_drain_perched                      , & ! Input:  [real(r8) (:)]  perched wt sub-surface runoff (mm H2O /s)         
   clandunit                           =>   col%landunit                                 , & ! Input:  [integer (:)]  column's landunit                                  
   ltype                               =>    lun%itype                                   , & ! Input:  [integer (:)]  landunit type                                      
   urbpoi                              =>    lun%urbpoi                                  , & ! Input:  [logical (:)]  true => landunit is an urban point                 
   ctype                               =>    col%itype                                   , & ! Input:  [integer (:)]  column type index                                  
   snl                                 =>    cps%snl                                     , & ! Input:  [integer (:)]  number of snow layers                              
   dz                                  =>    cps%dz                                      , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   bsw                                 =>    cps%bsw                                     , & ! Input:  [real(r8) (:,:)]  Clapp and Hornberger "b"                        
   t_soisno                            =>    ces%t_soisno                                , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)                       
   hksat                               =>    cps%hksat                                   , & ! Input:  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s)
   sucsat                              =>    cps%sucsat                                  , & ! Input:  [real(r8) (:,:)]  minimum soil suction (mm)                       
   z                                   =>    cps%z                                       , & ! Input:  [real(r8) (:,:)]  layer depth (m)                                 
   zi                                  =>    cps%zi                                      , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)           
   watsat                              =>    cps%watsat                                  , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
   hkdepth                             =>    cps%hkdepth                                 , & ! Input:  [real(r8) (:)]  decay factor (m)                                  
   icefrac                             =>    cps%icefrac                                 , & ! Input:  [real(r8) (:,:)] fraction of ice in layer                         
   zwt                                 =>    cws%zwt                                     , & ! Input:  [real(r8) (:)]  water table depth (m)                             
   wa                                  =>    cws%wa                                      , & ! Input:  [real(r8) (:)]  water in the unconfined aquifer (mm)              
   qcharge                             =>    cws%qcharge                                 , & ! Input:  [real(r8) (:)]  aquifer recharge rate (mm/s)                      
   eff_porosity                        =>    cps%eff_porosity                            , & ! Input:  [real(r8) (:,:)]  effective porosity = porosity - vol_ice         
   qflx_snwcp_liq                      =>    pwf_a%qflx_snwcp_liq                        , & ! Input:  [real(r8) (:)]  excess rainfall due to snow capping (mm H2O /s) [+]
   qflx_snwcp_ice                      =>    pwf_a%qflx_snwcp_ice                        , & ! Input:  [real(r8) (:)]  excess snowfall due to snow capping (mm H2O /s) [+]
   qflx_dew_grnd                       =>    pwf_a%qflx_dew_grnd                         , & ! Input:  [real(r8) (:)]  ground surface dew formation (mm H2O /s) [+]      
   qflx_dew_snow                       =>    pwf_a%qflx_dew_snow                         , & ! Input:  [real(r8) (:)]  surface dew added to snow pack (mm H2O /s) [+]    
   qflx_sub_snow                       =>    pwf_a%qflx_sub_snow                         , & ! Input:  [real(r8) (:)]  sublimation rate from snow pack (mm H2O /s) [+]   
   qflx_drain                          =>    cwf%qflx_drain                              , & ! Output: [real(r8) (:)]  sub-surface runoff (mm H2O /s)                    
   qflx_qrgwl                          =>    cwf%qflx_qrgwl                              , & ! Output: [real(r8) (:)]  qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
   qflx_rsub_sat                       =>    cwf%qflx_rsub_sat                           , & ! Output: [real(r8) (:)]  soil saturation excess [mm h2o/s]                 
   eflx_impsoil                        =>    cef%eflx_impsoil                            , & ! Output: [real(r8) (:)]  implicit evaporation for soil temperature equation
   Dsmax                               =>    cps%dsmax                                   , & ! Input:  [real(r8) (:)] max. velocity of baseflow (mm/day)
   max_moist                           =>    cps%max_moist                               , & ! Input:  [real(r8) (:,:)] maximum soil moisture (ice + liq)
   Ds                                  =>    cps%ds                                      , & ! Input:  [real(r8) (:)] fracton of Dsmax where non-linear baseflow begins
   Wsvic                               =>    cps%Wsvic                                   , & ! Input:  [real(r8) (:)] fraction of maximum soil moisutre where non-liear base flow occurs
   c_param                             =>    cps%c_param                                 , & ! Input:  [real(r8) (:)] baseflow exponent (Qb)                             
   depth                               =>    cps%depth                                   , & ! Input:  [real(r8) (:,:)] VIC soil depth                                   
   moist                               =>    cws%moist                                   , & ! Input:  [real(r8) (:,:)] soil layer moisture (mm)                         
   ice                                 =>    cws%ice                                     , & ! Input:  [real(r8) (:,:)] soil layer moisture (mm)                         
   hk_l                                =>    cws%hk_l                                    , & ! Input:  [real(r8) (:,:)] hydraulic conductivity (mm/s)                    
   h2osoi_liq                          =>    cws%h2osoi_liq                              , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)                            
   h2osoi_ice                          =>    cws%h2osoi_ice                                & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)                                
   )

    ! Get time step

    dtime = get_step_size()

    ! Convert layer thicknesses from m to mm

    do j = 1,nlevsoi
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
       qflx_drain(c)    = 0._r8 
       rsub_bot(c)      = 0._r8
       qflx_rsub_sat(c) = 0._r8
       rsub_top(c)      = 0._r8
       fracice_rsub(c)  = 0._r8
    end do

    ! The layer index of the first unsaturated layer, i.e., the layer right above
    ! the water table

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       jwt(c) = nlevsoi
       ! allow jwt to equal zero when zwt is in top layer
       do j = 1,nlevsoi
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
        q_perch_max = 1.e-5_r8 * sin(topo_slope(c) * (rpi/180._r8))

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
             rsub_top_max = 10._r8 * sin((rpi/180.) * topo_slope(c))
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
          rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*zwt(c))
       end if

       ! use analytical expression for aquifer specific yield
       rous = watsat(c,nlevsoi) &
            * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevsoi))**(-1./bsw(c,nlevsoi)))
       rous=max(rous,0.02_r8)

!--  water table is below the soil column  --------------------------------------
       if(jwt(c) == nlevsoi) then             
          wa(c)  = wa(c) - rsub_top(c) * dtime
          zwt(c)     = zwt(c) + (rsub_top(c) * dtime)/1000._r8/rous
          h2osoi_liq(c,nlevsoi) = h2osoi_liq(c,nlevsoi) + max(0._r8,(wa(c)-5000._r8))
          wa(c)  = min(wa(c), 5000._r8)
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
                do j = (nlvic(1)+nlvic(2)+1), nlevsoi
                   wtsub_vic = wtsub_vic + hk_l(c,j)*dzmm(c,j)
                end do
                
                do j = (nlvic(1)+nlvic(2)+1), nlevsoi
                   rsub_top_layer=max(rsub_top_tot, rsub_top_tot*hk_l(c,j)*dzmm(c,j)/wtsub_vic)
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                   rsub_top_tot = rsub_top_tot - rsub_top_layer
                end do
             else
                do j = jwt(c)+1, nlevsoi
                   ! use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)
                   
                   rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
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
            zwt(c) = zwt(c) - rsub_top_tot/1000._r8/rous
            wa(c) = wa(c) + rsub_top_tot
          endif

!-- recompute jwt  ---------------------------------------------------------
          ! allow jwt to equal zero when zwt is in top layer
          jwt(c) = nlevsoi
          do j = 1,nlevsoi
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

    do j = nlevsoi,2,-1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
          h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
          h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
       end do
    end do

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)

       !scs: watmin addition to fix water balance errors
       xs1(c)          = max(max(h2osoi_liq(c,1)-watmin,0._r8)- &
          max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)-watmin)),0._r8)
       h2osoi_liq(c,1) = h2osoi_liq(c,1) - xs1(c)

       if (urbpoi(clandunit(c))) then
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
       ! add in ice check
       xs1(c)          = max(max(h2osoi_ice(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1))),0._r8)
       h2osoi_ice(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1)), h2osoi_ice(c,1))
       qflx_snwcp_ice(c) = qflx_snwcp_ice(c) + xs1(c) / dtime
    end do

    ! Limit h2osoi_liq to be greater than or equal to watmin.
    ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water

    do j = 1, nlevsoi-1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
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
    j = nlevsoi
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (h2osoi_liq(c,j) < watmin) then
          xs(c) = watmin-h2osoi_liq(c,j)
          searchforwater: do i = nlevsoi-1, 1, -1
            available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._r8)
            if (available_h2osoi_liq .ge. xs(c)) then
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

       ! Implicit evaporation term is now zero

       eflx_impsoil(c) = 0._r8


    end do

    ! No drainage for urban columns (except for pervious road as computed above)

    do fc = 1, num_urbanc
       c = filter_urbanc(fc)
       if (ctype(c) /= icol_road_perv) then
         qflx_drain(c) = 0._r8
         ! This must be done for roofs and impervious road (walls will be zero)
         qflx_qrgwl(c) = qflx_snwcp_liq(c)
         eflx_impsoil(c) = 0._r8
       end if



    end do

    end associate 
   end subroutine Drainage

end module SoilHydrologyMod
