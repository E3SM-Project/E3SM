module SoilMoistStressMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates soil moisture stress for plant gpp and transpiration
  !
  ! After discussion with other developers, I have now removed all functions that
  ! return array, and decalared all variables that will be modified as intent(inout).
  ! The initialization will be done whenever the variable is initialized. This avoids
  ! code crash when initialization is not done appropriately, and make the code safer
  ! during the long-term maintenance
  ! Created by Jinyun Tang, Feb., 2014 
  !
  use ColumnDataType   , only : col_es, col_ws
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calc_root_moist_stress
  public :: calc_effective_soilporosity
  public :: calc_effective_snowporosity
  public :: calc_volumetric_h2oliq
  public :: set_perchroot_opt
  public :: init_root_moist_stress
  !
  ! !PRIVATE DATA MEMBERS:
  integer ::   root_moist_stress_method
  integer, parameter :: moist_stress_clm_default  = 0  !default method for calculating root moisture stress
  logical,  private :: perchroot     = .false.  ! true => btran is based only on unfrozen soil levels
  logical,  private :: perchroot_alt = .false.  ! true => btran is based on active layer (defined over two years); 
  !--------------------------------------------------------------------------------

contains

  !--------------------------------------------------------------------------------
  subroutine init_root_moist_stress()
    !
    !DESCRIPTION
    !specify the method to compute root soil moisture stress
    !
    implicit none

    root_moist_stress_method = moist_stress_clm_default   
  end subroutine init_root_moist_stress

  !--------------------------------------------------------------------------------
  subroutine set_perchroot_opt(perchroot_global, perchroot_alt_global)
    !
    !DESCRIPTIONS
    !set up local perchroot logical switches, in the future, this wil be
    !read in as namelist
    !
    ! !ARGUMENTS:
    implicit none
    logical, intent(in) :: perchroot_global
    logical, intent(in) :: perchroot_alt_global
    !------------------------------------------------------------------------------

    perchroot = perchroot_global
    perchroot_alt = perchroot_alt_global

  end subroutine set_perchroot_opt

  !--------------------------------------------------------------------------------
  subroutine calc_effective_soilporosity(bounds, ubj, numf, filter, &
       watsat, h2osoi_ice, denice, eff_por)
    !
    ! !DESCRIPTIONS
    ! compute the effective soil porosity
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use decompMod      , only : bounds_type
    use ColumnType     , only : col_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                          ! bounds
    integer           , intent(in)    :: ubj                             ! lbinning level indices
    integer           , intent(in)    :: numf                            ! filter dimension
    integer           , intent(in)    :: filter(:)                       ! filter  
    real(r8)          , intent(in)    :: watsat( bounds%begc: , 1: )     ! soil porosity
    real(r8)          , intent(in)    :: h2osoi_ice( bounds%begc: , 1: ) ! ice water content, kg H2o/m2
    real(r8)          , intent(in)    :: denice                          ! ice density, kg/m3
    real(r8)          , intent(inout) :: eff_por( bounds%begc: ,1: )     ! effective porosity
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, fc                                             !indices
    real(r8):: vol_ice    !volumetric ice
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(watsat)     == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osoi_ice) == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(eff_por)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))

    !main calculation loop
    !it assumes the soil layers start from 1
    do j = 1, ubj
       do fc = 1, numf
          c = filter(fc)
          !compute the volumetric ice content
          vol_ice=min(watsat(c,j), h2osoi_ice(c,j)/(denice*col_pp%dz(c,j)))

          !compute the maximum soil space to fill liquid water and air
          eff_por(c,j) = watsat(c,j) - vol_ice
       enddo
    enddo
  end subroutine calc_effective_soilporosity

  !--------------------------------------------------------------------------------
  subroutine calc_effective_snowporosity(bounds, lbj, jtop, numf, filter, &
       h2osoi_ice, denice, eff_por)
    !
    ! !DESCRIPTIONS
    ! compute the effective porosity snow
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use decompMod      , only : bounds_type
    use shr_log_mod    , only : errMsg => shr_log_errMsg    
    use ColumnType     , only : col_pp
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds                            !bounds
    integer           , intent(in)    :: lbj                               !ubing level indices
    integer           , intent(in)    :: jtop( bounds%begc: )              !top level for each column [col]    
    integer           , intent(in)    :: numf                              !filter dimension
    integer           , intent(in)    :: filter(:)                         !filter  
    real(r8)          , intent(in)    :: h2osoi_ice( bounds%begc: , lbj: ) !ice water content, kg H2o/m2
    real(r8)          , intent(in)    :: denice                            !ice density, kg/m3
    real(r8)          , intent(inout) :: eff_por( bounds%begc: ,lbj: )     !returning effective porosity
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc    !indices
    integer  :: ubj
    real(r8) :: vol_ice     !volumetric ice
    !------------------------------------------------------------------------------

    ubj = 0

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(jtop)       == (/bounds%endc/))     , errMsg(__FILE__, __LINE__)) 
    SHR_ASSERT_ALL((ubound(h2osoi_ice) == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(eff_por)    == (/bounds%endc,0/))   , errMsg(__FILE__, __LINE__))

    !main calculation loop

    !it assumes snow layer ends at 0
    do j = lbj,0
       do fc = 1, numf
          c = filter(fc)
          if (j>=jtop(c)) then
             !compute the volumetric ice content       
             vol_ice=min(1._r8, h2osoi_ice(c,j)/(denice*col_pp%dz(c,j)))

             !compute the maximum snow void space to fill liquid water and air         
             eff_por(c,j) = 1._r8 - vol_ice
          endif
       enddo
    enddo

  end subroutine calc_effective_snowporosity

  !--------------------------------------------------------------------------------
  subroutine calc_volumetric_h2oliq(bounds, jtop, lbj, ubj, numf, filter,&
       eff_porosity, h2osoi_liq, denh2o, vol_liq)
    !
    ! !DESCRIPTIONS
    ! compute the volumetric liquid water content
    !
    !
    ! !USES
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_log_mod    , only : errMsg => shr_log_errMsg  
    use decompMod      , only : bounds_type
    use ColumnType     , only : col_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                             ! bounds
    integer           , intent(in)    :: jtop( bounds%begc: )               ! top level for each column [col]  
    integer           , intent(in)    :: lbj, ubj                           ! lbinning and ubing level indices
    integer           , intent(in)    :: numf                               ! filter dimension
    integer           , intent(in)    :: filter(:)                          ! filter    
    real(r8)          , intent(in)    :: eff_porosity(bounds%begc: , lbj: ) ! effective soil porosity
    real(r8)          , intent(in)    :: h2osoi_liq(bounds%begc: , lbj: )   ! liquid water content [kg H2o/m2]
    real(r8)          , intent(in)    :: denh2o                             ! water density [kg/m3]
    real(r8)          , intent(inout) :: vol_liq(bounds%begc: , lbj: )      ! volumetric liquid water content  
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, fc  ! indices  
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes  
    SHR_ASSERT_ALL((ubound(jtop)         == (/bounds%endc/))     , errMsg(__FILE__, __LINE__)) 
    SHR_ASSERT_ALL((ubound(h2osoi_liq)   == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(eff_porosity) == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(vol_liq)      == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))  

    !main calculation loop
    do j = lbj, ubj
       do fc = 1, numf
          c = filter(fc)
          if(j>=jtop(c))then
             !volume of liquid is no greater than effective void space
             vol_liq(c,j) = min(eff_porosity(c,j), h2osoi_liq(c,j)/(col_pp%dz(c,j)*denh2o))
          endif
       enddo
    enddo

  end subroutine calc_volumetric_h2oliq

  !--------------------------------------------------------------------------------
  subroutine normalize_unfrozen_rootfr(bounds, ubj, fn, filterp, &
       canopystate_vars, soilstate_vars, temperature_vars, rootfr_unf)
    !
    ! !DESCRIPTIONS
    ! normalize root fraction for total unfrozen depth 
    !
    ! !USES
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use elm_varcon      , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment 
    use decompMod       , only : bounds_type
    use CanopyStateType , only : canopystate_type
    use EnergyFluxType  , only : energyflux_type
    use TemperatureType , only : temperature_type
    use SoilStateType   , only : soilstate_type
    use WaterSTateType  , only : waterstate_type
    use SimpleMathMod   , only : array_normalization
    use VegetationType       , only : veg_pp
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds                                     !bounds
    integer                , intent(in)    :: ubj                                        !ubinning level indices
    integer                , intent(in)    :: fn                                         !filter dimension
    integer                , intent(in)    :: filterp(:)                                 !filter
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(soilstate_type)   , intent(in)    :: soilstate_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    real(r8)               , intent(inout) :: rootfr_unf(bounds%begp:bounds%endp, 1:ubj) !normalized root fraction in unfrozen layers
    !
    ! !LOCAL VARIABLES:
    !real(r8) :: rootsum(bounds%begp:bounds%endp)  
    integer :: p, c, j, f  !indices  
    !------------------------------------------------------------------------------

    associate(                                                               &
         rootfr               => soilstate_vars%rootfr_patch               , & ! Input:  [real(r8)  (:,:) ]  fraction of roots in each soil layer

         t_soisno             => col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)                    

         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [real(r8) (:)   ]  prior year maximum annual depth of thaw                               
         altmax_indx          => canopystate_vars%altmax_indx_col            & ! Input:  [real(r8) (:)   ]  maximum annual depth of thaw                                          
         )

      ! main calculation loop  
      ! Initialize rootfr_unf to zero.
      ! I found it necessary to ensure the pgi compiler not
      ! to complain with float point exception. However, it raises a question how
      ! to make sure those values that are initialized with nan or spval are not reset
      ! to zero within similar coding style. Jinyun Tang, May 23, 2014.
      
      ! Define rootfraction for unfrozen soil only
      if (perchroot .or. perchroot_alt) then
         if (perchroot_alt) then
            ! use total active layer (defined ass max thaw depth for current and prior year)
            do j = 1, ubj
               do f = 1, fn
                  p = filterp(f)
                  c = veg_pp%column(p)

                  if ( j <= max(altmax_lastyear_indx(c), altmax_indx(c), 1) )then
                     rootfr_unf(p,j) = rootfr(p,j)
                  else
                     rootfr_unf(p,j) = 0._r8
                  end if
               end do
            end do
         else
            ! use instantaneous temperature
            do j = 1, ubj
               do f = 1, fn
                  p = filterp(f)
                  c = veg_pp%column(p)

                  if (t_soisno(c,j) >= tfrz) then
                     rootfr_unf(p,j) = rootfr(p,j)
                  else
                     rootfr_unf(p,j) = 0._r8
                  end if
               end do
            end do

         end if ! perchroot_alt          
      end if ! perchroot

      !normalize the root fraction for each pft
      call array_normalization(bounds%begp, bounds%endp, 1, ubj, &
           fn, filterp, rootfr_unf(bounds%begp:bounds%endp, 1:ubj))

    end associate        

  end subroutine normalize_unfrozen_rootfr
  
  !--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress_clm45default(bounds, &
       nlevgrnd, fn, filterp, rootfr_unf, &
       temperature_vars, soilstate_vars, energyflux_vars, waterstate_vars, &
       soil_water_retention_curve) 
    !
    ! DESCRIPTIONS
    ! compute the root water stress using the default clm45 approach
    !
    ! USES
    use shr_kind_mod         , only : r8 => shr_kind_r8  
    use shr_log_mod          , only : errMsg => shr_log_errMsg
    use decompMod            , only : bounds_type
    use elm_varcon           , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment
    use VegetationPropertiesType     , only : veg_vp
    use TemperatureType      , only : temperature_type
    use SoilStateType        , only : soilstate_type
    use EnergyFluxType       , only : energyflux_type
    use WaterSTateType       , only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use VegetationType            , only : veg_pp
    use clm_varctl       , only : use_hydrstress
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds                         !bounds
    integer                , intent(in)    :: nlevgrnd                       !number of vertical layers
    integer                , intent(in)    :: fn                             !number of filters
    integer                , intent(in)    :: filterp(:)                     !filter array          
    real(r8)               , intent(in)    :: rootfr_unf(bounds%begp: , 1: ) 
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8) :: smp_node, s_node  !temporary variables
    real(r8) :: smp_node_lf       !temporary variable
    integer :: p, f, j, c, l      !indices
    !------------------------------------------------------------------------------

    ! Enforce expected array sizes   
    SHR_ASSERT_ALL((ubound(rootfr_unf) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))  

    associate(                                                &
         smpso         => veg_vp%smpso                  , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal opening (mm)                    
         smpsc         => veg_vp%smpsc                  , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal closure (mm)                    
         tc_stress     => veg_vp%tc_stress              , & ! Input:  [real(r8)       ]  critical soil temperature for soil water stress (C)    
         t_soisno      => col_es%t_soisno     , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)                    

         watsat        => soilstate_vars%watsat_col         , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)                     
         sucsat        => soilstate_vars%sucsat_col         , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                        (constant)                                        
         bsw           => soilstate_vars%bsw_col            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                         (constant)                                        
         eff_porosity  => soilstate_vars%eff_porosity_col   , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
         rootfr        => soilstate_vars%rootfr_patch       , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer
         rootr         => soilstate_vars%rootr_patch        , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer                      

         btran         => energyflux_vars%btran_patch       , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1) (integrated soil water stress)
         btran2        => energyflux_vars%btran2_patch      , & ! Output: [real(r8) (:)   ]  integrated soil water stress square
         rresis        => energyflux_vars%rresis_patch      , & ! Output: [real(r8) (:,:) ]  root soil water stress (resistance) by layer (0-1)  (nlevgrnd)                          

         h2osoi_vol    => col_ws%h2osoi_vol    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liqvol => col_ws%h2osoi_liqvol   & ! Output: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR
         ) 

      do j = 1,nlevgrnd
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            l = veg_pp%landunit(p)

            ! Root resistance factors
            ! rootr effectively defines the active root fraction in each layer      
            if (h2osoi_liqvol(c,j) .le. 0._r8 .or. t_soisno(c,j) .le. tfrz + tc_stress) then
               rootr(p,j) = 0._r8
            else
               s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)

               !smp_node = max(smpsc(veg_pp%itype(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))
               call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node)
               smp_node = max(smpsc(veg_pp%itype(p)), smp_node)

               rresis(p,j) = min( (eff_porosity(c,j)/watsat(c,j))* &
                    (smp_node - smpsc(veg_pp%itype(p))) / (smpso(veg_pp%itype(p)) - smpsc(veg_pp%itype(p))), 1._r8)


               if (.not. (perchroot .or. perchroot_alt) ) then
                  rootr(p,j) = rootfr(p,j)*rresis(p,j)
               else
                  rootr(p,j) = rootfr_unf(p,j)*rresis(p,j)
               end if

               !it is possible to further separate out a btran function, but I will leave it for the moment, jyt
               if( .not. use_hydrstress ) then
                 btran(p)    = btran(p) + max(rootr(p,j),0._r8)
               endif

               !smp_node_lf = max(smpsc(veg_pp%itype(p)), -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j)))
               s_node = h2osoi_vol(c,j)/watsat(c,j)

               call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_node_lf)

               !smp_node_lf =  -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j))
               smp_node_lf = max(smpsc(veg_pp%itype(p)), smp_node_lf) 
               btran2(p)   = btran2(p) +rootfr(p,j)*min((smp_node_lf - smpsc(veg_pp%itype(p))) / &
                    (smpso(veg_pp%itype(p)) - smpsc(veg_pp%itype(p))), 1._r8)
            endif
         end do
      end do

      ! Normalize root resistances to get layer contribution to ET
      do j = 1,nlevgrnd
         do f = 1, fn
            p = filterp(f)
            if (btran(p) > btran0) then
               rootr(p,j) = rootr(p,j)/btran(p)
            else
               rootr(p,j) = 0._r8
            end if
         end do
      end do
    end associate

  end subroutine calc_root_moist_stress_clm45default

  !--------------------------------------------------------------------------------
  subroutine calc_root_moist_stress(bounds, nlevgrnd, fn, filterp, &
       canopystate_vars, energyflux_vars,  soilstate_vars, temperature_vars, &
       waterstate_vars, soil_water_retention_curve)
    !
    ! DESCRIPTIONS
    ! compute the root water stress using different approaches
    !
    ! USES
    use shr_kind_mod    , only : r8 => shr_kind_r8  
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use elm_varcon      , only : tfrz      !temperature where water freezes [K], this is taken as constant at the moment 
    use decompMod       , only : bounds_type
    use CanopyStateType , only : canopystate_type
    use EnergyFluxType  , only : energyflux_type
    use TemperatureType , only : temperature_type
    use SoilStateType   , only : soilstate_type
    use WaterSTateType  , only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use abortutils      , only : endrun       
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in)    :: bounds   !bounds
    integer                , intent(in)    :: nlevgrnd
    integer                , intent(in)    :: fn
    integer                , intent(in)    :: filterp(:)
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer :: p, f, j, c, l                                   ! indices
    real(r8) :: smp_node, s_node                               ! temporary variables
    real(r8) :: rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd) ! Rootfraction defined for unfrozen layers only.
    character(len=32) :: subname = 'calc_root_moist_stress'    ! subroutine name
    !------------------------------------------------------------------------------

    !define normalized rootfraction for unfrozen soil
    !define normalized rootfraction for unfrozen soil
    rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd) = 0._r8

    call normalize_unfrozen_rootfr(bounds,  &
         ubj = nlevgrnd,                    &
         fn = fn,                           &
         filterp = filterp,                 &
         canopystate_vars=canopystate_vars, &
         soilstate_vars=soilstate_vars,     &
         temperature_vars=temperature_vars, & 
         rootfr_unf=rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd))

    !suppose h2osoi_liq, eff_porosity are already computed somewhere else

    select case (root_moist_stress_method)
       !add other methods later
    case (moist_stress_clm_default)

       call calc_root_moist_stress_clm45default(bounds, &
            nlevgrnd = nlevgrnd,                        &
            fn = fn,                                    &
            filterp = filterp,                          &
            energyflux_vars=energyflux_vars,            &
            temperature_vars=temperature_vars,          &
            soilstate_vars=soilstate_vars,              &
            waterstate_vars=waterstate_vars,            &
            rootfr_unf=rootfr_unf(bounds%begp:bounds%endp,1:nlevgrnd), &
            soil_water_retention_curve=soil_water_retention_curve)

    case default
       call endrun(subname // ':: a root moisture stress function must be specified!')     
    end select

  end subroutine calc_root_moist_stress

end module SoilMoistStressMod
