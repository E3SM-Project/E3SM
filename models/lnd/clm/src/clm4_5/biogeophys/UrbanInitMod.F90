module UrbanInitMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Initialize urban data
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use shr_sys_mod , only : shr_sys_flush 
  use abortutils  , only : endrun  
  use clm_varctl  , only : iulog, use_vancouver, use_mexicocity
  use UrbanMod,     only : urban_traffic, urban_hac, urban_hac_off
  use decompMod   , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initTimeConstUrban ! Initialize urban time constant variables
  public :: initColdUrban      ! cold start initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initTimeConstUrban(bounds)
    !
    ! !DESCRIPTION: 
    ! Initialize urban time-constant variables
    !
    ! - Calculate view factors for road and one wall
    !                                                        WALL    |
    !                  ROAD                                          |
    !                                                         wall   |
    !          -----\          /-----   -             -  |\----------/
    !              | \  vsr   / |       |         r   |  | \  vww   /   s
    !              |  \      /  |       h         o   w  |  \      /    k
    !        wall  |   \    /   | wall  |         a   |  |   \    /     y
    !              |vwr \  / vwr|       |         d   |  |vrw \  / vsw 
    !              ------\/------       -             -  |-----\/-----
    !                   road                                  wall   |
    !              <----- w ---->                                    |
    !                                                    <---- h --->|
    !
    !    vsr = view factor of sky for road          vrw = view factor of road for wall
    !    vwr = view factor of one wall for road     vww = view factor of opposing wall for wall
    !                                               vsw = view factor of sky for wall
    !    vsr + vwr + vwr = 1                        vrw + vww + vsw = 1
    !
    ! Source: Masson, V. (2000) A physically-based scheme for the urban energy budget in 
    ! atmospheric models. Boundary-Layer Meteorology 94:357-397
    !
    ! - Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
    ! Grimmond and Oke (1999)
    !
    ! !USES:
    use clmtype      , only : col, cps  
    use clmtype      , only : lun, lps, lef, namel
    use clm_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, vkc, spval
    use clm_varcon   , only : icol_road_perv, icol_road_imperv, isturb_MIN
    use UrbanInputMod, only : urbinp, UrbanInput
    use UrbanMod     , only : UrbanParamInit
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: j,l,c,p,g        ! indices
    integer               :: nc,fl,ib         ! indices 
    integer               :: dindx            ! urban density type index
    integer               :: ier              ! error status
    real(r8)              :: sumvf            ! sum of view factors for wall or road
    real(r8), parameter   :: alpha = 4.43_r8  ! coefficient used to calculate z_d_town
    real(r8), parameter   :: beta = 1.0_r8    ! coefficient used to calculate z_d_town
    real(r8), parameter   :: C_d = 1.2_r8     ! drag coefficient as used in Grimmond and Oke (1999)
    real(r8)              :: plan_ai          ! plan area index - ratio building area to plan area (-)
    real(r8)              :: frontal_ai       ! frontal area index of buildings (-)
    real(r8)              :: build_lw_ratio   ! building short/long side ratio (-)
    !-----------------------------------------------------------------------

    ! cps%emg                 Output: [real(r8) (:)   ]  ground emissivity                       
    ! lun%z_0_town            Output: [real(r8) (:)   ]  urban landunit momentum roughness length (m)
    ! lun%z_d_town            Output: [real(r8) (:)   ]  urban landunit displacement height (m)  
    ! lun%ht_roof             Output: [real(r8) (:)   ]  height of urban roof (m)                
    ! lun%canyon_hwr          Output: [real(r8) (:)   ]  ratio of building height to street width (-)
    ! lun%wtroad_perv         Output: [real(r8) (:)   ]  weight of pervious column to total road 
    ! lun%ht_roof             Output: [real(r8) (:)   ]  height of urban roof (m)                
    ! lun%wtlunit_roof        Output: [real(r8) (:)   ]  weight of roof with respect to landunit 
    ! lun%canyon_hwr          Output: [real(r8) (:)   ]  urban canyon height to width ratio      
    ! lps%t_building_max      Output: [real(r8) (:)   ]  maximum internal building temperature (K)
    ! lps%t_building_min      Output: [real(r8) (:)   ]  minimum internal building temperature (K)
    ! lps%tk_wall             Output: [real(r8) (:,:) ]  thermal conductivity of urban wall (W/m/K)
    ! lps%tk_roof             Output: [real(r8) (:,:) ]  thermal conductivity of urban roof (W/m/K)
    ! lps%tk_improad          Output: [real(r8) (:,:) ]  thermal conductivity of urban impervious road (W/m/K)
    ! lps%cv_wall             Output: [real(r8) (:,:) ]  thermal conductivity of urban wall (J/m^3/K)
    ! lps%cv_roof             Output: [real(r8) (:,:) ]  thermal conductivity of urban roof (J/m^3/K)
    ! lps%cv_improad          Output: [real(r8) (:,:) ]  thermal conductivity of urban impervious road (J/m^3/K)
    ! lps%thick_wall          Output: [real(r8) (:)   ]  thickness of urban wall (m)             
    ! lps%thick_roof          Output: [real(r8) (:)   ]  thickness of urban roof (m)             
    ! lps%nlev_improad        Output: [integer (:)    ]  number of impervious road layers (-)     
    ! lef%eflx_traffic_factor Output: [real(r8) (:)   ]  multiplicative factor for sensible heat flux from urban traffic
    ! lps%vf_sr               Output: [real(r8) (:)   ]  view factor of sky for road                       
    ! lps%vf_wr               Output: [real(r8) (:)   ]  view factor of one wall for road                  
    ! lps%vf_sw               Output: [real(r8) (:)   ]  view factor of sky for one wall                   
    ! lps%vf_rw               Output: [real(r8) (:)   ]  view factor of road for one wall                  
    ! lps%vf_ww               Output: [real(r8) (:)   ]  view factor of opposing wall for one wall         
    ! lps%sabs_roof_dir       Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
    ! lps%sabs_roof_dif       Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
    ! lps%sabs_sunwall_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
    ! lps%sabs_sunwall_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
    ! lps%sabs_shadewall_dir  Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
    ! lps%sabs_shadewall_dif  Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
    ! lps%sabs_improad_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
    ! lps%sabs_improad_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
    ! lps%sabs_perroad_dir    Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
    ! lps%sabs_perroad_dif    Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux

    ! Initialization of urban parameter data structure

    call UrbanParamInit(bounds)

    ! Initialize time constant urban variables

    do l = bounds%begl,bounds%endl

       ! "0" refers to urban wall/roof surface and "nlevsoi" refers to urban wall/roof bottom
       if (lun%urbpoi(l)) then

          g = lun%gridcell(l)
          dindx = lun%itype(l) - isturb_MIN + 1

          ! Landunit level initialization for urban wall and roof layers and interfaces
          lun%canyon_hwr(l)         = urbinp%canyon_hwr(g,dindx)
          lun%wtroad_perv(l)        = urbinp%wtroad_perv(g,dindx)
          lun%ht_roof(l)            = urbinp%ht_roof(g,dindx)
          lun%wtlunit_roof(l)       = urbinp%wtlunit_roof(g,dindx)
          lps%tk_wall(l,:)          = urbinp%tk_wall(g,dindx,:)
          lps%tk_roof(l,:)          = urbinp%tk_roof(g,dindx,:)
          lps%tk_improad(l,:)       = urbinp%tk_improad(g,dindx,:)
          lps%cv_wall(l,:)          = urbinp%cv_wall(g,dindx,:)
          lps%cv_roof(l,:)          = urbinp%cv_roof(g,dindx,:)
          lps%cv_improad(l,:)       = urbinp%cv_improad(g,dindx,:)
          lps%thick_wall(l)         = urbinp%thick_wall(g,dindx)
          lps%thick_roof(l)         = urbinp%thick_roof(g,dindx)
          lps%nlev_improad(l)       = urbinp%nlev_improad(g,dindx)
          lps%t_building_min(l)     = urbinp%t_building_min(g,dindx)
          lps%t_building_max(l)     = urbinp%t_building_max(g,dindx)

          do c = lun%coli(l),lun%colf(l)
             if (col%itype(c) == icol_roof       ) cps%emg(c) = urbinp%em_roof(g,dindx)
             if (col%itype(c) == icol_sunwall    ) cps%emg(c) = urbinp%em_wall(g,dindx)
             if (col%itype(c) == icol_shadewall  ) cps%emg(c) = urbinp%em_wall(g,dindx)
             if (col%itype(c) == icol_road_imperv) cps%emg(c) = urbinp%em_improad(g,dindx)
             if (col%itype(c) == icol_road_perv  ) cps%emg(c) = urbinp%em_perroad(g,dindx)
          end do

          ! Inferred from Sailor and Lu 2004
          if (urban_traffic) then
             lef%eflx_traffic_factor(l) = 3.6_r8 * (lun%canyon_hwr(l)-0.5_r8) + 1.0_r8
          else
             lef%eflx_traffic_factor(l) = 0.0_r8
          end if

          if (use_vancouver .or. use_mexicocity) then
             ! Freely evolving
             lps%t_building_max(l) = 380.00_r8
             lps%t_building_min(l) = 200.00_r8
          else
             if (urban_hac == urban_hac_off) then
                ! Overwrite values read in from urbinp by freely evolving values
                lps%t_building_max(l) = 380.00_r8
                lps%t_building_min(l) = 200.00_r8
             end if
          end if

          !----------------------------------------------------------------------------------
          ! View factors for road and one wall in urban canyon (depends only on canyon_hwr)
          !----------------------------------------------------------------------------------

          ! road -- sky view factor -> 1 as building height -> 0 
          ! and -> 0 as building height -> infinity

          lps%vf_sr(l) = sqrt(lun%canyon_hwr(l)**2 + 1._r8) - lun%canyon_hwr(l)
          lps%vf_wr(l) = 0.5_r8 * (1._r8 - lps%vf_sr(l))

          ! one wall -- sky view factor -> 0.5 as building height -> 0 
          ! and -> 0 as building height -> infinity

          lps%vf_sw(l) = 0.5_r8 * (lun%canyon_hwr(l) + 1._r8 - sqrt(lun%canyon_hwr(l)**2+1._r8)) / lun%canyon_hwr(l)
          lps%vf_rw(l) = lps%vf_sw(l)
          lps%vf_ww(l) = 1._r8 - lps%vf_sw(l) - lps%vf_rw(l)

          ! error check -- make sure view factor sums to one for road and wall
          sumvf = lps%vf_sr(l) + 2._r8*lps%vf_wr(l)
          if (abs(sumvf-1._r8) > 1.e-06_r8 ) then
             write (iulog,*) 'urban road view factor error',sumvf
             write (iulog,*) 'clm model is stopping'
             call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
          endif
          sumvf = lps%vf_sw(l) + lps%vf_rw(l) + lps%vf_ww(l)
          if (abs(sumvf-1._r8) > 1.e-06_r8 ) then
             write (iulog,*) 'urban wall view factor error',sumvf
             write (iulog,*) 'clm model is stopping'
             call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
          endif

          !----------------------------------------------------------------------------------
          ! Calculate urban land unit aerodynamic constants using Macdonald (1998) as used in
          ! Grimmond and Oke (1999)
          !----------------------------------------------------------------------------------

          ! Calculate plan area index 
          plan_ai = lun%canyon_hwr(l)/(lun%canyon_hwr(l) + 1._r8)

          ! Building shape shortside/longside ratio (e.g. 1 = square )
          ! This assumes the building occupies the entire canyon length
          build_lw_ratio = plan_ai

          ! Calculate frontal area index
          frontal_ai = (1._r8 - plan_ai) * lun%canyon_hwr(l)

          ! Adjust frontal area index for different building configuration
          frontal_ai = frontal_ai * sqrt(1/build_lw_ratio) * sqrt(plan_ai)

          ! Calculate displacement height
          if (use_vancouver) then
             lun%z_d_town(l) = 3.5_r8
          else if (use_mexicocity) then
             lun%z_d_town(l) = 10.9_r8
          else
             lun%z_d_town(l) = (1._r8 + alpha**(-plan_ai) * (plan_ai - 1._r8)) * lun%ht_roof(l)
          end if

          ! Calculate the roughness length
          if (use_vancouver) then
             lun%z_0_town(l) = 0.35_r8
          else if (use_mexicocity) then
             lun%z_0_town(l) = 2.2_r8
          else
             lun%z_0_town(l) = lun%ht_roof(l) * (1._r8 - lun%z_d_town(l) / lun%ht_roof(l)) * &
                  exp(-1.0_r8 * (0.5_r8 * beta * C_d / vkc**2 * &
                  (1 - lun%z_d_town(l) / lun%ht_roof(l)) * frontal_ai)**(-0.5_r8))
          end if

       else

          lef%eflx_traffic_factor(l) = spval
          lps%t_building_max(l) = spval
          lps%t_building_min(l) = spval

          lps%vf_sr(l) = spval
          lps%vf_wr(l) = spval
          lps%vf_sw(l) = spval
          lps%vf_rw(l) = spval
          lps%vf_ww(l) = spval

       end if
    end do

    ! Deallocate memory for urbinp datatype

    call UrbanInput(bounds%begg, bounds%endg, mode='finalize')

  end subroutine initTimeConstUrban

  !-----------------------------------------------------------------------
  subroutine initColdUrban(bounds)
    !
    ! !DESCRIPTION: 
    ! Initialize urban time-varying variables
    !
    ! !USES:
    use clmtype   , only : pft, pes, pes, pef
    use clmtype   , only : col, cws, ces, cef, cwf
    use clmtype   , only : lun, lps, lef 
    use clm_varcon, only : spval, icol_road_perv
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: l,g,c,p       ! indices
    !-----------------------------------------------------------------------

    ! lps%taf                   Output: [real(r8) (:)]  urban canopy air temperature (K)        
    ! lps%qaf                   Output: [real(r8) (:)]  urban canopy air specific humidity (kg/kg)
    ! lps%t_building            Output: [real(r8) (:)]  internal building temperature (K)       
    ! lef%eflx_traffic          Output: [real(r8) (:)]  traffic sensible heat flux (W/m**2)     
    ! lef%eflx_wasteheat        Output: [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
    ! ces%t_grnd_u              Output: [real(r8) (:)]  Urban ground temperature (Kelvin)       
    ! cef%eflx_building_heat    Output: [real(r8) (:)]  heat flux from urban building interior to walls, roof (W/m**2)
    ! cef%eflx_urban_ac         Output: [real(r8) (:)]  urban air conditioning flux (W/m**2)    
    ! cef%eflx_urban_heat       Output: [real(r8) (:)]  urban heating flux (W/m**2)             
    ! cef%eflx_snomelt_u        Output: [real(r8) (:)]  Urban snow melt heat flux (W/m**2)      
    ! cws%fcov                  Output: [real(r8) (:)]  fractional impermeable area             
    ! cws%fsat                  Output: [real(r8) (:)]  fractional area with water table at surface
    ! cws%qcharge               Output: [real(r8) (:)]  aquifer recharge rate (mm/s)            
    ! cwf%qflx_runoff_u         Output: [real(r8) (:)]  Urban total runoff (qflx_drain+qflx_surf) (mm H2O /s)
    ! pes%t_ref2m_u             Output: [real(r8) (:)]  Urban 2 m height surface air temperature (Kelvin)
    ! pes%t_ref2m_min_u         Output: [real(r8) (:)]  Urban daily minimum of average 2 m height surface air temperature (K)
    ! pes%t_ref2m_max_u         Output: [real(r8) (:)]  Urban daily maximum of average 2 m height surface air temperature (K)
    ! pes%rh_ref2m_u            Output: [real(r8) (:)]  Urban 2 m height surface relative humidity (%)
    ! pef%eflx_wasteheat_pft    Output: [real(r8) (:)]  sensible heat flux from urban heating/cooling sources of waste heat at pft level (W/m**2)
    ! pef%eflx_heat_from_ac_pft Output: [real(r8) (:)]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
    ! pef%eflx_traffic_pft      Output: [real(r8) (:)]  sensible heat flux from traffic (W/m**2)
    ! pef%eflx_anthro           Output: [real(r8) (:)]  total anthropogenic heat flux (W/m**2)  
    ! pef%fsa_u                 Output: [real(r8) (:)]  Urban absorbed solar radiation (W/m**2) 
    ! pef%eflx_lwrad_net_u      Output: [real(r8) (:)]  Urban net longwave radiation (W/m**2)   
    ! pef%eflx_lwrad_out_u      Output: [real(r8) (:)]  Urban emitted longwave radiation (W/m**2)
    ! pef%eflx_lh_tot_u         Output: [real(r8) (:)]  Urban latent heat flux (W/m**2)         
    ! pef%eflx_sh_tot_u         Output: [real(r8) (:)]  Urban sensible heat flux (W/m**2)       
    ! pef%eflx_soil_grnd_u      Output: [real(r8) (:)]  Urban ground heat flux (W/m**2)         

    ! Cold start initialization for landunits: 
    do l = bounds%begl, bounds%endl 
       g = lun%gridcell(l)
       if (lun%urbpoi(l)) then
          if (use_vancouver) then
             lps%taf(l) = 297.56_r8
             lps%qaf(l) = 0.0111_r8
          else if (use_mexicocity) then
             lps%taf(l) = 289.46_r8
             lps%qaf(l) = 0.00248_r8
          else
             lps%taf(l) = 283._r8
             ! Arbitrary set since forc_q is not yet available
             lps%qaf(l) = 1.e-4_r8
          end if
       else
          lps%t_building(l)     = spval
          lef%eflx_traffic(l)   = spval
          lef%eflx_wasteheat(l) = spval
       end if
    end do

    ! Cold start initialization for columns: 
    do c = bounds%begc, bounds%endc 
       l = col%landunit(c)
       if (lun%urbpoi(l)) then
          cef%eflx_building_heat(c) = 0._r8
          cef%eflx_urban_ac(c)      = 0._r8
          cef%eflx_urban_heat(c)    = 0._r8

          ! Set hydrology variables for urban to spvalue 
          ! -- as only valid for pervious road
          if (col%itype(c) /= icol_road_perv  )then
             cws%fcov(c)    = spval
             cws%fsat(c)    = spval
             cws%qcharge(c) = spval
          end if
       else
          cef%eflx_building_heat(c) = 0._r8
          cef%eflx_urban_ac(c)      = 0._r8
          cef%eflx_urban_heat(c)    = 0._r8
          ces%t_grnd_u(c)           = spval
          cef%eflx_snomelt_u(c)     = spval
          cwf%qflx_runoff_u(c)      = spval
       end if
    end do

    ! Cold start initialization for pfts:
    do p = bounds%begp, bounds%endp 
       l = pft%landunit(p)
       if (.not. lun%urbpoi(l)) then
          pes%t_ref2m_u(p)             = spval
          pes%t_ref2m_min_u(p)         = spval
          pes%t_ref2m_max_u(p)         = spval
          pes%rh_ref2m_u(p)            = spval
          pef%eflx_wasteheat_pft(p)    = 0._r8
          pef%eflx_heat_from_ac_pft(p) = 0._r8
          pef%eflx_traffic_pft(p)      = 0._r8
          pef%eflx_anthro(p)           = 0._r8
          pef%eflx_lwrad_net_u(p)      = spval
          pef%eflx_lwrad_out_u(p)      = spval
          pef%eflx_lh_tot_u(p)         = spval
          pef%eflx_sh_tot_u(p)         = spval
          pef%eflx_soil_grnd_u(p)      = spval
          pef%fsa_u(p)                 = spval 
       end if
    end do

    ! Initialize urban albedos and absorbed fluxes to reasonable values
    lps%sabs_roof_dir(bounds%begl:bounds%endl, :)      = 0._r8    
    lps%sabs_roof_dif(bounds%begl:bounds%endl, :)      = 0._r8    
    lps%sabs_sunwall_dir(bounds%begl:bounds%endl, :)   = 0._r8
    lps%sabs_sunwall_dif(bounds%begl:bounds%endl, :)   = 0._r8
    lps%sabs_shadewall_dir(bounds%begl:bounds%endl, :) = 0._r8
    lps%sabs_shadewall_dif(bounds%begl:bounds%endl, :) = 0._r8
    lps%sabs_improad_dir(bounds%begl:bounds%endl, :)   = 0._r8
    lps%sabs_improad_dif(bounds%begl:bounds%endl, :)   = 0._r8
    lps%sabs_perroad_dir(bounds%begl:bounds%endl, :)   = 0._r8
    lps%sabs_perroad_dif(bounds%begl:bounds%endl, :)   = 0._r8

  end subroutine initColdUrban

end module UrbanInitMod
