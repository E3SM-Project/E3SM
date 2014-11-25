module UrbanRadiationMod

#include "shr_assert.h"

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_flush 
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_varpar        , only : numrad
  use clm_varcon        , only : isecspday, degpsec, namel
  use clm_varctl        , only : iulog
  use abortutils        , only : endrun  
  use UrbanParamsType   , only : urbanparams_type
  use atm2lndType       , only : atm2lnd_type
  use WaterStateType    , only : waterstate_type
  use TemperatureType   , only : temperature_type
  use SolarAbsorbedType , only : solarabs_type 
  use SurfaceAlbedoType , only : surfalb_type
  use UrbanParamsType   , only : urbanparams_type
  use EnergyFluxType    , only : energyflux_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanRadiation    ! Urban physics - radiative fluxes
  !
  ! PRIVATE MEMBER FUNCTIONS
  private :: net_longwave     ! Net longwave radiation for road and both walls in urban canyon 
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanRadiation (bounds                                   , &
       num_nourbanl, filter_nourbanl                                  , &
       num_urbanl, filter_urbanl                                      , &
       num_urbanc, filter_urbanc                                      , &
       num_urbanp, filter_urbanp                                      , &
       atm2lnd_vars, waterstate_vars, temperature_vars, urbanparams_vars, &
       solarabs_vars, surfalb_vars, energyflux_vars)
    !
    ! !DESCRIPTION: 
    ! Solar fluxes absorbed and reflected by roof and canyon (walls, road).
    ! Also net and upward longwave fluxes.

    ! !USES:
    use clm_varcon          , only : spval, sb, tfrz
    use column_varcon       , only : icol_road_perv, icol_road_imperv
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_time_manager    , only : get_curr_date, get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds    
    integer                , intent(in)    :: num_nourbanl       ! number of non-urban landunits in clump
    integer                , intent(in)    :: filter_nourbanl(:) ! non-urban landunit filter
    integer                , intent(in)    :: num_urbanl         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)   ! urban landunit filter
    integer                , intent(in)    :: num_urbanc         ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:)   ! urban column filter
    integer                , intent(in)    :: num_urbanp         ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:)   ! urban pft filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(waterstate_type)  , intent(in)    :: waterstate_vars
    type(temperature_type) , intent(in)    :: temperature_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(surfalb_type)     , intent(in)    :: surfalb_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,fl,p,c,l,g              ! indices
    integer  :: local_secp1                ! seconds into current date in local time
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: year,month,day             ! temporaries (not used)
    integer  :: secs                       ! seconds into current date

    real(r8), parameter :: mpe    = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: snoem  = 0.97_r8   ! snow emissivity (should use value from Biogeophysics1)

    real(r8) :: lwnet_roof(bounds%begl:bounds%endl)     ! net (outgoing-incoming) longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwnet_improad(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwnet_perroad(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwnet_sunwall(bounds%begl:bounds%endl)  ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwnet_shadewall(bounds%begl:bounds%endl)! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwnet_canyon(bounds%begl:bounds%endl)   ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: lwup_roof(bounds%begl:bounds%endl)      ! upward longwave radiation (per unit ground area), roof (W/m**2)
    real(r8) :: lwup_improad(bounds%begl:bounds%endl)   ! upward longwave radiation (per unit ground area), impervious road (W/m**2)
    real(r8) :: lwup_perroad(bounds%begl:bounds%endl)   ! upward longwave radiation (per unit ground area), pervious road (W/m**2)
    real(r8) :: lwup_sunwall(bounds%begl:bounds%endl)   ! upward longwave radiation, (per unit wall area), sunlit wall (W/m**2)
    real(r8) :: lwup_shadewall(bounds%begl:bounds%endl) ! upward longwave radiation, (per unit wall area), shaded wall (W/m**2)
    real(r8) :: lwup_canyon(bounds%begl:bounds%endl)    ! upward longwave radiation for canyon, per unit ground area (W/m**2)
    real(r8) :: t_roof(bounds%begl:bounds%endl)         ! roof temperature (K)
    real(r8) :: t_improad(bounds%begl:bounds%endl)      ! imppervious road temperature (K)
    real(r8) :: t_perroad(bounds%begl:bounds%endl)      ! pervious road temperature (K)
    real(r8) :: t_sunwall(bounds%begl:bounds%endl)      ! sunlit wall temperature (K)
    real(r8) :: t_shadewall(bounds%begl:bounds%endl)    ! shaded wall temperature (K)
    real(r8) :: lwdown(bounds%begl:bounds%endl)         ! atmospheric downward longwave radiation (W/m**2)
    real(r8) :: em_roof_s(bounds%begl:bounds%endl)      ! roof emissivity with snow effects
    real(r8) :: em_improad_s(bounds%begl:bounds%endl)   ! impervious road emissivity with snow effects
    real(r8) :: em_perroad_s(bounds%begl:bounds%endl)   ! pervious road emissivity with snow effects
    !-----------------------------------------------------------------------

    associate(                                                                 & 
         ctype              =>    col%itype                                  , & ! Input:  [integer (:)    ]  column type                                        
         coli               =>    lun%coli                                   , & ! Input:  [integer (:)    ]  beginning column index for landunit                
         colf               =>    lun%colf                                   , & ! Input:  [integer (:)    ]  ending column index for landunit                   
         canyon_hwr         =>    lun%canyon_hwr                             , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtroad_perv        =>    lun%wtroad_perv                            , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road            

         forc_solad         =>    atm2lnd_vars%forc_solad_grc                , & ! Input:  [real(r8) (:,:) ]  direct beam radiation  (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solai         =>    atm2lnd_vars%forc_solai_grc                , & ! Input:  [real(r8) (:,:) ]  diffuse beam radiation (vis=forc_sols , nir=forc_soll ) (W/m**2)
         forc_solar         =>    atm2lnd_vars%forc_solar_grc                , & ! Input:  [real(r8) (:)   ]  incident solar radiation (W/m**2)                 
         forc_lwrad         =>    atm2lnd_vars%forc_lwrad_not_downscaled_grc , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)   

         frac_sno           =>    waterstate_vars%frac_sno_col               , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       

         t_ref2m            =>    temperature_vars%t_ref2m_patch             , & ! Input:  [real(r8) (:)   ]  2 m height surface air temperature (K)            
         t_grnd             =>    temperature_vars%t_grnd_col                , & ! Input:  [real(r8) (:)   ]  ground temperature (K)                            

         em_roof            =>    urbanparams_vars%em_roof                   , & ! Input:  [real(r8) (:)   ]  roof emissivity                                   
         em_improad         =>    urbanparams_vars%em_improad                , & ! Input:  [real(r8) (:)   ]  impervious road emissivity                        
         em_perroad         =>    urbanparams_vars%em_perroad                , & ! Input:  [real(r8) (:)   ]  pervious road emissivity                          
         em_wall            =>    urbanparams_vars%em_wall                   , & ! Input:  [real(r8) (:)   ]  wall emissivity                                   

         albd               =>    surfalb_vars%albd_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (direct)                         
         albi               =>    surfalb_vars%albi_patch                    , & ! Input:  [real(r8) (:,:) ] pft surface albedo (diffuse)                        
         
         sabs_roof_dir      =>    solarabs_vars%sabs_roof_dir_lun            , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by roof per unit ground area per unit incident flux
         sabs_roof_dif      =>    solarabs_vars%sabs_roof_dif_lun            , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by roof per unit ground area per unit incident flux
         sabs_sunwall_dir   =>    solarabs_vars%sabs_sunwall_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_sunwall_dif   =>    solarabs_vars%sabs_sunwall_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by sunwall per unit wall area per unit incident flux
         sabs_shadewall_dir =>    solarabs_vars%sabs_shadewall_dir_lun       , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_shadewall_dif =>    solarabs_vars%sabs_shadewall_dif_lun       , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by shadewall per unit wall area per unit incident flux
         sabs_improad_dir   =>    solarabs_vars%sabs_improad_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_improad_dif   =>    solarabs_vars%sabs_improad_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by impervious road per unit ground area per unit incident flux
         sabs_perroad_dir   =>    solarabs_vars%sabs_perroad_dir_lun         , & ! Output: [real(r8) (:,:) ]  direct  solar absorbed  by pervious road per unit ground area per unit incident flux
         sabs_perroad_dif   =>    solarabs_vars%sabs_perroad_dif_lun         , & ! Output: [real(r8) (:,:) ]  diffuse solar absorbed  by pervious road per unit ground area per unit incident flux
         sabg               =>    solarabs_vars%sabg_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by ground (W/m**2)       
         sabv               =>    solarabs_vars%sabv_patch                   , & ! Output: [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)   
         fsa                =>    solarabs_vars%fsa_patch                    , & ! Output: [real(r8) (:)   ]  solar radiation absorbed (total) (W/m**2)         
         fsa_u              =>    solarabs_vars%fsa_u_patch                  , & ! Output: [real(r8) (:)   ]  urban solar radiation absorbed (total) (W/m**2)   

         eflx_lwrad_out     =>    energyflux_vars%eflx_lwrad_out_patch       , & ! Output: [real(r8) (:)   ]  emitted infrared (longwave) radiation (W/m**2)    
         eflx_lwrad_net     =>    energyflux_vars%eflx_lwrad_net_patch       , & ! Output: [real(r8) (:)   ]  net infrared (longwave) rad (W/m**2) [+ = to atm] 
         eflx_lwrad_net_u   =>    energyflux_vars%eflx_lwrad_net_u_patch     , & ! Output: [real(r8) (:)   ]  urban net infrared (longwave) rad (W/m**2) [+ = to atm]

         begl               =>    bounds%begl                                , &
         endl               =>    bounds%endl                                  &
         )

      ! Define fields that appear on the restart file for non-urban landunits 
      
      do fl = 1,num_nourbanl
         l = filter_nourbanl(fl)
         sabs_roof_dir(l,:)      = spval
         sabs_roof_dif(l,:)      = spval
         sabs_sunwall_dir(l,:)   = spval
         sabs_sunwall_dif(l,:)   = spval
         sabs_shadewall_dir(l,:) = spval
         sabs_shadewall_dif(l,:) = spval
         sabs_improad_dir(l,:)   = spval
         sabs_improad_dif(l,:)   = spval
         sabs_perroad_dir(l,:)   = spval
         sabs_perroad_dif(l,:)   = spval
      end do

      ! Set input forcing fields
      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l) 

         ! Need to set the following temperatures to some defined value even if it
         ! does not appear in the urban landunit for the net_longwave computation

         t_roof(l)      = 19._r8 + tfrz
         t_sunwall(l)   = 19._r8 + tfrz
         t_shadewall(l) = 19._r8 + tfrz
         t_improad(l)   = 19._r8 + tfrz
         t_perroad(l)   = 19._r8 + tfrz

         ! Initial assignment of emissivity
         em_roof_s(l)    = em_roof(l)
         em_improad_s(l) = em_improad(l)
         em_perroad_s(l) = em_perroad(l)

         ! Set urban temperatures and emissivity including snow effects.
         do c = coli(l),colf(l)
            if (ctype(c) == icol_roof       )  then
               t_roof(l)      = t_grnd(c)
               em_roof_s(l) = em_roof(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_imperv) then 
               t_improad(l)   = t_grnd(c)
               em_improad_s(l) = em_improad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_road_perv  ) then
               t_perroad(l)   = t_grnd(c)
               em_perroad_s(l) = em_perroad(l)*(1._r8-frac_sno(c)) + snoem*frac_sno(c)
            else if (ctype(c) == icol_sunwall    ) then
               t_sunwall(l)   = t_grnd(c)
            else if (ctype(c) == icol_shadewall  ) then
               t_shadewall(l) = t_grnd(c)
            end if
         end do
         lwdown(l) = forc_lwrad(g)
      end do

      ! Net longwave radiation for road and both walls in urban canyon allowing for multiple re-emission

      if (num_urbanl > 0) then
         call net_longwave (bounds,       &
              num_urbanl, filter_urbanl,  &
              canyon_hwr(begl:endl),      &
              wtroad_perv(begl:endl),     &
              lwdown(begl:endl),          &
              em_roof_s(begl:endl),       &
              em_improad_s(begl:endl),    &
              em_perroad_s(begl:endl),    &
              em_wall(begl:endl),         &
              t_roof(begl:endl),          &
              t_improad(begl:endl),       &
              t_perroad(begl:endl),       &
              t_sunwall(begl:endl),       &
              t_shadewall(begl:endl),     &
              lwnet_roof(begl:endl),      &
              lwnet_improad(begl:endl),   &
              lwnet_perroad(begl:endl),   &
              lwnet_sunwall(begl:endl),   &
              lwnet_shadewall(begl:endl), &
              lwnet_canyon(begl:endl),    &
              lwup_roof(begl:endl),       &
              lwup_improad(begl:endl),    &
              lwup_perroad(begl:endl),    &
              lwup_sunwall(begl:endl),    &
              lwup_shadewall(begl:endl),  &
              lwup_canyon(begl:endl),     &
              urbanparams_vars)
      end if

      dtime = get_step_size()
      call get_curr_date (year, month, day, secs)

      ! Determine variables needed for history output and communication with atm
      ! Loop over urban patches in clump

      do fp = 1,num_urbanp
         p = filter_urbanp(fp)
         c = pft%column(p)
         l = pft%landunit(p)
         g = pft%gridcell(p)

         ! Solar absorbed and longwave out and net
         ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)
         ! Each urban pft has its own column - this is used in the logic below

         if (ctype(c) == icol_roof) then   
            eflx_lwrad_out(p) = lwup_roof(l)
            eflx_lwrad_net(p) = lwnet_roof(l)
            eflx_lwrad_net_u(p) = lwnet_roof(l)
            sabg(p) = sabs_roof_dir(l,1)*forc_solad(g,1) + &
                 sabs_roof_dif(l,1)*forc_solai(g,1) + &
                 sabs_roof_dir(l,2)*forc_solad(g,2) + &
                 sabs_roof_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_sunwall) then   
            eflx_lwrad_out(p)   = lwup_sunwall(l)
            eflx_lwrad_net(p)   = lwnet_sunwall(l)
            eflx_lwrad_net_u(p) = lwnet_sunwall(l)
            sabg(p) = sabs_sunwall_dir(l,1)*forc_solad(g,1) + &
                 sabs_sunwall_dif(l,1)*forc_solai(g,1) + &
                 sabs_sunwall_dir(l,2)*forc_solad(g,2) + &
                 sabs_sunwall_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_shadewall) then   
            eflx_lwrad_out(p)   = lwup_shadewall(l)
            eflx_lwrad_net(p)   = lwnet_shadewall(l)
            eflx_lwrad_net_u(p) = lwnet_shadewall(l)
            sabg(p) = sabs_shadewall_dir(l,1)*forc_solad(g,1) + &
                 sabs_shadewall_dif(l,1)*forc_solai(g,1) + &
                 sabs_shadewall_dir(l,2)*forc_solad(g,2) + &
                 sabs_shadewall_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_road_perv) then       
            eflx_lwrad_out(p)   = lwup_perroad(l)
            eflx_lwrad_net(p)   = lwnet_perroad(l)
            eflx_lwrad_net_u(p) = lwnet_perroad(l)
            sabg(p) = sabs_perroad_dir(l,1)*forc_solad(g,1) + &
                 sabs_perroad_dif(l,1)*forc_solai(g,1) + &
                 sabs_perroad_dir(l,2)*forc_solad(g,2) + &
                 sabs_perroad_dif(l,2)*forc_solai(g,2) 

         else if (ctype(c) == icol_road_imperv) then       
            eflx_lwrad_out(p)   = lwup_improad(l)
            eflx_lwrad_net(p)   = lwnet_improad(l)
            eflx_lwrad_net_u(p) = lwnet_improad(l)
            sabg(p) = sabs_improad_dir(l,1)*forc_solad(g,1) + &
                 sabs_improad_dif(l,1)*forc_solai(g,1) + &
                 sabs_improad_dir(l,2)*forc_solad(g,2) + &
                 sabs_improad_dif(l,2)*forc_solai(g,2) 
         end if

         sabv(p)   = 0._r8
         fsa(p)    = sabv(p) + sabg(p)
         fsa_u(p)  = fsa(p)

      end do ! end loop over urban patches

    end associate

  end subroutine UrbanRadiation

  !-----------------------------------------------------------------------
  subroutine net_longwave (bounds                                                             , &
       num_urbanl, filter_urbanl, canyon_hwr, wtroad_perv                                     , &
       lwdown, em_roof, em_improad, em_perroad, em_wall                                       , &
       t_roof,  t_improad, t_perroad, t_sunwall, t_shadewall                                  , &
       lwnet_roof, lwnet_improad, lwnet_perroad, lwnet_sunwall, lwnet_shadewall, lwnet_canyon , &
       lwup_roof, lwup_improad, lwup_perroad, lwup_sunwall, lwup_shadewall, lwup_canyon, &
       urbanparams_vars)
    !
    ! !DESCRIPTION: 
    ! Net longwave radiation for road and both walls in urban canyon allowing for 
    ! multiple reflection. Also net longwave radiation for urban roof. 
    !
    ! !USES:
    use clm_varcon , only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                  
    integer , intent(in)  :: num_urbanl                      ! number of urban landunits
    integer , intent(in)  :: filter_urbanl(:)                ! urban landunit filter
    real(r8), intent(in)  :: canyon_hwr( bounds%begl: )      ! ratio of building height to street width [landunit]
    real(r8), intent(in)  :: wtroad_perv( bounds%begl: )     ! weight of pervious road wrt total road [landunit]

    real(r8), intent(in)  :: lwdown( bounds%begl: )          ! atmospheric longwave radiation (W/m**2) [landunit]
    real(r8), intent(in)  :: em_roof( bounds%begl: )         ! roof emissivity [landunit]
    real(r8), intent(in)  :: em_improad( bounds%begl: )      ! impervious road emissivity [landunit]
    real(r8), intent(in)  :: em_perroad( bounds%begl: )      ! pervious road emissivity [landunit]
    real(r8), intent(in)  :: em_wall( bounds%begl: )         ! wall emissivity [landunit]

    real(r8), intent(in)  :: t_roof( bounds%begl: )          ! roof temperature (K) [landunit]
    real(r8), intent(in)  :: t_improad( bounds%begl: )       ! impervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_perroad( bounds%begl: )       ! ervious road temperature (K) [landunit]
    real(r8), intent(in)  :: t_sunwall( bounds%begl: )       ! sunlit wall temperature (K) [landunit]
    real(r8), intent(in)  :: t_shadewall( bounds%begl: )     ! shaded wall temperature (K) [landunit]

    real(r8), intent(out) :: lwnet_roof( bounds%begl: )      ! net (outgoing-incoming) longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_improad( bounds%begl: )   ! net (outgoing-incoming) longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_perroad( bounds%begl: )   ! net (outgoing-incoming) longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_sunwall( bounds%begl: )   ! net (outgoing-incoming) longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_shadewall( bounds%begl: ) ! net (outgoing-incoming) longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwnet_canyon( bounds%begl: )    ! net (outgoing-incoming) longwave radiation for canyon, per unit ground area (W/m**2) [landunit]

    real(r8), intent(out) :: lwup_roof( bounds%begl: )       ! upward longwave radiation, roof (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_improad( bounds%begl: )    ! upward longwave radiation, impervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_perroad( bounds%begl: )    ! upward longwave radiation, pervious road (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_sunwall( bounds%begl: )    ! upward longwave radiation (per unit wall area), sunlit wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_shadewall( bounds%begl: )  ! upward longwave radiation (per unit wall area), shaded wall (W/m**2) [landunit]
    real(r8), intent(out) :: lwup_canyon( bounds%begl: )     ! upward longwave radiation for canyon, per unit ground area (W/m**2) [landunit]
    !
    type(urbanparams_type) , intent(in) :: urbanparams_vars
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lwdown_road(bounds%begl:bounds%endl)         ! atmospheric longwave radiation for total road (W/m**2)
    real(r8) :: lwdown_sunwall(bounds%begl:bounds%endl)      ! atmospheric longwave radiation (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: lwdown_shadewall(bounds%begl:bounds%endl)    ! atmospheric longwave radiation (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: lwtot(bounds%begl:bounds%endl)               ! incoming longwave radiation (W/m**2)

    real(r8) :: improad_a(bounds%begl:bounds%endl)           ! absorbed longwave for improad (W/m**2)
    real(r8) :: improad_r(bounds%begl:bounds%endl)           ! reflected longwave for improad (W/m**2)
    real(r8) :: improad_r_sky(bounds%begl:bounds%endl)       ! improad_r to sky (W/m**2)
    real(r8) :: improad_r_sunwall(bounds%begl:bounds%endl)   ! improad_r to sunlit wall (W/m**2)
    real(r8) :: improad_r_shadewall(bounds%begl:bounds%endl) ! improad_r to shaded wall (W/m**2)
    real(r8) :: improad_e(bounds%begl:bounds%endl)           ! emitted longwave for improad (W/m**2)
    real(r8) :: improad_e_sky(bounds%begl:bounds%endl)       ! improad_e to sky (W/m**2)
    real(r8) :: improad_e_sunwall(bounds%begl:bounds%endl)   ! improad_e to sunlit wall (W/m**2)
    real(r8) :: improad_e_shadewall(bounds%begl:bounds%endl) ! improad_e to shaded wall (W/m**2)

    real(r8) :: perroad_a(bounds%begl:bounds%endl)           ! absorbed longwave for perroad (W/m**2)
    real(r8) :: perroad_r(bounds%begl:bounds%endl)           ! reflected longwave for perroad (W/m**2)
    real(r8) :: perroad_r_sky(bounds%begl:bounds%endl)       ! perroad_r to sky (W/m**2)
    real(r8) :: perroad_r_sunwall(bounds%begl:bounds%endl)   ! perroad_r to sunlit wall (W/m**2)
    real(r8) :: perroad_r_shadewall(bounds%begl:bounds%endl) ! perroad_r to shaded wall (W/m**2)
    real(r8) :: perroad_e(bounds%begl:bounds%endl)           ! emitted longwave for perroad (W/m**2)
    real(r8) :: perroad_e_sky(bounds%begl:bounds%endl)       ! perroad_e to sky (W/m**2)
    real(r8) :: perroad_e_sunwall(bounds%begl:bounds%endl)   ! perroad_e to sunlit wall (W/m**2)
    real(r8) :: perroad_e_shadewall(bounds%begl:bounds%endl) ! perroad_e to shaded wall (W/m**2)

    real(r8) :: road_a(bounds%begl:bounds%endl)              ! absorbed longwave for total road (W/m**2)
    real(r8) :: road_r(bounds%begl:bounds%endl)              ! reflected longwave for total road (W/m**2)
    real(r8) :: road_r_sky(bounds%begl:bounds%endl)          ! total road_r to sky (W/m**2)
    real(r8) :: road_r_sunwall(bounds%begl:bounds%endl)      ! total road_r to sunlit wall (W/m**2)
    real(r8) :: road_r_shadewall(bounds%begl:bounds%endl)    ! total road_r to shaded wall (W/m**2)
    real(r8) :: road_e(bounds%begl:bounds%endl)              ! emitted longwave for total road (W/m**2)
    real(r8) :: road_e_sky(bounds%begl:bounds%endl)          ! total road_e to sky (W/m**2)
    real(r8) :: road_e_sunwall(bounds%begl:bounds%endl)      ! total road_e to sunlit wall (W/m**2)
    real(r8) :: road_e_shadewall(bounds%begl:bounds%endl)    ! total road_e to shaded wall (W/m**2)

    real(r8) :: sunwall_a(bounds%begl:bounds%endl)           ! absorbed longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r(bounds%begl:bounds%endl)           ! reflected longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_r_sky(bounds%begl:bounds%endl)       ! sunwall_r to sky (W/m**2)
    real(r8) :: sunwall_r_road(bounds%begl:bounds%endl)      ! sunwall_r to road (W/m**2)
    real(r8) :: sunwall_r_shadewall(bounds%begl:bounds%endl) ! sunwall_r to opposing (shaded) wall (W/m**2)
    real(r8) :: sunwall_e(bounds%begl:bounds%endl)           ! emitted longwave (per unit wall area) for sunlit wall (W/m**2)
    real(r8) :: sunwall_e_sky(bounds%begl:bounds%endl)       ! sunwall_e to sky (W/m**2)
    real(r8) :: sunwall_e_road(bounds%begl:bounds%endl)      ! sunwall_e to road (W/m**2)
    real(r8) :: sunwall_e_shadewall(bounds%begl:bounds%endl) ! sunwall_e to opposing (shaded) wall (W/m**2)

    real(r8) :: shadewall_a(bounds%begl:bounds%endl)         ! absorbed longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r(bounds%begl:bounds%endl)         ! reflected longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_r_sky(bounds%begl:bounds%endl)     ! shadewall_r to sky (W/m**2)
    real(r8) :: shadewall_r_road(bounds%begl:bounds%endl)    ! shadewall_r to road (W/m**2)
    real(r8) :: shadewall_r_sunwall(bounds%begl:bounds%endl) ! shadewall_r to opposing (sunlit) wall (W/m**2)
    real(r8) :: shadewall_e(bounds%begl:bounds%endl)         ! emitted longwave (per unit wall area) for shaded wall (W/m**2)
    real(r8) :: shadewall_e_sky(bounds%begl:bounds%endl)     ! shadewall_e to sky (W/m**2)
    real(r8) :: shadewall_e_road(bounds%begl:bounds%endl)    ! shadewall_e to road (W/m**2)
    real(r8) :: shadewall_e_sunwall(bounds%begl:bounds%endl) ! shadewall_e to opposing (sunlit) wall (W/m**2)
    integer  :: l,fl,iter                    ! indices
    integer, parameter  :: n = 50            ! number of interations
    real(r8) :: crit                         ! convergence criterion (W/m**2)
    real(r8) :: err                          ! energy conservation error (W/m**2)
    real(r8) :: wtroad_imperv(bounds%begl:bounds%endl)       ! weight of impervious road wrt total road
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes

    SHR_ASSERT_ALL((ubound(canyon_hwr)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wtroad_perv)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwdown)          == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_roof)         == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_improad)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_perroad)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(em_wall)         == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_roof)          == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_improad)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_perroad)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_sunwall)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_shadewall)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_roof)      == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_improad)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_perroad)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_sunwall)   == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_shadewall) == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwnet_canyon)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_roof)       == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_improad)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_perroad)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_sunwall)    == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_shadewall)  == (/bounds%endl/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(lwup_canyon)     == (/bounds%endl/)), errMsg(__FILE__, __LINE__))

    associate(                             & 
         vf_sr => urbanparams_vars%vf_sr , & ! Input:  [real(r8) (:)]  view factor of sky for road                       
         vf_wr => urbanparams_vars%vf_wr , & ! Input:  [real(r8) (:)]  view factor of one wall for road                  
         vf_sw => urbanparams_vars%vf_sw , & ! Input:  [real(r8) (:)]  view factor of sky for one wall                   
         vf_rw => urbanparams_vars%vf_rw , & ! Input:  [real(r8) (:)]  view factor of road for one wall                  
         vf_ww => urbanparams_vars%vf_ww   & ! Input:  [real(r8) (:)]  view factor of opposing wall for one wall         
         )

      ! Calculate impervious road

      do fl = 1,num_urbanl 
         l = filter_urbanl(fl)
         wtroad_imperv(l) = 1._r8 - wtroad_perv(l)
      end do

      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         ! atmospheric longwave radiation incident on walls and road in urban canyon.
         ! check for conservation (need to convert wall fluxes to ground area).
         ! lwdown (from atmosphere) = lwdown_road + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

         lwdown_road(l)      = lwdown(l) * vf_sr(l) 
         lwdown_sunwall(l)   = lwdown(l) * vf_sw(l)
         lwdown_shadewall(l) = lwdown(l) * vf_sw(l) 

         err = lwdown(l) - (lwdown_road(l) + (lwdown_shadewall(l) + lwdown_sunwall(l))*canyon_hwr(l))
         if (abs(err) > 0.10_r8 ) then
            write(iulog,*) 'urban incident atmospheric longwave radiation balance error',err
            write(iulog,*) 'l          = ',l
            write(iulog,*) 'lwdown     = ',lwdown(l)
            write(iulog,*) 'vf_sr      = ',vf_sr(l)
            write(iulog,*) 'vf_sw      = ',vf_sw(l)
            write(iulog,*) 'canyon_hwr = ',canyon_hwr(l)
            write(iulog,*) 'clm model is stopping'
            call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         endif
      end do

      do fl = 1,num_urbanl
         l = filter_urbanl(fl)

         ! initial absorption, reflection, and emission for road and both walls. 
         ! distribute reflected and emitted radiation to sky, road, and walls according 
         ! to appropriate view factor. radiation reflected to road and walls will
         ! undergo multiple reflections within the canyon.

         road_a(l)              = 0.0_r8
         road_r(l)              = 0.0_r8
         road_e(l)              = 0.0_r8
         improad_a(l)           =     em_improad(l)  * lwdown_road(l) 
         improad_r(l)           = (1._r8-em_improad(l)) * lwdown_road(l) 
         improad_r_sky(l)       = improad_r(l) * vf_sr(l)
         improad_r_sunwall(l)   = improad_r(l) * vf_wr(l)
         improad_r_shadewall(l) = improad_r(l) * vf_wr(l)
         improad_e(l)           = em_improad(l) * sb * (t_improad(l)**4) 
         improad_e_sky(l)       = improad_e(l) * vf_sr(l)
         improad_e_sunwall(l)   = improad_e(l) * vf_wr(l)
         improad_e_shadewall(l) = improad_e(l) * vf_wr(l)
         road_a(l)              = road_a(l) + improad_a(l)*wtroad_imperv(l)
         road_r(l)              = road_r(l) + improad_r(l)*wtroad_imperv(l)
         road_e(l)              = road_e(l) + improad_e(l)*wtroad_imperv(l)

         perroad_a(l)           =     em_perroad(l)  * lwdown_road(l)
         perroad_r(l)           = (1._r8-em_perroad(l)) * lwdown_road(l)
         perroad_r_sky(l)       = perroad_r(l) * vf_sr(l)
         perroad_r_sunwall(l)   = perroad_r(l) * vf_wr(l)
         perroad_r_shadewall(l) = perroad_r(l) * vf_wr(l)
         perroad_e(l)           = em_perroad(l) * sb * (t_perroad(l)**4) 
         perroad_e_sky(l)       = perroad_e(l) * vf_sr(l)
         perroad_e_sunwall(l)   = perroad_e(l) * vf_wr(l)
         perroad_e_shadewall(l) = perroad_e(l) * vf_wr(l)
         road_a(l)              = road_a(l) + perroad_a(l)*wtroad_perv(l)
         road_r(l)              = road_r(l) + perroad_r(l)*wtroad_perv(l)
         road_e(l)              = road_e(l) + perroad_e(l)*wtroad_perv(l)

         road_r_sky(l)          = road_r(l) * vf_sr(l)
         road_r_sunwall(l)      = road_r(l) * vf_wr(l)
         road_r_shadewall(l)    = road_r(l) * vf_wr(l)
         road_e_sky(l)          = road_e(l) * vf_sr(l)
         road_e_sunwall(l)      = road_e(l) * vf_wr(l)
         road_e_shadewall(l)    = road_e(l) * vf_wr(l)

         sunwall_a(l)           = em_wall(l) * lwdown_sunwall(l)
         sunwall_r(l)           = (1._r8-em_wall(l)) * lwdown_sunwall(l)
         sunwall_r_sky(l)       = sunwall_r(l) * vf_sw(l)
         sunwall_r_road(l)      = sunwall_r(l) * vf_rw(l)
         sunwall_r_shadewall(l) = sunwall_r(l) * vf_ww(l)
         sunwall_e(l)           = em_wall(l) * sb * (t_sunwall(l)**4) 
         sunwall_e_sky(l)       = sunwall_e(l) * vf_sw(l)
         sunwall_e_road(l)      = sunwall_e(l) * vf_rw(l)
         sunwall_e_shadewall(l) = sunwall_e(l) * vf_ww(l)

         shadewall_a(l)         = em_wall(l) * lwdown_shadewall(l)
         shadewall_r(l)         = (1._r8-em_wall(l)) * lwdown_shadewall(l)
         shadewall_r_sky(l)     = shadewall_r(l) * vf_sw(l)
         shadewall_r_road(l)    = shadewall_r(l) * vf_rw(l)
         shadewall_r_sunwall(l) = shadewall_r(l) * vf_ww(l)
         shadewall_e(l)         = em_wall(l) * sb * (t_shadewall(l)**4) 
         shadewall_e_sky(l)     = shadewall_e(l) * vf_sw(l)
         shadewall_e_road(l)    = shadewall_e(l) * vf_rw(l)
         shadewall_e_sunwall(l) = shadewall_e(l) * vf_ww(l)

         ! initialize sum of net and upward longwave radiation for road and both walls

         lwnet_improad(l)   = improad_e(l)   - improad_a(l)
         lwnet_perroad(l)   = perroad_e(l)   - perroad_a(l)
         lwnet_sunwall(l)   = sunwall_e(l)   - sunwall_a(l)
         lwnet_shadewall(l) = shadewall_e(l) - shadewall_a(l)

         lwup_improad(l)   = improad_r_sky(l)   + improad_e_sky(l)
         lwup_perroad(l)   = perroad_r_sky(l)   + perroad_e_sky(l)
         lwup_sunwall(l)   = sunwall_r_sky(l)   + sunwall_e_sky(l)
         lwup_shadewall(l) = shadewall_r_sky(l) + shadewall_e_sky(l)

      end do

      ! now account for absorption and reflection within canyon of fluxes from road and walls 
      ! allowing for multiple reflections
      !
      ! (1) absorption and reflection. note: emission from road and walls absorbed by walls and roads
      !     only occurs in first iteration. zero out for later iterations.
      !
      !     road: fluxes from walls need to be projected to ground area
      !     wall: fluxes from road need to be projected to wall area
      !
      ! (2) add net longwave for ith reflection to total net longwave
      !
      ! (3) distribute reflected radiation to sky, road, and walls according to view factors
      !
      ! (4) add upward longwave radiation to sky from road and walls for ith reflection to total
      !
      ! (5) stop iteration when absorption for ith reflection is less than some nominal amount. 
      !     small convergence criteria is required to ensure radiation is conserved

      do fl = 1,num_urbanl
         l = filter_urbanl(fl)

         do iter = 1, n
            ! step (1)

            lwtot(l) =  (sunwall_r_road(l) + sunwall_e_road(l)  &
                 + shadewall_r_road(l) + shadewall_e_road(l))*canyon_hwr(l)
            road_a(l)    = 0.0_r8
            road_r(l)    = 0.0_r8
            improad_r(l) = (1._r8-em_improad(l)) * lwtot(l)
            improad_a(l) =     em_improad(l)  * lwtot(l)
            road_a(l)    = road_a(l) + improad_a(l)*wtroad_imperv(l)
            road_r(l)    = road_r(l) + improad_r(l)*wtroad_imperv(l)
            perroad_r(l) = (1._r8-em_perroad(l)) * lwtot(l) 
            perroad_a(l) =     em_perroad(l)  * lwtot(l) 
            road_a(l)    = road_a(l) + perroad_a(l)*wtroad_perv(l)
            road_r(l)    = road_r(l) + perroad_r(l)*wtroad_perv(l)

            lwtot(l) = (road_r_sunwall(l) + road_e_sunwall(l))/canyon_hwr(l) &
                 + (shadewall_r_sunwall(l) + shadewall_e_sunwall(l))
            sunwall_a(l) =     em_wall(l)  * lwtot(l)
            sunwall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

            lwtot(l) = (road_r_shadewall(l) + road_e_shadewall(l))/canyon_hwr(l) &
                 + (sunwall_r_shadewall(l) + sunwall_e_shadewall(l))
            shadewall_a(l) =     em_wall(l)  * lwtot(l)
            shadewall_r(l) = (1._r8-em_wall(l)) * lwtot(l)

            sunwall_e_road(l)      = 0._r8
            shadewall_e_road(l)    = 0._r8
            road_e_sunwall(l)      = 0._r8
            shadewall_e_sunwall(l) = 0._r8
            road_e_shadewall(l)    = 0._r8
            sunwall_e_shadewall(l) = 0._r8

            ! step (2)

            lwnet_improad(l)   = lwnet_improad(l)   - improad_a(l)
            lwnet_perroad(l)   = lwnet_perroad(l)   - perroad_a(l)
            lwnet_sunwall(l)   = lwnet_sunwall(l)   - sunwall_a(l)
            lwnet_shadewall(l) = lwnet_shadewall(l) - shadewall_a(l)

            ! step (3)

            improad_r_sky(l)       = improad_r(l) * vf_sr(l)
            improad_r_sunwall(l)   = improad_r(l) * vf_wr(l)
            improad_r_shadewall(l) = improad_r(l) * vf_wr(l)

            perroad_r_sky(l)       = perroad_r(l) * vf_sr(l)
            perroad_r_sunwall(l)   = perroad_r(l) * vf_wr(l)
            perroad_r_shadewall(l) = perroad_r(l) * vf_wr(l)

            road_r_sky(l)          = road_r(l) * vf_sr(l)
            road_r_sunwall(l)      = road_r(l) * vf_wr(l)
            road_r_shadewall(l)    = road_r(l) * vf_wr(l)

            sunwall_r_sky(l)       = sunwall_r(l) * vf_sw(l)
            sunwall_r_road(l)      = sunwall_r(l) * vf_rw(l)
            sunwall_r_shadewall(l) = sunwall_r(l) * vf_ww(l)

            shadewall_r_sky(l)     = shadewall_r(l) * vf_sw(l)
            shadewall_r_road(l)    = shadewall_r(l) * vf_rw(l)
            shadewall_r_sunwall(l) = shadewall_r(l) * vf_ww(l)

            ! step (4)

            lwup_improad(l)   = lwup_improad(l)   + improad_r_sky(l)
            lwup_perroad(l)   = lwup_perroad(l)   + perroad_r_sky(l)
            lwup_sunwall(l)   = lwup_sunwall(l)   + sunwall_r_sky(l)
            lwup_shadewall(l) = lwup_shadewall(l) + shadewall_r_sky(l)

            ! step (5)

            crit = max(road_a(l), sunwall_a(l), shadewall_a(l))
            if (crit < .001_r8) exit
         end do
         if (iter >= n) then
            write (iulog,*) 'urban net longwave radiation error: no convergence'
            write (iulog,*) 'clm model is stopping'
            call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         endif

         ! total net longwave radiation for canyon. project wall fluxes to horizontal surface

         lwnet_canyon(l) = 0.0_r8
         lwnet_canyon(l) = lwnet_canyon(l) + lwnet_improad(l)*wtroad_imperv(l)
         lwnet_canyon(l) = lwnet_canyon(l) + lwnet_perroad(l)*wtroad_perv(l)
         lwnet_canyon(l) = lwnet_canyon(l) + (lwnet_sunwall(l) + lwnet_shadewall(l))*canyon_hwr(l)

         ! total emitted longwave for canyon. project wall fluxes to horizontal

         lwup_canyon(l) = 0.0_r8
         lwup_canyon(l) = lwup_canyon(l) + lwup_improad(l)*wtroad_imperv(l)
         lwup_canyon(l) = lwup_canyon(l) + lwup_perroad(l)*wtroad_perv(l)
         lwup_canyon(l) = lwup_canyon(l) + (lwup_sunwall(l) + lwup_shadewall(l))*canyon_hwr(l)

         ! conservation check. note: previous conservation check confirms partioning of incident
         ! atmospheric longwave radiation to road and walls is conserved as
         ! lwdown (from atmosphere) = lwdown_improad + lwdown_perroad + (lwdown_sunwall + lwdown_shadewall)*canyon_hwr

         err = lwnet_canyon(l) - (lwup_canyon(l) - lwdown(l))
         if (abs(err) > .10_r8 ) then
            write (iulog,*) 'urban net longwave radiation balance error',err
            write (iulog,*) 'clm model is stopping'
            call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if

      end do

      ! Net longwave radiation for roof

      do fl = 1,num_urbanl
         l = filter_urbanl(fl)
         lwup_roof(l) = em_roof(l)*sb*(t_roof(l)**4) + (1._r8-em_roof(l))*lwdown(l)
         lwnet_roof(l) = lwup_roof(l) - lwdown(l)
      end do

    end associate

  end subroutine net_longwave

end module UrbanRadiationMod
