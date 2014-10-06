module dynConsBiogeophysMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeophysical quantities (water & energy) with dynamic land
  ! cover.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use UrbanParamsType   , only : urbanparams_type
  use EnergyFluxType    , only : energyflux_type
  use LakeStateType     , only : lakestate_type
  use SoilHydrologyType , only : soilhydrology_type  
  use SoilStateType     , only : soilstate_type
  use TemperatureType   , only : temperature_type
  use WaterfluxType     , only : waterflux_type
  use WaterstateType    , only : waterstate_type
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft                
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dyn_hwcontent_init            ! compute grid-level heat and water content, before land cover change
  public :: dyn_hwcontent_final           ! compute grid-level heat and water content, after land cover change; also compute dynbal fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS
  private :: dyn_hwcontent                ! do the actual computation of grid-level heat and water content
  !---------------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_init(bounds, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Initialize variables used for dyn_hwcontent, and compute grid cell-level heat
    ! and water content before land cover change
    !
    ! Should be called BEFORE any subgrid weight updates this time step
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g   ! grid cell index
    !-------------------------------------------------------------------------------
    
    ! initialize heat and water content and dynamic balance fields to zero
    do g = bounds%begg, bounds%endg
       waterstate_vars%liq2_grc(g)           = 0._r8
       waterstate_vars%liq1_grc(g)           = 0._r8
       waterstate_vars%ice2_grc(g)           = 0._r8 
       waterstate_vars%ice1_grc(g)           = 0._r8

       waterflux_vars%qflx_liq_dynbal_grc(g) = 0._r8
       waterflux_vars%qflx_ice_dynbal_grc(g) = 0._r8

       temperature_vars%heat2_grc(g)         = 0._r8
       temperature_vars%heat1_grc(g)         = 0._r8

       energyflux_vars%eflx_dynbal_grc(g)    = 0._r8
    enddo

    call dyn_hwcontent( bounds, &
         waterstate_vars%liq1_grc(bounds%begg:bounds%endg), &
         waterstate_vars%ice1_grc(bounds%begg:bounds%endg), &
         temperature_vars%heat1_grc(bounds%begg:bounds%endg) , &
         urbanparams_vars, soilstate_vars, soilhydrology_vars, &
         temperature_vars, waterstate_vars, lakestate_vars)

  end subroutine dyn_hwcontent_init

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_final(bounds, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Compute grid cell-level heat and water content after land cover change, and compute
    ! the dynbal fluxes
    !
    ! Should be called AFTER all subgrid weight updates this time step
    !
    ! !USES:
    use clm_time_manager    , only : get_step_size
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g     ! grid cell index
    real(r8) :: dtime ! land model time step (sec)
    !---------------------------------------------------------------------------

    call dyn_hwcontent( bounds,                                &
         waterstate_vars%liq2_grc(bounds%begg:bounds%endg),    &
         waterstate_vars%ice2_grc(bounds%begg:bounds%endg),    &
         temperature_vars%heat2_grc(bounds%begg:bounds%endg) , &
         urbanparams_vars, soilstate_vars, soilhydrology_vars, &
         temperature_vars, waterstate_vars, lakestate_vars)

    dtime = get_step_size()
    do g = bounds%begg, bounds%endg
       waterflux_vars%qflx_liq_dynbal_grc (g) = (waterstate_vars%liq2_grc  (g) - waterstate_vars%liq1_grc  (g))/dtime
       waterflux_vars%qflx_ice_dynbal_grc (g) = (waterstate_vars%ice2_grc  (g) - waterstate_vars%ice1_grc  (g))/dtime
       energyflux_vars%eflx_dynbal_grc    (g) = (temperature_vars%heat2_grc(g) - temperature_vars%heat1_grc(g))/dtime
    end do

  end subroutine dyn_hwcontent_final

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent(bounds, gcell_liq, gcell_ice, gcell_heat, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, &
       temperature_vars, waterstate_vars, lakestate_vars)

    ! !DESCRIPTION:
    ! Compute grid-level heat and water content to track conservation with respect to
    ! dynamic land cover.

    ! !USES:
    use landunit_varcon , only : istsoil, istice, istwet, istdlak, istice_mec, istcrop
    use column_varcon   , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
    use clm_varcon      , only : cpice,  cpliq, denh2o
    use clm_varpar      , only : nlevsno, nlevgrnd, nlevurb, nlevlak
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)  :: bounds  
    real(r8)                 , intent(out) :: gcell_liq ( bounds%begg: )   ! [gridcell]
    real(r8)                 , intent(out) :: gcell_ice ( bounds%begg: )   ! [gridcell]
    real(r8)                 , intent(out) :: gcell_heat( bounds%begg: )   ! [gridcell]
    type(urbanparams_type)   , intent(in)  :: urbanparams_vars
    type(soilstate_type)     , intent(in)  :: soilstate_vars
    type(soilhydrology_type) , intent(in)  :: soilhydrology_vars
    type(temperature_type)   , intent(in)  :: temperature_vars
    type(waterstate_type)    , intent(in)  :: waterstate_vars
    type(lakestate_type)     , intent(in)  :: lakestate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: li,lf         ! loop initial/final indicies
    integer  :: ci,cf         ! loop initial/final indicies
    integer  :: pi,pf         ! loop initial/final indicies

    integer  :: g,l,c,p,k     ! loop indicies (grid,lunit,column,pft,vertical level)

    real(r8) :: wtgcell       ! weight relative to grid cell
    real(r8) :: wtcol         ! weight relative to column
    real(r8) :: liq           ! sum of liquid water at column level
    real(r8) :: ice           ! sum of frozen water at column level
    real(r8) :: heat          ! sum of heat content at column level
    real(r8) :: cv            ! heat capacity [J/(m^2 K)]

    integer ,pointer :: nlev_improad(:)  ! number of impervious road layers
    real(r8),pointer :: cv_wall(:,:)     ! thermal conductivity of urban wall
    real(r8),pointer :: cv_roof(:,:)     ! thermal conductivity of urban roof
    real(r8),pointer :: cv_improad(:,:)  ! thermal conductivity of urban impervious road
    integer ,pointer :: snl(:)           ! number of snow layers
    real(r8),pointer :: t_soisno(:,:)    ! soil temperature (Kelvin)
    real(r8),pointer :: h2osno(:)        ! snow water (mm H2O)
    real(r8),pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8),pointer :: h2osoi_ice(:,:)  ! frozen water (kg/m2)
    real(r8),pointer :: watsat(:,:)      ! volumetric soil water at saturation (porosity)
    real(r8),pointer :: csol(:,:)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(r8),pointer :: wa_col(:)        ! water in the unconfined aquifer (mm)
    real(r8),pointer :: dz(:,:)          ! layer depth (m)
    !-------------------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(gcell_liq)  == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(gcell_ice)  == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(gcell_heat) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    snl          => col%snl
    dz           => col%dz
    nlev_improad => urbanparams_vars%nlev_improad
    cv_wall      => urbanparams_vars%cv_wall
    cv_roof      => urbanparams_vars%cv_roof
    cv_improad   => urbanparams_vars%cv_improad
    watsat       => soilstate_vars%watsat_col
    csol         => soilstate_vars%csol_col
    wa_col       => soilhydrology_vars%wa_col
    t_soisno     => temperature_vars%t_soisno_col
    h2osoi_liq   => waterstate_vars%h2osoi_liq_col
    h2osoi_ice   => waterstate_vars%h2osoi_ice_col
    h2osno       => waterstate_vars%h2osno_col

    ! Get relevant sizes

    do g = bounds%begg,bounds%endg ! loop over grid cells
       gcell_liq  (g) = 0.0_r8   ! sum for one grid cell
       gcell_ice  (g) = 0.0_r8   ! sum for one grid cell
       gcell_heat (g) = 0.0_r8   ! sum for one grid cell
    end do

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       ci = lun%coli(l)
       cf = lun%colf(l)
       do c = ci,cf   ! loop over columns

          liq   = 0.0_r8 ! sum for one column
          ice   = 0.0_r8
          heat  = 0.0_r8

          !--- water & ice, above ground only ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop          )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_roof       )  &
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_imperv)  &
               .or. (lun%itype(l) == istdlak                                  )  &
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then

             if ( snl(c) < 0 ) then
                do k = snl(c)+1,0 ! loop over snow layers
                   liq   = liq   + h2osoi_liq(c,k)
                   ice   = ice   + h2osoi_ice(c,k)
                end do
             else                 ! no snow layers exist
                ice = ice + waterstate_vars%h2osno_col(c)
             end if
          end if

          !--- water & ice, below ground only ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop          )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istdlak                                  )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then
             do k = 1,nlevgrnd
                liq   = liq   + h2osoi_liq(c,k)
                ice   = ice   + h2osoi_ice(c,k)
             end do
          end if

          !--- water & ice, below ground, for lakes ---
          if ( lun%itype(l) == istdlak ) then
             do k = 1,nlevlak
                liq   = liq   + (1 - lakestate_vars%lake_icefrac_col(c,k)) * col%dz_lake(c,k) * denh2o
                ice   = ice   +      lakestate_vars%lake_icefrac_col(c,k)  * col%dz_lake(c,k) * denh2o
                ! lake layers do not change thickness when freezing, so denh2o should be used
                ! (thermal properties are appropriately adjusted; see LakeTemperatureMod)
             end do
          end if

          !--- water in aquifer ---
          if ( (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop          )  &
               .or. (lun%itype(l) == istwet                                   )  &
               .or. (lun%itype(l) == istice                                   )  &
               .or. (lun%itype(l) == istice_mec                               )  &           
               .or. (lun%urbpoi(l)          .and. col%itype(c) == icol_road_perv  )) then
             liq = liq + soilhydrology_vars%wa_col(c)
          end if

          !--- water in canopy (at pft level) ---
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then   ! note: soil specified at LU level
             do p = col%pfti(c),col%pftf(c) ! loop over patches
                if (pft%active(p)) then
                   liq = liq + waterstate_vars%h2ocan_patch(p) * pft%wtcol(p)
                end if
             end do
          end if

          !--- heat content, below ground only ---
          if (nlevurb > 0) then
             do k = 1,nlevurb
                if (col%itype(c)==icol_sunwall .OR. col%itype(c)==icol_shadewall) then
                   cv = cv_wall(l,k) * dz(c,k)
                   heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                else if (col%itype(c) == icol_roof) then
                   cv = cv_roof(l,k) * dz(c,k)
                   heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
                end if
             end do
          end if
          do k = 1,nlevgrnd
             if (col%itype(c) /= icol_sunwall .and. col%itype(c) /= icol_shadewall &
                  .and. col%itype(c) /= icol_roof) then
                if (col%itype(c) == icol_road_imperv .and. k >= 1 .and. k <= nlev_improad(l)) then
                   cv = cv_improad(l,k) * dz(c,k)
                else if (lun%itype(l) /= istwet .AND. lun%itype(l) /= istice .AND. lun%itype(l) /= istice_mec) then
                   cv = csol(c,k)*(1-watsat(c,k))*dz(c,k) + (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                else
                   cv = (h2osoi_ice(c,k)*cpice + h2osoi_liq(c,k)*cpliq)
                endif
                heat = heat + cv*t_soisno(c,k) / 1.e6_r8 
             end if
          end do

          !--- heat content, below ground in lake water, for lakes ---
          do k = 1,nlevlak
             if (lun%itype(l) == istdlak) then
                cv = denh2o*col%dz_lake(c,k)*( lakestate_vars%lake_icefrac_col(c,k)*cpice + &
                     (1 - lakestate_vars%lake_icefrac_col(c,k))*cpliq )
                heat = heat + cv*temperature_vars%t_lake_col(c,k) / 1.e6_r8
             end if
          end do

          !--- heat content, above ground only ---
          if ( snl(c) < 0 ) then
             do k = snl(c)+1,0 ! loop over snow layers
                cv = cpliq*h2osoi_liq(c,k) + cpice*h2osoi_ice(c,k)
                heat = heat + cv*t_soisno(c,k) / 1.e6_r8
             end do
          else if ( h2osno(c) > 0.0_r8 .and. lun%itype(l) /= istdlak) then
             ! the heat capacity (not latent heat) of snow without snow layers
             ! is currently ignored in LakeTemperature, so it should be ignored here
             k = 1
             cv = cpice*h2osno(c)
             heat = heat + cv*t_soisno(c,k) / 1.e6_r8
          end if

          !--- scale x/m^2 column-level values into x/m^2 gridcell-level values ---
          gcell_liq  (g) = gcell_liq  (g) + liq   * col%wtgcell(c)
          gcell_ice  (g) = gcell_ice  (g) + ice   * col%wtgcell(c)
          gcell_heat (g) = gcell_heat (g) + heat  * col%wtgcell(c)

       end do ! column loop      
    end do ! landunit loop

  end subroutine dyn_hwcontent

end module dynConsBiogeophysMod
