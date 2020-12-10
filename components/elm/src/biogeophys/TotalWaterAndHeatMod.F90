module TotalWaterAndHeatMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Routines for computing total column water and heat contents
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use elm_varcon         , only : cpice, cpliq, denh2o, tfrz, hfus, aquifer_water_baseline
  use elm_varpar         , only : nlevgrnd, nlevsoi, nlevurb, nlevlak
  use subgridAveMod      , only : p2c
  use SoilHydrologyType  , only : soilhydrology_type
  use WaterstateType     , only : waterstate_type
  use UrbanParamsType    , only : urbanparams_type
  use SoilStateType      , only : soilstate_type
  use TemperatureType    , only : temperature_type
  use LakeStateType      , only : lakestate_type
  use column_varcon      , only : icol_roof, icol_sunwall, icol_shadewall
  use column_varcon      , only : icol_road_perv, icol_road_imperv
  use landunit_varcon    , only : istdlak, istsoil,istcrop,istwet,istice,istice_mec
  use LandunitType       , only : lun_pp
  use ColumnType         , only : col_pp
  use ColumnDataType     , only : col_es, col_ws
  use VegetationType     , only : veg_pp
  use VegetationDataType , only : veg_ws 
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  ! For water (ComputeWaterMass* / ComputeLiqIceMass*): We have separate routines for lake
  ! vs. non-lake because parts of the code call just one or the other.
  !
  ! For heat (ComputeHeat*): We use separate routines for lake vs. non-lake to keep these
  ! routines parallel with the water routines.
  public :: ComputeWaterMassNonLake  ! Compute total water mass of non-lake columns
  public :: ComputeWaterMassLake     ! Compute total water mass of lake columns
  public :: ComputeLiqIceMassNonLake ! Compute total water mass of non-lake columns, separated into liquid and ice
  public :: ComputeLiqIceMassLake    ! Compute total water mass of lake columns, separated into liquid and ice
  public :: ComputeHeatNonLake       ! Compute heat content of non-lake columns
  public :: ComputeHeatLake          ! Compute heat content of lake columns
  public :: AdjustDeltaHeatForDeltaLiq ! Adjusts the change in gridcell heat content due to land cover change to account for the implicit heat flux associated with delta_liq
  public :: LiquidWaterHeat          ! Get the total heat content of some mass of liquid water at a given temperature

  !
  ! !PUBLIC MEMBER DATA:

  ! While some parts of the code work just fine with any heat_base_temp, other parts
  ! currently wouldn't work right if we changed this value. This is all related to the
  ! fact that we don't currently track temperature explicitly for all components of the
  ! system. Specifically:
  !
  ! (1) For liquid water pools that don't have an explicit temperature, we assume a
  !     temperature of heat_base_temp. This is not a terrible assumption for
  !     heat_base_temp = tfrz, but would be a terrible assumption for (e.g.)
  !     heat_base_temp = 0.
  !
  ! (2) In AdjustDeltaHeatForDeltaLiq, we currently don't account for the energy
  !     associated with delta_ice (as we do for delta_liq). This amounts to implicitly
  !     assuming that this ice runoff is at heat_base_temp (this is tied in with the fact
  !     that we don't explicitly track the temperature of runoff). This makes sense for
  !     heat_base_temp = tfrz, but wouldn't make sense for other values of heat_base_temp.
  real(r8), parameter, public :: heat_base_temp = tfrz  ! Base temperature for heat sums [K]

  ! ------------------------------------------------------------------------
  ! The following are public just to support unit testing; they shouldn't be used by other code
  ! ------------------------------------------------------------------------

  ! Minimum and maximum temperatures for the water temperature used by AdjustDeltaHeatForDeltaLiq
  real(r8), parameter :: DeltaLiqMinTemp = tfrz  ! [K]
  real(r8), parameter :: DeltaLiqMaxTemp = tfrz + 35._r8  ! [K]

  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: AccumulateLiquidWaterHeat ! For use by ComputeHeat* routines: accumulate quantities that we need to count for liquid water, for a single column
  private :: TempToHeat                ! For use by ComputeHeat* routines: convert temperature to heat content

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine ComputeWaterMassNonLake(bounds, num_nolakec, filter_nolakec, &
       soilhydrology_inst, waterstate_inst, water_mass)
    !
    ! !DESCRIPTION:
    ! Compute total water mass for all non-lake columns
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds     
    integer                  , intent(in)    :: num_nolakec                ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)          ! column filter for non-lake points
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterstate_type)    , intent(in)    :: waterstate_inst
    real(r8)                 , intent(inout) :: water_mass( bounds%begc: ) ! computed water mass (kg m-2)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: liquid_mass(bounds%begc:bounds%endc)  ! kg m-2
    real(r8) :: ice_mass(bounds%begc:bounds%endc)     ! kg m-2
    integer  :: fc, c

    character(len=*), parameter :: subname = 'ComputeWaterMassNonLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(water_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    call ComputeLiqIceMassNonLake( &
         bounds = bounds, &
         num_nolakec = num_nolakec, &
         filter_nolakec = filter_nolakec, &
         soilhydrology_inst = soilhydrology_inst, &
         waterstate_inst = waterstate_inst, &
         liquid_mass = liquid_mass(bounds%begc:bounds%endc), &
         ice_mass = ice_mass(bounds%begc:bounds%endc))

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       water_mass(c) = liquid_mass(c) + ice_mass(c)
    end do

  end subroutine ComputeWaterMassNonLake

  !-----------------------------------------------------------------------
  subroutine ComputeWaterMassLake(bounds, num_lakec, filter_lakec, &
       waterstate_inst, lakestate_vars, water_mass)
    !
    ! !DESCRIPTION:
    ! Compute total water mass for all lake columns
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds     
    integer                  , intent(in)    :: num_lakec                  ! number of column lake points in column filter
    integer                  , intent(in)    :: filter_lakec(:)            ! column filter for lake points
    type(waterstate_type)    , intent(in)    :: waterstate_inst
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    real(r8)                 , intent(inout) :: water_mass( bounds%begc: ) ! computed water mass (kg m-2)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: liquid_mass(bounds%begc:bounds%endc)  ! kg m-2
    real(r8) :: ice_mass(bounds%begc:bounds%endc)     ! kg m-2
    integer  :: fc, c

    character(len=*), parameter :: subname = 'ComputeWaterMassLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(water_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    call ComputeLiqIceMassLake( &
         bounds = bounds, &
         num_lakec = num_lakec, &
         filter_lakec = filter_lakec, &
         waterstate_inst = waterstate_inst, &
         lakestate_vars  = lakestate_vars, &
         liquid_mass = liquid_mass(bounds%begc:bounds%endc), &
         ice_mass = ice_mass(bounds%begc:bounds%endc))

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       water_mass(c) = liquid_mass(c) + ice_mass(c)
    end do

  end subroutine ComputeWaterMassLake


  !-----------------------------------------------------------------------
  subroutine ComputeLiqIceMassNonLake(bounds, num_nolakec, filter_nolakec, &
       soilhydrology_inst, waterstate_inst, liquid_mass, ice_mass)
    !
    ! !DESCRIPTION:
    ! Compute total water mass for all non-lake columns, separated into liquid and ice
    !
    ! Note: Changes to this routine should generally be accompanied by similar changes
    ! to ComputeHeatNonLake
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds     
    integer                  , intent(in)    :: num_nolakec                 ! number of column non-lake points in column filter
    integer                  , intent(in)    :: filter_nolakec(:)           ! column filter for non-lake points
    type(soilhydrology_type) , intent(in)    :: soilhydrology_inst
    type(waterstate_type)    , intent(in)    :: waterstate_inst
    real(r8)                 , intent(inout) :: liquid_mass( bounds%begc: ) ! computed liquid water mass (kg m-2)
    real(r8)                 , intent(inout) :: ice_mass( bounds%begc: )    ! computed ice mass (kg m-2)
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc, l, p                  ! indices
    logical  :: has_h2o  ! whether this point potentially has water to add
    real(r8) :: h2ocan_col(bounds%begc:bounds%endc)  ! canopy water (mm H2O)
    real(r8) :: snocan_col(bounds%begc:bounds%endc)  ! canopy snow water (mm H2O)
    real(r8) :: liqcan                               ! canopy liquid water (mm H2O)

    character(len=*), parameter :: subname = 'ComputeLiqIceMassNonLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(liquid_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ice_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl          =>    col_pp%snl                        , & ! Input:  [integer  (:)   ]  negative number of snow layers
         
         h2osfc       =>    col_ws%h2osfc     , & ! Input:  [real(r8) (:)   ]  surface water (mm)
         h2osno       =>    col_ws%h2osno     , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2ocan_patch =>    veg_ws%h2ocan   , & ! Input:  [real(r8) (:)   ]  canopy water (mm H2O)
!         snocan_patch =>    waterstate_inst%snocan_patch   , & ! Input:  [real(r8) (:)   ]  canopy snow water (mm H2O)
         h2osoi_ice   =>    col_ws%h2osoi_ice , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq   =>    col_ws%h2osoi_liq , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         total_plant_stored_h2o => col_ws%total_plant_stored_h2o, & 
                                                               ! Input:  [real(r8) (:,:) ] plant internal stored water (mm H2O)
         wa           =>    soilhydrology_inst%wa_col        & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
         )

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       liquid_mass(c) = 0._r8
       ice_mass(c) = 0._r8
    end do

    call p2c(bounds, num_nolakec, filter_nolakec, &
         h2ocan_patch(bounds%begp:bounds%endp), &
         h2ocan_col(bounds%begc:bounds%endc))

    !call p2c(bounds, num_nolakec, filter_nolakec, &
    !     snocan_patch(bounds%begp:bounds%endp), &
    !     snocan_col(bounds%begc:bounds%endc))
    snocan_col(:) = 0._r8

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)

       ! waterstate_inst%snocan_patch and waterstate_inst%liqcan_patch are only set if
       ! we're using snow-on-veg; otherwise they are 0. However, we can rely on
       ! h2ocan_patch being set in all cases, so we can always determine the liquid mass
       ! as (h2ocan - snocan).
       ! Note the difference between liqcan and total_plant_stored_h2o.  The prior
       ! is that which exists on the vegetation canopy surface, the latter is
       ! that which exists within the plant xylems and tissues.  In cases
       ! where FATES hydraulics is not turned on, this total_plant_stored_h2o is
       ! non-changing, and is set to 0 for a trivial solution.

       if (snl(c) < 0) then
          ! Loop over snow layers
          do j = snl(c)+1,0
             liquid_mass(c) = liquid_mass(c) + h2osoi_liq(c,j)
             ice_mass(c) = ice_mass(c) + h2osoi_ice(c,j)
          end do
       else if (h2osno(c) /= 0._r8) then
          ! No explicit snow layers, but there may still be some ice in h2osno (there is
          ! no liquid water in this case)
          ice_mass(c) = ice_mass(c) + h2osno(c)
       end if

    end do

    ! Soil water content
    do j = 1, nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall .or. &
              col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_road_imperv) then
             has_h2o = .false.
          else
             has_h2o = .true.
          end if

          if (has_h2o) then
             liquid_mass(c) = liquid_mass(c) + h2osoi_liq(c,j)
             ice_mass(c) = ice_mass(c) + h2osoi_ice(c,j)
         end if
       end do
    end do

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       l = col_pp%landunit(c)
       if ( (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop          )  &
            .or. (lun_pp%itype(l) == istwet                                   )  &
            .or. (lun_pp%itype(l) == istice                                   )  &
            .or. (lun_pp%itype(l) == istice_mec                               )  &
            .or. (lun_pp%urbpoi(l)          .and. col_pp%itype(c) == icol_road_perv  )) then
          liquid_mass(c) = liquid_mass(c) + wa(c)
       end if
       l = col_pp%landunit(c)

       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then   ! note: soil specified at LU level
          do p = col_pp%pfti(c),col_pp%pftf(c) ! loop over patches
             if (veg_pp%active(p)) then
                liquid_mass(c) = liquid_mass(c) + h2ocan_patch(p) * veg_pp%wtcol(p)
            end if
          end do
       end if
       !liqcan = h2ocan_col(c) - snocan_col(c)
       !liquid_mass(c) = liquid_mass(c) + liqcan + total_plant_stored_h2o(c)
       ice_mass(c)    = ice_mass(c) + snocan_col(c)

       if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_sunwall &
            .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_road_imperv) then
          ! Nothing more to add in this case
       else
          !liquid_mass(c) = liquid_mass(c) +  h2osfc(c)
       end if
    end do
  end associate

  end subroutine ComputeLiqIceMassNonLake

  !-----------------------------------------------------------------------
  subroutine ComputeLiqIceMassLake(bounds, num_lakec, filter_lakec, &
       waterstate_inst, lakestate_vars, liquid_mass, ice_mass)
    !
    ! !DESCRIPTION:
    ! Compute total water mass for all lake columns, separated into liquid and ice
    !
    ! Note: Changes to this routine should generally be accompanied by similar changes
    ! to ComputeHeatLake
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds     
    integer               , intent(in)    :: num_lakec                   ! number of column lake points in column filter
    integer               , intent(in)    :: filter_lakec(:)             ! column filter for lake points
    type(waterstate_type) , intent(in)    :: waterstate_inst
    type(lakestate_type)  , intent(in)    :: lakestate_vars
    real(r8)              , intent(inout) :: liquid_mass( bounds%begc: ) ! computed liquid water mass (kg m-2)
    real(r8)              , intent(inout) :: ice_mass( bounds%begc: )    ! computed ice mass (kg m-2)
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, fc                  ! indices

    character(len=*), parameter :: subname = 'ComputeLiqIceMassLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(liquid_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(ice_mass) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl          =>    col_pp%snl                        , & ! Input:  [integer  (:)   ]  negative number of snow layers
         
         h2osno       =>    col_ws%h2osno     , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         h2osoi_ice   =>    col_ws%h2osoi_ice , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq   =>    col_ws%h2osoi_liq   & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         )

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       liquid_mass(c) = 0._r8
       ice_mass(c) = 0._r8
    end do

    ! Snow water content
    do fc = 1, num_lakec
       c = filter_lakec(fc)
       if (snl(c) < 0) then
          ! Loop over snow layers
          do j = snl(c)+1,0
             liquid_mass(c) = liquid_mass(c) + h2osoi_liq(c,j)
             ice_mass(c) = ice_mass(c) + h2osoi_ice(c,j)
          end do
       else if (h2osno(c) /= 0._r8) then
          ! No explicit snow layers, but there may still be some ice in h2osno (there is
          ! no liquid water in this case)
          ice_mass(c) = ice_mass(c) + h2osno(c)
       end if
    end do

    ! Soil water content of the soil under the lake
    do j = 1, nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)
          liquid_mass(c) = liquid_mass(c) + h2osoi_liq(c,j)
          ice_mass(c) = ice_mass(c) + h2osoi_ice(c,j)
       end do
    end do

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       do j = 1,nlevlak
          liquid_mass(c) = liquid_mass(c) + (1 - lakestate_vars%lake_icefrac_col(c,j)) * col_pp%dz_lake(c,j) * denh2o
          ice_mass(c)    = ice_mass(c)    +      lakestate_vars%lake_icefrac_col(c,j)  * col_pp%dz_lake(c,j) * denh2o
          ! lake layers do not change thickness when freezing, so denh2o should be used
          ! (thermal properties are appropriately adjusted; see LakeTemperatureMod)
       end do
    end do
  end associate

  end subroutine ComputeLiqIceMassLake

  !-----------------------------------------------------------------------
  subroutine ComputeHeatNonLake(bounds, num_nolakec, filter_nolakec, &
       urbanparams_inst, soilstate_inst, &
       temperature_inst, waterstate_inst, soilhydrology_inst, &
       heat, heat_liquid, cv_liquid)
    !
    ! !DESCRIPTION:
    ! Compute total heat content for all non-lake columns.
    !
    ! Optionally, also return the total heat content just of liquid water for each column
    ! (excluding latent heat), and/or the total heat capacity just of liquid water for
    ! each column. Together, these can be used by the caller to compute the weighted
    ! average liquid water temperature (with weightings done by the water mass).
    !
    ! Note: Changes to this routine - for anything involving liquid or ice - should
    ! generally be accompanied by similar changes to ComputeLiqIceMassNonLake
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)  :: bounds
    integer                  , intent(in)  :: num_nolakec
    integer                  , intent(in)  :: filter_nolakec(:)
    type(urbanparams_type)   , intent(in)  :: urbanparams_inst
    type(soilstate_type)     , intent(in)  :: soilstate_inst
    type(temperature_type)   , intent(in)  :: temperature_inst
    type(waterstate_type)    , intent(in)  :: waterstate_inst
    type(soilhydrology_type) , intent(in)  :: soilhydrology_inst

    real(r8) , intent(inout) :: heat( bounds%begc: )        ! sum of heat content for all columns [J/m^2]
    real(r8) , intent(inout) :: heat_liquid( bounds%begc: ) ! sum of heat content for all columns: liquid water, excluding latent heat [J/m^2]
    real(r8) , intent(inout) :: cv_liquid( bounds%begc: )   ! sum of liquid heat capacity for all columns [J/(m^2 K)]
    !
    ! !LOCAL VARIABLES:
    integer :: fc
    integer :: l,c,j

    logical  :: has_h2o  ! whether this point potentially has water to add

    real(r8) :: h2ocan_col(bounds%begc:bounds%endc)  ! canopy water (mm H2O)
    real(r8) :: snocan_col(bounds%begc:bounds%endc)  ! canopy snow water (mm H2O)
    real(r8) :: liqcan        ! canopy liquid water (mm H2O)

    real(r8) :: heat_dry_mass(bounds%begc:bounds%endc) ! sum of heat content: dry mass [J/m^2]
    real(r8) :: heat_ice(bounds%begc:bounds%endc)      ! sum of heat content: ice [J/m^2]
    real(r8) :: latent_heat_liquid(bounds%begc:bounds%endc) ! sum of latent heat content of liquid water [J/m^2]

    character(len=*), parameter :: subname = 'ComputeHeatNonLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(heat) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(heat_liquid) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(cv_liquid) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl          => col_pp%snl, & ! number of snow layers
         dz           => col_pp%dz, &  ! layer depth (m)
         nlev_improad => urbanparams_inst%nlev_improad, & ! number of impervious road layers
         cv_wall      => urbanparams_inst%cv_wall, & ! heat capacity of urban wall (J/m^3/K)
         cv_roof      => urbanparams_inst%cv_roof, & ! heat capacity of urban roof (J/m^3/K)
         cv_improad   => urbanparams_inst%cv_improad, & ! heat capacity of urban impervious road (J/m^3/K)
         watsat       => soilstate_inst%watsat_col, & ! volumetric soil water at saturation (porosity)
         csol         => soilstate_inst%csol_col, & ! heat capacity, soil solids (J/m**3/Kelvin)
         t_soisno     => col_es%t_soisno, & ! soil temperature (Kelvin)
         t_h2osfc     => col_es%t_h2osfc, & ! surface water temperature (Kelvin)
         h2osoi_liq   => col_ws%h2osoi_liq, & ! liquid water (kg/m2)
         h2osoi_ice   => col_ws%h2osoi_ice, & ! frozen water (kg/m2)
         h2osno       => col_ws%h2osno, & ! snow water (mm H2O)
         h2osfc       => col_ws%h2osfc, & ! surface water (mm H2O)
         h2ocan_patch => veg_ws%h2ocan, & ! canopy water (mm H2O)
!         snocan_patch => waterstate_inst%snocan_patch, & ! canopy snow water (mm H2O)
         total_plant_stored_h2o_col => col_ws%total_plant_stored_h2o, & ! Input: [real(r8) (:)   ]  water mass in plant tissues (kg m-2)
         wa           => soilhydrology_inst%wa_col & ! water in the unconfined aquifer (mm)
         )

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)

       heat_liquid(c) = 0._r8
       cv_liquid(c) = 0._r8
       heat_dry_mass(c) = 0._r8
       heat_ice(c) = 0._r8
       latent_heat_liquid(c) = 0._r8
    end do

    call p2c(bounds, &
         parr = h2ocan_patch(bounds%begp:bounds%endp), &
         carr = h2ocan_col(bounds%begc:bounds%endc), &
         p2c_scale_type = 'unity')

    !call p2c(bounds, &
    !     parr = snocan_patch(bounds%begp:bounds%endp), &
    !     carr = snocan_col(bounds%begc:bounds%endc), &
    !     p2c_scale_type = 'unity')
    snocan_col(bounds%begc:bounds%endc) = 0._r8

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)

       !--- canopy water ---
       !
       ! TODO(wjs, 2017-03-11) Canopy water currently doesn't have an explicit
       ! temperature; thus, we only add its latent heat of fusion. Eventually, we should
       ! probably track its temperature explicitly - or at least give it an implicit
       ! temperature for the sake of these energy calculations (I think that's needed for
       ! full conservation).
       !
       ! However, we still call the standard AccumulateLiquidWaterHeat routine, so that we
       ! average in a heat_base_temp value in heat_liquid. I think this will generally
       ! lead to less of a sensible heat flux adjustment needed by the dynbal energy
       ! conservation code. (But I went back and forth on whether to do this, so could
       ! be convinced otherwise.)

       ! snocan and liqcan are only set if we're using snow-on-veg; otherwise they are 0.
       ! However, we can rely on h2ocan being set in all cases, so we can always
       ! determine the liquid mass as (h2ocan - snocan).
       
       ! Note (rgk 04-2017): added total_plant_stored_h2o_col(c), which is the
       ! water inside the plant, which is zero for all non-dynamic models. FATES hydraulics
       ! is the only one with dynamic storage atm.
       ! Commentary (rgk 04-2017): water has moved from the soil to the plant tissues,
       ! and the two pools have different temperatures associated with them. However,
       ! we are not accounting for or conserving the flux of energy between the two
       ! pools.  The energy in the plant water should "bring with it" the internal
       ! energy of the soil-to-root water flux.

       liqcan = (h2ocan_col(c) - snocan_col(c) + total_plant_stored_h2o_col(c))*0._r8
       call AccumulateLiquidWaterHeat( &
            temp = heat_base_temp, &
            h2o = liqcan, &
            cv_liquid = cv_liquid(c), &
            heat_liquid = heat_liquid(c), &
            latent_heat_liquid = latent_heat_liquid(c))

       !--- snow ---
       if ( snl(c) < 0 ) then
          ! Loop over snow layers
          do j = snl(c)+1,0
             call AccumulateLiquidWaterHeat( &
                  temp = t_soisno(c,j), &
                  h2o = h2osoi_liq(c,j), &
                  cv_liquid = cv_liquid(c), &
                  heat_liquid = heat_liquid(c), &
                  latent_heat_liquid = latent_heat_liquid(c))
             heat_ice(c) = heat_ice(c) + &
                  TempToHeat(temp = t_soisno(c,j), cv = (h2osoi_ice(c,j)*cpice))
          end do
       else if (h2osno(c) /= 0._r8) then
          ! No explicit snow layers, but there may still be some ice in h2osno (there is
          ! no liquid water in this case)
          j = 1
          heat_ice(c) = heat_ice(c) + &
               TempToHeat(temp = t_soisno(c,j), cv = (h2osno(c)*cpice))
       end if

       if (col_pp%hydrologically_active(c)) then
          ! NOTE(wjs, 2017-03-23) Water in the unconfined aquifer currently doesn't have
          ! an explicit temperature; thus, we only add its latent heat of
          ! fusion. However, we still call the standard AccumulateLiquidWaterHeat routine, so
          ! that we average in a heat_base_temp value in heat_liquid. I think this will
          ! generally lead to less of a sensible heat flux adjustment needed by the
          ! dynbal energy conservation code. (But I went back and forth on whether to do
          ! this, so could be convinced otherwise.) In the default CLM5 configuration,
          ! this should all be irrelevant, because (wa(c) - aquifer_water_baseline)
          ! should be fixed at 0 for all hydrologically-active points.

          call AccumulateLiquidWaterHeat( &
               temp = heat_base_temp, &
               h2o = (wa(c) - aquifer_water_baseline), &
               cv_liquid = cv_liquid(c), &
               heat_liquid = heat_liquid(c), &
               latent_heat_liquid = latent_heat_liquid(c))
       end if

       if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_sunwall &
            .or. col_pp%itype(c) == icol_shadewall .or. col_pp%itype(c) == icol_road_imperv) then
          ! Nothing more to add in this case
       else
          !--- surface water ---
          call AccumulateLiquidWaterHeat( &
               temp = t_h2osfc(c), &
               h2o = h2osfc(c), &
               cv_liquid = cv_liquid(c), &
               heat_liquid = heat_liquid(c), &
               latent_heat_liquid = latent_heat_liquid(c))
       end if

    end do


    !--- below ground (soil & soil water) and related urban columns
    do j = 1, nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)
          l = col_pp%landunit(c)

          if (col_pp%itype(c)==icol_sunwall .or. col_pp%itype(c)==icol_shadewall) then
             has_h2o = .false.
             if (j <= nlevurb) then
                heat_dry_mass(c) = heat_dry_mass(c) + &
                     TempToHeat(temp = t_soisno(c,j), cv = (cv_wall(l,j) * dz(c,j)))
             end if

          else if (col_pp%itype(c) == icol_roof) then
             if (j <= nlevurb) then
                has_h2o = .true.
                heat_dry_mass(c) = heat_dry_mass(c) + &
                     TempToHeat(temp = t_soisno(c,j), cv = (cv_roof(l,j) * dz(c,j)))
             else
                has_h2o = .false.
             end if

          else
             has_h2o = .true.

             if (col_pp%itype(c) == icol_road_imperv .and. j <= nlev_improad(l)) then
                heat_dry_mass(c) = heat_dry_mass(c) + &
                     TempToHeat(temp = t_soisno(c,j), cv = (cv_improad(l,j) * dz(c,j)))
             else if (lun_pp%itype(l) /= istwet .and. lun_pp%itype(l) /= istice .and. lun_pp%itype(l) /= istice_mec) then
                ! Note that this also includes impervious roads below nlev_improad (where
                ! we have soil)
                heat_dry_mass(c) = heat_dry_mass(c) + &
                     TempToHeat(temp = t_soisno(c,j), cv = (csol(c,j)*(1-watsat(c,j))*dz(c,j)))
             end if
          end if

          if (has_h2o) then
             call AccumulateLiquidWaterHeat( &
                  temp = t_soisno(c,j), &
                  h2o = h2osoi_liq(c,j), &
                  cv_liquid = cv_liquid(c), &
                  heat_liquid = heat_liquid(c), &
                  latent_heat_liquid = latent_heat_liquid(c))
             heat_ice(c) = heat_ice(c) + &
                  TempToHeat(temp = t_soisno(c,j), cv = (h2osoi_ice(c,j)*cpice))
          end if
       end do
    end do

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       heat(c) = heat_dry_mass(c) + heat_ice(c) + heat_liquid(c) + latent_heat_liquid(c)
    end do

    end associate

  end subroutine ComputeHeatNonLake

  !-----------------------------------------------------------------------
  subroutine ComputeHeatLake(bounds, num_lakec, filter_lakec, &
       soilstate_inst, temperature_inst, waterstate_inst, &
       heat, heat_liquid, cv_liquid)
    !
    ! !DESCRIPTION:
    ! Compute total heat content for all lake columns
    !
    ! Optionally, also return the total heat content just of liquid water for each column
    ! (excluding latent heat), and/or the total heat capacity just of liquid water for
    ! each column. Together, these can be used by the caller to compute the weighted
    ! average liquid water temperature (with weightings done by the water mass).
    !
    ! Note: Changes to this routine - for anything involving liquid or ice - should
    ! generally be accompanied by similar changes to ComputeLiqIceMassLake
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)  :: bounds
    integer                  , intent(in)  :: num_lakec
    integer                  , intent(in)  :: filter_lakec(:)
    type(soilstate_type)     , intent(in)  :: soilstate_inst
    type(temperature_type)   , intent(in)  :: temperature_inst
    type(waterstate_type)    , intent(in)  :: waterstate_inst

    real(r8) , intent(inout) :: heat( bounds%begc: )        ! sum of heat content for all columns [J/m^2]
    real(r8) , intent(inout) :: heat_liquid( bounds%begc: ) ! sum of heat content for all columns: liquid water, excluding latent heat [J/m^2]
    real(r8) , intent(inout) :: cv_liquid( bounds%begc: )   ! sum of liquid heat capacity for all columns [J/(m^2 K)]
    !
    ! !LOCAL VARIABLES:
    integer :: fc
    integer :: c,j

    real(r8) :: heat_dry_mass(bounds%begc:bounds%endc) ! sum of heat content: dry mass [J/m^2]
    real(r8) :: heat_ice(bounds%begc:bounds%endc)      ! sum of heat content: ice [J/m^2]
    real(r8) :: latent_heat_liquid(bounds%begc:bounds%endc) ! sum of latent heat content of liquid water [J/m^2]

    character(len=*), parameter :: subname = 'ComputeHeatLake'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(heat) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(heat_liquid) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(cv_liquid) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         snl          => col_pp%snl, & ! number of snow layers
         dz           => col_pp%dz, &  ! layer depth (m)
         watsat       => soilstate_inst%watsat_col, & ! volumetric soil water at saturation (porosity)
         csol         => soilstate_inst%csol_col, & ! heat capacity, soil solids (J/m**3/Kelvin)
         t_soisno     => col_es%t_soisno, & ! soil temperature (Kelvin)
         h2osoi_liq   => col_ws%h2osoi_liq, & ! liquid water (kg/m2)
         h2osoi_ice   => col_ws%h2osoi_ice, & ! frozen water (kg/m2)
         h2osno       => col_ws%h2osno & ! snow water (mm H2O)
         )

    do fc = 1, num_lakec
       c = filter_lakec(fc)

       heat_liquid(c) = 0._r8
       cv_liquid(c) = 0._r8
       heat_dry_mass(c) = 0._r8
       heat_ice(c) = 0._r8
       latent_heat_liquid(c) = 0._r8
    end do

    ! Snow heat content
    do fc = 1, num_lakec
       c = filter_lakec(fc)
       if ( snl(c) < 0 ) then
          ! Loop over snow layers
          do j = snl(c)+1,0
             call AccumulateLiquidWaterHeat( &
                  temp = t_soisno(c,j), &
                  h2o = h2osoi_liq(c,j), &
                  cv_liquid = cv_liquid(c), &
                  heat_liquid = heat_liquid(c), &
                  latent_heat_liquid = latent_heat_liquid(c))
             heat_ice(c) = heat_ice(c) + &
                  TempToHeat(temp = t_soisno(c,j), cv = (h2osoi_ice(c,j)*cpice))
          end do
       else if (h2osno(c) /= 0._r8) then
          ! TODO(wjs, 2017-03-16) (Copying this note from old code... I'm not positive
          ! it's still true.) The heat capacity (not latent heat) of snow without snow
          ! layers is currently ignored in LakeTemperature, so it should be ignored here.
          ! Eventually we should consider this.
       end if
    end do

    ! Soil water content of the soil under the lake
    do j = 1,nlevgrnd
       do fc = 1, num_lakec
          c = filter_lakec(fc)

          heat_dry_mass(c) = heat_dry_mass(c) + &
               TempToHeat(temp = t_soisno(c,j), cv = (csol(c,j)*(1-watsat(c,j))*dz(c,j)))
          call AccumulateLiquidWaterHeat( &
               temp = t_soisno(c,j), &
               h2o = h2osoi_liq(c,j), &
               cv_liquid = cv_liquid(c), &
               heat_liquid = heat_liquid(c), &
               latent_heat_liquid = latent_heat_liquid(c))
          heat_ice(c) = heat_ice(c) + &
               TempToHeat(temp = t_soisno(c,j), cv = (h2osoi_ice(c,j)*cpice))
       end do
    end do

    ! TODO(wjs, 2017-03-11) Include heat content of water in lakes, once we include
    ! lake water as an explicit water state (https://github.com/NCAR/CLM/issues/2)

    do fc = 1, num_lakec
       c = filter_lakec(fc)
       heat(c) = heat_dry_mass(c) + heat_ice(c) + heat_liquid(c) + latent_heat_liquid(c)
    end do

    end associate

  end subroutine ComputeHeatLake

  !-----------------------------------------------------------------------
  subroutine AdjustDeltaHeatForDeltaLiq(bounds, delta_liq, &
       liquid_water_temp1, liquid_water_temp2, &
       delta_heat)
    !
    ! !DESCRIPTION:
    ! Adjusts delta_heat (the change in gridcell heat content due to land cover change
    ! that needs to be accounted for via a heat flux) to account for the implicit heat
    ! flux associated with delta_liq.
    !
    ! Note that, throughout CLM, we don't explicitly track the temperature or heat content
    ! of runoff. Furthermore, we currently cannot compute the exact heat content of
    ! delta_liq (the dynamic landcover adjustment), because we aren't summing the liquid
    ! water heat content on a pool-by-pool (and layer-by-layer) basis, but rather on a
    ! bulk basis across each column. Thus, the formulation in this routine is currently
    ! using a rough approximation of the temperature of delta_liq - assuming it is at the
    ! average temperature of the liquid water in the grid cell. This can be a poor
    ! assumption in some cases (e.g., if the grid cell is 90% glacier, 5% natural veg and
    ! 5% crop, and the only transitions are between natural veg and crop - then the
    ! glacier's liquid water temperature factors into the average liquid water
    ! temperature, even though it doesn't contribute at all to delta_liq).
    !
    ! Also note that we don't account for delta_ice here. This implicitly assumes that
    ! ice runoff is at heat_base_temp (which is reasonable as long as heat_base_temp =
    ! tfrz).
    !
    ! Eventually, if we begin to explicitly account for the temperature / heat content of
    ! liquid and ice runoff in CLM, then this routine should be reworked to use the true
    ! heat contents of both liquid and ice runoff.
    !
    ! Sign convention: delta_liq and delta_heat are positive if the post-landcover change
    ! value is greater than the pre-landcover change value.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: delta_liq( bounds%begg: )  ! change in gridcell h2o liq content [kg/m^2]
    real(r8), intent(in) :: liquid_water_temp1( bounds%begg: ) ! average liquid water temperature before land cover change [K]
    real(r8), intent(in) :: liquid_water_temp2( bounds%begg: ) ! average liquid water temperature after land cover change [K]
    real(r8), intent(inout) :: delta_heat( bounds%begg: ) ! change in gridcell heat content [J/m^2]
    !
    ! !LOCAL VARIABLES:
    integer :: g
    real(r8) :: water_temperature  ! [K]
    real(r8) :: total_liquid_heat ! [J/m^2]
    real(r8) :: heat_liquid ! [J/m^2]
    real(r8) :: latent_heat_liquid ! [J/m^2]
    real(r8) :: cv ! heat capacity [J/(m^2 K)]

    character(len=*), parameter :: subname = 'AdjustDeltaHeatForDeltaLiq'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(delta_liq) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(liquid_water_temp1) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(liquid_water_temp2) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(delta_heat) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))

    do g = bounds%begg, bounds%endg
       if (delta_liq(g) /= 0._r8) then
          if (delta_liq(g) < 0._r8) then
             ! There was more water in the initial state than in the final state. We'll
             ! generate a positive runoff. We assume that the runoff has a temperature equal
             ! to the average temperature of liquid water in the initial state.
             water_temperature = liquid_water_temp1(g)
          else
             ! There is more water in the final state than in the initial state. We'll
             ! generate a negative runoff. We assume that we're sucking water out of the
             ! ocean at a temperature equal to the average temperature of liquid water in
             ! the final state.
             water_temperature = liquid_water_temp2(g)
          end if

          ! Since we're not trying to completely conserve energy here, it's better to
          ! ensure that the estimated water temperature is in some reasonable bounds.
          ! This protects against getting bad temperatures as a result of something like
          ! catastrophic cancellation, or the weirdness that can arise from having
          ! negative water volumes included in the averages.
          water_temperature = max(water_temperature, DeltaLiqMinTemp)
          water_temperature = min(water_temperature, DeltaLiqMaxTemp)

          total_liquid_heat = LiquidWaterHeat( &
               temp = water_temperature, &
               h2o = delta_liq(g))

          ! For delta_liq < 0 (liq2 < liq1): We'll generate a positive runoff; we want to
          ! effectively include some positive heat from that positive runoff in the heat2
          ! state, which means adding a positive term to delta_heat. Since the above heat
          ! quantities will be negative, we need to subtract them. The reverse is true
          ! for delta_liq > 0; again, we need to subtract the heat quantities.
          delta_heat(g) = delta_heat(g) - total_liquid_heat

       end if
    end do

  end subroutine AdjustDeltaHeatForDeltaLiq

  !-----------------------------------------------------------------------
  function LiquidWaterHeat(temp, h2o) result(heat)
    !
    ! !DESCRIPTION:
    ! Get the total heat content (including latent heat) of some mass of liquid water at
    ! a given temperature, using a base temperature of heat_base_temp.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: heat  ! function result
    real(r8), intent(in) :: temp  ! temperature [K]
    real(r8), intent(in) :: h2o   ! water mass [kg/m^2]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: heat_liquid ! heat content of liquid water, excluding latent heat [J/m^2]
    real(r8) :: latent_heat_liquid ! latent heat content of liquid water [J/m^2]

    character(len=*), parameter :: subname = 'LiquidWaterHeat'
    !-----------------------------------------------------------------------

    heat_liquid = 0._r8
    latent_heat_liquid = 0._r8
    call AccumulateLiquidWaterHeat(temp = temp, h2o = h2o, &
         heat_liquid = heat_liquid, latent_heat_liquid = latent_heat_liquid)

    heat = heat_liquid + latent_heat_liquid

  end function LiquidWaterHeat


  !-----------------------------------------------------------------------
  subroutine AccumulateLiquidWaterHeat(temp, h2o, &
       heat_liquid, latent_heat_liquid, cv_liquid)
    !
    ! !DESCRIPTION:
    ! In the course of accumulating heat contents: Accumulate quantities that we need to
    ! count for liquid water, for a single column
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: temp  ! temperature [K]
    real(r8), intent(in) :: h2o   ! water mass [kg/m^2]

    real(r8), intent(inout) :: heat_liquid        ! accumulated total heat content of liquid water for this column, excluding latent heat [J/m^2]
    real(r8), intent(inout) :: latent_heat_liquid ! accumulated total latent heat content of liquid water for this column [J/m^2]
    real(r8), intent(inout), optional :: cv_liquid ! accumulated total liquid heat capacity for this column [J/(m^2 K)]
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cv  ! heat capacity [J/(m^2 K)]

    character(len=*), parameter :: subname = 'AccumulateLiquidWaterHeat'
    !-----------------------------------------------------------------------

    cv = h2o*cpliq
    if (present(cv_liquid)) then
       cv_liquid = cv_liquid + cv
    end if
    heat_liquid = heat_liquid + TempToHeat(temp = temp, cv = cv)
    latent_heat_liquid = latent_heat_liquid + h2o*hfus
  end subroutine AccumulateLiquidWaterHeat

  !-----------------------------------------------------------------------
  pure function TempToHeat(temp, cv) result(heat)
    !
    ! !DESCRIPTION:
    ! Convert temperature to heat content
    !
    ! !ARGUMENTS:
    real(r8) :: heat  ! function result: heat in J/m^2
    real(r8), intent(in) :: temp  ! temperature [K]
    real(r8), intent(in) :: cv    ! heat capacity [J/(m^2 K)]
    !-----------------------------------------------------------------------

    heat = cv*(temp - heat_base_temp)

  end function TempToHeat

end module TotalWaterAndHeatMod
