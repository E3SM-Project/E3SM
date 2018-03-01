module BeTRSimulationCLM
  !
  ! !DESCRIPTION:
  !  CLM-BeTR interface
  !
  ! !USES:
  !
#include "shr_assert.h"

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use EcophysConType      , only : ecophyscon_type
  use BeTRSimulation      , only : betr_simulation_type
  use betr_decompMod      , only : betr_bounds_type
  use BeTR_TimeMod        , only : betr_time_type
  use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_clm_type


     ! NOTE(bja, 201603) CLM stubb types go here!

   contains
     procedure, public :: InitOnline                => CLMInit
     procedure, public :: Init                      => CLMInitOffline
     procedure, public :: StepWithoutDrainage       => CLMStepWithoutDrainage
     procedure, public :: StepWithDrainage          => CLMStepWithDrainage
     procedure, public :: SetBiophysForcing         => CLMSetBiophysForcing
     !clm unique subroutines
     procedure, public :: ConsistencyCheck          => clm_h2oiso_consistency_check
     procedure, public :: DiagnoseDtracerFreezeThaw => CLMDiagnoseDtracerFreezeThaw
     procedure, public :: CalcDewSubFlux            => CLMCalcDewSubFlux
     procedure, public :: SoilFluxStateRecv         => CLMBetrSoilFluxStateRecv
     procedure, public :: CalcSmpL                  => CLMCalcSmpL
  end type betr_simulation_clm_type

  public :: create_betr_simulation_clm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_clm()
    !DESCRIPTION
    ! constructor
    implicit none

    class(betr_simulation_clm_type), pointer :: create_betr_simulation_clm
    class(betr_simulation_clm_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_clm => simulation

  end function create_betr_simulation_clm

  !-------------------------------------------------------------------------------

  subroutine CLMInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, masterproc)
    !DESCRIPTION
    !initialize interface
    !
    !USES
    use betr_constants      , only : betr_namelist_buffer_size
    use betr_constants      , only : betr_filename_length
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use LandunitType        , only : landunit_type
    use ColumnType          , only : column_type
    use PatchType           , only : patch_type
    use LandunitType        , only : landunit_type
    use pftvarcon           , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType      , only : waterstate_type
    use landunit_varcon     , only : istcrop, istice, istsoil
    use BeTR_landvarconType , only : betr_landvarcon
    use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar          , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type)          , intent(inout) :: this
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
    type(patch_type)                         , intent(in) :: pft
    type(waterstate_type)                    , intent(inout) :: waterstate
    logical,                        optional , intent(in)    :: masterproc
    betr_nlevsoi                       = nlevsoi
    betr_nlevsno                       = nlevsno
    betr_nlevtrc_soil                  = nlevtrc_soil


    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice


    ! allocate the reaction types that may only be known to this
    ! simulation type.
    ! now call the base simulation init to continue initialization
    if(present(masterproc))then
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, masterproc=masterproc)
    else
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer)
    endif
  end subroutine CLMInit

  !-------------------------------------------------------------------------------

  subroutine CLMInitOffline(this, bounds, lun, col, pft, waterstate, namelist_buffer,base_filename)
    !DESCRIPTION
    !initialize interface
    !
    !USES
    use betr_constants      , only : betr_namelist_buffer_size
    use betr_constants      , only : betr_filename_length
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use LandunitType        , only : landunit_type
    use ColumnType          , only : column_type
    use PatchType           , only : patch_type
    use LandunitType        , only : landunit_type
    use pftvarcon           , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType      , only : waterstate_type
    use landunit_varcon     , only : istcrop, istice, istsoil
    use BeTR_landvarconType , only : betr_landvarcon
    use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar          , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type)          , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
    type(patch_type)                         , intent(in) :: pft
    type(waterstate_type)                    , intent(inout) :: waterstate

    betr_nlevsoi                       = nlevsoi
    betr_nlevsno                       = nlevsno
    betr_nlevtrc_soil                  = nlevtrc_soil


    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice


    ! allocate the reaction types that may only be known to this
    ! simulation type.
    ! now call the base simulation init to continue initialization
    call this%BeTRInit(bounds, lun, col, pft, waterstate,namelist_buffer,base_filename )

  end subroutine CLMInitOffline
  !---------------------------------------------------------------------------------
  subroutine CLMStepWithoutDrainage(this, bounds, col, pft)
   !DESCRIPTION
   !march one step without drainage
   !
   !USES
    use ColumnType        , only : column_type
    use tracer_varcon     , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use PatchType         , only : patch_type
    use LandunitType      , only : landunit_type
    use clm_varpar        , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    ! !ARGUMENTS :
    class(betr_simulation_clm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds ! bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(patch_type)                , intent(in)    :: pft
    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer                :: c

    !pass necessary data for correct subroutine call
    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%bsimstatus%reset()

    call this%BeTRSetcps(bounds, col, pft)
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds,  this%betr_col(c), &
         this%betr_pft(c), this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=this%bsimstatus%print_msg())
  end subroutine CLMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine CLMStepWithDrainage(this, bounds, col)
   !DESCRIPTION
   !march one step with drainage
  !
  !USES
    use ColumnType     , only : column_type
    use MathfuncMod    , only : safe_div
    use WaterFluxType  , only : waterflux_type
    use betr_decompMod , only : betr_bounds_type
    use LandunitType   , only : landunit_type
    use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use clm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_clm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col ! column type

    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer           :: c
    call this%bsimstatus%reset()

    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)

    do c = bounds%begc, bounds%endc
       if(.not. this%active_col(c))cycle
       call this%betr(c)%step_with_drainage(betr_bounds, this%betr_col(c),   &
         this%num_soilc, this%filter_soilc, this%jtops, &
         this%biogeo_flux(c), this%bstatus(c))
       if(this%bstatus(c)%check_status())then
         call this%bsimstatus%setcol(c)
         call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
         exit
       endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=this%bsimstatus%print_msg())
  end subroutine CLMStepWithDrainage

  !------------------------------------------------------------------------

  subroutine CLMDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec, col, lun)
    !
    ! DESCRIPTION
    ! aqueous tracer partition based on freeze-thaw
    !
    ! USES
    use ColumnType      , only : column_type
    use LandunitType    , only : landunit_type
    use landunit_varcon , only : istsoil
    use betr_decompMod  , only : betr_bounds_type

    implicit none
    !!ARGUMENTS
    class(betr_simulation_clm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_nolakec ! number of column non-lake points in column filter
    integer                         , intent(in)    :: filter_nolakec(:) ! column filter for non-lake points
    type(column_type)               , intent(in)    :: col ! column type
    type(landunit_type)             , intent(in)    :: lun

    !temporary variables
    type(betr_bounds_type)     :: betr_bounds
    integer :: fc, c

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)
   do fc = 1, num_nolakec
     c = filter_nolakec(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_dtracer_freeze_thaw(betr_bounds,  &
      this%num_soilc, this%filter_soilc,  this%biophys_forc(c))
  enddo
  end subroutine CLMDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine CLMCalcDewSubFlux(this,  &
       bounds, col, num_hydrologyc, filter_soilc_hydrologyc)
    !DESCRIPTION
    ! External interface called by CLM
    !
    !USES
    use LandunitType    , only : landunit_type
    use clm_varcon      , only : denh2o,spval
    use landunit_varcon , only : istsoil, istcrop
    use betr_decompMod  , only : betr_bounds_type
    use BeTR_TimeMod    , only : betr_time_type
    use ColumnType        , only : column_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col
    integer                         , intent(in)    :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer                         , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    integer :: fc, c

    type(betr_bounds_type)     :: betr_bounds

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)
    do fc = 1, num_hydrologyc
      c = filter_soilc_hydrologyc(fc)
      if(.not. this%active_col(c))cycle
      call this%betr(c)%calc_dew_sub_flux(this%betr_time, betr_bounds, this%betr_col(c), &
       this%num_soilc, this%filter_soilc, this%biophys_forc(c), this%betr(c)%tracers, &
       this%betr(c)%tracerfluxes, this%betr(c)%tracerstates)
    enddo

  end subroutine CLMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine CLMBetrSoilFluxStateRecv(this,  num_soilc, filter_soilc)
   !DESCRIPTION
   !this is to expaneded
   !
   !USES
    use betr_decompMod    , only : betr_bounds_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type) , intent(inout) :: this
    integer                         , intent(in)    :: num_soilc
    integer                         , intent(in)    :: filter_soilc(:)

    integer :: fc, c
    type(betr_bounds_type)     :: betr_bounds

    call this%BeTRSetBounds(betr_bounds)

!    do fc = 1, num_soilc
!      c = filter_soilc(fc)
!      if(.not. this%active_col(c))cycle
!      call this%betr(c)%bgc_reaction%lsm_betr_flux_state_receive(betr_bounds, &
!         this%num_soilc, this%filter_soilc,                                   &
!         this%betr(c)%tracerstates, this%betr(c)%tracerfluxes,  this%betr(c)%tracers)
!    enddo
  end subroutine CLMBetrSoilFluxStateRecv

  !------------------------------------------------------------------------
  subroutine CLMCalcSmpL(this, bounds, lbj, ubj, &
       numf, filter, t_soisno, soilstate_vars, waterstate_vars, soil_water_retention_curve)
   !DESCRIPTION
   ! calculate water suction potential
   !
   !USES
    use SoilStateType              , only : soilstate_type
    use WaterStateType             , only : waterstate_type
    use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
    use betr_decompMod             , only : betr_bounds_type
    use clm_varcon                 , only : grav,hfus,tfrz
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type)        , intent(inout)        :: this
    type(bounds_type)                      , intent(in)           :: bounds  ! bounds
    integer                                , intent(in)           :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer                                , intent(in)           :: numf                                              ! number of columns in column filter
    integer                                , intent(in)           :: filter(:)                                         ! column filter
    real(r8)                               , intent(in)           :: t_soisno(bounds%begc:, lbj: )                    ! soil temperature
    type(soilstate_type)                   , intent(in)           :: soilstate_vars
    type(waterstate_type)                  , intent(inout)        :: waterstate_vars
    class(soil_water_retention_curve_type) , intent(in), optional :: soil_water_retention_curve

    !temporary variables
    real(r8) :: s_node
    integer :: fc, c, j

    SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

    ! humor the compiler about unused variables
    if (this%num_soilc > 0)                  continue
    if (present(soil_water_retention_curve)) continue

    associate(                                            & !
         h2osoi_vol =>    waterstate_vars%h2osoi_vol_col, & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
         smp_l      =>    waterstate_vars%smp_l_col,      & ! Output: [real(r8) (:,:) ]  soil suction (mm)
         bsw        =>    soilstate_vars%bsw_col,         & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         watsat     =>    soilstate_vars%watsat_col,      & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         sucsat     =>    soilstate_vars%sucsat_col       & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         )

      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(.not. this%active_col(c))cycle
            if(j>=1)then
               if(t_soisno(c,j)<tfrz)then
                  smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
               else
                  s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                  ! FIXME(bja, 201603) this depends on the SWRC
                  ! implemnted by the driving LSM. Doesn't agree with
                  ! stub version...
                  ! the following call is CLM specific, jyt, Mar 13, 2016
                  !Xcall soil_water_retention_curve%soil_suction(c, j, s_node, soilstate_vars, smp_l(c,j))
               endif

            endif
         enddo
      enddo
    end associate
  end subroutine CLMCalcSmpL

!X!  !------------------------------------------------------------------------
!X!  subroutine betr_clm_readParams(this, ncid)
!X!
!X!    use ncdio_pio, only : file_desc_t
!X!
!X!    implicit none
!X!
!X!    class(betr_simulation_clm_type), intent(inout) :: this
!X!    type(file_desc_t), intent(inout) :: ncid  ! pio netCDF file id
!X!    call this%betr%bgc_reaction%readParams(ncid, this%betr%tracers)
!X!  end subroutine betr_clm_readParams

  !------------------------------------------------------------------------
  subroutine clm_h2oiso_consistency_check(this, &
       bounds, ubj, num_soilc, filter_soilc, waterstate_vars)
    !DESCRIPTION
    ! check the overall water mass consistency between betr and clm
    !
    !USES
    use MathfuncMod      , only : dot_sum
    use clm_varcon       , only : denh2o
    use WaterStateType   , only : waterstate_type
    use clm_time_manager , only : get_nstep
    implicit none
    !ARGUMENTS
    class(betr_simulation_clm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds  ! bounds
    integer                         , intent(in)    :: ubj
    integer                         , intent(in)    :: num_soilc
    integer                         , intent(in)    :: filter_soilc(:)
    type(waterstate_type)           , intent(in)    :: waterstate_vars
    !TEMPORARY VARIABLES
    real(r8), allocatable :: eyev(:)
    integer               :: fc, c
    real(r8)              :: totwater, err

    call this%bsimstatus%reset()
    ! humor the compiler about unused variables
    if (bounds%begc > 0) continue

    c = 1
    associate(                                                                    &
         h2osoi_ice           => waterstate_vars%h2osoi_ice_col,                  & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq           => waterstate_vars%h2osoi_liq_col,                  & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
         end_tracer_molarmass => this%betr(c)%tracerstates%end_tracer_molarmass_col, &
         id_trc_o18_h2o       => this%betr(c)%tracers%id_trc_o18_h2o                 &
         )


      allocate(eyev(1:ubj))
      eyev=1._r8

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         if(.not. this%active_col(c))cycle
         totwater=dot_sum(h2osoi_ice(c,1:ubj),eyev,this%bstatus(c)) + dot_sum(h2osoi_liq(c,1:ubj),eyev,this%bstatus(c))
         if(this%bstatus(c)%check_status())then
           call this%bsimstatus%setcol(c)
           call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
           exit
         endif
         err = totwater-end_tracer_molarmass(c,id_trc_o18_h2o)
         print*,get_nstep(),'diff',c, totwater, end_tracer_molarmass(c,id_trc_o18_h2o),err, err/totwater
      enddo
      deallocate(eyev)
      if(this%bsimstatus%check_status()) &
         call endrun(msg=this%bsimstatus%print_msg())
    end associate
  end subroutine clm_h2oiso_consistency_check

  !------------------------------------------------------------------------
  subroutine CLMSetBiophysForcing(this, bounds, col, pft, carbonflux_vars, waterstate_vars, &
    waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use SoilStateType     , only : soilstate_type
  use WaterStateType    , only : Waterstate_Type
  use TemperatureType   , only : temperature_type
  use ChemStateType     , only : chemstate_type
  use WaterfluxType     , only : waterflux_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CNCarbonFluxType  , only : carbonflux_type
  use CanopyStateType   , only : canopystate_type
  use clm_varpar        , only : nlevsno, nlevsoi
  use ColumnType        , only : column_type
  use PatchType         , only : patch_type
  implicit none
  !ARGUMENTS
  class(betr_simulation_clm_type) , intent(inout)        :: this
  type(bounds_type)               , intent(in)           :: bounds
  type(patch_type)                , intent(in) :: pft
  type(column_type)               , intent(in)    :: col ! column type
  type(carbonflux_type)           , optional, intent(in) :: carbonflux_vars
  type(Waterstate_Type)           , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)            , optional, intent(in) :: waterflux_vars
  type(temperature_type)          , optional, intent(in) :: temperature_vars
  type(soilhydrology_type)        , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)              , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)          , optional, intent(in) :: canopystate_vars
  type(chemstate_type)            , optional, intent(in) :: chemstate_vars
  type(soilstate_type)            , optional, intent(in) :: soilstate_vars

  call this%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, carbonflux_vars, waterstate_vars, &
      waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
      chemstate_vars, soilstate_vars)

  !the following will be CLM specific
  end subroutine CLMSetBiophysForcing
end module BeTRSimulationCLM
