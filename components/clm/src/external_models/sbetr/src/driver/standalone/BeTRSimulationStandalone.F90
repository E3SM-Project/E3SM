module BeTRSimulationStandalone
  !
  ! !DESCRIPTION:
  !  BeTR standalone simulation class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures. The standalone class use BeTR data structures
  !  natively, so is mapping BeTR to BeTR. This class shouldn't be
  !  doing much.
  !
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use betr_decompMod      , only : betr_bounds_type
  use decompMod           , only : bounds_type
  use BeTRSimulation      , only : betr_simulation_type
  use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use EcophysConType      , only : ecophyscon_type
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_standalone_type
     ! NOTE(bja, 201603) LSM specific types here!

     ! NOTE(bja, 201603) most (all?) BeTR types go into the base
     ! class.

   contains
     procedure, public :: Init                => StandaloneInitOffline
     procedure, public :: StepWithoutDrainage => StandaloneStepWithoutDrainage
     procedure, public :: StepWithDrainage    => StandaloneStepWithDrainage
     procedure, public :: SetBiophysForcing   => StandaloneSetBiophysForcing
  end type betr_simulation_standalone_type

  public :: create_betr_simulation_standalone

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_standalone()
  ! DESCRIPTION
  ! constructor
    implicit none

    class(betr_simulation_standalone_type), pointer :: create_betr_simulation_standalone
    class(betr_simulation_standalone_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_standalone => simulation

  end function create_betr_simulation_standalone

  !-------------------------------------------------------------------------------

  subroutine StandaloneInitOffline(this, bounds, lun, col, pft, waterstate, namelist_buffer, base_filename)

    !DESCRIPTION
    !initialize standalone betr
    !
    !USES
    use BeTRSimulation      , only : BeTRSimulationInit
    use betr_constants      , only : betr_namelist_buffer_size
    use betr_constants      , only : betr_filename_length
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use PatchType           , only : patch_type
    use ColumnType          , only : column_type
    use LandunitType        , only : landunit_type
    use pftvarcon           , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType      , only : waterstate_type
    use CNStateType         , only : cnstate_type
    use landunit_varcon     , only : istcrop, istice, istsoil
    use BeTR_landvarconType , only : betr_landvarcon
    use clm_varpar          , only : nlevsno, nlevsoi, nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_standalone_type)   , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in)    :: col
    type(patch_type)                         , intent(in)    :: pft
    type(waterstate_type)                    , intent(inout) :: waterstate

    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    !pass necessary data for correct subroutine call

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice

    ! allocate the reaction types that may only be known to this
    ! simulation type by deallocating and overriding methods created
    ! in betr%Init().

    ! now call the base simulation init to continue initialization
    call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, base_filename)

    !pass necessary data

  end subroutine StandaloneInitOffline


  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithoutDrainage(this, bounds, col, pft)
    !DESCRIPTION
    !march one step without drainage
    !
    !USES
    use ColumnType        , only : column_type
    use PatchType         , only : patch_type
    use BeTR_TimeMod      , only : betr_time_type
    use LandunitType      , only : lun
    use pftvarcon         , only : crop
    implicit none
    !ARGUMENTS
    class(betr_simulation_standalone_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds ! bounds
    type(column_type)                      , intent(in)    :: col ! column type
    type(patch_type)                       , intent(in)    :: pft
    !temporary variables
    type(betr_bounds_type)     :: betr_bounds

    integer  :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer  :: c
    !pass necessary data for correct subroutine call
    !set lbj and ubj

    call this%BeTRSetBounds(betr_bounds)

    call this%bsimstatus%reset()

    call this%BeTRSetcps(bounds, col, pft)
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine StandaloneStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine StandaloneStepWithDrainage(this, bounds, col)
    !DESCRIPTION
    !march one step with drainage
    !
    !USES
    use ColumnType        , only : column_type
    use MathfuncMod       , only : safe_div
    use WaterFluxType     , only : waterflux_type
    use LandunitType      , only : lun
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_standalone_type) , intent(inout) :: this
    type(bounds_type)                      , intent(in)    :: bounds
    type(column_type)                      , intent(in)    :: col ! column type

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer :: c

    !pass necessary data for correct subroutine call
    !set lbj and ubj

    call this%BeTRSetBounds(betr_bounds)

    call this%bsimstatus%reset()

    call this%BeTRSetcps(bounds, col)
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%step_with_drainage(betr_bounds, this%betr_col(c), &
         this%num_soilc, this%filter_soilc,        &
         this%jtops, this%biogeo_flux(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
  end subroutine StandaloneStepWithDrainage

  !------------------------------------------------------------------------
  subroutine StandaloneSetBiophysForcing(this, bounds, col, pft,  carbonflux_vars, waterstate_vars, &
    waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars, cnstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use ColumnType        , only : column_type
  use PatchType         , only : patch_type
  use SoilStateType     , only : soilstate_type
  use WaterStateType    , only : Waterstate_Type
  use TemperatureType   , only : temperature_type
  use ChemStateType     , only : chemstate_type
  use WaterfluxType     , only : waterflux_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CNStateType       , only : cnstate_type
  use CNCarbonFluxType  , only : carbonflux_type
  use CanopyStateType   , only : canopystate_type
  use clm_varpar        , only : nlevsno, nlevsoi
  implicit none
  !ARGUMENTS
  class(betr_simulation_standalone_type) , intent(inout)        :: this
  type(bounds_type)                      , intent(in)           :: bounds
  type(column_type)                      , intent(in)    :: col ! column type
  type(patch_type)                       , intent(in) :: pft
  type(cnstate_type)                     , optional, intent(in) :: cnstate_vars
  type(carbonflux_type)                  , optional, intent(in) :: carbonflux_vars
  type(Waterstate_Type)                  , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)                   , optional, intent(in) :: waterflux_vars
  type(temperature_type)                 , optional, intent(in) :: temperature_vars
  type(soilhydrology_type)               , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)                     , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)                 , optional, intent(in) :: canopystate_vars
  type(chemstate_type)                   , optional, intent(in) :: chemstate_vars
  type(soilstate_type)                   , optional, intent(in) :: soilstate_vars

  call this%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, carbonflux_vars, waterstate_vars, &
      waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
      chemstate_vars, soilstate_vars)

  !the following will be standalone bgc specific
  end subroutine StandaloneSetBiophysForcing
end module BeTRSimulationStandalone
