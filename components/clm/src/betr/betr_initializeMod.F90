module betr_initializeMod
  !
  ! !DESCRIPTION:
  !  subroutines to initialize betr based on namelist

  ! !USES:
  use BeTRTracerType            , only : BeTRtracer_type
  use TracerCoeffType           , only : TracerCoeff_type
  use TracerFluxType            , only : TracerFlux_type
  use TracerStateType           , only : TracerState_type
  use tracerboundarycondType    , only : tracerboundarycond_type
  use BGCReactionsMod           , only : bgc_reaction_type
  use PlantSoilnutrientFluxType , only : plantsoilnutrientflux_type
  use clm_varctl                , only : iulog
  use abortutils                , only : endrun
  use shr_log_mod               , only : errMsg => shr_log_errMsg
  implicit none
  save
  private   ! By default everything is public

  public :: betr_initialize
  public :: betr_readNL
  character(len=32) :: bgc_method='mock_run'

  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  type(BeTRtracer_type)                , public :: betrtracer_vars
  type(TracerCoeff_type)               , public :: tracercoeff_vars
  type(TracerFlux_type)                , public :: tracerflux_vars
  type(TracerState_type)               , public :: tracerState_vars
  type(tracerboundarycond_type)        , public :: tracerboundarycond_vars
  type(plantsoilnutrientflux_type)     , public :: plantsoilnutrientflux_vars
  class(bgc_reaction_type),allocatable , public :: bgc_reaction

contains

  !-------------------------------------------------------------------------------
  subroutine betr_readNL(NLFilename)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename              ! Namelist filename
                                                            !
                                                            ! !LOCAL VARIABLES:
    integer                      :: ierr                    ! error code
    integer                      :: unitn                   ! unit for namelist file
    character(len=32)            :: subname = 'betr_readNL' ! subroutine name
    !-----------------------------------------------------------------------

    namelist / betr_inparm / bgc_method

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in betr_inparm  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_CanopyHydrology_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, betr_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading betr_inparm namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if
    ! Broadcast namelist variables read in
    call shr_mpi_bcast(bgc_method, mpicom)

  end subroutine betr_readNL

  !-------------------------------------------------------------------------------
  subroutine betr_initialize(bounds, lbj, ubj, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Initialize BeTR
    !
    ! !USES:
    use decompMod             , only : bounds_type
    use BGCReactionsFactoryMod, only : ctreate_bgc_reaction_type
    use BetrBGCMod            , only : betrbgc_init
    use TransportMod          , only : init_transportmod
    use TracerParamsMod       , only : tracer_param_init
    use WaterstateType        , only : waterstate_type

    implicit none
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in) :: bounds
    integer              , intent(in) :: lbj, ubj
    type(waterstate_type), intent(in) :: waterstate_vars

    call betrtracer_vars%init_scalars()

    allocate(bgc_reaction, source=ctreate_bgc_reaction_type(bgc_method))

    call bgc_reaction%Init_betrbgc(bounds, lbj, ubj, betrtracer_vars)

    call init_transportmod

    call tracerState_vars%Init(bounds, lbj, ubj, betrtracer_vars)

    call tracerflux_vars%Init(bounds,  lbj, ubj, betrtracer_vars)

    call tracercoeff_vars%Init(bounds, lbj, ubj, betrtracer_vars)

    call tracerboundarycond_vars%Init(bounds, betrtracer_vars)

    call plantsoilnutrientflux_vars%Init(bounds, lbj, ubj)

    !initialize state variable
    call bgc_reaction%initCold(bounds,  betrtracer_vars, waterstate_vars, tracerstate_vars)

    !initialize boundary condition type
    call bgc_reaction%init_boundary_condition_type(bounds, betrtracer_vars, tracerboundarycond_vars)

    !initialize the betr parameterization module
    call tracer_param_init(bounds)

    !initialize the betrBGC module
    call betrbgc_init(bounds, betrtracer_vars)

  end subroutine betr_initialize
  !---------------------------------------------------------------------------------

end module betr_initializeMod
