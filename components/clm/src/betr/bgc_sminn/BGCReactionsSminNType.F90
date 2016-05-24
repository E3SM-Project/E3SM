module BGCReactionsSminNType


#include "shr_assert.h"
  !
  ! !DESCRIPTION:
  ! do NO3 leaching using betr and the default century bgc
  !
  ! HISTORY:
  ! Created by Jinyun Tang, Oct 2nd, 2014
  ! !USES:
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use BGCReactionsMod , only : bgc_reaction_type
  use tracer_varcon   , only : bndcond_as_conc, bndcond_as_flux
  use LandunitType    , only : lun
  use ColumnType      , only : col
  implicit none

  save
  private
  !
  ! !PUBLIC TYPES:
  public :: bgc_reaction_sminn_type

  type, extends(bgc_reaction_type) :: &
       bgc_reaction_sminn_type
  private
contains
  procedure :: Init_betrbgc                    ! initialize betr bgc
  procedure :: set_boundary_conditions         ! set top/bottom boundary conditions for various tracers
  procedure :: calc_bgc_reaction               ! doing bgc calculation
  procedure :: init_boundary_condition_type    ! initialize type of top boundary conditions
  procedure :: do_tracer_equilibration         ! do equilibrium tracer chemistry
  procedure :: InitCold                        ! do cold initialization
  procedure :: readParams                      ! read in parameters
  procedure :: betr_alm_flux_statevar_feedback !
  procedure :: init_betr_alm_bgc_coupler
end type bgc_reaction_sminn_type

     interface bgc_reaction_sminn_type
module procedure constructor

end interface bgc_reaction_sminn_type

contains
  !-------------------------------------------------------------------------------
  type(bgc_reaction_sminn_type) function constructor()
    !
    ! ! DESCRIPTION:
    ! create an object of type bgc_reaction_sminn_type.
    ! Right now it is purposely empty

  end function constructor

  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! !DESCRIPTION:
    ! initialize boundary condition types
    ! !USES:
    use BeTRTracerType        , only : betrtracer_type
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use BeTRTracerType        , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_sminn_type), intent(in) :: this
    type(BeTRtracer_type ),            intent(in) :: betrtracer_vars
    type(bounds_type),                 intent(in) :: bounds
    type(tracerboundarycond_type),     intent(in) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:
    integer,  parameter :: bndcond_as_conc = 1   !top boundary conditions as tracer concentration
    integer,  parameter :: bndcond_as_flx=2      !top boundary condition as tracer flux
    integer :: c

    tracerboundarycond_vars%topbc_type(:) = bndcond_as_flx
  end subroutine init_boundary_condition_type

  !-------------------------------------------------------------------------------
  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! initialize the betrbgc
    ! !USES:
    use BeTRTracerType        , only : betrtracer_type
    use MathfuncMod           , only : addone

    ! !ARGUMENTS:
    class(bgc_reaction_sminn_type), intent(in) :: this
    type(bounds_type),                 intent(in) :: bounds
    integer,                           intent(in) :: lbj, ubj
    type(BeTRtracer_type ),         intent(inout) :: betrtracer_vars

    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname ='Init_betrbgc'
    integer :: itemp_gwm
    integer :: itemp_gwm_grp
    integer :: item_grp
    integer :: dum
    integer :: itemp_grp

    itemp_gwm = 0; itemp_gwm_grp = 0

    betrtracer_vars%id_trc_nh3x = addone(itemp_gwm); dum = addone(itemp_gwm_grp)
    betrtracer_vars%id_trc_no3x = addone(itemp_gwm); dum = addone(itemp_gwm_grp)

    betrtracer_vars%ngwmobile_tracers      = itemp_gwm;   betrtracer_vars%ngwmobile_tracer_groups= itemp_gwm_grp
    betrtracer_vars%nmem_max               = 1

    call betrtracer_vars%Init()
    itemp_grp = 0    !group id

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_nh3x, trc_name='NH3x',    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = addone(itemp_grp), &
         trc_group_mem = 1, is_trc_volatile=.false.)

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_no3x, trc_name='NO3x',    &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.false.,trc_vtrans_scal=0._r8)

  end subroutine Init_betrbgc

  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
       waterflux_vars, tracerboundarycond_vars)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    ! !USES:
    use clm_varctl            , only : iulog
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use abortutils            , only : endrun
    use shr_log_mod           , only : errMsg => shr_log_errMsg
    use BeTRTracerType        , only : betrtracer_type
    use WaterfluxType         , only : waterflux_type

    ! !ARGUMENTS:
    class(bgc_reaction_sminn_type), intent(in)    :: this
    type(bounds_type),              intent(in)    :: bounds
    integer,                        intent(in)    :: num_soilc               ! number of columns in column filter_soilc
    integer,                        intent(in)    :: filter_soilc(:)         ! column filter_soilc
    type(betrtracer_type),          intent(in)    :: betrtracer_vars
    real(r8),                       intent(in)    :: dz_top(bounds%begc: )
    type(waterflux_type),           intent(in)    :: waterflux_vars
    type(tracerboundarycond_type),  intent(inout) :: tracerboundarycond_vars !

    ! !LOCAL VARIABLES:
    integer                                       :: fc, c
    character(len=255)                            :: subname = 'set_boundary_conditions'

    SHR_ASSERT_ALL((ubound(dz_top)                == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))

    associate(                                     &
         groupid  => betrtracer_vars%groupid       &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)

         !eventually, the following code will be implemented using polymorphism
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_no3x) = 0._r8
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_nh3x) = 0._r8

         tracerboundarycond_vars%bot_concflux_col(c,1,:) = 0._r8                                  !zero flux boundary condition

      enddo
    end associate
  end subroutine set_boundary_conditions

  !-------------------------------------------------------------------------------
  subroutine calc_bgc_reaction(this, bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, jtops,     &
       dtime, betrtracer_vars, tracercoeff_vars, waterstate_vars, temperature_vars, soilstate_vars, chemstate_vars, &
       cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, tracerstate_vars,    &
       tracerflux_vars, plantsoilnutrientflux_vars)
    !
    ! !DESCRIPTION:
    ! do bgc reaction
    ! !USES:
    !
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use tracercoeffType          , only : tracercoeff_type
    use BetrTracerType           , only : betrtracer_type
    use WaterStateType           , only : Waterstate_Type
    use TemperatureType          , only : temperature_type
    use SoilStatetype            , only : soilstate_type
    use ChemStateType            , only : chemstate_type
    use CanopyStateType          , only : canopystate_type
    use CNStateType              , only : cnstate_type
    use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
    use CNCarbonFluxType         , only : carbonflux_type
    use CNCarbonStateType        , only : carbonstate_type
    use CNNitrogenFluxType       , only : nitrogenflux_type
    use CNNitrogenStateType      , only : nitrogenstate_type
    use clm_varpar               , only : nlevtrc_soil
    use clm_varcon               , only : natomw
    !ARGUMENTS
    class(bgc_reaction_sminn_type)   , intent(in)    :: this
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    integer                          , intent(in)    :: num_soilp
    integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
    integer                          , intent(in)    :: jtops(bounds%begc: )       ! top index of each column
    integer                          , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0
    real(r8)                         , intent(in)    :: dtime                      ! model time step
    type(Waterstate_Type)            , intent(in)    :: waterstate_vars            ! water state variables
    type(temperature_type)           , intent(in)    :: temperature_vars           ! energy state variable
    type(soilstate_type)             , intent(in)    :: soilstate_vars
    type(cnstate_type)               , intent(inout) :: cnstate_vars
    type(carbonstate_type)           , intent(in)    :: carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(chemstate_type)             , intent(in)    :: chemstate_vars
    type(betrtracer_type)            , intent(in)    :: betrtracer_vars            ! betr configuration information
    type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(plantsoilnutrientflux_type) , intent(inout) :: plantsoilnutrientflux_vars !
    character(len=*), parameter :: subname ='calc_bgc_reaction'

    !do nothing here
  end subroutine calc_bgc_reaction

  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, &
       filter_soilc, betrtracer_vars, tracercoeff_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! requilibrate tracers that has solid and mobile phases
    ! using the theory of mass action.

    ! !USES:
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use BeTRTracerType        , only : betrtracer_type
    use BeTRTracerType        , only : betrtracer_type
    class(bgc_reaction_sminn_type),    intent(in) :: this

    ! !ARGUMENTS:
    type(bounds_type),      intent(in) :: bounds
    integer,                intent(in) :: lbj, ubj
    integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
    integer,                intent(in) :: num_soilc
    integer,                intent(in) :: filter_soilc(:)
    type(betrtracer_type),  intent(in) :: betrtracer_vars
    type(tracercoeff_type), intent(in) :: tracercoeff_vars
    type(tracerstate_type), intent(inout) :: tracerstate_vars
    character(len=255) :: subname = 'do_tracer_equilibration'


    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))

    !depending on the simulation type, an implementation of aqueous chemistry will be
    !employed to separate out the adsorbed phase
    !It should be noted that this formulation excludes the use of linear isotherm, which
    !can be integrated through the retardation factor


  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, betrtracer_vars, waterstate_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    !  do cold initialization
    !
    ! !USES:
    use BeTRTracerType           , only : BeTRTracer_Type
    use tracerstatetype          , only : tracerstate_type
    use WaterstateType           , only : waterstate_type
    use PatchType                , only : pft
    use clm_varcon               , only : spval, ispval
    use landunit_varcon          , only : istsoil, istcrop

    ! !ARGUMENTS:
    class(bgc_reaction_sminn_type) , intent(in)    :: this
    type(bounds_type)                 , intent(in)    :: bounds
    type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
    type(waterstate_type)             , intent(in)    :: waterstate_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter_soilc index
    integer               :: begc, endc
    integer               :: begg, endg
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------


    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          if(betrtracer_vars%ngwmobile_tracers>0)then
             tracerstate_vars%tracer_conc_mobile_col(c,:,:)        = spval
             tracerstate_vars%tracer_conc_surfwater_col(c,:)       = spval
             tracerstate_vars%tracer_conc_aquifer_col(c,:)         = spval
             tracerstate_vars%tracer_conc_grndwater_col(c,:)       = spval
          endif
          if(betrtracer_vars%ntracers > betrtracer_vars%ngwmobile_tracers)then
             tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = spval
          endif
          if(betrtracer_vars%nsolid_equil_tracers>0)then
             tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = spval
          endif
       endif
       tracerstate_vars%tracer_soi_molarmass_col(c,:)            = spval

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          !dual phase tracers

          tracerstate_vars%tracer_conc_mobile_col(c,:, :)          = 0._r8
          tracerstate_vars%tracer_conc_surfwater_col(c,:)          = 0._r8
          tracerstate_vars%tracer_conc_aquifer_col(c,:)            = 0._r8
          tracerstate_vars%tracer_conc_grndwater_col(c,:)          = 0._r8


          !solid tracers
          if(betrtracer_vars%ngwmobile_tracers < betrtracer_vars%ntracers)then
             tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = 0._r8
          endif

          if(betrtracer_vars%nsolid_equil_tracers>0)then
             tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
          endif
          tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8
       endif
    enddo

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine readParams(this, ncid, betrtracer_vars)

    ! !DESCRIPTION:
    ! read in parameters
    ! !USES:
    use ncdio_pio                , only : file_desc_t
    use BeTRTracerType           , only : BeTRTracer_Type

    ! !ARGUMENTS:
    class(bgc_reaction_sminn_type)    , intent(in)    :: this
    type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars
    type(file_desc_t)                 , intent(inout) :: ncid  ! pio netCDF file id

    !do nothing here
  end subroutine readParams

  !-------------------------------------------------------------------------------
  subroutine betr_alm_flux_statevar_feedback(this, bounds, num_soilc, filter_soilc, &
       carbonstate_vars, nitrogenstate_vars, nitrogenflux_vars, tracerstate_vars,   &
       tracerflux_vars,  betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! do state variable and flux variable exchange between betr and alm
    ! !USES:
    use shr_kind_mod               , only          : r8 => shr_kind_r8
    use tracerfluxType             , only          : tracerflux_type
    use tracerstatetype            , only          : tracerstate_type
    use decompMod                  , only          : bounds_type
    use BeTRTracerType             , only          : BeTRTracer_Type
    use CNCarbonStateType          , only          : carbonstate_type
    use CNNitrogenStateType        , only          : nitrogenstate_type
    use CNNitrogenFluxType         , only          : nitrogenflux_type

    ! !ARGUMENTS                                   :
    class(bgc_reaction_sminn_type) , intent(in)    :: this
    type(bounds_type)              , intent(in)    :: bounds             ! bounds
    integer                        , intent(in)    :: num_soilc          ! number of columns in column filter
    integer                        , intent(in)    :: filter_soilc(:)    ! column filter
    type(betrtracer_type)          , intent(in)    :: betrtracer_vars    ! betr configuration information
    type(tracerstate_type)         , intent(in)    :: tracerstate_vars
    type(tracerflux_type)          , intent(in)    :: tracerflux_vars
    type(carbonstate_type)         , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)        , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)       , intent(inout) :: nitrogenstate_vars !

    call assign_sminn_pools(bounds, num_soilc, filter_soilc,  carbonstate_vars, &
         nitrogenstate_vars, tracerstate_vars, betrtracer_vars)

    call assign_minnitrogen_hydroloss(bounds, num_soilc, filter_soilc, &
         tracerflux_vars, nitrogenflux_vars, betrtracer_vars)

  end subroutine betr_alm_flux_statevar_feedback

  !-------------------------------------------------------------------------------
  
  subroutine init_betr_alm_bgc_coupler(this, bounds, carbonstate_vars, &
       nitrogenstate_vars, betrtracer_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! initialize state and flux variable exchange between betr and alm
    ! !USES:
    use clm_varcon                 , only          : natomw, catomw
    use clm_varpar                 , only          : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    use CNCarbonStateType          , only          : carbonstate_type
    use CNNitrogenStateType        , only          : nitrogenstate_type
    use tracerstatetype            , only          : tracerstate_type
    use BetrTracerType             , only          : betrtracer_type
    use clm_varpar                 , only          : nlevtrc_soil
    use landunit_varcon            , only          : istsoil, istcrop

    ! !ARGUMENTS :
    class(bgc_reaction_sminn_type) , intent(in)    :: this
    type(bounds_type)              , intent(in)    :: bounds
    type(tracerstate_type)         , intent(inout) :: tracerstate_vars
    type(betrtracer_type)          , intent(in)    :: betrtracer_vars                    ! betr configuration information
    type(carbonstate_type)         , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type)       , intent(in)    :: nitrogenstate_vars

    ! !LOCAL VARIABLES:
    integer :: j, fc, c, l
    associate(                                                          &
         id_trc_no3x        => betrtracer_vars%id_trc_no3x            , &
         id_trc_nh3x        => betrtracer_vars%id_trc_nh3x            , &
         smin_no3_vr_col    => nitrogenstate_vars%smin_no3_vr_col     , &
         smin_nh4_vr_col    => nitrogenstate_vars%smin_nh4_vr_col     , &
         tracer_conc_mobile => tracerstate_vars%tracer_conc_mobile_col  &
         )

      do j = 1, nlevtrc_soil
         do c = bounds%begc, bounds%endc
            l = col%landunit(c)
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

               tracer_conc_mobile(c,j,id_trc_no3x) = smin_no3_vr_col(c,j) / natomw
               tracer_conc_mobile(c,j,id_trc_nh3x) = smin_nh4_vr_col(c,j) /natomw
            endif
         enddo
      enddo

    end associate
  end subroutine init_betr_alm_bgc_coupler
  !-------------------------------------------------------------------------------
  subroutine assign_minnitrogen_hydroloss(bounds, num_soilc, filter_soilc, tracerflux_vars, &
       nitrogenflux_vars, betrtracer_vars)

    !
    ! !DESCRIPTION:
    ! feedback the nitrogen hydrological fluxes, this comes after tracer mass balance, so the flux is with the unit of st/m2/s
    ! !USES:
    use tracerfluxType      , only : tracerflux_type
    use BetrTracerType      , only : betrtracer_type
    use CNNitrogenFluxType  , only : nitrogenflux_type
    use clm_varcon          , only : natomw
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds          ! bounds
    integer                 , intent(in)    :: num_soilc       ! number of columns in column filter
    integer                 , intent(in)    :: filter_soilc(:) ! column filter
    type(tracerflux_type)   , intent(in)    :: tracerflux_vars
    type(betrtracer_type)   , intent(in)    :: betrtracer_vars ! betr configuration information
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars

    ! !LOCAL VARIABLES:
    integer :: fc, c
    !get nitrogen leaching, and loss through surface runoff

    associate(                                      &
         id_trc_no3x => betrtracer_vars%id_trc_no3x &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         nitrogenflux_vars%smin_no3_leached_col(c) = tracerflux_vars%tracer_flx_totleached_col(c,id_trc_no3x)*natomw
         nitrogenflux_vars%smin_no3_runoff_col(c)  = tracerflux_vars%tracer_flx_surfrun_col(c,id_trc_no3x)*natomw
      enddo

    end associate
  end subroutine assign_minnitrogen_hydroloss

  !-------------------------------------------------------------------------------
  subroutine assign_sminn_pools(bounds, num_soilc, filter_soilc, carbonstate_vars, &
       nitrogenstate_vars, tracerstate_vars, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! update mineral nitrogen pool
    ! !USES:
    use clm_varcon           , only : natomw, catomw
    use clm_varpar           , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    use CNCarbonStateType    , only : carbonstate_type
    use CNNitrogenStateType  , only : nitrogenstate_type
    use tracerstatetype      , only : tracerstate_type
    use BetrTracerType       , only : betrtracer_type
    use clm_varpar           , only : nlevtrc_soil
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds             ! bounds
    integer                  , intent(in)    :: num_soilc          ! number of columns in column filter
    integer                  , intent(in)    :: filter_soilc(:)    ! column filter
    type(tracerstate_type)   , intent(in)    :: tracerstate_vars
    type(betrtracer_type)    , intent(in)    :: betrtracer_vars    ! betr configuration information
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars !

    ! !LOCAL VARIABLES:
    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7

    integer :: fc, c, j, k
    associate(                                                          &
         id_trc_no3x        => betrtracer_vars%id_trc_no3x            , &
         id_trc_nh3x        => betrtracer_vars%id_trc_nh3x            , &
         smin_no3_vr_col    => nitrogenstate_vars%smin_no3_vr_col     , &
         smin_nh4_vr_col    => nitrogenstate_vars%smin_nh4_vr_col     , &
         sminn_vr_col       => nitrogenstate_vars%sminn_vr_col        , &
         tracer_conc_mobile => tracerstate_vars%tracer_conc_mobile_col  &
         )

      do j = 1, nlevtrc_soil
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            smin_no3_vr_col(c,j) = tracer_conc_mobile(c,j,id_trc_no3x)*natomw
            smin_nh4_vr_col(c,j) = tracer_conc_mobile(c,j,id_trc_nh3x)*natomw
            sminn_vr_col   (c,j) = smin_no3_vr_col(c,j) + smin_nh4_vr_col(c,j)

         enddo
      enddo

    end associate
  end subroutine assign_sminn_pools

end module BGCReactionsSminNType
