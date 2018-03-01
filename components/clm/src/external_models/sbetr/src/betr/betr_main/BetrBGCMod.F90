module BetrBGCMod

#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  !  subroutines for betr application
  !
  !  !USES:
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
  use betr_ctrl                , only : iulog  => biulog
  use betr_constants           , only : betr_string_length
  use betr_varcon              , only : denh2o => bdenh2o
  use MathfuncMod              , only : dot_sum
  use BGCReactionsMod          , only : bgc_reaction_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use BeTRTracerType           , only : betrtracer_type
  use TracerCoeffType          , only : TracerCoeff_type
  use TracerFluxType           , only : TracerFlux_type
  use TracerStateType          , only : TracerState_type
  use tracerboundarycondType   , only : tracerboundarycond_type
  use BeTR_aerocondType        , only : betr_aerecond_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use betr_biogeoStateType     , only : betr_biogeo_state_type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use tracer_varcon            , only : nlevsno => betr_nlevsno
  implicit none

  private

  integer,  parameter :: diffusion_scheme = 1         ! do diffusive transport
  integer,  parameter :: advection_scheme = 2         ! do advective transport, aquesou phase only

  real(r8), parameter :: tiny_val         = 1.e-20_r8 !very small value, for tracer concentration etc.
  real(r8), parameter :: dtime_min        = 1._r8     !minimum time step 1 second
  real(r8), parameter :: err_tol_transp   = 1.e-8_r8  !error tolerance for tracer transport
  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: stage_tracer_transport
  public :: surface_tracer_hydropath_update
  public :: tracer_gws_transport
  public :: calc_ebullition
  public :: diagnose_gas_pressure
  public :: tracer_copy_a2b_div
  public :: tracer_copy_a2b_comb
  public :: tracer_col_mapping_div
  public :: tracer_col_mapping_comb
contains


  subroutine surface_tracer_hydropath_update(betr_time, bounds, col, &
     num_soilc, filter_soilc, &
     biophysforc, betrtracer_vars, tracerstate_vars, &
     tracercoeff_vars,  tracerflux_vars, betr_status)
    !
    ! !DESCRIPTIONS
    ! do tracer update through, residual snow, and surface runoff
    ! !USES:
    use tracerfluxType  , only : tracerflux_type
    use tracerstatetype , only : tracerstate_type
    use tracercoeffType , only : tracercoeff_type
    use BetrTracerType  , only : betrtracer_type
    use BeTR_TimeMod    , only : betr_time_type
    use BetrStatusType  , only : betr_status_type
    use betr_columnType          , only : betr_column_type
    !
    ! Arguments
    implicit none
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    type(betr_column_type)           , intent(in)    :: col
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars            ! betr configuration information
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(betr_status_type)           , intent(out)   :: betr_status

    integer :: lbj, ubj

    lbj = bounds%lbj; ubj = bounds%ubj
    !calculate flux from merging topsoil with surface ponding water and snow
    call calc_tracer_h2osfc_snow_residual_combine(betr_time, &
       bounds, col, num_soilc, filter_soilc,                 &
       biophysforc,                                          &
       betrtracer_vars,                                      &
       tracerstate_vars,                                     &
       tracerflux_vars, betr_status)
    if(betr_status%check_status())return

    ! do tracer wash with surface runoff
    call calc_tracer_surface_runoff(betr_time, bounds, lbj, ubj, &
       num_soilc,                                                &
       filter_soilc,                                             &
       biophysforc%fracice_col(bounds%begc:bounds%endc,1),       &
       col%dz(bounds%begc:bounds%endc, 1:ubj),                   &
       biophysforc,                                              &
       betrtracer_vars,                                          &
       tracerstate_vars,                                         &
       tracercoeff_vars,                                         &
       tracerflux_vars, betr_status)

  end subroutine surface_tracer_hydropath_update
  !-------------------------------------------------------------------------------

  subroutine stage_tracer_transport(betr_time, bounds, col, pft, num_soilc, &
       filter_soilc, num_soilp, filter_soilp, &
       biophysforc, biogeo_state, biogeo_flux, aerecond_vars,    &
       betrtracer_vars, tracercoeff_vars,                        &
       tracerboundarycond_vars, tracerstate_vars, tracerflux_vars, bgc_reaction, Rfactor, &
       advection_on, betr_status)
    !DESCRIPTION
    !set relevant coefficients for tracer transport
    ! !USES:
    use betr_ctrl              , only : betr_use_cn
    use tracer_varcon          , only : nlevsoi  => betr_nlevsoi
    use tracerfluxType         , only : tracerflux_type
    use tracerstatetype        , only : tracerstate_type
    use tracercoeffType        , only : tracercoeff_type
    use TracerBoundaryCondType , only : TracerBoundaryCond_type
    use BetrTracerType         , only : betrtracer_type
    use TracerParamsMod        , only : set_phase_convert_coeff, set_multi_phase_diffusion, calc_tracer_infiltration
    use TracerParamsMod        , only : get_zwt, calc_aerecond, betr_annualupdate
    use betr_columnType        , only : betr_column_type
    use BeTR_aerocondType      , only : betr_aerecond_type
    use BGCReactionsMod        , only : bgc_reaction_type
    use BeTR_TimeMod           , only : betr_time_type
    use BetrStatusType         , only : betr_status_type
    use BeTR_PatchType         , only : betr_patch_type
    implicit none
    !arguments
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    type(betr_column_type)           , intent(in)    :: col
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    integer                          , intent(in)    :: num_soilp
    integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_biogeo_state_type)     , intent(inout) :: biogeo_state
    type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
    class(betrtracer_type)           , intent(inout) :: betrtracer_vars            ! betr configuration information
    class(bgc_reaction_type)         , intent(inout) :: bgc_reaction
    type(betr_aerecond_type)         , intent(inout) :: aerecond_vars
    type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(TracerState_type)           , intent(in)    :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    real(r8)                         , intent(out)   :: Rfactor(bounds%begc: , bounds%lbj: ,1: ) !retardation factor
    logical                          , intent(in)    :: advection_on
    type(betr_status_type)           , intent(out)   :: betr_status

    !temporary variables
    integer :: jwt(bounds%begc:bounds%endc)
    integer :: lbj, ubj

    call betr_status%reset()

    SHR_ASSERT_ALL((ubound(Rfactor,1) == bounds%endc), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,2) == bounds%ubj ), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,3) == betrtracer_vars%ngwmobile_tracer_groups), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return

    lbj = bounds%lbj; ubj = bounds%ubj

    !set up jtops
    tracerboundarycond_vars%jtops_col(:)=1

    if(betr_use_cn)then
       !update npp for aerenchyma calculation
       call betr_annualupdate(betr_time, bounds, pft, num_soilc, filter_soilc, &
            num_soilp, filter_soilp, biophysforc, aerecond_vars, &
            tracercoeff_vars, betr_status)
       if(betr_status%check_status())return
    endif

    !obtain water table depth
    call get_zwt (bounds, num_soilc, filter_soilc,       &
         col%zi(bounds%begc:bounds%endc, 0:nlevsoi),     &
         biophysforc,                                    &
         biogeo_state%zwts_col(bounds%begc:bounds%endc), &
         jwt(bounds%begc:bounds%endc), betr_status)
    if(betr_status%check_status())return

    !calculate arenchyma conductance
    call calc_aerecond(bounds, col, pft, num_soilp, filter_soilp, &
         jwt(bounds%begc:bounds%endc),                  &
         biophysforc,                                   &
         betrtracer_vars,                               &
         aerecond_vars,                                 &
         tracercoeff_vars, betr_status)
    if(betr_status%check_status())return

    biophysforc%soil_pH(bounds%begc:bounds%endc,1:ubj)=7._r8

    call set_phase_convert_coeff(bounds, lbj, ubj, &
         tracerboundarycond_vars%jtops_col       , &
         num_soilc                               , &
         filter_soilc                            , &
         col%dz(bounds%begc:bounds%endc, lbj:ubj), &
         biophysforc = biophysforc               , &
         betrtracer_vars=betrtracer_vars         , &
         tracerstate_vars=tracerstate_vars       , &
         tracercoeff_vars=tracercoeff_vars       , &
         betr_status = betr_status)
    if(betr_status%check_status())return

    call set_multi_phase_diffusion(bounds, col, lbj, ubj, &
         tracerboundarycond_vars%jtops_col         , &
         num_soilc                                 , &
         filter_soilc                              , &
         biophysforc = biophysforc                 , &
         betrtracer_vars=betrtracer_vars           , &
         tracercoeff_vars=tracercoeff_vars         , &
         betr_status=betr_status)
    if(betr_status%check_status())return

    call bgc_reaction%set_boundary_conditions(bounds, num_soilc, filter_soilc, &
         col%dz(bounds%begc:bounds%endc,1),                                    &
         betrtracer_vars,                                                      &
         biophysforc                ,                                          &
         biogeo_flux                ,                                          &
         tracerboundarycond_vars, betr_status)
    if(betr_status%check_status())return

    if(advection_on) &
    call calc_tracer_infiltration(bounds,                        &
         tracerboundarycond_vars%jtops_col,                      &
         num_soilc,                                              &
         filter_soilc,                                           &
         tracercoeff_vars%bunsencef_col(bounds%begc:bounds%endc, &
         1,                                                      &
         1:betrtracer_vars%nvolatile_tracer_groups),             &
         betrtracer_vars,                                        &
         tracerboundarycond_vars,                                &
         biogeo_flux  ,                                          &
         tracerflux_vars%tracer_flx_infl_col, betr_status)
    if(betr_status%check_status())return

    call set_gwdif_Rfactor(bounds, lbj, ubj, &
         tracerboundarycond_vars%jtops_col,  &
         num_soilc,                          &
         filter_soilc,                       &
         tracercoeff_vars,                   &
         betrtracer_vars,                    &
         Rfactor, betr_status)

  end subroutine stage_tracer_transport

  !-------------------------------------------------------------------------------

  subroutine tracer_gws_transport(betr_time, bounds, col, pft, num_soilc, &
    filter_soilc, Rfactor, biophysforc, biogeo_flux, &
    betrtracer_vars, tracerboundarycond_vars, tracercoeff_vars, &
    tracerstate_vars, tracerflux_vars, bgc_reaction, advection_on, &
    diffusion_on, betr_status)
    !DESCRIPTION
    ! do multiphase transprot, diffusion, advection and turbation
    ! !USES:
    use tracerfluxType         , only : tracerflux_type
    use tracerstatetype        , only : tracerstate_type
    use tracercoeffType        , only : tracercoeff_type
    use TracerBoundaryCondType , only : TracerBoundaryCond_type
    use BetrTracerType         , only : betrtracer_type
    use BGCReactionsMod        , only : bgc_reaction_type
    use BeTR_TimeMod           , only : betr_time_type
    use BetrStatusType         , only : betr_status_type
    use betr_columnType        , only : betr_column_type
    use BeTR_PatchType         , only : betr_patch_type
    implicit none
    !ARGUMENTS
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    type(betr_column_type)           , intent(in)    :: col
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    real(r8)                         , intent(in)    :: Rfactor(bounds%begc: , bounds%lbj: ,1: ) !retardation factor
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
    class(betrtracer_type)           , intent(inout) :: betrtracer_vars            ! betr configuration information
    class(bgc_reaction_type)         , intent(inout) :: bgc_reaction
    logical                          , intent(in)    :: advection_on
    logical                          , intent(in)    :: diffusion_on
    type(tracerboundarycond_type)    , intent(inout) :: tracerboundarycond_vars
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(betr_status_type)           , intent(out)   :: betr_status
    integer :: lbj, ubj
    real(r8):: dtime

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(Rfactor,1) == bounds%endc), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,2) == bounds%ubj ), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,3) == betrtracer_vars%ngwmobile_tracer_groups), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return

    lbj = bounds%lbj; ubj = bounds%ubj
    dtime = betr_time%get_step_size()

    !do gas+aqueous diffusion and advection
    call tracer_gw_transport(betr_time, bounds, lbj, ubj, &
         col, pft,                                        &
         tracerboundarycond_vars%jtops_col,               &
         num_soilc,                                       &
         filter_soilc,                                    &
         Rfactor,                                         &
         col%dz(bounds%begc:bounds%endc, lbj:ubj),        &
         col%zi(bounds%begc:bounds%endc,lbj-1:ubj),       &
         (/advection_scheme, diffusion_scheme/),          &
         advection_on, diffusion_on,                      &
         betrtracer_vars,                                 &
         tracerboundarycond_vars,                         &
         tracercoeff_vars,                                &
         biophysforc,                                     &
         biogeo_flux,                                     &
         bgc_reaction,                                    &
         tracerstate_vars,                                &
         tracerflux_vars, betr_status)

    if(betr_status%check_status())return
    !do solid phase turbation
    call tracer_solid_transport(betr_time, bounds, 1, ubj,                         &
         num_soilc,                                                                &
         filter_soilc,                                                             &
         tracercoeff_vars%hmconductance_col(bounds%begc:bounds%endc, 1:ubj-1, : ), &
         col%dz(bounds%begc:bounds%endc, 1:ubj),                                   &
         betrtracer_vars,                                                          &
         tracerboundarycond_vars,                                                  &
         tracerflux_vars,                                                          &
         tracerstate_vars, betr_status)

  end subroutine tracer_gws_transport
  !-------------------------------------------------------------------------------
  subroutine tracer_solid_transport(betr_time,                           &
       bounds, lbj, ubj, num_soilc, filter_soilc, hmconductance_col, dz, &
       betrtracer_vars, tracerboundarycond_vars, tracerflux_vars,        &
       tracerstate_vars, betr_status)
    !
    ! !DESCRIPTION:
    !
    ! do solid phase tracer transport, due to various turbation processes,
    ! which are parameterized as diffusion
    ! the surface flux of solid tracer is zero
    !
    ! !USES:
    use tracerstateType        , only : tracerstate_type
    use tracerfluxType         , only : tracerflux_type
    use tracerboundarycondtype , only : tracerboundarycond_type
    use TransportMod           , only : DiffusTransp
    use BeTR_TimeMod           , only : betr_time_type
    use BetrStatusType         , only : betr_status_type
    use betr_constants         , only : betr_errmsg_len
    implicit none

    ! !ARGUMENTS:
    class(betr_time_type)         , intent(in)    :: betr_time
    type(bounds_type)             , intent(in)    :: bounds
    integer                       , intent(in)    :: lbj, ubj
    integer                       , intent(in)    :: num_soilc                                   ! number of columns in column filter_soilc
    integer                       , intent(in)    :: filter_soilc(:)                             ! column filter_soilc
    class(betrtracer_type)        , intent(inout) :: betrtracer_vars
    real(r8)                      , intent(in)    :: hmconductance_col(bounds%begc: , lbj: ,1: ) !weighted bulk conductance
    real(r8)                      , intent(in)    :: dz(bounds%begc: , lbj: )
    type(tracerboundarycond_type) , intent(in)    :: tracerboundarycond_vars
    type(tracerflux_type)         , intent(in)    :: tracerflux_vars
    type(tracerstate_type)        , intent(inout) :: tracerstate_vars
    type(betr_status_type)        , intent(out)   :: betr_status
    ! !LOCAL VARIABLES:
    integer              ::  kk, j, fc, c, l, ntrcs, trcid, k
    real(r8)             :: dtime_loc(bounds%begc:bounds%endc)
    real(r8)             :: time_remain(bounds%begc:bounds%endc)
    integer              :: jtops(bounds%begc:bounds%endc)
    logical              :: update_col(bounds%begc:bounds%endc)
    logical              :: lnegative_tracer
    logical              :: lexit_loop
    integer, allocatable :: difs_trc_group(:)
    real(r8), pointer    :: dtracer(:, :, :)
    real(r8), pointer    :: err_tracer(:, :)
    real(r8), pointer    :: local_source(:, :, :)
    real(r8), parameter  :: err_min_solid=1.e-12_r8
    character(len=255)   :: subname = 'tracer_solid_transport'
    integer              :: ntracer_groups
    character(len=betr_errmsg_len) :: msg

    call betr_status%reset()
    ! remove compiler warnings for unused dummy args
    if (size(tracerboundarycond_vars%jtops_col) > 0) continue

    ntracer_groups = betrtracer_vars%ntracer_groups
    if(.not. betrtracer_vars%is_solidtransport())return
    associate(                                                                                  &
         tracernames                   =>  betrtracer_vars%tracernames                        , &
         nmem_max                      =>  betrtracer_vars%nmem_max                           , &
         ngwmobile_tracers             =>  betrtracer_vars%ngwmobile_tracers                  , &
         tracer_group_memid            =>  betrtracer_vars%tracer_group_memid                 , &
         is_mobile                     =>  betrtracer_vars%is_mobile                          , &
         tracer_flx_netpro_vr          =>  tracerflux_vars%tracer_flx_netpro_vr_col           , & !
         tracer_conc_mobile_col        =>  tracerstate_vars%tracer_conc_mobile_col              &
         )

      SHR_ASSERT_ALL((ubound(hmconductance_col,1) == bounds%endc), errMsg(mod_filename,__LINE__),betr_status)
      if(betr_status%check_status())return
      SHR_ASSERT_ALL((ubound(hmconductance_col,2) == ubj-1), errMsg(mod_filename,__LINE__),betr_status)
      if(betr_status%check_status())return
      SHR_ASSERT_ALL((ubound(hmconductance_col,3) == ntracer_groups), errMsg(mod_filename,__LINE__),betr_status)
      if(betr_status%check_status())return
      SHR_ASSERT_ALL((ubound(dz,1)   == bounds%endc), errMsg(mod_filename,__LINE__),betr_status)
      if(betr_status%check_status())return
      SHR_ASSERT_ALL((ubound(dz,2)   == ubj), errMsg(mod_filename,__LINE__),betr_status)
      if(betr_status%check_status())return

      allocate (difs_trc_group (nmem_max                                   ))
      allocate (dtracer        (bounds%begc:bounds%endc, lbj:ubj, nmem_max ))
      allocate (err_tracer     (bounds%begc:bounds%endc, nmem_max          ))
      allocate (local_source   (bounds%begc:bounds%endc,lbj:ubj,nmem_max   ))

      jtops(:)            =1
      local_source(:,:,:) = 0._r8

      do j = betrtracer_vars%ngwmobile_tracer_groups+1, betrtracer_vars%ntracer_groups
         ntrcs = 0
         difs_trc_group(:) = 0
         do k = 1, nmem_max
            trcid = tracer_group_memid(j, k)
            if(trcid>0)then
               if(is_mobile(tracer_group_memid(j, k)))then
                  ntrcs = ntrcs + 1
                  difs_trc_group(ntrcs) = trcid
               endif
            endif
         enddo

         if (ntrcs==0) cycle

         !adaptive time stepping for solid phase transport
         kk = j - betrtracer_vars%ngwmobile_tracer_groups
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            dtime_loc(c)   = betr_time%get_step_size()
            time_remain(c) = betr_time%get_step_size()
            update_col(c) = .true.
         enddo

         do
            !do diffusive transport
            call DiffusTransp(bounds, betr_status, lbj, ubj, jtops, num_soilc, &
                 filter_soilc, ntrcs, tracer_conc_mobile_col(:,:,difs_trc_group(1:ntrcs)),  &
                 hmconductance_col(:,:,j),                                             &
                 dtime_loc,                                                            &
                 dz,                                                                   &
                 source=local_source(:,:,1:ntrcs),                                     &
                 update_col=update_col,                                                &
                 dtracer=dtracer(:,:,1:ntrcs))
            if(betr_status%check_status())return
            !do tracer update
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if(update_col(c))then
                  !do negative tracer screening
                  lnegative_tracer = .false.
                  !loop through layers
                  do k = 1, ntrcs
                     trcid = difs_trc_group(k)
                     do l = jtops(c), ubj
                        if(tracer_conc_mobile_col(c,l,trcid)<-dtracer(c,l,k))then
                           !if the tracer update is very tinty, then set it to zero
                           if(abs(dtracer(c,l,k))<tiny_val)dtracer(c,l,k) = 0._r8

                           if(tracer_conc_mobile_col(c,l,trcid)<0._r8)then
                              write(msg,*)'nstep=',betr_time%get_nstep(),'col=',c,'l=',l,'trcid=',trcid, &
                                   new_line('A')//'dtime=',dtime_loc(c), &
                                   new_line('A')//'trc=',tracer_conc_mobile_col(c,l,trcid),&
                                   'dtracer=',dtracer(c,l,k), &
                                   new_line('A')//'stopped for negative tracer '//&
                                   trim(betrtracer_vars%tracernames(j))//' '//trim(subname)//errMsg(mod_filename, __LINE__)
                              call betr_status%set_msg(msg=msg,err=-1)
                              return
                           endif

                           !if tracer concentration change is zero, then goto next layer
                           if(dtracer(c,l,k)==0._r8)cycle
                           !now negative tracer occurs, decrease the timestep and exit the layer loop
                           lnegative_tracer = .true.
                           dtime_loc(c) = dtime_loc(c)*0.5_r8
                           exit
                        endif
                     enddo
                     !negative tracer, ramp out the loop
                     if(lnegative_tracer)exit

                     !do error budget for good calculation
                     call daxpy(ubj-jtops(c)+1, 1._r8, dtracer(c,jtops(c):ubj,k), 1, &
                          tracer_conc_mobile_col(c,jtops(c):ubj,trcid),1)

                     err_tracer(c, k) = dot_sum(dtracer(c,jtops(c):ubj, k), dz(c,jtops(c):ubj),betr_status)
                     if(betr_status%check_status())return
                     if(abs(err_tracer(c,k))>=err_min_solid)then
                        call betr_status%set_msg('mass balance error for tracer ' &
                           //tracernames(trcid)//' in '//trim(subname)//errMsg(mod_filename, __LINE__), err=-1)
                        return
                     endif
                  enddo
                  !if negative tracer concentration is found, go to the next column
                  if(lnegative_tracer)cycle

                  !when everything is OK, update the remaining time to be evolved.
                  time_remain(c) = time_remain(c)-dtime_loc(c)
                  dtime_loc(c)   = max(dtime_loc(c),dtime_min)
                  dtime_loc(c)   = min(dtime_loc(c), time_remain(c))
               endif

            enddo

            !test for loop exit
            lexit_loop=exit_loop_by_threshold(bounds%begc, bounds%endc, time_remain, &
                 dtime_min, num_soilc, filter_soilc, update_col)
            if(lexit_loop)exit
         enddo
      enddo
      deallocate(difs_trc_group)
      deallocate(dtracer)
      deallocate(err_tracer)
      deallocate(local_source)
    end associate

  end subroutine tracer_solid_transport

  !-------------------------------------------------------------------------------
  subroutine tracer_gw_transport(betr_time, bounds, lbj, ubj, col, pft, jtops, num_soilc, filter_soilc,&
       Rfactor, dz, zi,  transp_pathway, advection_on, diffusion_on, betrtracer_vars, &
       tracerboundarycond_vars, tracercoeff_vars, biophysforc, biogeo_flux, bgc_reaction, &
       tracerstate_vars, tracerflux_vars, bstatus)
    !
    ! !DESCRIPTION:
    ! do dual-phase (gas+aqueous) vertical tracer transport
    !
    ! !USES:
    use tracerstateType        , only : tracerstate_type
    use tracerboundarycondtype , only : tracerboundarycond_type
    use tracerfluxtype         , only : tracerflux_type
    use tracercoeffType        , only : tracercoeff_type
    use BGCReactionsMod        , only : bgc_reaction_type
    use BeTR_TimeMod           , only : betr_time_type
    use BetrStatusType         , only : betr_status_type
    use betr_columnType        , only : betr_column_type
    use BeTR_PatchType         , only : betr_patch_type
    implicit none
    ! !ARGUMENTS:
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: lbj, ubj
    type(betr_column_type)           , intent(in)    :: col
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc                           ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)                     ! column filter_soilc
    integer                          , intent(in)    :: jtops(bounds%begc: )                ! top label of each column
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars
    real(r8)                         , intent(in)    :: dz(bounds%begc: ,lbj: )             !
    real(r8)                         , intent(in)    :: zi(bounds%begc: ,lbj-1: )           !
    real(r8)                         , intent(in)    :: Rfactor(bounds%begc: ,lbj: ,1: )    !rfactor for dual diffusive transport
    integer                          , intent(in)    :: transp_pathway(2)                   !the pathway vector
    logical                          , intent(in)    :: advection_on
    logical                          , intent(in)    :: diffusion_on
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
    type(tracerboundarycond_type)    , intent(in)    :: tracerboundarycond_vars
    class(bgc_reaction_type)         , intent(inout) :: bgc_reaction
    type(tracercoeff_type)           , intent(inout) :: tracercoeff_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(betr_status_type)           , intent(out)   :: bstatus
    ! !LOCAL VARIABLES:
    integer ::  kk
    integer :: jtops0(bounds%begc:bounds%endc)
    real(r8) :: dtime ! model time step
    character(len=255) :: subname = 'tracer_gw_transport'

    call bstatus%reset()
    SHR_ASSERT_ALL((ubound(jtops)     == (/bounds%endc/)),                          errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,1)      == bounds%endc),                              errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,2)      == ubj),                                      errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(zi,1)      == bounds%endc),                              errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(zi,2)      == ubj),                                      errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,1) == bounds%endc),                              errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,2) ==  ubj ),                                    errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,3) ==  betrtracer_vars%ngwmobile_tracer_groups), errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return

    dtime = betr_time%get_step_size()
    !
    !Exclude solid phase tracers, by doing tracer equilibration
    !This is equivalent to do aqueous chemistry without biological production/consumption
    !The reason for doing this is to account for change in phase
    !partitioning due to change in hydrological status.

    call bgc_reaction%do_tracer_equilibration(bounds, lbj, ubj, &
         jtops,                                                 &
         num_soilc,                                             &
         filter_soilc,                                          &
         betrtracer_vars,                                       &
         tracercoeff_vars, tracerstate_vars,bstatus)
    if(bstatus%check_status())return
    !do diffusive and advective transport, assuming aqueous and gaseous phase are in equilbrium
    do kk = 1 , 2
       if (transp_pathway(kk) == diffusion_scheme .and. diffusion_on) then
          call do_tracer_gw_diffusion(bounds, lbj, ubj,                                    &
               jtops,                                                                      &
               num_soilc,                                                                  &
               filter_soilc,                                                               &
               betrtracer_vars,                                                            &
               tracerboundarycond_vars,                                                    &
               Rfactor,                                                                    &
               tracercoeff_vars%hmconductance_col(bounds%begc:bounds%endc, lbj:ubj-1, : ), &
               dz,                                                                         &
               dtime,                                                                      &
               tracerstate_vars,                                                           &
               tracerflux_vars, bstatus)
          if(bstatus%check_status())return
       elseif (transp_pathway(kk) == advection_scheme .and. advection_on)then
          jtops0(:) = 1
          call do_tracer_advection(betr_time, bounds, lbj, ubj, &
               col, pft, &
               jtops0,                                          &
               num_soilc,                                       &
               filter_soilc,                                    &
               betrtracer_vars,                                 &
               dz,                                              &
               zi,                                              &
               dtime,                                           &
               biophysforc,                                     &
               biogeo_flux,                                     &
               tracercoeff_vars,                                &
               tracerstate_vars,                                &
               tracerflux_vars, bstatus)
           if(bstatus%check_status())return
       endif
    enddo

  end subroutine tracer_gw_transport

  !-------------------------------------------------------------------------------
  subroutine do_tracer_advection(betr_time, bounds, lbj, ubj, col, pft, &
       jtops, num_soilc, filter_soilc, betrtracer_vars, dz, zi, dtime,  &
       biophysforc, biogeo_flux, tracercoeff_vars, tracerstate_vars, tracerflux_vars, bstatus)
    !
    ! !DESCRIPTION:
    ! do aqueous advection for dissolved tracers, the advection of gasesous phase is done through pressure
    ! adjustment for ebullition
    ! the aquesou advection is formulated as
    ! \frac{\partial{p(vsm*C_aq)}}{\partial t} = - \frac{\partial u*C_aq}{\partial z} + S_{root vs soil}
    !
    ! now the code transport tracers by groups specified by phase conversion coefficients
    ! !USES:
    use tracerstateType , only : tracerstate_type
    use tracerfluxtype  , only : tracerflux_type
    use TracerCoeffType , only : tracercoeff_type
    use TransportMod    , only : semi_lagrange_adv_backward, set_debug_transp
    use MathfuncMod     , only : safe_div
    use BeTR_TimeMod    , only : betr_time_type
    use BetrStatusType  , only : betr_status_type
    use betr_constants  , only : betr_errmsg_len
    use betr_columnType , only : betr_column_type
    use BeTR_PatchType  , only : betr_patch_type
    implicit none
    !ARGUMENTS
    class(betr_time_type)            , intent(in)    :: betr_time!
    type(bounds_type)                , intent(in)    :: bounds
    type(betr_column_type)           , intent(in)    :: col
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: lbj, ubj
    integer                          , intent(in)    :: num_soilc                           ! number of columns in column filter_soilc
    integer                          , intent(in)    :: jtops(bounds%begc: )
    integer                          , intent(in)    :: filter_soilc(:)                     ! column filter_soilc
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars
    real(r8)                         , intent(in)    :: dz(bounds%begc: ,lbj: )
    real(r8)                         , intent(in)    :: zi(bounds%begc: ,lbj-1: )
    real(r8)                         , intent(in)    :: dtime                               !model time step
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)      , intent(in)    :: biogeo_flux
    type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(betr_status_type)           , intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    logical              :: update_col(bounds%begc:bounds%endc)
    real(r8)             :: time_remain(bounds%begc:bounds%endc)
    real(r8)             :: dtime_loc(bounds%begc:bounds%endc)
    real(r8)             :: denum, num
    real(r8)             :: qflx_adv_local(bounds%begc:bounds%endc,lbj-1:ubj)
    real(r8)             :: qflx_rootsoi_local(bounds%begc:bounds%endc,lbj:ubj) !
    integer, allocatable :: adv_trc_group( : )
    real(r8), pointer    :: err_tracer( : , : )
    real(r8), pointer    :: transp_mass_vr( : , : , : )
    real(r8), pointer    :: transp_mass(:, :)
    real(r8), pointer    :: leaching_mass( : , : )
    real(r8), pointer    :: inflx_top( : , : )
    real(r8), pointer    :: inflx_bot( : , : )
    real(r8), pointer    :: dmass( : , : )
    real(r8), pointer    :: trc_bot(  :, : )
    real(r8), pointer    :: trc_conc_out(:,:,:)
    real(r8), pointer    :: seep_mass(: , :)
    logical              :: halfdt_col(bounds%begc:bounds%endc)
    real(r8)             :: err_relative
    real(r8)             :: c_courant
    integer              :: num_loops                                           !number of loops as determined by the courant number condition
    logical              :: lexit_loop
    integer              :: c, p, fc, j, l, k, ntrcs, trcid,kk
    integer              :: ngwmobile_tracers
    logical              :: lshock
    real(r8), parameter  :: err_relative_threshold=1.e-2_r8                     !relative error threshold
    real(r8), parameter  :: err_adv_min=1.e-10_r8
    real(r8), parameter  :: loc_eps = 1.e-12_r8                                 !smoothing factor to avoid advection velocity spikes, dimension less
    real(r8)             :: mass0
    real(r8)             :: alpha
    character(len=255)   :: subname = 'do_tracer_advection'
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,1)  == bounds%endc),           errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,2)  == ubj),                   errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(zi,1)  == bounds%endc),           errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(zi,2)  == ubj),                   errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return

    associate(                                                                         &
         qflx_adv                  => biogeo_flux%qflx_adv_col                       , & !real(r8) (:,:)[intent(in)], advective velocity defined at layer interfatemperature_vars
         qflx_rootsoi              => biophysforc%qflx_rootsoi_col                   , & !real(r8) (:,:)[intent(in)], water flux between plant and soil at different layers
         qflx_rootsoi_frac_patch   => biophysforc%qflx_rootsoi_frac_patch            , &
         h2osoi_liqvol             => biophysforc%h2osoi_liqvol_col                  , & !real(r8) (:,:)[intent(in)]
         is_advective              => betrtracer_vars%is_advective                   , & !logical(:) [intent(in)], indicator whether the tracer undergoes advection
         is_mobile                 => betrtracer_vars%is_mobile                      , & !
         is_h2o                    => betrtracer_vars%is_h2o                         , & !logical(:) [intent(in)], indicator whether the tracer is h2o
         vtrans_scal               => betrtracer_vars%vtrans_scal                    , & !real(r8) (:) [intent(in)], transport scalar for tracer exchaning between root and soil
         ngwmobile_tracer_groups   => betrtracer_vars%ngwmobile_tracer_groups        , & !integer [intent(in)], number of mobile tracers undergoing dual phase transport
         nmem_max                  => betrtracer_vars%nmem_max                       , & !
         tracer_group_memid        => betrtracer_vars%tracer_group_memid             , & !
         tracernames               => betrtracer_vars%tracernames                    , & !
         tracer_conc_mobile_col    => tracerstate_vars%tracer_conc_mobile_col        , & !
         tracer_conc_grndwater_col => tracerstate_vars%tracer_conc_grndwater_col     , & !
         aqu2bulkcef_mobile_col    => tracercoeff_vars%aqu2bulkcef_mobile_col        , & !
         tracer_flx_leaching       => tracerflux_vars%tracer_flx_leaching_col        , & !
         tracer_flx_surfrun        => tracerflux_vars%tracer_flx_surfrun_col         , & !
         tracer_flx_vtrans         => tracerflux_vars%tracer_flx_vtrans_col          , & !
         tracer_flx_vtrans_patch   => tracerflux_vars%tracer_flx_vtrans_patch        , & !
         tracer_flx_vtrans_vr      => tracerflux_vars%tracer_flx_vtrans_vr_col       , & !
         tracer_flx_infl           => tracerflux_vars%tracer_flx_infl_col              & !
         )
      !allocate memories
      allocate (adv_trc_group (nmem_max                                    ))
      allocate (err_tracer    (bounds%begc:bounds%endc ,nmem_max           ))
      allocate (transp_mass_vr(bounds%begc:bounds%endc, lbj:ubj, nmem_max  ))
      allocate (transp_mass   (bounds%begc:bounds%endc, nmem_max           ))
      allocate (leaching_mass (bounds%begc:bounds%endc, nmem_max           ))
      allocate (inflx_top     (bounds%begc:bounds%endc, nmem_max           ))
      allocate (inflx_bot     (bounds%begc:bounds%endc, nmem_max           ))
      allocate (trc_bot       (bounds%begc:bounds%endc, nmem_max           ))
      allocate (seep_mass     (bounds%begc:bounds%endc, nmem_max           ))
      allocate (dmass         (bounds%begc:bounds%endc,nmem_max            ))
      allocate (trc_conc_out  (bounds%begc:bounds%endc,lbj:ubj, 1:nmem_max ))

      !initialize local variables
      update_col  (:) = .true.
      time_remain (:) = 0._r8
      dtime_loc   (:) = 0._r8
      !loop over all tracers
      do j = 1, ngwmobile_tracer_groups
         ntrcs = 0
         adv_trc_group(:) = 0
         do k = 1, nmem_max
            trcid = tracer_group_memid(j,k)
            if(trcid>0)then
               if(is_mobile(trcid) .and. is_advective(trcid)) then
                  ntrcs = ntrcs + 1
                  adv_trc_group(ntrcs) = trcid
               endif
            endif
         enddo
         !convert bulk mobile phase into aqueous phase'
         if(ntrcs==0)cycle
         do k = 1, ntrcs
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               inflx_top(c, k) = tracer_flx_infl(c,adv_trc_group(k))
               !set to 0 to ensure outgoing boundary condition is imposed, this may not be correct for water isotopes
               inflx_bot(c,k) = 0._r8
               trc_bot(c,k) = tracer_conc_grndwater_col(c,adv_trc_group(k))
            enddo
         enddo

         !obtain advective velocity for the tracer group
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            qflx_adv_local(c,jtops(c)-1) = qflx_adv(c,jtops(c)-1)
            do l = jtops(c), ubj
              if(l<ubj)then
                if(qflx_adv(c,l) > 0._r8)then
                  qflx_adv_local(c,l) = safe_div(qflx_adv(c,l),aqu2bulkcef_mobile_col(c,l,j),eps=loc_eps)
                else
                  qflx_adv_local(c,l) = safe_div(qflx_adv(c,l),aqu2bulkcef_mobile_col(c,l+1,j),eps=loc_eps)
                endif
              else
                if(qflx_adv(c,l) > 0._r8)then
                  qflx_adv_local(c,l) = safe_div(qflx_adv(c,l),aqu2bulkcef_mobile_col(c,l,j),eps=loc_eps)
                else
                  qflx_adv_local(c,l) = qflx_adv(c,l)
                endif
               endif
               if(is_h2o(j))then
                 qflx_rootsoi_local(c,l) = qflx_rootsoi(c,l)
               else
                 qflx_rootsoi_local(c,l) = safe_div(qflx_rootsoi(c,l),aqu2bulkcef_mobile_col(c,l,j),eps=loc_eps)
               endif
            enddo
         enddo

         !dertmine the local advection time step based on the existence of convergence grid cell, i.e.
         ! grid cells with ul * ur < 0
         ! note qflx_adv(c,jtops(c)-1) is defined with infiltration
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            dtime_loc(c)=dtime              !local advective time step

            !initialize the time keeper and make sure all columns are updated initially
            update_col(c)=.true.
            time_remain(c) = dtime
         enddo

         do
            !zero leaching flux, leaching is outgoing only.
            leaching_mass=0._r8

            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if(update_col(c))then
                  do k = 1, ntrcs
                     trcid = adv_trc_group(k)
                     dmass(c, k) = dot_sum(tracer_conc_mobile_col(c, jtops(c):ubj, trcid), &
                      dz(c, jtops(c):ubj),bstatus)
                     if(bstatus%check_status())return
                  enddo
               endif
            enddo

            ! do semi-lagrangian tracer transport
            call semi_lagrange_adv_backward(bounds, bstatus,lbj, ubj,                             &
                 jtops,                                                                           &
                 num_soilc,                                                                       &
                 filter_soilc,                                                                    &
                 ntrcs,                                                                           &
                 dtime_loc,                                                                       &
                 dz,                                                                              &
                 zi,                                                                              &
                 qflx_adv_local(bounds%begc:bounds%endc,lbj-1:ubj),                               &
                 inflx_top(bounds%begc:bounds%endc, 1:ntrcs),                                     &
                 inflx_bot(bounds%begc:bounds%endc, 1:ntrcs),                                     &
                 trc_bot(bounds%begc:bounds%endc,   1:ntrcs),                                     &
                 update_col,                                                                      &
                 halfdt_col,                                                                      &
                 tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj,adv_trc_group(1:ntrcs)), &
                 trc_conc_out(:,:,1:ntrcs),                                                       &
                 leaching_mass(bounds%begc:bounds%endc,1:ntrcs), seep_mass(bounds%begc:bounds%endc, 1:ntrcs))

            if(bstatus%check_status())return
            !do soil-root tracer exchange
            do k = 1, ntrcs
               trcid = adv_trc_group(k)

               do l = lbj, ubj
                  do fc = 1, num_soilc
                     c = filter_soilc(fc)
                     if(update_col(c) .and. (.not. halfdt_col(c)) .and. l>=jtops(c))then
                        tracer_conc_mobile_col(c,l,trcid)=trc_conc_out(c,l,k)
                     endif
                  enddo
               enddo
               transp_mass_vr(:,:, k) = 0._r8
               transp_mass(:, k) = 0._r8

               if(vtrans_scal(trcid)>0._r8)then
                  call calc_root_uptake_as_perfect_sink(bounds, lbj, ubj, num_soilc,   &
                       filter_soilc,                                                   &
                       dtime_loc,                                                      &
                       dz,                                                             &
                       qflx_rootsoi_local,                                             &
                       vtrans_scal(trcid),                                             &
                       is_h2o(trcid),                                                  &
                       update_col,                                                     &
                       halfdt_col,                                                     &
                       tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj,trcid), &
                       transp_mass_vr(bounds%begc:bounds%endc, lbj:ubj, k)           , &
                       transp_mass(bounds%begc:bounds%endc, k),bstatus)
                 if(bstatus%check_status())return
               endif
            enddo

            do k = 1, ntrcs
               trcid = adv_trc_group(k)
               do fc = 1, num_soilc
                  c = filter_soilc(fc)

                  if(update_col(c) .and. (.not. halfdt_col(c)))then
                     mass0   = dmass(c, k)
                     dmass(c, k) =  dot_sum(tracer_conc_mobile_col(c,jtops(c):ubj,trcid), &
                        dz(c,jtops(c):ubj),bstatus)- dmass(c, k)
                     if(bstatus%check_status())return
                     err_tracer(c, k) = dmass(c, k) - inflx_top(c,k) * dtime_loc(c) + leaching_mass(c,k) + &
                          transp_mass(c, k) + seep_mass(c,k)
                     if(abs(err_tracer(c,k))<err_adv_min .or. abs(err_tracer(c,k))/(mass0+1.e-10_r8) < 1.e-10_r8)then
                        !when the absolute value is too small, set relative error to
                        err_relative = err_relative_threshold*0.999_r8
                     else
                        err_relative = err_tracer(c,k)/maxval((/abs(inflx_top(c,k)*dtime_loc(c)), abs(leaching_mass(c,k)),&
                          abs(dmass(c,k)),tiny_val/))
                     endif
                     if(abs(err_relative)<err_relative_threshold)then
                        leaching_mass(c,k) = leaching_mass(c,k) - err_tracer(c,k)
                     else
                        write(msg,'(2(A,1X,I8),5X,A,7(5X,A,5X,E18.10))')'nstep=', betr_time%get_nstep(), ', col=',c, &
                             tracernames(trcid),' err=',err_tracer(c,k),&
                             ' transp=',transp_mass(c,k),' lech=',&
                             leaching_mass(c,k),' infl=',inflx_top(c,k),' dmass=',dmass(c,k), ' mass0=', &
                             mass0,'err_rel=',err_relative

                        msg=trim(msg)//new_line('A')//'advection mass balance error for tracer '//tracernames(j) &
                          //new_line('A')//errMsg(mod_filename, __LINE__)
                        call bstatus%set_msg(msg, err=-1)
                        return
                     endif
                     tracer_flx_vtrans(c, trcid)  = tracer_flx_vtrans(c,trcid) + transp_mass(c,k)
                     tracer_flx_leaching(c,trcid) = tracer_flx_leaching(c, trcid) + leaching_mass(c,k)
                     tracer_flx_surfrun(c,trcid)  = tracer_flx_surfrun(c, trcid) + seep_mass(c,k)
                     tracer_flx_vtrans_vr(c, lbj:ubj, trcid) = transp_mass_vr(c,lbj:ubj,k)

                     !the following implementation assumes there is only one column
                     do p = 1, pft%npfts
                       tracer_flx_vtrans_patch(p,trcid) = 0._r8
                       do l=lbj, ubj
                         tracer_flx_vtrans_patch(p,trcid) = tracer_flx_vtrans_patch(p,trcid) + &
                            tracer_flx_vtrans_vr(c, l, trcid) * qflx_rootsoi_frac_patch(p,l)
                       enddo
                     enddo
                  endif
               enddo
            enddo

            do fc = 1, num_soilc
               c = filter_soilc(fc)
               if(update_col(c))then
                  if(halfdt_col(c))then
                     dtime_loc(c) = max(dtime_loc(c)*0.5_r8,dtime_min)
                     dtime_loc(c) = min(dtime_loc(c), time_remain(c))
                  else
                     time_remain(c) = time_remain(c) - dtime_loc(c)
                  endif
               endif
            enddo
            ! do loop control test
            lexit_loop=exit_loop_by_threshold(bounds%begc, bounds%endc, time_remain, &
                 dtime_min, num_soilc, filter_soilc, update_col)

            if(lexit_loop)exit
         enddo
      enddo
      deallocate(adv_trc_group)
      deallocate(err_tracer   )
      deallocate(transp_mass  )
      deallocate(transp_mass_vr)
      deallocate(leaching_mass)
      deallocate(inflx_top    )
      deallocate(inflx_bot    )
      deallocate(dmass        )
      deallocate(trc_conc_out )
      deallocate(trc_bot)
      deallocate(seep_mass)
    end associate
  end subroutine do_tracer_advection

  !-------------------------------------------------------------------------------
  subroutine do_tracer_gw_diffusion(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, ttracerboundarycond_vars, Rfactor,                            &
       hmconductance_col, dz, dtime, tracerstate_vars, tracerflux_vars, bstatus)
    !
    !  !DESCRIPTION:
    !
    ! do diffusive tracer transport
    ! using a modified version from Tang and Riley (BG, 2014) for the dual diffusive transport. In this version
    ! the bulk concentration is used as state variable, rfactor is integrated into the solver
    ! !USES:
    use tracerstateType        , only : tracerstate_type
    use tracerboundarycondtype , only : tracerboundarycond_type
    use tracerfluxtype         , only : tracerflux_type
    use TransportMod           , only : DiffusTransp, get_cntheta
    use tracer_varcon          , only : bndcond_as_conc
    use betr_constants         , only : betr_errmsg_len
    use BetrStatusType         , only : betr_status_type
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)             , intent(in)    :: bounds
    integer                       , intent(in)    :: lbj, ubj
    integer                       , intent(in)    :: jtops(bounds%begc: )                        ! top label of each column
    integer                       , intent(in)    :: num_soilc                                   ! number of columns in column filter_soilc
    integer                       , intent(in)    :: filter_soilc(:)                             ! column filter_soilc
    class(betrtracer_type)        , intent(in)    :: betrtracer_vars
    real(r8)                      , intent(in)    :: hmconductance_col(bounds%begc: , lbj: ,1: ) !weighted bulk conductance
    real(r8)                      , intent(in)    :: Rfactor(bounds%begc: ,lbj:  ,1: )           !rfactor for dual diffusive transport
    real(r8)                      , intent(in)    :: dz(bounds%begc: ,lbj: )
    real(r8)                      , intent(in)    :: dtime                                       !model time step
    type(tracerboundarycond_type) , intent(in)    :: ttracerboundarycond_vars
    type(tracerstate_type)        , intent(inout) :: tracerstate_vars
    type(tracerflux_type)         , intent(inout) :: tracerflux_vars
    type(betr_status_type)        , intent(out)   :: bstatus
    !
    ! !LOCAL VARIABLES:
    character(len=255)    :: subname = 'do_tracer_gw_diffusion'
    integer               :: j, fc, c, l, ntrcs, k, trcid
    logical               :: update_col(bounds%begc:bounds%endc)   !logical switch for whether or not update a column
    real(r8)              :: time_remain(bounds%begc:bounds%endc)  !remaining time to evolve
    real(r8)              :: dtime_loc(bounds%begc:bounds%endc)    !local time step
    logical               :: lnegative_tracer                      !when true, negative tracer occurs
    logical               :: lexit_loop
    real(r8)              :: err_relative
    integer , allocatable :: dif_trc_group(:)
    real(r8), allocatable :: err_tracer(:, :)
    real(r8), allocatable :: diff_surf(:, :)
    real(r8), allocatable :: dtracer(:, :, :)
    real(r8), allocatable :: dmass(:, : )
    real(r8), allocatable :: local_source(:, :, :)
    real(r8)              :: mass0,mass1
    real(r8), parameter   :: err_relative_threshold=1.e-2_r8 !relative error threshold
    real(r8), parameter   :: err_dif_min = 1.e-12_r8  !minimum absolute error
    integer :: ntracer_groups
    character(len=betr_errmsg_len) :: msg, msg1


    call bstatus%reset()
    ntracer_groups = betrtracer_vars%ntracer_groups

    SHR_ASSERT_ALL((ubound(jtops)               == (/bounds%endc/)),     errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,1)                == bounds%endc),      errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dz,2)                == ubj),              errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(hmconductance_col,1) == bounds%endc),    errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(hmconductance_col,2) == ubj-1),           errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(hmconductance_col,3) == ntracer_groups),    errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,1)           == bounds%endc),     errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,2)           == ubj),           errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,3) == betrtracer_vars%ngwmobile_tracer_groups), errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return

    associate(                                                                                 &
         is_volatile              =>  betrtracer_vars%is_volatile                            , & !
         is_mobile                =>  betrtracer_vars%is_mobile                              , & !
         is_diffusive             =>  betrtracer_vars%is_diffusive                           , & !
         volatileid               =>  betrtracer_vars%volatileid                             , & !
         tracernames              =>  betrtracer_vars%tracernames                            , & !
         nmem_max                 =>  betrtracer_vars%nmem_max                               , & !
         ngwmobile_tracers        =>  betrtracer_vars%ngwmobile_tracers                      , & !
         ngwmobile_tracer_groups  => betrtracer_vars%ngwmobile_tracer_groups                 , & !
         tracer_group_memid       =>  betrtracer_vars%tracer_group_memid                     , & !
         tracer_flx_dif           =>  tracerflux_vars%tracer_flx_dif_col                     , & !
         tracer_flx_netpro_vr     => tracerflux_vars%tracer_flx_netpro_vr_col                , & !
         tracer_gwdif_concflux_top=>  ttracerboundarycond_vars%tracer_gwdif_concflux_top_col , & !
         condc_toplay             =>  ttracerboundarycond_vars%condc_toplay_col              , & !
         topbc_type               =>  ttracerboundarycond_vars%topbc_type                    , & !
         bot_concflux             =>  ttracerboundarycond_vars%bot_concflux_col              , & !
         tracer_conc_mobile_col   => tracerstate_vars%tracer_conc_mobile_col                   &
         )

      allocate (err_tracer    (bounds%begc:bounds%endc, nmem_max         ))
      allocate (diff_surf     (bounds%begc:bounds%endc, nmem_max         ))
      allocate (dtracer       (bounds%begc:bounds%endc,lbj:ubj, nmem_max ))
      allocate (dmass         (bounds%begc:bounds%endc, nmem_max         ))
      allocate (local_source  (bounds%begc:bounds%endc,lbj:ubj, nmem_max ))
      allocate (dif_trc_group (nmem_max                                  ))

      update_col(:)       = .true.
      time_remain(:)      = 0._r8
      dtime_loc(:)        = 0._r8
      local_source(:,:,:) = 0._r8

      do j = 1, ngwmobile_tracer_groups

         !assemable the tracer group for diffusion
         ntrcs = 0
         dif_trc_group(:) = 0
         do k = 1, nmem_max
            trcid = tracer_group_memid(j,k)
            if(trcid>0)then
               if(is_mobile(trcid) .and. is_diffusive(trcid)) then
                  ntrcs = ntrcs + 1
                  dif_trc_group(ntrcs) = trcid
               endif
            endif
         enddo
         if(ntrcs==0)cycle
         !initialize the time keeper
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            time_remain(c) = dtime
            dtime_loc(c) = dtime
            update_col(c) = .true.
         enddo

         !Do adpative time stepping to avoid negative tracer
         do
            call DiffusTransp(bounds, bstatus, lbj, ubj, jtops,                              &
                 num_soilc,                                                                  &
                 filter_soilc, ntrcs,                                                        &
                 tracer_conc_mobile_col( : , lbj:ubj ,dif_trc_group(1:ntrcs)),               &
                 Rfactor( : , lbj:ubj ,j),                                                   &
                 hmconductance_col( : , lbj:ubj-1 ,j), dtime_loc, dz(:,lbj:ubj),             &
                 local_source(:,lbj:ubj, 1:ntrcs),                                           &
                 tracer_gwdif_concflux_top( : , : ,dif_trc_group(1:ntrcs)),                  &
                 condc_toplay( : ,j),                                                        &
                 topbc_type(j),                                                              &
                 bot_concflux( : , : ,dif_trc_group(1:ntrcs)),                               &
                 update_col,                                                                 &
                 dtracer(:,lbj:ubj,1:ntrcs))

            if(bstatus%check_status())return
            !do tracer update
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               !update the column
               if(update_col(c))then
                  !do negative tracer screening
                  lnegative_tracer = .false.

                  do k = 1, ntrcs
                     trcid = dif_trc_group(k)
                     !loop through all layers
                     do l = jtops(c), ubj
                        if(tracer_conc_mobile_col(c,l,trcid)<-dtracer(c,l,k))then
                           !if the tracer update is very tinty, then set it to zero
                           if(abs(dtracer(c,l,k))<tiny_val)dtracer(c,l,k) = 0._r8
                           if(tracer_conc_mobile_col(c,l,trcid)<0._r8)then
                              !write error message and stop
                              write(msg,*) 'trcname=',tracernames(trcid),'c=',c,'l=',l, &
                                 new_line('A')//'trc=',tracer_conc_mobile_col(c,l,trcid),'dtracer=',dtracer(c,l,k),&
                                 'dtime=',dtime_loc(c)
                              msg=trim(msg)//'stopped '//trim(subname)//errMsg(mod_filename, __LINE__)
                              call bstatus%set_msg(msg=msg,err=-1)
                              return
                           endif

                           !if tracer concentration change is zero, then goto next layer
                           if(dtracer(c,l,k)==0._r8)cycle
                           !now negative tracer occurs, decrease the timestep and exit the layer loop
                           lnegative_tracer = .true.
                           dtime_loc(c) = dtime_loc(c)*0.5_r8
                           exit
                        endif
                     enddo
                     !negative tracer, ramp out the loop
                     if(lnegative_tracer)exit
                  enddo

                  !time stepping screening
                  if(dtime_loc(c)<1.e-3_r8)then
                     write(msg,*)'diffusion time step < 1.e-3_r8', dtime_loc(c), 'col ',c
                     do k = 1, ntrcs
                        write(msg1,*)'tracer '//tracernames(trcid),get_cntheta(), &
                           (l,tracer_conc_mobile_col(c,l,trcid),l=jtops(c),ubj)
                        msg=trim(msg)//new_line('A')//trim(msg1)
                        write(msg1,*)'dtracer', (l,dtracer(c,l,k),l=jtops(c),ubj)
                        msg=trim(msg)//new_line('A')//trim(msg1)
                     enddo
                     msg = trim(msg)//new_line('A')//'stopped '//trim(subname)//errMsg(mod_filename, __LINE__)
                     call bstatus%set_msg(msg=msg, err=-1)
                     return
                  endif

                  !if negative tracer concentration is found, go to the next column
                  if(lnegative_tracer)cycle

                  !update surface efflux for volatile tracers, positive for into the atmosphere
                  do k = 1, ntrcs
                     trcid = dif_trc_group(k)
                     if(topbc_type(j)==bndcond_as_conc)then
                        diff_surf(c, k) = -condc_toplay(c,j) *                                 &
                             (tracer_conc_mobile_col(c,jtops(c),trcid)/Rfactor(c,jtops(c),j) - &
                              tracer_gwdif_concflux_top(c,1,trcid) +                           &
                              get_cntheta()*(dtracer(c,jtops(c),k)/Rfactor(c,jtops(c),j) +     &
                              tracer_gwdif_concflux_top(c,1,trcid)-                            &
                              tracer_gwdif_concflux_top(c,2,trcid)))
                     else
                        diff_surf(c, k) = 0.5_r8*(tracer_gwdif_concflux_top(c,1,trcid)+tracer_gwdif_concflux_top(c,2,trcid))
                     endif

                     !do error budget for good calculations

                     mass0=dot_sum(x=tracer_conc_mobile_col(c,jtops(c):ubj,trcid),y=dz(c,jtops(c):ubj),bstatus=bstatus)
                     if(bstatus%check_status())return
                     call daxpy(ubj-jtops(c)+1, 1._r8, dtracer(c,jtops(c):ubj, k), 1, &
                          tracer_conc_mobile_col(c,jtops(c):ubj,trcid), 1)

                     dmass(c, k) = dot_sum(x=dtracer(c,jtops(c):ubj, k),y=dz(c,jtops(c):ubj),bstatus=bstatus)
                     if(bstatus%check_status())return
                     err_tracer(c, k) = dmass(c, k)-diff_surf(c, k) *dtime_loc(c)

                     mass1=dot_sum(x=tracer_conc_mobile_col(c,jtops(c):ubj,trcid),y=dz(c,jtops(c):ubj),bstatus=bstatus)
                     if(bstatus%check_status())return
                     !calculate relative error, defined as the ratio between absolute error with respect to surface flux
                     if(abs(err_tracer(c, k))<err_dif_min .or.  abs(err_tracer(c, k))/(mass1+1.e-10_r8) < 1.e-10_r8)then
                        !when the absolute value is too small, set relative error to
                        err_relative = err_relative_threshold*0.999_r8
                     else
                        err_relative = err_tracer(c, k)/max(abs(diff_surf(c, k)),tiny_val)
                     endif

                     if(abs(err_relative)<err_relative_threshold)then
                        !the calculation is good, use the error to correct the diffusive flux for volatile tracer
                        if(is_volatile(trcid))then
                           diff_surf(c,k) = diff_surf(c,k)+err_tracer(c,k)/dtime_loc(c)
                           !accumulate the diffusive flux at the given time step, + into the atmosphere
                           tracer_flx_dif(c,volatileid(trcid)) = tracer_flx_dif(c,volatileid(trcid)) - &
                                diff_surf(c,k) * dtime_loc(c)
                        endif
                     else
                        write(msg,*) 'mass bal error dif '//trim(tracernames(trcid))//new_line('A'), 'mass 0/1',mass0, &
                             mass1,'col=',c, &
                             new_line('A')//'err=', err_tracer(c,k), 'dmass=',dmass(c,k), ' dif=', diff_surf(c,k)*dtime_loc(c), &
                             ' prod=',dot_sum(x=local_source(c,jtops(c):ubj,k),y=dz(c,jtops(c):ubj),bstatus=bstatus)*dtime_loc(c)

                        if(bstatus%check_status())return
                        msg=trim(msg)//new_line('A')//'mass balance error for tracer '//trim(tracernames(trcid))//' in ' &
                           //trim(subname)//new_line('A')//errMsg(mod_filename, __LINE__)
                        call bstatus%set_msg(msg=msg, err=-1)
                        return
                     endif
                  enddo
               endif
               !when everything is OK, update the remaining time to be evolved.
               time_remain(c) = time_remain(c)-dtime_loc(c)
               dtime_loc(c)=max(dtime_loc(c),dtime_min)
               dtime_loc(c)=min(dtime_loc(c), time_remain(c))
            enddo

            !test for loop exit
            lexit_loop=exit_loop_by_threshold(bounds%begc, bounds%endc, time_remain, &
                 dtime_min, num_soilc, filter_soilc, update_col)
            if(lexit_loop)exit
         enddo
      enddo
      !
      deallocate(err_tracer)
      deallocate(diff_surf)
      deallocate(dtracer)
      deallocate(dmass)
      deallocate(local_source)
      deallocate(dif_trc_group)
    end associate
  end subroutine do_tracer_gw_diffusion

  !-------------------------------------------------------------------------------
  function exit_loop_by_threshold(beg,end, datain, threshold, num_soilc, filter_soilc, update_col)result(lexit_loop)
    !
    ! !DESCRIPTION:
    !
    ! decide if the loop should be terminated based on threshold

    ! !ARGUMENTS:
    integer  , intent(in)    :: beg, end
    real(r8) , intent(in)    :: datain(beg:end)
    real(r8) , intent(in)    :: threshold
    integer  , intent(in)    :: num_soilc
    integer  , intent(in)    :: filter_soilc(:)
    logical  , intent(inout) :: update_col(beg:end)

    ! !LOCAL VARIABLES:
    logical :: lexit_loop
    integer :: fc, c
    character(len=255) :: subname='exit_loop_by_threshold'

    lexit_loop = .true.
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if(datain(c)<=threshold)then
          update_col(c)=.false.
       else
          lexit_loop =.false.
       endif
    enddo

  end function exit_loop_by_threshold

  !-------------------------------------------------------------------------------
  subroutine set_gwdif_Rfactor(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
      tracercoeff_vars, betrtracer_vars, Rfactor, betr_status)
    !
    ! !DESCRIPTION:
    ! set up the retardation factor
    ! !USES:
    use tracercoeffType       , only : tracercoeff_type
    use BetrStatusType        , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: lbj, ubj
    integer                , intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer                , intent(in)    :: num_soilc                 ! number of columns in column filter_soilc
    integer                , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    class(betrtracer_type) , intent(in)    :: betrtracer_vars
    type(tracercoeff_type) , intent(in)    :: tracercoeff_vars
    real(r8)               , intent(inout) :: Rfactor(bounds%begc: ,lbj:  ,1: )  !rfactor for dual diffusive transport
    type(betr_status_type) , intent(out)   :: betr_status
    ! !LOCAL VARIABLES:
    integer            :: j, fc, c, k, kk, trcid
    character(len=255) :: subname = 'set_gwdif_Rfactor'

    call betr_status%reset()
    associate(                                                                       &  !
         ngwmobile_tracer_groups =>    betrtracer_vars%ngwmobile_tracer_groups     , &  !
         tracer_group_memid      =>    betrtracer_vars%tracer_group_memid          , &  !
         is_volatile             =>    betrtracer_vars%is_volatile                 , &  !
         is_h2o                  =>    betrtracer_vars%is_h2o                      , &  !
         volatilegroupid         =>    betrtracer_vars%volatilegroupid             , &  !
         gas2bulkcef_mobile      =>    tracercoeff_vars%gas2bulkcef_mobile_col     , &  !
         aqu2bulkcef_mobile      =>    tracercoeff_vars%aqu2bulkcef_mobile_col       &  !
         )

    SHR_ASSERT_ALL((ubound(jtops)   == (/bounds%endc/)) , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,1) == bounds%endc)   , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,2) == ubj)           , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(Rfactor,3) ==  ngwmobile_tracer_groups), errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return


      do j = lbj, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if(j>=jtops(c))then
               do k = 1, ngwmobile_tracer_groups
                  trcid = tracer_group_memid(k,1)
                  if(is_volatile(trcid))then
                     kk = volatilegroupid(trcid)
                     Rfactor(c,j, k) = gas2bulkcef_mobile(c,j,kk)
                  else
                     Rfactor(c,j, k) = aqu2bulkcef_mobile(c,j,k)
                  endif
               enddo
            endif
         enddo
      enddo
    end associate
  end subroutine set_gwdif_Rfactor

  !-------------------------------------------------------------------------------
  subroutine calc_ebullition(bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       forc_psrf, zi, dz, dtime, fracice, zwt, betrtracer_vars,                &
       tracercoeff_vars, tracerstate_vars, tracer_flx_ebu_col,ebullition_on, betr_status)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use tracercoeffType       , only : tracercoeff_type
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use betr_varcon           , only : grav => bgrav, oneatm => boneatm
    use BetrStatusType        , only : betr_status_type
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: lbj, ubj
    integer                , intent(in)    :: jtops(bounds%begc: )                                                             ! top label of each column
    integer                , intent(in)    :: num_soilc                                                                        ! number of columns in column filter_soilc
    integer                , intent(in)    :: filter_soilc(:)                                                                  ! column filter_soilc
    real(r8)               , intent(in)    :: dtime
    real(r8)               , intent(in)    :: dz(bounds%begc: , lbj: )
    real(r8)               , intent(in)    :: zi(bounds%begc: , lbj-1: )
    real(r8)               , intent(in)    :: forc_psrf(bounds%begc: )                                                         ! atmospheric pressure, [Pa]
    real(r8)               , intent(in)    :: fracice(bounds%begc: , lbj: )                                                    ! fraction of ice in the soil layer, 0-1
    real(r8)               , intent(in)    :: zwt(bounds%begc: )                                                               ! water table depth [m]
    class(betrtracer_type) , intent(in)    :: betrtracer_vars                                                                  ! tracer info data structure
    type(tracercoeff_type) , intent(in)    :: tracercoeff_vars                                                                 ! tracer phase conversion coefficients
    type(tracerstate_type) , intent(inout) :: tracerstate_vars                                                                 ! tracer state variables data structure
    real(r8)               , intent(inout) :: tracer_flx_ebu_col(bounds%begc:bounds%endc, 1:betrtracer_vars%nvolatile_tracers) ! tracer ebullition
    logical                , intent(in)    :: ebullition_on
    type(betr_status_type) , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    real(r8), parameter :: icefrac_sealed=0.99_r8                         !set the sealing up ice fraction
    real(r8)            :: bubble_flux(betrtracer_vars%nvolatile_tracers) !bubble flux, mol/m2/s
    real(r8)            :: press_hydro
    real(r8)            :: n2_pressure
    real(r8)            :: o2_pressure
    real(r8)            :: ar_pressure
    real(r8)            :: co2_pressure
    real(r8)            :: ch4_pressure
    real(r8)            :: total_pressure
    real(r8)            :: frac
    integer             :: vid
    integer             :: fc, c, j, kk

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops)     == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(forc_psrf) == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(fracice,1) == bounds%endc) ,          errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(fracice,2) == ubj) ,                  errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(dz,1)      == bounds%endc),           errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(dz,2)      == ubj),                   errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return
    SHR_ASSERT_ALL((ubound(zwt)       == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),betr_status)
    if(betr_status%check_status())return

    ! remove compiler warnings about unused dummy args
    if (dtime > 0) continue

    associate(                                                                &
         tracer_conc_mobile_col   => tracerstate_vars%tracer_conc_mobile_col, &
         aqu2bulkcef_mobile_col   => tracercoeff_vars%aqu2bulkcef_mobile_col, &
         henrycef_col             => tracercoeff_vars%henrycef_col          , &
         volatileid               => betrtracer_vars%volatileid             , &
         groupid                  => betrtracer_vars%groupid                , &
         volatilegroupid          => betrtracer_vars%volatilegroupid        , &
         ngwmobile_tracers        => betrtracer_vars%ngwmobile_tracers      , &
         nvolatile_tracers        => betrtracer_vars%nvolatile_tracers      , &
         is_mobile                => betrtracer_vars%is_mobile              , &
         is_volatile              => betrtracer_vars%is_volatile            , &
         is_h2o                   => betrtracer_vars%is_h2o                 , &
         id_trc_n2                => betrtracer_vars%id_trc_n2              , &
         id_trc_o2                => betrtracer_vars%id_trc_o2              , &
         id_trc_ar                => betrtracer_vars%id_trc_ar              , &
         id_trc_ch4               => betrtracer_vars%id_trc_ch4             , &
         id_trc_co2x              => betrtracer_vars%id_trc_co2x              &
         )

      if(.not. ebullition_on)return
      if (.not. all(is_mobile((/id_trc_n2,id_trc_o2,id_trc_ar,id_trc_ch4,id_trc_co2x/)))) return
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         !initialize bubble flux to zero
         tracer_flx_ebu_col(c,:) = 0._r8
         !do not do anything if the soil is ice sealed.
         if(fracice(c,1)>=icefrac_sealed)cycle
         !initialize temporary bubble collecting vector
         bubble_flux(:) = 0._r8
         do j = ubj, 1, -1
            !calculate the imposed atmospheric pressure plus hydrostatic pressure from water
            press_hydro= max(zwt(c)-zi(c,j-1),0._r8)*denh2o*grav + forc_psrf(c)
            !convert Pa into atm
            press_hydro=press_hydro/oneatm
            !calculate the total gas pressure
            n2_pressure=calc_gas_pressure(tracer_conc_mobile_col(c,j,id_trc_n2),    &
                 aqu2bulkcef_mobile_col(c,j,groupid(id_trc_n2)),                    &
                 henrycef_col(c, j, volatilegroupid(id_trc_n2)))

            o2_pressure=calc_gas_pressure(tracer_conc_mobile_col(c,j,id_trc_o2),    &
                 aqu2bulkcef_mobile_col(c,j,groupid(id_trc_o2)),                    &
                 henrycef_col(c, j, volatilegroupid(id_trc_o2)))

            ar_pressure=calc_gas_pressure(tracer_conc_mobile_col(c,j,id_trc_ar),    &
                 aqu2bulkcef_mobile_col(c,j,groupid(id_trc_ar)),                    &
                 henrycef_col(c, j, volatilegroupid(id_trc_ar)))

            co2_pressure=calc_gas_pressure(tracer_conc_mobile_col(c,j,id_trc_co2x), &
                 aqu2bulkcef_mobile_col(c,j,groupid(id_trc_co2x)),                  &
                 henrycef_col(c, j, volatilegroupid(id_trc_co2x)))

            ch4_pressure=calc_gas_pressure(tracer_conc_mobile_col(c,j,id_trc_ch4),  &
                 aqu2bulkcef_mobile_col(c,j,groupid(id_trc_ch4)),                   &
                 henrycef_col(c, j, volatilegroupid(id_trc_ch4)))

            total_pressure=n2_pressure+o2_pressure+ar_pressure+co2_pressure+ch4_pressure

            if(total_pressure>press_hydro)then
               !calculate the fraction of gas to be released as bubble
               frac=(total_pressure-press_hydro)/total_pressure
               !note because there exisiting a relationship gas_conc*gas2bulkcef=bulk_con
               !a fraction of amount frac will be released as bubbles to move upward
               do kk = 1, ngwmobile_tracers
                  if(is_volatile(kk) .and. (.not. is_h2o(kk)))then
                     vid = volatileid(kk)
                     !the upward moving bubble flux
                     bubble_flux(vid) = tracer_conc_mobile_col(c,j,kk)*frac*dz(c,j)
                     !the readjusted tracer concentration
                     tracer_conc_mobile_col(c,j,kk) = tracer_conc_mobile_col(c,j,kk)*(1._r8-frac)
                     !add the bubbles to next layer above
                     if(j>1)then
                        tracer_conc_mobile_col(c,j-1,kk) = tracer_conc_mobile_col(c,j-1,kk) + bubble_flux(vid)/dz(c,j-1)
                        bubble_flux(vid) = 0._r8
                     endif
                  endif
               enddo
            endif
         enddo
         do kk = 1, nvolatile_tracers
            !+ into the atmosphere
            tracer_flx_ebu_col(c,kk) = bubble_flux(kk)
         enddo

      enddo
    end associate
  end subroutine calc_ebullition

  !-------------------------------------------------------------------------------
  function calc_gas_pressure(tracer_conc, aqu2bulkcef, henrycef) result(pres_atm)
    !
    ! !DESCRIPTION:
    ! Calculate gas pressure using given conditions

    ! !ARGUMENTS:
    real(r8), intent(in) :: tracer_conc   !tracer concentrations [mol/m3]
    real(r8), intent(in) :: aqu2bulkcef   !conversion parameter between aqueous and bulk tracer concentrations
    real(r8), intent(in) :: henrycef      !henry's law constant

    ! !LOCAL VARIABLES:
    real(r8) :: aqucon
    real(r8) :: pres_atm   !gas pressure [atm]
    !compuate aqueous concentration, mol/m3
    aqucon = tracer_conc / aqu2bulkcef
    !compuate partial pressure, atm
    pres_atm = aqucon * 1.e-3_r8 / henrycef

  end function calc_gas_pressure

  !-------------------------------------------------------------------------------
  subroutine calc_root_uptake_as_perfect_sink(bounds, lbj, ubj,  num_soilc, filter_soilc, &
       dtime_loc, dz, qflx_rootsoi, vtrans_scal, is_h2o, update_col, halfdt_col, tracer_conc, transp_mass_vr , &
       transp_mass,bstatus)
    !
    ! !DESCRIPTION:
    ! calculate plant aqueous tracer uptake through transpiration into xylem
    use BetrStatusType        , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: lbj, ubj
    integer           , intent(in)    :: num_soilc                           ! number of columns in column filter_soilc
    integer           , intent(in)    :: filter_soilc(:)                     ! column filter_soilc
    real(r8)          , intent(in)    :: dz(bounds%begc: , lbj: )            ! layer thickness
    real(r8)          , intent(in)    :: dtime_loc(bounds%begc: )
    real(r8)          , intent(in)    :: qflx_rootsoi(bounds%begc: , lbj: )
    real(r8)          , intent(in)    :: vtrans_scal
    logical           , intent(in)    :: is_h2o
    logical           , intent(in)    :: update_col(bounds%begc:bounds%endc) ! logical switch for active col update
    logical           , intent(in)    :: halfdt_col(bounds%begc:bounds%endc)
    real(r8)          , intent(inout) :: tracer_conc(bounds%begc: , lbj: )   ! incoming tracer concentration
    real(r8)          , intent(out)   :: transp_mass_vr(bounds%begc: , lbj: )!
    real(r8)          , intent(out)   :: transp_mass(bounds%begc: )
    type(betr_status_type) , intent(out)   :: bstatus
    ! !LOCAL VARIABLES:
    real(r8) :: tracer_conc_new
    integer  :: fc, c, j

    call bstatus%reset()
    SHR_ASSERT_ALL((ubound(dz)               == (/bounds%endc, ubj/)) , errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(update_col)       == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(dtime_loc)        == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(tracer_conc,1)    == bounds%endc) ,          errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(tracer_conc,2)    ==  ubj) ,                 errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(qflx_rootsoi,1)   == bounds%endc) ,          errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(qflx_rootsoi,2)   == ubj) ,                  errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(transp_mass_vr,1) == bounds%endc)      ,     errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(transp_mass_vr,2) == ubj)      ,             errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(transp_mass)      == (/bounds%endc/))      , errMsg(mod_filename,__LINE__),bstatus)
    if(bstatus%check_status())return

    transp_mass(:) = 0._r8
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       if(update_col(c) .and. (.not. halfdt_col(c)))then
          do j = 1, ubj
             if(is_h2o)then
               transp_mass_vr(c,j)   = vtrans_scal*qflx_rootsoi(c,j)*denh2o*dtime_loc(c)
               tracer_conc_new = tracer_conc(c,j) - transp_mass_vr(c,j)/dz(c,j)
             else
               tracer_conc_new  = tracer_conc(c,j) * exp(-max(qflx_rootsoi(c,j)*vtrans_scal/dz(c,j),0._r8)*dtime_loc(c))
               transp_mass_vr(c,j)   = (tracer_conc(c,j)-tracer_conc_new)*dz(c,j)
             endif
             transp_mass(c) = transp_mass(c) + transp_mass_vr(c,j)
             tracer_conc(c,j) = tracer_conc_new
          enddo
       endif
    enddo

  end subroutine calc_root_uptake_as_perfect_sink


  !--------------------------------------------------------------------------------
  subroutine calc_tracer_surface_runoff(betr_time, bounds, lbj, ubj, num_soilc, filter_soilc, &
       fracice_top, dz_top2, biophysforc,  betrtracer_vars,    &
       tracerstate_vars, tracercoeff_vars, tracerflux_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! calculate tracer loss through surface water runoff
    !
    ! !USES:
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use MathfuncMod           , only : safe_div
    use BeTR_TimeMod          , only : betr_time_type
    use BetrStatusType        , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: lbj, ubj
    integer                          , intent(in)    :: num_soilc                               ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)                         ! column filter_soilc
    real(r8)                         , intent(in)    :: fracice_top(bounds%begc:bounds%endc)    ! ice fraction of topsoil
    real(r8)                         , intent(in)    :: dz_top2(bounds%begc:bounds%endc, 1:ubj) ! node depth of the first 2 soil layers
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars                         ! betr configuration information
    type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars                        ! tracer phase conversion coefficients
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars                         ! tracer flux
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars                        ! tracer state variables data structure
    type(betr_status_type)           , intent(out)   :: betr_status
    ! !LOCAL VARIABLES:
    integer  :: fc, c, j, k
    real(r8) :: scal
    real(r8) :: fracc(2)
    real(r8) :: h2o_srun           ! total amount of water lost as surface runoff
    real(r8) :: trc_srun           !
    real(r8) :: dloss
    real(r8) :: dtime
    real(r8) :: total
    real(r8) :: frac1
    real(r8) :: dtmp

    call betr_status%reset()
    ! remove compiler warnings about unused dummy args
    if (lbj > 0) continue

    associate(                                                                  &
         tracernames           => betrtracer_vars%tracernames                 , &
         ngwmobile_tracers     => betrtracer_vars%ngwmobile_tracers           , & !
         groupid               => betrtracer_vars%groupid                     , & !
         is_advective          => betrtracer_vars%is_advective                , & !
         is_h2o                => betrtracer_vars%is_h2o                      , & !Input [logical (:)] indicator whether it is a H2O tracer
         h2osoi_liqvol         => biophysforc%h2osoi_liqvol_col           ,     & !
         qflx_surf             => biophysforc%qflx_surf_col                ,    & !Input [real(r8) (:)]   surface runoff [mm H2O/s]
         tracer_conc_surfwater => tracerstate_vars%tracer_conc_surfwater_col  , & !Inout [real(r8) (:,:)] tracer concentration in surface water
         tracer_conc_mobile    => tracerstate_vars%tracer_conc_mobile_col     , & !
         aqu2bulkcef_mobile    => tracercoeff_vars%aqu2bulkcef_mobile_col     , & !
         tracer_flx_surfrun    => tracerflux_vars%tracer_flx_surfrun_col        & !Output[real(r8) (:,:)] tracer loss through surface runoff
         )

      dtime = betr_time%get_step_size()
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         !it is assumed the surface runoff water mixes perfectly with that of the first two soil nodes, so that a proportion goes off with surface runoff
         !Obtain the total volume
         if(abs(qflx_surf(c))<1.e-50_r8)cycle
         !volume of water coming from surface runoff
         h2o_srun = qflx_surf(c) * dtime / denh2o
         !total volume of water
         total = h2o_srun+ h2osoi_liqvol(c,1) * dz_top2(c,1) + h2osoi_liqvol(c,2) * dz_top2(c,2) * (1._r8-fracice_top(c))
         !fraction lost through liquid water surface runoff
         frac1 = h2o_srun/total

         do j = 1, ngwmobile_tracers

            if(.not. is_advective(j))cycle
            !Do not do this for water tracer, because surface runoff comes as the residual of infiltration.
            if(is_h2o(j))cycle
            tracer_conc_surfwater(c,j)  = 0._r8   !at this moment it is set to zero, however, when tracer is tracked in snow, it should be non-zero
            trc_srun = tracer_conc_surfwater(c,j) * h2o_srun   !total tracer mass in runoff water before mixing
            total = trc_srun
            do k = 1, 2
               if(k==1)then
                  scal = 1._r8
               else
                  scal=1._r8-fracice_top(c)      !reduce the water flush due to ice forst in layer 1
               endif
               fracc(k) = safe_div(tracer_conc_mobile(c,k,j), aqu2bulkcef_mobile(c,k,groupid(j))) * &
                    h2osoi_liqvol(c,k) * dz_top2(c,k) * scal
               total = total + fracc(k)        !total mass
            enddo
            !assume perfect mix and obtain the net loss through surface runoff
            dloss = total * frac1
            !total export through runoff
            tracer_flx_surfrun(c,j) = tracer_flx_surfrun(c,j) + dloss
            !increase of tracer in the surface runoff
            dloss = dloss - trc_srun
            dtmp = fracc(1)+fracc(2)
            tracer_conc_mobile(c,1,j) = tracer_conc_mobile(c,1,j) - dloss*safe_div(fracc(1),dtmp)/dz_top2(c,1)
            tracer_conc_mobile(c,2,j) = tracer_conc_mobile(c,2,j) - dloss*safe_div(fracc(2),dtmp)/dz_top2(c,2)
            tracer_conc_surfwater(c,j) = tracer_flx_surfrun(c,j)/h2o_srun       !revise the tracer concentration in runoff
         enddo
      enddo

    end associate
  end subroutine calc_tracer_surface_runoff


  !--------------------------------------------------------------------------------
  subroutine calc_tracer_h2osfc_snow_residual_combine(betr_time, &
       bounds, col, num_soilc, filter_soilc, biophysforc, betrtracer_vars, &
       tracerstate_vars, tracerflux_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! apply tracer flux from combining residual snow and ponding water
    ! !USES:
    use tracerfluxType  , only : tracerflux_type
    use tracerstatetype , only : tracerstate_type
    use BeTR_TimeMod    , only : betr_time_type
    use BetrStatusType  , only : betr_status_type
    use betr_columnType , only : betr_column_type
    ! !ARGUMENTS:
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds
    type(betr_column_type)           , intent(in)    :: col
    integer                          , intent(in)    :: num_soilc        ! number of column soil points in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)  ! column filter_soilc for soil points
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars  ! betr configuration information
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars  ! tracer flux
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars ! tracer state variables data structure
    type(betr_status_type)           , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    real(r8) :: dtime
    integer :: fc, c, j

    call betr_status%reset()
    ! remove compiler warnings about unused dummy args
    if (bounds%begc > 0) continue

    associate(                                                                                   & !
         dz                 =>    col%dz                                                       , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         qflx_snow2topsoi   =>    biophysforc%qflx_snow2topsoi_col                          ,    & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_h2osfc2topsoi =>    biophysforc%qflx_h2osfc2topsoi_col                        ,    & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s
         tracer_flx_h2osfc_snow_residual => tracerflux_vars%tracer_flx_h2osfc_snow_residual_col, & !
         tracer_conc_mobile =>    tracerstate_vars%tracer_conc_mobile_col                      , & !
         is_h2o             =>    betrtracer_vars%is_h2o                                       , & !
         tracernames        =>    betrtracer_vars%tracernames                                  , & !
         ngwmobile_tracers  =>    betrtracer_vars%ngwmobile_tracers                              & !
         )

      dtime = betr_time%get_step_size()
      do j = 1, ngwmobile_tracers

         if(.not. is_h2o(j))cycle
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            ! kg/m2/(kg/m3) = m
            tracer_flx_h2osfc_snow_residual(c,j) = (qflx_snow2topsoi(c) + qflx_h2osfc2topsoi(c))*dtime
            tracer_conc_mobile(c,1,j) = tracer_conc_mobile(c,1,j) + tracer_flx_h2osfc_snow_residual(c,j) /dz(c,1)
         enddo
      enddo

    end associate
  end subroutine calc_tracer_h2osfc_snow_residual_combine

  !--------------------------------------------------------------------------------
  subroutine diagnose_gas_pressure(bounds, lbj, ubj, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! diagnose gas pressure

    ! !USES:
    use tracercoeffType       , only : tracercoeff_type
    use tracerfluxType        , only : tracerflux_type
    use tracerstatetype       , only : tracerstate_type
    use MathfuncMod           , only : safe_div
    use betr_varcon           , only : oneatm => boneatm
    use BetrStatusType        , only : betr_status_type
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: lbj, ubj
    integer                , intent(in)    :: num_soilc        ! number of columns in column filter_soilc
    integer                , intent(in)    :: filter_soilc(:)  ! column filter_soilc
    class(betrtracer_type) , intent(in)    :: betrtracer_vars  ! tracer info data structure
    type(tracercoeff_type) , intent(in)    :: tracercoeff_vars ! tracer phase conversion coefficients
    type(tracerstate_type) , intent(inout) :: tracerstate_vars ! tracer state variables data structure
    type(betr_status_type) , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    integer  :: j, fc, c, jj
    real(r8) :: total_pres

    call betr_status%reset()
    ! remove compiler warnings about unused dummy args
    if (bounds%begc > 0) continue
    if (lbj > 0) continue

  associate(                                                                  &
         tracer_conc_mobile_col   => tracerstate_vars%tracer_conc_mobile_col, &
         tracer_P_gas_frac_col    => tracerstate_vars%tracer_P_gas_frac_col , &
         tracer_P_gas_col         => tracerstate_vars%tracer_P_gas_col      , &
         aqu2bulkcef_mobile_col   => tracercoeff_vars%aqu2bulkcef_mobile_col, &
         henrycef_col             => tracercoeff_vars%henrycef_col          , &
         volatilegroupid          => betrtracer_vars%volatilegroupid        , &
         volatileid               => betrtracer_vars%volatileid             , &
         groupid                  => betrtracer_vars%groupid                , &
         ngwmobile_tracers        => betrtracer_vars%ngwmobile_tracers      , &
         is_volatile              => betrtracer_vars%is_volatile            , &
         is_isotope               => betrtracer_vars%is_isotope             , &
         is_h2o                   => betrtracer_vars%is_h2o                   &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         do j = 1, ubj
            !calculate the total gas pressure
            total_pres=0._r8
            do jj = 1, ngwmobile_tracers
               if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
                  tracer_P_gas_frac_col(c,j, volatileid(jj))  = calc_gas_pressure(tracer_conc_mobile_col(c,j,jj), &
                       aqu2bulkcef_mobile_col(c,j,groupid(jj)), henrycef_col(c, j, volatilegroupid(jj)))
                  total_pres=total_pres + tracer_P_gas_frac_col(c,j, volatileid(jj))
               endif
            enddo
            tracer_P_gas_col(c,j) = total_pres * oneatm
            do jj = 1, ngwmobile_tracers
               if(is_volatile(jj) .and. (.not. is_h2o(jj)) .and. (.not. is_isotope(jj)))then
                  tracer_P_gas_frac_col(c,j, volatileid(jj)) = safe_div(tracer_P_gas_frac_col(c,j, volatileid(jj)), &
                       total_pres)
               endif
            enddo
         enddo
      enddo
    end associate
  end subroutine diagnose_gas_pressure


   !------------------------------------------------------------------------

   subroutine tracer_copy_a2b_div(bounds, num_snowc, filter_snowc, snl, tracer_in, tracer_copy)
   implicit none
   type(bounds_type) , intent(in)  :: bounds               ! bounds
   integer           , intent(in)  :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)  :: filter_snowc(:) ! column filter for soil points
   integer           , intent(in)  :: snl(bounds%begc:bounds%endc)
   real(r8)          , intent(in)  :: tracer_in(bounds%begc:bounds%endc, -nlevsno+1:0)
   real(r8)          , intent(out) :: tracer_copy(bounds%begc:bounds%endc, 1:nlevsno)

   integer :: fc, c, j

   tracer_copy(:,:) = 0._r8
   do j = 1, nlevsno
     do fc = 1, num_snowc
       c = filter_snowc(fc)
       if(j <= abs(snl(c)))then
         tracer_copy(c,j) = tracer_in(c,j+snl(c))
       endif
     enddo
   enddo

   end subroutine tracer_copy_a2b_div
   !------------------------------------------------------------------------

   subroutine tracer_copy_a2b_comb(bounds, num_snowc, filter_snowc, snl, tracer_in, tracer_copy)
   implicit none
   type(bounds_type) , intent(in)  :: bounds               ! bounds
   integer           , intent(in)  :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)  :: filter_snowc(:) ! column filter for soil points
   integer           , intent(in)  :: snl(bounds%begc:bounds%endc)
   real(r8)          , intent(in)  :: tracer_in(bounds%begc:bounds%endc, -nlevsno+1:1)
   real(r8)          , intent(out) :: tracer_copy(bounds%begc:bounds%endc,-nlevsno+1:1)

   integer :: fc, c, j

   tracer_copy(:,:) = 0._r8
   do j = -nlevsno+1, 1
     do fc = 1, num_snowc
       c = filter_snowc(fc)
       if(j >= snl(c)+1)then
         tracer_copy(c,j) = tracer_in(c,j)
       endif
     enddo
   enddo
   end subroutine tracer_copy_a2b_comb
   !------------------------------------------------------------------------

   subroutine tracer_col_mapping_div(bounds, num_snowc, filter_snowc, snl, transition_matrix, &
      tracer_copy, tracer_out, bstatus)

   use BetrStatusType, only : betr_status_type
   implicit none
   type(bounds_type) , intent(in)  :: bounds               ! bounds
   integer           , intent(in)  :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)  :: filter_snowc(:) ! column filter for soil points
   integer           , intent(in)  :: snl(bounds%begc:bounds%endc)
   real(r8)          , intent(in)  :: transition_matrix(bounds%begc:bounds%endc,1:nlevsno,1:nlevsno)
   real(r8)          , intent(in)  :: tracer_copy(bounds%begc:bounds%endc, 1:nlevsno)
   real(r8)          , intent(out) :: tracer_out(bounds%begc:bounds%endc, -nlevsno+1:0)
   type(betr_status_type), intent(out) :: bstatus
   integer :: fc, c, j

   call bstatus%reset()
   do j = -nlevsno + 1, 0
     do fc = 1, num_snowc
       c = filter_snowc(fc)
       if(j>= snl(c)+1)then
         tracer_out(c,j) = dot_sum(transition_matrix(c,j-snl(c),:),tracer_copy(c,:),bstatus)
         if(bstatus%check_status())return
       endif
     enddo
   enddo
   end subroutine tracer_col_mapping_div
   !------------------------------------------------------------------------

   subroutine tracer_col_mapping_comb(bounds, num_snowc, filter_snowc, snl, transition_matrix, &
      tracer_copy, tracer_out, bstatus)

    use BetrStatusType, only : betr_status_type
   implicit none
   type(bounds_type) , intent(in)  :: bounds               ! bounds
   integer           , intent(in)  :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)  :: filter_snowc(:) ! column filter for soil points
   integer           , intent(in)  :: snl(bounds%begc:bounds%endc)
   real(r8)          , intent(in)  :: transition_matrix(bounds%begc:bounds%endc,-nlevsno+1:1,-nlevsno+1:1)
   real(r8)          , intent(in)  :: tracer_copy(bounds%begc:bounds%endc,-nlevsno+1:1)
   real(r8)          , intent(out) :: tracer_out(bounds%begc:bounds%endc, -nlevsno+1:1)
   type(betr_status_type), intent(out) :: bstatus
   integer :: fc, c, j

   call bstatus%reset()
   do j = -nlevsno + 1, 1
     do fc = 1, num_snowc
       c = filter_snowc(fc)
       if(j>= snl(c)+1)then
         tracer_out(c,j) = dot_sum(transition_matrix(c,j,:),tracer_copy(c,:), bstatus)
         if(bstatus%check_status())return
       endif
     enddo
   enddo
   end subroutine tracer_col_mapping_comb


end module BetrBGCMod
