module TracerBalanceMod

!
! !DESCRIPTION:
! module contains subroutines to do
! tracer mass balance check

  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use BeTRTracerType     , only : betrtracer_type
  use ColumnType         , only : col
  use clm_time_manager   , only : get_nstep
  use clm_varctl         , only : iulog
implicit none
  save
  private


  public :: begin_betr_tracer_massbalance
  public :: betr_tracer_massbalance_check

  contains


    
    !--------------------------------------------------------------------------------
    subroutine begin_betr_tracer_massbalance(bounds, lbj, ubj, numf, filter, &
         betrtracer_vars, tracerstate_vars, tracerflux_vars)
      !
      ! !DESCRIPTION:
      ! Preparing for tracer mass balance check
      !
      ! !USES:
      use tracerstatetype       , only : tracerstate_type
      use clm_varpar            , only : nlevtrc_soil
      use tracerfluxType        , only : tracerflux_type

      implicit none
      ! !ARGUMENTS:
      type(bounds_type),      intent(in)    :: bounds
      integer,                intent(in)    :: lbj, ubj
      integer,                intent(in)    :: numf                        ! number of columns in column filter
      integer,                intent(in)    :: filter(:)                   ! column filter
      type(betrtracer_type),  intent(in)    :: betrtracer_vars             ! betr configuration information
      type(tracerstate_type), intent(inout) :: tracerstate_vars            ! tracer state variables data structure
      type(tracerflux_type) , intent(inout) :: tracerflux_vars

      ! !LOCAL VARIABLES:
      character(len=256) :: subname='begin_betr_tracer_massbalance'
      integer :: fc, c

      call tracerflux_vars%Reset(bounds, numf, filter)

      call betr_tracer_mass_summary(bounds, lbj, ubj, numf, filter, betrtracer_vars, tracerstate_vars, &
           tracerstate_vars%beg_tracer_molarmass_col(bounds%begc:bounds%endc, 1:betrtracer_vars%ntracers))

    end subroutine begin_betr_tracer_massbalance

    !--------------------------------------------------------------------------------
    subroutine betr_tracer_massbalance_check(bounds, lbj, ubj, numf, filter, betrtracer_vars, tracerstate_vars,  tracerflux_vars)
      !
      ! !DESCRIPTION:
      ! do mass balance check for betr tracers
      !
      ! for solid phase tracers, the only source/sink is biogeochemical production/consumption
      ! and it is currently assumed no solid phase input from atmospheric precipitation (either dry or wet)
      ! the equilibrium fraction is always associated with the (dual)-phase mobile tracer.
      ! However the situation is different for water isotopes, because ice is also part of the
      ! mass budget, and by assuming equilibrium partitioning, the chemical source/sink for ice is not tracked explicitly.
      !
      ! !USES:
      use tracerfluxType        , only : tracerflux_type
      use tracerstatetype       , only : tracerstate_type
      use abortutils            , only : endrun
      use clm_varctl            , only : iulog
      use clm_time_manager      , only : get_step_size,get_nstep
      use clm_varcon            , only : namec,catomw,natomw
      implicit none

      ! !ARGUMENTS:
      type(bounds_type),      intent(in)    :: bounds
      integer,                intent(in)    :: lbj, ubj
      integer,                intent(in)    :: numf             ! number of columns in column filter
      integer,                intent(in)    :: filter(:)        ! column filter
      type(betrtracer_type),  intent(in)    :: betrtracer_vars  ! betr configuration information
      type(tracerflux_type),  intent(inout) :: tracerflux_vars
      type(tracerstate_type), intent(inout) :: tracerstate_vars ! tracer state variables data structure
      ! !LOCAL VARIABLES:
      integer  :: jj, fc, c, kk
      real(r8) :: dtime
      real(r8) :: atw
      real(r8) :: err_rel, bal_beg, bal_end, bal_flx
      real(r8), parameter :: err_min = 1.e-8_r8
      real(r8), parameter :: err_min_rel=1.e-3_r8
      associate(                                                                            &
           beg_tracer_molarmass      => tracerstate_vars%beg_tracer_molarmass_col         , &
           end_tracer_molarmass      => tracerstate_vars%end_tracer_molarmass_col         , &
           tracer_flx_infl           => tracerflux_vars%tracer_flx_infl_col               , &
           tracer_flx_netpro         => tracerflux_vars%tracer_flx_netpro_col             , &
           tracer_flx_netphyloss     => tracerflux_vars%tracer_flx_netphyloss_col         , &
           is_mobile                 => betrtracer_vars%is_mobile                         , &
           errtracer                 => tracerstate_vars%errtracer_col                    , &
           ngwmobile_tracers         => betrtracer_vars%ngwmobile_tracers                 , &
           tracernames               => betrtracer_vars%tracernames                       , &
           ntracers                  => betrtracer_vars%ntracers                            &
           
           )

        call betr_tracer_mass_summary(bounds, lbj, ubj, numf, filter, betrtracer_vars, tracerstate_vars, &
             end_tracer_molarmass(bounds%begc:bounds%endc, 1:betrtracer_vars%ntracers))

        dtime=get_step_size()

        do fc = 1, numf
           c = filter(fc)
           !summarize the fluxes
           call tracerflux_vars%flux_summary(c, betrtracer_vars)

           do kk = 1, ngwmobile_tracers
              errtracer(c,kk) = beg_tracer_molarmass(c,kk)-end_tracer_molarmass(c,kk)  
              errtracer(c,kk) = errtracer(c,kk) + tracer_flx_netpro(c,kk)-tracer_flx_netphyloss(c,kk)
              if(abs(errtracer(c,kk))<err_min)then
                 err_rel=1.e-4_r8
              else
                 err_rel = errtracer(c,kk)/max(abs(beg_tracer_molarmass(c,kk)),abs(end_tracer_molarmass(c,kk)))
              endif

              if(abs(err_rel)>err_min_rel)then
                 write(iulog,*)'error exceeds the tolerance for tracer '//tracernames(kk), ' err=',errtracer(c,kk), ' col=',c
                 write(iulog,*)'nstep=',get_nstep()
                 write(iulog,'(4(A,X,E20.10))')'netpro=',tracer_flx_netpro(c,kk),' netphyloss=',tracer_flx_netphyloss(c,kk),&
                      ' begm=',beg_tracer_molarmass(c,kk), &
                      ' endm=',end_tracer_molarmass(c,kk)
                 call tracerflux_vars%flux_display(c,kk,betrtracer_vars)
                 call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
              endif
           enddo
           bal_beg=0._r8
           bal_end=0._r8
           bal_flx=0._r8
           do kk = ngwmobile_tracers+1, ntracers
              errtracer(c,kk) = beg_tracer_molarmass(c,kk)-end_tracer_molarmass(c,kk) + tracer_flx_netpro(c,kk)
              if(abs(errtracer(c,kk))>err_min)then
                 write(iulog,*)'error exceeds the tolerance for tracer '//tracernames(kk), 'err=',errtracer(c,kk), 'col=',c
                 write(iulog,*)get_nstep(),is_mobile(kk)
                 write(iulog,*)'begmss=', beg_tracer_molarmass(c,kk), 'endmass=',end_tracer_molarmass(c,kk),' netpro=',tracer_flx_netpro(c,kk)
                 call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
              endif
           enddo

           call tracerflux_vars%Temporal_average(c,dtime)
        enddo

      end associate

    end subroutine betr_tracer_massbalance_check


    !--------------------------------------------------------------------------------
    
    subroutine betr_tracer_mass_summary(bounds, lbj, ubj, numf, filter, betrtracer_vars, tracerstate_vars, tracer_molarmass_col)
      !
      ! !DESCRIPTION:
      ! summarize the column tracer mass
      !
      ! !USES:
      use tracerstatetype       , only : tracerstate_type
      use clm_varpar            , only : nlevtrc_soil
      use MathfuncMod           , only : dot_sum

      implicit none
      ! !ARGUMENTS:
      type(bounds_type),      intent(in)    :: bounds
      integer,                intent(in)    :: lbj, ubj
      integer,                intent(in)    :: numf                        ! number of columns in column filter
      integer,                intent(in)    :: filter(:)                   ! column filter
      type(betrtracer_type) , intent(in)    :: betrtracer_vars             ! betr configuration information
      type(tracerstate_type), intent(inout) :: tracerstate_vars            ! tracer state variables data structure
      real(r8)              , intent(inout) :: tracer_molarmass_col(bounds%begc:bounds%endc, 1:betrtracer_vars%ntracers)
      ! !LOCAL VARIABLES:
      integer :: jj, fc, c, kk

      associate(                                                                            &
           tracer_conc_mobile        => tracerstate_vars%tracer_conc_mobile_col           , &
           tracer_conc_solid_equil   => tracerstate_vars%tracer_conc_solid_equil_col      , &
           tracer_conc_solid_passive => tracerstate_vars%tracer_conc_solid_passive_col    , &
           dz                        => col%dz                                            , &
           ngwmobile_tracers         => betrtracer_vars%ngwmobile_tracers                 , &
           ntracers                  => betrtracer_vars%ntracers                          , &
           is_adsorb                 => betrtracer_vars%is_adsorb                         , &
           nsolid_passive_tracers    => betrtracer_vars%nsolid_passive_tracers            , &
           adsorbid                  => betrtracer_vars%adsorbid                            &
           
           )

        do jj = 1,   ngwmobile_tracers
           do fc = 1, numf
              c = filter(fc)

              tracer_molarmass_col(c,jj) = dot_sum(tracer_conc_mobile(c,1:nlevtrc_soil,jj), dz(c,1:nlevtrc_soil))

              if(is_adsorb(jj))then
                 tracer_molarmass_col(c,jj) = tracer_molarmass_col(c,jj) + &
                      dot_sum(tracer_conc_solid_equil(c,1:nlevtrc_soil,adsorbid(jj)),dz(c,1:nlevtrc_soil))
              endif
           enddo
        enddo

        do jj = 1, nsolid_passive_tracers
           kk = jj + ngwmobile_tracers
           do fc = 1, numf
              c = filter(fc)
              tracer_molarmass_col(c,kk) = dot_sum(tracer_conc_solid_passive(c,1:nlevtrc_soil,jj), dz(c,1:nlevtrc_soil))
           enddo
        enddo
      end associate
    end subroutine betr_tracer_mass_summary

end module TracerBalanceMod
