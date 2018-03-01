module TracerFluxType
  !!DESCRIPTION:
  ! tracer flux type
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod      , only : bounds_type  => betr_bounds_type
  use betr_varcon         , only : spval => bspval, ispval => bispval
  use tracer_varcon       , only : nlevtrc_soil => betr_nlevtrc_soil
  use BeTR_landvarconType , only : landvarcon => betr_landvarcon
  use betr_ctrl           , only : iulog => biulog
  use TracerBaseType      , only : tracerbase_type
  !
  ! !PUBLIC TYPES:

  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  !
  ! !PUBLIC DATA:
  !

  type, public, extends(tracerbase_type) :: TracerFlux_type

     !tracer flux defined at the column level
     real(r8), pointer :: tracer_flx_top_soil_col(:,:)  => null()  !tracer fluxes available for infiltration+runoff
     real(r8), pointer :: tracer_flx_can_loss_col(:,:)  => null()  !tracer loss from canopy
     real(r8), pointer :: tracer_flx_snowmelt_col(:,:) => null()   !tracer loss from snow melting
     real(r8), pointer :: tracer_flx_infl_col(:,:)     => null()   !tracer fluxes available for infiltration
     real(r8), pointer :: tracer_flx_netphyloss_col(:,:) => null() !total tracer loos through all possible physical pathways: drainage (+ runoff), leaching, ebullition, diffusion, minus precipitation/infiltration
     real(r8), pointer :: tracer_flx_netpro_col(:,:)   => null()   !total tracer production through chemical processes
     real(r8), pointer :: tracer_flx_dstor_col(:,:)    => null()   !net storage of tracer due to input-output, ideally, dstor=netpro-netloss at various scales
     real(r8), pointer :: tracer_flx_ebu_col(:,:)     => null()    !tracer emitted as bubbles, mol, lake, volatile
     real(r8), pointer :: tracer_flx_prec_col(:,:)   => null()     !tracer added to a column from precipitation, mol
     real(r8), pointer :: tracer_flx_dif_col(:,:)    => null()     !tracer emitted through diffusion, unsat, volatile

     real(r8), pointer :: tracer_flx_drain_col(:,:)   => null()    !tracer removal through subface drainage
     real(r8), pointer :: tracer_flx_surfemi_col(:,:)  => null()   !total emitted tracer fluxes at surface, volatile, including ebullition, diffusion, arenchyma transport
     real(r8), pointer :: tracer_flx_leaching_col(:,:)  => null()  !leaching fluxes
     real(r8), pointer :: tracer_flx_surfrun_col(:,:)   => null()  !tracer loss thru runoff, mol tracer / second
     real(r8), pointer :: tracer_flx_netpro_vr_col(:,:,:) => null()!total source strength for the tracers, chemical production, root exudation, excludes incoming root transport (by exchange with air) and (infiltration?)
     real(r8), pointer :: tracer_flx_tparchm_col(:,:)    => null() !total tracer flux through plant aerenchyma transport, for volatile species only, mol/m^2/s
     real(r8), pointer :: tracer_flx_parchm_vr_col(:,:,:) => null()!vertical resolved tracer flux through aerenchyma transport, for volatile species only, mol/m^3/s
     real(r8), pointer :: tracer_flx_totleached_col(:,:) => null() !total leaching flux, vertical + lateral leaching

     real(r8), pointer :: tracer_flx_vtrans_col(:,:)     => null() !column level tracer flux through transpiration
     real(r8), pointer :: tracer_flx_vtrans_vr_col(:,:,:) => null()!
     !real(r8), pointer :: tracer_flx_snowloss_col(:,:)  => null()  !tracer flux lost from snow dynamics, place holder

     !tracer fluxes defined at the pft level
     real(r8), pointer :: tracer_flx_vtrans_patch(:,:)      => null()       !tracer goes to the pathway of plant transpiration, currently not released, if it is nutrient, assumed it is taken by plants completely
     real(r8), pointer :: tracer_flx_snowfall_grnd_patch(:,:)=> null()
     real(r8), pointer :: tracer_flx_rainfall_grnd_patch(:,:)=> null()
     real(r8), pointer :: tracer_flx_prec_intr_patch(:,:) => null()         !interception of tracer from wet deposition [mol/s]
     real(r8), pointer :: tracer_flx_prec_grnd_patch(:,:) => null()         !tracer onto ground including from canopy runoff [mol /s]
     real(r8), pointer :: tracer_flx_snwcp_liq_patch(:,:)  => null()        !excess rainfall tracer due to snow capping [mol /s]
     real(r8), pointer :: tracer_flx_snwcp_ice_patch(:,:)  => null()        !excess snowfall tracer due to snow capping [mol /s], this is used for aerosol type and water type tracer input
     real(r8), pointer :: tracer_flx_dew_grnd_col(:,:)     => null()        !tracer flux to ground coming from dew formation
     real(r8), pointer :: tracer_flx_dew_snow_col(:,:)     => null()        !tracer flux to snow coming from dew formation
     real(r8), pointer :: tracer_flx_sub_snow_col(:,:)      => null()       !tracer flux loss from snow sublimation
     real(r8), pointer :: tracer_flx_h2osfc_snow_residual_col(:,:) => null() !tracer flux coming from residual standing water and residual snow

   contains
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: Reset
     procedure, public  :: Temporal_average
     procedure, public  :: Flux_summary
     procedure, public  :: Flux_display
     procedure, public  :: retrieve_hist
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type TracerFlux_type
contains
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, lbj, ubj, betrtracer_vars)
    ! !DESCRIPTION:
    ! initialize data type
    !
    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    implicit none
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout)  :: this
    type(bounds_type)     , intent(in) :: bounds
    integer               , intent(in) :: lbj, ubj
    type(BeTRTracer_Type) , intent(in) :: betrtracer_vars

    call this%InitAllocate(bounds, lbj, ubj, betrtracer_vars)
    call this%tracer_base_init()
    call this%InitHistory (bounds, betrtracer_vars)
    call this%InitCold    (bounds)
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, lbj, ubj, betrtracer_vars)
    ! !DESCRIPTION:
    ! memory allocation
    !
    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    implicit none
    !
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout) :: this
    type(bounds_type)     , intent(in) :: bounds
    integer               , intent(in) :: lbj, ubj
    type(BeTRTracer_Type) , intent(in) :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    integer :: begp, endp
    integer :: ngwmobile_tracers
    integer :: nvolatile_tracers
    integer :: ntracers
    integer :: nsolid_passive_tracers
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp
    ngwmobile_tracers      = betrtracer_vars%ngwmobile_tracers
    ntracers               = betrtracer_vars%ntracers
    nvolatile_tracers      = betrtracer_vars%nvolatile_tracers
    nsolid_passive_tracers = betrtracer_vars%nsolid_passive_tracers

    allocate(this%tracer_flx_prec_col             (begc:endc, 1:ntracers)); this%tracer_flx_prec_col           (:,:) = spval
    allocate(this%tracer_flx_snowfall_grnd_patch  (begp:endp, 1:ntracers)); this%tracer_flx_snowfall_grnd_patch(:,:) = spval
    allocate(this%tracer_flx_rainfall_grnd_patch  (begp:endp, 1:ntracers)); this%tracer_flx_rainfall_grnd_patch(:,:) = spval
    allocate(this%tracer_flx_prec_intr_patch      (begp:endp, 1:ntracers)); this%tracer_flx_prec_intr_patch    (:,:) = spval
    allocate(this%tracer_flx_prec_grnd_patch(begp:endp, 1:ntracers)); this%tracer_flx_prec_grnd_patch(:,:) = spval
    allocate(this%tracer_flx_snwcp_liq_patch      (begp:endp, 1:ntracers)); this%tracer_flx_snwcp_liq_patch    (:,:) = spval
    allocate(this%tracer_flx_snwcp_ice_patch      (begp:endp, 1:ntracers)); this%tracer_flx_snwcp_ice_patch    (:,:) = spval

    if(ngwmobile_tracers>0)then
       allocate(this%tracer_flx_drain_col       (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_drain_col    (:,:) = spval
       allocate(this%tracer_flx_top_soil_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_top_soil_col (:,:) = spval
       allocate(this%tracer_flx_can_loss_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_can_loss_col (:,:) = spval
       allocate(this%tracer_flx_snowmelt_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_snowmelt_col (:,:) = spval
       allocate(this%tracer_flx_infl_col        (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_infl_col     (:,:) = spval
       allocate(this%tracer_flx_leaching_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_leaching_col (:,:) = spval
       allocate(this%tracer_flx_surfrun_col     (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_surfrun_col  (:,:) = spval
       allocate(this%tracer_flx_vtrans_col      (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_vtrans_col   (:,:) = spval
       allocate(this%tracer_flx_dew_grnd_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_dew_grnd_col (:,:) = spval
       allocate(this%tracer_flx_dew_snow_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_dew_snow_col (:,:) = spval
       allocate(this%tracer_flx_sub_snow_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_sub_snow_col (:,:) = spval
       allocate(this%tracer_flx_vtrans_patch    (begp:endp, 1:ngwmobile_tracers)); this%tracer_flx_vtrans_patch(:,:) = spval

       allocate(this%tracer_flx_h2osfc_snow_residual_col(begc:endc, 1:ngwmobile_tracers))
       this%tracer_flx_h2osfc_snow_residual_col(:,:) = spval
       allocate(this%tracer_flx_totleached_col(begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_totleached_col(:,:) = spval
       allocate(this%tracer_flx_vtrans_vr_col   (begc:endc, lbj:ubj, 1:ngwmobile_tracers))
       this%tracer_flx_vtrans_vr_col   (:,:,:) = spval
    endif
    if(nvolatile_tracers>0)then
       allocate(this%tracer_flx_ebu_col         (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_ebu_col      (:,:) = spval
       allocate(this%tracer_flx_dif_col         (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_dif_col      (:,:) = spval
       allocate(this%tracer_flx_tparchm_col     (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_tparchm_col  (:,:) = spval
       allocate(this%tracer_flx_surfemi_col     (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_surfemi_col  (:,:) = spval
       allocate(this%tracer_flx_parchm_vr_col   (begc:endc, lbj:ubj, 1:nvolatile_tracers))
       this%tracer_flx_parchm_vr_col(:,:,:) = spval
    endif

    allocate(this%tracer_flx_netpro_vr_col  (begc:endc, lbj:ubj, 1:ntracers)); this%tracer_flx_netpro_vr_col (:,:,:) = spval
    allocate(this%tracer_flx_netphyloss_col (begc:endc, 1:ntracers)); this%tracer_flx_netphyloss_col(:,:)            = spval
    allocate(this%tracer_flx_netpro_col     (begc:endc, 1:ntracers)); this%tracer_flx_netpro_col(:,:)                = spval
    allocate(this%tracer_flx_dstor_col      (begc:endc, 1:ntracers)); this%tracer_flx_dstor_col(:,:)                 = spval

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! History fields initialization
    !
    ! !USES:
    !use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use BeTRTracerType, only: BeTRTracer_Type
    !
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout)  :: this
    type(bounds_type)     , intent(in) :: bounds
    type(BeTRTracer_Type) , intent(in) :: betrtracer_vars

    !
    ! !LOCAL VARIABLES:
    integer :: ntracers
    integer :: ngwmobile_tracers
    integer :: nsolid_passive_tracers
    integer :: jj, kk
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays
    integer :: num1d, num2d, it
    associate(                                                   &
      ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers   , &
      ntracers          => betrtracer_vars%ntracers            , &
      is_volatile       => betrtracer_vars%is_volatile         , &
      volatileid        => betrtracer_vars%volatileid          , &
      tracernames       => betrtracer_vars%tracernames           &
    )

    num2d = 0; num1d= 0
    do it = 1, 2
      do jj = 1, ntracers
        if(jj<= ngwmobile_tracers) then

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_DEW_GRND', units='mol/m2/s', &
           avgflag='A', long_name='incoming dew flux to ground for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_DEW_SNOW', units='mol/m2/s', &
           avgflag='A', long_name='incoming dew flux to snow from '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_H2OSFC_SNOW_RES', &
           units='mol/m2/s', avgflag='A', &
           long_name='incoming flux to topsoi from snow and h2osfc residual for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_SUB_SNOW', units='mol/m2/s',   &
           avgflag='A', long_name='sublimation flux from snow for '//trim(tracernames(jj)),      &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_TOPSOIL', units='mol/m2/s',    &
           avgflag='A', long_name='incoming flux at top of the soil for '//trim(tracernames(jj)),   &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_CAN_LOSS', units='mol/m2/s',     &
           avgflag='A', long_name='loss from canopy for '//trim(tracernames(jj)),       &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_SNOWMELT', units='mol/m2/s',     &
           avgflag='A', long_name='loss from snowmelt for '//trim(tracernames(jj)),      &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_INFIL', units='mol/m2/s',     &
           avgflag='A', long_name='infiltration for '//trim(tracernames(jj)),            &
           default='inactive')

          call this%add_hist_var2d (it, num2d, fname=trim(tracernames(jj))//'_FLX_NETPRO_vr',&
           units='mol/m3/s', type2d='levtrc', avgflag='A', &
           long_name='net production for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_LEACHING', units='mol/m2/s', &
           avgflag='A', long_name='bottom of soil leaching for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_SRUNOFF', units='mol/m2/s', &
           avgflag='A', long_name='loss from surface runoff for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_VTRANS', units='mol/m2/s', &
           avgflag='A', long_name='transport through transpiration for '//trim(tracernames(jj)), &
           default='inactive')

          call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_TLEACH', units='mol/m2/s',   &
           avgflag='A', long_name='transport through leaching for '//trim(tracernames(jj)),    &
           default='inactive')

          if(is_volatile(jj))then
            kk = volatileid(jj)

            call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_EBU', &
              units='mol/m2/s',  avgflag='A', &
              long_name='loss through ebullition (+ into atmosphere) for '//trim(tracernames(jj)), &
              default='inactive')

            call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_DIF', &
              units='mol/m2/s',  avgflag='A', &
              long_name='loss through diffusion (+ into atmosphere) for '//trim(tracernames(jj)),  &
              default='inactive')

            call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_ARCHM', units='mol/m2/s',    &
             avgflag='A', long_name='loss from aerenchyma transport for '//trim(tracernames(jj)),    &
             default='inactive')

            call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_SURFEMI', units='mol/m2/s',    &
             avgflag='A', long_name='loss from surface emission for '//trim(tracernames(jj)))
           endif

           call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_DRAIN', units='mol/m2/s', &
            avgflag='A', long_name='loss from drainage for '//trim(tracernames(jj)),        &
            default='inactive')
        endif

        call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_NETLOSS', units='mol/m2/s',  &
          avgflag='A', long_name='net loss for '//trim(tracernames(jj)),                    &
          default='inactive')

        call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_NETPRO', units='mol/m2/s',   &
          avgflag='A', long_name='net production for '//trim(tracernames(jj)),              &
          default='inactive')

        call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_DSTOR', units='mol/m2/s',    &
          avgflag='A', long_name='total concentration change for '//trim(tracernames(jj)),  &
          default='inactive')

        call this%add_hist_var1d (it, num1d, fname=trim(tracernames(jj))//'_FLX_PREC', units='mol/m2/s',     &
          avgflag='A', long_name='incoming from precipitation for '//trim(tracernames(jj)), &
          default='inactive')

      enddo
      if(it==1)call this%alloc_hist_list(num1d, num2d)
      num2d = 0; num1d= 0
    enddo
    end associate

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! cold initialization
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: c, j, p, l       ! index
    integer :: begp, endp
    integer :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    do p = begp,endp
      this%tracer_flx_vtrans_patch(p,:)         = 0._r8
      this%tracer_flx_snowfall_grnd_patch(p,:)  = 0._r8
      this%tracer_flx_rainfall_grnd_patch(p,:)  = 0._r8
      this%tracer_flx_prec_intr_patch(p,:)      = 0._r8
      this%tracer_flx_prec_grnd_patch(p,:)      = 0._r8
      this%tracer_flx_snwcp_liq_patch(p,:)      = 0._r8
      this%tracer_flx_snwcp_ice_patch(p,:)      = 0._r8
    enddo

    do c = begc, endc
      this%tracer_flx_top_soil_col(c,:)    = 0._r8
      this%tracer_flx_can_loss_col(c,:)    = 0._r8
      this%tracer_flx_snowmelt_col (c,:)   = 0._r8
      this%tracer_flx_infl_col(c,:)        = 0._r8
      this%tracer_flx_netphyloss_col(c,:)  = 0._r8
      this%tracer_flx_netpro_col(c,:)      = 0._r8
      this%tracer_flx_dstor_col(c,:)       = 0._r8
      this%tracer_flx_ebu_col(c,:)         = 0._r8
      this%tracer_flx_prec_col(c,:)        = 0._r8
      this%tracer_flx_dif_col(c,:)         = 0._r8
      this%tracer_flx_drain_col(c,:)       = 0._r8
      this%tracer_flx_surfemi_col(c,:)     = 0._r8
      this%tracer_flx_leaching_col(c,:)    = 0._r8
      this%tracer_flx_surfrun_col(c,:)     = 0._r8
      this%tracer_flx_tparchm_col(c,:)     = 0._r8
      this%tracer_flx_parchm_vr_col(c,:,:) = 0._r8
      this%tracer_flx_vtrans_col(c,:)      = 0._r8
      this%tracer_flx_dew_grnd_col   (c,:) = 0._r8
      this%tracer_flx_dew_snow_col   (c,:) = 0._r8
      this%tracer_flx_sub_snow_col   (c,:) = 0._r8
      this%tracer_flx_h2osfc_snow_residual_col(c,:) = 0._r8
      this%tracer_flx_totleached_col (c,:) = 0._r8
    enddo


  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, flag, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! Now it is purposely empty, but will be potentially useful in the future
    ! !USES:
    use BetrTracerType , only : betrtracer_type
    !
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout) :: this
    type(bounds_type)     , intent(in)    :: bounds
    character(len=*)      , intent(in)    :: flag                                         ! 'read' or 'write'
    type(BeTRTracer_Type) , intent(in)    :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file

    ! remove compiler warnings for unused dummy args
    if (size(this%tracer_flx_top_soil_col) > 0) continue
    if (bounds%begc > 0)                        continue
    if (len(flag) > 0)                          continue
    if (len(betrtracer_vars%betr_simname) > 0)  continue


   end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, bounds, numf, filter)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(TracerFlux_type), intent(inout) :: this
    type(bounds_type)    , intent(in)  :: bounds
    integer              , intent(in)  :: numf
    integer              , intent(in)  :: filter(:)
    !-----------------------------------------------------------------------

    integer :: fc, column

    ! remove compiler warnings for unused dummy args
    if (bounds%begc > 0) continue

    do fc = 1, numf
      column = filter(fc)
      this%tracer_flx_top_soil_col   (column,:)   = 0._r8
      this%tracer_flx_can_loss_col   (column,:)   = 0._r8
      this%tracer_flx_snowmelt_col   (column,:)   = 0._r8
      this%tracer_flx_infl_col       (column,:)   = 0._r8
      this%tracer_flx_netphyloss_col (column,:)   = 0._r8
      this%tracer_flx_netpro_col     (column,:)   = 0._r8
      this%tracer_flx_dstor_col      (column,:)   = 0._r8
      this%tracer_flx_ebu_col        (column,:)   = 0._r8
      this%tracer_flx_prec_col       (column,:)   = 0._r8
      this%tracer_flx_dif_col        (column,:)   = 0._r8
      this%tracer_flx_drain_col      (column,:)   = 0._r8
      this%tracer_flx_surfemi_col    (column,:)   = 0._r8
      this%tracer_flx_leaching_col   (column,:)   = 0._r8
      this%tracer_flx_surfrun_col    (column,:)   = 0._r8
      this%tracer_flx_tparchm_col    (column,:)   = 0._r8
      this%tracer_flx_parchm_vr_col  (column,:,:) = 0._r8
      this%tracer_flx_vtrans_col     (column,:)   = 0._r8
      this%tracer_flx_dew_grnd_col   (column,:)   = 0._r8
      this%tracer_flx_dew_snow_col   (column,:)   = 0._r8
      this%tracer_flx_sub_snow_col   (column,:)   = 0._r8
      this%tracer_flx_h2osfc_snow_residual_col(column,:)   = 0._r8
      this%tracer_flx_netpro_vr_col  (column,:,:)   = 0._r8
      this%tracer_flx_totleached_col (column,:)   = 0._r8
    enddo

  end subroutine Reset

!----------------------------------------
  subroutine Temporal_average(this, column, dtime)
    !
    ! !DESCRIPTION
    ! do temporal average for different fluxes

    !!ARGUMENTS:
    class(TracerFlux_type), intent(inout) :: this
    integer              , intent(in)  :: column     ! column index
    real(r8)             , intent(in)  :: dtime


    this%tracer_flx_top_soil_col   (column,:)   = this%tracer_flx_top_soil_col   (column,:)/dtime
    this%tracer_flx_can_loss_col   (column,:)   = this%tracer_flx_can_loss_col   (column,:)/dtime
    this%tracer_flx_snowmelt_col   (column,:)   = this%tracer_flx_snowmelt_col   (column,:)/dtime
    this%tracer_flx_infl_col       (column,:)   = this%tracer_flx_infl_col       (column,:)/dtime
    this%tracer_flx_netphyloss_col(column,:)    = this%tracer_flx_netphyloss_col (column,:)/dtime
    this%tracer_flx_netpro_col     (column,:)   = this%tracer_flx_netpro_col     (column,:)/dtime
    this%tracer_flx_ebu_col        (column,:)   = this%tracer_flx_ebu_col        (column,:)/dtime
    this%tracer_flx_prec_col       (column,:)   = this%tracer_flx_prec_col       (column,:)/dtime
    this%tracer_flx_dif_col        (column,:)   = this%tracer_flx_dif_col        (column,:)/dtime
    this%tracer_flx_drain_col      (column,:)   = this%tracer_flx_drain_col      (column,:)/dtime
    this%tracer_flx_surfemi_col    (column,:)   = this%tracer_flx_surfemi_col    (column,:)/dtime
    this%tracer_flx_leaching_col   (column,:)   = this%tracer_flx_leaching_col   (column,:)/dtime
    this%tracer_flx_surfrun_col    (column,:)   = this%tracer_flx_surfrun_col    (column,:)/dtime
    this%tracer_flx_tparchm_col    (column,:)   = this%tracer_flx_tparchm_col    (column,:)/dtime
    this%tracer_flx_vtrans_col     (column,:)   = this%tracer_flx_vtrans_col     (column,:)/dtime
    this%tracer_flx_dew_grnd_col   (column,:)   = this%tracer_flx_dew_grnd_col   (column,:)/dtime
    this%tracer_flx_dew_snow_col   (column,:)   = this%tracer_flx_dew_snow_col   (column,:)/dtime
    this%tracer_flx_sub_snow_col   (column,:)   = this%tracer_flx_sub_snow_col   (column,:)/dtime
    this%tracer_flx_h2osfc_snow_residual_col(column,:) =  this%tracer_flx_h2osfc_snow_residual_col(column,:)/dtime

    this%tracer_flx_totleached_col(column,:) = this%tracer_flx_drain_col(column,:) + this%tracer_flx_leaching_col(column,:)
  end subroutine temporal_average

  !----------------------------------------------------------------
  subroutine Flux_summary(this, col, betr_time, c, betrtracer_vars, bstatus)
    !
    ! aggregate fluxes for mass balance check
    ! USES
    use BetrTracerType , only : betrtracer_type
    use tracer_varcon  , only : nlevtrc_soil => betr_nlevtrc_soil
    use MathfuncMod    , only : dot_sum
    use BeTR_TimeMod   , only : betr_time_type
    use BetrStatusType, only : betr_status_type
    use betr_columnType  , only : betr_column_type
    implicit none
    class(TracerFlux_type) , intent(inout) :: this
    type(betr_column_type) , intent(in) :: col
    class(betr_time_type)  , intent(in) :: betr_time
    type(BeTRTracer_Type)  , intent(in) :: betrtracer_vars
    integer                , intent(in) :: c     ! column index
    type(betr_status_type) , intent(out):: bstatus

    !local variables
    integer :: jj, kk
    real(r8):: dtime
    call bstatus%reset()
    associate(                                                             &
         ntracers               => betrtracer_vars%ntracers              , &
         nvolatile_tracers      => betrtracer_vars%nvolatile_tracers     , &
         ngwmobile_tracers      => betrtracer_vars%ngwmobile_tracers     , &
         is_volatile            => betrtracer_vars%is_volatile           , &
         tracernames            => betrtracer_vars%tracernames           , &
         volatileid             => betrtracer_vars%volatileid              &
         )
      dtime = betr_time%get_step_size()

      do jj = 1, ngwmobile_tracers
         !the total net physical loss currently includes infiltration, surface runoff, transpiration aided transport,
         !lateral drainage, vertical leaching
         !for volatile tracers, this includes surface emission surface three different pathways
         this%tracer_flx_infl_col(c,jj) = this%tracer_flx_infl_col(c,jj)*dtime

         this%tracer_flx_netphyloss_col(c,jj) = - this%tracer_flx_infl_col(c,jj) - this%tracer_flx_dew_grnd_col(c,jj) &
              - this%tracer_flx_dew_snow_col(c,jj) - this%tracer_flx_h2osfc_snow_residual_col(c,jj) &
              + this%tracer_flx_sub_snow_col(c,jj) + this%tracer_flx_drain_col(c,jj) + &
              this%tracer_flx_surfrun_col(c,jj) + this%tracer_flx_vtrans_col(c,jj) + this%tracer_flx_leaching_col(c,jj)

         if(is_volatile(jj))then
            kk = volatileid(jj)
            this%tracer_flx_tparchm_col(c,kk) = dot_sum(x=this%tracer_flx_parchm_vr_col(c,1:nlevtrc_soil,kk), &
                 y=col%dz(c,1:nlevtrc_soil), bstatus=bstatus)
            if(bstatus%check_status())return
            this%tracer_flx_surfemi_col(c,kk) = this%tracer_flx_tparchm_col(c,kk) + this%tracer_flx_dif_col(c,kk) + &
                 this%tracer_flx_ebu_col(c,kk)

            this%tracer_flx_netphyloss_col(c,jj) = this%tracer_flx_netphyloss_col(c,jj)  +  this%tracer_flx_surfemi_col(c,kk)

         endif
      enddo

      do jj = 1, ntracers
         this%tracer_flx_netpro_col(c,jj) = dot_sum(x=this%tracer_flx_netpro_vr_col(c,1:nlevtrc_soil,jj),&
             y=col%dz(c,1:nlevtrc_soil), bstatus=bstatus)
         if(bstatus%check_status())return
         if(jj<=ngwmobile_tracers)then
            if(is_volatile(jj))then
               kk = volatileid(jj)
               this%tracer_flx_netpro_col(c,jj) = this%tracer_flx_netpro_col(c,jj) + this%tracer_flx_tparchm_col(c,kk)
            endif
         endif
      enddo
    end associate
  end subroutine Flux_summary


  !----------------------------------------------------------------
  subroutine Flux_display(this, c, jj, betrtracer_vars, msg)
    !
    ! aggregate fluxes for mass balance check

    use BetrTracerType        , only : betrtracer_type
    use betr_constants        , only : betr_errmsg_len

    class(TracerFlux_type), intent(inout)  :: this
    type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
    integer                , intent(in)  :: c     ! column index
    integer                , intent(in)  :: jj
    character(len=betr_errmsg_len), intent(out) :: msg
    character(len=betr_errmsg_len) :: msg1
    !local variables
    integer :: kk

    associate(                                                             &
         ntracers               => betrtracer_vars%ntracers              , &
         nvolatile_tracers      => betrtracer_vars%nvolatile_tracers     , &
         ngwmobile_tracers      => betrtracer_vars%ngwmobile_tracers     , &
         is_volatile            => betrtracer_vars%is_volatile           , &
         tracernames            => betrtracer_vars%tracernames           , &
         volatileid             => betrtracer_vars%volatileid              &
         )
      !the total net physical loss currently includes infiltration, surface runoff, transpiration aided transport,
      !lateral drainage, vertical leaching
      !for volatile tracers, this includes surface emission surface three different pathways
      write(msg,*)tracernames(jj), new_line('A')//'infl=',this%tracer_flx_infl_col(c,jj),&
           ',drain=',  this%tracer_flx_drain_col(c,jj),    &
           ',surfrun=',this%tracer_flx_surfrun_col(c,jj),',vtrans=', this%tracer_flx_vtrans_col(c,jj),&
           ',leaching=', this%tracer_flx_leaching_col(c,jj),',dew_grnd=',this%tracer_flx_dew_grnd_col(c, jj),&
           ',dew_snow=', this%tracer_flx_dew_snow_col(c, jj),',sub_snow=',this%tracer_flx_sub_snow_col(c,jj)

      if(is_volatile(jj))then
         kk = volatileid(jj)
         write(msg1,*) ',tpartm=', this%tracer_flx_tparchm_col(c,kk),',dif=', this%tracer_flx_dif_col(c,kk),  &
              ',ebu=',this%tracer_flx_ebu_col(c,kk)
      endif
      msg = trim(msg)//new_line('A')//trim(msg1)

    end associate
  end subroutine Flux_display

  !----------------------------------------------------------------
  subroutine retrieve_hist(this, bounds, lbj, ubj, flux_2d, flux_1d, betrtracer_vars)
  !DESCRIPTION
  !retrieve variable for history output
  use MathfuncMod, only : addone
  use BeTRTracerType , only : BeTRTracer_Type
  implicit none
  class(TracerFlux_type), intent(inout) :: this
  type(bounds_type)    , intent(in)  :: bounds
  integer, intent(in) :: lbj, ubj
  real(r8), intent(inout) :: flux_2d(bounds%begc:bounds%endc, lbj:ubj,1:this%num_hist2d)
  real(r8), intent(inout) :: flux_1d(bounds%begc:bounds%endc, 1:this%num_hist1d)
  type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
  integer :: begc, endc
  integer :: jj, kk, id
  integer :: idtemp1d, idtemp2d

   associate(                                                    &
      ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers   , &
      ntracers          => betrtracer_vars%ntracers            , &
      is_volatile       => betrtracer_vars%is_volatile         , &
      volatileid        => betrtracer_vars%volatileid          , &
      tracernames       => betrtracer_vars%tracernames           &
    )

    idtemp1d = 0; idtemp2d = 0
    begc=bounds%begc; endc=bounds%endc

    do jj = 1, ntracers
      if(jj<= ngwmobile_tracers) then

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_dew_grnd_col (begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_dew_snow_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_h2osfc_snow_residual_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_sub_snow_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_top_soil_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_can_loss_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_snowmelt_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_infl_col(begc:endc, jj)

        id=addone(idtemp2d); flux_2d(begc:endc,lbj:ubj, id) = this%tracer_flx_netpro_vr_col(begc:endc, lbj:ubj, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_leaching_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_surfrun_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_vtrans_col(begc:endc, jj)

        id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_totleached_col(begc:endc, jj)

        if(is_volatile(jj))then
          kk = volatileid(jj)
          id=addone(idtemp1d);flux_1d(begc:endc,id) = this%tracer_flx_ebu_col(begc:endc, kk)

          id=addone(idtemp1d);flux_1d(begc:endc,id) = this%tracer_flx_dif_col(begc:endc, kk)

          id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_tparchm_col(begc:endc, kk)

          id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_surfemi_col(begc:endc, kk)
         endif

         id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_drain_col(begc:endc, jj)
      endif
      id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_netphyloss_col(begc:endc, jj)

      id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_netpro_col(begc:endc, jj)

      id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_dstor_col(begc:endc, jj)

      id=addone(idtemp1d); flux_1d(begc:endc,id) = this%tracer_flx_prec_col(begc:endc, jj)
    enddo
  end associate
  end subroutine retrieve_hist
end module TracerFluxType
