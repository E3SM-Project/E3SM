module TracerFluxType
  !!DESCRIPTION:
  ! tracer flux type
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : pft
  use clm_varcon     , only : spval, ispval
  use clm_varpar     , only : nlevtrc_soil
  use landunit_varcon, only : istsoil, istcrop
  use clm_varctl     , only : iulog
  !
  ! !PUBLIC TYPES:

  implicit none
  save
  private
  !
  ! !PUBLIC DATA:
  !

  type, public :: TracerFlux_type

     !tracer flux defined at the column level
     real(r8), pointer :: tracer_flx_top_soil_col(:,:)    !tracer fluxes available for infiltration+runoff
     real(r8), pointer :: tracer_flx_can_loss_col(:,:)    !tracer loss from canopy
     real(r8), pointer :: tracer_flx_snowmelt_col(:,:)    !tracer loss from snow melting
     real(r8), pointer :: tracer_flx_infl_col(:,:)        !tracer fluxes available for infiltration
     real(r8), pointer :: tracer_flx_netphyloss_col(:,:)  !total tracer loos through all possible physical pathways: drainage (+ runoff), leaching, ebullition, diffusion, minus precipitation/infiltration
     real(r8), pointer :: tracer_flx_netpro_col(:,:)      !total tracer production through chemical processes
     real(r8), pointer :: tracer_flx_dstor_col(:,:)       !net storage of tracer due to input-output, ideally, dstor=netpro-netloss at various scales
     real(r8), pointer :: tracer_flx_ebu_col(:,:)         !tracer emitted as bubbles, mol, lake, volatile
     real(r8), pointer :: tracer_flx_prec_col(:,:)        !tracer added to a column from precipitation, mol
     real(r8), pointer :: tracer_flx_dif_col(:,:)         !tracer emitted through diffusion, unsat, volatile

     real(r8), pointer :: tracer_flx_drain_col(:,:)       !tracer removal through subface drainage
     real(r8), pointer :: tracer_flx_surfemi_col(:,:)     !total emitted tracer fluxes at surface, volatile, including ebullition, diffusion, arenchyma transport
     real(r8), pointer :: tracer_flx_leaching_col(:,:)    !leaching fluxes
     real(r8), pointer :: tracer_flx_surfrun_col(:,:)     !tracer loss thru runoff, mol tracer / second
     real(r8), pointer :: tracer_flx_netpro_vr_col(:,:,:) !total source strength for the tracers, chemical production, root exudation, excludes incoming root transport (by exchange with air) and (infiltration?)
     real(r8), pointer :: tracer_flx_tparchm_col(:,:)     !total tracer flux through plant aerenchyma transport, for volatile species only, mol/m^2/s
     real(r8), pointer :: tracer_flx_parchm_vr_col(:,:,:) !vertical resolved tracer flux through aerenchyma transport, for volatile species only, mol/m^3/s
     real(r8), pointer :: tracer_flx_totleached_col(:,:)  !total leaching flux, vertical + lateral leaching

     real(r8), pointer :: tracer_flx_vtrans_col(:,:)      !column level tracer flux through transpiration
     !real(r8), pointer :: tracer_flx_snowloss_col(:,:)    !tracer flux lost from snow dynamics, place holder

     !tracer fluxes defined at the pft level
     real(r8), pointer :: tracer_flx_vtrans_patch(:,:)             !tracer goes to the pathway of plant transpiration, currently not released, if it is nutrient, assumed it is taken by plants completely
     real(r8), pointer :: tracer_flx_snowfall_grnd_patch(:,:)
     real(r8), pointer :: tracer_flx_rainfall_grnd_patch(:,:)
     real(r8), pointer :: tracer_flx_prec_intr_patch(:,:)          !interception of tracer from wet deposition [mol/s]
     real(r8), pointer :: tracer_flx_prec_grnd_patch(:,:)          !tracer onto ground including from canopy runoff [mol /s]
     real(r8), pointer :: tracer_flx_snwcp_liq_patch(:,:)          !excess rainfall tracer due to snow capping [mol /s]
     real(r8), pointer :: tracer_flx_snwcp_ice_patch(:,:)          !excess snowfall tracer due to snow capping [mol /s], this is used for aerosol type and water type tracer input
     real(r8), pointer :: tracer_flx_dew_grnd_col(:,:)             !tracer flux to ground coming from dew formation
     real(r8), pointer :: tracer_flx_dew_snow_col(:,:)             !tracer flux to snow coming from dew formation
     real(r8), pointer :: tracer_flx_sub_snow_col(:,:)             !tracer flux loss from snow sublimation
     real(r8), pointer :: tracer_flx_h2osfc_snow_residual_col(:,:) !tracer flux coming from residual standing water and residual snow

   contains
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: Reset
     procedure, public  :: Temporal_average
     procedure, public  :: Flux_summary
     procedure, public  :: Flux_display
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
    class(TracerFlux_type)            :: this
    type(bounds_type), intent(in)     :: bounds
    integer, intent(in)               :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars


    call this%InitAllocate(bounds, lbj, ubj, betrtracer_vars)
    call this%InitHistory(bounds, betrtracer_vars)
    call this%InitCold(bounds)

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
    class(TracerFlux_type)            :: this
    type(bounds_type), intent(in)     :: bounds
    integer, intent(in)               :: lbj, ubj
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars
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

    allocate(this%tracer_flx_prec_col             (begc:endc, 1:ntracers)); this%tracer_flx_prec_col           (:,:) = nan
    allocate(this%tracer_flx_vtrans_patch         (begp:endp, 1:ntracers)); this%tracer_flx_vtrans_patch       (:,:) = nan
    allocate(this%tracer_flx_snowfall_grnd_patch  (begp:endp, 1:ntracers)); this%tracer_flx_snowfall_grnd_patch(:,:) = nan
    allocate(this%tracer_flx_rainfall_grnd_patch  (begp:endp, 1:ntracers)); this%tracer_flx_rainfall_grnd_patch(:,:) = nan
    allocate(this%tracer_flx_prec_intr_patch      (begp:endp, 1:ntracers)); this%tracer_flx_prec_intr_patch    (:,:) = nan
    allocate(this%tracer_flx_prec_grnd_patch(begp:endp, 1:ntracers)); this%tracer_flx_prec_grnd_patch(:,:) = nan
    allocate(this%tracer_flx_snwcp_liq_patch      (begp:endp, 1:ntracers)); this%tracer_flx_snwcp_liq_patch    (:,:) = nan
    allocate(this%tracer_flx_snwcp_ice_patch      (begp:endp, 1:ntracers)); this%tracer_flx_snwcp_ice_patch    (:,:) = nan

    if(ngwmobile_tracers>0)then
       allocate(this%tracer_flx_drain_col       (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_drain_col    (:,:) = nan
       allocate(this%tracer_flx_top_soil_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_top_soil_col (:,:) = nan
       allocate(this%tracer_flx_can_loss_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_can_loss_col (:,:) = nan
       allocate(this%tracer_flx_snowmelt_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_snowmelt_col (:,:) = nan
       allocate(this%tracer_flx_infl_col        (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_infl_col     (:,:) = nan
       allocate(this%tracer_flx_leaching_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_leaching_col (:,:) = nan
       allocate(this%tracer_flx_surfrun_col     (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_surfrun_col  (:,:) = nan
       allocate(this%tracer_flx_vtrans_col      (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_vtrans_col   (:,:) = nan
       allocate(this%tracer_flx_dew_grnd_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_dew_grnd_col (:,:) = nan
       allocate(this%tracer_flx_dew_snow_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_dew_snow_col (:,:) = nan
       allocate(this%tracer_flx_sub_snow_col    (begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_sub_snow_col (:,:) = nan

       allocate(this%tracer_flx_h2osfc_snow_residual_col(begc:endc, 1:ngwmobile_tracers));this%tracer_flx_h2osfc_snow_residual_col(:,:) = nan
       allocate(this%tracer_flx_totleached_col(begc:endc, 1:ngwmobile_tracers)); this%tracer_flx_totleached_col(:,:) = nan
    endif
    if(nvolatile_tracers>0)then
       allocate(this%tracer_flx_ebu_col         (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_ebu_col      (:,:) = nan
       allocate(this%tracer_flx_dif_col         (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_dif_col      (:,:) = nan
       allocate(this%tracer_flx_tparchm_col     (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_tparchm_col  (:,:) = nan
       allocate(this%tracer_flx_surfemi_col     (begc:endc, 1:nvolatile_tracers)); this%tracer_flx_surfemi_col  (:,:) = nan
       allocate(this%tracer_flx_parchm_vr_col   (begc:endc, lbj:ubj, 1:nvolatile_tracers)); this%tracer_flx_parchm_vr_col(:,:,:) = nan
    endif

    allocate(this%tracer_flx_netpro_vr_col  (begc:endc, lbj:ubj, 1:ntracers)); this%tracer_flx_netpro_vr_col (:,:,:) = nan
    allocate(this%tracer_flx_netphyloss_col (begc:endc, 1:ntracers)); this%tracer_flx_netphyloss_col(:,:)            = nan
    allocate(this%tracer_flx_netpro_col     (begc:endc, 1:ntracers)); this%tracer_flx_netpro_col(:,:)                = nan
    allocate(this%tracer_flx_dstor_col      (begc:endc, 1:ntracers)); this%tracer_flx_dstor_col(:,:)                 = nan


  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! History fields initialization
    !
    ! !USES:
    !use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only: spval
    use clm_varpar    , only: nlevsno
    use histFileMod   , only: hist_addfld1d, hist_addfld2d
    use histFileMod   , only: no_snow_normal, no_snow_zero
    use BeTRTracerType, only: BeTRTracer_Type
    !
    ! !ARGUMENTS:
    class(TracerFlux_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars

    !
    ! !LOCAL VARIABLES:
    integer :: ntracers
    integer :: ngwmobile_tracers
    integer :: nsolid_passive_tracers
    integer :: jj, kk
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

    associate(                                                   &
      ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers   , &
      ntracers          => betrtracer_vars%ntracers            , &
      is_volatile       => betrtracer_vars%is_volatile         , &
      volatileid        => betrtracer_vars%volatileid          , &
      tracernames       => betrtracer_vars%tracernames           &
    )
    begc=bounds%begc; endc=bounds%endc
    do jj = 1, ntracers
      if(jj<= ngwmobile_tracers) then

        this%tracer_flx_dew_grnd_col (begc:endc, jj) = spval
        data1dptr => this%tracer_flx_dew_grnd_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_DEW_GRND', units='none', &
         avgflag='A', long_name='incoming dew flux to ground for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')

        this%tracer_flx_dew_snow_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_dew_snow_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_DEW_SNOW', units='none', &
         avgflag='A', long_name='incoming dew flux to snow from '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')

        this%tracer_flx_h2osfc_snow_residual_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_h2osfc_snow_residual_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_H2OSFC_SNOW_RES', units='none', &
         avgflag='A', long_name='incoming flux to topsoi from snow and h2osfc residual for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')

        this%tracer_flx_sub_snow_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_sub_snow_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_SUB_SNOW', units='none', &
         avgflag='A', long_name='sublimation flux from snow for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')


        this%tracer_flx_top_soil_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_top_soil_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_TOPSOIL', units='none', &
         avgflag='A', long_name='incoming flux at top of the soil for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')


        this%tracer_flx_can_loss_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_can_loss_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_CAN_LOSS', units='none', &
         avgflag='A', long_name='loss from canopy for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')

        this%tracer_flx_snowmelt_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_snowmelt_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_SNOWMELT', units='none', &
         avgflag='A', long_name='loss from snowmelt for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')

        this%tracer_flx_infl_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_infl_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_INFIL', units='none', &
         avgflag='A', long_name='infiltration for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')


        this%tracer_flx_netpro_vr_col(begc:endc, :, jj)  = spval
        data2dptr =>this%tracer_flx_netpro_vr_col(:,:,jj)
        call hist_addfld2d (fname=trim(tracernames(jj))//'_FLX_NETPRO_vr', units='none', type2d='levtrc',  &
         avgflag='A', long_name='net production for '//trim(tracernames(jj)), &
         ptr_col=data2dptr, default='inactive')

        this%tracer_flx_leaching_col(begc:endc, jj) = spval
        data1dptr =>  this%tracer_flx_leaching_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_LEACHING', units='none', &
         avgflag='A', long_name='bottom of soil leaching for '//trim(tracernames(jj)), &
         ptr_col=data1dptr,  default='inactive')

        this%tracer_flx_surfrun_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_surfrun_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_SRUNOFF', units='none', &
         avgflag='A', long_name='loss from surface runoff for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')


        this%tracer_flx_vtrans_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_vtrans_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_VTRANS', units='none', &
         avgflag='A', long_name='transport through transpiration for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')

        this%tracer_flx_totleached_col(begc:endc, jj) = spval
        data1dptr => this%tracer_flx_totleached_col(:, jj)
        call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_TLEACH', units='none', &
         avgflag='A', long_name='transport through leaching for '//trim(tracernames(jj)), &
         ptr_col=data1dptr, default='inactive')

        if(is_volatile(jj))then
          kk = volatileid(jj)
          this%tracer_flx_ebu_col(begc:endc, kk) = spval
          data1dptr => this%tracer_flx_ebu_col(:, kk)
          call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_EBU', units='none', &
            avgflag='A', long_name='loss through ebullition (+ into atmosphere) for '//trim(tracernames(jj)), &
            ptr_col=data1dptr,  default='inactive')

          this%tracer_flx_dif_col(begc:endc, kk) = spval
          data1dptr => this%tracer_flx_dif_col(:, kk)
          call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_DIF', units='none', &
            avgflag='A', long_name='loss through diffusion (+ into atmosphere) for '//trim(tracernames(jj)), &
            ptr_col=data1dptr,  default='inactive')

          this%tracer_flx_tparchm_col(begc:endc, kk) = spval
          data1dptr => this%tracer_flx_tparchm_col(:, jj)
          call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_ARCHM', units='none', &
           avgflag='A', long_name='loss from aerenchyma transport for '//trim(tracernames(jj)), &
           ptr_col=data1dptr, default='inactive')

         endif

         this%tracer_flx_drain_col(begc:endc, jj) = spval
         data1dptr => this%tracer_flx_drain_col(:, jj)
         call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_DRAIN', units='none', &
          avgflag='A', long_name='loss from drainage for '//trim(tracernames(jj)), &
          ptr_col=data1dptr,  default='inactive')
      endif
      this%tracer_flx_netphyloss_col(begc:endc, jj) = spval
      data1dptr => this%tracer_flx_netphyloss_col(:, jj)
      call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_NETLOSS', units='none', &
        avgflag='A', long_name='net loss for '//trim(tracernames(jj)), &
        ptr_col=data1dptr, default='inactive')

      this%tracer_flx_netpro_col(begc:endc, jj) = spval
      data1dptr => this%tracer_flx_netpro_col(:, jj)
      call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_NETPRO', units='none', &
        avgflag='A', long_name='net production for '//trim(tracernames(jj)), &
        ptr_col=data1dptr, default='inactive')

      this%tracer_flx_dstor_col(begc:endc, jj)  = spval
      data1dptr => this%tracer_flx_dstor_col(:, jj)
      call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_DSTOR', units='none', &
        avgflag='A', long_name='total concentration change for '//trim(tracernames(jj)), &
        ptr_col=data1dptr,  default='inactive')

      this%tracer_flx_prec_col(begc:endc, jj) = spval
      data1dptr => this%tracer_flx_prec_col(:, jj)
      call hist_addfld1d (fname=trim(tracernames(jj))//'_FLX_PREC', units='none', &
        avgflag='A', long_name='incoming from precipitation for '//trim(tracernames(jj)), &
        ptr_col=data1dptr,  default='inactive')

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
    class(TracerFlux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, p, l       ! index
    integer               :: begp, endp
    integer               :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%ifspecial(l)) then
         this%tracer_flx_vtrans_patch(p,:)         = spval
         this%tracer_flx_snowfall_grnd_patch(p,:)  = spval
         this%tracer_flx_rainfall_grnd_patch(p,:)  = spval
         this%tracer_flx_prec_intr_patch(p,:)      = spval
         this%tracer_flx_prec_grnd_patch(p,:)      = spval
         this%tracer_flx_snwcp_liq_patch(p,:)      = spval
         this%tracer_flx_snwcp_ice_patch(p,:)      = spval
       endif
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
         this%tracer_flx_vtrans_patch(p,:)         = 0._r8
         this%tracer_flx_snowfall_grnd_patch(p,:)  = 0._r8
         this%tracer_flx_rainfall_grnd_patch(p,:)  = 0._r8
         this%tracer_flx_prec_intr_patch(p,:)      = 0._r8
         this%tracer_flx_prec_grnd_patch(p,:)      = 0._r8
         this%tracer_flx_snwcp_liq_patch(p,:)      = 0._r8
         this%tracer_flx_snwcp_ice_patch(p,:)      = 0._r8
       endif
    enddo
    do c = begc, endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
         this%tracer_flx_top_soil_col(c,:)    = spval
         this%tracer_flx_can_loss_col(c,:)    = spval
         this%tracer_flx_snowmelt_col (c,:)    = spval
         this%tracer_flx_infl_col(c,:)        = spval
         this%tracer_flx_netphyloss_col(c,:)  = spval
         this%tracer_flx_netpro_col(c,:)      = spval
         this%tracer_flx_dstor_col(c,:)       = spval
         this%tracer_flx_ebu_col(c,:)         = spval
         this%tracer_flx_prec_col(c,:)        = spval
         this%tracer_flx_dif_col(c,:)         = spval
         this%tracer_flx_drain_col(c,:)       = spval
         this%tracer_flx_surfemi_col(c,:)     = spval
         this%tracer_flx_leaching_col(c,:)    = spval
         this%tracer_flx_surfrun_col(c,:)     = spval
         this%tracer_flx_tparchm_col(c,:)     = spval
         this%tracer_flx_parchm_vr_col(c,:,:) = spval
         this%tracer_flx_vtrans_col(c,:)      = spval
         this%tracer_flx_dew_grnd_col   (c,:) = spval
         this%tracer_flx_dew_snow_col   (c,:) = spval
         this%tracer_flx_sub_snow_col   (c,:) = spval
         this%tracer_flx_h2osfc_snow_residual_col(c,:) = spval
         this%tracer_flx_totleached_col(c,:)  = spval
       endif

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
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
       endif
    enddo


  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! Now it is purposely empty, but will be potentially useful in the future
    ! !USES:
    use BetrTracerType        , only : betrtracer_type
    use clm_varpar , only : nlevsno, nlevsoi
    use clm_varcon , only : spval
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(TracerFlux_type) :: this
    type(bounds_type)    , intent(in)    :: bounds
    type(file_desc_t)    , intent(inout) :: ncid                                         ! netcdf id
    character(len=*)     , intent(in)    :: flag                                         ! 'read' or 'write'
    type(BeTRTracer_Type), intent(in)    :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file


   end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, bounds, numf, filter)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(TracerFlux_type)             :: this
    type(bounds_type)    , intent(in)  :: bounds
    integer              , intent(in)  :: numf
    integer              , intent(in)  :: filter(:)
    !-----------------------------------------------------------------------

    integer :: fc, column

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
    class(TracerFlux_type)             :: this
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
  subroutine Flux_summary(this, c, betrtracer_vars)
    !
    ! aggregate fluxes for mass balance check

    use BetrTracerType        , only : betrtracer_type
    use clm_time_manager      , only : get_step_size
    use clm_varpar            , only : nlevtrc_soil
    use MathfuncMod           , only : dot_sum
    class(TracerFlux_type)               :: this
    type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
    integer                , intent(in)  :: c     ! column index

    !local variables
    integer :: jj, kk
    real(r8):: dtime
    associate(                                                          &
         ntracers               => betrtracer_vars%ntracers              , &
         nvolatile_tracers      => betrtracer_vars%nvolatile_tracers     , &
         ngwmobile_tracers      => betrtracer_vars%ngwmobile_tracers     , &
         is_volatile            => betrtracer_vars%is_volatile           , &
         tracernames            => betrtracer_vars%tracernames           , &
         volatileid             => betrtracer_vars%volatileid              &
         )
      dtime = get_step_size()
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
            this%tracer_flx_tparchm_col(c,kk) = dot_sum(x=this%tracer_flx_parchm_vr_col(c,1:nlevtrc_soil,kk), y=col%dz(c,1:nlevtrc_soil))

            this%tracer_flx_surfemi_col(c,kk) = this%tracer_flx_tparchm_col(c,kk) + this%tracer_flx_dif_col(c,kk) + &
                 this%tracer_flx_ebu_col(c,kk)

            this%tracer_flx_netphyloss_col(c,jj) = this%tracer_flx_netphyloss_col(c,jj)  +  this%tracer_flx_surfemi_col(c,kk)

         endif
      enddo

      do jj = 1, ntracers
         this%tracer_flx_netpro_col(c,jj) = dot_sum(x=this%tracer_flx_netpro_vr_col(c,1:nlevtrc_soil,jj),y=col%dz(c,1:nlevtrc_soil))
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
  subroutine Flux_display(this, c, jj, betrtracer_vars)
    !
    ! aggregate fluxes for mass balance check

    use BetrTracerType        , only : betrtracer_type

    class(TracerFlux_type)               :: this
    type(BeTRTracer_Type)  , intent(in)  :: betrtracer_vars
    integer                , intent(in)  :: c     ! column index
    integer                , intent(in)  :: jj
    !local variables
    integer :: kk

    associate(                                                          &
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
      write(iulog,*)tracernames(jj)
      write(iulog,*)'infl=',this%tracer_flx_infl_col(c,jj),' drain=',  this%tracer_flx_drain_col(c,jj),    &
           ' surfrun=',this%tracer_flx_surfrun_col(c,jj),' vtrans=', this%tracer_flx_vtrans_col(c,jj),&
           ' leaching=', this%tracer_flx_leaching_col(c,jj)

      if(is_volatile(jj))then
         kk = volatileid(jj)
         write(iulog,*)'tpartm=', this%tracer_flx_tparchm_col(c,kk),' dif=', this%tracer_flx_dif_col(c,kk),  &
              ' ebu=',this%tracer_flx_ebu_col(c,kk)
      endif


    end associate
  end subroutine Flux_display

end module TracerFluxType
