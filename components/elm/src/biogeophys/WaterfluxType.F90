module WaterfluxType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use decompMod    , only : bounds_type, get_proc_global
  use elm_varcon   , only : spval
  use LandunitType , only : lun_pp                
  use ColumnType   , only : col_pp                
  use VegetationType    , only : veg_pp                
  use AnnualFluxDribbler, only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterflux_type

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_evap_grnd_patch     (:)   ! patch ground surface evaporation rate (mm H2O/s) [+]
     real(r8), pointer :: qflx_evap_grnd_col       (:)   ! col ground surface evaporation rate (mm H2O/s) [+]
     real(r8), pointer :: qflx_phs_neg_col         (:)   ! col sum of negative hydraulic redistribution fluxes (mm H2O/s) [+]
     real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_dew_snow_patch      (:)   ! patch surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_snow_col        (:)   ! col surface dew added to snow pack (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_patch      (:)   ! patch ground surface dew formation (mm H2O /s) [+]
     real(r8), pointer :: qflx_dew_grnd_col        (:)   ! col ground surface dew formation (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)


     real(r8), pointer :: dwb_col                  (:)   ! coll water mass change [+ increase] [mm H2O/s] 
     real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
     real(r8), pointer :: qflx_surf_col            (:)   ! col surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff
     real(r8), pointer :: qflx_snow_melt_col       (:)   ! col snow melt (net)
     real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)



     ! Dynamic land cover change

     ! Irrigation
     real(r8), pointer :: qflx_irrig_patch         (:)   ! patch irrigation flux (mm H2O/s)
     real(r8), pointer :: qflx_real_irrig_patch    (:)   ! patch real irrigation flux (mm H2O/s) 
     real(r8), pointer :: qflx_grnd_irrig_patch    (:)   ! groundwater irrigation (mm H2O/s) 
     real(r8), pointer :: qflx_surf_irrig_patch    (:)   ! surface water irrigation(mm H2O/s) 
     real(r8), pointer :: qflx_supply_patch        (:)   ! patch supply flux (mm H2O/s) 
     real(r8), pointer :: qflx_over_supply_patch   (:)   ! over supplied irrigation 

     ! For VSFM


   contains
 
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: Reset     
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type waterflux_type
  !------------------------------------------------------------------------
  type(waterflux_type)  ,public   :: waterflux_vars

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(waterflux_type) :: this
    type(bounds_type), intent(in)    :: bounds  

    call this%InitAllocate(bounds) ! same as "call initAllocate_type(hydro, bounds)"
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use elm_varpar     , only : nlevsno, nlevgrnd, nlevsoi
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: ncells
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg


    allocate(this%qflx_dew_grnd_patch      (begp:endp))              ; this%qflx_dew_grnd_patch      (:)   = nan
    allocate(this%qflx_dew_snow_patch      (begp:endp))              ; this%qflx_dew_snow_patch      (:)   = nan

    allocate(this%qflx_tran_veg_col        (begc:endc))              ; this%qflx_tran_veg_col        (:)   = nan
    allocate(this%qflx_evap_grnd_col       (begc:endc))              ; this%qflx_evap_grnd_col       (:)   = nan
    allocate(this%qflx_phs_neg_col         (begc:endc))              ; this%qflx_phs_neg_col       (:)   = nan
    allocate(this%qflx_dew_grnd_col        (begc:endc))              ; this%qflx_dew_grnd_col        (:)   = nan
    allocate(this%qflx_dew_snow_col        (begc:endc))              ; this%qflx_dew_snow_col        (:)   = nan
    allocate(this%qflx_evap_grnd_patch     (begp:endp))              ; this%qflx_evap_grnd_patch     (:)   = nan

    allocate(this%dwb_col                  (begc:endc))              ; this%dwb_col                  (:)   = nan

    allocate(this%qflx_drain_vr_col        (begc:endc,1:nlevgrnd))   ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_infl_col            (begc:endc))              ; this%qflx_infl_col            (:)   = nan
    allocate(this%qflx_surf_col            (begc:endc))              ; this%qflx_surf_col            (:)   = nan
    allocate(this%qflx_drain_col           (begc:endc))              ; this%qflx_drain_col           (:)   = nan
    allocate(this%qflx_top_soil_col        (begc:endc))              ; this%qflx_top_soil_col        (:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_snow_melt_col       (begc:endc))              ; this%qflx_snow_melt_col       (:)   = nan


    allocate(this%qflx_irrig_patch         (begp:endp))              ; this%qflx_irrig_patch         (:)   = nan 
    allocate(this%qflx_real_irrig_patch    (begp:endp))              ; this%qflx_real_irrig_patch    (:)   = nan 
    allocate(this%qflx_grnd_irrig_patch    (begp:endp))              ; this%qflx_grnd_irrig_patch    (:)   = nan
    allocate(this%qflx_surf_irrig_patch    (begp:endp))              ; this%qflx_surf_irrig_patch    (:)   = nan
    allocate(this%qflx_supply_patch        (begp:endp))              ; this%qflx_supply_patch        (:)   = nan 
    allocate(this%qflx_over_supply_patch   (begp:endp))              ; this%qflx_over_supply_patch   (:)   = nan 

    

    ncells = endc - begc + 1

    ncells = (endc - begc + 1)*nlevgrnd


  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use elm_varctl     , only : create_glacier_mec_landunit, use_cn, use_lch4
    use elm_varpar     , only : nlevsno, crop_prog, nlevsoi 
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    character(10)     :: active
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg


 
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use landunit_varcon, only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    !-----------------------------------------------------------------------

    this%qflx_evap_grnd_patch(bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_grnd_patch (bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_dew_snow_patch (bounds%begp:bounds%endp) = 0.0_r8

    this%qflx_evap_grnd_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_grnd_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_dew_snow_col (bounds%begc:bounds%endc) = 0.0_r8
    
    this%qflx_phs_neg_col(bounds%begc:bounds%endc)   = 0.0_r8

    this%qflx_h2osfc_surf_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_snow_melt_col(bounds%begc:bounds%endc)   = 0._r8

    this%dwb_col(bounds%begc:bounds%endc) = 0._r8
    ! needed for NitrogenLeaching 
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (col_pp%is_soil(c) .or. col_pp%is_crop(c)) then
          this%qflx_drain_col(c) = 0._r8
          this%qflx_surf_col(c)  = 0._r8
       end if
    end do

    do p = bounds%begp, bounds%endp
       l = veg_pp%landunit(p)
       
       if (veg_pp%is_on_soil_col(p)) then
       end if
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use spmdMod          , only : masterproc
    use elm_varcon       , only : denice, denh2o, pondmx, watmin
    use landunit_varcon  , only : istcrop, istdlak, istsoil 
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use elm_varpar       , only : nlevgrnd, nlevurb, nlevsno   
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterflux_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of pfts, globally
    integer :: err_code     ! error code
    logical :: do_io
    !-----------------------------------------------------------------------


  end subroutine Restart

  
  subroutine Reset(this, bounds, numf, filter)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(waterflux_type)              :: this
    type(bounds_type)    , intent(in)  :: bounds
    integer              , intent(in)  :: numf
    integer              , intent(in)  :: filter(:)
    !-----------------------------------------------------------------------
    
    integer :: fc, column
    
    do fc = 1, numf
      column = filter(fc)
    enddo  
  end subroutine Reset    
  
end module WaterfluxType
