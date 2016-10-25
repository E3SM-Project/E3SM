module ColumnEnergyFluxType

  #include "shr_assert.h"

  !------------------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use clm_varcon     , only : spval
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : pft                
  !-------------------------------------------------------------------------------
  implicit none
  save
  private


  type, public :: soilcol_energy_flux
   real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_r_col      (:)   ! col rural snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_snomelt_u_col      (:)   ! col urban snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
   real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
   real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
   real(r8), pointer :: eflx_building_heat_col  (:)   ! col heat flux from urban building interior to urban walls, roof (W/m**2)
   real(r8), pointer :: eflx_urban_ac_col       (:)   ! col urban air conditioning flux (W/m**2)
   real(r8), pointer :: eflx_urban_heat_col     (:)   ! col urban heating flux (W/m**2)
   real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)

   ! Latent heat
   real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

   ! Balance Checks
   real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
   real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
   real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)

  contains
    procedure, public :: Init => init_col_ef
    procedure, public :: Restart => restart_col_ef
    procedure, public :: Clean => clean_col_ef
    procedure, private :: InitAcclocate => initallocate_col_ef
    procedure, private :: InitHistory => inithistory_col_ef
    procedure, private :: InitCold => initcold_col_ef

  end type soilcol_energy_flux
  ! ------------------------------------------------------------------------

  ! Declare public api
  type(soilcol_energy_flux), public, target :: col_ef


  ! Begin Subroutines
  ! --------------------------------------------------------------------

  subroutine init_col_ef(this, bounds, t_grnd_col)

    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds, t_grnd_col ) 

  end subroutine init_col_ef


  subroutine initallocate_col_ef(this, begc, endc)
    class(soilcol_energy_flux) :: this
    integer, intent(in) :: begc   ! beginning soil column index
    integer, intent(in) :: endc   ! ending soil column index    

    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
    allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = nan
    allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = nan
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan
    allocate( this%eflx_building_heat_col  (begc:endc))             ; this%eflx_building_heat_col  (:)   = nan
    allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = nan
    allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = nan
    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan

  end subroutine init_col_ef


  subroutine inithistory_col_ef(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog 
    use clm_varctl     , only : use_cn
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
    !
    ! !ARGUMENTS:
    class(soilcol_energy_flux) :: this
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
    
    this%eflx_snomelt_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM',  units='W/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=this%eflx_snomelt_col, c2l_scale_type='urbanf')

    this%eflx_snomelt_r_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_R',  units='W/m^2',  &
         avgflag='A', long_name='Rural snow melt heat flux', &
         ptr_col=this%eflx_snomelt_r_col, set_spec=spval)

    this%eflx_snomelt_u_col(begc:endc) = spval
    call hist_addfld1d (fname='FSM_U',  units='W/m^2',  &
         avgflag='A', long_name='Urban snow melt heat flux', &
         ptr_col=this%eflx_snomelt_u_col, c2l_scale_type='urbanf', set_nourb=spval)

    this%eflx_building_heat_col(begc:endc) = spval
    call hist_addfld1d (fname='BUILDHEAT', units='W/m^2',  &
         avgflag='A', long_name='heat flux from urban building interior to walls and roof', &
         ptr_col=this%eflx_building_heat_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_ac_col(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
         avgflag='A', long_name='urban air conditioning flux', &
         ptr_col=this%eflx_urban_ac_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_urban_heat_col(begc:endc) = spval
    call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
         avgflag='A', long_name='urban heating flux', &
         ptr_col=this%eflx_urban_heat_col, set_nourb=0._r8, c2l_scale_type='urbanf')

    this%eflx_fgr12_col(begc:endc) = spval
    call hist_addfld1d (fname='FGR12',  units='W/m^2',  &
         avgflag='A', long_name='heat flux between soil layers 1 and 2', &
         ptr_col=this%eflx_fgr12_col, set_lake=spval)

    this%eflx_fgr_col(begc:endc,:) = spval
    call hist_addfld2d (fname='FGR_SOIL_R', units='watt/m^2', type2d='levgrnd', &
         avgflag='A', long_name='Rural downward heat flux at interface below each soil layer', &
         ptr_col=this%eflx_fgr_col, set_spec=spval, default='inactive')

    this%errsoi_col(begc:endc) = spval
    call hist_addfld1d (fname='ERRSOI',  units='W/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=this%errsoi_col)

  end subroutine inithistory_col_ef


 !-----------------------------------------------------------------------
  subroutine initcold_col_ef(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevlak, nlevurb
    use clm_varcon      , only : denice, denh2o, sb
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istice_mec
    use column_varcon   , only : icol_road_imperv, icol_roof, icol_sunwall
    use column_varcon   , only : icol_shadewall, icol_road_perv
    use clm_varctl      , only : iulog, use_vancouver, use_mexicocity
    !
    ! !ARGUMENTS:
    class(soilcol_energy_flux)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: j,l,c,p,levs,lev
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    associate(snl => col%snl) ! Output: [integer (:)    ]  number of snow layers   

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)

         if (lun%urbpoi(l)) then
            this%eflx_building_heat_col(c) = 0._r8
            this%eflx_urban_ac_col(c)      = 0._r8
            this%eflx_urban_heat_col(c)    = 0._r8
         else
            this%eflx_building_heat_col(c) = 0._r8
            this%eflx_urban_ac_col(c)      = 0._r8
            this%eflx_urban_heat_col(c)    = 0._r8
         end if

    end do      

    end associate

  end subroutine initcold_col_ef

  !------------------------------------------------------------------------


  subroutine restart_col_ef(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use shr_log_mod, only : errMsg => shr_log_errMsg
    use spmdMod    , only : masterproc
    use abortutils , only : endrun
    use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(soilcol_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double,  dim1name='column', &
         long_name='urban air conditioning flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac_col)

    call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double, dim1name='column', &
         long_name='urban heating flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat_col)

  end subroutine restart_col_ef


  subroutine clean_col_ef(this)
      class(soilcol_energy_flux) :: this

      deallocate( this%eflx_bot_col)
      deallocate( this%eflx_snomelt_col )
      deallocate( this%eflx_snomelt_r_col )
      deallocate( this%eflx_snomelt_u_col )
      deallocate( this%eflx_fgr12_col )
      deallocate( this%eflx_fgr_col )
      deallocate( this%eflx_building_heat_col )
      deallocate( this%eflx_urban_ac_col )
      deallocate( this%eflx_urban_heat_col )
      deallocate( this%htvp_col )
      deallocate( this%errsoi_col )
      deallocate( this%errseb_col )
      deallocate( this%errsol_col )
      deallocate( this%errlon_col )

  end subroutine clean_col_ef

! --------------------------------------------------------------------
! End Subroutines

end module ColumnEnergyFluxType 
