module GridcellDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gridcell data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varcon     , only : spval, ispval
  use histFileMod    , only : hist_addfld1d, hist_addfld2d, no_snow_normal
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use ColumnType     , only : col_pp
  use LandunitType   , only : lun_pp
  use GridcellType   , only : grc_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_energy_state
    ! temperature variables
    real(r8), pointer :: heat1                (:)   ! initial gridcell total heat content
    real(r8), pointer :: heat2                (:)   ! post land cover change total heat content
    real(r8), pointer :: liquid_water_temp1   (:)   ! initial weighted average liquid water temperature (K)
    real(r8), pointer :: liquid_water_temp2   (:)   ! post land cover change weighted average liquid water temperature (K)

  contains
    procedure, public :: Init    => grc_es_init
    procedure, public :: Clean   => grc_es_clean
  end type gridcell_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the gridcell level.
  !-----------------------------------------------------------------------
  type, public :: gridcell_energy_flux
    ! temperature variables
    real(r8), pointer :: eflx_dynbal           (:)   ! dynamic land cover change conversion energy flux (W/m**2)

  contains
    procedure, public :: Init    => grc_ef_init
    procedure, public :: Restart => grc_ef_restart
    procedure, public :: Clean   => grc_ef_clean
  end type gridcell_energy_flux
  
  !-----------------------------------------------------------------------
  ! declare the public instances of gridcell-level data types
  !-----------------------------------------------------------------------
  type(gridcell_energy_state)          , public, target :: grc_es    ! gridcell energy state
  type(gridcell_energy_flux)           , public, target :: grc_ef    ! gridcell energy flux
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell energy state data structure
  !------------------------------------------------------------------------
  subroutine grc_es_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_state) :: this
    integer, intent(in) :: begg,endg
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_es
    !-----------------------------------------------------------------------
    allocate(this%heat1                (begg:endg))                      ; this%heat1                (:)   = nan
    allocate(this%heat2                (begg:endg))                      ; this%heat2                (:)   = nan
    allocate(this%liquid_water_temp1   (begg:endg))                      ; this%liquid_water_temp1   (:)   = nan
    allocate(this%liquid_water_temp2   (begg:endg))                      ; this%liquid_water_temp2   (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_es
    !-----------------------------------------------------------------------
    this%heat1(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT1',  units='J/m^2',  &
         avgflag='A', long_name='initial gridcell total heat content', &
         ptr_lnd=this%heat1)

    this%heat2(begg:endg) = spval
    call hist_addfld1d (fname='GC_HEAT2',  units='J/m^2',  &
         avgflag='A', long_name='post land cover change total heat content', &
         ptr_lnd=this%heat2, default='inactive')  

    this%liquid_water_temp1(begg:endg) = spval
    call hist_addfld1d (fname='LIQUID_WATER_TEMP1', units='K', &
         avgflag='A', long_name='initial gridcell weighted average liquid water temperature', &
         ptr_lnd=this%liquid_water_temp1, default='inactive')

  end subroutine grc_es_init

  !------------------------------------------------------------------------
  subroutine grc_es_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%heat1)
    deallocate(this%heat2)
    deallocate(this%liquid_water_temp1)
    deallocate(this%liquid_water_temp2)
  end subroutine grc_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean gridcell energy flux data structure
  !------------------------------------------------------------------------
  subroutine grc_ef_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    integer, intent(in) :: begg,endg
    !------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! allocate for each member of grc_ef
    !-----------------------------------------------------------------------
    allocate(this%eflx_dynbal              (begg:endg))                  ; this%eflx_dynbal             (:)   = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of grc_ef
    !-----------------------------------------------------------------------
    this%eflx_dynbal(begg:endg) = spval 
    call hist_addfld1d (fname='EFLX_DYNBAL',  units='W/m^2',  &
         avgflag='A', long_name='dynamic land cover change conversion energy flux', &
         ptr_lnd=this%eflx_dynbal)

  end subroutine grc_ef_init

  !------------------------------------------------------------------------
  subroutine grc_ef_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write gridcell energy flux information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------
  end subroutine grc_ef_restart

  !------------------------------------------------------------------------
  subroutine grc_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_energy_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%eflx_dynbal)
  end subroutine grc_ef_clean
  

end module GridcellDataType

  
    