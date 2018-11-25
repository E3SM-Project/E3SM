module ColumnDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
  use clm_varcon     , only : spval, ispval
  use histFileMod    , only : hist_addfld1d
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_state
    real(r8), pointer :: t_h2osfc      (:)   => null() ! surface water temperature (K)
    real(r8), pointer :: t_h2osfc_bef  (:)   => null() ! surface water temperature at start of time step (K)
    real(r8), pointer :: t_ssbef       (:,:) => null() ! col soil/snow temperature before update (-nlevsno+1:nlevgrnd) 
  contains
    procedure, public :: Init    => col_es_init
    procedure, public :: Restart => col_es_restart
    procedure, public :: Clean   => col_es_clean
  end type column_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ws
    procedure, public :: Clean => clean_col_ws
  end type column_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cs
    procedure, public :: Clean => clean_col_cs
  end type column_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ns
    procedure, public :: Clean => clean_col_ns
  end type column_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ps
    procedure, public :: Clean => clean_col_ps
  end type column_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_ef
    procedure, public :: Clean => clean_col_ef
  end type column_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_water_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_wf
    procedure, public :: Clean => clean_col_wf
  end type column_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_cf
    procedure, public :: Clean => clean_col_cf
  end type column_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_nf
    procedure, public :: Clean => clean_col_nf
  end type column_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_col_pf
    procedure, public :: Clean => clean_col_pf
  end type column_phosphorus_flux

   
  !-----------------------------------------------------------------------
  ! declare the public instances of column-level data types
  !-----------------------------------------------------------------------
  type(column_energy_state)          , public, target :: col_es    ! column energy state
  type(column_water_state)           , public, target :: col_ws    ! column water state
  type(column_carbon_state)          , public, target :: col_cs    ! column carbon state
  type(column_nitrogen_state)        , public, target :: col_ns    ! column nitrogen state
  type(column_phosphorus_state)      , public, target :: col_ps    ! column phosphorus state
  type(column_energy_flux)           , public, target :: col_ef    ! column energy flux
  type(column_water_flux)            , public, target :: col_wf    ! column water flux
  type(column_carbon_flux)           , public, target :: col_cf    ! column carbon flux
  type(column_nitrogen_flux)         , public, target :: col_nf    ! column nitrogen flux
  type(column_phosphorus_flux)       , public, target :: col_pf    ! column phosphorus flux

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy state data structure
  !------------------------------------------------------------------------
  subroutine col_es_init(this, begc, endc)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
    !-----------------------------------------------------------------------
    ! allocate for each member of col_es
    !-----------------------------------------------------------------------
    allocate(this%t_h2osfc         (begc:endc))                     ; this%t_h2osfc           (:)   = nan
    allocate(this%t_h2osfc_bef     (begc:endc))                     ; this%t_h2osfc_bef       (:)   = nan
    allocate(this%t_ssbef          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_ssbef            (:,:) = nan

    !-----------------------------------------------------------------------
    ! initialize history fields for each member of col_es
    !-----------------------------------------------------------------------
    this%t_h2osfc(begc:endc) = spval
    call hist_addfld1d (fname='TH2OSFC',  units='K',  &
         avgflag='A', long_name='surface water temperature', &
         ptr_col=this%t_h2osfc)

    !-----------------------------------------------------------------------
    ! set cold-start initial values for each member of col_es
    !-----------------------------------------------------------------------
    this%t_h2osfc = 274._r8
    
  end subroutine col_es_init
    
  !------------------------------------------------------------------------
  subroutine col_es_restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write column energy state information to/from restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar   ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=this%t_h2osfc)
    if (flag=='read' .and. .not. readvar) then
       this%t_h2osfc(bounds%begc:bounds%endc) = 274.0_r8
    end if

  end subroutine col_es_restart

  !------------------------------------------------------------------------
  subroutine col_es_clean(this)
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_h2osfc)
  end subroutine col_es_clean
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ws(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ws
    
  !------------------------------------------------------------------------
  subroutine clean_col_ws(this)
    !
    ! !ARGUMENTS:
    class(column_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ws
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon state data structure
  !------------------------------------------------------------------------
  subroutine init_col_cs(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cs
    
  !------------------------------------------------------------------------
  subroutine clean_col_cs(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cs
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ns(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ns
    
  !------------------------------------------------------------------------
  subroutine clean_col_ns(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ns
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine init_col_ps(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ps
    
  !------------------------------------------------------------------------
  subroutine clean_col_ps(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ps

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_ef(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_ef
    
  !------------------------------------------------------------------------
  subroutine clean_col_ef(this)
    !
    ! !ARGUMENTS:
    class(column_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_ef
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column water flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_wf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_wf
    
  !------------------------------------------------------------------------
  subroutine clean_col_wf(this)
    !
    ! !ARGUMENTS:
    class(column_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_wf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column carbon flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_cf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_cf
    
  !------------------------------------------------------------------------
  subroutine clean_col_cf(this)
    !
    ! !ARGUMENTS:
    class(column_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_cf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_nf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_nf
    
  !------------------------------------------------------------------------
  subroutine clean_col_nf(this)
    !
    ! !ARGUMENTS:
    class(column_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_nf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine init_col_pf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_col_pf
    
  !------------------------------------------------------------------------
  subroutine clean_col_pf(this)
    !
    ! !ARGUMENTS:
    class(column_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_col_pf
  
    !------------------------------------------------------------------------
    
end module ColumnDataType

