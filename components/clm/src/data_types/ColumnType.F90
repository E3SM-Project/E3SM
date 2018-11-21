module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !DW  Converted from ColumnType
  !DW  Change the old function into PhysicalPropertiesType
  ! Column data type allocation and initialization
  ! -------------------------------------------------------- 
  ! column types can have values of
  ! -------------------------------------------------------- 
  !   1  => (istsoil)          soil (vegetated or bare soil)
  !   2  => (istcrop)          crop (only for crop configuration)
  !   3  => (istice)           land ice
  !   4  => (istice_mec)       land ice (multiple elevation classes)   
  !   5  => (istdlak)          deep lake
  !   6  => (istwet)           wetland
  !   71 => (icol_roof)        urban roof
  !   72 => (icol_sunwall)     urban sunwall
  !   73 => (icol_shadewall)   urban shadewall
  !   74 => (icol_road_imperv) urban impervious road
  !   75 => (icol_road_perv)   urban pervious road
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak
  use clm_varcon     , only : spval, ispval
  use clm_varctl     , only : use_cn
  use column_varcon  , only : is_hydrologically_active
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use spmdMod        , only : masterproc
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use histFieldsMod  , only : col_es_init_hist
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds physical property information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_physical_properties
     ! indices and weights for higher subgrid levels (landunit, topounit, gridcell)
     integer , pointer :: gridcell     (:) => null() ! index into gridcell level quantities
     real(r8), pointer :: wtgcell      (:) => null() ! weight (relative to gridcell)
     integer , pointer :: topounit     (:) => null() ! index into topounit level quantities
     real(r8), pointer :: wttopounit   (:) => null() ! weight (relative to topounit)
     integer , pointer :: landunit     (:) => null() ! index into landunit level quantities
     real(r8), pointer :: wtlunit      (:) => null() ! weight (relative to landunit)

     ! Starting and ending indices for subgrid types below the column level
     integer , pointer :: pfti         (:) => null() ! beginning pft index for each column
     integer , pointer :: pftf         (:) => null() ! ending pft index for each column
     integer , pointer :: npfts        (:) => null() ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype        (:) => null() ! column type
     logical , pointer :: active       (:) => null() ! true=>do computations on this column 

     ! topography
     real(r8), pointer :: glc_topo     (:) => null() ! surface elevation (m)
     real(r8), pointer :: micro_sigma  (:) => null() ! microtopography pdf sigma (m)
     real(r8), pointer :: n_melt       (:) => null() ! SCA shape parameter
     real(r8), pointer :: topo_slope   (:) => null() ! gridcell topographic slope
     real(r8), pointer :: topo_std     (:) => null() ! gridcell elevation standard deviation
     integer, pointer  :: nlevbed      (:) => null() ! number of layers to bedrock
     real(r8), pointer :: zibed        (:) => null() ! bedrock depth in model (interface level at nlevbed)

     ! vertical levels
     integer , pointer :: snl          (:)   => null() ! number of snow layers
     real(r8), pointer :: dz           (:,:) => null() ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z            (:,:) => null() ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi           (:,:) => null() ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii          (:)   => null() ! convective boundary height [m]
     real(r8), pointer :: dz_lake      (:,:) => null() ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake       (:,:) => null() ! layer depth for lake (m)
     real(r8), pointer :: lakedepth    (:)   => null() ! variable lake depth (m)                             

     ! other column characteristics
     logical , pointer :: hydrologically_active(:)   ! true if this column is a hydrologically active type

   contains

     procedure, public :: Init => col_pp_init
     procedure, public :: Clean => col_pp_clean

  end type column_physical_properties

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the column level.
  !-----------------------------------------------------------------------
  type, public :: column_energy_state
    real(r8), pointer :: t_h2osfc      (:) => null() ! surface water temperature (K)
  contains
    procedure, public :: Init    => init_col_es
    procedure, public :: Restart => restart_col_es
    procedure, public :: Clean   => clean_col_es
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
  type(column_physical_properties)   , public, target :: col_pp    ! column physical properties
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
  subroutine col_pp_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_physical_properties)  :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------

    ! The following is set in initGridCellsMod
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
    allocate(this%topounit    (begc:endc))                     ; this%topounit    (:)   = ispval
    allocate(this%wttopounit  (begc:endc))                     ; this%wttopounit  (:)   = nan
    allocate(this%landunit    (begc:endc))                     ; this%landunit    (:)   = ispval
    allocate(this%wtlunit     (begc:endc))                     ; this%wtlunit     (:)   = nan
    allocate(this%pfti        (begc:endc))                     ; this%pfti        (:)   = ispval
    allocate(this%pftf        (begc:endc))                     ; this%pftf        (:)   = ispval
    allocate(this%npfts       (begc:endc))                     ; this%npfts       (:)   = ispval
    allocate(this%itype       (begc:endc))                     ; this%itype       (:)   = ispval
    allocate(this%active      (begc:endc))                     ; this%active      (:)   = .false.

    ! The following is set in initVerticalMod
    allocate(this%snl         (begc:endc))                     ; this%snl         (:)   = ispval  !* cannot be averaged up
    allocate(this%dz          (begc:endc,-nlevsno+1:nlevgrnd)) ; this%dz          (:,:) = nan
    allocate(this%z           (begc:endc,-nlevsno+1:nlevgrnd)) ; this%z           (:,:) = nan
    allocate(this%zi          (begc:endc,-nlevsno+0:nlevgrnd)) ; this%zi          (:,:) = nan
    allocate(this%zii         (begc:endc))                     ; this%zii         (:)   = nan
    allocate(this%lakedepth   (begc:endc))                     ; this%lakedepth   (:)   = spval  
    allocate(this%dz_lake     (begc:endc,nlevlak))             ; this%dz_lake     (:,:) = nan
    allocate(this%z_lake      (begc:endc,nlevlak))             ; this%z_lake      (:,:) = nan

    allocate(this%glc_topo    (begc:endc))                     ; this%glc_topo    (:)   = nan
    allocate(this%micro_sigma (begc:endc))                     ; this%micro_sigma (:)   = nan
    allocate(this%n_melt      (begc:endc))                     ; this%n_melt      (:)   = nan 
    allocate(this%topo_slope  (begc:endc))                     ; this%topo_slope  (:)   = nan
    allocate(this%topo_std    (begc:endc))                     ; this%topo_std    (:)   = nan
    allocate(this%nlevbed     (begc:endc))                     ; this%nlevbed     (:)   = ispval
    allocate(this%zibed       (begc:endc))                     ; this%zibed       (:)   = nan

    allocate(this%hydrologically_active(begc:endc))            ; this%hydrologically_active(:) = .false.

  end subroutine col_pp_init

  !------------------------------------------------------------------------
  subroutine col_pp_clean(this)
    !
    ! !ARGUMENTS:
    class(column_physical_properties) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell   )
    deallocate(this%wtgcell    )
    deallocate(this%topounit   )
    deallocate(this%wttopounit )
    deallocate(this%landunit   )
    deallocate(this%wtlunit    )
    deallocate(this%pfti       )
    deallocate(this%pftf       )
    deallocate(this%npfts      )
    deallocate(this%itype      )
    deallocate(this%active     )
    deallocate(this%snl        )
    deallocate(this%dz         )
    deallocate(this%z          )
    deallocate(this%zi         )
    deallocate(this%zii        )
    deallocate(this%lakedepth  )
    deallocate(this%dz_lake    )
    deallocate(this%z_lake     )
    deallocate(this%glc_topo   )
    deallocate(this%micro_sigma)
    deallocate(this%n_melt     )
    deallocate(this%topo_slope )
    deallocate(this%topo_std   )
    deallocate(this%nlevbed    )
    deallocate(this%zibed      )
    deallocate(this%hydrologically_active)

  end subroutine col_pp_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean column energy state data structure
  !------------------------------------------------------------------------
  subroutine init_col_es(this, begc, endc)
    !
    ! !USES:
    !use histFieldsMod  , only : col_es_hist_init
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    allocate(this%t_h2osfc             (begc:endc))              ; this%t_h2osfc           (:)   = nan

    call col_es_init_hist(begc, endc)
    
    ! set cold-start initial values
    this%t_h2osfc = 274._r8
    
  end subroutine init_col_es
    
  !------------------------------------------------------------------------
  subroutine restart_col_es(this, bounds, ncid, flag)
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

  end subroutine restart_col_es

  !------------------------------------------------------------------------
  subroutine clean_col_es(this)
    !
    ! !ARGUMENTS:
    class(column_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_h2osfc)
  end subroutine clean_col_es
  
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
    
end module ColumnType
