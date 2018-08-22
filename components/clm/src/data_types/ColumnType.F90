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
  use column_varcon  , only : is_hydrologically_active
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: column_physical_properties_type
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

  end type column_physical_properties_type

  type(column_physical_properties_type), public, target :: col_pp !column data structure (soil/snow/canopy columns)
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine col_pp_init(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(column_physical_properties_type)  :: this
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
    class(column_physical_properties_type) :: this
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


end module ColumnType
