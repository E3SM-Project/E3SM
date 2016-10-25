module ColumnPhysicalPropertiesType

  ! -------------------------------------------------------- 
  ! ALM sub-grid hierarchy:
  ! Define column physical properties data type 
  ! -------------------------------------------------------- 
  ! 10 Oct 2016, DW  

  ! Moved from ColumnType.F90
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
  use clm_varcon     , only : ispval, spval
  use clm_varpar     , only : nlevsno, nlevgrnd, nlevlak


! moved from ChemStateType
  !------------------------------------------------------------------------------
  ! !USES:
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use clm_varpar      , only : nlevsoi

!  moved from CNDecompCascadeConType.F90

  use decompMod      , only : bounds_type
  use clm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools


  implicit none
  save
  private

  ! sub-grid geospatial and physical properties defined at the soil_column level
  ! migrate variables list from ColumnType.F90

  type, public :: column_phyisical_properties
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: landunit             (:)   ! index into landunit level quantities
     real(r8), pointer :: wtlunit              (:)   ! weight (relative to landunit)
     integer , pointer :: gridcell             (:)   ! index into gridcell level quantities
     real(r8), pointer :: wtgcell              (:)   ! weight (relative to gridcell)
     integer , pointer :: pfti                 (:)   ! beginning pft index for each column
     integer , pointer :: pftf                 (:)   ! ending pft index for each column
     integer , pointer :: npfts                (:)   ! number of patches for each column

     ! topological mapping functionality
     integer , pointer :: itype                (:)   ! column type
     logical , pointer :: active               (:)   ! true=>do computations on this column 

     ! topography
     real(r8), pointer :: glc_topo             (:)   ! surface elevation (m)
     real(r8), pointer :: micro_sigma          (:)   ! microtopography pdf sigma (m)
     real(r8), pointer :: n_melt               (:)   ! SCA shape parameter
     real(r8), pointer :: topo_slope           (:)   ! gridcell topographic slope
     real(r8), pointer :: topo_std             (:)   ! gridcell elevation standard deviation

     ! vertical levels
     integer , pointer :: snl                  (:)   ! number of snow layers
     real(r8), pointer :: dz                   (:,:) ! layer thickness (m)  (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: z                    (:,:) ! layer depth (m) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: zi                   (:,:) ! interface level below a "z" level (m) (-nlevsno+0:nlevgrnd) 
     real(r8), pointer :: zii                  (:)   ! convective boundary height [m]
     real(r8), pointer :: dz_lake              (:,:) ! lake layer thickness (m)  (1:nlevlak)
     real(r8), pointer :: z_lake               (:,:) ! layer depth for lake (m)
     real(r8), pointer :: lakedepth            (:)   ! variable lake depth (m)          
!  moved from chemStateType
     real(r8), pointer :: soil_pH              (:,:) ! soil pH (-nlevsno+1:nlevgrnd)                   

!  moved from CNDecompCascadeConType.F90
     !-- properties of each pathway along decomposition cascade 
     character(len=8)  , pointer :: cascade_step_name(:)               ! name of transition
     integer           , pointer :: cascade_donor_pool(:)              ! which pool is C taken from for a given decomposition step
     integer           , pointer :: cascade_receiver_pool(:)           ! which pool is C added to for a given decomposition step

     !-- properties of each decomposing pool
     logical           , pointer  :: floating_cn_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     logical           , pointer  :: floating_cp_ratio_decomp_pools(:) ! TRUE => pool has fixed C:N ratio
     character(len=8)  , pointer  :: decomp_pool_name_restart(:)       ! name of pool for restart files
     character(len=8)  , pointer  :: decomp_pool_name_history(:)       ! name of pool for history files
     character(len=20) , pointer  :: decomp_pool_name_long(:)          ! name of pool for netcdf long names
     character(len=8)  , pointer  :: decomp_pool_name_short(:)         ! name of pool for netcdf short names
     logical           , pointer  :: is_litter(:)                      ! TRUE => pool is a litter pool
     logical           , pointer  :: is_soil(:)                        ! TRUE => pool is a soil pool
     logical           , pointer  :: is_cwd(:)                         ! TRUE => pool is a cwd pool
     real(r8)          , pointer  :: initial_cn_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_cp_ratio(:)               ! c:n ratio for initialization of pools
     real(r8)          , pointer  :: initial_stock(:)                  ! initial concentration for seeding at spinup
     logical           , pointer  :: is_metabolic(:)                   ! TRUE => pool is metabolic material
     logical           , pointer  :: is_cellulose(:)                   ! TRUE => pool is cellulose
     logical           , pointer  :: is_lignin(:)                      ! TRUE => pool is lignin
     real(r8)          , pointer  :: spinup_factor(:)                  ! factor by which to scale AD and relevant processes by

   contains

     procedure, public :: Init => init_col_pp
     procedure, public :: Clean => clean_col_pp

  end type column_physical_properties

  ! declare the public instances of soilcolumn types
  type(soilcol_properties) , public, target :: col_pp

  contains 
!!!!***********DW***********change to (this, bounds) based on ChemStateType.F90
!  subroutine init_col_pp(this, begc, endc)
!    class(column_properties) :: this
!    integer, intent(in) :: begc   ! beginning soil column index
!    integer, intent(in) :: endc   ! ending soil column index

!  
   subroutine init_col_pp(this, bounds)
!   moved from ChemStateType.f90

!    use clm_varpar            , only : nlevsoi

    class(column_properties)         :: this
    type(bounds_type), intent(in) :: bounds 


    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj

    begc = bounds%begc;
    endc = bounds%endc
    lbj  = 1;
    ubj  = nlevsoi

    ! The following is set in initGridCellsMod  (***DW***: Need to considere Topounit?)
    allocate(this%gridcell    (begc:endc))                     ; this%gridcell    (:)   = ispval
    allocate(this%wtgcell     (begc:endc))                     ; this%wtgcell     (:)   = nan
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
!   moved from ChemStateType.F90
    allocate(this%soil_pH(begc:endc, lbj:ubj))

!  moved from CNDecompCascadeConType.F90
    ! !DESCRIPTION:
    ! Initialize decomposition cascade state
    !------------------------------------------------------------------------

    !-- properties of each pathway along decomposition cascade 
    allocate(decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions))
    allocate(decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions))

    !-- properties of each decomposing pool
    allocate(decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%floating_cp_ratio_decomp_pools(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools))
    allocate(decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_litter(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_soil(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cwd(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_cp_ratio(0:ndecomp_pools))
    allocate(decomp_cascade_con%initial_stock(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_metabolic(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_cellulose(0:ndecomp_pools))
    allocate(decomp_cascade_con%is_lignin(0:ndecomp_pools))
    allocate(decomp_cascade_con%spinup_factor(0:ndecomp_pools))

    !-- properties of each pathway along decomposition cascade 
    decomp_cascade_con%cascade_step_name(1:ndecomp_cascade_transitions) = ''
    decomp_cascade_con%cascade_donor_pool(1:ndecomp_cascade_transitions) = 0
    decomp_cascade_con%cascade_receiver_pool(1:ndecomp_cascade_transitions) = 0

    !-- first initialization of properties of each decomposing pool
    decomp_cascade_con%floating_cn_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%floating_cp_ratio_decomp_pools(0:ndecomp_pools) = .false.
    decomp_cascade_con%decomp_pool_name_history(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_restart(0:ndecomp_pools)       = ''
    decomp_cascade_con%decomp_pool_name_long(0:ndecomp_pools)          = ''
    decomp_cascade_con%decomp_pool_name_short(0:ndecomp_pools)         = ''
    decomp_cascade_con%is_litter(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%is_soil(0:ndecomp_pools)                        = .false.
    decomp_cascade_con%is_cwd(0:ndecomp_pools)                         = .false.
    decomp_cascade_con%initial_cn_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_cp_ratio(0:ndecomp_pools)               = nan
    decomp_cascade_con%initial_stock(0:ndecomp_pools)                  = nan
    decomp_cascade_con%is_metabolic(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_cellulose(0:ndecomp_pools)                   = .false.
    decomp_cascade_con%is_lignin(0:ndecomp_pools)                      = .false.
    decomp_cascade_con%spinup_factor(0:ndecomp_pools)                  = nan
  
  subroutine clean_col_pp(this)
    class(column_physical_properties) :: this
  
! !ARGUMENTS:
!    class(column_physical_properties) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell   )
    deallocate(this%wtgcell    )
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
    deallocate(this%soil_pH   )

  end subroutine clean_col_pp

end module ColumnPhysicalPropertiesType

