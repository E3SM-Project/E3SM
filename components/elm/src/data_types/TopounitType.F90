module TopounitType

  ! --------------------------------------------------------
  ! ELM sub-grid hierarchy:
  ! Define topographic unit data types, with Init and Clean for each
  ! Init() calls allocate memory and set history fields
  ! --------------------------------------------------------
  ! 3 Aug 2015, PET
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use landunit_varcon, only : max_lunit
  use elm_varcon     , only : ispval, spval
  use elm_varpar     , only : numrad
  use decompMod      , only : bounds_type

  implicit none
  save
  private

  !-----------------------------------------------------------------------
  ! sub-grid topology and physical properties defined at the topographic unit level
  type, public :: topounit_physical_properties

    ! indices and weights for higher subgrid level (gridcell)
    integer , pointer :: gridcell    (:) => null() ! index into gridcell level quantities
    integer , pointer :: topo_grc_ind(:) => null() ! index of topounit in the grid
    real(r8), pointer :: wtgcell     (:) => null() ! weight (relative to gridcell)

    ! connections between topounits in a single gridcell, for IM2 hillslope hydrology
    ! downhill_ti value of -1 indicates that there is no downhill topounit
    ! uphill_wt is a temporary store for the sum of column weights receiving water from uphill
    integer , pointer :: downhill_ti (:) => null() ! topounit index for downhill topounit on this gridcell
    real(r8), pointer :: uphill_wt   (:) => null() ! sum of column weights on topounit in the hydrologyc filter

    ! Starting and ending indices for all subgrid types below the landunit level
    integer , pointer :: lndi        (:) => null() ! beginning landunit index for each topounit
    integer , pointer :: lndf        (:) => null() ! ending landunit index for each topounit
    integer , pointer :: nlandunits  (:) => null() ! number of landunits for each topounit
    integer , pointer :: coli        (:) => null() ! beginning column index per landunit
    integer , pointer :: colf        (:) => null() ! ending column index for each landunit
    integer , pointer :: ncolumns    (:) => null() ! number of columns for each landunit
    integer , pointer :: pfti        (:) => null() ! beginning pft index for each landunit
    integer , pointer :: pftf        (:) => null() ! ending pft index for each landunit
    integer , pointer :: npfts       (:) => null() ! number of patches for each landunit

    ! indices into landunit-level arrays for landunits in this topounit (ispval implies
    ! this landunit doesn't exist on this topounit) [1:max_lunit, begt:endt]
    ! (note that the spatial dimension is last here, in contrast to most 2-d variables;
    ! this is for efficiency, since most loops will go over t in the outer loop, and
    ! landunit type in the inner loop)
    integer , pointer :: landunit_indices (:,:) => null() 
    logical , pointer :: active       (:) => null() ! true=>do computations on this topounit 
    ! physical properties
    real(r8), pointer :: area       (:) => null() ! land area (km^2)
    real(r8), pointer :: lat        (:) => null() ! mean latitude (radians)
    real(r8), pointer :: lon        (:) => null() ! mean longitude (radians)
    real(r8), pointer :: elevation  (:) => null() ! mean soil surface elevation, above mean sea level (m)
    real(r8), pointer :: slope      (:) => null() ! mean slope angle (radians)
    integer , pointer :: aspect     (:) => null() ! mean aspect angle, measured clockwise from north (radians)
    real(r8), pointer :: emissivity (:) => null() ! mean surface emissivity
    real(r8), pointer :: surfalb_dir(:,:) => null() ! (topunit,numrad) mean surface albedo (direct)
    real(r8), pointer :: surfalb_dif(:,:) => null() ! (topunit,numrad) mean surface albedo (diffuse)
  contains
    procedure, public :: Init  => init_top_pp
    procedure, public :: Clean => clean_top_pp
  end type topounit_physical_properties

  !-----------------------------------------------------------------------
  ! declare the public instances of topounit physical property type
  type(topounit_physical_properties),  public, target :: top_pp
  !$acc declare create(top_pp)
  contains

  !-----------------------------------------------------------------------
  subroutine init_top_pp(this, begt, endt)
    class(topounit_physical_properties) :: this
    integer, intent(in) :: begt   ! beginning topographic unit index
    integer, intent(in) :: endt   ! ending topographic unit index

    allocate(this%gridcell    (begt:endt)) ; this%gridcell    (:) = ispval
    allocate(this%topo_grc_ind(begt:endt)) ; this%topo_grc_ind(:) = ispval
    allocate(this%wtgcell     (begt:endt)) ; this%wtgcell     (:) = spval
    allocate(this%downhill_ti (begt:endt)) ; this%downhill_ti (:) = ispval
    allocate(this%uphill_wt   (begt:endt)) ; this%uphill_wt   (:) = ispval
    allocate(this%lndi        (begt:endt)) ; this%lndi        (:) = ispval
    allocate(this%lndf        (begt:endt)) ; this%lndf        (:) = ispval
    allocate(this%nlandunits  (begt:endt)) ; this%nlandunits  (:) = ispval
    allocate(this%coli        (begt:endt)) ; this%coli        (:) = ispval
    allocate(this%colf        (begt:endt)) ; this%colf        (:) = ispval
    allocate(this%ncolumns    (begt:endt)) ; this%ncolumns    (:) = ispval
    allocate(this%pfti        (begt:endt)) ; this%pfti        (:) = ispval
    allocate(this%pftf        (begt:endt)) ; this%pftf        (:) = ispval
    allocate(this%npfts       (begt:endt)) ; this%npfts       (:) = ispval

    allocate(this%landunit_indices(1:max_lunit, begt:endt)); this%landunit_indices(:,:) = ispval
    allocate(this%active      (begt:endt))                     ; this%active      (:)   = .false.

    allocate(this%area        (begt:endt)) ; this%area        (:) = spval 
    allocate(this%lat         (begt:endt)) ; this%lat         (:) = spval
    allocate(this%lon         (begt:endt)) ; this%lon         (:) = spval
    allocate(this%elevation   (begt:endt)) ; this%elevation   (:) = spval
    allocate(this%slope       (begt:endt)) ; this%slope       (:) = spval
    allocate(this%aspect      (begt:endt)) ; this%aspect      (:) = ispval
    allocate(this%emissivity  (begt:endt)) ; this%emissivity  (:) = spval
    allocate(this%surfalb_dir (begt:endt,1:numrad)) ; this%surfalb_dir(:,:) = spval
    allocate(this%surfalb_dif (begt:endt,1:numrad)) ; this%surfalb_dif(:,:) = spval 
  end subroutine init_top_pp

  !-----------------------------------------------------------------------
  subroutine clean_top_pp(this)
    class(topounit_physical_properties) :: this

    deallocate(this%gridcell    )
    deallocate(this%topo_grc_ind    )
    deallocate(this%wtgcell     )
    deallocate(this%downhill_ti )
    deallocate(this%uphill_wt   )
    deallocate(this%lndi        )
    deallocate(this%lndf        )
    deallocate(this%nlandunits  )
    deallocate(this%coli        )
    deallocate(this%colf        )
    deallocate(this%ncolumns    )
    deallocate(this%pfti        )
    deallocate(this%pftf        )
    deallocate(this%npfts       )
    deallocate(this%landunit_indices )
    deallocate(this%active     )

    deallocate(this%area        )
    deallocate(this%lat         )
    deallocate(this%lon         )
    deallocate(this%elevation   )
    deallocate(this%slope       )
    deallocate(this%aspect      )
    deallocate(this%emissivity  )
    deallocate(this%surfalb_dir )
    deallocate(this%surfalb_dif )
  end subroutine clean_top_pp

end module TopounitType
