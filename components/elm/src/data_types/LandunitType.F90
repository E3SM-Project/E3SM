module LandunitType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Landunit data type allocation
  ! --------------------------------------------------------
  ! landunits types can have values of (see landunit_varcon.F90)
  ! --------------------------------------------------------
  !   1  => (istsoil)    soil (vegetated or bare soil landunit)
  !   2  => (istcrop)    crop (only for crop configuration)
  !   3  => (istice)     land ice
  !   4  => (istice_mec) land ice (multiple elevation classes)
  !   5  => (istdlak)    deep lake
  !   6  => (istwet)     wetland
  !   7  => (isturb_tbd) urban tbd
  !   8  => (isturb_hd)  urban hd
  !   9  => (isturb_md)  urban md
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use elm_varcon     , only : ispval, spval
  use elm_varpar     , only : nh3dc_per_lunit
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: landunit_physical_properties

     ! indices and weights for higher subgrid levels (topounit, gridcell)
     integer , pointer :: gridcell     (:) => null() ! index into gridcell level quantities
     real(r8), pointer :: wtgcell      (:) => null() ! weight (relative to gridcell)
     integer , pointer :: topounit     (:) => null() ! index into topounit level quantities
     real(r8), pointer :: wttopounit   (:) => null() ! weight (relative to topounit)

     ! Starting and ending indices for all subgrid types below the landunit level
     integer , pointer :: coli         (:) => null() ! beginning column index per landunit
     integer , pointer :: colf         (:) => null() ! ending column index for each landunit
     integer , pointer :: ncolumns     (:) => null() ! number of columns for each landunit
     integer , pointer :: pfti         (:) => null() ! beginning pft index for each landunit
     integer , pointer :: pftf         (:) => null() ! ending pft index for each landunit
     integer , pointer :: npfts        (:) => null() ! number of patches for each landunit

     ! topological mapping functionality
     integer , pointer :: itype        (:) => null() ! landunit type
     logical , pointer :: ifspecial    (:) => null() ! true=>landunit is not vegetated
     logical , pointer :: lakpoi       (:) => null() ! true=>lake point
     logical , pointer :: urbpoi       (:) => null() ! true=>urban point
     logical , pointer :: glcmecpoi    (:) => null() ! true=>glacier_mec point
     logical , pointer :: active       (:) => null() ! true=>do computations on this landunit
     logical , pointer :: ispolygon    (:) => null() ! true=>is polygonal tundra
     integer , pointer :: polygontype  (:) => null() ! ice wedge polygon initial condition (LCP, FCP, or HCP)

     ! urban properties
     real(r8), pointer :: canyon_hwr   (:) => null() ! urban landunit canyon height to width ratio (-)
     real(r8), pointer :: wtroad_perv  (:) => null() ! urban landunit weight of pervious road column to total road (-)
     real(r8), pointer :: wtlunit_roof (:) => null() ! weight of roof with respect to urban landunit (-)
     real(r8), pointer :: ht_roof      (:) => null() ! height of urban roof (m)
     real(r8), pointer :: z_0_town     (:) => null() ! urban landunit momentum roughness length (m)
     real(r8), pointer :: z_d_town     (:) => null() ! urban landunit displacement height (m)

     real(r8), pointer :: hs_w_itf     (:,:) => null() ! hillslope width function defined at h3D soil column interface (N+1 values for N soil columns) (m)
     real(r8), pointer :: hs_x_itf     (:,:) => null() ! hillslope width function defined at h3D soil column interface (N+1 values for N soil columns) (m)
     real(r8), pointer :: hs_w_nod     (:,:) => null() ! hillslope width function defined at h3D soil column node (N+1 values for N soil columns) (m)
     real(r8), pointer :: hs_x_nod     (:,:) => null() ! hillslope width function defined at h3D soil column node (N+1 values for N soil columns) (m)

     real(r8), pointer :: hs_dx        (:,:) => null() ! dx of h3D soil column (m)
     real(r8), pointer :: hs_dx_node   (:,:) => null() ! dx between 2 adjunct h3D nodes (m)
     real(r8), pointer :: hs_dA        (:,:) => null() ! area of h3D soil column (m)
     real(r8), pointer :: hs_area      (:) => null() ! area of h3D hillslope (m)

   contains

     procedure, public :: Init => lun_pp_init
     procedure, public :: Clean => lun_pp_clean

  end type landunit_physical_properties

  type(landunit_physical_properties), public, target :: lun_pp  !geomorphological landunits
  !$acc declare create(lun_pp)
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lun_pp_Init(this, begl, endl)
    !
    ! !ARGUMENTS:
    class(landunit_physical_properties) :: this
    integer, intent(in) :: begl,endl
    !------------------------------------------------------------------------

    ! The following is set in InitGridCellsMod
    allocate(this%gridcell     (begl:endl)); this%gridcell  (:) = ispval
    allocate(this%wtgcell      (begl:endl)); this%wtgcell   (:) = spval
    allocate(this%topounit     (begl:endl)); this%topounit  (:) = ispval
    allocate(this%wttopounit   (begl:endl)); this%wttopounit(:) = spval
    allocate(this%coli         (begl:endl)); this%coli      (:) = ispval
    allocate(this%colf         (begl:endl)); this%colf      (:) = ispval
    allocate(this%ncolumns     (begl:endl)); this%ncolumns  (:) = ispval
    allocate(this%pfti         (begl:endl)); this%pfti      (:) = ispval
    allocate(this%pftf         (begl:endl)); this%pftf      (:) = ispval
    allocate(this%npfts        (begl:endl)); this%npfts     (:) = ispval
    allocate(this%itype        (begl:endl)); this%itype     (:) = ispval
    allocate(this%ifspecial    (begl:endl)); this%ifspecial (:) = .false.
    allocate(this%lakpoi       (begl:endl)); this%lakpoi    (:) = .false.
    allocate(this%urbpoi       (begl:endl)); this%urbpoi    (:) = .false.
    allocate(this%glcmecpoi    (begl:endl)); this%glcmecpoi (:) = .false.
    allocate(this%ispolygon    (begl:endl)); this%ispolygon (:) = .false.
    allocate(this%polygontype  (begl:endl)); this%polygontype(:) = ispval

    ! The following is initialized in routine setActive in module reweightMod
    allocate(this%active       (begl:endl))

    ! The following is set in routine urbanparams_vars%Init in module UrbanParamsMod
    allocate(this%canyon_hwr   (begl:endl)); this%canyon_hwr   (:) = spval
    allocate(this%wtroad_perv  (begl:endl)); this%wtroad_perv  (:) = spval
    allocate(this%wtlunit_roof (begl:endl)); this%wtlunit_roof (:) = spval
    allocate(this%ht_roof      (begl:endl)); this%ht_roof      (:) = spval
    allocate(this%z_0_town     (begl:endl)); this%z_0_town     (:) = spval
    allocate(this%z_d_town     (begl:endl)); this%z_d_town     (:) = spval

    ! The following is set in initVerticalMod
    allocate(this%hs_x_itf     (begl:endl,1:nh3dc_per_lunit+1)); this%hs_x_itf(:,:) = spval
    allocate(this%hs_w_itf     (begl:endl,1:nh3dc_per_lunit+1)); this%hs_w_itf(:,:) = spval
    allocate(this%hs_x_nod     (begl:endl,1:nh3dc_per_lunit  )); this%hs_x_nod(:,:) = spval
    allocate(this%hs_w_nod     (begl:endl,1:nh3dc_per_lunit  )); this%hs_w_nod(:,:) = spval
    allocate(this%hs_dx        (begl:endl,1:nh3dc_per_lunit  )); this%hs_dx   (:,:) = spval
    allocate(this%hs_dx_node   (begl:endl,1:nh3dc_per_lunit  )); this%hs_dx_node(:,:) = spval
    allocate(this%hs_dA        (begl:endl,1:nh3dc_per_lunit  )); this%hs_dA   (:,:) = spval
    allocate(this%hs_area      (begl:endl                    )); this%hs_area (:) = spval

  end subroutine lun_pp_init

  !------------------------------------------------------------------------
  subroutine lun_pp_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_physical_properties) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell     )
    deallocate(this%wtgcell      )
    deallocate(this%topounit     )
    deallocate(this%wttopounit   )
    deallocate(this%coli         )
    deallocate(this%colf         )
    deallocate(this%ncolumns     )
    deallocate(this%pfti         )
    deallocate(this%pftf         )
    deallocate(this%npfts        )
    deallocate(this%itype        )
    deallocate(this%ifspecial    )
    deallocate(this%lakpoi       )
    deallocate(this%urbpoi       )
    deallocate(this%glcmecpoi    )
    deallocate(this%ispolygon    )
    deallocate(this%polygontype  )
    deallocate(this%active       )
    deallocate(this%canyon_hwr   )
    deallocate(this%wtroad_perv  )
    deallocate(this%ht_roof      )
    deallocate(this%wtlunit_roof )
    deallocate(this%z_0_town     )
    deallocate(this%z_d_town     )
    deallocate(this%hs_x_itf     )
    deallocate(this%hs_w_itf     )
    deallocate(this%hs_x_nod     )
    deallocate(this%hs_dx        )
    deallocate(this%hs_dx_node   )
    deallocate(this%hs_dA        )
    deallocate(this%hs_area      )

  end subroutine lun_pp_clean

end module LandunitType
