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
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varcon     , only : ispval
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

     ! urban properties
     real(r8), pointer :: canyon_hwr   (:) => null() ! urban landunit canyon height to width ratio (-)   
     real(r8), pointer :: wtroad_perv  (:) => null() ! urban landunit weight of pervious road column to total road (-)
     real(r8), pointer :: wtlunit_roof (:) => null() ! weight of roof with respect to urban landunit (-)
     real(r8), pointer :: ht_roof      (:) => null() ! height of urban roof (m)
     real(r8), pointer :: z_0_town     (:) => null() ! urban landunit momentum roughness length (m)
     real(r8), pointer :: z_d_town     (:) => null() ! urban landunit displacement height (m)

   contains

     procedure, public :: Init => lun_pp_init
     procedure, public :: Clean => lun_pp_clean
     
  end type landunit_physical_properties
  type(landunit_physical_properties), public, target :: lun_pp  !geomorphological landunits
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
    allocate(this%wtgcell      (begl:endl)); this%wtgcell   (:) = nan
    allocate(this%topounit     (begl:endl)); this%topounit  (:) = ispval
    allocate(this%wttopounit   (begl:endl)); this%wttopounit(:) = nan
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

    ! The following is initialized in routine setActive in module reweightMod
    allocate(this%active       (begl:endl))

    ! The following is set in routine urbanparams_vars%Init in module UrbanParamsMod
    allocate(this%canyon_hwr   (begl:endl)); this%canyon_hwr   (:) = nan
    allocate(this%wtroad_perv  (begl:endl)); this%wtroad_perv  (:) = nan
    allocate(this%wtlunit_roof (begl:endl)); this%wtlunit_roof (:) = nan
    allocate(this%ht_roof      (begl:endl)); this%ht_roof      (:) = nan
    allocate(this%z_0_town     (begl:endl)); this%z_0_town     (:) = nan
    allocate(this%z_d_town     (begl:endl)); this%z_d_town     (:) = nan

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
    deallocate(this%active       )
    deallocate(this%canyon_hwr   )
    deallocate(this%wtroad_perv  )
    deallocate(this%ht_roof      )
    deallocate(this%wtlunit_roof )
    deallocate(this%z_0_town     )
    deallocate(this%z_d_town     )

  end subroutine lun_pp_clean

end module LandunitType
