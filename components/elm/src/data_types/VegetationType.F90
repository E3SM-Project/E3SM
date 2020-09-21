module VegetationType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Vegetation data type allocation 
  ! -------------------------------------------------------- 
  ! Vegetation types can have values of
  ! -------------------------------------------------------- 
  !   0  => not vegetated
  !   1  => needleleaf evergreen temperate tree
  !   2  => needleleaf evergreen boreal tree
  !   3  => needleleaf deciduous boreal tree
  !   4  => broadleaf evergreen tropical tree
  !   5  => broadleaf evergreen temperate tree
  !   6  => broadleaf deciduous tropical tree
  !   7  => broadleaf deciduous temperate tree
  !   8  => broadleaf deciduous boreal tree
  !   9  => broadleaf evergreen shrub
  !   10 => broadleaf deciduous temperate shrub
  !   11 => broadleaf deciduous boreal shrub
  !   12 => c3 arctic grass
  !   13 => c3 non-arctic grass
  !   14 => c4 grass
  !   15 => c3_crop
  !   16 => c3_irrigated
  !   17 => corn
  !   18 => irrigated corn
  !   19 => spring temperate cereal
  !   20 => irrigated spring temperate cereal
  !   21 => winter temperate cereal
  !   22 => irrigated winter temperate cereal
  !   23 => soybean
  !   24 => irrigated soybean
  ! -------------------------------------------------------- 
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varcon     , only : ispval
  use elm_varctl     , only : use_fates
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !-----------------------------------------------------------------------
  ! Define the data structure that holds physical property information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_physical_properties
     ! indices and weights for higher subgrid levels (column, landunit, topounit, gridcell)
     integer , pointer :: gridcell      (:) => null() ! index into gridcell level quantities
     real(r8), pointer :: wtgcell       (:) => null() ! weight (relative to gridcell) 
     integer , pointer :: topounit      (:) => null() ! index into topounit level quantities
     real(r8), pointer :: wttopounit    (:) => null() ! weight (relative to topounit)
     integer , pointer :: landunit      (:) => null() ! index into landunit level quantities
     real(r8), pointer :: wtlunit       (:) => null() ! weight (relative to landunit) 
     integer , pointer :: column        (:) => null() ! index into column level quantities
     real(r8), pointer :: wtcol         (:) => null() ! weight (relative to column) 

     ! topological mapping functionality
     integer , pointer :: itype         (:) => null() ! patch vegetation 
     integer , pointer :: mxy           (:) => null() ! m index for laixy(i,j,m),etc. (undefined for special landunits)
     logical , pointer :: active        (:) => null() ! true=>do computations on this patch

     ! Fates relevant types
     logical , pointer :: is_veg        (:) => null() ! This is an ACTIVE fates patch
     logical , pointer :: is_bareground (:) => null() ! ?
     real(r8), pointer :: wt_ed         (:) => null() ! TODO mv ? can this be removed
     logical , pointer :: is_fates      (:) => null() ! true for patch vector space reserved
                                                      ! for FATES.
                                                      ! this is static and is true for all 
                                                      ! patches within fates jurisdiction
                                                      ! including patches which are not currently
                                                      ! associated with a FATES linked-list patch
   contains

     procedure, public :: Init => veg_pp_init
     procedure, public :: Clean => veg_pp_clean
     
  end type vegetation_physical_properties

  !-----------------------------------------------------------------------
  ! declare the public instance of vegetation physical property data types
  !-----------------------------------------------------------------------
  type(vegetation_physical_properties)   , public, target :: veg_pp    ! vegetation physical properties

  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine veg_pp_init(this, begp, endp)
    !
    ! !ARGUMENTS:
    class(vegetation_physical_properties)   :: this
    integer, intent(in) :: begp,endp
    !
    ! LOCAL VARAIBLES:
    !------------------------------------------------------------------------

    ! The following is set in InitGridCells
    allocate(this%gridcell  (begp:endp)); this%gridcell    (:) = ispval
    allocate(this%wtgcell   (begp:endp)); this%wtgcell     (:) = nan
    allocate(this%topounit  (begp:endp)); this%topounit    (:) = ispval
    allocate(this%wttopounit(begp:endp)); this%wttopounit  (:) = nan
    allocate(this%landunit  (begp:endp)); this%landunit    (:) = ispval
    allocate(this%wtlunit   (begp:endp)); this%wtlunit     (:) = nan
    allocate(this%column    (begp:endp)); this%column      (:) = ispval
    allocate(this%wtcol     (begp:endp)); this%wtcol       (:) = nan
    allocate(this%itype     (begp:endp)); this%itype       (:) = ispval
    allocate(this%mxy       (begp:endp)); this%mxy         (:) = ispval
    allocate(this%active    (begp:endp)); this%active      (:) = .false.

    allocate(this%is_fates   (begp:endp)); this%is_fates   (:) = .false.
    if (use_fates) then
       allocate(this%is_veg  (begp:endp)); this%is_veg  (:) = .false.
       allocate(this%is_bareground (begp:endp)); this%is_bareground (:) = .false.
       allocate(this%wt_ed      (begp:endp)); this%wt_ed      (:) = nan 
    end if

	end subroutine veg_pp_init

  !------------------------------------------------------------------------
  subroutine veg_pp_clean(this)
    !
    ! !ARGUMENTS:
    class(vegetation_physical_properties) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell  )
    deallocate(this%wtgcell   )
    deallocate(this%topounit  )
    deallocate(this%wttopounit)
    deallocate(this%landunit  )
    deallocate(this%wtlunit   )
    deallocate(this%column    )
    deallocate(this%wtcol     )
    deallocate(this%itype     )
    deallocate(this%mxy       )
    deallocate(this%active    )

	 deallocate(this%is_fates)
    if (use_fates) then
       deallocate(this%is_veg)
       deallocate(this%is_bareground)
       deallocate(this%wt_ed)
    end if

  end subroutine veg_pp_clean

end module VegetationType
