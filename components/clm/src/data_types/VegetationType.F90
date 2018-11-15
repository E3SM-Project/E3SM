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
  use clm_varcon     , only : ispval
  use clm_varctl     , only : use_fates
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
  ! Define the data structure that holds energy state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_es
    procedure, public :: Clean => clean_veg_es
  end type vegetation_energy_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ws
    procedure, public :: Clean => clean_veg_ws
  end type vegetation_water_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_cs
    procedure, public :: Clean => clean_veg_cs
  end type vegetation_carbon_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ns
    procedure, public :: Clean => clean_veg_ns
  end type vegetation_nitrogen_state
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus state information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_state
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ps
    procedure, public :: Clean => clean_veg_ps
  end type vegetation_phosphorus_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_energy_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_ef
    procedure, public :: Clean => clean_veg_ef
  end type vegetation_energy_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_water_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_wf
    procedure, public :: Clean => clean_veg_wf
  end type vegetation_water_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds carbon flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_carbon_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_cf
    procedure, public :: Clean => clean_veg_cf
  end type vegetation_carbon_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds nitrogen flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_nitrogen_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_nf
    procedure, public :: Clean => clean_veg_nf
  end type vegetation_nitrogen_flux
  
  !-----------------------------------------------------------------------
  ! Define the data structure that holds phosphorus flux information at the vegetation level.
  !-----------------------------------------------------------------------
  type, public :: vegetation_phosphorus_flux
    real(r8), pointer :: xxx      (:) => null() ! xxx (xxx)
  contains
    procedure, public :: Init  => init_veg_pf
    procedure, public :: Clean => clean_veg_pf
  end type vegetation_phosphorus_flux
  
  !-----------------------------------------------------------------------
  ! declare the public instances of vegetation-level data types
  !-----------------------------------------------------------------------
  type(vegetation_physical_properties)   , public, target :: veg_pp    ! vegetation physical properties
  type(vegetation_energy_state)          , public, target :: veg_es    ! vegetation energy state
  type(vegetation_water_state)           , public, target :: veg_ws    ! vegetation water state
  type(vegetation_carbon_state)          , public, target :: veg_cs    ! vegetation carbon state
  type(vegetation_nitrogen_state)        , public, target :: veg_ns    ! vegetation nitrogen state
  type(vegetation_phosphorus_state)      , public, target :: veg_ps    ! vegetation phosphorus state
  type(vegetation_energy_flux)           , public, target :: veg_ef    ! vegetation energy flux
  type(vegetation_water_flux)            , public, target :: veg_wf    ! vegetation water flux
  type(vegetation_carbon_flux)           , public, target :: veg_cf    ! vegetation carbon flux
  type(vegetation_nitrogen_flux)         , public, target :: veg_nf    ! vegetation nitrogen flux
  type(vegetation_phosphorus_flux)       , public, target :: veg_pf    ! vegetation phosphorus flux

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

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_es(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_es
    
  !------------------------------------------------------------------------
  subroutine clean_veg_es(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_es
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ws(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ws
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ws(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ws
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_cs(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_cs
    
  !------------------------------------------------------------------------
  subroutine clean_veg_cs(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_cs
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ns(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ns
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ns(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ns
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus state data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ps(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ps
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ps(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_state) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ps

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation energy flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_ef(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_ef
    
  !------------------------------------------------------------------------
  subroutine clean_veg_ef(this)
    !
    ! !ARGUMENTS:
    class(vegetation_energy_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_ef
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation water flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_wf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_wf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_wf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_water_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_wf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation carbon flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_cf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_cf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_cf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_carbon_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_cf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation nitrogen flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_nf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_nf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_nf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_nitrogen_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_nf
  
  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean vegetation phosphorus flux data structure
  !------------------------------------------------------------------------
  subroutine init_veg_pf(this, begc, endc)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    integer, intent(in) :: begc,endc
    !------------------------------------------------------------------------
    
  end subroutine init_veg_pf
    
  !------------------------------------------------------------------------
  subroutine clean_veg_pf(this)
    !
    ! !ARGUMENTS:
    class(vegetation_phosphorus_flux) :: this
    !------------------------------------------------------------------------
    
  end subroutine clean_veg_pf
  

end module VegetationType
