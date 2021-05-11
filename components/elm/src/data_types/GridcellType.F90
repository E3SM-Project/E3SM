module GridcellType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gridcell data type allocation 
  ! -------------------------------------------------------- 
  ! gridcell types can have values of 
  ! -------------------------------------------------------- 
  !   1 => default
  !
  ! PET: 9 Feb 2015: Preparing to change the sub-grid hierarchy to include
  ! 	 topographic units between gridcell and landunit.
  !	 
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use landunit_varcon, only : max_lunit
  use elm_varcon     , only : ispval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: gridcell_physical_properties_type

     ! topological mapping functionality, local 1d gdc arrays
     integer , pointer :: gindex       (:) => null() ! global index
     real(r8), pointer :: area         (:) => null() ! total land area, gridcell (km^2)
     real(r8), pointer :: lat          (:) => null() ! latitude (radians)
     real(r8), pointer :: lon          (:) => null() ! longitude (radians)
     real(r8), pointer :: latdeg       (:) => null() ! latitude (degrees)
     real(r8), pointer :: londeg       (:) => null() ! longitude (degrees)

     ! Starting and ending indices for all subgrid types below the gridcell level
     integer , pointer :: topi         (:) => null() ! beginning topographic unit index for each gridcell
     integer , pointer :: topf         (:) => null() ! ending topographic unit index for each gridcell
     integer , pointer :: ntopounits   (:) => null() ! number of topographic units for each gridcell
     integer , pointer :: lndi         (:) => null() ! beginning landunit index for each gridcell
     integer , pointer :: lndf         (:) => null() ! ending landunit index for each gridcell
     integer , pointer :: nlandunits   (:) => null() ! number of landunits for each gridcell
     integer , pointer :: coli         (:) => null() ! beginning column index for each gridcell
     integer , pointer :: colf         (:) => null() ! ending column index for each gridcell
     integer , pointer :: ncolumns     (:) => null() ! number of columns for each gridcell
     integer , pointer :: pfti         (:) => null() ! beginning pft index for each gridcell
     integer , pointer :: pftf         (:) => null() ! ending pft index for each gridcell
     integer , pointer :: npfts        (:) => null() ! number of patches for each gridcell

     ! Daylength
     real(r8) , pointer :: max_dayl    (:) => null() ! maximum daylength for this grid cell (s)
     real(r8) , pointer :: dayl        (:) => null() ! daylength (seconds)
     real(r8) , pointer :: prev_dayl   (:) => null() ! daylength from previous timestep (seconds)

     ! indices into landunit-level arrays for landunits in this grid cell (ispval implies
     ! this landunit doesn't exist on this grid cell) [1:max_lunit, begg:endg]
     ! (note that the spatial dimension is last here, in contrast to most 2-d variables;
     ! this is for efficiency, since most loops will go over g in the outer loop, and
     ! landunit type in the inner loop)
     integer , pointer :: landunit_indices (:,:) => null() 

   contains

     procedure, public :: Init => grc_pp_init
     procedure, public :: Clean => grc_pp_clean
     
  end type gridcell_physical_properties_type
  type(gridcell_physical_properties_type), public, target :: grc_pp    !gridcell data structure
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine grc_pp_init(this, begg, endg)
    !
    ! !ARGUMENTS:
    class(gridcell_physical_properties_type) :: this
    integer, intent(in)  :: begg, endg
    !------------------------------------------------------------------------

    ! The following is set in InitGridCells
    allocate(this%gindex    (begg:endg)) ; this%gindex    (:) = ispval
    allocate(this%area      (begg:endg)) ; this%area      (:) = nan
    allocate(this%lat       (begg:endg)) ; this%lat       (:) = nan
    allocate(this%lon       (begg:endg)) ; this%lon       (:) = nan
    allocate(this%latdeg    (begg:endg)) ; this%latdeg    (:) = nan
    allocate(this%londeg    (begg:endg)) ; this%londeg    (:) = nan
    
    allocate(this%topi      (begg:endg)) ; this%topi      (:) = ispval
    allocate(this%topf      (begg:endg)) ; this%topf      (:) = ispval
    allocate(this%ntopounits(begg:endg)) ; this%ntopounits(:) = ispval
    allocate(this%lndi      (begg:endg)) ; this%lndi      (:) = ispval
    allocate(this%lndf      (begg:endg)) ; this%lndf      (:) = ispval
    allocate(this%nlandunits(begg:endg)) ; this%nlandunits(:) = ispval
    allocate(this%coli      (begg:endg)) ; this%coli      (:) = ispval
    allocate(this%colf      (begg:endg)) ; this%colf      (:) = ispval
    allocate(this%ncolumns  (begg:endg)) ; this%ncolumns  (:) = ispval
    allocate(this%pfti      (begg:endg)) ; this%pfti      (:) = ispval
    allocate(this%pftf      (begg:endg)) ; this%pftf      (:) = ispval
    allocate(this%npfts     (begg:endg)) ; this%npfts     (:) = ispval

    ! This is initiailized in module DayLength
    allocate(this%max_dayl  (begg:endg)) ; this%max_dayl  (:) = nan
    allocate(this%dayl      (begg:endg)) ; this%dayl      (:) = nan
    allocate(this%prev_dayl (begg:endg)) ; this%prev_dayl (:) = nan

    allocate(this%landunit_indices(1:max_lunit, begg:endg)); this%landunit_indices(:,:) = ispval

  end subroutine grc_pp_init

  !------------------------------------------------------------------------
  subroutine grc_pp_clean(this)
    !
    ! !ARGUMENTS:
    class(gridcell_physical_properties_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gindex           )
    deallocate(this%area             )
    deallocate(this%lat              )
    deallocate(this%lon              )
    deallocate(this%latdeg           )
    deallocate(this%londeg           )
    deallocate(this%topi             )
    deallocate(this%topf             )
    deallocate(this%ntopounits       )
    deallocate(this%lndi             )
    deallocate(this%lndf             )
    deallocate(this%nlandunits       )
    deallocate(this%coli             )
    deallocate(this%colf             )
    deallocate(this%ncolumns         )
    deallocate(this%pfti             )
    deallocate(this%pftf             )
    deallocate(this%npfts            )
    deallocate(this%max_dayl         )
    deallocate(this%dayl             )
    deallocate(this%prev_dayl        )
    deallocate(this%landunit_indices )
    
  end subroutine grc_pp_clean

end module GridcellType
