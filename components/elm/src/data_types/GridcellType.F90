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
  !    topographic units between gridcell and landunit.
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use landunit_varcon, only : max_lunit
  use elm_varcon     , only : ispval, spval
  use topounit_varcon, only : max_topounits
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
     integer , pointer :: ntopounits   (:) => null() ! number of topographic units for each gridcell calculated from topi and topf
     integer , pointer :: lndi         (:) => null() ! beginning landunit index for each gridcell
     integer , pointer :: lndf         (:) => null() ! ending landunit index for each gridcell
     integer , pointer :: nlandunits   (:) => null() ! number of landunits for each gridcell
     integer , pointer :: coli         (:) => null() ! beginning column index for each gridcell
     integer , pointer :: colf         (:) => null() ! ending column index for each gridcell
     integer , pointer :: ncolumns     (:) => null() ! number of columns for each gridcell
     integer , pointer :: pfti         (:) => null() ! beginning pft index for each gridcell
     integer , pointer :: pftf         (:) => null() ! ending pft index for each gridcell
     integer , pointer :: npfts        (:) => null() ! number of patches for each gridcell

     real(r8), pointer :: stdev_elev   (:) => null()     ! standard deviation of elevation within a gridcell
     real(r8), pointer :: sky_view     (:) => null()     ! mean of (sky view factor / cos(slope))
     real(r8), pointer :: terrain_config (:) => null()   ! mean of (terrain configuration factor / cos(slope))
     real(r8), pointer :: sinsl_cosas  (:) => null()     ! sin(slope)*cos(aspect) / cos(slope)
     real(r8), pointer :: sinsl_sinas  (:) => null()     ! sin(slope)*sin(aspect) / cos(slope)
     
     ! Daylength
     real(r8) , pointer :: max_dayl    (:) => null() ! maximum daylength for this grid cell (seconds)
     real(r8) , pointer :: dayl        (:) => null() ! daylength (seconds)
     real(r8) , pointer :: prev_dayl   (:) => null() ! daylength from previous timestep (seconds)
     real(r8) , pointer :: elevation   (:) => null() ! mean soil surface elevation, above mean sea level (m)
     real(r8) , pointer :: froudenum   (:) => null() ! Froude number (dimensionless)
     real(r8) , pointer :: MaxElevation   (:) => null() ! Maximum soil surface elevation, above mean sea level (meter) needed for precipitation downscaling
 
     !integer , pointer :: topounit_indices (:,:) => null()

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
  !$acc declare create(grc_pp)
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
    allocate(this%area      (begg:endg)) ; this%area      (:) = spval
    allocate(this%lat       (begg:endg)) ; this%lat       (:) = spval
    allocate(this%lon       (begg:endg)) ; this%lon       (:) = spval
    allocate(this%latdeg    (begg:endg)) ; this%latdeg    (:) = spval
    allocate(this%londeg    (begg:endg)) ; this%londeg    (:) = spval

    allocate(this%topi      (begg:endg)) ; this%topi      (:) = ispval
    allocate(this%topf      (begg:endg)) ; this%topf      (:) = ispval
    allocate(this%ntopounits(begg:endg)) ; this%ntopounits(:) = 1      ! Default number of topounits per grid is 1
    allocate(this%lndi      (begg:endg)) ; this%lndi      (:) = ispval
    allocate(this%lndf      (begg:endg)) ; this%lndf      (:) = ispval
    allocate(this%nlandunits(begg:endg)) ; this%nlandunits(:) = ispval
    allocate(this%coli      (begg:endg)) ; this%coli      (:) = ispval
    allocate(this%colf      (begg:endg)) ; this%colf      (:) = ispval
    allocate(this%ncolumns  (begg:endg)) ; this%ncolumns  (:) = ispval
    allocate(this%pfti      (begg:endg)) ; this%pfti      (:) = ispval
    allocate(this%pftf      (begg:endg)) ; this%pftf      (:) = ispval
    allocate(this%npfts     (begg:endg)) ; this%npfts     (:) = ispval

    allocate(this%stdev_elev(begg:endg)) ; this%stdev_elev(:) = ispval        ! standard deviation of elevation within a gridcell
    allocate(this%sky_view  (begg:endg)) ; this%sky_view  (:) = ispval        ! mean of (sky view factor / cos(slope))
    allocate(this%terrain_config(begg:endg)) ; this%terrain_config(:) = ispval! mean of (terrain configuration factor / cos(slope))
    allocate(this%sinsl_cosas(begg:endg)) ; this%sinsl_cosas(:) = ispval      ! sin(slope)*cos(aspect) / cos(slope)
    allocate(this%sinsl_sinas(begg:endg)) ; this%sinsl_sinas(:) = ispval      ! sin(slope)*sin(aspect) / cos(slope)
    
    ! This is initiailized in module DayLength
    allocate(this%max_dayl  (begg:endg)) ; this%max_dayl  (:) = spval
    allocate(this%dayl      (begg:endg)) ; this%dayl      (:) = spval
    allocate(this%prev_dayl (begg:endg)) ; this%prev_dayl (:) = spval
    
    allocate(this%elevation (begg:endg)) ; this%elevation (:) = spval
    allocate(this%froudenum (begg:endg)) ; this%froudenum (:) = spval 
    allocate(this%MaxElevation (begg:endg)) ; this%MaxElevation (:) = spval

    allocate(this%landunit_indices(1:max_lunit, begg:endg)); this%landunit_indices(:,:) = ispval
   
   ! allocate(this%topounit_indices (begg:endg,1:max_topounits)) ; this%topounit_indices (:,:) = ispval

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
    deallocate(this%elevation        )
    deallocate(this%froudenum        )
    deallocate(this%MaxElevation     )
    deallocate(this%landunit_indices )
    deallocate(this%stdev_elev       ) 
    deallocate(this%sky_view         ) 
    deallocate(this%terrain_config   ) 
    deallocate(this%sinsl_cosas      )
    deallocate(this%sinsl_sinas      )
    
  end subroutine grc_pp_clean

end module GridcellType
