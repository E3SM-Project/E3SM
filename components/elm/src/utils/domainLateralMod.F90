module domainLateralMod

#include "shr_assert.h"


#ifdef MOAB_LATERAL

  !-----------------------------------------------------------------------
  ! This is a stub for the case when PETSc is unavailable
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use elm_varctl  , only : iulog
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils  , only : endrun
  use MOABGridType, only : moab_gcell, mlndghostid
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !

  type, public :: oneD_int_data_for_moab
     integer               :: moab_app_id    ! ID of MAOB app
     character(len=1024)   :: tag_name       ! MOAB tag name
     integer               :: tag_type       ! type of MOAB tag: 0 = dense, int; 1 = dense, double
     integer               :: num_tags       ! Number of tags
     integer               :: entity_type(1) ! vertex or element based type
     integer               :: tag_index(1)   ! Index of tag after it is registered in MOAB
     integer               :: num_comp       ! number of components
     integer               :: ngcells        ! number of grid cells
     integer, allocatable  :: values(:)      ! data
  end type oneD_int_data_for_moab

  type, public :: twoD_real_data_for_moab
     integer               :: moab_app_id    ! ID of MAOB app
     character(len=1024)   :: tag_name       ! MOAB tag name
     integer               :: tag_type       ! type of MOAB tag: 0 = dense, int; 1 = dense, double
     integer               :: num_tags       ! Number of tags
     integer               :: entity_type(1) ! vertex or element based type
     integer               :: tag_index(1)   ! Index of tag after it is registered in MOAB
     integer               :: num_comp       ! number of components
     integer               :: ngcells        ! number of grid cells
     real(r8), allocatable :: values(:,:)    ! data
  end type twoD_real_data_for_moab

  type, public :: domainlateral_type
     type(oneD_int_data_for_moab)  :: grid_level_count
  end type domainlateral_type

  type(domainlateral_type)    , public :: ldomain_lateral
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public domainlateral_init                 ! initializes
  public setup_twoD_real_data_for_moab
  public GridLevelIntegerDataHaloExchange
  public GridLevelRealDataHaloExchange
  !
  !EOP
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine setup_oneD_int_data_for_moab(moab_app_id, tag_name, num_cells_ghosted, data)
    !
    ! DESCRIPTION:
    ! Sets up 1D integer-type data structure for performing MOAB-based halo exchange
    ! of data. This supports a single value per grid cell to be exchanged.
    !
    use iso_c_binding
    use iMOAB, only : iMOAB_DefineTagStorage
    !
    implicit none
    !
    ! ARGUMENTS:
    integer                      , intent(in)  :: moab_app_id
    character(len=*)             , intent(in)  :: tag_name
    integer                      , intent(in)  :: num_cells_ghosted
    type(oneD_int_data_for_moab) , intent(out) :: data
    !
    ! LOCAL VARIABLES:
    integer                                :: ierr

    data%moab_app_id    = moab_app_id
    data%tag_name       = trim(tag_name) // C_NULL_CHAR ! name
    data%tag_type       = 0                             ! 0 = dense, int
    data%num_tags       = 1                             ! a single tag
    data%entity_type(1) = 1                             ! element (== cell) based data
    data%num_comp       = 1                             ! number components in the tag
    data%ngcells        = num_cells_ghosted             ! number of grid cells

    ! allocate memory
    allocate(data%values(data%ngcells))

    ! define the tag in MOAB
    ierr = iMOAB_DefineTagStorage(data%moab_app_id, data%tag_name, data%tag_type, data%num_comp, data%tag_index(1))

  end subroutine setup_oneD_int_data_for_moab

  !------------------------------------------------------------------------------
  subroutine setup_twoD_real_data_for_moab(moab_app_id, tag_name, num_comp, num_cells_ghosted, data)
    !
    ! DESCRIPTION:
    ! Sets up 2D real-type data structure for performing MOAB-based halo exchange
    ! of data. This supports 'num_comp' values per grid cell to be exchanged.
    !
    use iso_c_binding
    use iMOAB, only : iMOAB_DefineTagStorage
    !
    ! ARGUMENT:
    integer                       , intent(in)  :: moab_app_id
    character(len=*)              , intent(in)  :: tag_name
    integer                       , intent(in)  :: num_comp
    integer                       , intent(in)  :: num_cells_ghosted
    type(twoD_real_data_for_moab) , intent(out) :: data
    !
    ! LOCAL VARIABLES:
    integer                                :: ierr

    data%moab_app_id    = moab_app_id
    data%tag_name       = trim(tag_name) // C_NULL_CHAR ! name
    data%tag_type       = 1                             ! 1 = dense, double
    data%num_tags       = 1                             ! a single tag
    data%entity_type(1) = 1                             ! element (== cell) based data
    data%num_comp       = num_comp                      ! number components in the tag
    data%ngcells        = num_cells_ghosted             ! number of grid cells

    ! allocate memory
    allocate(data%values(data%num_comp, data%ngcells))

    ! define the tag in MOAB
    ierr = iMOAB_DefineTagStorage(data%moab_app_id, data%tag_name, data%tag_type, data%num_comp, data%tag_index(1))

  end subroutine setup_twoD_real_data_for_moab

  !------------------------------------------------------------------------------
  subroutine domainlateral_init(domain_l)
    !
    ! DESCRIPTION:
    ! Creates data structure for doing halo exchanges using MOAB
    !
    implicit none
    !
    ! ARGUMENTS:
    type(domainlateral_type) :: domain_l ! domain datatype

    ! creates the MOAB tag 1D data at grid level
    call setup_oneD_int_data_for_moab(mlndghostid, 'grid_level_count', moab_gcell%num_ghosted, domain_l%grid_level_count)

  end subroutine domainlateral_init

  !------------------------------------------------------------------------------
  subroutine do_haloexchange_oneD_integer_data_for_moab(data)
    !
    ! DESCRIPTION:
    ! Perform MOAB-based halo exchange
    !
    use iMOAB, only : iMOAB_SetIntTagStorage, iMOAB_GetIntTagStorage, iMOAB_SynchronizeTags
    !
    ! ARGUMENT:
    type(oneD_int_data_for_moab) , intent(inout) :: data
    !
    ! LOCAL VARIABLE:
    integer :: ierr

    ! set the data in MOAB tag
    ierr = iMOAB_SetIntTagStorage(data%moab_app_id, data%tag_name, data%ngcells * data%num_comp, data%entity_type(1), data%values)
    if (ierr > 0) call endrun('Error: setting values in MOAB tag failed.')

    ! do the halo-exchange
    ierr = iMOAB_SynchronizeTags(data%moab_app_id, data%num_tags, data%tag_index(1), data%entity_type(1))
    if (ierr > 0) call endrun('Error: synchronization of MOAB tag failed.')

    ! get the data from MOAB tag
    ierr = iMOAB_GetIntTagStorage(data%moab_app_id, data%tag_name, data%ngcells * data%num_comp, data%entity_type(1), data%values)
    if (ierr > 0) call endrun('Error: getting values in MOAB tag failed.')

  end subroutine do_haloexchange_oneD_integer_data_for_moab

  !------------------------------------------------------------------------------
  subroutine do_haloexchange_twoD_real_data_for_moab(data)
    !
    ! DESCRIPTION:
    ! Perform MOAB-based halo exchange
    !
    use iMOAB, only : iMOAB_SetDoubleTagStorage, iMOAB_GetDoubleTagStorage, iMOAB_SynchronizeTags
    !
    ! INPUT ARGUMENT:
    type(twoD_real_data_for_moab) , intent(inout) :: data
    !
    ! LOCAL VARIABLE:
    integer :: ierr

    ! set the data in MOAB tag
    ierr = iMOAB_SetDoubleTagStorage(data%moab_app_id, data%tag_name, data%ngcells * data%num_comp, data%entity_type(1), data%values)
    if (ierr > 0) call endrun('Error: setting values in MOAB tag failed.')

    ! do the halo-exchange
    ierr = iMOAB_SynchronizeTags(data%moab_app_id, data%num_tags, data%tag_index(1), data%entity_type(1))
    if (ierr > 0) call endrun('Error: synchronization of MOAB tag failed.')

    ! get the data from MOAB tag
    ierr = iMOAB_GetDoubleTagStorage(data%moab_app_id, data%tag_name, data%ngcells * data%num_comp, data%entity_type(1), data%values)
    if (ierr > 0) call endrun('Error: getting values in MOAB tag failed.')

  end subroutine do_haloexchange_twoD_real_data_for_moab

  !------------------------------------------------------------------------------
  subroutine GridLevelIntegerDataHaloExchange(domain_l, begg, endg_owned, endg_all, elm_data)
    !
    ! DESCRIPTION:
    ! Performs halo exchange of integer data. It is assumed that there is only one value
    ! per grid cell. elm_data has data in ELM-format such that owned grid cells at the beginning
    ! followed by ghost grid cells. After MOAB-based halo exchange values are filled in
    ! elm_data corresponding to ghost cells.
    !
    implicit none
    !
    ! ARGUMENTS:
    type(domainlateral_type)         :: domain_l    ! domain datatype
    integer, intent(in)              :: begg        ! beginning index of grid cell
    integer, intent(in)              :: endg_owned  ! ending index for owned grid cells
    integer, intent(in)              :: endg_all    ! ending index for all (owned + ghost) grid cells
    integer, intent(inout) , pointer :: elm_data(:) ! data packed in ELM's format
    !
    ! LOCAL VARAIBLES:
    integer :: g, j, idx
    integer :: ierr

    ! convert data from ELM format to MOAB format
    do g = begg, endg_owned
       idx = moab_gcell%elm2moab(g)
       domain_l%grid_level_count%values(idx) = elm_data(g)
    end do

    ! perform halo exchange
    call do_haloexchange_oneD_integer_data_for_moab(domain_l%grid_level_count)

    ! convert data from MOAB format to ELM format
    do idx = 1, moab_gcell%num_ghosted
       if (.not.moab_gcell%is_owned(idx)) then
          g = moab_gcell%moab2elm(idx)
          elm_data(g) = domain_l%grid_level_count%values(idx)
       end if
    end do

  end subroutine GridLevelIntegerDataHaloExchange

  !------------------------------------------------------------------------------
  subroutine GridLevelRealDataHaloExchange(data_moab, begg, endg_owned, endg_all, elm_data)
    !
    ! DESCRIPTION:
    ! Performs halo exchange of real data. elm_data has data in ELM-format such that
    ! owned grid cells are at the beginning followed by ghost grid cells. After
    ! MOAB-based halo exchange, values in elm_data are filled for ghost cells.
    ! data_moab must be pre-set up by the caller via setup_twoD_real_data_for_moab.
    !
    implicit none
    !
    type(twoD_real_data_for_moab) , intent(inout) :: data_moab     ! caller-owned MOAB tag struct
    integer, intent(in)                           :: begg          ! beginning index of grid cell
    integer, intent(in)                           :: endg_owned    ! ending index for owned grid cells
    integer, intent(in)                           :: endg_all      ! ending index for all (owned + ghost) grid cells
    real(r8), intent(inout), pointer              :: elm_data(:,:) ! data packed in ELM's format
    !
    integer :: g, j, idx

    ! convert data from ELM format to MOAB format
    do g = begg, endg_owned
       idx = moab_gcell%elm2moab(g)
       do j = 1, data_moab%num_comp
          data_moab%values(j, idx) = elm_data(g, j)
       end do
    end do

    ! perform halo exchange
    call do_haloexchange_twoD_real_data_for_moab(data_moab)

    ! convert data from MOAB format to ELM format
    do idx = 1, moab_gcell%num_ghosted
       if (.not.moab_gcell%is_owned(idx)) then
          g = moab_gcell%moab2elm(idx)
          do j = 1, data_moab%num_comp
             elm_data(g, j) = data_moab%values(j, idx)
          end do
       end if
    end do

  end subroutine GridLevelRealDataHaloExchange

#else

  !-----------------------------------------------------------------------
  ! This is a stub for the case when PETSc is unavailable
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_abort
  use spmdMod     , only : masterproc
  use elm_varctl  , only : iulog
  use spmdMod     , only : masterproc, iam, npes, mpicom, comp_id
  use abortutils  , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
     
  type, public :: domainlateral_type
     integer :: dummy
  end type domainlateral_type

  type(domainlateral_type)    , public :: ldomain_lateral
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public domainlateral_init          ! allocates/nans domain types
  !
  !EOP
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: domainlateral_init
  !
  ! !INTERFACE:
  subroutine domainlateral_init(domain_l, cellsOnCell_old, edgesOnCell_old, &
       nEdgesOnCell_old, areaCell_old, dcEdge_old, dvEdge_old, &
       nCells_loc_old, nEdges_loc_old, maxEdges)
    !
    ! !ARGUMENTS:
    implicit none
    !
    !
    type(domainlateral_type) :: domain_l                     ! domain datatype
    integer , intent(in)     :: cellsOnCell_old(:,:)         ! grid cell level connectivity information
    integer , intent(in)     :: edgesOnCell_old(:,:)         ! index to determine distance between neighbors from dcEdge [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdgesOnCell_old(:)          ! number of edges                                           [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dcEdge_old(:)                ! distance between neighbors                                [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: dvEdge_old(:)                ! distance between vertices                                 [in natural order prior to domain decomposition]
    real(r8), intent(in)     :: areaCell_old(:)              ! area of grid cell                                         [in natural order prior to domain decomposition]
    integer , intent(in)     :: nCells_loc_old               ! number of local cell-to-cell connections                  [in natural order prior to domain decomposition]
    integer , intent(in)     :: nEdges_loc_old               ! number of edges                                           [in natural order prior to domain decomposition]
    integer , intent(in)     :: maxEdges                     ! max number of edges/neighbors

    character(len=*), parameter :: subname = 'domainlateral_init'

    call endrun(msg='ERROR ' // trim(subname) //': Requires '//&
         'PETSc, but the code was compiled without -DUSE_PETSC_LIB')

  end subroutine domainlateral_init


#endif

end module domainLateralMod
