module create_mapper_mod
  ! This module supports building mappers for use in unit tests

  use mct_mod
  use mct_wrapper_mod, only : mct_communicator, mct_compid
  use shr_kind_mod, only : r8 => shr_kind_r8
  use seq_map_type_mod, only : seq_map
  use simple_map_mod, only : simple_map_type

  implicit none
  private
  save

  public :: create_mapper ! create a simple mapper
  public :: clean_mapper  ! deallocate memory associated with a mapper

contains

  !-----------------------------------------------------------------------
  subroutine create_mapper(mapper, simple_map)
    !
    ! !DESCRIPTION:
    ! Create a simple mapper
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(seq_map), intent(out) :: mapper
    class(simple_map_type), intent(in) :: simple_map
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_mapper'
    !-----------------------------------------------------------------------

    mapper%copy_only = .false.
    mapper%rearrange_only = .false.
    mapper%esmf_map = .false.
    mapper%mapfile = ' '

    ! The strategy just relates to whether the mapping is done on the source or
    ! destination decomposition, which is irrelevant for a single-processor unit test
    mapper%strategy = 'X'

    ! May need to make this more sophisticated if it causes problems to use 0 for all
    ! mappers
    mapper%counter = 0

    allocate(mapper%gsmap_s)
    call create_gsmap(mapper%gsmap_s, simple_map%get_n_source_points())
    allocate(mapper%gsmap_d)
    call create_gsmap(mapper%gsmap_d, simple_map%get_n_dest_points())

    call mct_rearr_init(mapper%gsmap_s, mapper%gsmap_d, mct_communicator, mapper%rearr)

    call create_sMatp(mapper%sMatp, simple_map, mapper%gsmap_s, mapper%gsmap_d)

  end subroutine create_mapper


  !-----------------------------------------------------------------------
  subroutine clean_mapper(mapper)
    !
    ! !DESCRIPTION:
    ! Deallocate memory associated with a mapper.
    !
    ! This currently only deallocates the memory used in all mappers, NOT the
    ! cart3d-specific memory.
    !
    ! This assumes that gsmaps were created specially for this mapper, as is done in the
    ! convenience functions in this module (as opposed to having the mapper's gsmap
    ! pointers simply pointing to existing gsmaps).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(seq_map), intent(inout) :: mapper
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'clean_mapper'
    !-----------------------------------------------------------------------

    call mct_rearr_clean(mapper%rearr)
    call mct_sMatP_clean(mapper%sMatp)

    call mct_gsMap_clean(mapper%gsmap_s)
    deallocate(mapper%gsmap_s)
    call mct_gsMap_clean(mapper%gsmap_d)
    deallocate(mapper%gsmap_d)

  end subroutine clean_mapper

  !-----------------------------------------------------------------------
  subroutine create_gsmap(gsmap, npts)
    !
    ! !DESCRIPTION:
    ! Creates a simple, single-processor gsmap
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_gsMap), intent(out) :: gsmap
    integer, intent(in) :: npts
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'create_gsmap'
    !-----------------------------------------------------------------------

    call mct_gsMap_init(GSMap = gsmap, &
         comp_id = mct_compid, &
         ngseg = 1, &
         gsize = npts, &
         start = [1], &
         length = [npts], &
         pe_loc = [0])

  end subroutine create_gsmap

  !-----------------------------------------------------------------------
  subroutine create_sMatp(sMatp, simple_map, gsmap_s, gsmap_d)
    !
    ! !DESCRIPTION:
    ! Creates an sMatp object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(mct_sMatp), intent(out) :: sMatp
    class(simple_map_type), intent(in) :: simple_map
    type(mct_gsMap), intent(in) :: gsmap_s
    type(mct_gsMap), intent(in) :: gsmap_d
    !
    ! !LOCAL VARIABLES:
    integer :: n_elements   ! number of elements in the sparse matrix
    type(mct_sMat) :: sMati ! non-parallel sparse matrix

    ! The following pointers are needed because the MCT routines want pointer inputs
    integer, pointer :: source_indices(:)
    integer, pointer :: dest_indices(:)
    real(r8), pointer :: matrix_elements(:)

    character(len=*), parameter :: subname = 'create_sMatp'
    !-----------------------------------------------------------------------

    n_elements = simple_map%get_n_overlaps()

    call mct_sMat_init(sMati, &
         nrows = simple_map%get_n_dest_points(), &
         ncols = simple_map%get_n_source_points(), &
         lsize = n_elements)

    allocate(source_indices(n_elements))
    source_indices = simple_map%get_source_indices()
    call mct_sMat_ImpGColI(sMati, source_indices, n_elements)
    deallocate(source_indices)

    allocate(dest_indices(n_elements))
    dest_indices = simple_map%get_dest_indices()
    call mct_sMat_ImpGRowI(sMati, dest_indices, n_elements)
    deallocate(dest_indices)

    allocate(matrix_elements(n_elements))
    matrix_elements = simple_map%get_overlap_weights()
    call mct_sMat_ImpMatrix(sMati, matrix_elements, n_elements)
    deallocate(matrix_elements)

    call mct_sMatP_Init(sMatP, sMati, gsmap_s, gsmap_d, 0, mct_communicator, gsmap_s%comp_id)

    call mct_sMat_Clean(sMati)
  end subroutine create_sMatp


end module create_mapper_mod
