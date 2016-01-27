module cam_map_utils
  use pio,                 only: iMap=>PIO_OFFSET_KIND
  use cam_abortutils,      only: endrun
  use cam_logfile,         only: iulog
!!XXgoldyXX: v
use spmd_utils, only: npes, iam, mpicom, masterproc
use shr_sys_mod, only: shr_sys_flush
!!XXgoldyXX: ^

  implicit none
  private

  public iMap

!!XXgoldyXX: v
logical, public, save :: goldy_debug = .false.
!!XXgoldyXX: ^
  integer, private, save      :: unique_map_index = 0
  integer, private, parameter :: max_srcs         = 2
  integer, private, parameter :: max_dests        = 2

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_t: Information for a CAM map (between data array and
  !                 NetCDF file)
  !
  !    The map targets one or two dimensions of the NetCDF file.
  !    The 2-D map is useful for blocks of data (2 dimensions of array which
  !      map to one or two dimensions in the NetCDF file). For a 1-1 mapping,
  !      the first dimension is of size 3 (instead of 4).
  !      map(1, i) = index value of first src (e.g., lon, ncol)
  !      map(2, i) = index value of second src (e.g., lat, chunk)
  !      map(3, i) = global offset of first dest (e.g., lon, ncol)
  !      map(4, i) = global offset of second dest (e.g., lat, NA)
  !    src is the array dimension position corresponding to the map dimension
  !      in normal (not permuted) arrays.
  !      A negative pos denotes counting backwards from the last array dimension
  !    dest entries are the NetCDF dimension positions (must be increasing)
  !    It is an error to have src(1) == src(2).
  !
  !---------------------------------------------------------------------------
  type, public :: cam_filemap_t
    private
    integer                    :: index              =  0
    integer(iMap),   pointer   :: map(:,:)           => NULL()
    integer                    :: src(max_srcs)      =  0
    integer                    :: dest(max_dests)    =  0
    integer(iMap)              :: limits(max_srcs,2) =  -1
    integer                    :: nmapped            = -1
    integer,         pointer   :: blcksz(:)          => NULL() !e.g.,ncol(lchnk)
  contains
    procedure          :: init          => cam_filemap_init
    procedure          :: get_index     => cam_filemap_getIndex
    procedure          :: new_index     => cam_filemap_newIndex
    procedure          :: clear         => cam_filemap_clear
    procedure          :: copy          => cam_filemap_copy
    procedure          :: copy_elem     => cam_filemap_copyElem
    procedure          :: num_elem      => cam_filemap_size
    procedure          :: is_mapped     => cam_filemap_isMapped
    procedure          :: num_mapped    => cam_filemap_numMapped
    procedure          :: map_val       => cam_filemap_mapVal
    procedure          :: coord_vals    => cam_filemap_coordVals
    procedure          :: coord_dests   => cam_filemap_coordDests
    procedure          :: get_filemap   => cam_filemap_get_filemap
    procedure          :: has_blocksize => cam_filemap_has_blocksize
    procedure          :: blocksize     => cam_filemap_get_blocksize
    procedure          :: array_bounds  => cam_filemap_get_array_bounds
    procedure          :: active_cols   => cam_filemap_get_active_cols
    procedure          :: columnize     => cam_filemap_columnize
    procedure          :: compact       => cam_filemap_compact
!!XXgoldyXX: Cleanup when working
!    procedure          :: init_latlon   => cam_filemap_init_latlon
!    procedure          :: init_unstruct => cam_filemap_init_unstruct
!    generic, public    :: init          => init_latlon, init_unstruct
  end type cam_filemap_t

  !---------------------------------------------------------------------------
  !
  !  END: types BEGIN: private interfaces
  !
  !---------------------------------------------------------------------------

contains

!!#######################################################################
!!
!! index sorting routines:
!!  XXgoldyXX: Move to generic location?
!!
!!#######################################################################

  subroutine index_sort_vector(data, indices, compressval, dups_ok_in)
    use spmd_utils,      only: mpi_integer, mpi_integer8, iam, mpicom, npes
    use shr_mpi_mod,     only: shr_mpi_chkerr
    use cam_abortutils,  only: endrun
    use m_MergeSorts,    only: IndexSet, IndexSort

    ! Dummy arguments
    integer(iMap), pointer,  intent(in)    :: data(:)
    integer,       pointer,  intent(inout) :: indices(:)
    integer(iMap), optional, intent(in)    :: compressval
    logical,       optional, intent(in)    :: dups_ok_in

    ! Local variables
    integer                                :: num_elem   ! # mapped elements
    integer                                :: num_active ! # mapped pes
    integer                                :: ierr
    integer                                :: i
    integer                                :: lb, ub
    integer                                :: mycnt, my_first_elem
    integer                                :: mpi_group_world
    integer                                :: mpi_sort_group
    integer                                :: mpi_sort_comm
    integer                                :: my_sort_rank
    integer,       allocatable             :: sort_pes(:)
    integer,       allocatable             :: displs(:)
    integer,       allocatable             :: recvcounts(:)
    integer(iMap), allocatable             :: elements(:)
    integer(iMap), allocatable             :: local_elem(:)
    integer,       allocatable             :: ind(:)
    integer,       allocatable             :: temp(:)
    logical                                :: dups_ok ! .true. iff duplicates OK
    character(len=*),           parameter  :: subname = 'INDEX_SORT_VECTOR'

    ! Allow duplicate values?
    if (present(dups_ok_in)) then
      dups_ok = dups_ok_in
      if ((.not. dups_ok) .and. (.not. present(compressval))) then
        call endrun(trim(subname)//': dups_ok=.false. requires a compressval')
      end if
    else
      dups_ok = .true.
    end if

    ! The patch mapped values are in the number space of the master grid.
    ! They need to be compressed to go from 1 to the number of elements in the
    ! patch.
    ! Figure out the mapped elements in my portion of the patch mask
    if (.not. associated(data)) then
      mycnt = 0
      allocate(local_elem(mycnt))
      allocate(ind(mycnt))
      num_elem = 0
      lb = 0
      ub = -1
    else if (present(compressval)) then
      mycnt = COUNT(data /= compressval)
      allocate(local_elem(mycnt))
      allocate(ind(mycnt))
      num_elem = 0
      lb = LBOUND(data, 1)
      ub = UBOUND(data, 1)
      do i = lb, ub
        if (data(i) /= compressval) then
          num_elem = num_elem + 1
          local_elem(num_elem) = data(i)
          ind(num_elem) = i
        end if
      end do
    else
      lb = LBOUND(data, 1)
      ub = UBOUND(data, 1)
      mycnt = size(data)
      allocate(local_elem(mycnt))
      local_elem(1:mycnt) = data(lb:ub)
      num_elem = mycnt
    end if

    ! Find the tasks which have elements in this patch
    ! temp used for # elements per PE
    allocate(temp(0:npes-1))
    call MPI_allgather(mycnt, 1, MPI_integer, temp, 1, MPI_integer,           &
         mpicom, ierr)
    call shr_mpi_chkerr(ierr, subname//': MPI_allgather elements')
    num_active = COUNT(temp > 0)
    if (num_active > 1) then
      allocate(sort_pes(num_active))
      allocate(recvcounts(num_active))
      allocate(displs(num_active))
      num_elem = 0
      displs(1) = 0
      ! Find the number of mapped elements and number of pes in this patch
      my_sort_rank = -1
      do i = 0, npes - 1
        if (temp(i) > 0) then
          num_elem = num_elem + 1
          if (num_elem > num_active) then
            call endrun(subname//": overrun of sort_pes array")
          end if
          sort_pes(num_elem) = i
          if (iam == i) then
            my_sort_rank = num_elem - 1
            my_first_elem = displs(num_elem) + 1
          end if
          recvcounts(num_elem) = temp(i)
          if (num_elem < num_active) then
            displs(num_elem + 1) = displs(num_elem) + recvcounts(num_elem)
          end if
        end if
      end do
      if (num_elem < num_active) then
        call endrun(subname//": underrun of sort_pes array")
      end if
      if (my_sort_rank >= 0) then
        num_elem = SUM(temp) ! Total number of elements to sort
      else
        num_elem = 0
      end if
      deallocate(temp)  ! Cleanup
      ! Make a group with the active PEs
      call MPI_comm_group(mpicom, mpi_group_world, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_comm_group mpi_group_world')
      call MPI_group_incl(mpi_group_world, num_active, sort_pes,              &
           mpi_sort_group, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_group_incl sort_pes')
      ! Make a new communicator with the active PEs
      call MPI_comm_create(mpicom, mpi_sort_group, mpi_sort_comm, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_comm_create mpi_sort_comm')
      ! Collect all the elements for sorting (only active tasks now)
      allocate(elements(num_elem))
      if (mycnt > 0) then
        call MPI_allgatherv(local_elem, mycnt, MPI_integer8,                  &
             elements, recvcounts, displs, MPI_integer8, mpi_sort_comm, ierr)
        call shr_mpi_chkerr(ierr, subname//': MPI_allgatherv')
        ! Clean up for active PEs only
        call MPI_comm_free(mpi_sort_comm, ierr)
      end if
      ! General clean up
      call MPI_group_free(mpi_sort_group, ierr)
      deallocate(recvcounts)
      deallocate(displs)
    else if (mycnt > 0) then
      ! We are the only PE with patch info
      num_elem = mycnt
      allocate(elements(size(local_elem)))
      elements = local_elem
      my_first_elem = 1
    end if
    ! At this point, num_elem should always be the local # of items to sort
    if (num_elem > 0) then
      ! Sanity check
      if (size(elements) < num_elem) then
        call endrun(trim(subname)//": size(elements) must be >= num_elem")
      end if
      !! Do the sort, recvcounts will be the temporary index array
      allocate(recvcounts(size(elements)))
      allocate(displs(size(recvcounts)))
      call IndexSet(num_elem, recvcounts)
      call IndexSort(num_elem, recvcounts, elements, descend=.false.)
      ! Compress recvcounts (repeat data values)
      displs = 0
      do i = 1, num_elem - 1
        if (elements(recvcounts(i)) == elements(recvcounts(i + 1))) then
          displs(i + 1) = displs(i) + 1
        else
          displs(i + 1) = displs(i)
        end if
      end do
      ! Unload recvcounts into indices. Assume indices is initialized w/ default
      do i = 1, num_elem
        if ( (recvcounts(i) >= my_first_elem) .and.                           &
             (recvcounts(i) < (my_first_elem + mycnt))) then
          if (.not. dups_ok) then
            ! Use our indirect access to set the correct indices location
            ! ind array already has lb offset included
            ! Eliminate duplicate values
            if ((i > 1) .and. (displs(i) > displs(MAX((i - 1),1)))) then
              indices(ind(recvcounts(i) - my_first_elem + 1)) = compressval
            else
              indices(ind(recvcounts(i) - my_first_elem + 1)) = i - displs(i)
            end if
          else if (allocated(ind)) then
            indices(ind(recvcounts(i) - my_first_elem + 1)) = i - displs(i)
          else
            ! recvcounts points directly at a local location
            ! NB: repeat data values all get same index
            indices(recvcounts(i) - my_first_elem + lb) = i - displs(i)
          end if
        end if
      end do
      deallocate(recvcounts)
      deallocate(displs)
    end if

    if (allocated(ind)) then
      deallocate(ind)
    end if
    if (allocated(elements)) then
      deallocate(elements)
    end if
    deallocate(local_elem)

  end subroutine index_sort_vector

!!#######################################################################
!!
!! CAM grid mapping functions
!!
!!#######################################################################

  integer function cam_filemap_getIndex(this)
    ! Dummy variable
    class(cam_filemap_t)                 :: this

    cam_filemap_getIndex = this%index

  end function cam_filemap_getIndex

  subroutine cam_filemap_newIndex(this)
    ! Dummy variable
    class(cam_filemap_t)                 :: this

    unique_map_index = unique_map_index + 1
    this%index = unique_map_index

  end subroutine cam_filemap_newIndex

  subroutine cam_filemap_init(this, pemap, unstruct, src, dest)
    ! Dummy arguments
    class(cam_filemap_t)                 :: this
    integer(iMap),    pointer            :: pemap(:,:) ! Map elem for this PE
    logical,                  intent(in) :: unstruct
    integer,                  intent(in) :: src(:)
    integer, optional,        intent(in) :: dest(:)

    ! Local variables
    integer                              :: i ! Loop index
    integer                              :: index

    ! This shouldn't happen but maybe we will decide to reuse these
    if (associated(this%map)) then
      deallocate(this%map)
      nullify(this%map)
    end if

    ! Check in case these ever change (because algorithm then must be modified)
    if ((max_srcs /= 2) .or. (max_dests /= 2)) then
      call endrun('cam_filemap_init: max_src or max_dest modified')
    end if

    ! Some items are simply copied
    if (associated(pemap)) then
      this%map => pemap
    else
      nullify(this%map)
    end if
    this%src =  src
    if (present(dest)) then
      ! Structred grids will likely always have dest = (1, 2) but maybe . . .
      this%dest = dest
    else if (unstruct) then
      ! Unstructured grids are (so far) always spread along the first dimension
      this%dest(1) = 1
      this%dest(2) = 0
    else
      this%dest(1) = 1
      this%dest(2) = 2
    end if
    ! We may have holes in the 'block' decomposition which is specified by
    ! having src(2) < 0. 
    ! NB: This is currently a special purpose hack in that it is purely
    !     convention that the last dimension specifies the block index and
    !     that those blocks may not be filled.
    !     The proper way to generalize this functionality is to allow
    !     src(1) to also be < 0 and to look for holes there as well
    if (associated(this%map)) then
      do i = 1, max_srcs
        if (ANY(this%map(i,:) > 0)) then
          ! Min of all src(i) values
          this%limits(i, 1) = MINVAL(this%map(i,:), mask=(this%map(i,:)>0))
          ! Max of all src(i) values
          this%limits(i, 2) = MAXVAL(this%map(i,:), mask=(this%map(i,:)>0))
        else
          this%limits(i,1) = 0
          this%limits(i,2) = -1
        end if
      end do

      this%nmapped = 0
      do index = 1, this%num_elem()
        if ((this%dest(1) > 0) .and. (this%map(max_srcs+1, index) > 0)) then
          if (this%dest(2) > 0) then
            ! Can't do this test unless we know the dim 2 is large enough
            if (this%map(max_srcs+2, index) > 0) then
              this%nmapped = this%nmapped + 1
            end if
          else
            this%nmapped = this%nmapped + 1
          end if
        end if
      end do
    else
      this%limits(:,1) = 0
      this%limits(:,2) = -1
      this%nmapped = 0
    end if
    if (src(max_srcs) < 0) then
      ! This shouldn't happen but maybe we will decide to reuse these
      if (associated(this%blcksz)) then
        deallocate(this%blcksz)
        nullify(this%blcksz)
      end if
      allocate(this%blcksz(this%limits(max_srcs,1):this%limits(max_srcs,2)))
      this%blcksz = 0
      do i = 1, this%num_elem()
        index = this%map(max_srcs, i)
        if (this%is_mapped(i)) then
          this%blcksz(index) = this%blcksz(index) + 1
        end if
      end do
    end if

    call this%new_index()

  end subroutine cam_filemap_init

  subroutine cam_filemap_clear(this)
    ! Dummy arguments
    class(cam_filemap_t)               :: this

    if (associated(this%map)) then
      this%map = 0
    end if
    this%limits(:,1) = 0
    this%limits(:,2) = -1
    this%nmapped = 0

    ! Update the index in case a decomp was made with the old values
    call this%new_index()

  end subroutine cam_filemap_clear

  subroutine cam_filemap_copy(this, map)
    ! Dummy arguments
    class(cam_filemap_t)               :: this
    type(cam_filemap_t),  intent(in)   :: map

    ! This shouldn't happen but maybe we will decide to reuse these
    if (associated(this%map)) then
      deallocate(this%map)
      nullify(this%map)
    end if

    if (associated(map%map)) then
      allocate(this%map(size(map%map, 1), size(map%map, 2)))
      if (map%num_elem() > 0) then
        this%map = map%map
      end if
    else
      nullify(this%map)
    end if
    this%src = map%src
    this%dest = map%dest
    this%limits = map%limits
    this%nmapped = map%nmapped

    ! This shouldn't happen but maybe we will decide to reuse these
    if (associated(this%blcksz)) then
      deallocate(this%blcksz)
      nullify(this%blcksz)
    end if
    if (associated(map%blcksz)) then
      allocate(this%blcksz(size(map%blcksz)))
      this%blcksz = map%blcksz
    end if

    ! Even a copy has to have a unique index
    call this%new_index()

  end subroutine cam_filemap_copy

  subroutine cam_filemap_copyElem(this, map, index)
    ! Dummy arguments
    class(cam_filemap_t)               :: this
    type(cam_filemap_t), intent(in)    :: map
    integer,             intent(in)    :: index

    if (this%is_mapped(index)) then
      this%nmapped = this%nmapped - 1
    end if
    this%map(:, index) = map%map(:, index)
    if (this%is_mapped(index)) then
      this%nmapped = this%nmapped + 1
    end if

  end subroutine cam_filemap_copyElem

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_size: Total number of elements in the map
  !
  !---------------------------------------------------------------------------
  integer function cam_filemap_size(this)
    ! Dummy variable
    class(cam_filemap_t)                 :: this

    if (associated(this%map)) then
      cam_filemap_size = size(this%map,2)
    else
      cam_filemap_size = 0
    end if
  end function cam_filemap_size

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_isMapped: Return .true. iff the index value is mapped
  !
  !---------------------------------------------------------------------------
  elemental logical function cam_filemap_isMapped(this, index)
    ! Dummy arguments
    class(cam_filemap_t),      intent(in)    :: this
    integer,                   intent(in)    :: index

    if (associated(this%map)) then
      cam_filemap_isMapped = (this%map(3, index) > 0)
      if ((size(this%map, 1) > 3) .and. cam_filemap_isMapped) then
        cam_filemap_isMapped = (this%map(4, index) > 0)
        ! No else needed
      end if
    else
      cam_filemap_isMapped = .false.
    end if

  end function cam_filemap_isMapped

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_numMapped: Number of elements in the map with non-zero entries
  !
  !---------------------------------------------------------------------------
  integer function cam_filemap_numMapped(this)
    ! Dummy variable
    class(cam_filemap_t)                 :: this

    if (associated(this%map)) then
      cam_filemap_numMapped = this%nmapped
    else
      cam_filemap_numMapped = 0
    end if
  end function cam_filemap_numMapped

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_mapVal: Calculate an offset value from a map
  !
  !---------------------------------------------------------------------------
  integer(iMap) function cam_filemap_mapVal(this, index, dsize, dest_in)
    ! Dummy arguments
    class(cam_filemap_t)                     :: this
    integer,                   intent(in)    :: index
    integer(iMap),             intent(in)    :: dsize(:)
    integer,       optional,   intent(in)    :: dest_in(:)

    ! Local variables
    integer                                  :: d(max_dests)

    if (associated(this%map)) then
      if (present(dest_in)) then
        d = dest_in
      else
        d = this%dest
      end if
      if (this%map(3, index) > 0) then
        cam_filemap_mapVal = ((this%map(3, index) - 1) * dsize(d(1))) + 1
        if (size(this%map, 1) > 3) then
          if (this%map(4, index) > 0) then
            cam_filemap_mapVal = ((this%map(4, index) - 1) * dsize(d(2))) +   &
                 cam_filemap_mapVal  ! No +1 because it is offset from map(3,:)
          else
            cam_filemap_mapVal = 0
          end if
        ! No else needed
        end if
      else
        cam_filemap_mapVal = 0
      end if
    else
      cam_filemap_mapVal = 0
    end if
  end function cam_filemap_mapVal

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_coordVals: Find coord indices matching map index
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_coordVals(this, index, lonIndex, latIndex, isMapped)
    ! Dummy arguments
    class(cam_filemap_t)                     :: this
    integer,                   intent(in)    :: index
    integer,                   intent(out)   :: lonIndex
    integer,                   intent(out)   :: latIndex
    logical, optional,         intent(out)   :: isMapped

    if (associated(this%map)) then
      if (size(this%map,1) > (max_srcs + 1)) then
        lonIndex = this%map(1, index)
        latIndex = this%map(2, index)
      else
        lonIndex = index
        latIndex = index
      end if
      if (present(isMapped)) then
        isMapped = this%is_mapped(index)
      end if
    else if (present(isMapped)) then
      lonIndex = 0
      latIndex = 0
      isMapped = .false.
    else
      call endrun("cam_filemap_coordVals: must have map or pass isMapped")
    end if
  end subroutine  cam_filemap_coordVals

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_coordDests: Find coord indices matching map index
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_coordDests(this, index, lonIndex, latIndex, isMapped)
    ! Dummy arguments
    class(cam_filemap_t)                     :: this
    integer,                   intent(in)    :: index
    integer,                   intent(out)   :: lonIndex
    integer,                   intent(out)   :: latIndex
    logical, optional,         intent(out)   :: isMapped

    if (associated(this%map)) then
      lonIndex = this%map(3, index)
      if (size(this%map,1) > (max_srcs + 1)) then
        latIndex = this%map(4, index)
      else
        latIndex = 0
      end if
      if (present(isMapped)) then
        isMapped = this%is_mapped(index)
      end if
    else if (present(isMapped)) then
      lonIndex = 0
      latIndex = 0
      isMapped = .false.
    else
      call endrun("cam_filemap_coordDests: must have map or pass isMapped")
    end if
  end subroutine  cam_filemap_coordDests

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_get_filemap: Create the mapping between an array and a file
  !             This is the workhorse function for creating a PIO decomp DOF
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_get_filemap(this, fieldlens, filelens, filemap,      &
       src_in, dest_in, permutation_in)

    ! Dummy arguments
    class(cam_filemap_t)                      :: this
    integer,                    intent(in)    :: fieldlens(:)
    integer,                    intent(in)    :: filelens(:)
    integer(iMap), pointer                    :: filemap(:)
    integer,          optional, intent(in)    :: src_in(:)
    integer,          optional, intent(in)    :: dest_in(:)
    integer,          optional, intent(in)    :: permutation_in(:)

    ! Local variables
    integer                       :: srclens(7)          ! field dim lens
    integer                       :: srccnt              ! Rank of fieldlens
    integer(iMap), allocatable    :: dsize(:)
    integer                       :: mapind(max_srcs)   ! Source index for map
    integer,       allocatable    :: src_ind(:)     ! Map src to file dims
    integer                       :: fmind, j
    integer(iMap)                 :: mapSize, mapPos, pos, fileSize
    integer                       :: mapcnt         ! Dimension count
    integer                       :: locsize        ! Total local # elements
    integer                       :: tind, tlen     ! Temporarys
    integer                       :: i1, i2, i3, i4, i5, i6, i7
    integer                       :: i(7)

    ! This shouldn't happen but, who knows what evil lurks in the hearts of SEs
    if (associated(filemap)) then
      deallocate(filemap)
      nullify(filemap)
    end if

    !
    fileSize = product(filelens)
    srccnt = size(fieldlens)
    srclens(1:srccnt) = fieldlens(1:srccnt)
    if (srccnt < 7) then
      srclens(srccnt+1:7) = 1
    end if

    ! Allocate output filemap (dof)
    allocate(filemap(product(fieldlens)))
    filemap = 0

    ! Find map source dimensions in input array
    mapPos = 1                                ! To compare against map size
    mapind = 0
    if (present(src_in)) then
      mapcnt = size(src_in)  ! Just used until end of loop below
      if (mapcnt > max_srcs) then
        call endrun('cam_filemap_get_filemap: src_in too large')
      end if
    end if
    do j = 1, max_srcs
      if (present(src_in)) then
        if (mapcnt >= j) then
          if (src_in(j) < 0) then
            mapind(j) = srccnt + src_in(j) + 1
          else
            mapind(j) = src_in(j)
          end if
        end if
      else
        ! A src == 0 means that map dimension is not used
        if (this%src(j) /= 0) then
          ! src < 0 is count from last array dim
          if (this%src(j) < 0) then
            mapind(j) = srccnt + this%src(j) + 1
          else
            mapind(j) = this%src(j)
          end if
        end if
      end if
      if (mapind(j) > 0) then
        mapPos = mapPos * srclens(mapind(j))  ! To compare against map size
      end if
    end do
    mapcnt = COUNT(mapind /= 0)

    ! Check that map matches dims
    ! Since patch maps are compressed, we can't do an equal compare but
    !   it is still an error if the map has more elements than the array
    mapSize = this%num_elem()
    if (mapPos < this%num_mapped()) then
      call endrun('cam_filemap_get_filemap: Map size too large for array dims')
    end if

    ! dsize is a global offset for each dimension
    allocate(dsize(size(filelens)))
    dsize(1) = 1
    do j = 1, size(filelens) - 1
      dsize(j + 1) = dsize(j) * filelens(j)
    end do

!!XXgoldyXX:
!goldy_debug = (iam == 1)
if (goldy_debug) then
   write(iulog, *) 'cam_filemap_get_filemap: mapind:',mapind,', srccnt =',srccnt
   call shr_sys_flush(iulog)
   write(iulog, *) 'cam_filemap_get_filemap: srclens:',srclens
   call shr_sys_flush(iulog)
   write(iulog, *) 'cam_filemap_get_filemap: filelens:',filelens
   call shr_sys_flush(iulog)
   write(iulog, *) 'cam_filemap_get_filemap: dsize:',dsize
   call shr_sys_flush(iulog)
   if (present(dest_in)) then
     write(iulog, *) 'cam_filemap_get_filemap: dest_in:',dest_in
   else
     write(iulog, *) 'cam_filemap_get_filemap: dest:',this%dest
   end if
   call shr_sys_flush(iulog)
   if (present(src_in)) then
     write(iulog, *) 'cam_filemap_get_filemap: src_in:',src_in
   else
     write(iulog, *) 'cam_filemap_get_filemap: src:',this%src
   end if
   call shr_sys_flush(iulog)
end if
!!XXgoldyXX:
    ! src_ind maps each source dimension to the corresponding file dimension
    allocate(src_ind(srccnt))
    if (present(permutation_in)) then
      if (size(permutation_in) /= size(src_ind)) then
        call endrun('cam_filemap_get_filemap: permutation_in must have same rank as fieldlens')
      end if
      src_ind = permutation_in
    else
      src_ind = 0
      fmind = 1      ! Here fmind is first available permutation slot
      do j = 1, srccnt
        ! We need to find the offset location for each non-mapped src dimension
        if (ANY(mapind == j)) then
          continue
        else
          ! Find the next available output dimension
          if (present(dest_in)) then
            do while (ANY(dest_in == fmind))
              fmind = fmind + 1
              if (fmind > size(dsize)) then
                call endrun('cam_filemap_get_filemap: permutation calculation dest_in error')
              end if
            end do
          else
            do while (ANY(this%dest == fmind))
              fmind = fmind + 1
              if (fmind > size(dsize)) then
                call endrun('cam_filemap_get_filemap: permutation calculation dest error')
              end if
            end do
          end if
          if (fmind > size(dsize)) then
            call endrun('cam_filemap_get_filemap: permutation calculation error')
          end if
          src_ind(j) = fmind
          fmind = fmind + 1
        end if
      end do
    end if
       
    ! Step through the map and fill in local positions for each entry
    fmind = 1
    do i7 = 1, srclens(7)
      i(7) = i7
      do i6 = 1, srclens(6)
        i(6) = i6
        do i5 = 1, srclens(5)
          i(5) = i5
          do i4 = 1, srclens(4)
            i(4) = i4
            do i3 = 1, srclens(3)
              i(3) = i3
              do i2 = 1, srclens(2)
                i(2) = i2
                do i1 = 1, srclens(1)
                  i(1) = i1
                  pos = 0  ! Offset from map pos so starts at zero
                  tind = 1 ! Offset into the map
                  tlen = 1
                  do j = 1, srccnt
                    if (ANY(mapind == j)) then
                      ! j is a distributed field dimension index
                      tind = tind + ((i(j) - 1) * tlen)
                      tlen = tlen * srclens(j)
                    else
                      ! j is a local field dimension index
                      pos = pos + ((i(j) - 1) * dsize(src_ind(j)))
                    end if
                  end do
                  if (tind > mapSize) then
                    call endrun('cam_filemap_get_filemap: internal error, tind')
                  end if
                  mapPos = this%map_val(tind, dsize, dest_in)
                  if ((mapPos > 0) .and. ((pos + mapPos) > fileSize)) then
                    call endrun('cam_filemap_get_filemap: internal error, pos')
                  end if
                  if ((pos + mapPos) < 0) then
                    call endrun('cam_filemap_get_filemap: internal error, mpos')
                  end if
                  if (mapPos > 0) then
                    filemap(fmind) = pos + mapPos
                  else
                    ! This element is not mapped
                    filemap(fmind) = 0
                  end if
                  fmind = fmind + 1
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    if ((fmind - 1) /= size(filemap)) then
      call endrun('cam_filemap_get_filemap: internal error, fmind')
    end if
    deallocate(dsize)
  end subroutine cam_filemap_get_filemap

  integer function cam_filemap_get_blocksize(this, block_id)

    ! Dummy arguments
    class(cam_filemap_t)                       :: this
    integer,                    intent(in)     :: block_id

    if (.not. this%has_blocksize()) then
      call endrun('cam_filemap_get_blocksize: filemap has no blocks')
    else if ((block_id < LBOUND(this%blcksz, 1)) .or.                         &
         (block_id > UBOUND(this%blcksz, 1))) then
      call endrun('cam_filemap_get_blocksize: block_id out of range')
    else
      cam_filemap_get_blocksize = this%blcksz(block_id)
    end if
  end function cam_filemap_get_blocksize

  logical function cam_filemap_has_blocksize(this, lbnd, ubnd)

    ! Dummy arguments
    class(cam_filemap_t)                       :: this
    integer,             optional, intent(out) :: lbnd
    integer,             optional, intent(out) :: ubnd

    cam_filemap_has_blocksize = associated(this%blcksz)
    if (present(lbnd)) then
      lbnd = LBOUND(this%blcksz, 1)
    end if
    if (present(ubnd)) then
      ubnd = UBOUND(this%blcksz, 1)
    end if

  end function cam_filemap_has_blocksize

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_get_array_bounds: Sets grid bounds for the relevant array
  !            Only modifies the dimensions corresponding to the map's src
  !            dims should be sized (rank,2) with the second dimension used
  !            to store lower(1) and upper(2) bounds
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_get_array_bounds(this, dims)

    ! Dummy arguments
    class(cam_filemap_t)                       :: this
    integer,                    intent(inout)  :: dims(:,:)

    ! Local variables
    integer                                    :: rank ! rank of target array
    integer                                    :: i    ! Loop variable

    rank = size(dims,1)
    if (size(dims,2) < 2) then
      call endrun('cam_filemap_get_array_bounds: second dim of dims must be 2')
    end if

    if (MAXVAL(this%limits(1:max_srcs, 1:2)) > HUGE(kind(dims))) then
      call endrun('cam_filemap_get_array_bounds: limits too large')
    end if
    do i = 1, max_srcs
      if (this%src(i) > 0) then
        if (this%src(i) > rank) then
          call endrun('cam_filemap_get_array_bounds: rank too small')
        else
!!XXgoldyXX: Maybe modify definition of this%limits?
!          dims(i, 1:2) = INT(this%limits(i, 1:2), kind=kind(dims))
          if (associated(this%map)) then
            if (size(this%map) > 0) then
              dims(i, 1) = MINVAL(this%map(i,:))
              dims(i, 2) = MAXVAL(this%map(i,:))
            else
              dims(i, 1) = 0
              dims(i, 2) = -1
            end if
          else
            dims(i, 1) = 0
            dims(i, 2) = -1
          end if
        end if
      else if (this%src(i) < 0) then
!!XXgoldyXX: Maybe modify definition of this%limits?
!        dims(rank + this%src(i) + 1, 1:2) = INT(this%limits(i, 1:2), kind=kind(dims))
        if (associated(this%map)) then
          if (size(this%map) > 0) then
            dims(rank + this%src(i) + 1, 1) = MINVAL(this%map(i,:))
            dims(rank + this%src(i) + 1, 2) = MAXVAL(this%map(i,:))
          else
            dims(rank + this%src(i) + 1, 1) = 0
            dims(rank + this%src(i) + 1, 2) = -1
          end if
        else
          dims(rank + this%src(i) + 1, 1) = 0
          dims(rank + this%src(i) + 1, 2) = -1
        end if
       ! No else (zero means unused position)
      end if
    end do
  end subroutine cam_filemap_get_array_bounds

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_get_active_cols: Find which columns are active in a dimension
  !        Because we normally decompose columns (blocks or chunks) in the
  !        last dimension, we default to that dimension
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_get_active_cols(this, colnum, active, srcdim_in)

    ! Dummy arguments
    class(cam_filemap_t)                       :: this
    integer,                    intent(in)     :: colnum
    logical,                    intent(out)    :: active(:)
    integer, optional,          intent(in)     :: srcdim_in

    ! Local variables
    integer                                    :: srcdim

    if (present(srcdim_in)) then
      srcdim = srcdim_in
    else
      srcdim = max_srcs
    end if

    ! Sanity checks
    if ((srcdim < 1) .or. (srcdim > max_srcs)) then
      call endrun('cam_filemap_get_active_cols: srcdim out of range')
    else if (this%src(srcdim) >= 0) then
      call endrun('cam_filemap_get_active_cols: Invalid srcdim')
    else if (size(active) < size(this%map, (3 - srcdim))) then
      call endrun('cam_filemap_get_active_cols: active too small')
    else if (colnum < LBOUND(this%map, srcdim)) then
      call endrun('cam_filemap_get_active_cols: colnum too small')
    else if (colnum > UBOUND(this%map, srcdim)) then
      call endrun('cam_filemap_get_active_cols: colnum too large')
    ! else OK
    end if

    active = .false.
!!XXgoldyXX: This is probably completely wrong. What we want is column info
    select case(srcdim)
      case (1)
        active(1:size(this%map, 2)) = (this%map(colnum,:) > 0)
      case (2)
        active(1:size(this%map, 1)) = (this%map(:,colnum) > 0)
      case default
        call endrun('cam_filemap_get_active_cols: Invalid srcdim?!?')
      end select
  end subroutine cam_filemap_get_active_cols

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_columnize: Convert lon/lat map to ncol
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_columnize(this)
    use spmd_utils,  only: mpi_sum, mpi_integer, mpicom
    use shr_mpi_mod, only: shr_mpi_chkerr

    ! Dummy argument
    class(cam_filemap_t)                      :: this

    ! Local variables
    integer                                   :: i, j
    integer                                   :: lmax(2)
    integer                                   :: maxind(2)
    integer                                   :: offset(2)
    integer                                   :: ierr
    integer(iMap), pointer                    :: newmap(:,:) => NULL()
    character(len=*),              parameter  :: subname = 'CAM_FILEMAP_COLUMNIZE'
    ! Create a new map with same size and ordering
    ! We need the max lon/lat for calculating global offsets
    lmax = 0
    if (associated(this%map)) then
      if (size(this%map, 1) == max_srcs) then
        call endrun(trim(subname)//': must have at least 1 destination coord')
      else if (size(this%map, 1) > (max_srcs + 2)) then
        call endrun(trim(subname)//': has more than 2 destination coords')
      end if
      if (size(this%map, 2) > 0) then
        lmax(1) = MAXVAL(this%map(max_srcs + 1, :))
        if (size(this%map, 1) == (max_srcs + 2)) then
          lmax(2) = MAXVAL(this%map(max_srcs + 2, :))
        end if
      end if
    end if
    call MPI_allreduce(lmax(1), maxind(1), 1, MPI_integer, mpi_sum, mpicom, ierr)
    call shr_mpi_chkerr(ierr, subname//': MPI_allreduce maxlon')
    call MPI_allreduce(lmax(2), maxind(2), 1, MPI_integer, mpi_sum, mpicom, ierr)
    call shr_mpi_chkerr(ierr, subname//': MPI_allreduce maxlat')
    if (associated(this%map)) then
      if (size(this%map, 1) == (max_srcs + 2))then
        ! Create the new map
        allocate(newmap(max_srcs + 1, size(this%map, 2)))
        ! Who's on first?
        if (ANY(this%dest(1:2) <= 0)) then
          call endrun(trim(subname)//': can only handle positive dest indices')
        else if (this%dest(1) < this%dest(2)) then
          offset(1) = 1
          offset(2) = maxind(1)
        else if (this%dest(1) > this%dest(2)) then
          offset(1) = maxind(2)
          offset(2) = 1
        else
          call endrun(trim(subname)//': dest indices cannot be equal')
        end if
        do i = 1, size(newmap, 2)
          newmap(1:max_srcs, i) = this%map(1:max_srcs, i)
          j = ((this%map(max_srcs+1, i) * offset(1)) +                        &
               (this%map(max_srcs+2, i) * offset(2)))
          newmap(max_srcs+1, i) = j
        end do
        ! Replace our map with the new one
        deallocate(this%map)
        this%map => newmap
        nullify(newmap)
        ! Fixup dest
        this%dest(1) = 1
        this%dest(2:) = 0
        ! no else (do nothing if already ncol)
      end if ! End if lon/lat map
    end if

  end subroutine cam_filemap_columnize

  !---------------------------------------------------------------------------
  !
  !  cam_filemap_compact: Pack all active elements into a 1-based file order
  !                       Also pack latitude and longitude maps
  !
  !---------------------------------------------------------------------------
  subroutine cam_filemap_compact(this, lonmap, latmap,                        &
       num_lons, num_lats, num_mapped, columnize, dups_ok_in)
    use spmd_utils,  only: mpi_sum, mpi_integer, mpicom
    use shr_mpi_mod, only: shr_mpi_chkerr

    ! Dummy arguments
    class(cam_filemap_t)                     :: this
    integer(iMap),  pointer                  :: lonmap(:)
    integer(iMap),  pointer                  :: latmap(:)
    integer,        optional, intent(out)    :: num_lons
    integer,        optional, intent(out)    :: num_lats
    integer,        optional, intent(out)    :: num_mapped
    logical,        optional, intent(in)     :: columnize  ! Convert to ncol
    logical,        optional, intent(in)     :: dups_ok_in ! Dup coords OK

    ! Local variables
    integer                                  :: i, j
    integer                                  :: ierr
    integer(iMap), pointer                   :: data(:) => NULL()
    integer,       pointer                   :: indices(:) => NULL()
    integer                                  :: tmp_size
    logical                                  :: dok
    character(len=*),              parameter :: subname = 'CAM_FILEMAP_COMPACT'

    !! Possibly convert lon/lat map to ncol
    if (present(columnize)) then
      if (columnize) then
        call this%columnize()
      end if
    end if
    !! Are duplicate coordinate indices (lat/lon) OK?
    if (present(dups_ok_in)) then
      dok = dups_ok_in
    else
      dok = .false.
    end if
    ! Get a global index sort of mapped elements.
    do i = 1, max_dests
      if (this%dest(i) > 0) then
        if (associated(this%map)) then
          if (size(this%map, 1) >= max_srcs + i) then
            data => this%map(max_srcs + i, :)
          else
            nullify(data)
          end if
        else
          nullify(data)
        end if
        ! Allocate indices if necessary
        if (associated(indices)) then
          deallocate(indices)
          nullify(indices)
        end if
        if (associated(data)) then
          allocate(indices(LBOUND(data, 1):UBOUND(data, 1)))
          indices = 0
        else
          allocate(indices(0))
        end if
      end if
      call index_sort_vector(data, indices, compressval=0_iMap)
      if (associated(data) .and. associated(indices)) then
        data = indices
      end if
    end do
    if (associated(indices)) then
      deallocate(indices)
      nullify(indices)
    end if
    ! Get a global index sort of lat and lon maps 
    !! Compress latmap
    if (associated(latmap)) then
      ! Allocate indices
      allocate(indices(LBOUND(latmap, 1):UBOUND(latmap, 1)))
      indices = 0
    end if
    call index_sort_vector(latmap, indices, compressval=0_iMap, dups_ok_in=dok)
    if (associated(latmap)) then
      latmap = indices
      deallocate(indices)
      nullify(indices)
    end if
    !! Compress lonmap
    if (associated(lonmap)) then
      allocate(indices(LBOUND(lonmap, 1):UBOUND(lonmap, 1)))
      indices = 0
    end if
    call index_sort_vector(lonmap, indices, compressval=0_iMap, dups_ok_in=dok)
    if (associated(lonmap)) then
      lonmap = indices
      deallocate(indices)
      nullify(indices)
    end if

    if (present(num_mapped)) then
      if (this%num_elem() > 0) then
        ! Total number of mapped elements
        tmp_size=0
        do i = 1, this%num_elem()
          if (this%is_mapped(i)) then
            tmp_size = tmp_size + 1
          end if
        end do
      else
        tmp_size = 0
      end if
      call MPI_allreduce(tmp_size, num_mapped, 1, MPI_integer,              &
           mpi_sum, mpicom, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_allreduce num_mapped')
      if (num_mapped <= 0) then
        call endrun(trim(subname)//': num_mapped <= 0')
      end if
    end if
    if (present(num_lons)) then
      if (associated(lonmap)) then
        tmp_size = COUNT(lonmap /= 0)
      else
        tmp_size = 0
      end if
      call MPI_allreduce(tmp_size, num_lons, 1, MPI_integer,                  &
           mpi_sum, mpicom, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_allreduce num_lons')
      if (num_lons <= 0) then
        call endrun(trim(subname)//': numlons <= 0')
      end if
    end if
    if (present(num_lats)) then
      if (associated(latmap)) then
        tmp_size = COUNT(latmap /= 0)
      else
        tmp_size = 0
      end if
      call MPI_allreduce(tmp_size, num_lats, 1, MPI_integer,                  &
           mpi_sum, mpicom, ierr)
      call shr_mpi_chkerr(ierr, subname//': MPI_allreduce num_lats')
      if (num_lats <= 0) then
        call endrun(trim(subname)//': numlats <= 0')
      end if
    end if

  end subroutine cam_filemap_compact

end module cam_map_utils
