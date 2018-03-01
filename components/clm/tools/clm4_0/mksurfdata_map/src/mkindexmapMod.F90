module mkindexmapMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkindexmapMod
!
! !DESCRIPTION:
! Module containing subroutines for making maps of index data.
!
! This includes a routine for making a map using the dominant type among the input grid
! cells making up a given output cell, as well as routines for using an index map as
! indices into a lookup table, to essentially paint-by-number some other field, and some
! other related routines
!
! WJS (2-1-12): There is a lookup_2d subroutine, but not a lookup_1d (or any other
! dimensionality). That is simply because I needed lookup_2d, but have not yet needed a
! routine of other dimensionalities. In the future, it would probably be helpful to at
! least have lookup_1d and lookup_1d_netcdf. If this is done, see my notes under the
! lookup_2d_netcdf routine for some thoughts on avoiding duplication.
!
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8
   use mkncdio, only : nf_max_name

   implicit none
   private

! !PUBLIC TYPES:
!
   ! dim_slice_type: stores information about dimensions that we use for slicing a multi-
   ! dimensional variable
   type dim_slice_type
      character(len=nf_max_name) :: name  ! name of this dimension
      integer                    :: val   ! index to use for the slice
   end type dim_slice_type
   public :: dim_slice_type
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: get_dominant_indices  ! make output map based on dominant type in each grid cell
   public :: filter_same           ! build a filter of overlaps where src_val == dst_val
   public :: lookup_2d             ! create map based on a 2-d lookup table
   public :: lookup_2d_netcdf      ! wrapper to lookup_2d; first read table from netcdf file
   public :: which_max             ! get index of the maximum value in an array
!
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!EOP
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_dominant_indices
!
! !INTERFACE:
subroutine get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata, filter)
!
! !DESCRIPTION:
! Fills an output array on the destination grid (dst_array) whose values are equal to the
! (weighted) dominant value in the source grid cells overlapping a given destination grid
! cell
!
! Ignores all values in src_array that are less than minval or greater than maxval (treats
! those values the same as if they had wt=0).  (Note: for memory-use efficiency, it is
! best if the indices are designed such that most values between minval and maxval are
! actually used, since an array is allocated of size (maxval - minval + 1)*gridmap%nb.)
!
! The filter argument can be used to exclude certain overlaps -- if provided, we only
! consider overlaps where filter is .true. For example, see mkurbanparDomMod, where we
! first determine the dominant density class in each output cell, then determine the
! dominant region, but in the latter only considering overlapping source points whose
! density matches the dominant density in the output cell. If not provided, filter is
! treated as being .true. everywhere.
! 
! Output grid cells with no contributing valid source points are given the nodata value
!
! !USES:
   use mkgridmapMod, only : gridmap_type
!
! !ARGUMENTS:
   implicit none
   type(gridmap_type), intent(in) :: gridmap       ! provides mapping from src -> dst
   integer           , intent(in) :: src_array(:)  ! input values; length gridmap%na
   integer           , intent(out):: dst_array(:)  ! output values; length gridmap%nb
   integer           , intent(in) :: minval        ! minimum valid value in src_array
   integer           , intent(in) :: maxval        ! maximum valid value in src_array
   integer           , intent(in) :: nodata        ! value to assign to dst_array where there are no valid source points

   logical, intent(in), optional :: filter(:)  ! only consider overlaps where filter is .true.; length gridmap%ns
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   logical, allocatable  :: lfilter(:)    ! local version of filter
   logical, allocatable  :: hasdata(:)    ! true if an output cell has any valid data;
   real(r8), allocatable :: weights(:,:)  ! summed weight of each index value for each output cell

   integer  :: n, ni, no
   integer  :: k
   integer  :: maxindex
   real(r8) :: wt
   real(r8) :: maxwt

   character(len=*), parameter :: subname = "get_dominant_indices"
!-----------------------------------------------------------------------
   
   ! Error-check inputs and initialize local variables
   
   if (size(src_array) /= gridmap%na .or. &
        size(dst_array) /= gridmap%nb) then
      write(6,*) subname//' ERROR: incorrect sizes of src_array or dst_array'
      write(6,*) 'size(src_array) = ', size(src_array)
      write(6,*) 'gridmap%na      = ', gridmap%na
      write(6,*) 'size(dst_array) = ', size(dst_array)
      write(6,*) 'gridmap%nb      = ', gridmap%nb
      call abort()
   end if

   allocate(lfilter(gridmap%ns))

   if (present(filter)) then
      if (size(filter) /= gridmap%ns) then
         write(6,*) subname//' ERROR: incorrect size of filter'
         write(6,*) 'size(filter) = ', size(filter)
         write(6,*) 'gridmap%ns   = ', gridmap%ns
         call abort()
      end if

      lfilter(:) = filter(:)
   else
      lfilter(:) = .true.
   end if

   allocate(hasdata(gridmap%nb))
   hasdata(:) = .false.
   allocate(weights(minval:maxval, gridmap%nb))
   weights(minval:maxval,:) = 0.

   ! Determine weight of each index value for each output (destination) cell

   do n = 1, gridmap%ns
      if (lfilter(n)) then
         ni = gridmap%src_indx(n)
         no = gridmap%dst_indx(n)
         wt = gridmap%wovr(n)
         k = src_array(ni)
         if (k >= minval .and. k <= maxval) then
            ! Note: if we were doing something like weighted sums, I think we would
            ! want to divide wt by gridmap%frac_dst(no), as is done in
            ! gridmap_areaave_default. But since all we care about is the relative
            ! values of weights for a given destination cell, this is unnecessary
            weights(k,no) = weights(k,no) + wt
            hasdata(no) = .true.
         end if
      end if
   end do

   ! Determine output values
   ! Note: if a given destination cell has no contributing source points (thus
   ! hasdata(no) = false), or the max weight of any index overlapping this destination
   ! cell is <= 0, then the output value there will be nodata.
   ! (I don't think this latter condition -- weight <= 0 -- is possible, but we handle
   ! it anyway)

   dst_array(:) = nodata
   do no = 1, gridmap%nb
      if (hasdata(no)) then
         call which_max(weights(:,no), maxwt, maxindex, lbound=minval)
         if (maxwt > 0.) then
            dst_array(no) = maxindex
         end if
      end if
   end do

   deallocate(lfilter, weights, hasdata)

end subroutine get_dominant_indices
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: filter_same
!
! !INTERFACE:
subroutine filter_same(gridmap, filter, src_array, dst_array, nodata)
!
! !DESCRIPTION:
! Creates a filter of overlaps where src_array == dst_array.
!
! More specifically: given a src_array (of size gridmap%na) and an already-created
! dst_array (of size gridmap%nb):
!
! Creates a logical filter array, of size gridmap%ns (i.e., number of overlaps),
! according to the following rules:
! (1) anywhere where filter was already .false., it will remain .false.
! (2) if nodata is present: for any overlap where the value in dst_array is nodata,
!     filter will be .false.
! (3) for any overlap where the value in the given src_array differs from the value
!     in the given dst_array, filter will be .false.
! (4) anywhere else, filter will be .true.
!
! !USES:
   use mkgridmapMod, only : gridmap_type
!
! !ARGUMENTS:
   implicit none
   type(gridmap_type), intent(in)   :: gridmap      ! provides mapping from src -> dst
   logical           , intent(inout):: filter(:)    ! length gridmap%ns
   integer           , intent(in)   :: src_array(:) ! length gridmap%na
   integer           , intent(in)   :: dst_array(:) ! length gridmap%nb

   integer, intent(in), optional :: nodata  ! wherever dst_array == nodata, filter will be false
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer :: n, ni, no

   character(len=*), parameter :: subname = "make_filter"
!-----------------------------------------------------------------------
   
   ! Error check inputs
   
   if (size(filter) /= gridmap%ns .or. &
        size(src_array) /= gridmap%na .or. &
        size(dst_array) /= gridmap%nb) then
      write(6,*) subname//' ERROR: incorrect array sizes'
      write(6,*) 'size(src_array) = ', size(src_array)
      write(6,*) 'gridmap%na      = ', gridmap%na
      write(6,*) 'size(dst_array) = ', size(dst_array)
      write(6,*) 'gridmap%nb      = ', gridmap%nb
      write(6,*) 'size(filter)    = ', size(filter)
      write(6,*) 'gridmap%ns      = ', gridmap%ns
      call abort()
   end if

   ! Create the filter

   do n = 1, gridmap%ns
      ni = gridmap%src_indx(n)
      no = gridmap%dst_indx(n)

      if (present(nodata)) then
         if (dst_array(no) == nodata) then
            filter(n) = .false.
         end if
      end if

      if (dst_array(no) /= src_array(ni)) then
         filter(n) = .false.
      end if
   end do

end subroutine filter_same
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lookup_2d
!
! !INTERFACE:
subroutine lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, &
                     nodata, valid_entries, invalid_okay)
!
! !DESCRIPTION:
! Creates a data array using a paint-by-number approach according to a lookup table
!
! This routine operates on a 2-d lookup table. There are therefore two index arrays
! (index1 and index2); these index arrays are on the same grid as the desired data array
! (thus, index1, index2 and data must all have the same length).  Each output point, n, is
! then generally determined as:
!
! data(n) = lookup_table(index1(n), index2(n))
!      
! fill_val: value to put in data array where either:
! (a) index1 or index2 are equal to nodata (if nodata is given)
!     Note that this condition does NOT result in ierr being set
! (b) valid_entries(index1(n), index2(n)) is false (if valid_entries is given)
!     Note that this condition also results in ierr being set, unless invalid_okay is
!     present and .true.
!     (If valid_entries is not given, it is treated as being .true. everywhere)
! (c) index1 or index2 out of range
!     Note that this condition also results in ierr being set
! 
! ierr: error return code (if non-0, indicates first error encountered):
!    0: no error
!    1: attempt to assign values from the lookup table that are invalid according
!       to valid_entries (note: this is not considered an error if invalid_okay is
!       present and .true.)
!    2: attempt to access an out-of-range index in lookup table
! WJS (2-2-12): My main reason for using ierr rather than aborting in case of error
! is to facilitate unit testing
!
! !ARGUMENTS:
   implicit none
   integer , intent(in) :: index1(:)  ! index into dim 1 of lookup_table
   integer , intent(in) :: index2(:)  ! index into dim 2 of lookup_table
   real(r8), intent(in) :: lookup_table(:,:)
   real(r8), intent(in) :: fill_val   ! value to put in data where we don't have a valid value (see above for details)
   real(r8), intent(out):: data(:)    ! output arary
   integer , intent(out):: ierr       ! error return code (0 = no error)

   ! nodata flag in index1 and index2 (see above for details):
   integer, intent(in), optional :: nodata

   ! which entries are considered valid (see above for details):
   logical, intent(in), optional :: valid_entries(:,:)

   ! invalid_okay: if true, then assigning fill_val because valid_entries is false does
   ! NOT raise an error flag (invalid_okay defaults to false, meaning an error is
   ! raised in this case):
   logical, intent(in), optional :: invalid_okay
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer :: n
   integer :: i1, i2
   integer :: data_size          ! size of index1, index2 and data arrays
   integer :: table_n1           ! size of dimension 1 of lookup table
   integer :: table_n2           ! size of dimension 2 of lookup table
   logical :: linvalid_okay      ! local version of invalid_okay
   logical, allocatable :: lvalid_entries(:,:)  ! local version of valid_entries

   character(len=*), parameter :: subname = 'lookup_2d'
!-----------------------------------------------------------------------
   
   ierr = 0

   ! Error-check array sizes

   data_size = size(data)
   if (size(index1) /= data_size .or. size(index2) /= data_size) then
      write(6,*) subname//' ERROR: data array sizes do not match'
      write(6,*) 'size(data)   = ', data_size
      write(6,*) 'size(index1) = ', size(index1)
      write(6,*) 'size(index2) = ', size(index2)
      call abort()
   end if

   table_n1 = size(lookup_table,1)
   table_n2 = size(lookup_table,2)
   if (present(valid_entries)) then
      if (size(valid_entries,1) /= table_n1 .or. size(valid_entries,2) /= table_n2) then
         write(6,*) subname//' ERROR: size of valid_entries does not match lookup_table'
         write(6,*) 'size(lookup_table)  = ', table_n1, table_n2
         write(6,*) 'size(valid_entries) = ', size(valid_entries,1), &
              size(valid_entries,2)
         call abort()
      end if
   end if

   ! Set local version of invalid_okay & valid_entries

   if (present(invalid_okay)) then
      linvalid_okay = invalid_okay
   else
      linvalid_okay = .false.
   end if

   allocate(lvalid_entries(table_n1, table_n2))
   if (present(valid_entries)) then
      lvalid_entries(:,:) = valid_entries(:,:)
   else
      lvalid_entries(:,:) = .true.
   end if

   ! Do the lookups

   do n = 1, data_size
      i1 = index1(n)
      i2 = index2(n)

      ! First handle special cases:

      ! index is nodata flag (this is NOT an error)
      if (present(nodata)) then
         if (i1 == nodata .or. i2 == nodata) then
            data(n) = fill_val
            cycle
         end if
      end if

      ! index out of range
      if (i1 <= 0 .or. i1 > table_n1 .or. &
           i2 <= 0 .or. i2 > table_n2) then
         data(n) = fill_val
         if (ierr == 0) ierr = 2
         cycle
      end if

      ! lookup table entry is invalid
      if (.not. lvalid_entries(i1, i2)) then
         data(n) = fill_val
         if (.not. linvalid_okay) then
            if (ierr == 0) ierr = 1
         end if
         cycle
      end if

      ! Finally, the "normal" case, if none of the special cases were triggered:
      data(n) = lookup_table(i1, i2)
   end do

   deallocate(lvalid_entries)

end subroutine lookup_2d
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lookup_2d_netcdf
!
! !INTERFACE:
subroutine lookup_2d_netcdf(ncid, tablename, lookup_has_invalid, &
                            dimname1, dimname2, n_extra_dims, &
                            index1, index2, fill_val, data, ierr, &
                            extra_dims, nodata, invalid_okay)
!
! !DESCRIPTION:
! Wrapper to lookup_2d that first reads the lookup table from a netcdf file
!
! If lookup_has_invalid is false, then we treat all lookup table entries as valid data
! (i.e., all valid_entries are true in the call to lookup_2d). If lookup_has_invalid is
! true, then we read the _FillValue attribute for the lookup table variable, and consider
! any table entry with value _FillValue to be an invalid entry, thus putting fill_val in
! these data locations (and raising an error flag unless invalid_okay is present and
! true).
!
! The dimension given by dimname1 -- with the associated indices given by index1 -- is the
! fastest-varying dimension in the lookup table. Dimension dimname2 (associated with
! index2) is the second-fastest-varying dimension. Similarly, extra_dims should be ordered
! from faster-varying to slowest-varying dimension.  (The first dimension in extra_dims is
! the third-fastest-varying dimension in the lookup table.)
!
! n_extra_dims gives the number of extra dimensions (in addition to the first two) in the
! lookup table. We take a single 2-d slice of the lookup table, by using a single value of
! each of these other dimensions. If n_extra_dims > 0, then extra_dims must be present,
! with at least n_extra_dims entries. Each entry in extra_dims gives the name of a
! dimension and the dimension index to use for the slice.
!
! If size(extra_dims) > n_extra_dims, then we use the first n_extra_dims entries in
! extra_dims. If n_extra_dims = 0, then extra_dims is ignored.
!
! Note that we ignore any coordinate variables associated with the dimensions of the
! lookup table; we simply treat the lookup table indices as 1,2,3,...
!
! See the lookup_2d documentation for documentation of some other arguments
!
! WJS (2-1-12): Some thoughts on avoiding duplication if we eventually want similar
! routines, lookup_1d_netcdf, lookup_3d_netcdf, etc.:
!
! Much of the code in lookup_2d_netcdf could then be pulled out to a shared subroutine
! (e.g., much of the error-checking code).
!
! Or, maybe better: we could try to make a single lookup_netcdf subroutine that handles
! 1-d, 2-d and any other dimensionality. To do that, we would (1) make a generic interface
! (of which lookup_1d and lookup_2d would be implementations); (2) change the repeated
! arguments in lookup_2d_netcdf (*1 and *2) to arrays -- maybe using an array of a derived
! type containing these arguments; (3) if possible, initially read the lookup table into a
! 1-d array (if the netcdf call allows reading a n-d array into a 1-d array) (if netcdf
! doesn't allow this, then I think we could achieve the same thing by reading 1-d slices
! of the lookup table in a loop, building the full lookup table as a long 1-d array); (4)
! in the call to the generic 'lookup' function, reshape the 1-d lookup table
! appropriately. (Note: I think it would be challenging to combine lookup_1d and lookup_2d
! (etc.)  into a single routine using a similar method.)
!
! !USES:
   use mkncdio
! !ARGUMENTS:
   implicit none
   integer         , intent(in) :: ncid         ! ID of an open netcdf file
   character(len=*), intent(in) :: tablename    ! name of the lookup table variable
   logical         , intent(in) :: lookup_has_invalid ! should we use _FillValue? (see above)
   character(len=*), intent(in) :: dimname1     ! name of the first (fastest-varying) dimension of the lookup table
   character(len=*), intent(in) :: dimname2     ! name of the second dimension of the lookup table
   integer         , intent(in) :: n_extra_dims ! number of extra dimensions in the lookup table
   ! The following arguments are passed directly to lookup_2d:
   integer         , intent(in) :: index1(:)    ! index into dim 1 of lookup table
   integer         , intent(in) :: index2(:)    ! index into dim 2 of lookup table
   real(r8)        , intent(in) :: fill_val     ! value to put in data where we don't have a valid value
   real(r8)        , intent(out):: data(:)      ! output array
   integer         , intent(out):: ierr         ! error return code from the call to lookup_2d

   ! slice to use if lookup table variable has more than 2 dimensions:
   type(dim_slice_type), intent(in), optional :: extra_dims(:)

   ! nodata flag in index1 and index2, passed directly to lookup_2d:
   integer             , intent(in), optional :: nodata

   ! flag for whether trying to use a lookup table value that is equal to the _FillValue
   ! should raise an error flag
   ! (irrelevant if lookup_has_invalid is .false.)
   ! (passed directly to lookup_2d - see the documentation there for more details)
   logical             , intent(in), optional :: invalid_okay
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer  :: varid           ! netcdf variable id of the lookup table
   integer  :: ndims           ! total number of dimensions of lookup table
   integer  :: ndims_expected  ! value we expect for ndims, for error checking
   integer  :: i
   real(r8) :: table_fillval   ! value of the _FillValue attribute for the lookup table
   character(len=nf_max_name), allocatable :: dimnames(:)  ! dimension names
   integer , allocatable :: dimids(:)  ! dimension ids
   integer , allocatable :: dimlens(:) ! dimension lengths
   integer , allocatable :: starts(:)  ! starting indices for reading lookup table
   integer , allocatable :: counts(:)  ! dimension counts for reading lookup table
   real(r8), allocatable :: lookup_table(:,:)
   logical , allocatable :: valid_entries(:,:)  ! which entries of the lookup table are considered valid

   character(len=*), parameter :: subname = 'lookup_2d_netcdf'
!-----------------------------------------------------------------------
   
   ! Error-check extra_dims
   
   if (n_extra_dims > 0) then
      if (.not. present(extra_dims)) then
         write(6,*) subname//' ERROR: extra_dims must be present for n_extra_dims > 0'
         call abort()
      end if

      if (size(extra_dims) < n_extra_dims) then
         write(6,*) subname//' ERROR: not enough extra dimensions given'
         write(6,*) 'n_extra_dims =     ', n_extra_dims
         write(6,*) 'size(extra_dims) = ', size(extra_dims)
         call abort()
      end if
   end if

   ! Determine number of expected dimensions in the table, and actual number of
   ! dimensions in the netcdf file

   ndims_expected = 2 + n_extra_dims

   call check_ret(nf_inq_varid (ncid, tablename, varid), subname)
   call check_ret(nf_inq_varndims (ncid, varid, ndims), subname)
   if (ndims /= ndims_expected) then
      write(6,*) subname//' ERROR: unexpected number of dimensions in ', &
           trim(tablename)
      write(6,*) 'ndims = ', ndims
      write(6,*) 'expected (based on n_extra_dims): ', ndims_expected
      call abort()
   end if

   ! Get dimension names & sizes, and error-check them

   allocate(dimids(ndims), dimlens(ndims), dimnames(ndims))
   call check_ret(nf_inq_vardimid (ncid, varid, dimids), subname)
   do i = 1, ndims
      call check_ret(nf_inq_dimname (ncid, dimids(i), dimnames(i)), subname)
      call check_ret(nf_inq_dimlen (ncid, dimids(i), dimlens(i)), subname)
   end do

   call check_dimname(dimnames(1), dimname1, 1)
   call check_dimname(dimnames(2), dimname2, 2)
   do i = 1, n_extra_dims
      call check_dimname(dimnames(2+i), extra_dims(i)%name, 2+i)
      call check_dimsize(dimlens(2+i), extra_dims(i)%val, 2+i)
   end do

   ! Read the lookup table; if the given variable has more than 2 dimensions, we read
   ! a single 2-d slice

   allocate(starts(ndims), counts(ndims))
   allocate(lookup_table(dimlens(1), dimlens(2)))
   starts(1:2) = 1
   counts(1:2) = dimlens(1:2)
   do i = 1, n_extra_dims
      starts(2+i) = extra_dims(i)%val
      counts(2+i) = 1
   end do
   call check_ret(nf_get_vara_double (ncid, varid, starts, counts, lookup_table), subname)

   ! Determine which entries are valid

   allocate(valid_entries(size(lookup_table, 1), size(lookup_table, 2)))
   valid_entries(:,:) = .true.
   if (lookup_has_invalid) then
      call check_ret(nf_get_att_double (ncid, varid, '_FillValue', table_fillval), subname)
      where (lookup_table == table_fillval)
         valid_entries = .false.
      end where
   end if

   ! Do the lookups

   call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, nodata=nodata, &
        valid_entries=valid_entries, invalid_okay=invalid_okay)

   deallocate(valid_entries)
   deallocate(lookup_table)
   deallocate(starts, counts)
   deallocate(dimids, dimlens, dimnames)

contains
!------------------------------------------------------------------------------
   subroutine check_dimname(actual, expected, i)
   ! Make sure names are equal; if not, stop with an error message

      character(len=*), intent(in) :: actual, expected
      integer         , intent(in) :: i  ! dimension number, for output purposes

      if (actual /= expected) then
         write(6,*) subname//' ERROR: unexpected dimension name in ', trim(tablename)
         write(6,*) 'dimension #', i
         write(6,*) 'actual:   ', trim(actual)
         write(6,*) 'expected: ', trim(expected)
         call abort()
      end if
   end subroutine check_dimname

!------------------------------------------------------------------------------
   subroutine check_dimsize(length, index, i)
   ! Make sure dimension length is long enough; if not, stop with an error message

      integer, intent(in) :: length, index
      integer, intent(in) :: i  ! dimension number, for output purposes

      if (index > length) then
         write(6,*) subname//' ERROR: desired index exceeds dimension length in ', &
              trim(tablename)
         write(6,*) 'dimension #', i
         write(6,*) 'index:  ', index
         write(6,*) 'length: ', length
         call abort()
      end if
   end subroutine check_dimsize
      
end subroutine lookup_2d_netcdf
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: which_max
!
! !INTERFACE:
subroutine which_max(arr, maxval, maxindex, lbound)
!
! !DESCRIPTION:
! Returns maximum value in arr along with the index of the maximum value
!
! If multiple values are tied, returns index of the first maximum
!
! !ARGUMENTS:
   implicit none
   real(r8), intent(in) :: arr(:)
   real(r8), intent(out):: maxval   ! maximum value in arr(:)
   integer , intent(out):: maxindex ! first index of maxval

   ! lower bound of indices of arr; if not supplied, assumed to be 1:
   integer , intent(in), optional :: lbound
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
   integer :: i
!-----------------------------------------------------------------------
   
   maxindex = 1
   maxval = arr(1)

   do i = 2, size(arr)
      if (arr(i) > maxval) then
         maxindex = i
         maxval = arr(i)
      end if
   end do

   if (present(lbound)) then
      maxindex = maxindex + (lbound - 1)
   end if
end subroutine which_max
!------------------------------------------------------------------------------

end module mkindexmapMod
