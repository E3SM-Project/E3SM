module test_mkindexmapMod
! Module for testing mkindexmapMod

   use mkindexmapMod
   use test_mod
   use shr_kind_mod, only : r8 => shr_kind_r8

   implicit none
   private

   public :: test_get_dominant_indices
   public :: test_filter_same
   public :: test_lookup_2d
   public :: test_lookup_2d_netcdf
   public :: test_which_max   

   character(len=*), parameter :: modname = 'test_mkindexmapMod'

contains
   
!------------------------------------------------------------------------------
   subroutine test_get_dominant_indices

      use mkgridmapMod, only : gridmap_type

      implicit none

      type(gridmap_type) :: gridmap
      character(len=128) :: testname

      integer, allocatable :: src_array(:)
      integer, allocatable :: dst_array(:)
      integer, allocatable :: dst_array_t(:)
      logical, allocatable :: filter(:)
      integer :: minval, maxval, nodata

      character(len=*), parameter :: subname = 'test_get_dominant_indices'

      ! Set up a gridmap that will be used for most tests, and allocate corresponding
      ! arrays:
      ! Note that, for most tests here, the test arrays are: (1) simple case, (2) the main
      ! case to test, (3) simple case. Thus, the main case in question is #2 of 3, and
      ! we're always basically just testing one scenario in each call to the subroutine
      ! (rather than doing a bunch of tests at once, which could make setting up the test
      ! arrays more error-prone).

      ! This gridmap will have 3 src cells, 9 dest cells, and:
      ! src 1: just overlaps with dst 1
      ! src 2: overlaps with dst 1 & dst 2
      ! src 3..7: just overlaps with dst 2
      ! src 8: overlaps with dst 2 & dst 3
      ! src 9: just overlaps with dst 3
      ! Note: I'm not setting some things that aren't used in get_dominant_indices
      gridmap%na = 9
      gridmap%nb = 3
      gridmap%ns = 11
      allocate(gridmap%src_indx(gridmap%ns), &
               gridmap%dst_indx(gridmap%ns), &
               gridmap%wovr    (gridmap%ns))
      gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
      gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
      gridmap%wovr     = (/0.75,0.25, &  ! weights of sources 1:2 on dest 1
                           0.1,0.1,0.1,0.3,0.2,0.2,0.2, &  ! weights of sources 2:8 on dest 2
                           0.25,0.75/)   ! weights of sources 8:9 on test 3
      allocate(src_array  (gridmap%na), &
               dst_array  (gridmap%nb), &
               dst_array_t(gridmap%nb), &
               filter     (gridmap%ns))

      testname = 'basic test, all unique'
      src_array = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
      minval = 1
      maxval = 9
      nodata = -1
      ! dst 2 takes its value from src 5 because it has the largest weight:
      dst_array_t = (/1, 5, 9/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))
      
      testname = 'basic test, some duplicates'
      src_array = (/1, 2, 3, 3, 4, 2, 2, 1, 1/)
      minval = 1
      maxval = 4
      nodata = -1
      dst_array_t = (/1, 2, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'minval not 1'
      src_array = (/3, 4, 5, 5, 6, 4, 4, 3, 3/)
      minval = 3
      maxval = 6
      nodata = -1
      dst_array_t = (/3, 4, 3/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'single non-zero source value'
      src_array = (/1, 0, 0, 0, 0, 2, 0, 0, 1/)
      minval = 1
      maxval = 2
      nodata = -1
      dst_array_t = (/1, 2, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'single value within given min-max range'
      src_array = (/1, 0, 9, 9, 0, 2, 9, 9, 1/)
      minval = 1
      maxval = 2
      nodata = -1
      dst_array_t = (/1, 2, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))
      
      testname = 'no valid values'
      src_array = (/1, 0, 9, 9, 0, 0, 9, 9, 1/)
      minval = 1
      maxval = 2
      nodata = -1
      dst_array_t = (/1, nodata, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'some filters false'
      src_array = (/1, 2, 3, 3, 4, 2, 2, 1, 1/)
      minval = 1
      maxval = 4
      nodata = -1
      filter = (/.true., .true., &
                 .false., .true., .true., .true., .false., .true., .true., &
                 .true., .true./)
      dst_array_t = (/1, 4, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata, filter=filter)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))
      
      testname = 'all filters false'
      src_array = (/1, 2, 3, 3, 4, 2, 2, 1, 1/)
      minval = 1
      maxval = 4
      nodata = -1
      filter = (/.true., .true., &
                 .false., .false., .false., .false., .false., .false., .false., &
                 .true., .true./)      
      dst_array_t = (/1, nodata, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata, filter=filter)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      ! Modify gridmap weights for the following test
      gridmap%wovr     = (/0.75,0.25, &  ! weights of sources 1:2 on dest 1
                           0.0,0.0,0.0,0.0,0.0,0.0,0.0, &  ! weights of sources 2:8 on dest 2
                           0.25,0.75/)   ! weights of sources 8:9 on test 3
      testname='all weights 0'
      src_array = (/1, 1, 1, 1, 1, 1, 1, 1, 1/)
      minval = 1
      maxval = 2
      nodata = -1
      dst_array_t = (/1, nodata, 1/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))
      
      ! Make a new gridmap for the following test;
      ! this involves more output cells and a more complex mapping from src to dst
      ! This gridmap will have:
      ! dst 1: from src 1, 4, 7
      ! dst 2: from src 2, 4, 6
      ! dst 3: from src 1
      ! dst 4: no overlapping src cells
      ! dst 5: from src 5, 7, 8
      ! note that src 3 & 9 do not overlap with any dst
      deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, & 
                 src_array, dst_array, dst_array_t, filter)
      gridmap%na = 9
      gridmap%nb = 5
      gridmap%ns = 10
      allocate(gridmap%src_indx(gridmap%ns), &
               gridmap%dst_indx(gridmap%ns), &
               gridmap%wovr    (gridmap%ns))
      gridmap%src_indx = (/1, 2, 4, 4, 7, 6, 1, 5, 7, 8/)
      gridmap%dst_indx = (/1, 2, 1, 2, 1, 2, 3, 5, 5, 5/)
      gridmap%wovr     = (/1, 1, 2, 2, 1, 3, 1, 2, 2, 3/)
      allocate(src_array  (gridmap%na), &
               dst_array  (gridmap%nb), &
               dst_array_t(gridmap%nb), &
               filter     (gridmap%ns))

      testname = 'more complex gridmap'
      ! src index:  1  2  3  4  5  6  7  8  9
      src_array = (/1, 2, 3, 1, 5, 6, 5, 8, 9/)
      minval = 1
      maxval = 9
      nodata = -1
      dst_array_t = (/1, 6, 1, nodata, 5/)
      call get_dominant_indices(gridmap, src_array, dst_array, minval, maxval, nodata)
      call test_is(dst_array, dst_array_t, modname//' -- '//subname//' -- '//trim(testname))

      deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, & 
                 src_array, dst_array_t, filter)

   end subroutine test_get_dominant_indices
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   subroutine test_filter_same
      
      use mkgridmapMod, only : gridmap_type

      implicit none

      type(gridmap_type) :: gridmap
      character(len=128) :: testname
      
      integer, allocatable :: src_array(:)
      integer, allocatable :: dst_array(:)
      logical, allocatable :: filter(:)
      logical, allocatable :: filter_t(:)
      integer :: nodata

      character(len=*), parameter :: subname = 'test_filter_same'
      
      ! Set up a gridmap that will be used for most tests, and allocate corresponding
      ! arrays:
      ! Note that, for most tests here, the test arrays are: (1) simple case, (2) the main
      ! case to test, (3) simple case. Thus, the main case in question is #2 of 3, and
      ! we're always basically just testing one scenario in each call to the subroutine
      ! (rather than doing a bunch of tests at once, which could make setting up the test
      ! arrays more error-prone).

      ! This gridmap will have 3 src cells, 9 dest cells, and:
      ! src 1: just overlaps with dst 1
      ! src 2: overlaps with dst 1 & dst 2
      ! src 3..7: just overlaps with dst 2
      ! src 8: overlaps with dst 2 & dst 3
      ! src 9: just overlaps with dst 3
      ! Note: I'm not setting some things that aren't used in filter_same
      gridmap%na = 9
      gridmap%nb = 3
      gridmap%ns = 11
      allocate(gridmap%src_indx(gridmap%ns), &
               gridmap%dst_indx(gridmap%ns))
      gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
      gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
      allocate(src_array  (gridmap%na), &
               dst_array  (gridmap%nb), &
               filter     (gridmap%ns), &
               filter_t   (gridmap%ns))

      testname = 'maintain false values in filter'
      src_array(:) = 1
      dst_array(:) = 1
      filter(:) = .true.
      filter(3) = .false.
      filter(5) = .false.
      filter_t(:) = .true.
      filter_t(3) = .false.
      filter_t(5) = .false.
      call filter_same(gridmap, filter, src_array, dst_array)
      call test_is(filter, filter_t, modname//' -- '//subname//' -- '//trim(testname))
      
      testname = 'dst_array = nodata in some places'
      nodata = -1
      src_array(:) = 1
      src_array(5) = nodata  ! make sure that even when src_array = dst_array = nodata,
                             ! we still end up with filter = false
      dst_array = (/1, nodata, 1/)
      filter(:) = .true.
      filter_t(:) = .true.
      filter_t(3:9) = .false.  ! false for all overlaps with dst #2
      call filter_same(gridmap, filter, src_array, dst_array, nodata=nodata)
      call test_is(filter, filter_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'src_array not equal to dst_array in some places, no nodata argument'
      src_array(:) = (/1, 1, 1, 1, 2, 3, 1, 3, 1/)
      dst_array(:) = (/1, 1, 1/)
      filter(:) = .true.
      ! src_array index: 1      2      2      3      4       5       6      7       8       8      9
      filter_t(:) = (/.true.,.true.,.true.,.true.,.true.,.false.,.false.,.true.,.false.,.false.,.true./)
      call filter_same(gridmap, filter, src_array, dst_array)
      call test_is(filter, filter_t, modname//' -- '//subname//' -- '//trim(testname))

      testname = 'src_array not equal to dst_array in some places, nodata never applies'
      nodata = -1
      src_array(:) = (/1, 1, 1, 1, 2, 3, 1, 3, 1/)
      dst_array(:) = (/1, 1, 1/)
      filter(:) = .true.
      ! src_array index: 1      2      2      3      4       5       6      7       8       8      9
      filter_t(:) = (/.true.,.true.,.true.,.true.,.true.,.false.,.false.,.true.,.false.,.false.,.true./)
      call filter_same(gridmap, filter, src_array, dst_array, nodata=nodata)
      call test_is(filter, filter_t, modname//' -- '//subname//' -- '//trim(testname))
      
      testname = 'combination of false filter, src_array not equal to dst_array, and nodata'
      nodata = -1
      src_array(:) = (/1, 2, 1, 2, 1, 2, 1, 2, 1/)
      dst_array(:) = (/nodata, 1, 1/)
      filter(:) = .true.
      filter(4) = .false.
      filter_t(:) = (/.false.,.false.,.false.,.false.,.false.,.true.,.false.,.true.,.false.,.false.,.true./)
      call filter_same(gridmap, filter, src_array, dst_array, nodata=nodata)
      call test_is(filter, filter_t, modname//' -- '//subname//' -- '//trim(testname))


      deallocate(gridmap%src_indx, gridmap%dst_indx, & 
                 src_array, dst_array, filter, filter_t)

   end subroutine test_filter_same
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   subroutine test_lookup_2d

      implicit none

      character(len=128)    :: testname
      real(r8), allocatable :: lookup_table(:,:)
      logical , allocatable :: valid_entries(:,:)
      integer , allocatable :: index1(:), index2(:)
      real(r8), allocatable :: data(:), data_t(:)
      real(r8) :: fill_val
      integer  :: nodata
      integer  :: ierr, ierr_t

      character(len=*), parameter :: subname = 'test_lookup_2d'

      ! Create lookup table for use in most tests
      allocate(lookup_table(2,3), valid_entries(2,3))
      lookup_table(1,:) = (/11.,12.,13./)
      lookup_table(2,:) = (/21.,22.,23./)

      testname = 'basic test; no nodata or valid_entries'
      allocate(index1(5), index2(5), data(5), data_t(5))
      index1 = (/1,2,1,2,2/)
      index2 = (/1,2,3,2,3/)
      fill_val = -1.
      data_t = (/11., 22., 13., 22., 23./)
      ierr_t = 0
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr)
      call check_results
      deallocate(index1, index2, data, data_t)

      testname = 'basic test but with index out of range'
      allocate(index1(5), index2(5), data(5), data_t(5))
      index1 = (/1,2,3,2,2/)
      index2 = (/1,2,1,2,4/)
      fill_val = -1.
      data_t = (/11._r8, 22._r8, fill_val, 22._r8, fill_val/)
      ierr_t = 2
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr)
      call check_results
      deallocate(index1, index2, data, data_t)
      
      testname = 'basic test but with nodata present, and a nodata value in input'
      allocate(index1(5), index2(5), data(5), data_t(5))
      nodata = -1
      index1 = (/nodata,2,1,2,nodata/)
      index2 = (/1,2,3,nodata,nodata/)
      fill_val = -1.
      data_t = (/fill_val, 22._r8, 13._r8, fill_val, fill_val/)
      ierr_t = 0
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, nodata=nodata)
      call check_results
      deallocate(index1, index2, data, data_t)

      testname = 'valid_entries'
      allocate(index1(5), index2(5), data(5), data_t(5))
      index1 = (/1,1,2,2,1/)
      index2 = (/1,2,1,2,3/)
      valid_entries(1,:) = (/.false.,.false.,.true./)
      valid_entries(2,:) = (/.true. ,.true. ,.true./)
      fill_val = -1.
      data_t = (/fill_val, fill_val, 21._r8, 22._r8, 13._r8/)
      ierr_t = 1
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, valid_entries=valid_entries)
      call check_results

      testname = 'valid_entries, invalid_okay'
      ! Note: this test reuses some setup from the previous test
      ierr_t = 0
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, &
           valid_entries=valid_entries, invalid_okay=.true.)
      call check_results      
      deallocate(index1, index2, data, data_t)      
      

      testname = 'valid_entries, together with index out of range'
      ! in addition to checking both valid_entries and index out of range, this also
      ! makes sure that we get the appropriate ierr value when we have both errors
      ! (because we encounter the valid_entries error first)
      allocate(index1(5), index2(5), data(5), data_t(5))
      index1 = (/1,1,3,2,2/)
      index2 = (/1,2,1,1,0/)
      valid_entries(1,:) = (/.false.,.false.,.true./)
      valid_entries(2,:) = (/.true. ,.true. ,.true./)
      fill_val = -1.
      data_t = (/fill_val, fill_val, fill_val, 21._r8, fill_val/)
      ierr_t = 1
      call lookup_2d(index1, index2, lookup_table, fill_val, data, ierr, valid_entries=valid_entries)
      call check_results
      deallocate(index1, index2, data, data_t)      


      deallocate(lookup_table, valid_entries)

   contains
      subroutine check_results
         call test_is(data, data_t, modname//' -- '//subname//' -- '//trim(testname)//' -- data')
         call test_is(ierr, ierr_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ierr')
      end subroutine check_results

   end subroutine test_lookup_2d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   subroutine test_lookup_2d_netcdf

      use mkncdio

      implicit none

      character(len=128) :: testname
      character(len=64)  :: tablename
      character(len=4)   :: dimname1, dimname2
      logical :: invalid_lookup
      integer :: n_extra_dims
      integer , allocatable :: index1(:), index2(:)
      real(r8), allocatable :: data(:), data_t(:)
      real(r8) :: fill_val
      integer  :: nodata
      integer  :: ierr, ierr_t
      type(dim_slice_type), allocatable :: extra_dims(:)
      
      integer :: ncid
      character(len=*), parameter :: filename = 'unit_testers/inputs/test_lookup_2d_netcdf.nc'

      ! flags to enable tests that we don't usually want to run, because they result in
      ! an abort, but we may occasionally want to run to make sure this error-handling is
      ! working properly
      logical, parameter :: test_abort1 = .false.
      logical, parameter :: test_abort2 = .false.
      logical, parameter :: test_abort3 = .false.
      
      character(len=*), parameter :: subname = 'test_lookup_2d_netcdf'

      ! Open netcdf file that will be used for most tests:
      ! Note that this file was created such that lookup4d(i,j,k,l) = 1000*i+100*j+10*k+l,
      ! and similarly for the other variables
      ! Also, lookup2d(1,2) is missing (i.e., equal to the _FillVal)
      call check_ret(nf_open(filename, 0, ncid), subname)

      testname = '2-d lookup table with _FillValue resulting in valid_entries false somewhere'
      allocate(index1(5), index2(5), data(5), data_t(5))
      tablename = 'lookup2d'
      invalid_lookup = .true.
      dimname1 = 'dim1'
      dimname2 = 'dim2'
      n_extra_dims = 0
      index1 = (/1,2,1,2,2/)
      index2 = (/1,2,2,1,3/)
      fill_val = -1.
      ! Note that the third value is fill_val because lookup2d(1,2) is missing (i.e.,
      ! equal to the _FillVal in the netcdf file)
      data_t = (/11._r8, 22._r8, fill_val, 21._r8, 23._r8/)
      ierr_t = 1
      call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
           n_extra_dims, index1, index2, fill_val, data, ierr)
      call check_results

      testname = '2-d lookup table with _FillValue resulting in valid_entries false somewhere, invalid_okay'
      ! Note: this test reuses some setup from the previous test
      ierr_t = 0
      call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
           n_extra_dims, index1, index2, fill_val, data, ierr, invalid_okay=.true.)
      call check_results
      deallocate(index1, index2, data, data_t)

      testname = '3-d lookup table with no _FillValue; nodata in index arrays'
      allocate(index1(5), index2(5), data(5), data_t(5))
      tablename = 'lookup3d'
      invalid_lookup = .false.
      dimname1 = 'dim1'
      dimname2 = 'dim2'
      n_extra_dims = 1
      allocate(extra_dims(n_extra_dims))
      extra_dims(1) = dim_slice_type('dim3', 2)
      nodata = -999
      index1 = (/nodata,2,1,2,2/)
      index2 = (/1,2,2,1,nodata/)
      fill_val = -1.
      data_t = (/fill_val, 222._r8, 122._r8, 212._r8, fill_val/)
      ierr_t = 0
      call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
           n_extra_dims, index1, index2, fill_val, data, ierr, extra_dims=extra_dims, &
           nodata=nodata)
      call check_results
      deallocate(index1, index2, data, data_t, extra_dims)
      
      testname = '4-d lookup table'
      allocate(index1(5), index2(5), data(5), data_t(5))
      tablename = 'lookup4d'
      invalid_lookup = .true.
      dimname1 = 'dim1'
      dimname2 = 'dim2'
      n_extra_dims = 2
      allocate(extra_dims(n_extra_dims))
      extra_dims(1) = dim_slice_type('dim3', 4)
      extra_dims(2) = dim_slice_type('dim4', 5)
      index1 = (/1,2,1,2,2/)
      index2 = (/1,2,2,1,3/)
      fill_val = -1.
      data_t = (/1145., 2245., 1245., 2145., 2345./)
      ierr_t = 0
      call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
           n_extra_dims, index1, index2, fill_val, data, ierr, extra_dims=extra_dims)
      call check_results
      deallocate(index1, index2, data, data_t, extra_dims)      

      ! The following tests should result in the code aborting with an error message.  
      !
      ! We don't usually want to run these tests, because they result in the code
      ! aborting, but we may want to run them occasionally to make sure this
      ! error-handling is working correctly.

      if (test_abort1) then
         testname = '2-d lookup table with incorrect dimname for dimension 2'
         allocate(index1(5), index2(5), data(5), data_t(5))
         tablename = 'lookup2d'
         invalid_lookup = .true.
         dimname1 = 'dim1'
         dimname2 = 'bad2'  ! this differs from the value in the file
         n_extra_dims = 0
         index1 = (/1,2,1,2,2/)
         index2 = (/1,2,2,1,3/)
         fill_val = -1.
         ! Note that the third value is fill_val because lookup2d(1,2) is missing (i.e.,
         ! equal to the _FillVal in the netcdf file)
         data_t = (/11._r8, 22._r8, fill_val, 21._r8, 23._r8/)
         ierr_t = 1
         call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
              n_extra_dims, index1, index2, fill_val, data, ierr)
         deallocate(index1, index2, data, data_t)  
      end if

      if (test_abort2) then
         testname = '3-d lookup table with incorrect dimname for dimension 3'
         allocate(index1(5), index2(5), data(5), data_t(5))
         tablename = 'lookup3d'
         invalid_lookup = .false.
         dimname1 = 'dim1'
         dimname2 = 'dim2'
         n_extra_dims = 1
         allocate(extra_dims(n_extra_dims))
         extra_dims(1) = dim_slice_type('bad3', 2)  ! this name differs from the value in the file
         nodata = -999
         index1 = (/nodata,2,1,2,2/)
         index2 = (/1,2,2,1,nodata/)
         fill_val = -1.
         data_t = (/fill_val, 222._r8, 122._r8, 212._r8, fill_val/)
         ierr_t = 0
         call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
              n_extra_dims, index1, index2, fill_val, data, ierr, extra_dims=extra_dims, &
              nodata=nodata)
         deallocate(index1, index2, data, data_t, extra_dims)
      end if

      if (test_abort3) then
         testname = '3-d lookup table, trying to access too large index for dimension 3'
         allocate(index1(5), index2(5), data(5), data_t(5))
         tablename = 'lookup3d'
         invalid_lookup = .false.
         dimname1 = 'dim1'
         dimname2 = 'dim2'
         n_extra_dims = 1
         allocate(extra_dims(n_extra_dims))
         extra_dims(1) = dim_slice_type('dim3', 5)  ! this index is out of bounds
         nodata = -999
         index1 = (/nodata,2,1,2,2/)
         index2 = (/1,2,2,1,nodata/)
         fill_val = -1.
         data_t = (/fill_val, 222._r8, 122._r8, 212._r8, fill_val/)
         ierr_t = 0
         call lookup_2d_netcdf(ncid, tablename, invalid_lookup, dimname1, dimname2, &
              n_extra_dims, index1, index2, fill_val, data, ierr, extra_dims=extra_dims, &
              nodata=nodata)
         deallocate(index1, index2, data, data_t, extra_dims)
      end if

      call check_ret(nf_close(ncid), subname)

   contains
      subroutine check_results
         call test_is(data, data_t, modname//' -- '//subname//' -- '//trim(testname)//' -- data')
         call test_is(ierr, ierr_t, modname//' -- '//subname//' -- '//trim(testname)//' -- ierr')
      end subroutine check_results

   end subroutine test_lookup_2d_netcdf
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   subroutine test_which_max

      implicit none

      real(r8), dimension(:), allocatable :: arr

      character(len=128) :: testname

      real(r8) :: maxval, maxval_t
      integer  :: maxindex, maxindex_t

      character(len=*), parameter :: subname = 'test_which_max'


      testname = 'length-1 array'
      allocate(arr(1))
      arr = (/3.0/)
      maxval_t = 3.0
      maxindex_t = 1
      call which_max(arr, maxval, maxindex)
      call check_results
      deallocate(arr)

      testname = 'max @ 1'
      allocate(arr(5))
      arr = (/5.0, 2.0, 3.0, 2.5, 1.5/)
      maxval_t = 5.0
      maxindex_t = 1
      call which_max(arr, maxval, maxindex)
      call check_results
      deallocate(arr)

      testname = 'max in middle'
      allocate(arr(5))
      arr = (/1.0, 2.0, 3.0, 2.5, 1.5/)
      maxval_t = 3.0
      maxindex_t = 3
      call which_max(arr, maxval, maxindex)
      call check_results
      deallocate(arr)

      testname = 'max at end'
      allocate(arr(5))
      arr = (/1.0, 2.0, 3.0, 2.5, 8.0/)
      maxval_t = 8.0
      maxindex_t = 5
      call which_max(arr, maxval, maxindex)
      call check_results
      deallocate(arr)

      testname = 'multiple tied max values'
      allocate(arr(5))
      arr = (/1.0, 3.0, 3.0, 2.5, 1.5/)
      maxval_t = 3.0
      maxindex_t = 2
      call which_max(arr, maxval, maxindex)
      call check_results
      deallocate(arr)

      testname = 'max in middle, with lbound present'
      allocate(arr(3:7))
      arr = (/1.0, 3.0, 10.0, 2.5, 8.0/)
      maxval_t = 10.0
      maxindex_t = 5
      call which_max(arr, maxval, maxindex, lbound=3)
      call check_results
      deallocate(arr)

   contains
      subroutine check_results
         call test_is(maxval, maxval_t, modname//' -- '//subname//' -- '//trim(testname)//' -- maxval')
         call test_is(maxindex, maxindex_t, modname//' -- '//subname//' -- '//trim(testname)//' -- maxindex')
      end subroutine check_results

   end subroutine test_which_max
!------------------------------------------------------------------------------

end module test_mkindexmapMod
   
