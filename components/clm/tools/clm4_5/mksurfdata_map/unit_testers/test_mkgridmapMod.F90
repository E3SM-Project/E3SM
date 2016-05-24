module test_mkgridmapMod
  ! Module for testing mkgridmapMod

  use mkgridmapMod
  use test_mod
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: test_gridmap_areastddev
  public :: test_gridmap_areaave_default
  public :: test_gridmap_areaave_srcmask
  public :: test_gridmap_areaave_srcmask2
  
  character(len=*), parameter :: modname = 'test_mkgridmapMod'

contains

  !------------------------------------------------------------------------------
  subroutine test_gridmap_areaave_default

    implicit none

    type(gridmap_type) :: gridmap
    character(len=128) :: testname

    real(r8), allocatable :: src_array(:)
    real(r8), allocatable :: dst_array(:)
    real(r8), allocatable :: dst_array_t(:)

    real(r8), parameter :: nodata = -1._r8
    real(r8), parameter :: eps = 1.e-13_r8

    character(len=*), parameter :: subname = 'test_gridmap_areaave_default'

    ! Note about the gridmaps for the tests here:
    ! For most tests here, the test arrays are: (1) simple case, (2) the main case to
    ! test, (3) simple case. Thus, the main case in question is #2 of 3, and we're always
    ! basically just testing one scenario in each call to the subroutine (rather than
    ! doing a bunch of tests at once, which could make setting up the test arrays more
    ! error-prone).

    ! Set up a gridmap with 0 weight of overlap on dest #2
    gridmap%na = 4
    gridmap%nb = 3
    gridmap%ns = 4
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4/)
    gridmap%dst_indx = (/1,1,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.25_r8,0.75_r8/)   ! weights of sources 3:4 on test 3
    gridmap%frac_dst = (/1.0, 0.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'no overlap'
    src_array = (/0.1_r8,0.2_r8,0.3_r8,0.4_r8/)
    dst_array_t = (/0.125_r8, nodata, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, dst_array, dst_array_t)

    ! Set up a gridmap with a single point overlapping dest #2
    gridmap%na = 5
    gridmap%nb = 3
    gridmap%ns = 5
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4,5/)
    gridmap%dst_indx = (/1,1,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         1.0_r8, &           ! weight of source 3 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 4:5 on test 3
    gridmap%frac_dst = (/1.0, 1.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'single overlap'
    src_array = (/0.1_r8,0.2_r8,0.5_r8,0.3_r8,0.4_r8/)
    dst_array_t = (/0.125_r8, 0.5_r8, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! Now change the overlap point to have weight=0
    testname = 'single overlap with 0 weight'
    gridmap%wovr(3) = 0.0_r8
    gridmap%frac_dst(2) = 0.0_r8
    dst_array_t(2) = nodata
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, dst_array, dst_array_t)

    ! Set up a gridmap for the remaining tests
    ! This gridmap will have 3 src cells, 9 dest cells, and:
    ! src 1: just overlaps with dst 1
    ! src 2: overlaps with dst 1 & dst 2
    ! src 3..7: just overlaps with dst 2
    ! src 8: overlaps with dst 2 & dst 3
    ! src 9: just overlaps with dst 3
    gridmap%na = 9
    gridmap%nb = 3
    gridmap%ns = 11
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
    gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.05_r8,0.05_r8,0.1_r8,0.3_r8,0.2_r8,0.15_r8,0.15_r8, &  ! weights of sources 2:8 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 8:9 on test 3
    gridmap%frac_dst = (/1.0_r8, 1.0_r8, 1.0_r8/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))

   
    testname='multiple overlaps, all the same value'
    src_array = (/0.1_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.6_r8/)
    dst_array_t = (/0.2_r8, 0.5_r8, 0.575_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! dividing the weights by 2 shouldn't affect the mean
    testname='weights divided by 2'
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    
    ! using frac_dst > 1 should be okay
    testname='frac_dst > 1'
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8

    deallocate(src_array, dst_array, dst_array_t)

  end subroutine test_gridmap_areaave_default

  !------------------------------------------------------------------------------
  subroutine test_gridmap_areaave_srcmask

    implicit none

    type(gridmap_type) :: gridmap
    character(len=128) :: testname

    real(r8), allocatable :: src_array(:)
    real(r8), allocatable :: mask_src(:)
    real(r8), allocatable :: dst_array(:)
    real(r8), allocatable :: dst_array_t(:)

    real(r8), parameter :: nodata = -1._r8
    real(r8), parameter :: eps = 1.e-13_r8

    character(len=*), parameter :: subname = 'test_gridmap_areaave_srcmask'

    ! Note about the gridmaps for the tests here:
    ! For most tests here, the test arrays are: (1) simple case, (2) the main case to
    ! test, (3) simple case. Thus, the main case in question is #2 of 3, and we're always
    ! basically just testing one scenario in each call to the subroutine (rather than
    ! doing a bunch of tests at once, which could make setting up the test arrays more
    ! error-prone).

    ! Set up a gridmap with 0 weight of overlap on dest #2
    gridmap%na = 4
    gridmap%nb = 3
    gridmap%ns = 4
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4/)
    gridmap%dst_indx = (/1,1,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.25_r8,0.75_r8/)   ! weights of sources 3:4 on test 3
    gridmap%frac_dst = (/1.0, 0.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'no overlap'
    src_array = (/0.1_r8,0.2_r8,0.3_r8,0.4_r8/)
    mask_src(:) = 1.0_r8
    dst_array_t = (/0.125_r8, nodata, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, mask_src, dst_array, dst_array_t)

    ! Set up a gridmap with a single point overlapping dest #2
    gridmap%na = 5
    gridmap%nb = 3
    gridmap%ns = 5
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4,5/)
    gridmap%dst_indx = (/1,1,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         1.0_r8, &           ! weight of source 3 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 4:5 on test 3
    gridmap%frac_dst = (/1.0, 1.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'single overlap'
    src_array = (/0.1_r8,0.2_r8,0.5_r8,0.3_r8,0.4_r8/)
    mask_src(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.5_r8, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! Now change the overlap point to have src_mask=0
    testname = 'single overlap with 0 src_mask'
    mask_src(3) = 0.0_r8
    dst_array_t(2) = nodata
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, mask_src, dst_array, dst_array_t)

    ! Set up a gridmap for the remaining tests
    ! This gridmap will have 3 src cells, 9 dest cells, and:
    ! src 1: just overlaps with dst 1
    ! src 2: overlaps with dst 1 & dst 2
    ! src 3..7: just overlaps with dst 2
    ! src 8: overlaps with dst 2 & dst 3
    ! src 9: just overlaps with dst 3
    gridmap%na = 9
    gridmap%nb = 3
    gridmap%ns = 11
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
    gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.05_r8,0.05_r8,0.1_r8,0.3_r8,0.2_r8,0.15_r8,0.15_r8, &  ! weights of sources 2:8 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 8:9 on test 3
    gridmap%frac_dst = (/1.0_r8, 1.0_r8, 1.0_r8/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))

   
    testname='multiple overlaps, all the same value'
    src_array = (/0.1_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.6_r8/)
    mask_src(:) = 1.0_r8
    dst_array_t = (/0.2_r8, 0.5_r8, 0.575_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values, srcmask'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = (/1.0_r8, 1.0_r8, 0.0_r8, 0.5_r8, 1.0_r8, 0.5_r8, 0.0_r8, 1.0_r8, 1.0_r8/)
    dst_array_t = (/0.125_r8, 0.923076923076923_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! dividing the weights by 2 and dividing mask_src by a constant shouldn't affect the mean
    testname='weights divided by 2'
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 0.25_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    
    ! using frac_dst > 1 should be okay
    testname='frac_dst > 1'
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 0.25_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8


    deallocate(src_array, mask_src, dst_array, dst_array_t)

  end subroutine test_gridmap_areaave_srcmask

  !------------------------------------------------------------------------------
  subroutine test_gridmap_areaave_srcmask2

    implicit none

    type(gridmap_type) :: gridmap
    character(len=128) :: testname

    real(r8), allocatable :: src_array(:)
    real(r8), allocatable :: mask_src(:)
    real(r8), allocatable :: dst_array(:)
    real(r8), allocatable :: mask_dst(:)
    real(r8), allocatable :: dst_array_t(:)

    real(r8), parameter :: mask_dst_min = 0.0_r8
    real(r8), parameter :: nodata = -1._r8
    real(r8), parameter :: eps = 1.e-13_r8

    character(len=*), parameter :: subname = 'test_gridmap_areaave_srcmask2'

    ! Note about the gridmaps for the tests here:
    ! For most tests here, the test arrays are: (1) simple case, (2) the main case to
    ! test, (3) simple case. Thus, the main case in question is #2 of 3, and we're always
    ! basically just testing one scenario in each call to the subroutine (rather than
    ! doing a bunch of tests at once, which could make setting up the test arrays more
    ! error-prone).

    ! Set up a gridmap with 0 weight of overlap on dest #2
    gridmap%na = 4
    gridmap%nb = 3
    gridmap%ns = 4
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4/)
    gridmap%dst_indx = (/1,1,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.25_r8,0.75_r8/)   ! weights of sources 3:4 on test 3
    gridmap%frac_dst = (/1.0, 0.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         mask_dst   (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'no overlap'
    src_array = (/0.1_r8,0.2_r8,0.3_r8,0.4_r8/)
    mask_src(:) = 1.0_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, nodata, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, mask_src, dst_array, mask_dst, dst_array_t)

    ! Set up a gridmap with a single point overlapping dest #2
    gridmap%na = 5
    gridmap%nb = 3
    gridmap%ns = 5
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4,5/)
    gridmap%dst_indx = (/1,1,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         1.0_r8, &           ! weight of source 3 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 4:5 on test 3
    gridmap%frac_dst = (/1.0, 1.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         mask_dst   (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'single overlap'
    src_array = (/0.1_r8,0.2_r8,0.5_r8,0.3_r8,0.4_r8/)
    mask_src(:) = 1.0_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.5_r8, 0.375_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! Now change the overlap point to have src_mask=0
    testname = 'single overlap with 0 src_mask'
    mask_src(3) = 0.0_r8
    mask_dst(:) = 1.0_r8
    dst_array_t(2) = nodata
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, mask_src, dst_array, mask_dst, dst_array_t)

    ! Set up a gridmap for the remaining tests
    ! This gridmap will have 3 src cells, 9 dest cells, and:
    ! src 1: just overlaps with dst 1
    ! src 2: overlaps with dst 1 & dst 2
    ! src 3..7: just overlaps with dst 2
    ! src 8: overlaps with dst 2 & dst 3
    ! src 9: just overlaps with dst 3
    gridmap%na = 9
    gridmap%nb = 3
    gridmap%ns = 11
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
    gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.05_r8,0.05_r8,0.1_r8,0.3_r8,0.2_r8,0.15_r8,0.15_r8, &  ! weights of sources 2:8 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 8:9 on test 3
    gridmap%frac_dst = (/1.0_r8, 1.0_r8, 1.0_r8/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         mask_src   (gridmap%na), &
         dst_array  (gridmap%nb), &
         mask_dst   (gridmap%nb), &
         dst_array_t(gridmap%nb))

   
    testname='multiple overlaps, all the same value'
    src_array = (/0.1_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.6_r8/)
    mask_src(:) = 1.0_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.2_r8, 0.5_r8, 0.575_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 1.0_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values, dst mask'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 1.0_r8
    mask_dst(:) = (/1.0_r8, 0.0_r8, 1.0_r8/)
    dst_array_t = (/0.125_r8, nodata, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values, srcmask'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = (/1.0_r8, 1.0_r8, 0.0_r8, 0.5_r8, 1.0_r8, 0.5_r8, 0.0_r8, 1.0_r8, 1.0_r8/)
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.923076923076923_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! dividing the weights by 2 and dividing mask_src by a constant shouldn't affect the mean
    testname='weights divided by 2'
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 0.25_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8

    ! using frac_dst > 1 should be okay
    testname='frac_dst > 1'
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    mask_src(:) = 0.25_r8
    mask_dst(:) = 1.0_r8
    dst_array_t = (/0.125_r8, 0.875_r8, 1.775_r8/)
    call gridmap_areaave(gridmap, src_array, dst_array, nodata, mask_src, mask_dst, mask_dst_min)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    
    deallocate(src_array, mask_src, dst_array, mask_dst, dst_array_t)

  end subroutine test_gridmap_areaave_srcmask2

  !------------------------------------------------------------------------------
  subroutine test_gridmap_areastddev

    implicit none

    type(gridmap_type) :: gridmap
    character(len=128) :: testname

    real(r8), allocatable :: src_array(:)
    real(r8), allocatable :: dst_array(:)
    real(r8), allocatable :: dst_array_t(:)

    real(r8), parameter :: nodata = -1._r8
    real(r8), parameter :: eps = 1.e-13_r8

    character(len=*), parameter :: subname = 'test_gridmap_areastddev'

    ! Note about the gridmaps for the tests here:
    ! For most tests here, the test arrays are: (1) simple case, (2) the main case to
    ! test, (3) simple case. Thus, the main case in question is #2 of 3, and we're always
    ! basically just testing one scenario in each call to the subroutine (rather than
    ! doing a bunch of tests at once, which could make setting up the test arrays more
    ! error-prone).

    ! Set up a gridmap with 0 weight of overlap on dest #2
    gridmap%na = 4
    gridmap%nb = 3
    gridmap%ns = 4
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4/)
    gridmap%dst_indx = (/1,1,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.25_r8,0.75_r8/)   ! weights of sources 3:4 on test 3
    gridmap%frac_dst = (/1.0, 0.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'no overlap'
    src_array = (/0.1_r8,0.2_r8,0.3_r8,0.4_r8/)
    dst_array_t = (/0.04330127018922193_r8, nodata, 0.04330127018922195_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, dst_array, dst_array_t)

    ! Set up a gridmap with a single point overlapping dest #2
    gridmap%na = 5
    gridmap%nb = 3
    gridmap%ns = 5
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,3,4,5/)
    gridmap%dst_indx = (/1,1,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         1.0_r8, &           ! weight of source 3 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 4:5 on test 3
    gridmap%frac_dst = (/1.0, 1.0, 1.0/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))
    testname = 'single overlap'
    src_array = (/0.1_r8,0.2_r8,0.5_r8,0.3_r8,0.4_r8/)
    dst_array_t = (/0.04330127018922193_r8, 0.0, 0.04330127018922195_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    deallocate(gridmap%src_indx, gridmap%dst_indx, gridmap%wovr, gridmap%frac_dst)
    deallocate(src_array, dst_array, dst_array_t)

    ! Set up a gridmap for the remaining tests
    ! This gridmap will have 3 src cells, 9 dest cells, and:
    ! src 1: just overlaps with dst 1
    ! src 2: overlaps with dst 1 & dst 2
    ! src 3..7: just overlaps with dst 2
    ! src 8: overlaps with dst 2 & dst 3
    ! src 9: just overlaps with dst 3
    gridmap%na = 9
    gridmap%nb = 3
    gridmap%ns = 11
    allocate(gridmap%src_indx(gridmap%ns), &
             gridmap%dst_indx(gridmap%ns), &
             gridmap%wovr    (gridmap%ns), &
             gridmap%frac_dst(gridmap%nb))
    gridmap%src_indx = (/1,2,2,3,4,5,6,7,8,8,9/)
    gridmap%dst_indx = (/1,1,2,2,2,2,2,2,2,3,3/)
    gridmap%wovr     = (/0.75_r8,0.25_r8, &  ! weights of sources 1:2 on dest 1
                         0.05_r8,0.05_r8,0.1_r8,0.3_r8,0.2_r8,0.15_r8,0.15_r8, &  ! weights of sources 2:8 on dest 2
                         0.25_r8,0.75_r8/)   ! weights of sources 8:9 on test 3
    gridmap%frac_dst = (/1.0_r8, 1.0_r8, 1.0_r8/)
    gridmap%set = 'gridmap_IsSet'
    allocate(src_array  (gridmap%na), &
         dst_array  (gridmap%nb), &
         dst_array_t(gridmap%nb))

   
    testname='multiple overlaps, all the same value'
    src_array = (/0.1_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.6_r8/)
    dst_array_t = (/0.1732050807568877_r8, 0.0, 0.04330127018922193_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    testname='multiple overlaps, different values'
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.04330127018922193_r8, 0.5346727971385864_r8, 0.04330127018922197_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))

    ! dividing the weights by 2 shouldn't affect the standard deviation
    testname='weights divided by 2'
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.04330127018922193_r8, 0.5346727971385864_r8, 0.04330127018922197_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8

    ! using frac_dst > 1 should be okay
    testname='frac_dst > 1'
    gridmap%wovr(:) = gridmap%wovr(:) * 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) * 2.0_r8
    src_array = (/0.1_r8, 0.2_r8, 0.3_r8, 0.7_r8, 0.5_r8, 1.5_r8, 0.5_r8, 1.7_r8, 1.8_r8/)
    dst_array_t = (/0.04330127018922193_r8, 0.5346727971385864_r8, 0.04330127018922197_r8/)
    call gridmap_areastddev(gridmap, src_array, dst_array, nodata)
    call test_close(dst_array, dst_array_t, eps, modname//' -- '//subname//' -- '//trim(testname))
    ! restore wovr & frac_dst
    gridmap%wovr(:) = gridmap%wovr(:) / 2.0_r8
    gridmap%frac_dst(:) = gridmap%frac_dst(:) / 2.0_r8
    
    deallocate(src_array, dst_array, dst_array_t)

  end subroutine test_gridmap_areastddev
end module test_mkgridmapMod
