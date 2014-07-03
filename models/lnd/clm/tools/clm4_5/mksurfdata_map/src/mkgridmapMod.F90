module mkgridmapMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkgridmapMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

! !PUBLIC TYPES:
  type gridmap_type
     character(len=32) :: set ! If set or not
     character(len=32) :: name
     integer :: na            ! size of source domain
     integer :: nb            ! size of destination domain
     integer :: ni            ! number of row in the matrix
     integer :: nj            ! number of col in the matrix
     integer :: ns            ! number of non-zero elements in matrix
     real(r8), pointer :: yc_src(:)     ! "degrees" 
     real(r8), pointer :: yc_dst(:)     ! "degrees" 
     real(r8), pointer :: xc_src(:)     ! "degrees" 
     real(r8), pointer :: xc_dst(:)     ! "degrees" 
     integer , pointer :: mask_src(:)   ! "unitless" 
     integer , pointer :: mask_dst(:)   ! "unitless" 
     real(R8), pointer :: area_src(:)   ! area of a grid in map (radians)
     real(R8), pointer :: area_dst(:)   ! area of b grid in map (radians)
     real(r8), pointer :: frac_src(:)   ! "unitless" 
     real(r8), pointer :: frac_dst(:)   ! "unitless" 
     integer , pointer :: src_indx(:)   ! correpsonding column index
     integer , pointer :: dst_indx(:)   ! correpsonding row    index
     real(r8), pointer :: wovr(:)       ! wt of overlap input cell
  end type gridmap_type
  public :: gridmap_type
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: gridmap_setptrs     ! Set pointers to gridmap data
  public :: gridmap_mapread     ! Read in gridmap
  public :: gridmap_check       ! Check validity of a gridmap
  public :: gridmap_areaave     ! do area average
  public :: gridmap_areastddev  ! do area-weighted standard deviation
  public :: gridmap_clean       ! Clean and deallocate a gridmap structure
!
!
! !REVISION HISTORY:
! Author Mariana Vertenstein

  interface gridmap_areaave
     module procedure gridmap_areaave_default
     module procedure gridmap_areaave_srcmask
     module procedure gridmap_areaave_srcmask2
  end interface

  ! questions - how does the reverse mapping occur 
  ! is mask_dst read in - and what happens if this is very different
  ! from frac_dst which is calculated by mapping frac_src?
  ! in frac - isn't grid1_frac always 1 or 0?

! !PRIVATE MEMBER FUNCTIONS:
  private :: gridmap_checkifset

  character(len=32), parameter :: isSet = "gridmap_IsSet"
  
!
!EOP
!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_setptrs
!
! !INTERFACE:
  subroutine gridmap_setptrs(gridmap, nsrc, ndst, ns, yc_src, yc_dst, &
                             xc_src, xc_dst, mask_src, mask_dst,      &
                             frac_src, frac_dst, src_indx, dst_indx )
!
! !DESCRIPTION:
! This subroutine assigns pointers to some of the map type data.
!
! !ARGUMENTS:
     implicit none
     type(gridmap_type), intent(in) :: gridmap   ! mapping data
     integer, optional :: nsrc                   ! size of source domain
     integer, optional :: ndst                   ! size of destination domain
     integer, optional :: ns                     ! number of non-zero elements in matrix
     integer,  optional, pointer :: dst_indx(:)  ! Destination index
     integer,  optional, pointer :: src_indx(:)  ! Destination index
     real(r8), optional, pointer :: yc_src(:)    ! "degrees" 
     real(r8), optional, pointer :: yc_dst(:)    ! "degrees" 
     real(r8), optional, pointer :: xc_src(:)    ! "degrees" 
     real(r8), optional, pointer :: xc_dst(:)    ! "degrees" 
     integer , optional, pointer :: mask_src(:)  ! "unitless" 
     integer , optional, pointer :: mask_dst(:)  ! "unitless" 
     real(r8), optional, pointer :: frac_src(:)  ! "unitless" 
     real(r8), optional, pointer :: frac_dst(:)  ! "unitless" 
!
! !REVISION HISTORY:
!   Created by Erik Kluzek
!
! !LOCAL VARIABLES:
!EOP
!------------------------------------------------------------------------------
     character(*),parameter :: subName = '(gridmap_setptrs) '

     call gridmap_checkifset( gridmap, subname )
     if ( present(nsrc)     ) nsrc     = gridmap%na
     if ( present(ndst)     ) ndst     = gridmap%nb
     if ( present(ns)       ) ns       = gridmap%ns
     if ( present(yc_src)   ) yc_src   => gridmap%yc_src
     if ( present(xc_src)   ) xc_src   => gridmap%xc_src
     if ( present(mask_src) ) mask_src => gridmap%mask_src
     if ( present(frac_src) ) frac_src => gridmap%frac_src
     if ( present(yc_dst)   ) yc_dst   => gridmap%yc_dst
     if ( present(xc_dst)   ) xc_dst   => gridmap%xc_dst
     if ( present(mask_dst) ) mask_dst => gridmap%mask_dst
     if ( present(frac_dst) ) frac_dst => gridmap%frac_dst
     if ( present(dst_indx) ) dst_indx => gridmap%dst_indx
     if ( present(src_indx) ) src_indx => gridmap%src_indx
  end subroutine gridmap_setptrs

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_mapread
!
! !INTERFACE:
  subroutine gridmap_mapread(gridmap, fileName)
!
! !DESCRIPTION:
! This subroutine reads in the map file
!
! !USES:
    use mkutilsMod, only : convert_latlon
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(gridmap_type), intent(out) :: gridmap   ! mapping data
    character(len=*)   , intent(in)  :: filename  ! netCDF file to read
!
! !REVISION HISTORY:
!   Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
    integer :: n       ! generic loop indicies
    integer :: na      ! size of source domain
    integer :: nb      ! size of destination domain
    integer :: igrow   ! aVect index for matrix row
    integer :: igcol   ! aVect index for matrix column
    integer :: iwgt    ! aVect index for matrix element
    integer :: iarea   ! aVect index for area


    character,allocatable :: str(:)  ! variable length char string
    character(len=256)    :: attstr  ! netCDF attribute name string
    integer               :: rcode   ! netCDF routine return code
    integer               :: fid     ! netCDF file      ID
    integer               :: vid     ! netCDF variable  ID
    integer               :: did     ! netCDF dimension ID
    integer               :: ns      ! size of array

    real(r8), parameter   :: tol = 1.0e-5_r8  ! tolerance for checking that mapping data
                                              ! are within expected bounds

    !--- formats ---
    character(*),parameter :: subName = '(gridmap_map_read) '
    character(*),parameter :: F00 = '("(gridmap_map_read) ",4a)'
    character(*),parameter :: F01 = '("(gridmap_map_read) ",2(a,i7))'
!EOP
!------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    write(6,F00) "reading mapping matrix data..."

    ! open & read the file
    write(6,F00) "* file name                  : ",trim(fileName)

    rcode = nf_open(filename ,NF_NOWRITE, fid)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    !--- allocate memory & get matrix data ----------
    rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
    rcode = nf_inq_dimlen(fid, did  , gridmap%ns)
    rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
    rcode = nf_inq_dimlen(fid, did  , gridmap%na)
    rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
    rcode = nf_inq_dimlen(fid, did  , gridmap%nb)

    write(6,*) "* matrix dimensions rows x cols :",gridmap%na,' x',gridmap%nb
    write(6,*) "* number of non-zero elements: ",gridmap%ns

    ns = gridmap%ns
    na = gridmap%na
    nb = gridmap%nb
    allocate(gridmap%wovr(ns)    , &
             gridmap%src_indx(ns), &
             gridmap%dst_indx(ns), &
             gridmap%mask_src(na), &
             gridmap%area_src(na), &
             gridmap%frac_src(na), &
             gridmap%area_dst(nb), &
             gridmap%frac_dst(nb), &
             gridmap%mask_dst(nb), &
             gridmap%xc_dst(nb),   &
             gridmap%yc_dst(nb),   &
             gridmap%xc_src(na),   &
             gridmap%yc_src(na), stat=rcode)
    if (rcode /= 0) then
       write(6,*) SubName//' ERROR: allocate gridmap'
       call abort()
    endif

    rcode = nf_inq_varid(fid,'S'  ,vid)
    rcode = nf_get_var_double(fid,vid  ,gridmap%wovr)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'row',vid)
    rcode = nf_get_var_int(fid, vid  ,gridmap%dst_indx)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'col',vid)
    rcode = nf_get_var_int(fid, vid, gridmap%src_indx)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'area_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%area_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'area_b',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%area_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'frac_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%frac_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    if ( any(gridmap%frac_src(:) < 0.0_r8 .or. gridmap%frac_src > (1.0_r8 + tol)) )then
       write(6,*) SubName//' ERROR: frac_src out of bounds'
       write(6,*) 'max = ', maxval(gridmap%frac_src), ' min = ', minval(gridmap%frac_src)
       call abort()
    end if

    rcode = nf_inq_varid(fid,'frac_b',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%frac_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    if ( any(gridmap%frac_dst(:) < 0.0_r8 .or. gridmap%frac_dst > (1.0_r8 + tol)) )then
       write(6,*) SubName//' ERROR: frac_dst out of bounds'
       write(6,*) 'max = ', maxval(gridmap%frac_dst), ' min = ', minval(gridmap%frac_dst)
       call abort()
    end if

    rcode = nf_inq_varid(fid,'mask_a',vid)
    rcode = nf_get_var_int(fid, vid, gridmap%mask_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    if ( any(gridmap%mask_src(:) < 0 .or. gridmap%mask_src > 1) )then
       write(6,*) SubName//' ERROR: mask_src out of bounds'
       call abort()
    end if

    rcode = nf_inq_varid(fid,'mask_b',vid)
    rcode = nf_get_var_int(fid, vid, gridmap%mask_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    if ( any(gridmap%mask_dst(:) < 0 .or. gridmap%mask_dst > 1) )then
       write(6,*) SubName//' ERROR: mask_dst out of bounds'
       call abort()
    end if

    rcode = nf_inq_varid(fid,'xc_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%xc_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    call convert_latlon(fid, 'xc_a', gridmap%xc_src)

    rcode = nf_inq_varid(fid,'yc_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%yc_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    call convert_latlon(fid, 'yc_a', gridmap%yc_src)

    rcode = nf_inq_varid(fid,'xc_b',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%xc_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    call convert_latlon(fid, 'xc_b', gridmap%xc_dst)

    rcode = nf_inq_varid(fid,'yc_b',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%yc_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
    call convert_latlon(fid, 'yc_b', gridmap%yc_dst)

    rcode = nf_close(fid)

    gridmap%set = IsSet

  end subroutine gridmap_mapread

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_check
!
! !INTERFACE:
  subroutine gridmap_check(gridmap, caller)
!
! !DESCRIPTION:
! Check validity of a gridmap
! Aborts if there are any errors
!
! !USES:
    use mkvarctl, only : mksrf_gridtype
    use mkvarpar, only : re
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type) , intent(in) :: gridmap   ! mapping data
    character(len=*)   , intent(in) :: caller    ! calling subroutine (used for error messages)
!
! !REVISION HISTORY:
!   Created by Bill Sacks
!
! !LOCAL VARIABLES:
    real(r8) :: sum_area_i        ! global sum of input area
    real(r8) :: sum_area_o        ! global sum of output area
    integer  :: ni,no,ns_i,ns_o   ! indices

    real(r8), parameter :: relerr = 0.00001        ! max error: sum overlap wts ne 1
    character(len=*), parameter :: subname = 'gridmap_check'
!EOP
!------------------------------------------------------------------------------

    ns_i = gridmap%na
    ns_o = gridmap%nb

    ! -----------------------------------------------------------------
    ! Error check prep
    ! Global sum of output area -- must multiply by fraction of
    ! output grid that is land as determined by input grid
    ! -----------------------------------------------------------------
    
    sum_area_i = 0.0_r8
    do ni = 1,ns_i
       sum_area_i = sum_area_i + gridmap%area_src(ni)*gridmap%frac_src(ni)*re**2
    enddo

    sum_area_o = 0.
    do no = 1,ns_o
       sum_area_o = sum_area_o + gridmap%area_dst(no)*gridmap%frac_dst(no)*re**2
    end do

    ! -----------------------------------------------------------------
    ! Error check1
    ! Compare global sum_area_i to global sum_area_o.
    ! -----------------------------------------------------------------

    if ( trim(mksrf_gridtype) == 'global' ) then
       if ( abs(sum_area_o/sum_area_i-1.) > relerr ) then
          write (6,*) subname//' ERROR from '//trim(caller)//': mapping areas not conserved'
          write (6,'(a30,e20.10)') 'global sum output field = ',sum_area_o
          write (6,'(a30,e20.10)') 'global sum input  field = ',sum_area_i
          stop
       end if
    end if

  end subroutine gridmap_check
    

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_areaave_default
!
! !INTERFACE:
  subroutine gridmap_areaave_default (gridmap, src_array, dst_array, nodata)
!
! !DESCRIPTION:
! This subroutine does a simple area average
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    real(r8), intent(in) :: nodata               ! value to apply where there are no input data
!
! !REVISION HISTORY:
!   Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
    integer :: n,ns,ni,no
    real(r8):: wt,frac
    real(r8), allocatable :: sum_weights(:)      ! sum of weights on the output grid
    character(*),parameter :: subName = '(gridmap_areaave_default) '
!EOP
!------------------------------------------------------------------------------
    call gridmap_checkifset( gridmap, subname )
    allocate(sum_weights(size(dst_array)))
    sum_weights = 0._r8
    dst_array = 0._r8

    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       frac = gridmap%frac_dst(no)
       if (frac > 0.) then  
          dst_array(no) = dst_array(no) + wt * src_array(ni)/frac
          sum_weights(no) = sum_weights(no) + wt
       end if
    end do

    where (sum_weights == 0._r8)
       dst_array = nodata
    end where

    deallocate(sum_weights)

  end subroutine gridmap_areaave_default

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_areaave_srcmask
!
! !INTERFACE:
  subroutine gridmap_areaave_srcmask (gridmap, src_array, dst_array, nodata, mask_src)
!
! !DESCRIPTION:
! This subroutine does an area average with the source mask
!
! !ARGUMENTS:
     implicit none
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    real(r8), intent(in) :: nodata               ! value to apply where there are no input data
    real(r8), intent(in) :: mask_src(:)
!
! !REVISION HISTORY:
!   Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
    integer :: n,ns,ni,no
    real(r8):: wt
    real(r8), allocatable :: wtnorm(:)
    character(*),parameter :: subName = '(gridmap_areaave_srcmask) '
!EOP
!------------------------------------------------------------------------------
    call gridmap_checkifset( gridmap, subname )
    ns = size(dst_array)
    allocate(wtnorm(ns)) 
    wtnorm(:) = 0._r8

    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (mask_src(ni) > 0) then
          wtnorm(no) = wtnorm(no) + wt*mask_src(ni)
       end if
    end do

    dst_array = 0._r8
    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (mask_src(ni) > 0) then 
          dst_array(no) = dst_array(no) + wt*mask_src(ni)*src_array(ni)/wtnorm(no)
       end if
    end do

    where (wtnorm == 0._r8)
       dst_array = nodata
    end where

    deallocate(wtnorm)

  end subroutine gridmap_areaave_srcmask

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_areaave_srcmask2
!
! !INTERFACE:
  subroutine gridmap_areaave_srcmask2 (gridmap, src_array, dst_array, nodata, mask_src, &
       mask_dst, mask_dst_min)
!
! !DESCRIPTION:
! This subroutine does an area average with the source mask and making sure the
! destination mask is valid as well.
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    real(r8), intent(in) :: nodata               ! value to apply where there are no input data
    real(r8), intent(in) :: mask_src(:)
    real(r8), intent(in) :: mask_dst(:)
    real(r8), intent(in) :: mask_dst_min
!
! !REVISION HISTORY:
!   Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
    integer :: n,ns,ni,no
    real(r8):: wt
    real(r8), allocatable :: wtnorm(:)
    character(*),parameter :: subName = '(gridmap_areaave_srcmask2) '
!EOP
!------------------------------------------------------------------------------

    call gridmap_checkifset( gridmap, subname )
    ns = size(dst_array)
    allocate(wtnorm(ns)) 
    wtnorm(:) = 0._r8

    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (mask_src(ni) > 0) then
          wtnorm(no) = wtnorm(no) + wt*mask_src(ni)
       end if
    end do

    dst_array = 0._r8
    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (mask_dst(no) > mask_dst_min) then
          if (mask_src(ni) > 0) then 
             dst_array(no) = dst_array(no) + wt*mask_src(ni)*src_array(ni)/wtnorm(no)
          end if
       end if
    end do

    where ((wtnorm == 0._r8) .or. (mask_dst <= mask_dst_min))
       dst_array = nodata
    end where

    deallocate(wtnorm)

  end subroutine gridmap_areaave_srcmask2

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_areastddev
!
! !INTERFACE:
  subroutine gridmap_areastddev (gridmap, src_array, dst_array, nodata)
!
! !DESCRIPTION:
! Computes area-weighted standard deviation
!
! We use the definition of standard deviation that applies if you measure the full
! population (as opposed to the unbiased standard deviation that should be used when
! sampling a subset of the full population). (This is equivalent to using 1/N rather than
! 1/(N-1).) This makes sense if we assume that the underlying values are constant
! throughout each source grid cell -- in that case, we know the full population as long as
! we know the values in all source grid cells, which is generally the case.
!
! The formula is from <http://en.wikipedia.org/wiki/Weighted_mean#Weighted_sample_variance>
! (accessed 3-4-13). 
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    real(r8), intent(in) :: nodata               ! value to apply where there are no input data
!
! !REVISION HISTORY:
!   Created by Bill Sacks
!
! !LOCAL VARIABLES:
    integer :: n,ni,no
    integer :: ns_o                                ! number of output points
    real(r8):: wt                                  ! weight of overlap
    real(r8), allocatable :: weighted_means(:)     ! weighted mean on the output grid
    real(r8), allocatable :: sum_weights(:)        ! sum of weights on the output grid
    character(*),parameter :: subName = '(gridmap_areastddev) '
!EOP
!------------------------------------------------------------------------------
    call gridmap_checkifset( gridmap, subname )

    ns_o = size(dst_array)
    allocate(weighted_means(ns_o))
    call gridmap_areaave(gridmap, src_array, weighted_means, nodata=0._r8)

    ! WJS (3-5-13): I believe that sum_weights should be the same as gridmap%frac_dst,
    ! but I'm not positive of this, so we compute it explicitly to be safe
    allocate(sum_weights(ns_o))
    sum_weights(:) = 0._r8
    dst_array(:)   = 0._r8
    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       ! The following accumulates the numerator of the weighted sigma-squared
       dst_array(no) = dst_array(no) + wt * (src_array(ni) - weighted_means(no))**2
       sum_weights(no) = sum_weights(no) + wt
    end do

    do no = 1,ns_o
       if (sum_weights(no) > 0._r8) then
          dst_array(no) = sqrt(dst_array(no)/sum_weights(no))
       else
          dst_array(no) = nodata
       end if
    end do

    deallocate(weighted_means, sum_weights)

  end subroutine gridmap_areastddev

!==========================================================================

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridmap_clean
!
! !INTERFACE:
  subroutine gridmap_clean(gridmap)
!
! !DESCRIPTION:
! This subroutine deallocates the gridmap type
!
! !ARGUMENTS:
    implicit none
    type(gridmap_type), intent(inout)       :: gridmap
!
! !REVISION HISTORY:
!   Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
    character(len=*), parameter :: subName = "gridmap_clean"
    integer ier    ! error flag
!EOP
!------------------------------------------------------------------------------
    if ( gridmap%set .eq. IsSet )then
       deallocate(gridmap%wovr    , &
                  gridmap%src_indx, &
                  gridmap%dst_indx, &
                  gridmap%mask_src, &
                  gridmap%mask_dst, &
                  gridmap%area_src, &
                  gridmap%area_dst, &
                  gridmap%frac_src, &
                  gridmap%frac_dst, &
                  gridmap%xc_src,   &
                  gridmap%yc_src, stat=ier)
       if (ier /= 0) then
          write(6,*) SubName//' ERROR: deallocate gridmap'
          call abort()
       endif
    else
       write(6,*) SubName//' Warning: calling '//trim(subName)//' on unallocated gridmap'
    end if
    gridmap%set = "NOT-set"

  end subroutine gridmap_clean

!==========================================================================

  subroutine gridmap_checkifset( gridmap, subname )

    implicit none
    type(gridmap_type), intent(in) :: gridmap
    character(len=*),   intent(in) :: subname

    if ( gridmap%set .ne. IsSet )then
       write(6,*) SubName//' ERROR: gridmap NOT set yet, run gridmap_mapread first'
       call abort()
    end if
  end subroutine gridmap_checkifset

end module mkgridmapMod


