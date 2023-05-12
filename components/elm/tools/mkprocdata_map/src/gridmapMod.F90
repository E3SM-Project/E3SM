module gridmapMod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use fileutils,    only : getfil

  implicit none
  private
  include 'netcdf.inc'

  type gridmap_type
     character(len=32) :: name
     integer :: na      ! size of source domain
     integer :: nb      ! size of destination domain
     integer :: ni      ! number of row in the matrix
     integer :: nj      ! number of col in the matrix
     integer :: ns      ! number of non-zero elements in matrix
     real(r8), pointer :: yc_src(:) 	! "degrees" 
     real(r8), pointer :: yc_dst(:) 	! "degrees" 
     real(r8), pointer :: xc_src(:) 	! "degrees" 
     real(r8), pointer :: xc_dst(:) 	! "degrees" 
     integer , pointer :: mask_src(:) 	! "unitless" 
     integer , pointer :: mask_dst(:) 	! "unitless" 
     real(R8), pointer :: area_src(:)   ! area of a grid in map (radians)
     real(R8), pointer :: area_dst(:)   ! area of b grid in map (radians)
     real(r8), pointer :: frac_src(:) 	! "unitless" 
     real(r8), pointer :: frac_dst(:) 	! "unitless" 
     integer , pointer :: src_indx(:)   ! correpsonding column index
     integer , pointer :: dst_indx(:)   ! correpsonding row    index
     real(r8), pointer :: wovr(:)       ! wt of overlap input cell
     real(r8), pointer :: scalepft_i(:) ! PFT wt of overlap input cell
  end type gridmap_type
  public :: gridmap_type

  public :: gridmap_mapread
  public :: gridmap_areaave
  public :: gridmap_clean

  interface gridmap_areaave
     module procedure gridmap_areaave_default
     module procedure gridmap_areaave_mask
  end interface

  ! questions - how does the reverse mapping occur 
  ! is mask_dst read in - and what happens if this is very different
  ! from frac_dst which is calculated by mapping frac_src?
  ! in frac - isn't grid1_frac always 1 or 0?
  
contains

  subroutine gridmap_mapread(gridmap, fileName)

    !--- input/output parameters ---
    type(gridmap_type), intent(out) :: gridmap   ! mapping data
    character(len=*)   , intent(in)  :: filename  ! netCDF file to read


    !--- local ---
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
    character(len=256)    :: locfn

    !--- formats ---
    character(*),parameter :: subName = '(gridmap_map_read) '
    character(*),parameter :: F00 = '("(gridmap_map_read) ",4a)'
    character(*),parameter :: F01 = '("(gridmap_map_read) ",2(a,i7))'

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    write(6,F00) "reading mapping matrix data..."

    ! open & read the file
    write(6,F00) "* file name                  : ",trim(fileName)

    call getfil (trim(filename), locfn, 0)
    rcode = nf_open(locfn ,NF_NOWRITE, fid)
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
             gridmap%xc_src(na),   &
             gridmap%yc_src(na), stat=rcode)
    if (rcode /= 0) then
       write(6,*) SubName//' ERROR: allocate gridmap'
       stop 1
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

    rcode = nf_inq_varid(fid,'frac_b',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%frac_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'mask_a',vid)
    rcode = nf_get_var_int(fid, vid, gridmap%mask_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'mask_b',vid)
    rcode = nf_get_var_int(fid, vid, gridmap%mask_dst)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'xc_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%xc_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_inq_varid(fid,'yc_a',vid)
    rcode = nf_get_var_double(fid, vid, gridmap%yc_src)
    if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

    rcode = nf_close(fid)

  end subroutine gridmap_mapread

!==========================================================================

  subroutine gridmap_areaave_default (gridmap, src_array, dst_array)

    !--- input/output parameters ---
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    
    !--- local ---
    integer :: n,ns,ni,no
    real(r8):: wt,frac
    !----------------------------------------------------------

    dst_array = 0._r8
    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       frac = gridmap%frac_dst(no)
       if (frac > 0.) then  
          dst_array(no) = dst_array(no) + wt * src_array(ni)/frac
       end if
    end do

  end subroutine gridmap_areaave_default

!==========================================================================

  subroutine gridmap_areaave_mask (gridmap, src_array, dst_array, src_mask, &
       spval)

    !--- input/output parameters ---
    type(gridmap_type) , intent(in) :: gridmap   ! gridmap data
    real(r8), intent(in) :: src_array(:)
    real(r8), intent(out):: dst_array(:)
    real(r8), intent(in) :: src_mask(:)
    real(r8), intent(in), optional :: spval
    
    !--- local ---
    integer :: n,ns,ni,no
    real(r8):: wt
    real(r8), allocatable :: wtnorm(:)
    !----------------------------------------------------------

    ns = size(dst_array)
    allocate(wtnorm(ns)) 

    wtnorm(:) = 0._r8
    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (src_mask(ni) > 0) then
          wtnorm(no) = wtnorm(no) + wt*src_mask(ni)
       end if
    end do

    do no = 1,ns
       if (wtnorm(no) > 0) then
          dst_array(no) = 0._r8
       else
          if (present(spval)) then
             dst_array(no) = spval
          else
             dst_array(no) = 1.e36
          end if
       end if
    end do

    do n = 1,gridmap%ns
       ni = gridmap%src_indx(n)
       no = gridmap%dst_indx(n)
       wt = gridmap%wovr(n)
       if (wtnorm(no) > 0) then 
          dst_array(no) = dst_array(no) + wt*src_mask(ni)*src_array(ni)/wtnorm(no)
       end if
    end do

    deallocate(wtnorm)

  end subroutine gridmap_areaave_mask

!==========================================================================

  subroutine gridmap_clean(gridmap)

    !--- input/output parameters ---
    implicit none
    type(gridmap_type), intent(inout)       :: gridmap


    !--- local ---
    character(len=*), parameter :: subName = "gridmap_clean"
    integer ier    ! error flag
    !----------------------------------------------------------

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
       stop 1
    endif

  end subroutine gridmap_clean

end module gridmapMod


