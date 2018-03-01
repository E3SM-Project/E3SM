module mkdomainMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domain1Mod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkvarpar    , only : re
  use nanMod      , only : nan, bigint
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  type domain_type
     character*16     :: set        ! flag to check if domain is set
     integer          :: ns         ! global size of domain
     integer          :: ni,nj      ! for 2d domains only
     real(r8)         :: edgen      ! lsmedge north
     real(r8)         :: edgee      ! lsmedge east
     real(r8)         :: edges      ! lsmedge south
     real(r8)         :: edgew      ! lsmedge west
     integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
     real(r8),pointer :: frac(:)    ! fractional land
     real(r8),pointer :: latc(:)    ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)    ! longitude of grid cell (deg)
     real(r8),pointer :: lats(:)    ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:)    ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:)    ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:)    ! grid cell longitude, E edge (deg)
     real(r8),pointer :: area(:)    ! grid cell area (km**2) (only used for output grid)
     logical          :: is_2d      ! if this is a 2-d domain
     logical          :: fracset    ! if frac is set
     logical          :: maskset    ! if mask is set
  end type domain_type

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_clean         
  public domain_check         
  public domain_read 
  public domain_read_dims  ! get dimensions from a domain file (only public for unit testing)
  public domain_read_map
  public domain_write         
  public domain_checksame
!
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
  character*16,parameter :: set   = 'domain_set      '
  character*16,parameter :: unset = 'NOdomain_unsetNO'

  real(r8) :: flandmin = 0.001            !minimum land frac for land cell
!
! !PRIVATE MEMBER FUNCTIONS:
  private domain_init          
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_init
!
! !INTERFACE:
  subroutine domain_init(domain,ns)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
    integer            :: ns            ! grid size, 2d
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
    integer nb,ne
!
!------------------------------------------------------------------------------

    nb = 1
    ne = ns
    
    if (domain%set == set) then
       call domain_clean(domain)
    endif
    
    allocate(domain%mask(ns), &
             domain%frac(ns), &
             domain%latc(ns), &
             domain%lonc(ns), &
             domain%lats(ns), &
             domain%latn(ns), &
             domain%lonw(ns), &
             domain%lone(ns), &
             domain%area(ns), stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate mask, frac, lat, lon, area '
    endif
    
    domain%ns       = ns
    domain%mask     = -9999
    domain%frac     = -1.0e36
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan
    domain%set      = set
    domain%fracset  = .false.
    domain%maskset  = .false.
    
  end subroutine domain_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_clean
!
! !INTERFACE:
  subroutine domain_clean(domain)
!
! !DESCRIPTION:
! This subroutine deallocates the domain type
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
!
!------------------------------------------------------------------------------

    if (domain%set == set) then
       write(6,*) 'domain_clean: cleaning ',domain%ns
       deallocate(domain%mask, &
                  domain%frac, &
                  domain%latc, &
                  domain%lonc, &
                  domain%lats, &
                  domain%latn, &
                  domain%lonw, &
                  domain%lone, &
                  domain%area, stat=ier)
       if (ier /= 0) then
          write(6,*) 'domain_clean ERROR: deallocate mask, frac, lat, lon, area '
          call abort()
       endif
    else
       write(6,*) 'domain_clean WARN: clean domain unecessary '
    endif

    domain%ns       = bigint
    domain%set      = unset
    domain%fracset  = .false.
    domain%maskset  = .false.

end subroutine domain_clean

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_check
!
! !INTERFACE:
  subroutine domain_check(domain)
!
! !DESCRIPTION:
! This subroutine write domain info
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

    write(6,*) '  domain_check set       = ',trim(domain%set)
    write(6,*) '  domain_check ns        = ',domain%ns
    write(6,*) '  domain_check lonc      = ',minval(domain%lonc),maxval(domain%lonc)
    write(6,*) '  domain_check latc      = ',minval(domain%latc),maxval(domain%latc)
    write(6,*) '  domain_check mask      = ',minval(domain%mask),maxval(domain%mask)
    write(6,*) '  domain_check frac      = ',minval(domain%frac),maxval(domain%frac)
    write(6,*) ' '

end subroutine domain_check

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_read_map
!
! !INTERFACE:
  logical function domain_read_map(domain, fname)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
    use mkutilsMod, only : convert_latlon
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)     :: fname    ! this assumes a SCRIP mapping file - look at destination grid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    include 'netcdf.inc'
    integer :: i,j,n                           ! indices
    integer :: grid_rank                       ! rank of domain grid 
    integer :: ns                              ! size of domain grid
    integer :: ncid                            ! netCDF file id
    integer :: dimid                           ! netCDF dimension id
    integer :: varid                           ! netCDF variable id
    integer :: ndims                           ! number of dims for variable
    integer :: ier                             ! error status
    real(r8), allocatable :: xv(:,:)           ! local array for corner lons
    real(r8), allocatable :: yv(:,:)           ! local array for corner lats
    integer :: grid_dims(2)
    character(len= 32) :: subname = 'domain_read'
!-----------------------------------------------------------------

    domain_read_map = .true.

    ! Read domain file and compute stuff as needed

    call check_ret(nf_open(fname, 0, ncid), subname)

    ! Assume unstructured grid

    domain%ni = -9999
    domain%nj = -9999
    domain%is_2d = .false.

    ier = nf_inq_dimid  (ncid, 'n_b', dimid)
    if ( ier /= NF_NOERR )then
       domain_read_map = .false.
    else
      call check_ret(nf_inq_dimlen (ncid, dimid, domain%ns), subname)
    
      call check_ret(nf_inq_dimid  (ncid, 'dst_grid_rank', dimid), subname)
      call check_ret(nf_inq_dimlen (ncid, dimid, grid_rank), subname)

      if (grid_rank == 2) then
         call check_ret(nf_inq_varid  (ncid, 'dst_grid_dims', varid), subname)
         call check_ret(nf_get_var_int (ncid, varid, grid_dims), subname)
         domain%ni = grid_dims(1)
         domain%nj = grid_dims(2)
         domain%is_2d = .true.
      end if

      call domain_init(domain, domain%ns)
      ns = domain%ns

      call check_ret(nf_inq_varid (ncid, 'xc_b', varid), subname)
      call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
      call convert_latlon(ncid, 'xc_b', domain%lonc)

      call check_ret(nf_inq_varid (ncid, 'yc_b', varid), subname)
      call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)
      call convert_latlon(ncid, 'yc_b', domain%latc)

      if (grid_rank == 2 ) then
         allocate(yv(4,ns), xv(4,ns))
         call check_ret(nf_inq_varid (ncid, 'yv_b', varid), subname)
         call check_ret(nf_get_var_double (ncid, varid, yv), subname)
         call check_ret(nf_inq_varid (ncid, 'xv_b', varid), subname)
         call check_ret(nf_get_var_double (ncid, varid, xv), subname)

         domain%lats(:) = yv(1,:)  
         call convert_latlon(ncid, 'yv_b', domain%lats(:))

         domain%latn(:) = yv(3,:)  
         call convert_latlon(ncid, 'yv_b', domain%latn(:))

         domain%lonw(:) = xv(1,:)
         call convert_latlon(ncid, 'xv_b', domain%lonw(:))

         domain%lone(:) = xv(2,:)
         call convert_latlon(ncid, 'xv_b', domain%lone(:))

         domain%edgen = maxval(domain%latn)
         domain%edgee = maxval(domain%lone)
         domain%edges = minval(domain%lats)
         domain%edgew = minval(domain%lonw)
         deallocate(yv,xv)
      end if
       
      call check_ret(nf_inq_varid (ncid, 'frac_b', varid), subname)
      call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)

      call check_ret(nf_inq_varid (ncid, 'mask_b', varid), subname)
      call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)

      call check_ret(nf_inq_varid (ncid, 'area_b', varid), subname)
      call check_ret(nf_get_var_double (ncid, varid, domain%area), subname)
      domain%area = domain%area * re**2
    end if
    domain%maskset = .true.
    domain%fracset = .true.

    call check_ret(nf_close(ncid), subname)

  end function domain_read_map

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_read
!
! !INTERFACE:
  subroutine domain_read(domain, fname, readmask)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
    use mkutilsMod, only : convert_latlon
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    character(len=*) ,intent(in)     :: fname
    logical,optional, intent(in)     :: readmask ! true => read mask instead of landmask for urban parameters
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    include 'netcdf.inc'
    integer :: i,j,n
    real(r8), allocatable :: lon1d(:)          ! local array for 1d lon
    real(r8), allocatable :: lat1d(:)          ! local array for 1d lat
    real(r8), allocatable :: xv(:,:)           ! local array for corner lons
    real(r8), allocatable :: yv(:,:)           ! local array for corner lats
    integer :: ncid                            ! netCDF file id
    integer :: varid                           ! netCDF variable id
    logical :: edgeNESWset                     ! local EDGE[NESW]
    logical :: lonlatset                       ! local lon(:,:), lat(:,:)
    logical :: llneswset                       ! local lat[ns],lon[we]
    logical :: landfracset                     ! local landfrac
    logical :: maskset                         ! local mask
    integer :: ndims                           ! number of dims for variable
    integer :: ier                             ! error status
    logical :: lreadmask                       ! local readmask
    character(len= 32) :: lonvar               ! name of 2-d longitude variable
    character(len= 32) :: latvar               ! name of 2-d latitude variable
    character(len= 32) :: subname = 'domain_read'
!-----------------------------------------------------------------

    lonlatset   = .false. 
    edgeNESWset = .false. 
    llneswset   = .false. 
    landfracset = .false. 
    maskset     = .false. 
    lreadmask   = .false.

    if (present(readmask)) then
       lreadmask = readmask
    end if

    call check_ret(nf_open(fname, 0, ncid), subname)

    call domain_read_dims(domain, ncid)
    call domain_init(domain, domain%ns)
    write(6,*) trim(subname),' initialized domain'

    ! ----- Set lat/lon variable ------

    lonvar = ' '
    latvar = ' '

    if (.not. lonlatset) then
       ier = nf_inq_varid (ncid, 'LONGXY', varid)
       if (ier == NF_NOERR) then
          lonvar = 'LONGXY'
          latvar = 'LATIXY'
          lonlatset = .true.
       end if
    end if

    if (.not. lonlatset) then
       ier = nf_inq_varid (ncid, 'lon', varid)
       if (ier == NF_NOERR) then
          lonvar = 'lon'
          latvar = 'lat'
          lonlatset = .true.
       end if
    end if

    if (.not. lonlatset) then
       ier = nf_inq_varid (ncid, 'LONGITUDE', varid)
       if (ier == NF_NOERR) then
          lonvar = 'LONGITUDE'
          latvar = 'LATITUDE'
          lonlatset = .true.
       end if
    end if

    if (.not. lonlatset) then
       write(6,*)'lon/lat values not set' 
       write(6,*)'currently assume either that lon/lat or LONGXY/LATIXY', &
            ' or LONGITUDE/LATITUDE variables are on input dataset'
       call abort()
    end if

    call check_ret(nf_inq_varid (ncid, lonvar, varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
    call convert_latlon(ncid, lonvar, domain%lonc)

    call check_ret(nf_inq_varid (ncid, latvar, varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)
    call convert_latlon(ncid, latvar, domain%latc)

    ! ----- Set landmask/landfrac  ------

    ier = nf_inq_varid (ncid, 'frac', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read frac'
       call check_ret(nf_inq_varid (ncid, 'frac', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    ier = nf_inq_varid (ncid, 'LANDFRAC', varid)
    if (ier == NF_NOERR) then
       if (landfracset) write(6,*) trim(subname),' WARNING, overwriting frac'
       landfracset = .true.
       write(6,*) trim(subname),' read LANDFRAC'
       call check_ret(nf_inq_varid (ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, domain%frac), subname)
    endif

    if (lreadmask) then
       ier = nf_inq_varid (ncid, 'mask', varid)
       if (ier == NF_NOERR) then
          if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
          maskset = .true.
          write(6,*) trim(subname),' read mask with lreadmask set'
          call check_ret(nf_inq_varid (ncid, 'mask', varid), subname)
          call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
       endif
    else
       ier = nf_inq_varid (ncid, 'mask', varid)
       if (ier == NF_NOERR) then
          if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
          maskset = .true.
          write(6,*) trim(subname),' read mask'
          call check_ret(nf_inq_varid (ncid, 'mask', varid), subname)
          call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
       endif
       ier = nf_inq_varid (ncid, 'LANDMASK', varid)
       if (ier == NF_NOERR) then
          if (maskset) write(6,*) trim(subname),' WARNING, overwriting mask'
          maskset = .true.
          write(6,*) trim(subname),' read LANDMASK'
          call check_ret(nf_inq_varid (ncid, 'LANDMASK', varid), subname)
          call check_ret(nf_get_var_int (ncid, varid, domain%mask), subname)
       endif
    end if

    call check_ret(nf_close(ncid), subname)

    ! ----- set derived variables ----

    if (.not.maskset.and.landfracset) then
       maskset = .true.
       where (domain%frac < flandmin)
          domain%mask = 0     !ocean
       elsewhere
          domain%mask = 1     !land
       endwhere
    endif

    if (.not.landfracset.and.maskset) then
       landfracset = .true.
       do n = 1,domain%ns
          if ( domain%mask(n) == 0 )then
             domain%frac(n) = 0._r8     !ocean
          else
             domain%frac(n) = 1._r8     !land
          end if
       end do
    endif
    domain%maskset = maskset
    domain%fracset = landfracset

  end subroutine domain_read

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_read_dims
!
! !INTERFACE:
  subroutine domain_read_dims(domain, ncid)
!
! !DESCRIPTION:
! get dimension size(s) from a domain file
! sets domain%ns, domain%is_2d; and (if 2-d) domain%ni and domain%nj
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    integer          ,intent(in)    :: ncid    ! ID of an open netcdf file
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
    logical :: dimset                          ! has dimension information been set?
    character(len= 32) :: subname = 'domain_read_dims'
!-----------------------------------------------------------------

    ! Assume unstructured grid
    domain%ni = -9999
    domain%nj = -9999
    domain%is_2d = .false.

    dimset = .false.

    ! Note: We use the first dimension that is found in the following list

    ! ----- First try to find 2-d info ------

    call domain_read_dims_2d(domain, dimset, ncid, 'lsmlon', 'lsmlat')
    call domain_read_dims_2d(domain, dimset, ncid, 'ni', 'nj')
    call domain_read_dims_2d(domain, dimset, ncid, 'lon', 'lat')

    ! ----- If we haven't found 2-d info, try to find 1-d info -----

    call domain_read_dims_1d(domain, dimset, ncid, 'num_pixels')

    ! ----- If we haven't found any info, abort -----

    if (.not. dimset) then
       write(6,*) trim(subname),' ERROR: dims not set'
       call abort()
    endif    

  contains

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_read_dims_2d
!
! !INTERFACE:
    subroutine domain_read_dims_2d(domain, dimset, ncid, lon_name, lat_name)
!
! !DESCRIPTION:
! Try to read 2-d dimension size information
!
! Checks whether the given lon_name is found in the netcdf file. If it is:
! (a) If dimset is already true, then it issues a warning and returns
! (b) If dimset is false, then this sets:
! - domain%ni
! - domain%nj
! - domain%ns
! - domain%is_2d
! - dimset = true
!
! If the given lon_name is not found, the above variables are left unchanged
!
! !ARGUMENTS:
      implicit none
      type(domain_type),intent(inout) :: domain
      logical          ,intent(inout) :: dimset    ! has dimension information been set?
      integer          ,intent(in)    :: ncid      ! ID of an open netCDF file
      character(len=*) ,intent(in)    :: lon_name
      character(len=*) ,intent(in)    :: lat_name
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
      include 'netcdf.inc'
      integer :: dimid                             ! netCDF dimension id
      integer :: nlon,nlat                         ! size
      integer :: ier                               ! error status
      
      character(len= 32) :: subname = 'domain_read_dims_2d'

!-----------------------------------------------------------------

      ier = nf_inq_dimid (ncid, lon_name, dimid)
      if (ier == NF_NOERR) then
         if (dimset) then
            write(6,*) trim(subname),' WARNING: dimension sizes already set; skipping ', &
                 lon_name, '/', lat_name
         else
            write(6,*) trim(subname),' read lon and lat dims from ', lon_name, '/', lat_name
            call check_ret(nf_inq_dimid  (ncid, lon_name, dimid), subname)
            call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
            call check_ret(nf_inq_dimid  (ncid, lat_name, dimid), subname)
            call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
            domain%ni = nlon
            domain%nj = nlat
            domain%ns = nlon * nlat
            domain%is_2d = .true.
            dimset = .true.
         end if
      endif
      
    end subroutine domain_read_dims_2d


!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_read_dims_1d
!
! !INTERFACE:
    subroutine domain_read_dims_1d(domain, dimset, ncid, dim_name)
!
! !DESCRIPTION:
! Try to read 1-d dimension size information
!
! Checks whether the given dim_name is found in the netcdf file. If it is:
! (a) If dimset is already true, then it issues a warning and returns
! (b) If dimset is false, then this sets:
! - domain%ns
! - domain%is_2d
! - dimset = true
!
! If the given dim_name is not found, the above variables are left unchanged
!
! !ARGUMENTS:
      implicit none
      type(domain_type),intent(inout) :: domain
      logical          ,intent(inout) :: dimset    ! has dimension information been set?
      integer          ,intent(in)    :: ncid      ! ID of an open netCDF file
      character(len=*) ,intent(in)    :: dim_name
!
! !REVISION HISTORY:
! Author: Bill Sacks
!
!
! !LOCAL VARIABLES:
!EOP
      include 'netcdf.inc'
      integer :: dimid                             ! netCDF dimension id
      integer :: npts                              ! size
      integer :: ier                               ! error status
      
      character(len= 32) :: subname = 'domain_read_dims_1d'

!-----------------------------------------------------------------

      ier = nf_inq_dimid (ncid, dim_name, dimid)
      if (ier == NF_NOERR) then
         if (dimset) then
            write(6,*) trim(subname),' WARNING: dimension sizes already set; skipping ', dim_name
         else
            write(6,*) trim(subname),' read 1-d length from ', dim_name
            call check_ret(nf_inq_dimid  (ncid, dim_name, dimid), subname)
            call check_ret(nf_inq_dimlen (ncid, dimid, npts), subname)
            domain%ns = npts
            domain%is_2d = .false.
            dimset = .true.
         end if
      endif
      
    end subroutine domain_read_dims_1d

  end subroutine domain_read_dims


!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_write
!
! !INTERFACE:
  subroutine domain_write(domain,fname)
!
! !DESCRIPTION:
! Write a domain to netcdf

! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain
    character(len=*)  ,intent(in)    :: fname
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: varid                          !netCDF variable id
    integer  :: ncid                           !netCDF file id
    integer  :: omode                          !netCDF output mode
    character(len= 32) :: subname = 'domain_write'
!-----------------------------------------------------------------

    call check_ret(nf_open(trim(fname), nf_write, ncid), subname)
    ! File will be in define mode. Set fill mode to "no fill" to optimize performance

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write domain fields 

    call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%area), subname)

    call check_ret(nf_inq_varid(ncid, 'LONGXY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%lonc), subname)

    call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%latc), subname)

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Close grid data dataset

    call check_ret(nf_close(ncid), subname)

  end subroutine domain_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_ret
!
! !INTERFACE:
  subroutine check_ret(ret, calling)
!
! !DESCRIPTION:
! Check return status from netcdf call
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer, intent(in) :: ret
    character(len=*)    :: calling
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error from ',trim(calling), ' rcode = ', ret, &
                 ' error = ', NF_STRERROR(ret)
       call abort()
    end if

  end subroutine check_ret

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_checksame
!
! !INTERFACE:
  subroutine domain_checksame( srcdomain, dstdomain, tgridmap )
!
! !DESCRIPTION:
! Check that the input domains agree with the input map
!
! USES:
    use mkgridmapMod, only : gridmap_type, gridmap_setptrs
! !ARGUMENTS:
    implicit none
    type(domain_type), intent(in) :: srcdomain ! input domain
    type(domain_type), intent(in) :: dstdomain ! output domain
    type(gridmap_type),intent(in) :: tgridmap  ! grid map
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
     integer :: na, nb, ns             ! gridmap sizes
     integer :: n, ni                  ! indices
     real(r8), pointer :: xc_src(:)    ! Source longitude
     real(r8), pointer :: yc_src(:)    ! Source latitude
     real(r8), pointer :: frac_src(:)  ! Source fraction
     integer,  pointer :: mask_src(:)  ! Source mask
     integer,  pointer :: src_indx(:)  ! Source index
     real(r8), pointer :: xc_dst(:)    ! Destination longitude
     real(r8), pointer :: yc_dst(:)    ! Destination latitude
     real(r8), pointer :: frac_dst(:)  ! Destination fraction
     integer,  pointer :: mask_dst(:)  ! Destination mask
     integer,  pointer :: dst_indx(:)  ! Destination index
     character(len= 32) :: subname = 'domain_checksame'

     ! tolerance for checking equality of lat & lon
     ! We allow for single-precision rounding-level differences (approx. 1.2e-7 relative
     ! error) For a value of 360 (max value for lat / lon), this means we can allow
     ! absolute errors of about 5e-5.
     real(r8), parameter :: eps = 5.e-5_r8


     if (srcdomain%set == unset) then
        write(6,*) trim(subname)//'ERROR: source domain is unset!'
        call abort()
     end if
     if (srcdomain%set == unset) then
        write(6,*) trim(subname)//'ERROR: destination domain is unset!'
        call abort()
     end if

     call gridmap_setptrs( tgridmap, nsrc=na, ndst=nb, ns=ns,    &
                           xc_src=xc_src, yc_src=yc_src,         &
                           xc_dst=xc_dst, yc_dst=yc_dst,         &
                           mask_src=mask_src, mask_dst=mask_dst, &
                           src_indx=src_indx, dst_indx=dst_indx  &
                         )
       
     if (srcdomain%ns /= na) then
        write(6,*) trim(subname)// &
              ' ERROR: input domain size and gridmap source size are not the same size'
        write(6,*)' domain size = ',srcdomain%ns
        write(6,*)' map src size= ',na
        call abort()
     end if
     if (dstdomain%ns /= nb) then
        write(6,*) trim(subname)// &              
           ' ERROR: output domain size and gridmap destination size are not the same size'
        write(6,*)' domain size = ',dstdomain%ns
        write(6,*)' map dst size= ',nb
        call abort()
     end if
     do n = 1,ns
        ni = src_indx(n)
        if ( srcdomain%maskset )then
           if (srcdomain%mask(ni) /= mask_src(ni)) then
              write(6,*) trim(subname)// &              
                 ' ERROR: input domain mask and gridmap mask are not the same at ni = ',ni
              write(6,*)' domain  mask= ',srcdomain%mask(ni)
              write(6,*)' gridmap mask= ',mask_src(ni)
              call abort()
           end if
        end if
        if (abs(srcdomain%lonc(ni) - xc_src(ni)) > eps) then
           write(6,*) trim(subname)// &
               ' ERROR: input domain lon and gridmap lon not the same at ni = ',ni
           write(6,*)' domain  lon= ',srcdomain%lonc(ni)
           write(6,*)' gridmap lon= ',xc_src(ni)
           call abort()
        end if
        if (abs(srcdomain%latc(ni) - yc_src(ni)) > eps) then
           write(6,*) trim(subname)// &               
               ' ERROR: input domain lat and gridmap lat not the same at ni = ',ni
           write(6,*)' domain  lat= ',srcdomain%latc(ni)
           write(6,*)' gridmap lat= ',yc_src(ni)
           call abort()
        end if
     end do
     do n = 1,ns
        ni = dst_indx(n)
        if ( dstdomain%maskset )then
           if (dstdomain%mask(ni) /= mask_dst(ni)) then
              write(6,*) trim(subname)// &                              
                  ' ERROR: output domain mask and gridmap mask are not the same at ni = ',ni
              write(6,*)' domain  mask= ',dstdomain%mask(ni)
              write(6,*)' gridmap mask= ',mask_dst(ni)
              call abort()
           end if
        end if
        if (abs(dstdomain%lonc(ni) - xc_dst(ni)) > eps) then
           write(6,*) trim(subname)// &
               ' ERROR: output domain lon and gridmap lon not the same at ni = ',ni
           write(6,*)' domain  lon= ',dstdomain%lonc(ni)
           write(6,*)' gridmap lon= ',xc_dst(ni)
           call abort()
        end if
        if (abs(dstdomain%latc(ni) - yc_dst(ni)) > eps) then
           write(6,*) trim(subname)// &
                ' ERROR: output domain lat and gridmap lat not the same at ni = ',ni
           write(6,*)' domain  lat= ',dstdomain%latc(ni)
           write(6,*)' gridmap lat= ',yc_dst(ni)
           call abort()
        end if
     end do

  end subroutine domain_checksame

end module mkdomainMod
