module mkdomainMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkdomainMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mkfileutils , only : getfil
  use mkvarpar    , only : re
  use mknanMod
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain1_type

  type domain1_type
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
  end type domain1_type

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain1_init          
  public domain1_clean         
  public domain1_check         
  public domain1_read 
  public domain1_read_map
  public domain1_write         
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
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_init
!
! !INTERFACE:
  subroutine domain1_init(domain,ns)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain1_type) :: domain        ! domain datatype
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
       call domain1_clean(domain)
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
       write(6,*) 'domain1_init ERROR: allocate mask, frac, lat, lon, area '
    endif
    
    domain%ns       = ns
    domain%mask     = -9999
    domain%frac     = -1.0e36
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan
    domain%set      = set
    
  end subroutine domain1_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_clean
!
! !INTERFACE:
  subroutine domain1_clean(domain)
!
! !DESCRIPTION:
! This subroutine deallocates the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain1_type) :: domain        ! domain datatype
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
       write(6,*) 'domain1_clean: cleaning ',domain%ns
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
          write(6,*) 'domain1_clean ERROR: deallocate mask, frac, lat, lon, area '
          call abort()
       endif
    else
       write(6,*) 'domain1_clean WARN: clean domain unecessary '
    endif

    domain%ns         = bigint
    domain%set        = unset

end subroutine domain1_clean

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_check
!
! !INTERFACE:
  subroutine domain1_check(domain)
!
! !DESCRIPTION:
! This subroutine write domain info
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain1_type),intent(in)  :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

    write(6,*) '  domain1_check set       = ',trim(domain%set)
    write(6,*) '  domain1_check ns        = ',domain%ns
    write(6,*) '  domain1_check lonc      = ',minval(domain%lonc),maxval(domain%lonc)
    write(6,*) '  domain1_check latc      = ',minval(domain%latc),maxval(domain%latc)
    write(6,*) '  domain1_check mask      = ',minval(domain%mask),maxval(domain%mask)
    write(6,*) '  domain1_check frac      = ',minval(domain%frac),maxval(domain%frac)
    write(6,*) ' '

end subroutine domain1_check

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_read_map
!
! !INTERFACE:
  subroutine domain1_read_map(domain, fname)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain1_type),intent(inout) :: domain
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
    character(len=256) :: locfn                ! local file name
    character(len= 32) :: subname = 'domain1_read'
!-----------------------------------------------------------------

    ! Read domain file and compute stuff as needed

    call getfil (fname, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    ! Assume unstructured grid

    domain%ni = -9999
    domain%nj = -9999

    call check_ret(nf_inq_dimid  (ncid, 'n_b', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, domain%ns), subname)
    
    call check_ret(nf_inq_dimid  (ncid, 'dst_grid_rank', dimid), subname)
    call check_ret(nf_inq_dimlen (ncid, dimid, grid_rank), subname)

    if (grid_rank == 2) then
       call check_ret(nf_inq_varid  (ncid, 'dst_grid_dims', varid), subname)
       call check_ret(nf_get_var_int (ncid, varid, grid_dims), subname)
       domain%ni = grid_dims(1)
       domain%nj = grid_dims(2)
    end if

    ns = domain%ns
    call domain1_init(domain, ns)

    call check_ret(nf_inq_varid (ncid, 'xc_b', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
    call check_ret(nf_inq_varid (ncid, 'yc_b', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)

    if (grid_rank == 2 ) then
       allocate(yv(4,ns), xv(4,ns))
       call check_ret(nf_inq_varid (ncid, 'yv_b', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, yv), subname)
       call check_ret(nf_inq_varid (ncid, 'xv_b', varid), subname)
       call check_ret(nf_get_var_double (ncid, varid, xv), subname)
       domain%lats(:) = yv(1,:)  
       domain%latn(:) = yv(3,:)  
       domain%lonw(:) = xv(1,:)
       domain%lone(:) = xv(2,:)
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

    call check_ret(nf_close(ncid), subname)

  end subroutine domain1_read_map

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_read
!
! !INTERFACE:
  subroutine domain1_read(domain, fname, readmask)
!
! !DESCRIPTION:
! Read a grid file
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain1_type),intent(inout) :: domain
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
    integer :: nlon,nlat                       ! size
    integer :: ns                              ! size of domain
    real(r8), allocatable :: lon1d(:)          ! local array for 1d lon
    real(r8), allocatable :: lat1d(:)          ! local array for 1d lat
    real(r8), allocatable :: xv(:,:)           ! local array for corner lons
    real(r8), allocatable :: yv(:,:)           ! local array for corner lats
    character(len=256) :: locfn                ! local file name
    integer :: ncid                            ! netCDF file id
    integer :: dimid                           ! netCDF dimension id
    integer :: varid                           ! netCDF variable id
    logical :: dimset                          ! local ni,nj
    logical :: edgeNESWset                     ! local EDGE[NESW]
    logical :: lonlatset                       ! local lon(:,:), lat(:,:)
    logical :: llneswset                       ! local lat[ns],lon[we]
    logical :: landfracset                     ! local landfrac
    logical :: maskset                         ! local mask
    integer :: ndims                           ! number of dims for variable
    integer :: ier                             ! error status
    logical :: lreadmask                       ! local readmask
    character(len= 32) :: subname = 'domain1_read'
!-----------------------------------------------------------------

    dimset      = .false. 
    lonlatset   = .false. 
    edgeNESWset = .false. 
    llneswset   = .false. 
    landfracset = .false. 
    maskset     = .false. 
    lreadmask   = .false.

    if (present(readmask)) then
       lreadmask = readmask
    end if

    ! Read domain file and compute stuff as needed

    call getfil (fname, locfn, 0)
    call check_ret(nf_open(locfn, 0, ncid), subname)

    ! Assume unstructured grid
    domain%ni = -9999
    domain%nj = -9999

    ! ----- Set lat/lon dimension ------

    ier = nf_inq_dimid (ncid, 'lon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lon and lat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
       domain%ni = nlon
       domain%nj = nlat
    endif

    ier = nf_inq_dimid (ncid, 'ni', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read ni and nj dims'
       call check_ret(nf_inq_dimid  (ncid, 'ni', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'nj', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
       domain%ni = nlon
       domain%nj = nlat
    endif

    ier = nf_inq_dimid (ncid, 'lsmlon', dimid)
    if (ier == NF_NOERR) then
       if (dimset) write(6,*) trim(subname),' WARNING, overwriting dims'
       dimset = .true.
       write(6,*) trim(subname),' read lsmlon and lsmlat dims'
       call check_ret(nf_inq_dimid  (ncid, 'lsmlon', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlon), subname)
       call check_ret(nf_inq_dimid  (ncid, 'lsmlat', dimid), subname)
       call check_ret(nf_inq_dimlen (ncid, dimid, nlat), subname)
       domain%ni = nlon
       domain%nj = nlat
    endif

    if (dimset) then
       write(6,*) trim(subname),' initialized domain'
       call domain1_init(domain,nlon*nlat)
    else
       write(6,*) trim(subname),' ERROR: dims not set for domain1_init'
       stop
    endif
    ns = domain%ns

    ! ----- Set lat/lon variable ------

    call check_ret(nf_inq_varid (ncid, 'LONGXY', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%lonc), subname)
    call check_ret(nf_inq_varid (ncid, 'LATIXY', varid), subname)
    call check_ret(nf_get_var_double (ncid, varid, domain%latc), subname)

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
       do n = 1,ns
          if ( domain%mask(n) == 0 )then
             domain%frac(n) = 0._r8     !ocean
          else
             domain%frac(n) = 1._r8     !land
          end if
       end do
    endif

  end subroutine domain1_read

!----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain1_write
!
! !INTERFACE:
  subroutine domain1_write(domain,fname)
!
! !USES:
!
! !DESCRIPTION:
! Write a domain to netcdf

! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain1_type),intent(inout) :: domain
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
    character(len= 32) :: subname = 'domain1_write'
!-----------------------------------------------------------------

    call check_ret(nf_open(trim(fname), nf_write, ncid), subname)
    ! File will be in define mode. Set fill mode to "no fill" to optimize performance

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write domain fields 

    call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%area), subname)

    call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%latn), subname)

    call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%lone), subname)

    call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%lats), subname)

    call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%lonw), subname)

    call check_ret(nf_inq_varid(ncid, 'LONGXY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%lonc), subname)

    call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, domain%latc), subname)

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! Close grid data dataset

    call check_ret(nf_close(ncid), subname)

  end subroutine domain1_write

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

end module mkdomainMod
