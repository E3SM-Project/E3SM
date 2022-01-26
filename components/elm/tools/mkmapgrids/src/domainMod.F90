
module domainMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domainMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use nanMod      , only : nan, bigint
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  character*16,parameter, private :: domain_unset = 'NOdomain_unsetNO'
  real(r8), parameter, private :: frac_init = -9999._r8
  integer,  parameter, private :: mask_init = bigint

  type domain_type
     integer          :: ns            ! global size of domain
     integer          :: ni,nj         ! size of global arrays (lsmlon,lsmlat)
     integer ,pointer :: mask(:)       ! land mask: 1 = land. 0 = ocean
     real(r8),pointer :: frac(:)       ! fractional land
     real(r8),pointer :: latc(:)       ! latitude of grid cell (deg)
     real(r8),pointer :: lonc(:)       ! longitude of grid cell (deg)
     real(r8),pointer :: lats(:)       ! grid cell latitude, S edge (deg)
     real(r8),pointer :: latn(:)       ! grid cell latitude, N edge (deg)
     real(r8),pointer :: lonw(:)       ! grid cell longitude, W edge (deg)
     real(r8),pointer :: lone(:)       ! grid cell longitude, E edge (deg)
  end type domain_type

!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_celledge_global
  public domain_celledge_regional
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
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
  subroutine domain_init(domain,ns,ni,nj)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
    integer           :: ns            ! grid size
    integer, optional :: ni
    integer, optional :: nj
!
! !REVISION HISTORY:
!   Created by T Craig
!
! !LOCAL VARIABLES:
    integer ier
!
!EOP
!------------------------------------------------------------------------------

    allocate(domain%mask(ns), &
             domain%frac(ns), &
             domain%latc(ns), &
             domain%lonc(ns), &
             domain%lats(ns), &
             domain%latn(ns), &
             domain%lonw(ns), &
             domain%lone(ns), stat=ier)
    if (ier /= 0) then
       write(6,*) 'domain_init ERROR: allocate lats, latn, lonw, lone'
       stop
    endif

    domain%ns       = ns
    if (present(ni)) domain%ni = ni
    if (present(nj)) domain%nj = nj
    domain%mask     = bigint
    domain%frac     = frac_init
    domain%latc     = nan
    domain%lonc     = nan
    domain%lats     = nan
    domain%latn     = nan
    domain%lonw     = nan
    domain%lone     = nan

end subroutine domain_init

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
! !USES:
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

    write(6,*) ' '
    write(6,*) 'domain_check ns    = ',domain%ns
    write(6,*) 'domain_check ni,nj = ',domain%ni,domain%nj
    write(6,*) 'domain_check lonc  = ',minval(domain%lonc),maxval(domain%lonc)
    write(6,*) 'domain_check latc  = ',minval(domain%latc),maxval(domain%latc)
    write(6,*) 'domain_check mask  = ',minval(domain%mask),maxval(domain%mask)
    write(6,*) 'domain_check frac  = ',minval(domain%frac),maxval(domain%frac)
    write(6,*) 'domain_check latn  = ',minval(domain%latn),maxval(domain%latn)
    write(6,*) 'domain_check lone  = ',minval(domain%lone),maxval(domain%lone)
    write(6,*) 'domain_check lats  = ',minval(domain%lats),maxval(domain%lats)
    write(6,*) 'domain_check lonw  = ',minval(domain%lonw),maxval(domain%lonw)

end subroutine domain_check

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_celledge_regional
!
! !INTERFACE:
  subroutine domain_celledge_regional (domain, edgen, edgee, edges, edgew)
!
! !DESCRIPTION:
! Southern and western edges of grid cells - regional grid
! (can become global as special case)
! Latitudes -- southern/northern edges for each latitude strip.
! For grids oriented South to North, the southern
! and northern edges of latitude strip [j] are:
!        southern = lats(j  )
!        northern = lats(j+1)
! For grids oriented North to South: the southern
! and northern edges of latitude strip [j] are:
!        northern = lats(j  )
!        southern = lats(j+1)
! In both cases, [lats] must be dimensioned lats(lat+1)
! Longitudes -- western edges. Longitudes for the western edge of the
! cells must increase continuously and span 360 degrees. Assume that
! grid starts at Dateline with western edge on Dateline Western edges
! correspond to [longxy] (longitude at center of cell) and range from
! -180 to 180 with negative longitudes west of Greenwich.
! Partial grids that do not span 360 degrees are allowed so long as they
! have the convention of Grid 1 with
!      western edge of grid: >= -180 and < 180
!      eastern edge of grid: > western edge  and <= 180
! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
! strip can have variable longitudinal resolution
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain
    real(r8), intent(in) :: edgen             !northern edge of grid (degrees)
    real(r8), intent(in) :: edgee             !eastern edge of grid (degrees)
    real(r8), intent(in) :: edges             !southern edge of grid (degrees)
    real(r8), intent(in) :: edgew             !western edge of grid (degrees)

!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005.11.20 Updated to domain datatype by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: longxy(:,:) 
    real(r8),pointer :: latixy(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j,n     
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    write(6,*) 'domain_celledge, using celledge_regional'

    nlon = domain%ni
    nlat = domain%nj
    allocate(longxy(nlon,nlat))
    allocate(latixy(nlon,nlat))
    allocate(lats(nlon,nlat))
    allocate(latn(nlon,nlat))
    allocate(lonw(nlon,nlat))
    allocate(lone(nlon,nlat))
    
    do j = 1,nlat
    do i = 1,nlon
       n = (j-1)*nlon + i
       longxy(i,j) = domain%lonc(n)
       latixy(i,j) = domain%latc(n)
    end do
    end do

    ! Latitudes
    ! Assumes lats are constant on an i line

    if (nlat == 1) then                      ! single latitude
       lats(:,1)    = edges
       latn(:,nlat) = edgen
    elseif (latixy(1,2) > latixy(1,1)) then  ! South to North grid
       lats(:,1)    = edges
       latn(:,nlat) = edgen
       do j = 2, nlat
          lats(:,j) = (latixy(1,j-1) + latixy(1,j)) / 2._r8
          latn(:,j-1) = lats(:,j)
       end do
    else                                     ! North to South grid
       latn(:,1)    = edgen
       lats(:,nlat) = edges
       do j = 2, nlat
          latn(:,j) = (latixy(1,j-1) + latixy(1,j)) / 2._r8
          lats(:,j-1) = latn(:,j)
       end do
    end if

    ! Longitudes
    ! Western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180. This is the same as [edgew].

    do j = 1, nlat
       lonw(1,j)    = edgew
       lone(nlon,j) = edgee
       dx = (edgee - edgew) / nlon
       do i = 2, nlon
          lonw(i,j)   = lonw(1,j) + (i-1)*dx
          lone(i-1,j) = lonw(i,j)
       end do
    end do

    do n = 1,domain%ns
       j = (n-1)/nlon + 1
       i = n - ((j-1)*nlon)
       domain%lonc(n) = longxy(i,j)
       domain%latc(n) = latixy(i,j)
       domain%lats(n) = lats(i,j)
       domain%latn(n) = latn(i,j)
       domain%lonw(n) = lonw(i,j)
       domain%lone(n) = lone(i,j)
    end do

    deallocate( longxy, latixy, lats, latn, lonw, lone )

  end subroutine domain_celledge_regional

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_celledge_global
!
! !INTERFACE:
  subroutine domain_celledge_global (domain)
!
! !DESCRIPTION:
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain

!
! !REVISION HISTORY:
! 2005.11.20 Updated to domain datatype by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: nlon              
    integer  :: nlat              
    real(r8),pointer :: longxy(:,:) 
    real(r8),pointer :: latixy(:,:) 
    real(r8),pointer :: lats(:,:)   
    real(r8),pointer :: latn(:,:)   
    real(r8),pointer :: lonw(:,:)   
    real(r8),pointer :: lone(:,:)   
    integer i,j,n           !indices
    real(r8) dx             !cell width
!------------------------------------------------------------------------

    write(6,*) 'domain_celledge, using celledge_global'

    nlon = domain%ni
    nlat = domain%nj
    write(6,*)'nlon= ',nlon,' nlat= ',nlat
    allocate(longxy(nlon,nlat))
    allocate(latixy(nlon,nlat))
    allocate(lats(nlon,nlat))
    allocate(latn(nlon,nlat))
    allocate(lonw(nlon,nlat))
    allocate(lone(nlon,nlat))

    do j = 1,nlat
    do i = 1,nlon
       n = (j-1)*nlon + i
       longxy(i,j) = domain%lonc(n)
       latixy(i,j) = domain%latc(n)
    end do
    end do

    ! Latitudes
    lats(:,1)    = -90._r8
    latn(:,nlat) =  90._r8
    do j = 2, nlat   
    do i = 1, nlon
       lats(i,j) = (latixy(i,j-1) + latixy(i,j)) / 2._r8
       latn(i,j-1) = lats(i,j)
    enddo
    enddo

    ! Longitudes

    do j = 1, nlat
      lonw(1,j)    = 1.5_r8*longxy(1,j) - 0.5_r8*longxy(2,j)
      lone(nlon,j) = lonw(1,j) + 360._r8
    enddo

    do j = 1, nlat   
    do i = 2, nlon
      lonw(i,j) = 0.5_r8*longxy(i-1,j) + 0.5_r8*longxy(i,j)
      lone(i-1,j) = lonw(i,j)
    enddo
    enddo

    do n = 1,domain%ns
       j = (n-1)/nlon + 1
       i = n - ((j-1)*nlon)
       domain%lonc(n) = longxy(i,j)
       domain%latc(n) = latixy(i,j)
       domain%lats(n) = lats(i,j)
       domain%latn(n) = latn(i,j)
       domain%lonw(n) = lonw(i,j)
       domain%lone(n) = lone(i,j)
    end do

    deallocate( longxy, latixy, lats, latn, lonw, lone )

  end subroutine domain_celledge_global

!-----------------------------------------------------------------------

end module domainMod





