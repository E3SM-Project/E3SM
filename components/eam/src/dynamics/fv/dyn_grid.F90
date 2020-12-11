module dyn_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of dynamics computational grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code. Global column index range
!         defined using full (unreduced) grid. 
! 
! Entry points:
!      get_block_bounds_d       get first and last indices in global 
!                               block ordering
!      get_block_gcol_d         get column indices for given block
!      get_block_gcol_cnt_d     get number of columns in given block
!      get_block_lvl_cnt_d      get number of vertical levels in column
!      get_block_levels_d       get vertical levels in column
!      get_gcol_block_d         get global block indices and local columns 
!                               index for given global column index
!      get_gcol_block_cnt_d     get number of blocks containing data
!                               from a given global column index
!      get_block_owner_d        get process "owning" given block
!      get_horiz_grid_d         get horizontal grid coordinates
!      get_horiz_grid_dim_d     get horizontal dimensions of dynamics grid
!      dyn_grid_get_pref        get reference pressures for the dynamics grid
!      dyn_grid_get_elem_coords get coordinates of a specified element (latitude) 
!                               of the dynamics grid (lat slice of the block)
!      dyn_grid_find_gcols      finds nearest column for given lat/lon
!      dyn_grid_get_colndx      get global lat and lon coordinate and MPI process indices 
!                               corresponding to a specified global column index
!
! Author: John Drake and Patrick Worley
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use pmgrid,           only: plev
   use cam_abortutils,   only: endrun
   use cam_logfile,      only: iulog

  implicit none

  integer, private   :: ngcols_d = 0     ! number of dynamics columns
! WS 2006.04.12:  moved ptimelevels here from prognostics, which is gone
  integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore

  real(r8), parameter ::  D0_5                    =   0.5_r8
  real(r8), parameter ::  D90_0                   =  90.0_r8
  real(r8), parameter ::  D360_0                  = 360.0_r8

  ! The FV dynamics grids
  integer, parameter, public :: dyn_decomp         = 101
  integer, parameter, public :: dyn_stagger_decomp = 102 !Backward compatibility
  integer, parameter, public :: dyn_ustag_decomp   = 102
  integer, parameter, public :: dyn_vstag_decomp   = 103

contains

  subroutine initgrid

    ! !USES:
    use shr_kind_mod,   only: r8 => shr_kind_r8
    use physconst,      only: pi, rair
    use pmgrid
    use pspect
    use commap,         only: w, clat, clon, w_staggered, clat
    use commap,         only: clat_staggered, clon, latdeg, londeg
    use commap,         only: latdeg_st, londeg_st
    use cam_abortutils, only: endrun
    use rgrid,          only: fullgrid
    use spmd_utils,     only: masterproc
    use cam_logfile,    only: iulog
    implicit none

    !-----------------------------------------------------------------------

    ! !DESCRIPTION:
    !
    !   Initialize Model commons.
    ! 
    ! !REVISION HISTORY:
    !
    !   92.06.01      Bath          Creation from CCM1
    !   96.02.01      Buja          Modifications
    !   01.01.19      Lin           incorporated Rasch's bug fix for the 
    !                               "Gaussian" weights
    !   02.04.04      Sawyer        Removed comspe
    !
    !EOP
    !-----------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:

    real(r8), parameter ::  D0_0                     =  0.0_r8
    real(r8), parameter ::  D1_0                     =  1.0_r8
    real(r8), parameter ::  D1_5                     =  1.5_r8
    real(r8), parameter ::  D2_0                     =  2.0_r8
    real(r8), parameter ::  D90_0                    =  90.0_r8
    real(r8), parameter ::  D180_0                   = 180.0_r8
    real(r8), parameter ::  D360_0                   = 360.0_r8
    real(r8), parameter ::  D1EM8                    =  1.0e-8_r8

    integer i           ! longitude index
    integer j           ! Latitude index
    integer k           ! Level index
    integer m           ! Index for legendre array
    integer irow        ! Latitude pair index
    integer lat         ! Latitude index

    real(r8) dp             ! Spacing between latitudes
    real(r8) sum
    !
    !-----------------------------------------------------------------------

    ! L-R dynamics uses a regular latitude distribution (not gausian).
    ! The algorithm below is a bastardized version of LSM: map.F.

    dp = D180_0/(plat-1)
    do lat = 1, plat
      latdeg(lat) = -D90_0 + (lat-1)*dp
      clat(lat) = latdeg(lat)*pi/D180_0
    end do

    ! Calculate latitudes for the staggered grid

    do lat = 1, plat-1
      clat_staggered(lat) = (clat(lat) + clat(lat+1)) / D2_0
      latdeg_st     (lat) = clat_staggered(lat)*D180_0/pi
    end do

    ! Weights are defined as cos(phi)*(delta-phi)
    ! For a sanity check, the sum of w across all lats should be 2, or 1 across
    ! half of the latitudes.

    do lat = 2, plat-1
      w(lat) = sin(clat_staggered(lat)) - sin(clat_staggered(lat-1))
    end do
    w(1) = sin(clat_staggered(1)) + D1_0
    w(plat) = w(1)

    sum = D0_0
    do lat=1,plat
      if (masterproc) write(iulog,*) 'initcom: lat, clat, w ', lat, clat(lat), w(lat)
      sum = sum + w(lat)
    end do

    if (abs(sum - D2_0) > D1EM8) then
      write(iulog,*) 'INITCOM 1: weights do not sum to 2. sum=',sum
      call endrun
    end if

    dp = pi / real(plat-1,r8)
    do lat = 1, plat-1
      w_staggered(lat) = sin(clat(lat+1)) - sin(clat(lat))
    end do

    sum = D0_0
    do lat=1,plat-1
      sum = sum + w_staggered(lat)
    end do

    if (abs(sum - D2_0) > D1EM8) then
      write(iulog,*) 'INITCOM 2: weights do not sum to 2. sum=',sum
      call endrun
    end if
    !
    ! Determine whether full or reduced grid
    !
    fullgrid = .true.
    if (masterproc) write(iulog,*) 'Number of longitudes per latitude = ', plon
    !
    ! Longitude array
    !
    do lat=1,plat
      do i=1,plon
        londeg(i,lat) = (i-1)*D360_0/plon
        clon(i,lat)   = (i-1)*D2_0*pi/plon
      end do
    end do
    do lat=1,plat
      do i=1,splon
        londeg_st(i,lat) = (i-D1_5)*D360_0/splon
      end do
    end do

    return

  end subroutine initgrid

  subroutine get_block_ldof_d( nlev, ldof)
    use pmgrid, only : beglonxy, endlonxy, beglatxy, endlatxy
    use pio, only : pio_offset_kind
    integer, intent(in) :: nlev
    integer(kind=pio_offset_kind), intent(out) :: ldof(:)

    integer:: i,j,k, lcnt
#if 0
    lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
    allocate(ldof(lcnt))
    lcnt=0
    ldof(:)=0	
    do k=1,nlev
      do j=beglatxy,endlatxy
        do i=beglonxy, endlonxy
          lcnt=lcnt+1
          if(j.eq.1.and.hdim2_d==plat-1) then
            ldof(lcnt)=0
          else if(hdim2_d>1) then
            if(fileorder.eq.'xyz') then
              ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d+(k-1)*hdim1_d*hdim2_d
            else 
              ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d*nlev+(k-1)*hdim1_d
            end if
          else  ! lon, lev decomp used for history nacs
            ldof(lcnt)=i+(k-1)*hdim1_d
          end if
        end do
      end do
    end do
#endif

  end subroutine get_block_ldof_d


!========================================================================
!
  subroutine get_block_bounds_d(block_first,block_last)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return first and last indices used in global block ordering
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid, only: plat, spmd_on, nprxy_x, nprxy_y

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
    block_first = 1
    if (spmd_on .eq. 1) then
! Assume 1 block per subdomain
      block_last  = nprxy_x*nprxy_y
    else
      block_last  = plat
    end if

    return
  end subroutine get_block_bounds_d
!
!========================================================================
!
  subroutine get_block_gcol_d(blockid,size,cdex)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return list of dynamics column indices in given block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid,   only: spmd_on, plon, nprxy_x
    use spmd_utils, only: iam
#if ( defined SPMD )
    use spmd_dyn, only: lonrangexy, latrangexy, npes_xy
#endif

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size

    integer, intent(out):: cdex(size)   ! global column indices
!---------------------------Local workspace-----------------------------
!
    integer i,j                        ! block coordinates
    integer blksiz                     ! block size
    integer k,l                        ! loop indices
    integer n                          ! column index
!-----------------------------------------------------------------------

    if (spmd_on .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      blksiz = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
               (latrangexy(2,j)-latrangexy(1,j)+1)
      if (size < blksiz) then
        write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
                      size,' < ',blksiz,' ) '
        call endrun
      else
        n = 0
        do k=latrangexy(1,j),latrangexy(2,j)
          do l=lonrangexy(1,i),lonrangexy(2,i)
            n = n + 1
            cdex(n) = l + (k-1)*plon
          end do
        end do
      end if
#endif
    else
      if (size < plon) then
        write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
             size,' < ',plon,' ) '
        call endrun
      else
        n = (blockid-1)*plon
        do i = 1,plon
          n = n + 1
          cdex(i) = n
        end do
      end if
    end if
!
    return
  end subroutine get_block_gcol_d
!
!========================================================================
!
  integer function get_block_gcol_cnt_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of dynamics columns in indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid, only: spmd_on, nprxy_x, plon
#if ( defined SPMD )
    use spmd_dyn, only: lonrangexy, latrangexy
#endif

    implicit none

!------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id

!---------------------------Local workspace-----------------------------
    integer i, j

!-----------------------------------------------------------------------
    if (spmd_on .eq. 1) then
      j = (blockid-1) / nprxy_x + 1
      i = blockid - (j-1) * nprxy_x
#if ( defined SPMD )
      get_block_gcol_cnt_d = (lonrangexy(2,i)-lonrangexy(1,i)+1) *       &
           (latrangexy(2,j)-latrangexy(1,j)+1)
#endif
    else
      get_block_gcol_cnt_d = plon
    end if

    return
  end function get_block_gcol_cnt_d

!
!========================================================================
!
  integer function get_block_lvl_cnt_d(blockid,bcid)

!-----------------------------------------------------------------------
!
!
! Purpose: Return number of levels in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid    ! column index within block

!-----------------------------------------------------------------------
!  latitude slice block
    get_block_lvl_cnt_d = plev + 1

    return
  end function get_block_lvl_cnt_d
!
!========================================================================
!
  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)

!-----------------------------------------------------------------------
!
!
! Purpose: Return level indices in indicated column. If column
!          includes surface fields, then it is defined to also
!          include level 0.
!
! Method:
!
! Author: Patrick Worley
!
!-----------------------------------------------------------------------

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid    ! column index within block
    integer, intent(in) :: lvlsiz   ! dimension of levels array

    integer, intent(out) :: levels(lvlsiz) ! levels indices for block

!---------------------------Local workspace-----------------------------
!
    integer k                      ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
    if (lvlsiz < plev + 1) then
      write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
                          lvlsiz,' < ',plev + 1,' ) '
      call endrun
    else
      do k=0,plev
        levels(k+1) = k
      end do
      do k=plev+2,lvlsiz
        levels(k) = -1
      end do
    end if

    return
  end subroutine get_block_levels_d

!
!========================================================================
!
  subroutine get_gcol_block_d(gcol,cnt,blockid,bcid,localblockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return global block index and local column index
!          for global column index
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid, only: plon, spmd_on, nprxy_x, nprxy_y
#if ( defined SPMD )
    use spmd_dyn, only: lonrangexy, latrangexy
#endif

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: gcol     ! global column index
    integer, intent(in) :: cnt      ! size of blockid and bcid arrays

    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)
!---------------------------Local workspace-----------------------------
!
    integer i,j,ii,jj                   ! loop indices
    integer glon, glat                  ! global longitude and latitude
                                        ! indices
    integer ddlon                       ! number of longitudes in block
!-----------------------------------------------------------------------

!  lon/lat block
    if (cnt < 1) then
      write(iulog,*)'GET_GCOL_BLOCK_D: arrays not large enough (', cnt,' < ',1,' )'
      call endrun
    else
      if (spmd_on .eq. 1) then
! Determine global latitude and longitude coordinate indices from
! global column index
        glat = (gcol-1)/plon + 1
        glon = gcol - ((glat-1)*plon)

! Determine block coordinates (ii,jj), where ii ranges from 1 to
! nprxy_x and jj ranges from 1 to nprxy_y.
#if ( defined SPMD )
        ii=0
        do i=1,nprxy_x
          if (lonrangexy(1,i) .le. glon .and. glon .le. lonrangexy(2,i)) ii=i
        end do
        jj=0
        do j=1,nprxy_y
          if (latrangexy(1,j) .le. glat .and. glat .le. latrangexy(2,j)) jj=j
        end do
        if (ii .eq. 0 .or. jj .eq. 0) then
          write(iulog,*)'GET_GCOL_BLOCK_D: could not find block indices for (', &
                        glon,',',glat,' ) '
          call endrun
        end if

! Global block index
        blockid(1) = (jj-1)*nprxy_x+ii

! Local coordinates in block
        j = glat-latrangexy(1,jj)+1
        i = glon-lonrangexy(1,ii)+1
        ddlon = lonrangexy(2,ii)-lonrangexy(1,ii)+1

! Local column index in block
        bcid(1) = (j-1)*ddlon+i
!
#endif
      else
        glat = (gcol-1)/plon + 1
        glon = gcol - ((glat-1)*plon)

        blockid(1) = glat
        bcid(1)    = glon
      end if
!
      do j=2,cnt
        blockid(j) = -1
        bcid(j)    = -1
      end do
!
    end if
!
    return
  end subroutine get_gcol_block_d

  subroutine get_gcol_lat(gcid, lat)
    use pmgrid,     only: plat, plon
    use commap,     only: clat
    integer, intent(in) :: gcid(:)
    real(r8), intent(out) :: lat(:)
    integer :: k, glen, ilat 

    glen = size(gcid)

    do k=1,glen
      ilat = (gcid(k)-1)/plon + 1 
      lat(k) = CLAT(ilat)
    end do

  end subroutine get_gcol_lat

  subroutine get_gcol_lon(gcid, lon)
    use pmgrid,     only: plat, plon
    use commap,     only: clon
    integer, intent(in) :: gcid(:)
    real(r8), intent(out) :: lon(:)

    integer :: k, glen, ilat

    glen = size(gcid)

    do k=1,glen
      ilat = (gcid(k)-1)/plon + 1
      lon(k) = CLON(gcid(k) - ((ilat-1)*plon),1)
    end do

  end subroutine get_gcol_lon

!
!========================================================================
!
  integer function get_gcol_block_cnt_d(gcol)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return number of blocks contain data for the vertical column
!          with the given global column index
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: gcol     ! global column index
!-----------------------------------------------------------------------
!  lon/lat block
    get_gcol_block_cnt_d = 1

    return
  end function get_gcol_block_cnt_d

!
!========================================================================
!
  integer function get_block_owner_d(blockid)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return id of processor that "owns" the indicated block
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid, only: spmd_on
#if ( defined SPMD )
    use spmd_dyn, only: proc
#endif

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
    if (spmd_on .eq. 1) then
      get_block_owner_d = blockid - 1
    else
      get_block_owner_d = proc(blockid)
    end if
#else
    get_block_owner_d = 0
#endif

    return
  end function get_block_owner_d

!
!========================================================================
!
  subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Returns declared horizontal dimensions of computational grid.
!          Note that global column ordering is assumed to be compatible
!          with the first dimension major ordering of the 2D array.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid,     only: plat, plon

!------------------------------Arguments--------------------------------
    integer, intent(out) :: hdim1_d           ! first horizontal dimension
    integer, intent(out) :: hdim2_d           ! second horizontal dimension
!-----------------------------------------------------------------------
    if (ngcols_d == 0) then
      ngcols_d = plat*plon
    end if
    hdim1_d = plon
    hdim2_d = plat

    return
  end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
  subroutine get_horiz_grid_d(size,clat_d_out,clon_d_out, &
       area_d_out, wght_d_out, lat_d_out, lon_d_out)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Return latitude and longitude (in radians), column surface
!          area (in radians squared) and surface integration weights
!          for global column indices that will be passed to/from physics
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
    use pmgrid,     only: plat, plon
    use commap,     only: clat, clon, latdeg, londeg, w
    use physconst,  only: pi, spval

    implicit none
!------------------------------Arguments--------------------------------
    integer, intent(in)   :: size             ! array sizes

    real(r8), intent(out), optional :: clat_d_out(size) ! column latitudes
    real(r8), intent(out), optional :: clon_d_out(size) ! column longitudes

    real(r8), intent(out), optional :: area_d_out(size) ! column surface 
                                                        !  area
    real(r8), intent(out), optional :: wght_d_out(size) ! column integration
                                                        !  weight
    real(r8), intent(out), optional :: lat_d_out(size)  ! column deg latitudes
    real(r8), intent(out), optional :: lon_d_out(size)  ! column deg longitudes
!---------------------------Local workspace-----------------------------
!
    integer i,j                 ! loop indices
    integer n                   ! column index
    real(r8) :: ns_vert(2,plon) ! latitude grid vertices
    real(r8) :: ew_vert(2,plon) ! longitude grid vertices
    real(r8) :: del_theta       ! difference in latitude at a grid point
    real(r8) :: del_phi         ! difference in longitude at a grid point
    real(r8),parameter :: degtorad=pi/180.0_r8 ! convert degrees to radians
    character(len=128) :: errormsg

!-----------------------------------------------------------------------
!   if (size < ngcols_d) then
!      write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
!                          size,' < ',ngcols_d,' ) '
!      call endrun
!   else
    if(present(clon_d_out)) then
      if(size == ngcols_d) then
        n = 0
        do j = 1,plat
          do i = 1,plon
            n = n + 1
            clon_d_out(n) = clon(i,j)
          end do
        enddo
      else if(size == plon) then
        clon_d_out(:) = clon(:,1)
      else
        write(errormsg, '(a,4(i0,a))')'clon_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plon,') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      end if
    end if
    if(present(clat_d_out)) then
      if(size == ngcols_d) then
        n = 0
        do j = 1,plat
          do i = 1,plon
            n = n + 1
            clat_d_out(n) = clat(j)
          end do
        end do
      else if(size == plat) then
        clat_d_out(:) = clat(:)
      else
        write(errormsg, '(a,4(i0,a))')'clat_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plat,') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      end if
    end if
    if(size==plat .and. present(wght_d_out)) then
      wght_d_out(:) = w(:)
      return
    end if

    if ( ( present(area_d_out) ) .or. ( present(wght_d_out) ) ) then
      if ((size < ngcols_d) .and. present(area_d_out)) then
        write(errormsg, '(a,2(i0,a))')'area_d_out array size incorrect (',  &
             size, ' /= ', ngcols_d, ') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      else if ((size < ngcols_d) .and. present(area_d_out)) then
        write(errormsg, '(a,2(i0,a))')'wght_d_out array size incorrect (',  &
             size, ' /= ', ngcols_d, ') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      end if
      n = 0
      do j = 1,plat

        ! First, determine vertices of each grid point. 
        ! Verticies are ordered as follows: 
        ! ns_vert: 1=lower left, 2 = upper left
        ! ew_vert: 1=lower left, 2 = lower right

        ! Latitude vertices
        ns_vert(:,:) = spval
        if (j .eq. 1) then
          ns_vert(1,:plon)    = -D90_0 + (latdeg(1) - latdeg(2))*D0_5
        else
          ns_vert(1,:plon)    = (latdeg(j) + latdeg(j-1) )*D0_5
        end if

        if (j .eq. plat) then
          ns_vert(2,:plon) =  D90_0 + (latdeg(plat) - latdeg(plat-1))*D0_5
        else
          ns_vert(2,:plon)    = (latdeg(j) + latdeg(j+1) )*D0_5
        end if

        ! Longitude vertices
        ew_vert(:,:) = spval
        ew_vert(1,1)          = (londeg(1,j) - D360_0 + londeg(plon,j))*D0_5
        ew_vert(1,2:plon)  = (londeg(1:plon-1,j)+ londeg(2:plon,j))*D0_5
        ew_vert(2,:plon-1) = ew_vert(1,2:plon)
        ew_vert(2,plon)    = (londeg(plon,j) + (D360_0 + londeg(1,j)))*D0_5

        do i = 1,plon
          n = n + 1

          if (j .eq. 1) then
            del_phi = -sin( latdeg(j)*degtorad )    + sin( ns_vert(2,i)*degtorad )
          else if (j .eq. plat) then
            del_phi =  sin( latdeg(j)*degtorad )    - sin( ns_vert(1,i)*degtorad )
          else
            del_phi =  sin( ns_vert(2,i)*degtorad ) - sin( ns_vert(1,i)*degtorad )
          end if

          del_theta = ( ew_vert(2,i) - ew_vert(1,i) )*degtorad

          if ( present(area_d_out) ) area_d_out(n) = del_theta*del_phi
          if (present(wght_d_out) ) wght_d_out(n) = del_theta*del_phi
        end do

      end do

    end if

    if(present(lon_d_out)) then
      if(size == ngcols_d) then
        n = 0
        do j = 1, plat
          do i = 1, plon
            n = n + 1
            lon_d_out(n) = londeg(i,j)
          end do
        end do
      else if(size == plon) then
        lon_d_out(:) = londeg(:,1)
      else
        write(errormsg, '(a,4(i0,a))')'lon_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plon,') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      end if
    end if
    if(present(lat_d_out)) then
      if(size == ngcols_d) then
        n = 0
        do j = 1, plat
          do i = 1, plon
            n = n +  1
            lat_d_out(n) = latdeg(j)
          end do
        end do
      else if(size == plat) then
        lat_d_out(:) = latdeg(:)
      else
        write(errormsg, '(a,4(i0,a))')'lat_d_out array size incorrect (',    &
             size, ' /= ', ngcols_d, ' .and. ', size, ' /= ', plat,') '
        call endrun('GET_HORIZ_GRID_D: '//errormsg)
      end if
    end if

    return
  end subroutine get_horiz_grid_d


!#######################################################################
    function get_dyn_grid_parm_real2d(name) result(rval)
      use commap, only : londeg, londeg_st, clon
      character(len=*), intent(in) :: name
      real(r8), pointer :: rval(:,:)

      if(name.eq.'clon') then
        rval => clon
      else if(name.eq.'londeg') then
        rval => londeg
      else if(name.eq.'londeg_st') then
        rval => londeg_st
      else
        nullify(rval)
      end if
    end function get_dyn_grid_parm_real2d

    !#######################################################################
    function get_dyn_grid_parm_real1d(name) result(rval)
      use commap, only : latdeg, latdeg_st, clat_staggered, w_staggered, clat, w
      character(len=*), intent(in) :: name
      real(r8), pointer :: rval(:)

      if(name.eq.'clat') then
        rval => clat
      else if(name.eq.'latdeg') then
        rval => latdeg
      else if(name.eq.'latdeg_st') then
        rval => latdeg_st
      else if(name.eq.'clatdeg_staggered') then
        rval => latdeg_st
      else if(name.eq.'w') then
        rval => w
      else if(name.eq.'w_staggered') then
        rval => w_staggered
      else
        nullify(rval)
      end if
    end function get_dyn_grid_parm_real1d


   integer function get_dyn_grid_parm(name) result(ival)
     use pmgrid, only : beglat, endlat, plat, plon, plev, plevp, &
          splon, beglev, endlev, endlevp, &
          beglonxy, endlonxy, beglatxy, endlatxy
     character(len=*), intent(in) :: name

     if(name.eq.'splon') then
        ival = splon
     else if(name.eq.'beglev') then
        ival = beglev
     else if(name.eq.'endlev') then
        ival = endlev
     else if(name.eq.'endlevp') then
        ival = endlevp
     else if(name.eq.'beglonxy') then
        ival = beglonxy
     else if(name.eq.'endlonxy') then
        ival = endlonxy
     else if(name.eq.'beglatxy') then
        ival = beglatxy
     else if(name.eq.'endlatxy') then
        ival = endlatxy
     else if(name.eq.'beglat') then
        ival = beglat
     else if(name.eq.'endlat') then
        ival = endlat
     else if(name.eq.'plat') then
        ival = plat
     else if(name.eq.'plon') then
        ival = plon
     else if(name.eq.'plev') then
        ival = plev
     else if(name.eq.'plevp') then
        ival = plevp
     else	
        ival = -1
     end if
     return
   end function get_dyn_grid_parm

!#######################################################################

subroutine dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)

   ! return reference pressures for the dynamics grid

   use hycoef, only: hypi, hypm, nprlev

   ! arguments
   real(r8), intent(out) :: pref_edge(:) ! reference pressure at layer edges (Pa)
   real(r8), intent(out) :: pref_mid(:)  ! reference pressure at layer midpoints (Pa)
   integer,  intent(out) :: num_pr_lev   ! number of top levels using pure pressure representation

   integer :: k
   !-----------------------------------------------------------------------

   do k = 1, plev
      pref_edge(k) = hypi(k)
      pref_mid(k)  = hypm(k)
   end do
   pref_edge(plev+1) = hypi(plev+1)
   
   num_pr_lev = nprlev

end subroutine dyn_grid_get_pref

!#######################################################################

!-------------------------------------------------------------------------------
! This returns the lat/lon information (and corresponding MPI task numbers (owners)) 
! of the global model grid columns nearest to the input satellite coordinate (lat,lon)
!-------------------------------------------------------------------------------
subroutine dyn_grid_find_gcols( lat, lon, nclosest, owners, indx, jndx, rlat, rlon, idyn_dists ) 
  use spmd_utils,     only: iam
  use shr_const_mod,  only: SHR_CONST_PI, SHR_CONST_REARTH
  use pmgrid,         only: plon, plat

  real(r8), intent(in) :: lat
  real(r8), intent(in) :: lon
  integer, intent(in)  :: nclosest
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: indx(nclosest)
  integer, intent(out) :: jndx(nclosest)

  real(r8),optional, intent(out) :: rlon(nclosest)
  real(r8),optional, intent(out) :: rlat(nclosest)
  real(r8),optional, intent(out) :: idyn_dists(nclosest)

  real(r8) :: dist            ! the distance (in radians**2 from lat, lon)
  real(r8) :: latr, lonr      ! lat, lon inputs converted to radians
  integer  :: ngcols
  integer  :: i, j

  integer :: blockid(1), bcid(1), lclblockid(1)

  real(r8), allocatable :: clat_d(:), clon_d(:), distmin(:)
  integer, allocatable :: igcol(:)
  real(r8), parameter :: rad2deg = 180._r8/SHR_CONST_PI

  latr = lat/rad2deg
  lonr = lon/rad2deg

  ngcols = plon*plat
  allocate( clat_d(1:ngcols) )
  allocate( clon_d(1:ngcols) )
  allocate( igcol(nclosest) )
  allocate( distmin(nclosest) )

  call get_horiz_grid_d(ngcols, clat_d_out=clat_d, clon_d_out=clon_d)

  igcol(:)    = -999
  distmin(:) = 1.e10_r8

  do i = 1,ngcols

     ! Use the Spherical Law of Cosines to find the great-circle distance.
     dist = acos(sin(latr) * sin(clat_d(i)) + cos(latr) * cos(clat_d(i)) * cos(clon_d(i) - lonr)) * SHR_CONST_REARTH
     do j = nclosest, 1, -1
        if (dist < distmin(j)) then
           
           if (j < nclosest) then
              distmin(j+1) = distmin(j)
              igcol(j+1)    = igcol(j)
           end if
           
           distmin(j) = dist
           igcol(j)    = i
        else
           exit
        end if
     enddo

  enddo

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        ! get global lat and lon coordinate indices from global column index
        ! -- plon is global number of longitude grid points
        jndx(i) = (igcol(i)-1)/plon + 1 
        indx(i) = igcol(i) - (jndx(i)-1)*plon
     else
        jndx(i) = -1
        indx(i) = -1
     endif
     
     if ( present(rlat) ) rlat(i) = clat_d(igcol(i)) * rad2deg
     if ( present(rlon) ) rlon(i) = clon_d(igcol(i)) * rad2deg
     
     if (present(idyn_dists)) then
        idyn_dists(i) = distmin(i)
     end if
     
  end do

  deallocate( clat_d )
  deallocate( clon_d )
  deallocate( igcol )
  deallocate( distmin )

end subroutine dyn_grid_find_gcols

!#######################################################################
subroutine dyn_grid_get_colndx( igcol, nclosest, owners, indx, jndx ) 
  use spmd_utils, only: iam
  use pmgrid,     only: plon

  integer, intent(in)  :: nclosest
  integer, intent(in)  :: igcol(nclosest)
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: indx(nclosest)
  integer, intent(out) :: jndx(nclosest)

  integer  :: i
  integer :: blockid(1), bcid(1), lclblockid(1)

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))

     if ( iam==owners(i) ) then
        ! get global lat and lon coordinate indices from global column index
        ! -- plon is global number of longitude grid points
        jndx(i) = (igcol(i)-1)/plon + 1 
        indx(i) = igcol(i) - (jndx(i)-1)*plon
     else
        jndx(i) = -1
        indx(i) = -1
     endif

  end do

end subroutine dyn_grid_get_colndx
!#######################################################################

! this returns coordinates of a latitude slice of the block corresponding 
! to latitude index latndx 
 
subroutine dyn_grid_get_elem_coords( latndx, rlon, rlat, cdex )
  use commap, only : clat, clon
  use pmgrid, only : beglonxy, endlonxy
  use pmgrid, only : plon

  integer, intent(in) :: latndx ! lat  index

  real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the latndx slice
  real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the latndx slice
  integer, optional, intent(out) :: cdex(:) ! global column index

  integer :: i,ii,j

  if (present(cdex)) cdex(:) = -1
  if (present(rlat)) rlat(:) = -999._r8
  if (present(rlon)) rlon(:) = -999._r8

  j = latndx
  ii=0
  do i = beglonxy,endlonxy
     ii = ii+1
     if (present(cdex)) cdex(ii) = i + (j-1)*plon
     if (present(rlat)) rlat(ii) = clat(j)
     if (present(rlon)) rlon(ii) = clon(i,1)
  enddo

end subroutine dyn_grid_get_elem_coords

!#######################################################################

subroutine define_cam_grids()
  use pspect,           only: ptrm
  use pmgrid,           only: plon, plat, beglonxy, endlonxy, beglatxy, endlatxy
  use spmd_utils,       only: iam
#ifdef SPMD
  use spmd_dyn,         only: npes_xy
#endif
  use commap,           only: londeg, londeg_st, latdeg, latdeg_st
  use commap,           only: w_staggered, w
  use cam_grid_support, only: horiz_coord_t, horiz_coord_create, iMap
  use cam_grid_support, only: cam_grid_register, cam_grid_attribute_register

  integer                      :: i, j, ind
  type(horiz_coord_t), pointer :: lat_coord
  type(horiz_coord_t), pointer :: lon_coord
  type(horiz_coord_t), pointer :: slat_coord
  type(horiz_coord_t), pointer :: slon_coord
  integer(iMap),       pointer :: coord_map(:)
  integer(iMap),       pointer :: grid_map(:,:)
  real(r8), allocatable        :: latvals(:)
  real(r8),            pointer :: rattval(:)
  integer,             pointer :: iattval(:)
#ifndef SPMD
  integer, parameter           :: npes_xy = 1
#endif

  ! Note: not using get_horiz_grid_dim_d or get_horiz_grid_d since those
 !       are deprecated ('cause I said so' -- goldy)

  nullify(lat_coord)
  nullify(lon_coord)
  nullify(slat_coord)
  nullify(slon_coord)
  nullify(coord_map)
  nullify(grid_map)
  nullify(rattval)
  nullify(iattval)
  ! Grid for cell centers
  ! Make a map
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      grid_map(4, ind) = i
    end do
  end do
  ! Cell-centered latitude coordinate
  allocate(coord_map(endlatxy - beglatxy + 1))
  coord_map = (/ (i, i = beglatxy, endlatxy) /)
  lat_coord => horiz_coord_create('lat', '', plat, 'latitude',                &
       'degrees_north', beglatxy, endlatxy, latdeg(beglatxy:endlatxy),        &
       map=coord_map)
  nullify(coord_map)
  ! Cell-centered longitude coordinate
  if ((beglonxy == 1) .and. (endlonxy == plon)) then
    ! Assume this is not a distributed dimension
    lon_coord => horiz_coord_create('lon', '', plon, 'longitude',             &
         'degrees_east', beglonxy, endlonxy, londeg(beglonxy:endlonxy,1))
  else
    allocate(coord_map(endlonxy - beglonxy + 1))
    coord_map = (/ (i, i = beglonxy, endlonxy) /)
    lon_coord => horiz_coord_create('lon', '', plon, 'longitude',             &
         'degrees_east', beglonxy, endlonxy, londeg(beglonxy:endlonxy,1),     &
         map=coord_map)
    deallocate(coord_map)
    nullify(coord_map)
  end if
  ! Cell-centered grid
  call cam_grid_register('fv_centers', dyn_decomp, lat_coord, lon_coord,      &
       grid_map, unstruct=.false.)
  allocate(rattval(size(w)))
  rattval = w
  call cam_grid_attribute_register('fv_centers', 'gw', 'latitude weights', 'lat', rattval)
  nullify(rattval)
  nullify(grid_map) ! Belongs to the grid
  ! Staggered grid for U_S
  ! Make a map
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      if ((i == beglatxy) .and. (beglatxy == 1)) then
        grid_map(4, ind) = 0
      else
        grid_map(4, ind) = i - 1
      end if
    end do
  end do

  ! Staggered latitudes 'skip' the first one so they are 'off by one'
  ! This means we always must have a coordinate map
  if (iam .lt. npes_xy) then
    allocate(coord_map(endlatxy - beglatxy + 1))
    ! NB: coord_map(1) == 0 when beglat == 1, that element is not output
    do i = 1, size(coord_map)
      coord_map(i) = i + beglatxy - 2
    end do
  else
    allocate(coord_map(0))
  end if
  allocate(latvals(beglatxy:endlatxy))
  if (beglatxy == 1) then
    latvals(1) = 0
    latvals(2:endlatxy) = latdeg_st(1:endlatxy-1)
  else
    i = beglatxy - 1 ! Stupid NAG 'error'
    latvals(beglatxy:endlatxy) = latdeg_st(i:endlatxy-1)
  end if
  slat_coord => horiz_coord_create('slat', '', (plat - 1),                  &
       'staggered latitude', 'degrees_north', beglatxy, endlatxy, latvals,  &
       map=coord_map)
  deallocate(coord_map)
  nullify(coord_map)

  call cam_grid_register('fv_u_stagger', dyn_ustag_decomp, slat_coord,        &
       lon_coord, grid_map, unstruct=.false.)
  nullify(grid_map) ! Belongs to the grid
  deallocate(latvals)

  ! Staggered grid for V_S
  ! Make a map (need to do this because lat indices are distributed)
  allocate(grid_map(4, ((endlonxy - beglonxy + 1) * (endlatxy - beglatxy + 1))))
  ind = 0
  do i = beglatxy, endlatxy
    do j = beglonxy, endlonxy
      ind = ind + 1
      grid_map(1, ind) = j
      grid_map(2, ind) = i
      grid_map(3, ind) = j
      grid_map(4, ind) = i
    end do
  end do
  ! Staggered longitude coordinate
  if ((beglonxy == 1) .and. (endlonxy == plon)) then
    ! Assume this is not a distributed dimension
    slon_coord => horiz_coord_create('slon', '', plon, 'staggered longitude', &
         'degrees_east', beglonxy, endlonxy, londeg_st(beglonxy:endlonxy,1))
  else
    allocate(coord_map(endlonxy - beglonxy + 1))
    coord_map = (/ (i, i = beglonxy, endlonxy) /)
    slon_coord => horiz_coord_create('slon', '', plon, 'staggered longitude', &
         'degrees_east', beglonxy, endlonxy, londeg_st(beglonxy:endlonxy,1),  &
         map=coord_map)
    deallocate(coord_map)
    nullify(coord_map)
  end if
  call cam_grid_register('fv_v_stagger', dyn_vstag_decomp, lat_coord,         &
       slon_coord, grid_map, unstruct=.false.)
  call cam_grid_attribute_register('fv_v_stagger', 'w_stag',                  &
       'staggered latitude weights', 'lat', w_staggered)
  nullify(grid_map) ! Belongs to the grid

end subroutine define_cam_grids

subroutine physgrid_copy_attributes_d(gridname, grid_attribute_names)
  use cam_grid_support, only: max_hcoordname_len

  ! Dummy arguments
  character(len=max_hcoordname_len),          intent(out) :: gridname
  character(len=max_hcoordname_len), pointer, intent(out) :: grid_attribute_names(:)

  gridname = 'fv_centers'
  allocate(grid_attribute_names(1))
  grid_attribute_names(1) = 'gw'

end subroutine physgrid_copy_attributes_d

!#######################################################################

end module dyn_grid
