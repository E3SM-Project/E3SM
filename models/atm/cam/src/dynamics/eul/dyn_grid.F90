module dyn_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of dynamics computational grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code. Global column index range
!         defined using full (unreduced) grid. 
!
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
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plev
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none
   save

   integer, private   :: ngcols_d = 0     ! number of dynamics columns

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore


contains

  subroutine get_block_ldof_d(nlev, ldof)
    use pio, only : pio_offset	    
    integer, intent(in) ::  nlev
    integer(kind=pio_offset), intent(out) :: ldof(:)
    
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
   use pmgrid, only: plat

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(out) :: block_first  ! first (global) index used for blocks
   integer, intent(out) :: block_last   ! last (global) index used for blocks

!-----------------------------------------------------------------------
!  latitude slice block
   block_first = 1
   block_last  = plat

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
   use pmgrid,     only: plat, plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! global block id
   integer, intent(in) :: size         ! array size

   integer, intent(out):: cdex(size)   ! global column indices
!---------------------------Local workspace-----------------------------
!
    integer i,j                            ! loop indices
    integer n                              ! column index
!-----------------------------------------------------------------------
! block == latitude slice
   if (size < plon) then
      write(iulog,*)'GET_BLOCK_GCOL_D: array not large enough (', &
                          size,' < ',plon,' ) '
      call endrun
   else
      n = (blockid-1)*plon
      do i = 1,plon
         n = n + 1
         cdex(i) = n
      enddo
   endif
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
   use pmgrid, only: plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
   get_block_gcol_cnt_d = plon

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
      enddo
      do k=plev+2,lvlsiz
         levels(k) = -1
      enddo
   endif

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
   use pmgrid,     only: plat, plon

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: gcol     ! global column index
   integer, intent(in) :: cnt      ! size of blockid and bcid arrays

   integer, intent(out) :: blockid(cnt) ! block index
   integer, intent(out) :: bcid(cnt)    ! column index within block
   integer, intent(out), optional :: localblockid(cnt)
!---------------------------Local workspace-----------------------------
!
    integer jb                     ! loop index
!-----------------------------------------------------------------------
!  latitude slice block
   if (cnt < 1) then
      write(iulog,*)'GET_GCOL_BLOCK_D: arrays not large enough (', &
                          cnt,' < ',1,' ) '
      call endrun
   else
      blockid(1) = (gcol-1)/plon + 1
      bcid(1)    = gcol - (blockid(1)-1)*plon
      do jb=2,cnt
         blockid(jb) = -1
         bcid(jb)    = -1
      enddo
   endif
!
   return
   end subroutine get_gcol_block_d
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
!  latitude slice block
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
#if ( defined SPMD )
   use spmd_dyn, only: proc
#endif

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid  ! global block id

!-----------------------------------------------------------------------
!  latitude slice block
#if (defined SPMD)
   get_block_owner_d = proc(blockid)
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
   endif
   hdim1_d = plon
   hdim2_d = plat

   return
   end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
   subroutine get_horiz_grid_d(size,clat_d_out,clon_d_out,area_d_out, &
                               wght_d_out)

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
   use pmgrid,        only: plat, plon
   use rgrid,         only: nlon
   use commap,        only: clat, clon, londeg, latdeg, w
   use physconst,     only: pi, spval
   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in)   :: size             ! array sizes

   real(r8), intent(out), optional :: clat_d_out(size) ! column latitudes
   real(r8), intent(out), optional :: clon_d_out(size) ! column longitudes
   real(r8), intent(out), optional :: area_d_out(size) ! column surface 
                                                       !  area
   real(r8), intent(out), optional :: wght_d_out(size) ! column integration
                                                       !  weight
!---------------------------Local workspace-----------------------------
!
    integer i,j                      ! loop indices
    integer n                        ! column index
    real(r8) :: ns_vert(2,plon)      ! latitude grid vertices
    real(r8) :: ew_vert(2,plon)      ! longitude grid vertices
    real(r8) :: del_theta            ! difference in latitude at a grid point
    real(r8) :: del_phi              ! difference in longitude at a grid point
    real(r8), parameter :: degtorad=pi/180_r8
!-----------------------------------------------------------------------
    if(present(clon_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1,nlon(j)
                n = n + 1
                clon_d_out(n) = clon(i,j)
             enddo
          enddo
       else if(size == plon) then
          clon_d_out(:) = clon(:,1)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if(present(clat_d_out)) then
       if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1,nlon(j)
                n = n + 1
                clat_d_out(n) = clat(j)
             enddo
          enddo
       else if(size == plat) then
          clat_d_out(:) = clat(:)
       else
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
       end if
    end if
    if ( ( present(wght_d_out) ) ) then

       if(size==plat) then
          wght_d_out(:) = (0.5_r8*w(:)/nlon(:))* (4.0_r8*pi)
       else if(size == ngcols_d) then
          n = 0
          do j = 1,plat
             do i = 1,nlon(j)
                n = n + 1
                wght_d_out(n) = ( 0.5_r8*w(j)/nlon(j) ) * (4.0_r8*pi)
             enddo
          enddo
       end if
    end if
    if ( present(area_d_out) ) then
       if(size < ngcols_d) then
          write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
               size,' < ',ngcols_d,' ) '
          call endrun
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
             ns_vert(1,:nlon(1)) = -90.0_r8
          else
             ns_vert(1,:nlon(j)) = (latdeg(j) + latdeg(j-1) )*0.5_r8
          endif
          
          if (j .eq. plat) then
             ns_vert(2,:nlon(plat)) =  90.0_r8
          else
             ns_vert(2,:nlon(j)) = (latdeg(j) + latdeg(j+1) )*0.5_r8
          endif

          ! Longitude vertices
          ew_vert(:,:) = spval
          ew_vert(1,1)          = (londeg(1,j) - 360.0_r8 + londeg(nlon(j),j))*0.5_r8
          ew_vert(1,2:nlon(j))  = (londeg(1:nlon(j)-1,j)+ londeg(2:nlon(j),j))*0.5_r8
          ew_vert(2,:nlon(j)-1) = ew_vert(1,2:nlon(j))
          ew_vert(2,nlon(j))    = (londeg(nlon(j),j) + (360.0_r8 + londeg(1,j)))*0.5_r8
          
          do i = 1,nlon(j)
             n = n + 1
             del_phi = sin( ns_vert(2,i)*degtorad ) - sin( ns_vert(1,i)*degtorad )
             del_theta = ( ew_vert(2,i) - ew_vert(1,i) )*degtorad
             area_d_out(n) = del_theta*del_phi
          enddo

       enddo
    endif
!
    return
  end subroutine get_horiz_grid_d


!#######################################################################
   function get_dyn_grid_parm_real2d(name) result(rval)
     use commap, only : londeg, clon
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:,:)

     if(name.eq.'clon') then
        rval => clon
     else if(name.eq.'londeg') then
        rval => londeg
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real2d

!#######################################################################
   function get_dyn_grid_parm_real1d(name) result(rval)
     use commap, only : latdeg, clat, w
     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:)

     if(name.eq.'clat') then
        rval => clat
     else if(name.eq.'latdeg') then
        rval => latdeg
     else if(name.eq.'w') then
        rval => w
     else
        nullify(rval)
     end if
   end function get_dyn_grid_parm_real1d




   integer function get_dyn_grid_parm(name) result(ival)
     use pmgrid, only : beglat, endlat, plat, plon, plev, plevp
     character(len=*), intent(in) :: name

     if(name.eq.'beglat' .or. name .eq. 'beglatxy') then
        ival = beglat
     else if(name.eq.'endlat' .or. name .eq. 'endlatxy') then
        ival = endlat
     else if(name.eq.'plat') then
        ival = plat
     else if(name.eq.'plon' .or. name .eq. 'endlonxy') then
        ival = plon
     else if(name.eq.'plev') then
        ival = plev
     else if(name.eq.'plevp') then
        ival = plevp
     else if(name .eq. 'beglonxy') then
	ival = 1
     else
        ival = -1
     end if


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
  do i = 1,plon
     ii = ii+1
     if (present(cdex)) cdex(ii) = i + (j-1)*plon
     if (present(rlat)) rlat(ii) = clat(j)
     if (present(rlon)) rlon(ii) = clon(i,1)
  enddo

endsubroutine dyn_grid_get_elem_coords

!#######################################################################

end module dyn_grid
