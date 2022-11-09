module ConnectionSetType
!This module is for creating the 2D grid connections for lateral GW flow.
!Han Qiu 2021.12
use shr_kind_mod   , only : r8 => shr_kind_r8
use shr_infnan_mod  , only : isnan => shr_infnan_isnan,nan => shr_infnan_nan, assignment(=)
use decompMod      , only : bounds_type
implicit none
save
public
  
type, public :: connection_set_type
    Integer, pointer :: nconn  => null()         ! number of connections
    Integer, pointer :: grid_id_up(:)  => null()    ! list of ids of upwind cells
    Integer, pointer :: grid_id_dn(:)  => null()    ! list of ids of downwind cells
    Real(r8), pointer :: dist(:)   => null()      ! list of distance vectors
    Real(r8), pointer :: area(:)   => null()      ! list of areas of faces normal to distance vectors
    Real(r8), pointer :: uparea(:)   => null()     ! list of up cell areas of horizaontal faces 
    Real(r8), pointer :: downarea(:)   => null()   ! list of down cell areas of horizaontal faces 
    Real(r8), pointer :: dzg(:)   => null()        ! list of areas of dz between downwind and upwind cells
    Real(r8), pointer :: facecos(:)   => null()    !dot product of the cell face normal vector and cell centroid vector
    Real(r8), pointer :: vertcos(:)   => null()    !dot product of the cell face normal vector and cell centroid vector for vertical flux, the rank for vertcos
                                                   !is from 1 to column size which is different from rank of lateral faces
    Real(r8), pointer :: facesin(:)   => null()    !the sinusoidal value of the angle correponding to facecos
    Real(r8), pointer :: slope(:)   => null()      !slope for each connection
    Real(r8), pointer :: aniso(:)   => null()      !anisotropic ratio
contains
    procedure, public :: Init => col_connect_init
end type connection_set_type

public :: get_natveg_column_id 
 
type (connection_set_type), public, target :: conn   ! connection type

contains

! ************************************************************************** !
subroutine col_connect_init(this, bounds)
   use ColumnType , only:col_pp
   type(bounds_type), intent(in)    :: bounds
   class(connection_set_type)       :: this 
   Integer                          :: n, nn, iconn,begc,endc,begg, endg,ii,jj 
   Integer                          :: g,nx,ny
   Real(r8)                         :: x(36,51),y(36,51),zhc(35,50),zh(36,51),dx(35,49),dy(34,50),dz(35,50),slopex(35,49),slopey(34,50),slopexx(35,50),slopeyy(35,50)
   !Real(r8)                         :: x(3,51),y(3,51),zh(3,51),dx(2,49),dy(1,50),dz(2,50),slopex(2,49),slopey(1,50),slopexx(2,50),slopeyy(2,50)
   begc = bounds%begc;  endc = bounds%endc
   begg = bounds%begg;  endg = bounds%endg
   nx = 50
   ny = 35
   iconn=0
   n = (nx-1)*ny+(ny-1)*nx
   !n = endg-begg
   ! number of connections in each layer
   allocate(this%nconn)             ;  this%nconn = n                   
   allocate(this%grid_id_up(n)) ;  this%grid_id_up(:) = 0
   allocate(this%grid_id_dn(n)) ;  this%grid_id_dn(:) = 0
   allocate(this%area(n))       ;  this%area(:) = 0
   allocate(this%uparea(n))       ;  this%uparea(:) = 0
   allocate(this%downarea(n))       ;  this%downarea(:) = 0
   allocate(this%dist(n))       ;  this%dist(:) = 0
   allocate(this%dzg(n))        ;  this%dzg(:) = 0
   allocate(this%slope(n))        ;  this%slope(:) = 0
   allocate(this%facecos(n))        ;  this%facecos(:) = 0
   allocate(this%facesin(n))        ;  this%facesin(:) = 0
   allocate(this%vertcos(nx*ny))        ;  this%vertcos(:) = 0
   allocate(this%aniso(nx*ny))        ;  this%aniso(:) = 0
   !dz = 1.0_r8
   do ii = 1, nx+1
     do jj = 1,ny+1
        x (jj , ii) = ii*1000-1000
        y (jj , ii) = jj*1000-1000         
     end do
   end do

   nn = 0
   do jj = 1, ny
   do ii = 1, nx
            nn = nn+1
            zhc(jj,ii) = col_pp%topo_ele(nn)  !zhc is cell centered elevation
     end do
   end do

   do jj = 2, ny
   do ii = 2, nx
            !zh is cell edge elevation calculated by averaging neighbouring cells
            zh(jj,ii) = (zhc(jj-1,ii-1) + zhc(jj,ii) + zhc(jj,ii-1) + zhc(jj-1,ii))/4._r8   
     end do
   end do

  !define boundary edge elevation
   do jj = 2,ny 
       zh(jj,1) = (zhc(jj-1, 1) + zhc(jj, 1)) / 2._r8
       zh(jj,nx+1) = (zhc(jj-1, nx) + zhc(jj, nx)) / 2._r8
   enddo

   do ii = 2,nx 
       zh(1,ii) = (zhc(1, ii-1) + zhc(1,ii)) / 2._r8
       zh(ny+1,ii) = (zhc(ny, ii-1) + zhc(ny, jj-1)) / 2._r8
   enddo

   zh(1, 1) = zhc (1,1)
   zh(ny+1,nx+1) = zhc (ny,nx)
   zh(1,nx+1) = zhc (1,nx)
   zh(ny+1,1) = zhc (ny,1)

   do ii = 1, ny           
     do jj = 1,nx-1
      slopex(ii,jj) = (zh(ii,jj+2)+zh(ii+1,jj+2)-zh(ii,jj)-zh(ii+1,jj))/(x(ii,jj+2)+x(ii+1,jj+2)-x(ii,jj)-x(ii+1,jj))
      dx(ii,jj) = abs((x(ii,jj+2)+x(ii+1,jj+2)-x(ii,jj)-x(ii+1,jj))/4)
     end do
   end do
   
   if(ny > 1) then
   do ii = 1, ny-1
     do jj = 1,nx
     slopey(ii,jj) = (zh(ii+2,jj)+zh(ii+2,jj+1)-zh(ii,jj)-zh(ii,jj+1))/(y(ii+2,jj)+y(ii+2,jj+1)-y(ii,jj)-y(ii,jj+1))
     dy(ii,jj) = abs((y(ii+2,jj)+y(ii+2,jj+1)-y(ii,jj)-y(ii,jj+1))/4)
     end do
   end do
  endif

! grid rank e.g. nx*ny = 5*4
! 1 2 3 4 5
! 6 7 8 9 10
! 11 12 13 14 15
! 16 17 18 19 20

   if(ny>1) then
     do ii= 1, ny
      do jj = 1, nx-1
       !g = begg+iconn-1
       ! assume grid is counted along nx, cross-sectional area is dy*dz, top and bottom area is , dy-bar*dx
       iconn = iconn+1
       this%grid_id_up(iconn) = (ii-1)*nx+jj       !  Step-2: Eventually will need to read from surface dataset
       this%grid_id_dn(iconn) = (ii-1)*nx+jj+1  !  There is already some code that we will be able to
       this%area(iconn)       = abs((y(ii+1,jj+1)-y(ii,jj+1)))      ! cross-sectional length, multiply by dz is area
       this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*dx(ii,jj))/2.0_r8      ! up cell grid area
       this%downarea(iconn)   = abs((y(ii+1,jj+1)-y(ii,jj+1)+y(ii+1,jj+2)-y(ii,jj+2))*dx(ii,jj))/2.0_r8       ! down cell grid area
       this%dist(iconn)       = sqrt((slopex(ii,jj)*dx(ii,jj))**2+dx(ii,jj)**2)  ! triangle law 
       this%dzg(iconn)        = slopex(ii,jj)*dx(ii,jj)   ! down cell elevation - up cell elevation, positive go up hill, negative go down hill
       this%slope(iconn)      = slopex(ii,jj)
       this%facecos(iconn)    = 1._r8/sqrt(1._r8 + slopex(ii,jj)**2)
       this%facesin(iconn)    = 1._r8/sqrt(1._r8 + (1._r8/slopex(ii,jj))**2)
       enddo
     enddo
     do ii= 1, ny-1
      do jj = 1, nx
       !g = begg+iconn-1
       ! assume grid is counted along ny, cross-sectional area is dx*dz
       iconn = iconn+1
       this%grid_id_up(iconn) = jj+ nx*(ii-1)  !  Step-2: Eventually will need to read from surface dataset
       this%grid_id_dn(iconn) = jj+ nx*ii !  There is already some code that we will be able to
       this%area(iconn)       = abs((x(ii+1,jj+1)-x(ii+1,jj)))      ! cross-sectional length, multiply by dz is area
       this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*(x(ii+1,jj+1)-x(ii+1,jj)))/2.0;      ! up cell vertical area
       this%downarea(iconn)   = abs((y(ii+2,jj)-y(ii+1,jj)+y(ii+2,jj+1)-y(ii+1,jj+1))*(x(ii+1,jj+1)-x(ii+1,jj)))/2.0;   ! down cell vertical area       
       this%dist(iconn)       = sqrt((slopey(ii,jj)*dy(ii,jj))**2 + dy(ii, jj)**2) ! triangle law 
       this%dzg(iconn)        = slopey(ii,jj)*dy(ii,jj)
       this%slope(iconn)      = slopey(ii,jj)
       this%facecos(iconn)    = 1._r8/sqrt(1._r8 + slopey(ii,jj)**2)
       this%facesin(iconn)    = 1._r8/sqrt(1._r8 + (1._r8/slopey(ii,jj))**2)
       enddo
     enddo 
   else
     do ii= 1, ny
      do jj = 1, nx-1
       !g = begg+iconn-1
       ! assume grid is counted along nx
       iconn = iconn+1
       this%grid_id_up(iconn) = jj    !  Step-2: Eventually will need to read from surface dataset
       this%grid_id_dn(iconn) = jj+1  !  There is already some code that we will be able to
       this%area(iconn)       = abs(dx(jj,1)*dz(jj,1))      !  use to fill this data structure
       this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*dx(ii,jj))/2.0_r8      ! up cell vertical area, trapzoid
       this%downarea(iconn)   = abs((y(ii+1,jj+1)-y(ii,jj+1)+y(ii+1,jj+2)-y(ii,jj+2))*dx(ii,jj))/2.0_r8       ! down cell vertical area
       this%facecos(iconn)    = 1/sqrt(1 + slopex(ii,jj)**2)
       this%facesin(iconn)    = 1._r8/sqrt(1 + (1._r8/slopex(ii,jj))**2)
       this%dist(iconn)       = sqrt((slopex(jj,1)*dx(jj,1))**2+dx(ii,jj)**2)     ! triangle law !may not used for now
       this%dzg(iconn)        = slopex(ii,jj)*dx(ii,jj) 
       enddo
      enddo
    endif
  iconn = 0
 if (ny > 1) then
  do ii = 1, ny
    do jj = 1, nx
      iconn = iconn+1
      slopexx(ii,jj) = (zh(ii,jj+1)+zh(ii+1,jj+1)-zh(ii,jj)-zh(ii+1,jj))/(x(ii,jj+1)+x(ii+1,jj+1)-x(ii,jj)-x(ii+1,jj))
      slopeyy(ii,jj) = (zh(ii+1,jj)+zh(ii+1,jj+1)-zh(ii,jj)-zh(ii,jj+1))/(y(ii+1,jj)+y(ii+1,jj+1)-y(ii,jj)-y(ii,jj+1))
      this%vertcos(iconn) = 1._r8/sqrt(1 + slopeyy(ii,jj)**2) * 1._r8/sqrt(1 + slopexx(ii,jj)**2)
      !this%aniso(iconn) ! will calculate based on clay sand ratio     
    enddo
  enddo
 else
    do jj = 1, nx
      iconn = iconn+1
      slopexx(1,jj) = (zh(1,jj+1)+zh(2,jj+1)-zh(1,jj)-zh(2,jj))/(x(1,jj+1)+x(2,jj+1)-x(1,jj)-x(2,jj))
      this%vertcos(iconn) = 1._r8/sqrt(1 + slopexx(1,jj)**2)
      !this%aniso(iconn) ! will calculate based on clay sand ratio 
    enddo
endif

end subroutine col_connect_init

function get_natveg_column_id(id, bounds) result(id_out)

  implicit none
  type(bounds_type), intent(in)    :: bounds
  integer, intent(in) :: id
  integer :: id_out,begc,endc,begg,endg
  begc = bounds%begc;  endc = bounds%endc
  begg = bounds%begg;  endg = bounds%endg

  
  id_out=begg*16-16+id ! for 2D transect 1 processor for 3d need beg or use filter_hydrologyc


end function get_natveg_column_id
end module ConnectionSetType
