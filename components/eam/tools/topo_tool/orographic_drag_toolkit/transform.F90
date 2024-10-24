Module transform
use shr_kind_mod, only: r8 => shr_kind_r8
contains
!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereABPFromRLL
!
! Description:
!   Determine the (alpha,beta,panel) coordinate of a point on the sphere from
!   a given regular lat lon coordinate.
!
! Parameters:
!   lon - Coordinate longitude
!   lat - Coordinate latitude
!   alpha (OUT) - Alpha coordinate
!   beta (OUT) - Beta coordinate
!   ipanel (OUT) - Face panel
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereABPFromRLL(lon, lat, alpha, beta, ipanel, ldetermine_panel)
  use shr_kind_mod, only: r8 => shr_kind_r8
  IMPLICIT NONE
  
  REAL    (R8), INTENT(IN)  :: lon, lat
  REAL    (R8), INTENT(OUT) :: alpha, beta
  INTEGER :: ipanel
  LOGICAL, INTENT(IN) :: ldetermine_panel
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  REAL    (r8), PARAMETER :: rotate_cube = 0.0
  
  ! Local variables
  REAL    (R8) :: xx, yy, zz, pm
  REAL    (R8) :: sx, sy, sz
  INTEGER  :: ix, iy, iz
  
  ! Translate to (x,y,z) space
  xx = COS(lon-rotate_cube) * COS(lat)
  yy = SIN(lon-rotate_cube) * COS(lat)
  zz = SIN(lat)
  
  pm = MAX(ABS(xx), ABS(yy), ABS(zz))
  
  ! Check maximality of the x coordinate
  IF (pm == ABS(xx)) THEN
    IF (xx > 0) THEN; ix = 1; ELSE; ix = -1; ENDIF
  ELSE
    ix = 0
  ENDIF

  ! Check maximality of the y coordinate
  IF (pm == ABS(yy)) THEN
    IF (yy > 0) THEN; iy = 1; ELSE; iy = -1; ENDIF
  ELSE
    iy = 0
  ENDIF
    
  ! Check maximality of the z coordinate
  IF (pm == ABS(zz)) THEN
    IF (zz > 0) THEN; iz = 1; ELSE; iz = -1; ENDIF
  ELSE
    iz = 0
  ENDIF
  
  ! Panel assignments
  IF (ldetermine_panel) THEN
    IF (iz  ==  1) THEN
      ipanel = 6; sx = yy; sy = -xx; sz = zz
      
    ELSEIF (iz  == -1) THEN
      ipanel = 5; sx = yy; sy = xx; sz = -zz
      
    ELSEIF ((ix == 1) .AND. (iy /= 1)) THEN
      ipanel = 1; sx = yy; sy = zz; sz = xx
      
    ELSEIF ((ix == -1) .AND. (iy /= -1)) THEN
      ipanel = 3; sx = -yy; sy = zz; sz = -xx
      
    ELSEIF ((iy == 1) .AND. (ix /= -1)) THEN
      ipanel = 2; sx = -xx; sy = zz; sz = yy
      
    ELSEIF ((iy == -1) .AND. (ix /=  1)) THEN
      ipanel = 4; sx = xx; sy = zz; sz = -yy
      
    ELSE
      WRITE(*,*) 'Fatal Error: CubedSphereABPFromRLL failed'
      WRITE(*,*) '(xx, yy, zz) = (', xx, ',', yy, ',', zz, ')'
      WRITE(*,*) 'pm =', pm, ' (ix, iy, iz) = (', ix, ',', iy, ',', iz, ')'
      STOP
    ENDIF
  ELSE
    IF (ipanel  ==  6) THEN
      sx = yy; sy = -xx; sz = zz
    ELSEIF (ipanel  == 5) THEN
      sx = yy; sy = xx; sz = -zz
    ELSEIF (ipanel == 1) THEN
      sx = yy; sy = zz; sz = xx        
    ELSEIF (ipanel == 3) THEN
      sx = -yy; sy = zz; sz = -xx
    ELSEIF (ipanel == 2) THEN
      sx = -xx; sy = zz; sz = yy
    ELSEIF (ipanel == 4) THEN
      sx = xx; sy = zz; sz = -yy
    ELSE
      WRITE(*,*) "ipanel out of range",ipanel
      STOP
    END IF
  END IF
  
  ! Use panel information to calculate (alpha, beta) coords
  alpha = ATAN(sx / sz)
  beta = ATAN(sy / sz)
  
END SUBROUTINE CubedSphereABPFromRLL

!------------------------------------------------------------------------------
! SUBROUTINE EquiangularAllAreas
!
! Description:
!   Compute the area of all cubed sphere grid cells, storing the results in
!   a two dimensional array.
!
! Parameters: 
!   icube - Resolution of the cubed sphere
!   dA (OUT) - Output array containing the area of all cubed sphere grid cells
!------------------------------------------------------------------------------
SUBROUTINE EquiangularAllAreas(icube, dA)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)                           :: icube
  REAL (r8), DIMENSION(icube,icube), INTENT(OUT) :: dA
  
  ! Local variables
  INTEGER                       :: k, k1, k2
  REAL (r8)                          :: a1, a2, a3, a4
  REAL (r8), DIMENSION(icube+1,icube+1)  :: ang
  REAL (r8), DIMENSION(icube+1)      :: gp
  
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq   = 0.25*pi
  
  
  !#ifdef DBG 
  REAL (r8)   :: dbg1 !DBG
  !#endif
  
  ! Recall that we are using equi-angular spherical gridding
  !   Compute the angle between equiangular cubed sphere projection grid lines.
  DO k = 1, icube+1
    gp(k) = -piq + (pi/DBLE(2*(icube))) * DBLE(k-1)
  ENDDO
  
  DO k2=1,icube+1
    DO k1=1,icube+1
      ang(k1,k2) =ACOS(-SIN(gp(k1)) * SIN(gp(k2)))
    ENDDO
  ENDDO
  
  DO k2=1,icube
    DO k1=1,icube
      a1 =      ang(k1  , k2  )
      a2 = pi - ang(k1+1, k2  )
      a3 = pi - ang(k1  , k2+1)
      a4 =      ang(k1+1, k2+1)      
      ! area = r*r*(-2*pi+sum(interior angles))
      DA(k1,k2) = -2.0*pi+a1+a2+a3+a4
    ENDDO
  ENDDO
  
  !#ifdef DBG 
  ! Only for debugging - test consistency
  dbg1 = 0.0                           !DBG
  DO k2=1,icube
    DO k1=1,icube
      dbg1 = dbg1 + DA(k1,k2)         !DBG
    ENDDO
  ENDDO
  write(*,*) 'DAcube consistency: ',dbg1-4.0*pi/6.0 !DBG
  !#endif
END SUBROUTINE EquiangularAllAreas


!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereRLLFromABP
!
! Description:
!   Determine the lat lon coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - Alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   lon (OUT) - Calculated longitude
!   lat (OUT) - Calculated latitude
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereRLLFromABP(alpha, beta, ipanel, lon, lat)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE        
  REAL    (r8), INTENT(IN)  :: alpha, beta
  INTEGER     , INTENT(IN)  :: ipanel
  REAL    (r8), INTENT(OUT) :: lon, lat        
  ! Local variables
  REAL    (r8) :: xx, yy, zz, rotate_cube
  REAL    (r8), PARAMETER :: pi   = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: piq  = 0.25*pi
  
  rotate_cube = 0.0
  ! Convert to cartesian coordinates
  CALL CubedSphereXYZFromABP(alpha, beta, ipanel, xx, yy, zz)        
  ! Convert back to lat lon
  lat = ASIN(zz)
  if (xx==0.0.and.yy==0.0) THEN
    lon = 0.0
  else
    lon = ATAN2(yy, xx) +rotate_cube 
    IF (lon<0.0) lon=lon+2.0*pi
    IF (lon>2.0*pi) lon=lon-2.0*pi
  end if
END SUBROUTINE CubedSphereRLLFromABP

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereXYZFromABP
!
! Description:
!   Determine the Cartesian coordinate of a point on a sphere given its
!   (alpha,beta,panel) coordinate.
!
! Parameters:
!   alpha - Alpha coordinate
!   beta - Beta coordinate
!   panel - Cubed sphere panel id
!   xx (OUT) - Calculated x coordinate
!   yy (OUT) - Calculated y coordinate
!   zz (OUT) - Calculated z coordinate
!------------------------------------------------------------------------------
SUBROUTINE CubedSphereXYZFromABP(alpha, beta, ipanel, xx, yy, zz)
  use shr_kind_mod, only: r8 => shr_kind_r8        
  IMPLICIT NONE
  
  REAL    (r8), INTENT(IN)  :: alpha, beta
  INTEGER     , INTENT(IN)  :: ipanel
  REAL    (r8), INTENT(OUT) :: xx, yy, zz        
  ! Local variables
  REAL    (r8) :: a1, b1, pm
  REAL    (r8) :: sx, sy, sz       
  
  ! Convert to Cartesian coordinates
  a1 = TAN(alpha)
  b1 = TAN(beta)
  
  sz = (1.0 + a1 * a1 + b1 * b1)**(-0.5)
  sx = sz * a1
  sy = sz * b1        
  ! Panel assignments
  IF (ipanel == 6) THEN
    yy = sx; xx = -sy; zz = sz          
  ELSEIF (ipanel == 5) THEN
    yy = sx; xx = sy; zz = -sz          
  ELSEIF (ipanel == 1) THEN
    yy = sx; zz = sy; xx = sz          
  ELSEIF (ipanel == 3) THEN
    yy = -sx; zz = sy; xx = -sz          
  ELSEIF (ipanel == 2) THEN
    xx = -sx; zz = sy; yy = sz          
  ELSEIF (ipanel == 4) THEN
    xx = sx; zz = sy; yy = -sz          
  ELSE
    WRITE(*,*) 'Fatal Error: Panel out of range in CubedSphereXYZFromABP'
    WRITE(*,*) '(alpha, beta, panel) = (', alpha, ',', beta, ',', ipanel, ')'
    STOP
  ENDIF
END SUBROUTINE CubedSphereXYZFromABP


SUBROUTINE remove_duplicates_integer(n_in,f_in,n_out,f_out)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  integer,dimension(n_in), intent(in) :: f_in
  integer, intent(out) :: n_out
  integer,dimension(n_in), intent(out) :: f_out
  !
  ! local work space
  !
  integer :: k,i,j
  !
  ! remove duplicates in ipanel_tmp
  !
  k = 1
  f_out(1) = f_in(1)
  outer: do i=2,n_in
    do j=1,k
      !            if (f_out(j) == f_in(i)) then
      if (ABS(f_out(j)-f_in(i))<1.0E-10) then
        ! Found a match so start looking again
        cycle outer
      end if
    end do
    ! No match found so add it to the output
    k = k + 1
    f_out(k) = f_in(i)
  end do outer
  n_out = k
END SUBROUTINE remove_duplicates_integer

SUBROUTINE remove_duplicates_latlon(n_in,lon_in,lat_in,n_out,lon_out,lat_out,tiny,ldbg)
  use shr_kind_mod, only: r8 => shr_kind_r8
  integer, intent(in) :: n_in
  real(r8),dimension(n_in), intent(inout) :: lon_in,lat_in
  real, intent(in) :: tiny
  integer, intent(out) :: n_out
  real(r8),dimension(n_in), intent(out) :: lon_out,lat_out
  logical :: ldbg
  !
  ! local work space
  !
  integer :: k,i,j
  REAL    (r8), PARAMETER :: pi        = 3.14159265358979323846264338327
  REAL    (r8), PARAMETER :: pih       = 0.50*pi
  !
  ! for pole points: make sure the longitudes are identical so that algorithm below works properly
  !
  do i=2,n_in
    if (abs(lat_in(i)-pih)<tiny.or.abs(lat_in(i)+pih)<tiny) then 
      lon_in(i) = lon_in(i-1)    
      write(*,*) "pole fix"
    end if
  end do

  lon_out = -9999999.9
  lat_out = -9999999.9
  !
  k = 1
  lon_out(1) = lon_in(1)
  lat_out(1) = lat_in(1)
  outer: do i=2,n_in
    do j=1,k
      if (ABS(lon_out(j)-lon_in(i))<tiny.AND.ABS(lat_out(j)-lat_in(i))<tiny) then
        ! Found a match so start looking again
        cycle outer
      end if
    end do
    ! No match found so add it to the output
    k = k + 1
    lon_out(k) = lon_in(i)
    lat_out(k) = lat_in(i)
  end do outer
  n_out = k
END SUBROUTINE remove_duplicates_latlon


end Module
