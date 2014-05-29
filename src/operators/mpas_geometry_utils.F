! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_geometry_utils

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_io_units

   contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{
   
      implicit none
   
      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
   
      real (kind=RKIND) :: a, b, c          ! Side lengths of spherical triangle ABC
   
      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
   
      real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC
   
      real (kind=RKIND) :: s                ! Semiperimeter of the triangle
      real (kind=RKIND) :: sin_angle
   
      a = acos(max(min(bx*cx + by*cy + bz*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (3)
      b = acos(max(min(ax*cx + ay*cy + az*cz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (2)
      c = acos(max(min(ax*bx + ay*by + az*bz,1.0_RKIND),-1.0_RKIND))      ! Eqn. (1)
   
      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
   
      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
   
      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)
   
      s = 0.5*(a + b + c)
!      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   ! Eqn. (28)
      sin_angle = sqrt(min(1.0_RKIND,max(0.0_RKIND,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   ! Eqn. (28)
   
      if ((Dx*ax + Dy*ay + Dz*az) >= 0.0) then
         mpas_sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      else
         mpas_sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      end if
   
   end function mpas_sphere_angle!}}}
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_PLANE_ANGLE
   !
   ! Computes the angle between vectors AB and AC, given points A, B, and C, and
   !   a vector (u,v,w) normal to the plane.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_plane_angle(ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w)!{{{
   
      implicit none
   
      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w
   
      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: mAB              ! The magnitude of AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=RKIND) :: mAC              ! The magnitude of AC
   
      real (kind=RKIND) :: Dx               ! The i-components of the cross product AB x AC
      real (kind=RKIND) :: Dy               ! The j-components of the cross product AB x AC
      real (kind=RKIND) :: Dz               ! The k-components of the cross product AB x AC
   
      real (kind=RKIND) :: cos_angle
   
      ABx = bx - ax
      ABy = by - ay
      ABz = bz - az
      mAB = sqrt(ABx**2.0 + ABy**2.0 + ABz**2.0)
   
      ACx = cx - ax
      ACy = cy - ay
      ACz = cz - az
      mAC = sqrt(ACx**2.0 + ACy**2.0 + ACz**2.0)
   
   
      Dx =   (ABy * ACz) - (ABz * ACy)
      Dy = -((ABx * ACz) - (ABz * ACx))
      Dz =   (ABx * ACy) - (ABy * ACx)
   
      cos_angle = (ABx*ACx + ABy*ACy + ABz*ACz) / (mAB * mAC)
   
      if ((Dx*u + Dy*v + Dz*w) >= 0.0) then
         mpas_plane_angle =  acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      else
         mpas_plane_angle = -acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      end if
   
   end function mpas_plane_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION MPAS_ARC_LENGTH
   !
   ! Returns the length of the great circle arc from A=(ax, ay, az) to 
   !    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   !    same sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function mpas_arc_length(ax, ay, az, bx, by, bz)!{{{
   
      implicit none
   
      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz
   
      real (kind=RKIND) :: r, c
      real (kind=RKIND) :: cx, cy, cz
   
      cx = bx - ax
      cy = by - ay
      cz = bz - az

!      r = ax*ax + ay*ay + az*az
!      c = cx*cx + cy*cy + cz*cz
!
!      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

      r = sqrt(ax*ax + ay*ay + az*az)
      c = sqrt(cx*cx + cy*cy + cz*cz)
!      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
      mpas_arc_length = r * 2.0 * asin(c/(2.0*r))

   end function mpas_arc_length!}}}
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTine mpas_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mpas_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{
   
      implicit none
   
      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz
      real (kind=RKIND), intent(out) :: cx, cy, cz
   
      real (kind=RKIND) :: r           ! Radius of the sphere
      real (kind=RKIND) :: d           
   
      r = sqrt(ax*ax + ay*ay + az*az)
   
      cx = 0.5*(ax + bx)
      cy = 0.5*(ay + by)
      cz = 0.5*(az + bz)
   
      if (cx == 0. .and. cy == 0. .and. cz == 0.) then
         write(stderrUnit,*) 'Error: arc_bisect: A and B are diametrically opposite'
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if
   
   end subroutine mpas_arc_bisect!}}}


   subroutine mpas_poly_fit_2(a_in,b_out,weights_in,m,n,ne)!{{{

      implicit none

      integer, intent(in) :: m,n,ne
      real (kind=RKIND), dimension(ne,ne), intent(in) :: a_in, weights_in
      real (kind=RKIND), dimension(ne,ne), intent(out) :: b_out
   
      ! local storage
   
      real (kind=RKIND), dimension(m,n)  :: a
      real (kind=RKIND), dimension(n,m)  :: b
      real (kind=RKIND), dimension(m,m)  :: w,wt,h
      real (kind=RKIND), dimension(n,m)  :: at, ath
!      real (kind=RKIND), dimension(n,n)  :: ata, ata_inv, atha, atha_inv
      real (kind=RKIND), dimension(n,n)  :: ata, atha, atha_inv
      integer, dimension(n) :: indx
!      integer :: i,j
   
      if ( (ne < n) .or. (ne < m) ) then
         write(stderrUnit,*) ' error in poly_fit_2 inversion ',m,n,ne
         call mpas_dmpar_global_abort('ERROR: in subroutine poly_fit_2()')
      end if
   
!      a(1:m,1:n) = a_in(1:n,1:m) 
      a(1:m,1:n) = a_in(1:m,1:n)
      w(1:m,1:m) = weights_in(1:m,1:m) 
      b_out(:,:) = 0.   

      wt = transpose(w)
      h = matmul(wt,w)
      at = transpose(a)
      ath = matmul(at,h)
      atha = matmul(ath,a)
      
      ata = matmul(at,a)

!      if (m == n) then
!         call mpas_migs(a,n,b,indx)
!      else

         call mpas_migs(atha,n,atha_inv,indx)

         b = matmul(atha_inv,ath)

!         call mpas_migs(ata,n,ata_inv,indx)
!         b = matmul(ata_inv,at)
!      end if
      b_out(1:n,1:m) = b(1:n,1:m)

!     do i=1,n
!        write(stdoutUnit,*) ' i, indx ',i,indx(i)
!     end do
!
!     write(stdoutUnit,*) ' '

   end subroutine mpas_poly_fit_2!}}}


! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine mpas_migs (a,n,x,indx)!{{{
!
! subroutine to invert matrix a(n,n) with the inverse stored
! in x(n,n) in the output.  copyright (c) tao pang 2001.
!
  implicit none
  integer, intent (in) :: n
  integer :: i,j,k
  integer, intent (out), dimension (n) :: indx
  real (kind=RKIND), intent (inout), dimension (n,n):: a
  real (kind=RKIND), intent (out), dimension (n,n):: x
  real (kind=RKIND), dimension (n,n) :: b
!
  do i = 1, n
    do j = 1, n
      b(i,j) = 0.0
    end do
  end do
  do i = 1, n
    b(i,i) = 1.0
  end do
!
  call mpas_elgs (a,n,indx)
!
  do i = 1, n-1
    do j = i+1, n
      do k = 1, n
        b(indx(j),k) = b(indx(j),k)-a(indx(j),i)*b(indx(i),k)
      end do
    end do
  end do
!
  do i = 1, n
    x(n,i) = b(indx(n),i)/a(indx(n),n)
    do j = n-1, 1, -1
      x(j,i) = b(indx(j),i)
      do k = j+1, n
        x(j,i) = x(j,i)-a(indx(j),k)*x(k,i)
      end do
      x(j,i) =  x(j,i)/a(indx(j),j)
    end do
  end do
end subroutine mpas_migs!}}}


subroutine mpas_elgs (a,n,indx)!{{{
!
! subroutine to perform the partial-pivoting gaussian elimination.
! a(n,n) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! indx(n) records the pivoting order.  copyright (c) tao pang 2001.
!
  implicit none
  integer, intent (in) :: n
  integer :: i,j,k,itmp
  integer, intent (out), dimension (n) :: indx
  real (kind=RKIND) :: c1,pi,pi1,pj
  real (kind=RKIND), intent (inout), dimension (n,n) :: a
  real (kind=RKIND), dimension (n) :: c
!
! initialize the index
!
  do i = 1, n
    indx(i) = i
  end do
!
! find the rescaling factors, one from each row
!
  do i = 1, n
    c1= 0.0
    do j = 1, n
      c1 = max(c1,abs(a(i,j)))
    end do
    c(i) = c1
  end do
!
! search the pivoting (largest) element from each column
!
  do j = 1, n-1
    pi1 = 0.0
    do i = j, n
      pi = abs(a(indx(i),j))/c(indx(i))
      if (pi.gt.pi1) then
        pi1 = pi
        k   = i
      endif
    end do
!
! interchange the rows via indx(n) to record pivoting order
!
    itmp    = indx(j)
    indx(j) = indx(k)
    indx(k) = itmp
    do i = j+1, n
      pj  = a(indx(i),j)/a(indx(j),j)
!
! record pivoting ratios below the diagonal
!
      a(indx(i),j) = pj
!
! modify other elements accordingly
!
      do k = j+1, n
        a(indx(i),k) = a(indx(i),k)-pj*a(indx(j),k)
      end do
    end do
  end do
!
end subroutine mpas_elgs!}}}

!-------------------------------------------------------------

end module mpas_geometry_utils
