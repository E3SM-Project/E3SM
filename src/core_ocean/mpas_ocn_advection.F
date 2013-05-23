! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_advection

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_constants


   contains


   subroutine ocn_initialize_advection_rk( grid, err )!{{{
                                      
!
! compute the cell coefficients for the polynomial fit.
! this is performed during setup for model integration.
! WCS, 31 August 2009
!
! Described in:
! Skamarock, W. C., & Gassmann, A. (2011). 
! Conservative Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for ODE-Based Time Integration. 
! Monthly Weather Review, 139(9), 2962-2975. doi:10.1175/MWR-D-10-05056.1
!
      implicit none

      type (mesh_type), intent(in) :: grid
      integer, intent(out) :: err

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:), pointer :: advCells

!  local variables

      real (kind=RKIND), dimension(2, grid % nEdges) :: thetae
      real (kind=RKIND), dimension(grid % nEdges) :: xe, ye
      real (kind=RKIND), dimension(grid % nCells) :: theta_abs

      real (kind=RKIND), dimension(25) :: xc, yc, zc ! cell center coordinates
      real (kind=RKIND), dimension(25) :: thetav, thetat, dl_sphere
      real (kind=RKIND) :: xm, ym, zm, dl, xec, yec, zec
      real (kind=RKIND) :: thetae_tmp, xe_tmp, ye_tmp
      real (kind=RKIND) :: xv1, xv2, yv1, yv2, zv1, zv2
      integer :: i, j, k, ip1, ip2, m, n, ip1a, ii
      integer :: iCell, iEdge
      real (kind=RKIND) :: pii
      real (kind=RKIND) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
      real (kind=RKIND) :: pdx1, pdx2, pdx3, pdy1, pdy2, pdy3, dx1, dx2, dy1, dy2
      real (kind=RKIND) :: angv1, angv2, dl1, dl2
      real (kind=RKIND), dimension(25) :: dxe, dye, x2v, y2v, xp, yp
      
      real (kind=RKIND) :: amatrix(25,25), bmatrix(25,25), wmatrix(25,25)
      real (kind=RKIND) :: length_scale
      integer :: ma,na, cell_add, mw, nn
      integer, dimension(25) :: cell_list


      integer :: cell1, cell2
      integer, parameter :: polynomial_order = 2
      logical, parameter :: debug = .false.
      logical, parameter :: least_squares = .true.
      logical :: add_the_cell, do_the_cell

      logical, parameter :: reset_poly = .true.

      real (kind=RKIND) :: rcell, cos2t, costsint, sin2t
      real (kind=RKIND), dimension(grid%maxEdges) :: angle_2d

!---    
      err = 0

      if(polynomial_order > 2) then
        write (*,*) 'Polynomial for second derivitave can only be 2'
        err = 1
        return
      end if

      pii = 2.*asin(1.0)

!     advCells => grid % advCells % array
      allocate(advCells(grid % maxEdges2))
      deriv_two => grid % deriv_two % array
      deriv_two(:,:,:) = 0.

      do iCell = 1, grid % nCells !  is this correct? - we need first halo cell also...

         cell_list(1) = iCell
         do i=2, grid % nEdgesOnCell % array(iCell)+1
            cell_list(i) = grid % CellsOnCell % array(i-1,iCell)
         end do
         n = grid % nEdgesOnCell % array(iCell) + 1

         if ( polynomial_order > 2 ) then
            do i=2,grid % nEdgesOnCell % array(iCell) + 1
               do j=1,grid % nEdgesOnCell % array ( cell_list(i) )
                  cell_add = grid % CellsOnCell % array (j,cell_list(i))
                  add_the_cell = .true.
                  do k=1,n
                     if ( cell_add == cell_list(k) ) add_the_cell = .false.
                  end do
                  if (add_the_cell) then
                     n = n+1
                     cell_list(n) = cell_add
                  end if
               end do
            end do
         end if
 
         advCells(1) = n

!  check to see if we are reaching outside the halo

         do_the_cell = .true.
         do i=1,n
            if (cell_list(i) > grid % nCells) do_the_cell = .false.
         end do


         if ( .not. do_the_cell ) cycle


!  compute poynomial fit for this cell if all needed neighbors exist
         if ( grid % on_a_sphere ) then

            do i=1,n
               advCells(i+1) = cell_list(i)
               xc(i) = grid % xCell % array(advCells(i+1))/grid % sphere_radius
               yc(i) = grid % yCell % array(advCells(i+1))/grid % sphere_radius
               zc(i) = grid % zCell % array(advCells(i+1))/grid % sphere_radius
            end do

            if ( zc(1) == 1.0_RKIND) then
               theta_abs(iCell) = pii/2.0_RKIND
            else
               theta_abs(iCell) =  pii/2.0_RKIND - sphere_angle( xc(1), yc(1), zc(1),  &
                                                                 xc(2), yc(2), zc(2),  &
                                                                 0.0_RKIND, 0.0_RKIND, 1.0_RKIND ) 
            end if

! angles from cell center to neighbor centers (thetav)

            do i=1,n-1
   
               ip2 = i+2
               if (ip2 > n) ip2 = 2
    
               thetav(i) = sphere_angle( xc(1),   yc(1),   zc(1),    &
                                         xc(i+1), yc(i+1), zc(i+1),  &
                                         xc(ip2), yc(ip2), zc(ip2)   )

               dl_sphere(i) = grid % sphere_radius*arc_length( xc(1),   yc(1),   zc(1),  &
                                            xc(i+1), yc(i+1), zc(i+1) )
            end do

            length_scale = 1.
            do i=1,n-1
               dl_sphere(i) = dl_sphere(i)/length_scale
            end do

!            thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
            thetat(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
            do i=2,n-1
               thetat(i) = thetat(i-1) + thetav(i-1)
            end do
   
            do i=1,n-1
               xp(i) = cos(thetat(i)) * dl_sphere(i)
               yp(i) = sin(thetat(i)) * dl_sphere(i)
            end do

         else     ! On an x-y plane

            do i=1,n-1

               angle_2d(i) = grid%angleEdge%array(grid % EdgesOnCell % array(i,iCell))
               iEdge = grid % EdgesOnCell % array(i,iCell)
               if ( iCell .ne. grid % CellsOnEdge % array(1,iEdge)) &
                  angle_2d(i) = angle_2d(i) - pii

!              xp(i) = grid % xCell % array(cell_list(i)) - grid % xCell % array(iCell)
!              yp(i) = grid % yCell % array(cell_list(i)) - grid % yCell % array(iCell)

               xp(i) = grid % dcEdge % array(grid % EdgesOnCell % array(i,iCell)) * cos(angle_2d(i))
               yp(i) = grid % dcEdge % array(grid % EdgesOnCell % array(i,iCell)) * sin(angle_2d(i))

            end do

         end if


         ma = n-1
         mw = grid % nEdgesOnCell % array (iCell)

         bmatrix = 0.
         amatrix = 0.
         wmatrix = 0.

         if (polynomial_order == 2) then
            na = 6
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               wmatrix(i,i) = 1.
            end do
 
         else if (polynomial_order == 3) then
            na = 10
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
   
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               amatrix(i,7) = xp(i-1)**3
               amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
               amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
               amatrix(i,10) = yp(i-1)**3
   
               wmatrix(i,i) = 1.
 
            end do

         else
            na = 15
            ma = ma+1
  
            amatrix(1,1) = 1.
            wmatrix(1,1) = 1.
            do i=2,ma
               amatrix(i,1) = 1.
               amatrix(i,2) = xp(i-1)
               amatrix(i,3) = yp(i-1)
   
               amatrix(i,4) = xp(i-1)**2
               amatrix(i,5) = xp(i-1) * yp(i-1)
               amatrix(i,6) = yp(i-1)**2
   
               amatrix(i,7) = xp(i-1)**3
               amatrix(i,8) = yp(i-1) * (xp(i-1)**2)
               amatrix(i,9) = xp(i-1) * (yp(i-1)**2)
               amatrix(i,10) = yp(i-1)**3
   
               amatrix(i,11) = xp(i-1)**4
               amatrix(i,12) = yp(i-1) * (xp(i-1)**3)
               amatrix(i,13) = (xp(i-1)**2)*(yp(i-1)**2)
               amatrix(i,14) = xp(i-1) * (yp(i-1)**3)
               amatrix(i,15) = yp(i-1)**4
   
               wmatrix(i,i) = 1.
  
            end do
 
            do i=1,mw
               wmatrix(i,i) = 1.
            end do
 
         end if
 
         call ocn_poly_fit_2( amatrix, bmatrix, wmatrix, ma, na, 25 )

         do i=1,grid % nEdgesOnCell % array (iCell)
            ip1 = i+1
            if (ip1 > n-1) ip1 = 1
  
            iEdge = grid % EdgesOnCell % array (i,iCell)

            if(grid % on_a_sphere) then
              xv1 = grid % xVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              yv1 = grid % yVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              zv1 = grid % zVertex % array(grid % verticesOnEdge % array (1,iedge))/grid % sphere_radius
              xv2 = grid % xVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
              yv2 = grid % yVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
              zv2 = grid % zVertex % array(grid % verticesOnEdge % array (2,iedge))/grid % sphere_radius
            else
              xv1 = grid % xVertex % array(grid % verticesOnEdge % array (1,iedge))
              yv1 = grid % yVertex % array(grid % verticesOnEdge % array (1,iedge))
              zv1 = grid % zVertex % array(grid % verticesOnEdge % array (1,iedge))
              xv2 = grid % xVertex % array(grid % verticesOnEdge % array (2,iedge))
              yv2 = grid % yVertex % array(grid % verticesOnEdge % array (2,iedge))
              zv2 = grid % zVertex % array(grid % verticesOnEdge % array (2,iedge))
            end if
  
            if ( grid % on_a_sphere ) then
               call ocn_arc_bisect( xv1, yv1, zv1,  &
                                xv2, yv2, zv2,  &
                                xec, yec, zec   )
  
               thetae_tmp = sphere_angle( xc(1),   yc(1),   zc(1),    &
                                          xc(i+1), yc(i+1), zc(i+1),  &
                                          xec,     yec,     zec       )
               thetae_tmp = thetae_tmp + thetat(i)
               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
                  thetae(1,grid % EdgesOnCell % array (i,iCell)) = thetae_tmp
               else
                  thetae(2,grid % EdgesOnCell % array (i,iCell)) = thetae_tmp
               end if
!            else
!
!               xe(grid % EdgesOnCell % array (i,iCell)) = 0.5 * (xv1 + xv2)
!               ye(grid % EdgesOnCell % array (i,iCell)) = 0.5 * (yv1 + yv2)

            end if
  
         end do

!  fill second derivative stencil for rk advection 

         do i=1, grid % nEdgesOnCell % array (iCell)
            iEdge = grid % EdgesOnCell % array (i,iCell)
  
  
            if ( grid % on_a_sphere ) then
               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
  
                  cos2t = cos(thetae(1,grid % EdgesOnCell % array (i,iCell)))
                  sin2t = sin(thetae(1,grid % EdgesOnCell % array (i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
   
                  do j=1,n
                     deriv_two(j,1,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               else
     
                  cos2t = cos(thetae(2,grid % EdgesOnCell % array (i,iCell)))
                  sin2t = sin(thetae(2,grid % EdgesOnCell % array (i,iCell)))
                  costsint = cos2t*sin2t
                  cos2t = cos2t**2
                  sin2t = sin2t**2
      
                  do j=1,n
                     deriv_two(j,2,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               end if

            else

               cos2t = cos(angle_2d(i))
               sin2t = sin(angle_2d(i))
               costsint = cos2t*sin2t
               cos2t = cos2t**2
               sin2t = sin2t**2

!               do j=1,n
!
!                  deriv_two(j,1,iEdge) =   2.*xe(iEdge)*xe(iEdge)*bmatrix(4,j)  &
!                                         + 2.*xe(iEdge)*ye(iEdge)*bmatrix(5,j)  &
!                                         + 2.*ye(iEdge)*ye(iEdge)*bmatrix(6,j)
!               end do

               if (iCell == grid % cellsOnEdge % array(1,iEdge)) then
                  do j=1,n
                     deriv_two(j,1,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               else
                  do j=1,n
                     deriv_two(j,2,iEdge) =   2.*cos2t*bmatrix(4,j)  &
                                            + 2.*costsint*bmatrix(5,j)  &
                                            + 2.*sin2t*bmatrix(6,j)
                  end do
               end if

            end if
         end do
 
      end do ! end of loop over cells

      if (debug) stop


!      write(0,*) ' check for deriv2 coefficients, iEdge 4 '
!
!      iEdge = 4
!      j = 1
!      iCell = grid % cellsOnEdge % array(1,iEdge)
!      write(0,*) ' j, icell, coef ',j,iCell,deriv_two(j,1,iEdge)
!      do j=2,7
!         write(0,*) ' j, icell, coef ',j,grid % CellsOnCell % array(j-1,iCell),deriv_two(j,1,iEdge)
!      end do
!
!      j = 1
!      iCell = grid % cellsOnEdge % array(2,iEdge)
!      write(0,*) ' j, icell, coef ',j,iCell,deriv_two(j,2,iEdge)
!      do j=2,7
!         write(0,*) ' j, icell, coef ',j,grid % CellsOnCell % array(j-1,iCell),deriv_two(j,2,iEdge)
!      end do
!      stop

   end subroutine ocn_initialize_advection_rk!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION SPHERE_ANGLE
   !
   ! Computes the angle between arcs AB and AC, given points A, B, and C
   ! Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{
   
      implicit none
   
      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
   
      real (kind=RKIND) :: a, b, c          ! Side lengths of spherical triangle ABC
   
      real (kind=RKIND) :: ABx, ABy, ABz    ! The components of the vector AB
      real (kind=RKIND) :: mAB              ! The magnitude of AB
      real (kind=RKIND) :: ACx, ACy, ACz    ! The components of the vector AC
      real (kind=RKIND) :: mAC              ! The magnitude of AC
   
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
         sphere_angle =  2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      else
         sphere_angle = -2.0 * asin(max(min(sin_angle,1.0_RKIND),-1.0_RKIND))
      end if
   
   end function sphere_angle!}}}
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION PLANE_ANGLE
   !
   ! Computes the angle between vectors AB and AC, given points A, B, and C, and
   !   a vector (u,v,w) normal to the plane.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function plane_angle(ax, ay, az, bx, by, bz, cx, cy, cz, u, v, w)!{{{
   
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
         plane_angle =  acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      else
         plane_angle = -acos(max(min(cos_angle,1.0_RKIND),-1.0_RKIND))
      end if
   
   end function plane_angle!}}}


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FUNCTION ARC_LENGTH
   !
   ! Returns the length of the great circle arc from A=(ax, ay, az) to 
   !    B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   !    same sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real (kind=RKIND) function arc_length(ax, ay, az, bx, by, bz)!{{{
   
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
      arc_length = r * 2.0 * asin(c/(2.0*r))

   end function arc_length!}}}
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTine ocn_arc_bisect
   !
   ! Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   !   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   !   surface of a sphere centered at the origin.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ocn_arc_bisect(ax, ay, az, bx, by, bz, cx, cy, cz)!{{{
   
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
         write(0,*) 'Error: arc_bisect: A and B are diametrically opposite'
      else
         d = sqrt(cx*cx + cy*cy + cz*cz)
         cx = r * cx / d
         cy = r * cy / d
         cz = r * cz / d
      end if
   
   end subroutine ocn_arc_bisect!}}}


   subroutine ocn_poly_fit_2(a_in,b_out,weights_in,m,n,ne)!{{{

      implicit none

      integer, intent(in) :: m,n,ne
      real (kind=RKIND), dimension(ne,ne), intent(in) :: a_in, weights_in
      real (kind=RKIND), dimension(ne,ne), intent(out) :: b_out
   
      ! local storage
   
      real (kind=RKIND), dimension(m,n)  :: a
      real (kind=RKIND), dimension(n,m)  :: b
      real (kind=RKIND), dimension(m,m)  :: w,wt,h
      real (kind=RKIND), dimension(n,m)  :: at, ath
      real (kind=RKIND), dimension(n,n)  :: ata, ata_inv, atha, atha_inv
      integer, dimension(n) :: indx
      integer :: i,j
   
      if ( (ne<n) .or. (ne<m) ) then
         write(6,*) ' error in poly_fit_2 inversion ',m,n,ne
         stop
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
!         call ocn_migs(a,n,b,indx)
!      else

         call ocn_migs(atha,n,atha_inv,indx)

         b = matmul(atha_inv,ath)

!         call ocn_migs(ata,n,ata_inv,indx)
!         b = matmul(ata_inv,at)
!      end if
      b_out(1:n,1:m) = b(1:n,1:m)

!     do i=1,n
!        write(6,*) ' i, indx ',i,indx(i)
!     end do
!
!     write(6,*) ' '

   end subroutine ocn_poly_fit_2!}}}


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
subroutine ocn_migs (a,n,x,indx)!{{{
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
  call ocn_elgs (a,n,indx)
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
end subroutine ocn_migs!}}}


subroutine ocn_elgs (a,n,indx)!{{{
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
end subroutine ocn_elgs!}}}

!-------------------------------------------------------------

   subroutine ocn_initialize_deformation_weights( grid )!{{{
                                      
!
! compute the cell coefficients for the deformation calculations
! WCS, 13 July 2010
!
      implicit none

      type (mesh_type), intent(in) :: grid
!      type (grid_meta), intent(in) :: grid

      real (kind=RKIND), dimension(:,:), pointer :: defc_a, defc_b
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell

!  local variables

      real (kind=RKIND), dimension(2, grid % nEdges) :: thetae
      real (kind=RKIND), dimension(grid % nEdges) :: xe, ye
      real (kind=RKIND), dimension(grid % nCells) :: theta_abs

      real (kind=RKIND), dimension(25) :: xc, yc, zc ! cell center coordinates
      real (kind=RKIND), dimension(25) :: thetav, thetat, dl_sphere
      real (kind=RKIND) :: xm, ym, zm, dl, xec, yec, zec
      real (kind=RKIND) :: thetae_tmp, xe_tmp, ye_tmp
      real (kind=RKIND) :: xv1, xv2, yv1, yv2, zv1, zv2
      integer :: i, j, k, ip1, ip2, m, n, ip1a, ii
      integer :: iCell, iEdge
      real (kind=RKIND) :: pii
      real (kind=RKIND) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
      real (kind=RKIND) :: pdx1, pdx2, pdx3, pdy1, pdy2, pdy3, dx1, dx2, dy1, dy2
      real (kind=RKIND) :: angv1, angv2, dl1, dl2
      real (kind=RKIND), dimension(25) :: dxe, dye, x2v, y2v, xp, yp, xpt, ypt
      
      real (kind=RKIND) :: length_scale
      integer :: ma,na, cell_add, mw, nn
      integer, dimension(25) :: cell_list

      integer :: cell1, cell2, iv
      logical :: do_the_cell
      real (kind=RKIND) :: area_cell, sint2, cost2, sint_cost, sumw1, sumw2, xptt, area_cellt

      logical, parameter :: debug = .false.

      if (debug) write(0,*) ' in def weight calc '

      defc_a => grid % defc_a % array
      defc_b => grid % defc_b % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgesOnCell => grid % edgesOnCell % array

      defc_a(:,:) = 0.
      defc_b(:,:) = 0.

      pii = 2.*asin(1.0)

      if (debug) write(0,*) ' beginning cell loop '

      do iCell = 1, grid % nCells

         if (debug) write(0,*) ' cell loop ', iCell

         cell_list(1) = iCell
         do i=2, grid % nEdgesOnCell % array(iCell)+1
            cell_list(i) = grid % CellsOnCell % array(i-1,iCell)
         end do
         n = grid % nEdgesOnCell % array(iCell) + 1

!  check to see if we are reaching outside the halo

         if (debug) write(0,*) ' points ', n

         do_the_cell = .true.
         do i=1,n
            if (cell_list(i) > grid % nCells) do_the_cell = .false.
         end do


         if (.not. do_the_cell) cycle


!  compute poynomial fit for this cell if all needed neighbors exist
         if (grid % on_a_sphere) then

            xc(1) = grid % xCell % array(iCell)/a
            yc(1) = grid % yCell % array(iCell)/a
            zc(1) = grid % zCell % array(iCell)/a


            do i=2,n
               iv = grid % verticesOnCell % array(i-1,iCell)
               xc(i) = grid % xVertex % array(iv)/a
               yc(i) = grid % yVertex % array(iv)/a
               zc(i) = grid % zVertex % array(iv)/a
            end do

            theta_abs(iCell) =  pii/2. - sphere_angle( xc(1), yc(1), zc(1),  &
                                                       xc(2), yc(2), zc(2),  &
                                                       0.0_RKIND, 0.0_RKIND, 1.0_RKIND ) 

! angles from cell center to neighbor centers (thetav)

            do i=1,n-1
   
               ip2 = i+2
               if (ip2 > n) ip2 = 2
    
               thetav(i) = sphere_angle( xc(1),   yc(1),   zc(1),    &
                                         xc(i+1), yc(i+1), zc(i+1),  &
                                         xc(ip2), yc(ip2), zc(ip2)   )

               dl_sphere(i) = a*arc_length( xc(1),   yc(1),   zc(1),  &
                                            xc(i+1), yc(i+1), zc(i+1) )
            end do

            length_scale = 1.
            do i=1,n-1
               dl_sphere(i) = dl_sphere(i)/length_scale
            end do

            thetat(1) = 0.  !  this defines the x direction, cell center 1 -> 
!            thetat(1) = theta_abs(iCell)  !  this defines the x direction, longitude line
            do i=2,n-1
               thetat(i) = thetat(i-1) + thetav(i-1)
            end do
   
            do i=1,n-1
               xp(i) = cos(thetat(i)) * dl_sphere(i)
               yp(i) = sin(thetat(i)) * dl_sphere(i)
            end do

         else     ! On an x-y plane

            xp(1) = grid % xCell % array(iCell)
            yp(1) = grid % yCell % array(iCell)


            do i=2,n
               iv = grid % verticesOnCell % array(i-1,iCell)
               xp(i) = grid % xVertex % array(iv)
               yp(i) = grid % yVertex % array(iv)
            end do

         end if

!         thetat(1) = 0.
         thetat(1) = theta_abs(iCell)
         do i=2,n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            thetat(i) = plane_angle( 0.0_RKIND, 0.0_RKIND, 0.0_RKIND,  &
                                     xp(i)-xp(i-1), yp(i)-yp(i-1), 0.0_RKIND,  &
                                     xp(ip1)-xp(i), yp(ip1)-yp(i), 0.0_RKIND,  &
                                     0.0_RKIND, 0.0_RKIND, 1.0_RKIND)
            thetat(i) = thetat(i) + thetat(i-1)
         end do

         area_cell = 0.
         area_cellt = 0.
         do i=1,n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            dl = sqrt((xp(ip1)-xp(i))**2 + (yp(ip1)-yp(i))**2)
            area_cell = area_cell + 0.25*(xp(i)+xp(ip1))*(yp(ip1)-yp(i)) - 0.25*(yp(i)+yp(ip1))*(xp(ip1)-xp(i))
            area_cellt = area_cellt + (0.25*(xp(i)+xp(ip1))*cos(thetat(i)) + 0.25*(yp(i)+yp(ip1))*sin(thetat(i)))*dl
         end do
         if (debug) write(0,*) ' area_cell, area_cellt ',area_cell, area_cellt,area_cell-area_cellt

         do i=1,n-1
            ip1 = i+1
            if (ip1 == n) ip1 = 1
            dl = sqrt((xp(ip1)-xp(i))**2 + (yp(ip1)-yp(i))**2)
            sint2 = (sin(thetat(i)))**2
            cost2 = (cos(thetat(i)))**2
            sint_cost = sin(thetat(i))*cos(thetat(i))
            defc_a(i,iCell) = dl*(cost2 - sint2)/area_cell
            defc_b(i,iCell) = dl*2.*sint_cost/area_cell
            if (cellsOnEdge(1,EdgesOnCell(i,iCell)) /= iCell) then
               defc_a(i,iCell) = - defc_a(i,iCell)
               defc_b(i,iCell) = - defc_b(i,iCell)
            end if
 
         end do

      end do

      if (debug) write(0,*) ' exiting def weight calc '

   end subroutine ocn_initialize_deformation_weights!}}}

end module ocn_advection
