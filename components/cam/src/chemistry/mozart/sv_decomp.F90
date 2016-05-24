!------------------------------------------------------------------------- 
! purpose: singular value decomposition
!
! method:
! given a matrix a(1:m,1:n), with physical dimensions mp by np,
! this routine computes its singular value decomposition,
! the matrix u replaces a on output. the
! diagonal matrix of singular values w is output as a vector
! w(1:n). the matrix v (not the transpose v^t) is output as
! v(1:n,1:n).
!
! author: a. maute dec 2003      
! (* copyright (c) 1985 numerical recipes software -- svdcmp *!
! from numerical recipes 1986 pp. 60 or can be find on web-sites
!------------------------------------------------------------------------- 

      module sv_decomp

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      private
      public :: svdcmp
      public :: svbksb

      integer, parameter :: nmax = 1600

      contains

      subroutine svdcmp( a, m, n, mp, np, w, v )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      integer, intent(in)     :: m
      integer, intent(in)     :: n
      integer, intent(in)     :: mp
      integer, intent(in)     :: np
      real(r8), intent(inout) :: a(mp,np)
      real(r8), intent(out)   :: v(np,np)
      real(r8), intent(out)   :: w(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, its, j, k, l, nm
      real(r8) :: anorm
      real(r8) :: c
      real(r8) :: f
      real(r8) :: g
      real(r8) :: h
      real(r8) :: s
      real(r8) :: scale
      real(r8) :: x, y, z
      real(r8) :: rv1(nmax)
      logical  :: cnd1
      logical  :: cnd2

      g     = 0.0_r8
      scale = 0.0_r8
      anorm = 0.0_r8

loop1 : &
      do i = 1,n
        l = i + 1
        rv1(i) = scale*g
        g     = 0.0_r8
        s     = 0.0_r8
        scale = 0.0_r8
        if( i <= m ) then
          do k = i,m
            scale = scale + abs(a(k,i))
          end do
          if( scale /= 0.0_r8 ) then
            do k = i,m
              a(k,i) = a(k,i)/scale
              s = s + a(k,i)*a(k,i)
            end do
            f = a(i,i)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,i) = f - g
            if( i /= n ) then
              do j = l,n
                s = 0.0_r8
                do k = i,m
                  s = s + a(k,i)*a(k,j)
                end do
                f = s/h
                do k = i,m
                  a(k,j) = a(k,j) + f*a(k,i)
                end do
              end do
            end if
            do k = i,m
              a(k,i) = scale*a(k,i)
            end do
          endif
        endif
        w(i) = scale *g
        g     = 0.0_r8
        s     = 0.0_r8
        scale = 0.0_r8
        if( i <= m .and. i /= n ) then
          do k = l,n
            scale = scale + abs(a(i,k))
          end do
          if( scale /= 0.0_r8 ) then
            do k = l,n
              a(i,k) = a(i,k)/scale
              s      = s + a(i,k)*a(i,k)
            end do
            f = a(i,l)
            g = -sign(sqrt(s),f)
            h = f*g - s
            a(i,l) = f - g
            do k = l,n
              rv1(k) = a(i,k)/h
            end do
            if( i /= m ) then
              do j = l,m
                s = 0.0_r8
                do k = l,n
                  s = s + a(j,k)*a(i,k)
                end do
                do k = l,n
                  a(j,k) = a(j,k) + s*rv1(k)
                end do
              end do
            end if
            do k = l,n
              a(i,k) = scale*a(i,k)
            end do
          end if
        end if
        anorm = max( anorm,(abs(w(i)) + abs(rv1(i))) )
      end do loop1

      do i = n,1,-1
        if( i < n ) then
          if( g /= 0.0_r8 ) then
            do j = l,n
              v(j,i) = (a(i,j)/a(i,l))/g
            end do
            do j = l,n
              s = 0.0_r8
              do k = l,n
                s = s + a(i,k)*v(k,j)
              end do
              do k = l,n
                v(k,j) = v(k,j) + s*v(k,i)
              end do
            end do
          end if
          do j = l,n
            v(i,j) = 0.0_r8
            v(j,i) = 0.0_r8
          end do
        end if
        v(i,i) = 1.0_r8
        g = rv1(i)
        l = i
      end do

      do i = n,1,-1
        l = i + 1
        g = w(i)
        if( i < n ) then
          do j = l,n
            a(i,j) = 0.0_r8
          end do
        end if
        if( g /= 0.0_r8 ) then
          g = 1.0_r8/g
          if( i /= n ) then
            do j = l,n
              s = 0.0_r8
              do k = l,m
                s = s + a(k,i)*a(k,j)
              end do
              f = (s/a(i,i))*g
              do k = i,m
                a(k,j) = a(k,j) + f*a(k,i)
              end do
            end do
          end if
          do j = i,m
            a(j,i) = a(j,i)*g
          end do
        else
          do j = i,m
            a(j,i) = 0.0_r8
          end do
        end if
        a(i,i) = a(i,i) + 1.0_r8
      end do

      do k = n,1,-1
loop2 : do its = 1,30
          do l = k,1,-1
            nm = l - 1
            cnd1 = abs( rv1(l) ) + anorm == anorm
            if( cnd1 ) then
              cnd2 = .false.
              exit
            end if
            cnd2 = abs( w(nm) ) + anorm == anorm
            if( cnd2 ) then
              cnd1 = .true.
              exit
            else if( l == 1 ) then
              cnd1 = .true.
              cnd2 = .true.
            end if
          end do

          if( cnd2 ) then
            c = 0.0_r8
            s = 1.0_r8
            do i = l,k
              f = s*rv1(i)
              if( (abs(f) + anorm) /= anorm ) then
                g = w(i)
                h = sqrt(f*f + g*g)
                w(i) = h
                h = 1.0_r8/h
                c = (g*h)
                s = -(f*h)
                do j = 1,m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
                end do
              end if
            end do
          end if

          if( cnd1 ) then
            z = w(k)
            if( l == k ) then
              if( z < 0.0_r8 ) then
                w(k) = -z
                do j = 1,n
                  v(j,k) = -v(j,k)
                end do
              end if
              exit loop2
            end if
          end if

          x = w(l)
          nm = k - 1
          y = w(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0_r8*h*y)
          g = sqrt( f*f + 1.0_r8 )
          f = ((x - z)*(x + z) + h*((y/(f + sign(g,f))) - h))/x
          c = 1.0_r8
          s = 1.0_r8
          do j = l,nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = sqrt( f*f + h*h )
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do nm = 1,n
              x = v(nm,j)
              z = v(nm,i)
              v(nm,j) = (x*c)+(z*s)
              v(nm,i) = -(x*s)+(z*c)
            end do
            z = sqrt( f*f + h*h )
            w(j) = z
            if( z /= 0.0_r8 ) then
              z = 1.0_r8/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do nm = 1,m
              y = a(nm,j)
              z = a(nm,i)
              a(nm,j) = (y*c)+(z*s)
              a(nm,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0_r8
          rv1(k) = f
          w(k)   = x
        end do loop2
      end do
      
      end subroutine svdcmp

!-------------------------------------------------------------------------      
! purpose: solves a*x = b
!
! method:     
! solves a*x = b for a vector x, where a is specified by the arrays
! u,w,v as returned by svdcmp. m and n
! are the logical dimensions of a, and will be equal for square matrices.
! mp and np are the physical dimensions of a. b(1:m) is the input right-hand 
! side. x(1:n) is the output solution vector. no input quantities are 
! destroyed, so the routine may be called sequentially with different b
!
! author:  a. maute dec 2002   
! (* copyright (c) 1985 numerical recipes software -- svbksb *!
! from numerical recipes 1986 pp. 57 or can be find on web-sites
!-------------------------------------------------------------------------      

      subroutine svbksb( u, w, v, m, n, mp, np, b, x )
!------------------------------------------------------------------------- 
!	... dummy arguments
!------------------------------------------------------------------------- 
      integer, intent(in)   :: m
      integer, intent(in)   :: n
      integer, intent(in)   :: mp
      integer, intent(in)   :: np
      real(r8), intent(in)  :: u(mp,np)
      real(r8), intent(in)  :: w(np)
      real(r8), intent(in)  :: v(np,np)
      real(r8), intent(in)  :: b(mp)
      real(r8), intent(out) :: x(np)

!------------------------------------------------------------------------- 
!	... local variables
!------------------------------------------------------------------------- 
      integer  :: i, j, jj
      real(r8) :: s
      real(r8) :: tmp(nmax)

      do j = 1,n
        s = 0._r8
        if( w(j) /= 0._r8 ) then
          do i = 1,m
            s = s + u(i,j)*b(i)
          end do
          s = s/w(j)
        endif
        tmp(j) = s
      end do

      do j = 1,n
        s = 0._r8
        do jj = 1,n
          s = s + v(j,jj)*tmp(jj)
        end do
        x(j) = s
      end do

      end subroutine svbksb

      end module sv_decomp
