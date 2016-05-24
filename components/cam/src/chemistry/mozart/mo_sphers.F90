
      module mo_sphers

      contains

      subroutine sphers( z, zen, dsdh, nid )
!-----------------------------------------------------------------------------
!   purpose:
!   calculate slant path over vertical depth ds/dh in spherical geometry.
!   calculation is based on:  a.dahlback, and k.stamnes, a new spheric model
!   for computing the radiation field available for photolysis and heating
!   at twilight, planet.space sci., v39, n5, pp. 671-683, 1991 (appendix b)
!-----------------------------------------------------------------------------
!   parameters:
!   nz      - integer, number of specified altitude levels in the working (i)
!             grid
!   z       - real, specified altitude working grid (km)                  (i)
!   zen     - real, solar zenith angle (degrees)                          (i)
!   dsdh    - real, slant path of direct beam through each layer crossed  (o)
!             when travelling from the top of the atmosphere to layer i;
!             dsdh(i,j), i = 0..nz-1, j = 1..nz-1
!   nid     - integer, number of layers crossed by the direct beam when   (o)
!             travelling from the top of the atmosphere to layer i;
!             nid(i), i = 0..nz-1
!-----------------------------------------------------------------------------

      use mo_params,    only : rearth
      use mo_constants, only : d2r
      use ppgrid,       only : pver, pverp
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)     :: zen
      real(r8), intent(in)     :: z(pverp)
      integer, intent(out) :: nid(0:pver)
      real(r8), intent(out)    :: dsdh(0:pver,pver)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: i, j, k
      integer :: id
      real(r8)    :: re, ze(pverp)
      real(r8)    :: zd(0:pver)
      real(r8)    :: zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
      real(r8)    :: radius

      radius = rearth*1.e-3_r8    ! rearth m -> km
      zenrad = zen*d2r
!-----------------------------------------------------------------------------
! 	... include the elevation above sea level to the radius of the earth:
!-----------------------------------------------------------------------------
      re = radius + z(1)
!-----------------------------------------------------------------------------
! 	... correspondingly z changed to the elevation above earth surface:
!-----------------------------------------------------------------------------
      ze(1:pverp) = z(1:pverp) - z(1)

!-----------------------------------------------------------------------------
! 	... inverse coordinate of z
!-----------------------------------------------------------------------------
      zd(0) = ze(pverp)
      do k = 1,pver
        zd(k) = ze(pverp - k)
      end do

!-----------------------------------------------------------------------------
! 	... initialize dsdh(i,j), nid(i)
!-----------------------------------------------------------------------------
      dsdh(0:pver,1:pver) = 0._r8
      nid(0:pver)         = 0

!-----------------------------------------------------------------------------
! 	... calculate ds/dh of every layer
!-----------------------------------------------------------------------------
      do i = 0,pver
        rpsinz = (re + zd(i)) * sin(zenrad)
        if( zen > 90._r8 .and. rpsinz < re ) then
           nid(i) = -1
        else
!-----------------------------------------------------------------------------
! 	... find index of layer in which the screening height lies
!-----------------------------------------------------------------------------
           id = i 
           if( zen > 90._r8 ) then
              do j = 1,pver
                 if( (rpsinz < ( zd(j-1) + re ) ) .and. (rpsinz >= ( zd(j) + re )) ) then
		    id = j
		 end if
              end do
           end if
 
           do j = 1,id
             if( j == id .and. id == i .and. zen > 90._r8) then
                sm = -1._r8
	     else
                sm = 1._r8
             end if 
             rj = re + zd(j-1)
             rjp1 = re + zd(j)
             dhj = zd(j-1) - zd(j)
             ga = max( 0._r8,rj*rj - rpsinz*rpsinz )
             gb = max( 0._r8,rjp1*rjp1 - rpsinz*rpsinz )
             if( id > i .and. j == id ) then
                dsj = sqrt( ga )
             else
                dsj = sqrt( ga ) - sm*sqrt( gb )
             end if
             dsdh(i,j) = dsj / dhj
           end do
           nid(i) = id
        end if
      end do

      end subroutine sphers

      end module mo_sphers
