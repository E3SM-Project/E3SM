      module mo_airmas

      private
      public :: airmas

      contains

      subroutine airmas( z, zen, dsdh, nid, cz, &
                         vcol, scol )
!-----------------------------------------------------------------------------
!   purpose:
!   calculate vertical and slant air columns, in spherical geometry, as a
!   function of altitude.
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
!   vcol    - real, output, vertical air column, molec cm-2, above level iz
!   scol    - real, output, slant air column in direction of sun, above iz
!             also in molec cm-2
!-----------------------------------------------------------------------------

      use mo_params, only : largest
      use ppgrid,    only : pverp, pver
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nid(0:pver)
      real(r8), intent(in)    :: z(pverp)
      real(r8), intent(in)    :: zen
      real(r8), intent(in)    :: dsdh(0:pver,pver)
      real(r8), intent(in)    :: cz(pverp)
      real(r8), intent(out)   :: vcol(pverp)
      real(r8), intent(out)   :: scol(pverp)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: id, j
      real(r8)    :: sum, ssum, vsum, ratio

!-----------------------------------------------------------------------------
! 	... calculate vertical and slant column from each level: work downward
!-----------------------------------------------------------------------------
      vsum = 0._r8
      ssum = 0._r8
      do id = 1,pver
         vsum = vsum + cz(pverp-id)
         vcol(pverp-id) = vsum
         sum = 0._r8
         if( nid(id) < 0 ) then
            sum = largest
         else
!-----------------------------------------------------------------------------
! 	... single pass layers:
!-----------------------------------------------------------------------------
            do j = 1,min( nid(id),id )
               sum = sum + cz(pverp-j)*dsdh(id,j)
            end do
!-----------------------------------------------------------------------------
! 	... double pass layers:
!-----------------------------------------------------------------------------
            do j = min( nid(id),id )+1,nid(id)
               sum = sum + 2._r8*cz(pverp-j)*dsdh(id,j)
            end do
         end if
         scol(pverp - id) = sum
      end do

!-----------------------------------------------------------------------------
! 	... special section to set scol(pverp)
!-----------------------------------------------------------------------------
      if( scol(pver-1) /= 0._r8 ) then
         ratio       = scol(pver)/scol(pver-1)
	 scol(pverp) = ratio * scol(pver)
      else
         scol(pverp) = 0._r8
      end if

      end subroutine airmas

      end module mo_airmas
