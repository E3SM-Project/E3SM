
      module mo_setozo

      use shr_kind_mod, only : r8 => shr_kind_r8

      contains

      subroutine setozo( z, nw, wl, tlay, dto3, &
                         to3, o3, airlev, o3top )
!-----------------------------------------------------------------------------
!   purpose:
!   set up an altitude profile of ozone, and corresponding absorption
!   optical depths.  subroutine includes a shape-conserving scaling method
!   that allows scaling of the entire profile to a given overhead ozone
!   column amount.
!-----------------------------------------------------------------------------
!   parameters:
!   dobnew - real, overhead ozone column amount (du) to which profile     (i)
!            should be scaled.  if dobnew < 0, no scaling is done
!   nz     - integer, number of specified altitude levels in the working  (i)
!            grid
!   z      - real, specified altitude working grid (km)
!   nw     - integer, number of specified intervals + 1 in working        (i)
!            wavelength grid
!   wl     - real, vector of lower limits of wavelength intervals in      (i)
!            working wavelength grid
!   xso3   - real, molecular absoprtion cross section (cm^2) of o3 at     (i)
!            each specified wavelength (wmo value at 273)
!   s226   - real, molecular absoprtion cross section (cm^2) of o3 at     (i)
!            each specified wavelength (value from molina and molina at 226k)
!   s263   - real, molecular absoprtion cross section (cm^2) of o3 at     (i)
!            each specified wavelength (value from molina and molina at 263k)
!   s298   - real, molecular absoprtion cross section (cm^2) of o3 at     (i)
!            each specified wavelength (value from molina and molina at 298k)
!   tlay   - real, temperature (k) at each specified altitude layer       (i)
!   dto3   - real, optical depth due to ozone absorption at each          (o)
!            specified altitude at each specified wavelength
!   to3    - real, totol ozone column                                     (o)
!-----------------------------------------------------------------------------

      use mo_params, only : kw
      use mo_waveo3, only : xso3, s226, s263, s298
      use ppgrid,    only : pver, pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nw
      real(r8), intent(in)    :: o3top
      real(r8), intent(in)    :: wl(kw)
      real(r8), intent(in)    :: z(pverp)
      real(r8), intent(in)    :: tlay(pver) 
      real(r8), intent(in)    :: airlev(pverp) 
      real(r8), intent(out)   :: dto3(pver,nw)     
      real(r8), intent(out)   :: to3(pverp)
      real(r8), intent(inout) :: o3(pverp)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      real(r8), parameter :: wave_lims(2) = (/ 240.5_r8, 350._r8 /)
      real(r8), parameter :: t_lims(3)    = (/ 226._r8, 263._r8, 298._r8 /)
      real(r8), parameter :: tfac1 = 1._r8/(t_lims(2) - t_lims(1))
      real(r8), parameter :: tfac2 = 1._r8/(t_lims(3) - t_lims(2))

      integer :: k, wn
      real(r8)    :: so3
      real(r8)    :: cz(pverp)

!-----------------------------------------------------------------------------
! 	... compute column increments
!-----------------------------------------------------------------------------
      o3(1:pverp) = o3(1:pverp)*airlev(1:pverp) 
      cz(1:pver)  = .5_r8*(o3(2:pverp) + o3(1:pver)) * 1.e5_r8 * (z(2:pverp) - z(1:pver))
      to3(pverp)  = o3top
      do k = pver,1,-1
         to3(k) = to3(k+1) + cz(k)
      end do

!-----------------------------------------------------------------------------
! 	... include exponential tail integral from infinity to 50 km,
!           fold tail integral into top layer
!           specify scale height near top of data.
!-----------------------------------------------------------------------------

      cz(pver) = cz(pver) + o3top
!-----------------------------------------------------------------------------
! 	... calculate ozone optical depth for each layer, with temperature 
!           correction.  output, dto3(kz,kw)
!-----------------------------------------------------------------------------
      do wn = 1,nw
         if( wl(wn) > wave_lims(1)  .and. wl(wn+1) < wave_lims(2) ) then
            do k = 1,pver
               if( tlay(k) < t_lims(2) ) then
                  so3 = s226(wn) + (s263(wn) - s226(wn)) * tfac1 * (tlay(k) - t_lims(1))
               else
                  so3 = s263(wn) + (s298(wn) - s263(wn)) * tfac2 * (tlay(k) - t_lims(2))
               end if
               dto3(k,wn) = cz(k)*so3
            end do
         else
            dto3(1:pver,wn) = cz(1:pver) * xso3(wn)
         end if
      end do

      end subroutine setozo

      end module mo_setozo

