
      module mo_lymana

      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

      integer, parameter :: nla = 2

      contains

      subroutine lymana( o2col, secchi, dto2la, xso2la )
!-----------------------------------------------------------------------------
!   purpose:
!   calculate the effective absorption cross section of o2 in the lyman-alpha
!   bands and an effective o2 optical depth at all altitudes.  parameterized
!   after:  chabrillat, s., and g. kockarts, simple parameterization of the
!   absorption of the solar lyman-alpha line, geophysical research letters,
!   vol.24, no.21, pp 2659-2662, 1997.
!-----------------------------------------------------------------------------
!   parameters:
!   nz      - integer, number of specified altitude levels in the working (i)
!             grid
!   o2col   - real, slant overhead o2 column (molec/cc) at each specified (i)
!             altitude
!   dto2la  - real, optical depth due to o2 absorption at each specified  (o)
!             vertical layer
!   xso2la  - real, molecular absorption cross section in la bands        (o)
!-----------------------------------------------------------------------------

      use mo_params
      use ppgrid, only: pver, pverp

      implicit none

!-----------------------------------------------------------------------------
!	... dummy arguments
!-----------------------------------------------------------------------------
      real(r8), intent(in)  :: o2col(pverp)
      real(r8), intent(in)  :: secchi(pverp)
      real(r8), intent(out) :: dto2la(pver,nla-1)
      real(r8), intent(out) :: xso2la(pverp,nla-1)

!-----------------------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------------------
      integer :: i, k, kp1
      real(r8), dimension(pverp) :: rm, ro2
      real(r8), save :: b(3), c(3), d(3), e(3)

      data b / 6.8431e-01_r8, 2.29841e-01_r8,  8.65412e-02_r8/, &
           c /8.22114e-21_r8, 1.77556e-20_r8,  8.22112e-21_r8/, &
           d / 6.0073e-21_r8, 4.28569e-21_r8,  1.28059e-20_r8/, &
           e /8.21666e-21_r8, 1.63296e-20_r8,  4.85121e-17_r8/

!-----------------------------------------------------------------------------
! 	... calculate reduction factors at every altitude
!-----------------------------------------------------------------------------
      rm(:)  = 0._r8
      ro2(:) = 0._r8
      do k = 1,pverp
        do i = 1,3
          rm(k)  = rm(k) + b(i) * exp( -c(i)*o2col(k) )
          ro2(k) = ro2(k) + d(i) * exp( -e(i)*o2col(k) )
        end do
      end do

!-----------------------------------------------------------------------------
! 	... calculate effective o2 optical depths and effective o2 cross sections
!-----------------------------------------------------------------------------
      do k = 1,pver
         if( rm(k) > 1.e-100_r8 ) then
            kp1 = k + 1
            if( rm(kp1) > 0._r8 ) then
               dto2la(k,1) = log( rm(kp1) )/secchi(kp1) - log( rm(k) )/secchi(k)
            else
               dto2la(k,1) = 1000._r8
            end if
         else
            dto2la(k,1) = 1000._r8
         end if
      end do
      do k = 1,pverp
         if( rm(k) > 1.e-100_r8 ) then
            if( ro2(k) > 1.e-100_r8 ) then
               xso2la(k,1) = ro2(k)/rm(k)
            else
               xso2la(k,1) = 0._r8          
            end if
         else
            xso2la(k,1) = 0._r8
         end if
      end do

      end subroutine lymana

      end module mo_lymana
