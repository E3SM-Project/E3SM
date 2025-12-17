MODULE shr_wave_mod

   use shr_kind_mod, only : R8 => shr_kind_r8

   implicit none
   public :: shr_calc_wave_freq

contains
   subroutine shr_calc_wave_freq(nfr,freq1,xfreq,freq,dfreq)
      !This subroutine calculates the wave frequencies as done in WW3

      integer, intent(in) :: nfr ! number of frequency bin to use (or wave numbers)
      real(R8), intent(in) :: freq1 ! first freq
      real(R8), intent(in) :: xfreq ! freq multiplication factor

      real(R8), dimension(:), intent(out) :: freq  ! FREQUENCIES     
      real(R8), dimension(:), intent(out) :: dfreq    !Frequency bandwidths     

      !local vars
      integer :: ik
      real(R8) ::  sigma
      real(R8) ::  sxfr

      sigma   = freq1 / xfreq**2 
      sxfr    = 0.5 * (xfreq-1./xfreq)

      DO ik=1, nfr+2
         sigma    = sigma * xfreq
         freq(ik)  = sigma
         dfreq(ik) = sigma * sxfr      
      ENDDO

   end subroutine shr_calc_wave_freq        

END MODULE shr_wave_mod
