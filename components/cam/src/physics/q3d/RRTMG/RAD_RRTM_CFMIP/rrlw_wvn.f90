      module rrlw_wvn

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : nbndlw, mg, ngptlw, maxinpx

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: For the lower atmosphere, the number of reference
!                    atmospheres that are stored for each spectral band
!                    per pressure level and temperature.  Each of these
!                    atmospheres has different relative amounts of the 
!                    key species for the band (i.e. different binary
!                    species parameters).
! nspb   :  integer: Same as nspa for the upper atmosphere
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
! totplnk:  real   : Integrated Planck value for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!
! ngc    :  integer: The number of new g-intervals in each band
! ngs    :  integer: The cumulative sum of new g-intervals for each band
! ngm    :  integer: The index of each new g-interval relative to the
!                    original 16 g-intervals in each band
! ngn    :  integer: The number of original g-intervals that are 
!                    combined to make each new g-intervals in each band
! ngb    :  integer: The band index for each new g-interval
! wt     :  real   : RRTM weights for the original 16 g-intervals
! rwgt   :  real   : Weights for combining original 16 g-intervals 
!                    (256 total) into reduced set of g-intervals 
!                    (140 total)
! nxmol  :  integer: Number of cross-section molecules
! ixindx :  integer: Flag for active cross-sections in calculation
!------------------------------------------------------------------

      integer(kind=im) :: ng(nbndlw)
      integer(kind=im) :: nspa(nbndlw)
      integer(kind=im) :: nspb(nbndlw)

      real(kind=rb) :: wavenum1(nbndlw)
      real(kind=rb) :: wavenum2(nbndlw)
      real(kind=rb) :: delwave(nbndlw)

      real(kind=rb) :: totplnk(181,nbndlw)
      real(kind=rb) :: totplk16(181)

      integer(kind=im) :: ngc(nbndlw)
      integer(kind=im) :: ngs(nbndlw)
      integer(kind=im) :: ngn(ngptlw)
      integer(kind=im) :: ngb(ngptlw)
      integer(kind=im) :: ngm(nbndlw*mg)

      real(kind=rb) :: wt(mg)
      real(kind=rb) :: rwgt(nbndlw*mg)

      integer(kind=im) :: nxmol
      integer(kind=im) :: ixindx(maxinpx)

      end module rrlw_wvn
