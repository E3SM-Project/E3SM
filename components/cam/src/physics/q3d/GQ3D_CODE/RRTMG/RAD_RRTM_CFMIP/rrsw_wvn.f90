      module rrsw_wvn

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : nbndsw, mg, ngptsw, jpb1, jpb2

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: 
! nspb   :  integer: 
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
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
!                    (224 total) into reduced set of g-intervals 
!                    (112 total)
!------------------------------------------------------------------

      integer(kind=im) :: ng(jpb1:jpb2)
      integer(kind=im) :: nspa(jpb1:jpb2)
      integer(kind=im) :: nspb(jpb1:jpb2)

      real(kind=rb) :: wavenum1(jpb1:jpb2)
      real(kind=rb) :: wavenum2(jpb1:jpb2)
      real(kind=rb) :: delwave(jpb1:jpb2)

      integer(kind=im) :: ngc(nbndsw)
      integer(kind=im) :: ngs(nbndsw)
      integer(kind=im) :: ngn(ngptsw)
      integer(kind=im) :: ngb(ngptsw)
      integer(kind=im) :: ngm(nbndsw*mg)

      real(kind=rb) :: wt(mg)
      real(kind=rb) :: rwgt(nbndsw*mg)

      end module rrsw_wvn
