      module rrlw_ref

      use parkind, only : im => kind_im, rb => kind_rb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw reference atmosphere 
! Based on standard mid-latitude summer profile
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! pref   :  real   : Reference pressure levels
! preflog:  real   : Reference pressure levels, ln(pref)
! tref   :  real   : Reference temperature levels for MLS profile
! chi_mls:  real   : 
!------------------------------------------------------------------

      real(kind=rb) , dimension(59) :: pref
      real(kind=rb) , dimension(59) :: preflog
      real(kind=rb) , dimension(59) :: tref
      real(kind=rb) :: chi_mls(7,59)

      end module rrlw_ref
