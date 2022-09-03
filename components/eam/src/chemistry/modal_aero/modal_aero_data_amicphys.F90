!IMPORTANT NOTE: This is a *temporary file* for testing a configuration
!variables declared here should be initialized in this file. It is very important
! that public variables here be declared as protected at the least!!

!If we decide to make this file stay/permanent, logic for computing
!igas_soag, igas_soagzz, nufi, mode_aging_optaa(max_mode), npca, iaer_pom, iaer_soa
!and mw_gas(max_gas) should be added to this file


module modal_aero_data_amicphys
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use modal_aero_data, only:  ntot_amode, nsoag, nsoa, npoa, nbc
  implicit none
  private
  save

#if ( defined MODAL_AERO_4MODE_SOA_MOM)
 logical, parameter, public :: is_soa_vbs = .true.
#else
 logical, parameter, public :: is_soa_vbs = .false.
#endif


#if ( defined MODAL_AERO_3MODE )
  integer, parameter, public :: max_gas = nsoag + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 3
#elif ( defined MODAL_AERO_4MODE )
  integer, parameter, public :: max_gas = nsoag + 1
  ! the +3 in max_aer are dst, ncl, so4
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 3
#elif ( ( defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_4MODE_SOA_MOM ) && ( defined MOSAIC_SPECIES ) )
  integer, parameter, public :: max_gas = nsoa + 10 !4 !QZR to match ngas=11
  ! the +9 in max_aer are dst, ncl, so4, mom, nh4, no3, cl, ca, co3
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 9
#elif ( defined MODAL_AERO_4MODE_MOM || defined MODAL_AERO_4MODE_SOA_MOM )
  integer, parameter, public :: max_gas = nsoag + 1
  !integer, parameter, public, public :: max_gas = nsoag + 1
  ! the +4 in max_aer are dst, ncl, so4, mom
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 4
#elif ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
  integer, parameter, public :: max_gas = nsoag + 4
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 8
#elif ( defined MODAL_AERO_7MODE )
  integer, parameter, public :: max_gas = nsoag + 2
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_8MODE )
  integer, parameter, public :: max_gas = nsoag + 2
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 4
#elif ( defined MODAL_AERO_9MODE )
  integer, parameter, public :: max_gas = nsoag + 2
  integer, parameter, public :: max_aer = nsoa + npoa + nbc + 4 + 5
#endif

#if (( defined MODAL_AERO_8MODE ) || ( defined MODAL_AERO_4MODE ) || ( defined MODAL_AERO_4MODE_MOM )|| (defined MODAL_AERO_4MODE_SOA_MOM))
  integer, parameter, public :: ntot_amode_extd = ntot_amode
#else
  integer, parameter, public :: ntot_amode_extd = ntot_amode + 1
#endif

  integer, parameter :: max_mode_fresh = 1

  integer, parameter, public:: max_mode = ntot_amode_extd + max_mode_fresh


  integer, public:: igas_soag, igas_soagzz, nufi, mode_aging_optaa(max_mode), &
       npca, iaer_pom, iaer_soa
  real(r8), public :: mw_gas(max_gas)



end module modal_aero_data_amicphys
