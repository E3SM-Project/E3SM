module SoilBiogeochemNitrogenUptakeMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the nitrogen uptake profile
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilBiogeochemNitrogenUptake
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemNitrogenUptake(bounds, nlevdecomp, num_soilc, filter_soilc, &
       sminn_vr, dzsoi_decomp, nfixation_prof, nuptake_prof)
    !
    ! DESCRIPTION
    ! Calculate the nitrogen uptake profile
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds  
    integer           , intent(in)    :: nlevdecomp                                          ! number of vertical layers
    integer           , intent(in)    :: num_soilc                                           ! number of soil columns in filter
    integer           , intent(in)    :: filter_soilc(:)                                     ! filter for soil columns
    real(r8)          , intent(in)    :: sminn_vr(bounds%begc: , 1: )                        ! soil mineral nitrogen profile
    real(r8)          , intent(in)    :: dzsoi_decomp(1: )                                   ! layer thickness
    real(r8)          , intent(in)    :: nfixation_prof(bounds%begc: , 1: )                  ! nitrogen fixation profile
    real(r8)          , intent(inout) :: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp) ! nitrogen uptake profile
    !
    ! !LOCAL VARIABLES:
    integer :: fc, j, c      ! indices
    real(r8):: sminn_tot(bounds%begc:bounds%endc)  !vertically integrated mineral nitrogen
    !-----------------------------------------------------------------------
  
    SHR_ASSERT_ALL((ubound(dzsoi_decomp)   == (/nlevdecomp/))              , errMsg(__FILE__, __LINE__))   
    SHR_ASSERT_ALL((ubound(sminn_vr)       == (/bounds%endc, nlevdecomp/)) , errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(nfixation_prof) == (/bounds%endc, nlevdecomp/)) , errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(nuptake_prof)   == (/bounds%endc, nlevdecomp/)) , errMsg(__FILE__, __LINE__))

    ! init sminn_tot
    do fc=1,num_soilc
       c = filter_soilc(fc)
       sminn_tot(c) = 0.
    end do

    do j = 1, nlevdecomp
       do fc=1,num_soilc
          c = filter_soilc(fc)
          sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
       end do
    end do

    do j = 1, nlevdecomp
       do fc=1,num_soilc
          c = filter_soilc(fc)      
          if (sminn_tot(c)  >  0.) then
             nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
          else
             nuptake_prof(c,j) = nfixation_prof(c,j)
          endif

       end do
    end do

  end subroutine SoilBiogeochemNitrogenUptake
  
end module SoilBiogeochemNitrogenUptakeMod
