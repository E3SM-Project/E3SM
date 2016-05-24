module RootBiophysMod

#include "shr_assert.h"
  !-------------------------------------------------------------------------------------- 
  ! DESCRIPTION:
  ! module contains subroutine for root biophysics
  !
  ! HISTORY
  ! created by Jinyun Tang, Mar 1st, 2014
  implicit none
  private
  public :: init_vegrootfr
  public :: init_rootprof
  integer, parameter :: zeng_2001_root = 0 !the zeng 2001 root profile function

  integer :: root_prof_method              !select the type of root profile parameterization   
  !-------------------------------------------------------------------------------------- 

contains

  !-------------------------------------------------------------------------------------- 
  subroutine init_rootprof()
    !
    !DESCRIPTION
    ! initialize methods for root profile calculation
    implicit none

    root_prof_method = zeng_2001_root

  end subroutine init_rootprof

  !-------------------------------------------------------------------------------------- 
  subroutine init_vegrootfr(bounds, nlevsoi, nlevgrnd, rootfr)
    !
    !DESCRIPTION
    !initialize plant root profiles
    !
    ! USES
    use shr_kind_mod   , only : r8 => shr_kind_r8   
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use decompMod      , only : bounds_type
    use abortutils     , only : endrun         
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                     ! bounds
    integer,           intent(in) :: nlevsoi                    ! number of hydactive layers
    integer,           intent(in) :: nlevgrnd                   ! number of soil layers
    real(r8),          intent(out):: rootfr(bounds%begp: , 1: ) !
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'init_vegrootfr'  ! subroutine name
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(rootfr) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    select case (root_prof_method)
    case (zeng_2001_root)
       rootfr(bounds%begp:bounds%endp, 1 : nlevsoi) = zeng2001_rootfr(bounds, nlevsoi)

       !case (jackson_1996_root)
       !jackson root, 1996, to be defined later
       !rootfr(bounds%begp:bounds%endp, 1 : ubj) = jackson1996_rootfr(bounds, ubj, pcolumn, ivt, zi)
       !case (schenk_jackson_2002_root)
       !schenk and Jackson root, 2002, to be defined later
       !rootfr(bounds%begp:bounds%endp, 1 : ubj) = schenk2002_rootfr(bounds, ubj, pcolumn, ivt, zi)        
    case default
       call endrun(subname // ':: a root fraction function must be specified!')   
    end select
    rootfr(bounds%begp:bounds%endp,nlevsoi+1:nlevgrnd)=0._r8   

  end subroutine init_vegrootfr

  !--------------------------------------------------------------------------------------   
  function zeng2001_rootfr(bounds, ubj) result(rootfr)
    !
    ! DESCRIPTION
    ! compute root profile for soil water uptake
    ! using equation from Zeng 2001, J. Hydrometeorology
    !
    ! USES
    use shr_kind_mod   , only : r8 => shr_kind_r8   
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg   
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, roota_par, rootb_par  !these pars shall be moved to here and set as private in the future
    use PatchType      , only : pft
    use ColumnType     , only : col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                  ! bounds
    integer           , intent(in)    :: ubj                     ! ubnd
    !
    ! !RESULT
    real(r8) :: rootfr(bounds%begp:bounds%endp , 1:ubj ) !
    !
    ! !LOCAL VARIABLES:
    integer :: p, lev, c
    !------------------------------------------------------------------------

    !(computing from surface, d is depth in meter):
    ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
    ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
    ! beta & d_obs given in Zeng et al. (1998).   

    do p = bounds%begp,bounds%endp   

       if (pft%itype(p) /= noveg) then
          c = pft%column(p)
          do lev = 1, ubj-1
             rootfr(p,lev) = .5_r8*( exp(-roota_par(pft%itype(p)) * col%zi(c,lev-1))  &
                  + exp(-rootb_par(pft%itype(p)) * col%zi(c,lev-1))  &
                  - exp(-roota_par(pft%itype(p)) * col%zi(c,lev  ))  &
                  - exp(-rootb_par(pft%itype(p)) * col%zi(c,lev  )) )
          end do
          rootfr(p,ubj) = .5_r8*( exp(-roota_par(pft%itype(p)) * col%zi(c,ubj-1))  &
               + exp(-rootb_par(pft%itype(p)) * col%zi(c,ubj-1)) )

       else
          rootfr(p,1:ubj) = 0._r8
       endif

    enddo
    return

  end function zeng2001_rootfr

end module RootBiophysMod
