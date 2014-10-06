
module pft2colMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pft2colMod
!
! !DESCRIPTION:
! Contains calls to methods to perfom averages over from pfts to columns
! for model variables.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use subgridAveMod
  use clmtype
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c  ! obtain column properties from average over column pfts
!
! !REVISION HISTORY:
! 03/09/08: Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pft2col
!
! !INTERFACE:
  subroutine pft2col (lbc, ubc, num_nolakec, filter_nolakec)
!
! !DESCRIPTION:
! Averages over all pfts for variables defined over both soil and lake
! to provide the column-level averages of state and flux variables
! defined at the pft level.
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
!
! !REVISION HISTORY:
! 03/09/08: Created by Mariana Vertenstein
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer :: c,fc                      ! indices
    integer :: num_allc                  ! number of total column points
    integer :: filter_allc(ubc-lbc+1)    ! filter for all column points
    real(r8), pointer :: ptrp(:)         ! pointer to input pft array
    real(r8), pointer :: ptrc(:)         ! pointer to output column array
! -----------------------------------------------------------------

    ! Set up a filter for all column points

    num_allc = ubc-lbc+1
    fc = 0
    do c = lbc,ubc
       fc = fc + 1
       filter_allc(fc) = c
    end do

    ! Note: lake points are excluded from many of the following averages. For some fields,
    ! this is because the field doesn't apply over lakes. However, for many others, this
    ! is because the field is computed in HydrologyLake, which is called after this
    ! routine; thus, for lakes, the column-level values of these fields are explicitly set
    ! in HydrologyLakeMod. (The fields that are included here for lakes are computed
    ! elsewhere, e.g., in BiogeophysicsLake.)

    ! Averaging for pft water state variables

    ptrp => pws%h2ocan
    ptrc => pws_a%h2ocan
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ! Averaging for pft water flux variables

    ptrp => pwf%qflx_evap_tot
    ptrc => pwf_a%qflx_evap_tot
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_rain_grnd
    ptrc => pwf_a%qflx_rain_grnd
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_snow_grnd
    ptrc => pwf_a%qflx_snow_grnd
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)
    
    ptrp => pwf%qflx_snwcp_liq
    ptrc => pwf_a%qflx_snwcp_liq
    call p2c (num_allc, filter_allc, ptrp, ptrc)

    ptrp => pwf%qflx_snwcp_ice
    ptrc => pwf_a%qflx_snwcp_ice
    ! For lakes, this field is initially set in BiogeophysicsLake (which is called before
    ! this routine; hence it is appropriate to include lake columns in this p2c call).
    ! However, it is later overwritten in HydrologyLake, both on the pft and the column
    ! level.
    call p2c (num_allc, filter_allc, ptrp, ptrc)

    ptrp => pwf%qflx_tran_veg
    ptrc => pwf_a%qflx_tran_veg
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_evap_grnd
    ptrc => pwf_a%qflx_evap_grnd
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_evap_soi
    ptrc => pwf_a%qflx_evap_soi
    call p2c (num_allc, filter_allc, ptrp, ptrc)

    ptrp => pwf%qflx_prec_grnd
    ptrc => pwf_a%qflx_prec_grnd
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_dew_grnd
    ptrc => pwf_a%qflx_dew_grnd
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_sub_snow
    ptrc => pwf_a%qflx_sub_snow
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

    ptrp => pwf%qflx_dew_snow
    ptrc => pwf_a%qflx_dew_snow
    call p2c (num_nolakec, filter_nolakec, ptrp, ptrc)

  end subroutine pft2col

end module pft2colMod
