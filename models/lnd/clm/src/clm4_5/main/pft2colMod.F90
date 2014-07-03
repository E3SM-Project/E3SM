module pft2colMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains calls to methods to perfom averages over from pfts to columns
  ! for model variables.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use subgridAveMod
  use clmtype
  use decompMod, only: bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c  ! obtain column properties from average over column pfts
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine pft2col (bounds, num_nolakec, filter_nolakec)
    !
    ! !DESCRIPTION:
    ! Averages over all pfts for variables defined over both soil and lake
    ! to provide the column-level averages of state and flux variables
    ! defined at the pft level.
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_nolakec       ! number of column non-lake points in column filter
    integer, intent(in) :: filter_nolakec(:) ! column filter for non-lake points
    !
    ! !LOCAL VARIABLES:
    integer :: c,fc              ! indices
    integer :: num_allc          ! number of active column points
    integer :: filter_allc(bounds%endp-bounds%begp+1)    ! filter for all active column points
    ! -----------------------------------------------------------------

    ! Set up a filter for all active column points

    fc = 0
    do c = bounds%begc,bounds%endc
       if (col%active(c)) then
          fc = fc + 1
          filter_allc(fc) = c
       end if
    end do
    num_allc = fc

    ! Note: lake points are excluded from many of the following averages. For some fields,
    ! this is because the field doesn't apply over lakes. However, for many others, this
    ! is because the field is computed in SLakeHydrologyMod, which is called after this
    ! routine; thus, for lakes, the column-level values of these fields are explicitly set
    ! in SLakeHydrologyMod. (The fields that are included here for lakes are computed
    ! elsewhere, e.g., in SLakeFluxesMod.)

    ! Averaging for pft water state variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pws%h2ocan(bounds%begp:bounds%endp), &
         pws_a%h2ocan(bounds%begc:bounds%endc))

    ! Averaging for pft water flux variables

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_ev_snow(bounds%begp:bounds%endp), &
         pwf_a%qflx_ev_snow(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_ev_soil(bounds%begp:bounds%endp), &
         pwf_a%qflx_ev_soil(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_ev_h2osfc(bounds%begp:bounds%endp), &
         pwf_a%qflx_ev_h2osfc(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_evap_soi(bounds%begp:bounds%endp), &
         pwf_a%qflx_evap_soi(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_evap_tot(bounds%begp:bounds%endp), &
         pwf_a%qflx_evap_tot(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_rain_grnd(bounds%begp:bounds%endp), &
         pwf_a%qflx_rain_grnd(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_snow_grnd(bounds%begp:bounds%endp), &
         pwf_a%qflx_snow_grnd(bounds%begc:bounds%endc))
    
    call p2c (bounds, num_allc, filter_allc, &
         pwf%qflx_snwcp_liq(bounds%begp:bounds%endp), &
         pwf_a%qflx_snwcp_liq(bounds%begc:bounds%endc))

    ! For lakes, this field is initially set in SLakeFluxesMod (which is called before this routine; 
    ! hence it is appropriate to include lake columns in this p2c call.
    ! However, it is later overwritten in SLakeHydrologyMod, both on the pft and the column level.

    call p2c (bounds, num_allc, filter_allc, &
         pwf%qflx_snwcp_ice(bounds%begp:bounds%endp), &
         pwf_a%qflx_snwcp_ice(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_tran_veg(bounds%begp:bounds%endp), &
         pwf_a%qflx_tran_veg(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_evap_grnd(bounds%begp:bounds%endp), &
         pwf_a%qflx_evap_grnd(bounds%begc:bounds%endc))

    call p2c (bounds, num_allc, filter_allc, &
         pwf%qflx_evap_soi(bounds%begp:bounds%endp), &
         pwf_a%qflx_evap_soi(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_prec_grnd(bounds%begp:bounds%endp), &
         pwf_a%qflx_prec_grnd(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_dew_grnd(bounds%begp:bounds%endp), &
         pwf_a%qflx_dew_grnd(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_sub_snow(bounds%begp:bounds%endp), &
         pwf_a%qflx_sub_snow(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_dew_snow(bounds%begp:bounds%endp), &
         pwf_a%qflx_dew_snow(bounds%begc:bounds%endc))

    call p2c (bounds, num_nolakec, filter_nolakec, &
         pwf%qflx_irrig(bounds%begp:bounds%endp), &
         pwf_a%qflx_irrig(bounds%begc:bounds%endc))

  end subroutine pft2col

end module pft2colMod
