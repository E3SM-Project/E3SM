module AnnualUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for updating annual summation variables
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type
  use CNStateType      , only : cnstate_type
  use ColumnDataType   , only : col_cf
  use VegetationDataType, only : veg_cf
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: AnnualUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine AnnualUpdate(bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, &
       cnstate_vars, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update annual summation variables
    !
    ! !USES:
    !#py use clm_time_manager, only: get_step_size, get_days_per_year
      !$acc routine seq
    use clm_varcon      , only: secspday
    use SubgridAveMod   , only: p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer               , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer               , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer               , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(cnstate_type)    , intent(inout) :: cnstate_vars
    real(r8) ,intent(in)  :: dt           ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: fp,fc        ! lake filter indices
    real(r8) :: dayspyr = 365.d0
    !-----------------------------------------------------------------------
    associate(&
      annsum_counter_col          => cnstate_vars%annsum_counter_col , &
      annsum_potential_gpp_patch  => cnstate_vars%annsum_potential_gpp_patch , &
      tempsum_potential_gpp_patch => cnstate_vars%tempsum_potential_gpp_patch , &
      annmax_retransn_patch       => cnstate_vars%annmax_retransn_patch  , &
      tempmax_retransn_patch      => cnstate_vars%tempmax_retransn_patch  , &
      annmax_retransp_patch       => cnstate_vars%annmax_retransp_patch  , &
      tempmax_retransp_patch      => cnstate_vars%tempmax_retransp_patch  , &
      annavg_t2m_patch            => cnstate_vars%annavg_t2m_patch    , &
      tempavg_t2m_patch           => cnstate_vars%tempavg_t2m_patch    , &
      annsum_npp                  => veg_cf%annsum_npp   , &
      tempsum_npp                 => veg_cf%tempsum_npp  , &
      annsum_npp_col              => col_cf%annsum_npp      , &
      annavg_t2m_col              => cnstate_vars%annavg_t2m_col &
      )
    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       annsum_counter_col(c) = annsum_counter_col(c) + dt
    end do

    if (num_soilc > 0) then
      if (annsum_counter_col(filter_soilc(1)) >= dayspyr * secspday) then
          ! patch loop
          do fp = 1,num_soilp
             p = filter_soilp(fp)

             ! update annual plant ndemand accumulator
             annsum_potential_gpp_patch(p)  = tempsum_potential_gpp_patch(p)
             tempsum_potential_gpp_patch(p) = 0._r8

             ! update annual total N retranslocation accumulator
             annmax_retransn_patch(p)  = tempmax_retransn_patch(p)
             tempmax_retransn_patch(p) = 0._r8

             ! update annual total P retranslocation accumulator
             annmax_retransp_patch(p)  = tempmax_retransp_patch(p)
             tempmax_retransp_patch(p) = 0._r8

             ! update annual average 2m air temperature accumulator
             annavg_t2m_patch(p)  = tempavg_t2m_patch(p)
             tempavg_t2m_patch(p) = 0._r8

             ! update annual NPP accumulator, convert to annual total
             annsum_npp(p) = tempsum_npp(p) * dt
             tempsum_npp(p) = 0._r8

          end do
          ! use p2c routine to get selected column-average pft-level fluxes and states

          call p2c(bounds, num_soilc, filter_soilc, &
               annsum_npp(bounds%begp:bounds%endp), &
               annsum_npp_col(bounds%begc:bounds%endc))

          call p2c(bounds, num_soilc, filter_soilc, &
               annavg_t2m_patch(bounds%begp:bounds%endp), &
               annavg_t2m_col(bounds%begc:bounds%endc))
       end if

    end if

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       if (annsum_counter_col(c) >= dayspyr * secspday) then
           annsum_counter_col(c) = 0._r8
       end if
    end do

  end associate

 end subroutine AnnualUpdate

end module AnnualUpdateMod
