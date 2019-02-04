module AnnualUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for updating annual summation variables
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type
  use CNCarbonFluxType , only : carbonflux_type
  use CNStateType      , only : cnstate_type
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
       cnstate_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update annual summation variables
    !
    ! !USES:
    use clm_time_manager, only: get_step_size, get_days_per_year
    use clm_varcon      , only: secspday
    use SubgridAveMod   , only: p2c
    use clm_varctl      , only: use_cndv 
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    integer               , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer               , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer               , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer               , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(cnstate_type)    , intent(inout) :: cnstate_vars
    type(carbonflux_type) , intent(inout) :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: fp,fc        ! lake filter indices
    real(r8):: dt           ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       cnstate_vars%annsum_counter_col(c) = cnstate_vars%annsum_counter_col(c) + dt
    end do

    if (num_soilc > 0) then

       if (cnstate_vars%annsum_counter_col(filter_soilc(1)) >= get_days_per_year() * secspday) then

          ! patch loop
          do fp = 1,num_soilp
             p = filter_soilp(fp)

             ! update annual plant ndemand accumulator
             cnstate_vars%annsum_potential_gpp_patch(p)  = cnstate_vars%tempsum_potential_gpp_patch(p)
             cnstate_vars%tempsum_potential_gpp_patch(p) = 0._r8

             ! update annual total N retranslocation accumulator
             cnstate_vars%annmax_retransn_patch(p)  = cnstate_vars%tempmax_retransn_patch(p)
             cnstate_vars%tempmax_retransn_patch(p) = 0._r8

             ! update annual total P retranslocation accumulator
             cnstate_vars%annmax_retransp_patch(p)  = cnstate_vars%tempmax_retransp_patch(p)
             cnstate_vars%tempmax_retransp_patch(p) = 0._r8

             ! update annual average 2m air temperature accumulator
             cnstate_vars%annavg_t2m_patch(p)  = cnstate_vars%tempavg_t2m_patch(p)
             cnstate_vars%tempavg_t2m_patch(p) = 0._r8

             ! update annual NPP accumulator, convert to annual total
             carbonflux_vars%annsum_npp_patch(p) = carbonflux_vars%tempsum_npp_patch(p) * dt
             carbonflux_vars%tempsum_npp_patch(p) = 0._r8

             if (use_cndv) then
                ! update annual litfall accumulator, convert to annual total
                carbonflux_vars%annsum_litfall_patch(p) = carbonflux_vars%tempsum_litfall_patch(p) * dt
                carbonflux_vars%tempsum_litfall_patch(p) = 0._r8
             end if
          end do

          ! use p2c routine to get selected column-average pft-level fluxes and states

          call p2c(bounds, num_soilc, filter_soilc, &
               carbonflux_vars%annsum_npp_patch(bounds%begp:bounds%endp), &
               carbonflux_vars%annsum_npp_col(bounds%begc:bounds%endc))

          call p2c(bounds, num_soilc, filter_soilc, &
               cnstate_vars%annavg_t2m_patch(bounds%begp:bounds%endp), &
               cnstate_vars%annavg_t2m_col(bounds%begc:bounds%endc))
       end if

    end if

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       if (cnstate_vars%annsum_counter_col(c) >= get_days_per_year() * secspday) then
          cnstate_vars%annsum_counter_col(c) = 0._r8
       end if
    end do

 end subroutine AnnualUpdate

end module AnnualUpdateMod
