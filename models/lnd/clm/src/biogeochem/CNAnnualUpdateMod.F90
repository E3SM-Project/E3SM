module CNAnnualUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for updating annual summation variables
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : bounds_type
  use CNVegCarbonFluxType , only : cnveg_carbonflux_type
  use CNvegStateType      , only : cnveg_state_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNAnnualUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNAnnualUpdate(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_state_inst, cnveg_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update annual summation variables
    !
    ! !USES:
    use clm_time_manager, only: get_step_size, get_days_per_year
    use clm_varcon      , only: secspday
    use SubgridAveMod   , only: p2c
    !
    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds  
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(cnveg_state_type)      , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type) , intent(inout) :: cnveg_carbonflux_inst
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
       cnveg_state_inst%annsum_counter_col(c) = cnveg_state_inst%annsum_counter_col(c) + dt
    end do

    if (num_soilc > 0) then

       if (cnveg_state_inst%annsum_counter_col(filter_soilc(1)) >= get_days_per_year() * secspday) then

          ! patch loop
          do fp = 1,num_soilp
             p = filter_soilp(fp)

             ! update annual plant ndemand accumulator
             cnveg_state_inst%annsum_potential_gpp_patch(p)  = cnveg_state_inst%tempsum_potential_gpp_patch(p)
             cnveg_state_inst%tempsum_potential_gpp_patch(p) = 0._r8

             ! update annual total N retranslocation accumulator
             cnveg_state_inst%annmax_retransn_patch(p)  = cnveg_state_inst%tempmax_retransn_patch(p)
             cnveg_state_inst%tempmax_retransn_patch(p) = 0._r8

             ! update annual average 2m air temperature accumulator
             cnveg_state_inst%annavg_t2m_patch(p)  = cnveg_state_inst%tempavg_t2m_patch(p)
             cnveg_state_inst%tempavg_t2m_patch(p) = 0._r8

             ! update annual NPP accumulator, convert to annual total
             cnveg_carbonflux_inst%annsum_npp_patch(p) = cnveg_carbonflux_inst%tempsum_npp_patch(p) * dt
             cnveg_carbonflux_inst%tempsum_npp_patch(p) = 0._r8

             ! update annual litfall accumulator, convert to annual total
             cnveg_carbonflux_inst%annsum_litfall_patch(p) = cnveg_carbonflux_inst%tempsum_litfall_patch(p) * dt
             cnveg_carbonflux_inst%tempsum_litfall_patch(p) = 0._r8
          end do

          ! use p2c routine to get selected column-average patch-level fluxes and states

          call p2c(bounds, num_soilc, filter_soilc, &
               cnveg_carbonflux_inst%annsum_npp_patch(bounds%begp:bounds%endp), &
               cnveg_carbonflux_inst%annsum_npp_col(bounds%begc:bounds%endc))

          call p2c(bounds, num_soilc, filter_soilc, &
               cnveg_state_inst%annavg_t2m_patch(bounds%begp:bounds%endp), &
               cnveg_state_inst%annavg_t2m_col(bounds%begc:bounds%endc))
       end if

    end if

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       if (cnveg_state_inst%annsum_counter_col(c) >= get_days_per_year() * secspday) then
          cnveg_state_inst%annsum_counter_col(c) = 0._r8
       end if
    end do

 end subroutine CNAnnualUpdate

end module CNAnnualUpdateMod
