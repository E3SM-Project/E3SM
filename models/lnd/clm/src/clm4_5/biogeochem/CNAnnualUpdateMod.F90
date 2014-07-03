module CNAnnualUpdateMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for updating annual summation variables
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use decompMod   , only: bounds_type
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNAnnualUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNAnnualUpdate(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update annual summation variables
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size, get_days_per_year
    use clm_varcon      , only: secspday
    use pft2colMod      , only: p2c
    use clm_varctl      , only: use_cndv 
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_soilc         ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
    integer, intent(in) :: num_soilp         ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: fp,fc        ! lake filter indices
    real(r8):: dt           ! radiation time step (seconds)
    !-----------------------------------------------------------------------

   associate(& 
   annsum_counter            => cps%annsum_counter         , & ! InOut:  [real(r8) (:)]  seconds since last annual accumulator turnover
   tempsum_potential_gpp     => pepv%tempsum_potential_gpp , & ! InOut:  [real(r8) (:)]  temporary annual sum of potential GPP   
   annsum_potential_gpp      => pepv%annsum_potential_gpp  , & ! InOut:  [real(r8) (:)]  annual sum of potential GPP             
   tempmax_retransn          => pepv%tempmax_retransn      , & ! InOut:  [real(r8) (:)]  temporary annual max of retranslocated N pool (gN/m2)
   annmax_retransn           => pepv%annmax_retransn       , & ! InOut:  [real(r8) (:)]  annual max of retranslocated N pool (gN/m2)
   tempavg_t2m               => pepv%tempavg_t2m           , & ! InOut:  [real(r8) (:)]  temporary average 2m air temperature (K)
   annavg_t2m                => pepv%annavg_t2m            , & ! InOut:  [real(r8) (:)]  annual average 2m air temperature (K)   
   tempsum_npp               => pepv%tempsum_npp           , & ! InOut:  [real(r8) (:)]  temporary sum NPP (gC/m2/yr)            
   annsum_npp                => pepv%annsum_npp            , & ! InOut:  [real(r8) (:)]  annual sum NPP (gC/m2/yr)               
   cannsum_npp               => cps%cannsum_npp            , & ! InOut:  [real(r8) (:)]  column annual sum NPP (gC/m2/yr)        
   cannavg_t2m               => cps%cannavg_t2m            , & ! InOut:  [real(r8) (:)] annual average of 2m air temperature, averaged from pft-level (K)
   tempsum_litfall           => pepv%tempsum_litfall       , & ! InOut:  [real(r8) (:)]  temporary sum litfall (gC/m2/yr)        
   annsum_litfall            => pepv%annsum_litfall        , & ! InOut:  [real(r8) (:)]  annual sum litfall (gC/m2/yr)           
   pcolumn                   => pft%column                   & ! Input:  [integer (:)]  index into column level                  
   )

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      annsum_counter(c) = annsum_counter(c) + dt
   end do

   if (num_soilc .gt. 0) then

   if (annsum_counter(filter_soilc(1)) >= get_days_per_year() * secspday) then
      ! pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
            ! update annual plant ndemand accumulator
            annsum_potential_gpp(p)  = tempsum_potential_gpp(p)
            tempsum_potential_gpp(p) = 0._r8

            ! update annual total N retranslocation accumulator
            annmax_retransn(p)  = tempmax_retransn(p)
            tempmax_retransn(p) = 0._r8

            ! update annual average 2m air temperature accumulator
            annavg_t2m(p)  = tempavg_t2m(p)
            tempavg_t2m(p) = 0._r8

            ! update annual NPP accumulator, convert to annual total
            annsum_npp(p) = tempsum_npp(p) * dt
            tempsum_npp(p) = 0._r8

            if (use_cndv) then
               ! update annual litfall accumulator, convert to annual total
               annsum_litfall(p) = tempsum_litfall(p) * dt
               tempsum_litfall(p) = 0._r8
            end if
      end do

      ! use p2c routine to get selected column-average pft-level fluxes and states
      call p2c(bounds, num_soilc, filter_soilc, &
           annsum_npp(bounds%begp:bounds%endp), &
           cannsum_npp(bounds%begc:bounds%endc))
      call p2c(bounds, num_soilc, filter_soilc, &
           annavg_t2m(bounds%begp:bounds%endp), &
           cannavg_t2m(bounds%begc:bounds%endc))
   end if

   end if

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      if (annsum_counter(c) >= get_days_per_year() * secspday) annsum_counter(c) = 0._r8
   end do

    end associate 
 end subroutine CNAnnualUpdate
!-----------------------------------------------------------------------

end module CNAnnualUpdateMod
