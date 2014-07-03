
module CNBalanceCheckMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNBalanceCheckMod
!
! !DESCRIPTION:
! Module for carbon mass balance checking.
!
! !USES:
    use abortutils  , only: endrun
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varctl  , only: iulog
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: BeginCBalance
    public :: BeginNBalance
    public :: CBalanceCheck
    public :: NBalanceCheck
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginCBalance
!
! !INTERFACE:
subroutine BeginCBalance(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning carbon balance for mass
! conservation checks.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_begcb(:)   ! carbon mass, beginning of time step (gC/m**2)

!
! !OTHER LOCAL VARIABLES:
   integer :: c     ! indices
   integer :: fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begcb                      => clm3%g%l%c%ccbal%begcb
   totcolc                        => clm3%g%l%c%ccs%totcolc

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level carbon balance,
      ! for mass conservation check
 
      col_begcb(c) = totcolc(c)

   end do ! end of columns loop
 

end subroutine BeginCBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BeginNBalance
!
! !INTERFACE:
subroutine BeginNBalance(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, calculate the beginning nitrogen balance for mass
! conservation checks.
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 2/4/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_begnb(:)   ! nitrogen mass, beginning of time step (gN/m**2)
!

! !OTHER LOCAL VARIABLES:
   integer :: c     ! indices
   integer :: fc   ! lake filter indices
!
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers at the column level
   col_begnb                      => clm3%g%l%c%cnbal%begnb
   totcoln                        => clm3%g%l%c%cns%totcoln

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level nitrogen balance,
      ! for mass conservation check
 
      col_begnb(c) = totcoln(c)

   end do ! end of columns loop
 
end subroutine BeginNBalance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CBalanceCheck
!
! !INTERFACE:
subroutine CBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, perform carbon mass conservation check for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: gpp(:)            ! (gC/m2/s) gross primary production 
   real(r8), pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
   real(r8), pointer :: col_hrv_xsmrpool_to_atm(:)  ! excess MR pool harvest mortality (gC/m2/s)
   real(r8), pointer :: dwt_closs(:)              ! (gC/m2/s) total carbon loss from product pools and conversion
   real(r8), pointer :: product_closs(:)      ! (gC/m2/s) total wood product carbon loss
   real(r8), pointer :: som_c_leached(:)    ! total SOM C loss from vertical transport (gC/m^2/s)
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_cinputs(:)  ! (gC/m2/s) total column-level carbon inputs (for balance check)
   real(r8), pointer :: col_coutputs(:) ! (gC/m2/s) total column-level carbon outputs (for balance check)
   real(r8), pointer :: col_begcb(:)    ! carbon mass, beginning of time step (gC/m**2)
   real(r8), pointer :: col_endcb(:)    ! carbon mass, end of time step (gC/m**2)
   real(r8), pointer :: col_errcb(:)    ! carbon balance error for the timestep (gC/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,err_index  ! indices
   integer :: fc          ! lake filter indices
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to column-level arrays
   totcolc                        => clm3%g%l%c%ccs%totcolc
   gpp                            => clm3%g%l%c%ccf%pcf_a%gpp
   er                             => clm3%g%l%c%ccf%er
   col_fire_closs                 => clm3%g%l%c%ccf%col_fire_closs
   col_hrv_xsmrpool_to_atm        => clm3%g%l%c%ccf%pcf_a%hrv_xsmrpool_to_atm
   dwt_closs                      => clm3%g%l%c%ccf%dwt_closs
   product_closs                  => clm3%g%l%c%ccf%product_closs
   
   col_cinputs                    => clm3%g%l%c%ccf%col_cinputs
   col_coutputs                   => clm3%g%l%c%ccf%col_coutputs
   col_begcb                      => clm3%g%l%c%ccbal%begcb
   col_endcb                      => clm3%g%l%c%ccbal%endcb
   col_errcb                      => clm3%g%l%c%ccbal%errcb
   som_c_leached                  => clm3%g%l%c%ccf%som_c_leached
   
   
   ! set time steps
   dt = real( get_step_size(), r8 )
   
   err_found = .false.
   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      
      ! calculate the total column-level carbon storage, for mass conservation check
      
      col_endcb(c) = totcolc(c)
      
      ! calculate total column-level inputs
      
      col_cinputs(c) = gpp(c)
      
      ! calculate total column-level outputs
      ! er = ar + hr, col_fire_closs includes pft-level fire losses
      
      col_coutputs(c) = er(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)
      
      ! subtract leaching flux
      col_coutputs(c) = col_coutputs(c) - som_c_leached(c)
      
      ! calculate the total column-level carbon balance error for this time step
      col_errcb(c) = (col_cinputs(c) - col_coutputs(c))*dt - &
           (col_endcb(c) - col_begcb(c))
      
      ! check for significant errors
      if (abs(col_errcb(c)) > 1e-8_r8) then
         err_found = .true.
         err_index = c
      end if
      
   end do ! end of columns loop
   
   if (err_found) then
      c = err_index
      write(iulog,*)'column cbalance error = ', col_errcb(c), c
      write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(clm3%g%l%c%gridcell(c)),clm3%g%londeg(clm3%g%l%c%gridcell(c))
      write(iulog,*)'begcb       = ',col_begcb(c)
      write(iulog,*)'endcb       = ',col_endcb(c)
      write(iulog,*)'delta store = ',col_endcb(c)-col_begcb(c)
      
      call endrun
   end if
   

end subroutine CBalanceCheck
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NBalanceCheck
!
! !INTERFACE:
subroutine NBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! On the radiation time step, perform nitrogen mass conservation check
! for column and pft
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use surfrdMod       , only: crop_prog
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
!
! !CALLED FROM:
! subroutine clm_driver1
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   real(r8), pointer :: totcoln(:)               ! (gN/m2) total column nitrogen, incl veg
   real(r8), pointer :: ndep_to_sminn(:)         ! atmospheric N deposition to soil mineral N (gN/m2/s)
   real(r8), pointer :: nfix_to_sminn(:)         ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
   real(r8), pointer :: fert_to_sminn(:)
   real(r8), pointer :: soyfixn_to_sminn(:)
   real(r8), pointer :: supplement_to_sminn(:)   ! supplemental N supply (gN/m2/s)
   real(r8), pointer :: denit(:)                 ! total rate of denitrification (gN/m2/s)
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_leached(:)         ! soil mineral N pool loss to leaching (gN/m2/s)
#else
   real(r8), pointer :: smin_no3_leached(:)      ! soil mineral NO3 pool loss to leaching (gN/m2/s)
   real(r8), pointer :: smin_no3_runoff(:)       ! soil mineral NO3 pool loss to runoff (gN/m2/s)
   real(r8), pointer :: f_n2o_nit(:)             ! flux of N2o from nitrification [gN/m^2/s]
#endif
   real(r8), pointer :: col_fire_nloss(:)        ! total column-level fire N loss (gN/m2/s)
   real(r8), pointer :: dwt_nloss(:)              ! (gN/m2/s) total nitrogen loss from product pools and conversion
   real(r8), pointer :: product_nloss(:)          ! (gN/m2/s) total wood product nitrogen loss
   real(r8), pointer :: som_n_leached(:)    ! total SOM N loss from vertical transport
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_ninputs(:)           ! column-level N inputs (gN/m2/s)
   real(r8), pointer :: col_noutputs(:)          ! column-level N outputs (gN/m2/s)
   real(r8), pointer :: col_begnb(:)             ! nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: col_endnb(:)             ! nitrogen mass, end of time step (gN/m**2)
   real(r8), pointer :: col_errnb(:)             ! nitrogen balance error for the timestep (gN/m**2)

! !OTHER LOCAL VARIABLES:
   integer :: c,err_index,j    ! indices
   integer :: fc             ! lake filter indices
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers to column-level arrays
   
   totcoln                        => clm3%g%l%c%cns%totcoln
   ndep_to_sminn                  => clm3%g%l%c%cnf%ndep_to_sminn
   nfix_to_sminn                  => clm3%g%l%c%cnf%nfix_to_sminn
   fert_to_sminn                  => clm3%g%l%c%cnf%fert_to_sminn
   soyfixn_to_sminn               => clm3%g%l%c%cnf%soyfixn_to_sminn
   supplement_to_sminn            => clm3%g%l%c%cnf%supplement_to_sminn
   denit                          => clm3%g%l%c%cnf%denit
#ifndef NITRIF_DENITRIF
   sminn_leached                  => clm3%g%l%c%cnf%sminn_leached
#else
   smin_no3_leached               => clm3%g%l%c%cnf%smin_no3_leached
   smin_no3_runoff                => clm3%g%l%c%cnf%smin_no3_runoff
   f_n2o_nit                      => clm3%g%l%c%cnf%f_n2o_nit
#endif
   col_fire_nloss                 => clm3%g%l%c%cnf%col_fire_nloss
   dwt_nloss                      => clm3%g%l%c%cnf%dwt_nloss
   product_nloss                  => clm3%g%l%c%cnf%product_nloss
   som_n_leached                  => clm3%g%l%c%cnf%som_n_leached
   
   col_ninputs                    => clm3%g%l%c%cnf%col_ninputs
   col_noutputs                   => clm3%g%l%c%cnf%col_noutputs
   col_begnb                      => clm3%g%l%c%cnbal%begnb
   col_endnb                      => clm3%g%l%c%cnbal%endnb
   col_errnb                      => clm3%g%l%c%cnbal%errnb
   
   ! set time steps
   dt = real( get_step_size(), r8 )
   
   err_found = .false.
   ! column loop
   do fc = 1,num_soilc
      c=filter_soilc(fc)
      
      ! calculate the total column-level nitrogen storage, for mass conservation check
      
      col_endnb(c) = totcoln(c)
      
      ! calculate total column-level inputs
      
      col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)
      if (crop_prog) col_ninputs(c) = col_ninputs(c) + &
                       fert_to_sminn(c) + soyfixn_to_sminn(c)
      
      ! calculate total column-level outputs
      
      col_noutputs(c) = denit(c) + col_fire_nloss(c) + dwt_nloss(c) + product_nloss(c)
      
#ifndef NITRIF_DENITRIF
      col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
#else
      col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)
      
      col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
#endif
      
      col_noutputs(c) = col_noutputs(c) - som_n_leached(c)
      
      ! calculate the total column-level nitrogen balance error for this time step
      
      col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
           (col_endnb(c) - col_begnb(c))
      
      if (abs(col_errnb(c)) > 1e-8_r8) then
         err_found = .true.
         err_index = c
      end if
      
   end do ! end of columns loop
   
   if (err_found) then
      c = err_index
      write(iulog,*)'column nbalance error = ', col_errnb(c), c
      write(iulog,*)'Latdeg,Londeg=',clm3%g%latdeg(clm3%g%l%c%gridcell(c)),clm3%g%londeg(clm3%g%l%c%gridcell(c))
      write(iulog,*)'begnb       = ',col_begnb(c)
      write(iulog,*)'endnb       = ',col_endnb(c)
      write(iulog,*)'delta store = ',col_endnb(c)-col_begnb(c)
      write(iulog,*)'input mass  = ',col_ninputs(c)*dt
      write(iulog,*)'output mass = ',col_noutputs(c)*dt
      write(iulog,*)'net flux    = ',(col_ninputs(c)-col_noutputs(c))*dt
      
      call endrun
   end if
   
 end subroutine NBalanceCheck
 !-----------------------------------------------------------------------
#endif
 
end module CNBalanceCheckMod
