module CNBalanceCheckMod

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
   col_begcb                      => ccbal%begcb
   totcolc                        => ccs%totcolc

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
   col_begnb                      => cnbal%begnb
   totcoln                        => cns%totcoln

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
	totcolc                        => ccs%totcolc
	gpp                            => pcf_a%gpp
	er                             => ccf%er
	col_fire_closs                 => ccf%col_fire_closs
	col_hrv_xsmrpool_to_atm        => pcf_a%hrv_xsmrpool_to_atm
	dwt_closs                      => ccf%dwt_closs
	product_closs                  => ccf%product_closs
	
    col_cinputs                    => ccf%col_cinputs
    col_coutputs                   => ccf%col_coutputs
    col_begcb                      => ccbal%begcb
    col_endcb                      => ccbal%endcb
    col_errcb                      => ccbal%errcb

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
      write(iulog,*)'begcb       = ',col_begcb(c)
      write(iulog,*)'endcb       = ',col_endcb(c)
      write(iulog,*)'delta store = ',col_endcb(c)-col_begcb(c)
      write(iulog,*)'input mass  = ',col_cinputs(c)*dt
      write(iulog,*)'output mass = ',col_coutputs(c)*dt
      write(iulog,*)'net flux    = ',(col_cinputs(c)-col_coutputs(c))*dt
	  write(iulog,*)'nee         = ',ccf%nee(c) * dt
	  write(iulog,*)'gpp         = ',gpp(c) * dt
	  write(iulog,*)'er          = ',er(c) * dt
	  write(iulog,*)'col_fire_closs         = ',col_fire_closs(c) * dt
     write(iulog,*)'col_hrv_xsmrpool_to_atm = ',col_hrv_xsmrpool_to_atm(c) * dt
	  write(iulog,*)'dwt_closs         = ',dwt_closs(c) * dt
	  write(iulog,*)'product_closs         = ',product_closs(c) * dt
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
   real(r8), pointer :: supplement_to_sminn(:)   ! supplemental N supply (gN/m2/s)
   real(r8), pointer :: denit(:)                 ! total rate of denitrification (gN/m2/s)
   real(r8), pointer :: sminn_leached(:)         ! soil mineral N pool loss to leaching (gN/m2/s)
   real(r8), pointer :: col_fire_nloss(:)        ! total column-level fire N loss (gN/m2/s)
   real(r8), pointer :: dwt_nloss(:)              ! (gN/m2/s) total nitrogen loss from product pools and conversion
   real(r8), pointer :: product_nloss(:)          ! (gN/m2/s) total wood product nitrogen loss
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
   real(r8), pointer :: col_ninputs(:)           ! column-level N inputs (gN/m2/s)
   real(r8), pointer :: col_noutputs(:)          ! column-level N outputs (gN/m2/s)
   real(r8), pointer :: col_begnb(:)             ! nitrogen mass, beginning of time step (gN/m**2)
   real(r8), pointer :: col_endnb(:)             ! nitrogen mass, end of time step (gN/m**2)
   real(r8), pointer :: col_errnb(:)             ! nitrogen balance error for the timestep (gN/m**2)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,err_index    ! indices
   integer :: fc             ! lake filter indices
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers to column-level arrays

    totcoln                        => cns%totcoln
	ndep_to_sminn                  => cnf%ndep_to_sminn
    nfix_to_sminn                  => cnf%nfix_to_sminn
    supplement_to_sminn            => cnf%supplement_to_sminn
    denit                          => cnf%denit
    sminn_leached                  => cnf%sminn_leached
    col_fire_nloss                 => cnf%col_fire_nloss
    dwt_nloss                      => cnf%dwt_nloss
    product_nloss                  => cnf%product_nloss

    col_ninputs                    => cnf%col_ninputs
    col_noutputs                   => cnf%col_noutputs
    col_begnb                      => cnbal%begnb
    col_endnb                      => cnbal%endnb
    col_errnb                      => cnbal%errnb

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

      ! calculate total column-level outputs

      col_noutputs(c) = denit(c) + sminn_leached(c) + col_fire_nloss(c) + dwt_nloss(c) + product_nloss(c)

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

end module CNBalanceCheckMod
