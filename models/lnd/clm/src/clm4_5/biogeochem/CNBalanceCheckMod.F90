module CNBalanceCheckMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon mass balance checking.
  !
  ! !USES:
  use abortutils  , only: endrun
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: iulog, use_nitrif_denitrif
  use decompMod   , only: bounds_type
  use shr_log_mod , only : errMsg => shr_log_errMsg
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginCBalance
  public :: BeginNBalance
  public :: CBalanceCheck
  public :: NBalanceCheck
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginCBalance(bounds, num_soilc, filter_soilc)
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
    type(bounds_type) , intent(in) :: bounds          ! bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

   associate(& 
   col_begcb =>  ccbal%begcb , & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
   totcolc   =>  ccs%totcolc   & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
   )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level carbon balance,
      ! for mass conservation check
 
      col_begcb(c) = totcolc(c)

   end do ! end of columns loop
 

   end associate
 end subroutine BeginCBalance

 !-----------------------------------------------------------------------
 subroutine BeginNBalance(bounds, num_soilc, filter_soilc)
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
   type(bounds_type), intent(in) :: bounds ! bounds
   integer, intent(in) :: num_soilc        ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   integer :: c     ! indices
   integer :: fc   ! lake filter indices
   !-----------------------------------------------------------------------
   associate(& 
   col_begnb => cnbal%begnb , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
   totcoln   => cns%totcoln   & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
   )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
 
      ! calculate beginning column-level nitrogen balance,
      ! for mass conservation check
 
      col_begnb(c) = totcoln(c)

   end do ! end of columns loop
 
   end associate
 end subroutine BeginNBalance

 !-----------------------------------------------------------------------
 subroutine CBalanceCheck(bounds, num_soilc, filter_soilc)
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
   type(bounds_type) , intent(in) :: bounds          ! bounds
   integer           , intent(in) :: num_soilc       ! number of soil columns in filter
   integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   integer :: c,err_index    ! indices
   integer :: fc             ! lake filter indices
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
   !-----------------------------------------------------------------------

   associate(& 
   totcolc                 =>    ccs%totcolc               , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
   gpp                     =>    pcf_a%gpp                 , & ! Input:  [real(r8) (:)]  (gC/m2/s) gross primary production      
   er                      =>    ccf%er                    , & ! Input:  [real(r8) (:)]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   col_fire_closs          =>    ccf%col_fire_closs        , & ! Input:  [real(r8) (:)]  (gC/m2/s) total column-level fire C loss
   col_hrv_xsmrpool_to_atm =>    pcf_a%hrv_xsmrpool_to_atm , & ! Input:  [real(r8) (:)]  excess MR pool harvest mortality (gC/m2/s)
   dwt_closs               =>    ccf%dwt_closs             , & ! Input:  [real(r8) (:)]  (gC/m2/s) total carbon loss from product pools and conversion
   product_closs           =>    ccf%product_closs         , & ! Input:  [real(r8) (:)]  (gC/m2/s) total wood product carbon loss
   col_cinputs             =>    ccf%col_cinputs           , & ! Output: [real(r8) (:)]  (gC/m2/s) total column-level carbon inputs (for balance check)
   col_coutputs            =>    ccf%col_coutputs          , & ! Output: [real(r8) (:)]  (gC/m2/s) total column-level carbon outputs (for balance check)
   col_begcb               =>    ccbal%begcb               , & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
   col_endcb               =>    ccbal%endcb               , & ! Output: [real(r8) (:)]  carbon mass, end of time step (gC/m**2) 
   col_errcb               =>    ccbal%errcb               , & ! Output: [real(r8) (:)]  carbon balance error for the timestep (gC/m**2)
   som_c_leached           =>    ccf%som_c_leached           & ! Input:  [real(r8) (:)]  total SOM C loss from vertical transport (gC/m^2/s)
   )
   
   
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
      write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
      write(iulog,*)'begcb       = ',col_begcb(c)
      write(iulog,*)'endcb       = ',col_endcb(c)
      write(iulog,*)'delta store = ',col_endcb(c)-col_begcb(c)
      call endrun(msg=errMsg(__FILE__, __LINE__))
   end if

 end associate
 end subroutine CBalanceCheck

!-----------------------------------------------------------------------
 subroutine NBalanceCheck(bounds, num_soilc, filter_soilc)
   !
   ! !DESCRIPTION:
   ! On the radiation time step, perform nitrogen mass conservation check
   ! for column and pft
   !
   ! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar      , only: crop_prog
   !
   ! !ARGUMENTS:
   implicit none
   type(bounds_type), intent(in) :: bounds  ! bounds
   integer, intent(in) :: num_soilc         ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   integer :: c,err_index,j    ! indices
   integer :: fc             ! lake filter indices
   logical :: err_found      ! error flag
   real(r8):: dt             ! radiation time step (seconds)
   !-----------------------------------------------------------------------
   
   associate(& 
   totcoln             =>    cns%totcoln             , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
   ndep_to_sminn       =>    cnf%ndep_to_sminn       , & ! Input:  [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
   nfix_to_sminn       =>    cnf%nfix_to_sminn       , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
   fert_to_sminn       =>    cnf%fert_to_sminn       , & ! Input:  [real(r8) (:)]                                          
   soyfixn_to_sminn    =>    cnf%soyfixn_to_sminn    , & ! Input:  [real(r8) (:)]                                          
   supplement_to_sminn =>    cnf%supplement_to_sminn , & ! Input:  [real(r8) (:)]  supplemental N supply (gN/m2/s)         
   denit               =>    cnf%denit               , & ! Input:  [real(r8) (:)]  total rate of denitrification (gN/m2/s) 
   sminn_leached       =>    cnf%sminn_leached       , & ! Input:  [real(r8) (:)]  soil mineral N pool loss to leaching (gN/m2/s)
   smin_no3_leached    =>    cnf%smin_no3_leached    , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to leaching (gN/m2/s)
   smin_no3_runoff     =>    cnf%smin_no3_runoff     , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to runoff (gN/m2/s)
   f_n2o_nit           =>    cnf%f_n2o_nit           , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m^2/s]
   col_fire_nloss      =>    cnf%col_fire_nloss      , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
   dwt_nloss           =>    cnf%dwt_nloss           , & ! Input:  [real(r8) (:)]  (gN/m2/s) total nitrogen loss from product pools and conversion
   product_nloss       =>    cnf%product_nloss       , & ! Input:  [real(r8) (:)]  (gN/m2/s) total wood product nitrogen loss
   som_n_leached       =>    cnf%som_n_leached       , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
   col_ninputs         =>    cnf%col_ninputs         , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)         
   col_noutputs        =>    cnf%col_noutputs        , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)        
   col_begnb           =>    cnbal%begnb             , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
   col_endnb           =>    cnbal%endnb             , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
   col_errnb           =>    cnbal%errnb               & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
   )
   
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
      
      if (.not. use_nitrif_denitrif) then
         col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
      else
         col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)
      
         col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
      end if
      
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
      write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
      write(iulog,*)'begnb       = ',col_begnb(c)
      write(iulog,*)'endnb       = ',col_endnb(c)
      write(iulog,*)'delta store = ',col_endnb(c)-col_begnb(c)
      write(iulog,*)'input mass  = ',col_ninputs(c)*dt
      write(iulog,*)'output mass = ',col_noutputs(c)*dt
      write(iulog,*)'net flux    = ',(col_ninputs(c)-col_noutputs(c))*dt
      call endrun(msg=errMsg(__FILE__, __LINE__))
   end if
   
 end associate
end subroutine NBalanceCheck
 
end module CNBalanceCheckMod
