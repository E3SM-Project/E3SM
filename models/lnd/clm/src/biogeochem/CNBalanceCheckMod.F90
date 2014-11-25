module CNBalanceCheckMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon mass balance checking.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog, use_nitrif_denitrif, use_ed
  use clm_time_manager    , only : get_step_size
  use clm_varpar          , only : crop_prog
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnType          , only : col                
  use GridcellType        , only : grc
  !
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
  subroutine BeginCBalance(bounds, num_soilc, filter_soilc, &
       carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning carbon balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                        & 
         totcolc   =>  carbonstate_vars%totcolc_col , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         col_begcb =>  carbonstate_vars%begcb_col     & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
         )

      ! calculate beginning column-level carbon balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begcb(c) = totcolc(c)
      end do

    end associate

  end subroutine BeginCBalance
 
  !-----------------------------------------------------------------------
  subroutine BeginNBalance(bounds, num_soilc, filter_soilc, &
       nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning nitrogen balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    integer                  , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                         & 
         totcoln   => nitrogenstate_vars%totcoln_col , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         col_begnb => nitrogenstate_vars%begnb_col     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         )

      ! calculate beginning column-level nitrogen balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begnb(c) = totcoln(c)
      end do

    end associate

  end subroutine BeginNBalance

  !-----------------------------------------------------------------------
  subroutine CBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       carbonstate_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,err_index    ! indices
    integer  :: fc             ! lake filter indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: col_cinputs
    real(r8) :: col_coutputs
    !-----------------------------------------------------------------------

    associate(                                                                   & 
         totcolc                 =>    carbonstate_vars%totcolc_col            , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         
         gpp                     =>    carbonflux_vars%gpp_col                 , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production      
         er                      =>    carbonflux_vars%er_col                  , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs          =>    carbonflux_vars%fire_closs_col      , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_hrv_xsmrpool_to_atm =>    carbonflux_vars%hrv_xsmrpool_to_atm_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality 
         dwt_closs               =>    carbonflux_vars%dwt_closs_col           , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total carbon loss from product pools and conversion
         product_closs           =>    carbonflux_vars%product_closs_col       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total wood product carbon loss
         som_c_leached           =>    carbonflux_vars%som_c_leached_col       , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss from vertical transport 
         
         col_begcb               =>    carbonstate_vars%begcb_col              , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         col_endcb               =>    carbonstate_vars%endcb_col              , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2) 
         col_errcb               =>    carbonstate_vars%errcb_col                & ! Output: [real(r8) (:) ]  carbon balance error for the timestep (gC/m**2)
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
         col_cinputs = gpp(c)

         ! calculate total column-level outputs
         ! er = ar + hr, col_fire_closs includes pft-level fire losses
         col_coutputs = er(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - (col_endcb(c) - col_begcb(c))

         ! check for significant errors
         if (abs(col_errcb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if

      end do ! end of columns loop

      if (.not. use_ed) then
         if (err_found) then
            c = err_index
            write(iulog,*)'column cbalance error = ', col_errcb(c), c
            write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
            write(iulog,*)'begcb       = ',col_begcb(c)
            write(iulog,*)'endcb       = ',col_endcb(c)
            write(iulog,*)'delta store = ',col_endcb(c)-col_begcb(c)
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end if !use_ed

    end associate

  end subroutine CBalanceCheck

  !-----------------------------------------------------------------------
  subroutine NBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform nitrogen mass conservation check
    ! for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                             & 
         totcoln             =>    nitrogenstate_vars%totcoln_col            , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg 
         ndep_to_sminn       =>    nitrogenflux_vars%ndep_to_sminn_col       , & ! Input:  [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         nfix_to_sminn       =>    nitrogenflux_vars%nfix_to_sminn_col       , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         fert_to_sminn       =>    nitrogenflux_vars%fert_to_sminn_col       , & ! Input:  [real(r8) (:)]                                          
         soyfixn_to_sminn    =>    nitrogenflux_vars%soyfixn_to_sminn_col    , & ! Input:  [real(r8) (:)]                                          
         supplement_to_sminn =>    nitrogenflux_vars%supplement_to_sminn_col , & ! Input:  [real(r8) (:)]  supplemental N supply (gN/m2/s)         
         denit               =>    nitrogenflux_vars%denit_col               , & ! Input:  [real(r8) (:)]  total rate of denitrification (gN/m2/s) 
         sminn_leached       =>    nitrogenflux_vars%sminn_leached_col       , & ! Input:  [real(r8) (:)]  soil mineral N pool loss to leaching (gN/m2/s)
         smin_no3_leached    =>    nitrogenflux_vars%smin_no3_leached_col    , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to leaching (gN/m2/s)
         smin_no3_runoff     =>    nitrogenflux_vars%smin_no3_runoff_col     , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to runoff (gN/m2/s)
         f_n2o_nit           =>    nitrogenflux_vars%f_n2o_nit_col           , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m^2/s]
         col_fire_nloss      =>    nitrogenflux_vars%fire_nloss_col      , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
         dwt_nloss           =>    nitrogenflux_vars%dwt_nloss_col           , & ! Input:  [real(r8) (:)]  (gN/m2/s) total nitrogen loss from product pools and conversion
         product_nloss       =>    nitrogenflux_vars%product_nloss_col       , & ! Input:  [real(r8) (:)]  (gN/m2/s) total wood product nitrogen loss
         som_n_leached       =>    nitrogenflux_vars%som_n_leached_col       , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
         
         col_ninputs         =>    nitrogenflux_vars%ninputs_col         , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)         
         col_noutputs        =>    nitrogenflux_vars%noutputs_col        , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)        
         col_begnb           =>    nitrogenstate_vars%begnb_col              , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         col_endnb           =>    nitrogenstate_vars%endnb_col              , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
         col_errnb           =>    nitrogenstate_vars%errnb_col                & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
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
