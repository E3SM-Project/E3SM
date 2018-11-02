module EcosystemBalanceCheckMod

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
  use clm_varctl          , only : iulog, use_nitrif_denitrif, use_fates
  use clm_time_manager    , only : get_step_size,get_nstep
  use clm_varpar          , only : crop_prog
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnType          , only : col_pp                
  use GridcellType        , only : grc_pp
  use clm_varpar          , only : nlevdecomp
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : nu_com
  use clm_varctl          , only : ECA_Pconst_RGspin

  use CNDecompCascadeConType , only : decomp_cascade_con
  use clm_varpar          , only: ndecomp_cascade_transitions
  use subgridAveMod       , only : p2c, c2g
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode, pf_hmode
  ! forest fertilization experiment
  use clm_time_manager    , only : get_curr_date
  use CNStateType         , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use clm_varctl          , only : forest_fert_exp
  use VegetationType      , only : veg_pp
  use pftvarcon           , only: noveg
  use clm_varctl          , only : NFIX_PTASE_plant

  !
  implicit none
  save
  private
  real(r8), parameter :: balance_check_tolerance = 1e-8_r8
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginColCBalance
  public :: BeginColNBalance
  public :: BeginColPBalance
  public :: ColCBalanceCheck
  public :: ColNBalanceCheck
  public :: ColPBalanceCheck
  public :: BeginGridCBalanceBeforeDynSubgridDriver
  public :: BeginGridNBalanceBeforeDynSubgridDriver
  public :: BeginGridPBalanceBeforeDynSubgridDriver
  public :: EndGridCBalanceAfterDynSubgridDriver
  public :: EndGridNBalanceAfterDynSubgridDriver
  public :: EndGridPBalanceAfterDynSubgridDriver
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BeginColCBalance(bounds, num_soilc, filter_soilc, &
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

  end subroutine BeginColCBalance
 
  !-----------------------------------------------------------------------
  subroutine BeginColNBalance(bounds, num_soilc, filter_soilc, &
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

  end subroutine BeginColNBalance

  !-----------------------------------------------------------------------
  subroutine BeginColPBalance(bounds, num_soilc, filter_soilc, &
       phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning phosphorus balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds          
    integer                  , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                           &
         totcolp   => phosphorusstate_vars%totcolp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         !X.YANG - checking P balance problem, starting from VEGP 
         totpftp   => phosphorusstate_vars%totpftp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp   => phosphorusstate_vars%totsomp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp   => phosphorusstate_vars%cwdp_col       , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp   => phosphorusstate_vars%totlitp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         sminp   => phosphorusstate_vars%sminp_col     , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
 
         col_begpb => phosphorusstate_vars%begpb_col     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         )

      ! calculate beginning column-level phosphorus balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begpb(c) = totcolp(c)
      end do

    end associate

  end subroutine BeginColPBalance

  !-----------------------------------------------------------------------
  subroutine ColCBalanceCheck(bounds, &
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

    associate(                                                                           &
         totcolc                   =>    carbonstate_vars%totcolc_col                  , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         gpp                       =>    carbonflux_vars%gpp_col                       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production
         er                        =>    carbonflux_vars%er_col                        , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs            =>    carbonflux_vars%fire_closs_col                , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_prod1c_loss           =>    carbonflux_vars%prod1c_loss_col               , & ! Input:  [real(r8) (:) ]  (gC/m2/s) crop leafc harvested
         col_hrv_xsmrpool_to_atm   =>    carbonflux_vars%hrv_xsmrpool_to_atm_col       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality
         som_c_leached             =>    carbonflux_vars%som_c_leached_col             , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss from vertical transport
         col_decompc_delta         =>    carbonflux_vars%externalc_to_decomp_delta_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) summarized net change of whole column C i/o to decomposing pool bwtn time-step
         hrv_deadstemc_to_prod10c  =>    carbonflux_vars%hrv_deadstemc_to_prod10c_col  , & ! Input:  [real(r8) (:) ]  (gC/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemc_to_prod100c =>    carbonflux_vars%hrv_deadstemc_to_prod100c_col , & ! Input:  [real(r8) (:) ]  (gC/m2/s) dead stem C harvest mortality to 100-year product pool
         col_begcb                 =>    carbonstate_vars%begcb_col                    , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         col_endcb                 =>    carbonstate_vars%endcb_col                    , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         col_errcb                 =>    carbonstate_vars%errcb_col                      & ! Output: [real(r8) (:) ]  carbon balance error for the timestep (gC/m**2)
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
         col_coutputs = er(c) + col_fire_closs(c) + col_hrv_xsmrpool_to_atm(c)

         ! Fluxes to product pools are included in column-level outputs: the product
         ! pools are not included in totcolc, so are outside the system with respect to
         ! these balance checks
         col_coutputs = col_coutputs + &
              hrv_deadstemc_to_prod10c(c) + hrv_deadstemc_to_prod100c(c)

         col_coutputs = col_coutputs + col_prod1c_loss(c)

         ! subtract leaching flux
         col_coutputs = col_coutputs - som_c_leached(c)

         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs - col_coutputs)*dt - (col_endcb(c) - col_begcb(c))

         ! adjusting the time-lag of org. C increments to decomposing pools when coupled with PFLOTRAN bgc
         ! (because PF bgc uses the extern C as sink (in, + ) at previous time-step,
         ! but note that it includes possible negative adding)
         if (use_pflotran .and. pf_cmode) then
            col_errcb(c) = col_errcb(c) - col_decompc_delta(c)*dt
            ! here is '-' adjustment. It says that the adding to PF decomp c pools was less.
         end if

         ! check for significant errors
         if (abs(col_errcb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if
      end do ! end of columns loop
      
      ! Consider adapting this check to be fates compliant (rgk 04-2017)
      if (.not. use_fates) then
         if (err_found) then
            c = err_index
            write(iulog,*)'column cbalance error = ', col_errcb(c), c
            write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
            write(iulog,*)'input                 = ',col_cinputs*dt
            write(iulog,*)'output                = ',col_coutputs*dt
            write(iulog,*)'er                    = ',er(c)*dt,carbonflux_vars%hr_col(c)*dt
            write(iulog,*)'fire                  = ',col_fire_closs(c)*dt
            write(iulog,*)'hrv_to_atm            = ',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'hrv_to_prod10         = ',hrv_deadstemc_to_prod10c(c)*dt
            write(iulog,*)'hrv_to_prod100        = ',hrv_deadstemc_to_prod100c(c)*dt
            write(iulog,*)'leach                 = ',som_c_leached(c)*dt
            write(iulog,*)'begcb                 = ',col_begcb(c)
            write(iulog,*)'endcb                 = ',col_endcb(c),carbonstate_vars%totsomc_col(c)
            write(iulog,*)'delta store           = ',col_endcb(c)-col_begcb(c)

            if (use_pflotran .and. pf_cmode) then
               write(iulog,*)'pf_delta_decompc      = ',col_decompc_delta(c)*dt
            end if

            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end if !use_fates

    end associate

  end subroutine ColCBalanceCheck

  !-----------------------------------------------------------------------
  subroutine ColNBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform nitrogen mass conservation check
    ! for column and pft
    !
    use tracer_varcon,  only : is_active_betr_bgc
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j,p  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day) 
    !-----------------------------------------------------------------------

    associate(                                                                             &
         totcoln                   =>    nitrogenstate_vars%totcoln_col                  , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         ndep_to_sminn             =>    nitrogenflux_vars%ndep_to_sminn_col             , & ! Input:  [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         nfix_to_sminn             =>    nitrogenflux_vars%nfix_to_sminn_col             , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         nfix_to_ecosysn           =>    nitrogenflux_vars%nfix_to_ecosysn_col           , &
         fert_to_sminn             =>    nitrogenflux_vars%fert_to_sminn_col             , & ! Input:  [real(r8) (:)]
         soyfixn_to_sminn          =>    nitrogenflux_vars%soyfixn_to_sminn_col          , & ! Input:  [real(r8) (:)]
         supplement_to_sminn       =>    nitrogenflux_vars%supplement_to_sminn_col       , & ! Input:  [real(r8) (:)]  supplemental N supply (gN/m2/s)
         denit                     =>    nitrogenflux_vars%denit_col                     , & ! Input:  [real(r8) (:)]  total rate of denitrification (gN/m2/s)
         sminn_leached             =>    nitrogenflux_vars%sminn_leached_col             , & ! Input:  [real(r8) (:)]  soil mineral N pool loss to leaching (gN/m2/s)
         smin_no3_leached          =>    nitrogenflux_vars%smin_no3_leached_col          , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to leaching (gN/m2/s)
         smin_no3_runoff           =>    nitrogenflux_vars%smin_no3_runoff_col           , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to runoff (gN/m2/s)
         f_n2o_nit                 =>    nitrogenflux_vars%f_n2o_nit_col                 , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m^2/s]
         col_prod1n_loss           =>    nitrogenflux_vars%prod1n_loss_col               , & ! Input:  [real(r8) (:) ]  (gN/m2/s) crop leafc harvested
         col_fire_nloss            =>    nitrogenflux_vars%fire_nloss_col                , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
         hrv_deadstemn_to_prod10n  =>    nitrogenflux_vars%hrv_deadstemn_to_prod10n_col  , & ! Input:  [real(r8) (:)]  (gN/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemn_to_prod100n =>    nitrogenflux_vars%hrv_deadstemn_to_prod100n_col , & ! Input:  [real(r8) (:)]  (gN/m2/s) dead stem C harvest mortality to 100-year product pool
         som_n_leached             =>    nitrogenflux_vars%som_n_leached_col             , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
         supplement_to_plantn      =>    nitrogenflux_vars%supplement_to_plantn          , &
         ! pflotran:
         col_decompn_delta         =>    nitrogenflux_vars%externaln_to_decomp_delta_col , & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column N i/o to decomposing pool bwtn time-step
         col_ninputs               =>    nitrogenflux_vars%ninputs_col                   , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)
         col_noutputs              =>    nitrogenflux_vars%noutputs_col                  , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)
         col_begnb                 =>    nitrogenstate_vars%begnb_col                    , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         col_endnb                 =>    nitrogenstate_vars%endnb_col                    , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
         col_errnb                 =>    nitrogenstate_vars%errnb_col                      & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      call get_curr_date(kyr, kmo, kda, mcsec)

      err_found = .false.
      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         ! calculate total column-level inputs
         if (NFIX_PTASE_plant) then
            col_ninputs(c) = ndep_to_sminn(c) + nfix_to_ecosysn(c) + supplement_to_sminn(c)
         else
            col_ninputs(c) = ndep_to_sminn(c) + nfix_to_sminn(c) + supplement_to_sminn(c)
         end if

         if (crop_prog) col_ninputs(c) = col_ninputs(c) + &
              fert_to_sminn(c) + soyfixn_to_sminn(c)

         do p = col_pp%pfti(c), col_pp%pftf(c)
            if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                col_ninputs(c) = col_ninputs(c) + supplement_to_plantn(p) * veg_pp%wtcol(p)
            end if
         end do

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 1 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_ninputs(c) = col_ninputs(c) + fert_dose(c,kmo)/dt
             end if
         end if

         ! calculate total column-level outputs
         col_noutputs(c) = denit(c) + col_fire_nloss(c)

         ! Fluxes to product pools are included in column-level outputs: the product
         ! pools are not included in totcolc, so are outside the system with respect to
         ! these balance checks
         col_noutputs(c) = col_noutputs(c) + &
              hrv_deadstemn_to_prod10n(c) + hrv_deadstemn_to_prod100n(c)

         if (is_active_betr_bgc)then
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         else
           if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
           else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            if(use_pflotran .and. pf_cmode) then
               ! inclusion of aq. NH4 transport by PFLOTRAN-bgc
               col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
            else

               col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)

            endif

           end if
         endif

         col_noutputs(c) = col_noutputs(c) + col_prod1n_loss(c)
         
         col_noutputs(c) = col_noutputs(c) - som_n_leached(c)

         ! calculate the total column-level nitrogen balance error for this time step
         col_errnb(c) = (col_ninputs(c) - col_noutputs(c))*dt - &
              (col_endnb(c) - col_begnb(c))

         ! adjusting the time-lag of org. N increments to decomposing pools when coupled with PFLOTRAN bgc
         ! (because PF bgc uses the extern N sink (in, +) at previous time-step,
         ! but note that it includes possible negative adding)
         if (use_pflotran .and. pf_cmode) then
            col_errnb(c) = col_errnb(c) - col_decompn_delta(c)*dt
            ! here is '-' adjustment. It says that the adding to PF decomp n pools was less.
         end if

         if (abs(col_errnb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if
      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column nbalance error = ',col_errnb(c), c, get_nstep()
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(iulog,*)'begnb                 = ',col_begnb(c)
         write(iulog,*)'endnb                 = ',col_endnb(c)
         write(iulog,*)'delta store           = ',col_endnb(c)-col_begnb(c)
         write(iulog,*)'input mass            = ',col_ninputs(c)*dt
         write(iulog,*)'output mass           = ',col_noutputs(c)*dt
         write(iulog,*)'net flux              = ',(col_ninputs(c)-col_noutputs(c))*dt
         write(iulog,*)'denit                 = ',denit(c)*dt
         write(iulog,*)'n2onit                = ',f_n2o_nit(c)*dt
         write(iulog,*)'no3 leach             = ',smin_no3_leached(c)*dt 
         write(iulog,*)'no3 runof             = ',smin_no3_runoff(c)*dt
         write(iulog,*)'ndep                  = ',ndep_to_sminn(c)*dt
         write(iulog,*)'nfix                  = ',nfix_to_sminn(c)*dt
         write(iulog,*)'nsup                  = ',supplement_to_sminn(c)*dt
         if(crop_prog) then
            write(iulog,*)'fertm                 = ',fert_to_sminn(c)*dt
            write(iulog,*)'soyfx                 = ',soyfixn_to_sminn(c)*dt
         endif
         write(iulog,*)'fire                  = ',col_fire_nloss(c)*dt

         if (use_pflotran .and. pf_cmode) then
            write(iulog,*)'pf_delta_decompn      = ',col_decompn_delta(c)*dt
         end if
         call endrun(msg=errMsg(__FILE__, __LINE__))


      end if

    end associate

  end subroutine ColNBalanceCheck


  !-----------------------------------------------------------------------
  subroutine ColPBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform phosphorus mass conservation check
    ! for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds          
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,err_index,j,k,p  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)

    real(r8) :: leafp_to_litter_col(bounds%begc:bounds%endc) 
    real(r8) :: frootp_to_litter_col(bounds%begc:bounds%endc) 
    real(r8):: flux_mineralization_col(bounds%begc:bounds%endc)   !  local temperary variable

    integer:: kyr                     ! current year 
    integer:: kmo                     ! month of year  (1, ..., 12)
    integer:: kda                     ! day of month   (1, ..., 31) 
    integer:: mcsec                   ! seconds of day (0, ..., seconds/day) 
    !-----------------------------------------------------------------------

    associate(                                                                            &
         totcolp                   => phosphorusstate_vars%totcolp_col                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         supplement_to_sminp       => phosphorusflux_vars%supplement_to_sminp_col       , & ! Input:  [real(r8) (:)]  supplemental P supply (gP/m2/s)
         sminp_leached             => phosphorusflux_vars%sminp_leached_col             , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         col_fire_ploss            => phosphorusflux_vars%fire_ploss_col                , & ! Input:  [real(r8) (:)]  total column-level fire P loss (gP/m2/s)
         hrv_deadstemp_to_prod10p  => phosphorusflux_vars%hrv_deadstemp_to_prod10p_col  , & ! Input:  [real(r8) (:)]  (gP/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemp_to_prod100p => phosphorusflux_vars%hrv_deadstemp_to_prod100p_col , & ! Input:  [real(r8) (:)]  (gP/m2/s) dead stem C harvest mortality to 100-year product pool
         primp_to_labilep          => phosphorusflux_vars%primp_to_labilep_col          , &
         secondp_to_occlp          => phosphorusflux_vars%secondp_to_occlp_col          , &
         fert_p_to_sminp           => phosphorusflux_vars%fert_p_to_sminp_col           , &
         supplement_to_plantp      => phosphorusflux_vars%supplement_to_plantp          , &
         col_prod1p_loss           => phosphorusflux_vars%prod1p_loss_col               , & ! Input:  [real(r8) (:) ]  (gP/m2/s) crop leafc harvested 
         col_pinputs               => phosphorusflux_vars%pinputs_col                   , & ! Output: [real(r8) (:)]  column-level P inputs (gP/m2/s)
         col_poutputs              => phosphorusflux_vars%poutputs_col                  , & ! Output: [real(r8) (:)]  column-level P outputs (gP/m2/s)
         col_begpb                 => phosphorusstate_vars%begpb_col                    , & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         col_endpb                 => phosphorusstate_vars%endpb_col                    , & ! Output: [real(r8) (:)]  phosphorus mass, end of time step (gP/m**2)
         col_errpb                 => phosphorusstate_vars%errpb_col                    , & ! Output: [real(r8) (:)]  phosphorus balance error for the timestep (gP/m**2)

         totpftp                   => phosphorusstate_vars%totpftp_col                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp                   => phosphorusstate_vars%totsomp_col                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp                      => phosphorusstate_vars%cwdp_col                     , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp                   => phosphorusstate_vars%totlitp_col                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         sminp                     => phosphorusstate_vars%sminp_col                    , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         leafp_to_litter           => phosphorusflux_vars%leafp_to_litter_patch         , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         frootp_to_litter          => phosphorusflux_vars%frootp_to_litter_patch        , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         sminp_to_plant            => phosphorusflux_vars%sminp_to_plant_col            , &
         cascade_receiver_pool     => decomp_cascade_con%cascade_receiver_pool          , &
         pf                        =>  phosphorusflux_vars                              , &
         ps                        =>  phosphorusstate_vars                               &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      call get_curr_date(kyr, kmo, kda, mcsec)

      err_found = .false.

      call p2c(bounds,num_soilc,filter_soilc, &
           leafp_to_litter(bounds%begp:bounds%endp), &
           leafp_to_litter_col(bounds%begc:bounds%endc))
      call p2c(bounds,num_soilc,filter_soilc, &
           frootp_to_litter(bounds%begp:bounds%endp), &
           frootp_to_litter_col(bounds%begc:bounds%endc))

      !! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         flux_mineralization_col(c) = 0._r8
      enddo

      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  flux_mineralization_col(c) = flux_mineralization_col(c) - &
                                               pf%decomp_cascade_sminp_flux_col(c,k)
               end do
         else
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                               pf%decomp_cascade_sminp_flux_col(c,k)

               end do
         endif
      end do

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                       pf%biochem_pmin_col(c)
      end do


      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level phosphorus storage, for mass conservation check
         col_endpb(c) = totcolp(c)


         ! calculate total column-level inputs
         col_pinputs(c) = primp_to_labilep(c) + supplement_to_sminp(c)

         do p = col_pp%pfti(c), col_pp%pftf(c)
            if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                col_pinputs(c) = col_pinputs(c) + supplement_to_plantp(p) * veg_pp%wtcol(p)
            end if
         end do

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 2 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_pinputs(c) = col_pinputs(c) + fert_dose(c,kmo)/dt
             end if
         end if

         col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + col_fire_ploss(c)

         ! Fluxes to product pools are included in column-level outputs: the product
         ! pools are not included in totcolc, so are outside the system with respect to
         ! these balance checks
         col_poutputs(c) = col_poutputs(c) + &
              hrv_deadstemp_to_prod10p(c) + hrv_deadstemp_to_prod100p(c)

         if ((nu_com .ne. 'RD') .and. ECA_Pconst_RGspin) then
            do j = 1, nlevdecomp               
               col_poutputs(c) = col_poutputs(c) + &
                  (ps%solutionp_vr_col_cur(c,j) -  ps%solutionp_vr_col_prev(c,j)  + &
                  ps%labilep_vr_col_cur(c,j) -  ps%labilep_vr_col_prev(c,j) + &
                  ps%secondp_vr_col_cur(c,j) - ps%secondp_vr_col_prev(c,j) ) * dzsoi_decomp(j)/dt
            end do 
         end if

         col_poutputs(c) = col_poutputs(c) + col_prod1p_loss(c)
         
         ! calculate the total column-level phosphorus balance error for this time step
         col_errpb(c) = (col_pinputs(c) - col_poutputs(c))*dt - &
              (col_endpb(c) - col_begpb(c))

         if (abs(col_errpb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if
      end do ! end of columns loop


      if (err_found) then
         c = err_index
         write(iulog,*)'column pbalance error = ', col_errpb(c), c
         write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(iulog,*)'begpb       = ',col_begpb(c)
         write(iulog,*)'endpb       = ',col_endpb(c)
         write(iulog,*)'delta store = ',col_endpb(c)-col_begpb(c)
         write(iulog,*)'input mass  = ',col_pinputs(c)*dt
         write(iulog,*)'output mass = ',col_poutputs(c)*dt
         write(iulog,*)'net flux    = ',(col_pinputs(c)-col_poutputs(c))*dt
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine ColPBalanceCheck

  !-----------------------------------------------------------------------
  subroutine BeginGridCBalanceBeforeDynSubgridDriver(bounds, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning carbon balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    !
    !-----------------------------------------------------------------------

    associate(                                                              &
         totcolc               =>  carbonstate_vars%totcolc_col           , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         begcb_grc             =>  carbonstate_vars%begcb_grc               & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
         )

      call c2g( bounds = bounds, &
           carr = totcolc(bounds%begc:bounds%endc), &
           garr = begcb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')

    end associate

  end subroutine BeginGridCBalanceBeforeDynSubgridDriver
 
  !-----------------------------------------------------------------------
  subroutine BeginGridNBalanceBeforeDynSubgridDriver(bounds, nitrogenstate_vars)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning nitrogen balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    !-----------------------------------------------------------------------

    associate(                                         &
         totcoln   => nitrogenstate_vars%totcoln_col , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         begnb_grc => nitrogenstate_vars%begnb_grc     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         )

      call c2g( bounds = bounds, &
           carr = totcoln(bounds%begc:bounds%endc), &
           garr = begnb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')

    end associate

  end subroutine BeginGridNBalanceBeforeDynSubgridDriver

  !-----------------------------------------------------------------------
  subroutine BeginGridPBalanceBeforeDynSubgridDriver(bounds, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning phosphorus balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    !-----------------------------------------------------------------------

    associate(                                           &
         totcolp   => phosphorusstate_vars%totcolp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         begpb_grc => phosphorusstate_vars%begpb_grc     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         )

      call c2g( bounds = bounds, &
           carr = totcolp(bounds%begc:bounds%endc), &
           garr = begpb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')


    end associate

  end subroutine BeginGridPBalanceBeforeDynSubgridDriver

  !-----------------------------------------------------------------------
  subroutine EndGridCBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, &
       carbonstate_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: grc_cinputs
    real(r8) :: grc_coutputs
    !-----------------------------------------------------------------------

    associate(                                                                       &
         totcolc                   =>    carbonstate_vars%totcolc_col              , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         dwt_prod10c_gain_grc      =>    carbonflux_vars%dwt_prod10c_gain_grc      , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_prod100c_gain_grc     =>    carbonflux_vars%dwt_prod100c_gain_grc     , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_conv_cflux_grc        =>    carbonflux_vars%dwt_conv_cflux_grc        , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_leaf_grc     =>    carbonflux_vars%dwt_seedc_to_leaf_grc     , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_deadstem_grc =>    carbonflux_vars%dwt_seedc_to_deadstem_grc , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         begcb_grc                 =>    carbonstate_vars%begcb_grc                , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         endcb_grc                 =>    carbonstate_vars%endcb_grc                , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         errcb_grc                 =>    carbonstate_vars%errcb_grc                  & ! Output: [real(r8) (:) ]  carbon balance error for the time step (gC/m**2)
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcolc(bounds%begc:bounds%endc), &
           garr = endcb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')

      do g = bounds%begg, bounds%endg
         endcb_grc(g) = endcb_grc(g)

         grc_cinputs = &
              dwt_seedc_to_leaf_grc(g)     + &
              dwt_seedc_to_deadstem_grc(g)

         grc_coutputs = &
              dwt_conv_cflux_grc(g)        + &
              dwt_prod10c_gain_grc(g)      + &
              dwt_prod100c_gain_grc(g)
              

         errcb_grc(g) = (grc_cinputs - grc_coutputs)*dt - (endcb_grc(g) - begcb_grc(g))

         ! check for significant errors
         if (abs(errcb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
         g = err_index
         write(iulog,*)'Grid cbalance error   = ',errcb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_cinputs*dt
         write(iulog,*)'output                = ',grc_coutputs*dt
         write(iulog,*)'error                 = ',errcb_grc(g)*dt
         write(iulog,*)'begcb                 = ',begcb_grc(g)
         write(iulog,*)'endcb                 = ',endcb_grc(g)
         write(iulog,*)'delta store           = ',endcb_grc(g)-begcb_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine EndGridCBalanceAfterDynSubgridDriver

  !-----------------------------------------------------------------------
  subroutine EndGridNBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, &
       nitrogenstate_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform nitrogen mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(in)    :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: grc_ninputs
    real(r8) :: grc_noutputs
    !-----------------------------------------------------------------------

    associate(                                                                       &
         totcoln                   =>    nitrogenstate_vars%totcoln_col              , & ! Input:  [real(r8) (:) ]  (gN/m22)   total column nitrogen, incl veg and cpool
         dwt_prod10n_gain_grc      =>    nitrogenflux_vars%dwt_prod10n_gain_grc      , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_prod100n_gain_grc     =>    nitrogenflux_vars%dwt_prod100n_gain_grc     , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_conv_nflux_grc        =>    nitrogenflux_vars%dwt_conv_nflux_grc        , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_leaf_grc     =>    nitrogenflux_vars%dwt_seedn_to_leaf_grc     , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_deadstem_grc =>    nitrogenflux_vars%dwt_seedn_to_deadstem_grc , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         begnb_grc                 =>    nitrogenstate_vars%begnb_grc                , & ! Output: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         endnb_grc                 =>    nitrogenstate_vars%endnb_grc                , & ! Output: [real(r8) (:) ]  nitrogen mass, end of time step (gN/m2**2)
         errnb_grc                 =>    nitrogenstate_vars%errnb_grc                  & ! Output: [real(r8) (:) ]  nitrogen balance error for the time step (gN/m2**2)
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcoln(bounds%begc:bounds%endc), &
           garr = endnb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')

      do g = bounds%begg, bounds%endg
         endnb_grc(g) = endnb_grc(g)

         grc_ninputs = &
              dwt_seedn_to_leaf_grc(g)     + &
              dwt_seedn_to_deadstem_grc(g)

         grc_noutputs = &
              dwt_conv_nflux_grc(g)        + &
              dwt_prod10n_gain_grc(g)      + &
              dwt_prod100n_gain_grc(g)
              

         errnb_grc(g) = (grc_ninputs - grc_noutputs)*dt - (endnb_grc(g) - begnb_grc(g))

         ! check for significant errors
         if (abs(errnb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
         g = err_index
         write(iulog,*)'Grid nbalance error   = ',errnb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_ninputs*dt
         write(iulog,*)'output                = ',grc_noutputs*dt
         write(iulog,*)'error                 = ',errnb_grc(g)*dt
         write(iulog,*)'begcb                 = ',begnb_grc(g)
         write(iulog,*)'endcb                 = ',endnb_grc(g)
         write(iulog,*)'delta store           = ',endnb_grc(g)-begnb_grc(g)
         write(iulog,*)''
         write(iulog,*)'dwt_conv                ',dwt_conv_nflux_grc(g)
         write(iulog,*)'dwt_prod10              ',dwt_prod10n_gain_grc(g)
         write(iulog,*)'dwt_prod100             ',dwt_prod100n_gain_grc(g)
         write(iulog,*)''
         write(iulog,*)'dwt_seedn_leaf          ',dwt_seedn_to_leaf_grc(g)
         write(iulog,*)'dwt_seedn_deadstem      ',dwt_seedn_to_deadstem_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine EndGridNBalanceAfterDynSubgridDriver

  !-----------------------------------------------------------------------

  subroutine EndGridPBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, &
       phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform phosphorus mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(in)    :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: grc_pinputs
    real(r8) :: grc_poutputs
    !-----------------------------------------------------------------------

    associate(                                                                       &
         totcolp                   =>    phosphorusstate_vars%totcolp_col              , & ! Input:  [real(r8) (:) ]  (gP/m2)   total column phosphorus, incl veg and cpool
         dwt_prod10p_gain_grc      =>    phosphorusflux_vars%dwt_prod10p_gain_grc      , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_prod100p_gain_grc     =>    phosphorusflux_vars%dwt_prod100p_gain_grc     , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_conv_pflux_grc        =>    phosphorusflux_vars%dwt_conv_pflux_grc        , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_leaf_grc     =>    phosphorusflux_vars%dwt_seedp_to_leaf_grc     , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_deadstem_grc =>    phosphorusflux_vars%dwt_seedp_to_deadstem_grc , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         begpb_grc                 =>    phosphorusstate_vars%begpb_grc                , & ! Output: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         endpb_grc                 =>    phosphorusstate_vars%endpb_grc                , & ! Output: [real(r8) (:) ]  phosphorus mass, end of time step (gP/m2**2)
         errpb_grc                 =>    phosphorusstate_vars%errpb_grc                  & ! Output: [real(r8) (:) ]  phosphorus balance error for the time step (gP/m2**2)
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcolp(bounds%begc:bounds%endc), &
           garr = endpb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', &
           l2g_scale_type = 'unity')

      do g = bounds%begg, bounds%endg
         endpb_grc(g) = endpb_grc(g)

         grc_pinputs = &
              dwt_seedp_to_leaf_grc(g)     + &
              dwt_seedp_to_deadstem_grc(g)

         grc_poutputs = &
              dwt_conv_pflux_grc(g)        + &
              dwt_prod10p_gain_grc(g)      + &
              dwt_prod100p_gain_grc(g)
              

         errpb_grc(g) = (grc_pinputs - grc_poutputs)*dt - (endpb_grc(g) - begpb_grc(g))

         ! check for significant errors
         if (abs(errpb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
         g = err_index
         write(iulog,*)'Grid pbalance error   = ',errpb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_pinputs*dt
         write(iulog,*)'output                = ',grc_poutputs*dt
         write(iulog,*)'error                 = ',errpb_grc(g)*dt
         write(iulog,*)'begcb                 = ',begpb_grc(g)
         write(iulog,*)'endcb                 = ',endpb_grc(g)
         write(iulog,*)'delta store           = ',endpb_grc(g)-begpb_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine EndGridPBalanceAfterDynSubgridDriver

end module EcosystemBalanceCheckMod
