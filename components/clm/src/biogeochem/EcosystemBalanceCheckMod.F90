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
  use clm_varpar          , only : nlevdecomp
  use clm_varcon          , only : dzsoi_decomp
  use clm_varctl          , only : nu_com
  use clm_varctl          , only : ECA_Pconst_RGspin

  use CNDecompCascadeConType , only : decomp_cascade_con
  use clm_varpar          , only: ndecomp_cascade_transitions
  use subgridAveMod       , only : p2c, c2g
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  ! soil erosion
  use clm_varctl          , only : use_erosion, ero_ccycle
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode, pf_hmode
  ! forest fertilization experiment
  use clm_time_manager    , only : get_curr_date
  use CNStateType         , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use clm_varctl          , only : forest_fert_exp
  use pftvarcon           , only: noveg
  use clm_varctl          , only : NFIX_PTASE_plant
  use GridcellType        , only : grc_pp
  use GridcellDataType    , only : gridcell_carbon_state, grc_cf
  use GridcellDataType    , only : grc_ns, grc_nf, grc_ps, grc_pf
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : column_carbon_state, col_cf
  use ColumnDataType      , only : col_ns, col_nf, col_ps, col_pf
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_cf, veg_nf, veg_pf


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
       col_cs)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning carbon balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_carbon_state) , intent(inout) :: col_cs
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                        &
         totcolc   =>  col_cs%totcolc , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         col_begcb =>  col_cs%begcb     & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
         )

      ! calculate beginning column-level carbon balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begcb(c) = totcolc(c)
         col_cs%totblgc_beg(c) = col_cs%totblgc(c)
         col_cs%totabgc_beg(c) = col_cs%totabgc(c)
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
         totcoln   => col_ns%totcoln , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         col_begnb => col_ns%begnb     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         )

      ! calculate beginning column-level nitrogen balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begnb(c) = totcoln(c)
         col_ns%totabgn_beg(c)=col_ns%totabgn(c)
         col_ns%totblgn_beg(c)=col_ns%totblgn(c)
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
         totcolp   => col_ps%totcolp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         !X.YANG - checking P balance problem, starting from VEGP
         totpftp   => col_ps%totpftp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totsomp   => col_ps%totsomp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         cwdp   => col_ps%cwdp       , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp   => col_ps%totlitp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         sminp   => col_ps%sminp     , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg

         col_begpb => col_ps%begpb     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
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
       col_cs, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_carbon_state) , intent(inout) :: col_cs
    type(carbonflux_type)     , intent(in)    :: carbonflux_vars
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
         totcolc                   =>    col_cs%totcolc                  , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         gpp                       =>    col_cf%gpp                       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production
         er                        =>    col_cf%er                        , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs            =>    col_cf%fire_closs                , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_prod1c_loss           =>    col_cf%prod1c_loss               , & ! Input:  [real(r8) (:) ]  (gC/m2/s) crop leafc harvested
         col_hrv_xsmrpool_to_atm   =>    col_cf%hrv_xsmrpool_to_atm       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality
         som_c_leached             =>    col_cf%som_c_leached             , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss from vertical transport
         som_c_yield               =>    col_cf%somc_yield                , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss by erosion
         col_decompc_delta         =>    col_cf%externalc_to_decomp_delta , & ! Input:  [real(r8) (:) ]  (gC/m2/s) summarized net change of whole column C i/o to decomposing pool bwtn time-step
         hrv_deadstemc_to_prod10c  =>    col_cf%hrv_deadstemc_to_prod10c  , & ! Input:  [real(r8) (:) ]  (gC/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemc_to_prod100c =>    col_cf%hrv_deadstemc_to_prod100c , & ! Input:  [real(r8) (:) ]  (gC/m2/s) dead stem C harvest mortality to 100-year product pool
         col_begcb                 =>    col_cs%begcb                    , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         col_endcb                 =>    col_cs%endcb                    , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         col_errcb                 =>    col_cs%errcb                      & ! Output: [real(r8) (:) ]  carbon balance error for the timestep (gC/m**2)
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

         ! add erosion flux
         if (ero_ccycle) then
            col_coutputs = col_coutputs + som_c_yield(c)
         end if

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
            write(iulog,*)'nstep                 = ', get_nstep()
            write(iulog,*)'column cbalance error = ', col_errcb(c), c
            write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
            write(iulog,*)'input                 = ',col_cinputs*dt
            write(iulog,*)'output                = ',col_coutputs*dt
            write(iulog,*)'er                    = ',er(c)*dt,col_cf%hr(c)*dt
            write(iulog,*)'fire                  = ',col_fire_closs(c)*dt
            write(iulog,*)'hrv_to_atm            = ',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'hrv_to_prod10         = ',hrv_deadstemc_to_prod10c(c)*dt
            write(iulog,*)'hrv_to_prod100        = ',hrv_deadstemc_to_prod100c(c)*dt
            write(iulog,*)'leach                 = ',som_c_leached(c)*dt
            write(iulog,*)'begcb                 = ',col_begcb(c)
            write(iulog,*)'endcb                 = ',col_endcb(c),col_cs%totsomc(c)
            write(iulog,*)'delta store           = ',col_endcb(c)-col_begcb(c)
            write(iulog,*)'litter fall           = ',col_cf%plant_to_litter_cflux(c)*dt
            write(iulog,*)'plant to cwd          = ',col_cf%plant_to_cwd_cflux(c)*dt
            write(iulog,*)'tot pft               = ',col_cs%totpftc(c)
            if (ero_ccycle) then
               write(iulog,*)'erosion               = ',som_c_yield(c)*dt
            end if

            if (use_pflotran .and. pf_cmode) then
               write(iulog,*)'pf_delta_decompc      = ',col_decompc_delta(c)*dt
            end if
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
         totcoln                   =>    col_ns%totcoln                  , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         ndep_to_sminn             =>    col_nf%ndep_to_sminn             , & ! Input:  [real(r8) (:)]  atmospheric N deposition to soil mineral N (gN/m2/s)
         nfix_to_sminn             =>    col_nf%nfix_to_sminn             , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
         nfix_to_ecosysn           =>    col_nf%nfix_to_ecosysn           , &
         fert_to_sminn             =>    col_nf%fert_to_sminn             , & ! Input:  [real(r8) (:)]
         soyfixn_to_sminn          =>    col_nf%soyfixn_to_sminn          , & ! Input:  [real(r8) (:)]
         supplement_to_sminn       =>    col_nf%supplement_to_sminn       , & ! Input:  [real(r8) (:)]  supplemental N supply (gN/m2/s)
         denit                     =>    col_nf%denit                     , & ! Input:  [real(r8) (:)]  total rate of denitrification (gN/m2/s)
         sminn_leached             =>    col_nf%sminn_leached             , & ! Input:  [real(r8) (:)]  soil mineral N pool loss to leaching (gN/m2/s)
         smin_no3_leached          =>    col_nf%smin_no3_leached          , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to leaching (gN/m2/s)
         smin_no3_runoff           =>    col_nf%smin_no3_runoff           , & ! Input:  [real(r8) (:)]  soil mineral NO3 pool loss to runoff (gN/m2/s)
         f_n2o_nit                 =>    col_nf%f_n2o_nit                 , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m^2/s]
         col_prod1n_loss           =>    col_nf%prod1n_loss               , & ! Input:  [real(r8) (:) ]  (gN/m2/s) crop leafc harvested
         col_fire_nloss            =>    col_nf%fire_nloss                , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
         hrv_deadstemn_to_prod10n  =>    col_nf%hrv_deadstemn_to_prod10n  , & ! Input:  [real(r8) (:)]  (gN/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemn_to_prod100n =>    col_nf%hrv_deadstemn_to_prod100n , & ! Input:  [real(r8) (:)]  (gN/m2/s) dead stem C harvest mortality to 100-year product pool
         som_n_leached             =>    col_nf%som_n_leached             , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
         som_n_yield               =>    col_nf%somn_yield                , & ! Input:  [real(r8) (:)]  total SOM N loss by erosion
         supplement_to_plantn      =>    col_nf%supplement_to_plantn      , &
         ! pflotran:
         col_decompn_delta         =>    col_nf%externaln_to_decomp_delta , & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column N i/o to decomposing pool bwtn time-step
         col_ninputs               =>    col_nf%ninputs                   , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)
         col_noutputs              =>    col_nf%noutputs                  , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)
         col_begnb                 =>    col_ns%begnb                    , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         col_endnb                 =>    col_ns%endnb                    , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
         col_errnb                 =>    col_ns%errnb                      & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
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

         col_ninputs(c) = col_ninputs(c) + supplement_to_plantn(c)

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

         col_noutputs(c) = col_noutputs(c) + col_prod1n_loss(c)

         col_noutputs(c) = col_noutputs(c) - som_n_leached(c)

         ! subtracted erosion flux
         if (ero_ccycle) then
            col_noutputs(c) = col_noutputs(c) + som_n_yield(c)
         end if

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
         if (NFIX_PTASE_plant) then
           write(iulog,*)'nfix                  = ',nfix_to_ecosysn(c)*dt,nfix_to_sminn(c)*dt
         else
           write(iulog,*)'nfix                  = ',nfix_to_sminn(c)*dt
         endif
         write(iulog,*)'nsup                  = ',supplement_to_sminn(c)*dt
         write(iulog,*)'nsup_plant            = ',supplement_to_plantn(c)*dt
         if(crop_prog) then
            write(iulog,*)'fertm                 = ',fert_to_sminn(c)*dt
            write(iulog,*)'soyfx                 = ',soyfixn_to_sminn(c)*dt
         endif
         write(iulog,*)'fire                  = ',col_fire_nloss(c)*dt

         if (ero_ccycle) then
            write(iulog,*)'erosion               = ',som_n_yield(c)*dt
         end if

         if (use_pflotran .and. pf_cmode) then
            write(iulog,*)'pf_delta_decompn      = ',col_decompn_delta(c)*dt
         end if

         write(iulog,*)'abgn=',col_ns%totabgn_beg(c),col_ns%totabgn(c)
         write(iulog,*)'blgn=',col_ns%totblgn_beg(c),col_ns%totblgn(c)
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
         totcolp                   => col_ps%totcolp                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         supplement_to_sminp       => col_pf%supplement_to_sminp       , & ! Input:  [real(r8) (:)]  supplemental P supply (gP/m2/s)
         sminp_leached             => col_pf%sminp_leached             , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         col_fire_ploss            => col_pf%fire_ploss                , & ! Input:  [real(r8) (:)]  total column-level fire P loss (gP/m2/s)
         hrv_deadstemp_to_prod10p  => col_pf%hrv_deadstemp_to_prod10p  , & ! Input:  [real(r8) (:)]  (gP/m2/s) dead stem C harvest mortality to 10-year product pool
         hrv_deadstemp_to_prod100p => col_pf%hrv_deadstemp_to_prod100p , & ! Input:  [real(r8) (:)]  (gP/m2/s) dead stem C harvest mortality to 100-year product pool
         primp_to_labilep          => col_pf%primp_to_labilep          , &
         secondp_to_occlp          => col_pf%secondp_to_occlp          , &
         fert_p_to_sminp           => col_pf%fert_p_to_sminp           , &
         som_p_yield               => col_pf%somp_yield                , & ! Input:  [real(r8) (:)]  SOM P pool loss by erosion (gP/m^2/s)
         labilep_yield             => col_pf%labilep_yield             , & ! Input:  [real(r8) (:)]  soil labile mineral P loss by erosion (gP/m^s/s)
         secondp_yield             => col_pf%secondp_yield             , & ! Input:  [real(r8) (:)]  soil secondary mineral P loss by erosion (gP/m^s/s)
         occlp_yield               => col_pf%occlp_yield               , & ! Input:  [real(r8) (:)]  soil occluded mineral P loss by erosion (gP/m^s/s)
         primp_yield               => col_pf%primp_yield               , & ! Input:  [real(r8) (:)]  soil primary mineral P loss by erosion (gP/m^s/s)
         supplement_to_plantp      => veg_pf%supplement_to_plantp          , &
         col_prod1p_loss           => col_pf%prod1p_loss               , & ! Input:  [real(r8) (:) ]  (gP/m2/s) crop leafc harvested
         col_pinputs               => col_pf%pinputs                   , & ! Output: [real(r8) (:)]  column-level P inputs (gP/m2/s)
         col_poutputs              => col_pf%poutputs                  , & ! Output: [real(r8) (:)]  column-level P outputs (gP/m2/s)
         col_begpb                 => col_ps%begpb                    , & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         col_endpb                 => col_ps%endpb                    , & ! Output: [real(r8) (:)]  phosphorus mass, end of time step (gP/m**2)
         col_errpb                 => col_ps%errpb                    , & ! Output: [real(r8) (:)]  phosphorus balance error for the timestep (gP/m**2)

         totpftp                   => col_ps%totpftp                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totsomp                   => col_ps%totsomp                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         cwdp                      => col_ps%cwdp                     , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp                   => col_ps%totlitp                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         sminp                     => col_ps%sminp                    , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         leafp_to_litter           => veg_pf%leafp_to_litter         , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         frootp_to_litter          => veg_pf%frootp_to_litter        , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         sminp_to_plant            => col_pf%sminp_to_plant            , &
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
                                               col_pf%decomp_cascade_sminp_flux(c,k)
               end do
         else
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                               col_pf%decomp_cascade_sminp_flux(c,k)

               end do
         endif
      end do

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                       col_pf%biochem_pmin(c)
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
                  (col_ps%solutionp_vr_cur(c,j) -  col_ps%solutionp_vr_prev(c,j)  + &
                  col_ps%labilep_vr_cur(c,j) -  col_ps%labilep_vr_prev(c,j) + &
                  col_ps%secondp_vr_cur(c,j) - col_ps%secondp_vr_prev(c,j) ) * dzsoi_decomp(j)/dt
            end do
         end if

         col_poutputs(c) = col_poutputs(c) + col_prod1p_loss(c)

         ! soil erosion
         if (ero_ccycle) then
            col_poutputs(c) = col_poutputs(c) + som_p_yield(c) + labilep_yield(c) + &
               secondp_yield(c) !+ occlp_yield(c) + primp_yield(c)
         end if

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
         write(iulog,*)'primp in    = ',primp_to_labilep(c)
         write(iulog,*)'supp        = ',supplement_to_sminp(c)
         if (ero_ccycle) then
            write(iulog,*)'SOP erosion = ',som_p_yield(c)*dt
            write(iulog,*)'SIP erosion = ',(labilep_yield(c)+secondp_yield(c)+occlp_yield(c)+primp_yield(c))*dt
         end if

         write(*,*)'column pbalance error = ', col_errpb(c), c
         write(*,*)'Latdeg,Londeg=',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(*,*)'begpb       = ',col_begpb(c)
         write(*,*)'endpb       = ',col_endpb(c)
         write(*,*)'delta store = ',col_endpb(c)-col_begpb(c)
         write(*,*)'input mass  = ',col_pinputs(c)*dt
         write(*,*)'output mass = ',col_poutputs(c)*dt
         write(*,*)'net flux    = ',(col_pinputs(c)-col_poutputs(c))*dt
         write(*,*)'primp in    = ',primp_to_labilep(c)
         write(*,*)'supp        = ',supplement_to_sminp(c)
         if (ero_ccycle) then
            write(iulog,*)'SOP erosion = ',som_p_yield(c)*dt
            write(iulog,*)'SIP erosion = ',(labilep_yield(c)+secondp_yield(c)+occlp_yield(c)+primp_yield(c))*dt
         end if
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine ColPBalanceCheck

  !-----------------------------------------------------------------------
  subroutine BeginGridCBalanceBeforeDynSubgridDriver(bounds, col_cs, grc_cs)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning carbon balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    type(column_carbon_state)  , intent(inout) :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    !
    !-----------------------------------------------------------------------

    associate(                                                              &
         totcolc           =>  col_cs%totcolc           , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool
         begcb_grc         =>  grc_cs%begcb               & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
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
         totcoln   => col_ns%totcoln , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         begnb_grc => grc_ns%begnb     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
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
         totcolp   => col_ps%totcolp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         begpb_grc => grc_ps%begpb     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
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
       col_cs, grc_cs, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    integer                    , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                    , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_carbon_state)  , intent(inout) :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    type(carbonflux_type)      , intent(in)    :: carbonflux_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    real(r8) :: grc_cinputs
    real(r8) :: grc_coutputs
    !-----------------------------------------------------------------------

    associate(                                                                       &
         totcolc                   =>    col_cs%totcolc              , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         dwt_prod10c_gain_grc      =>    grc_cf%dwt_prod10c_gain      , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_prod100c_gain_grc     =>    grc_cf%dwt_prod100c_gain     , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_conv_cflux_grc        =>    grc_cf%dwt_conv_cflux        , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_leaf_grc     =>    grc_cf%dwt_seedc_to_leaf     , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_deadstem_grc =>    grc_cf%dwt_seedc_to_deadstem , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         begcb_grc                 =>    grc_cs%begcb                , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         endcb_grc                 =>    grc_cs%endcb                , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         errcb_grc                 =>    grc_cs%errcb                  & ! Output: [real(r8) (:) ]  carbon balance error for the time step (gC/m**2)
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
         totcoln                   =>    col_ns%totcoln              , & ! Input:  [real(r8) (:) ]  (gN/m22)   total column nitrogen, incl veg and cpool
         dwt_prod10n_gain_grc      =>    grc_nf%dwt_prod10n_gain      , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_prod100n_gain_grc     =>    grc_nf%dwt_prod100n_gain     , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_conv_nflux_grc        =>    grc_nf%dwt_conv_nflux        , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_leaf_grc     =>    grc_nf%dwt_seedn_to_leaf     , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_deadstem_grc =>    grc_nf%dwt_seedn_to_deadstem , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         begnb_grc                 =>    grc_ns%begnb                , & ! Output: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         endnb_grc                 =>    grc_ns%endnb                , & ! Output: [real(r8) (:) ]  nitrogen mass, end of time step (gN/m2**2)
         errnb_grc                 =>    grc_ns%errnb                  & ! Output: [real(r8) (:) ]  nitrogen balance error for the time step (gN/m2**2)
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
         totcolp                   =>    col_ps%totcolp              , & ! Input:  [real(r8) (:) ]  (gP/m2)   total column phosphorus, incl veg and cpool
         dwt_prod10p_gain_grc      =>    grc_pf%dwt_prod10p_gain      , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_prod100p_gain_grc     =>    grc_pf%dwt_prod100p_gain     , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_conv_pflux_grc        =>    grc_pf%dwt_conv_pflux        , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_leaf_grc     =>    grc_pf%dwt_seedp_to_leaf     , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_deadstem_grc =>    grc_pf%dwt_seedp_to_deadstem , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         begpb_grc                 =>    grc_ps%begpb                , & ! Output: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         endpb_grc                 =>    grc_ps%endpb                , & ! Output: [real(r8) (:) ]  phosphorus mass, end of time step (gP/m2**2)
         errpb_grc                 =>    grc_ps%errpb                  & ! Output: [real(r8) (:) ]  phosphorus balance error for the time step (gP/m2**2)
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
