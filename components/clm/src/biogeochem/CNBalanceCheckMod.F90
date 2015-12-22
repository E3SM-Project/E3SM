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
  use clm_time_manager    , only : get_step_size,get_nstep
  use clm_varpar          , only : crop_prog
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnType          , only : col                
  use GridcellType        , only : grc

  use CNDecompCascadeConType , only : decomp_cascade_con
  use clm_varpar          , only: ndecomp_cascade_transitions
  use subgridAveMod       , only : p2c
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  ! bgc interface & pflotran:
  use clm_varctl          , only : use_pflotran, pf_cmode, pf_hmode


  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BeginCBalance
  public :: BeginNBalance
  public :: BeginPBalance
  public :: CBalanceCheck
  public :: NBalanceCheck
  public :: PBalanceCheck
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
  subroutine BeginPBalance(bounds, num_soilc, filter_soilc, &
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

    associate(                                         & 
         totcolp   => phosphorusstate_vars%totcolp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         !X.YANG - checking P balance problem, starting from VEGP 
         totpftp   => phosphorusstate_vars%totpftp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp   => phosphorusstate_vars%totsomp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp   => phosphorusstate_vars%cwdp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totlitp   => phosphorusstate_vars%totlitp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         sminp   => phosphorusstate_vars%sminp_col , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
 
         col_begpb => phosphorusstate_vars%begpb_col     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         )

      ! calculate beginning column-level phosphorus balance, for mass conservation check
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_begpb(c) = totcolp(c)
!         col_begpb(c) = totpftp(c)
!         col_begpb(c) = sminp(c)
!         col_begpb(c) = totsomp(c)+totlitp(c)+cwdp(c)
      end do

    end associate

  end subroutine BeginPBalance

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
         col_decompc_delta       =>    carbonflux_vars%externalc_to_decomp_delta_col, & ! Input: [real(r8) (:) ] (gC/m2/s) summarized net change of whole column C i/o to decomposing pool bwtn time-step
         
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
      
      if (.not. use_ed) then
         if (err_found) then
            c = err_index
            write(iulog,*)'column cbalance error = ', col_errcb(c), c
            write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
            write(iulog,*)'input=',col_cinputs*dt
            write(iulog,*)'output=',col_coutputs*dt
            write(iulog,*)'er=',er(c)*dt,carbonflux_vars%hr_col(c)*dt
            write(iulog,*)'fire=',col_fire_closs(c)*dt
            write(iulog,*)'dwt=',dwt_closs(c)*dt
            write(iulog,*)'product=',product_closs(c)*dt
            write(iulog,*)'hrv=',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'leach=',som_c_leached(c)*dt
            write(iulog,*)'begcb       = ',col_begcb(c)
            write(iulog,*)'endcb       = ',col_endcb(c),carbonstate_vars%totsomc_col(c)
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
    use tracer_varcon,  only : is_active_betr_bgc
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
         ! pflotran:
         col_decompn_delta   =>    nitrogenflux_vars%externaln_to_decomp_delta_col, & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column N i/o to decomposing pool bwtn time-step
         col_no3_delta       =>    nitrogenflux_vars%no3_net_transport_delta_col  , & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column NO3 leaching bwtn time-step

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

         if(is_active_betr_bgc)then
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         else
           if (.not. use_nitrif_denitrif) then
            col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
           else
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
           end if
         endif

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

            ! if not hydrology-coupled, NO3 leaching/runoff at previous time-step used as 'source' (out, -) to PFLOTRAN bgc
            if (.not.pf_hmode) col_errnb(c) = col_errnb(c) + col_no3_delta(c)*dt
            ! here is '+' adjustment. It says that the taking to PF no3 pools was more.
         end if

         if (abs(col_errnb(c)) > 1e-8_r8) then
            err_found = .true.
            err_index = c
         end if
      end do ! end of columns loop

      if (err_found) then
         c = err_index
         write(iulog,*)'column nbalance error = ', col_errnb(c), c, get_nstep()
         write(iulog,*)'Latdeg,Londeg         = ',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begnb                 = ',col_begnb(c)
         write(iulog,*)'endnb                 = ',col_endnb(c)
         write(iulog,*)'delta store           = ',col_endnb(c)-col_begnb(c)
         write(iulog,*)'input mass            = ',col_ninputs(c)*dt
         write(iulog,*)'output mass           = ',col_noutputs(c)*dt
         write(iulog,*)'net flux              = ',(col_ninputs(c)-col_noutputs(c))*dt
         write(iulog,*)'denit                 = ',denit(c)*dt
         write(iulog,*)'n2onit                = ',f_n2o_nit(c)*dt
         write(iulog,*)'no3 leach             = ', smin_no3_leached(c)*dt 
         write(iulog,*)'no3 runof             = ', smin_no3_runoff(c)*dt
         write(iulog,*)'ndep                  = ',ndep_to_sminn(c)*dt
         write(iulog,*)'nfix                  = ', nfix_to_sminn(c)*dt
         write(iulog,*)'nsup                  = ',supplement_to_sminn(c)*dt
         if(crop_prog) then
            write(iulog,*)'fertm                 = ',fert_to_sminn(c)*dt
            write(iulog,*)'soyfx                 = ',soyfixn_to_sminn(c)*dt
         endif
         write(iulog,*)'fire                  = ',col_fire_nloss(c)*dt
         write(iulog,*)'dwt                   = ',dwt_nloss(c)*dt 
         write(iulog,*)'prod                  = ',product_nloss(c)*dt
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine NBalanceCheck


  !-----------------------------------------------------------------------
  subroutine PBalanceCheck(bounds, &
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
    integer :: c,err_index,j,k  ! indices
    integer :: fc             ! lake filter indices
    logical :: err_found      ! error flag
    real(r8):: dt             ! radiation time step (seconds)

    real(r8) :: leafp_to_litter_col(bounds%begc:bounds%endc) 
    real(r8) :: frootp_to_litter_col(bounds%begc:bounds%endc) 
   real(r8):: flux_mineralization_col(bounds%begc:bounds%endc)   !!  local temperary variable
    !-----------------------------------------------------------------------

    associate(                                                             & 
         totcolp             =>    phosphorusstate_vars%totcolp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         supplement_to_sminp =>    phosphorusflux_vars%supplement_to_sminp_col , & ! Input:  [real(r8) (:)]  supplemental P supply (gP/m2/s)         
         sminp_leached       =>    phosphorusflux_vars%sminp_leached_col       , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         col_fire_ploss      =>    phosphorusflux_vars%fire_ploss_col      , & ! Input:  [real(r8) (:)]  total column-level fire P loss (gP/m2/s)
         dwt_ploss           =>    phosphorusflux_vars%dwt_ploss_col           , & ! Input:  [real(r8) (:)]  (gP/m2/s) total phosphorus loss from product pools and conversion
         product_ploss       =>    phosphorusflux_vars%product_ploss_col       , & ! Input:  [real(r8) (:)]  (gP/m2/s) total wood product phosphorus loss
         primp_to_labilep=>    phosphorusflux_vars%primp_to_labilep_col    , &
         secondp_to_occlp=>    phosphorusflux_vars%secondp_to_occlp_col    , &
         fert_p_to_sminp =>    phosphorusflux_vars%fert_p_to_sminp_col     , &
 
         col_pinputs         =>    phosphorusflux_vars%pinputs_col         , & ! Output: [real(r8) (:)]  column-level P inputs (gP/m2/s)         
         col_poutputs        =>    phosphorusflux_vars%poutputs_col        , & ! Output: [real(r8) (:)]  column-level P outputs (gP/m2/s)        
         col_begpb           =>    phosphorusstate_vars%begpb_col              , & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         col_endpb           =>    phosphorusstate_vars%endpb_col              , & ! Output: [real(r8) (:)]  phosphorus mass, end of time step (gP/m**2)
         col_errpb           =>    phosphorusstate_vars%errpb_col              ,  & ! Output: [real(r8) (:)]  phosphorus balance error for the timestep (gP/m**2)

         !X.YANG testing P Balance, from VEGP
         totpftp             =>    phosphorusstate_vars%totpftp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp             =>    phosphorusstate_vars%totsomp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp             =>    phosphorusstate_vars%cwdp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totlitp             =>    phosphorusstate_vars%totlitp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         sminp             =>    phosphorusstate_vars%sminp_col            , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         leafp_to_litter       =>    phosphorusflux_vars%leafp_to_litter_patch       , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         frootp_to_litter       =>   phosphorusflux_vars%frootp_to_litter_patch      , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         sminp_to_plant         =>   phosphorusflux_vars%sminp_to_plant_col          ,&
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool ,     &
         pf => phosphorusflux_vars                                               &
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      err_found = .false.

      call p2c(bounds,num_soilc,filter_soilc, &
           leafp_to_litter(bounds%begp:bounds%endp), &
           leafp_to_litter_col(bounds%begc:bounds%endc))
      call p2c(bounds,num_soilc,filter_soilc, &
           frootp_to_litter(bounds%begp:bounds%endp), &
           frootp_to_litter_col(bounds%begc:bounds%endc))

      !! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes
      !! - X.YANG
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
                                               pf%decomp_cascade_sminp_flux_col(c,k)*dt
               end do
         else
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                    flux_mineralization_col(c) = flux_mineralization_col(c) + &
                                               pf%decomp_cascade_sminp_flux_col(c,k)*dt

               end do
         endif
      end do



      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level phosphorus storage, for mass conservation check
!         col_endpb(c) = totpftp(c)
!         col_endpb(c) = sminp(c)
         col_endpb(c) = totcolp(c)
!         col_endpb(c) = totsomp(c)+cwdp(c)+totlitp(c)


         ! calculate total column-level inputs
         col_pinputs(c) = primp_to_labilep(c) + supplement_to_sminp(c)
!         col_pinputs(c) = leafp_to_litter_col(c)+frootp_to_litter_col(c)
!         col_pinputs(c) = sminp_to_plant(c)
!         col_pinputs(c) =  flux_mineralization_col(c)/dt
!         if (crop_prog) col_pinputs(c) = col_pinputs(c) + fert_p_to_sminp(c) 

         col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + col_fire_ploss(c) + dwt_ploss(c) + product_ploss(c)
!         col_poutputs(c) = leafp_to_litter_col(c)+frootp_to_litter_col(c)
!         col_poutputs(c) =  flux_mineralization_col(c)/dt
!          col_poutputs(c) = sminp_to_plant(c)

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
         write(iulog,*)'Latdeg,Londeg=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
         write(iulog,*)'begpb       = ',col_begpb(c)
         write(iulog,*)'endpb       = ',col_endpb(c)
         write(iulog,*)'delta store = ',col_endpb(c)-col_begpb(c)
         write(iulog,*)'input mass  = ',col_pinputs(c)*dt
         write(iulog,*)'output mass = ',col_poutputs(c)*dt
         write(iulog,*)'net flux    = ',(col_pinputs(c)-col_poutputs(c))*dt
         call endrun(msg=errMsg(__FILE__, __LINE__))
      end if

    end associate

  end subroutine PBalanceCheck

end module CNBalanceCheckMod
