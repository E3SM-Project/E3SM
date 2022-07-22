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
  use elm_varctl          , only : iulog, use_fates
  use clm_time_manager    , only : get_step_size,get_nstep
  use elm_varpar          , only : crop_prog
  use elm_varpar          , only : nlevdecomp
  use elm_varcon          , only : dzsoi_decomp
  use elm_varctl          , only : nu_com
  use elm_varctl          , only : ECA_Pconst_RGspin
  use spmdMod             , only : masterproc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use elm_varpar          , only: ndecomp_cascade_transitions
  use subgridAveMod       , only : p2c, c2g, unity
  ! soil erosion
  use elm_varctl          , only : use_erosion, ero_ccycle
  ! bgc interface & pflotran:
  use elm_varctl          , only : use_pflotran, pf_cmode, pf_hmode
  ! forest fertilization experiment
  use clm_time_manager    , only : get_curr_date
  use CNStateType         , only : fert_type , fert_continue, fert_dose, fert_start, fert_end
  use elm_varctl          , only : forest_fert_exp
  use pftvarcon           , only: noveg
  use elm_varctl          , only : NFIX_PTASE_plant
  use GridcellType        , only : grc_pp
  use GridcellDataType    , only : gridcell_carbon_state, gridcell_carbon_flux
  use GridcellDataType    , only : gridcell_nitrogen_state, gridcell_nitrogen_flux
  use GridcellDataType    , only : gridcell_phosphorus_state, gridcell_phosphorus_flux
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : column_carbon_state, column_carbon_flux
  use ColumnDataType      , only : column_nitrogen_state, column_nitrogen_flux
  use ColumnDataType      , only : column_phosphorus_state, column_phosphorus_flux
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_cf, veg_nf, veg_pf

  use timeinfoMod

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
  public :: BeginGridCBalance
  public :: GridCBalanceCheck
  public :: BeginGridNBalance
  public :: BeginGridPBalance
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
      end do

    end associate

  end subroutine BeginColCBalance

  !-----------------------------------------------------------------------
  subroutine BeginColNBalance(bounds, num_soilc, filter_soilc, &
       col_ns)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning nitrogen balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds          
    integer           , intent(in)    :: num_soilc       ! number of soil columns filter
    integer           , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_nitrogen_state) , intent(inout) :: col_ns
    !
    ! !LOCAL VARIABLES:
    integer :: c    ! indices
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
      end do

    end associate

  end subroutine BeginColNBalance

  !-----------------------------------------------------------------------
  subroutine BeginColPBalance(bounds, num_soilc, filter_soilc, &
       col_ps)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the beginning phosphorus balance for mass
    ! conservation checks.
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_phosphorus_state) , intent(inout) :: col_ps
    !
    ! !LOCAL VARIABLES:
    integer :: c     ! indices
    integer :: fc   ! lake filter indices
    !-----------------------------------------------------------------------

    associate(                                           &
         totcolp   => col_ps%totcolp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totpftp   => col_ps%totpftp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         totsomp   => col_ps%totsomp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         cwdp      => col_ps%cwdp    , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         totlitp   => col_ps%totlitp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg 
         sminp     => col_ps%sminp   , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
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
       col_cs, col_cf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform carbon mass conservation check for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_carbon_state) , intent(inout) :: col_cs
    type(column_carbon_flux)  , intent(in)    :: col_cf
    !
    ! !LOCAL VARIABLES:
    integer  :: c,err_index    ! indices
    integer  :: fc             ! lake filter indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    integer  :: nstep
    !-----------------------------------------------------------------------

    associate(                                                                           &
         totcolc                   =>    col_cs%totcolc                  , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         gpp                       =>    col_cf%gpp                       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) gross primary production
         er                        =>    col_cf%er                        , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs            =>    col_cf%fire_closs                , & ! Input:  [real(r8) (:) ]  (gC/m2/s) total column-level fire C loss
         col_prod1c_loss           =>    col_cf%prod1c_loss               , & ! Input:  [real(r8) (:) ]  (gC/m2/s) crop leafc harvested
         col_prod10c_loss          =>    col_cf%prod10c_loss              , & ! Input:  [real(r8) (:) ]  (gC/m2/s) 10-year wood C harvested
         col_prod100c_loss         =>    col_cf%prod100c_loss             , & ! Input:  [real(r8) (:) ]  (gC/m2/s) 100-year wood C harvested
         col_hrv_xsmrpool_to_atm   =>    col_cf%hrv_xsmrpool_to_atm       , & ! Input:  [real(r8) (:) ]  (gC/m2/s) excess MR pool harvest mortality
         som_c_leached             =>    col_cf%som_c_leached             , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss from vertical transport
         som_c_yield               =>    col_cf%somc_yield                , & ! Input:  [real(r8) (:) ]  (gC/m^2/s)total SOM C loss by erosion
         col_decompc_delta         =>    col_cf%externalc_to_decomp_delta , & ! Input:  [real(r8) (:) ]  (gC/m2/s) summarized net change of whole column C i/o to decomposing pool bwtn time-step
         col_cinputs               =>    col_cf%cinputs                   , & ! Output: [real(r8) (:)]  column-level C inputs (gC/m2/s)
         col_coutputs              =>    col_cf%coutputs                  , & ! Output: [real(r8) (:)]  column-level C outputs (gC/m2/s)
         col_begcb                 =>    col_cs%begcb                    , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         col_endcb                 =>    col_cs%endcb                    , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         col_errcb                 =>    col_cs%errcb                    ,  & ! Output: [real(r8) (:) ]  carbon balance error for the timestep (gC/m**2)
         hr                        =>    col_cf%hr                       , &  ! Input: heterotrophic respiration flux (gC/m2/s)
         litfall                   =>    col_cf%litfall )                     ! Input: carbon flux from litterfall (particularly fates) ( gc/m2/s)


      ! set time steps
      dt = real( get_step_size(), r8 )
      nstep = get_nstep()

      err_found = .false.
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate the total column-level carbon storage, for mass conservation check
         col_endcb(c) = totcolc(c)

         ! FATES also checks to see if input fluxes match
         ! a change in the total stock. So hwere we assume that
         ! the fates stocks are 0 for simplicity, and are only
         ! concerned that the changes in the soil stocks (including
         ! fragmented litter, match the fluxes in and out of fates
         if(use_fates) then

            col_cinputs  = litfall(c)

            col_coutputs = er(c)

         else

            ! calculate total column-level inputs
            col_cinputs(c) = gpp(c)

            ! calculate total column-level outputs
            ! er = ar + hr, col_fire_closs includes pft-level fire losses
            col_coutputs(c) = er(c) + col_fire_closs(c) + col_hrv_xsmrpool_to_atm(c)
         end if

         ! Wood product losses and crop export losses
         col_coutputs(c) = col_coutputs(c) + &
                 col_prod1c_loss(c) + col_prod10c_loss(c) + col_prod100c_loss(c)

         ! subtract leaching flux
         col_coutputs(c) = col_coutputs(c) - som_c_leached(c)


         ! add erosion flux
         if (ero_ccycle) then
            col_coutputs(c) = col_coutputs(c) + som_c_yield(c)
         end if


         ! calculate the total column-level carbon balance error for this time step
         col_errcb(c) = (col_cinputs(c) - col_coutputs(c))*dt - (col_endcb(c) - col_begcb(c))

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
#ifndef _OPENACC
         if (err_found .and. nstep > 1) then
            c = err_index
            write(iulog,*)'column cbalance error = ', col_errcb(c), c
            write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
            write(iulog,*)'input                 = ',col_cinputs(c)*dt
            write(iulog,*)'output                = ',col_coutputs(c)*dt
            write(iulog,*)'er                    = ',er(c)*dt,col_cf%hr(c)*dt
            write(iulog,*)'fire                  = ',col_fire_closs(c)*dt
            write(iulog,*)'hrv_to_atm            = ',col_hrv_xsmrpool_to_atm(c)*dt
            write(iulog,*)'leach                 = ',som_c_leached(c)*dt
            write(iulog,*)'begcb                 = ',col_begcb(c)
            write(iulog,*)'endcb                 = ',col_endcb(c),col_cs%totsomc(c)
            write(iulog,*)'totsomc               = ',col_cs%totsomc(c)
            write(iulog,*)'delta store           = ',col_endcb(c)-col_begcb(c)

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
         else
             if (masterproc .and. nstep < 2) then
                write(iulog,*) '--WARNING-- skipping CN balance check for first timestep'
             end if
         end if
#endif
      end if


    end associate

  end subroutine ColCBalanceCheck

  !-----------------------------------------------------------------------
  subroutine ColNBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       col_ns, col_nf)
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
    type(column_nitrogen_state) , intent(inout) :: col_ns
    type(column_nitrogen_flux)  , intent(inout) :: col_nf
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
    integer :: nstep
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
         plant_to_litter_nflux     =>    col_nf%plant_to_litter_nflux     , & ! Input                   flux of N from FATES litter into ELM
                                                                              !                         litter (gP/m2/s)
         col_prod1n_loss           =>    col_nf%prod1n_loss               , & ! Input:  [real(r8) (:)]  crop leafc harvested [gN/m2/s]
         col_prod10n_loss          =>    col_nf%prod10n_loss              , & ! Input:  [real(r8) (:)]  10-year wood product harvested [gN/m2/s]
         col_prod100n_loss         =>    col_nf%prod100n_loss             , & ! Input:  [real(r8) (:)]  100-year wood product harvestd [gN/m2/s]
         col_fire_nloss            =>    col_nf%fire_nloss                , & ! Input:  [real(r8) (:)]  total column-level fire N loss (gN/m2/s)
         som_n_leached             =>    col_nf%som_n_leached             , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport
         som_n_yield               =>    col_nf%somn_yield                , & ! Input:  [real(r8) (:)]  total SOM N loss by erosion
         supplement_to_plantn      =>    veg_nf%supplement_to_plantn      , &
         ! pflotran:
         col_decompn_delta         =>    col_nf%externaln_to_decomp_delta , & ! Input: [real(r8) (:) ] (gN/m2/s) summarized net change of whole column N i/o to decomposing pool bwtn time-step
         col_ninputs               =>    col_nf%ninputs                   , & ! Output: [real(r8) (:)]  column-level N inputs (gN/m2/s)
         col_noutputs              =>    col_nf%noutputs                  , & ! Output: [real(r8) (:)]  column-level N outputs (gN/m2/s)
         col_begnb                 =>    col_ns%begnb                    , & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         col_endnb                 =>    col_ns%endnb                    , & ! Output: [real(r8) (:)]  nitrogen mass, end of time step (gN/m**2)
         col_errnb                 =>    col_ns%errnb                    ,  & ! Output: [real(r8) (:)]  nitrogen balance error for the timestep (gN/m**2)
         sminn_to_plant            =>    col_nf%sminn_to_plant           )  ! nitrogen flux to plants [gN/m2/s]


      ! set time steps
      dt = dtime_mod
      nstep = nstep_mod
      kyr = year_curr; kmo = mon_curr; kda = day_curr; mcsec = secs_curr;

      err_found = .false.
      ! column loop
      do fc = 1,num_soilc
         c=filter_soilc(fc)

         ! calculate the total column-level nitrogen storage, for mass conservation check
         col_endnb(c) = totcoln(c)

         if(use_fates) then

            ! calculate total column-level inputs
            ! In FATES the nfix sent to the plant, if any, is not contained in
            ! sminn_to_plant(c), therefore, in neither case do we account for nfix_to_ecosysn
            if (NFIX_PTASE_plant) then
               col_ninputs(c) = ndep_to_sminn(c) + supplement_to_sminn(c) + nfix_to_sminn(c)
            else
               col_ninputs(c) = ndep_to_sminn(c) + supplement_to_sminn(c) + nfix_to_sminn(c)
            end if

            ! plant_to_litter_nflux is used by FATES to store
            ! column level fragmentation fluxes of nitrogen FATES litter to
            ! ELM litter

            col_ninputs(c) = col_ninputs(c) + plant_to_litter_nflux(c)

            ! calculate total column-level outputs
            col_noutputs(c) = denit(c) +  sminn_to_plant(c)

         else

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

            ! calculate total column-level outputs
            col_noutputs(c) = denit(c) + col_fire_nloss(c)

         end if

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
                 .and. fert_type(c) == 1 &
                 .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_ninputs(c) = col_ninputs(c) + fert_dose(c,kmo)/dt
            end if
         end if

         if (is_active_betr_bgc)then
            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)

            col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
         else

            col_noutputs(c) = col_noutputs(c) + f_n2o_nit(c)
            
            if(use_pflotran .and. pf_cmode) then
               ! inclusion of aq. NH4 transport by PFLOTRAN-bgc
               col_noutputs(c) = col_noutputs(c) + sminn_leached(c)
            else
               col_noutputs(c) = col_noutputs(c) + smin_no3_leached(c) + smin_no3_runoff(c)
            endif
            

         endif

         col_noutputs(c) = col_noutputs(c) + &
               col_prod1n_loss(c) + col_prod10n_loss(c) + col_prod100n_loss(c)

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

      if (err_found .and. nstep > 1) then
#ifndef _OPENACC
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
         write(iulog,*)'n_to_plant            = ',sminn_to_plant(c)*dt
         write(iulog,*)'plant_to_litter       = ',plant_to_litter_nflux(c)*dt
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
         call endrun(msg=errMsg(__FILE__, __LINE__))
#endif

      end if

    end associate

  end subroutine ColNBalanceCheck


  !-----------------------------------------------------------------------
  subroutine ColPBalanceCheck(bounds, &
       num_soilc, filter_soilc, &
       col_ps, col_pf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform phosphorus mass conservation check
    ! for column and pft
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_phosphorus_state) , intent(inout) :: col_ps
    type(column_phosphorus_flux)  , intent(inout) :: col_pf
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
    integer :: nstep
    !-----------------------------------------------------------------------

    associate(                                                                            &
         totcolp                   => col_ps%totcolp                  , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         supplement_to_sminp       => col_pf%supplement_to_sminp       , & ! Input:  [real(r8) (:)]  supplemental P supply (gP/m2/s)
         sminp_leached             => col_pf%sminp_leached             , & ! Input:  [real(r8) (:)]  soil mineral P pool loss to leaching (gP/m2/s)
         col_fire_ploss            => col_pf%fire_ploss                , & ! Input:  [real(r8) (:)]  total column-level fire P loss (gP/m2/s)
         biochem_pmin_to_plant     => col_pf%biochem_pmin_to_plant     , & ! Input:  total column level biochemical p flux to plant (gP/m2/s)
         primp_to_labilep          => col_pf%primp_to_labilep          , &
         secondp_to_occlp          => col_pf%secondp_to_occlp          , &
         fert_p_to_sminp           => col_pf%fert_p_to_sminp           , &
         som_p_yield               => col_pf%somp_yield                , & ! Input:  [real(r8) (:)]  SOM P pool loss by erosion (gP/m^2/s)
         labilep_yield             => col_pf%labilep_yield             , & ! Input:  [real(r8) (:)]  soil labile mineral P loss by erosion (gP/m^s/s)
         secondp_yield             => col_pf%secondp_yield             , & ! Input:  [real(r8) (:)]  soil secondary mineral P loss by erosion (gP/m^s/s)
         occlp_yield               => col_pf%occlp_yield               , & ! Input:  [real(r8) (:)]  soil occluded mineral P loss by erosion (gP/m^s/s)
         primp_yield               => col_pf%primp_yield               , & ! Input:  [real(r8) (:)]  soil primary mineral P loss by erosion (gP/m^s/s)
         supplement_to_plantp      => veg_pf%supplement_to_plantp          , &
         plant_to_litter_pflux     => col_pf%plant_to_litter_pflux     , & ! Input                   flux of P from FATES litter into ELM
                                                                           !                         litter (gP/m2/s)
         col_prod1p_loss           => col_pf%prod1p_loss               , & ! Input:  [real(r8) (:) ]  crop leafc harvested (gP/m2/s)
         col_prod10p_loss          => col_pf%prod10p_loss              , & ! Input:  [real(r8) (:) ]  10-yr wood product harvested (gP/m2/s)
         col_prod100p_loss         => col_pf%prod100p_loss             , & ! Input:  [real(r8) (:) ]  100-yr wood product harvested (gP/m2/s)
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
         cascade_receiver_pool     => decomp_cascade_con%cascade_receiver_pool   &
         )

      ! set time steps
      dt = dtime_mod
      kyr = year_curr; kmo = mon_curr; kda = day_curr; mcsec = secs_curr;

      err_found = .false.

      if(.not.use_fates)then
         call p2c(bounds,num_soilc,filter_soilc, &
              leafp_to_litter(bounds%begp:bounds%endp), &
              leafp_to_litter_col(bounds%begc:bounds%endc))
         call p2c(bounds,num_soilc,filter_soilc, &
              frootp_to_litter(bounds%begp:bounds%endp), &
              frootp_to_litter_col(bounds%begc:bounds%endc))
      end if

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

         if(use_fates) then

            col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + sminp_to_plant(c) + biochem_pmin_to_plant(c)
            ! plant_to_litter_nflux is used by FATES to store
            ! column level fragmentation fluxes of nitrogen FATES litter to
            ! ELM litter

            col_pinputs(c) = col_pinputs(c) + plant_to_litter_pflux(c)

         else
            do p = col_pp%pfti(c), col_pp%pftf(c)
               if (veg_pp%active(p) .and. (veg_pp%itype(p) .ne. noveg)) then
                  col_pinputs(c) = col_pinputs(c) + supplement_to_plantp(p) * veg_pp%wtcol(p)
               end if
            end do

            col_poutputs(c) = secondp_to_occlp(c) + sminp_leached(c) + col_fire_ploss(c)

         end if

         ! forest fertilization
         if (forest_fert_exp) then
            if ( ((fert_continue(c) == 1 .and. kyr > fert_start(c) .and. kyr <= fert_end(c)) .or.  kyr == fert_start(c)) &
               .and. fert_type(c) == 2 &
               .and. kda == 1  .and. mcsec == 1800) then ! fertilization assumed to occur at the begnining of each month
               col_pinputs(c) = col_pinputs(c) + fert_dose(c,kmo)/dt
             end if
         end if

         if ((nu_com .ne. 'RD') .and. ECA_Pconst_RGspin) then
            do j = 1, nlevdecomp
               col_poutputs(c) = col_poutputs(c) + &
                  (col_ps%solutionp_vr_cur(c,j) -  col_ps%solutionp_vr_prev(c,j)  + &
                  col_ps%labilep_vr_cur(c,j) -  col_ps%labilep_vr_prev(c,j) + &
                  col_ps%secondp_vr_cur(c,j) - col_ps%secondp_vr_prev(c,j) ) * dzsoi_decomp(j)/dt
            end do
         end if

         col_poutputs(c) = col_poutputs(c) + &
            col_prod1p_loss(c) + col_prod10p_loss(c) + col_prod100p_loss(c)

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


      if (err_found .and. nstep>1) then
#ifndef _OPENACC
         c = err_index
         write(iulog,*)'column pbalance error = ', col_errpb(c), c
         write(iulog,*)'Latdeg,Londeg=',grc_pp%latdeg(col_pp%gridcell(c)),grc_pp%londeg(col_pp%gridcell(c))
         write(iulog,*)'begpb       = ',col_begpb(c)
         write(iulog,*)'endpb       = ',col_endpb(c)
         write(iulog,*)'delta store = ',col_endpb(c)-col_begpb(c)
         write(iulog,*)'input mass  = ',col_pinputs(c)*dt
         write(iulog,*)'output mass = ',col_poutputs(c)*dt
         write(iulog,*)'net flux    = ',(col_pinputs(c)-col_poutputs(c))*dt
         write(iulog,*)'ptol_pflux: = ',plant_to_litter_pflux(c)*dt
         write(iulog,*)'sminp_to_plant: = ',sminp_to_plant(c)*dt
         write(iulog,*)'primp_to_labilep = ',primp_to_labilep(c)*dt
         write(iulog,*)'supplement_to_sminp = ',supplement_to_sminp(c)*dt
         write(iulog,*)'secondp_to_occlp = ',secondp_to_occlp(c)*dt
         write(iulog,*)'sminp_leached = ',sminp_leached(c)*dt
         if(use_fates)then
            write(iulog,*) 'plant_to_litter_flux = ',plant_to_litter_pflux(c)*dt
            write(iulog,*) 'biochem_pmin_to_plant = ',biochem_pmin_to_plant(c)*dt
         end if
         if (ero_ccycle) then
            write(iulog,*)'SOP erosion = ',som_p_yield(c)*dt
            write(iulog,*)'SIP erosion = ',(labilep_yield(c)+secondp_yield(c)+occlp_yield(c)+primp_yield(c))*dt
         end if
         call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
      end if

    end associate

  end subroutine ColPBalanceCheck

  !-----------------------------------------------------------------------
  subroutine BeginGridCBalance(bounds, col_cs, grc_cs)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning carbon balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)  :: bounds
    type(column_carbon_state)  , intent(in)  :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    !-----------------------------------------------------------------------

    associate(                                                    &
         totcolc               =>  col_cs%totcolc               , & ! Input:  [real(r8) (:)] (gC/m2) total column carbon, incl veg and cpool
         totpftc               =>  col_cs%totpftc               , & ! Input:  [real(r8) (:)] (gC/m2) patch-level carbon aggregated to column level, incl veg and cpool
         cwdc                  =>  col_cs%cwdc                  , & ! Input:  [real(r8) (:)] (gC/m2) total column coarse woody debris carbon,
         totsomc               =>  col_cs%totsomc               , & ! Input:  [real(r8) (:)] (gC/m2) total column soil organic matter carbon
         totlitc               =>  col_cs%totlitc               , & ! Input:  [real(r8) (:)] (gC/m2) total column litter carbon
         totprodc              =>  col_cs%totprodc              , & ! Input:  [real(r8) (:)] (gC/m2) total column wood product carbon
         ctrunc                =>  col_cs%ctrunc                , & ! Input:  [real(r8) (:)] (gC/m2) total column truncation carbon sink
         cropseedc_deficit     =>  col_cs%cropseedc_deficit     , & ! Input:  [real(r8) (:)] (gC/m2) column carbon pool for seeding new growth
         begcb_grc             =>  grc_cs%begcb                 , & ! Output: [real(r8) (:)]  carbon mass, beginning of time step (gC/m**2)
         beg_totc              =>  grc_cs%beg_totc              , & ! Output: [real(r8) (:)] (gC/m2) total column carbon, incl veg and cpool
         beg_totpftc           =>  grc_cs%beg_totpftc           , & ! Output: [real(r8) (:)] (gC/m2) patch-level carbon aggregated to column level, incl veg and cpool
         beg_cwdc              =>  grc_cs%beg_cwdc              , & ! Output: [real(r8) (:)] (gC/m2) total column coarse woody debris carbon
         beg_totsomc           =>  grc_cs%beg_totsomc           , & ! Output: [real(r8) (:)] (gC/m2) total column soil organic matter carbon
         beg_totlitc           =>  grc_cs%beg_totlitc           , & ! Output: [real(r8) (:)] (gC/m2) total column litter carbon
         beg_totprodc          =>  grc_cs%beg_totprodc          , & ! Output: [real(r8) (:)] (gC/m2) total column wood product carbon
         beg_ctrunc            =>  grc_cs%beg_ctrunc            , & ! Output: [real(r8) (:)] (gC/m2) total column truncation carbon sink
         beg_cropseedc_deficit =>  grc_cs%beg_cropseedc_deficit   & ! Output: [real(r8) (:)] (gC/m2) column carbon pool for seeding new growth
         )

      call c2g(bounds, totcolc(bounds%begc:bounds%endc), begcb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, totcolc(bounds%begc:bounds%endc), beg_totc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, totpftc(bounds%begc:bounds%endc), beg_totpftc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, cwdc(bounds%begc:bounds%endc), beg_cwdc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, totlitc(bounds%begc:bounds%endc), beg_totlitc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, totsomc(bounds%begc:bounds%endc), beg_totsomc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, totprodc(bounds%begc:bounds%endc), beg_totprodc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, ctrunc(bounds%begc:bounds%endc), beg_ctrunc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      call c2g(bounds, cropseedc_deficit(bounds%begc:bounds%endc), beg_cropseedc_deficit(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

    end associate

  end subroutine BeginGridCBalance

  !-----------------------------------------------------------------------
  subroutine GridCBalanceCheck(bounds, col_cs, col_cf, grc_cs, grc_cf)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning carbon balance for mass conservation checks
    ! at grid cell level
    !
    use GridcellDataType, only : gridcell_carbon_flux
    use ColumnDataType  , only : column_carbon_flux
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    type(column_carbon_state)  , intent(in)    :: col_cs
    type(column_carbon_flux)   , intent(in)    :: col_cf
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    type(gridcell_carbon_flux) , intent(inout) :: grc_cf
    !
    integer             :: g, nstep
    real(r8)            :: dt
    real(r8), parameter :: error_tol = 1.e-8_r8
    !-----------------------------------------------------------------------

    associate(                                                       &
         col_totc                  => col_cs%totcolc               , & ! Input:  [real(r8) (:) ] (gC/m2) total column carbon, incl veg and cpool
         col_totpftc               => col_cs%totpftc               , & ! Input:  [real(r8) (:) ] (gC/m2) patch-level carbon aggregated to column level, incl veg and cpool
         col_cwdc                  => col_cs%cwdc                  , & ! Input:  [real(r8) (:) ] (gC/m2) total column coarse woody debris carbon,
         col_totsomc               => col_cs%totsomc               , & ! Input:  [real(r8) (:) ] (gC/m2) total column soil organic matter carbon
         col_totlitc               => col_cs%totlitc               , & ! Input:  [real(r8) (:) ] (gC/m2) total column litter carbon
         col_totprodc              => col_cs%totprodc              , & ! Input:  [real(r8) (:) ] (gC/m2) total column wood product carbon
         col_ctrunc                => col_cs%ctrunc                , & ! Input:  [real(r8) (:) ] (gC/m2) total column truncation carbon sink
         col_cropseedc_deficit     => col_cs%cropseedc_deficit     , & ! Input:  [real(r8) (:) ] (gC/m2) column carbon pool for seeding new growth
         col_gpp                   => col_cf%gpp                   , & ! Input:  [real(r8) (:) ] (gC/m2/s) gross primary production
         col_er                    => col_cf%er                    , & ! Input:  [real(r8) (:) ] (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         col_fire_closs            => col_cf%fire_closs            , & ! Input:  [real(r8) (:) ] (gC/m2/s) total column-level fire C loss
         col_prod1c_loss           => col_cf%prod1c_loss           , & ! Input:  [real(r8) (:) ] (gC/m2/s) crop leafc harvested
         col_prod10c_loss          => col_cf%prod10c_loss          , & ! Input:  [real(r8) (:) ] (gC/m2/s) 10-year wood C harvested
         col_prod100c_loss         => col_cf%prod100c_loss         , & ! Input:  [real(r8) (:) ] (gC/m2/s) 100-year wood C harvested 
         col_hrv_xsmrpool_to_atm   => col_cf%hrv_xsmrpool_to_atm   , & ! Input:  [real(r8) (:) ] (gC/m2/s) excess MR pool harvest mortality
         col_som_c_leached         => col_cf%som_c_leached         , & ! Input:  [real(r8) (:) ] (gC/m^2/s)total SOM C loss from vertical transport
         col_som_c_yield           => col_cf%somc_yield            , & ! Input:  [real(r8) (:) ] (gC/m^2/s)total SOM C loss by erosion
         col_cinputs               => col_cf%cinputs               , & ! Input:  [real(r8) (:) ] (gC/m2/s) column-level C inputs
         col_coutputs              => col_cf%coutputs              , & ! Input:  [real(r8) (:) ] (gC/m2/s) column-level C outputs
         beg_totc                  => grc_cs%beg_totc              , & ! Input:  [real(r8) (:) ] (gC/m2) total column carbon, incl veg and cpool
         end_totc                  => grc_cs%end_totc              , & ! Output: [real(r8) (:) ] (gC/m2) total column carbon, incl veg and cpool
         end_totpftc               => grc_cs%end_totpftc           , & ! Output: [real(r8) (:) ] (gC/m2) total column carbon, incl veg and cpool
         end_cwdc                  => grc_cs%end_cwdc              , & ! Output: [real(r8) (:) ] (gC/m2) total column coarse woody debris carbon
         end_totsomc               => grc_cs%end_totsomc           , & ! Output: [real(r8) (:) ] (gC/m2) total column soil organic matter carbon
         end_totlitc               => grc_cs%end_totlitc           , & ! Output: [real(r8) (:) ] (gC/m2) total column litter carbon
         end_totprodc              => grc_cs%end_totprodc          , & ! Output: [real(r8) (:) ] (gC/m2) total column wood product carbon
         end_ctrunc                => grc_cs%end_ctrunc            , & ! Output: [real(r8) (:) ] (gC/m2) total column truncation carbon sink
         end_cropseedc_deficit     => grc_cs%end_cropseedc_deficit , & ! Output: [real(r8) (:) ] (gC/m2) column carbon pool for seeding new growth
         grc_errcb                 => grc_cs%errcb                 , & ! Output: [real(r8) (:) ] carbon balance error for the timestep (gC/m**2)
         grc_dwt_conv_cflux        => grc_cf%dwt_conv_cflux        , & !Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         grc_dwt_seedc_to_leaf     => grc_cf%dwt_seedc_to_leaf     , & !Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         grc_dwt_seedc_to_deadstem => grc_cf%dwt_seedc_to_deadstem , & !Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         grc_gpp                   => grc_cf%gpp                   , & ! Output: [real(r8) (:) ] (gC/m2/s) gross primary production
         grc_er                    => grc_cf%er                    , & ! Output: [real(r8) (:) ] (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         grc_fire_closs            => grc_cf%fire_closs            , & ! Output: [real(r8) (:) ] (gC/m2/s) total column-level fire C loss
         grc_prod1c_loss           => grc_cf%prod1_loss            , & ! Output: [real(r8) (:) ] (gC/m2/s) crop leafc harvested
         grc_prod10c_loss          => grc_cf%prod10_loss           , & ! Output: [real(r8) (:) ] (gC/m2/s) 10-year wood C harvested
         grc_prod100c_loss         => grc_cf%prod100_loss          , & ! Output: [real(r8) (:) ] (gC/m2/s) 100-year wood C harvested 
         grc_hrv_xsmrpool_to_atm   => grc_cf%hrv_xsmrpool_to_atm   , & ! Output: [real(r8) (:) ] (gC/m2/s) excess MR pool harvest mortality
         grc_som_c_leached         => grc_cf%som_c_leached         , & ! Output: [real(r8) (:) ] (gC/m^2/s)total SOM C loss from vertical transport
         grc_som_c_yield           => grc_cf%somc_yield            , & ! Output: [real(r8) (:) ] (gC/m^2/s)total SOM C loss by erosion
         grc_cinputs               => grc_cf%cinputs               , & ! Output: [real(r8) (:) ] (gC/m2/s) column-level C inputs
         grc_coutputs              => grc_cf%coutputs                & ! Output: [real(r8) (:) ] (gC/m2/s) column-level C outputs
         )

      ! c2g states
      call c2g(bounds, col_totc(bounds%begc:bounds%endc), end_totc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_totpftc(bounds%begc:bounds%endc), end_totpftc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_cwdc(bounds%begc:bounds%endc), end_cwdc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_totlitc(bounds%begc:bounds%endc), end_totlitc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_totsomc(bounds%begc:bounds%endc), end_totsomc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_totprodc(bounds%begc:bounds%endc), end_totprodc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_ctrunc(bounds%begc:bounds%endc), end_ctrunc(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_cropseedc_deficit(bounds%begc:bounds%endc), end_cropseedc_deficit(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      ! c2g fluxes
      call c2g(bounds, col_gpp(bounds%begc:bounds%endc), grc_gpp(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_er(bounds%begc:bounds%endc), grc_er(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_fire_closs(bounds%begc:bounds%endc), grc_fire_closs(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_prod1c_loss(bounds%begc:bounds%endc), grc_prod1c_loss(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_prod10c_loss(bounds%begc:bounds%endc), grc_prod10c_loss(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_prod100c_loss(bounds%begc:bounds%endc), grc_prod100c_loss(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_hrv_xsmrpool_to_atm(bounds%begc:bounds%endc), grc_hrv_xsmrpool_to_atm(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_som_c_leached(bounds%begc:bounds%endc), grc_som_c_leached(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')
      call c2g(bounds, col_som_c_yield(bounds%begc:bounds%endc), grc_som_c_yield(bounds%begg:bounds%endg), &
               c2l_scale_type = 'unity', l2g_scale_type = 'unity')

      dt = real( get_step_size(), r8 )
      nstep = get_nstep()

      do g = bounds%begg, bounds%endg
         grc_cinputs(g) = grc_gpp(g) + grc_dwt_seedc_to_leaf(g) + grc_dwt_seedc_to_deadstem(g)

         grc_coutputs(g) = grc_er(g) + grc_fire_closs(g) + grc_hrv_xsmrpool_to_atm(g) + &
              grc_prod1c_loss(g) + grc_prod10c_loss(g) + grc_prod100c_loss(g) - grc_som_c_leached(g) + &
              grc_dwt_conv_cflux(g)

         if (ero_ccycle) then
            grc_coutputs(g) = grc_coutputs(g) + grc_som_c_yield(g)
         end if

         grc_errcb(g) = (grc_cinputs(g) - grc_coutputs(g))*dt - (end_totc(g) - beg_totc(g))

         if (grc_errcb(g) > error_tol .and. nstep > 1) then
            write(iulog,*)'grid cbalance error = ', grc_errcb(g), g
            write(iulog,*)'Latdeg,Londeg       = ', grc_pp%latdeg(g), grc_pp%londeg(g)
            write(iulog,*)'input               = ', grc_cinputs(g)*dt
            write(iulog,*)'output              = ', grc_coutputs(g)*dt
            write(iulog,*)'er                  = ', grc_er(g)*dt
            write(iulog,*)'fire                = ', grc_fire_closs(g)*dt
            write(iulog,*)'hrv_to_atm          = ', grc_hrv_xsmrpool_to_atm(g)*dt
            write(iulog,*)'leach               = ', grc_som_c_leached(g)*dt

            if (ero_ccycle) then
               write(iulog,*)'erosion             = ',grc_som_c_yield(g)*dt
            end if

            write(iulog,*)'begcb               = ', beg_totc(g)
            write(iulog,*)'endcb               = ', end_totc(g)
            write(iulog,*)'delta store         = ', end_totc(g) - beg_totc(g)

            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if
      end do

    end associate

  end subroutine GridCBalanceCheck
 
  !-----------------------------------------------------------------------
  subroutine BeginGridNBalance(bounds, col_ns, grc_ns)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning nitrogen balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)             , intent(in)    :: bounds
    type(column_nitrogen_state)   , intent(in)    :: col_ns
    type(gridcell_nitrogen_state) , intent(inout) :: grc_ns
    !-----------------------------------------------------------------------

    associate(                                         &
         totcoln   => col_ns%totcoln , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg
         begnb_grc => grc_ns%begnb     & ! Output: [real(r8) (:)]  nitrogen mass, beginning of time step (gN/m**2)
         )

      call c2g(bounds, totcoln(bounds%begc:bounds%endc), begnb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', l2g_scale_type = 'unity')

    end associate

  end subroutine BeginGridNBalance

  !-----------------------------------------------------------------------
  subroutine BeginGridPBalance(bounds, col_ps, grc_ps)
    !
    ! !DESCRIPTION:
    ! Calculate the beginning phosphorus balance for mass conservation checks
    ! at grid cell level
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    type(column_phosphorus_state)   , intent(in)    :: col_ps
    type(gridcell_phosphorus_state) , intent(inout) :: grc_ps
    !
    !-----------------------------------------------------------------------

    associate(                                           &
         totcolp   => col_ps%totcolp , & ! Input:  [real(r8) (:)]  (gP/m2) total column phosphorus, incl veg
         begpb_grc => grc_ps%begpb     & ! Output: [real(r8) (:)]  phosphorus mass, beginning of time step (gP/m**2)
         )

      call c2g(bounds, totcolp(bounds%begc:bounds%endc), begpb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = 'unity', l2g_scale_type = 'unity')


    end associate

  end subroutine BeginGridPBalance

  !-----------------------------------------------------------------------
  subroutine EndGridCBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, col_cs, grc_cs, grc_cf)
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
    type(gridcell_carbon_flux) , intent(inout) :: grc_cf
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                                       &
         totcolc                   =>    col_cs%totcolc              , & ! Input:  [real(r8) (:) ]  (gC/m2)   total column carbon, incl veg and cpool
         dwt_conv_cflux_grc        =>    grc_cf%dwt_conv_cflux        , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_leaf_grc     =>    grc_cf%dwt_seedc_to_leaf     , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         dwt_seedc_to_deadstem_grc =>    grc_cf%dwt_seedc_to_deadstem , & ! Input: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         grc_cinputs               =>    grc_cf%cinputs               , & ! Output: [real(r8) (:)]  grid-level C inputs (gC/m2/s)
         grc_coutputs              =>    grc_cf%coutputs              , & ! Output: [real(r8) (:)]  grid-level C outputs (gC/m2/s)
         begcb_grc                 =>    grc_cs%begcb                , & ! Output: [real(r8) (:) ]  carbon mass, beginning of time step (gC/m**2)
         endcb_grc                 =>    grc_cs%endcb                , & ! Output: [real(r8) (:) ]  carbon mass, end of time step (gC/m**2)
         errcb_grc                 =>    grc_cs%errcb                  & ! Output: [real(r8) (:) ]  carbon balance error for the time step (gC/m**2)
         )

      ! set time steps
      dt = dtime_mod

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcolc(bounds%begc:bounds%endc), &
           garr = endcb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = unity, &
           l2g_scale_type = unity)

      do g = bounds%begg, bounds%endg
         endcb_grc(g) = endcb_grc(g)

         grc_cinputs(g) = &
              dwt_seedc_to_leaf_grc(g)     + &
              dwt_seedc_to_deadstem_grc(g)

         grc_coutputs(g) = &
              dwt_conv_cflux_grc(g)


         errcb_grc(g) = (grc_cinputs(g) - grc_coutputs(g))*dt - (endcb_grc(g) - begcb_grc(g))

         ! check for significant errors
         if (abs(errcb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
#ifndef _OPENACC
         g = err_index
         write(iulog,*)'Grid cbalance error   = ',errcb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_cinputs(g)*dt
         write(iulog,*)'output                = ',grc_coutputs(g)*dt
         write(iulog,*)'error                 = ',errcb_grc(g)*dt
         write(iulog,*)'begcb                 = ',begcb_grc(g)
         write(iulog,*)'endcb                 = ',endcb_grc(g)
         write(iulog,*)'delta store           = ',endcb_grc(g)-begcb_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
      end if

    end associate

  end subroutine EndGridCBalanceAfterDynSubgridDriver

  !-----------------------------------------------------------------------
  subroutine EndGridNBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, col_ns, grc_ns, grc_nf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform nitrogen mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)             , intent(in)    :: bounds          
    integer                       , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                       , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_nitrogen_state)   , intent(inout) :: col_ns
    type(gridcell_nitrogen_state) , intent(inout) :: grc_ns
    type(gridcell_nitrogen_flux)  , intent(inout) :: grc_nf
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                          &
         totcoln                   =>    col_ns%totcoln               , & ! Input:  [real(r8) (:) ]  (gN/m22)   total column nitrogen, incl veg and cpool
         dwt_conv_nflux_grc        =>    grc_nf%dwt_conv_nflux        , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_leaf_grc     =>    grc_nf%dwt_seedn_to_leaf     , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         dwt_seedn_to_deadstem_grc =>    grc_nf%dwt_seedn_to_deadstem , & ! Input: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         grc_ninputs               =>    grc_nf%ninputs               , & ! Output: [real(r8) (:)]  grid-level N inputs (gN/m2/s)
         grc_noutputs              =>    grc_nf%noutputs              , & ! Output: [real(r8) (:)]  grid-level N outputs (gN/m2/s)
         begnb_grc                 =>    grc_ns%begnb                 , & ! Output: [real(r8) (:) ]  nitrogen mass, beginning of time step (gN/m2**2)
         endnb_grc                 =>    grc_ns%endnb                 , & ! Output: [real(r8) (:) ]  nitrogen mass, end of time step (gN/m2**2)
         errnb_grc                 =>    grc_ns%errnb                   & ! Output: [real(r8) (:) ]  nitrogen balance error for the time step (gN/m2**2)
         )

      ! set time steps
      dt = dtime_mod

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcoln(bounds%begc:bounds%endc), &
           garr = endnb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = unity, &
           l2g_scale_type = unity)

      do g = bounds%begg, bounds%endg
         endnb_grc(g) = endnb_grc(g)

         grc_ninputs(g) = &
              dwt_seedn_to_leaf_grc(g)     + &
              dwt_seedn_to_deadstem_grc(g)

         grc_noutputs(g) = &
              dwt_conv_nflux_grc(g)


         errnb_grc(g) = (grc_ninputs(g) - grc_noutputs(g))*dt - (endnb_grc(g) - begnb_grc(g))

         ! check for significant errors
         if (abs(errnb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
#ifndef _OPENACC
         g = err_index
         write(iulog,*)'Grid nbalance error   = ',errnb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_ninputs(g)*dt
         write(iulog,*)'output                = ',grc_noutputs(g)*dt
         write(iulog,*)'error                 = ',errnb_grc(g)*dt
         write(iulog,*)'begnb                 = ',begnb_grc(g)
         write(iulog,*)'endnb                 = ',endnb_grc(g)
         write(iulog,*)'delta store           = ',endnb_grc(g)-begnb_grc(g)
         write(iulog,*)''
         write(iulog,*)'dwt_conv                ',dwt_conv_nflux_grc(g)
         write(iulog,*)''
         write(iulog,*)'dwt_seedn_leaf          ',dwt_seedn_to_leaf_grc(g)
         write(iulog,*)'dwt_seedn_deadstem      ',dwt_seedn_to_deadstem_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
      end if

    end associate

  end subroutine EndGridNBalanceAfterDynSubgridDriver

  !-----------------------------------------------------------------------

  subroutine EndGridPBalanceAfterDynSubgridDriver(bounds, &
       num_soilc, filter_soilc, col_ps, grc_ps, grc_pf)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform phosphorus mass conservation check
    ! at grid level after dynamic subgrid driver has been called
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds          
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_phosphorus_state)   , intent(inout) :: col_ps
    type(gridcell_phosphorus_state) , intent(inout) :: grc_ps
    type(gridcell_phosphorus_flux)  , intent(inout) :: grc_pf
    !
    ! !LOCAL VARIABLES:
    integer  :: g,err_index    ! indices
    logical  :: err_found      ! error flag
    real(r8) :: dt             ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                                          &
         totcolp                   =>    col_ps%totcolp               , & ! Input:  [real(r8) (:) ]  (gP/m2)   total column phosphorus, incl veg and cpool
         dwt_conv_pflux_grc        =>    grc_pf%dwt_conv_pflux        , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_leaf_grc     =>    grc_pf%dwt_seedp_to_leaf     , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         dwt_seedp_to_deadstem_grc =>    grc_pf%dwt_seedp_to_deadstem , & ! Input: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         grc_pinputs               =>    grc_pf%pinputs               , & ! Output: [real(r8) (:)]  grid-level P inputs (gP/m2/s)
         grc_poutputs              =>    grc_pf%poutputs              , & ! Output: [real(r8) (:)]  grid-level P outputs (gP/m2/s)
         begpb_grc                 =>    grc_ps%begpb                 , & ! Output: [real(r8) (:) ]  phosphorus mass, beginning of time step (gP/m2**2)
         endpb_grc                 =>    grc_ps%endpb                 , & ! Output: [real(r8) (:) ]  phosphorus mass, end of time step (gP/m2**2)
         errpb_grc                 =>    grc_ps%errpb                   & ! Output: [real(r8) (:) ]  phosphorus balance error for the time step (gP/m2**2)
         )

      ! set time steps
      dt = dtime_mod

      err_found = .false.

      call c2g( bounds = bounds, &
           carr = totcolp(bounds%begc:bounds%endc), &
           garr = endpb_grc(bounds%begg:bounds%endg), &
           c2l_scale_type = unity, &
           l2g_scale_type = unity)

      do g = bounds%begg, bounds%endg
         endpb_grc(g) = endpb_grc(g)

         grc_pinputs(g) = &
              dwt_seedp_to_leaf_grc(g)     + &
              dwt_seedp_to_deadstem_grc(g)

         grc_poutputs(g) = &
              dwt_conv_pflux_grc(g)


         errpb_grc(g) = (grc_pinputs(g) - grc_poutputs(g))*dt - (endpb_grc(g) - begpb_grc(g))

         ! check for significant errors
         if (abs(errpb_grc(g)) > balance_check_tolerance) then
            err_found = .true.
            err_index = g
         end if

      end do

      if (err_found) then
#ifndef _OPENACC
         g = err_index
         write(iulog,*)'Grid pbalance error   = ',errpb_grc(g), g
         write(iulog,*)'Latdeg,Londeg         = ',grc_pp%latdeg(g),grc_pp%londeg(g)
         write(iulog,*)'input                 = ',grc_pinputs(g)*dt
         write(iulog,*)'output                = ',grc_poutputs(g)*dt
         write(iulog,*)'error                 = ',errpb_grc(g)*dt
         write(iulog,*)'begpb                 = ',begpb_grc(g)
         write(iulog,*)'endpb                 = ',endpb_grc(g)
         write(iulog,*)'delta store           = ',endpb_grc(g)-begpb_grc(g)
         call endrun(msg=errMsg(__FILE__, __LINE__))
#endif
      end if

    end associate

  end subroutine EndGridPBalanceAfterDynSubgridDriver

end module EcosystemBalanceCheckMod
