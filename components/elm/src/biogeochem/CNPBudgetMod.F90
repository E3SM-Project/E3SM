module CNPBudgetMod
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_sys_mod         , only : shr_sys_abort
  use decompMod           , only : bounds_type
  use abortutils          , only : endrun
  use elm_varctl          , only : iulog
  use atm2lndType         , only : atm2lnd_type
  use lnd2atmType         , only : lnd2atm_type
  use spmdMod             , only : masterproc
  use GridcellDataType    , only : grc_ws
  use ColumnDataType      , only : col_ws
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use ColumnDataType      , only : column_carbon_state, col_cf 
  use ColumnDataType      , only : col_ns, col_nf, col_ps, col_pf 
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use GridcellDataType    , only : gridcell_carbon_state, gridcell_carbon_flux

  implicit none
  save
  private

  public :: CNPBudget_Reset
  public :: CNPBudget_Run
  public :: CNPBudget_Accum
  public :: CNPBudget_Print
  public :: CNPBudget_Restart
  public :: CNPBudget_SetBeginningMonthlyStates
  public :: CNPBudget_SetEndingMonthlyStates

  integer, parameter :: carbon_budget     = 1
  integer, parameter :: nitrogen_budget   = 2
  integer, parameter :: phosphorus_budget = 3

  !--- F for flux ---

  ! C inputs
  integer, parameter :: f_gpp                  = 1

  ! C outputs
  integer, parameter :: f_er                   =  2
  integer, parameter :: f_fire_closs           =  3
  integer, parameter :: f_hrv_xsmrpool_to_atm  =  4
  integer, parameter :: f_prod1c_loss          =  5
  integer, parameter :: f_prod10c_loss         =  6
  integer, parameter :: f_prod100c_loss        =  7
  integer, parameter :: f_som_c_leached        =  8
  integer, parameter :: f_som_c_yield          =  9
  integer, parameter :: f_dwt_conv_cflux       = 10
  integer, parameter :: f_dwt_seedc_to_leaf    = 11
  integer, parameter :: f_dwt_seedc_to_deadstem= 12

  integer, parameter, public :: c_f_size = f_dwt_seedc_to_deadstem

  character(len=51), parameter :: c_f_name(c_f_size) = &
       (/&
       '                              gross primary product', &
       '                              ecosystem respiration', &
       '                                        fire C loss', &
       '                   excess MR pool harvest mortality', &
       '        decomposition loss from 1-year product pool', &
       '       decomposition loss from 10-year product pool', &
       '      decomposition loss from 100-year product pool', &
       '                 SOM C loss from vertical transport', &
       '                                         SOM C loss', &
       '          flux to atmosphere due to dynamic weights', &
       '         seed source to leaf due to dynamic weights', &
       '   seed source to dead steam due to dynamic weights'  &
       /)
       
       
  ! N inputs
  integer, parameter :: f_ndep_to_sminn        = 13
  integer, parameter :: f_nfix_to_ecosysn      = 14
  integer, parameter :: f_nfix_to_sminn        = 15
  integer, parameter :: f_supplement_to_sminn  = 16
  integer, parameter :: f_fert_to_sminn        = 17
  integer, parameter :: f_soyfixn_to_sminn     = 18
  integer, parameter :: f_supplement_to_plantn = 19
  integer, parameter :: f_nfert_dose           = 20

  ! N outputs
  integer, parameter :: f_denit                = 21
  integer, parameter :: f_fire_ploss           = 22
  integer, parameter :: f_n2o_nit              = 23
  integer, parameter :: f_smin_no3_leached     = 24
  integer, parameter :: f_smin_no3_runoff      = 25
  integer, parameter :: f_sminn_leached        = 26
  integer, parameter :: f_col_prod1n_loss      = 27
  integer, parameter :: f_col_prod10n_loss     = 28
  integer, parameter :: f_col_prod100n_loss    = 29
  integer, parameter :: f_som_n_leached        = 30
  integer, parameter :: f_som_n_yield          = 31

  integer, parameter, public :: n_f_size = f_som_n_yield - f_dwt_seedc_to_deadstem

  character(len=39), parameter :: n_f_name(n_f_size) = &
       (/&
       '                      atm. N deposition', &
       '                total nitrogen fixation', &
       '        symbiotic/asymbiotic N fixation', &
       '         vert-int supplemental N supply', &
       '                           fertilizer N', &
       '                       soybean fixation', &
       '         supplementary P flux for plant', &
       '                    forest N fertiziler', &
       '                rate of denitrification', &
       '                            fire N loss', &
       '            N2O flux from nitrification', &
       ' soil mineral NO3 pool loss to leaching', &
       '   soil mineral NO3 pool loss to runoff', &
       '   soil mineral N pool loss to leaching', &
       '          1-year wood product harvested', &
       '         10-year wood product harvested', &
       '        100-year wood product harvested', &
       '     SOM N loss from vertical transport', &
       '                  SOM N loss by erosion' &
       /)

  ! P inputs
  integer, parameter :: f_primp_to_labilep     = 32
  integer, parameter :: f_supplement_to_sminp  = 33
  integer, parameter :: f_supplement_to_plantp = 34
  integer, parameter :: f_pfert_dose           = 35

  ! P outputs
  integer, parameter :: f_secondp_to_occlp     = 36
  integer, parameter :: f_sminp_leached        = 37
  integer, parameter :: f_col_fire_ploss       = 38
  integer, parameter :: f_solutionp            = 39
  integer, parameter :: f_labilep              = 40
  integer, parameter :: f_secondp              = 41
  integer, parameter :: f_col_prod1p_loss      = 42
  integer, parameter :: f_col_prod10p_loss     = 43
  integer, parameter :: f_col_prod100p_loss    = 44
  integer, parameter :: f_som_p_yield          = 45
  integer, parameter :: f_labilep_yield        = 46
  integer, parameter :: f_secondp_yield        = 47

  integer, parameter, public :: p_f_size = f_secondp_yield - f_som_n_yield

  character(len=44), parameter :: p_f_name(p_f_size) = &
       (/&
       '                   primary mineral to labile', &
       '                         supplemental supply', &
       '                     supplementary for plant', &
       '                           forest fertiziler', &
       '   desorption of secondary mineral to labile', &
       '          soil mineral pool loss to leaching', &
       '                                   fire loss', &
       '                               solution_???_', &
       '                                 labile_???_', &
       '                               scondary_???_', &
       '   decomposition loss from 1-yr crop product', &
       '  decomposition loss from 10-yr crop product', &
       ' decomposition loss from 100-yr crop product', &
       '                         erosion loss of SOM', &
       '                 erosion loss of soil labile', &
       '      erosion loss of soil secondary mineral' &
       /)

  !--- S for state ---

  ! C
  integer, parameter :: s_totc_beg              =  1
  integer, parameter :: s_totpftc_beg           =  2
  integer, parameter :: s_cwdc_beg              =  3
  integer, parameter :: s_totlitc_beg           =  4
  integer, parameter :: s_totsomc_beg           =  5
  integer, parameter :: s_totprodc_beg          =  6
  integer, parameter :: s_ctrunc_beg            =  7
  integer, parameter :: s_cropseedc_deficit_beg =  8

  integer, parameter :: s_totc_end              =  9
  integer, parameter :: s_totpftc_end           = 10
  integer, parameter :: s_cwdc_end              = 11
  integer, parameter :: s_totlitc_end           = 12
  integer, parameter :: s_totsomc_end           = 13
  integer, parameter :: s_totprodc_end          = 14
  integer, parameter :: s_ctrunc_end            = 15
  integer, parameter :: s_cropseedc_deficit_end = 16

  integer, parameter :: s_c_error               = 17

  integer, parameter, public :: c_s_size = s_c_error

  integer, parameter, public :: c_s_name_size = 8
  character(len=25),parameter :: c_s_name(c_s_name_size) = &
       (/&
       '                PFT', &
       '                CWD', &
       '       Total litter', &
       '          Total SOM', &
       ' Total wood product', &
       '    Truncation sink', &
       '  Crop seed deficit', &
       '     Grid-level Err'  &
       /)

  ! N
  integer, parameter :: s_totpftn_beg           = 18
  integer, parameter :: s_cwdn_beg              = 19
  integer, parameter :: s_totlitn_beg           = 20
  integer, parameter :: s_totsomn_beg           = 21
  integer, parameter :: s_sminn_beg             = 22
  integer, parameter :: s_totprodn_beg          = 23
  integer, parameter :: s_plant_n_buffer_beg    = 24
  integer, parameter :: s_ntrunc_beg            = 25
  integer, parameter :: s_cropseedn_deficit_beg = 26

  integer, parameter :: s_totpftn_end           = 27
  integer, parameter :: s_cwdn_end              = 28
  integer, parameter :: s_totlitn_end           = 29
  integer, parameter :: s_totsomn_end           = 30
  integer, parameter :: s_sminn_end             = 31
  integer, parameter :: s_totprodn_end          = 32
  integer, parameter :: s_plant_n_buffer_end    = 33
  integer, parameter :: s_ntrunc_end            = 34
  integer, parameter :: s_cropseedn_deficit_end = 35

  integer, parameter :: s_n_error               = 36

  integer, parameter, public :: n_s_size = s_n_error - s_c_error
  
  character(len=25),parameter :: n_s_name(n_s_size) = &
       (/&
       '              total_n_beg', &
       '              total_n_end', &
       'coarse woody_debris_n_beg', &
       'coarse woody_debris_n_end', &
       '       total_litter_n_beg', &
       '       total_litter_n_end', &
       '          total_som_n_beg', &
       '          total_som_n_end', &
       '       soil_mineral_n_beg', &
       '       soil_mineral_n_end', &
       ' total_wood_product_n_beg', &
       ' total_wood_product_n_end', &
       '    truncation_sink_n_beg', &
       '    truncation_sink_n_end', &
       '   abstract_storage_n_beg', &
       '   abstract_storage_n_end', &
       '  crop_seed_deficit_n_beg', &
       '  crop_seed_deficit_n_end', &
       '                  error c'  &
       /)

  ! P
  integer, parameter :: s_totpftp_beg           = 37
  integer, parameter :: s_cwdp_beg              = 38
  integer, parameter :: s_totlitp_beg           = 39
  integer, parameter :: s_totsomp_beg           = 40
  integer, parameter :: s_totprodp_beg          = 41
  integer, parameter :: s_ptrunc_beg            = 42
  integer, parameter :: s_solutionp_beg         = 43
  integer, parameter :: s_labilep_beg           = 44
  integer, parameter :: s_secondp_beg           = 45
  integer, parameter :: s_cropseedp_deficit_beg = 46

  integer, parameter :: s_totpftp_end           = 47
  integer, parameter :: s_cwd_endp              = 48
  integer, parameter :: s_totlitp_end           = 49
  integer, parameter :: s_totsomp_end           = 50
  integer, parameter :: s_totprodp_end          = 51
  integer, parameter :: s_ptrunc_end            = 52
  integer, parameter :: s_solutionp_end         = 53
  integer, parameter :: s_labilep_end           = 54
  integer, parameter :: s_secondp_end           = 55
  integer, parameter :: s_cropseedp_deficit_end = 56

  integer, parameter :: s_p_error               = 57

  integer, parameter, public :: p_s_size = s_p_error - s_n_error

  character(len=25),parameter :: s_name_p(p_s_size) = &
       (/&
       '              total_n_beg', &
       '              total_n_end', &
       'coarse woody_debris_n_beg', &
       'coarse woody_debris_n_end', &
       '       total_litter_n_beg', &
       '       total_litter_n_end', &
       '          total_som_n_beg', &
       '          total_som_n_end', &
       ' total_wood_product_n_beg', &
       ' total_wood_product_n_end', &
       '       solution_???_n_beg', &
       '       solution_???_n_end', &
       '         labile_???_n_beg', &
       '         labile_???_n_end', &
       '      secondary_???_n_beg', &
       '      secondary_???_n_end', &
       '    truncation_sink_n_beg', &
       '    truncation_sink_n_end', &
       '  crop_seed_deficit_n_beg', &
       '  crop_seed_deficit_n_end', &
       '                  error c'  &
       /)

  !--- P for period ---

  integer, parameter :: p_inst = 1
  integer, parameter :: p_day  = 2
  integer, parameter :: p_mon  = 3
  integer, parameter :: p_ann  = 4
  integer, parameter :: p_inf  = 5

  integer, parameter, public :: p_size = p_inf

  character(len=8),parameter :: pname(p_size) = &
       (/'    inst','   daily',' monthly','  annual','all_time' /)

  real(r8) :: c_budg_fluxL(c_f_size, p_size) ! local sum, valid on all pes
  real(r8) :: n_budg_fluxL(n_f_size, p_size) ! local sum, valid on all pes
  real(r8) :: p_budg_fluxL(p_f_size, p_size) ! local sum, valid on all pes

  real(r8) :: c_budg_fluxG(c_f_size, p_size) ! global sum, valid only on root pe
  real(r8) :: n_budg_fluxG(n_f_size, p_size) ! global sum, valid only on root pe
  real(r8) :: p_budg_fluxG(p_f_size, p_size) ! global sum, valid only on root pe

  real(r8) :: c_budg_fluxN(c_f_size, p_size) ! counter, valid only on root pe
  real(r8) :: n_budg_fluxN(n_f_size, p_size) ! counter, valid only on root pe
  real(r8) :: p_budg_fluxN(p_f_size, p_size) ! counter, valid only on root pe

  real(r8) :: c_budg_stateL(c_s_size, p_size) ! local sum, valid on all pes
  real(r8) :: n_budg_stateL(n_s_size, p_size) ! local sum, valid on all pes
  real(r8) :: p_budg_stateL(p_s_size, p_size) ! local sum, valid on all pes

  real(r8), public :: c_budg_stateG(c_s_size, p_size) ! global sum, valid only on root pe
  real(r8), public :: n_budg_stateG(n_s_size, p_size) ! global sum, valid only on root pe
  real(r8), public :: p_budg_stateG(p_s_size, p_size) ! global sum, valid only on root pe

  !----- formats -----
  character(*),parameter :: C_FA0  = "('    ',12x,(42x,a10,2x),' | ',(3x,a10,2x))"
  character(*),parameter :: C_FF = "('    ',a51,f15.8,' | ',f18.2)"
  character(*),parameter :: FF2= "('    ',a12,a15,' | ',f18.2)"
  character(*),parameter :: C_FS = "('    ',a12,7(f18.2),26x,' | ',(f18.2))"
  character(*),parameter :: C_FS0= "('    ',12x,8(a19),' | ',(a19))"
  character(*),parameter :: C_FS2= "('    ',a12,67x,f18.2,67x,' | ',f18.2)"
  character(*),parameter :: C_FS3= "('    ',a12,8(f18.2),8x,' | ',(f18.2))"
  character(*),parameter :: C_FS_2 = "('    ',a25,f15.2,5x,f15.2,5x,' | ',f18.2)"
  character(*),parameter :: C_SA0  = "('    ',33x,2(5x,a3,8x),' | ',(8x,a12,2x))"
  character(*),parameter :: C_FS2_2= "('    ',a12,17x,f18.2,18x,' | ',f18.2)"
  character(*),parameter :: C_SA0_2= "('    ',31x,2(5x,a3,9x),' |',(8x,a12,2x))"
  character(*),parameter :: C_FS3_3= "('    ',a12,53x,' | ',(f18.2))"

contains

  !-----------------------------------------------------------------------
  subroutine CNPBudget_Reset(mode)

    !
    implicit none
    !
    character(len=*), intent(in),optional :: mode
    !
    call Reset(mode, c_budg_fluxL, c_budg_fluxG, c_budg_fluxN, c_budg_stateL, c_budg_stateG)
    call Reset(mode, n_budg_fluxL, n_budg_fluxG, n_budg_fluxN, n_budg_stateL, n_budg_stateG)
    call Reset(mode, p_budg_fluxL, p_budg_fluxG, p_budg_fluxN, p_budg_stateL, p_budg_stateG)

  end subroutine CNPBudget_Reset


  !-----------------------------------------------------------------------
  subroutine Reset(mode, budg_fluxL, budg_fluxG, budg_fluxN, budg_stateL, budg_stateG)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep
    !
    implicit none
    !
    character(len=*), intent(in), optional :: mode
    real(r8) :: budg_fluxL(:,:), budg_fluxG(:,:), budg_fluxN(:,:), budg_stateL(:,:), budg_stateG(:,:)
    !
    integer :: year, mon, day, sec
    integer :: ip
    character(*),parameter :: subName = '(CNPBudget_Reset) '

    if (.not.present(mode)) then
       call get_curr_date(year, mon, day, sec)
       call get_prev_date(year, mon, day, sec)

       do ip = 1,p_size
          if (ip == p_inst) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
             budg_stateL(:,ip) = 0.0_r8
             budg_stateG(:,ip) = 0.0_r8
          endif
          if (ip==p_inf .and. get_nstep()==1) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
       enddo

    else

       if (trim(mode) == 'inst') then
          budg_fluxL  (:,p_inst)   = 0.0_r8
          budg_fluxG  (:,p_inst)   = 0.0_r8
          budg_fluxN  (:,p_inst)   = 0.0_r8
          budg_stateL (:,p_inst)   = 0.0_r8
          budg_stateG (:,p_inst)   = 0.0_r8
       elseif (trim(mode) == 'day') then
          budg_fluxL  (:,p_day)    = 0.0_r8
          budg_fluxG  (:,p_day)    = 0.0_r8
          budg_fluxN  (:,p_day)    = 0.0_r8
          budg_stateL (:,p_day)    = 0.0_r8
          budg_stateG (:,p_day)    = 0.0_r8
       elseif (trim(mode) == 'mon') then
          budg_fluxL  (:,p_mon)    = 0.0_r8
          budg_fluxG  (:,p_mon)    = 0.0_r8
          budg_fluxN  (:,p_mon)    = 0.0_r8
          budg_stateL (:,p_mon)    = 0.0_r8
          budg_stateG (:,p_mon)    = 0.0_r8
       elseif (trim(mode) == 'ann') then
          budg_fluxL  (:,p_ann)    = 0.0_r8
          budg_fluxG  (:,p_ann)    = 0.0_r8
          budg_fluxN  (:,p_ann)    = 0.0_r8
          budg_stateL (:,p_ann)    = 0.0_r8
          budg_stateG (:,p_ann)    = 0.0_r8
       elseif (trim(mode) == 'inf') then
          budg_fluxL  (:,p_inf)    = 0.0_r8
          budg_fluxG  (:,p_inf)    = 0.0_r8
          budg_fluxN  (:,p_inf)    = 0.0_r8
          budg_stateL (:,p_inf)    = 0.0_r8
          budg_stateG (:,p_inf)    = 0.0_r8
       elseif (trim(mode) == 'all') then
          budg_fluxL  (:,:)        = 0.0_r8
          budg_fluxG  (:,:)        = 0.0_r8
          budg_fluxN  (:,:)        = 0.0_r8
          budg_stateL (:,:)        = 0.0_r8
          budg_stateG (:,:)        = 0.0_r8
       else
          call shr_sys_abort(subname//' ERROR in mode '//trim(mode))
       endif
    endif

  end subroutine Reset

  !-----------------------------------------------------------------------
  subroutine CNPBudget_Restart(bounds, ncid, flag)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio, only : ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    character(len=*),parameter :: subname = 'WaterBudget_Restart'

    select case (trim(flag))
    case ('define')
       call Restart_Define(bounds, ncid, 'C')
       call Restart_Define(bounds, ncid, 'N')
       call Restart_Define(bounds, ncid, 'P')

    case ('write')
       call Restart_Write(bounds, ncid, flag, 'C', c_f_size, c_s_size, &
            c_budg_fluxL, c_budg_fluxG, c_budg_fluxN, c_budg_stateL, c_budg_stateG)
       call Restart_Write(bounds, ncid, flag, 'N', n_f_size, n_s_size, &
            n_budg_fluxL, n_budg_fluxG, n_budg_fluxN, n_budg_stateL, n_budg_stateG)
       call Restart_Write(bounds, ncid, flag, 'P', p_f_size, p_s_size, &
            p_budg_fluxL, p_budg_fluxG, p_budg_fluxN, p_budg_stateL, p_budg_stateG)

    case ('read')
       call Restart_Read(bounds, ncid, flag, 'C', c_f_size, c_s_size, &
            c_budg_fluxG, c_budg_fluxN, c_budg_stateL)
       call Restart_Read(bounds, ncid, flag, 'N', n_f_size, n_s_size, &
            n_budg_fluxG, n_budg_fluxN, n_budg_stateL)
       call Restart_Read(bounds, ncid, flag, 'P', p_f_size, p_s_size, &
            p_budg_fluxG, p_budg_fluxN, n_budg_stateL)

    case default
       write(iulog,*) trim(subname),' ERROR: unknown flag = ',flag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CNPBudget_Restart

!-----------------------------------------------------------------------
  subroutine CNPBudget_Accum()
    !
    implicit none

    call Accum(c_s_size,c_budg_fluxN, c_budg_fluxL, c_budg_stateL, c_budg_stateG)
    !call Accum(n_s_size,n_budg_fluxN, n_budg_fluxL, n_budg_stateL, n_budg_stateG)
    !call Accum(p_s_size,p_budg_fluxN, p_budg_fluxL, p_budg_stateL, p_budg_stateG)

  end subroutine CNPBudget_Accum
!-----------------------------------------------------------------------
  subroutine Accum(s_size, budg_fluxN, budg_fluxL, budg_stateL, budg_stateG)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep
    !
    implicit none
    !
    integer                :: s_size
    real(r8)               :: budg_fluxN(:,:), budg_fluxL(:,:), budg_stateL(:,:), budg_stateG(:,:)
    !
    integer                :: ip, is
    integer                :: year_prev, month_prev, day_prev, sec_prev
    integer                :: year_curr, month_curr, day_curr, sec_curr
    character(*),parameter :: subName = '(CarbonBudget_Accum)'
    logical                :: update_state_beg, update_state_end

    call get_prev_date(year_prev, month_prev, day_prev, sec_prev)
    call get_curr_date(year_curr, month_curr, day_curr, sec_curr)

    do ip = p_inst+1, p_size
       budg_fluxL(:,ip) = budg_fluxL(:,ip) + budg_fluxL(:,p_inst)
       update_state_beg = .false.
       update_state_end = .false.

       select case (ip)
       case (p_day)
          if (sec_prev == 0) update_state_beg = .true.
          if (sec_curr == 0) update_state_end = .true.
       case (p_mon)
          if (sec_prev == 0 .and. day_prev == 1) update_state_beg = .true.
          if (sec_curr == 0 .and. day_curr == 1) update_state_end = .true.
       case (p_ann)
          if (sec_prev == 0 .and. day_prev == 1 .and. month_prev == 1) update_state_beg = .true.
          if (sec_curr == 0 .and. day_curr == 1 .and. month_curr == 1) update_state_end = .true.
       case (p_inf)
          if (get_nstep() == 1) update_state_beg = .true.
          update_state_end = .true.
       end select

       if (update_state_beg) then
          do is = 1, s_size/2 - 1
             budg_stateL(is,ip) = budg_stateL(is, p_inst)
          end do
       endif

       if (update_state_end) then
          do is = s_size/2 + 1, s_size-1
             budg_stateL(is,ip) = budg_stateL(is, p_inst)
          end do
       endif
       is= s_size; budg_stateL(is,ip) = budg_stateL(is,ip) + budg_stateL(is, p_inst)
    end do
    budg_fluxN(:,:) = budg_fluxN(:,:) + 1._r8

  end subroutine Accum

  !-----------------------------------------------------------------------
  subroutine CNPBudget_Run(bounds, atm2lnd_vars, lnd2atm_vars, grc_cs, grc_cf)
    !
    ! !DESCRIPTION:
    !
    use domainMod, only : ldomain
    use elm_varcon, only : re
    !
    implicit none

    type(bounds_type)           , intent(in) :: bounds
    type(atm2lnd_type)          , intent(in) :: atm2lnd_vars
    type(lnd2atm_type)          , intent(in) :: lnd2atm_vars
    type(gridcell_carbon_state) , intent(in) :: grc_cs
    type(gridcell_carbon_flux)  , intent(in) :: grc_cf

    call CBudget_Run(bounds, atm2lnd_vars, lnd2atm_vars, grc_cs, grc_cf, c_budg_fluxL, c_budg_stateL)

  end subroutine CNPBudget_Run
    
  !-----------------------------------------------------------------------
  subroutine CBudget_Run(bounds, atm2lnd_vars, lnd2atm_vars, grc_cs, grc_cf, budg_fluxL, budg_stateL)
    !
    ! !DESCRIPTION:
    !
    use domainMod, only : ldomain
    use elm_varcon, only : re
    !
    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(lnd2atm_type)          , intent(in)    :: lnd2atm_vars
    type(gridcell_carbon_state) , intent(in)    :: grc_cs
    type(gridcell_carbon_flux)  , intent(in)    :: grc_cf
    real(r8)                    , intent(inout) :: budg_fluxL(:,:), budg_stateL(:,:)
    !
    ! !LOCAL VARIABLES:
    integer  :: g, nf, ns, ip
    real(r8) :: af, one_over_re2

    associate(                                                       &
         beg_totc                  => grc_cs%beg_totc              , & ! Input: [real(r8) (:)] (gC/m2) total column carbon, incl veg and cpool
         beg_totpftc               => grc_cs%beg_totpftc           , & ! Input: [real(r8) (:)] (gC/m2) patch-level carbon aggregated to column level, incl veg and cpool
         beg_cwdc                  => grc_cs%beg_cwdc              , & ! Input: [real(r8) (:)] (gC/m2) total column coarse woody debris carbon
         beg_totsomc               => grc_cs%beg_totsomc           , & ! Input: [real(r8) (:)] (gC/m2) total column soil organic matter carbon
         beg_totlitc               => grc_cs%beg_totlitc           , & ! Input: [real(r8) (:)] (gC/m2) total column litter carbon
         beg_totprodc              => grc_cs%beg_totprodc          , & ! Input: [real(r8) (:)] (gC/m2) total column wood product carbon
         beg_ctrunc                => grc_cs%beg_ctrunc            , & ! Input: [real(r8) (:)] (gC/m2) total column truncation carbon sink
         beg_cropseedc_deficit     => grc_cs%beg_cropseedc_deficit , & ! Input: [real(r8) (:)] (gC/m2) column carbon pool for seeding new growth
         end_totc                  => grc_cs%end_totc              , & ! Input: [real(r8) (:)] (gC/m2) total column carbon, incl veg and cpool
         end_totpftc               => grc_cs%end_totpftc           , & ! Input: [real(r8) (:)] (gC/m2) patch-level carbon aggregated to column level, incl veg and cpool
         end_cwdc                  => grc_cs%end_cwdc              , & ! Input: [real(r8) (:)] (gC/m2) total column coarse woody debris carbon
         end_totsomc               => grc_cs%end_totsomc           , & ! Input: [real(r8) (:)] (gC/m2) total column soil organic matter carbon
         end_totlitc               => grc_cs%end_totlitc           , & ! Input: [real(r8) (:)] (gC/m2) total column litter carbon
         end_totprodc              => grc_cs%end_totprodc          , & ! Input: [real(r8) (:)] (gC/m2) total column wood product carbon
         end_ctrunc                => grc_cs%end_ctrunc            , & ! Input: [real(r8) (:)] (gC/m2) total column truncation carbon sink
         end_cropseedc_deficit     => grc_cs%end_cropseedc_deficit , & ! Input: [real(r8) (:)] (gC/m2) column carbon pool for seeding new growth
         errcb                     => grc_cs%errcb                 , & ! Input: [real(r8) (:)] (gC/m^2) carbon mass balance error
         gpp                       => grc_cf%gpp                   , & ! Input: [real(r8) (:)] (gC/m2/s) gross primary production
         er                        => grc_cf%er                    , & ! Input: [real(r8) (:)] (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
         fire_closs                => grc_cf%fire_closs            , & ! Input: [real(r8) (:)] (gC/m2/s) total column-level fire C loss
         prod1c_loss               => grc_cf%prod1_loss            , & ! Input: [real(r8) (:)] (gC/m2/s) crop leafc harvested
         prod10c_loss              => grc_cf%prod10_loss           , & ! Input: [real(r8) (:)] (gC/m2/s) 10-year wood C harvested
         prod100c_loss             => grc_cf%prod100_loss          , & ! Input: [real(r8) (:)] (gC/m2/s) 100-year wood C harvested
         hrv_xsmrpool_to_atm       => grc_cf%hrv_xsmrpool_to_atm   , & ! Input: [real(r8) (:)] (gC/m2/s) excess MR pool harvest mortality
         som_c_leached             => grc_cf%som_c_leached         , & ! Input: [real(r8) (:)] (gC/m^2/s)total SOM C loss from vertical transport
         som_c_yield               => grc_cf%somc_yield            , & ! Input: [real(r8) (:)] (gC/m^2/s)total SOM C loss by erosion
         grc_dwt_conv_cflux        => grc_cf%dwt_conv_cflux        , & ! Input: [real(r8) (:)] (gC/m^2/s) flux to atmosphere during transient run
         grc_dwt_seedc_to_leaf     => grc_cf%dwt_seedc_to_leaf     , & ! Input: [real(r8) (:)] (gC/m^2/s) seed to leaf flux during transient run
         grc_dwt_seedc_to_deadstem => grc_cf%dwt_seedc_to_deadstem   & !Input: [real(r8) (:)] (gC/m^2/s) seed to dead stem flux during transient run
         )

      ip = p_inst
      budg_stateL(:,ip) = 0._r8
      one_over_re2 = 1._r8/(re**2._r8)

      do g = bounds%begg, bounds%endg
         af = (ldomain%area(g) * one_over_re2) * & ! area (converting km**2 to radians**2)
              ldomain%frac(g)                      ! land fraction

         ! fluxes
         nf = f_gpp                   ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + gpp(g)                       *af
         nf = f_er                    ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - er(g)                        *af
         nf = f_fire_closs            ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - fire_closs(g)                *af
         nf = f_hrv_xsmrpool_to_atm   ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - hrv_xsmrpool_to_atm(g)       *af
         nf = f_prod1c_loss           ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - prod1c_loss(g)               *af
         nf = f_prod10c_loss          ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - prod10c_loss(g)              *af
         nf = f_prod100c_loss         ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - prod100c_loss(g)             *af
         nf = f_som_c_leached         ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - som_c_leached(g)             *af
         nf = f_som_c_yield           ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - som_c_yield(g)               *af
         nf = f_dwt_conv_cflux        ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) - grc_dwt_conv_cflux(g)        *af
         nf = f_dwt_seedc_to_leaf     ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + grc_dwt_seedc_to_leaf(g)     *af
         nf = f_dwt_seedc_to_deadstem ; budg_fluxL(nf,ip) = budg_fluxL(nf,ip) + grc_dwt_seedc_to_deadstem(g) *af

         ! states
         ns = s_totc_beg              ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_totc(g)              *af
         ns = s_totpftc_beg           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_totpftc(g)           *af
         ns = s_cwdc_beg              ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_cwdc(g)              *af
         ns = s_totlitc_beg           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_totlitc(g)           *af
         ns = s_totsomc_beg           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_totsomc(g)           *af
         ns = s_totprodc_beg          ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_totprodc(g)          *af
         ns = s_ctrunc_beg            ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_ctrunc(g)            *af
         ns = s_cropseedc_deficit_beg ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + beg_cropseedc_deficit(g) *af

         ns = s_totc_end              ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_totc(g)              *af
         ns = s_totpftc_end           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_totpftc(g)           *af
         ns = s_cwdc_end              ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_cwdc(g)              *af
         ns = s_totlitc_end           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_totlitc(g)           *af
         ns = s_totsomc_end           ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_totsomc(g)           *af
         ns = s_totprodc_end          ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_totprodc(g)          *af
         ns = s_ctrunc_end            ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_ctrunc(g)            *af
         ns = s_cropseedc_deficit_end ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + end_cropseedc_deficit(g) *af

         ns = s_c_error               ; budg_stateL(ns,ip) = budg_stateL(ns,ip) + errcb(g)                 *af
      end do

    end associate

  end subroutine CBudget_Run

  !-----------------------------------------------------------------------
  subroutine CNPBudget_Sum0()
    !
    implicit none

    call Sum0(c_f_size, c_s_size, c_budg_fluxL, c_budg_fluxG, c_budg_stateL, c_budg_stateG)
    call Sum0(n_f_size, n_s_size, n_budg_fluxL, n_budg_fluxG, n_budg_stateL, n_budg_stateG)
    call Sum0(p_f_size, p_s_size, p_budg_fluxL, p_budg_fluxG, p_budg_stateL, p_budg_stateG)

  end subroutine CNPBudget_Sum0

  !-----------------------------------------------------------------------
  subroutine Sum0(f_size, s_size, budg_fluxL, budg_fluxG, budg_stateL, budg_stateG)
    !
    use spmdMod    , only : mpicom
    use shr_mpi_mod, only : shr_mpi_sum
    !
    implicit none
    !
    integer                :: f_size, s_size
    real(r8)               :: budg_fluxL(:,:), budg_fluxG(:,:)
    real(r8)               :: budg_stateL(:,:), budg_stateG(:,:)
    !
    real(r8)               :: budg_fluxGtmp(f_size,p_size) ! temporary sum
    real(r8)               :: budg_stateGtmp(s_size,p_size) ! temporary sum
    character(*),parameter :: subName = '(Sum0)'

    budg_fluxGtmp = 0._r8
    budg_stateGtmp = 0._r8

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName)
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName)


    budg_fluxG  = budg_fluxG + budg_fluxGtmp
    budg_stateG = budg_stateGtmp

    budg_fluxL(:,:)       = 0._r8 ! reset all fluxes
    budg_stateL(:,p_inst) = 0._r8 ! only reset instantaneous states

  end subroutine Sum0

  !-----------------------------------------------------------------------

  subroutine CNPBudget_Print(budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend)
    !
    implicit none
    !
    integer , intent(in) :: budg_print_inst
    integer , intent(in) :: budg_print_daily
    integer , intent(in) :: budg_print_month
    integer , intent(in) :: budg_print_ann
    integer , intent(in) :: budg_print_ltann
    integer , intent(in) :: budg_print_ltend

    call Print( carbon_budget, c_f_size, c_s_size            , &
         c_budg_fluxL, c_budg_fluxG, c_budg_fluxN            , &
         c_budg_stateL, c_budg_stateG                        , &
         budg_print_inst, budg_print_daily, budg_print_month , &
         budg_print_ann, budg_print_ltann, budg_print_ltend)

  end subroutine CNPBudget_Print

  !-----------------------------------------------------------------------

  subroutine Print(budget_type, f_size, s_size, budg_fluxL, budg_fluxG, budg_fluxN, &
       budg_stateL, budg_stateG, &
       budg_print_inst,  budg_print_daily,  budg_print_month,  &
       budg_print_ann,  budg_print_ltann,  budg_print_ltend)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
    use shr_const_mod   , only : shr_const_pi
    !
    implicit none
    !
    integer              :: budget_type, f_size, s_size
    real(r8)             :: budg_fluxL(:,:), budg_fluxG(:,:), budg_fluxN(:,:)
    real(r8)             :: budg_stateL(:,:), budg_stateG(:,:)
    integer , intent(in) :: budg_print_inst
    integer , intent(in) :: budg_print_daily
    integer , intent(in) :: budg_print_month
    integer , intent(in) :: budg_print_ann
    integer , intent(in) :: budg_print_ltann
    integer , intent(in) :: budg_print_ltend
    !
    ! !LOCAL VARIABLES:
    integer :: s,f,ic,nf,ip,is ! data array indicies
    integer :: plev        ! print level
    integer :: year, mon, day, sec
    integer :: cdate
    logical :: sumdone
    real(r8) :: unit_conversion
    real(r8) :: budg_fluxGpr (f_size,p_size) ! values to print, scaled and such

    sumdone = .false.

    if (get_nstep() <= 1) then
       call get_prev_date(year, mon, day, sec);
    else
       call get_curr_date(year, mon, day, sec);
    end if

    cdate = year*10000 + mon*100 + day

    do ip = 1,p_size
       plev = 0
       if (ip == p_inst) then
          plev = max(plev,budg_print_inst)
       endif
       if (ip==p_day .and. sec==0) then
          plev = max(plev,budg_print_daily)
       endif
       if (ip==p_mon .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_month)
       endif
       if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ann)
       endif
       if (ip==p_inf .and. mon==1 .and. day==1 .and. sec==0) then
          plev = max(plev,budg_print_ltann)
       endif

       if (plev > 0) then
          unit_conversion = 1._r8/(4.0_r8*shr_const_pi)*1.0e9_r8
          if (.not.sumdone) then
             sumdone = .true.
             call Sum0(f_size, s_size, budg_fluxL, budg_fluxG, budg_stateL, budg_stateG)
             budg_fluxGpr = budg_fluxG
             budg_fluxGpr = budg_fluxGpr*unit_conversion
             budg_fluxGpr = budg_fluxGpr/budg_fluxN
          end if

          if (ip == p_day .and. get_nstep() == 1) cycle
          if (ip == p_mon .and. get_nstep() == 1) cycle
          if (ip == p_ann .and. get_nstep() == 1) cycle
          if (ip == p_inf .and. get_nstep() == 1) cycle

          if (masterproc) then
             select case (budget_type)
             case (carbon_budget)
                call CarbonBudget_Message(ip, cdate, sec, f_size, s_size, &
                     budg_stateG, budg_fluxG, budg_fluxGpr, unit_conversion)

             case (nitrogen_budget)
             case (phosphorus_budget)
             case default
                call endrun(msg='ERROR: Unknown budget type'//&
                     errMsg(__FILE__, __LINE__))
             end select
          end if
       end if
    end do

  end subroutine Print

  !-----------------------------------------------------------------------
  subroutine CarbonBudget_Message(ip, cdate, sec, f_size, s_size, budg_stateG, budg_fluxG, budg_fluxGpr, unit_conversion)
    !
    use clm_time_manager, only : get_curr_date, get_prev_date, get_nstep, get_step_size
    use shr_const_mod   , only : shr_const_pi
    !
    implicit none
    !
    integer  , intent(in) :: ip, cdate, sec, f_size, s_size
    real(r8) , intent(in) :: budg_stateG(:,:), budg_fluxG(:,:), budg_fluxGpr (:,:)
    real(r8) , intent(in) :: unit_conversion
    !
    !
    ! !LOCAL VARIABLES:
    integer :: f, s, s_beg, s_end ! data array indicies
    real(r8) :: time_integrated_flux, state_net_change
    real(r8) :: relative_error
    real(r8), parameter :: error_tol = 0.01_r8
    real(r8), parameter :: relative_error_tol = 1.e-10_r8 ! [%]

    write(iulog,*   )''
    write(iulog,*   )'NET CARBON FLUXES : period ',trim(pname(ip)),': date = ',cdate,sec
    write(iulog,C_FA0 )'  Time  ','  Time    '
    write(iulog,C_FA0 )'averaged','integrated'
    write(iulog,C_FA0 )'kgC/m2/s*1e6','kgC/m2*1e6'

    write(iulog,'(71("-"),"|",20("-"))')

    do f = 1, f_size
       write(iulog,C_FF)c_f_name(f),budg_fluxGpr(f,ip),budg_fluxG(f,ip)*unit_conversion*get_step_size()
    end do

    write(iulog,'(71("-"),"|",20("-"))')

    write(iulog,C_FF)'   *SUM*', &
         sum(budg_fluxGpr(:,ip)), sum(budg_fluxG(:,ip))*unit_conversion*get_step_size()

    time_integrated_flux = sum(budg_fluxG(:,ip))*unit_conversion*get_step_size()

    write(iulog,'(71("-"),"|",20("-"))')

    write(iulog,*)''
    write(iulog,*)'CARBON STATES (kgC/m2*1e6): period ',trim(pname(ip)),': date = ',cdate,sec

    write(iulog,*)''
    write(iulog,C_SA0_2)'beg','end','*NET CHANGE*'
    do s = 1,c_s_name_size-1
       s_beg = s_totc_beg + s
       s_end = s_totc_end + s
       write(iulog,C_FS_2)c_s_name(s),&
            budg_stateG(s_beg,ip)*unit_conversion, &
            budg_stateG(s_end,ip)*unit_conversion, &
            (budg_stateG(s_end,ip) - budg_stateG(s_beg,ip))*unit_conversion
    end do
    write(iulog,C_FS_2)c_s_name(c_s_name_size),0._r8, budg_stateG(s_c_error,ip) *unit_conversion, &
       budg_stateG(s_c_error,ip) *unit_conversion


    write(iulog,'(70("-"),"|",23("-"))')

    write(iulog,C_FS2_2)'       *SUM*', &
         (budg_stateG(s_totpftc_end           ,ip) - budg_stateG(s_totpftc_beg          ,ip))*unit_conversion + &
         (budg_stateG(s_cwdc_end              ,ip) - budg_stateG(s_cwdc_beg             ,ip))*unit_conversion + &
         (budg_stateG(s_totlitc_end           ,ip) - budg_stateG(s_totlitc_beg          ,ip))*unit_conversion + &
         (budg_stateG(s_totsomc_end           ,ip) - budg_stateG(s_totsomc_beg          ,ip))*unit_conversion + &
         (budg_stateG(s_totprodc_end          ,ip) - budg_stateG(s_totprodc_beg         ,ip))*unit_conversion + &
         (budg_stateG(s_ctrunc_end            ,ip) - budg_stateG(s_ctrunc_beg           ,ip))*unit_conversion + &
         (budg_stateG(s_cropseedc_deficit_end ,ip) - budg_stateG(s_cropseedc_deficit_beg,ip))*unit_conversion + &
         budg_stateG(s_c_error,ip) *unit_conversion, &
         (budg_stateG(s_totc_end     ,ip) - budg_stateG(s_totc_beg     ,ip))*unit_conversion + &
         budg_stateG(s_c_error,ip) *unit_conversion

    state_net_change = (budg_stateG(s_totc_end, ip) - budg_stateG(s_totc_beg, ip))*unit_conversion + &
         budg_stateG(s_c_error,ip) *unit_conversion

    relative_error = abs(time_integrated_flux - state_net_change)/(budg_stateG(s_totc_end, ip)*unit_conversion) * 100._r8

    if (relative_error > relative_error_tol) then
       write(iulog,*)'time integrated flux = ',time_integrated_flux
       write(iulog,*)'net change in state  = ',state_net_change
       write(iulog,*)'current state        = ',budg_stateG(s_totc_end, ip)
       write(iulog,*)'relative error [%]   = ',relative_error
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    write(iulog,'(70("-"),"|",23("-"))')

  end subroutine CarbonBudget_Message

  !-----------------------------------------------------------------------
  subroutine Restart_Define(bounds, ncid, name)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*)                 :: name

    call ncd_defvar(varname=trim(name)//'_budg_fluxG', xtype=ncd_double, &
         dim1name=trim(name)//'_budg_flux', &
         long_name=trim(name)//'_budg_fluxG', units='mm', ncid=ncid)

    call ncd_defvar(varname=trim(name)//'_budg_fluxN', xtype=ncd_double, &
         dim1name=trim(name)//'_budg_flux', &
         long_name=trim(name)//'_budg_fluxN', units='-', ncid=ncid)

    call ncd_defvar(varname=trim(name)//'_budg_stateG', xtype=ncd_double, &
         dim1name=trim(name)//'_budg_state', &
         long_name=trim(name)//'_budg_stateG', units='-', ncid=ncid)

  end subroutine Restart_Define

  !-----------------------------------------------------------------------
  subroutine Restart_Write(bounds, ncid, flag, name, f_size, s_size, &
       budg_fluxL, budg_fluxG, budg_fluxN, budg_stateL, budg_stateG)
    !
    use ncdio_pio   , only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio   , only : ncd_defvar
    use spmdMod     , only : mpicom
    use shr_mpi_mod , only : shr_mpi_sum
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid                               ! netcdf id
    character(len=*) , intent(in)    :: flag                               ! 'read' or 'write'
    character(len=*), intent(in)     :: name                               ! name of bgc
    integer                          :: f_size, s_size                     ! size of fluxes and states
    real(r8)                         :: budg_fluxL(:,:), budg_fluxG(:,:)   ! fluxes local and global
    real(r8)                         :: budg_fluxN(:,:)                    ! count
    real(r8)                         :: budg_stateL(:,:), budg_stateG(:,:) ! states local and global
    !
    ! !LOCAL VARIABLES:
    real(r8) :: budg_fluxGtmp(f_size,p_size) ! temporary sum
    real(r8) :: budg_stateGtmp(s_size,p_size) ! temporary sum
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    integer  :: f, s, p, count
    character(*),parameter :: subName = '(Restart_Write) '

    call shr_mpi_sum(budg_fluxL, budg_fluxGtmp, mpicom, subName)
    call shr_mpi_sum(budg_stateL, budg_stateGtmp, mpicom, subName)

    ! Copy data from 2D into 1D array
    count = 0
    do f = 1, f_size
       do p = 1, p_size
          count = count + 1
          budg_fluxG_1D(count) = budg_fluxG(f,p) + budg_fluxGtmp(f,p)
          budg_fluxN_1D(count) = budg_fluxN(f,p)
       end do
    end do

    ! Copy data from 2D into 1D array
    count = 0
    do s = 1, s_size
       do p = 1, p_size
          count = count + 1
          budg_stateG_1D(count) = budg_stateGtmp(s,p)
       end do
    end do

    call ncd_io(flag=flag, varname=trim(name)//'_budg_fluxG', data=budg_fluxG_1D, ncid=ncid)

    call ncd_io(flag=flag, varname=trim(name)//'_budg_fluxN', data=budg_fluxN_1D, ncid=ncid)

    call ncd_io(flag=flag, varname=trim(name)//'_budg_stateG', data=budg_stateG_1D, ncid=ncid)

  end subroutine Restart_Write

  !-----------------------------------------------------------------------
  subroutine Restart_Read(bounds, ncid, flag, name, f_size, s_size, budg_fluxG, budg_fluxN, budg_stateL)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_int
    use ncdio_pio, only : ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    character(len=*) , intent(in)     :: name
    integer                          :: f_size, s_size
    real(r8)                         :: budg_fluxG(:,:)
    real(r8)                         :: budg_fluxN(:,:)
    real(r8)                         :: budg_stateL(:,:)
    character(len=20)                :: var_name
    !
    ! !LOCAL VARIABLES:
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    integer  :: f, s, p, count

    call ncd_io(flag=flag, varname=trim(name)//'_budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=trim(name)//'_budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=trim(name)//'_budg_stateG', data=budg_stateG_1D, ncid=ncid)

    ! Copy data from 1D into 2D array
    count = 0
    do f = 1, f_size
       do p = 1, p_size
          count = count + 1
          budg_fluxG(f,p) = budg_fluxG_1D(count)
          budg_fluxN(f,p) = budg_fluxN_1D(count)
       end do
    end do

    ! Copy data from 1D into 2D array
    if (masterproc) then
       count = 0
       do s = 1, s_size
          do p = 1, p_size
             count = count + 1
             budg_stateL(s,p) = budg_stateG_1D(count)
          end do
       end do
    end if

  end subroutine Restart_Read

  !-----------------------------------------------------------------------

  subroutine CNPBudget_SetBeginningMonthlyStates(bounds, col_cs, grc_cs)
    !
    use GridcellDataType, only : gridcell_carbon_state
    use ColumnDataType  , only : column_carbon_state
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    type(column_carbon_state)  , intent(in)    :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs

    call CBudget_SetBeginningMonthlyStates(bounds, col_cs, grc_cs)

  end subroutine CNPBudget_SetBeginningMonthlyStates


  !-----------------------------------------------------------------------

  subroutine CNPBudget_SetEndingMonthlyStates(bounds, col_cs, grc_cs)
    !
    use GridcellDataType, only : gridcell_carbon_state
    use ColumnDataType  , only : column_carbon_state
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)          , intent(in)    :: bounds
    type(column_carbon_state)  , intent(in)    :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs

    call CBudget_SetEndingMonthlyStates(bounds, col_cs, grc_cs)

  end subroutine CNPBudget_SetEndingMonthlyStates

  !-----------------------------------------------------------------------
  subroutine CBudget_SetBeginningMonthlyStates(bounds, col_cs, grc_cs)
    !
    ! !DESCRIPTION:
    ! Set grid-level carbon states at the beginning of a month
    !
    ! !USES:
    use subgridAveMod    , only : p2c, c2g
    use elm_varpar       , only : nlevgrnd, nlevsoi, nlevurb
    use elm_varcon       , only : spval
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon    , only : icol_road_perv, icol_road_imperv
    use clm_time_manager , only : get_curr_date, get_prev_date, get_nstep
    use GridcellDataType , only : gridcell_carbon_state
    use ColumnDataType   , only : column_carbon_state
    !
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    type(column_carbon_state)  , intent(in)    :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    !
    ! !LOCAL VARIABLES:
    integer :: year_prev, month_prev, day_prev, sec_prev
    integer :: year_curr, month_curr, day_curr, sec_curr
    !-----------------------------------------------------------------------

    associate(                                                       &
         begcb             =>    col_cs%begcb         , & ! Input : [real(r8) (:)   ]  carbon mass begining of the time step
         endcb             =>    col_cs%endcb         , & ! Input : [real(r8) (:)   ]  carbon mass begining of the time step
         tcs_month_beg_grc =>    grc_cs%tcs_month_beg   & ! Output: [real(r8) (:)   ]  grid-level carbon mass at the begining of a month
         )

      ! Get current and previous dates to determine if a new month started
      call get_prev_date(year_curr, month_curr, day_curr, sec_curr);
      call get_prev_date(year_prev, month_prev, day_prev, sec_prev)

      ! If at the beginning of a simulation, save grid-level TCS based on
      ! 'begcb' from the current time step
      if ( day_curr == 1 .and. sec_curr == 0 .and. get_nstep() <= 1 ) then
         call c2g( bounds, &
              begcb(bounds%begc:bounds%endc), &
              tcs_month_beg_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )
      endif

      ! If multiple steps into a simulation and the last time step was the
      ! end of a month, save grid-level TCS based on 'endcb' from the last
      ! time step
      if (get_nstep() > 1 .and. day_prev == 1 .and. sec_prev == 0) then
         call c2g( bounds, &
              endcb(bounds%begc:bounds%endc), &
              tcs_month_beg_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )
      endif

    end associate

  end subroutine CBudget_SetBeginningMonthlyStates

  !-----------------------------------------------------------------------
  subroutine CBudget_SetEndingMonthlyStates(bounds, col_cs, grc_cs)
    !
    ! !DESCRIPTION:
    ! Set grid-level carbon states at the beginning of a month
    !
    ! !USES:
    use subgridAveMod    , only : p2c, c2g
    use elm_varpar       , only : nlevgrnd, nlevsoi, nlevurb
    use elm_varcon       , only : spval
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon    , only : icol_road_perv, icol_road_imperv
    use clm_time_manager , only : get_curr_date, get_prev_date, get_nstep
    use GridcellDataType , only : gridcell_carbon_state
    use ColumnDataType   , only : column_carbon_state
    !
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds
    type(column_carbon_state)  , intent(in)    :: col_cs
    type(gridcell_carbon_state), intent(inout) :: grc_cs
    !
    ! !LOCAL VARIABLES:
    integer :: year_prev, month_prev, day_prev, sec_prev
    integer :: year_curr, month_curr, day_curr, sec_curr
    !-----------------------------------------------------------------------

    associate(                                                       &
         begcb             =>    col_cs%begcb         , & ! Input : [real(r8) (:)   ]  carbon mass begining of the time step
         endcb             =>    col_cs%endcb         , & ! Input : [real(r8) (:)   ]  carbon mass begining of the time step
         tcs_month_end_grc =>    grc_cs%tcs_month_beg   & ! Output: [real(r8) (:)   ]  grid-level carbon mass at the begining of a month
         )

      ! Get current and previous dates to determine if a new month started
      call get_prev_date(year_curr, month_curr, day_curr, sec_curr);
      call get_prev_date(year_prev, month_prev, day_prev, sec_prev)

      ! If at the beginning of a simulation, save grid-level TCS based on
      ! 'begcb' from the current time step
      if ( day_curr == 1 .and. sec_curr == 0 .and. get_nstep() <= 1 ) then
         call c2g( bounds, &
              begcb(bounds%begc:bounds%endc), &
              tcs_month_end_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )
      endif

      ! If multiple steps into a simulation and the last time step was the
      ! end of a month, save grid-level TCS based on 'endcb' from the last
      ! time step
      if (get_nstep() > 1 .and. day_prev == 1 .and. sec_prev == 0) then
         call c2g( bounds, &
              endcb(bounds%begc:bounds%endc), &
              tcs_month_end_grc(bounds%begg:bounds%endg), &
              c2l_scale_type= 'unity', l2g_scale_type='unity' )
      endif

    end associate

  end subroutine CBudget_SetEndingMonthlyStates

end module CNPBudgetMod
