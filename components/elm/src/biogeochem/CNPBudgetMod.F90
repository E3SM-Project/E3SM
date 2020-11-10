module CNPBudgetMod
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_abort
  use decompMod         , only : bounds_type
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use atm2lndType       , only : atm2lnd_type
  use lnd2atmType       , only : lnd2atm_type
  use spmdMod           , only : masterproc
  use GridcellDataType  , only : grc_ws
  use ColumnDataType    , only : col_ws

  implicit none
  save
  private

  public :: CNPBudget_Reset

  !--- F for flux ---

  ! C inputs
  integer, parameter :: f_gpp                  = 1

  ! C outputs
  integer, parameter :: f_ar                   =  2
  integer, parameter :: f_hr                   =  3
  integer, parameter :: f_fire_closs           =  4
  integer, parameter :: f_hrv_xsmrpool_to_atm  =  5
  integer, parameter :: f_prod1c_loss          =  6
  integer, parameter :: f_prod10c_loss         =  7
  integer, parameter :: f_prod100c_loss        =  8
  integer, parameter :: f_som_c_leached        =  9
  integer, parameter :: f_som_c_yield          = 10

  integer, parameter, public :: c_f_size = f_som_c_yield

  character(len=51), parameter :: c_f_name(c_f_size) = &
       (/&
       '                              gross primary product', &
       '                            autotrophic respiration', &
       '                          heterotrophic respiration', &
       '                                        fire C loss', &
       '                   excess MR pool harvest mortality', &
       '        decomposition loss from 1-year product pool', &
       '       decomposition loss from 10-year product pool', &
       '      decomposition loss from 100-year product pool', &
       '                 SOM C loss from vertical transport', &
       '                                         SOM C loss'  &
       /)
       
       
  ! N inputs
  integer, parameter :: f_ndep_to_sminn        = 11
  integer, parameter :: f_nfix_to_ecosysn      = 12
  integer, parameter :: f_nfix_to_sminn        = 13
  integer, parameter :: f_supplement_to_sminn  = 14
  integer, parameter :: f_fert_to_sminn        = 15
  integer, parameter :: f_soyfixn_to_sminn     = 16
  integer, parameter :: f_supplement_to_plantn = 17
  integer, parameter :: f_nfert_dose           = 18

  ! N outputs
  integer, parameter :: f_denit                = 19
  integer, parameter :: f_fire_ploss           = 20
  integer, parameter :: f_n2o_nit            = 21
  integer, parameter :: f_smin_no3_leached     = 22
  integer, parameter :: f_smin_no3_runoff      = 23
  integer, parameter :: f_sminn_leached        = 24
  integer, parameter :: f_col_prod1n_loss      = 25
  integer, parameter :: f_col_prod10n_loss     = 26
  integer, parameter :: f_col_prod100n_loss    = 27
  integer, parameter :: f_som_n_leached        = 28
  integer, parameter :: f_som_n_yield          = 29

  integer, parameter, public :: n_f_size = f_som_n_yield     - f_som_c_yield

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
  integer, parameter :: f_primp_to_labilep     = 30
  integer, parameter :: f_supplement_to_sminp  = 31
  integer, parameter :: f_supplement_to_plantp = 32
  integer, parameter :: f_pfert_dose           = 33

  ! P outputs
  integer, parameter :: f_secondp_to_occlp     = 34
  integer, parameter :: f_sminp_leached        = 35
  integer, parameter :: f_col_fire_ploss       = 36
  integer, parameter :: f_solutionp            = 37
  integer, parameter :: f_labilep              = 38
  integer, parameter :: f_secondp              = 39
  integer, parameter :: f_col_prod1p_loss      = 40
  integer, parameter :: f_col_prod10p_loss     = 41
  integer, parameter :: f_col_prod100p_loss    = 42
  integer, parameter :: f_som_p_yield          = 43
  integer, parameter :: f_labilep_yield        = 44
  integer, parameter :: f_secondp_yield        = 45

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
       '                                solution ???', &
       '                                  labile ???', &
       '                                scondary ???', &
       '   decomposition loss from 1-yr crop product', &
       '  decomposition loss from 10-yr crop product', &
       ' decomposition loss from 100-yr crop product', &
       '                         erosion loss of SOM', &
       '                 erosion loss of soil labile', &
       '      erosion loss of soil secondary mineral' &
       /)

  !--- S for state ---

  ! C
  integer, parameter :: s_totpftc_beg           =  1
  integer, parameter :: s_totpftc_end           =  2
  integer, parameter :: s_cwdc_beg              =  3
  integer, parameter :: s_cwdc_end              =  4
  integer, parameter :: s_totlitc_beg           =  5
  integer, parameter :: s_totlitc_end           =  6
  integer, parameter :: s_totsomc_beg           =  7
  integer, parameter :: s_totsomc_end           =  8
  integer, parameter :: s_totprodc_beg          =  9
  integer, parameter :: s_totprodc_end          = 10
  integer, parameter :: s_ctrunc_beg            = 11
  integer, parameter :: s_ctrunc_end            = 12
  integer, parameter :: s_cropseedc_deficit_beg = 13
  integer, parameter :: s_cropseedc_deficit_end = 14
  integer, parameter :: s_c_error               = 15

  integer, parameter, public :: c_s_size = s_c_error

  character(len=25),parameter :: c_s_name(c_s_size) = &
       (/&
       '              total_c_beg', &
       '              total_c_end', &
       'coarse woody debris_c_beg', &
       'coarse woody debris_c_end', &
       '       total litter_c_beg', &
       '       total litter_c_end', &
       '          total_som_c_beg', &
       '          total_som_c_end', &
       ' total_wood_product_c_beg', &
       ' total_wood_product_c_end', &
       '    truncation_sink_c_beg', &
       '    truncation_sink_c_end', &
       '  crop_seed_deficit_c_beg', &
       '  crop_seed_deficit_c_end', &
       '                  error c'  &
       /)

  ! N
  integer, parameter :: s_totpftn_beg           = 16
  integer, parameter :: s_totpftn_end           = 17
  integer, parameter :: s_cwdn_beg              = 18
  integer, parameter :: s_cwdn_end              = 19
  integer, parameter :: s_totlitn_beg           = 20
  integer, parameter :: s_totlitn_end           = 21
  integer, parameter :: s_totsomn_beg           = 22
  integer, parameter :: s_totsomn_end           = 23
  integer, parameter :: s_sminn_beg             = 24
  integer, parameter :: s_sminn_end             = 25
  integer, parameter :: s_totprodn_beg          = 26
  integer, parameter :: s_totprodn_end          = 27
  integer, parameter :: s_ntrunc_beg            = 28
  integer, parameter :: s_ntrunc_end            = 29
  integer, parameter :: s_plant_n_buffer_beg    = 30
  integer, parameter :: s_plant_n_buffer_end    = 31
  integer, parameter :: s_cropseedn_deficit_beg = 32
  integer, parameter :: s_cropseedn_deficit_end = 33
  integer, parameter :: s_n_error               = 34

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
  integer, parameter :: s_totpftp_beg           = 35
  integer, parameter :: s_totpftp_end           = 36
  integer, parameter :: s_cwdp_beg              = 37
  integer, parameter :: s_cwd_endp              = 38
  integer, parameter :: s_totlitp_beg           = 39
  integer, parameter :: s_totlitp_end           = 40
  integer, parameter :: s_totsomp_beg           = 41
  integer, parameter :: s_totsomp_end           = 42
  integer, parameter :: s_totprodp_beg          = 43
  integer, parameter :: s_totprodp_end          = 44
  integer, parameter :: s_ptrunc_beg            = 45
  integer, parameter :: s_ptrunc_end            = 46
  integer, parameter :: s_solutionp_beg         = 47
  integer, parameter :: s_solutionp_end         = 48
  integer, parameter :: s_labilep_beg           = 49
  integer, parameter :: s_labilep_end           = 50
  integer, parameter :: s_secondp_beg           = 51
  integer, parameter :: s_secondp_end           = 52
  integer, parameter :: s_cropseedp_deficit_beg = 53
  integer, parameter :: s_cropseedp_deficit_end = 54
  integer, parameter :: s_p_error               = 55

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
  real(r8) :: p_budg_stateL(n_s_size, p_size) ! local sum, valid on all pes

  real(r8), public :: c_budg_stateG(c_s_size, p_size) ! global sum, valid only on root pe
  real(r8), public :: n_budg_stateG(n_s_size, p_size) ! global sum, valid only on root pe
  real(r8), public :: p_budg_stateG(p_s_size, p_size) ! global sum, valid only on root pe

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
    use clm_time_manager, only : get_curr_date, get_prev_date
    !
    implicit none
    !
    character(len=*), intent(in), optional :: mode
    real(r8) :: budg_fluxL(:,:), budg_fluxG(:,:), budg_fluxN(:,:), budg_stateL(:,:), budg_stateG(:,:)
    !
    integer :: year, mon, day, sec
    integer :: ip
    character(*),parameter :: subName = '(WaterBudget_Reset) '

    if (.not.present(mode)) then
       call get_curr_date(year, mon, day, sec)
       call get_prev_date(year, mon, day, sec)

       do ip = 1,p_size
          if (ip == p_inst) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_day .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_mon .and. day==1 .and. sec==0) then
             budg_fluxL(:,ip)  = 0.0_r8
             budg_fluxG(:,ip)  = 0.0_r8
             budg_fluxN(:,ip)  = 0.0_r8
          endif
          if (ip==p_ann .and. mon==1 .and. day==1 .and. sec==0) then
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
  subroutine CNP_Restart(bounds, ncid, flag)
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
       call Restart_Write(bounds, ncid, flag, 'P', n_f_size, n_s_size, &
            n_budg_fluxL, n_budg_fluxG, n_budg_fluxN, n_budg_stateL, n_budg_stateG)

    case ('read')
       call Restart_Read(bounds, ncid, flag, 'C', c_f_size, c_s_size, &
            c_budg_fluxG, c_budg_fluxN, c_budg_stateL)
       call Restart_Read(bounds, ncid, flag, 'N', n_f_size, n_s_size, &
            n_budg_fluxG, n_budg_fluxN, n_budg_stateL)
       call Restart_Read(bounds, ncid, flag, 'P', n_f_size, n_s_size, &
            n_budg_fluxG, n_budg_fluxN, n_budg_stateL)

    case default
       write(iulog,*) trim(subname),' ERROR: unknown flag = ',flag
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CNP_Restart

  !-----------------------------------------------------------------------
  subroutine Restart_Define(bounds, ncid, name)
    !
    use ncdio_pio, only : file_desc_t, ncd_io, ncd_double, ncd_defvar
    !
    implicit none
    !
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=1)                 :: name

    call ncd_defvar(varname=name//'_budg_fluxG', xtype=ncd_double, &
         dim1name=name//'_budg_flux', &
         long_name=name//'_budg_fluxG', units='mm', ncid=ncid)

    call ncd_defvar(varname=name//'_budg_fluxN', xtype=ncd_double, &
         dim1name=name//'_budg_flux', &
         long_name=name//'_budg_fluxN', units='-', ncid=ncid)

    call ncd_defvar(varname=name//'_budg_stateG', xtype=ncd_double, &
         dim1name=name//'_budg_state', &
         long_name=name//'_budg_stateG', units='mm', ncid=ncid)

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
    type(file_desc_t), intent(inout) :: ncid                                               ! netcdf id
    character(len=*) , intent(in)    :: flag                                               ! 'read' or 'write'
    character(len=1), intent(in)     :: name
    integer                          :: f_size, s_size
    real(r8)                         :: budg_fluxL(:,:), budg_fluxG(:,:)
    real(r8)                         :: budg_fluxN(:,:)
    real(r8)                         :: budg_stateL(:,:), budg_stateG(:,:)
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

    call ncd_io(flag=flag, varname=name//'_budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=name//'_budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=name//'_budg_stateG', data=budg_stateG_1D, ncid=ncid)

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
    character(len=1), intent(in)     :: name
    integer                          :: f_size, s_size
    real(r8)                         :: budg_fluxG(:,:)
    real(r8)                         :: budg_fluxN(:,:)
    real(r8)                         :: budg_stateL(:,:)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: budg_fluxG_1D (f_size*p_size)
    real(r8) :: budg_fluxN_1D (f_size*p_size)
    real(r8) :: budg_stateG_1D(s_size*p_size)
    integer  :: f, s, p, count

    call ncd_io(flag=flag, varname=name//'_budg_fluxG', data=budg_fluxG_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=name//'_budg_fluxN', data=budg_fluxN_1D, ncid=ncid)
    call ncd_io(flag=flag, varname=name//'_budg_stateG', data=budg_stateG_1D, ncid=ncid)

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

end module CNPBudgetMod
