module ace_comp_mod

  use esmf
  use eatmMod
  use eatmIO
  use mct_mod
  use seq_timemgr_mod, only: seq_timemgr_EClockGetData
  use shr_const_mod
  use shr_kind_mod, only: R4=>SHR_KIND_R4, R8=>SHR_KIND_R8, CL=>SHR_KIND_CL, IN=>SHR_KIND_IN
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort

  use ftorch, only: &
    torch_kCPU, &
    torch_model, &
    torch_tensor, &
    torch_delete, &
    torch_kFloat32, &
    torch_model_load, &
    torch_tensor_print, &
    torch_model_forward, &
    torch_tensor_from_array, &
    torch_tensor_from_blob

  use, intrinsic :: iso_c_binding, only: c_loc, c_int64_t, c_int

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: ace_comp_init
  public :: ace_comp_run
  public :: ace_comp_finalize

  !--------------------------------------------------------------------------
  ! Public module data
  !--------------------------------------------------------------------------
  integer, public :: n_input_channels=39  ! number of input channels to emulator
  integer, public :: n_output_channels=44 ! number of output channels to emulator
  integer, public :: eatm_idt=6 * 60 * 60 ! eatm timestep (6hr) in seconds

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  real(R8), parameter :: rdair  = SHR_CONST_RDAIR  ! dry air gas constant   ~ J/K/kg
  real(R8), parameter :: tKFrz  = SHR_CONST_TKFRZ

  ! Set up Torch data structures
  type(torch_model) :: ace_model
  type(torch_tensor), dimension(1) :: input_tensor
  type(torch_tensor), dimension(1) :: output_tensor

  integer(c_int64_t) :: input_tensor_shape(4)
  integer(c_int64_t) :: output_tensor_shape(4)

  integer(c_int) :: tensor_layout(4)

  ! TODO (AN): Parse from namelist
  character(len=*), parameter :: torchscript_file="/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/ace2_EAMv3_ckpt_traced.tar"
  character(len=*), parameter :: norm_file="/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/ace2_EAMv3_normalize.nc"
  character(len=*), parameter :: denorm_file="/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/ace2_EAMv3_denormalize.nc"
  save

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

  subroutine ace_comp_init(ggrid)

    implicit none

    type(mct_gGrid), intent(in), pointer :: ggrid
    integer     :: i, j, k           ! loop indicies

    input_tensor_shape = [ &
      int(1, kind=c_int), &
      int(n_input_channels, kind=c_int), &
      int(lsize_y, kind=c_int), &
      int(lsize_x, kind=c_int) &
    ]

    output_tensor_shape = [ &
      int(1, kind=c_int), &
      int(n_output_channels, kind=c_int), &
      int(lsize_y, kind=c_int), &
      int(lsize_x, kind=c_int) &
    ]

    tensor_layout = [1_c_int, 2_c_int, 4_c_int, 3_c_int]

    ! load the initial condition data into input array
    call ace_read_ic()

    ! load the traced model
    call torch_model_load(ace_model, torchscript_file, torch_kCPU)

    ! load the restart data (i.e. prognostic state) associated with IC
    call ace_read_restart()

    ! fill both time levels of intrp struct with restart data
    do k = 1, n_output_channels
      do j = 1, lsize_y
        do i = 1, lsize_x
          eatm_intrp%t_im1(k, i, j) = net_outputs(1, k, i, j)
          eatm_intrp%t_ip1(k, i, j) = net_outputs(1, k, i, j)
        end do
      end do
    end do

    ! using restart data from ACE set the fields passed to the coupler
    call ace_eatm_export(ggrid)

    call init_normalizer(normalizer, norm_file, n_input_channels)
    call init_normalizer(denormalizer, denorm_file, n_output_channels)

  end subroutine ace_comp_init

  subroutine ace_comp_run(EClock, ggrid)
    ! !DESCRIPTION: run method for ace model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock), intent(in) :: EClock
    type(mct_gGrid), intent(in), pointer :: ggrid

    !--- local ---
    integer     :: i, j, k           ! loop indicies
    real(R8)    :: t_frac            ! frac of cpl_t / eatm_t
    integer     :: t_modulo          ! frac of cpl_t / eatm_t
    real(R8)    :: cpl_dt            ! timestep
    integer(in) :: cpl_idt           ! integer timestep
    integer(in) :: stepno            ! step number
    integer(in) :: CurrentYMD        ! model date
    integer(in) :: CurrentTOD        ! model sec into model date
    logical     :: call_inference

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
    call seq_timemgr_EClockGetData( EClock, stepno=stepno, dtime=cpl_idt)

    cpl_dt = real(cpl_idt, kind=r8)

    write(logunit_atm, *) "stepno: ", stepno
    write(logunit_atm, *) "cpl_dt: ", cpl_dt
    write(logunit_atm, *) "cpl_idt: ", cpl_idt
    write(logunit_atm, *) "eatm_idt: ", eatm_idt
    write(logunit_atm, *) "CurrentYMD: ", CurrentYMD
    write(logunit_atm, *) "CurrentTOD: ", CurrentTOD
    call shr_sys_flush(logunit_atm)

    ! integer remained
    t_modulo = mod(CurrentTOD, eatm_idt)

    if (t_modulo .eq. 0) then

      ! populate net_imports array with IC/restart data passed to coupler for initializtion
      call ace_eatm_import()

      ! normalize, can probably happen after tensor is made becuase it's a pointer
      call normalizer%normalize(net_inputs)

      ! create input/output tensors based off net input/output arrays
      call torch_tensor_from_blob(&
        input_tensor(1), &
        c_loc(net_inputs), &
        4_c_int, &
        input_tensor_shape, &
        tensor_layout, &
        torch_kFloat32, &
        torch_kCPU &
      )
      call torch_tensor_from_blob(&
        output_tensor(1), &
        c_loc(net_outputs), &
        4_c_int, &
        output_tensor_shape, &
        tensor_layout, &
        torch_kFloat32, &
        torch_kCPU &
      )

      ! run inference
      call torch_model_forward(ace_model, input_tensor, output_tensor)

      ! denormalize
      call denormalizer%denormalize(net_outputs)

      ! advance the time levels
      do k = 1, n_output_channels
        do j = 1, lsize_y
          do i = 1, lsize_x
            eatm_intrp%t_im1(k, i, j) = eatm_intrp%t_ip1(k, i, j)
            eatm_intrp%t_ip1(k, i, j) = net_outputs(1, k, i, j)
          end do
        end do
      end do

    end if

    t_frac = real(t_modulo, kind=r8)/real(eatm_idt, kind=r8)

    ! time interpolate the results
    do k = 1, n_output_channels
      do j = 1, lsize_y
        do i = 1, lsize_x
          net_outputs(1, k, i, j) = eatm_intrp%t_im1(k, i, j) + t_frac * (eatm_intrp%t_ip1(k, i, j) - eatm_intrp%t_im1(k, i, j))
        end do
      end do
    end do

    call ace_eatm_export(ggrid)
  end subroutine ace_comp_run

  subroutine ace_comp_finalize()
    call torch_delete(ace_model)
    call torch_delete(input_tensor)
    call torch_delete(output_tensor)

    call finalize_normalizer(normalizer)
    call finalize_normalizer(denormalizer)
  end subroutine ace_comp_finalize

  subroutine ace_read_ic()
    ! AN TODO: add this to the namelist
    character(CL)     :: ic_filename = "/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/initial_conditions/1971010100_time_1.nc"
    type(file_desc_t) :: ncid    ! netcdf file id

    call ncd_pio_openfile(ncid, trim(ic_filename), 0)

    ! populate the input tensor with data from initial condition file
    call ace_read_netcdf(varname='LANDFRAC', data=net_inputs(1, 1, :, :), ncid=ncid)
    call ace_read_netcdf(varname='OCNFRAC', data=net_inputs(1, 2, :, :), ncid=ncid)
    call ace_read_netcdf(varname='ICEFRAC', data=net_inputs(1, 3, :, :), ncid=ncid)
    call ace_read_netcdf(varname='PHIS', data=net_inputs(1, 4, :, :), ncid=ncid)
    call ace_read_netcdf(varname='SOLIN', data=net_inputs(1, 5, :, :), ncid=ncid)
    call ace_read_netcdf(varname='PS', data=net_inputs(1, 6, :, :), ncid=ncid)
    call ace_read_netcdf(varname='TS', data=net_inputs(1, 7, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_0', data=net_inputs(1, 8, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_1', data=net_inputs(1, 9, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_2', data=net_inputs(1, 10, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_3', data=net_inputs(1, 11, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_4', data=net_inputs(1, 12, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_5', data=net_inputs(1, 13, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_6', data=net_inputs(1, 14, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_7', data=net_inputs(1, 15, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_0', data=net_inputs(1, 16, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_1', data=net_inputs(1, 17, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_2', data=net_inputs(1, 18, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_3', data=net_inputs(1, 19, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_4', data=net_inputs(1, 20, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_5', data=net_inputs(1, 21, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_6', data=net_inputs(1, 22, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_7', data=net_inputs(1, 23, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_0', data=net_inputs(1, 24, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_1', data=net_inputs(1, 25, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_2', data=net_inputs(1, 26, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_3', data=net_inputs(1, 27, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_4', data=net_inputs(1, 28, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_5', data=net_inputs(1, 29, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_6', data=net_inputs(1, 30, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_7', data=net_inputs(1, 31, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_0', data=net_inputs(1, 32, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_1', data=net_inputs(1, 33, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_2', data=net_inputs(1, 34, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_3', data=net_inputs(1, 35, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_4', data=net_inputs(1, 36, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_5', data=net_inputs(1, 37, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_6', data=net_inputs(1, 38, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_7', data=net_inputs(1, 39, :, :), ncid=ncid)

    call ncd_pio_closefile(ncid)
  end subroutine ace_read_ic

  subroutine ace_read_restart()
    ! AN TODO: add this to the namelist
    character(CL)     :: restart_filename = "/global/cfs/cdirs/e3sm/anolan/ACE2-E3SMv3/inference_output/ace_eatm_restart.nc"
    type(file_desc_t) :: ncid    ! netcdf file id

    call ncd_pio_openfile(ncid, trim(restart_filename), 0)

    ! probably only need to read the subset used in the export to cpl
    call ace_read_netcdf(varname='PS', data=net_outputs(1, 1, :, :), ncid=ncid)
    call ace_read_netcdf(varname='TS', data=net_outputs(1, 2, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_0', data=net_outputs(1, 3, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_1', data=net_outputs(1, 4, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_2', data=net_outputs(1, 5, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_3', data=net_outputs(1, 6, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_4', data=net_outputs(1, 7, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_5', data=net_outputs(1, 8, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_6', data=net_outputs(1, 9, :, :), ncid=ncid)
    call ace_read_netcdf(varname='T_7', data=net_outputs(1, 10, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_0', data=net_outputs(1, 11, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_1', data=net_outputs(1, 12, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_2', data=net_outputs(1, 13, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_3', data=net_outputs(1, 14, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_4', data=net_outputs(1, 15, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_5', data=net_outputs(1, 16, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_6', data=net_outputs(1, 17, :, :), ncid=ncid)
    call ace_read_netcdf(varname='specific_total_water_7', data=net_outputs(1, 18, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_0', data=net_outputs(1, 19, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_1', data=net_outputs(1, 20, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_2', data=net_outputs(1, 21, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_3', data=net_outputs(1, 22, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_4', data=net_outputs(1, 23, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_5', data=net_outputs(1, 24, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_6', data=net_outputs(1, 25, :, :), ncid=ncid)
    call ace_read_netcdf(varname='U_7', data=net_outputs(1, 26, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_0', data=net_outputs(1, 27, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_1', data=net_outputs(1, 28, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_2', data=net_outputs(1, 29, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_3', data=net_outputs(1, 30, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_4', data=net_outputs(1, 31, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_5', data=net_outputs(1, 32, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_6', data=net_outputs(1, 33, :, :), ncid=ncid)
    call ace_read_netcdf(varname='V_7', data=net_outputs(1, 34, :, :), ncid=ncid)
    call ace_read_netcdf(varname='LHFLX', data=net_outputs(1, 35, :, :), ncid=ncid)
    call ace_read_netcdf(varname='SHFLX', data=net_outputs(1, 36, :, :), ncid=ncid)
    call ace_read_netcdf(varname='surface_precipitation_rate', data=net_outputs(1, 37, :, :), ncid=ncid)
    call ace_read_netcdf(varname='surface_upward_longwave_flux', data=net_outputs(1, 38, :, :), ncid=ncid)
    call ace_read_netcdf(varname='FLUT', data=net_outputs(1, 39, :, :), ncid=ncid)
    call ace_read_netcdf(varname='FLDS', data=net_outputs(1, 40, :, :), ncid=ncid)
    call ace_read_netcdf(varname='FSDS', data=net_outputs(1, 41, :, :), ncid=ncid)
    call ace_read_netcdf(varname='surface_upward_shortwave_flux', data=net_outputs(1, 42, :, :), ncid=ncid)
    call ace_read_netcdf(varname='top_of_atmos_upward_shortwave_flux', data=net_outputs(1, 43, :, :), ncid=ncid)
    call ace_read_netcdf(varname='tendency_of_total_water_path_due_to_advection', data=net_outputs(1, 44, :, :), ncid=ncid)

    call ncd_pio_closefile(ncid)
  end subroutine ace_read_restart

  subroutine ace_eatm_import()
    ! !LOCAL VARIABLES:
    integer :: i, j

    do j = 1, lsize_y
      do i = 1, lsize_x
        ! net_inputs(1,  1, i, j) = lndfrac(i, j)            ! ACE2-EAMv3: LANDFRAC
        ! net_inputs(1,  2, i, j) = ocnfrac(i, j)            ! ACE2-EAMv3: OCNFRAC
        ! net_inputs(1,  3, i, j) = icefrac(i, j)            ! ACE2-EAMv3: ICEFRAC
        net_inputs(1,  1, i, j) = net_inputs(1,  1, i, j)  ! ACE2-EAMv3: LANDFRAC
        net_inputs(1,  2, i, j) = net_inputs(1,  2, i, j)  ! ACE2-EAMv3: OCNFRAC
        net_inputs(1,  3, i, j) = net_inputs(1,  3, i, j)  ! ACE2-EAMv3: ICEFRAC
        net_inputs(1,  4, i, j) = net_inputs(1, 4, i, j)   ! ACE2-EAMv3: PHIS
        ! -----------------------------------------------------------------------
        ! TODO (AN): Evolve `SOLIN`, `PS`, and TS fileds intime
        ! -----------------------------------------------------------------------
        net_inputs(1,  5, i, j) = net_inputs(1, 5, i, j)   ! ACE2-EAMv3: SOLIN
        net_inputs(1,  6, i, j) = net_outputs(1, 1, i, j)  ! ACE2-EAMv3: PS
        ! use landfrac as weights...
        net_inputs(1,  7, i, j) = net_outputs(1, 2, i, j)  ! ACE2-EAMv3: TS
        ! For 3D fields just advance through with time
        net_inputs(1,  8, i, j) = net_outputs(1,  3, i, j) ! ACE2-EAMv3: T_0
        net_inputs(1,  9, i, j) = net_outputs(1,  4, i, j) ! ACE2-EAMv3: T_1
        net_inputs(1, 10, i, j) = net_outputs(1,  5, i, j) ! ACE2-EAMv3: T_2
        net_inputs(1, 11, i, j) = net_outputs(1,  6, i, j) ! ACE2-EAMv3: T_3
        net_inputs(1, 12, i, j) = net_outputs(1,  7, i, j) ! ACE2-EAMv3: T_4
        net_inputs(1, 13, i, j) = net_outputs(1,  8, i, j) ! ACE2-EAMv3: T_5
        net_inputs(1, 14, i, j) = net_outputs(1,  9, i, j) ! ACE2-EAMv3: T_6
        net_inputs(1, 15, i, j) = net_outputs(1, 10, i, j) ! ACE2-EAMv3: T_7
        net_inputs(1, 16, i, j) = net_outputs(1, 11, i, j) ! ACE2-EAMv3: specific_total_water_0
        net_inputs(1, 17, i, j) = net_outputs(1, 12, i, j) ! ACE2-EAMv3: specific_total_water_1
        net_inputs(1, 18, i, j) = net_outputs(1, 13, i, j) ! ACE2-EAMv3: specific_total_water_2
        net_inputs(1, 19, i, j) = net_outputs(1, 14, i, j) ! ACE2-EAMv3: specific_total_water_3
        net_inputs(1, 20, i, j) = net_outputs(1, 15, i, j) ! ACE2-EAMv3: specific_total_water_4
        net_inputs(1, 21, i, j) = net_outputs(1, 16, i, j) ! ACE2-EAMv3: specific_total_water_5
        net_inputs(1, 22, i, j) = net_outputs(1, 17, i, j) ! ACE2-EAMv3: specific_total_water_6
        net_inputs(1, 23, i, j) = net_outputs(1, 18, i, j) ! ACE2-EAMv3: specific_total_water_7
        net_inputs(1, 24, i, j) = net_outputs(1, 19, i, j) ! ACE2-EAMv3: U_0
        net_inputs(1, 25, i, j) = net_outputs(1, 20, i, j) ! ACE2-EAMv3: U_1
        net_inputs(1, 26, i, j) = net_outputs(1, 21, i, j) ! ACE2-EAMv3: U_2
        net_inputs(1, 27, i, j) = net_outputs(1, 22, i, j) ! ACE2-EAMv3: U_3
        net_inputs(1, 28, i, j) = net_outputs(1, 23, i, j) ! ACE2-EAMv3: U_4
        net_inputs(1, 29, i, j) = net_outputs(1, 24, i, j) ! ACE2-EAMv3: U_5
        net_inputs(1, 30, i, j) = net_outputs(1, 25, i, j) ! ACE2-EAMv3: U_6
        net_inputs(1, 31, i, j) = net_outputs(1, 26, i, j) ! ACE2-EAMv3: U_7
        net_inputs(1, 32, i, j) = net_outputs(1, 27, i, j) ! ACE2-EAMv3: V_0
        net_inputs(1, 33, i, j) = net_outputs(1, 28, i, j) ! ACE2-EAMv3: V_1
        net_inputs(1, 34, i, j) = net_outputs(1, 29, i, j) ! ACE2-EAMv3: V_2
        net_inputs(1, 35, i, j) = net_outputs(1, 30, i, j) ! ACE2-EAMv3: V_3
        net_inputs(1, 36, i, j) = net_outputs(1, 31, i, j) ! ACE2-EAMv3: V_4
        net_inputs(1, 37, i, j) = net_outputs(1, 32, i, j) ! ACE2-EAMv3: V_5
        net_inputs(1, 38, i, j) = net_outputs(1, 33, i, j) ! ACE2-EAMv3: V_6
        net_inputs(1, 39, i, j) = net_outputs(1, 34, i, j) ! ACE2-EAMv3: V_7
      enddo
    enddo

    write(logunit_atm, *) "----------------------------------------------------------------"
    write(logunit_atm, *) "ace_eatm_import"
    write(logunit_atm, *) "----------------------------------------------------------------"
    write(logunit_atm, *) "ts  (min, max):   ( ", minval(ts(:, :)),  maxval(ts(:, :)), " )"
    call shr_sys_flush(logunit_atm)

  end subroutine ace_eatm_import

  subroutine ace_eatm_export(ggrid)
    ! !LOCAL VARIABLES:
    type(mct_gGrid), pointer :: ggrid
    real(R8),        pointer :: yc(:)

    integer(IN) :: klat
    integer     :: i, j, n
    real(R8)    :: e, avg_alb
    real(R8), parameter :: degtorad = SHR_CONST_PI/180.0_R8

    allocate(yc(lsize))

    klat = mct_aVect_indexRA(ggrid%data,'lat')
    yc(:) = ggrid%data%rAttr(klat,:)

    n = 0
    do j = 1, lsize_y
      do i = 1, lsize_x
        n = n + 1

        zbot(i, j) = 10.0_R8
        ubot(i, j) = net_outputs(1, 26, i, j) ! U_7
        vbot(i, j) = net_outputs(1, 34, i, j) ! V_7
        tbot(i, j) = net_outputs(1, 10, i, j) ! T_7
        ptem(i, j) = net_outputs(1, 10, i, j) ! T_7
        ! do we want to make this correspond to level 0 by using ak_i and bk_i
        pbot(i, j) = net_outputs(1,  1, i, j) ! PS (Surface pressure)
        pslv(i, j) = net_outputs(1,  1, i, j) ! PS (Surface pressure)
        lwdn(i, j) = net_outputs(1, 40, i, j) ! FLDS (Downwelling longwave flux at surface)

        !--- saturation vapor pressure ---
        e = datm_shr_esat(tbot(i, j), tbot(i, j))
        !--- specific humidity ---
        shum(i, j) = (0.622_R8 * e)/(pbot(i, j) - 0.378_R8 * e)
        !--- density ---
        dens(i, j) = pbot(i, j)  / (rdair * tbot(i, j) * (1 + 0.608_R8 * shum(i, j)))

        snowc(i, j) = 0.0_R8
        rainc(i, j) = 0.0_R8
        if (tbot(i, j) < tKFrz) then
          rainl(i, j) = 0.0_R8
          snowl(i, j) = max(net_outputs(1, 37, i, j), 0.0_R8)
        else
          rainl(i, j) = max(net_outputs(1, 37, i, j), 0.0_R8)
          snowl(i, j) = 0.0_R8
        endif

        ! Downwelling solar flux at surface
        swnet(i, j) = net_outputs(1, 41, i, j)
        !--- fabricate required sw[n,v]d[r,f] components from swnet ---
        swvdr(i, j) = swnet(i, j) * 0.28_R8
        swndr(i, j) = swnet(i, j) * 0.31_R8
        swvdf(i, j) = swnet(i, j) * 0.24_R8
        swndf(i, j) = swnet(i, j) * 0.17_R8

        ! avg_alb = ( 0.069 - 0.011*cos(2.0_R8*yc(n)*degtorad ) )
        ! swnet(i, j) = swnet(i, j) * (1.0_R4 - REAL(avg_alb, R4))
      enddo
    enddo

    deallocate(yc)

    write(logunit_atm, *) "----------------------------------------------------------------"
    write(logunit_atm, *) "ace_eatm_export"
    write(logunit_atm, *) "----------------------------------------------------------------"
    write(logunit_atm, *) "zbot  (min, max):   ( ", minval(zbot(:, :)),  maxval(zbot(:, :)), " )"
    write(logunit_atm, *) "tbot   (min, max):  ( ", minval(tbot(:, :)),  maxval(tbot(:, :)), " )"
    write(logunit_atm, *) "pbot   (min, max):  ( ", minval(pbot(:, :)),  maxval(pbot(:, :)), " )"
    write(logunit_atm, *) "ubot   (min, max):  ( ", minval(ubot(:, :)),  maxval(ubot(:, :)), " )"
    write(logunit_atm, *) "vbot   (min, max):  ( ", minval(vbot(:, :)),  maxval(vbot(:, :)), " )"
    write(logunit_atm, *) "swnet  (min, max):  ( ", minval(swnet(:, :)), maxval(swnet(:, :)), " )"
    write(logunit_atm, *) "rainl (min, max):   ( ", minval(rainl(:, :)), maxval(rainl(:, :)), " )"
    write(logunit_atm, *) "snowl (min, max):   ( ", minval(snowl(:, :)), maxval(snowl(:, :)), " )"
    call shr_sys_flush(logunit_atm)

  end subroutine ace_eatm_export

  subroutine ace_read_netcdf(varname, data, ncid)
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: varname             ! variable name
    real(R4)        , intent(inout) :: data(:,:)           ! raw data
    ! !LOCAL VARIABLES:
    logical                    :: found
    character(len=*),parameter :: subname = '(ace_read_netcdf) '

    call ncd_io(varname=varname, data=data, flag='read', ncid=ncid, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: reading ' //trim(varname))
  end subroutine ace_read_netcdf

  !===============================================================================
  real(R8) function datm_shr_eSat(tK,tKbot)

    !--- arguments ---
    real(R8),intent(in) :: tK    ! temp used in polynomial calculation
    real(R8),intent(in) :: tKbot ! bottom atm temp

    !--- local ---
    real(R8)           :: t     ! tK converted to Celcius

    !--- coefficients for esat over water ---
    real(R8),parameter :: a0=6.107799961_R8
    real(R8),parameter :: a1=4.436518521e-01_R8
    real(R8),parameter :: a2=1.428945805e-02_R8
    real(R8),parameter :: a3=2.650648471e-04_R8
    real(R8),parameter :: a4=3.031240396e-06_R8
    real(R8),parameter :: a5=2.034080948e-08_R8
    real(R8),parameter :: a6=6.136820929e-11_R8

    !--- coefficients for esat over ice ---
    real(R8),parameter :: b0=6.109177956_R8
    real(R8),parameter :: b1=5.034698970e-01_R8
    real(R8),parameter :: b2=1.886013408e-02_R8
    real(R8),parameter :: b3=4.176223716e-04_R8
    real(R8),parameter :: b4=5.824720280e-06_R8
    real(R8),parameter :: b5=4.838803174e-08_R8
    real(R8),parameter :: b6=1.838826904e-10_R8

    !----------------------------------------------------------------------------
    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    ! required to convert relative humidity to specific humidity
    !----------------------------------------------------------------------------

    t = min( 50.0_R8, max(-50.0_R8,(tK-tKfrz)) )
    if ( tKbot < tKfrz) then
       datm_shr_eSat = 100.0_R8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    else
       datm_shr_eSat = 100.0_R8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    end if

  end function datm_shr_eSat

  ! Define here becasue we need the eatmIO mod and trying to avoid circular imports
  subroutine init_normalizer(norm, norm_file, n)
    implicit none
    class(t_normalization_struct), intent(out) :: norm
    character(len=*),   intent(in)  :: norm_file
    integer,            intent (in) :: n  ! number of variables
    ! !LOCAL VARIABLES:
    type(file_desc_t)          :: ncid    ! netcdf file id
    logical                    :: found
    character(len=*),parameter :: subname = '(init_normalizer) '

    allocate(norm%stds(n))
    allocate(norm%means(n))

    call ncd_pio_openfile(ncid, trim(norm_file), 0)

    call ncd_io(varname='stds', data=norm%stds, flag='read', ncid=ncid, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: reading -- stds -- from ' // trim(norm_file))

    call ncd_io(varname='means', data=norm%means, flag='read', ncid=ncid, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: reading -- means -- from ' // trim(norm_file))

    call ncd_pio_closefile(ncid)
  end subroutine init_normalizer

  subroutine finalize_normalizer(norm)
    implicit none
    class(t_normalization_struct), intent(inout) :: norm

    deallocate(norm%stds)
    deallocate(norm%means)

  end subroutine finalize_normalizer
end module ace_comp_mod
