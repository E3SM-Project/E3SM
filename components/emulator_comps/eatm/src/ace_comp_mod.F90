module ace_comp_mod

  use eatmMod
  use eatmIO
  use shr_const_mod
  use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod,  only: shr_sys_flush, shr_sys_abort

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------
  public :: ace_comp_init

  !--------------------------------------------------------------------------
  ! Public module data
  !--------------------------------------------------------------------------
  integer, public :: n_input_channels=39  ! number of input channels to emulator
  integer, public :: n_output_channels=44 ! number of input channels to emulator

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  real(R8), parameter :: tKFrz  = SHR_CONST_TKFRZ

  save

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS

  subroutine ace_comp_init()

    ! load the initial condition data into input array
    call ace_read_ic()
    ! TODO: create a tensor from the input array

    ! load the restart data (i.e. prognostic state) associated with IC
    call ace_read_restart()
    ! TODO: create a tensor from the output array

    ! using restart data from ACE set the fields passed to the coupler
    call ace_eatm_export()

  end subroutine ace_comp_init

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

  subroutine ace_eatm_export()
    ! !LOCAL VARIABLES:
    integer :: i, j

    do i = 1, lsize_x
      do j = 1, lsize_y
        zbot(i, j) = 10.0_R8
        ubot(i, j) = net_outputs(1, 26, i, j) ! U_7
        vbot(i, j) = net_outputs(1, 34, i, j) ! V_7
        tbot(i, j) = net_outputs(1, 10, i, j) ! T_7
        ptem(i, j) = net_outputs(1, 10, i, j) ! T_7
        ! do we want to make this correspond to level 0 by using ak_i and bk_i
        pbot(i, j) = net_outputs(1,  1, i, j) ! PS (Surface pressure)
        pslv(i, j) = net_outputs(1,  1, i, j) ! PS (Surface pressure)
        lwdn(i, j) = net_outputs(1, 40, i, j) ! FLDS (Downwelling longwave flux at surface)

        snowc(i, j) = 0.0_R8
        rainc(i, j) = 0.0_R8
        if (tbot(i, j) < tKFrz) then
          rainl(i, j) = 0.0_R8
          snowl(i, j) = net_outputs(1, 37, i, j)
        else
          rainl(i, j) = net_outputs(1, 37, i, j)
          snowl(i, j) = 0.0_R8
        endif

        ! swnet = (Downwelling solar flux at surface) - (surface upward shortwave flux)
        swnet(i, j) = net_outputs(1, 41, i, j) - net_outputs(1, 42, i, j)
        !--- fabricate required sw[n,v]d[r,f] components from swnet ---
        swvdr(i, j) = swnet(i, j) * 0.28_R8
        swndr(i, j) = swnet(i, j) * 0.31_R8
        swvdf(i, j) = swnet(i, j) * 0.24_R8
        swndf(i, j) = swnet(i, j) * 0.17_R8
      enddo
    enddo
  end subroutine ace_eatm_export

  subroutine ace_read_netcdf(varname, data, ncid)
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid                ! netcdf file id
    character(len=*), intent(in)    :: varname             ! variable name
    real(R8)        , intent(inout) :: data(:,:)           ! raw data
    ! !LOCAL VARIABLES:
    logical                    :: found
    character(len=*),parameter :: subname = '(ace_read_netcdf) '

    call ncd_io(varname=varname, data=data, flag='read', ncid=ncid, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: reading ' //trim(varname))
  end subroutine ace_read_netcdf

end module ace_comp_mod
