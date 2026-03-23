module eatm_restart_mod

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_abort

  use eatmIO
  use eatmMod

  implicit none
  save

  public :: eatm_restart_file_write
  public :: eatm_restart_file_read
  public :: eatm_initial_condition_file_read

  private

  contains
    subroutine eatm_restart_file_write( file, rdate, stepno )

      !----------------------------------------------------------------
      ! !DESCRIPTION:
      ! ...
      implicit none
      ! !ARGUMENTS
      character(len=*), intent(in) :: file
      character(len=*), intent(in) :: rdate
      integer, intent(in) :: stepno

      ! !LOCAL VARIABLES:
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i ! index
      !----------------------------------------------------------------

      write(logunit_atm,'(72a1)') ("-",i=1,60)
      write(logunit_atm, *) 'restart_file_open: writing EATM restart dataset '
      write(logunit_atm, *)

      ! Define dimensions and variables
      call ncd_pio_createfile(ncid, trim(file))
      call set_restart_file_dimensions(ncid)
      call eatm_restart(ncid, 'define')
      call ncd_enddef(ncid)

      ! Write restart file variables
      call eatm_restart(ncid, 'write')
      call ncd_pio_closefile(ncid)

      write(logunit_atm, *)
      write(logunit_atm, *) 'Successfully wrote out restart data at nstep = ',stepno
      write(logunit_atm,'(72a1)') ("-",i=1,60)

    end subroutine eatm_restart_file_write

    subroutine eatm_restart_file_read( file )
      !----------------------------------------------------------------
      ! !DESCRIPTION:
      ! ...
      implicit none
      ! !ARGUMENTS
      character(len=*), intent(in) :: file

      ! !LOCAL VARIABLES:
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i ! index
      !----------------------------------------------------------------

      write(logunit_atm, *) 'Reading restart dataset'
      call ncd_pio_openfile(ncid, trim(file), 0)
      call eatm_restart(ncid, 'read')
      call ncd_pio_closefile(ncid)

      write(logunit_atm, *)
      write(logunit_atm, *) 'Successfully read restart data for restart run'
      write(logunit_atm,'(72a1)') ("-",i=1,60)

    end subroutine eatm_restart_file_read

    subroutine eatm_initial_condition_file_read( file )
      !----------------------------------------------------------------
      ! !DESCRIPTION:
      ! ...
      implicit none
      ! !ARGUMENTS
      character(len=*), intent(in) :: file

      ! !LOCAL VARIABLES:
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i ! index
      !----------------------------------------------------------------

      write(logunit_atm, *) 'Reading initial condition dataset'
      call ncd_pio_openfile(ncid, trim(file), 0)
      call eatm_initial_condition(ncid, 'read')
      call ncd_pio_closefile(ncid)

      write(logunit_atm, *)
      write(logunit_atm, *) 'Successfully read initial condition data for startup run'
      write(logunit_atm,'(72a1)') ("-",i=1,60)

    end subroutine eatm_initial_condition_file_read

    subroutine set_restart_file_dimensions( ncid )

      !----------------------------------------------------------------
      ! !DESCRIPTION:
      ! ...
      implicit none
      ! !ARGUMENTS
      type(file_desc_t) :: ncid ! netcdf id

      ! !LOCAL VARIABLES:
      integer :: dimid               ! netCDF dimension id
      integer :: ier                 ! error status
      character(len=  8) :: curdate  ! current date
      character(len=  8) :: curtime  ! current time
      character(len=256) :: str
      !----------------------------------------------------------------

      ! Define dimensions
      call ncd_defdim(ncid, 'lon', lsize_x, dimid)
      call ncd_defdim(ncid, 'lat', lsize_y, dimid)
      call ncd_defdim(ncid, 'time', 2, dimid)

      ! Define global attributes
    end subroutine set_restart_file_dimensions

    subroutine eatm_restart( ncid, flag )
      !-----------------------------------------------------------------------
      ! DESCRIPTION:
      ! define/read/write eatm restart data.

      ! ARGUMENTS:
      implicit none
      type(file_desc_t), intent(inout) :: ncid ! netcdf id
      character(len=*) , intent(in)    :: flag ! 'define' or 'read' or 'write'

      ! LOCAL VARIABLES:
      integer :: c ! channel index
      logical :: readvar ! determine if variable is read
      character(len=1)   :: level
      character(len=32)  :: vname
      character(len=32)  :: uname
      character(len=256) :: lname

      do c = 1, 44
        if (c == 1) then
          vname = "PS"
          lname = "Surface pressure"
          uname = "Pa"
        elseif (c == 2) then
          vname = "TS"
          lname = "Surface temperature (radiative)"
          uname = "K"
        elseif ((c > 2) .and. (c < 11)) then
          level = achar( iachar('0') + mod(c - 3, 8) )
          vname = "T_"//level
          lname = "Temperature level-"//level
          uname = "K"
        elseif ((c > 10) .and. (c < 19)) then
          level = achar( iachar('0') + mod(c - 3, 8) )
          vname = "specific_total_water_"//level
          lname = "specific total water level-"//level
          uname = "kg/kg"
        elseif ((c > 18) .and. (c < 27)) then
          level = achar( iachar('0') + mod(c - 3, 8) )
          vname = "U_"//level
          lname = "Zonal wind level-"//level
          uname = "m/s"
        elseif ((c > 26) .and. (c < 35)) then
          level = achar( iachar('0') + mod(c - 3, 8) )
          vname = "V_"//level
          lname = "Meridional wind level-"//level
          uname = "m/s"
        elseif (c == 35) then
          vname = "LHFLX"
          lname = "Surface latent heat flux"
          uname = "W/m2"
        elseif (c == 36) then
          vname = "SHFLX"
          lname = "Surface sensible heat flux"
          uname = "W/m2"
        elseif (c == 37) then
          vname = "surface_precipitation_rate"
          lname = "surface precipitation rate (all phases)"
          uname = "kg/m2/s"
        elseif (c == 38) then
          vname = "surface_upward_longwave_flux"
          lname = "Upward longwave radiative flux at surface"
          uname = "W/m2"
        elseif (c == 39) then
          vname = "FLUT"
          lname = "Upward longwave radiative flux at TOA"
          uname = "W/m2"
        elseif (c == 40) then
          vname = "FLDS"
          lname = "Downward longwave radiative flux at surface"
          uname = "W/m2"
        elseif (c == 41) then
          vname = "FSDS"
          lname = "Downward shortwave radiative flux at surface"
          uname = "W/m2"
        elseif (c == 42) then
          vname = "surface_upward_shortwave_flux"
          lname = "Upward shortwave radiative flux at surface"
          uname = "W/m2"
        elseif (c == 43) then
          vname = "top_of_atmos_upward_shortwave_flux"
          lname = "Upward shortwave radiative flux at TOA"
          uname = "W/m2"
        elseif (c == 44) then
          vname = "tendency_of_total_water_path_due_to_advection"
          lname = "Tendency of total water path form advection"
          uname = "kg/m2/s"
        else
          print *, "Something is wrong"
        endif

        if ( flag == 'define' ) then
          call ncd_defvar(&
            ncid=ncid, &
            varname=trim(vname), &
            xtype=ncd_double, &
            dim1name='lon', &
            dim2name='lat', &
            dim3name='time', &
            long_name=trim(lname), &
            units=trim(uname) &
        )
        elseif (flag == 'read' .or. flag == 'write') then
          call ncd_io(&
            varname=trim(vname), &
            data=eatm_intrp%t_im1(c, :, :), &
            ncid=ncid, &
            flag=flag, &
            nt=1, &
            readvar=readvar &
          )
          call ncd_io(&
            varname=trim(vname), &
            data=eatm_intrp%t_ip1(c, :, :), &
            ncid=ncid, &
            flag=flag, &
            nt=2, &
            readvar=readvar &
          )
        endif
      enddo

      if ( flag == 'define' ) then
        call ncd_defvar(&
          ncid=ncid, &
          varname="PHIS", &
          xtype=ncd_double, &
          dim1name="lon", &
          dim2name="lat", &
          long_name="Surface geopotential", &
          units="m2/s2" &
        )
        call ncd_defvar(&
          ncid=ncid, &
          varname="SOLIN", &
          xtype=ncd_double, &
          dim1name="lon", &
          dim2name="lat", &
          long_name="Solar insolation", &
          units="W/m2" &
        )
      elseif (flag == 'read' .or. flag == 'write') then
        call ncd_io(&
          varname="PHIS", &
          data=net_inputs(1, 4, :, :), &
          ncid=ncid, &
          flag=flag, &
          readvar=readvar &
        )
        call ncd_io(&
          varname="SOLIN", &
          data=net_inputs(1, 5, :, :), &
          ncid=ncid, &
          flag=flag, &
          readvar=readvar &
        )
      endif
    end subroutine eatm_restart

    subroutine eatm_initial_condition( ncid, flag )
      !-----------------------------------------------------------------------
      ! DESCRIPTION:
      ! read eatm initial condition data.

      ! ARGUMENTS:
      implicit none
      type(file_desc_t), intent(inout) :: ncid ! netcdf id
      character(len=*) , intent(in)    :: flag ! 'define' or 'read' or 'write'

      ! LOCAL VARIABLES:
      integer :: c ! channel index
      logical :: readvar ! determine if variable is read
      character(len=1)   :: level
      character(len=32)  :: vname
      character(len=32)  :: uname
      character(len=256) :: lname

      integer :: oc               ! output channel index
      character(len=32) :: oname  ! output variable name

      do c = 1, 39
        if (c == 1) then
          vname = "LANDFRAC"
          lname = "Fraction of sfc area covered by land"
          uname = "1"
        elseif (c == 2) then
          vname = "OCNFRAC"
          lname = "Fraction of sfc area covered by ocean"
          uname = "1"
        elseif (c == 3) then
          vname = "ICEFRAC"
          lname = "Fraction of sfc area covered by sea-ice"
          uname = "1"
        elseif (c == 4) then
          vname = "PHIS"
          lname = "Surface geopotential"
          uname = "m2/s2"
        elseif (c == 5) then
          vname = "SOLIN"
          lname = "Solar insolation"
          uname = "W/m2"
        elseif (c == 6) then
          vname = "PS"
          lname = "Surface pressure"
          uname = "Pa"
        elseif (c == 7) then
          vname = "TS"
          lname = "Surface temperature (radiative)"
          uname = "K"
        elseif (c > 7 .and. c < 16) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "T_"//level
          lname = "Temperature level-"//level
          uname = "K"
        elseif (c > 15 .and. c < 24) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "specific_total_water_"//level
          lname = "specific total water level-"//level
          uname = "kg/kg"
        elseif (c > 23 .and. c < 32) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "U_"//level
          lname = "Zonal wind level-"//level
          uname = "m/s"
        elseif (c > 31 .and. c < 40) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "V_"//level
          lname = "Meridional wind level-"//level
          uname = "m/s"
        else
          print *, "Something is wrong"
        endif

        if (flag == 'read') then
          call ncd_io(&
            varname=trim(vname), &
            data=net_inputs(1, c, :, :), &
            ncid=ncid, &
            flag=flag, &
            readvar=readvar &
          )
        else
            call shr_sys_abort()
        endif
      enddo
    end subroutine eatm_initial_condition
end module eatm_restart_mod
