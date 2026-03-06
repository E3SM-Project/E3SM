module eatm_restart_mod

  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_sys_mod       , only : shr_sys_abort

  use eatmIO
  use eatmMod

  implicit none
  save

  public :: eatm_restart_file_write
  public :: eatm_restart_file_read

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
        elseif (c .gt. 7 .and. c .lt. 16) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "T_"//level
          lname = "Temperature level-"//level
          uname = "K"
        elseif (c .gt. 15 .and. c .lt. 24) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "specific_total_water_"//level
          lname = "specific total water level-"//level
          uname = "kg/kg"
        elseif (c .gt. 23 .and. c .lt. 32) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "U_"//level
          lname = "Zonal wind level-"//level
          uname = "m/s"
        elseif (c .gt. 31 .and. c .lt. 40) then
          level = achar( iachar('0') + mod(c, 8) )
          vname = "V_"//level
          lname = "Meridional wind level-"//level
          uname = "m/s"
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
            long_name=trim(lname), &
            units=trim(uname) &
        )
        elseif (flag == 'read' .or. flag == 'write') then
          call ncd_io(&
            varname=trim(vname), &
            data=net_inputs(1, c, :, :), &
            ncid=ncid, &
            flag=flag, &
            readvar=readvar &
          )
        endif

        !! TODO: abort is we were supposed to read and we didn't
        !if (flag=='read' .and. .not. readvar .and. nsrest == nsrContinue) then
        !    call shr_sys_abort()
        !endif
      enddo
    end subroutine eatm_restart
end module eatm_restart_mod
