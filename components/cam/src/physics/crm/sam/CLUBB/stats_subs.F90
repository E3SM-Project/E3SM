!-----------------------------------------------------------------------
!  $Id: stats_subs.F90 6146 2013-04-05 18:02:22Z raut@uwm.edu $
module stats_subs

  implicit none

  private ! Set Default Scope

  public :: stats_init, stats_begin_timestep, stats_end_timestep, & 
    stats_accumulate, stats_finalize, stats_accumulate_hydromet, &
    stats_accumulate_LH_tend

  private :: stats_zero, stats_avg

  contains

  !-----------------------------------------------------------------------
  subroutine stats_init( iunit, fname_prefix, fdir, l_stats_in, &
                         stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                         nzmax, gzt, gzm, nnrad_zt, &
                         grad_zt, nnrad_zm, grad_zm, day, month, year, &
                         rlat, rlon, time_current, delt )
    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
      zt,      & ! Variables
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      ztscr07, & 
      ztscr08, & 
      ztscr09, & 
      ztscr10, & 
      ztscr11, & 
      ztscr12, & 
      ztscr13, & 
      ztscr14, & 
      ztscr15, & 
      ztscr16, & 
      ztscr17, & 
      ztscr18, & 
      ztscr19, & 
      ztscr20, & 
      ztscr21

    use stats_variables, only: & 
      LH_zt, &  ! Variable(s)
      LH_sfc

    use stats_variables, only: & 
      zm,      & ! Variables
      zmscr01, & 
      zmscr02, & 
      zmscr03, & 
      zmscr04, & 
      zmscr05, & 
      zmscr06, & 
      zmscr07, & 
      zmscr08, & 
      zmscr09, & 
      zmscr10, & 
      zmscr11, & 
      zmscr12, & 
      zmscr13, & 
      zmscr14, & 
      zmscr15, &
      zmscr16, &
      zmscr17, &
      rad_zt

    use stats_variables, only: &
      rad_zm,  &
      sfc,     & 
      l_stats, &
      l_output_rad_files, & 
      stats_tsamp,   & 
      stats_tout,    & 
      l_stats_samp,  & 
      l_stats_last, & 
      fname_zt, & 
      fname_LH_zt, & 
      fname_LH_sfc, & 
      fname_zm, &
      fname_rad_zt, &
      fname_rad_zm, & 
      fname_sfc, & 
      l_netcdf, & 
      l_grads

    use clubb_precision, only: & 
      time_precision, & ! Constant(s)
      core_rknd

    use output_grads, only: & 
      open_grads ! Procedure

#ifdef NETCDF
    use output_netcdf, only: & 
      open_netcdf     ! Procedure
#endif

    use stats_zm, only: &
      nvarmax_zm, & ! Constant(s) 
      stats_init_zm ! Procedure(s)

    use stats_zt, only: & 
      nvarmax_zt, & ! Constant(s)
      stats_init_zt ! Procedure(s)

    use stats_LH_zt, only: & 
      nvarmax_LH_zt, & ! Constant(s)
      stats_init_LH_zt ! Procedure(s)

    use stats_LH_sfc, only: & 
      nvarmax_LH_sfc, & ! Constant(s)
      stats_init_LH_sfc ! Procedure(s)

    use stats_rad_zt, only: & 
      nvarmax_rad_zt, & ! Constant(s)
      stats_init_rad_zt ! Procedure(s)

    use stats_rad_zm, only: & 
      nvarmax_rad_zm, & ! Constant(s)
      stats_init_rad_zm ! Procedure(s)

    use stats_sfc, only: &
      nvarmax_sfc, & ! Constant(s)
      stats_init_sfc ! Procedure(s)

    use error_code, only: &
      clubb_at_least_debug_level ! Function

    use constants_clubb, only: &
      fstdout, fstderr, var_length ! Constants

    use parameters_microphys, only: &
      LH_microphys_disabled, & ! Constant
      LH_microphys_type ! Variable

    implicit none

    ! Input Variables

    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: nzmax ! Grid points in the vertical [count]

    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  & 
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(1), intent(in) ::  & 
      rlat, rlon   ! Latitude and Longitude             [Degrees N/E]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Model time                         [s]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]


    ! Local Variables
    logical :: l_error

    character(len=200) :: fname

    integer :: i, ntot, read_status

    ! Namelist Variables

    character(len=10) :: stats_fmt  ! File storage convention

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_LH_zt) ::  & 
      vars_LH_zt  ! Latin Hypercube variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_LH_sfc) ::  & 
      vars_LH_sfc  ! Latin Hypercube variables at the surface

    character(len=var_length), dimension(nvarmax_zm) ::  & 
      vars_zm  ! Variables on the momentum levels

    character(len=var_length), dimension(nvarmax_rad_zt) ::  & 
      vars_rad_zt  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_rad_zm) ::  & 
      vars_rad_zm  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_sfc) ::  &
      vars_sfc ! Variables at the model surface

    namelist /statsnl/ & 
      vars_zt, & 
      vars_zm, &
      vars_LH_zt, &
      vars_LH_sfc, &
      vars_rad_zt, &
      vars_rad_zm, & 
      vars_sfc

    ! ---- Begin Code ----

    ! Initialize
    l_error = .false.

    ! Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    stats_tsamp = stats_tsamp_in
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in
    stats_fmt   = trim( stats_fmt_in )

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_last  = .false.
      return
    end if

    ! Initialize namelist variables

    vars_zt  = ''
    vars_zm  = ''
    vars_LH_zt = ''
    vars_LH_sfc = ''
    vars_rad_zt = ''
    vars_rad_zm = ''
    vars_sfc = ''

    ! Reads list of variables that should be output to GrADS/NetCDF (namelist &statsnl)

    open(unit=iunit, file=fnamelist)
    read(unit=iunit, nml=statsnl, iostat=read_status, end=100)
    if ( read_status /= 0 ) then
      if ( read_status > 0 ) then
        write(fstderr,*) "Error reading stats namelist in file ",  &
                         trim( fnamelist )
      else ! Read status < 0
        write(fstderr,*) "End of file marker reached while reading stats namelist in file ", &
          trim( fnamelist )
      end if
      write(fstderr,*) "One cause is having more statistical variables ",  &
                       "listed in the namelist for var_zt, var_zm, or ",  &
                       "var_sfc than allowed by nvarmax_zt, nvarmax_zm, ",  &
                       "or nvarmax_sfc, respectively."
      write(fstderr,*) "Maximum variables allowed for var_zt = ", nvarmax_zt
      write(fstderr,*) "Maximum variables allowed for var_zm = ", nvarmax_zm
      write(fstderr,*) "Maximum variables allowed for var_rad_zt = ", nvarmax_rad_zt
      write(fstderr,*) "Maximum variables allowed for var_rad_zm = ", nvarmax_rad_zm
      write(fstderr,*) "Maximum variables allowed for var_sfc = ", nvarmax_sfc
      stop "stats_init: Error reading stats namelist."
    end if ! read_status /= 0

    close(unit=iunit)

    if ( clubb_at_least_debug_level( 1 ) ) then
      write(fstdout,*) "--------------------------------------------------"

      write(fstdout,*) "Statistics"

      write(fstdout,*) "--------------------------------------------------"
      write(fstdout,*) "vars_zt = "
      i = 1
      do while ( vars_zt(i) /= '' )
        write(fstdout,*) vars_zt(i)
        i = i + 1
      end do

      write(fstdout,*) "vars_zm = "
      i = 1
      do while ( vars_zm(i) /= '' )
        write(fstdout,*) vars_zm(i)
        i = i + 1
      end do

      if ( LH_microphys_type /= LH_microphys_disabled ) then
        write(fstdout,*) "vars_LH_zt = "
        i = 1
        do while ( vars_LH_zt(i) /= '' )
          write(fstdout,*) vars_LH_zt(i)
          i = i + 1
        end do

        write(fstdout,*) "vars_LH_sfc = "
        i = 1
        do while ( vars_LH_sfc(i) /= '' )
          write(fstdout,*) vars_LH_sfc(i)
          i = i + 1
        end do
      end if ! LH_microphys_type /= LH_microphys_disabled

      if ( l_output_rad_files ) then
        write(fstdout,*) "vars_rad_zt = "
        i = 1
        do while ( vars_rad_zt(i) /= '' )
          write(fstdout,*) vars_rad_zt(i)
          i = i + 1
        end do

        write(fstdout,*) "vars_rad_zm = "
        i = 1
        do while ( vars_rad_zm(i) /= '' )
          write(fstdout,*) vars_rad_zm(i)
          i = i + 1
        end do
      end if ! l_output_rad_files

      write(fstdout,*) "vars_sfc = "
      i = 1
      do while ( vars_sfc(i) /= '' )
        write(fstdout,*) vars_sfc(i)
        i = i + 1
      end do

      write(fstdout,*) "--------------------------------------------------"
    end if ! clubb_at_least_debug_level 1

    ! Determine file names for GrADS or NetCDF files
    fname_zt  = trim( fname_prefix )//"_zt"
    fname_zm  = trim( fname_prefix )//"_zm"
    fname_LH_zt  = trim( fname_prefix )//"_LH_zt"
    fname_LH_sfc  = trim( fname_prefix )//"_LH_sfc"
    fname_rad_zt  = trim( fname_prefix )//"_rad_zt"
    fname_rad_zm  = trim( fname_prefix )//"_rad_zm"
    fname_sfc = trim( fname_prefix )//"_sfc"

    ! Parse the file type for stats output.  Currently only GrADS and
    ! netCDF > version 3.5 are supported by this code.
    select case ( trim( stats_fmt ) )
    case ( "GrADS", "grads", "gr" )
      l_netcdf = .false.
      l_grads  = .true.

    case ( "NetCDF", "netcdf", "nc" )
      l_netcdf = .true.
      l_grads  = .false.

    case default
      write(fstderr,*) "In module stats_subs subroutine stats_init: "
      write(fstderr,*) "Invalid stats output format "//trim( stats_fmt )
      stop "Fatal error"

    end select

    ! Check sampling and output frequencies

    ! The model time step length, delt (which is dt_main), should multiply
    ! evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - real( floor( stats_tsamp/delt ), kind=time_precision ) )  & 
           > 1.e-8_time_precision ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                       'delt (which is dt_main).  Check the appropriate ',  &
                       'model.in file.'
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
      write(fstderr,*) 'delt = ', delt
    end if

    ! The statistical sampling time step length, stats_tsamp, should multiply
    ! evenly into the statistical output time step length, stats_tout.
    if ( abs( stats_tout/stats_tsamp &
           - real( floor( stats_tout/stats_tsamp ), kind=time_precision ) ) & 
         > 1.e-8_time_precision ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tout should be an even multiple of ',  &
                       'stats_tsamp.  Check the appropriate model.in file.'
      write(fstderr,*) 'stats_tout = ', stats_tout
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
    end if

    ! Initialize zt (mass points)

    i = 1
    do while ( ichar(vars_zt(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zt(i)) /= 0 & 
               .and. i <= nvarmax_zt )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      stop "stats_init:  number of zt statistical variables exceeds limit"
    end if

    zt%nn = ntot
    zt%kk = nzmax

    allocate( zt%z( zt%kk ) )
    zt%z = gzt

    allocate( zt%x( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%n( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%l_in_update( 1, 1, zt%kk, zt%nn ) )
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )

    allocate( zt%f%var( zt%nn ) )
    allocate( zt%f%z( zt%kk ) )

    ! Allocate scratch space

    allocate( ztscr01(zt%kk) )
    allocate( ztscr02(zt%kk) )
    allocate( ztscr03(zt%kk) )
    allocate( ztscr04(zt%kk) )
    allocate( ztscr05(zt%kk) )
    allocate( ztscr06(zt%kk) )
    allocate( ztscr07(zt%kk) )
    allocate( ztscr08(zt%kk) )
    allocate( ztscr09(zt%kk) )
    allocate( ztscr10(zt%kk) )
    allocate( ztscr11(zt%kk) )
    allocate( ztscr12(zt%kk) )
    allocate( ztscr13(zt%kk) )
    allocate( ztscr14(zt%kk) )
    allocate( ztscr15(zt%kk) )
    allocate( ztscr16(zt%kk) )
    allocate( ztscr17(zt%kk) )
    allocate( ztscr18(zt%kk) )
    allocate( ztscr19(zt%kk) )
    allocate( ztscr20(zt%kk) )
    allocate( ztscr21(zt%kk) )

    ztscr01 = 0.0_core_rknd
    ztscr02 = 0.0_core_rknd
    ztscr03 = 0.0_core_rknd
    ztscr04 = 0.0_core_rknd
    ztscr05 = 0.0_core_rknd
    ztscr06 = 0.0_core_rknd
    ztscr07 = 0.0_core_rknd
    ztscr08 = 0.0_core_rknd
    ztscr09 = 0.0_core_rknd
    ztscr10 = 0.0_core_rknd
    ztscr11 = 0.0_core_rknd
    ztscr12 = 0.0_core_rknd
    ztscr13 = 0.0_core_rknd
    ztscr14 = 0.0_core_rknd
    ztscr15 = 0.0_core_rknd
    ztscr16 = 0.0_core_rknd
    ztscr17 = 0.0_core_rknd
    ztscr18 = 0.0_core_rknd
    ztscr19 = 0.0_core_rknd
    ztscr20 = 0.0_core_rknd
    ztscr21 = 0.0_core_rknd

    fname = trim( fname_zt )

    if ( l_grads ) then

      ! Open GrADS file
      call open_grads( iunit, fdir, fname,  & 
                       1, zt%kk, zt%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zt%nn, zt%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, zt%kk, zt%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, zt%nn, &  ! In
                        zt%f ) ! InOut
#else
      stop "This CLUBB program was not compiled with netCDF support."
#endif

    end if

    ! Default initialization for array indices for zt

    call stats_init_zt( vars_zt, l_error )


    ! Setup output file for LH_zt (Latin Hypercube stats)

    if ( LH_microphys_type /= LH_microphys_disabled ) then

      i = 1
      do while ( ichar(vars_LH_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_LH_zt(i)) /= 0 & 
                 .and. i <= nvarmax_LH_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_LH_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_zt than allowed for by nvarmax_LH_zt."
        write(fstderr,*) "Check the number of variables listed for vars_LH_zt ",  &
                         "in the stats namelist, or change nvarmax_LH_zt."
        write(fstderr,*) "nvarmax_LH_zt = ", nvarmax_LH_zt
        stop "stats_init:  number of LH_zt statistical variables exceeds limit"
      end if

      LH_zt%nn = ntot
      LH_zt%kk = nzmax

      allocate( LH_zt%z( LH_zt%kk ) )
      LH_zt%z = gzt

      allocate( LH_zt%x( 1, 1, LH_zt%kk, LH_zt%nn ) )
      allocate( LH_zt%n( 1, 1, LH_zt%kk, LH_zt%nn ) )
      allocate( LH_zt%l_in_update( 1, 1, LH_zt%kk, LH_zt%nn ) )
      call stats_zero( LH_zt%kk, LH_zt%nn, LH_zt%x, LH_zt%n, LH_zt%l_in_update )

      allocate( LH_zt%f%var( LH_zt%nn ) )
      allocate( LH_zt%f%z( LH_zt%kk ) )


      fname = trim( fname_LH_zt )

      if ( l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & 
                         1, LH_zt%kk, LH_zt%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         LH_zt%nn, LH_zt%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname, 1, LH_zt%kk, LH_zt%z, &  ! In
                          day, month, year, rlat, rlon, &  ! In
                          time_current+stats_tout, stats_tout, LH_zt%nn, &  ! In
                          LH_zt%f ) ! InOut
#else
        stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

      call stats_init_LH_zt( vars_LH_zt, l_error )

      i = 1
      do while ( ichar(vars_LH_sfc(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_LH_sfc(i)) /= 0 & 
                 .and. i <= nvarmax_LH_sfc )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_LH_sfc ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_zt than allowed for by nvarmax_LH_sfc."
        write(fstderr,*) "Check the number of variables listed for vars_LH_sfc ",  &
                         "in the stats namelist, or change nvarmax_LH_sfc."
        write(fstderr,*) "nvarmax_LH_sfc = ", nvarmax_LH_sfc
        stop "stats_init:  number of LH_sfc statistical variables exceeds limit"
      end if

      LH_sfc%nn = ntot
      LH_sfc%kk = 1

      allocate( LH_sfc%z( LH_sfc%kk ) )
      LH_sfc%z = gzm(1)

      allocate( LH_sfc%x( 1, 1, LH_sfc%kk, LH_sfc%nn ) )
      allocate( LH_sfc%n( 1, 1, LH_sfc%kk, LH_sfc%nn ) )
      allocate( LH_sfc%l_in_update( 1, 1, LH_sfc%kk, LH_sfc%nn ) )

      call stats_zero( LH_sfc%kk, LH_sfc%nn, LH_sfc%x, LH_sfc%n, LH_sfc%l_in_update )

      allocate( LH_sfc%f%var( LH_sfc%nn ) )
      allocate( LH_sfc%f%z( LH_sfc%kk ) )

      fname = trim( fname_LH_sfc )

      if ( l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & 
                         1, LH_sfc%kk, LH_sfc%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         LH_sfc%nn, LH_sfc%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname, 1, LH_sfc%kk, LH_sfc%z, &  ! In
                          day, month, year, rlat, rlon, &  ! In
                          time_current+stats_tout, stats_tout, LH_sfc%nn, &  ! In
                          LH_sfc%f ) ! InOut
#else
        stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

      call stats_init_LH_sfc( vars_LH_sfc, l_error )

    end if ! LH_microphys_type /= LH_microphys_disabled

    ! Initialize zm (momentum points)

    i = 1
    do while ( ichar(vars_zm(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zm(i)) /= 0 & 
               .and. i <= nvarmax_zm )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      stop "stats_init:  number of zm statistical variables exceeds limit"
    end if

    zm%nn = ntot
    zm%kk = nzmax

    allocate( zm%z( zm%kk ) )
    zm%z = gzm

    allocate( zm%x( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%n( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%l_in_update( 1, 1, zm%kk, zm%nn ) )

    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )

    allocate( zm%f%var( zm%nn ) )
    allocate( zm%f%z( zm%kk ) )

    ! Allocate scratch space

    allocate( zmscr01(zm%kk) )
    allocate( zmscr02(zm%kk) )
    allocate( zmscr03(zm%kk) )
    allocate( zmscr04(zm%kk) )
    allocate( zmscr05(zm%kk) )
    allocate( zmscr06(zm%kk) )
    allocate( zmscr07(zm%kk) )
    allocate( zmscr08(zm%kk) )
    allocate( zmscr09(zm%kk) )
    allocate( zmscr10(zm%kk) )
    allocate( zmscr11(zm%kk) )
    allocate( zmscr12(zm%kk) )
    allocate( zmscr13(zm%kk) )
    allocate( zmscr14(zm%kk) )
    allocate( zmscr15(zm%kk) )
    allocate( zmscr16(zm%kk) )
    allocate( zmscr17(zm%kk) )

    ! Initialize to 0
    zmscr01 = 0.0_core_rknd
    zmscr02 = 0.0_core_rknd
    zmscr03 = 0.0_core_rknd
    zmscr04 = 0.0_core_rknd
    zmscr05 = 0.0_core_rknd
    zmscr06 = 0.0_core_rknd
    zmscr07 = 0.0_core_rknd
    zmscr08 = 0.0_core_rknd
    zmscr09 = 0.0_core_rknd
    zmscr10 = 0.0_core_rknd
    zmscr11 = 0.0_core_rknd
    zmscr12 = 0.0_core_rknd
    zmscr13 = 0.0_core_rknd
    zmscr14 = 0.0_core_rknd
    zmscr15 = 0.0_core_rknd
    zmscr16 = 0.0_core_rknd
    zmscr17 = 0.0_core_rknd


    fname = trim( fname_zm )
    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, zm%kk, zm%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zm%nn, zm%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, zm%kk, zm%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, zm%nn, &  ! In
                        zm%f ) ! InOut

#else
      stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    call stats_init_zm( vars_zm, l_error )

    ! Initialize rad_zt (radiation points)

    if (l_output_rad_files) then

      i = 1
      do while ( ichar(vars_rad_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zt(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        stop "stats_init:  number of rad_zt statistical variables exceeds limit"
      end if

      rad_zt%nn = ntot
      rad_zt%kk = nnrad_zt

      allocate( rad_zt%z( rad_zt%kk ) )
      rad_zt%z = grad_zt

      allocate( rad_zt%x( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%n( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%l_in_update( 1, 1, rad_zt%kk, rad_zt%nn ) )

      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )

      allocate( rad_zt%f%var( rad_zt%nn ) )
      allocate( rad_zt%f%z( rad_zt%kk ) )

      ! Allocate scratch space

      !allocate( radscr01(rad%kk) )
      !allocate( radscr02(rad%kk) )
      !allocate( radscr03(rad%kk) )
      !allocate( radscr04(rad%kk) )
      !allocate( radscr05(rad%kk) )
      !allocate( radscr06(rad%kk) )
      !allocate( radscr07(rad%kk) )
      !allocate( radscr08(rad%kk) )
      !allocate( radscr09(rad%kk) )
      !allocate( radscr10(rad%kk) )
      !allocate( radscr11(rad%kk) )
      !allocate( radscr12(rad%kk) )
      !allocate( radscr13(rad%kk) )
      !allocate( radscr14(rad%kk) )
      !allocate( radscr15(rad%kk) )
      !allocate( radscr16(rad%kk) )
      !allocate( radscr17(rad%kk) )

      !radscr01 = 0.0_core_rknd
      !radscr02 = 0.0_core_rknd
      !radscr03 = 0.0_core_rknd
      !radscr04 = 0.0_core_rknd
      !radscr05 = 0.0_core_rknd
      !radscr06 = 0.0_core_rknd
      !radscr07 = 0.0_core_rknd
      !radscr08 = 0.0_core_rknd
      !radscr09 = 0.0_core_rknd
      !radscr10 = 0.0_core_rknd
      !radscr11 = 0.0_core_rknd
      !radscr12 = 0.0_core_rknd
      !radscr13 = 0.0_core_rknd
      !radscr14 = 0.0_core_rknd
      !radscr15 = 0.0_core_rknd
      !radscr16 = 0.0_core_rknd
      !radscr17 = 0.0_core_rknd


      fname = trim( fname_rad_zt )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & 
                         1, rad_zt%kk, rad_zt%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         rad_zt%nn, rad_zt%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname,  & 
                          1, rad_zt%kk, rad_zt%z, & 
                          day, month, year, rlat, rlon, & 
                          time_current+stats_tout, stats_tout, & 
                          rad_zt%nn, rad_zt%f )

#else
        stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

      call stats_init_rad_zt( vars_rad_zt, l_error )

      ! Initialize rad_zm (radiation points)

      i = 1
      do while ( ichar(vars_rad_zm(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zm(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zm )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        stop "stats_init:  number of rad_zm statistical variables exceeds limit"
      end if

      rad_zm%nn = ntot
      rad_zm%kk = nnrad_zm

      allocate( rad_zm%z( rad_zm%kk ) )
      rad_zm%z = grad_zm

      allocate( rad_zm%x( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%n( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%l_in_update( 1, 1, rad_zm%kk, rad_zm%nn ) )

      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )

      allocate( rad_zm%f%var( rad_zm%nn ) )
      allocate( rad_zm%f%z( rad_zm%kk ) )

      ! Allocate scratch space

      !allocate( radscr01(rad%kk) )
      !allocate( radscr02(rad%kk) )
      !allocate( radscr03(rad%kk) )
      !allocate( radscr04(rad%kk) )
      !allocate( radscr05(rad%kk) )
      !allocate( radscr06(rad%kk) )
      !allocate( radscr07(rad%kk) )
      !allocate( radscr08(rad%kk) )
      !allocate( radscr09(rad%kk) )
      !allocate( radscr10(rad%kk) )
      !allocate( radscr11(rad%kk) )
      !allocate( radscr12(rad%kk) )
      !allocate( radscr13(rad%kk) )
      !allocate( radscr14(rad%kk) )
      !allocate( radscr15(rad%kk) )
      !allocate( radscr16(rad%kk) )
      !allocate( radscr17(rad%kk) )

      !radscr01 = 0.0_core_rknd
      !radscr02 = 0.0_core_rknd
      !radscr03 = 0.0_core_rknd
      !radscr04 = 0.0_core_rknd
      !radscr05 = 0.0_core_rknd
      !radscr06 = 0.0_core_rknd
      !radscr07 = 0.0_core_rknd
      !radscr08 = 0.0_core_rknd
      !radscr09 = 0.0_core_rknd
      !radscr10 = 0.0_core_rknd
      !radscr11 = 0.0_core_rknd
      !radscr12 = 0.0_core_rknd
      !radscr13 = 0.0_core_rknd
      !radscr14 = 0.0_core_rknd
      !radscr15 = 0.0_core_rknd
      !radscr16 = 0.0_core_rknd
      !radscr17 = 0.0_core_rknd


      fname = trim( fname_rad_zm )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & 
                         1, rad_zm%kk, rad_zm%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         rad_zm%nn, rad_zm%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname,  & 
                          1, rad_zm%kk, rad_zm%z, & 
                          day, month, year, rlat, rlon, & 
                          time_current+stats_tout, stats_tout, & 
                          rad_zm%nn, rad_zm%f )

#else
        stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

      call stats_init_rad_zm( vars_rad_zm, l_error )
    end if ! l_output_rad_files


    ! Initialize sfc (surface point)

    i = 1
    do while ( ichar(vars_sfc(i)(1:1)) /= 0  & 
               .and. len_trim(vars_sfc(i)) /= 0 & 
               .and. i <= nvarmax_sfc )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      stop "stats_init:  number of sfc statistical variables exceeds limit"
    end if

    sfc%nn = ntot
    sfc%kk = 1

    allocate( sfc%z( sfc%kk ) )
    sfc%z = gzm(1)

    allocate( sfc%x( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%n( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%l_in_update( 1, 1, sfc%kk, sfc%nn ) )

    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    allocate( sfc%f%var( sfc%nn ) )
    allocate( sfc%f%z( sfc%kk ) )

    fname = trim( fname_sfc )

    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, sfc%kk, sfc%z, & 
                       day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         sfc%nn, sfc%f )

    else ! Open NetCDF files
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, sfc%kk, sfc%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, sfc%nn, &  ! In
                        sfc%f ) ! InOut

#else
      stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    call stats_init_sfc( vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      stop "Fatal error"
    endif

    return

    ! If namelist was not found in input file, turn off statistics

    100 continue
    write(fstderr,*) 'Error with statsnl, statistics is turned off'
    l_stats       = .false.
    l_stats_samp  = .false.
    l_stats_last  = .false.

    return
  end subroutine stats_init
  !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    ! Description:
    !   Initialize stats to zero
    ! References:
    !   None
    !-----------------------------------------------------------------------
    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: kk, nn

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical, dimension(1,1,kk,nn), intent(out) :: l_in_update

    ! Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0_stat_rknd
      n(:,:,:,:) = 0_stat_nknd
      l_in_update(:,:,:,:) = .false.
    end if

    return
  end subroutine stats_zero

  !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    ! Description:
    !   Compute the average of stats fields
    ! References:
    !   None
    !-----------------------------------------------------------------------
    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    ! External
    intrinsic :: real

    ! Input Variable(s)
    integer, intent(in) :: &
      kk, & ! Number of levels in vertical (i.e. Z) dimension
      nn    ! Number of variables being sampled in x

    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: &
      n ! The variable n is the number of samples per x per kk

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout) :: &
      x ! The variable x is a set of nn variables being averaged over n

    ! ---- Begin Code ----

    ! Compute averages
    where ( n(1,1,1:kk,1:nn) > 0 )
      x(1,1,1:kk,1:nn) = x(1,1,1:kk,1:nn) / real( n(1,1,1:kk,1:nn), kind=stat_rknd )
    end where

    return
  end subroutine stats_avg

  !-----------------------------------------------------------------------
  subroutine stats_begin_timestep( time_elapsed )

    !     Description:
    !       Given the elapsed time, set flags determining specifics such as
    !       if this time set should be sampled or if this is the first or
    !       last time step.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        l_stats,  & ! Variable(s)
        l_stats_samp, & 
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout

    use clubb_precision, only: & 
        time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    real(kind=time_precision), intent(in) ::  & 
      time_elapsed ! Elapsed model time       [s]

    if ( .not. l_stats ) return

    ! Only sample time steps that are multiples of "stats_tsamp"
    ! in a case's "model.in" file to shorten length of run
    if ( mod( time_elapsed, stats_tsamp ) < 1.e-8_time_precision ) then
      l_stats_samp = .true.
    else
      l_stats_samp = .false.
    end if

    ! Indicates the end of the sampling time period. Signals to start writing to the file
    if ( mod( time_elapsed, stats_tout ) < 1.e-8_time_precision ) then
      l_stats_last = .true.
    else
      l_stats_last = .false.
    end if

    return

  end subroutine stats_begin_timestep

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep( )

    ! Description: 
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        LH_zt, &
        LH_sfc, &
        zm, & 
        rad_zt, &
        rad_zm, &
        sfc, & 
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout, &
        l_output_rad_files, & 
        l_grads

    use clubb_precision, only: & 
        time_precision ! Variable(s)

    use output_grads, only: & 
        write_grads ! Procedure(s)

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use parameters_microphys, only: &
      LH_microphys_disabled  ! Constant

    use parameters_microphys, only: &
      LH_microphys_type, & ! Variable(s)
      LH_microphys_calls

#ifdef NETCDF
    use output_netcdf, only: & 
        write_netcdf ! Procedure(s)
#endif

    implicit none

    ! External
    intrinsic :: floor

    ! Local Variables

    integer :: i, k

    logical :: l_error

    ! ---- Begin Code ----

    ! Check if it is time to write to file

    if ( .not. l_stats_last ) return

    ! Initialize
    l_error = .false.

    ! Look for errors by checking the number of sampling points
    ! for each variable in the zt statistics at each vertical level.
    do i = 1, zt%nn
      do k = 1, zt%kk

        if ( zt%n(1,1,k,i) /= 0 .and.  &
             zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zt%f%var(i)%name), ' in zt ',  &
                             'at k = ', k,  &
                             '; zt%n(',k,',',i,') = ', zt%n(1,1,k,i)
          end if ! clubb_at_lest_debug_level 1

        end if ! n /= 0 and n /= stats_tout/stats_tsamp

      end do ! k = 1 .. zt%kk
    end do ! i = 1 .. zt%nn

    ! Look for errors by checking the number of sampling points
    ! for each variable in the zm statistics at each vertical level.
    do i = 1, zm%nn
      do k = 1, zm%kk

        if ( zm%n(1,1,k,i) /= 0 .and.  &
             zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zm%f%var(i)%name), ' in zm ',  &
                             'at k = ', k,  &
                             '; zm%n(',k,',',i,') = ', zm%n(1,1,k,i)
          end if ! clubb_at_least_debug_level 1

        end if ! n /= 0 and n /= stats_tout/stats_tsamp

      end do ! k = 1 .. zm%kk
    end do ! i = 1 .. zm%nn

    if ( LH_microphys_type /= LH_microphys_disabled ) then
      ! Look for errors by checking the number of sampling points
      ! for each variable in the LH_zt statistics at each vertical level.
      do i = 1, LH_zt%nn
        do k = 1, LH_zt%kk

          if ( LH_zt%n(1,1,k,i) /= 0 .and.  &
               LH_zt%n(1,1,k,i) /= floor( stats_tout/stats_tsamp ) .and. &
               LH_zt%n(1,1,k,i) /= LH_microphys_calls * floor( stats_tout/stats_tsamp ) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                trim(LH_zt%f%var(i)%name), ' in LH_zt ',  &
                'at k = ', k,  &
                '; LH_zt%n(',k,',',i,') = ', LH_zt%n(1,1,k,i)
            end if ! clubb_at_lest_debug_level 1

          end if ! n /= 0 and n /= LH_microphys_calls * stats_tout/stats_tsamp

        end do ! k = 1 .. LH_zt%kk
      end do ! i = 1 .. LH_zt%nn

      ! Look for errors by checking the number of sampling points
      ! for each variable in the LH_zt statistics at each vertical level.
      do i = 1, LH_sfc%nn
        do k = 1, LH_sfc%kk

          if ( LH_sfc%n(1,1,k,i) /= 0 .and.  &
               LH_sfc%n(1,1,k,i) /= floor( stats_tout/stats_tsamp ) .and. &
               LH_sfc%n(1,1,k,i) /= LH_microphys_calls * floor( stats_tout/stats_tsamp ) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                trim(LH_sfc%f%var(i)%name), ' in LH_sfc ',  &
                'at k = ', k,  &
                '; LH_sfc%n(',k,',',i,') = ', LH_sfc%n(1,1,k,i)
            end if ! clubb_at_lest_debug_level 1

          end if ! n /= 0 and n /= LH_microphys_calls * stats_tout/stats_tsamp

        end do ! k = 1 .. LH_sfc%kk
      end do ! i = 1 .. LH_sfc%nn
    end if ! LH_microphys_type /= LH_microphys_disabled


    if ( l_output_rad_files ) then
      ! Look for errors by checking the number of sampling points
      ! for each variable in the rad_zt statistics at each vertical level.
      do i = 1, rad_zt%nn
        do k = 1, rad_zt%kk

          if ( rad_zt%n(1,1,k,i) /= 0 .and.  &
               rad_zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zt%f%var(i)%name), ' in rad_zt ',  &
                               'at k = ', k,  &
                               '; rad_zt%n(',k,',',i,') = ', rad_zt%n(1,1,k,i)
            end if ! clubb_at_lest_debug_level 1

          end if ! n /= 0 and n /= stats_tout/stats_tsamp

        end do ! k = 1 .. rad_zt%kk
      end do !  i = 1 .. rad_zt%nn

      ! Look for errors by checking the number of sampling points
      ! for each variable in the rad_zm statistics at each vertical level.
      do i = 1, rad_zm%nn
        do k = 1, rad_zm%kk

          if ( rad_zm%n(1,1,k,i) /= 0 .and.  &
               rad_zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zm%f%var(i)%name), ' in rad_zm ',  &
                               'at k = ', k,  &
                               '; rad_zm%n(',k,',',i,') = ', rad_zm%n(1,1,k,i)
            end if ! clubb_at_lest_debug_level 1

          end if ! n /= 0 and n /= stats_tout/stats_tsamp

        end do ! k = 1 .. rad_zm%kk
      end do !  i = 1 .. rad_zm%nn

    end if ! l_output_rad_files

    ! Look for errors by checking the number of sampling points
    ! for each variable in the sfc statistics at each vertical level.
    do i = 1, sfc%nn
      do k = 1, sfc%kk

        if ( sfc%n(1,1,k,i) /= 0 .and.  &
             sfc%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(sfc%f%var(i)%name), ' in sfc ',  &
                             'at k = ', k,  &
                             '; sfc%n(',k,',',i,') = ', sfc%n(1,1,k,i)
          end if ! clubb_at_lest_debug_level 1

        end if ! n /= 0 and n /= stats_tout/stats_tsamp

      end do ! k = 1 .. sfc%kk
    end do !  i = 1 .. sfc%nn

    ! Stop the run if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      stop 'stats_end_timestep:  error(s) found'
    end if ! l_error

    ! Compute averages
    call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
    call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
    if ( LH_microphys_type /= LH_microphys_disabled ) then
      call stats_avg( LH_zt%kk, LH_zt%nn, LH_zt%x, LH_zt%n )
      call stats_avg( LH_sfc%kk, LH_sfc%nn, LH_sfc%x, LH_sfc%n )
    end if
    if ( l_output_rad_files ) then
      call stats_avg( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n )
      call stats_avg( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n )
    end if
    call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

    ! Write to file
    if ( l_grads ) then
      call write_grads( zt%f  )
      call write_grads( zm%f  )
      if ( LH_microphys_type /= LH_microphys_disabled ) then
        call write_grads( LH_zt%f  )
        call write_grads( LH_sfc%f  )
      end if
      if ( l_output_rad_files ) then
        call write_grads( rad_zt%f  )
        call write_grads( rad_zm%f  )
      end if
      call write_grads( sfc%f  )
    else ! l_netcdf
#ifdef NETCDF
      call write_netcdf( zt%f  )
      call write_netcdf( zm%f  )
      if ( LH_microphys_type /= LH_microphys_disabled ) then
        call write_netcdf( LH_zt%f  )
        call write_netcdf( LH_sfc%f  )
      end if
      if ( l_output_rad_files ) then
        call write_netcdf( rad_zt%f  )
        call write_netcdf( rad_zm%f  )
      end if
      call write_netcdf( sfc%f  )
#else
      stop "This program was not compiled with netCDF support"
#endif /* NETCDF */
    end if ! l_grads

    ! Reset sample fields
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )
    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )
    if ( LH_microphys_type /= LH_microphys_disabled ) then
      call stats_zero( LH_zt%kk, LH_zt%nn, LH_zt%x, LH_zt%n, LH_zt%l_in_update )
      call stats_zero( LH_sfc%kk, LH_sfc%nn, LH_sfc%x, LH_sfc%n, LH_sfc%l_in_update )
    end if
    if ( l_output_rad_files ) then
      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )
      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )
    end if
    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )


    return
  end subroutine stats_end_timestep

  !----------------------------------------------------------------------
  subroutine stats_accumulate & 
                   ( um, vm, upwp, vpwp, up2, vp2, &
                     thlm, rtm, wprtp, wpthlp, &
                     wp2, wp3, rtp2, thlp2, rtpthlp, &
                     p_in_Pa, exner, rho, rho_zm, &
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, &
                     thv_ds_zt, wm_zt, wm_zm, rcm, wprcp, rc_coef, &
                     rcm_zm, rtm_zm, thlm_zm, cloud_frac, ice_supersat_frac, &
                     cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
                     cloud_cover, sigma_sqd_w, pdf_params, &
                     sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
                     wpsclrp, edsclrm, edsclrm_forcing )

    ! Description:
    !   Accumulate those stats variables that are preserved in CLUBB from timestep to
    !   timestep, but not those stats that are not, (e.g. budget terms, longwave and
    !   shortwave components, etc.)
    !
    ! References:
    !   None
    !----------------------------------------------------------------------

    use stats_variables, only: & 
        zt, & ! Variables
        zm, & 
        sfc, & 
        l_stats_samp, & 
        ithlm, & 
        iT_in_K, & 
        ithvm, & 
        irtm, & 
        ircm, & 
        ium, & 
        ivm, & 
        iwm_zt, & 
        iwm_zm, & 
        iug, & 
        ivg, & 
        icloud_frac, &
        iice_supersat_frac, & 
        ircm_in_layer, &
        icloud_cover

    use stats_variables, only: &
        ip_in_Pa, & 
        iexner, & 
        irho_ds_zt, &
        ithv_ds_zt, &
        iLscale, & 
        iwp3, & 
        iwp3_zm, & 
        iwpthlp2, & 
        iwp2thlp,  & 
        iwprtp2, & 
        iwp2rtp, & 
        iLscale_up, & 
        iLscale_down, & 
        itau_zt, & 
        iKh_zt

    use stats_variables, only: & 
        iwp2thvp, &  ! Variable(s)
        iwp2rcp, & 
        iwprtpthlp, & 
        isigma_sqd_w_zt, & 
        irho, & 
        irsat, & 
        irsati

    use stats_variables, only: & 
        imixt_frac, &  ! Variable(s)
        iw1, & 
        iw2, & 
        ivarnce_w1, & 
        ivarnce_w2, & 
        ithl1, & 
        ithl2, & 
        ivarnce_thl1, & 
        ivarnce_thl2, & 
        irt1, & 
        irt2, & 
        ivarnce_rt1, & 
        ivarnce_rt2, & 
        irc1, & 
        irc2, & 
        irsl1, & 
        irsl2, & 
        icloud_frac1, & 
        icloud_frac2

    use stats_variables, only: & 
        is1, & 
        is2, & 
        istdev_s1, & 
        istdev_s2, &
        istdev_t1, &
        istdev_t2, &
        icovar_st_1, &
        icovar_st_2, &
        icorr_st_1, &
        icorr_st_2, &
        icrt1, &
        icrt2, &
        icthl1, &
        icthl2, &
        irrtthl, &
        is_mellor

    use stats_variables, only: & 
        iwp2_zt, &  ! Variable(s)
        ithlp2_zt, & 
        iwpthlp_zt, & 
        iwprtp_zt, & 
        irtp2_zt, & 
        irtpthlp_zt, &
        iup2_zt, &
        ivp2_zt, &
        iupwp_zt, &
        ivpwp_zt, & 
        iwp2, & 
        irtp2, & 
        ithlp2, & 
        irtpthlp, & 
        iwprtp,  & 
        iwpthlp, & 
        iwp4,  & 
        iwpthvp, & 
        irtpthvp

    use stats_variables, only: & 
        ithlpthvp, & 
        itau_zm, & 
        iKh_zm, & 
        iwprcp, & 
        irc_coef, &
        ithlprcp, & 
        irtprcp, & 
        ircp2, & 
        iupwp, & 
        ivpwp, & 
        iup2, & 
        ivp2, & 
        irho_zm, & 
        isigma_sqd_w, &
        irho_ds_zm, &
        ithv_ds_zm, &
        iem

    use stats_variables, only: & 
        ishear, &  ! Variable(s)
        iFrad, & 
        icc, & 
        iz_cloud_base, & 
        ilwp, &
        ivwp, &
        ithlm_vert_avg, &
        irtm_vert_avg, &
        ium_vert_avg, &
        ivm_vert_avg, &
        iwp2_vert_avg, &
        iup2_vert_avg, &
        ivp2_vert_avg, &
        irtp2_vert_avg, &
        ithlp2_vert_avg

    use stats_variables, only: & 
        isclrm, &  ! Variable(s)
        isclrm_f, & 
        iedsclrm, & 
        iedsclrm_f, & 
        isclrprtp, & 
        isclrp2, & 
        isclrpthvp, & 
        isclrpthlp, & 
        isclrprcp, & 
        iwpsclrp, & 
        iwp2sclrp, & 
        iwpsclrp2, & 
        iwpsclrprtp, & 
        iwpsclrpthlp, & 
        iwpedsclrp

    use stats_variables, only: &
      icloud_frac_zm, &
      iice_supersat_frac_zm, &
      ircm_zm, &
      irtm_zm, &
      ithlm_zm

    use stats_variables, only: &
      iwp3_on_wp2, &
      iwp3_on_wp2_zt, &
      iSkw_velocity

    use stats_variables, only: &
      ia3_coef, & ! Variables
      ia3_coef_zt

    use grid_class, only: & 
        gr ! Variable

    use grid_class, only: & 
      zt2zm ! Procedure(s)

    use variables_diagnostic_module, only: & 
        thvm, & ! Variable(s)
        ug, & 
        vg, & 
        Lscale, & 
        wpthlp2, & 
        wp2thlp, & 
        wprtp2, & 
        wp2rtp, & 
        Lscale_up, & 
        Lscale_down, & 
        tau_zt, &
        Kh_zt, & 
        wp2thvp, & 
        wp2rcp, & 
        wprtpthlp, & 
        sigma_sqd_w_zt, & 
        rsat

    use variables_diagnostic_module, only: & 
        wp2_zt, &  ! Variable(s)
        thlp2_zt, & 
        wpthlp_zt, & 
        wprtp_zt, & 
        rtp2_zt, & 
        rtpthlp_zt, &
        up2_zt, &
        vp2_zt, &
        upwp_zt, &
        vpwp_zt, & 
        wp4, & 
        rtpthvp, & 
        thlpthvp, & 
        wpthvp, &
        tau_zm, &
        Kh_zm, & 
        thlprcp, & 
        rtprcp, & 
        rcp2, & 
        em, & 
        Frad, & 
        sclrpthvp, & 
        sclrprcp, & 
        wp2sclrp, & 
        wpsclrp2, & 
        wpsclrprtp, & 
        wpsclrpthlp, & 
        wpedsclrp

    use variables_diagnostic_module, only: & 
      a3_coef, & ! Variable(s)
      a3_coef_zt, &
      wp3_zm, &
      wp3_on_wp2, &
      wp3_on_wp2_zt, &
      Skw_velocity

    use pdf_parameter_module, only: & 
      pdf_parameter ! Type

    use T_in_K_module, only: & 
        thlm2T_in_K ! Procedure

    use constants_clubb, only: & 
        rc_tol, & ! Constant(s)
        w_tol_sqd

    use parameters_model, only: & 
        sclr_dim,  &        ! Variable(s)
        edsclr_dim

    use stats_type, only: & 
        stat_update_var,  & ! Procedure(s)
        stat_update_var_pt

    use fill_holes, only: &
        vertical_avg, &     ! Procedure(s)
        vertical_integral

    use interpolation, only: & 
        lin_int             ! Procedure

    use saturation, only: &
      sat_mixrat_ice ! Procedure

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      um,      & ! u wind                        [m/s]
      vm,      & ! v wind                        [m/s]
      upwp,    & ! vertical u momentum flux      [m^2/s^2]
      vpwp,    & ! vertical v momentum flux      [m^2/s^2]
      up2,     & ! u'^2                          [m^2/s^2]
      vp2,     & ! v'^2                          [m^2/s^2]
      thlm,    & ! liquid potential temperature  [K]
      rtm,     & ! total water mixing ratio      [kg/kg]
      wprtp,   & ! w'rt'                         [(kg/kg) m/s]
      wpthlp,  & ! w'thl'                        [m K /s]
      wp2,     & ! w'^2                          [m^2/s^2]
      wp3,     & ! w'^3                          [m^3/s^3]
      rtp2,    & ! rt'^2                         [(kg/kg)^2]
      thlp2,   & ! thl'^2                        [K^2]
      rtpthlp    ! rt'thl'                       [kg/kg K]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
      exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
      rho,          & ! Density                                  [kg/m^3]
      rho_zm,       & ! Density                                  [kg/m^3]
      rho_ds_zm,    & ! Dry, static density (momentum levels)    [kg/m^3]
      rho_ds_zt,    & ! Dry, static density (thermo. levs.)      [kg/m^3]
      thv_ds_zm,    & ! Dry, base-state theta_v (momentum levs.) [K]
      thv_ds_zt,    & ! Dry, base-state theta_v (thermo. levs.)  [K]
      wm_zt,        & ! w on thermodynamic levels                [m/s]
      wm_zm           ! w on momentum levels                     [m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      rcm_zm,               & ! Total water mixing ratio                 [kg/kg]
      rtm_zm,               & ! Total water mixing ratio                 [kg/kg]
      thlm_zm,              & ! Liquid potential temperature             [K]
      rcm,                  & ! Cloud water mixing ratio                 [kg/kg]
      wprcp,                & ! w'rc'                                    [(kg/kg) m/s]
      rc_coef,              & ! Coefficient of X' R_l' in Eq. (34)       [-]
      cloud_frac,           & ! Cloud fraction                           [-]
      ice_supersat_frac,    & ! Ice cloud fracion                        [-]
      cloud_frac_zm,        & ! Cloud fraction on zm levels              [-]
      ice_supersat_frac_zm, & ! Ice cloud fraction on zm levels          [-]
      rcm_in_layer,         & ! Cloud water mixing ratio in cloud layer  [kg/kg]
      cloud_cover             ! Cloud cover                              [-]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      sigma_sqd_w    ! PDF width parameter (momentum levels)    [-]

    type(pdf_parameter), dimension(gr%nz), intent(in) :: & 
      pdf_params ! PDF parameters [units vary]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: & 
      sclrm,           & ! High-order passive scalar            [units vary]
      sclrp2,          & ! High-order passive scalar variance   [units^2]
      sclrprtp,        & ! High-order passive scalar covariance [units kg/kg]
      sclrpthlp,       & ! High-order passive scalar covariance [units K]
      sclrm_forcing,   & ! Large-scale forcing of scalar        [units/s]
      wpsclrp            ! w'sclr'                              [units m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: & 
      edsclrm,         & ! Eddy-diff passive scalar      [units vary] 
      edsclrm_forcing    ! Large-scale forcing of edscalar  [units vary]

    ! Local Variables

    integer :: i, k

    real( kind = core_rknd ), dimension(gr%nz) :: &
      T_in_K, &  ! Absolute temperature         [K]
      rsati,  &  ! Saturation w.r.t ice         [kg/kg]
      shear,  &  ! Wind shear production term   [m^2/s^3]
      s_mellor   ! Mellor's 's'                 [kg/kg]

    real( kind = core_rknd ) :: xtmp

    ! ---- Begin Code ----

    ! Sample fields

    if ( l_stats_samp ) then

      ! zt variables


      if ( iT_in_K > 0 .or. irsati > 0 ) then
        T_in_K = thlm2T_in_K( thlm, exner, rcm )
      else
        T_in_K = -999._core_rknd
      end if

      call stat_update_var( iT_in_K, T_in_K, zt )

      call stat_update_var( ithlm, thlm, zt )
      call stat_update_var( ithvm, thvm, zt )
      call stat_update_var( irtm, rtm, zt )
      call stat_update_var( ircm, rcm, zt )
      call stat_update_var( ium, um, zt )
      call stat_update_var( ivm, vm, zt )
      call stat_update_var( iwm_zt, wm_zt, zt )
      call stat_update_var( iwm_zm, wm_zm, zm )
      call stat_update_var( iug, ug, zt )
      call stat_update_var( ivg, vg, zt )
      call stat_update_var( icloud_frac, cloud_frac, zt )
      call stat_update_var( iice_supersat_frac, ice_supersat_frac, zt)
      call stat_update_var( ircm_in_layer, rcm_in_layer, zt )
      call stat_update_var( icloud_cover, cloud_cover, zt )
      call stat_update_var( ip_in_Pa, p_in_Pa, zt )
      call stat_update_var( iexner, exner, zt )
      call stat_update_var( irho_ds_zt, rho_ds_zt, zt )
      call stat_update_var( ithv_ds_zt, thv_ds_zt, zt )
      call stat_update_var( iLscale, Lscale, zt )
      call stat_update_var( iwp3, wp3, zt )
      call stat_update_var( iwpthlp2, wpthlp2, zt )
      call stat_update_var( iwp2thlp, wp2thlp, zt )
      call stat_update_var( iwprtp2, wprtp2, zt )
      call stat_update_var( iwp2rtp, wp2rtp, zt )
      call stat_update_var( iLscale_up, Lscale_up, zt )
      call stat_update_var( iLscale_down, Lscale_down, zt )
      call stat_update_var( itau_zt, tau_zt, zt )
      call stat_update_var( iKh_zt, Kh_zt, zt )
      call stat_update_var( iwp2thvp, wp2thvp, zt )
      call stat_update_var( iwp2rcp, wp2rcp, zt )
      call stat_update_var( iwprtpthlp, wprtpthlp, zt )
      call stat_update_var( isigma_sqd_w_zt, sigma_sqd_w_zt, zt )
      call stat_update_var( irho, rho, zt )
      call stat_update_var( irsat, rsat, zt )
      if ( irsati > 0 ) then
        rsati = sat_mixrat_ice( p_in_Pa, T_in_K )
        call stat_update_var( irsati, rsati, zt )
      end if

      call stat_update_var( imixt_frac, pdf_params%mixt_frac, zt )
      call stat_update_var( iw1, pdf_params%w1, zt )
      call stat_update_var( iw2, pdf_params%w2, zt )
      call stat_update_var( ivarnce_w1, pdf_params%varnce_w1, zt )
      call stat_update_var( ivarnce_w2, pdf_params%varnce_w2, zt )
      call stat_update_var( ithl1, pdf_params%thl1, zt )
      call stat_update_var( ithl2, pdf_params%thl2, zt )
      call stat_update_var( ivarnce_thl1, pdf_params%varnce_thl1, zt )
      call stat_update_var( ivarnce_thl2, pdf_params%varnce_thl2, zt )
      call stat_update_var( irt1, pdf_params%rt1, zt )
      call stat_update_var( irt2, pdf_params%rt2, zt )
      call stat_update_var( ivarnce_rt1, pdf_params%varnce_rt1, zt )
      call stat_update_var( ivarnce_rt2, pdf_params%varnce_rt2, zt )
      call stat_update_var( irc1, pdf_params%rc1, zt )
      call stat_update_var( irc2, pdf_params%rc2, zt )
      call stat_update_var( irsl1, pdf_params%rsl1, zt )
      call stat_update_var( irsl2, pdf_params%rsl2, zt )
      call stat_update_var( icloud_frac1, pdf_params%cloud_frac1, zt )
      call stat_update_var( icloud_frac2, pdf_params%cloud_frac2, zt )
      call stat_update_var( is1, pdf_params%s1, zt )
      call stat_update_var( is2, pdf_params%s2, zt )
      call stat_update_var( istdev_s1, pdf_params%stdev_s1, zt )
      call stat_update_var( istdev_s2, pdf_params%stdev_s2, zt )
      call stat_update_var( istdev_t1, pdf_params%stdev_t1, zt )
      call stat_update_var( istdev_t2, pdf_params%stdev_t2, zt )
      call stat_update_var( icovar_st_1, pdf_params%covar_st_1, zt )
      call stat_update_var( icovar_st_2, pdf_params%covar_st_2, zt )
      call stat_update_var( icorr_st_1, pdf_params%corr_st_1, zt )
      call stat_update_var( icorr_st_2, pdf_params%corr_st_2, zt )
      call stat_update_var( irrtthl, pdf_params%rrtthl, zt )
      call stat_update_var( icrt1, pdf_params%crt1, zt )
      call stat_update_var( icrt2, pdf_params%crt2, zt )
      call stat_update_var( icthl1, pdf_params%cthl1, zt )
      call stat_update_var( icthl2, pdf_params%cthl2, zt )
      call stat_update_var( iwp2_zt, wp2_zt, zt )
      call stat_update_var( ithlp2_zt, thlp2_zt, zt )
      call stat_update_var( iwpthlp_zt, wpthlp_zt, zt )
      call stat_update_var( iwprtp_zt, wprtp_zt, zt )
      call stat_update_var( irtp2_zt, rtp2_zt, zt )
      call stat_update_var( irtpthlp_zt, rtpthlp_zt, zt )
      call stat_update_var( iup2_zt, up2_zt, zt )
      call stat_update_var( ivp2_zt, vp2_zt, zt )
      call stat_update_var( iupwp_zt, upwp_zt, zt )
      call stat_update_var( ivpwp_zt, vpwp_zt, zt )
      call stat_update_var( ia3_coef_zt, a3_coef_zt, zt )
      call stat_update_var( iwp3_on_wp2_zt, wp3_on_wp2_zt, zt )

      if ( is_mellor > 0 ) then
        ! Determine 's' from Mellor (1977) (extended liquid water)
        s_mellor(:) = pdf_params%mixt_frac * pdf_params%s1 &
                    + (1.0_core_rknd-pdf_params%mixt_frac) * pdf_params%s2
        call stat_update_var( is_mellor, s_mellor, zt )
      end if

      if ( sclr_dim > 0 ) then
        do i=1, sclr_dim
          call stat_update_var( isclrm(i), sclrm(:,i), zt )
          call stat_update_var( isclrm_f(i), sclrm_forcing(:,i),  zt )
        end do
      end if

      if ( edsclr_dim > 0 ) then
        do i=1, edsclr_dim
          call stat_update_var( iedsclrm(i), edsclrm(:,i), zt )
          call stat_update_var( iedsclrm_f(i), edsclrm_forcing(:,i), zt )
        end do
      end if

      ! zm variables

      call stat_update_var( iwp2, wp2, zm )
      call stat_update_var( iwp3_zm, wp3_zm, zm )
      call stat_update_var( irtp2, rtp2, zm )
      call stat_update_var( ithlp2, thlp2, zm )
      call stat_update_var( irtpthlp, rtpthlp, zm )
      call stat_update_var( iwprtp, wprtp, zm )
      call stat_update_var( iwpthlp, wpthlp, zm )
      call stat_update_var( iwp4, wp4, zm )
      call stat_update_var( iwpthvp, wpthvp, zm )
      call stat_update_var( irtpthvp, rtpthvp, zm )
      call stat_update_var( ithlpthvp, thlpthvp, zm )
      call stat_update_var( itau_zm, tau_zm, zm )
      call stat_update_var( iKh_zm, Kh_zm, zm )
      call stat_update_var( iwprcp, wprcp, zm )
      call stat_update_var( irc_coef, rc_coef, zm )
      call stat_update_var( ithlprcp, thlprcp, zm )
      call stat_update_var( irtprcp, rtprcp, zm )
      call stat_update_var( ircp2, rcp2, zm )
      call stat_update_var( iupwp, upwp, zm )
      call stat_update_var( ivpwp, vpwp, zm )
      call stat_update_var( ivp2, vp2, zm )
      call stat_update_var( iup2, up2, zm )
      call stat_update_var( irho_zm, rho_zm, zm )
      call stat_update_var( isigma_sqd_w, sigma_sqd_w, zm )
      call stat_update_var( irho_ds_zm, rho_ds_zm, zm )
      call stat_update_var( ithv_ds_zm, thv_ds_zm, zm )
      call stat_update_var( iem, em, zm )
      call stat_update_var( iFrad, Frad, zm )

      call stat_update_var( iSkw_velocity, Skw_velocity, zm )
      call stat_update_var( ia3_coef, a3_coef, zm )
      call stat_update_var( iwp3_on_wp2, wp3_on_wp2, zm )

      call stat_update_var( icloud_frac_zm, cloud_frac_zm, zm )
      call stat_update_var( iice_supersat_frac_zm, ice_supersat_frac_zm, zm )
      call stat_update_var( ircm_zm, rcm_zm, zm )
      call stat_update_var( irtm_zm, rtm_zm, zm )
      call stat_update_var( ithlm_zm, thlm_zm, zm )

      if ( sclr_dim > 0 ) then
        do i=1, sclr_dim
          call stat_update_var( isclrp2(i), sclrp2(:,i), zm )
          call stat_update_var( isclrprtp(i), sclrprtp(:,i), zm )
          call stat_update_var( isclrpthvp(i), sclrpthvp(:,i), zm )
          call stat_update_var( isclrpthlp(i), sclrpthlp(:,i), zm )
          call stat_update_var( isclrprcp(i), sclrprcp(:,i), zm )
          call stat_update_var( iwpsclrp(i), wpsclrp(:,i), zm )
          call stat_update_var( iwp2sclrp(i), wp2sclrp(:,i), zm )
          call stat_update_var( iwpsclrp2(i), wpsclrp2(:,i), zm )
          call stat_update_var( iwpsclrprtp(i), wpsclrprtp(:,i), zm )
          call stat_update_var( iwpsclrpthlp(i), wpsclrpthlp(:,i), zm )
        end do
      end if
      if ( edsclr_dim > 0 ) then
        do i=1, edsclr_dim
          call stat_update_var( iwpedsclrp(i), wpedsclrp(:,i), zm )
        end do
      end if

      ! Calculate shear production
      if ( ishear > 0 ) then
        do k = 1, gr%nz-1, 1
          shear(k) = - upwp(k) * ( um(k+1) - um(k) ) * gr%invrs_dzm(k)  &
                     - vpwp(k) * ( vm(k+1) - vm(k) ) * gr%invrs_dzm(k)
        enddo
        shear(gr%nz) = 0.0_core_rknd
      end if
      call stat_update_var( ishear, shear, zm )

      ! sfc variables

      ! Cloud cover
      call stat_update_var_pt( icc, 1, maxval( cloud_frac(1:gr%nz) ), sfc )

      ! Cloud base
      if ( iz_cloud_base > 0 ) then

        k = 1
        do while ( rcm(k) < rc_tol .and. k < gr%nz )
          k = k + 1
        enddo

        if ( k > 1 .and. k < gr%nz) then

          ! Use linear interpolation to find the exact height of the
          ! rc_tol kg/kg level.  Brian.
          call stat_update_var_pt( iz_cloud_base, 1, lin_int( rc_tol, rcm(k),  &
                                   rcm(k-1), gr%zt(k), gr%zt(k-1) ), sfc )

        else

          ! Set the cloud base output to -10m, if it's clear. 
          call stat_update_var_pt( iz_cloud_base, 1, -10.0_core_rknd , sfc ) ! Known magic number
 
        end if ! k > 1 and k < gr%nz

      end if ! iz_cloud_base > 0

      ! Liquid Water Path
      if ( ilwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 rcm(2:gr%nz), gr%invrs_dzt(2:gr%nz) )

        call stat_update_var_pt( ilwp, 1, xtmp, sfc )

      end if

      ! Vapor Water Path (Preciptable Water)
      if ( ivwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 ( rtm(2:gr%nz) - rcm(2:gr%nz) ), gr%invrs_dzt(2:gr%nz) )

        call stat_update_var_pt( ivwp, 1, xtmp, sfc )

      end if


      ! Vertical average of thermodynamic level variables.

      ! Find the vertical average of thermodynamic level variables, averaged from
      ! level 2 (the first thermodynamic level above model surface) through
      ! level gr%nz (the top of the model).  Use the vertical averaging function
      ! found in fill_holes.F90.

      ! Vertical average of thlm.
      call stat_update_var_pt( ithlm_vert_avg, 1,  &
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), &
                         thlm(2:gr%nz), gr%invrs_dzt(2:gr%nz) ), &
                               sfc )

      ! Vertical average of rtm.
      call stat_update_var_pt( irtm_vert_avg, 1,  &
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), &
                         rtm(2:gr%nz), gr%invrs_dzt(2:gr%nz) ), &
                               sfc )

      ! Vertical average of um.
      call stat_update_var_pt( ium_vert_avg, 1,  &
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), &
                         um(2:gr%nz), gr%invrs_dzt(2:gr%nz) ), &
                               sfc )

      ! Vertical average of vm.
      call stat_update_var_pt( ivm_vert_avg, 1,  &
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), &
                         vm(2:gr%nz), gr%invrs_dzt(2:gr%nz) ), &
                               sfc )

      ! Vertical average of momentum level variables.

      ! Find the vertical average of momentum level variables, averaged over the
      ! entire vertical profile (level 1 through level gr%nz).  Use the vertical
      ! averaging function found in fill_holes.F90.

      ! Vertical average of wp2.
      call stat_update_var_pt( iwp2_vert_avg, 1,  &
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), &
                         wp2(1:gr%nz), gr%invrs_dzm(1:gr%nz) ), &
                               sfc )

      ! Vertical average of up2.
      call stat_update_var_pt( iup2_vert_avg, 1,  &
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), &
                         up2(1:gr%nz), gr%invrs_dzm(1:gr%nz) ), &
                               sfc )

      ! Vertical average of vp2.
      call stat_update_var_pt( ivp2_vert_avg, 1,  &
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), &
                         vp2(1:gr%nz), gr%invrs_dzm(1:gr%nz) ), &
                               sfc )

      ! Vertical average of rtp2.
      call stat_update_var_pt( irtp2_vert_avg, 1,  &
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), &
                         rtp2(1:gr%nz), gr%invrs_dzm(1:gr%nz) ), &
                               sfc )

      ! Vertical average of thlp2.
      call stat_update_var_pt( ithlp2_vert_avg, 1,  &
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), &
                         thlp2(1:gr%nz), gr%invrs_dzm(1:gr%nz) ), &
                               sfc )


    end if ! l_stats_samp


    return
  end subroutine stats_accumulate
!------------------------------------------------------------------------------
  subroutine stats_accumulate_hydromet( hydromet, rho_ds_zt )
! Description:
!   Compute stats related the hydrometeors

! References:
!   None
!------------------------------------------------------------------------------
    use parameters_model, only: &
      hydromet_dim ! Variable(s)

    use grid_class, only: &
      gr ! Variable(s)

    use array_index, only:  & 
      iirrainm, iirsnowm, iiricem, iirgraupelm, & ! Variable(s)
      iiNrm, iiNsnowm, iiNim, iiNgraupelm

    use stats_variables, only: &
      sfc, & ! Variable(s)
      irrainm, & 
      irsnowm, & 
      iricem, & 
      irgraupelm, & 
      iNim, & 
      iNrm, & 
      iNsnowm, &
      iNgraupelm, &
      iswp, &
      irwp, &
      iiwp

    use fill_holes, only: &
      vertical_integral ! Procedure(s)

    use stats_type, only: & 
      stat_update_var, & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: &
      zt, & ! Variables
      l_stats_samp

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    ! Local Variables
    real(kind=core_rknd) :: xtmp

    ! ---- Begin Code ----

    if ( l_stats_samp ) then

      if ( iirrainm > 0 ) then
        call stat_update_var( irrainm, hydromet(:,iirrainm), zt )
      end if

      if ( iirsnowm > 0 ) then
        call stat_update_var( irsnowm, hydromet(:,iirsnowm), zt )
      end if

      if ( iiricem > 0 ) then
        call stat_update_var( iricem, hydromet(:,iiricem), zt )
      end if

      if ( iirgraupelm > 0 ) then
        call stat_update_var( irgraupelm,  & 
                              hydromet(:,iirgraupelm), zt )
      end if

      if ( iiNim > 0 ) then
        call stat_update_var( iNim, hydromet(:,iiNim), zt )
      end if

      if ( iiNrm > 0 ) then
        call stat_update_var( iNrm, hydromet(:,iiNrm), zt )
      end if

      if ( iiNsnowm > 0 ) then
        call stat_update_var( iNsnowm, hydromet(:,iiNsnowm), zt )
      end if

      if ( iiNgraupelm > 0 ) then
        call stat_update_var( iNgraupelm, hydromet(:,iiNgraupelm), zt )
      end if

      ! Snow Water Path
      if ( iswp > 0 .and. iirsnowm > 0 ) then

        ! Calculate snow water path
        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iirsnowm), gr%invrs_dzt(2:gr%nz) )

        call stat_update_var_pt( iswp, 1, xtmp, sfc )

      end if ! iswp > 0 .and. iirsnowm > 0

      ! Ice Water Path
      if ( iiwp > 0 .and. iiricem > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iiricem), gr%invrs_dzt(2:gr%nz) )

        call stat_update_var_pt( iiwp, 1, xtmp, sfc )

      end if

      ! Rain Water Path
      if ( irwp > 0 .and. iirrainm > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iirrainm), gr%invrs_dzt(2:gr%nz) )

        call stat_update_var_pt( irwp, 1, xtmp, sfc )

      end if ! irwp > 0 .and. irrainm > 0
    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_hydromet
!------------------------------------------------------------------------------
  subroutine stats_accumulate_LH_tend( LH_hydromet_mc, LH_thlm_mc, LH_rvm_mc, LH_rcm_mc )
! Description:
!   Compute stats for the tendency of latin hypercube sample points.

! References:
!   None
!------------------------------------------------------------------------------
    use parameters_model, only: &
      hydromet_dim ! Variable(s)

    use grid_class, only: &
      gr ! Variable(s)

    use array_index, only:  & 
      iirrainm, iirsnowm, iiricem, iirgraupelm, & ! Variable(s)
      iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm

    use stats_variables, only: &
      iLH_rrainm_mc, & ! Variable(s)
      iLH_rsnowm_mc, & 
      iLH_ricem_mc, & 
      iLH_rgraupelm_mc, & 
      iLH_Ncm_mc, &
      iLH_Nim_mc, & 
      iLH_Nrm_mc, & 
      iLH_Nsnowm_mc, &
      iLH_Ngraupelm_mc, &
      iLH_rcm_mc, &
      iLH_rvm_mc, &
      iLH_thlm_mc

    use stats_variables, only: &
      iAKstd, & ! Variable(s)
      iAKstd_cld, &
      iAKm_rcm, &
      iAKm_rcc, &
      iAKm, & 
      iLH_AKm, &
      iLH_rcm_avg

     use variables_diagnostic_module, only: &
      AKm, & ! Variable(s)
      lh_AKm, &
      AKstd, & 
      lh_rcm_avg, &
      AKstd_cld, &
      AKm_rcm, &
      AKm_rcc

    use stats_type, only: & 
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
      LH_zt, & ! Variables
      l_stats_samp

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      LH_hydromet_mc ! Tendency of hydrometeors except for rvm/rcm  [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      LH_thlm_mc, & ! Tendency of liquid potential temperature [kg/kg/s]
      LH_rcm_mc,  & ! Tendency of cloud water                  [kg/kg/s]
      LH_rvm_mc     ! Tendency of vapor                        [kg/kg/s]

    if ( l_stats_samp ) then

      call stat_update_var( iLH_thlm_mc, LH_thlm_mc, LH_zt )
      call stat_update_var( iLH_rcm_mc, LH_rcm_mc, LH_zt )
      call stat_update_var( iLH_rvm_mc, LH_rvm_mc, LH_zt )

      if ( iiNcm > 0 ) then
        call stat_update_var( iLH_Ncm_mc, LH_hydromet_mc(:,iiNcm), LH_zt )
      end if

      if ( iirrainm > 0 ) then
        call stat_update_var( iLH_rrainm_mc, LH_hydromet_mc(:,iirrainm), LH_zt )
      end if

      if ( iirsnowm > 0 ) then
        call stat_update_var( iLH_rsnowm_mc, LH_hydromet_mc(:,iirsnowm), LH_zt )
      end if

      if ( iiricem > 0 ) then
        call stat_update_var( iLH_ricem_mc, LH_hydromet_mc(:,iiricem), LH_zt )
      end if

      if ( iirgraupelm > 0 ) then
        call stat_update_var( iLH_rgraupelm_mc, LH_hydromet_mc(:,iirgraupelm), LH_zt )
      end if

      if ( iiNim > 0 ) then
        call stat_update_var( iLH_Nim_mc, LH_hydromet_mc(:,iiNim), LH_zt )
      end if

      if ( iiNrm > 0 ) then
        call stat_update_var( iLH_Nrm_mc, LH_hydromet_mc(:,iiNrm), LH_zt )
      end if

      if ( iiNsnowm > 0 ) then
        call stat_update_var( iLH_Nsnowm_mc, LH_hydromet_mc(:,iiNsnowm), LH_zt )
      end if

      if ( iiNgraupelm > 0 ) then
        call stat_update_var( iLH_Ngraupelm_mc, LH_hydromet_mc(:,iiNgraupelm), LH_zt )
      end if

      call stat_update_var( iAKm, AKm, LH_zt )
      call stat_update_var( iLH_AKm, lh_AKm, LH_zt)
      call stat_update_var( iLH_rcm_avg, lh_rcm_avg, LH_zt )
      call stat_update_var( iAKstd, AKstd, LH_zt )
      call stat_update_var( iAKstd_cld, AKstd_cld, LH_zt )

      call stat_update_var( iAKm_rcm, AKm_rcm, LH_zt)
      call stat_update_var( iAKm_rcc, AKm_rcc, LH_zt )

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_LH_tend

  !-----------------------------------------------------------------------
  subroutine stats_finalize( )

    !     Description:
    !     Close NetCDF files and deallocate scratch space and
    !     stats file structures.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        LH_zt, &
        LH_sfc, &
        zm, &
        rad_zt, &
        rad_zm, & 
        sfc, & 
        l_netcdf, & 
        l_stats, &
        l_output_rad_files

    use stats_variables, only: & 
        ztscr01, &  ! Variable(s)
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        ztscr10, & 
        ztscr11, & 
        ztscr12, & 
        ztscr13, & 
        ztscr14, & 
        ztscr15, & 
        ztscr16, & 
        ztscr17, & 
        ztscr18, & 
        ztscr19, & 
        ztscr20, & 
        ztscr21

    use stats_variables, only: & 
        zmscr01, &  ! Variable(s)
        zmscr02, & 
        zmscr03, & 
        zmscr04, & 
        zmscr05, & 
        zmscr06, & 
        zmscr07, & 
        zmscr08, & 
        zmscr09, & 
        zmscr10, & 
        zmscr11, & 
        zmscr12, & 
        zmscr13, & 
        zmscr14, & 
        zmscr15, & 
        zmscr16, & 
        zmscr17

    !use stats_variables, only: &
    !    radscr01, &  ! Variable(s)
    !    radscr02, &
    !    radscr03, &
    !    radscr04, &
    !    radscr05, &
    !    radscr06, &
    !    radscr07, &
    !    radscr08, &
    !    radscr09, &
    !    radscr10, &
    !    radscr11, &
    !    radscr12, &
    !    radscr13, &
    !    radscr14, &
    !    radscr15, &
    !    radscr16, &
    !    radscr17

    use stats_variables, only: & 
      isclrm, & 
      isclrm_f, & 
      iedsclrm, & 
      iedsclrm_f, & 
      isclrprtp, & 
      isclrp2, & 
      isclrpthvp, & 
      isclrpthlp, & 
      isclrprcp, & 
      iwpsclrp, & 
      iwp2sclrp, & 
      iwpsclrp2, & 
      iwpsclrprtp, & 
      iwpsclrpthlp, & 
      iwpedsclrp

    use parameters_microphys, only: &
      LH_microphys_disabled ! Constant(s)

    use parameters_microphys, only: &
      LH_microphys_type ! Variable(s)

#ifdef NETCDF
    use output_netcdf, only:  & 
        close_netcdf ! Procedure
#endif

    implicit none

    if ( l_stats .and. l_netcdf ) then
#ifdef NETCDF
      call close_netcdf( zt%f )
      call close_netcdf( LH_zt%f )
      call close_netcdf( LH_sfc%f )
      call close_netcdf( zm%f )
      call close_netcdf( rad_zt%f )
      call close_netcdf( rad_zm%f )
      call close_netcdf( sfc%f )
#else
      stop "This program was not compiled with netCDF support"
#endif
    end if

    if ( l_stats ) then
      ! De-allocate all zt variables
      deallocate( zt%z )

      deallocate( zt%x )

      deallocate( zt%n )
      deallocate( zt%l_in_update )


      deallocate( zt%f%var )
      deallocate( zt%f%z )
      deallocate( zt%f%rlat )
      deallocate( zt%f%rlon )

      deallocate ( ztscr01 )
      deallocate ( ztscr02 )
      deallocate ( ztscr03 )
      deallocate ( ztscr04 )
      deallocate ( ztscr05 )
      deallocate ( ztscr06 )
      deallocate ( ztscr07 )
      deallocate ( ztscr08 )
      deallocate ( ztscr09 )
      deallocate ( ztscr10 )
      deallocate ( ztscr11 )
      deallocate ( ztscr12 )
      deallocate ( ztscr13 )
      deallocate ( ztscr14 )
      deallocate ( ztscr15 )
      deallocate ( ztscr16 )
      deallocate ( ztscr17 )
      deallocate ( ztscr18 )
      deallocate ( ztscr19 )
      deallocate ( ztscr20 )
      deallocate ( ztscr21 )

      if ( LH_microphys_type /= LH_microphys_disabled ) then
        ! De-allocate all LH_zt variables
        deallocate( LH_zt%z )

        deallocate( LH_zt%x )

        deallocate( LH_zt%n )
        deallocate( LH_zt%l_in_update )


        deallocate( LH_zt%f%var )
        deallocate( LH_zt%f%z )
        deallocate( LH_zt%f%rlat )
        deallocate( LH_zt%f%rlon )

        ! De-allocate all LH_sfc variables
        deallocate( LH_sfc%z )

        deallocate( LH_sfc%x )

        deallocate( LH_sfc%n )
        deallocate( LH_sfc%l_in_update )


        deallocate( LH_sfc%f%var )
        deallocate( LH_sfc%f%z )
        deallocate( LH_sfc%f%rlat )
        deallocate( LH_sfc%f%rlon )
      end if

      ! De-allocate all zm variables
      deallocate( zm%z )

      deallocate( zm%x )
      deallocate( zm%n )

      deallocate( zm%f%var )
      deallocate( zm%f%z )
      deallocate( zm%f%rlat )
      deallocate( zm%f%rlon )
      deallocate( zm%l_in_update )

      deallocate ( zmscr01 )
      deallocate ( zmscr02 )
      deallocate ( zmscr03 )
      deallocate ( zmscr04 )
      deallocate ( zmscr05 )
      deallocate ( zmscr06 )
      deallocate ( zmscr07 )
      deallocate ( zmscr08 )
      deallocate ( zmscr09 )
      deallocate ( zmscr10 )
      deallocate ( zmscr11 )
      deallocate ( zmscr12 )
      deallocate ( zmscr13 )
      deallocate ( zmscr14 )
      deallocate ( zmscr15 )
      deallocate ( zmscr16 )
      deallocate ( zmscr17 )

      if (l_output_rad_files) then
        ! De-allocate all rad_zt variables
        deallocate( rad_zt%z )

        deallocate( rad_zt%x )
        deallocate( rad_zt%n )

        deallocate( rad_zt%f%var )
        deallocate( rad_zt%f%z )
        deallocate( rad_zt%f%rlat )
        deallocate( rad_zt%f%rlon )
        deallocate( rad_zt%l_in_update )

        ! De-allocate all rad_zm variables
        deallocate( rad_zm%z )

        deallocate( rad_zm%x )
        deallocate( rad_zm%n )

        deallocate( rad_zm%f%var )
        deallocate( rad_zm%f%z )
        deallocate( rad_zm%l_in_update )

        !deallocate ( radscr01 )
        !deallocate ( radscr02 )
        !deallocate ( radscr03 )
        !deallocate ( radscr04 )
        !deallocate ( radscr05 )
        !deallocate ( radscr06 )
        !deallocate ( radscr07 )
        !deallocate ( radscr08 )
        !deallocate ( radscr09 )
        !deallocate ( radscr10 )
        !deallocate ( radscr11 )
        !deallocate ( radscr12 )
        !deallocate ( radscr13 )
        !deallocate ( radscr14 )
        !deallocate ( radscr15 )
        !deallocate ( radscr16 )
        !deallocate ( radscr17 )
      end if ! l_output_rad_files

      ! De-allocate all sfc variables
      deallocate( sfc%z )

      deallocate( sfc%x )
      deallocate( sfc%n )
      deallocate( sfc%l_in_update )

      deallocate( sfc%f%var )
      deallocate( sfc%f%z )
      deallocate( sfc%f%rlat )
      deallocate( sfc%f%rlon )

      ! De-allocate scalar indices
      deallocate( isclrm )
      deallocate( isclrm_f )
      deallocate( iedsclrm )
      deallocate( iedsclrm_f )
      deallocate( isclrprtp )
      deallocate( isclrp2 )
      deallocate( isclrpthvp )
      deallocate( isclrpthlp )
      deallocate( isclrprcp )
      deallocate( iwpsclrp )
      deallocate( iwp2sclrp )
      deallocate( iwpsclrp2 )
      deallocate( iwpsclrprtp )
      deallocate( iwpsclrpthlp )
      deallocate( iwpedsclrp )

    end if ! l_stats


    return
  end subroutine stats_finalize

!===============================================================================

end module stats_subs
