!-----------------------------------------------------------------------
!  $Id$
!===============================================================================
module stats_clubb_utilities

  implicit none

  private ! Set Default Scope

  public :: stats_init, stats_begin_timestep, stats_end_timestep, & 
    stats_accumulate, stats_finalize, stats_accumulate_hydromet, &
    stats_accumulate_lh_tend

  private :: stats_zero, stats_avg, stats_check_num_samples

  contains

  !-----------------------------------------------------------------------
  subroutine stats_init( iunit, fname_prefix, fdir, l_stats_in, &
                         stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                         nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
                         grad_zt, nnrad_zm, grad_zm, day, month, year, &
                         lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
                         stats_zt, stats_zm, stats_sfc, &
                         stats_lh_zt, stats_lh_sfc, &
                         stats_rad_zt, stats_rad_zm )
    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
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
        l_silhs_out ! Variable(s)

    use stats_variables, only: &
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
        zmscr17

    use stats_variables, only: &
        l_stats, &
        l_output_rad_files, &
        stats_tsamp,   &
        stats_tout,    &
        l_stats_samp,  &
        l_stats_last, &
        fname_zt, &
        fname_lh_zt, &
        fname_lh_sfc, &
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
        open_netcdf_for_writing     ! Procedure
#endif

    use stats_zm_module, only: &
        nvarmax_zm, & ! Constant(s)
        stats_init_zm ! Procedure(s)

    use stats_zt_module, only: &
        nvarmax_zt, & ! Constant(s)
        stats_init_zt ! Procedure(s)

    use stats_lh_zt_module, only: &
        nvarmax_lh_zt, & ! Constant(s)
        stats_init_lh_zt ! Procedure(s)

    use stats_lh_sfc_module, only: &
        nvarmax_lh_sfc, & ! Constant(s)
        stats_init_lh_sfc ! Procedure(s)

    use stats_rad_zt_module, only: &
        nvarmax_rad_zt, & ! Constant(s)
        stats_init_rad_zt ! Procedure(s)

    use stats_rad_zm_module, only: &
        nvarmax_rad_zm, & ! Constant(s)
        stats_init_rad_zm ! Procedure(s)

    use stats_sfc_module, only: &
        nvarmax_sfc, & ! Constant(s)
        stats_init_sfc ! Procedure(s)

    use constants_clubb, only: &
        fstdout, fstderr, var_length ! Constants

    use parameters_model, only: &
        hydromet_dim, &  ! Variable(s)
        sclr_dim, &
        edsclr_dim

    use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        err_code, &                     ! Error Indicator
        clubb_fatal_error               ! Constant

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    ! Local Constants
    integer, parameter :: &
      silhs_num_importance_categories = 8

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real( kind = core_rknd ), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      nlon, & ! Number of points in the X direction [-]
      nlat, & ! Number of points in the Y direction [-]
      nzmax   ! Grid points in the vertical         [-]

    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  & 
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  & 
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  & 
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = time_precision ), intent(in) ::  & 
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt, stats_lh_sfc)  [boolean]

    ! Local Variables
    logical :: l_error

    character(len=200) :: fname

    integer :: ivar, ntot, read_status

    ! Namelist Variables

    character(len=10) :: stats_fmt  ! File storage convention

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_lh_zt) ::  & 
      vars_lh_zt  ! Latin Hypercube variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_lh_sfc) ::  & 
      vars_lh_sfc  ! Latin Hypercube variables at the surface

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
      vars_lh_zt, &
      vars_lh_sfc, &
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
    l_silhs_out = l_silhs_out_in

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_last  = .false.
      return
    end if

    ! Initialize namelist variables

    vars_zt  = ''
    vars_zm  = ''
    vars_lh_zt = ''
    vars_lh_sfc = ''
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
      write(fstderr,*) "stats_init: Error reading stats namelist."
      err_code = clubb_fatal_error
      close(unit=iunit)
      return
    end if ! read_status /= 0

    close(unit=iunit)

    if ( clubb_at_least_debug_level( 1 ) ) then
      write(fstdout,*) "--------------------------------------------------"

      write(fstdout,*) "Statistics"

      write(fstdout,*) "--------------------------------------------------"
      write(fstdout,*) "vars_zt = "
      ivar = 1
      do while ( vars_zt(ivar) /= '' )
        write(fstdout,*) vars_zt(ivar)
        ivar = ivar + 1
      end do

      write(fstdout,*) "vars_zm = "
      ivar = 1
      do while ( vars_zm(ivar) /= '' )
        write(fstdout,*) vars_zm(ivar)
        ivar = ivar + 1
      end do

      if ( l_silhs_out ) then
        write(fstdout,*) "vars_lh_zt = "
        ivar = 1
        do while ( vars_lh_zt(ivar) /= '' )
          write(fstdout,*) vars_lh_zt(ivar)
          ivar = ivar + 1
        end do

        write(fstdout,*) "vars_lh_sfc = "
        ivar = 1
        do while ( vars_lh_sfc(ivar) /= '' )
          write(fstdout,*) vars_lh_sfc(ivar)
          ivar = ivar + 1
        end do
      end if ! l_silhs_out

      if ( l_output_rad_files ) then
        write(fstdout,*) "vars_rad_zt = "
        ivar = 1
        do while ( vars_rad_zt(ivar) /= '' )
          write(fstdout,*) vars_rad_zt(ivar)
          ivar = ivar + 1
        end do

        write(fstdout,*) "vars_rad_zm = "
        ivar = 1
        do while ( vars_rad_zm(ivar) /= '' )
          write(fstdout,*) vars_rad_zm(ivar)
          ivar = ivar + 1
        end do
      end if ! l_output_rad_files

      write(fstdout,*) "vars_sfc = "
      ivar = 1
      do while ( vars_sfc(ivar) /= '' )
        write(fstdout,*) vars_sfc(ivar)
        ivar = ivar + 1
      end do

      write(fstdout,*) "--------------------------------------------------"
    end if ! clubb_at_least_debug_level 1

    ! Determine file names for GrADS or NetCDF files
    fname_zt  = trim( fname_prefix )//"_zt"
    fname_zm  = trim( fname_prefix )//"_zm"
    fname_lh_zt  = trim( fname_prefix )//"_lh_zt"
    fname_lh_sfc  = trim( fname_prefix )//"_lh_sfc"
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
      write(fstderr,*) "In module stats_clubb_utilities subroutine stats_init: "
      write(fstderr,*) "Invalid stats output format "//trim( stats_fmt )
      err_code = clubb_fatal_error
      return

    end select

    ! Check sampling and output frequencies

    ! The model time step length, delt (which is dt_main), should multiply
    ! evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - real( floor( stats_tsamp/delt ), kind=core_rknd ) )  & 
           > 1.e-8_core_rknd) then
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
           - real( floor( stats_tout/stats_tsamp ), kind=core_rknd) ) & 
         > 1.e-8_core_rknd) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tout should be an even multiple of ',  &
                       'stats_tsamp.  Check the appropriate model.in file.'
      write(fstderr,*) 'stats_tout = ', stats_tout
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
    end if

    ! Initialize zt (mass points)

    ivar = 1
    do while ( ichar(vars_zt(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_zt(ivar)) /= 0 & 
               .and. ivar <= nvarmax_zt )
      ivar = ivar + 1
    end do
    ntot = ivar - 1

    if ( any( vars_zt == "hm_i" ) ) then
       ! Correct for number of variables found under "hm_i".
       ! Subtract "hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_hm_i" ) ) then
       ! Correct for number of variables found under "mu_hm_i".
       ! Subtract "mu_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_Ncn_i" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i".
       ! Subtract "mu_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "mu_hm_i_n" ) ) then
       ! Correct for number of variables found under "mu_hm_i_n".
       ! Subtract "mu_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "mu_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i_n".
       ! Subtract "mu_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "sigma_hm_i" ) ) then
       ! Correct for number of variables found under "sigma_hm_i".
       ! Subtract "sigma_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "sigma_Ncn_i" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i".
       ! Subtract "sigma_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "sigma_hm_i_n" ) ) then
       ! Correct for number of variables found under "sigma_hm_i_n".
       ! Subtract "sigma_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "sigma_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i_n".
       ! Subtract "sigma_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif

    if ( any( vars_zt == "corr_w_hm_i" ) ) then
       ! Correct for number of variables found under "corr_w_hm_i".
       ! Subtract "corr_w_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_w_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_w_Ncn_i".
       ! Subtract "corr_w_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_chi_hm_i" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i".
       ! Subtract "corr_chi_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i".
       ! Subtract "corr_chi_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_eta_hm_i" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i".
       ! Subtract "corr_eta_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i".
       ! Subtract "corr_eta_Ncn_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i" ) ) then
       ! Correct for number of variables found under "corr_Ncn_hm_i".
       ! Subtract "corr_Ncn_hm_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i" ) ) then
       ! Correct for number of variables found under "corr_hmx_hmy_i".
       ! Subtract "corr_hmx_hmy_i" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of correlations of two hydrometeors, which is found by:
       ! (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zt statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 )
    endif

    if ( any( vars_zt == "corr_w_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_w_hm_i_n".
       ! Subtract "corr_w_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_w_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_w_Ncn_i_n".
       ! Subtract "corr_w_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_chi_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i_n".
       ! Subtract "corr_chi_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i_n".
       ! Subtract "corr_chi_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_eta_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i_n".
       ! Subtract "corr_eta_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i_n".
       ! Subtract "corr_eta_Ncn_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) to the number of zt
       ! statistical variables.
       ntot = ntot + 2
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_Ncn_hm_i_n".
       ! Subtract "corr_Ncn_hm_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) for each hydrometeor
       ! to the number of zt statistical variables.
       ntot = ntot + 2 * hydromet_dim
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i_n" ) ) then
       ! Correct for number of variables found under "corr_hmx_hmy_i_n".
       ! Subtract "corr_hmx_hmy_i_n" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of normal space correlations of two hydrometeors, which is
       ! found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zt statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 )
    endif

    if ( any( vars_zt == "hmp2_zt" ) ) then
       ! Correct for number of variables found under "hmp2_zt".
       ! Subtract "hmp2_zt" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zt statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zt == "wp2hmp" ) ) then
       ! Correct for number of variables found under "wp2hmp".
       ! Subtract "wp2hmp" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zt statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zt == "sclrm" ) ) then
       ! Correct for number of variables found under "sclrm".
       ! Subtract "sclrm" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + sclr_dim
    endif   

    if ( any( vars_zt == "sclrm_f" ) ) then
       ! Correct for number of variables found under "sclrm_f".
       ! Subtract "sclrm_f" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + sclr_dim
    endif

    if ( any( vars_zt == "edsclrm" ) ) then
       ! Correct for number of variables found under "edsclrm".
       ! Subtract "edsclrm" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + edsclr_dim
    endif

    if ( any( vars_zt == "edsclrm_f" ) ) then
       ! Correct for number of variables found under "edsclrm_f".
       ! Subtract "edsclrm_f" from the number of zt statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zt statistical variables.
       ntot = ntot + edsclr_dim
    endif

    if ( ntot >= nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      write(fstderr,*) "number of variables in vars_zt = ", ntot
      write(fstderr,*) "stats_init:  number of zt statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return
    end if

    stats_zt%num_output_fields = ntot
    stats_zt%kk = nzmax
    stats_zt%ii = nlon
    stats_zt%jj = nlat

    allocate( stats_zt%z( stats_zt%kk ) )
    stats_zt%z = gzt

    allocate( stats_zt%accum_field_values( stats_zt%ii, stats_zt%jj, &
      stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%accum_num_samples( stats_zt%ii, stats_zt%jj, &
      stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%l_in_update( stats_zt%ii, stats_zt%jj, stats_zt%kk, &
      stats_zt%num_output_fields ) )
    call stats_zero( stats_zt%ii, stats_zt%jj, stats_zt%kk, stats_zt%num_output_fields, & ! In
      stats_zt%accum_field_values, stats_zt%accum_num_samples, stats_zt%l_in_update ) ! Out

    allocate( stats_zt%file%grid_avg_var( stats_zt%num_output_fields ) )
    allocate( stats_zt%file%z( stats_zt%kk ) )

    ! Allocate scratch space

    allocate( ztscr01(stats_zt%kk) )
    allocate( ztscr02(stats_zt%kk) )
    allocate( ztscr03(stats_zt%kk) )
    allocate( ztscr04(stats_zt%kk) )
    allocate( ztscr05(stats_zt%kk) )
    allocate( ztscr06(stats_zt%kk) )
    allocate( ztscr07(stats_zt%kk) )
    allocate( ztscr08(stats_zt%kk) )
    allocate( ztscr09(stats_zt%kk) )
    allocate( ztscr10(stats_zt%kk) )
    allocate( ztscr11(stats_zt%kk) )
    allocate( ztscr12(stats_zt%kk) )
    allocate( ztscr13(stats_zt%kk) )
    allocate( ztscr14(stats_zt%kk) )
    allocate( ztscr15(stats_zt%kk) )
    allocate( ztscr16(stats_zt%kk) )
    allocate( ztscr17(stats_zt%kk) )
    allocate( ztscr18(stats_zt%kk) )
    allocate( ztscr19(stats_zt%kk) )
    allocate( ztscr20(stats_zt%kk) )
    allocate( ztscr21(stats_zt%kk) )

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
      call open_grads( iunit, fdir, fname,  &  ! intent(in)
                       1, stats_zt%kk, nlat, nlon, stats_zt%z, & ! intent(in) 
                       day, month, year, lat_vals, lon_vals, &  ! intent(in)
                       time_current+real(stats_tout,kind=time_precision), stats_tout, &!intent(in) 
                       stats_zt%num_output_fields, & ! intent(in)
                       stats_zt%file ) ! intent(inout)

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_zt%kk, stats_zt%z, &  ! In
                        day, month, year, lat_vals, lon_vals, &  ! In
                        time_current, stats_tout, stats_zt%num_output_fields, &  ! In
                        stats_zt%file ) ! InOut

      if ( err_code == clubb_fatal_error ) return
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif

    end if

    ! Default initialization for array indices for zt

    call stats_init_zt( vars_zt, l_error, & !intent(in)
                        stats_zt ) ! intent(inout)


    ! Setup output file for stats_lh_zt (Latin Hypercube stats)

    if ( l_silhs_out ) then

      ivar = 1
      do while ( ichar(vars_lh_zt(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_lh_zt(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_lh_zt )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( any( vars_lh_zt == "silhs_variance_category" ) ) then
        ! Correct for number of variables found under "silhs_variance_category".
        ! Subtract "silhs_variance_category" from the number of lh_zt statistical
        ! variables.
        ntot = ntot - 1
        ! Add 1 for each SILHS category to the number of lh_zt statistical variables
        ntot = ntot + silhs_num_importance_categories
      end if

      if ( any( vars_lh_zt == "lh_samp_frac_category" ) ) then
        ! Correct for number of variables found under "lh_samp_frac_category".
        ! Subtract "lh_samp_frac_category" from the number of lh_zt statistical
        ! variables.
        ntot = ntot - 1
        ! Add 1 for each SILHS category to the number of lh_zt statistical variables
        ntot = ntot + silhs_num_importance_categories
      end if

      if ( ntot == nvarmax_lh_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_lh_zt than allowed for by nvarmax_lh_zt."
        write(fstderr,*) "Check the number of variables listed for vars_lh_zt ",  &
                         "in the stats namelist, or change nvarmax_lh_zt."
        write(fstderr,*) "nvarmax_lh_zt = ", nvarmax_lh_zt
        write(fstderr,*) "number of variables in vars_lh_zt = ", ntot
        write(fstderr,*) "stats_init:  number of lh_zt statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      stats_lh_zt%num_output_fields = ntot
      stats_lh_zt%kk = nzmax
      stats_lh_zt%ii = nlon
      stats_lh_zt%jj = nlat

      allocate( stats_lh_zt%z( stats_lh_zt%kk ) )
      stats_lh_zt%z = gzt

      allocate( stats_lh_zt%accum_field_values( stats_lh_zt%ii, stats_lh_zt%jj, &
        stats_lh_zt%kk, stats_lh_zt%num_output_fields ) )
      allocate( stats_lh_zt%accum_num_samples( stats_lh_zt%ii, stats_lh_zt%jj, &
        stats_lh_zt%kk, stats_lh_zt%num_output_fields ) )
      allocate( stats_lh_zt%l_in_update( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, &
        stats_lh_zt%num_output_fields ) )
      call stats_zero( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, & ! intent(in)
        stats_lh_zt%num_output_fields, & ! intent(in)
        stats_lh_zt%accum_field_values, stats_lh_zt%accum_num_samples, & ! intent(out)
        stats_lh_zt%l_in_update ) ! intent(out)

      allocate( stats_lh_zt%file%grid_avg_var( stats_lh_zt%num_output_fields ) )
      allocate( stats_lh_zt%file%z( stats_lh_zt%kk ) )


      fname = trim( fname_lh_zt )

      if ( l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & ! intent(in)
                         1, stats_lh_zt%kk, nlat, nlon, stats_lh_zt%z, & ! intent(in)
                         day, month, year, lat_vals, lon_vals, &  ! intent(in)
                         time_current+real(stats_tout,kind=time_precision), stats_tout, & ! In
                         stats_lh_zt%num_output_fields, & ! intent(in)
                         stats_lh_zt%file ) ! intent(inout)

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_lh_zt%kk, &  ! In
                          stats_lh_zt%z, day, month, year, lat_vals, lon_vals, &  ! In
                          time_current, stats_tout, stats_lh_zt%num_output_fields, &  ! In
                          stats_lh_zt%file ) ! InOut

        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

      call stats_init_lh_zt( vars_lh_zt, l_error, & !intent(in)
                             stats_lh_zt ) ! intent(inout)

      ivar = 1
      do while ( ichar(vars_lh_sfc(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_lh_sfc(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_lh_sfc )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_lh_sfc ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_lh_sfc than allowed for by nvarmax_lh_sfc."
        write(fstderr,*) "Check the number of variables listed for vars_lh_sfc ",  &
                         "in the stats namelist, or change nvarmax_lh_sfc."
        write(fstderr,*) "nvarmax_lh_sfc = ", nvarmax_lh_sfc
        write(fstderr,*) "number of variables in vars_lh_sfc = ", ntot
        write(fstderr,*) "stats_init:  number of lh_sfc statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      stats_lh_sfc%num_output_fields = ntot
      stats_lh_sfc%kk = 1
      stats_lh_sfc%ii = nlon
      stats_lh_sfc%jj = nlat

      allocate( stats_lh_sfc%z( stats_lh_sfc%kk ) )
      stats_lh_sfc%z = gzm(1)

      allocate( stats_lh_sfc%accum_field_values( stats_lh_sfc%ii, stats_lh_sfc%jj, &
        stats_lh_sfc%kk, stats_lh_sfc%num_output_fields ) )
      allocate( stats_lh_sfc%accum_num_samples( stats_lh_sfc%ii, stats_lh_sfc%jj, &
        stats_lh_sfc%kk, stats_lh_sfc%num_output_fields ) )
      allocate( stats_lh_sfc%l_in_update( stats_lh_sfc%ii, stats_lh_sfc%jj, &
        stats_lh_sfc%kk, stats_lh_sfc%num_output_fields ) )

      call stats_zero( stats_lh_sfc%ii, stats_lh_sfc%jj, stats_lh_sfc%kk, & ! intent(in)
          stats_lh_sfc%num_output_fields, & ! intent(in)
          stats_lh_sfc%accum_field_values, & ! intent(out)
          stats_lh_sfc%accum_num_samples, stats_lh_sfc%l_in_update ) ! intent(out)

      allocate( stats_lh_sfc%file%grid_avg_var( stats_lh_sfc%num_output_fields ) )
      allocate( stats_lh_sfc%file%z( stats_lh_sfc%kk ) )

      fname = trim( fname_lh_sfc )

      if ( l_grads ) then

        ! Open GrADS file
        call open_grads( iunit, fdir, fname,  & ! intent(in)
                         1, stats_lh_sfc%kk, nlat, nlon, stats_lh_sfc%z, & ! intent(in) 
                         day, month, year, lat_vals, lon_vals, &  ! intent(in)
                         time_current+real(stats_tout,kind=time_precision), stats_tout, &  ! In
                         stats_lh_sfc%num_output_fields, & ! intent(in)
                         stats_lh_sfc%file ) ! intent(inout)

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_lh_sfc%kk, &  ! In
                          stats_lh_sfc%z, day, month, year, lat_vals, lon_vals, &  ! In
                          time_current, stats_tout, stats_lh_sfc%num_output_fields, &  ! In
                          stats_lh_sfc%file ) ! InOut

        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif

      end if

      call stats_init_lh_sfc( vars_lh_sfc, l_error, & !intent(in)
                              stats_lh_sfc ) ! intent(inout)

    end if ! l_silhs_out

    ! Initialize stats_zm (momentum points)

    ivar = 1
    do while ( ichar(vars_zm(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_zm(ivar)) /= 0 & 
               .and. ivar <= nvarmax_zm )
      ivar = ivar + 1
    end do
    ntot = ivar - 1

    if ( any( vars_zm == "hydrometp2" ) ) then
       ! Correct for number of variables found under "hydrometp2".
       ! Subtract "hydrometp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "wphydrometp" ) ) then
       ! Correct for number of variables found under "wphydrometp".
       ! Subtract "wphydrometp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "rtphmp" ) ) then
       ! Correct for number of variables found under "rtphmp".
       ! Subtract "rtphmp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "thlphmp" ) ) then
       ! Correct for number of variables found under "thlphmp".
       ! Subtract "thlphmp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "hmxphmyp" ) ) then
       ! Correct for number of variables found under "hmxphmyp".
       ! Subtract "hmxphmyp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add the number of overall covariances of two hydrometeors, which is
       ! found by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! to the number of zm statistical variables.
       ntot = ntot + hydromet_dim * ( hydromet_dim - 1 ) / 2
    endif

    if ( any( vars_zm == "K_hm" ) ) then
       ! Correct for number of variables found under "K_hm".
       ! Subtract "K_hm" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each hydrometeor to the number of zm statistical variables.
       ntot = ntot + hydromet_dim
    endif

    if ( any( vars_zm == "sclrprtp" ) ) then
       ! Correct for number of variables found under "sclrprtp".
       ! Subtract "sclrprtp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif

    if ( any( vars_zm == "sclrp2" ) ) then
       ! Correct for number of variables found under "sclrp2".
       ! Subtract "sclrp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrpthvp" ) ) then
       ! Correct for number of variables found under "sclrpthvp".
       ! Subtract "sclrpthvp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrpthlp" ) ) then
       ! Correct for number of variables found under "sclrpthlp".
       ! Subtract "sclrpthlp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "sclrprcp" ) ) then
       ! Correct for number of variables found under "sclrprcp".
       ! Subtract "sclrprcp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrp" ) ) then
       ! Correct for number of variables found under "wpsclrp".
       ! Subtract "wpsclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrp2" ) ) then
       ! Correct for number of variables found under "wpsclrp2".
       ! Subtract "wpsclrp2" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wp2sclrp" ) ) then
       ! Correct for number of variables found under "wp2sclrp".
       ! Subtract "wp2sclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrprtp" ) ) then
       ! Correct for number of variables found under "wpsclrprtp".
       ! Subtract "wpsclrprtp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpsclrpthlp" ) ) then
       ! Correct for number of variables found under "wpsclrpthlp".
       ! Subtract "wpsclrpthlp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + sclr_dim
    endif


    if ( any( vars_zm == "wpedsclrp" ) ) then
       ! Correct for number of variables found under "wpedsclrp".
       ! Subtract "wpedsclrp" from the number of zm statistical variables.
       ntot = ntot - 1
       ! Add 1 for each scalar to the number of zm statistical variables.
       ntot = ntot + edsclr_dim
    endif



    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      write(fstderr,*) "number of variables in vars_zm = ", ntot
      write(fstderr,*) "stats_init:  number of zm statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return
    end if

    stats_zm%num_output_fields = ntot
    stats_zm%kk = nzmax
    stats_zm%ii = nlon
    stats_zm%jj = nlat

    allocate( stats_zm%z( stats_zm%kk ) )
    stats_zm%z = gzm

    allocate( stats_zm%accum_field_values( stats_zm%ii, stats_zm%jj, &
      stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%accum_num_samples( stats_zm%ii, stats_zm%jj, &
      stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%l_in_update( stats_zm%ii, stats_zm%jj, stats_zm%kk, &
      stats_zm%num_output_fields ) )

    call stats_zero( stats_zm%ii, stats_zm%jj, stats_zm%kk, stats_zm%num_output_fields, & ! In
      stats_zm%accum_field_values, stats_zm%accum_num_samples, stats_zm%l_in_update ) ! intent(out)

    allocate( stats_zm%file%grid_avg_var( stats_zm%num_output_fields ) )
    allocate( stats_zm%file%z( stats_zm%kk ) )

    ! Allocate scratch space

    allocate( zmscr01(stats_zm%kk) )
    allocate( zmscr02(stats_zm%kk) )
    allocate( zmscr03(stats_zm%kk) )
    allocate( zmscr04(stats_zm%kk) )
    allocate( zmscr05(stats_zm%kk) )
    allocate( zmscr06(stats_zm%kk) )
    allocate( zmscr07(stats_zm%kk) )
    allocate( zmscr08(stats_zm%kk) )
    allocate( zmscr09(stats_zm%kk) )
    allocate( zmscr10(stats_zm%kk) )
    allocate( zmscr11(stats_zm%kk) )
    allocate( zmscr12(stats_zm%kk) )
    allocate( zmscr13(stats_zm%kk) )
    allocate( zmscr14(stats_zm%kk) )
    allocate( zmscr15(stats_zm%kk) )
    allocate( zmscr16(stats_zm%kk) )
    allocate( zmscr17(stats_zm%kk) )

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
      call open_grads( iunit, fdir, fname,  & ! intent(in)
                       1, stats_zm%kk, nlat, nlon, stats_zm%z, & ! intent(in)
                       day, month, year, lat_vals, lon_vals, & ! intent(in)
                       time_current+real(stats_tout,kind=time_precision), stats_tout, & !intent(in)
                       stats_zm%num_output_fields, & ! intent(in)
                       stats_zm%file ) ! intent(inout)

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_zm%kk, stats_zm%z, &  ! In
                        day, month, year, lat_vals, lon_vals, &  ! In
                        time_current, stats_tout, stats_zm%num_output_fields, &  ! In
                        stats_zm%file ) ! InOut

      if ( err_code == clubb_fatal_error ) return
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    call stats_init_zm( vars_zm, l_error, & !intent(in)
                        stats_zm ) ! intent(inout)

    ! Initialize stats_rad_zt (radiation points)

    if (l_output_rad_files) then

      ivar = 1
      do while ( ichar(vars_rad_zt(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zt(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_rad_zt )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        write(fstderr,*) "number of variables in vars_rad_zt = ", ntot
        write(fstderr,*) "stats_init:  number of rad_zt statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      stats_rad_zt%num_output_fields = ntot
      stats_rad_zt%kk = nnrad_zt
      stats_rad_zt%ii = nlon
      stats_rad_zt%jj = nlat
      allocate( stats_rad_zt%z( stats_rad_zt%kk ) )
      stats_rad_zt%z = grad_zt

      allocate( stats_rad_zt%accum_field_values( stats_rad_zt%ii, stats_rad_zt%jj, &
        stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%accum_num_samples( stats_rad_zt%ii, stats_rad_zt%jj, &
        stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%l_in_update( stats_rad_zt%ii, stats_rad_zt%jj, &
        stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )

      call stats_zero( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zt%kk, & ! intent(in)
                       stats_rad_zt%num_output_fields, & ! intent(in)
                       stats_rad_zt%accum_field_values, & ! intent(out)
                       stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )! intent(out)

      allocate( stats_rad_zt%file%grid_avg_var( stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%file%z( stats_rad_zt%kk ) )

      fname = trim( fname_rad_zt )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & ! intent(in)
                         1, stats_rad_zt%kk, nlat, nlon, stats_rad_zt%z, & ! intent(in)
                         day, month, year, lat_vals, lon_vals, & 
                         time_current+real(stats_tout, kind=time_precision), stats_tout, & ! In
                         stats_rad_zt%num_output_fields, & ! intent(in)
                         stats_rad_zt%file ) ! intent(inout)

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname,  & ! intent(in)
                          1, stats_rad_zt%kk, stats_rad_zt%z, & ! intent(in)
                          day, month, year, lat_vals, lon_vals, & ! intent(in)
                          time_current, stats_tout, & ! intent(in)
                          stats_rad_zt%num_output_fields, & ! intent(in)
                          stats_rad_zt%file ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

      call stats_init_rad_zt( vars_rad_zt, l_error, & !intent(in)
                              stats_rad_zt ) ! intent(inout)

      ! Initialize stats_rad_zm (radiation points)

      ivar = 1
      do while ( ichar(vars_rad_zm(ivar)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zm(ivar)) /= 0 & 
                 .and. ivar <= nvarmax_rad_zm )
        ivar = ivar + 1
      end do
      ntot = ivar - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        write(fstderr,*) "number of variables in vars_rad_zm = ", ntot
        write(fstderr,*) "stats_init:  number of rad_zm statistical variables exceeds limit"
        err_code = clubb_fatal_error
        return
      end if

      stats_rad_zm%num_output_fields = ntot
      stats_rad_zm%kk = nnrad_zm
      stats_rad_zm%ii = nlon
      stats_rad_zm%jj = nlat

      allocate( stats_rad_zm%z( stats_rad_zm%kk ) )
      stats_rad_zm%z = grad_zm

      allocate( stats_rad_zm%accum_field_values( stats_rad_zm%ii, stats_rad_zm%jj, &
        stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
      allocate( stats_rad_zm%accum_num_samples( stats_rad_zm%ii, stats_rad_zm%jj, &
        stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
      allocate( stats_rad_zm%l_in_update( stats_rad_zm%ii, stats_rad_zm%jj, &
        stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )

      call stats_zero( stats_rad_zm%ii, stats_rad_zm%jj, stats_rad_zm%kk, & ! intent(in)
                       stats_rad_zm%num_output_fields, & ! intent(in)
                       stats_rad_zm%accum_field_values, & ! intent(out)
                       stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update ) ! intent(out)

      allocate( stats_rad_zm%file%grid_avg_var( stats_rad_zm%num_output_fields ) )
      allocate( stats_rad_zm%file%z( stats_rad_zm%kk ) )

      fname = trim( fname_rad_zm )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & ! intent(in)
                         1, stats_rad_zm%kk, nlat, nlon, stats_rad_zm%z, & ! intent(in)
                         day, month, year, lat_vals, lon_vals, & 
                         time_current+real(stats_tout,kind=time_precision), stats_tout, & ! In
                         stats_rad_zm%num_output_fields, & ! intent(in)
                         stats_rad_zm%file ) ! intent(inout)

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf_for_writing( nlat, nlon, fdir, fname,  & ! intent(in)
                          1, stats_rad_zm%kk, stats_rad_zm%z, & ! intent(in)
                          day, month, year, lat_vals, lon_vals, & ! intent(in)
                          time_current, stats_tout, & ! intent(in)
                          stats_rad_zm%num_output_fields, & ! intent(in)
                          stats_rad_zm%file ) ! intent(inout)

        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This CLUBB program was not compiled with netCDF support."
#endif
      end if

      call stats_init_rad_zm( vars_rad_zm, l_error, & !intent(in)
                              stats_rad_zm ) ! intent(inout)
    end if ! l_output_rad_files


    ! Initialize stats_sfc (surface point)

    ivar = 1
    do while ( ichar(vars_sfc(ivar)(1:1)) /= 0  & 
               .and. len_trim(vars_sfc(ivar)) /= 0 & 
               .and. ivar <= nvarmax_sfc )
      ivar = ivar + 1
    end do
    ntot = ivar - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      write(fstderr,*) "number of variables in vars_sfc = ", ntot
      write(fstderr,*) "stats_init:  number of sfc statistical variables exceeds limit"
      err_code = clubb_fatal_error
      return

    end if

    stats_sfc%num_output_fields = ntot
    stats_sfc%kk = 1
    stats_sfc%ii = nlon
    stats_sfc%jj = nlat

    allocate( stats_sfc%z( stats_sfc%kk ) )
    stats_sfc%z = gzm(1)

    allocate( stats_sfc%accum_field_values( stats_sfc%ii, stats_sfc%jj, &
      stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%accum_num_samples( stats_sfc%ii, stats_sfc%jj, &
      stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%l_in_update( stats_sfc%ii, stats_sfc%jj, &
      stats_sfc%kk, stats_sfc%num_output_fields ) )

    call stats_zero( stats_sfc%ii, stats_sfc%jj, stats_sfc%kk, stats_sfc%num_output_fields, & ! In
      stats_sfc%accum_field_values, stats_sfc%accum_num_samples, stats_sfc%l_in_update ) ! out

    allocate( stats_sfc%file%grid_avg_var( stats_sfc%num_output_fields ) )
    allocate( stats_sfc%file%z( stats_sfc%kk ) )

    fname = trim( fname_sfc )

    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & ! intent(in)
                       1, stats_sfc%kk, nlat, nlon, stats_sfc%z, & ! intent(in)
                       day, month, year, lat_vals, lon_vals, & ! intent(in)
                       time_current+real(stats_tout,kind=time_precision), stats_tout, & !intent(in)
                       stats_sfc%num_output_fields, & ! intent(in)
                       stats_sfc%file ) ! intent(inout)

    else ! Open NetCDF files
#ifdef NETCDF
      call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, stats_sfc%kk, stats_sfc%z, &  ! In
                        day, month, year, lat_vals, lon_vals, &  ! In
                        time_current, stats_tout, stats_sfc%num_output_fields, &  ! In
                        stats_sfc%file ) ! InOut

      if ( err_code == clubb_fatal_error ) return
#else
      error stop "This CLUBB program was not compiled with netCDF support."
#endif
    end if

    call stats_init_sfc( vars_sfc, l_error, & !intent(in)
                         stats_sfc ) ! intent(inout)

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      err_code = clubb_fatal_error
      return
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
  subroutine stats_zero( ii, jj, kk, nn, &
                         x, n, l_in_update )

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
    integer, intent(in) :: ii, jj, kk, nn

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(ii,jj,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(ii,jj,kk,nn), intent(out) :: n
    logical, dimension(ii,jj,kk,nn), intent(out) :: l_in_update

    ! Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0_stat_rknd
      n(:,:,:,:) = 0_stat_nknd
      l_in_update(:,:,:,:) = .false.
    end if

    return
  end subroutine stats_zero

  !-----------------------------------------------------------------------
  subroutine stats_avg( ii, jj, kk, nn, n, &
                        x )

    ! Description:
    !   Compute the average of stats fields
    ! References:
    !   None
    !-----------------------------------------------------------------------
    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    use stat_file_module, only: &
        clubb_i, clubb_j ! Variable(s)

    implicit none

    ! External
    intrinsic :: real

    ! Input Variable(s)
    integer, intent(in) :: &
      ii, & ! Number of points in X (i.e. latitude) dimension
      jj, & ! Number of points in Y (i.e. longitude) dimension
      kk, & ! Number of levels in vertical (i.e. Z) dimension
      nn    ! Number of variables being output to disk (e.g. cloud_frac, rain rate, etc.)

    integer(kind=stat_nknd), dimension(ii,jj,kk,nn), intent(in) :: &
      n ! n is the number of samples for each of the nn fields 
        ! and each of the kk vertical levels

    ! Output Variable(s)
    real(kind=stat_rknd), dimension(ii,jj,kk,nn), intent(inout) :: &
      x ! The variable x contains the cumulative sums of n sample values of each of
        ! the nn output fields (e.g. the sum of the sampled rain rate values)

    ! ---- Begin Code ----

    ! Compute averages
    where ( n(1,1,1:kk,1:nn) > 0 )
      x(clubb_i,clubb_j,1:kk,1:nn) = x(clubb_i,clubb_j,1:kk,1:nn) &
         / real( n(clubb_i,clubb_j,1:kk,1:nn), kind=stat_rknd )
    end where

    return
  end subroutine stats_avg

  !-----------------------------------------------------------------------
  subroutine stats_begin_timestep( itime, stats_nsamp, stats_nout)

    !     Description:
    !       Given the elapsed time, set flags determining specifics such as
    !       if this time set should be sampled or if this is the first or
    !       last time step.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        l_stats,  & ! Variable(s)
        l_stats_samp, & 
        l_stats_last 


    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) ::  & 
      itime, &       ! Elapsed model time       [timestep]
      stats_nsamp, & ! Stats sampling interval  [timestep]
      stats_nout     ! Stats output interval    [timestep]

    if ( .not. l_stats ) return

    ! Only sample time steps that are multiples of "stats_tsamp"
    ! in a case's "model.in" file to shorten length of run
    if ( mod( itime, stats_nsamp ) == 0 ) then
      l_stats_samp = .true.
    else
      l_stats_samp = .false.
    end if

    ! Indicates the end of the sampling time period. Signals to start writing to the file
    if ( mod( itime, stats_nout ) == 0 ) then
      l_stats_last = .true.
    else
      l_stats_last = .false.
    end if
   
    return

  end subroutine stats_begin_timestep

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep( clubb_params, &              ! intent(in)
                                 stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                                 stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                                 stats_rad_zt, stats_rad_zm & ! intent(inout)
#ifdef NETCDF
                                 , l_uv_nudge, &
                                 l_tke_aniso, &
                                 l_standard_term_ta &
#endif
                                  )

    ! Description: 
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        l_stats_last, &
        l_output_rad_files, & 
        l_grads, &
        l_silhs_out

    use output_grads, only: &
        write_grads ! Procedure(s)

    use stat_file_module, only: &
        clubb_i, & ! Variable(s)
        clubb_j

    use parameter_indices, only: &
        nparams    ! Variable(s)

#ifdef NETCDF
    use output_netcdf, only: & 
        write_netcdf ! Procedure(s)
#endif

    use error_code, only : &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

    use stats_type, only: stats ! Type

    implicit none

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    ! External
    intrinsic :: floor

#ifdef NETCDF
    ! Input Variables
    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging.
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e.
                            ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.
#endif

    ! Local Variables

    logical :: l_error

    ! ---- Begin Code ----

    ! Check if it is time to write to file

    if ( .not. l_stats_last ) return

    ! Initialize
    l_error = .false.

    call stats_check_num_samples( stats_zt, & ! intent(in)
                                  l_error ) ! intent(inout)
    call stats_check_num_samples( stats_zm, & ! intent(in)
                                  l_error ) ! intent(inout)
    call stats_check_num_samples( stats_sfc, & ! intent(in)
                                  l_error ) ! intent(inout)
    if ( l_silhs_out ) then
      call stats_check_num_samples( stats_lh_zt, & ! intent(in)
                                    l_error ) ! intent(inout)
      call stats_check_num_samples( stats_lh_sfc, & ! intent(in)
                                    l_error ) ! intent(inout)
    end if
    if ( l_output_rad_files ) then
      call stats_check_num_samples( stats_rad_zt, & ! intent(in)
                                    l_error ) ! intent(inout)
      call stats_check_num_samples( stats_rad_zm, & ! intent(in)
                                    l_error ) ! intent(inout)
    end if

    ! Return if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      write(fstderr,*) 'stats_end_timestep:  error(s) found'
      err_code = clubb_fatal_error
      return
    end if ! l_error

    ! Compute averages
    call stats_avg( stats_zt%ii, stats_zt%jj, stats_zt%kk, stats_zt%num_output_fields, & ! In
                    stats_zt%accum_num_samples, & ! intent(in)
                    stats_zt%accum_field_values ) ! intent(inout)
    call stats_avg( stats_zm%ii, stats_zm%jj, stats_zm%kk, stats_zm%num_output_fields, & ! In
                    stats_zm%accum_num_samples, & ! intent(in)
                    stats_zm%accum_field_values ) ! intent(inout)
    if ( l_silhs_out ) then
      call stats_avg( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, & ! intent(in)
         stats_lh_zt%num_output_fields, stats_lh_zt%accum_num_samples, & ! intent(in)
         stats_lh_zt%accum_field_values ) ! intent(inout)
      call stats_avg( stats_lh_sfc%ii, stats_lh_sfc%jj, stats_lh_sfc%kk, & ! intent(in)
        stats_lh_sfc%num_output_fields, stats_lh_sfc%accum_num_samples, & ! intent(in)
        stats_lh_sfc%accum_field_values ) ! intent(inout)
    end if
    if ( l_output_rad_files ) then
      call stats_avg( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zt%kk, & ! intent(in)
        stats_rad_zt%num_output_fields, & ! intent(in)
        stats_rad_zt%accum_num_samples, & ! intent(in)
        stats_rad_zt%accum_field_values ) ! intent(inout)
      call stats_avg( stats_rad_zm%ii, stats_rad_zm%jj, stats_rad_zm%kk, & ! intent(in)
        stats_rad_zm%num_output_fields, & ! intent(in)
        stats_rad_zm%accum_num_samples, & ! intent(in)
        stats_rad_zm%accum_field_values ) ! intent(inout)
    end if
    call stats_avg( stats_sfc%ii, stats_sfc%jj, stats_sfc%kk, stats_sfc%num_output_fields, & ! In
        stats_sfc%accum_num_samples, & ! intent(in)
        stats_sfc%accum_field_values ) ! intent(inout)

    ! Only write to the file and zero out the stats fields if we've reach the horizontal
    ! limits of the domain (this is always true in the single-column case because it's 1x1).
    if ( clubb_i == stats_zt%ii .and. clubb_j == stats_zt%jj ) then
      ! Write to file
      if ( l_grads ) then
        call write_grads( stats_zt%file  ) ! intent(inout)
        call write_grads( stats_zm%file  ) ! intent(inout)
        if ( l_silhs_out ) then
          call write_grads( stats_lh_zt%file  ) ! intent(inout)
          call write_grads( stats_lh_sfc%file  ) ! intent(inout)
        end if
        if ( l_output_rad_files ) then
          call write_grads( stats_rad_zt%file  ) ! intent(inout)
          call write_grads( stats_rad_zm%file  ) ! intent(inout)
        end if
        call write_grads( stats_sfc%file  ) ! intent(inout)
      else ! l_netcdf

#ifdef NETCDF
        call write_netcdf( clubb_params, & ! intent(in)
                           l_uv_nudge, & ! intent(in)
                           l_tke_aniso, & ! intent(in)
                           l_standard_term_ta, & ! intent(in)
                           stats_zt%file  ) ! intent(inout)
        call write_netcdf( clubb_params, & ! intent(in)
                           l_uv_nudge, & ! intent(in)
                           l_tke_aniso, & ! intent(in)
                           l_standard_term_ta, & ! intent(in)
                           stats_zm%file  ) ! intent(inout)
        if ( l_silhs_out ) then
          call write_netcdf( clubb_params, & ! intent(in)
                             l_uv_nudge, & ! intent(in)
                             l_tke_aniso, & ! intent(in)
                             l_standard_term_ta, & ! intent(in)
                             stats_lh_zt%file  ) ! intent(inout)
          call write_netcdf( clubb_params, & ! intent(in)
                             l_uv_nudge, & ! intent(in)
                             l_tke_aniso, & ! intent(in)
                             l_standard_term_ta, & ! intent(in)
                             stats_lh_sfc%file  ) ! intent(inout)
        end if
        if ( l_output_rad_files ) then
          call write_netcdf( clubb_params, & ! intent(in)
                             l_uv_nudge, & ! intent(in)
                             l_tke_aniso, & ! intent(in)
                             l_standard_term_ta, & ! intent(in)
                             stats_rad_zt%file  ) ! intent(inout)
          call write_netcdf( clubb_params, & ! intent(in)
                             l_uv_nudge, & ! intent(in)
                             l_tke_aniso, & ! intent(in)
                             l_standard_term_ta, & ! intent(in)
                             stats_rad_zm%file  ) ! intent(inout)
        end if
        call write_netcdf( clubb_params, & ! intent(in)
                           l_uv_nudge, & ! intent(in)
                           l_tke_aniso, & ! intent(in)
                           l_standard_term_ta, & ! intent(in)
                           stats_sfc%file  ) ! intent(inout)
            
        if ( err_code == clubb_fatal_error ) return
#else
        error stop "This program was not compiled with netCDF support"
#endif /* NETCDF */
      end if ! l_grads

      ! Reset sample fields
      call stats_zero( stats_zt%ii, stats_zt%jj, stats_zt%kk, stats_zt%num_output_fields, & ! In
      stats_zt%accum_field_values, stats_zt%accum_num_samples, stats_zt%l_in_update ) ! out
      call stats_zero( stats_zm%ii, stats_zm%jj, stats_zm%kk, stats_zm%num_output_fields, & ! In
        stats_zm%accum_field_values, stats_zm%accum_num_samples, stats_zm%l_in_update ) ! Out
      if ( l_silhs_out ) then
        call stats_zero( stats_lh_zt%ii, stats_lh_zt%jj, stats_lh_zt%kk, & ! intent(in)
          stats_lh_zt%num_output_fields, & ! intent(in)
          stats_lh_zt%accum_field_values, & ! intent(out)
          stats_lh_zt%accum_num_samples, stats_lh_zt%l_in_update ) ! intent(out)
        call stats_zero( stats_lh_sfc%ii, stats_lh_sfc%jj, stats_lh_sfc%kk, & ! intent(in)
          stats_lh_sfc%num_output_fields, & ! intent(in)
          stats_lh_sfc%accum_field_values, & ! intent(out)
          stats_lh_sfc%accum_num_samples, stats_lh_sfc%l_in_update ) ! intent(out)
      end if
      if ( l_output_rad_files ) then
        call stats_zero( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zt%kk, & ! intent(in)
          stats_rad_zt%num_output_fields, & ! intent(in)
          stats_rad_zt%accum_field_values, & ! intent(out)
          stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update ) ! intent(out)
        call stats_zero( stats_rad_zt%ii, stats_rad_zt%jj, stats_rad_zm%kk, & ! intent(in)
          stats_rad_zm%num_output_fields, & ! intent(in)
          stats_rad_zm%accum_field_values, & ! intent(out)
          stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update ) ! intent(out)
      end if
      call stats_zero( stats_sfc%ii, stats_sfc%jj, stats_sfc%kk, stats_sfc%num_output_fields, & !IN
        stats_sfc%accum_field_values, & ! intent(out)
        stats_sfc%accum_num_samples, stats_sfc%l_in_update ) ! intent(out)

    end if ! clubb_i = stats_zt%ii .and. clubb_j == stats_zt%jj


    return
  end subroutine stats_end_timestep

  !----------------------------------------------------------------------
  subroutine stats_accumulate & 
                   ( gr, um, vm, upwp, vpwp, up2, vp2, &
                     thlm, rtm, wprtp, wpthlp, &
                     wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
                     wpthvp, wp2thvp, rtpthvp, thlpthvp, &
                     p_in_Pa, exner, rho, rho_zm, &
                     rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
                     wm_zt, wm_zm, rcm, wprcp, rc_coef, rc_coef_zm, &
                     rcm_zm, rtm_zm, thlm_zm, cloud_frac, ice_supersat_frac, &
                     cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
                     cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
                     thvm, ug, vg, Lscale, wpthlp2, wp2thlp, wprtp2, wp2rtp, &
                     Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
                     wprtpthlp, sigma_sqd_w_zt, rsat, wp2_zt, thlp2_zt, &
                     wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, up2_zt, &
                     vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, & 
                     wp2up2, wp2vp2, wp4, &
                     tau_zm, Kh_zm, thlprcp, &
                     rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
                     wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
                     w_up_in_cloud, pdf_params, pdf_params_zm, sclrm, sclrp2, &
                     sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp, &
                     wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, wpsclrprtp, &
                     wpsclrpthlp, wpedsclrp, edsclrm, edsclrm_forcing, &
                     stats_zt, stats_zm, stats_sfc )

    ! Description:
    !   Accumulate those stats variables that are preserved in CLUBB from timestep to
    !   timestep, but not those stats that are not, (e.g. budget terms, longwave and
    !   shortwave components, etc.)
    !
    ! References:
    !   None
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        cloud_frac_min, &  ! Constant
        eps


    use pdf_utilities, only: &
        compute_variance_binormal    ! Procedure

    use stats_variables, only: & 
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
        icloud_cover, &
        ircm_supersat_adj

    use stats_variables, only: &
        ip_in_Pa, & ! Variable(s)
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
        iw_up_in_cloud, &
        iwprtpthlp, &
        irc_coef, &
        isigma_sqd_w_zt, & 
        irho, & 
        irsat, & 
        irsati

    use stats_variables, only: & 
        imixt_frac, &  ! Variable(s)
        iw_1, & 
        iw_2, & 
        ivarnce_w_1, & 
        ivarnce_w_2, & 
        ithl_1, & 
        ithl_2, & 
        ivarnce_thl_1, & 
        ivarnce_thl_2, & 
        irt_1, & 
        irt_2, & 
        ivarnce_rt_1, & 
        ivarnce_rt_2, & 
        irc_1, & 
        irc_2, & 
        irsatl_1, & 
        irsatl_2, & 
        icloud_frac_1, & 
        icloud_frac_2

    use stats_variables, only: & 
        ichi_1, & ! Variable(s)
        ichi_2, &
        istdev_chi_1, &
        istdev_chi_2, &
        ichip2, &
        istdev_eta_1, &
        istdev_eta_2, &
        icovar_chi_eta_1, &
        icovar_chi_eta_2, &
        icorr_w_chi_1, &
        icorr_w_chi_2, &
        icorr_w_eta_1, &
        icorr_w_eta_2, &
        icorr_chi_eta_1, &
        icorr_chi_eta_2, &
        icrt_1, &
        icrt_2, &
        icthl_1, &
        icthl_2, &
        icorr_w_rt_1, &
        icorr_w_rt_2, &
        icorr_w_thl_1, &
        icorr_w_thl_2, &
        icorr_rt_thl_1, &
        icorr_rt_thl_2, &
        ichi

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
        irtp3, &
        ithlp2, &
        ithlp3, &
        irtpthlp, &
        iwprtp,  &
        iwpthlp, &
        iwpup2, &
        iwpvp2, &
        iwp2up2, &
        iwp2vp2, &
        iwp4,  &
        iwpthvp, &
        irtpthvp

    use stats_variables, only: &
        ithlpthvp, & ! Variable(s)
        itau_zm, &
        iKh_zm, &
        iwprcp, &
        irc_coef_zm, &
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
        iw_1_zm, & ! Variable(s)
        iw_2_zm, &
        ivarnce_w_1_zm, &
        ivarnce_w_2_zm, &
        imixt_frac_zm

    use stats_variables, only: &
        iwp3_on_wp2, &
        iwp3_on_wp2_zt, &
        iSkw_velocity

    use stats_variables, only: &
        ia3_coef, & ! Variables
        ia3_coef_zt, &
        ircm_in_cloud
        
    use stats_variables, only: &
        itot_vartn_normlzd_rtm, &
        itot_vartn_normlzd_thlm, &
        itot_vartn_normlzd_wprtp

    use grid_class, only: & 
        grid ! Type

    use grid_class, only: & 
        zt2zm ! Procedure(s)

    use pdf_parameter_module, only: & 
        pdf_parameter ! Type

    use T_in_K_module, only: & 
        thlm2T_in_K ! Procedure

    use constants_clubb, only: & 
        rc_tol, fstderr    ! Constant(s)

    use parameters_model, only: & 
        sclr_dim,  &        ! Variable(s)
        edsclr_dim

    use stats_type_utilities, only: & 
        stat_update_var,  & ! Procedure(s)
        stat_update_var_pt

    use fill_holes, only: &
        vertical_avg, &     ! Procedure(s)
        vertical_integral

    use interpolation, only: & 
        lin_interpolate_two_points             ! Procedure

    use saturation, only: &
        sat_mixrat_ice ! Procedure

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type (grid), target, intent(in) :: gr

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      um,       & ! u wind (thermodynamic levels)          [m/s]
      vm,       & ! v wind (thermodynamic levels)          [m/s]
      upwp,     & ! vertical u momentum flux (m levs.)     [m^2/s^2]
      vpwp,     & ! vertical v momentum flux (m levs.)     [m^2/s^2]
      up2,      & ! < u'^2 > (momentum levels)             [m^2/s^2]
      vp2,      & ! < v'^2 > (momentum levels)             [m^2/s^2]
      thlm,     & ! liquid potential temperature (t levs.) [K]
      rtm,      & ! total water mixing ratio (t levs.)     [kg/kg]
      wprtp,    & ! < w' r_t' > (momentum levels)          [m/s kg/kg]
      wpthlp,   & ! < w' th_l' > (momentum levels)         [m/s K]
      wp2,      & ! < w'^2 > (momentum levels)             [m^2/s^2]
      wp3,      & ! < w'^3 > (thermodynamic levels)        [m^3/s^3]
      rtp2,     & ! < r_t'^2 > (momentum levels)           [(kg/kg)^2]
      rtp3,     & ! < r_t'^3 > (thermodynamic levels)      [(kg/kg)^3]
      thlp2,    & ! < th_l'^2 > (momentum levels)          [K^2]
      thlp3,    & ! < th_l'^3 > (thermodynamic levels)     [K^3]
      rtpthlp,  & ! < r_t' th_l' > (momentum levels)       [kg/kg K]
      wpthvp,   & ! < w' th_v' > (momentum levels)         [kg/kg K]
      wp2thvp,  & ! < w'^2 th_v' > (thermodynamic levels)  [m^2/s^2 K]
      rtpthvp,  & ! < r_t' th_v' > (momentum levels)       [kg/kg K]
      thlpthvp    ! < th_l' th_v' > (momentum levels)      [K^2]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
      exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
      rho,          & ! Density (thermodynamic levels)           [kg/m^3]
      rho_zm,       & ! Density on momentum levels               [kg/m^3]
      rho_ds_zm,    & ! Dry, static density (momentum levels)    [kg/m^3]
      rho_ds_zt,    & ! Dry, static density (thermo. levs.)      [kg/m^3]
      thv_ds_zm,    & ! Dry, base-state theta_v (momentum levs.) [K]
      thv_ds_zt,    & ! Dry, base-state theta_v (thermo. levs.)  [K]
      wm_zt,        & ! w on thermodynamic levels                [m/s]
      wm_zm           ! w on momentum levels                     [m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      rcm_zm,               & ! Cloud water mixing ratio on m levs.      [kg/kg]
      rtm_zm,               & ! Total water mixing ratio on m levs.      [kg/kg]
      thlm_zm,              & ! Liquid potential temperature on m levs.  [K]
      rcm,                  & ! Cloud water mixing ratio (t levs.)       [kg/kg]
      wprcp,                & ! < w' r_c' > (momentum levels)            [m/s kg/kg]
      rc_coef,              & ! Coefficient of X'r_c' (t-levs.)      [K/(kg/kg)]
      rc_coef_zm,           & ! Coefficient of X'r_c' on m-levs.     [K/(kg/kg)]
      cloud_frac,           & ! Cloud fraction (thermodynamic levels)    [-]
      ice_supersat_frac,    & ! Ice cloud fracion (thermodynamic levels) [-]
      cloud_frac_zm,        & ! Cloud fraction on zm levels              [-]
      ice_supersat_frac_zm, & ! Ice cloud fraction on zm levels          [-]
      rcm_in_layer,         & ! Cloud water mixing ratio in cloud layer  [kg/kg]
      cloud_cover,          & ! Cloud cover                              [-]
      rcm_supersat_adj        ! rcm adjustment due to supersaturation    [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      sigma_sqd_w    ! PDF width parameter (momentum levels)    [-]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
        thvm,           & ! Virtual potential temperature        [K]
        ug,             & ! u geostrophic wind                   [m/s]
        vg,             & ! v geostrophic wind                   [m/s]
        Lscale,         & ! Length scale                         [m]
        wpthlp2,        & ! w'thl'^2                             [m K^2/s]
        wp2thlp,        & ! w'^2 thl'                            [m^2 K/s^2]
        wprtp2,         & ! w'rt'^2                              [m/s kg^2/kg^2]
        wp2rtp,         & ! w'^2rt'                              [m^2/s^2 kg/kg]
        Lscale_up,      & ! Length scale (upwards component)     [m]
        Lscale_down,    & ! Length scale (downwards component)   [m]
        tau_zt,         & ! Eddy diss. time scale; thermo. levs. [s]
        Kh_zt,          & ! Eddy diff. coef. on thermo. levels   [m^2/s]
        wp2rcp,         & ! w'^2 rc'                             [m^2/s^2 kg/kg]
        wprtpthlp,      & ! w'rt'thl'                            [m/s kg/kg K]
        sigma_sqd_w_zt, & ! PDF width parameter (thermo. levels) [-]
        rsat              ! Saturation mixing ratio              [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
        wp2_zt,        & ! w'^2 on thermo. grid                  [m^2/s^2]
        thlp2_zt,      & ! thl'^2 on thermo. grid                [K^2]
        wpthlp_zt,     & ! w'thl' on thermo. grid                [m K/s]
        wprtp_zt,      & ! w'rt' on thermo. grid                 [m kg/(kg s)]
        rtp2_zt,       & ! rt'^2 on therm. grid                  [(kg/kg)^2]
        rtpthlp_zt,    & ! rt'thl' on thermo. grid               [kg K/kg]
        up2_zt,        & ! u'^2 on thermo. grid                  [m^2/s^2]
        vp2_zt,        & ! v'^2 on thermo. grid                  [m^2/s^2]
        upwp_zt,       & ! u'w' on thermo. grid                  [m^2/s^2]
        vpwp_zt,       & ! v'w' on thermo. grid                  [m^2/s^2]
        wpup2,         & ! w'u'^2 (thermodynamic levels)         [m^3/s^3]
        wpvp2,         & ! w'v'^2 (thermodynamic levels)         [m^3/s^3]
        wp2up2,        & ! < w'^2u'^2 > (momentum levels)        [m^4/s^4]
        wp2vp2,        & ! < w'^2v'^2 > (momentum levels)        [m^4/s^4]
        wp4,           & ! < w'^4 > (momentum levels)            [m^4/s^4]
        tau_zm,        & ! Eddy diss. time scale; momentum levs. [s]
        Kh_zm,         & ! Eddy diff. coef. on momentum levels   [m^2/s]
        thlprcp,       & ! thl'rc'                               [K kg/kg]
        rtprcp,        & ! rt'rc'                                [kg^2/kg^2]
        rcp2,          & ! rc'^2                                 [kg^2/kg^2]
        em,            & ! Turbulent Kinetic Energy (TKE)        [m^2/s^2]
        a3_coef,       & ! The a3 coefficient from CLUBB eqns    [-]
        a3_coef_zt,    & ! The a3 coef. interp. to the zt grid   [-]
        wp3_zm,        & ! w'^3 interpolated to momentum levels  [m^3/s^3]
        wp3_on_wp2,    & ! w'^3 / w'^2 on the zm grid            [m/s]
        wp3_on_wp2_zt, & ! w'^3 / w'^2 on the zt grid            [m/s]
        Skw_velocity,  & ! Skewness velocity                     [m/s]
        w_up_in_cloud    ! Upward velocity inside Clouds         [m/s]

    type(pdf_parameter), intent(in) :: & 
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: & 
      sclrm,           & ! High-order passive scalar            [units vary]
      sclrp2,          & ! High-order passive scalar variance   [units^2]
      sclrprtp,        & ! High-order passive scalar covariance [units kg/kg]
      sclrpthlp,       & ! High-order passive scalar covariance [units K]
      sclrm_forcing,   & ! Large-scale forcing of scalar        [units/s]
      sclrpthvp,       & ! High-order passive scalar covariance [units K]
      wpsclrp            ! w'sclr'                              [units m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: & 
      sclrprcp,    & ! sclr'rc'     [units vary]
      wp2sclrp,    & ! w'^2 sclr'   [units vary]
      wpsclrp2,    & ! w'sclr'^2    [units vary]
      wpsclrprtp,  & ! w'sclr'rt'   [units vary]
      wpsclrpthlp    ! w'sclr'thl'  [units vary]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: & 
      wpedsclrp,       & ! w'edsclr'                        [units vary]
      edsclrm,         & ! Eddy-diff passive scalar         [units vary] 
      edsclrm_forcing    ! Large-scale forcing of edscalar  [units vary]

    ! Local Variables

    integer :: isclr, k
    integer :: grid_level = 1  ! grid level for stats where there is only one sensible level (eg timeseries)

    real( kind = core_rknd ), dimension(gr%nz) :: &
      T_in_K,      &  ! Absolute temperature         [K]
      rsati,       &  ! Saturation w.r.t ice         [kg/kg]
      shear,       &  ! Wind shear production term   [m^2/s^3]
      chi,         &  ! Mellor's 's'                 [kg/kg]
      chip2,         &  ! Variance of Mellor's 's'     [kg/kg]
      rcm_in_cloud    ! rcm in cloud                 [kg/kg]

    real( kind = core_rknd ) :: xtmp

    ! ---- Begin Code ----

    ! Sample fields

    if ( l_stats_samp ) then

      ! stats_zt variables


      if ( iT_in_K > 0 .or. irsati > 0 ) then
        T_in_K = thlm2T_in_K( thlm, exner, rcm )
      else
        T_in_K = -999._core_rknd
      end if

      call stat_update_var( iT_in_K, T_in_K, & ! intent(in)
                            stats_zt ) ! intent(inout)
 
      call stat_update_var( ithlm, thlm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithvm, thvm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irtm, rtm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ircm, rcm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ium, um, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivm, vm, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwm_zt, wm_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwm_zm, wm_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iug, ug, & ! intent(in) 
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivg, vg, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icloud_frac, cloud_frac, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iice_supersat_frac, ice_supersat_frac, & ! intent(in)
                            stats_zt) ! intent(inout)
      call stat_update_var( ircm_in_layer, rcm_in_layer, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icloud_cover, cloud_cover, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ircm_supersat_adj, rcm_supersat_adj, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ip_in_Pa, p_in_Pa, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iexner, exner, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irho_ds_zt, rho_ds_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithv_ds_zt, thv_ds_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iLscale, Lscale, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwpup2, wpup2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwpvp2, wpvp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp3, wp3, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwpthlp2, wpthlp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp2thlp, wp2thlp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwprtp2, wprtp2, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp2rtp, wp2rtp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iLscale_up, Lscale_up, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iLscale_down, Lscale_down, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( itau_zt, tau_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iKh_zt, Kh_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp2thvp, wp2thvp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp2rcp, wp2rcp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iw_up_in_cloud, w_up_in_cloud, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwprtpthlp, wprtpthlp, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irc_coef, rc_coef, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( isigma_sqd_w_zt, sigma_sqd_w_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irho, rho, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irsat, rsat, & ! intent(in)
                            stats_zt ) ! intent(inout)
      if ( irsati > 0 ) then
        rsati = sat_mixrat_ice( p_in_Pa, T_in_K )
        call stat_update_var( irsati, rsati, & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      call stat_update_var( imixt_frac, pdf_params%mixt_frac(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iw_1, pdf_params%w_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iw_2, pdf_params%w_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_w_1, pdf_params%varnce_w_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_w_2, pdf_params%varnce_w_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithl_1, pdf_params%thl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithl_2, pdf_params%thl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_thl_1, pdf_params%varnce_thl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_thl_2, pdf_params%varnce_thl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irt_1, pdf_params%rt_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irt_2, pdf_params%rt_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_rt_1, pdf_params%varnce_rt_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivarnce_rt_2, pdf_params%varnce_rt_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout )
      call stat_update_var( irc_1, pdf_params%rc_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irc_2, pdf_params%rc_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irsatl_1, pdf_params%rsatl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irsatl_2, pdf_params%rsatl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icloud_frac_1, pdf_params%cloud_frac_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icloud_frac_2, pdf_params%cloud_frac_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ichi_1, pdf_params%chi_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ichi_2, pdf_params%chi_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( istdev_chi_1, pdf_params%stdev_chi_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( istdev_chi_2, pdf_params%stdev_chi_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( istdev_eta_1, pdf_params%stdev_eta_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( istdev_eta_2, pdf_params%stdev_eta_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icovar_chi_eta_1, pdf_params%covar_chi_eta_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icovar_chi_eta_2, pdf_params%covar_chi_eta_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_chi_1, pdf_params%corr_w_chi_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_chi_2, pdf_params%corr_w_chi_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_eta_1, pdf_params%corr_w_eta_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_eta_2, pdf_params%corr_w_eta_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_chi_eta_1, pdf_params%corr_chi_eta_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_chi_eta_2, pdf_params%corr_chi_eta_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_rt_1, pdf_params%corr_w_rt_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_rt_2, pdf_params%corr_w_rt_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_thl_1, pdf_params%corr_w_thl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_w_thl_2, pdf_params%corr_w_thl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_rt_thl_1, pdf_params%corr_rt_thl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icorr_rt_thl_2, pdf_params%corr_rt_thl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icrt_1, pdf_params%crt_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icrt_2, pdf_params%crt_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icthl_1, pdf_params%cthl_1(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( icthl_2, pdf_params%cthl_2(1,:), & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp2_zt, wp2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithlp2_zt, thlp2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ithlp3, thlp3, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwpthlp_zt, wpthlp_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwprtp_zt, wprtp_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irtp2_zt, rtp2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irtp3, rtp3, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( irtpthlp_zt, rtpthlp_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iup2_zt, up2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivp2_zt, vp2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iupwp_zt, upwp_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ivpwp_zt, vpwp_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( ia3_coef_zt, a3_coef_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)
      call stat_update_var( iwp3_on_wp2_zt, wp3_on_wp2_zt, & ! intent(in)
                            stats_zt ) ! intent(inout)

      if ( ichi > 0 ) then
        ! Determine 's' from Mellor (1977) (extended liquid water)
        chi(:) = pdf_params%mixt_frac(1,:) * pdf_params%chi_1(1,:) &
                    + (1.0_core_rknd-pdf_params%mixt_frac(1,:)) * pdf_params%chi_2(1,:)
        call stat_update_var( ichi, chi, & ! intent(in)
                             stats_zt ) ! intent(inout)
      end if 

      ! Calculate variance of chi
      if ( ichip2 > 0 ) then
        chip2 = compute_variance_binormal( chi, pdf_params%chi_1(1,:), pdf_params%chi_2(1,:), &
                                         pdf_params%stdev_chi_1(1,:), pdf_params%stdev_chi_2(1,:), &
                                         pdf_params%mixt_frac(1,:) )
        call stat_update_var( ichip2, chip2, & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( sclr_dim > 0 ) then
        do isclr=1, sclr_dim
          call stat_update_var( isclrm(isclr), sclrm(:,isclr), & ! intent(in)
                                stats_zt ) ! intent(inout)
          call stat_update_var( isclrm_f(isclr), sclrm_forcing(:,isclr),  & ! intent(in)
                                stats_zt ) ! intent(inout)
        end do
      end if

      if ( edsclr_dim > 0 ) then
        do isclr = 1, edsclr_dim
          call stat_update_var( iedsclrm(isclr), edsclrm(:,isclr), & ! intent(in)
                                stats_zt ) ! intent(inout)
          call stat_update_var( iedsclrm_f(isclr), edsclrm_forcing(:,isclr), & ! intent(in)
                                stats_zt ) ! intent(inout)
        end do
      end if

      ! Calculate rcm in cloud
      if ( ircm_in_cloud > 0 ) then
        where ( cloud_frac(:) > cloud_frac_min )
            rcm_in_cloud(:) = rcm / cloud_frac
        elsewhere
            rcm_in_cloud(:) = rcm
        endwhere

        call stat_update_var( ircm_in_cloud, rcm_in_cloud, & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      ! stats_zm variables

      call stat_update_var( iwp2, wp2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwp3_zm, wp3_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irtp2, rtp2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ithlp2, thlp2, & ! intent(in)
                            stats_zm ) ! intent(inout) 
      call stat_update_var( irtpthlp, rtpthlp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwprtp, wprtp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwpthlp, wpthlp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwp2up2, wp2up2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwp2vp2, wp2vp2, &  ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwp4, wp4, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwpthvp, wpthvp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irtpthvp, rtpthvp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ithlpthvp, thlpthvp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( itau_zm, tau_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iKh_zm, Kh_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwprcp, wprcp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irc_coef_zm, rc_coef_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ithlprcp, thlprcp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irtprcp, rtprcp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ircp2, rcp2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iupwp, upwp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ivpwp, vpwp, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ivp2, vp2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iup2, up2, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irho_zm, rho_zm, & ! intent(in)
                            stats_zm ) ! intent(inout) 
      call stat_update_var( isigma_sqd_w, sigma_sqd_w, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irho_ds_zm, rho_ds_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ithv_ds_zm, thv_ds_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iem, em, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iSkw_velocity, Skw_velocity, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ia3_coef, a3_coef, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iwp3_on_wp2, wp3_on_wp2, & ! intent(in)
                            stats_zm ) ! intent(inout)

      call stat_update_var( icloud_frac_zm, cloud_frac_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iice_supersat_frac_zm, ice_supersat_frac_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ircm_zm, rcm_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( irtm_zm, rtm_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ithlm_zm, thlm_zm, & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iw_1_zm, pdf_params_zm%w_1(1,:), & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( iw_2_zm, pdf_params_zm%w_2(1,:), & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ivarnce_w_1_zm, pdf_params_zm%varnce_w_1(1,:), & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( ivarnce_w_2_zm, pdf_params_zm%varnce_w_2(1,:), & ! intent(in)
                            stats_zm ) ! intent(inout)
      call stat_update_var( imixt_frac_zm, pdf_params_zm%mixt_frac(1,:), & ! intent(in)
                            stats_zm ) ! intent(inout)

      if ( sclr_dim > 0 ) then
        do isclr=1, sclr_dim
          call stat_update_var( isclrp2(isclr), sclrp2(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( isclrprtp(isclr), sclrprtp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( isclrpthvp(isclr), sclrpthvp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( isclrpthlp(isclr), sclrpthlp(:,isclr), & ! intent(in)
                                 stats_zm ) ! intent(inout)
          call stat_update_var( isclrprcp(isclr), sclrprcp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( iwpsclrp(isclr), wpsclrp(:,isclr), & ! intent(in)
                               stats_zm ) ! intent(inout)
          call stat_update_var( iwp2sclrp(isclr), wp2sclrp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( iwpsclrp2(isclr), wpsclrp2(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( iwpsclrprtp(isclr), wpsclrprtp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
          call stat_update_var( iwpsclrpthlp(isclr), wpsclrpthlp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
        end do
      end if
      if ( edsclr_dim > 0 ) then
        do isclr = 1, edsclr_dim
          call stat_update_var( iwpedsclrp(isclr), wpedsclrp(:,isclr), & ! intent(in)
                                stats_zm ) ! intent(inout)
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
      call stat_update_var( ishear, shear, & ! intent(in)
                            stats_zm ) ! intent(inout)

      ! stats_sfc variables

      ! Cloud cover
      call stat_update_var_pt( icc, grid_level, maxval( cloud_frac(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Cloud base
      if ( iz_cloud_base > 0 ) then

        k = 1
        do while ( rcm(k) < rc_tol .and. k < gr%nz )
          k = k + 1
        enddo

        if ( k > 1 .and. k < gr%nz) then

          ! Use linear interpolation to find the exact height of the
          ! rc_tol kg/kg level.  Brian.
          call stat_update_var_pt( iz_cloud_base, grid_level, & ! intent(in)
                                   lin_interpolate_two_points( rc_tol, rcm(k), & ! intent(in)
                                   rcm(k-1), gr%zt(k), gr%zt(k-1) ), & ! intent(in)
                                   stats_sfc ) ! intent(inout)

        else

          ! Set the cloud base output to -10m, if it's clear. 
          ! Known magic number
          call stat_update_var_pt( iz_cloud_base, grid_level, -10.0_core_rknd , & ! intent(in)
                                   stats_sfc ) ! intent(inout)
 
        end if ! k > 1 and k < gr%nz

      end if ! iz_cloud_base > 0

      ! Liquid Water Path
      if ( ilwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 rcm(2:gr%nz), gr%dzt(2:gr%nz) )

        call stat_update_var_pt( ilwp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if

      ! Vapor Water Path (Precipitable Water)
      if ( ivwp > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 ( rtm(2:gr%nz) - rcm(2:gr%nz) ), gr%dzt(2:gr%nz) )

        call stat_update_var_pt( ivwp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if


      ! Vertical average of thermodynamic level variables.

      ! Find the vertical average of thermodynamic level variables, averaged from
      ! level 2 (the first thermodynamic level above model surface) through
      ! level gr%nz (the top of the model).  Use the vertical averaging function
      ! found in fill_holes.F90.

      ! Vertical average of thlm.
      call stat_update_var_pt( ithlm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), & ! intent(in)
                         thlm(2:gr%nz), gr%dzt(2:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of rtm.
      call stat_update_var_pt( irtm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), & ! intent(in)
                         rtm(2:gr%nz), gr%dzt(2:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of um.
      call stat_update_var_pt( ium_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), & ! intent(in)
                         um(2:gr%nz), gr%dzt(2:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of vm.
      call stat_update_var_pt( ivm_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-2+1), rho_ds_zt(2:gr%nz), & ! intent(in)
                         vm(2:gr%nz), gr%dzt(2:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of momentum level variables.

      ! Find the vertical average of momentum level variables, averaged over the
      ! entire vertical profile (level 1 through level gr%nz).  Use the vertical
      ! averaging function found in fill_holes.F90.

      ! Vertical average of wp2.
      call stat_update_var_pt( iwp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), & ! intent(in)
                         wp2(1:gr%nz), gr%dzm(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of up2.
      call stat_update_var_pt( iup2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), & ! intent(in)
                         up2(1:gr%nz), gr%dzm(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of vp2.
      call stat_update_var_pt( ivp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), & ! intent(in)
                         vp2(1:gr%nz), gr%dzm(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of rtp2.
      call stat_update_var_pt( irtp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), & ! intent(in)
                         rtp2(1:gr%nz), gr%dzm(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)

      ! Vertical average of thlp2.
      call stat_update_var_pt( ithlp2_vert_avg, grid_level,  & ! intent(in)
           vertical_avg( (gr%nz-1+1), rho_ds_zm(1:gr%nz), & ! intent(in)
                         thlp2(1:gr%nz), gr%dzm(1:gr%nz) ), & ! intent(in)
                               stats_sfc ) ! intent(inout)
      
      
      if (itot_vartn_normlzd_rtm > 0) then
        if (abs(rtm(gr%nz) - rtm(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_rtm tried to divide by zero denominator ", &
                            "(surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(rtm(2 : gr%nz) - rtm(1 : gr%nz-1)) / abs(rtm(gr%nz) - rtm(1)))
        end if
        
        call stat_update_var_pt( itot_vartn_normlzd_rtm, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
     
      if (itot_vartn_normlzd_thlm > 0) then
        if (abs(thlm(gr%nz) - thlm(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_thlm tried to divide by zero denominator ", &
                            "(surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(thlm(2 : gr%nz) - thlm(1 : gr%nz-1)) / abs(thlm(gr%nz) - thlm(1)))
        end if
        
        call stat_update_var_pt( itot_vartn_normlzd_thlm, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
     
      if (itot_vartn_normlzd_wprtp > 0) then
        if (abs(wprtp(gr%nz) - wprtp(1)) < eps) then
          write(fstderr, *) "Warning: tot_vartn_normlzd_wprtp tried to divide by zero denominator ", &
                            "(surface level value was equal to top level value)"
          xtmp = -999_core_rknd  ! workaround to signify zero denominator 
        else
          xtmp = sum(abs(wprtp(2 : gr%nz) - wprtp(1 : gr%nz-1)) / abs(wprtp(gr%nz) - wprtp(1)))
        end if
        
        call stat_update_var_pt( itot_vartn_normlzd_wprtp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
      end if
    end if ! l_stats_samp


    return
  end subroutine stats_accumulate
!------------------------------------------------------------------------------
  subroutine stats_accumulate_hydromet( gr, hydromet, rho_ds_zt, & !intent(in)
                                        stats_zt, stats_sfc ) ! intent(inout)
! Description:
!   Compute stats related the hydrometeors

! References:
!   None
!------------------------------------------------------------------------------
    use parameters_model, only: &
        hydromet_dim ! Variable(s)

    use grid_class, only: &
        grid ! Type
        
    use array_index, only:  & 
        iirr, iirs, iiri, iirg, & ! Variable(s)
        iiNr, iiNs, iiNi, iiNg

    use stats_variables, only: &
        irrm, & 
        irsm, & 
        irim, & 
        irgm, & 
        iNim, & 
        iNrm, & 
        iNsm, &
        iNgm, &
        iswp, &
        irwp, &
        iiwp

    use fill_holes, only: &
        vertical_integral ! Procedure(s)

    use stats_type_utilities, only: & 
        stat_update_var, & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        l_stats_samp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_sfc

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    ! Local Variables
    real(kind=core_rknd) :: xtmp
    
    integer :: grid_level = 1

    ! ---- Begin Code ----

    if ( l_stats_samp ) then

      if ( iirr > 0 ) then
        call stat_update_var( irrm, hydromet(:,iirr), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iirs > 0 ) then
        call stat_update_var( irsm, hydromet(:,iirs), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if 

      if ( iiri > 0 ) then 
        call stat_update_var( irim, hydromet(:,iiri), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iirg > 0 ) then
        call stat_update_var( irgm,  &  ! intent(in)
                              hydromet(:,iirg), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iiNi > 0 ) then
        call stat_update_var( iNim, hydromet(:,iiNi), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iiNr > 0 ) then
        call stat_update_var( iNrm, hydromet(:,iiNr), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iiNs > 0 ) then
        call stat_update_var( iNsm, hydromet(:,iiNs), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      if ( iiNg > 0 ) then
        call stat_update_var( iNgm, hydromet(:,iiNg), & ! intent(in)
                              stats_zt ) ! intent(inout)
      end if

      ! Snow Water Path
      if ( iswp > 0 .and. iirs > 0 ) then

        ! Calculate snow water path
        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iirs), gr%dzt(2:gr%nz) )

        call stat_update_var_pt( iswp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if ! iswp > 0 .and. iirs > 0

      ! Ice Water Path
      if ( iiwp > 0 .and. iiri > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iiri), gr%dzt(2:gr%nz) )

        call stat_update_var_pt( iiwp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)

      end if

      ! Rain Water Path
      if ( irwp > 0 .and. iirr > 0 ) then

        xtmp &
        = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 hydromet(2:gr%nz,iirr), gr%dzt(2:gr%nz) )

        call stat_update_var_pt( irwp, grid_level, xtmp, & ! intent(in)
                                 stats_sfc ) ! intent(inout)
 
      end if ! irwp > 0 .and. irrm > 0
    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_hydromet
!------------------------------------------------------------------------------
  subroutine stats_accumulate_lh_tend( gr, lh_hydromet_mc, lh_Ncm_mc, &
                                       lh_thlm_mc, lh_rvm_mc, lh_rcm_mc, &
                                       lh_AKm, AKm, AKstd, AKstd_cld, &
                                       lh_rcm_avg, AKm_rcm, AKm_rcc, &
                                       stats_lh_zt )

! Description:
!   Compute stats for the tendency of latin hypercube sample points.

! References:
!   None
!------------------------------------------------------------------------------

    use parameters_model, only: &
        hydromet_dim ! Variable(s)

    use grid_class, only: &
        grid ! Type

    use array_index, only:  & 
        iirr, iirs, iiri, iirg, & ! Variable(s)
        iiNr, iiNs, iiNi, iiNg

    use stats_variables, only: &
        ilh_rrm_mc, & ! Variable(s)
        ilh_rsm_mc, & 
        ilh_rim_mc, & 
        ilh_rgm_mc, & 
        ilh_Ncm_mc, &
        ilh_Nim_mc, & 
        ilh_Nrm_mc, & 
        ilh_Nsm_mc, &
        ilh_Ngm_mc, &
        ilh_rcm_mc, &
        ilh_rvm_mc, &
        ilh_thlm_mc

    use stats_variables, only: &
        iAKstd, & ! Variable(s)
        iAKstd_cld, &
        iAKm_rcm, &
        iAKm_rcc, &
        iAKm, & 
        ilh_AKm, &
        ilh_rcm_avg

    use stats_type_utilities, only: & 
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        l_stats_samp

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_lh_zt

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      lh_hydromet_mc ! Tendency of hydrometeors except for rvm/rcm  [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      lh_Ncm_mc,  & ! Tendency of cloud droplet concentration  [num/kg/s]
      lh_thlm_mc, & ! Tendency of liquid potential temperature [kg/kg/s]
      lh_rcm_mc,  & ! Tendency of cloud water                  [kg/kg/s]
      lh_rvm_mc     ! Tendency of vapor                        [kg/kg/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      lh_AKm,     & ! Kessler ac estimate                 [kg/kg/s]
      AKm,        & ! Exact Kessler ac                    [kg/kg/s]
      AKstd,      & ! St dev of exact Kessler ac          [kg/kg/s]
      AKstd_cld,  & ! Stdev of exact w/in cloud ac        [kg/kg/s]
      lh_rcm_avg, & ! Monte Carlo rcm estimate            [kg/kg]
      AKm_rcm,    & ! Kessler ac based on rcm             [kg/kg/s]
      AKm_rcc       ! Kessler ac based on rcm/cloud_frac  [kg/kg/s]

    if ( l_stats_samp ) then

      call stat_update_var( ilh_thlm_mc, lh_thlm_mc, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)
      call stat_update_var( ilh_rcm_mc, lh_rcm_mc, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)
      call stat_update_var( ilh_rvm_mc, lh_rvm_mc, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)

      call stat_update_var( ilh_Ncm_mc, lh_Ncm_mc, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)

      if ( iirr > 0 ) then
        call stat_update_var( ilh_rrm_mc, lh_hydromet_mc(:,iirr), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iirs > 0 ) then
        call stat_update_var( ilh_rsm_mc, lh_hydromet_mc(:,iirs), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if 

      if ( iiri > 0 ) then
        call stat_update_var( ilh_rim_mc, lh_hydromet_mc(:,iiri), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iirg > 0 ) then
        call stat_update_var( ilh_rgm_mc, lh_hydromet_mc(:,iirg), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iiNi > 0 ) then
        call stat_update_var( ilh_Nim_mc, lh_hydromet_mc(:,iiNi), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iiNr > 0 ) then
        call stat_update_var( ilh_Nrm_mc, lh_hydromet_mc(:,iiNr), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iiNs > 0 ) then
        call stat_update_var( ilh_Nsm_mc, lh_hydromet_mc(:,iiNs), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if

      if ( iiNg > 0 ) then
        call stat_update_var( ilh_Ngm_mc, lh_hydromet_mc(:,iiNg), & ! intent(in)
                              stats_lh_zt ) ! intent(inout)
      end if 

      call stat_update_var( iAKm, AKm, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)
      call stat_update_var( ilh_AKm, lh_AKm, & ! intent(in)
                            stats_lh_zt) ! intent(inout)
      call stat_update_var( ilh_rcm_avg, lh_rcm_avg, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)
      call stat_update_var( iAKstd, AKstd, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)
      call stat_update_var( iAKstd_cld, AKstd_cld, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)

      call stat_update_var( iAKm_rcm, AKm_rcm, & ! intent(in)
                            stats_lh_zt) ! intent(inout)
      call stat_update_var( iAKm_rcc, AKm_rcc, & ! intent(in)
                            stats_lh_zt ) ! intent(inout)

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_lh_tend
    
  !-----------------------------------------------------------------------
  subroutine stats_finalize( stats_zt, stats_zm, stats_sfc, &
                             stats_lh_zt, stats_lh_sfc, &
                             stats_rad_zt, stats_rad_zm ) ! intent(inout)

    !     Description:
    !     Close NetCDF files and deallocate scratch space and
    !     stats file structures.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        l_netcdf, & 
        l_stats, &
        l_output_rad_files, &
        l_silhs_out

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

    use stats_variables, only: &
        ihm_1, &
        ihm_2, &
        imu_hm_1, &
        imu_hm_2, &
        imu_hm_1_n, &
        imu_hm_2_n, &
        isigma_hm_1, &
        isigma_hm_2, &
        isigma_hm_1_n, &
        isigma_hm_2_n, &
        icorr_w_hm_1, &
        icorr_w_hm_2, &
        icorr_chi_hm_1, &
        icorr_chi_hm_2, &
        icorr_eta_hm_1, &
        icorr_eta_hm_2, &
        icorr_Ncn_hm_1, &
        icorr_Ncn_hm_2, &
        icorr_hmx_hmy_1, &
        icorr_hmx_hmy_2, &
        icorr_w_hm_1_n, &
        icorr_w_hm_2_n, &
        icorr_chi_hm_1_n, &
        icorr_chi_hm_2_n, &
        icorr_eta_hm_1_n, &
        icorr_eta_hm_2_n, &
        icorr_Ncn_hm_1_n, &
        icorr_Ncn_hm_2_n, &
        icorr_hmx_hmy_1_n, &
        icorr_hmx_hmy_2_n, &
        ihmp2_zt, &
        iwp2hmp, &
        ihydrometp2, &
        iwphydrometp, &
        iK_hm, &
        irtphmp, &
        ithlphmp, &
        ihmxphmyp

    use stats_variables, only: &
        isilhs_variance_category, & ! Variable(s)
        ilh_samp_frac_category

#ifdef NETCDF
    use output_netcdf, only:  & 
        close_netcdf ! Procedure
#endif

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    if ( l_stats .and. l_netcdf ) then
#ifdef NETCDF
      call close_netcdf( stats_zt%file ) ! intent(inout)
      call close_netcdf( stats_lh_zt%file ) ! intent(inout)
      call close_netcdf( stats_lh_sfc%file ) ! intent(inout)
      call close_netcdf( stats_zm%file ) ! intent(inout)
      call close_netcdf( stats_rad_zt%file ) ! intent(inout)
      call close_netcdf( stats_rad_zm%file ) ! intent(inout)
      call close_netcdf( stats_sfc%file ) ! intent(inout)
#else
      error stop "This program was not compiled with netCDF support"
#endif
    end if

    if ( l_stats ) then
      ! De-allocate all stats_zt variables
      if (allocated(stats_zt%z)) then
        deallocate( stats_zt%z )

        deallocate( stats_zt%accum_field_values )
        deallocate( stats_zt%accum_num_samples )
        deallocate( stats_zt%l_in_update )

        deallocate( stats_zt%file%grid_avg_var )
        deallocate( stats_zt%file%z )
               
        ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
        if ( allocated( stats_zt%file%lat_vals ) ) then
          deallocate( stats_zt%file%lat_vals )
        end if
        if ( allocated( stats_zt%file%lon_vals ) ) then
          deallocate( stats_zt%file%lon_vals )
        end if

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
      end if

      if ( l_silhs_out .and. allocated(stats_lh_zt%z) ) then
        ! De-allocate all stats_lh_zt variables
        deallocate( stats_lh_zt%z )

        deallocate( stats_lh_zt%accum_field_values )
        deallocate( stats_lh_zt%accum_num_samples )
        deallocate( stats_lh_zt%l_in_update )

        deallocate( stats_lh_zt%file%grid_avg_var )
        deallocate( stats_lh_zt%file%z )
        
        ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
        if ( allocated(stats_lh_zt%file%lat_vals ) ) then
          deallocate( stats_lh_zt%file%lat_vals )
        end if
        if ( allocated(stats_lh_zt%file%lon_vals ) ) then
          deallocate( stats_lh_zt%file%lon_vals )
        end if

        ! De-allocate all stats_lh_sfc variables
        deallocate( stats_lh_sfc%z )

        deallocate( stats_lh_sfc%accum_field_values )
        deallocate( stats_lh_sfc%accum_num_samples )
        deallocate( stats_lh_sfc%l_in_update )

        deallocate( stats_lh_sfc%file%grid_avg_var )
        deallocate( stats_lh_sfc%file%z )
             
        ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
        if ( allocated( stats_lh_sfc%file%lat_vals ) ) then
          deallocate( stats_lh_sfc%file%lat_vals )
        end if
        if ( allocated( stats_lh_sfc%file%lon_vals ) ) then
          deallocate( stats_lh_sfc%file%lon_vals )
        end if
      end if ! l_silhs_out

      ! De-allocate all stats_zm variables
      if (allocated(stats_zm%z)) then
        deallocate( stats_zm%z )

        deallocate( stats_zm%accum_field_values )
        deallocate( stats_zm%accum_num_samples )
        deallocate( stats_zm%l_in_update )

        deallocate( stats_zm%file%grid_avg_var )
        deallocate( stats_zm%file%z )
             
        ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
        if ( allocated( stats_zm%file%lat_vals ) ) then
          deallocate( stats_zm%file%lat_vals )
        end if
        if ( allocated( stats_zm%file%lon_vals ) ) then
          deallocate( stats_zm%file%lon_vals )
        end if

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
      end if

      if ( l_output_rad_files ) then
        ! De-allocate all stats_rad_zt variables
        if (allocated(stats_rad_zt%z)) then
          deallocate( stats_rad_zt%z )

          deallocate( stats_rad_zt%accum_field_values )
          deallocate( stats_rad_zt%accum_num_samples )
          deallocate( stats_rad_zt%l_in_update )

          deallocate( stats_rad_zt%file%grid_avg_var )
          deallocate( stats_rad_zt%file%z )
               
          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated( stats_rad_zt%file%lat_vals ) ) then
            deallocate( stats_rad_zt%file%lat_vals )
          end if
          if ( allocated( stats_rad_zt%file%lon_vals ) ) then
            deallocate( stats_rad_zt%file%lon_vals )
          end if

          ! De-allocate all stats_rad_zm variables
          deallocate( stats_rad_zm%z )

          deallocate( stats_rad_zm%accum_field_values )
          deallocate( stats_rad_zm%accum_num_samples )
          deallocate( stats_rad_zm%l_in_update )

          deallocate( stats_rad_zm%file%grid_avg_var )
          deallocate( stats_rad_zm%file%z )

          ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
          if ( allocated( stats_rad_zm%file%lat_vals ) ) then
            deallocate( stats_rad_zm%file%lat_vals )
          end if
          if ( allocated( stats_rad_zm%file%lon_vals ) ) then
            deallocate( stats_rad_zm%file%lon_vals )
          end if

        end if

      end if ! l_output_rad_files

      ! De-allocate all stats_sfc variables
      if (allocated(stats_sfc%z)) then
        deallocate( stats_sfc%z )

        deallocate( stats_sfc%accum_field_values )
        deallocate( stats_sfc%accum_num_samples )
        deallocate( stats_sfc%l_in_update )

        deallocate( stats_sfc%file%grid_avg_var )
        deallocate( stats_sfc%file%z )
      end if
             
      ! Check if pointer is allocated to prevent crash in netcdf (ticket 765)
      if ( allocated( stats_sfc%file%lat_vals ) ) then
        deallocate( stats_sfc%file%lat_vals )
      end if
      if ( allocated( stats_sfc%file%lon_vals ) ) then
        deallocate( stats_sfc%file%lon_vals )
      end if

      ! De-allocate scalar indices
      if (allocated(isclrm)) then
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
      end if

      ! De-allocate hyderometeor statistical variables
      if (allocated(ihm_1)) then
        deallocate( ihm_1 )
        deallocate( ihm_2 )
        deallocate( imu_hm_1 )
        deallocate( imu_hm_2 )
        deallocate( imu_hm_1_n )
        deallocate( imu_hm_2_n )
        deallocate( isigma_hm_1 )
        deallocate( isigma_hm_2 )
        deallocate( isigma_hm_1_n )
        deallocate( isigma_hm_2_n )
        deallocate( icorr_w_hm_1 )
        deallocate( icorr_w_hm_2 )
        deallocate( icorr_chi_hm_1 )
        deallocate( icorr_chi_hm_2 )
        deallocate( icorr_eta_hm_1 )
        deallocate( icorr_eta_hm_2 )
        deallocate( icorr_Ncn_hm_1 )
        deallocate( icorr_Ncn_hm_2 )
        deallocate( icorr_hmx_hmy_1 )
        deallocate( icorr_hmx_hmy_2 )
        deallocate( icorr_w_hm_1_n )
        deallocate( icorr_w_hm_2_n )
        deallocate( icorr_chi_hm_1_n )
        deallocate( icorr_chi_hm_2_n )
        deallocate( icorr_eta_hm_1_n )
        deallocate( icorr_eta_hm_2_n )
        deallocate( icorr_Ncn_hm_1_n )
        deallocate( icorr_Ncn_hm_2_n )
        deallocate( icorr_hmx_hmy_1_n )
        deallocate( icorr_hmx_hmy_2_n )
        deallocate( ihmp2_zt )
        deallocate( iwp2hmp )
        deallocate( ihydrometp2 )
        deallocate( iwphydrometp )
        deallocate( irtphmp )
        deallocate( ithlphmp )
        deallocate( ihmxphmyp )
        deallocate( iK_hm )
      end if

      if ( allocated( isilhs_variance_category ) ) then
        deallocate( isilhs_variance_category )
        deallocate( ilh_samp_frac_category )
      end if

    end if ! l_stats

    return
  end subroutine stats_finalize

!===============================================================================

!-----------------------------------------------------------------------
subroutine stats_check_num_samples( stats_grid, &
                                    l_error )

! Description:
!   Ensures that each variable in a stats grid is sampled the correct
!   number of times.
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_tsamp, & ! Variable(s)
        stats_tout

    use error_code, only: &
        clubb_at_least_debug_level   ! Procedure

    implicit none

  ! Input Variables
    type (stats), intent(in) :: &
      stats_grid               ! Grid type              [grid]

  ! Input/Output Variables
    logical, intent(inout) :: &
      l_error                  ! Indicates an error     [boolean]

  ! Local Variables
    integer :: ivar, kvar      ! Loop variable          [index]

    logical :: l_proper_sample

!-----------------------------------------------------------------------

  !----- Begin Code -----

  ! Look for errors by checking the number of sampling points
  ! for each variable in the statistics grid at each vertical level.
  do ivar = 1, stats_grid%num_output_fields
    do kvar = 1, stats_grid%kk

      l_proper_sample = ( stats_grid%accum_num_samples(1,1,kvar,ivar) == 0 .or. &
                          stats_grid%accum_num_samples(1,1,kvar,ivar) == &
                            floor(stats_tout/stats_tsamp) )

      if ( .not. l_proper_sample ) then

        l_error = .true.  ! This will stop the run

        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) 'Possible sampling error for variable ',  &
                           trim(stats_grid%file%grid_avg_var(ivar)%name), ' in stats_grid ',  &
                           'at k = ', kvar,  &
                           '; stats_grid%accum_num_samples(',kvar,',',ivar,') = ', &
                            stats_grid%accum_num_samples(1,1,kvar,ivar)
        end if ! clubb_at_lest_debug_level 1


      end if ! .not. l_proper_sample

    end do ! kvar = 1 .. stats_grid%kk
  end do ! ivar = 1 .. stats_grid%num_output_fields

  return
end subroutine stats_check_num_samples
!-----------------------------------------------------------------------

end module stats_clubb_utilities
