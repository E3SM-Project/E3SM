! $Id: stat_clubb.F90 1070 2013-04-19 20:05:10Z minghuai.wang@pnl.gov $
module stat_clubb
#ifdef CLUBB_CRM
  use grid, only: nx, ny, nz, nzm
  implicit none

  public :: stats_clubb_update

#ifdef CLUBB_LH
  public stats_clubb_silhs_update
#endif

  public :: stats_end_timestep_clubb, stats_init_clubb
#ifndef CRM
  public ::  hbuf_stats_init_clubb
#endif

  !  Output arrays for CLUBB statistics    
  real, allocatable, dimension(:,:,:,:) :: out_zt, out_zm, out_rad_zt, out_rad_zm,  &
        out_sfc, out_LH_zt, out_LH_sfc

  private

  contains 
!---------------------------------------------------------------------------------------------------
  subroutine stats_clubb_update( upwp, vpwp, up2, vp2, wprtp, wpthlp, &
    wp2, wp3, rtp2, thlp2, rtpthlp, cloud_frac, rcm, um, vm, t_tndcy, & 
    qc_tndcy, qv_tndcy,u_tndcy,v_tndcy )

! Description:
!   Update statistics for CLUBB variables
!
! References:
!   None
!---------------------------------------------------------------------------------------------------
  use grid, only: nx, ny, nzm, nz, dimx1_s, dimx2_s, dimy1_s, dimy2_s

#ifndef CRM
  use hbuffer, only: hbuf_put, hbuf_avg_put
#endif

  ! Modules from CLUBB
  use clubb_precision, only: core_rknd ! Constant

  use interpolation, only: lin_int ! Procedure(s)

  use grid_class, only: gr

  use clubbvars, only: tndcy_precision, l_stats_samgrid

  implicit none

  real(kind=core_rknd), dimension(nx, ny, nz), intent(in) :: &
    upwp,        &! u'w'                          [m^2/s^2]
    vpwp,        &! u'w'                          [m^2/s^2]
    up2,         &! u'^2                          [m^2/s^2]
    vp2,         &! v'^2                          [m^2/s^2]
    wprtp,       &! w' r_t'                       [(m kg)/(s kg)]
    wpthlp,      &! w' th_l'                      [(m K)/s]
    wp2,         &! w'^2                          [m^2/s^2]
    rtp2,        &! r_t'^2                        [(kg/kg)^2]
    thlp2,       &! th_l'^2                       [K^2]
    rtpthlp,     &! r_t' th_l'                    [(kg K)/kg]
    cloud_frac,  &! Cloud Fraction                [-]
    rcm           ! Cloud water                   [kg/kg]

  ! w'^3 is requires additional ghost points on the x and y dimension
  real(kind=core_rknd), dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nz), intent(in) :: &
    wp3,&    ! w'^3                       [m^3/s^3]
    um, &    ! x-wind                     [m/s]
    vm       ! y-wind                     [m/s]

  real(tndcy_precision), dimension(nx, ny, nzm), intent(in) :: &
    t_tndcy,  & ! CLUBB contribution to moist static energy  [K/s]
    qc_tndcy, & ! CLUBB contribution to liquid water         [kg/kg/s]
    qv_tndcy, & ! CLUBB contribution to vapor water          [kg/kg/s]
    u_tndcy,  & ! CLUBB contribution to x-wind               [m/s^2]
    v_tndcy     ! CLUBB contribution to y-wind               [m/s^2]

  ! Local variables
  real, dimension(nzm) :: &
    upwp_avg,   &
    vpwp_avg,   &
    up2_avg,    &
    vp2_avg,    &
    wprtp_avg,  &
    wpthlp_avg, &
    wp2_avg,    &
    thlp2_avg,  &
    rtp2_avg,   &
    rtpthlp_avg,&
    sigma_sqd_w_avg, &
    Kh_zt_avg,  &
    tau_zm_avg

  real :: factor_xy

  integer :: i, j, k

  !---------------------------------------------------------
  ! CLUBB variables
  ! Notes: The variables located on the vertical velocity levels 
  ! must be interpolated for the stats grid, which is on the pressure levels.
  ! -dschanen 21 Jul 2008
  factor_xy = 1. / real( nx*ny )

  upwp_avg   = 0.0
  vpwp_avg   = 0.0
  vp2_avg    = 0.0
  up2_avg    = 0.0
  wprtp_avg  = 0.0
  wpthlp_avg = 0.0
  wp2_avg    = 0.0

  thlp2_avg   = 0.0
  rtp2_avg    = 0.0
  rtpthlp_avg = 0.0

  ! Here we omit the ghost point, since the SAM stats don't have one
  do i = 1, nx
    do j = 1, ny
      do k = 1, nzm
        upwp_avg(k) = upwp_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), upwp(i,j,k+1), upwp(i,j,k) )
        vpwp_avg(k) = vpwp_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), vpwp(i,j,k+1), vpwp(i,j,k) )
        vp2_avg(k) = vp2_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), vp2(i,j,k+1), vp2(i,j,k) )
        up2_avg(k) = up2_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), up2(i,j,k+1), up2(i,j,k) )
        wprtp_avg(k) = wprtp_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), wprtp(i,j,k+1), wprtp(i,j,k) )
        wpthlp_avg(k) = wpthlp_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), wpthlp(i,j,k+1), wpthlp(i,j,k) )
        wp2_avg(k) = wp2_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), wp2(i,j,k+1), wp2(i,j,k) )
        rtp2_avg(k) = rtp2_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), rtp2(i,j,k+1), rtp2(i,j,k) )
        thlp2_avg(k) = thlp2_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), thlp2(i,j,k+1), thlp2(i,j,k) )
        rtpthlp_avg(k) = rtpthlp_avg(k) &
          + lin_int( gr%zt(k+1), gr%zm(k+1), gr%zm(k), rtpthlp(i,j,k+1), rtpthlp(i,j,k) )
      end do ! k = 1..nzm
    end do ! j = 1..ny
  end do ! i = 1..nx

#ifndef CRM
  ! Velocity grid variables
  call hbuf_put('UPWP', upwp_avg, factor_xy)
  call hbuf_put('VPWP', vpwp_avg, factor_xy)
  call hbuf_put('VP2', vp2_avg, factor_xy)
  call hbuf_put('UP2', up2_avg, factor_xy)
  call hbuf_put('WPRTP', wprtp_avg, factor_xy)
  call hbuf_put('WPTHLP', wpthlp_avg, factor_xy)
  call hbuf_put('WP2', wp2_avg, factor_xy)
  call hbuf_put('RTP2', rtp2_avg, factor_xy)
  call hbuf_put('THLP2', thlp2_avg, factor_xy)
  call hbuf_put('RTPTHLP', rtpthlp_avg, factor_xy)

  ! CLUBB thermodynamic grid varibles (SAM pressure levels + ghost point)
  call hbuf_avg_put('CLD_FRAC', real( cloud_frac(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('RCM', real( rcm(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('UM', real( um(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('VM', real( vm(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('WP3', real( wp3(1:nx,1:ny,2:nz) ), 1,nx, 1,ny, nzm, 1.)

  ! CLUBB tendency of state variables
  call hbuf_avg_put('T_TNDCY', real(t_tndcy(1:nx,1:ny,1:nzm)), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('QC_TNDCY', real(qc_tndcy(1:nx,1:ny,1:nzm)), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('QV_TNDCY', real(qv_tndcy(1:nx,1:ny,1:nzm)), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('U_TNDCY', real(U_tndcy(1:nx,1:ny,1:nzm)), 1,nx, 1,ny, nzm, 1.)
  call hbuf_avg_put('V_TNDCY', real(V_tndcy(1:nx,1:ny,1:nzm)), 1,nx, 1,ny, nzm, 1.)

  if(l_stats_samgrid) then  !output clubb statistics in SAM
    call hbuf_clubb_output ()
  end if
#endif

  return
  end subroutine stats_clubb_update

#ifdef CLUBB_LH
!---------------------------------------------------------------------------------------------------
  subroutine stats_clubb_silhs_update( )

! Description:
!   Update statistics for CLUBB SILHS variables
!
! References:
!   None
!---------------------------------------------------------------------------------------------------
    use grid, only: nx, ny, nzm, nz

    use hbuffer, only: hbuf_put, hbuf_avg_put

    use microphysics, only: &
      nmicro_fields, mkname, index_water_vapor

    ! Modules from CLUBB
    use clubb_precision, only: core_rknd ! Constant

    use interpolation, only: lin_int ! Procedure(s)

    use grid_class, only: gr

    use clubb_silhs_vars, only: &
      LH_rt, LH_t, X_nl_all_levs, LH_sample_point_weights, LH_t_avg_tndcy, &
      LH_micro_field_avg_tndcy

    use latin_hypercube_arrays, only: &
      d_variables

    use parameters_microphys, only: &
      LH_microphys_calls

    use corr_matrix_module, only: &
      iiLH_s_mellor, iiLH_w, &
      iiLH_rrain, iiLH_rsnow, iiLH_rice, &
      iiLH_Nr, iiLH_Nsnow, iiLH_Ni, iiLH_Nc

    use array_index, only: & 
      iirrainm, iiNrm, iirsnowm, iiricem, & ! Variables
      iiNcm, iiNsnowm, iiNim

    implicit none

    ! Local Variables
    real, dimension(nx,ny,nzm) :: &
      LH_rt_weighted, &
      LH_t_weighted

    real, dimension(nx,ny,nzm,d_variables) :: &
      X_nl_all_levs_weighted

    character(len=8) :: stat_name
    integer :: indx, ivar, k

    ! ---- Begin Code ----

    ! Determine cloud weighted sample averages
    LH_rt_weighted   = 0.
    LH_t_weighted    = 0.
    X_nl_all_levs_weighted = 0.

    do indx = 1, LH_microphys_calls
      do k = 1, nzm
        LH_rt_weighted(:,:,k) = LH_rt_weighted(:,:,k) &
          + LH_rt(:,:,k,indx) * LH_sample_point_weights(:,:,indx)
        LH_t_weighted(:,:,k)  = LH_t_weighted(:,:,k) &
          + LH_t(:,:,k,indx) * LH_sample_point_weights(:,:,indx)

        do ivar = 1, d_variables
          X_nl_all_levs_weighted(:,:,k,ivar) = X_nl_all_levs_weighted(:,:,k,ivar) &
            + X_nl_all_levs(:,:,k,indx,ivar) * LH_sample_point_weights(:,:,indx)
        end do

      end do ! k = 1..nzm
    end do ! indx = 1..LH_microphys_calls

    LH_rt_weighted = LH_rt_weighted / real( LH_microphys_calls )
    LH_t_weighted = LH_t_weighted / real( LH_microphys_calls )
    X_nl_all_levs_weighted = X_nl_all_levs_weighted / real( LH_microphys_calls )

    call hbuf_avg_put( 'LH_RT', LH_rt_weighted, 1,nx, 1,ny, nzm, 1. )
    call hbuf_avg_put( 'LH_TL',  LH_t_weighted, 1,nx, 1,ny, nzm, 1. )

    do ivar = 1, d_variables
      if ( ivar  == iiLH_s_mellor ) then
        stat_name = "LH_S_MEL"
      else if ( ivar == iiLH_w ) then
        stat_name = "LH_W"
      else if ( ivar == iiLH_rrain ) then
        stat_name = "LH_RRAIN"
      else if ( ivar == iiLH_rsnow ) then
        stat_name = "LH_RSNOW"
      else if ( ivar == iiLH_rice ) then
        stat_name = "LH_RICE"
      else if ( ivar == iiLH_Nr ) then
        stat_name = "LH_NR"
      else if ( ivar == iiLH_Nsnow ) then
        stat_name = "LH_NSNOW"
      else if ( ivar == iiLH_Ni ) then
        stat_name = "LH_NI"
      else if ( ivar == iiLH_Nc ) then
        stat_name = "LH_NC"
      end if ! ivar

      call hbuf_avg_put( stat_name, X_nl_all_levs_weighted(:,:,:,ivar), 1,nx, 1,ny, nzm, 1. )
    end do

    ! Tendency averages

    call hbuf_avg_put( 'LH_TL_MC', real( LH_t_avg_tndcy ), &
                       1,nx, 1,ny, nzm, 1. )

    do ivar = 1, nmicro_fields
      if ( ivar == index_water_vapor ) then
        stat_name = 'LH_RT_MC'
      else if ( ivar == iirrainm ) then
        stat_name = 'LH_RR_MC'
      else if ( ivar == iirsnowm ) then
        stat_name = 'LH_RS_MC'
      else if ( ivar == iiricem ) then
        stat_name = 'LH_RI_MC'
      else if ( ivar == iiNim ) then
        stat_name = 'LH_NI_MC'
      else if ( ivar == iiNrm ) then
        stat_name = 'LH_NR_MC'
      else if ( ivar == iiNsnowm ) then
        stat_name = 'LH_NS_MC'
      else
        stat_name = ''
      end if
      if ( stat_name /= '' ) then
        call hbuf_avg_put( stat_name, &
                           real( LH_micro_field_avg_tndcy(:,:,:,ivar) ), &
                           1,nx, 1,ny, nzm, 1. )
      end if
    end do

    return
  end subroutine stats_clubb_silhs_update
#endif /* CLUBB_LH */

subroutine stats_init_clubb( l_stats_in, l_output_rad_files_in, stats_tsamp_in, stats_tout_in, &
                         nzmax, nnrad_zt,nnrad_zm, time_current, delt )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of SAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with SAM output. This is adopted from clubb_intr.F90 in CAM5.2.   
    
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
      fname_rad_zt, &
      fname_rad_zm, & 
      fname_sfc, & 
      l_netcdf, & 
      l_grads

    use clubb_precision, only: & 
      time_precision, &   ! Constant(s)
      core_rknd

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

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    logical, intent(in) :: l_output_rad_files_in ! Rad Stats on? T/F

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nzmax ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count] 
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Model time                         [s]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Timestep (dt_main in CLUBB)         [s]


    !  Local Variables

    !  Namelist Variables

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      clubb_vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_LH_zt) ::  & 
      clubb_vars_LH_zt  ! Latin Hypercube variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_LH_sfc) ::  & 
      clubb_vars_LH_sfc  ! Latin Hypercube variables at the surface

    character(len=var_length), dimension(nvarmax_zm) ::  & 
      clubb_vars_zm  ! Variables on the momentum levels

    character(len=var_length), dimension(nvarmax_rad_zt) ::  & 
      clubb_vars_rad_zt  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_rad_zm) ::  & 
      clubb_vars_rad_zm  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_sfc) ::  &
      clubb_vars_sfc ! Variables at the model surface

    namelist /clubb_stats_nl/ & 
      clubb_vars_zt, & 
      clubb_vars_zm, &
      clubb_vars_LH_zt, &
      clubb_vars_LH_sfc, &
      clubb_vars_rad_zt, &
      clubb_vars_rad_zm, & 
      clubb_vars_sfc

    !  Local Variables

    logical :: l_error

    character(len=200) :: fname, temp1, sub

    integer :: i, ntot, read_status
    integer :: iunit

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    l_output_rad_files = l_output_rad_files_in
    
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_last  = .false.
      return
    end if

    !  Initialize namelist variables

    clubb_vars_zt     = ''
    clubb_vars_zm     = ''
    clubb_vars_LH_zt = ''
    clubb_vars_LH_sfc = ''
    clubb_vars_rad_zt = ''
    clubb_vars_rad_zm = ''
    clubb_vars_sfc    = ''

    !  Read variables to compute from the namelist    
    ! in SAM, namelist is read on every MPI task, so no need for mpibcast
!    if (masterproc) then
      iunit= 55 
      open(unit=iunit,file="clubb_stats_sam")
      read(unit=iunit, nml=clubb_stats_nl, iostat=read_status)
      if (read_status /= 0) then
        stop 'stats_init_clubb:  error reading namelist'
      end if
      close(unit=iunit)
!    end if

!#ifdef SPMD
      ! Broadcast namelist variables
!      call mpibcast(clubb_vars_zt,      var_length*nvarmax_zt,       mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_zm,      var_length*nvarmax_zm,       mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_LH_zt,      var_length*nvarmax_LH_zt,       mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_LH_sfc,      var_length*nvarmax_LH_sfc,       mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpichar,   0, mpicom)
!      call mpibcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpichar,   0, mpicom)
!#endif

    !  Hardcode these for use in SAM-CLUBB, don't want either
    l_netcdf = .false.
    l_grads  = .false.

    !  Check sampling and output frequencies

    !  The model time step length, delt (which is dtmain), should multiply
    !  evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - real(floor(stats_tsamp/delt), kind=time_precision ) )  & 
           > 1.e-8_time_precision ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                       'delt (which is dtmain).  Check the appropriate ',  &
                       'model.in file.'
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
      write(fstderr,*) 'delt = ', delt
    endif

    !  Initialize zt (mass points)

    i = 1
    do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_zt(i)) /= 0 & 
               .and. i <= nvarmax_zt )
       i = i + 1
       write(2001, *) 'i=', i-1, ' clubb_vars_zt  ', trim(clubb_vars_zt(i))
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      stop "stats_init_clubb:  number of zt statistical variables exceeds limit"
    endif

    zt%nn = ntot
    zt%kk = nzmax

    allocate( zt%z( zt%kk ) )

    allocate( zt%x( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%n( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%l_in_update( 1, 1, zt%kk, zt%nn ) )
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )

    allocate( zt%f%var( zt%nn ) )
    allocate( zt%f%z( zt%kk ) )

    !  Allocate scratch space

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

    !  Default initialization for array indices for zt

    call stats_init_zt( clubb_vars_zt, l_error )

    ! Setup output file for LH_zt (Latin Hypercube stats)

    if ( LH_microphys_type /= LH_microphys_disabled ) then

      i = 1
      do while ( ichar(clubb_vars_LH_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_LH_zt(i)) /= 0 & 
                 .and. i <= nvarmax_LH_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_LH_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_zt than allowed for by nvarmax_LH_zt."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_LH_zt ",  &
                         "in the stats namelist, or change nvarmax_LH_zt."
        write(fstderr,*) "nvarmax_LH_zt = ", nvarmax_LH_zt
        stop "stats_init:  number of LH_zt statistical variables exceeds limit"
      end if

      LH_zt%nn = ntot
      LH_zt%kk = nzmax

      allocate( LH_zt%z( LH_zt%kk ) )
!      LH_zt%z = gzt

      allocate( LH_zt%x( 1, 1, LH_zt%kk, LH_zt%nn ) )
      allocate( LH_zt%n( 1, 1, LH_zt%kk, LH_zt%nn ) )
      allocate( LH_zt%l_in_update( 1, 1, LH_zt%kk, LH_zt%nn ) )
      call stats_zero( LH_zt%kk, LH_zt%nn, LH_zt%x, LH_zt%n, LH_zt%l_in_update )

      allocate( LH_zt%f%var( LH_zt%nn ) )
      allocate( LH_zt%f%z( LH_zt%kk ) )

      call stats_init_LH_zt( clubb_vars_LH_zt, l_error )

      i = 1
      do while ( ichar(clubb_vars_LH_sfc(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_LH_sfc(i)) /= 0 & 
                 .and. i <= nvarmax_LH_sfc )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_LH_sfc ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_zt than allowed for by nvarmax_LH_sfc."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_LH_sfc ",  &
                         "in the stats namelist, or change nvarmax_LH_sfc."
        write(fstderr,*) "nvarmax_LH_sfc = ", nvarmax_LH_sfc
        stop "stats_init:  number of LH_sfc statistical variables exceeds limit"
      end if

      LH_sfc%nn = ntot
      LH_sfc%kk = 1

      allocate( LH_sfc%z( LH_sfc%kk ) )

      allocate( LH_sfc%x( 1, 1, LH_sfc%kk, LH_sfc%nn ) )
      allocate( LH_sfc%n( 1, 1, LH_sfc%kk, LH_sfc%nn ) )
      allocate( LH_sfc%l_in_update( 1, 1, LH_sfc%kk, LH_sfc%nn ) )

      call stats_zero( LH_sfc%kk, LH_sfc%nn, LH_sfc%x, LH_sfc%n, LH_sfc%l_in_update )

      allocate( LH_sfc%f%var( LH_sfc%nn ) )
      allocate( LH_sfc%f%z( LH_sfc%kk ) )

      call stats_init_LH_sfc( clubb_vars_LH_sfc, l_error )

    end if ! LH_microphys_type /= LH_microphys_disabled

    !  Initialize zm (momentum points)

    i = 1
    do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_zm(i)) /= 0 & 
               .and. i <= nvarmax_zm )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      stop "stats_init_clubb:  number of zm statistical variables exceeds limit"
    endif

    zm%nn = ntot
    zm%kk = nzmax

    allocate( zm%z( zm%kk ) )

    allocate( zm%x( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%n( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%l_in_update( 1, 1, zm%kk, zm%nn ) )

    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )

    allocate( zm%f%var( zm%nn ) )
    allocate( zm%f%z( zm%kk ) )

    !  Allocate scratch space

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

    call stats_init_zm( clubb_vars_zm, l_error )

    !  Initialize rad_zt (radiation points)

    if (l_output_rad_files) then
    
      i = 1
      do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_rad_zt(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        stop "stats_init_clubb:  number of rad_zt statistical variables exceeds limit"
      endif

      rad_zt%nn = ntot
      rad_zt%kk = nnrad_zt

      allocate( rad_zt%z( rad_zt%kk ) )

      allocate( rad_zt%x( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%n( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%l_in_update( 1, 1, rad_zt%kk, rad_zt%nn ) )

      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )

      allocate( rad_zt%f%var( rad_zt%nn ) )
      allocate( rad_zt%f%z( rad_zt%kk ) )


      call stats_init_rad_zt( clubb_vars_rad_zt, l_error )

      !  Initialize rad_zm (radiation points)

      i = 1
      do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0  & 
                 .and. len_trim(clubb_vars_rad_zm(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zm )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        stop "stats_init_clubb:  number of rad_zm statistical variables exceeds limit"
      endif

      rad_zm%nn = ntot
      rad_zm%kk = nnrad_zm

      allocate( rad_zm%z( rad_zm%kk ) )

      allocate( rad_zm%x( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%n( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%l_in_update( 1, 1, rad_zm%kk, rad_zm%nn ) )

      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )

      allocate( rad_zm%f%var( rad_zm%nn ) )
      allocate( rad_zm%f%z( rad_zm%kk ) )


      call stats_init_rad_zm( clubb_vars_rad_zm, l_error )
    end if ! l_output_rad_files


    !  Initialize sfc (surface point)

    i = 1
    do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0  & 
               .and. len_trim(clubb_vars_sfc(i)) /= 0 & 
               .and. i <= nvarmax_sfc )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "clubb_vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      stop "stats_init_clubb:  number of sfc statistical variables exceeds limit"
    endif

    sfc%nn = ntot
    sfc%kk = 1

    allocate( sfc%z( sfc%kk ) )

    allocate( sfc%x( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%n( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%l_in_update( 1, 1, sfc%kk, sfc%nn ) )

    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    allocate( sfc%f%var( sfc%nn ) )
    allocate( sfc%f%z( sfc%kk ) )

    call stats_init_sfc( clubb_vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      stop
    endif

    allocate(out_zt(nx, ny, nz, zt%nn))
    allocate(out_zm(nx, ny, nz, zm%nn))
    allocate(out_sfc(nx, ny, nz, sfc%nn))

    if(l_output_rad_files) then
      allocate(out_rad_zt(nx, ny, nz, rad_zt%nn))
      allocate(out_rad_zm(nx, ny, nz, rad_zm%nn))
    end if

    if(LH_microphys_type /= LH_microphys_disabled ) then
      allocate(out_LH_zt(nx, ny, nz, LH_zt%nn))
      allocate(out_LH_sfc(nx, ny, nz, LH_sfc%nn))
    end if

    return

  end subroutine stats_init_clubb
!================================================================================== !
!                                                                                   !
!================================================================================== !
#ifndef CRM  
  subroutine hbuf_stats_init_clubb(namelist,deflist,unitlist,status,average_type,count,clubbcount)

   use stats_variables, only: &
       zt, LH_zt, zm, rad_zm, rad_zt, sfc, LH_sfc, l_output_rad_files
    use parameters_microphys, only: &
      LH_microphys_disabled, & ! Constant
      LH_microphys_type ! Variable

   implicit none

   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count, clubbcount, n, ii, jj, ncond

   character*8 name
   character*80 longname
   character*10 units

!  Local variables
   integer  :: i
   character*100  temp1, sub

   clubbcount = 0

!   Now call add fields
    do i = 1, zt%nn
    
      temp1 = trim(zt%f%var(i)%name)
      sub   = temp1
!      if (len(temp1) > 16) sub = temp1(1:16)
     
!       call addfld(trim(sub),trim(zt%f%var(i)%units),nnzp,&
!            'A',trim(zt%f%var(i)%description),phys_decomp)
       call add_to_namelist(count, clubbcount, trim(sub), trim(zt%f%var(i)%description),  &
            trim(zt%f%var(i)%units), 0)
    enddo
    
    do i = 1, zm%nn
    
      temp1 = trim(zm%f%var(i)%name)
      sub   = temp1
!      if (len(temp1) > 16) sub = temp1(1:16)
    
!      call addfld(trim(sub),trim(zm%f%var(i)%units),nnzp,&
!           'A',trim(zm%f%var(i)%description),phys_decomp)
       call add_to_namelist(count, clubbcount, trim(sub), trim(zm%f%var(i)%description),  & 
            trim(zm%f%var(i)%units), 0)
    enddo

    if (l_output_rad_files) then     
      do i = 1, rad_zt%nn
!        call addfld(trim(rad_zt%f%var(i)%name),trim(rad_zt%f%var(i)%units),nnzp,&
!           'A',trim(rad_zt%f%var(i)%description),phys_decomp)
       call add_to_namelist(count, clubbcount, trim(rad_zt%f%var(i)%name),  & 
            trim(rad_zt%f%var(i)%description), trim(rad_zt%f%var(i)%units), 0)
      enddo
    
      do i = 1, rad_zm%nn
!        call addfld(trim(rad_zm%f%var(i)%name),trim(rad_zm%f%var(i)%units),nnzp,&
!           'A',trim(rad_zm%f%var(i)%description),phys_decomp)
       call add_to_namelist(count, clubbcount, trim(rad_zm%f%var(i)%name),  & 
            trim(rad_zm%f%var(i)%description), trim(rad_zm%f%var(i)%units), 0)
      enddo
    endif 

    if ( LH_microphys_type /= LH_microphys_disabled ) then
       do i=1, LH_zt%nn
       call add_to_namelist(count, clubbcount, trim(LH_zt%f%var(i)%name),  & 
            trim(LH_zt%f%var(i)%description), trim(LH_zt%f%var(i)%units), 0)
       end do       
       do i=1, LH_sfc%nn
       call add_to_namelist(count, clubbcount, trim(LH_sfc%f%var(i)%name),  &
            trim(LH_sfc%f%var(i)%description), trim(LH_sfc%f%var(i)%units), 0)
       end do
    endif
    
    do i = 1, sfc%nn
       call add_to_namelist(count, clubbcount, trim(sfc%f%var(i)%name),  &
            trim(sfc%f%var(i)%description), trim(sfc%f%var(i)%units), 0)
    enddo

   return

  end subroutine hbuf_stats_init_clubb
  !================================================================================

  subroutine hbuf_clubb_output()

    use stats_variables, only: &
       zt, LH_zt, zm, rad_zm, rad_zt, sfc, LH_sfc, l_output_rad_files
    use parameters_microphys, only: &
      LH_microphys_disabled, & ! Constant
      LH_microphys_type ! Variable
    use hbuffer, only: hbuf_avg_put

    implicit none
    
    ! locale variables 
    integer :: i
    character*100  temp1, sub     

    do i = 1, zt%nn
      call hbuf_avg_put(trim(zt%f%var(i)%name), out_zt(1:nx, 1:ny, 2:nz, i), 1, nx, 1, ny, nzm, 1.)
    enddo

    do i = 1, zm%nn
      !Velocity level. Here we just simplely put the last nz-1 onto the pressure level. 
      call hbuf_avg_put(trim(zm%f%var(i)%name), out_zm(1:nx, 1:ny, 1:(nz-1), i),  &
           1, nx, 1, ny, nzm, 1.)
    enddo

    if (l_output_rad_files) then
      do i = 1, rad_zt%nn
        call hbuf_avg_put(trim(rad_zt%f%var(i)%name), &
             out_rad_zt(1:nx, 1:ny, 2:nz, i), 1, nx, 1, ny, nzm, 1.)
      enddo

      do i = 1, rad_zm%nn
        call hbuf_avg_put(trim(rad_zm%f%var(i)%name),  & 
             out_rad_zm(1:nx, 1:ny, 1:(nz-1), i), 1, nx, 1, ny, nzm, 1.)
      enddo
    endif

    if ( LH_microphys_type /= LH_microphys_disabled ) then
       do i=1, LH_zt%nn
         call hbuf_avg_put(trim(LH_zt%f%var(i)%name),  & 
              out_LH_zt(1:nx, 1:ny, 2:nz, i), 1, nx, 1, ny, nzm, 1.)
       end do

       do i=1, LH_sfc%nn
         ! For simplicity, hbuf_avg_put is also called for surface varialbes. 
         ! so zeroout values from level 2 to nz
         out_LH_sfc(:, :, 2:nz, i) = 0.0 
         call hbuf_avg_put(trim(LH_sfc%f%var(i)%name),  &
              out_LH_sfc(1:nx, 1:ny, 1:(nz-1), i), 1, nx, 1, ny, nzm, 1.)
       end do
    end if

    do i = 1, sfc%nn
      ! For simplicity, hbuf_avg_put is also called for surface varialbes. 
      ! so zeroout values from level 2 to nz
      out_sfc(:, :, 2:nz, i) = 0.0
      call hbuf_avg_put(trim(sfc%f%var(i)%name),  &
           out_sfc(1:nx, 1:ny, 1:(nz-1), i), 1, nx, 1, ny, nzm, 1.)
    enddo

    return 

  end subroutine hbuf_clubb_output
#endif /*CRM*/
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    !-----------------------------------------------------------------------
  subroutine stats_end_timestep_clubb(ix, jy)

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
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
        l_output_rad_files

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use parameters_microphys, only: &
      LH_microphys_disabled  ! Constant

    use parameters_microphys, only: &
      LH_microphys_type, & ! Variable(s)
      LH_microphys_calls


    implicit none


    integer, intent(in) :: ix
    integer, intent(in) :: jy

    ! Local Variables

    integer :: i, k
    logical :: l_error

    ! ---- Begin Code ----

    ! Check if it is time to write to file

    if ( .not. l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zt statistics at each vertical level.
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
          endif

        endif

      enddo
    enddo

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zm statistics at each vertical level.
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
          endif

        endif

      enddo
    enddo

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


    if (l_output_rad_files) then
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zt statistics at each vertical level.
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
            endif

          endif

        enddo
      enddo
    
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zm statistics at each vertical level.
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
            endif

          endif

        enddo
      enddo
    end if ! l_output_rad_files

    !  Look for errors by checking the number of sampling points
    !  for each variable in the sfc statistics at each vertical level.
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
          endif

        endif

      enddo
    enddo

    !  Stop the run if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      stop 'stats_end_timestep:  error(s) found'
    endif

    !  Compute averages
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

   !  Here we are not outputting the data, rather reading the stats into 
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.    
    do i = 1, zt%nn
      do k = 1, zt%kk 
         out_zt(ix,jy,k,i) = zt%x(1,1,k,i)
	 if(out_zt(ix,jy,k,i) /= out_zt(ix,jy,k,i)) out_zt(ix,jy,k,i) = 0.0
      enddo   
    enddo
    
    do i = 1, zm%nn
      do k = 1, zt%kk 
         out_zm(ix,jy,k,i) = zm%x(1,1,k,i)
	 if(out_zm(ix,jy,k,i) /= out_zm(ix,jy,k,i)) out_zm(ix,jy,k,i) = 0.0
      enddo   
    enddo

    if (l_output_rad_files) then 
      do i = 1, rad_zt%nn
        do k = 1, rad_zt%kk 
          out_rad_zt(ix,jy,k,i) = rad_zt%x(1,1,k,i)
	  if(out_rad_zt(ix,jy,k,i) /= out_rad_zt(ix,jy,k,i)) out_rad_zt(ix,jy,k,i) = 0.0
        enddo   
      enddo
    
      do i = 1, rad_zm%nn
        do k = 1, rad_zm%kk 
          out_rad_zm(ix,jy,k,i) = rad_zm%x(1,1,k,i)
	  if(out_rad_zm(ix,jy,k,i) /= out_rad_zm(ix,jy,k,i)) out_rad_zm(ix,jy,k,i) = 0.0
        enddo   
      enddo
    endif

    if ( LH_microphys_type /= LH_microphys_disabled ) then
      do i=1, LH_zt%nn
        do k=1, LH_zt%kk
          out_LH_zt(ix,jy,k,i) = LH_zt%x(1,1,k,i)
          if(out_LH_zt(ix,jy,k,i) /= out_LH_zt(ix,jy,k,i)) out_LH_zt(ix,jy,k,i) = 0.0 
        enddo
      enddo
      
      out_LH_sfc(ix,jy,:,:) = 0.0
      do i=1, LH_sfc%nn
        out_LH_sfc(ix,jy,1,i) = LH_sfc%x(1,1,1,i)
        if(out_LH_sfc(ix,jy,1,i) /= out_LH_sfc(ix,jy,1,i)) out_LH_sfc(ix,jy,1,i) = 0.0
      end do 
    endif
    
    out_sfc(ix, jy, :, :) = 0.0
    do i = 1, sfc%nn
      out_sfc(ix,jy,1,i) = sfc%x(1,1,1,i)   
      if(out_sfc(ix,jy,1,i) /= out_sfc(ix,jy,1,i)) out_sfc(ix,jy,1,i) = 0.0
    enddo

    !  Reset sample fields
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )
    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )
      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )
    end if
    if ( LH_microphys_type /= LH_microphys_disabled) then
      call stats_zero( LH_zt%kk, LH_zt%nn, LH_zt%x, LH_zt%n, LH_zt%l_in_update )
      call stats_zero( LH_sfc%kk, LH_sfc%nn, LH_sfc%x, LH_sfc%n, LH_sfc%l_in_update )
    end if
    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    return

  end subroutine stats_end_timestep_clubb
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, nn

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical, dimension(1,1,kk,nn), intent(out)                 :: l_in_update

    !  Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0_stat_rknd
      n(:,:,:,:) = 0_stat_nknd
      l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: nn, kk
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m=1,nn
      do k=1,kk

        if ( n(1,1,k,m) > 0 ) then
          x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m), kind=stat_rknd )
        end if

      end do
    end do

    return

  end subroutine stats_avg
#endif /* CLUBB_CRM*/
end module stat_clubb
