
module qbo
!--------------------------------------------------------------------
! This module performes a relaxation towards either a cyclic idealized
! QBO sequence derived from observations (here:28months) or towards
! observed equatorial wind data
!
! Author: Katja Matthes
! Date: February 2005
!
! implementation into WACCM standard version, K. Matthes, October 2005 
! 3 possibilities:
! default: no QBO relaxation
! options: 1. cyclic QBO time series, currently a 28month sequence available: qbocyclic28months.nc
!          2. observed QBO time series, currently CCMVal QBO series, extended: qboseries_ext.nc 
!          3. fixed QBOe or QBOw phase, currently : qboeast.nc, qbowest.nc
!
!--------------------------------------------------------------------
! modified by Anne Smith, December 2009
!
! Now accepts alternative input in the form of Fourier coefficients.
! In this case, the QBO wind is calculated from the expression
! u_qbo(k,n) = ubar_qbo(k) + sum_n[fcos_qbo(k,n) * cos(ffreq_qbo(n)*(cday-cday_ref)) 
!                                + fsin_qbo(k,n) * sin(ffreq_qbo(n)*(cday-cday_ref))]
! where cday is day of model run
!       cday_ref is a reference date so that historical data line up correctly
!       k is level index
!       n is coefficient index
!       sum_n is the sum over n
!---------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use pmgrid,         only: plon, plev, plevp
  use ppgrid,         only: pcols, pver
  use time_manager,   only: get_curr_date, get_curr_calday
  use phys_grid,      only: get_rlat_p
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use cam_history
  use error_messages, only: alloc_err
  use abortutils,     only: endrun  
  use cam_logfile,    only: iulog

  implicit none

  private
  save

!---------------------------------------------------------------------
! Public methods
!---------------------------------------------------------------------
  public               :: qbo_init               ! initialize qbo package
  public               :: qbo_timestep_init      ! interpolate to current time
  public               :: qbo_relax              ! relax zonal mean wind
  public               :: qbo_defaultopts        ! set default values of namelist variables
  public               :: qbo_setopts            ! get namelist input

!---------------------------------------------------------------------
! Private module data
!---------------------------------------------------------------------
  integer,parameter    :: qbo_dypm = 30          ! days in a cyclic qbo "month"

  integer              :: nm, np                 ! Indices for previous, next month
  integer              :: np1 = 0                ! tempory for next time index of dataset
  integer              :: ktop                   ! Model layers within qbo region
  integer              :: kbot                   ! Model layers within qbo region
  integer              :: timesiz                ! size of time dimension on dataset
  integer              :: levsiz                 ! size of lev  dimension on dataset
  integer              :: qbo_mons               ! length of cyclic qbo in months
  integer              :: delt                   ! time step in seconds
  real(r8)             :: qbo_days               ! length of cyclic qbo in days
  integer, allocatable :: date_qbo(:)            ! Date on qbo dataset (YYYYMMDD) (timesiz)
  integer, allocatable :: secd_qbo(:)            ! seconds of date (0-86399) (timesiz)

  real(r8)             :: cdaym                  ! dataset calendar day previous month
  real(r8)             :: cdayp                  ! dataset calendar day next month
  real(r8),allocatable :: u_qbo(:,:)             ! qbo winds(pver,timesiz)
  real(r8),allocatable :: tauz(:)
  real(r8)             :: u_tstep(pver)          ! qbo winds for this time step

  integer              :: coefsiz                ! size of coefficient dimension on fft dataset
  real(r8)             :: cday_ref               ! reference day for fft input
  real(r8),allocatable :: ubar_qbo(:)            ! qbo mean winds(pver)
  real(r8),allocatable :: fcos_qbo(:,:)          ! qbo cosine coefficients(pver,coefsiz)
  real(r8),allocatable :: fsin_qbo(:,:)          ! qbo sine coefficients(pver,coefsiz)
  real(r8),allocatable :: ffreq_qbo(:)           ! frequencies for expanding qbo coefficients (coefsiz)

!
! logistic for controling QBO relaxation 
!
  character(len=256)   :: qbo_forcing_file = 'NO_QBO_FILE'
  logical,public       :: qbo_use_forcing  = .FALSE.        ! .TRUE. => this package is active
  logical              :: qbo_cyclic       = .FALSE.        ! .TRUE. => assume cyclic qbo data

  logical              :: has_monthly_data = .TRUE.         ! .TRUE. => data file has monthly winds
                                                            ! .FALSE.=> data file has fft coefficients


contains

!================================================================================================

    subroutine qbo_defaultopts( &
         qbo_forcing_file_out,  &
         qbo_use_forcing_out,   &
         qbo_cyclic_out        )

      character(len=*), intent(out), optional :: qbo_forcing_file_out
      logical,          intent(out), optional :: qbo_use_forcing_out
      logical,          intent(out), optional :: qbo_cyclic_out


      if ( present(qbo_forcing_file_out) ) then
         qbo_forcing_file_out = qbo_forcing_file
      endif
      if ( present(qbo_use_forcing_out) ) then
         qbo_use_forcing_out = qbo_use_forcing
      endif
      if ( present(qbo_cyclic_out) ) then
         qbo_cyclic_out = qbo_cyclic
      endif

    end subroutine qbo_defaultopts

!================================================================================================

    subroutine qbo_setopts(    &
         qbo_forcing_file_in,  &
         qbo_use_forcing_in,   &
         qbo_cyclic_in        )

      character(len=*), intent(in), optional :: qbo_forcing_file_in
      logical,          intent(in), optional :: qbo_use_forcing_in
      logical,          intent(in), optional :: qbo_cyclic_in

      if ( present(qbo_forcing_file_in) ) then
         qbo_forcing_file = qbo_forcing_file_in
      endif
      if ( present(qbo_use_forcing_in) ) then
         qbo_use_forcing = qbo_use_forcing_in
      endif
      if ( present(qbo_cyclic_in) ) then
         qbo_cyclic = qbo_cyclic_in
      endif

    end subroutine qbo_setopts

!================================================================================================


  subroutine qbo_init
!---------------------------------------------------------------------
! initialize qbo module
!---------------------------------------------------------------------
   use time_manager, only : get_step_size
   use ref_pres,     only : pref_mid
   use wrap_nf
#if (defined SPMD )
   use mpishorthand
#endif

!---------------------------------------------------------------------
! Local workspace
!---------------------------------------------------------------------
    real(r8), parameter :: hPa2Pa = 100._r8

    integer  :: i, m, y, n         ! testing indexes
    integer  :: yr, mon, day       ! components of a date
    integer  :: ncdate             ! current date in integer format [yyyymmdd]
    integer  :: ncsec              ! current time of day [seconds]
    integer  :: ncid               ! netcdf ID for input dataset 
    integer  :: astat              ! allocate status

    integer  :: k,kk               ! vertical interpolation indexes
    integer  :: levdimid           ! netcdf id for level dimension
    integer  :: timdimid           ! netcdf id for time  dimension
    integer  :: levid              ! netcdf id for level variable
    integer  :: dateid             ! netcdf id for date variable
    integer  :: secdid             ! netcdf id for seconds variable
    integer  :: uqboid             ! netcdf id for qbo wind variable

    real(r8) :: fac1, fac2         ! interpolation factors
    real(r8) :: calday             ! current calendar day
    real(r8) :: cday               ! current day (in cycle, if cycling)
    real(r8), allocatable :: p_inp(:)    ! qbo pressure levels(levsiz)
    real(r8), allocatable :: u_inp(:,:)  ! input qbo winds(levsiz,timesiz) have to be sorted
                                         ! from top to bottom like winds in model

    integer  :: fftdimid           ! netcdf id for fft coefficient dimension
    integer  :: ubarqboid          ! netcdf id for mean wind variable
    integer  :: cosqboid           ! netcdf id for cosine coefficient variable
    integer  :: sinqboid           ! netcdf id for sine coefficient variable
    integer  :: freqqboid          ! netcdf id for fft frequency variable
    integer  :: refdatid           ! netcdf id for reference date

    integer               :: ref_date       ! date corresponding to beginning of QBO data
    integer               :: yr_ref         ! year for reference day
    real(r8), allocatable :: ubar_inp(:)    ! input time-mean winds for fft input (levsiz)
    real(r8), allocatable :: cosq_inp(:,:)  ! input cosine coefficients for fft input (levsiz,coefsiz)
    real(r8), allocatable :: sinq_inp(:,:)  ! input sine coefficients for fft input (levsiz,coefsiz)

    logical  :: found

    if( .not. qbo_use_forcing ) then
       return
    end if
!---------------------------------------------------------------------
! SPMD: Master does all the work of reading files.  Sends needed info to slaves
!---------------------------------------------------------------------
    if( masterproc ) then
!---------------------------------------------------------------------
! Get qbo file
!---------------------------------------------------------------------
       call wrap_open( qbo_forcing_file, NF90_NOERR, ncid )
       write(iulog,*) 'qbo_init: successfully opened ',trim(qbo_forcing_file)
!---------------------------------------------------------------------
! Figure out if the file contains qbo winds by month or fft coefficients 
!   by looking for the variable DATE
!---------------------------------------------------------------------
       call wrap_inq_varid( ncid, 'date' , dateid, abort=.false.  )

       if (dateid > 0) then
          has_monthly_data=.true.
       else
          has_monthly_data=.false.
       end if
       write(iulog,*) 'has_monthly_data=', has_monthly_data

!---------------------------------------------------------------------
! Get and check dimension info
!---------------------------------------------------------------------
       if (has_monthly_data) then
          call wrap_inq_dimid( ncid, 'level' , levdimid  )
          call wrap_inq_dimid( ncid, 'time', timdimid  )

          call wrap_inq_dimlen( ncid, levdimid, levsiz   )
          call wrap_inq_dimlen( ncid, timdimid, timesiz  )

          call wrap_inq_varid( ncid, 'date' , dateid  )
          call wrap_inq_varid( ncid, 'secs' , secdid  )
          call wrap_inq_varid( ncid, 'qbo'  , uqboid  )
          call wrap_inq_varid( ncid, 'level', levid   )
       else
          call wrap_inq_dimid( ncid, 'level' , levdimid  )
          call wrap_inq_dimid( ncid, 'ncoef', fftdimid  )

          call wrap_inq_dimlen( ncid, levdimid, levsiz   )
          call wrap_inq_dimlen( ncid, fftdimid, coefsiz  )

          call wrap_inq_varid( ncid, 'ref_date', refdatid  )

          call wrap_inq_varid( ncid, 'ubar'   , ubarqboid  )
          call wrap_inq_varid( ncid, 'cosqbo' , cosqboid  )
          call wrap_inq_varid( ncid, 'sinqbo' , sinqboid  )
          call wrap_inq_varid( ncid, 'freqqbo', freqqboid   )
          call wrap_inq_varid( ncid, 'level'  , levid   )
       end if
    end if


!---------------------------------------------------------------------
! Broadcast the logical flag has_monthly_data
!---------------------------------------------------------------------
#if (defined SPMD )
    call mpibcast( has_monthly_data, 1, mpilog, 0, mpicom )
#endif

    if (.not.has_monthly_data) then
!---------------------------------------------------------------------
! Broadcast coefsiz
!---------------------------------------------------------------------
#if (defined SPMD )
       call mpibcast( coefsiz, 1, mpiint, 0, mpicom )
#endif
    end if

    if (has_monthly_data) then
#if (defined SPMD )
!---------------------------------------------------------------------
! Broadcast the time size to all tasks for case with monthly input
!---------------------------------------------------------------------
          call mpibcast( timesiz, 1, mpiint, 0, mpicom )
#endif

       delt = get_step_size()
       if( qbo_cyclic ) then
          qbo_mons = timesiz 
          qbo_days = qbo_mons*qbo_dypm
       end if

       allocate( date_qbo(timesiz), secd_qbo(timesiz), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'qbo_init: failed to allocate date_qbo ... secd_qbo; error = ',astat
          call endrun
       end if
    end if
!---------------------------------------------------------------------
! Get input variables
!---------------------------------------------------------------------
master_proc : &
    if( masterproc ) then
!---------------------------------------------------------------------
! Allocate arrays, depending on type of input data: monthly or fft
!---------------------------------------------------------------------
       if (has_monthly_data) then
          allocate( u_inp(levsiz,timesiz), p_inp(levsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate u_inp ... p_inp; error = ',astat
             call endrun
          end if
          astat = nf90_get_var( ncid, uqboid, u_inp )
          if (astat/=NF90_NOERR) then
             write(iulog,*) "QBO: NF90_GET_VAR: error reading varid =", uqboid
             call endrun
          end if
          call wrap_get_var_realx( ncid, levid , p_inp )
       else
          allocate( p_inp(levsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate p_inp; error = ',astat
             call endrun
          end if
          allocate( ubar_inp(levsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate ubar_inp; error = ',astat
             call endrun
          end if
          allocate( cosq_inp(levsiz,coefsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate cosq_inp; error = ',astat
             call endrun
          end if
          allocate( sinq_inp(levsiz,coefsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate sinq_inp; error = ',astat
             call endrun
          end if
          allocate( ffreq_qbo(coefsiz), stat=astat )
          if( astat /= 0 ) then
             write(iulog,*) 'qbo_init: failed to allocate ffreq_qbo; error = ',astat
             call endrun
          end if
          call wrap_get_var_realx( ncid, levid ,    p_inp )
          call wrap_get_var_realx( ncid, ubarqboid, ubar_inp )
          astat = nf90_get_var( ncid, cosqboid, cosq_inp )
          if (astat/=NF90_NOERR) then
             write(iulog,*) "QBO: NF90_GET_VAR: error reading varid =", cosqboid
             call endrun
          end if
          astat = nf90_get_var( ncid, sinqboid, sinq_inp )
          if (astat/=NF90_NOERR) then
             write(iulog,*) "QBO: NF90_GET_VAR: error reading varid =", sinqboid
             call endrun
          end if
          call wrap_get_var_realx( ncid, freqqboid, ffreq_qbo )
          call wrap_get_scalar_int(ncid, refdatid,  ref_date )

          call bnddyi( ref_date, 0, cday_ref )
          yr_ref = ref_date/10000
          cday_ref = cday_ref + yr_ref*365._r8

       end if

!---------------------------------------------------------------------
! Convert from millibars to pascals
!---------------------------------------------------------------------
       p_inp(:) = p_inp(:)*hPa2Pa
       write(iulog,*) 'qbo_init: p_inp', p_inp/hPa2Pa

       if (has_monthly_data) then
          call wrap_get_var_int( ncid, dateid, date_qbo )
          call wrap_get_var_int( ncid, secdid, secd_qbo )

          write(iulog,*) 'qbo_init: u_inp', u_inp(:,1)
       end if

!---------------------------------------------------------------------
! Find first model layer within qbo range
!---------------------------------------------------------------------
       do ktop = 1,pver
          if( pref_mid(ktop) >= p_inp(1) ) then
             exit
          end if
       end do
       write(iulog,*) 'qbo_init: ktop = ', ktop, pref_mid(ktop)/hPa2Pa

!---------------------------------------------------------------------
! Find last model layer within qbo range
!---------------------------------------------------------------------
       do kbot = pver,ktop,-1
          if( pref_mid(kbot) <= p_inp(levsiz) ) then
             exit
          end if
       end do
       write(iulog,*) 'qbo_init: kbot = ', kbot, pref_mid(kbot)/hPa2Pa

       if (has_monthly_data) then
           allocate( u_qbo(ktop:kbot,timesiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate u_qbo; error = ',astat
              call endrun
           end if
        else
           allocate( ubar_qbo(ktop:kbot), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate ubar_qbo; error = ',astat
              call endrun
           end if
           allocate( fcos_qbo(ktop:kbot,coefsiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate fcos_qbo; error = ',astat
              call endrun
           end if
           allocate( fsin_qbo(ktop:kbot,coefsiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate fsin_qbo; error = ',astat
              call endrun
           end if
        end if
!---------------------------------------------------------------------
! Vertically interpolate input winds to model reference pressures
!---------------------------------------------------------------------
        do k = ktop,kbot
           do kk = 1,levsiz-1
              if( p_inp(kk+1) > pref_mid(k) ) then
                 exit
              end if
           end do
           fac1 = (pref_mid(k) - p_inp(kk+1)) / (p_inp(kk) - p_inp(kk+1))
           fac2 = (p_inp(kk) - pref_mid(k)  ) / (p_inp(kk) - p_inp(kk+1))
           if (has_monthly_data) then
              u_qbo(k,:) = u_inp(kk,:)*fac1 + u_inp(kk+1,:)*fac2
           else
              ubar_qbo(k) = ubar_inp(kk)*fac1 + ubar_inp(kk+1)*fac2
              fcos_qbo(k,:) = cosq_inp(kk,:)*fac1 + cosq_inp(kk+1,:)*fac2
              fsin_qbo(k,:) = sinq_inp(kk,:)*fac1 + sinq_inp(kk+1,:)*fac2
           endif
        end do

        if (has_monthly_data) then
           deallocate( u_inp, p_inp )
           write(iulog,*) 'qbo_init: u', u_qbo(ktop:kbot,1)
        else
           deallocate( p_inp )
           deallocate( ubar_inp )
           deallocate( cosq_inp )
           deallocate( sinq_inp )
           write(iulog,*) 'qbo_init: fcos_qbo', fcos_qbo(ktop:kbot,1)
        end if


!---------------------------------------------------------------------
! Dates are not used for cyclic QBO
!---------------------------------------------------------------------
        if( has_monthly_data .and. qbo_cyclic ) then
           secd_qbo(:) = 0
           date_qbo(:) = 0
       end if
    end if master_proc

#if (defined SPMD )
!---------------------------------------------------------------------
! Broadcast the vertical limits
!---------------------------------------------------------------------
    call mpibcast( ktop, 1, mpiint, 0, mpicom )
    call mpibcast( kbot, 1, mpiint, 0, mpicom )
#endif

    if( .not. masterproc ) then
       if (has_monthly_data) then
           allocate( u_qbo(ktop:kbot,timesiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate u_qbo; error = ',astat
              call endrun
           end if
        else
           allocate( ubar_qbo(ktop:kbot), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate ubar_qbo; error = ',astat
              call endrun
           end if
           allocate( fcos_qbo(ktop:kbot,coefsiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate fcos_qbo; error = ',astat
              call endrun
           end if
           allocate( fsin_qbo(ktop:kbot,coefsiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate fsin_qbo; error = ',astat
              call endrun
           end if
           allocate( ffreq_qbo(coefsiz), stat=astat )
           if( astat /= 0 ) then
              write(iulog,*) 'qbo_init: failed to allocate ffreq_qbo; error = ',astat
              call endrun
           end if
        end if
    end if

!---------------------------------------------------------------------
! Broadcast input data  to all tasks
!---------------------------------------------------------------------
    kk = kbot-ktop+1
    if (has_monthly_data) then
#if (defined SPMD )
          call mpibcast( u_qbo   , (kbot-ktop+1)*timesiz, mpir8, 0, mpicom )
          call mpibcast( date_qbo, timesiz, mpiint, 0, mpicom )
          call mpibcast( secd_qbo, timesiz, mpiint, 0, mpicom )
#endif
    else
#if (defined SPMD )
          call mpibcast( cday_ref  , 1,                     mpir8, 0, mpicom )
#endif
#if (defined SPMD )
          call mpibcast( ubar_qbo  , (kbot-ktop+1),         mpir8, 0, mpicom )
          call mpibcast( fcos_qbo  , (kbot-ktop+1)*coefsiz, mpir8, 0, mpicom )
          call mpibcast( fsin_qbo  , (kbot-ktop+1)*coefsiz, mpir8, 0, mpicom )
#endif
#if (defined SPMD )
          call mpibcast( ffreq_qbo , coefsiz,               mpir8, 0, mpicom )
#endif
    endif

!---------------------------------------------------------------------
! setup vertical factor
!---------------------------------------------------------------------
    allocate( tauz(ktop-1:kbot+1), stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'qbo_init: failed to allocate tauz; error = ',astat
       call endrun
    end if

    tauz(ktop-1)    = 2._r8
    tauz(kbot+1)    = 2._r8
    tauz(ktop:kbot) = 1._r8

!---------------------------------------------------------------------
! Get current day of year and date
!---------------------------------------------------------------------
    calday = get_curr_calday()
    call get_curr_date( yr, mon, day, ncsec )
    write(iulog,*) 'qbo_init: get_curr_date = ', yr,mon,day,ncsec
    ncdate = yr*10000 + mon*100 + day

!---------------------------------------------------------------------
! Set current day in run or in qbo cycle
!---------------------------------------------------------------------
    cday = calday + yr*365._r8
    write(iulog,*) 'qbo_init: cday = ', cday
    if( has_monthly_data .and. qbo_cyclic ) then
       cday = mod( cday,qbo_days )
    end if
    write(iulog,*) 'qbo_init: cday = ', cday

    if( has_monthly_data ) then
       if( qbo_cyclic ) then
!---------------------------------------------------------------------
! Set up cyclic qbo
!---------------------------------------------------------------------
! Set past and future month indexes
!---------------------------------------------------------------------
          if( cday == 0._r8 ) then
             nm  = qbo_mons - 1
             np  = qbo_mons
          else
             nm  = cday / qbo_dypm
             np  = mod( nm,qbo_mons ) + 1
             if( nm == 0 ) then
                nm = qbo_mons
             end if
          end if
          write(iulog,*) 'qbo_init: nm,np = ', nm,np
!---------------------------------------------------------------------
! Set past and future days for data, generate day for cyclic data
!---------------------------------------------------------------------
          cdayp = np * qbo_dypm
          cdaym = nm * qbo_dypm
       else
!---------------------------------------------------------------------
! Set up noncyclic qbo
!---------------------------------------------------------------------
          found = .false.
          do nm = 1, timesiz-1
             np = nm+1
             call bnddyi( date_qbo(nm), secd_qbo(nm), cdaym )
             call bnddyi( date_qbo(np), secd_qbo(np), cdayp )
             yr    = date_qbo(nm)/10000
             cdaym = cdaym + yr*365._r8
             yr    = date_qbo(np)/10000
             cdayp = cdayp + yr*365._r8
             write(iulog,*) 'qbo_init: nm, date_qbo(nm), date_qbo(np), cdaym, cdayp = ', &
                  nm, date_qbo(nm), date_qbo(np), cdaym, cdayp
             if( cday >= cdaym .and. cday < cdayp ) then
                found = .true.
                exit
             end if
          end do
          if( .not. found ) then
             call endrun( 'QBO_INIT: failed to find bracketing dates' )
          end if
       end if
       write(iulog,*) 'qbo_init: cdaym,cdayp = ', cdaym,cdayp
    endif

!---------------------------------------------------------------------
! Initialize output buffer for two fields: QBO forcing wind and wind tendency of qbo relaxation
!----------------------------------------------------------------------
! 
!----------------------------------------------------------------------
    call addfld ('QBOTEND ','M/S/S   ',plev, 'A','Wind tendency from QBO relaxation',phys_decomp)
    call add_default ('QBOTEND', 1, ' ')

    call addfld ('QBO_U0 ','M/S   ',plev, 'A','Specified wind used for QBO',phys_decomp)
    call add_default ('QBO_U0', 1, ' ')


    write(iulog,*) 'end of qbo_init'

  end subroutine qbo_init

!================================================================================================

  subroutine qbo_timestep_init
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate QBO zonal wind to current time
! 
! Method: Linear interpolation between dates on QBO file, 
!         vertically and horizontally
! 
! Author: 
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
    integer  :: k                   ! level index
    integer  :: yr, mon, day        ! components of a date
    integer  :: yrl, monl, dayl     ! components of a date
    integer  :: ncdate              ! current date in integer format [yyyymmdd]
    integer  :: ncsec               ! current time of day [seconds]
    integer  :: ncsecl              ! current time of day [seconds]

    real(r8) :: fact1, fact2        ! time interpolation factors
    real(r8) :: calday              ! day of year at end of present time step
    real(r8) :: caldayl             ! day of year at begining of present time step
    real(r8) :: cday                ! day within qbo period at end of present time step
    real(r8) :: cdayl               ! day within qbo period at begining of present time step
    real(r8) :: deltat              ! time difference (days) between cdaym and cdayp

    integer  :: n                   ! coefficient index
    real(r8) :: ccc                 ! cosine term for expanding coefficients
    real(r8) :: sss                 ! sine term for expanding coefficients

    logical  :: new_interval        ! flag for new time interval

has_qbo_forcing : &
    if( qbo_use_forcing ) then
!-----------------------------------------------------------------------
! Get current day of year and date
!-----------------------------------------------------------------------
       caldayl = get_curr_calday( -delt )
       call get_curr_date( yrl, monl, dayl, ncsecl, -delt )
       calday  = get_curr_calday()
       call get_curr_date( yr, mon, day, ncsec )
       ncdate = yr*10000 + mon*100 + day
#ifdef QBO_DIAGS
       write(iulog,*) 'qbo_timestep_init: calday = ', calday
       write(iulog,*) 'qbo_timestep_init: ncdate = ', ncdate
#endif

!-----------------------------------------------------------------------
! Set current day in run or in qbo cycle
!-----------------------------------------------------------------------
       cday  = calday + yr*365._r8
       cdayl = caldayl + yrl*365._r8
#ifdef QBO_DIAGS
       write(iulog,*) 'qbo_timestep_init: cday = ', cday
#endif

       if( has_monthly_data ) then
!-----------------------------------------------------------------------
! Time interpolation for cases with monthly input data
!-----------------------------------------------------------------------
          if( .not. qbo_cyclic ) then
             new_interval = cday > cdayp
          else
             cday  = mod( cday,qbo_days )
             cdayl = mod( cdayl,qbo_days )
             if( cday > cdayl ) then
                new_interval = cday > cdayp
             else
                new_interval = .true.
             end if
          end if
#ifdef QBO_DIAGS
          write(iulog,*) 'qbo_timestep_init: cday = ', cday
#endif
!-----------------------------------------------------------------------
! If model time is past current forward timeslice, then switch to next one.
! If qbo_cyclic = .true. interpolation between end and beginning of data (np == 1).  
! Note that np is never 1 when qbo_cyclic is .false.
!-----------------------------------------------------------------------
next_interval : &
          if( new_interval ) then

!-----------------------------------------------------------------------
! Increment index of future time sample
!-----------------------------------------------------------------------
             if( qbo_cyclic ) then
                np1 = mod( np,qbo_mons ) + 1
             else
                np1 = np + 1
             end if
             if( np1 > timesiz ) then
                call endrun ('QBO_TIMESTEP_INIT: Attempt to go past end of QBO data')
             end if
!-----------------------------------------------------------------------
! Set indexes into u table
!-----------------------------------------------------------------------
             nm = np
             np = np1
#ifdef QBO_DIAGS
             write(iulog,*) 'qbo_timestep_init: nm,np = ', nm,np
             write(iulog,*) 'qbo_timestep_init: date_qbo(np), secd_qbo(np) = ', date_qbo(np), secd_qbo(np)
#endif

!-----------------------------------------------------------------------
! Set past and future days for data, generate day for cyclic data
!-----------------------------------------------------------------------
             cdaym = cdayp
             if( qbo_cyclic ) then
                cdayp = np * qbo_dypm
             else
                call bnddyi( date_qbo(np), secd_qbo(np), cdayp )
                yr = date_qbo(np)/10000
                cdayp = cdayp + yr*365._r8
             end if
#ifdef QBO_DIAGS
             write(iulog,*) 'qbo_timestep_init: cdaym,cdayp = ', cdaym,cdayp
#endif
             if( np /= 1 .and. cday > cdayp ) then
                write(iulog,*) 'qbo_timestep_init: Input qbo for date',date_qbo(np),' sec ',secd_qbo(np), &
                  'does not exceed model date',ncdate,' sec ',ncsec,' Stopping.'
                call endrun
             end if
          end if next_interval

!-----------------------------------------------------------------------
! Determine time interpolation factors.  Account for December-January 
! interpolation if dataset is being cycled.
!-----------------------------------------------------------------------
          if( qbo_cyclic .and. np == 1 ) then                   ! Dec-Jan interpolation
             deltat = cdayp + qbo_days - cdaym
             if (cday > cdayp) then                         ! We are in December
                fact1 = (cdayp + qbo_days - cday)/deltat
                fact2 = (cday - cdaym)/deltat
             else                                           ! We are in January
                fact1 = (cdayp - cday)/deltat
                fact2 = (cday + qbo_days - cdaym)/deltat
             end if
          else
             deltat = cdayp - cdaym
             fact1 = (cdayp - cday )/deltat
             fact2 = (cday  - cdaym)/deltat
          end if
#ifdef QBO_DIAGS
          write(iulog,*) 'qbo_timestep_init: fact1,fact2 = ', fact1, fact2
#endif

!-----------------------------------------------------------------------
! Time interpolation
!-----------------------------------------------------------------------
          do k = ktop, kbot
             u_tstep(k) = u_qbo(k,nm)*fact1 + u_qbo(k,np)*fact2
          end do
          if( ktop > 1 ) then
             u_tstep(ktop-1) = u_tstep(ktop)
          end if
          if( kbot < pver ) then
             u_tstep(kbot+1) = u_tstep(kbot)
          end if

       else

!-----------------------------------------------------------------------
! Wind at this timestep for fft input data
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Set past and future days for data, generate winds for current day from fft data
!-----------------------------------------------------------------------

          do k = ktop, kbot
             u_tstep(k) =  ubar_qbo(k)
          end do
          do n=1,coefsiz
             ccc = cos(ffreq_qbo(n)*(cday-cday_ref))
             sss = sin(ffreq_qbo(n)*(cday-cday_ref))
             do k = ktop, kbot
                u_tstep(k) =  u_tstep(k) + fcos_qbo(k,n)*ccc + fsin_qbo(k,n)*sss
             end do
          end do
          if( ktop > 1 ) then
             u_tstep(ktop-1) = u_tstep(ktop)
          end if
          if( kbot < pver ) then
             u_tstep(kbot+1) = u_tstep(kbot)
          end if
       end if

#ifdef QBO_DIAGS
       write(iulog,*) 'qbo_timestep_init: u_tstep ', u_tstep(ktop:kbot)
#endif

    end if has_qbo_forcing

  end subroutine qbo_timestep_init

!================================================================================================

    subroutine qbo_relax( state, ptend )
      use shr_assert_mod, only: shr_assert
!------------------------------------------------------------------------
! relax zonal mean wind towards qbo sequence 
!------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------------------
    type(physics_state), intent(in)    :: state                ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend                ! individual parameterization tendencies

!--------------------------------------------------------------------------------
! Local variables
!--------------------------------------------------------------------------------
    integer                            :: lchnk                ! chunk identifier
    integer                            :: ncol                 ! number of atmospheric columns
    integer                            :: i, k                 ! loop indicies
    integer                            :: kl, ku               ! loop indicies

    real(r8)                           :: tauxi(pcols)         ! latitudes in radians for present chunck
    real(r8)                           :: tauzz
    real(r8)                           :: u
    real(r8)                           :: rlat                 ! latitudes in radians for present chunck
    real(r8)                           :: crelax               ! relaxation constant

    real(r8), parameter                :: tconst = 10._r8      ! relaxation time constant in days     
    real(r8), parameter                :: tconst1 = tconst * 86400._r8
    real(r8)                           :: qbo_u0(pcols,pver)   ! QBO wind used for driving parameterization

   lchnk = state%lchnk
   ncol  = state%ncol

   call physics_ptend_init(ptend, state%psetcols, 'qbo', lu=.true.)

has_qbo_forcing : &
   if( qbo_use_forcing ) then

   ! Make sure uzm is initialized if qbo is on.
   call shr_assert(allocated(state%uzm), &
        msg="Dycore has not allocated state%uzm before qbo_relax is called!")

    kl          = max( 1,ktop-1 )
    ku          = min( plev,kbot+1 )
!--------------------------------------------------------------------------------
! get latitude in radians for present chunk
!--------------------------------------------------------------------------------
    do i = 1,ncol
       rlat     = get_rlat_p( lchnk, i )
       tauxi(i) = tconst1*taux( rlat )
    end do

    qbo_u0(:,:) = 0._r8

    do k = kl,ku
      tauzz = tauz(k)
      u     = u_tstep(k)
      do i = 1,ncol
!--------------------------------------------------------------------------------
! determine relaxation constant 
!--------------------------------------------------------------------------------
         crelax = tauxi(i)*tauzz
         if( crelax /= 0._r8 ) then
            crelax = 1._r8 / crelax
!--------------------------------------------------------------------------------
! do relaxation of zonal mean wind
!--------------------------------------------------------------------------------

         if(u < 50.0_r8) then
            ptend%u(i,k) = crelax * (u - state%uzm(i,k))
         end if
         end if
!--------------------------------------------------------------------------------
! variable representing QBO wind
!--------------------------------------------------------------------------------
         if((u < 50.0_r8) .and. (crelax /= 0._r8)) then
            qbo_u0(i,k) = u/tauzz/tauxi(i)*tconst1
         end if
      end do
    end do

!--------------------------------------------------------------------------------
!output tendency of relaxation to monthly ('h1') output file
!--------------------------------------------------------------------------------
    call outfld( 'QBOTEND', ptend%u(:,:), pcols, lchnk )

!--------------------------------------------------------------------------------
!output specified QBO wind to h0 output file
!--------------------------------------------------------------------------------
    call outfld( 'QBO_U0', qbo_u0, pcols, lchnk )
   end if has_qbo_forcing

   end subroutine qbo_relax

!================================================================================================

  function taux( rlat )
!------------------------------------------------------------------------
! calculates relaxation constant in latitude
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!       ... dummy arguments
!------------------------------------------------------------------------
    real(r8), intent(in)  :: rlat    ! latitude in radians for present chunk

!------------------------------------------------------------------------
!       ... local variables
!------------------------------------------------------------------------
    real(r8), parameter   :: factor = 1._r8/(2._r8*0.174532925_r8*0.174532925_r8)
    real(r8)              :: alat    ! abs rlat

!------------------------------------------------------------------------
!       ... function declaration
!------------------------------------------------------------------------
    real(r8)              :: taux    ! relaxation constant in latitude
    

    alat = abs( rlat )
    if( alat <= .035_r8 ) then
!------------------------------------------------------------------------
! rlat=0.035 (latitude in radians): rlat*180/pi=2 degrees, around equator full relaxation
!------------------------------------------------------------------------
       taux = 1._r8
    else if( alat <= .384_r8)  then
!------------------------------------------------------------------------
! from 6 to 22 degrees latitude weakening of relaxation with Gaussian distribution
! half width=10° =>  in radians: 0.174532925
!------------------------------------------------------------------------
       taux = exp( rlat*rlat*factor )
    else
!------------------------------------------------------------------------
! other latitudes no relaxation
!------------------------------------------------------------------------
       taux = 0._r8
    end if
 
  end function taux

end module qbo
