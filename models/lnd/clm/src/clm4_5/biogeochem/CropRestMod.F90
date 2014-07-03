module CropRestMod
  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Read/Write to/from Crop info to CLM restart file. 
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CropRest        ! Restart prognostic crop model
  public :: CropRestYear    ! Get the number of years crop has spunup
  public :: CropRestIncYear ! Increment the crop spinup years
  !
  ! !PRIVATE DATA MEMBERS:
  integer :: restyear = 0 ! Restart year from the initial conditions file, 
                          ! incremented as time elapses
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CropRest ( bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write Crop restart data
    !
    ! !USES:
    use clmtype
    use restUtilMod  
    use pio,       only: file_desc_t
    use ncdio_pio, only: ncd_int, ncd_double, ncd_log
    use decompMod, only: bounds_type, get_proc_global
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j                      ! indices 
    real(r8):: m                          ! multiplier for the exit_spinup code
    logical :: readvar                    ! determine if variable is on initial file
    integer :: ier                        ! error status
    integer, pointer :: temp1d(:)         ! temporary
    !-----------------------------------------------------------------------

    ! Prognostic crop restart year 
    call restartvar(ncid=ncid, flag=flag,  varname='restyear', xtype=ncd_int,  &
         long_name='Number of years prognostic crop ran', units="years", &
         interpinic_flag='copy', readvar=readvar, data=restyear)
    if (flag=='read' .and. readvar)  call checkDates( )

    ! peaklai interpinic 
    call restartvar(ncid=ncid, flag=flag,  varname='peaklai', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='Flag if at max allowed LAI or not', &
         flag_values=(/0,1/), nvalid_range=(/0,1/), &
         flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) , &
         interpinic_flag='interp', readvar=readvar, data=pps%peaklai)

    ! idop 
    call restartvar(ncid=ncid, flag=flag,  varname='idop', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='Date of planting', units='jday', nvalid_range=(/1,366/), & 
         interpinic_flag='interp', readvar=readvar, data=pps%idop)

    ! aleaf
    call restartvar(ncid=ncid, flag=flag,  varname='aleaf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='leaf allocation coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%aleaf)

    ! aleafi
    call restartvar(ncid=ncid, flag=flag,  varname='aleafi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Saved leaf allocation coefficient from phase 2', &
         units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%aleafi)

    ! astem
    call restartvar(ncid=ncid, flag=flag,  varname='astem', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='stem allocation coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%astem)

    ! astemi
    call restartvar(ncid=ncid, flag=flag,  varname='astemi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Saved stem allocation coefficient from phase 2',&
         units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%astemi)

    ! htmx 
    call restartvar(ncid=ncid, flag=flag,  varname='htmx', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='max height attained by a crop during year',&
         units='m', &
         interpinic_flag='interp', readvar=readvar, data=pps%htmx)

    ! hdidx
    call restartvar(ncid=ncid, flag=flag,  varname='hdidx', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='cold hardening index', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%hdidx)

    ! vf
    call restartvar(ncid=ncid, flag=flag,  varname='vf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vernalization factor', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%vf)

    ! cumvd 
    call restartvar(ncid=ncid, flag=flag,  varname='cumvd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='cumulative vernalization d', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%cumvd)

    ! croplive 
    allocate(temp1d(bounds%begp:bounds%endp))
    if (flag == 'write') then 
       do p= bounds%begp,bounds%endp
          if (pps%croplive(p)) then
             temp1d(p) = 1
          else
             temp1d(p) = 0
          end if
       end do
    end if
    call restartvar(ncid=ncid, flag=flag,  varname='croplive', xtype=ncd_log,  &
         dim1name='pft', &
         long_name='Flag that crop is alive, but not harvested', &
         interpinic_flag='interp', readvar=readvar, data=temp1d)
    if (flag == 'read') then 
       do p= bounds%begp,bounds%endp
          if (temp1d(p) == 1) then
             pps%croplive(p) = .true.
          else
             pps%croplive(p) = .false.
          end if
       end do
    end if
    deallocate(temp1d)

    ! cropplant 
    allocate(temp1d(bounds%begp:bounds%endp))
    if (flag == 'write') then 
       do p= bounds%begp,bounds%endp
          if (pps%cropplant(p)) then
             temp1d(p) = 1
          else
             temp1d(p) = 0
          end if
       end do
    end if
    call restartvar(ncid=ncid, flag=flag,  varname='cropplant', xtype=ncd_log,  &
         dim1name='pft', &
         long_name='Flag that crop is planted, but not harvested' , &
         interpinic_flag='interp', readvar=readvar, data=temp1d)
    if (flag == 'read') then 
       do p= bounds%begp,bounds%endp
          if (temp1d(p) == 1) then
             pps%cropplant(p) = .true.
          else
             pps%cropplant(p) = .false.
          end if
       end do
    end if
    deallocate(temp1d)

    ! harvdate 
    call restartvar(ncid=ncid, flag=flag,  varname='harvdate', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='harvest date', units='jday', nvalid_range=(/1,366/), & 
         interpinic_flag='interp', readvar=readvar, data=pps%harvdate)

    ! gdd1020
    call restartvar(ncid=ncid, flag=flag,  varname='gdd1020', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='20 year average of growing degree-days base 10C from planting', units='ddays', &
         interpinic_flag='interp', readvar=readvar, data=pps%gdd1020)

    ! gdd820
    call restartvar(ncid=ncid, flag=flag,  varname='gdd820', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='20 year average of growing degree-days base 8C from planting', units='ddays', &
         interpinic_flag='interp', readvar=readvar, data=pps%gdd820)

    ! gdd020
    call restartvar(ncid=ncid, flag=flag,  varname='gdd020', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='20 year average of growing degree-days base 0C from planting', units='ddays', &
         interpinic_flag='interp', readvar=readvar, data=pps%gdd020)

    ! gddmaturity
    call restartvar(ncid=ncid, flag=flag,  varname='gddmaturity', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Growing degree days needed to harvest', units='ddays', &
         interpinic_flag='interp', readvar=readvar, data=pps%gddmaturity)

    ! huileaf
    call restartvar(ncid=ncid, flag=flag,  varname='huileaf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='heat unit index needed from planting to leaf emergence', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%huileaf)

    ! huigrain
    call restartvar(ncid=ncid, flag=flag,  varname='huigrain', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='heat unit index needed to reach vegetative maturity', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%huigrain)

    ! grainc
    call restartvar(ncid=ncid, flag=flag,  varname='grainc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C', units='gC/m2', &
         interpinic_flag='interp', readvar=readvar, data=pcs%grainc)

    ! grainc_storage
    call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C storage', units='gC/m2', &
         interpinic_flag='interp', readvar=readvar, data=pcs%grainc_storage)

    ! grainc_xfer
    call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C transfer', units='gC/m2', &
         interpinic_flag='interp', readvar=readvar, data=pcs%grainc_xfer)

    ! grainn
    call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N', units='gN/m2', &
         interpinic_flag='interp', readvar=readvar, data=pns%grainn)

    ! grainn_storage
    call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N storage', units='gN/m2', &
         interpinic_flag='interp', readvar=readvar, data=pns%grainn_storage)

    ! grainn_xfer 
    call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N transfer', units='gN/m2', &
         interpinic_flag='interp', readvar=readvar, data=pns%grainn_xfer)

    ! grainc_xfer_to_grainc
    call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer_to_grainc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C growth from storage', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%grainc_xfer_to_grainc)

    ! livestemc_to_litter
    call restartvar(ncid=ncid, flag=flag,  varname='livestemc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='live stem C litterfall', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%livestemc_to_litter)

    ! grainc_to_food
    call restartvar(ncid=ncid, flag=flag,  varname='grainc_to_food', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C to food', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%grainc_to_food)

    ! grainn_xfer_to_grainn
    call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer_to_grainn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N growth from storage', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%grainn_xfer_to_grainn)

    ! livestemn_to_litter
    call restartvar(ncid=ncid, flag=flag,  varname='livestemn_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='livestem N to litter', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%livestemn_to_litter)

    ! grainn_to_food
    call restartvar(ncid=ncid, flag=flag,  varname='grainn_to_food', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N to food', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%grainn_to_food)

    ! cpool_to_grainc
    call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='allocation to grain C', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%cpool_to_grainc)

    ! cpool_to_grainc_storage
    call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='allocation to grain C storage', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%cpool_to_grainc_storage)

    ! npool_to_grainn
    call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='allocation to grain N', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%npool_to_grainn)

    ! npool_to_grainn_storage
    call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='allocation to grain N storage', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%npool_to_grainn_storage)

    ! cpool_grain_gr
    call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_gr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain growth respiration', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%cpool_grain_gr)

    ! cpool_grain_storage_gr
    call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_storage_gr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain growth respiration to storage', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%cpool_grain_storage_gr)

    ! transfer_grain_gr
    call restartvar(ncid=ncid, flag=flag,  varname='transfer_grain_gr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain growth respiration from storage', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%transfer_grain_gr)

    ! grainc_storage_to_xfer
    call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage_to_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain C shift storage to transfer', units='gC/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pcf%grainc_storage_to_xfer)

    ! grainn_storage_to_xfer
    call restartvar(ncid=ncid, flag=flag, varname='grainn_storage_to_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='grain N shift storage to transfer', units='gN/m2/s', &
         interpinic_flag='interp', readvar=readvar, data=pnf%grainn_storage_to_xfer)

  end subroutine CropRest

  !-----------------------------------------------------------------------
  integer function CropRestYear ( )

    ! !DESCRIPTION: 
    ! Return the restart year for prognostic crop

     CropRestYear = restyear
  end function CropRestYear

  !-----------------------------------------------------------------------
  subroutine CropRestIncYear ()
    !
    ! !DESCRIPTION: 
    ! Increment the crop restart year, if appropriate
    !
    ! This routine should be called every time step, but only once per clump (to avoid
    ! inadvertently updating nyrs multiple times)
    !
    ! !USES:
    use clm_varpar       , only : crop_prog
    use clm_time_manager , only : get_curr_date, is_first_step
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer kyr   ! current year
    integer kmo   ! month of year  (1, ..., 12)
    integer kda   ! day of month   (1, ..., 31)
    integer mcsec ! seconds of day (0, ..., seconds/day)
    !-----------------------------------------------------------------------

    ! Update restyear only when running with prognostic crop
    if ( crop_prog )then
       ! Update restyear when it's the start of a new year - but don't do that at the
       ! very start of the run
       call get_curr_date (   kyr, kmo, kda, mcsec)
       if ((kmo == 1 .and. kda == 1 .and. mcsec == 0) .and. .not. is_first_step()) then
          restyear = restyear + 1
       end if
    end if

  end subroutine CropRestIncYear

  !-----------------------------------------------------------------------
  subroutine checkDates( )
    !
    ! !DESCRIPTION: 
    ! Make sure the dates are compatible. The date given to startup the model
    ! and the date on the restart file must be the same although years can be
    ! different. The dates need to be checked when the restart file is being
    ! read in for a startup or branch case (they are NOT allowed to be different
    ! for a restart case).
    !
    ! For the prognostic crop model the date of planting is tracked and growing
    ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
    ! messes up these bits of saved information.
    !
    ! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
    !
    ! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !-----------------------------------------------------------------------
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( masterproc ) &
            write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
               (stmon_day - stmon_day/100)
          call endrun(msg=' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
               ' and the date on the restart file needs to match (years can be different)'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

  end subroutine checkDates

end module CropRestMod
