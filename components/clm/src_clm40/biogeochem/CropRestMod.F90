module CropRestMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CropRestMod
! 
! !DESCRIPTION: 
! Read/Write to/from Crop info to CLM restart file. 
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
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
! !REVISION HISTORY:
! Module created by slevis following CNRestMod by Peter Thornton
!

! !PRIVATE DATA MEMBERS:
  integer :: restyear = 0         ! Restart year from the initial conditions file, incremented as time elapses

!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRest
!
! !INTERFACE:
  subroutine CropRest ( ncid, flag )
!
! !DESCRIPTION: 
! Read/write Crop restart data
!
! !USES:
    use clmtype
    use clm_atmlnd      , only : clm_a2l
    use clm_varpar      , only : numrad
    use decompMod       , only : get_proc_bounds
    use clm_time_manager, only : is_restart
    use ncdio_pio
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t)  :: ncid             ! netcdf id
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: slevis
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,p,j                      ! indices 
    integer :: begp, endp                 ! per-proc beginning and ending pft indices
    integer :: begc, endc                 ! per-proc beginning and ending column indices 
    integer :: begl, endl                 ! per-proc beginning and ending landunit indices
    integer :: begg, endg                 ! per-proc gridcell ending gridcell indices
    real(r8):: m                          ! multiplier for the exit_spinup code
    logical :: readvar                    ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    integer , pointer :: iptemp(:)        ! pointer to memory to be allocated
    integer :: ier                        ! error status
!-----------------------------------------------------------------------

    ! Prognostic crop restart year
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='restyear', xtype=ncd_int,  &
            long_name='Number of years prognostic crop ran', units="years")
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='restyear', data=restyear, &
            ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' )then
           if ( readvar ) then
              call checkDates( )
           else
              if ( is_restart()) call endrun
           end if
       end if       
    end if

    ! Set pointers into derived type

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !--------------------------------
    ! pft physical state variables 
    !--------------------------------

    ! peaklai
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='peaklai', xtype=ncd_int,  &
            dim1name='pft',long_name='Flag if at max allowed LAI or not', &
            flag_values=(/0,1/), nvalid_range=(/0,1/),                    &
            flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='peaklai', data=pps%peaklai, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! idop
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='idop', xtype=ncd_int,  &
            dim1name='pft',long_name='Date of planting',units='jday', &
            nvalid_range=(/1,366/) )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='idop', data=pps%idop, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! aleaf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='aleaf', xtype=ncd_double,  &
            dim1name='pft',long_name='leaf allocation coefficient',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='aleaf', data=pps%aleaf, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! aleafi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='aleafi', xtype=ncd_double,  &
            dim1name='pft',long_name='Saved leaf allocation coefficient from phase 2', &
            units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='aleafi', data=pps%aleafi, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! astem
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='astem', xtype=ncd_double,  &
            dim1name='pft',long_name='stem allocation coefficient',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='astem', data=pps%astem, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! astemi
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='astemi', xtype=ncd_double,  &
            dim1name='pft',long_name='Saved stem allocation coefficient from phase 2',&
            units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='astemi', data=pps%astemi, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! htmx 
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='htmx', xtype=ncd_double,  &
            dim1name='pft',long_name='max height attained by a crop during year',&
            units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='htmx', data=pps%htmx, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! hdidx
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='hdidx', xtype=ncd_double,  &
            dim1name='pft',long_name='cold hardening index',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='hdidx', data=pps%hdidx, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! vf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='vf', xtype=ncd_double,  &
            dim1name='pft',long_name='vernalization factor',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='vf', data=pps%vf, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cumvd
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cumvd', xtype=ncd_double,  &
            dim1name='pft',long_name='cumulative vernalization d',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cumvd', data=pps%cumvd, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! croplive
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='croplive', xtype=ncd_log,  &
            dim1name='pft',long_name='Flag that crop is alive, but not harvested')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='croplive', data=pps%croplive, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cropplant
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cropplant', xtype=ncd_log,  &
            dim1name='pft',long_name='Flag that crop is planted, but not harvested' )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cropplant', data=pps%cropplant, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! harvdate
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='harvdate', xtype=ncd_int,  &
            dim1name='pft',long_name='harvest date',units='jday', &
            nvalid_range=(/1,366/) )
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='harvdate', data=pps%harvdate, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! gdd1020
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gdd1020', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='20 year average of growing degree-days base 10C from planting', &
            units='ddays')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gdd1020', data=pps%gdd1020, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! gdd820
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gdd820', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='20 year average of growing degree-days base 8C from planting', &
            units='ddays')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gdd820', data=pps%gdd820, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! gdd020
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gdd020', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='20 year average of growing degree-days base 0C from planting', &
            units='ddays')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gdd020', data=pps%gdd020, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! gddmaturity
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='gddmaturity', xtype=ncd_double,  &
            dim1name='pft',long_name='Growing degree days needed to harvest',units='ddays')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='gddmaturity', data=pps%gddmaturity, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! huileaf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='huileaf', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='heat unit index needed from planting to leaf emergence',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='huileaf', data=pps%huileaf, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! huigrain
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='huigrain', xtype=ncd_double,  &
            dim1name='pft',long_name='heat unit index needed to reach vegetative maturity', &
            units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='huigrain', data=pps%huigrain, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C',units='gC/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc', data=pcs%grainc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C storage',units='gC/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc_storage', data=pcs%grainc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C transfer',units='gC/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc_xfer', data=pcs%grainc_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N',units='gN/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn', data=pns%grainn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N storage',units='gN/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn_storage', data=pns%grainn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N transfer',units='gN/m2')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn_xfer', data=pns%grainn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc_xfer_to_grainc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc_xfer_to_grainc', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C growth from storage',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc_xfer_to_grainc', data=pcf%grainc_xfer_to_grainc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! livestemc_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemc_to_litter', xtype=ncd_double,  &
            dim1name='pft',long_name='live stem C litterfall',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemc_to_litter', data=pcf%livestemc_to_litter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc_to_food
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc_to_food', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C to food',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc_to_food', data=pcf%grainc_to_food, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn_xfer_to_grainn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn_xfer_to_grainn', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N growth from storage',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn_xfer_to_grainn', data=pnf%grainn_xfer_to_grainn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! livestemn_to_litter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_to_litter', xtype=ncd_double,  &
            dim1name='pft',long_name='livestem N to litter',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn_to_litter', data=pnf%livestemn_to_litter, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn_to_food
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn_to_food', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N to food',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn_to_food', data=pnf%grainn_to_food, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cpool_to_grainc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_grainc', xtype=ncd_double,  &
            dim1name='pft',long_name='allocation to grain C',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cpool_to_grainc', data=pcf%cpool_to_grainc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cpool_to_grainc_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_to_grainc_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='allocation to grain C storage',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cpool_to_grainc_storage', data=pcf%cpool_to_grainc_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! npool_to_grainn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_grainn', xtype=ncd_double,  &
            dim1name='pft',long_name='allocation to grain N',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='npool_to_grainn', data=pnf%npool_to_grainn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! npool_to_grainn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool_to_grainn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='allocation to grain N storage',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='npool_to_grainn_storage', data=pnf%npool_to_grainn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cpool_grain_gr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_grain_gr', xtype=ncd_double,  &
            dim1name='pft',long_name='grain growth respiration',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cpool_grain_gr', data=pcf%cpool_grain_gr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! cpool_grain_storage_gr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cpool_grain_storage_gr', xtype=ncd_double,  &
            dim1name='pft',long_name='grain growth respiration to storage',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cpool_grain_storage_gr', data=pcf%cpool_grain_storage_gr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! transfer_grain_gr
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='transfer_grain_gr', xtype=ncd_double,  &
            dim1name='pft',long_name='grain growth respiration from storage',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='transfer_grain_gr', data=pcf%transfer_grain_gr, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainc_storage_to_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainc_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='grain C shift storage to transfer',units='gC/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainc_storage_to_xfer', data=pcf%grainc_storage_to_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

    ! grainn_storage_to_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grainn_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='grain N shift storage to transfer',units='gN/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='grainn_storage_to_xfer', data=pnf%grainn_storage_to_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar) 
       if (flag=='read' .and. .not. readvar) then
           if (is_restart()) call endrun
       end if       
    end if

  end subroutine CropRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestYear
!
! !INTERFACE:
  integer function CropRestYear ( )
!
! !DESCRIPTION: 
! Return the restart year for prognostic crop
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
     CropRestYear = restyear
  end function CropRestYear

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestIncYear
!
! !INTERFACE:
  subroutine CropRestIncYear ()
!
! !DESCRIPTION: 
! Increment the crop restart year, if appropriate
!
! This routine should be called every time step, but only once per clump (to avoid
! inadvertently updating nyrs multiple times)
!
! !USES:
    use surfrdMod        , only : crop_prog
    use clm_time_manager , only : get_curr_date, is_first_step
    implicit none
!
! !LOCAL VARIABLES:
    integer kyr                     ! current year
    integer kmo                     !         month of year  (1, ..., 12)
    integer kda                     !         day of month   (1, ..., 31)
    integer mcsec                   !         seconds of day (0, ..., seconds/day)
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
!BOP
!
! !IROUTINE: checkDates
!
! !INTERFACE:
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
! !USES:
!
! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
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
          call endrun( trim(subname)// &
          ' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
          ' and the date on the restart file needs to match (years can be different)' )
       end if
    end if

  end subroutine checkDates

end module CropRestMod

