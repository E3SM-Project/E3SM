module SLakeRestMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes restart data
  !
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  ! save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeRest
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SLakeRest( ncid, flag )
    !
    ! USES
    use clmtype
    use pio,         only: file_desc_t
    use ncdio_pio,   only: ncd_double 
    use restUtilMod
    !
    ! DESCRIPTION:
    ! Read/Write biogeophysics information to/from restart file.
    !
    ! arguments:
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'read' or 'write'
    !
    ! local variables:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    ! Note t_lake is already in BiogeophysRest.

    ! column water state variable - lake_icefrac
    call restartvar(ncid=ncid, flag=flag, varname='LAKE_ICEFRAC', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake layer ice fraction', units='kg/kg', &
         interpinic_flag='interp', readvar=readvar, data=cws%lake_icefrac)

    ! column physical state variable - savedtke1
    call restartvar(ncid=ncid, flag=flag, varname='SAVEDTKE1', xtype=ncd_double,  &
         dim1name='column', &
         long_name='top lake layer eddy conductivity', units='W/(m K)', &
         interpinic_flag='interp', readvar=readvar, data=cps%savedtke1)

    ! column physical state variable - ust_lake
    call restartvar(ncid=ncid, flag=flag, varname='USTLAKE', xtype=ncd_double,  &
         dim1name='column', &
         long_name='friction velocity for lakes', units='m/s', &
         interpinic_flag='interp', readvar=readvar, data=cps%ust_lake)

    ! column physical state variable - z0mg
    call restartvar(ncid=ncid, flag=flag, varname='Z0MG', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground momentum roughness length', units='m', &
         interpinic_flag='interp', readvar=readvar, data=cps%z0mg)

  end subroutine SLakeRest

end module SLakeRestMod
