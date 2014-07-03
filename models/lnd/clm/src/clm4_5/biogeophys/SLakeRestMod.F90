module SLakeRestMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SLakeRestMod
!
! !DESCRIPTION:
! Reads from or writes restart data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeRest
!
! !REVISION HISTORY:
! 2009, June: Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SLakeRest
!
! !INTERFACE:
  subroutine SLakeRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio_pio
    use decompMod     , only : get_proc_bounds
    use clm_time_manager , only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 12/11/2003, Peter Thornton: Added cps%coszen, pps%gdir, and pps%omega
!   for new sunlit/shaded canopy algorithm (in SUNSHA ifdef block)
! 4/25/2005, Peter Thornton: Removed the SUNSHA ifdefs, since this is now the
!   default code behavior.
! 6/12/2005, Moved to netcdf format and renamed file
! 6/2009, Zack Subin: Adapted for S Lake physics.
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Note t_lake is already in BiogeophysRest.

    ! column water state variable - lake_icefrac

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='LAKE_ICEFRAC', xtype=ncd_double, &
            dim1name='column', dim2name='levlak', switchdim=.true., &
            long_name='lake layer ice fraction', units='kg/kg')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='LAKE_ICEFRAC', data=cptr%cws%lake_icefrac, &
            dim1name='column', switchdim=.true., &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - savedtke1

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='SAVEDTKE1', xtype=ncd_double,  &
            dim1name='column', &
            long_name='top lake layer eddy conductivity', units='W/(m K)')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='SAVEDTKE1', data=cptr%cps%savedtke1, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - ust_lake

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='USTLAKE', xtype=ncd_double,  &
            dim1name='column', &
            long_name='friction velocity for lakes', units='m/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='USTLAKE', data=cptr%cps%ust_lake, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! column physical state variable - z0mg

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='Z0MG', xtype=ncd_double,  &
            dim1name='column', &
            long_name='ground momentum roughness length', units='m')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='Z0MG', data=cptr%cps%z0mg, &
            dim1name='column', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag == 'read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

  end subroutine SLakeRest

end module SLakeRestMod
