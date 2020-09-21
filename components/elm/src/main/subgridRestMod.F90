module subgridRestMod

#include "shr_assert.h"

  !------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use decompMod          , only : bounds_type, BOUNDS_LEVEL_PROC, ldecomp
  use domainMod          , only : ldomain
  use clm_time_manager   , only : get_curr_date
  use elm_varcon         , only : nameg, namel, namec, namep
  use elm_varpar         , only : nlevsno
  use pio                , only : file_desc_t
  use ncdio_pio          , only : ncd_int, ncd_double
  use GetGlobalValuesMod , only : GetGlobalIndexArray
  use GridcellType       , only : grc_pp
  use LandunitType       , only : lun_pp
  use ColumnType         , only : col_pp
  use VegetationType     , only : veg_pp
  use perf_mod           , only : t_startf, t_stopf
  use restUtilMod

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: subgridRest                   ! handle restart of subgrid variables
  public :: subgridRest_check_consistency ! check consistency of variables read by subgridRest
  public :: subgridRest_read_cleanup      ! do cleanup of variables allocated when reading the restart file; should be called after subgridRest and subgridRest_check_consistency are complete

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: subgridRest_write_only     ! handle restart of subgrid variables that only need to be written, not read
  private :: subgridRest_write_and_read ! handle restart of subgrid variables that need to be read as well as written
  private :: save_old_weights

  ! !PRIVATE TYPES:
  real(r8), allocatable :: pft_wtlunit_before_rest_read(:)  ! veg_pp%wtlunit weights - saved values from before the restart read
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine subgridRest( bounds, ncid, flag )
    !
    ! !DESCRIPTION:
    ! Handle restart of subgrid variables
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='SubgridRest' ! subroutine name
    !------------------------------------------------------------------------

    if (flag /= 'read') then
       call t_startf('subgridRest_write')
       call subgridRest_write_only(bounds, ncid, flag)
       call t_stopf('subgridRest_write')
    end if

    call t_startf('subgridRest_write-read')
    call subgridRest_write_and_read(bounds, ncid, flag)
    call t_stopf('subgridRest_write-read')

  end subroutine subgridRest

  !-----------------------------------------------------------------------
  subroutine subgridRest_write_only(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart for variables that only need to be written, not read. This applies
    ! to variables that are time-constant and are only put on the restart file for the
    ! sake of having some additional metadata there.
    !
    ! Note that 'active' flags appear in this routine: they don't need to be read because
    ! they can be computed using other info on the restart file (particularly subgrid
    ! weights).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p,i             ! indices
    logical :: readvar               ! temporary
    real(r8), pointer :: rgarr(:)    ! temporary
    real(r8), pointer :: rlarr(:)    ! temporary
    real(r8), pointer :: rcarr(:)    ! temporary
    real(r8), pointer :: rparr(:)    ! temporary
    integer , pointer :: igarr(:)    ! temporary
    integer , pointer :: ilarr(:)    ! temporary
    integer , pointer :: icarr(:)    ! temporary
    integer , pointer :: iparr(:)    ! temporary
    
    character(len=*), parameter :: subname = 'subgridRest_write_only'
    !-----------------------------------------------------------------------
    
    !------------------------------------------------------------------
    ! Write gridcell info
    !------------------------------------------------------------------

    allocate(rgarr(bounds%begg:bounds%endg), igarr(bounds%begg:bounds%endg))

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lon', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell longitude', units='degrees_east',         &
         interpinic_flag='skip', readvar=readvar, data=grc_pp%londeg)

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lat', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell latitude', units='degrees_north',         &
         interpinic_flag='skip', readvar=readvar, data=grc_pp%latdeg)

    do g=bounds%begg,bounds%endg
       igarr(g)= mod(ldecomp%gdc2glo(g)-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_ixy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d longitude index of corresponding gridcell',     &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    do g=bounds%begg,bounds%endg
       igarr(g)= (ldecomp%gdc2glo(g) - 1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_jxy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d latitude index of corresponding gridcell',      &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    deallocate(rgarr,igarr)

    !------------------------------------------------------------------
    ! Write landunit info
    !------------------------------------------------------------------

    allocate(rlarr(bounds%begl:bounds%endl), ilarr(bounds%begl:bounds%endl))

    do l=bounds%begl,bounds%endl
       rlarr(l) = grc_pp%londeg(lun_pp%gridcell(l))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_lon', xtype=ncd_double,  &
         dim1name='landunit',                                                      &
         long_name='landunit longitude', units='degrees_east',                     &
         interpinic_flag='skip', readvar=readvar, data=rlarr)
 
    do l=bounds%begl,bounds%endl
       rlarr(l) = grc_pp%latdeg(lun_pp%gridcell(l))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_lat', xtype=ncd_double,  &
         dim1name='landunit',                                                      &
         long_name='landunit latitude', units='degrees_north',                     &
         interpinic_flag='skip', readvar=readvar, data=rlarr)

    do l=bounds%begl,bounds%endl
       ilarr(l) = mod(ldecomp%gdc2glo(lun_pp%gridcell(l))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_ixy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d longitude index of corresponding landunit',                 &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

     do l=bounds%begl,bounds%endl
        ilarr(l) = (ldecomp%gdc2glo(lun_pp%gridcell(l))-1)/ldomain%ni + 1
     enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_jxy', xtype=ncd_int,     &
         dim1name='landunit',                                                      &
         long_name='2d latitude index of corresponding landunit',                  &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    ilarr = GetGlobalIndexArray(lun_pp%gridcell(bounds%begl:bounds%endl), bounds%begl, bounds%endl, clmlevel=nameg)
    call restartvar(ncid=ncid, flag=flag, varname='land1d_gridcell_index', xtype=ncd_int, &
         dim1name='landunit',                                                             &
         long_name='gridcell index of corresponding landunit',                            &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    call restartvar(ncid=ncid, flag=flag, varname='land1d_ityplun', xtype=ncd_int, &
         dim1name='landunit',                                                      &
         long_name='landunit type (see global attributes)', units=' ',             &
         interpinic_flag='skip', readvar=readvar, data=lun_pp%itype)

    do l=bounds%begl,bounds%endl
       if (lun_pp%active(l)) then
          ilarr(l) = 1
       else
          ilarr(l) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='land1d_active', xtype=ncd_int,  &
         dim1name='landunit',                                                      &
         long_name='landunit active flag (1=active, 0=inactive)',                  &
         interpinic_flag='skip', readvar=readvar, data=ilarr)

    deallocate(rlarr, ilarr)

    !------------------------------------------------------------------
    ! Write column info
    !------------------------------------------------------------------

    allocate(rcarr(bounds%begc:bounds%endc), icarr(bounds%begc:bounds%endc))

    do c= bounds%begc, bounds%endc
       rcarr(c) = grc_pp%londeg(col_pp%gridcell(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_lon', xtype=ncd_double,   &
         dim1name='column',                                                         &
         long_name='column longitude', units='degrees_east',                        &
         interpinic_flag='skip', readvar=readvar, data=rcarr)

    do c= bounds%begc, bounds%endc
       rcarr(c) = grc_pp%latdeg(col_pp%gridcell(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_lat', xtype=ncd_double,   &
         dim1name='column',                                                         &
         long_name='column latitude', units='degrees_north',                        &
         interpinic_flag='skip', readvar=readvar, data=rcarr)

    do c= bounds%begc, bounds%endc
       icarr(c) = mod(ldecomp%gdc2glo(col_pp%gridcell(c))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ixy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d longitude index of corresponding column', units=' ',         &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    do c= bounds%begc, bounds%endc
       icarr(c) = (ldecomp%gdc2glo(col_pp%gridcell(c))-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_jxy', xtype=ncd_int,      &
         dim1name='column',                                                         &
         long_name='2d latitude index of corresponding column', units=' ',          &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    icarr = GetGlobalIndexArray(col_pp%gridcell(bounds%begc:bounds%endc), bounds%begc, bounds%endc, clmlevel=nameg)
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_gridcell_index', xtype=ncd_int, &
         dim1name='column',                                                               &
         long_name='gridcell index of corresponding column',                              &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    icarr = GetGlobalIndexArray(col_pp%landunit(bounds%begc:bounds%endc), bounds%begc, bounds%endc, clmlevel=namel)
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_landunit_index', xtype=ncd_int, &
         dim1name='column',                                                               &
         long_name='landunit index of corresponding column',                              &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    do c= bounds%begc, bounds%endc
       icarr(c) = lun_pp%itype(col_pp%landunit(c))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ityplun', xtype=ncd_int,  &
         dim1name='column',                                                         &
         long_name='column landunit type (see global attributes)', units=' ',       &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_ityp', xtype=ncd_int,     &
         dim1name='column',                                                         &
         long_name='column type (see global attributes)', units=' ',                &
         interpinic_flag='skip', readvar=readvar, data=col_pp%itype)

    do c=bounds%begc,bounds%endc
       if (col_pp%active(c)) then 
          icarr(c) = 1
       else
          icarr(c) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='cols1d_active', xtype=ncd_int,   &
         dim1name='column',                                                         &
         long_name='column active flag (1=active, 0=inactive)', units=' ',          &
         interpinic_flag='skip', readvar=readvar, data=icarr)

    deallocate(rcarr, icarr)

    !------------------------------------------------------------------
    ! Write pft info
    !------------------------------------------------------------------

    allocate(rparr(bounds%begp:bounds%endp), iparr(bounds%begp:bounds%endp))

    do p=bounds%begp,bounds%endp
       rparr(p) = grc_pp%londeg(veg_pp%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lon', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft longitude', units='degrees_east',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       rparr(p) = grc_pp%latdeg(veg_pp%gridcell(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_lat', xtype=ncd_double, &
         dim1name='pft',                                                          &
         long_name='pft latitude', units='degrees_north',                         &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = mod(ldecomp%gdc2glo(veg_pp%gridcell(p))-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ixy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d longitude index of corresponding pft', units='',        &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = (ldecomp%gdc2glo(veg_pp%gridcell(p))-1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_jxy', xtype=ncd_int, &
         dim1name='pft',                                                       &
         long_name='2d latitude index of corresponding pft', units='',         &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = GetGlobalIndexArray(veg_pp%gridcell(bounds%begp:bounds%endp), bounds%begp, bounds%endp, clmlevel=nameg)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_gridcell_index', xtype=ncd_int, &
         dim1name='pft',                                                                  &
         long_name='gridcell index of corresponding pft',                                 &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = GetGlobalIndexArray(veg_pp%landunit(bounds%begp:bounds%endp), bounds%begp, bounds%endp, clmlevel=namel)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_landunit_index', xtype=ncd_int, &
         dim1name='pft',                                                                  &
         long_name='landunit index of corresponding pft',                                 &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    iparr = GetGlobalIndexArray(veg_pp%column(bounds%begp:bounds%endp), bounds%begp, bounds%endp, clmlevel=namec)
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_column_index', xtype=ncd_int,   &
         dim1name='pft',                                                                  &
         long_name='column index of corresponding pft',                                   &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypveg', xtype=ncd_int,  &
         dim1name='pft',                                                            &
         long_name='pft vegetation type', units='',                                 &
         interpinic_flag='skip', readvar=readvar, data=veg_pp%itype)

    do p=bounds%begp,bounds%endp
       iparr(p) = col_pp%itype(veg_pp%column(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_itypcol', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft column type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       iparr(p) = lun_pp%itype(veg_pp%landunit(p))
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_ityplun', xtype=ncd_int, &
         dim1name='pft',                                                           &
         long_name='pft landunit type (see global attributes)', units='',          &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       if (veg_pp%active(p)) then
          iparr(p) = 1
       else
          iparr(p) = 0
       end if
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_active', xtype=ncd_int, &
         dim1name='pft',                                                          &
         long_name='pft active flag (1=active, 0=inactive)', units='',            &
         interpinic_flag='skip', readvar=readvar, data=iparr)

    do p=bounds%begp,bounds%endp
       c = veg_pp%column(p)
       rparr(p) = col_pp%glc_topo(c)
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='skip', readvar=readvar, data=rparr)

    deallocate(rparr, iparr)

  end subroutine subgridRest_write_only

  !-----------------------------------------------------------------------
  subroutine subgridRest_write_and_read(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    logical :: readvar              ! temporary
    real(r8), pointer :: temp2d(:,:) ! temporary for sno column variables
    
    character(len=*), parameter :: subname = 'subgridRest_write_and_read'
    !-----------------------------------------------------------------------
    
    if (flag == 'read') then
       call save_old_weights(bounds)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='land1d_wtxy', xtype=ncd_double, &
         dim1name='landunit',                                                      &
         long_name='landunit weight relative to corresponding gridcell',           &
         interpinic_flag='skip', readvar=readvar, data=lun_pp%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtxy', xtype=ncd_double,  &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding gridcell', units=' ',   &
         interpinic_flag='skip', readvar=readvar, data=col_pp%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_wtlnd', xtype=ncd_double, &
         dim1name='column',                                                         &
         long_name='column weight relative to corresponding landunit', units=' ',   &
         interpinic_flag='skip', readvar=readvar, data=col_pp%wtlunit)

    call restartvar(ncid=ncid, flag=flag, varname='cols1d_topoglc', xtype=ncd_double,   &
         dim1name='column',                                                             &
         long_name='mean elevation on glacier elevation classes', units='m',            &
         interpinic_flag='skip', readvar=readvar, data=col_pp%glc_topo)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtxy', xtype=ncd_double,  &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding gridcell', units='',       &  
         interpinic_flag='skip', readvar=readvar, data=veg_pp%wtgcell)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtlnd', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding landunit', units='',       & 
         interpinic_flag='skip', readvar=readvar, data=veg_pp%wtlunit)

    call restartvar(ncid=ncid, flag=flag, varname='pfts1d_wtcol', xtype=ncd_double, &
         dim1name='pft',                                                            &
         long_name='pft weight relative to corresponding column', units='',         &
         interpinic_flag='skip', readvar=readvar, data=veg_pp%wtcol)

    ! Snow column variables

    call restartvar(ncid=ncid, flag=flag, varname='SNLSNO', xtype=ncd_int,  & 
         dim1name='column', &
         long_name='number of snow layers', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=col_pp%snl)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = col_pp%dz(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='DZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer thickness', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col_pp%dz(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = col_pp%z(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col_pp%z(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno:-1))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno:-1) = col_pp%zi(bounds%begc:bounds%endc,-nlevsno:-1)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZISNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno, upperb2=-1, &
         long_name='snow interface depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       col_pp%zi(bounds%begc:bounds%endc,-nlevsno:-1) = temp2d(bounds%begc:bounds%endc,-nlevsno:-1) 
    end if
    deallocate(temp2d)

  end subroutine subgridRest_write_and_read

  !-----------------------------------------------------------------------
  subroutine save_old_weights(bounds)
    !
    ! !DESCRIPTION:
    ! Save old weights, from before the restart read, for later consistency checks.
    !
    ! !USES:
    type(bounds_type), intent(in)    :: bounds ! bounds (expected to be proc-level)
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'save_old_weights'
    !-----------------------------------------------------------------------
    
    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname//' ERROR: expect proc-level bounds')

    allocate(pft_wtlunit_before_rest_read(bounds%begp:bounds%endp))
    pft_wtlunit_before_rest_read(bounds%begp:bounds%endp) = veg_pp%wtlunit(bounds%begp:bounds%endp)

  end subroutine save_old_weights


  !-----------------------------------------------------------------------
  subroutine subgridRest_check_consistency(bounds)
    !
    ! !DESCRIPTION:
    ! Check consistency of variables read by subgridRest.
    !
    ! This should be called AFTER subgridRest is called to read the restart file.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgridRest_check_consistency'
    !-----------------------------------------------------------------------

    if (do_check_weights()) then
       call check_weights(bounds)
    end if

  contains

    !-----------------------------------------------------------------------
    logical function do_check_weights()
      !
      ! !DESCRIPTION:
      ! Return true if we should check weights
      !
      ! !USES:
      use elm_varctl          , only : nsrest, nsrContinue, use_fates
      use dynSubgridControlMod, only : get_do_transient_pfts
      !
      ! !ARGUMENTS:
      !
      ! !LOCAL VARIABLES:
      
      character(len=*), parameter :: subname = 'do_check_weights'
      !-----------------------------------------------------------------------
      
      if (get_do_transient_pfts()) then
         ! Don't check weights for a transient PFT case, because it's harder to come up with the
         ! correct weights to check against
         do_check_weights = .false.
      else if (nsrest == nsrContinue) then
         ! Don't check weights for a restart run
         !
         ! WJS (3-25-14): I'm not sure why we don't do the check in this case, but I'm
         ! maintaining the logic that used to be in BiogeophysRestMod regarding these
         ! weight checks
         do_check_weights = .false.
      else if (use_fates) then
         ! Don't check weights for a ed case, because the weights will almost certainly
         ! differ from the surface dataset in this case
         do_check_weights = .false.
      else
         do_check_weights = .true.
      end if

    end function do_check_weights

    !-----------------------------------------------------------------------
    subroutine check_weights(bounds)
      !
      ! !DESCRIPTION:
      ! Make sure that pft weights on the landunit agree with the weights read from the
      ! surface dataset, for the natural veg landunit.
      !
      ! Note that we do NOT do a more general check of all subgrid weights, because it's
      ! possible that some other subgrid weights have changed relative to the surface
      ! dataset, e.g., due to dynamic landunits. It would probably be possible to do more
      ! checking than is done here, but the check here should be sufficient to catch major
      ! inconsistencies between the restart file and the surface dataset.
      !
      ! !USES:
      use landunit_varcon, only : istsoil
      use elm_varctl, only : iulog
      !
      ! !ARGUMENTS:
      type(bounds_type), intent(in)    :: bounds ! bounds
      !
      ! !LOCAL VARIABLES:
      integer  :: p, l ! indices
      real(r8) :: diff ! difference in weights

      real(r8), parameter :: tol = 5.e-3  ! tolerance for checking weights
      
      character(len=*), parameter :: subname = 'check_weights'
      !-----------------------------------------------------------------------
      
      do p = bounds%begp, bounds%endp
         l = veg_pp%landunit(p)
         if (lun_pp%itype(l) == istsoil) then
            diff = abs(veg_pp%wtlunit(p) - pft_wtlunit_before_rest_read(p))
            if (diff > tol) then
               write(iulog,*) 'ERROR: PFT weights are SIGNIFICANTLY different between the restart (finidat) file'
               write(iulog,*) 'and the surface dataset (fsurdat).'
               write(iulog,*) 'Maximum allowed difference: ', tol
               write(iulog,*) 'Difference found: ', diff
               write(iulog,*) 'This match is a requirement for non-transient runs'
               write(iulog,*)
               write(iulog,*) 'Possible solutions to this problem:'
               write(iulog,*) '(1) Make sure you are using the intended finidat and fsurdat files'
               write(iulog,*) '(2) If you are running a present-day simulation, then make sure that your'
               write(iulog,*) '    initial conditions file is from the END of a 20th century transient run'
               write(iulog,*) '(3) If you are confident that you are using the correct finidat and fsurdat files,'
               write(iulog,*) '    yet are still experiencing this error, then you can bypass this check by setting:'
               write(iulog,*) '      check_finidat_pct_consistency = .false.'
               write(iulog,*) '    in user_nl_clm'
               write(iulog,*) '    In this case, CLM will take the weights from the initial conditions file.'
               write(iulog,*) ' '
               call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
            end if
         end if
      end do

    end subroutine check_weights

  end subroutine subgridRest_check_consistency


  !-----------------------------------------------------------------------
  subroutine subgridRest_read_cleanup
    !
    ! !DESCRIPTION:
    ! Do cleanup of variables allocated when reading the restart file
    !
    ! Should be called after subgridRest and subgridRest_check_consistency are complete.
    ! Note that this must be called after subgridRest is called to read the restart file,
    ! in order to avoid a memory leak.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'subgridRest_read_cleanup'
    !-----------------------------------------------------------------------
    
    deallocate(pft_wtlunit_before_rest_read)

  end subroutine subgridRest_read_cleanup


end module subgridRestMod
