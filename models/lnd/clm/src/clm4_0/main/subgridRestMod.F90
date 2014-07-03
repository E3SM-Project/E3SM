module subgridRestMod

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subgridRest
!
! !INTERFACE:
  subroutine subgridRest( ncid, flag )

    use shr_kind_mod , only : r8 => shr_kind_r8
    use clmtype
    use ncdio_pio           
    use decompMod        , only : get_proc_bounds, ldecomp
    use domainMod        , only : ldomain
    use clm_time_manager , only : get_curr_date
    use abortutils       , only : endrun
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid  ! netCDF dataset id
    character(len=*) , intent(in)    :: flag  ! flag to determine if define, write or read data
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: g,l,c,p,j,i         ! indices
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    integer :: begp, endp          ! per-proc beg/end pft indices
    integer :: begc, endc          ! per-proc beg/end column indices
    integer :: begl, endl          ! per-proc beg/end landunit indices
    integer :: begg, endg          ! per-proc beg/end gridcell indices
    integer :: ier                 ! error status
    real(r8),pointer :: rgarr(:)   ! temporary
    real(r8),pointer :: rlarr(:)   ! temporary
    real(r8),pointer :: rcarr(:)   ! temporary
    real(r8),pointer :: rparr(:)   ! temporary
    integer ,pointer :: igarr(:)   ! temporary
    integer ,pointer :: ilarr(:)   ! temporary
    integer ,pointer :: icarr(:)   ! temporary
    integer ,pointer :: iparr(:)   ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    character(len=32) :: subname='SubgridRest' ! subroutine name
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => grc
    lptr => lun
    cptr => col
    pptr => pft

    ! Get relevant sizes

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    ! Allocate dynamic memory

    if (flag == 'write') then
       allocate(rgarr(begg:endg),rlarr(begl:endl),rcarr(begc:endc),rparr(begp:endp),stat=ier)
       if (ier /= 0) call endrun('allocation error from inicfile_fields rarrs')
       allocate(igarr(begg:endg),ilarr(begl:endl),icarr(begc:endc),iparr(begp:endp),stat=ier)
       if (ier /= 0) call endrun('allocation error from inicfile_fields iarrs')
    end if

    ! Write output data (first write current date and seconds of current date)

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='mcdate', xtype=ncd_int, &
            long_name='current date as 8 digit integer (YYYYMMDD)')
       call ncd_defvar(ncid=ncid, varname='mcsec', xtype=ncd_int,  &
            long_name='current seconds of current date', units='s')
    else if (flag == 'write') then
       call get_curr_date (yr, mon, day, mcsec)
       mcdate = yr*10000 + mon*100 + day
       !TODO - add this to the file - get this to work
!DEBUG  call ncd_io(varname='mcdate', data=mcdate, ncid=ncid, flag=flag)
!DEBUG  call ncd_io(varname='mcsec' , data=mcsec , ncid=ncid, flag=flag)
    end if

    ! Write gridcell info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='grid1d_lon', xtype=ncd_double,  &
            dim1name='gridcell', long_name='gridcell longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='grid1d_lat', xtype=ncd_double,  &
            dim1name='gridcell', long_name='gridcell latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='grid1d_ixy', xtype=ncd_int,  &
            dim1name='gridcell', long_name='2d longitude index of corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='grid1d_jxy', xtype=ncd_int,  &
            dim1name='gridcell', long_name='2d latitude index of corresponding gridcell')
    else if (flag == 'write') then
       do g=begg,endg
          igarr(g)= mod(ldecomp%gdc2glo(g)-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='grid1d_ixy', data=igarr      , dim1name=nameg, ncid=ncid, flag=flag)
       do g=begg,endg
          igarr(g)= (ldecomp%gdc2glo(g) - 1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='grid1d_jxy', data=igarr      , dim1name=nameg, ncid=ncid, flag=flag)
       call ncd_io(varname='grid1d_lon', data=gptr%londeg, dim1name=nameg, ncid=ncid, flag=flag)
       call ncd_io(varname='grid1d_lat', data=gptr%latdeg, dim1name=nameg, ncid=ncid, flag=flag)
    end if

    ! Write landunit info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='land1d_lon', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='land1d_lat', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='land1d_ixy', xtype=ncd_int,  &
            dim1name='landunit', long_name='2d longitude index of corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='land1d_jxy', xtype=ncd_int,  &
            dim1name='landunit', long_name='2d latitude index of corresponding landunit')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_defvar(ncid=ncid, varname='land1d_gi', xtype=ncd_int,  &
       !     dim1name='landunit', long_name='1d grid index of corresponding landunit')
       ! ----------------------------------------------------------------
       call ncd_defvar(ncid=ncid, varname='land1d_wtxy', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='land1d_ityplun', xtype=ncd_int,  &
            dim1name='landunit', long_name='landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       do l=begl,endl
          rlarr(l) = gptr%londeg(lptr%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lon'    , data=rlarr        , dim1name=namel, ncid=ncid, flag=flag)
       do l=begl,endl
          rlarr(l) = gptr%latdeg(lptr%gridcell(l))
       enddo
       call ncd_io(varname='land1d_lat'    , data=rlarr        , dim1name=namel, ncid=ncid, flag=flag)
       do l=begl,endl
          ilarr(l) = mod(ldecomp%gdc2glo(lptr%gridcell(l))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='land1d_ixy'    , data=ilarr        , dim1name=namel, ncid=ncid, flag=flag)
       do l=begl,endl
          ilarr(l) = (ldecomp%gdc2glo(lptr%gridcell(l))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='land1d_jxy'    , data=ilarr        , dim1name=namel, ncid=ncid, flag=flag)
       call ncd_io(varname='land1d_wtxy'   , data=lptr%wtgcell , dim1name=namel, ncid=ncid, flag=flag)
       call ncd_io(varname='land1d_ityplun', data=lptr%itype   , dim1name=namel, ncid=ncid, flag=flag)
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_io(varname='land1d_gi'     , data=lptr%gridcell, dim1name=namel, ncid=ncid, flag=flag)
       ! ----------------------------------------------------------------
    end if

    ! Write column info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cols1d_lon', xtype=ncd_double,  &
            dim1name='column', long_name='column longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='cols1d_lat', xtype=ncd_double,  &
            dim1name='column', long_name='column latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='cols1d_ixy', xtype=ncd_int,   &
            dim1name='column', long_name='2d longitude index of corresponding column')
       call ncd_defvar(ncid=ncid, varname='cols1d_jxy', xtype=ncd_int,   &
            dim1name='column', long_name='2d latitude index of corresponding column')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_defvar(ncid=ncid, varname='cols1d_gi', xtype=ncd_int,   &
       !     dim1name='column', long_name='1d grid index of corresponding column')
       !call ncd_defvar(ncid=ncid, varname='cols1d_li', xtype=ncd_int,   &
       !     dim1name='column', long_name='1d landunit index of corresponding column')
       ! ----------------------------------------------------------------
       call ncd_defvar(ncid=ncid, varname='cols1d_wtxy', xtype=ncd_double,   &
            dim1name='column', long_name='column weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='cols1d_wtlnd', xtype=ncd_double,   &
            dim1name='column', long_name='column weight relative to corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='cols1d_ityplun', xtype=ncd_int,   &
            dim1name='column', long_name='column landunit type (vegetated,urban,lake,wetland or glacier)')
       call ncd_defvar(ncid=ncid, varname='cols1d_ityp', xtype=ncd_int,   &
            dim1name='column', long_name=&
           'column type (61-roof,62-sunwall,63-shadewall,64-impervious road,65-pervious road,1-all other columns)')
    else if (flag == 'write') then
       do c=begc,endc
          rcarr(c) = gptr%londeg(cptr%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lon'  , data=rcarr        , dim1name=namec, ncid=ncid, flag=flag)
       do c=begc,endc
          rcarr(c) = gptr%latdeg(cptr%gridcell(c))
       enddo
       call ncd_io(varname='cols1d_lat'  , data=rcarr        , dim1name=namec, ncid=ncid, flag=flag)
       do c=begc,endc
          icarr(c) = mod(ldecomp%gdc2glo(cptr%gridcell(c))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='cols1d_ixy'  , data=icarr        , dim1name=namec, ncid=ncid, flag=flag)
       do c=begc,endc
          icarr(c) = (ldecomp%gdc2glo(cptr%gridcell(c))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='cols1d_jxy'  , data=icarr        , dim1name=namec, ncid=ncid, flag=flag)
       call ncd_io(varname='cols1d_wtxy' , data=cptr%wtgcell , dim1name=namec, ncid=ncid, flag=flag)
       call ncd_io(varname='cols1d_wtlnd', data=cptr%wtlunit , dim1name=namec, ncid=ncid, flag=flag)
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_io(varname='cols1d_gi'   , data=cptr%gridcell, dim1name=namec, ncid=ncid, flag=flag)
       !call ncd_io(varname='cols1d_li'   , data=cptr%landunit, dim1name=namec, ncid=ncid, flag=flag)
       ! ----------------------------------------------------------------
       do c=begc,endc
          icarr(c) = lptr%itype(cptr%landunit(c))
       enddo
       call ncd_io(varname='cols1d_ityplun', data=icarr      , dim1name=namec, ncid=ncid, flag=flag)
       do c=begc,endc
          icarr(c) = cptr%itype((c))
       enddo
       call ncd_io(varname='cols1d_ityp', data=icarr      , dim1name=namec, ncid=ncid, flag=flag)
    end if

    ! Write pft info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pfts1d_lon', xtype=ncd_double,  &
            dim1name='pft', long_name='pft longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='pfts1d_lat', xtype=ncd_double,  &
            dim1name='pft', long_name='pft latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ixy', xtype=ncd_int,  &
            dim1name='pft', long_name='2d longitude index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_jxy', xtype=ncd_int,  &
            dim1name='pft', long_name='2d latitude index of corresponding pft')
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_defvar(ncid=ncid, varname='pfts1d_gi', xtype=ncd_int,  &
       !     dim1name='pft', long_name='1d grid index of corresponding pft')
       !call ncd_defvar(ncid=ncid, varname='pfts1d_li', xtype=ncd_int,  &
       !     dim1name='pft', long_name='1d landunit index of corresponding pft')
       ! ----------------------------------------------------------------
       call ncd_defvar(ncid=ncid, varname='pfts1d_ci', xtype=ncd_int,  &
            dim1name='pft', long_name='1d column index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtxy', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtlnd', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtcol', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding column')
       call ncd_defvar(ncid=ncid, varname='pfts1d_itypveg', xtype=ncd_int,  &
            dim1name='pft', long_name='pft vegetation type')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ityplun', xtype=ncd_int,  &
            dim1name='pft', long_name='pft landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       do p=begp,endp
          rparr(p) = gptr%londeg(pptr%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lon'    , data=rparr        , dim1name=namep, ncid=ncid, flag=flag)
       do p=begp,endp
          rparr(p) = gptr%latdeg(pptr%gridcell(p))
       enddo
       call ncd_io(varname='pfts1d_lat'    , data=rparr        , dim1name=namep, ncid=ncid, flag=flag)
       do p=begp,endp
          iparr(p) = mod(ldecomp%gdc2glo(pptr%gridcell(p))-1,ldomain%ni) + 1
       enddo
       call ncd_io(varname='pfts1d_ixy'    , data=iparr        , dim1name=namep, ncid=ncid, flag=flag)
       do p=begp,endp
          iparr(p) = (ldecomp%gdc2glo(pptr%gridcell(p))-1)/ldomain%ni + 1
       enddo
       call ncd_io(varname='pfts1d_jxy'    , data=iparr        , dim1name=namep, ncid=ncid, flag=flag)
       call ncd_io(varname='pfts1d_wtxy'   , data=pptr%wtgcell , dim1name=namep, ncid=ncid, flag=flag)
       call ncd_io(varname='pfts1d_wtlnd'  , data=pptr%wtlunit , dim1name=namep, ncid=ncid, flag=flag)
       call ncd_io(varname='pfts1d_wtcol'  , data=pptr%wtcol   , dim1name=namep, ncid=ncid, flag=flag)
       call ncd_io(varname='pfts1d_itypveg', data=pptr%itype   , dim1name=namep, ncid=ncid, flag=flag)
       ! --- EBK Do NOT write out indices that are incorrect 4/1/2011 --- Bug 1310
       !call ncd_io(varname='pfts1d_gi'     , data=pptr%gridcell, dim1name=namep, ncid=ncid, flag=flag)
       !call ncd_io(varname='pfts1d_li'     , data=pptr%landunit, dim1name=namep, ncid=ncid, flag=flag)
       !call ncd_io(varname='pfts1d_ci'     , data=pptr%column  , dim1name=namep, ncid=ncid, flag=flag)
       ! ----------------------------------------------------------------
       do p=begp,endp
          iparr(p) = lptr%itype(pptr%landunit(p))
       enddo
       call ncd_io(varname='pfts1d_ityplun', data=iparr      , dim1name=namep, ncid=ncid, flag=flag)
    end if

    if (flag == 'write') then
       deallocate(rgarr,rlarr,rcarr,rparr)
       deallocate(igarr,ilarr,icarr,iparr)
    end if

  end subroutine subgridRest

end module subgridRestMod
