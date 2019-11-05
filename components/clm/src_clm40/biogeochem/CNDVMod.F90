module CNDVMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDVMod
!
! !DESCRIPTION:
! Module containing routines to drive the annual dynamic vegetation
! that works with CN, reset related variables,
! and initialize/reset time invariant variables
!
! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use abortutils          , only : endrun
  use CNVegStructUpdateMod, only : CNVegStructUpdate
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public dv                 ! Drives the annual dynamic vegetation that
                            ! works with CN
  public histCNDV           ! Output CNDV history file
!
! !REVISION HISTORY:
! Module modified by Sam Levis from similar module DGVMMod
! created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dv
!
! !INTERFACE:
  subroutine dv(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp, kyr)
!
! !DESCRIPTION:
! Drives the annual dynamic vegetation that works with CN
!
! !USES:
    use clmtype
    use CNDVLightMod        , only : Light
    use CNDVEstablishmentMod, only : Establishment
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg       ! gridcell bounds
    integer, intent(in) :: lbp, ubp       ! pft bounds
    integer, intent(inout) :: num_natvegp ! number of naturally-vegetated
                                          ! pfts in filter
    integer, intent(inout) :: filter_natvegp(ubp-lbp+1) ! filter for
                                          ! naturally-vegetated pfts
    integer, intent(in) :: kyr            ! used in routine climate20 below
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: mxy(:)         ! pft m index (for laixy(i,j,m),etc.)
   integer , pointer :: pgridcell(:)   ! gridcell of corresponding pft
   real(r8), pointer :: fpcgrid(:)     ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: agdd(:)        ! accumulated growing degree days above 5
   real(r8), pointer :: t_mo_min(:)    ! annual min of t_mo (Kelvin)
!
! local pointers to implicit inout arguments
!
   real(r8), pointer :: tmomin20(:)         ! 20-yr running mean of tmomin
   real(r8), pointer :: agdd20(:)           ! 20-yr running mean of agdd
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p                    ! indices
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    agdd20    => gdgvs%agdd20
    tmomin20  => gdgvs%tmomin20

    ! Assign local pointers to derived type members (pft-level)

    mxy       => pft%mxy
    pgridcell => pft%gridcell
    fpcgrid   => pdgvs%fpcgrid
    t_mo_min  => pdgvs%t_mo_min
    agdd      => pdgvs%agdd

    ! *************************************************************************
    ! S. Levis version of LPJ's routine climate20: 'Returns' tmomin20 & agdd20
    ! for use in routine bioclim, which I have placed in routine Establishment
    ! Instead of 20-yr running mean of coldest monthly temperature,
    ! use 20-yr running mean of minimum 10-day running mean
    ! *************************************************************************

    do p = lbp,ubp
       g = pgridcell(p)
       if (kyr == 2) then ! slevis: add ".and. start_type==arb_ic" here?
          tmomin20(g) = t_mo_min(p) ! NO, b/c want to be able to start dgvm
          agdd20(g) = agdd(p)       ! w/ clmi file from non-dgvm simulation
       end if
       tmomin20(g) = (19._r8 * tmomin20(g) + t_mo_min(p)) / 20._r8
       agdd20(g)   = (19._r8 * agdd20(g)   + agdd(p)    ) / 20._r8
    end do

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns fpcgrid and nind

    call Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated fpcgrid, nind, crownarea, and present. Due to updated
    ! present, we do not use the natveg filter in this subroutine.

    call Establishment(lbg, ubg, lbp, ubp)

    ! Reset dgvm variables needed in next yr (too few to keep subr. dvreset)

    do p = lbp,ubp
       pcs%leafcmax(p) = 0._r8
       pdgvs%t_mo_min(p) = 1.0e+36_r8
    end do
  end subroutine dv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: histCNDV
!
! !INTERFACE:
  subroutine histCNDV()
!
! !DESCRIPTION:
! Create CNDV history file
!
! !USES:
    use clmtype
    use decompMod       , only : get_proc_bounds, get_proc_global
    use clm_varpar      , only : maxpatch_pft
    use domainMod       , only : ldomain
    use clm_varctl      , only : caseid, ctitle, finidat, fsurdat, fpftcon, iulog
    use clm_varcon      , only : spval
    use clm_time_manager, only : get_ref_date, get_nstep, get_curr_date, get_curr_time
    use fileutils       , only : get_filename
    use shr_sys_mod     , only : shr_sys_getenv
    use spmdMod         , only : masterproc
    use shr_const_mod   , only : SHR_CONST_CDAY
    use ncdio_pio
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   logical , pointer :: ifspecial(:)        ! true=>landunit is not vegetated (landunit-level)
   integer , pointer :: pgridcell(:)        ! gridcell index of corresponding pft (pft-level)
   integer , pointer :: plandunit(:)        ! landunit index of corresponding pft (pft-level)
   integer , pointer :: mxy(:)              ! pft m index (for laixy(i,j,m),etc.)
   real(r8), pointer :: fpcgrid(:)          ! foliar projective cover on gridcell (fraction)
   real(r8), pointer :: nind(:)             ! number of individuals (#/m**2)
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: dgvm_fn      ! dgvm history filename
    type(file_desc_t)  :: ncid         ! netcdf file id
    integer :: ncprec                  ! output precision
    integer :: g,p,l                   ! indices
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: ier                     ! error status
    integer :: mdcur, mscur, mcdate    ! outputs from get_curr_time
    integer :: yr,mon,day,mcsec        ! outputs from get_curr_date
    integer :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer :: nstep                   ! time step
    integer :: nbsec                   ! seconds components of a date
    integer :: dimid                   ! dimension, variable id
    real(r8):: time                    ! current time
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    real(r8), pointer :: rbuf2dg(:,:)  ! temporary
    character(len=32) :: subname='histCNDV'
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (gridcell-level)

    ! NONE

    ! Assign local pointers to derived type members (landunit-level)

    ifspecial  => lun%ifspecial

    ! Assign local pointers to derived subtypes components (pft-level)

    mxy       => pft%mxy
    pgridcell => pft%gridcell
    plandunit => pft%landunit
    fpcgrid   => pdgvs%fpcgrid
    nind      => pdgvs%nind

    ! Determine subgrid bounds for this processor and allocate dynamic memory

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    allocate(rbuf2dg(begg:endg,maxpatch_pft), stat=ier)
    if (ier /= 0) call endrun('histCNDV: allocation error for rbuf2dg')

    ! Set output precision

    ncprec = ncd_double

    ! -----------------------------------------------------------------------
    ! Create new netCDF file. File will be in define mode
    ! -----------------------------------------------------------------------

    dgvm_fn = set_dgvm_filename()
    call ncd_pio_createfile(ncid, trim(dgvm_fn))

    ! -----------------------------------------------------------------------
    ! Create global attributes.
    ! -----------------------------------------------------------------------
    
    str = 'CF1.0'
    call ncd_putatt (ncid, ncd_global, 'conventions', trim(str))
    
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, ncd_global,'history', trim(str))
    
    call shr_sys_getenv('LOGNAME', str, ier)
    if (ier /= 0) call endrun('error: LOGNAME environment variable not defined')
       
    call ncd_putatt (ncid, ncd_global, 'logname', trim(str))
       
    call shr_sys_getenv('HOST', str, ier)
    call ncd_putatt (ncid, ncd_global, 'host', trim(str))
       
    str = 'Community Land Model: CLM3'
    call ncd_putatt (ncid, ncd_global, 'source',  trim(str))
       
    str = '$Name$'
    call ncd_putatt (ncid, ncd_global, 'version', trim(str))
       
    str = '$Id$'
    call ncd_putatt (ncid, ncd_global, 'revision_id',  trim(str))

    str = ctitle
    call ncd_putatt (ncid, ncd_global, 'case_title', trim(str))

    str = caseid
    call ncd_putatt (ncid, ncd_global, 'case_id',  trim(str))

    str = get_filename(fsurdat)
    call ncd_putatt(ncid, ncd_global, 'Surface_dataset',  trim(str))

    str = 'arbitrary initialization'
    if (finidat /= ' ') str = get_filename(finidat)
    call ncd_putatt(ncid, ncd_global, 'Initial_conditions_dataset',  trim(str))

    str = get_filename(fpftcon)
    call ncd_putatt(ncid, ncd_global, 'PFT_physiological_constants_dataset', trim(str))

    ! -----------------------------------------------------------------------
    ! Define dimensions.
    ! -----------------------------------------------------------------------
    
    if (ldomain%isgrid2d) then
       call ncd_defdim (ncid, 'lon' ,ldomain%ni, dimid)
       call ncd_defdim (ncid, 'lat' ,ldomain%nj, dimid)
    else
       call ncd_defdim (ncid, 'gridcell', ldomain%ns, dimid)
    end if
    call ncd_defdim (ncid, 'pft' , maxpatch_pft , dimid)
    call ncd_defdim (ncid, 'time', ncd_unlimited, dimid)
    call ncd_defdim (ncid, 'string_length', 80  , dimid)
    
    ! -----------------------------------------------------------------------
    ! Define variables
    ! -----------------------------------------------------------------------
    
    ! Define coordinate variables (including time)
    
    if (ldomain%isgrid2d) then
       call ncd_defvar(ncid=ncid, varname='lon', xtype=ncprec, dim1name='lon', &
            long_name='coordinate longitude', units='degrees_east')
    
       call ncd_defvar(ncid=ncid, varname='lat', xtype=ncprec, dim1name='lat', &
            long_name='coordinate latitude', units='degrees_north')
    end if
    
    call get_curr_time(mdcur, mscur)
    call get_ref_date(yr, mon, day, nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,80) yr,mon,day
80  format(i4.4,'-',i2.2,'-',i2.2)
    write(basesec ,90) hours, minutes, secs
90  format(i2.2,':',i2.2,':',i2.2)
    str = 'days since ' // basedate // " " // basesec
    time = mdcur + mscur/SHR_CONST_CDAY
    
    call ncd_defvar(ncid=ncid, varname='time', xtype=ncd_double, dim1name='time', &
         long_name='time', units=str)
       
    ! Define surface grid (coordinate variables, latitude, longitude, surface type).
    
    if (ldomain%isgrid2d) then
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, &
            dim1name='lon', dim2name='lat', &
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, &
            dim1name='lon', dim2name='lat', &
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=ncd_int, &
            dim1name='lon', dim2name='lat', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
    else
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, &
            dim1name='gridcell',&
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, &
            dim1name='gridcell',&
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=ncd_int, &
            dim1name='gridcell', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
    end if

    ! Define time information

    call ncd_defvar(ncid=ncid, varname='mcdate', xtype=ncd_int, dim1name='time',&
         long_name='current date (YYYYMMDD)')

    call ncd_defvar(ncid=ncid, varname='mcsec', xtype=ncd_int, dim1name='time',&
         long_name='current seconds of current date', units='s')

    call ncd_defvar(ncid=ncid, varname='mdcur', xtype=ncd_int, dim1name='time',&
         long_name='current day (from base day)')

    call ncd_defvar(ncid=ncid, varname='mscur', xtype=ncd_int, dim1name='time',&
         long_name='current seconds of current day', units='s')

    call ncd_defvar(ncid=ncid, varname='nstep', xtype=ncd_int, dim1name='time',&
         long_name='time step', units='s')

    ! Define time dependent variables

    if (ldomain%isgrid2d) then
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)
       
       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
    else 
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='gridcell', dim2name='pft', dim3name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)
       
       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='gridcell', dim2name='pft', dim3name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
    end if

    call ncd_enddef(ncid)

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    ! Write surface grid (coordinate variables, latitude, longitude, surface type).

    call ncd_io(ncid=ncid, varname='longxy'  , data=ldomain%lonc, flag='write', &
         dim1name=grlnd)
    call ncd_io(ncid=ncid, varname='latixy'  , data=ldomain%latc, flag='write', &
         dim1name=grlnd)
    call ncd_io(ncid=ncid, varname='landmask', data=ldomain%mask, flag='write', &
         dim1name=grlnd)

    ! Write current date, current seconds, current day, current nstep

    call get_curr_date(yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day
    nstep = get_nstep()

    call ncd_io(ncid=ncid, varname='mcdate', data=mcdate, nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mcsec' , data=mcsec , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mdcur' , data=mdcur , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mscur' , data=mcsec , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='nstep' , data=nstep , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='time'  , data=time  , nt=1, flag='write')

    ! Write time dependent variables to CNDV history file

    ! The if .not. ifspecial statment below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    rbuf2dg(:,:) = 0._r8
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = fpcgrid(p)*100._r8
    end do
    call ncd_io(ncid=ncid, varname='FPCGRID', dim1name=grlnd, data=rbuf2dg, &
         nt=1, flag='write')

    rbuf2dg(:,:) = 0._r8
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = nind(p)
    end do
    call ncd_io(ncid=ncid, varname='NIND', dim1name=grlnd, data=rbuf2dg, &
         nt=1, flag='write')

    ! Deallocate dynamic memory

    deallocate(rbuf2dg)

    !------------------------------------------------------------------
    ! Close and archive netcdf CNDV history file
    !------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*)'(histCNDV): Finished writing CNDV history dataset ',&
            trim(dgvm_fn), 'at nstep = ',get_nstep()
    end if

  end subroutine histCNDV

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dgvm_filename
!
! !INTERFACE:
  character(len=256) function set_dgvm_filename ()
!
! !DESCRIPTION:
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl       , only : caseid, inst_suffix
    use clm_time_manager , only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_dgvm_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
                        ".hv."//trim(cdate)//".nc"

  end function set_dgvm_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildNatVegFilter
!
! !INTERFACE:
  subroutine BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbp, ubp                   ! pft bounds
    integer, intent(out) :: num_natvegp                ! number of pfts in naturally-vegetated filter
    integer, intent(out) :: filter_natvegp(ubp-lbp+1)  ! pft filter for naturally-vegetated points
!
! !CALLED FROM:
! subroutine lpj in this module
!
! !REVISION HISTORY:
! Author: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    logical , pointer :: present(:)     ! whether this pft present in patch
!EOP
!
! !LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type members (pft-level)
    present   => pdgvs%present

    num_natvegp = 0
    do p = lbp,ubp
       if (present(p)) then
          num_natvegp = num_natvegp + 1
          filter_natvegp(num_natvegp) = p
       end if
    end do

  end subroutine BuildNatVegFilter

end module CNDVMod
