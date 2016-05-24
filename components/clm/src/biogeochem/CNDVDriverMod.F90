module CNDVDriverMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Note that this module was created simply to contain the subroutine dv
  ! which cannot cannot be in CNDVMod due to circular dependencies 
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use abortutils        , only : endrun
  use decompMod         , only : bounds_type
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use atm2lndType       , only : atm2lnd_type
  use CNDVType          , only : dgvs_type
  use CNCarbonStateType , only : carbonstate_type
  use CNCarbonFluxType  , only : carbonflux_type
  use clm_varcon        , only : grlnd
  use LandunitType      , only : lun                
  use PatchType         , only : pft                
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNDVDriver
  public :: CNDVHist
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private set_dgvm_filename
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNDVDriver(bounds, &
       num_natvegp, filter_natvegp, kyr, &
       atm2lnd_vars, carbonflux_vars, carbonstate_vars, dgvs_vars)
    !
    ! !DESCRIPTION:
    ! Drives the annual dynamic vegetation that works with CN
    !
    ! !USES:
    use CNDVLightMod         , only : Light
    use CNDVEstablishmentMod , only : Establishment
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds                  
    integer                , intent(inout) :: num_natvegp             ! number of naturally-vegetated patches in filter
    integer                , intent(inout) :: filter_natvegp(:)       ! filter for naturally-vegetated patches
    integer                , intent(in)    :: kyr                     ! used in routine climate20 below
    type(atm2lnd_type)     , intent(inout) :: atm2lnd_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(carbonstate_type) , intent(inout) :: carbonstate_vars
    type(dgvs_type)        , intent(inout) :: dgvs_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: p                    ! patch index
    !-----------------------------------------------------------------------

    associate(                                            & 
         fpcgrid     => dgvs_vars%fpcgrid_patch         , & ! Input:  [real(r8) (:) ]  foliar projective cover on gridcell (fraction)    
         agdd20      => dgvs_vars%agdd20_patch          , & ! Output: [real(r8) (:) ]  20-yr running mean of agdd                        
         tmomin20    => dgvs_vars%tmomin20_patch        , & ! Output: [real(r8) (:) ]  20-yr running mean of tmomin                      

         t_mo_min    => atm2lnd_vars%t_mo_min_patch     , & ! Output: [real(r8) (:) ]  annual min of t_mo (Kelvin)                       

         agdd        => dgvs_vars%agdd_patch            , & ! Input:  [real(r8) (:) ]  accumulated growing degree days above 5           
         leafcmax    => carbonstate_vars%leafcmax_patch   & ! Output: [real(r8) (:) ]  (gC/m2) ann max leaf C 
         )

      ! *************************************************************************
      ! S. Levis version of LPJ's routine climate20: 'Returns' tmomin20 & agdd20
      ! for use in routine bioclim, which I have placed in routine Establishment
      ! Instead of 20-yr running mean of coldest monthly temperature,
      ! use 20-yr running mean of minimum 10-day running mean
      ! *************************************************************************
      
      do p = bounds%begp, bounds%endp
         if (kyr == 2) then ! slevis: add ".and. start_type==arb_ic" here?
            tmomin20(p) = t_mo_min(p) ! NO, b/c want to be able to start dgvm
            agdd20(p) = agdd(p)       ! w/ clmi file from non-dgvm simulation
         end if
         tmomin20(p) = (19._r8 * tmomin20(p) + t_mo_min(p)) / 20._r8
         agdd20(p)   = (19._r8 * agdd20(p)   + agdd(p)    ) / 20._r8
      end do

      ! Rebuild filter of present natually-vegetated patches after Kill()

      num_natvegp = 0
      do p = bounds%begp,bounds%endp
         if (dgvs_vars%present_patch(p)) then
            num_natvegp = num_natvegp + 1
            filter_natvegp(num_natvegp) = p
         end if
      end do

      ! Returns fpcgrid and nind

      call Light(bounds, num_natvegp, filter_natvegp, &
           carbonstate_vars, dgvs_vars)

      ! Returns updated fpcgrid, nind, crownarea, and present. Due to updated
      ! present, we do not use the natveg filter in this subroutine.

      call Establishment(bounds, &
           atm2lnd_vars, carbonflux_vars, carbonstate_vars, dgvs_vars)

      ! Reset dgvm variables needed in next yr (too few to keep subr. dvreset)

      do p = bounds%begp,bounds%endp
         carbonstate_vars%leafcmax_patch(p) = 0._r8
         atm2lnd_vars%t_mo_min_patch(p) = 1.0e+36_r8
      end do

    end associate 

  end subroutine CNDVDriver

  !-----------------------------------------------------------------------
  subroutine CNDVHist(bounds, dgvs_vars) 
    !
    ! !DESCRIPTION:
    ! Write CNDV history file
    !
    ! !USES:
    use shr_const_mod   , only : SHR_CONST_CDAY
    use shr_sys_mod     , only : shr_sys_getenv
    use clm_varpar      , only : maxpatch_pft
    use clm_varctl      , only : caseid, ctitle, finidat, fsurdat, paramfile, iulog
    use clm_varcon      , only : spval
    use clm_time_manager, only : get_ref_date, get_nstep, get_curr_date, get_curr_time
    use domainMod       , only : ldomain
    use fileutils       , only : get_filename
    use spmdMod         , only : masterproc
    use ncdio_pio
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  
    type(dgvs_type)  , intent(in) :: dgvs_vars
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: dgvm_fn      ! dgvm history filename
    type(file_desc_t)  :: ncid         ! netcdf file id
    integer :: ncprec                  ! output precision
    integer :: g,p,l                   ! indices
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
    real(r8) , pointer :: rbuf2dg (:,:)   ! Input:  [real(r8) (:,:)]  temporary 
    !-----------------------------------------------------------------------

    associate(& 
         fpcgrid => dgvs_vars%fpcgrid_patch , & ! Input:  [real(r8) (:)]  foliar projective cover on gridcell (fraction)    
         nind    => dgvs_vars%nind_patch      & ! Input:  [real(r8) (:)]  number of individuals (#/m**2)                    
         )

      allocate(rbuf2dg(bounds%begg:bounds%endg,maxpatch_pft), stat=ier)
      if (ier /= 0) call endrun(msg='histCNDV: allocation error for rbuf2dg'//&
           errMsg(__FILE__, __LINE__))

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
      if (ier /= 0) call endrun(msg='error: LOGNAME environment variable not defined'//&
           errMsg(__FILE__, __LINE__))

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

      str = get_filename(paramfile)
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
80    format(i4.4,'-',i2.2,'-',i2.2)
      write(basesec ,90) hours, minutes, secs
90    format(i2.2,':',i2.2,':',i2.2)
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

      rbuf2dg(bounds%begg : bounds%endg, :) = 0._r8
      do p = bounds%begp,bounds%endp
         g = pft%gridcell(p)
         l = pft%landunit(p)
         if (.not. lun%ifspecial(l)) rbuf2dg(g,pft%mxy(p)) = fpcgrid(p)*100._r8
      end do
      call ncd_io(ncid=ncid, varname='FPCGRID', dim1name=grlnd, data=rbuf2dg, &
           nt=1, flag='write')

      rbuf2dg(bounds%begg : bounds%endg, :) = 0._r8
      do p = bounds%begp,bounds%endp
         g = pft%gridcell(p)
         l = pft%landunit(p)
         if (.not. lun%ifspecial(l)) rbuf2dg(g,pft%mxy(p)) = nind(p)
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
      
    end associate

  end subroutine CNDVHist

  !-----------------------------------------------------------------------
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
  subroutine BuildNatVegFilter(bounds, num_natvegp, filter_natvegp, dgvs_vars)
    !
    ! !DESCRIPTION:
    ! Reconstruct a filter of naturally-vegetated Patches for use in DGVM
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)  :: bounds   
    integer           , intent(out) :: num_natvegp       ! number of patches in naturally-vegetated filter
    integer           , intent(out) :: filter_natvegp(:) ! patch filter for naturally-vegetated points
    type(dgvs_type)   , intent(in)  :: dgvs_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p
    !-----------------------------------------------------------------------

    num_natvegp = 0
    do p = bounds%begp,bounds%endp
       if (dgvs_vars%present_patch(p)) then
          num_natvegp = num_natvegp + 1
          filter_natvegp(num_natvegp) = p
       end if
    end do

  end subroutine BuildNatVegFilter

end module CNDVDriverMod
