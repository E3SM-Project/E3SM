module RtmRestFile

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or writes to/ the RTM restart file.
!
! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_sys_mod   , only : shr_sys_abort
  use RtmSpmd       , only : masterproc 
  use RtmVar        , only : rtmlon, rtmlat, iulog, inst_suffix, rpntfil, &
                             caseid, nsrest, brnch_retain_casename, &
                             finidat_rtm, nrevsn_rtm, wrmflag, &
                             nsrContinue, nsrBranch, nsrStartup, &
                             ctitle, version, username, hostname, conventions, source
  use RtmHistFile   , only : RtmHistRestart
  use RtmFileUtils  , only : relavu, getavu, opnfil, getfil
  use RtmTimeManager, only : timemgr_restart, get_nstep, get_curr_date, is_last_step
  use RunoffMod     , only : rtmCTL
  use RtmIO       
  use RtmDateTime
#ifdef INCLUDE_WRM
  use WRM_type_mod  , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_subw_io_mod, only : WRM_computeRelease
#endif
  use rof_cpl_indices , only : nt_rtm, rtm_tracers 
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: RtmRestFileName  
  public :: RtmRestFileRead
  public :: RtmRestFileWrite
  public :: RtmRestGetfile
  public :: RtmRestTimeManager
  public :: RtmRestart
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_dimset
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
! !PRIVATE TYPES: None
  private

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!=======================================================================

  subroutine RtmRestFileWrite( file, rdate )

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Read/write RTM restart file.

    ! !ARGUMENTS:
    implicit none
    character(len=*) , intent(in) :: file            ! output netcdf restart file
    character(len=*) , intent(in) :: rdate           ! restart file time stamp for name

    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    ! Define dimensions and variables

    if (masterproc) then	
       write(iulog,*)
       write(iulog,*)'restFile_open: writing RTM restart dataset '
       write(iulog,*)
    end if
    call ncd_pio_createfile(ncid, trim(file))
    call restFile_dimset( ncid )
    call RtmRestart( ncid, flag='define' )
    call RtmHistRestart ( ncid, flag='define', rdate=rdate )
    call timemgr_restart( ncid, flag='define' )
    call ncd_enddef(ncid)

    ! Write restart file variables
    call RtmRestart( ncid, flag='write' )
    call RtmHistRestart ( ncid, flag='write' )
    call timemgr_restart( ncid, flag='write' )
    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
    end if

    ! Write restart pointer file
    call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine RtmRestFileWrite

!-----------------------------------------------------------------------

  subroutine RtmRestFileRead( file )

    ! !DESCRIPTION:
    ! Read a RTM restart file.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
    !-------------------------------------

    ! Read file
    if (masterproc) write(iulog,*) 'Reading restart dataset'
    call ncd_pio_openfile (ncid, trim(file), 0)
    call RtmRestart( ncid, flag='read' )
    call RtmHistRestart(ncid, flag='read')
    call ncd_pio_closefile(ncid)

    ! Write out diagnostic info
    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine RtmRestFileRead

!-----------------------------------------------------------------------

  subroutine RtmRestTimeManager( file )

    ! !DESCRIPTION:
    ! Read a RTM restart file.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
    !-------------------------------------

    ! Read file
    if (masterproc) write(iulog,*) 'Reading restart Timemanger'
    call ncd_pio_openfile (ncid, trim(file), 0)
    call timemgr_restart(ncid, flag='read')
    call ncd_pio_closefile(ncid)

    ! Write out diagnostic info
    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine RtmRestTimeManager

!-----------------------------------------------------------------------

  subroutine RtmRestGetfile( file, path )

    !---------------------------------------------------
    ! DESCRIPTION:
    ! Determine and obtain netcdf restart file

    ! ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file

    ! LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !---------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 
    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if
       
    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn_rtm"
    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn_rtm)
       if (nrevsn_rtm(length-2:length) == '.nc') then
          path = trim(nrevsn_rtm) 
       else
          path = trim(nrevsn_rtm) // '.nc'
       end if
       call getfil( path, file, 0 )
       
       ! Check case name consistency (case name must be different 
       ! for branch run, unless brnch_retain_casename is set)
       ctest = 'xx.'//trim(caseid)//'.mosart'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          write(iulog,*) 'Must change case name on branch run if ',&
               'brnch_retain_casename namelist is not set'
          write(iulog,*) 'previous case filename= ',trim(file),&
               ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
               ' ftest = ',trim(ftest)
          call shr_sys_abort()
       end if
    end if

    ! Initial run 
    if (nsrest==nsrStartup) then
       call getfil( finidat_rtm, file, 0 )
    end if
    
  end subroutine RtmRestGetfile

!-----------------------------------------------------------------------

  subroutine restFile_read_pfile( pnamer )

    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks

    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of restart file

    ! !LOCAL VARIABLES:
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !--------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn_rtm] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.
       
    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = './'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

!-----------------------------------------------------------------------

  subroutine restFile_write_pfile( fnamer )

    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name

    if (masterproc) then
       nio = getavu()
       filename= './'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )
       
       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile


!-----------------------------------------------------------------------

  character(len=256) function RtmRestFileName( rdate )

    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name 

    RtmRestFileName = "./"//trim(caseid)//".mosart"//trim(inst_suffix)//".r."//trim(rdate)//".nc"
    if (masterproc) then
       write(iulog,*)'writing restart file ',trim(RtmRestFileName),' for model date = ',rdate
    end if
 
  end function RtmRestFileName

!------------------------------------------------------------------------

  subroutine restFile_dimset( ncid )

    !----------------------------------------------------------------
    ! !DESCRIPTION:
    ! Read/Write initial data from/to netCDF instantaneous initial data file

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid

    ! !LOCAL VARIABLES:
    integer :: dimid               ! netCDF dimension id
    integer :: ier                 ! error status
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len=*),parameter :: subname='restFile_dimset' ! subroutine name
    !----------------------------------------------------------------

    ! Define dimensions
    
    call ncd_defdim(ncid, 'rtmlon'  , rtmlon         , dimid)
    call ncd_defdim(ncid, 'rtmlat'  , rtmlat         , dimid)
    call ncd_defdim(ncid, 'string_length', 64        , dimid)
       
    ! Define global attributes
    
    call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
    call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
    call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
    call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
    str = '$Id: restFileMod.F90 35676 2012-03-22 21:48:04Z tcraig $'
    call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
    call ncd_putatt(ncid, NCD_GLOBAL, 'title', &
          'RTM Restart information, required to continue a simulation' )

  end subroutine restFile_dimset
  
!-----------------------------------------------------------------------

  subroutine RtmRestart(ncid, flag)

    !-----------------------------------------------------------------------
    ! DESCRIPTION:
    ! Read/write RTM restart data.
    !
    ! ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout)  :: ncid ! netcdf id
    character(len=*) , intent(in) :: flag   ! 'read' or 'write'
    ! LOCAL VARIABLES:
    logical :: readvar          ! determine if variable is on initial file
    integer :: nt,nv,n          ! indices
    integer :: nvmax            ! number of variables to write
    integer :: ig, idam         ! indices
    logical :: varok            ! variable ok flag for reading/writing
    real(r8) , pointer :: dfld(:) ! temporary array
    logical :: compute_release  ! if release or stormthstop not read
    logical :: storage_read     ! check if storage on restart file
    logical :: release_read     ! check if release on restart file
    logical :: stormth_read     ! check if stormthstop on restart file
    character(len=32)  :: vname,uname
    character(len=128) :: lname
    !-----------------------------------------------------------------------

    compute_release = .false.
    storage_read = .true.
    release_read = .true.
    stormth_read = .true.
#ifdef INCLUDE_WRM
    nvmax = 10
#else
    nvmax = 7
#endif

    do nv = 1,nvmax
    do nt = 1,nt_rtm

       varok = .true.

       if (nv == 1) then
          vname = 'RTM_VOLR_'//trim(rtm_tracers(nt))
          lname = 'water volume in cell (volr)'
          uname = 'm3'
          dfld  => rtmCTL%volr(:,nt)
       elseif (nv == 2) then
          vname = 'RTM_RUNOFF_'//trim(rtm_tracers(nt))
          lname = 'runoff (runoff)'
          uname = 'm3/s'
          dfld  => rtmCTL%runoff(:,nt)
       elseif (nv == 3) then
          vname = 'RTM_DVOLRDT_'//trim(rtm_tracers(nt))
          lname = 'water volume change in cell (dvolrdt)'
          uname = 'mm/s'
          dfld  => rtmCTL%dvolrdt(:,nt)
       elseif (nv == 4) then
          vname = 'RTM_WH_'//trim(rtm_tracers(nt))
          lname = 'surface water storage at hillslopes in cell'
          uname = 'm'
          dfld  => rtmCTL%wh(:,nt)
       elseif (nv == 5) then
          vname = 'RTM_WT_'//trim(rtm_tracers(nt))
          lname = 'water storage in tributary channels in cell'
          uname = 'm3'
          dfld  => rtmCTL%wt(:,nt)
       elseif (nv == 6) then
          vname = 'RTM_WR_'//trim(rtm_tracers(nt))
          lname = 'water storage in main channel in cell'
          uname = 'm3'
          dfld  => rtmCTL%wr(:,nt)
       elseif (nv == 7) then
          vname = 'RTM_EROUT_'//trim(rtm_tracers(nt))
          lname = 'instataneous flow out of main channel in cell'
          uname = 'm3/s'
          dfld  => rtmCTL%erout(:,nt)
#ifdef INCLUDE_WRM
       elseif (nv == 8 .and. trim(rtm_tracers(nt)) == 'LIQ') then
          varok = .false.
          if (wrmflag) then
             varok = .true.
             StorWater%storageG = 0._r8
             if (flag == 'write') then
                do idam = 1, ctlSubwWRM%localNumDam
                   ig = WRMUnit%icell(idam)
                   StorWater%storageG(ig) = StorWater%storage(idam)
                enddo
             endif
             vname = 'DAM_STORAGE_'//trim(rtm_tracers(nt))
             lname = 'dam storage'
             uname = 'm3'
             dfld  => StorWater%storageG(:)
          endif
       elseif (nv == 9 .and. trim(rtm_tracers(nt)) == 'LIQ') then
          varok = .false.
          if (wrmflag) then
             varok = .true.
             StorWater%releaseG = 0._r8
             if (flag == 'write') then
                do idam = 1, ctlSubwWRM%localNumDam
                   ig = WRMUnit%icell(idam)
                   StorWater%releaseG(ig) = StorWater%release(idam)
                enddo
             endif
             vname = 'DAM_RELEASE_'//trim(rtm_tracers(nt))
             lname = 'dam release'
             uname = 'm3'
             dfld  => StorWater%releaseG(:)
          endif
       elseif (nv == 10 .and. trim(rtm_tracers(nt)) == 'LIQ') then
          varok = .false.
          if (wrmflag) then
             varok = .true.
             WRMUnit%StorMthStOpG = 0._r8
             if (flag == 'write') then
                do idam = 1, ctlSubwWRM%localNumDam
                   ig = WRMUnit%icell(idam)
                   WRMUnit%StorMthStOpG(ig) = WRMUnit%StorMthStOp(idam)
                enddo
             endif
             vname = 'DAM_STORMTHSTOP_'//trim(rtm_tracers(nt))
             lname = 'dam StorMthStOp'
             uname = 'm3'
             dfld  => WRMUnit%StorMthStOpG(:)
          endif
#endif
       else
          varok = .false.
       endif

       if (varok) then
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=ncd_double,  dim1name='rtmlon', dim2name='rtmlat', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(vname), data=dfld, dim1name='allrof', &
               ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (nsrest == nsrContinue) then
                call shr_sys_abort()
             else
                dfld = 0._r8
                if (vname == 'DAM_STORAGE_LIQ') then
                   storage_read = .false.
                elseif (vname == 'DAM_RELEASE_LIQ') then
                   release_read = .false.
                   compute_release = .true.
                elseif (vname == 'DAM_STORMTHSTOP_LIQ') then
                   stormth_read = .false.
                   compute_release = .true.
                endif
             end if
          end if
       end if
       end if  ! varok

    enddo
    enddo

    if (flag == 'read') then
       do n = rtmCTL%begr,rtmCTL%endr
          do nt = 1,nt_rtm
             if (abs(rtmCTL%volr(n,nt))    > 1.e30) rtmCTL%volr(n,nt) = 0.
             if (abs(rtmCTL%runoff(n,nt))  > 1.e30) rtmCTL%runoff(n,nt) = 0.
             if (abs(rtmCTL%dvolrdt(n,nt)) > 1.e30) rtmCTL%dvolrdt(n,nt) = 0.
             if (abs(rtmCTL%wh(n,nt))      > 1.e30) rtmCTL%wh(n,nt) = 0.
             if (abs(rtmCTL%wt(n,nt))      > 1.e30) rtmCTL%wt(n,nt) = 0.
             if (abs(rtmCTL%wr(n,nt))      > 1.e30) rtmCTL%wr(n,nt) = 0.
             if (abs(rtmCTL%erout(n,nt))   > 1.e30) rtmCTL%erout(n,nt) = 0.
          end do  ! nt
#ifdef INCLUDE_WRM
          if (wrmflag) then
             if (abs(storWater%storageG(n)) > 1.e30) storWater%storageG(n) = 0.
             if (abs(storWater%releaseG(n)) > 1.e30) storWater%releaseG(n) = 0.
             if (abs(WRMUnit%StorMthStOpG(n)) > 1.e30) WRMUnit%StorMthStOpG(n) = 0.
          endif
#endif
          if (rtmCTL%mask(n) == 1) then
             do nt = 1,nt_rtm
                rtmCTL%runofflnd(n,nt) = rtmCTL%runoff(n,nt)
                rtmCTL%dvolrdtlnd(n,nt)= rtmCTL%dvolrdt(n,nt)
             end do
          elseif (rtmCTL%mask(n) >= 2) then
             do nt = 1,nt_rtm
                rtmCTL%runoffocn(n,nt) = rtmCTL%runoff(n,nt)
                rtmCTL%dvolrdtocn(n,nt)= rtmCTL%dvolrdt(n,nt)
             enddo
          endif
       enddo  ! n

#ifdef INCLUDE_WRM
       ! only overwrite fields that have been read, otherwise, use initial values
       if (wrmflag) then
          do idam = 1, ctlSubwWRM%localNumDam
             ig = WRMUnit%icell(idam)
             if (storage_read) StorWater%storage(idam) = StorWater%storageG(ig)
             if (release_read) StorWater%release(idam) = StorWater%releaseG(ig)
             if (stormth_read) WRMUnit%StorMthStOp(idam) = WRMUnit%StorMthStOpG(ig)
          enddo
          if (compute_release) then
             call WRM_computeRelease()
          endif
       endif
#endif
    endif  ! read

  end subroutine RtmRestart

end module RtmRestFile



