module RtmMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! Mosart Routing Model
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_PI, SHR_CONST_CDAY
  use rof_cpl_indices , only : nt_rtm, rtm_tracers 
  use RtmSpmd         , only : masterproc, npes, iam, mpicom_rof, ROFID, mastertask, &
                               MPI_REAL8,MPI_INTEGER,MPI_CHARACTER,MPI_LOGICAL,MPI_MAX
  use RtmVar          , only : re, spval, rtmlon, rtmlat, iulog, ice_runoff, &
                               frivinp_rtm, finidat_rtm, nrevsn_rtm, &
                               nsrContinue, nsrBranch, nsrStartup, nsrest, &
                               inst_index, inst_suffix, inst_name, wrmflag, &
                               smat_option, decomp_option, barrier_timers
  use RtmFileUtils    , only : getfil, getavu, relavu
  use RtmTimeManager  , only : timemgr_init, get_nstep, get_curr_date, advance_timestep
  use RtmHistFlds     , only : RtmHistFldsInit, RtmHistFldsSet 
  use RtmHistFile     , only : RtmHistUpdateHbuf, RtmHistHtapesWrapup, RtmHistHtapesBuild, &
                               rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq,     &
                               rtmhist_avgflag_pertape, rtmhist_avgflag_pertape, & 
                               rtmhist_fincl1, rtmhist_fincl2, rtmhist_fincl3,   &
                               rtmhist_fexcl1, rtmhist_fexcl2, rtmhist_fexcl3,   &
                               max_tapes, max_namlen
  use RtmRestFile     , only : RtmRestTimeManager, RtmRestGetFile, RtmRestFileRead, &
                               RtmRestFileWrite, RtmRestFileName
  use RunoffMod       , only : RunoffInit, rtmCTL, Tctl, Tunit, TRunoff, Tpara, &
                               gsmap_r, &
                               SMatP_dnstrm, avsrc_dnstrm, avdst_dnstrm, &
                               SMatP_direct, avsrc_direct, avdst_direct, &
                               SMatP_eroutUp, avsrc_eroutUp, avdst_eroutUp
  use MOSART_physics_mod, only : Euler
  use MOSART_physics_mod, only : updatestate_hillslope, updatestate_subnetwork, &
                                 updatestate_mainchannel
#ifdef INCLUDE_WRM
  use WRM_type_mod    , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_subw_IO_mod , only : WRM_init, WRM_computeRelease
#endif
  use RtmIO
  use mct_mod
  use perf_mod
  use pio
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini          ! Initialize MOSART grid
  public Rtmrun          ! River routing model
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: RtmFloodInit

! !PRIVATE TYPES:

! MOSART tracers
  character(len=256) :: rtm_trstr   ! tracer string

! MOSART namelists
  integer, save :: coupling_period   ! mosart coupling period
  integer, save :: delt_mosart       ! mosart internal timestep (->nsub)

! MOSART constants
  real(r8) :: cfl_scale = 1.0_r8    ! cfl scale factor, must be <= 1.0
  real(r8) :: river_depth_minimum = 1.e-4 ! gridcell average minimum river depth [m]

!global (glo)
  integer , pointer :: ID0_global(:)  ! local ID index
  integer , pointer :: dnID_global(:) ! downstream ID based on ID0
  real(r8), pointer :: area_global(:) ! area
  integer , pointer :: IDkey(:)       ! translation key from ID to gindex

!local (gdc)
  real(r8), save, pointer :: evel(:,:)       ! effective tracer velocity (m/s)
  real(r8), save, pointer :: flow(:,:)       ! mosart flow (m3/s)
  real(r8), save, pointer :: eroup_lagi(:,:) ! erout previous timestep (m3/s)
  real(r8), save, pointer :: eroup_lagf(:,:) ! erout current timestep (m3/s)
  real(r8), save, pointer :: erowm_regi(:,:) ! erout previous timestep (m3/s)
  real(r8), save, pointer :: erowm_regf(:,:) ! erout current timestep (m3/s)
  real(r8), save, pointer :: eroutup_avg(:,:)! eroutup average over coupling period (m3/s)
  real(r8), save, pointer :: erlat_avg(:,:)  ! erlateral average over coupling period (m3/s)

! global MOSART grid
  real(r8),pointer :: rlatc(:)    ! latitude of 1d grid cell (deg)
  real(r8),pointer :: rlonc(:)    ! longitude of 1d grid cell (deg)
  real(r8),pointer :: rlats(:)    ! latitude of 1d south grid cell edge (deg)
  real(r8),pointer :: rlatn(:)    ! latitude of 1d north grid cell edge (deg)
  real(r8),pointer :: rlonw(:)    ! longitude of 1d west grid cell edge (deg)
  real(r8),pointer :: rlone(:)    ! longitude of 1d east grid cell edge (deg)

  logical :: do_rtmflood
  logical :: do_rtm

!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmini
!
! !INTERFACE:
  subroutine Rtmini(rtm_active,flood_active)
!
! !DESCRIPTION:
! Initialize MOSART grid, mask, decomp
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    logical, intent(out) :: rtm_active
    logical, intent(out) :: flood_active
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
! Update: T Craig, Dec 2006
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: effvel0 = 10.0_r8             ! default velocity (m/s)
    real(r8) :: effvel(nt_rtm)                ! downstream velocity (m/s)
    real(r8) :: edgen                         ! North edge of the direction file
    real(r8) :: edgee                         ! East edge of the direction file
    real(r8) :: edges                         ! South edge of the direction file
    real(r8) :: edgew                         ! West edge of the direction file
    integer  :: i,j,k,n,ng,g,n2,nt,nn         ! loop indices
    integer  :: i1,j1,i2,j2
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8) :: lrtmarea                      ! tmp local sum of area
    real(r8),allocatable :: tempr(:,:)        ! temporary buffer
    integer ,allocatable :: itempr(:,:)       ! temporary buffer
    integer ,allocatable :: idxocn(:)         ! downstream ocean outlet cell
    integer ,allocatable :: nupstrm(:)        ! number of upstream cells including own cell
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of gridcells on a pe
    integer ,allocatable :: nba(:)            ! number of basins on each pe
    integer ,allocatable :: nrs(:)            ! begr on each pe
    integer ,allocatable :: basin(:)          ! basin to mosart mapping
    integer  :: nmos,nmos_chk                 ! number of mosart points
    integer  :: nout,nout_chk                 ! number of basin with outlets
    integer  :: nbas,nbas_chk                 ! number of basin/ocean points
    integer  :: nrof,nrof_chk                 ! num of active mosart points
    integer  :: baspe                         ! pe with min number of mosart cells
    integer  :: maxrtm                        ! max num of rtms per pe for decomp
    integer  :: minbas,maxbas                 ! used for decomp search
    integer  :: nl,nloops                     ! used for decomp search
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: numr                          ! tot num of roff pts on all pes
    real(r8) :: dtover,dtovermax              ! ts calc temporaries
    type(file_desc_t) :: ncid                 ! netcdf file id
    integer  :: dimid                         ! netcdf dimension identifier
    integer  :: nroflnd                       ! local number of land runoff 
    integer  :: nrofocn                       ! local number of ocn runoff
    integer  :: pid,np,npmin,npmax,npint      ! log loop control
    integer  :: na,nb,ns                      ! mct sizes
    integer  :: ni,no,go                      ! tmps
    integer ,pointer  :: rgdc2glo(:)          ! temporary for initialization
    integer ,pointer  :: rglo2gdc(:)          ! temporary for initialization
    integer ,pointer  :: gmask(:)             ! global mask
    logical           :: found                ! flag
    character(len=256):: fnamer               ! name of netcdf restart file 
    character(len=256):: pnamer               ! full pathname of netcdf restart file
    character(len=256):: locfn                ! local file name
    character(len=16384) :: rList             ! list of fields for SM multiply
    integer           :: unitn                ! unit for namelist file
    integer,parameter :: dbug = 3             ! 0 = none, 1=normal, 2=much, 3=max
    logical :: lexist                         ! File exists
    character(len= 7) :: runtyp(4)            ! run type
    integer ,allocatable :: gindex(:)         ! global index
    integer           :: cnt, lsize, gsize    ! counter
    integer           :: igrow,igcol,iwgt     ! mct field indices
    character(len=256):: nlfilename_rof       ! namelist filename
    type(mct_avect)   :: avtmp, avtmpG        ! temporary avects
    type(mct_sMat)    :: sMat                 ! temporary sparse matrix, needed for sMatP
    character(len=*),parameter :: subname = '(Rtmini) '
!-----------------------------------------------------------------------

    !-------------------------------------------------------
    ! Read in mosart namelist
    !-------------------------------------------------------

    namelist /mosart_inparm / ice_runoff, do_rtm, do_rtmflood, &
         frivinp_rtm, finidat_rtm, nrevsn_rtm, coupling_period, &
         rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq, &
         rtmhist_fincl1,  rtmhist_fincl2, rtmhist_fincl3, &
         rtmhist_fexcl1,  rtmhist_fexcl2, rtmhist_fexcl3, &
         rtmhist_avgflag_pertape, decomp_option, wrmflag, &
         smat_option, delt_mosart, barrier_timers

    ! Preset values
    do_rtm      = .true.
    do_rtmflood = .false.
    ice_runoff  = .true.
    wrmflag     = .false.
    barrier_timers = .false.
    finidat_rtm = ' '
    nrevsn_rtm  = ' '
    coupling_period   = -1
    delt_mosart = 3600
    decomp_option = 'basin'
    smat_option = 'opt'

    nlfilename_rof = "mosart_in" // trim(inst_suffix)
    inquire (file = trim(nlfilename_rof), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_rof does NOT exist:'&
            //trim(nlfilename_rof)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_rof does not exist')
    end if
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in mosart_inparm namelist from: ', trim(nlfilename_rof)
       open( unitn, file=trim(nlfilename_rof), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, mosart_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on mosart_inparm read' )
          endif
       end do
       call relavu( unitn )
    end if

    call mpi_bcast (coupling_period,   1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (delt_mosart    ,   1, MPI_INTEGER, 0, mpicom_rof, ier)

    call mpi_bcast (finidat_rtm  , len(finidat_rtm)  , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (frivinp_rtm  , len(frivinp_rtm)  , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (nrevsn_rtm   , len(nrevsn_rtm)   , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (decomp_option, len(decomp_option), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (smat_option  , len(smat_option)  , MPI_CHARACTER, 0, mpicom_rof, ier)

    call mpi_bcast (do_rtm,         1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (do_rtmflood,    1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (ice_runoff,     1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (wrmflag,        1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (barrier_timers, 1, MPI_LOGICAL, 0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_nhtfrq, size(rtmhist_nhtfrq), MPI_INTEGER,   0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_mfilt , size(rtmhist_mfilt) , MPI_INTEGER,   0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_ndens , size(rtmhist_ndens) , MPI_INTEGER,   0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_fexcl1, (max_namlen+2)*size(rtmhist_fexcl1), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fexcl2, (max_namlen+2)*size(rtmhist_fexcl2), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fexcl3, (max_namlen+2)*size(rtmhist_fexcl3), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl1, (max_namlen+2)*size(rtmhist_fincl1), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl2, (max_namlen+2)*size(rtmhist_fincl2), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl3, (max_namlen+2)*size(rtmhist_fincl3), MPI_CHARACTER, 0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_avgflag_pertape, size(rtmhist_avgflag_pertape), MPI_CHARACTER, 0, mpicom_rof, ier)

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    if (masterproc) then
       write(iulog,*) 'define run:'
       write(iulog,*) '   run type              = ',runtyp(nsrest+1)
      !write(iulog,*) '   case title            = ',trim(ctitle)
      !write(iulog,*) '   username              = ',trim(username)
      !write(iulog,*) '   hostname              = ',trim(hostname)
       write(iulog,*) '   coupling_period       = ',coupling_period
       write(iulog,*) '   delt_mosart           = ',delt_mosart
       write(iulog,*) '   decomp_option         = ',trim(decomp_option)
       write(iulog,*) '   smat_option           = ',trim(smat_option)
       write(iulog,*) '   wrmflag               = ',wrmflag
       write(iulog,*) '   barrier_timers        = ',barrier_timers
       if (nsrest == nsrStartup .and. finidat_rtm /= ' ') then
          write(iulog,*) '   MOSART initial data   = ',trim(finidat_rtm)
       end if
    endif

    rtm_active = do_rtm
    flood_active = do_rtmflood
    
    if (do_rtm) then
       if (frivinp_rtm == ' ') then
          call shr_sys_abort( subname//' ERROR: do_rtm TRUE, but frivinp_rtm NOT set' )
       else
          if (masterproc) then
             write(iulog,*) '   MOSART river data       = ',trim(frivinp_rtm)
          endif
       end if
    else
       if (masterproc) then
          write(iulog,*)'MOSART will not be active '
       endif
       RETURN
    end if

    if (coupling_period <= 0) then
       write(iulog,*) subname,' ERROR MOSART coupling_period invalid',coupling_period
       call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
    endif

    if (delt_mosart <= 0) then
       write(iulog,*) subname,' ERROR MOSART delt_mosart invalid',delt_mosart
       call shr_sys_abort( subname//' ERROR: delt_mosart invalid' )
    endif
       
    do i = 1, max_tapes
       if (rtmhist_nhtfrq(i) == 0) then
          rtmhist_mfilt(i) = 1
       else if (rtmhist_nhtfrq(i) < 0) then
          rtmhist_nhtfrq(i) = nint(-rtmhist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*coupling_period))
       endif
    end do

    !-------------------------------------------------------
    ! Initialize MOSART time manager 
    !-------------------------------------------------------

    ! Intiialize MOSART pio
    call ncd_pio_init()

    ! Obtain restart file if appropriate
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
        (nsrest == nsrBranch  )) then
       call RtmRestGetfile( file=fnamer, path=pnamer )
    endif       

    ! Initialize time manager
    if (nsrest == nsrStartup) then  
       call timemgr_init(dtime_in=coupling_period)
    else
       call RtmRestTimeManager(file=fnamer)
    end if

    !-------------------------------------------------------
    ! Initialize rtm_trstr
    !-------------------------------------------------------

    rtm_trstr = trim(rtm_tracers(1))
    do n = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(n))
    enddo
    if (masterproc) then
       write(iulog,*)'MOSART tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !-------------------------------------------------------
    ! Read input data (river direction file)
    !-------------------------------------------------------

    ! Useful constants and initial values
    deg2rad = SHR_CONST_PI / 180._r8

    call t_startf('mosarti_grid')

    call getfil(frivinp_rtm, locfn, 0 )
    if (masterproc) then
       write(iulog,*) 'Read in MOSART file name: ',trim(frivinp_rtm)
       call shr_sys_flush(iulog)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'lon',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlon)
    call ncd_inqdid(ncid,'lat',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlat)

    if (masterproc) then
       write(iulog,*) 'Values for rtmlon/rtmlat: ',rtmlon,rtmlat
       write(iulog,*) 'Successfully read MOSART dimensions'
       call shr_sys_flush(iulog)
    endif

    ! Allocate variables
    allocate(rlonc(rtmlon), rlatc(rtmlat), &
             rlonw(rtmlon), rlone(rtmlon), &
             rlats(rtmlat), rlatn(rtmlat), &
             rtmCTL%rlon(rtmlon),          &
             rtmCTL%rlat(rtmlat),          &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' : Allocation ERROR for rlon'
       call shr_sys_abort(subname//' ERROR alloc for rlon')
    end if

    ! reading the routing parameters
    allocate ( &
              ID0_global(rtmlon*rtmlat), area_global(rtmlon*rtmlat), &
              dnID_global(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname, ' : Allocation error for ID0_global'
       call shr_sys_abort(subname//' ERROR alloc for ID0')
    end if

    allocate(tempr(rtmlon,rtmlat))  
    allocate(itempr(rtmlon,rtmlat))  

    call ncd_io(ncid=ncid, varname='longxy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART longitudes')
    if (masterproc) write(iulog,*) 'Read longxy ',minval(tempr),maxval(tempr)
    do i=1,rtmlon
       rtmCTL%rlon(i) = tempr(i,1)
       rlonc(i) = tempr(i,1)
    enddo
    if (masterproc) write(iulog,*) 'rlonc ',minval(rlonc),maxval(rlonc)

    call ncd_io(ncid=ncid, varname='latixy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART latitudes')
    if (masterproc) write(iulog,*) 'Read latixy ',minval(tempr),maxval(tempr)
    do j=1,rtmlat
       rtmCTL%rlat(j) = tempr(1,j)
       rlatc(j) = tempr(1,j)
    end do
    if (masterproc) write(iulog,*) 'rlatc ',minval(rlatc),maxval(rlatc)

    call ncd_io(ncid=ncid, varname='area', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART area')
    if (masterproc) write(iulog,*) 'Read area ',minval(tempr),maxval(tempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       area_global(n) = tempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'area ',minval(tempr),maxval(tempr)

    call ncd_io(ncid=ncid, varname='ID', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART ID')
    if (masterproc) write(iulog,*) 'Read ID ',minval(itempr),maxval(itempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       ID0_global(n) = itempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'ID ',minval(itempr),maxval(itempr)

    call ncd_io(ncid=ncid, varname='dnID', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART dnID')
    if (masterproc) write(iulog,*) 'Read dnID ',minval(itempr),maxval(itempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       dnID_global(n) = itempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'dnID ',minval(itempr),maxval(itempr)

    deallocate(tempr)
    deallocate(itempr)             

    call ncd_pio_closefile(ncid)

    !-------------------------------------------------------
    ! RESET dnID indices based on ID0
    ! rename the dnID values to be consistent with global grid indexing.
    ! where 1 = lower left of grid and rtmlon*rtmlat is upper right.
    ! ID0 is the "key", modify dnID based on that.  keep the IDkey around
    ! for as long as needed.  This is a key that translates the ID0 value
    ! to the gindex value.  compute the key, then apply the key to dnID_global.
    ! As part of this, check that each value of ID0 is unique and within
    ! the range of 1 to rtmlon*rtmlat.
    !-------------------------------------------------------

    allocate(IDkey(rtmlon*rtmlat))
    IDkey = 0
    do n=1,rtmlon*rtmlat
       if (ID0_global(n) < 0 .or. ID0_global(n) > rtmlon*rtmlat) then
          write(iulog,*) subname,' ERROR ID0 out of range',n,ID0_global(n)
          call shr_sys_abort(subname//' ERROR error ID0 out of range')
       endif
       if (IDkey(ID0_global(n)) /= 0) then
          write(iulog,*) subname,' ERROR ID0 value occurs twice',n,ID0_global(n)
          call shr_sys_abort(subname//' ERROR ID0 value occurs twice')
       endif
       IDkey(ID0_global(n)) = n
    enddo
    if (minval(IDkey) < 1) then
       write(iulog,*) subname,' ERROR IDkey incomplete'
       call shr_sys_abort(subname//' ERROR IDkey incomplete')
    endif
    do n=1,rtmlon*rtmlat
       if (dnID_global(n) > 0 .and. dnID_global(n) <= rtmlon*rtmlat) then
          if (IDkey(dnID_global(n)) > 0 .and. IDkey(dnID_global(n)) <= rtmlon*rtmlat) then
             dnID_global(n) = IDkey(dnID_global(n))
          else
             write(iulog,*) subname,' ERROR bad IDkey',n,dnID_global(n),IDkey(dnID_global(n))
             call shr_sys_abort(subname//' ERROR bad IDkey')
          endif
       endif
    enddo
    deallocate(ID0_global)

    !-------------------------------------------------------
    ! Derive gridbox edges
    !-------------------------------------------------------

    ! assuming equispaced grid, calculate edges from rtmlat/rtmlon
    ! w/o assuming a global grid
    edgen = maxval(rlatc) + 0.5*abs(rlatc(1) - rlatc(2))
    edges = minval(rlatc) - 0.5*abs(rlatc(1) - rlatc(2))
    edgee = maxval(rlonc) + 0.5*abs(rlonc(1) - rlonc(2))
    edgew = minval(rlonc) - 0.5*abs(rlonc(1) - rlonc(2))

    if ( edgen .ne.  90._r8 )then
       if (masterproc) write(iulog,*) 'Regional grid: edgen = ', edgen
    end if
    if ( edges .ne. -90._r8 )then
       if (masterproc) write(iulog,*) 'Regional grid: edges = ', edges
    end if
    if ( edgee .ne. 180._r8 )then
       if (masterproc) write(iulog,*) 'Regional grid: edgee = ', edgee
    end if
    if ( edgew .ne.-180._r8 )then
       if (masterproc) write(iulog,*) 'Regional grid: edgew = ', edgew
    end if

    ! Set edge latitudes (assumes latitudes are constant for a given longitude)
    rlats(:) = edges
    rlatn(:) = edgen
    do j = 2, rtmlat
       if (rlatc(2) > rlatc(1)) then ! South to North grid
          rlats(j)   = (rlatc(j-1) + rlatc(j)) / 2._r8
          rlatn(j-1) = rlats(j)
       else  ! North to South grid
          rlatn(j)   = (rlatc(j-1) + rlatc(j)) / 2._r8
          rlats(j-1) = rlatn(j)
       end if
    end do

!    if (masterproc) then
!       write(iulog,*) 'tcx rlats = ',rlats
!       write(iulog,*) 'tcx rlatn = ',rlatn
!    endif

    ! Set edge longitudes
    rlonw(:) = edgew
    rlone(:) = edgee
    dx = (edgee - edgew) / rtmlon
    do i = 2, rtmlon
       rlonw(i)   = rlonw(i) + (i-1)*dx
       rlone(i-1) = rlonw(i)
    end do

!    if (masterproc) then
!       write(iulog,*) 'tcx rlonw = ',rlonw
!       write(iulog,*) 'tcx rlone = ',rlone
!    endif

    call t_stopf('mosarti_grid')

    !-------------------------------------------------------
    ! Determine mosart ocn/land mask (global, all procs)
    !-------------------------------------------------------

    call t_startf('mosarti_decomp')

    allocate (gmask(rtmlon*rtmlat), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname, ' : Allocation ERROR for gmask'
       call shr_sys_abort(subname//' ERROR alloc for gmask')
    end if

    !  1=land, 
    !  2=ocean,
    !  3=ocean outlet from land

    gmask = 2    ! assume ocean point
    do n=1,rtmlon*rtmlat         ! mark all downstream points as outlet
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= rtmlon*rtmlat)) then
          gmask(nr) = 3          ! <- nr
       end if
    enddo
    do n=1,rtmlon*rtmlat         ! now mark all points with downstream points as land
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= rtmlon*rtmlat)) then
          gmask(n) = 1           ! <- n
       end if
    enddo

    !-------------------------------------------------------
    ! Compute total number of basins and runoff points
    !-------------------------------------------------------

    nbas = 0
    nrof = 0
    nout = 0
    nmos = 0
    do nr=1,rtmlon*rtmlat
       if (gmask(nr) == 3) then
          nout = nout + 1
          nbas = nbas + 1
          nmos = nmos + 1
          nrof = nrof + 1
       elseif (gmask(nr) == 2) then
          nbas = nbas + 1
          nrof = nrof + 1
       elseif (gmask(nr) == 1) then
          nmos = nmos + 1
          nrof = nrof + 1
       endif
    enddo
    if (masterproc) then
       write(iulog,*) 'Number of outlet basins = ',nout
       write(iulog,*) 'Number of total  basins = ',nbas
       write(iulog,*) 'Number of mosart points = ',nmos
       write(iulog,*) 'Number of runoff points = ',nrof
    endif

    !-------------------------------------------------------
    ! Compute river basins, actually compute ocean outlet gridcell
    !-------------------------------------------------------

    ! idxocn = final downstream cell, index is global 1d ocean gridcell
    ! nupstrm = number of source gridcells upstream including self

    allocate(idxocn(rtmlon*rtmlat),nupstrm(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' : Allocation ERROR for ',&
            'idxocn,nupstrm'
       call shr_sys_abort(subname//' ERROR alloc for idxocn nupstrm')
    end if

    call t_startf('mosarti_dec_basins')
    idxocn  = 0
    nupstrm = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             nupstrm(n) = nupstrm(n) + 1
             n = dnID_global(n)
             g = g + 1
          end do
          if (gmask(n) == 3) then           ! found ocean outlet 
             nupstrm(n) = nupstrm(n) + 1    ! one more land cell for n
             idxocn(nr) = n                 ! set ocean outlet or nr to n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) subname,' ERROR closed basin found', &
               g,nr,gmask(nr),dnID_global(nr), &
               n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR closed basin found')
          elseif (gmask(n) == 2) then
             write(iulog,*) subname,' ERROR found invalid ocean cell ',nr
             call shr_sys_abort(subname//' ERROR found invalid ocean cell')
          else 
             write(iulog,*) subname,' ERROR downstream cell is unknown', &
               g,nr,gmask(nr),dnID_global(nr), &
               n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR downstream cell is unknown')
          endif
       elseif (gmask(n) >= 2) then  ! ocean, give to self
          nupstrm(n) = nupstrm(n) + 1
          idxocn(nr) = n
       endif
    enddo
    call t_stopf('mosarti_dec_basins')

    ! check

    nbas_chk = 0
    nrof_chk = 0
    do nr=1,rtmlon*rtmlat
!      !if (masterproc) write(iulog,*) 'nupstrm check ',nr,gmask(nr),nupstrm(nr),idxocn(nr)
       if (gmask(nr) >= 2 .and. nupstrm(nr) > 0) then
          nbas_chk = nbas_chk + 1
          nrof_chk = nrof_chk + nupstrm(nr)
       endif
    enddo

    if (nbas_chk /= nbas .or. nrof_chk /= nrof) then
       write(iulog,*) subname,' ERROR nbas nrof check',nbas,nbas_chk,nrof,nrof_chk
       call shr_sys_abort(subname//' ERROR nbas nrof check')
    endif

    !-------------------------------------------------------
    !--- Now allocate those basins to pes
    !-------------------------------------------------------

    call t_startf('mosarti_dec_distr')

    !--- this is the heart of the decomp, need to set pocn and nop by the end of this
    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of mosart cells/pe 

    allocate(pocn(rtmlon*rtmlat),     &  !global mosart array
             nop(0:npes-1), &
             nba(0:npes-1))

    pocn = -99
    nop = 0
    nba = 0

    if (trim(decomp_option) == 'basin') then
       baspe = 0
       maxrtm = int(float(nrof)/float(npes)*0.445) + 1
       nloops = 3
       minbas = nrof
       do nl=1,nloops
          maxbas = minbas - 1
          minbas = maxval(nupstrm)/(2**nl)
          if (nl == nloops) minbas = min(minbas,1)
          do nr=1,rtmlon*rtmlat
             if (gmask(nr) >= 2 .and. nupstrm(nr) > 0 .and. nupstrm(nr) >= minbas .and. nupstrm(nr) <= maxbas) then
                ! Decomp options
                !   find min pe (implemented but scales poorly)
                !   use increasing thresholds (implemented, ok load balance for l2r or calc)
                !   distribute basins using above methods but work from max to min basin size
                !
                !--------------
                ! find min pe
                !             baspe = 0
                !             do n = 1,npes-1
                !                if (nop(n) < nop(baspe)) baspe = n
                !             enddo
                !--------------
                ! find next pe below maxrtm threshhold and increment
                do while (nop(baspe) > maxrtm)
                   baspe = baspe + 1
                   if (baspe > npes-1) then
                      baspe = 0
                      maxrtm = max(maxrtm*1.5, maxrtm+1.0)   ! 3 loop, .445 and 1.5 chosen carefully
                   endif
                enddo
                !--------------
                if (baspe > npes-1 .or. baspe < 0) then
                   write(iulog,*) 'ERROR in decomp for MOSART ',nr,npes,baspe
                   call shr_sys_abort('ERROR mosart decomp')
                endif
                nop(baspe) = nop(baspe) + nupstrm(nr)
                nba(baspe) = nba(baspe) + 1
                pocn(nr) = baspe
             endif
          enddo ! nr
       enddo ! nl

       ! set pocn for land cells, was set for ocean above
       do nr=1,rtmlon*rtmlat
          if (idxocn(nr) > 0) then
             pocn(nr) = pocn(idxocn(nr))
             if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
                write(iulog,*) subname,' ERROR pocn lnd setting ',&
                   nr,idxocn(nr),idxocn(idxocn(nr)),pocn(idxocn(nr)),pocn(nr),npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    elseif (trim(decomp_option) == '1d') then
       ! distribute active points in 1d fashion to pes
       ! baspe is the pe assignment
       ! maxrtm is the maximum number of points to assign to each pe
       baspe = 0
       maxrtm = (nrof-1)/npes + 1
       do nr=1,rtmlon*rtmlat
          if (gmask(nr) >= 1) then
             pocn(nr) = baspe
             nop(baspe) = nop(baspe) + 1
             if (nop(baspe) >= maxrtm) then
                baspe = (mod(baspe+1,npes))
                if (baspe < 0 .or. baspe > npes-1) then
                   write(iulog,*) subname,' ERROR basepe ',baspe,npes
                   call shr_sys_abort(subname//' ERROR pocn lnd')
                endif
             endif
          endif
       enddo

    elseif (trim(decomp_option) == 'roundrobin') then
       ! distribute active points in roundrobin fashion to pes
       ! baspe is the pe assignment
       ! maxrtm is the maximum number of points to assign to each pe
       baspe = 0
       do nr=1,rtmlon*rtmlat
          if (gmask(nr) >= 1) then
             pocn(nr) = baspe
             nop(baspe) = nop(baspe) + 1
             baspe = (mod(baspe+1,npes))
             if (baspe < 0 .or. baspe > npes-1) then
                write(iulog,*) subname,' ERROR basepe ',baspe,npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    else
       write(iulog,*) subname,' ERROR decomp option unknown ',trim(decomp_option)
       call shr_sys_abort(subname//' ERROR pocn lnd')
    endif  ! decomp_option

    if (masterproc) then
       write(iulog,*) 'MOSART cells and basins total  = ',nrof,nbas
       write(iulog,*) 'MOSART cells per basin avg/max = ',nrof/nbas,maxval(nupstrm)
       write(iulog,*) 'MOSART cells per pe    min/max = ',minval(nop),maxval(nop)
       write(iulog,*) 'MOSART basins per pe   min/max = ',minval(nba),maxval(nba)
    endif

    deallocate(nupstrm)

    !-------------------------------------------------------
    !--- Count and distribute cells to rglo2gdc
    !-------------------------------------------------------

    rtmCTL%numr   = 0
    rtmCTL%lnumr  = 0

    do n = 0,npes-1
       if (iam == n) then
          rtmCTL%begr  = rtmCTL%numr  + 1
       endif
       rtmCTL%numr  = rtmCTL%numr  + nop(n)
       if (iam == n) then
          rtmCTL%lnumr = rtmCTL%lnumr + nop(n)
          rtmCTL%endr  = rtmCTL%begr  + rtmCTL%lnumr  - 1
       endif
    enddo

    allocate(rglo2gdc(rtmlon*rtmlat), &  !global mosart array
             nrs(0:npes-1))
    nrs = 0
    rglo2gdc = 0

    ! nrs is begr on each pe
    nrs(0) = 1
    do n = 1,npes-1
       nrs(n) = nrs(n-1) + nop(n-1)
    enddo

    ! reuse nba for nop-like counter here
    ! pocn -99 is unused cell
    nba = 0
    do nr = 1,rtmlon*rtmlat
       if (pocn(nr) >= 0) then
          rglo2gdc(nr) = nrs(pocn(nr)) + nba(pocn(nr))
          nba(pocn(nr)) = nba(pocn(nr)) + 1          
       endif
    enddo
    do n = 0,npes-1
       if (nba(n) /= nop(n)) then
          write(iulog,*) subname,' ERROR mosart cell count ',n,nba(n),nop(n)
          call shr_sys_abort(subname//' ERROR mosart cell count')
       endif
    enddo

    deallocate(nop,nba,nrs)
    deallocate(pocn)
    call t_stopf('mosarti_dec_distr')

    !-------------------------------------------------------
    !--- adjust area estimation from DRT algorithm for those outlet grids
    !--- useful for grid-based representation only
    !--- need to compute areas where they are not defined in input file
    !-------------------------------------------------------

    do n=1,rtmlon*rtmlat
       if (area_global(n) <= 0._r8) then
          i = mod(n-1,rtmlon) + 1
          j = (n-1)/rtmlon + 1
          dx = (rlone(i) - rlonw(i)) * deg2rad
          dy = sin(rlatn(j)*deg2rad) - sin(rlats(j)*deg2rad)
          area_global(n) = abs(1.e6_r8 * dx*dy*re*re)
          if (masterproc .and. area_global(n) <= 0) then
             write(iulog,*) 'Warning! Zero area for unit ', n, area_global(n),dx,dy,re
          end if
       end if
    end do

    call t_stopf('mosarti_decomp')

    !-------------------------------------------------------
    !--- Write per-processor runoff bounds depending on dbug level
    !-------------------------------------------------------

    call t_startf('mosarti_print')

    call shr_sys_flush(iulog)
    if (masterproc) then
       write(iulog,*) 'total runoff cells numr  = ',rtmCTL%numr
    endif
    call shr_sys_flush(iulog)
    call mpi_barrier(mpicom_rof,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    elseif (dbug == 3) then
       npint = 1
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)
       if (iam == pid) then
          write(iulog,'(2a,i9,a,i9,a,i9,a,i9)') &
             'MOSART decomp info',' proc = ',iam, &
             ' begr = ',rtmCTL%begr,&
             ' endr = ',rtmCTL%endr, &
             ' numr = ',rtmCTL%lnumr
       endif
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom_rof,ier)
    enddo

    call t_stopf('mosarti_print')

    !-------------------------------------------------------
    ! Allocate local flux variables
    !-------------------------------------------------------

    call t_startf('mosarti_vars')

    allocate (evel    (rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              flow    (rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroup_lagi(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroup_lagf(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erowm_regi(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erowm_regf(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroutup_avg(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erlat_avg(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' Allocation ERROR for flow'
       call shr_sys_abort(subname//' Allocationt ERROR flow')
    end if
    flow(:,:)        = 0._r8
    eroup_lagi(:,:)  = 0._r8
    eroup_lagf(:,:)  = 0._r8
    erowm_regi(:,:)  = 0._r8
    erowm_regf(:,:)  = 0._r8
    eroutup_avg(:,:) = 0._r8
    erlat_avg(:,:)   = 0._r8

    !-------------------------------------------------------
    ! Allocate runoff datatype 
    !-------------------------------------------------------

    call RunoffInit(rtmCTL%begr, rtmCTL%endr, rtmCTL%numr)

    !-------------------------------------------------------
    ! Initialize mosart flood - rtmCTL%fthresh and evel
    !-------------------------------------------------------

    if (do_rtmflood) then
       write(iulog,*) subname,' Flood not validated in this version, abort'
       call shr_sys_abort(subname//' Flood feature unavailable')
       call RtmFloodInit (frivinp_rtm, rtmCTL%begr, rtmCTL%endr, rtmCTL%fthresh, evel)
    else
       effvel(:) = effvel0  ! downstream velocity (m/s)
       rtmCTL%fthresh(:) = abs(spval)
       do nt = 1,nt_rtm
          do nr = rtmCTL%begr,rtmCTL%endr
             evel(nr,nt) = effvel(nt)
          enddo
       enddo
    end if

    !-------------------------------------------------------
    ! Initialize runoff data type
    !-------------------------------------------------------

    allocate(rgdc2glo(rtmCTL%numr), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' ERROR allocation of rgdc2glo'
       call shr_sys_abort(subname//' ERROR allocate of rgdc2glo')
    end if

    ! Set map from local to global index space
    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = rglo2gdc(n)
       if (nr > 0) then
          numr = numr + 1
          rgdc2glo(nr) = n
       endif
    end do
    end do
    if (numr /= rtmCTL%numr) then
       write(iulog,*) subname,'ERROR numr and rtmCTL%numr are different ',numr,rtmCTL%numr
       call shr_sys_abort(subname//' ERROR numr')
    endif

    ! Determine runoff datatype variables
    lrtmarea = 0.0_r8
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       rtmCTL%gindex(nr) = rgdc2glo(nr)
       rtmCTL%mask(nr) = gmask(rgdc2glo(nr))
       n = rgdc2glo(nr)
       i = mod(n-1,rtmlon) + 1
       j = (n-1)/rtmlon + 1
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(iulog,*) subname,' ERROR gdc2glo, nr,ng= ',nr,n
          call shr_sys_abort(subname//' ERROR gdc2glo values')
       endif
       rtmCTL%lonc(nr) = rtmCTL%rlon(i)
       rtmCTL%latc(nr) = rtmCTL%rlat(j)

       rtmCTL%outletg(nr) = idxocn(n)
       rtmCTL%area(nr) = area_global(n)
       lrtmarea = lrtmarea + rtmCTL%area(nr)
       if (dnID_global(n) <= 0) then
          rtmCTL%dsig(nr) = 0
       else
          if (rglo2gdc(dnID_global(n)) == 0) then
             write(iulog,*) subname,' ERROR glo2gdc dnID_global ',&
                  nr,n,dnID_global(n),rglo2gdc(dnID_global(n))
             call shr_sys_abort(subname//' ERROT glo2gdc dnID_global')
          endif
          cnt = cnt + 1
          rtmCTL%dsig(nr) = dnID_global(n)
       endif
    enddo
    deallocate(gmask)
    deallocate(rglo2gdc)
    deallocate(rgdc2glo)
    deallocate (dnID_global,area_global)
    deallocate(idxocn)
    call shr_mpi_sum(lrtmarea,rtmCTL%totarea,mpicom_rof,'mosart totarea',all=.true.)
    if (masterproc) write(iulog,*) subname,'  earth area ',4.0_r8*shr_const_pi*1.0e6_r8*re*re
    if (masterproc) write(iulog,*) subname,' MOSART area ',rtmCTL%totarea
    if (minval(rtmCTL%mask) < 1) then
       write(iulog,*) subname,'ERROR rtmCTL mask lt 1 ',minval(rtmCTL%mask),maxval(rtmCTL%mask)
       call shr_sys_abort(subname//' ERROR rtmCTL mask')
    endif

    !-------------------------------------------------------
    ! Compute Sparse Matrix for downstream advection
    !-------------------------------------------------------

    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat
    allocate(gindex(lsize))
    do nr = rtmCTL%begr,rtmCTL%endr
       gindex(nr-rtmCTL%begr+1) = rtmCTL%gindex(nr)
    enddo
    call mct_gsMap_init( gsMap_r, gindex, mpicom_rof, ROFID, lsize, gsize )
!   write(iulog,*) subname,' gsMap_r lsize = ',iam,mct_gsMap_lsize(gsMap_r,mpicom_rof)
    deallocate(gindex)

    if (smat_option == 'opt') then
       ! distributed smat initialization
       ! mct_sMat_init must be given the number of rows and columns that
       ! would be in the full matrix.  Nrows= size of output vector=nb.
       ! Ncols = size of input vector = na.

       cnt = 0
       do nr=rtmCTL%begr,rtmCTL%endr
          if(rtmCTL%dsig(nr) > 0) cnt = cnt + 1
       enddo

       call mct_sMat_init(sMat, gsize, gsize, cnt)
       igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
       igrow = mct_sMat_indexIA(sMat,'grow')   ! dst
       iwgt  = mct_sMat_indexRA(sMat,'weight')
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%dsig(nr) > 0) then
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%dsig(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          endif
       enddo

       call mct_sMatP_Init(sMatP_dnstrm, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

    elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then

       ! root initialization

       call mct_aVect_init(avtmp,rList='f1:f2',lsize=lsize)
       call mct_aVect_zero(avtmp)
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avtmp%rAttr(1,cnt) = rtmCTL%gindex(nr)
          avtmp%rAttr(2,cnt) = rtmCTL%dsig(nr)
       enddo
       call mct_avect_gather(avtmp,avtmpG,gsmap_r,mastertask,mpicom_rof)
       if (masterproc) then
          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
             endif
          enddo

          call mct_sMat_init(sMat, gsize, gsize, cnt)
          igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
          igrow = mct_sMat_indexIA(sMat,'grow')
          iwgt  = mct_sMat_indexRA(sMat,'weight')

          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(2,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             endif
          enddo
          call mct_avect_clean(avtmpG)
       else
          call mct_sMat_init(sMat,1,1,1)
       endif
       call mct_avect_clean(avtmp)

       call mct_sMatP_Init(sMatP_dnstrm, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

    else

       write(iulog,*) trim(subname),' MOSART ERROR: invalid smat_option '//trim(smat_option)
       call shr_sys_abort(trim(subname)//' ERROR invald smat option')

    endif

    ! initialize the AVs to go with sMatP
    write(rList,'(a,i3.3)') 'tr',1
    do nt = 2,nt_rtm
       write(rList,'(a,i3.3)') trim(rList)//':tr',nt
    enddo
!    write(iulog,*) trim(subname),' MOSART initialize avect ',trim(rList)
    call mct_aVect_init(avsrc_dnstrm,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avdst_dnstrm,rList=rList,lsize=rtmCTL%lnumr)
!    write(iulog,*) subname,' avsrc_dnstrm lsize = ',iam,mct_aVect_lsize(avsrc_dnstrm)
!    write(iulog,*) subname,' avdst_dnstrm lsize = ',iam,mct_aVect_lsize(avdst_dnstrm)

    lsize = mct_smat_gNumEl(sMatP_dnstrm%Matrix,mpicom_rof)
    if (masterproc) write(iulog,*) subname," Done initializing SmatP_dnstrm, nElements = ",lsize

    ! keep only sMatP
    call mct_sMat_clean(sMat)

    !-------------------------------------------------------
    ! Compute Sparse Matrix for direct to outlet transfer
    ! reuse gsmap_r
    !-------------------------------------------------------

    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat

    if (smat_option == 'opt') then
       ! distributed smat initialization
       ! mct_sMat_init must be given the number of rows and columns that
       ! would be in the full matrix.  Nrows= size of output vector=nb.
       ! Ncols = size of input vector = na.

       cnt = rtmCTL%endr - rtmCTL%begr + 1

       call mct_sMat_init(sMat, gsize, gsize, cnt)
       igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
       igrow = mct_sMat_indexIA(sMat,'grow')
       iwgt  = mct_sMat_indexRA(sMat,'weight')
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%outletg(nr) > 0) then
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%outletg(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          else
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%gindex(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          endif
       enddo
       if (cnt /= rtmCTL%endr - rtmCTL%begr + 1) then
          write(iulog,*) trim(subname),' MOSART ERROR: smat cnt1 ',cnt,rtmCTL%endr-rtmCTL%begr+1
          call shr_sys_abort(trim(subname)//' ERROR smat cnt1')
       endif

       call mct_sMatP_Init(sMatP_direct, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

    elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then

       ! root initialization

       call mct_aVect_init(avtmp,rList='f1:f2',lsize=lsize)
       call mct_aVect_zero(avtmp)
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avtmp%rAttr(1,cnt) = rtmCTL%gindex(nr)
          avtmp%rAttr(2,cnt) = rtmCTL%outletg(nr)
       enddo
       call mct_avect_gather(avtmp,avtmpG,gsmap_r,mastertask,mpicom_rof)
       if (masterproc) then

          cnt = rtmlon*rtmlat

          call mct_sMat_init(sMat, gsize, gsize, cnt)
          igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
          igrow = mct_sMat_indexIA(sMat,'grow')
          iwgt  = mct_sMat_indexRA(sMat,'weight')

          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(2,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             else
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(1,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             endif
          enddo
          if (cnt /= rtmlon*rtmlat) then
             write(iulog,*) trim(subname),' MOSART ERROR: smat cnt2 ',cnt,rtmlon*rtmlat
             call shr_sys_abort(trim(subname)//' ERROR smat cnt2')
          endif
          call mct_avect_clean(avtmpG)
       else
          call mct_sMat_init(sMat,1,1,1)
       endif
       call mct_avect_clean(avtmp)

       call mct_sMatP_Init(sMatP_direct, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

    else

       write(iulog,*) trim(subname),' MOSART ERROR: invalid smat_option '//trim(smat_option)
       call shr_sys_abort(trim(subname)//' ERROR invald smat option')

    endif

    ! initialize the AVs to go with sMatP
    write(rList,'(a,i3.3)') 'tr',1
    do nt = 2,nt_rtm
       write(rList,'(a,i3.3)') trim(rList)//':tr',nt
    enddo
!    write(iulog,*) trim(subname),' MOSART initialize avect ',trim(rList)
    call mct_aVect_init(avsrc_direct,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avdst_direct,rList=rList,lsize=rtmCTL%lnumr)
!    write(iulog,*) subname,' avsrc_direct lsize = ',iam,mct_aVect_lsize(avsrc_direct)
!    write(iulog,*) subname,' avdst_direct lsize = ',iam,mct_aVect_lsize(avdst_direct)

    lsize = mct_smat_gNumEl(sMatP_direct%Matrix,mpicom_rof)
    if (masterproc) write(iulog,*) subname," Done initializing SmatP_direct, nElements = ",lsize

    ! keep only sMatP
    call mct_sMat_clean(sMat)

    !-------------------------------------------------------
    ! Compute timestep and subcycling number
    !-------------------------------------------------------

! tcraig, old code based on cfl
!    dtover = 0._r8
!    dtovermax = 0._r8
!!    write(iulog,*) "tcx ddist ",minval(ddist),maxval(ddist)
!!    write(iulog,*) "tcx evel  ",minval(evel),maxval(evel)
!    do nt=1,nt_rtm
!       do nr=rtmCTL%begr,rtmCTL%endr
!          if (ddist(nr) /= 0._r8) then
!             dtover = evel(nr,nt)/ddist(nr)
!          else
!             dtover = 0._r8
!          endif
!          dtovermax = max(dtovermax,dtover)
!       enddo
!    enddo
!    dtover = dtovermax
!    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom_rof,ier)
!
!    write(iulog,*) "tcx dtover ",dtover,dtovermax
!
!    if (dtovermax > 0._r8) then
!       delt_mosart = (1.0_r8/dtovermax)*cfl_scale
!    else
!       write(iulog,*) subname,' ERROR in delt_mosart ',delt_mosart,dtover
!       call shr_sys_abort(subname//' ERROR delt_mosart')
!    endif
!
!    if (masterproc) write(iulog,*) 'mosart max timestep = ',delt_mosart,' (sec) for cfl_scale = ',cfl_scale
!    if (masterproc) call shr_sys_flush(iulog)
!
!    delt_mosart = 600._r8  ! here set the time interval for routing as 10 mins, which is sufficient for 1/8th degree resolution and coarser.
!    if (masterproc) write(iulog,*) 'mosart max timestep hardwired to = ',delt_mosart
!    if (masterproc) call shr_sys_flush(iulog)

    call t_stopf('mosarti_vars')

    !-------------------------------------------------------
    ! Initialize mosart
    !-------------------------------------------------------

    call t_startf('mosarti_mosart_init')

    !=== initialize MOSART related variables
!    if (masterproc) write(iulog,*) ' call mosart_init'
!    if (masterproc) call shr_sys_flush(iulog)
    call MOSART_init()

    call t_stopf('mosarti_mosart_init')

#ifdef INCLUDE_WRM
    call t_startf('mosarti_wrm_init')
    if (wrmflag) then
       call WRM_init()
    endif
    call t_startf('mosarti_wrm_init')
#endif

    !-------------------------------------------------------
    ! Read restart/initial info
    !-------------------------------------------------------

    call t_startf('mosarti_restart')

!    if (masterproc) write(iulog,*) ' call RtmRestFileRead'
!    if (masterproc) call shr_sys_flush(iulog)

    ! The call below opens and closes the file
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
        (nsrest == nsrBranch  )) then
       call RtmRestFileRead( file=fnamer )
       !write(iulog,*) ' MOSART init file is read'
       TRunoff%wh   = rtmCTL%wh
       TRunoff%wt   = rtmCTL%wt
       TRunoff%wr   = rtmCTL%wr
       TRunoff%erout= rtmCTL%erout
    else
!       do nt = 1,nt_rtm
!       do nr = rtmCTL%begr,rtmCTL%endr
!          TRunoff%wh(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 1.e-10_r8
!          TRunoff%wt(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 1.e-8_r8
!          TRunoff%wr(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 10._r8
!       enddo
!       enddo
       call WRM_computeRelease()
    endif

    do nt = 1,nt_rtm
    do nr = rtmCTL%begr,rtmCTL%endr
       call UpdateState_hillslope(nr,nt)
       call UpdateState_subnetwork(nr,nt)
       call UpdateState_mainchannel(nr,nt)
       rtmCTL%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + &
                             TRunoff%wh(nr,nt)*rtmCTL%area(nr))
    enddo
    enddo

    call t_stopf('mosarti_restart')

    !-------------------------------------------------------
    ! Initialize mosart history handler and fields
    !-------------------------------------------------------

    call t_startf('mosarti_histinit')

!    if (masterproc) write(iulog,*) ' call RtmHistFldsInit'
!    if (masterproc) call shr_sys_flush(iulog)

    call RtmHistFldsInit()
    if (nsrest==nsrStartup .or. nsrest==nsrBranch) then
       call RtmHistHtapesBuild()
    end if
    call RtmHistFldsSet()

    if (masterproc) write(iulog,*) subname,' done'
    if (masterproc) call shr_sys_flush(iulog)

    call t_stopf('mosarti_histinit')

  end subroutine Rtmini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmrun
!
! !INTERFACE:
  subroutine Rtmrun(rstwr,nlend,rdate)
!
! !DESCRIPTION:
! River routing model
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    logical ,         intent(in) :: rstwr          ! true => write restart file this step)
    logical ,         intent(in) :: nlend          ! true => end of run on this step
    character(len=*), intent(in) :: rdate          ! restart file time stamp for name
!
! !CALLED FROM:
! subroutine RtmMap in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i, j, n, nr, ns, nt, n2, nf, idam ! indices
    integer, parameter :: budget_terms_total = 80
    logical  :: output_all_budget_terms = .false.   ! output flag
    real(r8) :: budget_terms (budget_terms_total,nt_rtm)    ! local budget sums
    real(r8) :: budget_global(budget_terms_total,nt_rtm)    ! global budget sums
    real(r8),save :: budget_accum(nt_rtm)   ! BUDGET accumulator over run
    integer ,save :: budget_accum_cnt       ! counter for budget_accum
    real(r8) :: budget_input, budget_output, budget_volume, budget_total, &
                budget_other
    logical  :: budget_check, budget_write  ! do global budget check

    ! BUDGET term ids
    ! budget computed in m3 over each coupling period
    ! use delt_coupling to convert from rates to volumes
    ! General equation is as follows
    !   Vf = Vi + input - output + other
    !     or
    !   Vf-Vi - input + output - other == 0 for conservation
    ! budget_volume = Vf-Vi
    ! budget_input = input
    ! budget_output = output
    ! budget_other = other terms
    ! budget_total = budget_volume - budget_input + budget_output - budget_other
    ! bv_ generally accumulates a volume (m3)
    ! br_ generally accumulates a rate (m3/s)

    ! Storage/Volume TERMS (volumes, m3)
    integer,parameter :: bv_volt_i = 1  ! initial total volume
    integer,parameter :: bv_volt_f = 2  ! final   total volume
    integer,parameter :: bv_wt_i   = 3  ! initial wt volume
    integer,parameter :: bv_wt_f   = 4  ! final   wt volume
    integer,parameter :: bv_wr_i   = 5  ! initial wr volume
    integer,parameter :: bv_wr_f   = 6  ! final   wr volume
    integer,parameter :: bv_wh_i   = 7  ! initial wh volume
    integer,parameter :: bv_wh_f   = 8  ! final   wh volume
    integer,parameter :: bv_dstor_i= 9  ! initial dam storage
    integer,parameter :: bv_dstor_f= 10 ! final   dam storage

    ! Input TERMS (rates, m3/s)
    integer,parameter :: br_qsur   = 20 ! input qsur
    integer,parameter :: br_qsub   = 21 ! input qsub
    integer,parameter :: br_qgwl   = 22 ! input qgwl
    integer,parameter :: br_qdto   = 23 ! input qdto

    ! Output TERMS (rates m3/s or volumes m3)
    integer,parameter :: br_ocnout = 40 ! runoff output to ocean
    integer,parameter :: br_lndout = 41 ! runoff output on non ocean points
    integer,parameter :: br_flood  = 42 ! flood term back to land
    integer,parameter :: br_direct = 43 ! direct output term
    integer,parameter :: bv_dsupp_i= 44 ! initial dam supply
    integer,parameter :: bv_dsupp_f= 45 ! final   dam supply

    ! Other Diagnostic TERMS (rates, m3/s)
    integer,parameter :: br_erolpo = 60 ! erout lag ocn previous
    integer,parameter :: br_erolco = 61 ! erout lag ocn current
    integer,parameter :: br_erorpo = 62 ! erout lag ocn previous
    integer,parameter :: br_erorco = 63 ! erout lag ocn current
    integer,parameter :: br_eroutup= 64 ! erout upstream average
    integer,parameter :: br_erolpn = 65 ! erout lag non-ocn previous
    integer,parameter :: br_erolcn = 66 ! erout lag non-ocn current
    integer,parameter :: br_erorpn = 67 ! erout lag non-ocn previous
    integer,parameter :: br_erorcn = 68 ! erout lag non-ocn current
    integer,parameter :: br_erlat  = 69 ! erlateral 

    ! Accumuluation TERMS
    integer,parameter :: bv_naccum = 80 ! accumulated net budget

    !   volume = 2 - 1 + bv_dstor_f - bv_dstor_i
    !   input  = br_qsur + br_qsub + br_qgwl + br_qdto
    !   output = br_ocnout + br_flood + br_direct + 42
    !   total  = volume - input + output
    !   erlag  = br_erolpn - br_erolcn

    real(r8) :: volr_init                   ! temporary storage to compute dvolrdt
    real(r8),parameter :: budget_tolerance = 1.0e-6   ! budget tolerance, m3/day
    logical  :: abort                       ! abort flag
    real(r8) :: sum1,sum2
    integer  :: yr, mon, day, ymd, tod      ! time information
    integer  :: nsub                        ! subcyling for cfl
    real(r8) :: delt                        ! delt associated with subcycling
    real(r8) :: delt_coupling               ! real value of coupling_period
    integer , save :: nsub_save             ! previous nsub
    real(r8), save :: delt_save = -1.0      ! previous delt
    logical , save :: first_call = .true.   ! first time flag (for backwards compatibility)
    character(len=256) :: filer             ! restart file name
    integer  :: cnt                         ! counter for gridcells
    integer  :: ier                         ! error code
    integer,parameter  :: dbug = 1          ! local debug flag
!scs
! parameters used in negative runoff partitioning algorithm
    real(r8) :: river_volume_minimum        ! gridcell area multiplied by average river_depth_minimum [m3]
    real(r8) :: qgwl_volume                 ! volume of runoff during time step [m3]
!scs
    character(len=*),parameter :: subname = '(Rtmrun) '
!-----------------------------------------------------------------------

    call t_startf('mosartr_tot')
    call shr_sys_flush(iulog)

    delt_coupling = coupling_period*1.0_r8
    budget_check = .true.   ! leave this on all the time for now
    if (first_call) then
       budget_accum = 0._r8
       budget_accum_cnt = 0
       if (masterproc) write(iulog,'(2a,g20.12)') trim(subname),' MOSART coupling period ',delt_coupling
    end if

    budget_terms = 0._r8

    flow = 0._r8
    eroup_lagi = 0._r8
    eroup_lagf = 0._r8
    erowm_regi = 0._r8
    erowm_regf = 0._r8
    eroutup_avg = 0._r8
    erlat_avg = 0._r8
    rtmCTL%runoff = 0._r8
    rtmCTL%direct = 0._r8
    rtmCTL%flood = 0._r8
    rtmCTL%runofflnd = spval
    rtmCTL%runoffocn = spval
    rtmCTL%dvolrdt = 0._r8
    rtmCTL%dvolrdtlnd = spval
    rtmCTL%dvolrdtocn = spval

    if (budget_check) then
       call t_startf('mosartr_budget')
       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(bv_volt_i,nt) = budget_terms( bv_volt_i,nt) + rtmCTL%volr(nr,nt)
          budget_terms(bv_wt_i,nt) = budget_terms(bv_wt_i,nt) + TRunoff%wt(nr,nt)
          budget_terms(bv_wr_i,nt) = budget_terms(bv_wr_i,nt) + TRunoff%wr(nr,nt)
          budget_terms(bv_wh_i,nt) = budget_terms(bv_wh_i,nt) + TRunoff%wh(nr,nt)*rtmCTL%area(nr)
          budget_terms(br_qsur,nt) = budget_terms(br_qsur,nt) + rtmCTL%qsur(nr,nt)*delt_coupling
          budget_terms(br_qsub,nt) = budget_terms(br_qsub,nt) + rtmCTL%qsub(nr,nt)*delt_coupling
          budget_terms(br_qgwl,nt) = budget_terms(br_qgwl,nt) + rtmCTL%qgwl(nr,nt)*delt_coupling
          budget_terms(br_qdto,nt) = budget_terms(br_qdto,nt) + rtmCTL%qdto(nr,nt)*delt_coupling
       enddo
       enddo

#ifdef INCLUDE_WRM
       if (wrmflag) then
          StorWater%supply = 0._r8
          nt = 1
          do nr = rtmCTL%begr,rtmCTL%endr
             budget_terms(bv_dsupp_i,nt) = budget_terms(bv_dsupp_i,nt) + StorWater%supply(nr)
          enddo
          do idam = 1,ctlSubwWRM%LocalNumDam
             budget_terms(bv_dstor_i,nt) = budget_terms(bv_dstor_i,nt) + StorWater%storage(idam)
          enddo
       endif
#endif
       call t_stopf('mosartr_budget')
    endif ! budget_check

    ! data for euler solver, in m3/s here
    do nr = rtmCTL%begr,rtmCTL%endr
    do nt = 1,nt_rtm
       TRunoff%qsur(nr,nt) = rtmCTL%qsur(nr,nt)
       TRunoff%qsub(nr,nt) = rtmCTL%qsub(nr,nt)
       TRunoff%qgwl(nr,nt) = rtmCTL%qgwl(nr,nt)
    enddo
    enddo

    !-----------------------------------
    ! Compute flood
    ! Remove water from mosart and send back to clm
    ! Just consider land points and only remove liquid water 
    ! rtmCTL%flood is m3/s here
    !-----------------------------------

    call t_startf('mosartr_flood')
    nt = 1 
    rtmCTL%flood = 0._r8
    do nr = rtmCTL%begr,rtmCTL%endr
       ! initialize rtmCTL%flood to zero
       if (rtmCTL%mask(nr) == 1) then
          if (rtmCTL%volr(nr,nt) > rtmCTL%fthresh(nr)) then 
             ! determine flux that is sent back to the land
             ! this is in m3/s
             rtmCTL%flood(nr) = &
                  (rtmCTL%volr(nr,nt)-rtmCTL%fthresh(nr)) / (delt_coupling)

             ! rtmCTL%flood will be sent back to land - so must subtract this 
             ! from the input runoff from land
             ! tcraig, comment - this seems like an odd approach, you
             !   might create negative forcing.  why not take it out of
             !   the volr directly?  it's also odd to compute this
             !   at the initial time of the time loop.  why not do
             !   it at the end or even during the run loop as the
             !   new volume is computed.  fluxout depends on volr, so
             !   how this is implemented does impact the solution.
             TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) - rtmCTL%flood(nr)
          endif
       endif
    enddo
    call t_stopf('mosartr_flood')

    !-----------------------------------
    ! DIRECT sMAT transfer to outlet point using sMat
    ! Remember to subract water from TRunoff forcing
    !-----------------------------------

    if (barrier_timers) then
       call t_startf('mosartr_SMdirect_barrier')
       call mpi_barrier(mpicom_rof,ier)
       call t_stopf ('mosartr_SMdirect_barrier')
    endif

    call t_startf('mosartr_SMdirect')

    !--- copy direct transfer fields to AV
    !--- convert kg/m2s to m3/s
    call mct_avect_zero(avsrc_direct)
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       cnt = cnt + 1
       do nt = 1,nt_rtm

          !-------------------------------
          ! This water is routed directly to the outlet and passed out
          ! Turn on and off terms as desired via commenting them out
          !-------------------------------

          !---- all dto water, none was going to TRunoff ---
          !---- *** DO NOT TURN THIS ONE OFF, conservation will fail *** ---
          avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + rtmCTL%qdto(nr,nt)

          !---- negative gwl water less than channel volume, remove from TRunoff ---
          !---- scs
          qgwl_volume = TRunoff%qgwl(nr,nt) * delt_mosart
          river_volume_minimum = river_depth_minimum * rtmCTL%area(nr)
          !---- if qgwl is negative, and adding it to the main channel would bring 
          !---- main channel storage below a threshold, send qgwl directly to ocean
          if (((qgwl_volume + TRunoff%wr(nr,nt)) < river_volume_minimum) &
             .and. (TRunoff%qgwl(nr,nt) < 0._r8)) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
             TRunoff%qgwl(nr,nt) = 0._r8
          endif
          !---- scs

          !---- negative qgwl water, remove from TRunoff ---
          if (TRunoff%qgwl(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
             TRunoff%qgwl(nr,nt) = 0._r8
          endif

          !---- all gwl water, remove from TRunoff ---
          avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
          TRunoff%qgwl(nr,nt) = 0._r8

          !---- negative qsub water, remove from TRunoff ---
          if (TRunoff%qsub(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qsub(nr,nt)
             TRunoff%qsub(nr,nt) = 0._r8
          endif

          !---- negative qsur water, remove from TRunoff ---
          if (TRunoff%qsur(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qsur(nr,nt)
             TRunoff%qsur(nr,nt) = 0._r8
          endif

          !---- water outside the basin ---
          !---- *** DO NOT TURN THIS ONE OFF, conservation will fail *** ---
          if (TUnit%mask(nr) > 0) then
             ! mosart euler
          else
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + &
                TRunoff%qsub(nr,nt) + &
                TRunoff%qsur(nr,nt) + &
                TRunoff%qgwl(nr,nt)
             TRunoff%qsub(nr,nt) = 0._r8
             TRunoff%qsur(nr,nt) = 0._r8
             TRunoff%qgwl(nr,nt) = 0._r8
          endif

          !---- all nt=2 water ---
          !---- can turn off euler_calc for this tracer ----
          if (nt == 2) then
             TUnit%euler_calc(nt) = .false.
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + &
                TRunoff%qsub(nr,nt) + &
                TRunoff%qsur(nr,nt) + &
                TRunoff%qgwl(nr,nt)
             TRunoff%qsub(nr,nt) = 0._r8
             TRunoff%qsur(nr,nt) = 0._r8
             TRunoff%qgwl(nr,nt) = 0._r8
          endif

       enddo
    enddo
    call mct_avect_zero(avdst_direct)

    call mct_sMat_avMult(avsrc_direct, sMatP_direct, avdst_direct)

    !--- copy direct transfer water from AV to output field ---
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       cnt = cnt + 1
       do nt = 1,nt_rtm
          rtmCTL%direct(nr,nt) = avdst_direct%rAttr(nt,cnt)
       enddo
    enddo
    call t_stopf('mosartr_SMdirect')

    !-----------------------------------
    ! MOSART Subcycling
    !-----------------------------------

    call t_startf('mosartr_subcycling')

    if (first_call) then
       do nt = 1,nt_rtm
          write(iulog,'(2a,i6,l4)') trim(subname),' euler_calc for nt = ',nt,TUnit%euler_calc(nt)
       enddo
    endif

    nsub = coupling_period/delt_mosart
    if (nsub*delt_mosart < coupling_period) then
       nsub = nsub + 1
    end if
    delt = delt_coupling/float(nsub)
    if (delt /= delt_save) then
       if (masterproc) then
          write(iulog,'(2a,2g20.12,2i12)') trim(subname),' MOSART delt update from/to',delt_save,delt,nsub_save,nsub
       end if
    endif

    nsub_save = nsub
    delt_save = delt
    Tctl%DeltaT = delt

    !-----------------------------------
    ! mosart euler solver
    ! --- convert TRunoff fields from m3/s to m/s before calling Euler
    !-----------------------------------

    do nt = 1,nt_rtm
    do nr = rtmCTL%begr,rtmCTL%endr
       TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) / rtmCTL%area(nr)
       TRunoff%qsub(nr,nt) = TRunoff%qsub(nr,nt) / rtmCTL%area(nr)
       TRunoff%qgwl(nr,nt) = TRunoff%qgwl(nr,nt) / rtmCTL%area(nr)
    enddo
    enddo

    do ns = 1,nsub

       ! this advances the model time one coupling period when ns=1
       ! we want the advance in the subcycling because of the new_month flag
       ! this should really be advancing delt per ns, not coupling_period per call, but
       ! need to modify clock initialization and restart for that to be supported
       ! and make sure a rational timestep doesn't lead to accumulating errors
       call advance_timestep(ns)

       call get_curr_date(yr, mon, day, tod)
       ymd = yr*10000 + mon*100 + day
       if (tod == 0 .and. masterproc) then
          write(iulog,*) ' '
          write(iulog,'(2a,i4,a,i10,i6)') trim(subname),' subcycling=',ns,': model date=',ymd,tod
       endif

       call t_startf('mosartr_euler')
       call Euler()
       call t_stopf('mosartr_euler')

! tcraig - NOT using this now, but leave it here in case it's useful in the future
!   for some runoff terms.
!       !-----------------------------------
!       ! downstream advection using sMat
!       !-----------------------------------
!
!       if (barrier_timers) then
!          call t_startf('mosartr_SMdnstrm_barrier')
!          call mpi_barrier(mpicom_rof,ier)
!          call t_stopf ('mosartr_SMdnstrm_barrier')
!       endif
!
!       call t_startf('mosartr_SMdnstrm')
!
!       !--- copy fluxout into avsrc_dnstrm ---
!       cnt = 0
!       do n = rtmCTL%begr,rtmCTL%endr
!          cnt = cnt + 1
!          do nt = 1,nt_rtm
!             avsrc_dnstrm%rAttr(nt,cnt) = fluxout(n,nt)
!          enddo
!       enddo
!       call mct_avect_zero(avdst_dnstrm)
!
!       call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
!
!       !--- add mapped fluxout to sfluxin ---
!       cnt = 0
!       sfluxin = 0._r8
!       do n = rtmCTL%begr,rtmCTL%endr
!          cnt = cnt + 1
!          do nt = 1,nt_rtm
!             sfluxin(n,nt) = sfluxin(n,nt) + avdst_dnstrm%rAttr(nt,cnt)
!          enddo
!       enddo
!       call t_stopf('mosartr_SMdnstrm')

       !-----------------------------------
       ! accumulate local flow field
       !-----------------------------------

       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          flow(nr,nt) = flow(nr,nt) + TRunoff%flow(nr,nt)
          eroup_lagi(nr,nt) = eroup_lagi(nr,nt) + TRunoff%eroup_lagi(nr,nt)
          eroup_lagf(nr,nt) = eroup_lagf(nr,nt) + TRunoff%eroup_lagf(nr,nt)
          erowm_regi(nr,nt) = erowm_regi(nr,nt) + TRunoff%erowm_regi(nr,nt)
          erowm_regf(nr,nt) = erowm_regf(nr,nt) + TRunoff%erowm_regf(nr,nt)
          eroutup_avg(nr,nt) = eroutup_avg(nr,nt) + TRunoff%eroutup_avg(nr,nt)
          erlat_avg(nr,nt) = erlat_avg(nr,nt) + TRunoff%erlat_avg(nr,nt)
       enddo
       enddo

    enddo ! nsub

    !-----------------------------------
    ! average flow over subcycling
    !-----------------------------------

    flow        = flow        / float(nsub)
    eroup_lagi  = eroup_lagi  / float(nsub)
    eroup_lagf  = eroup_lagf  / float(nsub)
    erowm_regi  = erowm_regi  / float(nsub)
    erowm_regf  = erowm_regf  / float(nsub)
    eroutup_avg = eroutup_avg / float(nsub)
    erlat_avg   = erlat_avg   / float(nsub)

    !-----------------------------------
    ! update states when subsycling completed
    !-----------------------------------

    rtmCTL%wh      = TRunoff%wh
    rtmCTL%wt      = TRunoff%wt
    rtmCTL%wr      = TRunoff%wr
    rtmCTL%erout   = TRunoff%erout

    do nt = 1,nt_rtm
    do nr = rtmCTL%begr,rtmCTL%endr
       volr_init = rtmCTL%volr(nr,nt)
       rtmCTL%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + &
                             TRunoff%wh(nr,nt)*rtmCTL%area(nr))
       rtmCTL%dvolrdt(nr,nt) = (rtmCTL%volr(nr,nt) - volr_init) / delt_coupling
       rtmCTL%runoff(nr,nt) = flow(nr,nt)

       rtmCTL%runofftot(nr,nt) = rtmCTL%direct(nr,nt)
       if (rtmCTL%mask(nr) == 1) then
          rtmCTL%runofflnd(nr,nt) = rtmCTL%runoff(nr,nt)
          rtmCTL%dvolrdtlnd(nr,nt)= rtmCTL%dvolrdt(nr,nt)
       elseif (rtmCTL%mask(nr) >= 2) then
          rtmCTL%runoffocn(nr,nt) = rtmCTL%runoff(nr,nt)
          rtmCTL%runofftot(nr,nt) = rtmCTL%runofftot(nr,nt) + rtmCTL%runoff(nr,nt)
          rtmCTL%dvolrdtocn(nr,nt)= rtmCTL%dvolrdt(nr,nt)
       endif
    enddo
    enddo

    call t_stopf('mosartr_subcycling')

    !-----------------------------------
    ! BUDGET
    !-----------------------------------

    budget_write = .false.
    if (day == 1 .and. mon == 1) budget_write = .true.
    if (tod == 0) budget_write = .true.

    if (budget_check) then
       call t_startf('mosartr_budget')
       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(bv_volt_f,nt) = budget_terms(bv_volt_f,nt) + rtmCTL%volr(nr,nt)
          budget_terms(bv_wt_f,nt) = budget_terms(bv_wt_f,nt) + TRunoff%wt(nr,nt)
          budget_terms(bv_wr_f,nt) = budget_terms(bv_wr_f,nt) + TRunoff%wr(nr,nt)
          budget_terms(bv_wh_f,nt) = budget_terms(bv_wh_f,nt) + TRunoff%wh(nr,nt)*rtmCTL%area(nr)
          budget_terms(br_direct,nt) = budget_terms(br_direct,nt) + rtmCTL%direct(nr,nt)*delt_coupling
          if (rtmCTL%mask(nr) >= 2) then
             budget_terms(br_ocnout,nt) = budget_terms(br_ocnout,nt) + rtmCTL%runoff(nr,nt)*delt_coupling
             budget_terms(br_erolpo,nt) = budget_terms(br_erolpo,nt) + eroup_lagi(nr,nt)*delt_coupling
             budget_terms(br_erolco,nt) = budget_terms(br_erolco,nt) + eroup_lagf(nr,nt)*delt_coupling
             budget_terms(br_erorpo,nt) = budget_terms(br_erorpo,nt) + erowm_regi(nr,nt)*delt_coupling
             budget_terms(br_erorco,nt) = budget_terms(br_erorco,nt) + erowm_regf(nr,nt)*delt_coupling
          else
             budget_terms(br_lndout,nt) = budget_terms(br_lndout,nt) + rtmCTL%runoff(nr,nt)*delt_coupling
             budget_terms(br_erolpn,nt) = budget_terms(br_erolpn,nt) + eroup_lagi(nr,nt)*delt_coupling
             budget_terms(br_erolcn,nt) = budget_terms(br_erolcn,nt) + eroup_lagf(nr,nt)*delt_coupling
             budget_terms(br_erorpn,nt) = budget_terms(br_erorpn,nt) + erowm_regi(nr,nt)*delt_coupling
             budget_terms(br_erorcn,nt) = budget_terms(br_erorcn,nt) + erowm_regf(nr,nt)*delt_coupling
          endif
          budget_terms(br_eroutup,nt) = budget_terms(br_eroutup,nt) - eroutup_avg(nr,nt)*delt_coupling
          budget_terms(br_erlat,nt) = budget_terms(br_erlat,nt) - erlat_avg(nr,nt)*delt_coupling
       enddo
       enddo
       nt = 1
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(br_flood,nt) = budget_terms(br_flood,nt) + rtmCTL%flood(nr)*delt_coupling
       enddo
#ifdef INCLUDE_WRM
       if (wrmflag) then
          nt = 1
          do nr = rtmCTL%begr,rtmCTL%endr
             budget_terms(bv_dsupp_f,nt) = budget_terms(bv_dsupp_f,nt) + StorWater%supply(nr)
          enddo
          do idam = 1,ctlSubwWRM%LocalNumDam
             budget_terms(bv_dstor_f,nt) = budget_terms(bv_dstor_f,nt) + StorWater%storage(idam)
          enddo
       endif
#endif

       ! accumulate the budget total over the run to make sure it's decreasing on avg
       budget_accum_cnt = budget_accum_cnt + 1
       do nt = 1,nt_rtm
          budget_volume =  budget_terms(bv_volt_f,nt) - budget_terms(bv_volt_i,nt) + &
                           budget_terms(bv_dstor_f,nt) - budget_terms(bv_dstor_i,nt)
          budget_input  =  budget_terms(br_qsur,nt) + budget_terms(br_qsub,nt) + &
                           budget_terms(br_qgwl,nt) + budget_terms(br_qdto,nt)
          budget_output =  budget_terms(br_ocnout,nt) + budget_terms(br_flood,nt) + &
                           budget_terms(br_direct,nt) + &
                           budget_terms(bv_dsupp_f,nt) - budget_terms(bv_dsupp_i,nt)
          budget_accum(nt) = budget_accum(nt) + budget_volume - budget_input + budget_output
          budget_terms(bv_naccum,nt) = budget_accum(nt)/budget_accum_cnt
       enddo
       call t_stopf('mosartr_budget')
    endif  ! budget_check

    if (budget_check .and. budget_write) then
       call t_startf('mosartr_budget')
       !--- check budget

       ! convert terms from m3 to million m3
       budget_terms(:,:) = budget_terms(:,:) * 1.0e-6_r8

       ! global sum
       call shr_mpi_sum(budget_terms,budget_global,mpicom_rof,'mosart global budget',all=.false.)

       ! write budget
       if (masterproc) then
          write(iulog,'(2a,i10,i6)') trim(subname),' MOSART BUDGET diagnostics (million m3) for ',ymd,tod
          do nt = 1,nt_rtm
            budget_volume = (budget_global(bv_volt_f,nt) - budget_global(bv_volt_i,nt) + &
                             budget_global(bv_dstor_f,nt) - budget_global(bv_dstor_i,nt))
            budget_input  = (budget_global(br_qsur,nt) + budget_global(br_qsub,nt) + &
                             budget_global(br_qgwl,nt) + budget_global(br_qdto,nt))
            budget_output = (budget_global(br_ocnout,nt) + budget_global(br_flood,nt) + &
                             budget_global(br_direct,nt) + &
                             budget_global(bv_dsupp_f,nt) - budget_global(bv_dsupp_i,nt))
            ! erout lag, need to remove current term and add in previous term, current term used in next timestep
            budget_other  = budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt) + &
                            budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)
            budget_total  = budget_volume - budget_input + budget_output - budget_other

            write(iulog,'(2a)') trim(subname),'-----------------------------------------------------------------'
            write(iulog,'(2a,i4)')        trim(subname),'  tracer = ',nt
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wh    = ',nt,budget_global(bv_wh_f,nt)-budget_global(bv_wh_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wt    = ',nt,budget_global(bv_wt_f,nt)-budget_global(bv_wt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wr    = ',nt,budget_global(bv_wr_f,nt)-budget_global(bv_wr_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume dstor = ',nt,budget_global(bv_dstor_f,nt)-budget_global(bv_dstor_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * dvolume total = ',nt,budget_volume
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x dvolume check = ',nt,budget_volume - &
                                                                             (budget_global(bv_wh_f,nt)-budget_global(bv_wh_i,nt) + &
                                                                              budget_global(bv_wt_f,nt)-budget_global(bv_wt_i,nt) + &
                                                                              budget_global(bv_wr_f,nt)-budget_global(bv_wr_i,nt) + &
                                                                              budget_global(bv_dstor_f,nt)-budget_global(bv_dstor_i,nt)),&
                                                                              ' (should be zero)'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volume   init = ',nt,budget_global(bv_volt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volume  final = ',nt,budget_global(bv_volt_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumeh  init = ',nt,budget_global(bv_wh_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumeh final = ',nt,budget_global(bv_wh_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumet  init = ',nt,budget_global(bv_wt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumet final = ',nt,budget_global(bv_wt_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumer  init = ',nt,budget_global(bv_wr_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumer final = ',nt,budget_global(bv_wr_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x storage  init = ',nt,budget_global(bv_dstor_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x storage final = ',nt,budget_global(bv_dstor_f,nt)
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input surface = ',nt,budget_global(br_qsur,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input subsurf = ',nt,budget_global(br_qsub,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input gwl     = ',nt,budget_global(br_qgwl,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input dto     = ',nt,budget_global(br_qdto,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * input total   = ',nt,budget_input
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x input check   = ',nt,budget_input - &
                                                                             (budget_global(br_qsur,nt)+budget_global(br_qsub,nt)+ &
                                                                              budget_global(br_qgwl,nt)+budget_global(br_qdto,nt)), &
                                                                             ' (should be zero)'
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output runoff = ',nt,budget_global(br_ocnout,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output direct = ',nt,budget_global(br_direct,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output flood  = ',nt,budget_global(br_flood,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output supply = ',nt,budget_global(bv_dsupp_f,nt)-budget_global(bv_dsupp_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * output total  = ',nt,budget_output
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x output check  = ',nt,budget_output - &
                                                                             (budget_global(br_ocnout,nt) + budget_global(br_direct,nt) + &
                                                                              budget_global(br_flood,nt) + &
                                                                              budget_global(bv_dsupp_f,nt)-budget_global(bv_dsupp_i,nt)), &
                                                                             ' (should be zero)'
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   other dwn lag = ',nt,budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   other reg lnd = ',nt,budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * other total   = ',nt,budget_other
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x other check   = ',nt,budget_other - &
                                                                            (budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt) + &
                                                                             budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)), &
                                                                            ' (should be zero)'
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum dvolume   = ',nt,budget_volume
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum input     = ',nt,budget_input
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum output    = ',nt,budget_output
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum other     = ',nt,budget_other 
          endif
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' * sum budget ** = ',nt,budget_total,' (should be zero, dv-in+out-oth)'
          if (output_all_budget_terms) then
            ! accum budget is just dv-i+o and should show that over time, the other terms go to zero (lag yes, reg land no)
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x accum budget  = ',nt,budget_global(bv_naccum,nt),' (should tend to zero over run, dv-in+out)'
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x runoff     ocn= ',nt,budget_global(br_ocnout,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagi ocn= ',nt,budget_global(br_erolpo,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagf ocn= ',nt,budget_global(br_erolco,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regi ocn= ',nt,budget_global(br_erorpo,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regf ocn= ',nt,budget_global(br_erorco,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x other reg ocn = ',nt,budget_global(br_erorpo,nt) - budget_global(br_erorco,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x runoff     lnd= ',nt,budget_global(br_lndout,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagi lnd= ',nt,budget_global(br_erolpn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagf lnd= ',nt,budget_global(br_erolcn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regi lnd= ',nt,budget_global(br_erorpn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regf lnd= ',nt,budget_global(br_erorcn,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroutup_avg   = ',nt,budget_global(br_eroutup,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erlateral     = ',nt,budget_global(br_erlat,nt)
            write(iulog,'(2a)') trim(subname),'----------------'
          endif

            if ((budget_total) > 1.0e-6) then
               write(iulog,'(2a,i4)') trim(subname),' ***** BUDGET WARNING error gt 1. m3 for nt = ',nt
            endif
          enddo
          write(iulog,'(a)') '----------------------------------- '
       endif

       call t_stopf('mosartr_budget')
    endif  ! budget_write

    !-----------------------------------
    ! Write out MOSART history file
    !-----------------------------------

    call t_startf('mosartr_hbuf')
    call RtmHistFldsSet()
    call RtmHistUpdateHbuf()
    call t_stopf('mosartr_hbuf')

    call t_startf('mosartr_htapes')
    call RtmHistHtapesWrapup( rstwr, nlend )
    call t_stopf('mosartr_htapes')

    !-----------------------------------
    ! Write out MOSART restart file
    !-----------------------------------

    if (rstwr) then
       call t_startf('mosartr_rest')
       filer = RtmRestFileName(rdate=rdate)
       call RtmRestFileWrite( filer, rdate=rdate )
       call t_stopf('mosartr_rest')
    end if

    !-----------------------------------
    ! Done
    !-----------------------------------

    first_call = .false.

    call shr_sys_flush(iulog)
    call t_stopf('mosartr_tot')

  end subroutine Rtmrun

!-----------------------------------------------------------------------

  subroutine RtmFloodInit(frivinp, begr, endr, fthresh, evel )

    !-----------------------------------------------------------------------
    ! Uses

    ! Input variables
    character(len=*), intent(in) :: frivinp
    integer , intent(in)  :: begr, endr
    real(r8), intent(out) :: fthresh(begr:endr)
    real(r8), intent(out) :: evel(begr:endr,nt_rtm) 

    ! Local variables
    real(r8) , pointer :: rslope(:)   
    real(r8) , pointer :: max_volr(:)
    integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
    integer :: nt,n,cnt              ! indices
    logical :: readvar               ! read variable in or not
    integer :: ier                   ! status variable
    integer :: dids(2)               ! variable dimension ids 
    type(file_desc_t)  :: ncid       ! pio file desc
    type(var_desc_t)   :: vardesc    ! pio variable desc 
    type(io_desc_t)    :: iodesc     ! pio io desc
    character(len=256) :: locfn      ! local file name

    !MOSART Flood variables for spatially varying celerity
    real(r8) :: effvel(nt_rtm) = 0.7_r8   ! downstream velocity (m/s)
    real(r8) :: min_ev(nt_rtm) = 0.35_r8  ! minimum downstream velocity (m/s)
    real(r8) :: fslope = 1.0_r8           ! maximum slope for which flooding can occur
    character(len=*),parameter :: subname = '(RtmFloodInit) '
    !-----------------------------------------------------------------------

    allocate(rslope(begr:endr), max_volr(begr:endr), stat=ier)
    if (ier /= 0) call shr_sys_abort(subname // ' allocation ERROR')

    ! Assume that if SLOPE is on river input dataset so is MAX_VOLR and that
    ! both have the same io descriptor

    call getfil(frivinp, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ier = pio_inq_varid(ncid, name='SLOPE', vardesc=vardesc)
    if (ier /= PIO_noerr) then
       if (masterproc) write(iulog,*) subname//' variable SLOPE is not on dataset'
       readvar = .false.
    else
       readvar = .true.
    end if
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    if (readvar) then
       ier = pio_inq_vardimid(ncid, vardesc, dids)
       allocate(compdof(rtmCTL%lnumr))
       cnt = 0
       do n = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          compDOF(cnt) = rtmCTL%gindex(n)
       enddo
       call pio_initdecomp(pio_subsystem, pio_double, dids, compDOF, iodesc)
       deallocate(compdof)
! tcraig, there ia bug here, shouldn't use same vardesc for two different variable
       call pio_read_darray(ncid, vardesc, iodesc, rslope, ier)
       call pio_read_darray(ncid, vardesc, iodesc, max_volr, ier)
       call pio_freedecomp(ncid, iodesc)
    else
       rslope(:)   = 1._r8
       max_volr(:) = spval
    end if
    call pio_closefile(ncid)

    do nt = 1,nt_rtm
       do n = rtmCTL%begr, rtmCTL%endr
          fthresh(n) = 0.95*max_volr(n)*max(1._r8,rslope(n))
          ! modify velocity based on gridcell average slope (Manning eqn)
          evel(n,nt) = max(min_ev(nt),effvel(nt_rtm)*sqrt(max(0._r8,rslope(n)))) 
       end do
    end do

    deallocate(rslope, max_volr)

  end subroutine RtmFloodInit 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
  subroutine MOSART_init
!
! !REVISION HISTORY:
! Author: Hongyi Li

! !DESCRIPTION:
! initialize MOSART variables
! 
! !USES:
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
  type(file_desc_t)  :: ncid       ! pio file desc
  type(var_desc_t)   :: vardesc    ! pio variable desc 
  type(io_desc_t)    :: iodesc_dbl ! pio io desc
  type(io_desc_t)    :: iodesc_int ! pio io desc
  integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
  integer :: dids(2)               ! variable dimension ids 
  integer :: dsizes(2)             ! variable dimension lengths
  integer :: ier                  ! error code
  integer :: begr, endr, iunit, nn, n, cnt, nr, nt
  integer :: numDT_r, numDT_t
  integer :: lsize, gsize
  integer :: igrow, igcol, iwgt
  type(mct_avect) :: avtmp, avtmpG ! temporary avects
  type(mct_sMat)  :: sMat          ! temporary sparse matrix, needed for sMatP
  real(r8):: areatot_prev, areatot_tmp, areatot_new
  integer :: tcnt
  character(len=16384) :: rList             ! list of fields for SM multiply
  character(len=1000) :: fname
  character(len=*),parameter :: subname = '(MOSART_init)'
  character(len=*),parameter :: FORMI = '(2A,2i10)'
  character(len=*),parameter :: FORMR = '(2A,2g15.7)'
 
  begr = rtmCTL%begr
  endr = rtmCTL%endr
  
  if(endr >= begr) then
     ! routing parameters
     call ncd_pio_openfile (ncid, trim(frivinp_rtm), 0)
     call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
     allocate(compdof(rtmCTL%lnumr))
     cnt = 0
     do n = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        compDOF(cnt) = rtmCTL%gindex(n)
     enddo

     ! setup iodesc based on frac dids
     ier = pio_inq_varid(ncid, name='frac', vardesc=vardesc)
     ier = pio_inq_vardimid(ncid, vardesc, dids)
     ier = pio_inq_dimlen(ncid, dids(1),dsizes(1))
     ier = pio_inq_dimlen(ncid, dids(2),dsizes(2))
     call pio_initdecomp(pio_subsystem, pio_double, dsizes, compDOF, iodesc_dbl)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes, compDOF, iodesc_int)
     deallocate(compdof)
     call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

     allocate(TUnit%euler_calc(nt_rtm))
     Tunit%euler_calc = .true.

     allocate(TUnit%frac(begr:endr))
     ier = pio_inq_varid(ncid, name='frac', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%frac, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read frac ',minval(Tunit%frac),maxval(Tunit%frac)
     call shr_sys_flush(iulog)

     ! read fdir, convert to mask
     ! fdir <0 ocean, 0=outlet, >0 land
     ! tunit mask is 0=ocean, 1=land, 2=outlet for mosart calcs

     allocate(TUnit%mask(begr:endr))  
     ier = pio_inq_varid(ncid, name='fdir', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%mask, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read fdir mask ',minval(Tunit%mask),maxval(Tunit%mask)
     call shr_sys_flush(iulog)

     do n = rtmCtl%begr, rtmCTL%endr
        if (Tunit%mask(n) < 0) then
           Tunit%mask(n) = 0
        elseif (Tunit%mask(n) == 0) then
           Tunit%mask(n) = 2
           if (abs(Tunit%frac(n)-1.0_r8)>1.0e-9) then
              write(iulog,*) subname,' ERROR frac ne 1.0',n,Tunit%frac(n)
              call shr_sys_abort(subname//' ERROR frac ne 1.0')
           endif
        elseif (Tunit%mask(n) > 0) then
           Tunit%mask(n) = 1
           if (abs(Tunit%frac(n)-1.0_r8)>1.0e-9) then
              write(iulog,*) subname,' ERROR frac ne 1.0',n,Tunit%frac(n)
              call shr_sys_abort(subname//' ERROR frac ne 1.0')
           endif
        else
           call shr_sys_abort(subname//' Tunit mask error')
        endif
     enddo

     allocate(TUnit%ID0(begr:endr))  
     ier = pio_inq_varid(ncid, name='ID', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%ID0, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read ID0 ',minval(Tunit%ID0),maxval(Tunit%ID0)
     call shr_sys_flush(iulog)

     allocate(TUnit%dnID(begr:endr))  
     ier = pio_inq_varid(ncid, name='dnID', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%dnID, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read dnID ',minval(Tunit%dnID),maxval(Tunit%dnID)
     call shr_sys_flush(iulog)

     !-------------------------------------------------------
     ! RESET ID0 and dnID indices using the IDkey to be consistent
     ! with standard gindex order to leverage gsmap_r
     !-------------------------------------------------------
     do n=rtmCtl%begr, rtmCTL%endr
        TUnit%ID0(n)  = IDkey(TUnit%ID0(n))
        if (Tunit%dnID(n) > 0 .and. TUnit%dnID(n) <= rtmlon*rtmlat) then
           if (IDkey(TUnit%dnID(n)) > 0 .and. IDkey(TUnit%dnID(n)) <= rtmlon*rtmlat) then
              TUnit%dnID(n) = IDkey(TUnit%dnID(n))
           else
              write(iulog,*) subname,' ERROR bad IDkey for TUnit%dnID',n,TUnit%dnID(n),IDkey(TUnit%dnID(n))
              call shr_sys_abort(subname//' ERROR bad IDkey for TUnit%dnID')
           endif
        endif
     enddo

     allocate(TUnit%area(begr:endr))  
     ier = pio_inq_varid(ncid, name='area', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%area, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read area ',minval(Tunit%area),maxval(Tunit%area)
     call shr_sys_flush(iulog)

     do n=rtmCtl%begr, rtmCTL%endr
        if (TUnit%area(n) <= 0._r8) TUnit%area(n) = rtmCTL%area(n)
        if (TUnit%area(n) /= rtmCTL%area(n)) then
           write(iulog,*) subname,' ERROR area mismatch',TUnit%area(n),rtmCTL%area(n)
           call shr_sys_abort(subname//' ERROR area mismatch')
        endif
     enddo

     allocate(TUnit%areaTotal(begr:endr))  
     ier = pio_inq_varid(ncid, name='areaTotal', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%areaTotal, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read areaTotal ',minval(Tunit%areaTotal),maxval(Tunit%areaTotal)
     call shr_sys_flush(iulog)

     allocate(TUnit%rlenTotal(begr:endr))
     TUnit%rlenTotal = 0._r8

     allocate(TUnit%nh(begr:endr))  
     ier = pio_inq_varid(ncid, name='nh', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nh, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read nh ',minval(Tunit%nh),maxval(Tunit%nh)
     call shr_sys_flush(iulog)

     allocate(TUnit%hslp(begr:endr))  
     ier = pio_inq_varid(ncid, name='hslp', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%hslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read hslp ',minval(Tunit%hslp),maxval(Tunit%hslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%hslpsqrt(begr:endr))  
     TUnit%hslpsqrt = 0._r8

     allocate(TUnit%gxr(begr:endr))  
     ier = pio_inq_varid(ncid, name='gxr', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%gxr, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read gxr ',minval(Tunit%gxr),maxval(Tunit%gxr)
     call shr_sys_flush(iulog)

     allocate(TUnit%hlen(begr:endr))
     TUnit%hlen = 0._r8

     allocate(TUnit%tslp(begr:endr))  
     ier = pio_inq_varid(ncid, name='tslp', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%tslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read tslp ',minval(Tunit%tslp),maxval(Tunit%tslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%tslpsqrt(begr:endr))  
     TUnit%tslpsqrt = 0._r8

     allocate(TUnit%tlen(begr:endr))
     TUnit%tlen = 0._r8

     allocate(TUnit%twidth(begr:endr))  
     ier = pio_inq_varid(ncid, name='twid', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%twidth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read twidth ',minval(Tunit%twidth),maxval(Tunit%twidth)
     call shr_sys_flush(iulog)

     allocate(TUnit%nt(begr:endr))  
     ier = pio_inq_varid(ncid, name='nt', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nt, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read nt ',minval(Tunit%nt),maxval(Tunit%nt)
     call shr_sys_flush(iulog)

     allocate(TUnit%rlen(begr:endr))  
     ier = pio_inq_varid(ncid, name='rlen', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rlen, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rlen ',minval(Tunit%rlen),maxval(Tunit%rlen)
     call shr_sys_flush(iulog)

     allocate(TUnit%rslp(begr:endr))  
     ier = pio_inq_varid(ncid, name='rslp', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rslp ',minval(Tunit%rslp),maxval(Tunit%rslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%rslpsqrt(begr:endr))  
     TUnit%rslpsqrt = 0._r8

     allocate(TUnit%rwidth(begr:endr))  
     ier = pio_inq_varid(ncid, name='rwid', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rwidth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rwidth ',minval(Tunit%rwidth),maxval(Tunit%rwidth)
     call shr_sys_flush(iulog)

     allocate(TUnit%rwidth0(begr:endr))  
     ier = pio_inq_varid(ncid, name='rwid0', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rwidth0, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rwidth0 ',minval(Tunit%rwidth0),maxval(Tunit%rwidth0)
     call shr_sys_flush(iulog)

     allocate(TUnit%rdepth(begr:endr))  
     ier = pio_inq_varid(ncid, name='rdep', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rdepth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rdepth ',minval(Tunit%rdepth),maxval(Tunit%rdepth)
     call shr_sys_flush(iulog)

     allocate(TUnit%nr(begr:endr))  
     ier = pio_inq_varid(ncid, name='nr', vardesc=vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nr, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read nr ',minval(Tunit%nr),maxval(Tunit%nr)
     call shr_sys_flush(iulog)

     allocate(TUnit%nUp(begr:endr))
     TUnit%nUp = 0

     allocate(TUnit%iUp(begr:endr,8))
     TUnit%iUp = 0

     allocate(TUnit%indexDown(begr:endr))
     TUnit%indexDown = 0

     ! initialize water states and fluxes
     allocate (TRunoff%wh(begr:endr,nt_rtm))
     TRunoff%wh = 0._r8

     allocate (TRunoff%dwh(begr:endr,nt_rtm))
     TRunoff%dwh = 0._r8

     allocate (TRunoff%yh(begr:endr,nt_rtm))
     TRunoff%yh = 0._r8

     allocate (TRunoff%qsur(begr:endr,nt_rtm))
     TRunoff%qsur = 0._r8

     allocate (TRunoff%qsub(begr:endr,nt_rtm))
     TRunoff%qsub = 0._r8

     allocate (TRunoff%qgwl(begr:endr,nt_rtm))
     TRunoff%qgwl = 0._r8

     allocate (TRunoff%ehout(begr:endr,nt_rtm))
     TRunoff%ehout = 0._r8

     allocate (TRunoff%tarea(begr:endr,nt_rtm))
     TRunoff%tarea = 0._r8

     allocate (TRunoff%wt(begr:endr,nt_rtm))
     TRunoff%wt= 0._r8

     allocate (TRunoff%dwt(begr:endr,nt_rtm))
     TRunoff%dwt = 0._r8

     allocate (TRunoff%yt(begr:endr,nt_rtm))
     TRunoff%yt = 0._r8

     allocate (TRunoff%mt(begr:endr,nt_rtm))
     TRunoff%mt = 0._r8

     allocate (TRunoff%rt(begr:endr,nt_rtm))
     TRunoff%rt = 0._r8

     allocate (TRunoff%pt(begr:endr,nt_rtm))
     TRunoff%pt = 0._r8

     allocate (TRunoff%vt(begr:endr,nt_rtm))
     TRunoff%vt = 0._r8

     allocate (TRunoff%tt(begr:endr,nt_rtm))
     TRunoff%tt = 0._r8

     allocate (TRunoff%etin(begr:endr,nt_rtm))
     TRunoff%etin = 0._r8

     allocate (TRunoff%etout(begr:endr,nt_rtm))
     TRunoff%etout = 0._r8

     allocate (TRunoff%rarea(begr:endr,nt_rtm))
     TRunoff%rarea = 0._r8

     allocate (TRunoff%wr(begr:endr,nt_rtm))
     TRunoff%wr = 0._r8

     allocate (TRunoff%dwr(begr:endr,nt_rtm))
     TRunoff%dwr = 0._r8

     allocate (TRunoff%yr(begr:endr,nt_rtm))
     TRunoff%yr = 0._r8

     allocate (TRunoff%mr(begr:endr,nt_rtm))
     TRunoff%mr = 0._r8

     allocate (TRunoff%rr(begr:endr,nt_rtm))
     TRunoff%rr = 0._r8

     allocate (TRunoff%pr(begr:endr,nt_rtm))
     TRunoff%pr = 0._r8

     allocate (TRunoff%vr(begr:endr,nt_rtm))
     TRunoff%vr = 0._r8

     allocate (TRunoff%tr(begr:endr,nt_rtm))
     TRunoff%tr = 0._r8

     allocate (TRunoff%erlg(begr:endr,nt_rtm))
     TRunoff%erlg = 0._r8

     allocate (TRunoff%erlateral(begr:endr,nt_rtm))
     TRunoff%erlateral = 0._r8

     allocate (TRunoff%erin(begr:endr,nt_rtm))
     TRunoff%erin = 0._r8

     allocate (TRunoff%erout(begr:endr,nt_rtm))
     TRunoff%erout = 0._r8

     allocate (TRunoff%eroup_lagi(begr:endr,nt_rtm))
     TRunoff%eroup_lagi = 0._r8

     allocate (TRunoff%eroup_lagf(begr:endr,nt_rtm))
     TRunoff%eroup_lagf = 0._r8

     allocate (TRunoff%erowm_regi(begr:endr,nt_rtm))
     TRunoff%erowm_regi = 0._r8

     allocate (TRunoff%erowm_regf(begr:endr,nt_rtm))
     TRunoff%erowm_regf = 0._r8

     allocate (TRunoff%eroutUp(begr:endr,nt_rtm))
     TRunoff%eroutUp = 0._r8

     allocate (TRunoff%eroutUp_avg(begr:endr,nt_rtm))
     TRunoff%eroutUp_avg = 0._r8

     allocate (TRunoff%erlat_avg(begr:endr,nt_rtm))
     TRunoff%erlat_avg = 0._r8

     allocate (TRunoff%ergwl(begr:endr,nt_rtm))
     TRunoff%ergwl = 0._r8

     allocate (TRunoff%flow(begr:endr,nt_rtm))
     TRunoff%flow = 0._r8
    
     allocate (TPara%c_twid(begr:endr))
     TPara%c_twid = 1.0_r8

     call pio_freedecomp(ncid, iodesc_dbl)
     call pio_freedecomp(ncid, iodesc_int)
     call pio_closefile(ncid)

   ! control parameters and some other derived parameters
   ! estimate derived input variables

     do iunit=rtmCTL%begr,rtmCTL%endr
        if(TUnit%Gxr(iunit) > 0._r8) then
           TUnit%rlenTotal(iunit) = TUnit%area(iunit)*TUnit%Gxr(iunit)
        end if
     end do

     do iunit=rtmCTL%begr,rtmCTL%endr
        if(TUnit%rlen(iunit) > TUnit%rlenTotal(iunit)) then
           TUnit%rlenTotal(iunit) = TUnit%rlen(iunit)
        end if
     end do     

     do iunit=rtmCTL%begr,rtmCTL%endr
      
        if(TUnit%rlen(iunit) > 0._r8) then
           TUnit%hlen(iunit) = TUnit%area(iunit) / TUnit%rlenTotal(iunit) / 2._r8
           if(TUnit%hlen(iunit) > 50000_r8) then
              TUnit%hlen(iunit) = 50000_r8   ! allievate the outlier in drainage density estimation. TO DO
           end if
           TUnit%tlen(iunit) = TUnit%area(iunit) / TUnit%rlen(iunit) / 2._r8 - TUnit%hlen(iunit)
           if(TUnit%twidth(iunit) < 0._r8) then
              TUnit%twidth(iunit) = 0._r8
           end if
           if(TUnit%tlen(iunit) > 0._r8 .and. (TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit) > 1._r8) then
              TUnit%twidth(iunit) = TPara%c_twid(iunit)*TUnit%twidth(iunit)*((TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit))
           end if
          
           if(TUnit%tlen(iunit) > 0._r8 .and. TUnit%twidth(iunit) <= 0._r8) then
              TUnit%twidth(iunit) = 0._r8
           end if
        else
           TUnit%hlen(iunit) = 0._r8
           TUnit%tlen(iunit) = 0._r8
           TUnit%twidth(iunit) = 0._r8
        end if
        
        if(TUnit%rslp(iunit) <= 0._r8) then
           TUnit%rslp(iunit) = 0.0001_r8
        end if
        if(TUnit%tslp(iunit) <= 0._r8) then
           TUnit%tslp(iunit) = 0.0001_r8
        end if
        if(TUnit%hslp(iunit) <= 0._r8) then
           TUnit%hslp(iunit) = 0.005_r8
        end if
        TUnit%rslpsqrt(iunit) = sqrt(Tunit%rslp(iunit))
        TUnit%tslpsqrt(iunit) = sqrt(Tunit%tslp(iunit))
        TUnit%hslpsqrt(iunit) = sqrt(Tunit%hslp(iunit))
     end do

     lsize = rtmCTL%lnumr
     gsize = rtmlon*rtmlat

     if (smat_option == 'opt') then
        ! distributed smat initialization
        ! mct_sMat_init must be given the number of rows and columns that
        ! would be in the full matrix.  Nrows= size of output vector=nb.
        ! Ncols = size of input vector = na.

        cnt = 0
        do iunit=rtmCTL%begr,rtmCTL%endr
           if(TUnit%dnID(iunit) > 0) cnt = cnt + 1
        enddo

        call mct_sMat_init(sMat, gsize, gsize, cnt)
        igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
        igrow = mct_sMat_indexIA(sMat,'grow')
        iwgt  = mct_sMat_indexRA(sMat,'weight')
        cnt = 0
        do iunit = rtmCTL%begr,rtmCTL%endr
           if (TUnit%dnID(iunit) > 0) then
              cnt = cnt + 1
              sMat%data%iAttr(igcol,cnt) = TUnit%ID0(iunit)
              sMat%data%iAttr(igrow,cnt) = TUnit%dnID(iunit)
              sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
           endif
        enddo

        call mct_sMatP_Init(sMatP_eroutUp, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

     elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then
        ! root initialization
        call mct_aVect_init(avtmp,rList='f1:f2',lsize=lsize)
        call mct_aVect_zero(avtmp)
        cnt = 0
        do iunit = rtmCTL%begr,rtmCTL%endr
           cnt = cnt + 1
           avtmp%rAttr(1,cnt) = TUnit%ID0(iunit)
           avtmp%rAttr(2,cnt) = TUnit%dnID(iunit)
        enddo
        call mct_avect_gather(avtmp,avtmpG,gsmap_r,mastertask,mpicom_rof)
        if (masterproc) then
           cnt = 0
           do n = 1,rtmlon*rtmlat
              if (avtmpG%rAttr(2,n) > 0) then
                 cnt = cnt + 1
              endif
           enddo

           call mct_sMat_init(sMat, gsize, gsize, cnt)
           igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
           igrow = mct_sMat_indexIA(sMat,'grow')
           iwgt  = mct_sMat_indexRA(sMat,'weight')

           cnt = 0
           do n = 1,rtmlon*rtmlat
              if (avtmpG%rAttr(2,n) > 0) then
                 cnt = cnt + 1
                 sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                 sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(2,n)
                 sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
              endif
           enddo
           call mct_avect_clean(avtmpG)
        else
           call mct_sMat_init(sMat,1,1,1)
        endif
        call mct_avect_clean(avtmp)

        call mct_sMatP_Init(sMatP_eroutUp, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

     else

        write(iulog,*) trim(subname),' MOSART ERROR: invalid smat_option '//trim(smat_option)
        call shr_sys_abort(trim(subname)//' ERROR invald smat option')

     endif

     ! initialize the AVs to go with sMatP
     write(rList,'(a,i3.3)') 'tr',1
     do nt = 2,nt_rtm
        write(rList,'(a,i3.3)') trim(rList)//':tr',nt
     enddo
!    write(iulog,*) trim(subname),' MOSART initialize avect ',trim(rList)
     call mct_aVect_init(avsrc_eroutUp,rList=rList,lsize=rtmCTL%lnumr)
     call mct_aVect_init(avdst_eroutUp,rList=rList,lsize=rtmCTL%lnumr)

     lsize = mct_smat_gNumEl(sMatP_eroutUp%Matrix,mpicom_rof)
     if (masterproc) write(iulog,*) subname," Done initializing SmatP_eroutUp, nElements = ",lsize

     ! keep only sMatP
     call mct_sMat_clean(sMat)

  end if  ! endr >= begr

  !--- compute areatot from area using dnID ---
  !--- this basically advects upstream areas downstream and
  !--- adds them up as it goes until all upstream areas are accounted for

  allocate(Tunit%areatotal2(rtmCTL%begr:rtmCTL%endr))
  Tunit%areatotal2 = 0._r8

  ! initialize avdst to local area and add that to areatotal2
  cnt = 0
  call mct_avect_zero(avdst_eroutUp)
  do nr = rtmCTL%begr,rtmCTL%endr
     cnt = cnt + 1
     avdst_eroutUp%rAttr(1,cnt) = rtmCTL%area(nr)
     Tunit%areatotal2(nr) = avdst_eroutUp%rAttr(1,cnt)
  enddo

  tcnt = 0
  areatot_prev = -99._r8
  areatot_new = -50._r8
  do while (areatot_new /= areatot_prev .and. tcnt < rtmlon*rtmlat)

     tcnt = tcnt + 1

     ! copy avdst to avsrc for next downstream step
     cnt = 0
     call mct_avect_zero(avsrc_eroutUp)
     do nr = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        avsrc_eroutUp%rAttr(1,cnt) = avdst_eroutUp%rAttr(1,cnt)
     enddo

     call mct_avect_zero(avdst_eroutUp)

     call mct_sMat_avMult(avsrc_eroutUp, sMatP_eroutUp, avdst_eroutUp)

     ! add avdst to areatot and compute new global sum
     cnt = 0
     areatot_prev = areatot_new
     areatot_tmp = 0._r8
     do nr = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        Tunit%areatotal2(nr) = Tunit%areatotal2(nr) + avdst_eroutUp%rAttr(1,cnt)
        areatot_tmp = areatot_tmp + Tunit%areatotal2(nr)
     enddo
     call shr_mpi_sum(areatot_tmp, areatot_new, mpicom_rof, 'areatot_new', all=.true.)

     !if (masterproc) then
     !   write(iulog,*) trim(subname),' areatot calc ',tcnt,areatot_new
     !endif

  enddo

  if (areatot_new /= areatot_prev) then
     write(iulog,*) trim(subname),' MOSART ERROR: areatot incorrect ',areatot_new, areatot_prev
     call shr_sys_abort(trim(subname)//' ERROR areatot incorrect')
  endif

!  do nr = rtmCTL%begr,rtmCTL%endr
!     if (TUnit%areatotal(nr) > 0._r8 .and. Tunit%areatotal2(nr) /= TUnit%areatotal(nr)) then
!        write(iulog,'(2a,i12,2e16.4,f16.4)') trim(subname),' areatot diff ',nr,TUnit%areatotal(nr),Tunit%areatota!l2(nr),&
!           abs(TUnit%areatotal(nr)-Tunit%areatotal2(nr))/(TUnit%areatotal(nr))
!     endif
!  enddo


  ! control parameters
  Tctl%RoutingMethod = 1
  !Tctl%DATAH = rtm_nsteps*get_step_size()
  !Tctl%DeltaT = 60._r8  !
  !   if(Tctl%DATAH > 0 .and. Tctl%DATAH < Tctl%DeltaT) then
  !       Tctl%DeltaT = Tctl%DATAH
  !   end if      
  Tctl%DLevelH2R = 5
  Tctl%DLevelR = 3
  call SubTimestep ! prepare for numerical computation

  call shr_mpi_max(maxval(Tunit%numDT_r),numDT_r,mpicom_rof,'numDT_r',all=.false.)
  call shr_mpi_max(maxval(Tunit%numDT_t),numDT_t,mpicom_rof,'numDT_t',all=.false.)
  if (masterproc) then
     write(iulog,*) subname,' DLevelH2R = ',Tctl%DlevelH2R
     write(iulog,*) subname,' numDT_r   = ',minval(Tunit%numDT_r),maxval(Tunit%numDT_r)
     write(iulog,*) subname,' numDT_r max  = ',numDT_r
     write(iulog,*) subname,' numDT_t   = ',minval(Tunit%numDT_t),maxval(Tunit%numDT_t)
     write(iulog,*) subname,' numDT_t max  = ',numDT_t
  endif
    
  !if(masterproc) then
  !    fname = '/lustre/liho745/DCLM_model/ccsm_hy/run/clm_MOSART_subw2/run/test.dat'
  !    call createFile(1111,fname)
  !end if

  end subroutine MOSART_init

!----------------------------------------------------------------------------

  subroutine SubTimestep
  ! !DESCRIPTION: predescribe the sub-time-steps for channel routing
    implicit none    
    integer :: iunit   !local index
    character(len=*),parameter :: subname = '(SubTimestep)'

    allocate(TUnit%numDT_r(rtmCTL%begr:rtmCTL%endr),TUnit%numDT_t(rtmCTL%begr:rtmCTL%endr))
    TUnit%numDT_r = 1
    TUnit%numDT_t = 1
    allocate(TUnit%phi_r(rtmCTL%begr:rtmCTL%endr),TUnit%phi_t(rtmCTL%begr:rtmCTL%endr))
    TUnit%phi_r = 0._r8
    TUnit%phi_t = 0._r8

    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0 .and. TUnit%rlen(iunit) > 0._r8) then
          TUnit%phi_r(iunit) = TUnit%areaTotal2(iunit)*sqrt(TUnit%rslp(iunit))/(TUnit%rlen(iunit)*TUnit%rwidth(iunit))
          if(TUnit%phi_r(iunit) >= 10._r8) then
             TUnit%numDT_r(iunit) = (TUnit%numDT_r(iunit)*log10(TUnit%phi_r(iunit))*Tctl%DLevelR) + 1
          else 
             TUnit%numDT_r(iunit) = TUnit%numDT_r(iunit)*1.0_r8*Tctl%DLevelR + 1
          end if
       end if
       if(TUnit%numDT_r(iunit) < 1) TUnit%numDT_r(iunit) = 1
      
       if(TUnit%tlen(iunit) > 0._r8) then
          TUnit%phi_t(iunit) =      TUnit%area(iunit)*sqrt(TUnit%tslp(iunit))/(TUnit%tlen(iunit)*TUnit%twidth(iunit))
          if(TUnit%phi_t(iunit) >= 10._r8) then 
             TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*log10(TUnit%phi_t(iunit))*Tctl%DLevelR) + 1
          else 
             TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*1.0*Tctl%DLevelR) + 1
          end if
       end if
       if(TUnit%numDT_t(iunit) < 1) TUnit%numDT_t(iunit) = 1
    end do
  end subroutine SubTimestep

!-----------------------------------------------------------------------

end module RtmMod

