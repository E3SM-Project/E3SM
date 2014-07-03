module RtmMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! River Routing Model (U. of Texas River Transport
! Model)~\cite{Branstetter:2001}
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_PI, SHR_CONST_CDAY
  use rtm_cpl_indices , only : nt_rtm, rtm_tracers 
  use RtmSpmd         , only : masterproc, npes, iam, mpicom_rof, &
                               MPI_REAL8,MPI_INTEGER,MPI_CHARACTER,MPI_LOGICAL,MPI_MAX
  use RtmVar          , only : re, spval, rtmlon, rtmlat, iulog, ice_runoff, &
                               frivinp_rtm, finidat_rtm, nrevsn_rtm, &
                               nsrContinue, nsrBranch, nsrStartup, nsrest, &
                               inst_index, inst_suffix, inst_name, &
                               rtm_active, flood_active, effvel_active
  use RtmFileUtils    , only : getfil, getavu, relavu
  use RtmTimeManager  , only : timemgr_init, get_nstep, get_curr_date
  use RtmHistFlds     , only : RtmHistFldsInit, RtmHistFldsSet 
  use RtmHistFile     , only : RtmHistUpdateHbuf, RtmHistHtapesWrapup, RtmHistHtapesBuild, &
                               rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq,     &
                               rtmhist_avgflag_pertape, rtmhist_avgflag_pertape, & 
                               rtmhist_fincl1, rtmhist_fincl2, rtmhist_fincl3,   &
                               rtmhist_fexcl1, rtmhist_fexcl2, rtmhist_fexcl3,   &
                               max_tapes, max_namlen
  use RtmRestFile     , only : RtmRestTimeManager, RtmRestGetFile, RtmRestFileRead, &
                               RtmRestFileWrite, RtmRestFileName
  use RunoffMod       , only : RunoffInit, runoff
  use RtmIO
  use mct_mod
  use perf_mod
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini          ! Initialize RTM grid
  public Rtmrun          ! River routing model (based on U. Texas code)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: RtmFloodInit

! !PRIVATE TYPES:

! RTM tracers
  character(len=256) :: rtm_trstr   ! tracer string

! RTM naemlists
  integer :: rtm_tstep                    ! RTM time step

! RTM constants
  real(r8),save :: delt_rtm_max              ! RTM max timestep
  real(r8) :: cfl_scale = 0.1_r8        ! cfl scale factor, must be <= 1.0

!global (glo)
  integer , pointer :: dwnstrm_index(:)! downstream index

!local (gdc)
  real(r8), save, pointer :: ddist(:)        ! downstream dist (m)
  real(r8), save, pointer :: evel(:,:)       ! effective tracer velocity (m/s)
  real(r8), save, pointer :: sfluxin(:,:)    ! cell tracer influx (m3/s)
  real(r8), save, pointer :: fluxout(:,:)    ! cell tracer outlflux (m3/s)

! global rtm grid
  real(r8),pointer :: rlatc(:)    ! latitude of 1d grid cell (deg)
  real(r8),pointer :: rlonc(:)    ! longitude of 1d grid cell (deg)
  real(r8),pointer :: rlats(:)    ! latitude of 1d south grid cell edge (deg)
  real(r8),pointer :: rlatn(:)    ! latitude of 1d north grid cell edge (deg)
  real(r8),pointer :: rlonw(:)    ! longitude of 1d west grid cell edge (deg)
  real(r8),pointer :: rlone(:)    ! longitude of 1d east grid cell edge (deg)

  character(len=256) :: flood_mode
  character(len=256) :: rtm_mode
  character(len=256) :: rtm_effvel    

  character(len=256) :: nlfilename_rof = 'rof_in' 
  character(len=256) :: nlfilename_lnd = 'lnd_in' 
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
  subroutine Rtmini()
!
! !DESCRIPTION:
! Initialize RTM grid, mask, decomp
!
! !USES:
!
! !ARGUMENTS:
    implicit none
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


    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !rdirc input to i
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !rdirc input to j
    real(r8) :: edgen =   90._r8
    real(r8) :: edgee =  180._r8
    real(r8) :: edges =  -90._r8
    real(r8) :: edgew = -180._r8
    integer  :: i,j,k,n,ng,g,n2,nt            ! loop indices
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8) :: lrtmarea                      ! tmp local sum of area
    real(r8),allocatable :: tempr(:,:)        ! temporary buffer
    integer ,allocatable :: rdirc(:)          ! temporary buffer
    integer ,allocatable :: iocn(:)           ! downstream ocean cell
    integer ,allocatable :: nocn(:)           ! number of rtm cells in basin
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of rtm cells on a pe
    integer ,allocatable :: nba(:)            ! number of basins on each pe
    integer ,allocatable :: nrs(:)            ! begr on each pe
    integer ,allocatable :: basin(:)          ! basin to rtm mapping
    integer  :: nbas                          ! number of basins
    integer  :: nrtm                          ! num of rtm points
    integer  :: baspe                         ! pe with min number of rtm cells
    integer  :: maxrtm                        ! max num of rtms per pe for decomp
    integer  :: minbas,maxbas                 ! used for decomp search
    integer  :: nl,nloops                     ! used for decomp search
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: begr,endr,numr                ! tot num of roff pts on all pes
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
    logical           :: found                ! if variable found on rdirc file
    character(len=256):: fnamer               ! name of netcdf restart file 
    character(len=256):: pnamer               ! full pathname of netcdf restart file
    character(len=256):: locfn                ! local file name
    integer , pointer :: num_rtm(:)           ! num of cells on each pe
    integer           :: begro,endro          ! local start/stop indices
    integer           :: begrl,endrl          ! local start/stop indices
    integer           :: lnumrl               ! rtm gdc local number of lnd cells
    integer           :: lnumro               ! rtm gdc local number of ocn cells
    integer           :: unitn                ! unit for namelist file
    integer,parameter :: dbug = 1             ! 0 = none, 1=normal, 2=much, 3=max
    logical :: lexist                         ! File exists
    character(len= 7) :: runtyp(4)            ! run type
    character(*),parameter :: subname = '(Rtmini) '
!-----------------------------------------------------------------------

    !-------------------------------------------------------
    ! Read in rtm namelist
    !-------------------------------------------------------

    namelist /rtm_inparm / &
         ice_runoff, rtm_mode, flood_mode, rtm_effvel, &
         frivinp_rtm, finidat_rtm, nrevsn_rtm, rtm_tstep, &
         rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq, &
         rtmhist_fincl1,  rtmhist_fincl2, rtmhist_fincl3, &
         rtmhist_fexcl1,  rtmhist_fexcl2, rtmhist_fexcl3, &
         rtmhist_avgflag_pertape

    ! Preset values
    rtm_effvel     = 'NULL'
    rtm_mode    = 'ACTIVE'
    flood_mode  = 'NULL'
    ice_runoff  = .true.
    finidat_rtm = ' '
    nrevsn_rtm  = ' '
    rtm_tstep   = -1

    nlfilename_rof = "rof_in" // trim(inst_suffix)
    inquire (file = trim(nlfilename_rof), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_rof does NOT exist:'&
            //trim(nlfilename_rof)
       call shr_sys_abort()
    end if
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in rtm_inparm namelist from: ', trim(nlfilename_rof)
       open( unitn, file=trim(nlfilename_rof), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, rtm_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on rtm_inparm read' )
          endif
       end do
       call relavu( unitn )
    end if

    call mpi_bcast (rtm_tstep,   1, MPI_INTEGER, 0, mpicom_rof, ier)

    call mpi_bcast (finidat_rtm,  len(finidat_rtm), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (frivinp_rtm,  len(frivinp_rtm), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (nrevsn_rtm ,  len(nrevsn_rtm) , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtm_mode,     len(rtm_mode)   , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (flood_mode,   len(flood_mode) , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtm_effvel,   len(rtm_effvel) , MPI_CHARACTER, 0, mpicom_rof, ier)

    call mpi_bcast (ice_runoff,  1, MPI_LOGICAL, 0, mpicom_rof, ier)

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
       if (nsrest == nsrStartup .and. finidat_rtm /= ' ') then
          write(iulog,*) '   RTM initial data   = ',trim(finidat_rtm)
       end if
    endif

    rtm_active = .true.
    flood_active = .true.
    effvel_active = .true.

    if (trim(rtm_mode) == 'NULL') then
       rtm_active = .false.
       flood_active = .false.
       effvel_active = .true.
    endif

    if (trim(flood_mode) == 'NULL') then
       flood_active = .false.
    endif

    if (trim(rtm_effvel) == 'NULL') then
       effvel_active = .false.
    endif

    if (masterproc) then
      if (flood_active) then
         write(iulog,*) '   RTM :: flooding is on '
      else
         write(iulog,*) '   RTM :: flooding is off '
      endif
      if (effvel_active) then
         write(iulog,*) '   RTM :: calculate effective velocity with slope (4.5) '
      else
         write(iulog,*) '   RTM :: use default effective velocity (4.0) '
      endif
    endif
    
    if (rtm_active) then
       if (frivinp_rtm == ' ') then
          call shr_sys_abort( subname//' ERROR: rtm_mode ACTIVE, but frivinp_rtm NOT set' )
       else
          if (masterproc) then
             write(iulog,*) '   RTM river data       = ',trim(frivinp_rtm)
          endif
       end if
    else
       if (masterproc) then
          write(iulog,*)'RTM will not be active '
       endif
       RETURN
    end if

    if (rtm_tstep <= 0) then
       write(iulog,*) subname,': ERROR rtm step invalid',rtm_tstep
       call shr_sys_abort( subname//' ERROR: rtm_tstep invalid' )
    endif
       
    do i = 1, max_tapes
       if (rtmhist_nhtfrq(i) == 0) then
          rtmhist_mfilt(i) = 1
       else if (rtmhist_nhtfrq(i) < 0) then
          rtmhist_nhtfrq(i) = nint(-rtmhist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*rtm_tstep))
       endif
    end do

    !-------------------------------------------------------
    ! Initialize rtm time manager 
    !-------------------------------------------------------

    ! Intiialize RTM pio
    call ncd_pio_init()

    ! Obtain restart file if appropriate
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
 	(nsrest == nsrBranch  )) then
       call RtmRestGetfile( file=fnamer, path=pnamer )
    endif       

    ! Initialize time manager
    if (nsrest == nsrStartup) then  
       call timemgr_init(dtime_in=rtm_tstep)
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
       write(iulog,*)'rtm tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !-------------------------------------------------------
    ! Read input data (river direction file)
    !-------------------------------------------------------

    ! Useful constants and initial values
    deg2rad = SHR_CONST_PI / 180._r8

    call t_startf('rtmi_grid')

    call getfil(frivinp_rtm, locfn, 0 )
    if (masterproc) then
       write(iulog,*)'Read in RTM file name: ',trim(frivinp_rtm)
       call shr_sys_flush(iulog)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqdid(ncid,'ni',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlon)
    call ncd_inqdid(ncid,'nj',dimid)
    call ncd_inqdlen(ncid,dimid,rtmlat)

    if (masterproc) then
       write(iulog,*) 'Values for rtmlon/rtmlat: ',rtmlon,rtmlat
       write(iulog,*) 'Successfully read RTM dimensions'
       call shr_sys_flush(iulog)
    endif

    ! Allocate variables
    allocate(rlonc(rtmlon), rlatc(rtmlat), &
             rlonw(rtmlon), rlone(rtmlon), &
             rlats(rtmlat), rlatn(rtmlat), &
             runoff%rlon(rtmlon),          &
             runoff%rlat(rtmlat),          &
             rdirc(rtmlon*rtmlat),         &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation ERROR for rdirc'
       call shr_sys_abort()
    end if

    allocate(tempr(rtmlon,rtmlat))  
    call ncd_io(ncid=ncid, varname='RTM_FLOW_DIRECTION', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: RTM_FLOW_DIRECTION NOT on rdirc file' )
    do j=1,rtmlat
      do i=1,rtmlon
         
         !-------------------------------------------------------
         ! Put in a check for a negative rdirc value and abort.
         !-------------------------------------------------------
         if ( tempr(i,j) .lt. 0 ) call shr_sys_abort( trim(subname)//' ERROR: Found a negative RTM_FLOW_DIRECTION. &
              This is currently not supported. ' )
         n = (j-1)*rtmlon + i
         rdirc(n) = nint(tempr(i,j))
      enddo
    enddo
    call ncd_io(ncid=ncid, varname='xc', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: RTM longitudes NOT on rdirc file' )
    do i=1,rtmlon
       runoff%rlon(i) = tempr(i,1)
       rlonc(i) = tempr(i,1)
    enddo
    call ncd_io(ncid=ncid, varname='yc', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: RTM latitudes NOT on rdirc file' )
    do j=1,rtmlat
       runoff%rlat(j) = tempr(1,j)
       rlatc(j) = tempr(1,j)
    end do
    deallocate(tempr)             

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*)'RTM netcdf river direction file successfully read '
       call shr_sys_flush(iulog)
    endif

    !-------------------------------------------------------
    ! Set dwnstrm_index from rdirc values
    !-------------------------------------------------------

    ! The following assumes that there is no runoff  
    ! south of j=1 or north of j=rtmlat
    ! This is true for rdirc.05
    ! Determine dwnstrmm_index from rtm river flow direction (0-8)

    allocate (dwnstrm_index(rtmlon*rtmlat), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation ERROR for dwnstrm_index'
       call shr_sys_abort()
    end if

    dwnstrm_index(:) = 0
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       if (rdirc(n) /= 0) then
          ir = i + ioff(rdirc(n))
          jr = j + joff(rdirc(n))
          if (ir < 1     ) ir = ir + rtmlon
          if (ir > rtmlon) ir = ir - rtmlon
          !--- check cross pole flow, etc
          if (jr < 1 .or. jr > rtmlat .or. ir < 1 .or. ir > rtmlon) then
             write(iulog,*) 'Rtmini ERROR ir jr bounds ',i,j,rdirc(n),ir,jr
             call shr_sys_abort()
          endif
          nr = (jr-1)*rtmlon + ir
          if (n == nr) then
             write(iulog,*) 'Rtmini ERROR dwnstrm_index ',i,j,n,rdirc(n),ir,jr,nr
             call shr_sys_abort()
          endif
          dwnstrm_index(n) = nr
       endif
    enddo
    enddo

    !-------------------------------------------------------
    ! Determine rtm ocn/land mask (global, all procs)
    !-------------------------------------------------------

    !  0=none, 
    !  1=land, 
    !  2=ocean outflow, 
    ! -1=reroute over ocean to ocean outflow points

    call t_startf('rtmi_decomp')

    allocate (gmask(rtmlon*rtmlat), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation ERROR for gmask'
       call shr_sys_abort()
    end if

    ! Initialize gmask to 2 everywhere, then compute land cells

    gmask(:) = 2

    do n=1,rtmlon*rtmlat         ! override downstream setting from local info
       nr = dwnstrm_index(n)
       if (nr /= 0) then         ! n is always land if dwnstrm_index exists
          if (rdirc(n) > 0) then
             gmask(n) = 1
          else if (rdirc(n) < 0) then 
             gmask(n) = -1
          end if
       end if
    enddo

    deallocate(rdirc)

    !-------------------------------------------------------
    ! Compute river basins, actually compute ocean outlet gridcell
    !-------------------------------------------------------

    ! iocn = final downstream cell, index is global 1d ocean gridcell
    ! nocn = number of source gridcells for ocean gridcell

    allocate(iocn(rtmlon*rtmlat),nocn(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation ERROR for ',&
            'iocn,nocn'
       call shr_sys_abort()
    end if

    call t_startf('rtmi_dec_basins')
    iocn = 0
    nocn = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             n = dwnstrm_index(n)
             g = g + 1
          end do
          if (gmask(n) == 2) then           ! found ocean outlet 
             iocn(nr) = n                   ! set ocean outlet or nr to n
             nocn(n) = nocn(n) + 1          ! one more land cell for n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) 'rtmini ERROR no downstream ocean cell ', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
             call shr_sys_abort()
          else 
             write(iulog,*) 'rtmini ERROR downstream cell is non-ocean,non-land', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
             call shr_sys_abort()
          endif
       elseif (gmask(n) == 2) then  ! ocean, give to self
          iocn(nr) = n
          nocn(n) = nocn(n) + 1
       endif
    enddo
    call t_stopf('rtmi_dec_basins')

    !-------------------------------------------------------
    !--- Now allocate those basins to pes
    !-------------------------------------------------------

    call t_startf('rtmi_dec_distr')

    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of rtm cells/pe 

    nbas = 0
    nrtm = 0
    do nr=1,rtmlon*rtmlat
       if (nocn(nr) > 0) then
          nbas = nbas + 1
          nrtm = nrtm + nocn(nr)
       endif
    enddo

    allocate(pocn(rtmlon*rtmlat),     &  !global rtm array
             rglo2gdc(rtmlon*rtmlat), &  !global rtm array
             nop(0:npes-1), &
             nba(0:npes-1), &
             nrs(0:npes-1), &
             num_rtm(0:npes-1))

    nop = 0
    nba = 0
    nrs = 0
    pocn = -99
    rglo2gdc = 0
    baspe = 0
    maxrtm = int(float(nrtm)/float(npes)*0.445) + 1
    nloops = 3
    minbas = nrtm
    do nl=1,nloops
       maxbas = minbas - 1
       minbas = maxval(nocn)/(2**nl)
       if (nl == nloops) minbas = min(minbas,1)
       do nr=1,rtmlon*rtmlat
          if (nocn(nr) > 0 .and. nocn(nr) >= minbas .and. nocn(nr) <= maxbas) then
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
                write(iulog,*) 'ERROR in decomp for rtm ',nr,npes,baspe
                call shr_sys_abort()
             endif
             nop(baspe) = nop(baspe) + nocn(nr)
             nba(baspe) = nba(baspe) + 1
             pocn(nr) = baspe
          endif
       enddo ! nr
    enddo ! nl

    ! set pocn for land cells, was set for ocean above
    do nr=1,rtmlon*rtmlat
       if (iocn(nr) > 0) then
          pocn(nr) = pocn(iocn(nr))
          if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
             write(iulog,*) 'Rtmini ERROR pocn lnd setting ',&
                  nr,iocn(nr),iocn(iocn(nr)),pocn(iocn(nr)),pocn(nr),npes
             call shr_sys_abort()
          endif
       endif
    enddo

    if (masterproc) then
       write(iulog,*) 'rtm cells and basins total  = ',nrtm,nbas
       write(iulog,*) 'rtm cells per basin avg/max = ',nrtm/nbas,maxval(nocn)
       write(iulog,*) 'rtm cells per pe    min/max = ',minval(nop),maxval(nop)
       write(iulog,*) 'basins    per pe    min/max = ',minval(nba),maxval(nba)
    end if

    !-------------------------------------------------------
    !--- Count and distribute cells to rglo2gdc
    !-------------------------------------------------------

    runoff%numr   = 0
    runoff%numro  = 0
    runoff%numrl  = 0
    runoff%lnumr  = 0
    lnumro = 0
    lnumrl = 0
    num_rtm = 0

    do n = 0,npes-1
       if (iam == n) then
          runoff%begr  = runoff%numr  + 1
          begrl = runoff%numrl + 1
          begro = runoff%numro + 1
       endif
       num_rtm(n)   = num_rtm(n)   + nop(n)
       runoff%numr  = runoff%numr  + nop(n)
       runoff%numro = runoff%numro + nba(n)
       runoff%numrl = runoff%numrl + nop(n) - nba(n)
       if (iam == n) then
          runoff%lnumr = runoff%lnumr + nop(n)
          runoff%endr  = runoff%begr  + runoff%lnumr  - 1
          lnumro = lnumro + nba(n)
          lnumrl = lnumrl + nop(n) - nba(n)
          endro  = begro + lnumro - 1
          endrl  = begrl + lnumrl - 1
       endif
    enddo

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
          write(iulog,*) 'Rtmini ERROR rtm cell count ',n,nba(n),nop(n)
          call shr_sys_abort()
       endif
    enddo

    deallocate(nop,nba,nrs)
    deallocate(iocn,nocn)
    deallocate(pocn)
    call t_stopf('rtmi_dec_distr')

    !--- set some local values

    nroflnd = runoff%numrl
    nrofocn = runoff%numro
    numr = nroflnd + nrofocn
    begr = runoff%begr
    endr = runoff%endr
    call t_stopf('rtmi_decomp')

    !--- Write per-processor runoff bounds depending on dbug level

    call t_startf('rtmi_print')

    call shr_sys_flush(iulog)
    if (masterproc) then
       write(iulog,*) 'total runoff cells numr  = ',runoff%numr
       write(iulog,*) 'total runoff cells numrl = ',runoff%numrl
       write(iulog,*) 'total runoff cells numro = ',runoff%numro
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
          write(iulog,*) 'rtm decomp info',' proc = ',iam, &
             ' begr = ',runoff%begr,&
             ' endr = ',runoff%endr, &
             ' numr = ',runoff%lnumr
          write(iulog,*) '               ',' proc = ',iam, &
             ' begrl= ',begrl,&
             ' endrl= ',endrl, &
             ' numrl= ',lnumrl
          write(iulog,*) '               ',' proc = ',iam, &
             ' begro= ',begro,&
             ' endro= ',endro, &
             ' numro= ',lnumro
       endif
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom_rof,ier)
    enddo

    call t_stopf('rtmi_print')

    !-------------------------------------------------------
    ! Allocate local flux variables
    !-------------------------------------------------------

    call t_startf('rtmi_vars')

    allocate (fluxout (begr:endr,nt_rtm), &
              ddist   (begr:endr),        &
              evel    (begr:endr,nt_rtm), &
              sfluxin (begr:endr,nt_rtm), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation ERROR for ',&
            'volr, fluxout, ddist'
       call shr_sys_abort()
    end if
    fluxout(:,:) = 0._r8
    ddist(:)     = 0._r8
    sfluxin(:,:) = 0._r8

    !-------------------------------------------------------
    ! Allocate runoff datatype 
    !-------------------------------------------------------

    call RunoffInit(begr, endr, numr)

    allocate(rgdc2glo(numr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of runoff%gcd2glo'
       call shr_sys_abort()
    end if

    ! Set map from local to global index space
    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = rglo2gdc(n)
       if (nr /= 0) then
          numr = numr + 1
          rgdc2glo(nr) = n         
          runoff%mask(nr) = gmask(n)
       endif
    end do
    end do


    !-------------------------------------------------------
    ! Initialize runoff data type
    !-------------------------------------------------------

    if (numr /= runoff%numr) then
       write(iulog,*) 'Rtmini ERROR numr and runoff%numr are different ',numr,runoff%numr
       call shr_sys_abort()
    endif
    deallocate(gmask)

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

    ! Set edge longitudes
    rlonw(:) = edgew
    rlone(:) = edgee
    dx = (edgee - edgew) / rtmlon
    do i = 2, rtmlon
       rlonw(i)   = rlonw(i) + (i-1)*dx
       rlone(i-1) = rlonw(i)
    end do

    ! Determine runoff datatype variables
    lrtmarea = 0.0_r8
    do nr = begr,endr
       runoff%gindex(nr) = rgdc2glo(nr)
       n = rgdc2glo(nr)
       i = mod(n-1,rtmlon) + 1
       j = (n-1)/rtmlon + 1
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(iulog,*) 'Rtmini ERROR gdc2glo, nr,ng= ',nr,n
          call shr_sys_abort()
       endif
       runoff%lonc(nr) = runoff%rlon(i)
       runoff%latc(nr) = runoff%rlat(j)

       if (runoff%mask(nr) == 1) then
          do nt = 1,nt_rtm
             runoff%runofflnd(nr,nt) = runoff%runoff(nr,nt)
             runoff%dvolrdtlnd(nr,nt)= runoff%dvolrdt(nr,nt)
             runoff%volrlnd(nr,nt)= runoff%volr(nr,nt)
          enddo
       elseif (runoff%mask(nr) == 2) then
          do nt = 1,nt_rtm
             runoff%runoffocn(nr,nt) = runoff%runoff(nr,nt)
             runoff%dvolrdtocn(nr,nt)= runoff%dvolrdt(nr,nt)
          enddo
       endif

       dx = (rlone(i) - rlonw(i)) * deg2rad
       dy = sin(rlatn(j)*deg2rad) - sin(rlats(j)*deg2rad)
       runoff%area(nr) = 1.e6_r8 * dx*dy*re*re
       lrtmarea = lrtmarea + runoff%area(nr)
       if (dwnstrm_index(n) == 0) then
          runoff%dsi(nr) = 0
       else
          if (rglo2gdc(dwnstrm_index(n)) == 0) then
             write(iulog,*) 'Rtmini ERROR glo2gdc dwnstrm ',&
                  nr,n,dwnstrm_index(n),rglo2gdc(dwnstrm_index(n))
             call shr_sys_abort()
          endif
          runoff%dsi(nr) = rglo2gdc(dwnstrm_index(n))
       endif

    enddo
    deallocate(dwnstrm_index)
    deallocate(rglo2gdc)
    deallocate(rgdc2glo)
    call shr_mpi_sum(lrtmarea,runoff%totarea,mpicom_rof,'rtm totarea',all=.true.)
    if (masterproc) write(iulog,*) 'Rtmini earth area',4.0_r8*shr_const_pi*1.0e6_r8*re*re
    if (masterproc) write(iulog,*) 'Rtmini rtm area  ',runoff%totarea

    !-------------------------------------------------------
    ! Determine downstream distance 
    !-------------------------------------------------------

    ! Instead of reading a distance file calculate the downstream distance
    do nr=begr,endr
       g = runoff%dsi(nr)
       if (g == 0) then
          ddist(nr) = 0._r8
       elseif (g < begr .or. g > endr) then
          write(iulog,*) 'Rtmini: ERROR in ddist calc ',nr,g,begr,endr
          call shr_sys_abort()
       else
          dy  = deg2rad * abs(runoff%latc(nr)-runoff%latc(g)) * re*1000._r8
          dx  = runoff%lonc(nr)-runoff%lonc(g)
          dx1 = abs(dx)
          dx2 = abs(dx+360._r8)
          dx3 = abs(dx-360._r8)
          dx  = min(dx1,dx2,dx3)
          dx  = deg2rad * dx * re*1000._r8 * &
                0.5_r8*(cos(runoff%latc(nr)*deg2rad)+ &
                        cos(runoff%latc(g)*deg2rad))
          ddist(nr) = sqrt(dx*dx + dy*dy)
       endif
    enddo

    !-------------------------------------------------------
    ! Initialize rtm flood - runoff%fthresh and evel
    ! tcraig - recommend RtmFloodInit be called always here and that
    !   the flood_active logic be inside rtmfloodinit, not here.
    !   this will place the fthresh and evel initialization in
    !   a consistent location for extensibility.
    !-------------------------------------------------------
    call RtmFloodInit (frivinp_rtm, begr, endr, runoff%fthresh, evel, &
        flood_active, effvel_active)

    !-------------------------------------------------------
    ! Compute timestep and subcycling number
    !-------------------------------------------------------

    dtover = 0._r8
    dtovermax = 0._r8
    do nt=1,nt_rtm
       do nr=begr,endr
          if (ddist(nr) /= 0._r8) then
             dtover = evel(nr,nt)/ddist(nr)
          else
             dtover = 0._r8
          endif
          dtovermax = max(dtovermax,dtover)
       enddo
    enddo
    dtover = dtovermax
    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom_rof,ier)
    if (dtovermax > 0._r8) then
       delt_rtm_max = (1.0_r8/dtovermax)*cfl_scale
    else
       write(iulog,*) 'rtmini ERROR in delt_rtm_max ',delt_rtm_max,dtover
       call shr_sys_abort()
    endif
    if (masterproc) write(iulog,*) 'rtm max timestep = ',delt_rtm_max,' (sec) for cfl_scale = ',cfl_scale

    call t_stopf('rtmi_vars')

    !-------------------------------------------------------
    ! Read restart/initial info
    !-------------------------------------------------------

    ! The call below opens and closes the file
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
 	(nsrest == nsrBranch  )) then
       call RtmRestFileRead( file=fnamer )
       fluxout(:,:) = runoff%fluxout(:,:)
    end if

    !-------------------------------------------------------
    ! Initialize rtm history handler and fields
    !-------------------------------------------------------

    call RtmHistFldsInit()
    if (nsrest==nsrStartup .or. nsrest==nsrBranch) then
       call RtmHistHtapesBuild()
    end if
    call RtmHistFldsSet()

    if (masterproc) write(iulog,*) subname //':: Success '

  end subroutine Rtmini
  !=======================================================================

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmrun
!
! !INTERFACE:
  subroutine Rtmrun(totrunin, rstwr, nlend, rdate)
!
! !DESCRIPTION:
! River routing model (based on U. Texas code).
! Input is totrunin
! Input/output is fluxout, volr.
! Outputs are dvolrdt\_r, dvolrdt\_lnd\_r, dvolrdt\_ocn\_r, flxocn\_r, flxlnd\_r.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8),         pointer    :: totrunin(:,:)  ! cell tracer lnd forcing on rtm grid (mm/s)
    logical ,         intent(in) :: rstwr          ! true => write restart file this step
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
    integer  :: i, j, n, nr, ns, nt         ! indices
    real(r8) :: dvolrdt                     ! change in storage in discharge units (m3/s)
    real(r8) :: sumfin(nt_rtm),sumfex(nt_rtm)
    real(r8) :: sumrin(nt_rtm),sumdvt(nt_rtm)
    real(r8) :: budget1(nt_rtm)             ! budget check in m3/s
    real(r8) :: budget_global(nt_rtm)       ! global budget sum
    logical  :: budget_check                ! do global budget check
    real(r8),parameter :: budget_tolerance = 1.0e-6   ! budget tolerance, m3/day
    logical  :: abort                       ! abort flag
    real(r8) :: sum1,sum2
    integer  :: yr, mon, day, ymd, tod      ! time information
    integer  :: nsub                        ! subcyling for cfl
    real(r8) :: delt                        ! delt associated with subcycling
    real(r8) :: delt_rtm                    ! real value of rtm_tstep
    integer , save :: nsub_save             ! previous nsub
    real(r8), save :: delt_save             ! previous delt
    logical , save :: first_time = .true.   ! first time flag (for backwards compatibility)
    character(len=256) :: filer             ! restart file name
    integer,parameter  :: dbug = 1          ! local debug flag
    character(*),parameter :: subname = '(Rtmrun) '
!-----------------------------------------------------------------------

    call t_startf('RTMrun')
    call shr_sys_flush(iulog)

    call get_curr_date(yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    if (tod == 0 .and. masterproc) then
       write(iulog,*) ' '
       write(iulog,*) trim(subname),' model date is',ymd,tod
    endif

    delt_rtm = rtm_tstep*1.0_r8
    if (first_time) then
       if (masterproc) write(iulog,*) trim(subname),': rtm act timestep ~ ',delt_rtm
       first_time = .false.
    end if

    budget_check = .false.
    if (day == 1 .and. mon == 1) budget_check = .true.

    ! BUDGET per water tracer, in m3/s = totrunin + (volr_i - volr_f)/dt - flood - runoff ~ 0.0
    ! BUDGET, totrunin input, add initial volume, final volume removed later
    if (budget_check) then
       call t_startf('RTMbudget')
       budget1 = 0.0_r8
       do nt = 1,nt_rtm
       do nr = runoff%begr,runoff%endr
          budget1(nt) = budget1(nt) + totrunin(nr,nt)*runoff%area(nr)*0.001_r8 &
                                    + runoff%volr(nr,nt)/delt_rtm
       enddo
       enddo
       call t_stopf('RTMbudget')
    endif

    ! Remove water from rtm and send back to clm
    ! Just consider land points and only remove liquid water 
    ! runoff%flood needs to be a flux - in units of mm/s
    ! totrunin is a flux (mm/s)

    call t_startf('RTMflood')
    nt = 1 
    do nr = runoff%begr,runoff%endr
       ! initialize runoff%flood to zero
       runoff%flood(nr) = 0._r8
       if (flood_active .and. runoff%mask(nr) == 1) then
          if (runoff%volr(nr,nt) > runoff%fthresh(nr)) then 
             ! determine flux that is sent back to the land
             ! need to convert to mm/s to be consistent with totrunin units
             runoff%flood(nr) = &
                  1000._r8*(runoff%volr(nr,nt)-runoff%fthresh(nr)) / &
                  (delt_rtm*runoff%area(nr))
             ! runoff%flood will be sent back to land - so must subtract this 
             ! from the input runoff from land
             ! tcraig, comment - this seems like an odd approach, you
             !   might create negative forcing.  why not take it out of
             !   the volr directly?  it's also odd to compute this
             !   at the initial time of the time loop.  why not do
             !   it at the end or even during the run loop as the
             !   new volume is computed.  fluxout depends on volr, so
             !   how this is implemented does impact the solution.
             totrunin(nr,nt)= totrunin(nr,nt) - runoff%flood(nr)
          endif
       endif
    enddo
    call t_stopf('RTMflood')
    
    ! BUDGET, flood out, just liquid water term
    if (budget_check) then
       call t_startf('RTMbudget')
       nt = 1
       do nr = runoff%begr,runoff%endr
          budget1(nt) = budget1(nt) - runoff%flood(nr)*runoff%area(nr)*0.001_r8
       enddo
       call t_stopf('RTMbudget')
    endif

    ! Subcycling

    call t_startf('RTMsubcycling')

    nsub = int(delt_rtm/delt_rtm_max) + 1
    delt = delt_rtm/float(nsub)
    if (delt /= delt_save) then
       if (masterproc) then
          write(iulog,*) trim(subname),': rtm delt update from/to',delt_save,delt,nsub_save,nsub
       end if
    endif

    nsub_save = nsub
    delt_save = delt

    sumfin = 0._r8
    sumfex = 0._r8
    sumrin = 0._r8
    sumdvt = 0._r8
    runoff%runoff = 0._r8
    runoff%runofflnd = spval
    runoff%runoffocn = spval
    runoff%dvolrdt = 0._r8
    runoff%dvolrdtlnd = spval
    runoff%dvolrdtocn = spval

    do ns = 1,nsub

       sfluxin = 0._r8
       do n = runoff%begr,runoff%endr
          nr = runoff%dsi(n)
          if (abs(runoff%mask(n)) == 1) then
             if (nr < runoff%begr .or. nr > runoff%endr) then
                write(iulog,*) trim(subname),':Rtm ERROR: non-local communication ',n,nr
                call shr_sys_abort()
             endif
             do nt = 1,nt_rtm
                sfluxin(nr,nt) = sfluxin(nr,nt) + fluxout(n,nt)
             enddo
          endif
       enddo

       if (dbug > 1) then
          do nt = 1,nt_rtm
             sum1 = 0._r8
             sum2 = 0._r8
             do n = runoff%begr,runoff%endr
                sum1 = sum1 + sfluxin(n,nt)
                sum2 = sum2 + fluxout(n,nt)
             enddo
             if (abs(sum1+sum2) > 0.0_r8) then
             if (abs(sum1-sum2)/(sum1+sum2) > 1.0e-12) then
                write(iulog,*) trim(subname),':RTM Warning: fluxin = ',sum1,&
                     ' not equal to fluxout = ',sum2,' for tracer ',nt
             endif
             endif
          enddo
       endif

       do nt = 1,nt_rtm
       do n = runoff%begr,runoff%endr
          dvolrdt = sfluxin(n,nt) + 0.001_r8*totrunin(n,nt)*runoff%area(n) - fluxout(n,nt)

          if (dbug > 1) then
             sumfin(nt) = sumfin(nt) + sfluxin(n,nt)
             sumfex(nt) = sumfex(nt) + fluxout(n,nt)
             sumrin(nt) = sumrin(nt) + 0.001_r8*totrunin(n,nt)*runoff%area(n)
             sumdvt(nt) = sumdvt(nt) + dvolrdt
          endif

          if (abs(runoff%mask(n)) == 1) then         ! land points
             runoff%volr(n,nt) = runoff%volr(n,nt) + dvolrdt*delt
             fluxout(n,nt) = runoff%volr(n,nt) * evel(n,nt)/ddist(n)
             ! --- old cfl constraint.  now use subcycling.  for reference only
             ! fluxout(n)  = min(fluxout(n), volr(n)/delt)
             ! --- this would stop upstream flow if volr/fluxout < 0
             ! --- negative flow largely a function of negative forcing
             ! --- can still have negative runoff where there is negative
             !     forcing over a mask=2 (ocn) point since forcing is put onto
             !     the ocean instantly at these points
             ! --- also, want to allow negative flow so it doesn't build up
             ! fluxout(n) = max(fluxout(n),0._r8)
          else
             runoff%volr(n,nt) = 0._r8
             fluxout(n,nt) = 0._r8
          endif

          if (abs(runoff%mask(n)) == 1) then
             runoff%runoff(n,nt) = runoff%runoff(n,nt) + fluxout(n,nt)
          elseif (runoff%mask(n) == 2) then
             runoff%runoff(n,nt) = runoff%runoff(n,nt) + dvolrdt
          elseif (dvolrdt /= 0.0_r8) then
             ! this water has no where to go.....
             write(iulog,*) subname,' runoff dvolrdt ERROR ',nt,n,runoff%mask(n),dvolrdt
             call shr_sys_abort(subname//' runoff dvolrdt ERROR')
          endif
          ! Convert local dvolrdt (in m3/s) to output dvolrdt (in mm/s)
          runoff%dvolrdt(n,nt) = runoff%dvolrdt(n,nt) + 1000._r8*dvolrdt/runoff%area(n)

       enddo
       enddo

    enddo

    ! average fluxes over subcycling
    runoff%runoff  = runoff%runoff / float(nsub)
    runoff%dvolrdt = runoff%dvolrdt / float(nsub)
    runoff%fluxout = fluxout

    do n = runoff%begr,runoff%endr
       if (runoff%mask(n) == 1) then
          do nt = 1,nt_rtm
             runoff%volrlnd(n,nt)= runoff%volr(n,nt)
             runoff%runofflnd(n,nt) = runoff%runoff(n,nt)
             runoff%dvolrdtlnd(n,nt)= runoff%dvolrdt(n,nt)
          enddo
       elseif (runoff%mask(n) == 2) then
          do nt = 1,nt_rtm
             runoff%runoffocn(n,nt) = runoff%runoff(n,nt)
             runoff%dvolrdtocn(n,nt)= runoff%dvolrdt(n,nt)
          enddo
       endif
    enddo

    call t_stopf('RTMsubcycling')

    ! BUDGET, runoff out (only ocean points), subtract final volume out
    if (budget_check) then
       call t_startf('RTMbudget')
       do nt = 1,nt_rtm
       do nr = runoff%begr,runoff%endr
          if (runoff%mask(nr) == 2) then
             budget1(nt) = budget1(nt) - runoff%runoff(nr,nt)
          endif
          budget1(nt) = budget1(nt) - runoff%volr(nr,nt)/delt_rtm
       enddo
       enddo
       call t_stopf('RTMbudget')
    endif

    ! BUDGET, compute global sums and check
    if (budget_check) then
       call t_startf('RTMbudget')
       call shr_mpi_sum(budget1,budget_global,mpicom_rof,'rtm global budget',all=.false.)
       if (masterproc) then
          write(iulog,*) ' '
          abort = .false.
          do nt = 1,nt_rtm
             budget_global(nt) = budget_global(nt) * 1000._r8 / runoff%totarea * 1.0e6_r8   ! convert from m3/s to kg/m2s*1e6
             write(iulog,'(2a,i10,i6,i4,g20.12)') trim(subname),':BUDGET check (kg/m2s*1e6)= ',ymd,tod,nt,budget_global(nt)
             if (abs(budget_global(nt)) > budget_tolerance) abort = .true.
          enddo
          if (abort) then
             write(iulog,*) trim(subname),':BUDGET abort balance too large ',budget_tolerance
             call shr_sys_abort(trim(subname)//':rtm budget ERROR')
          endif
       endif
       call t_stopf('RTMbudget')
    endif

    ! Write out RTM history and restart file
    call t_startf('RTMhbuf')
    call RtmHistFldsSet()
    call RtmHistUpdateHbuf()
    call t_stopf('RTMhbuf')

    call t_startf('RTMhtapes')
    call RtmHistHtapesWrapup( rstwr, nlend )
    call t_stopf('RTMhtapes')

    if (rstwr) then
       call t_startf('RTMrest')
       filer = RtmRestFileName(rdate=rdate)
       call RtmRestFileWrite( filer, rdate=rdate )
       call t_stopf('RTMrest')
    end if

    ! Global water balance calculation and ERROR check
    if (dbug > 1) then
       do nt = 1,nt_rtm
       if (abs(sumdvt(nt)+sumrin(nt)) > 0.0_r8) then
       if (abs((sumdvt(nt)-sumrin(nt))/(sumdvt(nt)+sumrin(nt))) > 1.0e-6) then
          write(iulog,*) trim(subname),':RTM Warning: water balance nt,dvt,rin,fin,fex = ', &
             nt,sumdvt(nt),sumrin(nt),sumfin(nt),sumfex(nt)
          !call shr_sys_abort()
       endif
       endif
       enddo
    endif
    call shr_sys_flush(iulog)
    call t_stopf('RTMrun')

  end subroutine Rtmrun

  !=======================================================================
  !
  !=======================================================================
  subroutine RtmFloodInit(frivinp, begr, endr, fthresh, evel, &
                          is_rtmflood_on, is_effvel_on )


    !-----------------------------------------------------------------------
    ! Uses
    use pio
    use RtmVar, only :  spval 

    ! Subroutine arguments 
    ! in mode arguments
    character(len=*), intent(in) :: frivinp
    integer ,         intent(in) :: begr, endr
    logical ,         intent(in) :: is_rtmflood_on  !control flooding
    logical ,         intent(in) :: is_effvel_on    !control eff. velocity 

    ! out mode arguments
    real(r8), intent(out) :: fthresh(begr:endr)
    real(r8), intent(out) :: evel(begr:endr,nt_rtm) 

    ! Local dynamically alloc'd variables
    real(r8) , allocatable :: rslope(:)   
    real(r8) , allocatable :: max_volr(:)
    real(r8) , allocatable :: tempr1(:,:),tempr2(:,:) ! temporary buffer for netcdf read

    integer(kind=pio_offset), pointer   :: compdof(:) ! computational degrees of freedom for pio 
    integer :: nt,n,cnt,nr           ! indices
    logical :: readvar               ! read variable in or not
    integer :: ier                   ! status variable
    integer :: dids(2)               ! variable dimension ids 
    integer :: dsizes(2)             ! variable global sizes
    type(file_desc_t)  :: ncid       ! pio file desc
    type(var_desc_t)   :: vardesc1   ! pio variable desc 
    type(var_desc_t)   :: vardesc2   ! pio variable desc 
    type(io_desc_t)    :: iodesc     ! pio io desc
    character(len=256) :: locfn      ! local file name

    !Rtm Flood constants for spatially varying celerity
    real(r8),parameter :: effvel4_0(nt_rtm) = 0.35_r8   ! downstream velocity (m/s)

    real(r8),parameter :: effvel4_5(nt_rtm) = 1.0_r8   ! downstream velocity (m/s)
    real(r8),parameter :: min_ev4_5(nt_rtm) = 0.05_r8  ! minimum downstream velocity (m/s)

    ! name of this subroutine for logging
    character(*),parameter :: subname = '(RtmFloodInit) '
    !-----------------------------------------------------------------------

    !----------------------
    ! if either is_rtmflood_on = .true. or is_effvel_on is .true. then do 
    ! read slope and max_volr out of rdric file.  Below we make the distinction
    ! between using SLOPE (only when is_effvel_on=.true.) and MAX_VOLR (which is
    ! always used when is_rtmflood_on is .true.).
    !----------------------
    if (is_rtmflood_on .or. is_effvel_on) then

       allocate(rslope(begr:endr), max_volr(begr:endr), stat=ier)
       if (ier /= 0) call shr_sys_abort(subname // ':: allocation ERROR')

       ! Get file
       call getfil(frivinp, locfn, 0 )
       if (masterproc) then
          write(iulog,*) subname//':: reading RTM file name for SLOPE and MAX_VOLR: ',&
               trim(frivinp_rtm)
          call shr_sys_flush(iulog)
       endif

       ! Open file and make sure reuqired variables are in file
       ! Assume that if SLOPE is on river input dataset so is MAX_VOLR and that
       ! both have the same io descriptor
       call ncd_pio_openfile (ncid, trim(locfn), 0)
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ier = pio_inq_varid(ncid, name='SLOPE', vardesc=vardesc1)
       if (ier /= PIO_noerr) then
          call shr_sys_abort( trim(subname)//':: ERROR SLOPE not on rdirc file' )
       end if

       ier = pio_inq_varid(ncid, name='MAX_VOLR', vardesc=vardesc2)
       if (ier /= PIO_noerr) then
          call shr_sys_abort( trim(subname)//':: ERROR MAX_VOLR not on rdirc file' )
       end if
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

       ! Set iodesc
       ier = pio_inq_vardimid(ncid, vardesc1, dids)
       ier = pio_inq_dimlen(ncid, dids(1),dsizes(1))
       ier = pio_inq_dimlen(ncid, dids(2),dsizes(2))
       allocate(compdof(runoff%lnumr))
       cnt = 0
       do n = runoff%begr,runoff%endr
          cnt = cnt + 1
          compdof(cnt) = runoff%gindex(n)
       enddo
       call pio_initdecomp(pio_subsystem, pio_double, dsizes, compdof, iodesc)
       deallocate(compdof)
   
       ! Read data
       call pio_read_darray(ncid, vardesc1, iodesc, rslope, ier)
       call pio_read_darray(ncid, vardesc2, iodesc, max_volr, ier)
   
       ! Cleanup and close file
       call pio_freedecomp(ncid, iodesc)
       call pio_closefile(ncid)

    endif
   
    ! done reading rdirc file, now set fthresh and effvel
    if (is_rtmflood_on) then
       do nt = 1,nt_rtm
          do n = begr, endr
             fthresh(n) = max_volr(n)*max(0.005_r8 , 3. * rslope(n))
          end do
       end do
    else
       runoff%fthresh(:) = spval
    endif

    if (is_effvel_on) then
       if (masterproc) write(iulog,*) subname //':: using effvel4_5 '
       do nt = 1,nt_rtm
          do n = begr, endr
             ! modify velocity based on gridcell average slope (Manning eqn)
             evel(n,nt) = max(min_ev4_5(nt),effvel4_5(nt)*sqrt(max(0._r8,rslope(n)))) 
          end do
       end do
    else
       if (masterproc) write(iulog,*) subname //':: using effvel4_0 '
       do nt = 1,nt_rtm
          do nr = begr,endr
             evel(nr,nt) = effvel4_0(nt)
          enddo
       enddo
    endif

    !clean up and exit subroutine
    if (is_rtmflood_on .or. is_effvel_on) then
       deallocate(rslope, max_volr)
    endif

    if (masterproc) write(iulog,*) subname //':: Success '

  end subroutine RtmFloodInit 

  !=======================================================================

end module RtmMod
