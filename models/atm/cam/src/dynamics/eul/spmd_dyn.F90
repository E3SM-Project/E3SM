module spmd_dyn

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM spectral Eulerian dynamics.
! 
! Author: CCM Core Group
! Modified: P. Worley, September 2002, November 2003, December 2003,
!                      November 2004, January 2005, April 2007
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use rgrid,        only: nlon
   use pmgrid,       only: plat, numlats, &
                           beglat, endlat, begirow, endirow, plev
   use spmd_utils,   only: iam, masterproc, npes, proc_smp_map
   use scanslt,      only: beglatex, endlatex, numbnd, numlatsex
   use mpishorthand, only: mpir8, mpicom
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   private
   save

   public spmdinit_dyn, compute_gsfactors, spmdbuf
   public spmd_readnl

   logical, public :: local_dp_map=.false. ! flag indicates that mapping between dynamics 
                                           !  and physics decompositions does not require 
                                           !  interprocess communication
   integer, public :: block_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                           !  in dynamics decomposition (including level 0)
   integer, public :: chunk_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                           !  in physics decomposition (including level 0)

   integer, public, allocatable ::        &
    cut(:,:),                             &! partition for MPI tasks
    cutex(:,:)                             ! extended partition 
   integer, public :: proc(plat)           ! MPI task id associated with a given lat.
   integer, public :: neighs               ! number of south neighbors to comm guardcells
   integer, public, allocatable :: neighs_proc(:)    ! sorted south process neighbors
   integer, public :: neighn               ! number of north neighbors to comm guardcells
   integer, public, allocatable :: neighn_proc(:)    ! sorted north process neighbors
   integer, public :: npessp               ! number of MPI tasks in spectral space
   integer, public :: maxlats              ! max number of lats on any MPI task
   integer, public :: maxcols              ! max number of columns on any MPI task
   integer, public, allocatable :: nlat_p(:)    ! number of latitudes per MPI task
   integer, public, allocatable :: ncol_p(:)    ! number of columns per MPI task
   integer, public :: realloc4_steps       ! number of swaps in realloc4 algorithms
   integer, public, allocatable :: realloc4_proc(:)
                                           ! swap partner in each step of 
                                           ! realloc4 algorithms
   integer, public, allocatable :: realloc4_step(:)
                                           ! step in realloc4 algorithms
                                           ! in which communicate with a given
                                           ! process
   integer, public :: allgather_steps      ! number of swaps in allgather algorithm
   integer, public, allocatable :: allgather_proc(:)
                                           ! swap partner in each step of 
                                           ! allgather (realloc5/7) algorithm
   integer, public, allocatable :: allgather_step(:)
                                           ! step in allgather (realloc5/7) algorithm
                                           ! in which communicate with a given
                                           ! process
!
   logical, private, parameter :: def_equi_by_col = .true.          ! default
   logical, private :: dyn_equi_by_col = def_equi_by_col 
                                           ! flag indicating whether to assign
                                           ! latitudes to equidistribute columns or
                                           ! latitudes. This only matters when using a
                                           ! reduced grid.
!
   logical, private, parameter :: def_mirror = .false.          ! default
   logical, private :: mirror = def_mirror ! flag indicating whether latitudes and their
                                           ! reflections across the equator should assigned 
                                           ! to consecutive processes
!
! Dynamics communication transpose algorithm option:
!  0: use mpi_alltoallv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
   integer, private, parameter :: min_alltoall = 0
   integer, private, parameter :: max_alltoall = 3
   integer, private, parameter :: def_alltoall = 0         ! default
   integer, public :: dyn_alltoall  = def_alltoall
!
! Dynamics communication allgather (realloc5/7) algorithm option:
!  0: use mpi_allgatherv
!  1: use point-to-point MPI-1 two-sided implementation
!  2: use point-to-point MPI-2 one-sided implementation if supported, 
!       otherwise use MPI-1 implementation
!  3: use Co-Array Fortran implementation if supported, 
!       otherwise use MPI-1 implementation
   integer, private, parameter :: min_allgather = 0
   integer, private, parameter :: max_allgather = 3
   integer, private, parameter :: def_allgather = 0         ! default
   integer, public :: dyn_allgather = def_allgather
!
! Dynamics dyn_npes option:
!  1 <=  dyn_npes <= min( 2*(npes/2), plat )
   integer, private, parameter :: min_npes = 1
   integer, private, parameter :: max_npes = plat
   integer, private, parameter :: def_npes = plat
   integer, public :: dyn_npes = def_npes
!
! Dynamics dyn_npes_stride option:
!  1 <=  dyn_npes_stride <= npes/dyn_npes
   integer, private, parameter :: min_npes_stride = 1
   integer, private, parameter :: max_npes_stride = plat
   integer, private, parameter :: def_npes_stride = 1
   integer, public :: dyn_npes_stride = def_npes_stride
!
! MPI communicator for active dynamics processes
!
   integer, public :: mpicom_dyn_active    
!
! Collective communication send/receive buffers
#if (defined CAF)
   real(r8), public, allocatable :: buf1(:)[:],buf2(:)[:] ! buffers for packing MPI msgs
#else
   real(r8), public, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs
#endif
   integer, public :: spmdbuf_siz = 0        ! buffer size (in r8s)
   integer, public :: buf1win                ! buf1 Window id
   integer, public :: buf2win                ! buf2 Window id

contains
   
!----------------------------------------------------------------------

  subroutine spmd_readnl(nlfilename)

    ! !USES:
    use units,           only: getunit, freeunit	
    use namelist_utils,  only: find_group_name	
    use spmd_utils,      only: npes, masterproc	
    use pmgrid,          only: plat, plev, plon	
    use mpishorthand	
     
    implicit none
     
     !
     ! !PARAMETERS:
   character(len=*), intent(in) :: nlfilename

! !DESCRIPTION: Read in EUL-specific namelist variables.  Must be 
!               performed before dyn\_init
!
! !REVISION HISTORY:
!   2010.05.15   Sawyer  Creation
!
!EOP
!=========================================================================
!BOC
! Local variables
   integer :: ierr           ! error code
   integer :: unitn          ! namelist unit number
   character(len=*), parameter ::  subname = "spmd_readnl"

   namelist /spmd_dyn_inparm/ dyn_alltoall, &
             dyn_allgather,  &
             dyn_equi_by_col,&
             dyn_npes,       &
             dyn_npes_stride 

   if (masterproc) then
      write(iulog,*) 'Read in spmd_dyn_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )
      
      ! Look for dyn_eul_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'spmd_dyn_inparm', status=ierr)
      if (ierr == 0) then  ! found spmd_dyn_inparm
         read(unitn, spmd_dyn_inparm, iostat=ierr)  ! read the spmd_dyn_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                 ' error condition for spmd_dyn_inparm' )
         end if
      end if
      close( unitn )
      call freeunit( unitn )
   endif

   call mpibcast (dyn_alltoall   ,1,mpiint,0,mpicom)
   call mpibcast (dyn_allgather  ,1,mpiint,0,mpicom)
   call mpibcast (dyn_equi_by_col,1,mpilog,0,mpicom)
   call mpibcast (dyn_npes       ,1,mpiint,0,mpicom)
   call mpibcast (dyn_npes_stride,1,mpiint,0,mpicom)

   if ((dyn_alltoall.lt.min_alltoall).or. &
        (dyn_alltoall.gt.max_alltoall)) then
      write(iulog,*)                                          &
           'spmd_readnl:  ERROR:  dyn_alltoall=', &
           dyn_alltoall,                              &
           '  is out of range.  It must be between ',        &
           min_alltoall,' and ',max_alltoall
      call endrun
   endif
   
   if ((dyn_allgather.lt.min_allgather).or. &
        (dyn_allgather.gt.max_allgather)) then
      write(iulog,*)                                          &
           'spmd_readnl:  ERROR:  dyn_allgather=', &
           dyn_allgather,                              &
           '  is out of range.  It must be between ',        &
           min_allgather,' and ',max_allgather
      call endrun
   endif
   !
   if ((dyn_npes.lt.min_npes).or. &
        (dyn_npes.gt.max_npes)) then
      write(iulog,*)                                          &
           'spmd_readnl:  ERROR:  dyn_npes=', &
           dyn_npes,                              &
           '  is out of range.  It must be between ',        &
           min_npes,' and ',max_npes
      call endrun
   endif
   !
   if ((dyn_npes_stride.lt.min_npes_stride).or. &
        (dyn_npes_stride.gt.max_npes_stride)) then
      write(iulog,*)                                          &
           'spmd_readnl:  ERROR:  dyn_npes_stride=', &
           dyn_npes_stride,                              &
           '  is out of range.  It must be between ',        &
           min_npes_stride,' and ',max_npes_stride
      call endrun
   endif
   
   
 end subroutine spmd_readnl


!========================================================================

   subroutine spmdinit_dyn ()
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processes
! 
! Method: Distribution is S->N for processes 0->dyn_npes
! 
! Author: CCM Core Group
! Modified: P. Worley, November 2003 to improve SMP load balance, and to
!           change distribution to 
!             S->E for processes 0,2,..,dyn_npes-2
!           and 
!             N->E for processes 1,3,..,dyn_npes-1
!           when mirror flag is set (at request of physics)
! Modified: P. Worley, November 2004 to improve load balance for 
!           reduced grid by equidistributing columns (not latitudes)
!           in latitude decomposition. Used when equi_by_col flag is set.
!           On by default, and gives identical decomposition as
!           equidistributing by latitude when using a full grid.
! Modified: P. Worley, April 2007 to support idle processes when
!           in the dynamics (dyn_npes < npes)
! 
!-----------------------------------------------------------------------
      use comspe, only: numm
      use spmd_utils
#if (defined MODCM_DP_TRANSPOSE)
      use parutilitiesmodule, only : parinit
#endif
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer i         ! loop index
      integer tot_cols  ! total number of columns in computational grid
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer tot_nx    ! total number of latitudes/columns in 
                        ! computational grid
      integer nx_base   ! approx. number of latitudes/columns per proc
      integer nx_p(0:npes-1)      ! number of latitudes/columns per process
      integer nx_smp(0:npes-1)    ! number of latitudes/columns per SMP
      integer nproc_smp(0:npes-1) ! number of MPI processes per SMP
      integer workleft  ! amount of work still to be parcelled out

      integer smpid     ! SMP id
      integer smpids    ! SMP id for SH process
      integer smpidn    ! SMP id for NH process
      integer procj     ! process offset loop index
      integer procid    ! process id
      integer procids   ! process id SH
      integer procidn   ! process id NH
      integer procid_s  ! strided process id
      integer procids_s ! strided process id SH
      integer procidn_s ! strided process id NH

      integer max_ncols ! maximum number of columns assigned to a process
      integer min_max_ncols    ! minmax number of columns assigned 
                               ! to a process over all latitude assignments
      integer ncol      ! number of columns assigned to current process
      integer ncol_curtot  ! current total number of columns assigned
      integer ncol_curgoal ! target number of columns to be assigned to process
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer neighn_minlat(plat)    ! minimum latitude in north neighbor
      integer neighs_maxlat(plat)    ! maximum latitude in south neighbor
      integer active_proc            ! +1 for active dynamics processes
      integer ierror                 ! MPI error return

      real(r8) avgnx_proc(0:npes-1) ! average number of latitudes/columns per 
                                    ! MPI process in a given SMP node
      real(r8) minavgnx_proc        ! minimum average number of 
                                    ! latitudes/columns per 
                                    ! MPI process over SMP nodes
      real(r8) alpha    ! slop factor in assigning latitudes to processes
      real(r8) opt_alpha! best slop factor in assigning latitudes to processes

      logical done      ! exit flag for latitude assignment loop
!
!-----------------------------------------------------------------------
!
! Initialize Pilgrim library
!
#if (defined MODCM_DP_TRANSPOSE)
      call parinit(mpicom)
#endif
!
! Initialize mirror flag
!
      mirror = phys_mirror_decomp_req
!
! Allocate memory for MPI task partition array
! and extended partition
!
      allocate (cut  (2,0:npes-1))
      cut(1,0:npes-1)  = 1
      cut(2,0:npes-1)  = 0
!
      allocate (cutex(2,0:npes-1))
      cutex(1,0:npes-1) = 1
      cutex(2,0:npes-1) = 0
!
! Allocate memory for number of lats per proc
!
      allocate (nlat_p (0:npes-1))
      nlat_p(0:npes-1) = 0
!
! Allocate memory for number of columns per proc
!
      allocate (ncol_p (0:npes-1))
      ncol_p(0:npes-1) = 0
!
! determine total number of columns
!
      tot_cols = 0
      do lat=1,plat
         tot_cols = tot_cols + nlon(lat)
      enddo
!
! Make sure number of PEs, latitudes, and columns are kosher
!
      call factor (plat, m2, m3, m5)

      if (m2 < 1) then
         call endrun ('SPMDINIT_DYN: Problem size is not divisible by 2')
      end if

      if (masterproc) then
         write(iulog,*) 'Problem factors: 2**',m2,' * 3**',m3,' * 5**',m5
      end if

      if (npes > 1) then
         if (dyn_npes > min( 2*(npes/2), plat ) ) then
            dyn_npes = min( 2*(npes/2), plat )
         endif
         if (dyn_npes_stride > npes/dyn_npes) then
            dyn_npes_stride = npes/dyn_npes
         endif
      else
         dyn_npes = 1
         dyn_npes_stride = 1
      endif

      if ((dyn_equi_by_col) .and. (mod(tot_cols,2) /= 0)) then
         write(iulog,*)'SPMDINIT_DYN: Total number of columns(', &
                   tot_cols,') must be a multiple of 2'
         call endrun
      end if
!
!  Initialization for inactive processes
!
      beglat  = 1
      endlat  = 0
      numlats = 0
      begirow = 1
      endirow = 0

      beglatex = 1
      endlatex = 0
      numlatsex = 0
!
! Special initialization for dyn_npes == 1 case
!
      if (dyn_npes .eq. 1) then
!
         nlat_p(0) = plat
         cut(1,0)  = 1
         cut(2,0)  = plat
!
         ncol_p(0) = 0
         do lat=1,plat
            ncol_p(0) = ncol_p(0) + nlon(lat)
         enddo
!
         if (iam .eq. 0) then
            beglat  = 1
            endlat  = plat
            numlats = plat
            begirow = 1
            endirow = plat/2
         endif
!
      else
!
! Determine approximate number of columns or latitudes per process
!
         if (dyn_equi_by_col) then
            tot_nx = tot_cols
         else
            tot_nx = plat
         endif
         nx_base = tot_nx/dyn_npes
         do procid=0,dyn_npes-1
            procid_s = dyn_npes_stride*procid
            nx_p(procid_s) = nx_base
         enddo
!
! Calculate initial distribution of columns or latitudes and 
! distribution of processes by SMP
!
         nx_smp(0:npes-1) = 0
         nproc_smp(0:npes-1) = 0
         do procid=0,dyn_npes-1
            procid_s = dyn_npes_stride*procid
            smpid = proc_smp_map(procid_s)
            nproc_smp(smpid) = nproc_smp(smpid) + 1
         enddo
!
         do smpid=0,nsmps-1
            nx_smp(smpid)     = nx_base*nproc_smp(smpid)
            avgnx_proc(smpid) = real(nx_base,r8)
         enddo
!
! Equi-distribute remaining columns or latitudes across SMPs
! without increasing per process imbalance beyond minimum
!
         workleft = tot_nx - dyn_npes*nx_base
         do while (workleft > 0)
!
! (a) Find minimun number of columns or latitudes assigned to an SMP
!
            minavgnx_proc = avgnx_proc(0)
            do smpid=1,nsmps-1
               if (minavgnx_proc > avgnx_proc(smpid)) then
                  minavgnx_proc = avgnx_proc(smpid)
               endif
            enddo
!
! (b) Assign an additional column or latitude to processes with 
!     nx_base latitudes/columns in SMPs with the minimum 
!     average number of latitudes/columns
!
            do procid=dyn_npes/2-1,0,-1
               if (mirror) then
                  procids = 2*procid
                  procidn = procids + 1
               else
                  procids = procid
                  procidn = dyn_npes - procids - 1
               endif
!
               procids_s = dyn_npes_stride*procids
               procidn_s = dyn_npes_stride*procidn
!
               smpids = proc_smp_map(procids_s)
               smpidn = proc_smp_map(procidn_s)
               if ((nx_p(procids_s) .eq. nx_base)  .and. &
                   ((avgnx_proc(smpids) .eq. minavgnx_proc) .or. &
                    (avgnx_proc(smpidn) .eq. minavgnx_proc)) .and. &
                   (workleft > 0)) then
!
                  nx_p(procids_s) = nx_p(procids_s) + 1
                  nx_smp(smpids) = nx_smp(smpids) + 1
                  avgnx_proc(smpids) = &
                     real(nx_smp(smpids),r8)/real(nproc_smp(smpids),r8)
!
                  nx_p(procidn_s) = nx_p(procids_s)
                  nx_smp(smpidn) = nx_smp(smpidn) + 1
                  avgnx_proc(smpidn) = &
                     real(nx_smp(smpidn),r8)/real(nproc_smp(smpidn),r8)
!
                  workleft = workleft - 2
               endif
            enddo
         end do
!
! Partition latitudes over processes, equidistributing either 
! a) columns, or
! b) latitudes
!
         if (dyn_equi_by_col) then
!
! Evaluate different latitude assignments
!
            min_max_ncols = tot_cols
            do i=0,10
               alpha = .05_r8*i
               max_ncols = 0
!
               iend = 0
               ncol_curtot  = 0
               ncol_curgoal = 0
               do procid=0,dyn_npes/2-1
                  if (mirror) then
                     procids = 2*procid
                  else
                     procids = procid
                  endif
                  procids_s = dyn_npes_stride*procids
                  ncol_curgoal = ncol_curgoal + nx_p(procids_s)
                  ncol = 0
!
                  done = .false.
!
! Add latitudes until near column per process goal for current process
!
                  do while ((.not. done) .and. &
                     (ncol_curtot < ncol_curgoal))
                     if (iend .ge. plat/2) then
                        write(iulog,*)'SPMDINIT_DYN: error in assigning latitudes to processes'
                        call endrun
                     endif
                     if (ncol_curtot + nlon(iend+1) .le. &
                         ncol_curgoal + alpha*nlon(iend+1)) then
                        iend = iend + 1
                        ncol = ncol + nlon(iend)
                        ncol_curtot = ncol_curtot + nlon(iend)
                     else
                        done = .true.
                     endif
                  enddo
                  if (ncol > max_ncols) max_ncols = ncol
!
               enddo
               if (max_ncols < min_max_ncols) then
                  min_max_ncols = max_ncols
                  opt_alpha = alpha
               endif
            enddo
!
! Determine latitude assignments when equidistributing columns
!
            iend = 0
            ncol_curtot = 0
            ncol_curgoal = 0
            do procid=0,dyn_npes/2-1
               if (mirror) then
                  procids = 2*procid
                  procidn = procids + 1
               else
                  procids = procid
                  procidn = dyn_npes - procids - 1
               endif
!
               procids_s = dyn_npes_stride*procids
               procidn_s = dyn_npes_stride*procidn
!
               ncol_curgoal = ncol_curgoal + nx_p(procids_s)
               ncol_p(procids_s) = 0
!
               cut(1,procids_s) = iend + 1
               cut(2,procids_s) = iend
               done = .false.
!
! Add latitudes until near column per process goal for current process
!
               do while ((.not. done) .and. &
                  (ncol_curtot < ncol_curgoal))
                  if (ncol_curtot + nlon(iend+1) .le. &
                     ncol_curgoal + opt_alpha*nlon(iend+1)) then
                     iend = iend + 1
                     cut(2,procids_s) = iend
                     ncol_p(procids_s) = ncol_p(procids_s) + nlon(iend)
                     ncol_curtot = ncol_curtot + nlon(iend)
                     nlat_p(procids_s) = nlat_p(procids_s) + 1
                  else
                     done = .true.
                  endif
               enddo
!
! Assign mirror latitudes
!
               cut(1,procidn_s) = plat - cut(2,procids_s) + 1
               cut(2,procidn_s) = plat - cut(1,procids_s) + 1
               ncol_p(procidn_s) = ncol_p(procids_s)
               nlat_p(procidn_s) = nlat_p(procids_s)
!
! Save local information
!
               if (iam == procids_s .or. iam == procidn_s) then
                  beglat = cut(1,iam)
                  endlat = cut(2,iam)
                  numlats = nlat_p(iam)
                  begirow = cut(1,procids_s)
                  endirow = cut(2,procids_s)
               end if
!
            enddo
!
         else
!
! Determine latitude assignments when
! equidistributing latitudes
!
            iend = 0
            do procid=0,dyn_npes/2-1
               if (mirror) then
                  procids = 2*procid
                  procidn = procids + 1
               else
                  procids = procid
                  procidn = dyn_npes - procids - 1
               endif
!
               procids_s = dyn_npes_stride*procids
               procidn_s = dyn_npes_stride*procidn
!
               nlat_p(procids_s) = nx_p(procids_s)
               cut(1,procids_s) = iend + 1
               cut(2,procids_s) = iend + nlat_p(procids_s)
               iend = iend + nlat_p(procids_s)
!
               ncol_p(procids_s) = 0
               do lat=cut(1,procids_s),cut(2,procids_s)
                  ncol_p(procids_s) = ncol_p(procids_s) + nlon(lat)
               enddo
!
! Assign mirror latitudes
!
               nlat_p(procidn_s) = nx_p(procidn_s)
               cut(1,procidn_s) = plat - cut(2,procids_s) + 1
               cut(2,procidn_s) = plat - cut(1,procids_s) + 1
!
               ncol_p(procidn_s) = 0
               do lat=cut(1,procidn_s),cut(2,procidn_s)
                  ncol_p(procidn_s) = ncol_p(procidn_s) + nlon(lat)
               enddo
!
! Save local information
!
               if (iam == procids_s .or. iam == procidn_s) then
                  beglat = cut(1,iam)
                  endlat = cut(2,iam)
                  numlats = nlat_p(iam)
                  begirow = cut(1,procids_s)
                  endirow = cut(2,procids_s)
               end if
!
            enddo
         endif
!
      endif
!
! Calculate maximum number of latitudes and columns assigned to a process
!
      maxlats = maxval(nlat_p)
      maxcols = maxval(ncol_p)
!
      do procid=0,dyn_npes-1
         procid_s = dyn_npes_stride*procid
         if (masterproc) then
            write(iulog,*)'procid ',procid_s,' assigned ', &
                      cut(2,procid_s)-cut(1,procid_s)+1,' latitude values from', &
                      cut(1,procid_s),' through ',cut(2,procid_s),' containing', &
                      ncol_p(procid_s),' vertical columns'
         end if
!
! Determine which process is responsible for the defined latitudes
!
         do lat=cut(1,procid_s),cut(2,procid_s)
            proc(lat) = procid_s
         end do
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
         cutex(1,procid_s) = cut(1,procid_s) - numbnd
         cutex(2,procid_s) = cut(2,procid_s) + numbnd
         if (iam == procid_s) then
            beglatex = cutex(1,procid_s) + numbnd
            endlatex = cutex(2,procid_s) + numbnd
            numlatsex = endlatex - beglatex + 1
         end if
      end do
!
! Determine neighbor processes needed for boundary communication.  
! North first.
!
      neighn = 0
      neighn_minlat(:) = -1
      do procid=0,dyn_npes-1
         procid_s = dyn_npes_stride*procid
         if (procid_s /= iam) then
            if ((cut(1,procid_s) > cut(2,iam)) .and. &
                (cut(1,procid_s) <= cut(2,iam)+numbnd)) then
               neighn_minlat(cut(1,procid_s)) = procid_s
               neighn = neighn + 1
            endif
         endif
      enddo
!
! Sort north processes by increasing latitude
!
      allocate (neighn_proc (neighn))
      neighn = 0
      do lat=1,plat
         if (neighn_minlat(lat) /= -1) then
            neighn = neighn + 1
            neighn_proc(neighn) = neighn_minlat(lat)
         endif
      enddo
!
! South next.
!
      neighs = 0
      neighs_maxlat(:) = -1
      do procid=0,dyn_npes-1
         procid_s = dyn_npes_stride*procid
         if (procid_s /= iam) then
            if ((cut(2,procid_s) < cut(1,iam)) .and. &
                (cut(2,procid_s) >= cut(1,iam)-numbnd)) then
               neighs_maxlat(cut(2,procid_s)) = procid_s
               neighs = neighs + 1
            endif
         endif
      enddo
!
! Sort south processes by decreasing latitude
!
      allocate (neighs_proc (neighs))
      neighs = 0
      do lat=plat,1,-1
         if (neighs_maxlat(lat) /= -1) then
            neighs = neighs + 1
            neighs_proc(neighs) = neighs_maxlat(lat)
         endif
      enddo
!
      if (masterproc) then
         write(iulog,*)'-----------------------------------------'
         write(iulog,*)'Number of lats passed north & south = ',numbnd
         write(iulog,*)'Node  Partition  Extended Partition'
         write(iulog,*)'-----------------------------------------'
         do procid=0,dyn_npes-1
            procid_s = dyn_npes_stride*procid
            write(iulog,200) procid_s,cut(1,procid_s),cut(2,procid_s) ,cutex(1,procid_s), &
                         cutex(2,procid_s)
200         format(i3,4x,i3,'-',i3,7x,i3,'-',i3)
         end do
      end if
!      write(iulog,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(iulog,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      call decomp_wavenumbers ()
!
! Make communicator for active dynamics processors (for use in realloc4a/4b)
      if (beglat <= endlat) then
         active_proc = 1
      else
         active_proc = 0
      endif
      call mpi_comm_split(mpicom, active_proc, iam, mpicom_dyn_active, ierror)
!
! Precompute swap partners and number of steps in realloc4 alltoall algorithm.
! First, determine number of swaps.
!
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            if (((numm(iam) > 0) .and. (nlat_p(procid) > 0)) .or. &
               ((numm(procid) > 0) .and. (numlats > 0))) then
               realloc4_steps = realloc4_steps + 1
            end if
         end if
      end do
!
! Second, determine swap partners.
!
      allocate( realloc4_proc(realloc4_steps) )
      allocate( realloc4_step(0:npes-1) )
      realloc4_step(:) = -1
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            if (((numm(iam) > 0) .and. (nlat_p(procid) > 0)) .or. &
               ((numm(procid) > 0) .and. (numlats > 0))) then
               realloc4_steps = realloc4_steps + 1
               realloc4_proc(realloc4_steps) = procid
               realloc4_step(procid) = realloc4_steps
            end if
         end if
      end do
!
! Precompute swap partners in realloc5/7 allgather algorithm.
      allocate( allgather_proc(npes-1) )
      allocate( allgather_step(0:npes-1) )
      allgather_step(:) = -1
      allgather_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            allgather_steps = allgather_steps + 1
            allgather_proc(allgather_steps) = procid
            allgather_step(procid) = allgather_steps
         end if
      end do
!
      return
   end subroutine spmdinit_dyn

!========================================================================

   subroutine factor (nitems, m2, m3, m5)
!----------------------------------------------------------------------- 
! 
! Purpose: Factor a given number into powers of 2,3,5
! 
! Method: Brute force application of "mod" function
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nitems      ! Number to be factored into powers of 2,3,5
      integer, intent(out) :: m2,m3,m5   ! Powers of 2, 3, and 5 respectively
!
! Local workspace
!
      integer num                        ! current number to be factored
!
!-----------------------------------------------------------------------
!
      num = nitems
      m2 = 0
      m3 = 0
      m5 = 0
      
2     if (mod(num,2) == 0) then
         m2 = m2 + 1
         num = num/2
         goto 2
      end if
      
3     if (mod(num,3) == 0) then
         m3 = m3 + 1
         num = num/3
         goto 3
      end if
      
5     if (mod(num,5) == 0) then
         m5 = m5 + 1
         num = num/5
         goto 5
      end if
      
      if (num /= 1) then
         write(iulog,*) 'FACTOR: ',nitems,' has a prime factor other than 2, 3, or 5.  Aborting...'
         call endrun
      end if
      
      return
   end subroutine factor

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processes
! 
! Method: Approximately equidistribute both the number of spectral 
!         coefficients and the number of wavenumbers assigned to each 
!         MPI task using a modified version of the mapping due to
!         Barros and Kauranne. 
! 
! Author: P. Worley, September 2002
! 
!-----------------------------------------------------------------------
      use pspect, only: pmmax
      use comspe, only: numm, maxm, locm, locrm, nlen, lpspt, lnstart
!
! Local workspace
!
      integer procid      ! process id
      integer procid_s    ! strided process id
      integer m, lm       ! global and local fourier wavenumber indices
      integer mstride     ! Stride over wavenumbers used in decomposition
      integer begm1       ! Starting Fourier wavenumbers owned by an MPI task
      integer begm2       !  when using Barros & Kauranne decomposition
      integer speccount(0:npes-1)
                          ! number of spectral coefficients assigned to
                          ! each MPI task
!-----------------------------------------------------------------------
!
! determine upper bound on number of wavenumbers to be assigned to each 
! process
      if (mod(pmmax,dyn_npes) .eq. 0) then
         maxm = pmmax/dyn_npes
      else
         maxm = (pmmax/dyn_npes) + 1
      endif
      allocate ( numm(0:npes-1) )
      allocate ( locm(1:maxm, 0:npes-1) )
      allocate ( locrm(1:2*maxm, 0:npes-1) )
!
! assign wavenumbers to approximately equidistribute the number 
! of spectral coefficients assigned to each process
      numm(:) = 0
      locm(:,:) = huge(1)
      locrm(:,:) = huge(1)
      speccount(:) = 0
      mstride = 2*dyn_npes
      npessp = 0
      do procid = 0,dyn_npes-1
         procid_s = dyn_npes_stride*procid
         begm1 = procid + 1
         begm2 = mstride - procid
         do m=begm1,pmmax,mstride
            numm(procid_s) = numm(procid_s) + 1
            locm(numm(procid_s),procid_s) = m
            speccount(procid_s) = speccount(procid_s) + nlen(m)
         enddo
         do m=begm2,pmmax,mstride
            numm(procid_s) = numm(procid_s) + 1
            locm(numm(procid_s),procid_s) = m
            speccount(procid_s) = speccount(procid_s) + nlen(m)
         enddo
!
         if (numm(procid_s) .gt. 0) then
            npessp = npessp + 1
         endif
!
      enddo
!
      do procid = 0,dyn_npes-1
         procid_s = dyn_npes_stride*procid
         if (masterproc) then
            write(iulog,*)'procid ',procid_s,' assigned ', speccount(procid_s), &
                      ' spectral coefficients and ', numm(procid_s), &
                      ' m values: ', (locm(lm,procid_s),lm=1,numm(procid_s))
         end if
         do lm=1,numm(procid_s)
            locrm(2*lm-1,procid_s) = 2*locm(lm,procid_s)-1
            locrm(2*lm  ,procid_s) = 2*locm(lm,procid_s)
         enddo
      enddo
!
! Calculate number of local spectral coefficients
      lpspt = 0
      do lm=1,numm(iam)
         lpspt = lpspt + nlen(locm(lm,iam))
      enddo
!
! Evaluate displacement info based on truncation params and
! wavenumber assignment
      allocate ( lnstart(1:maxm) )
      lnstart(1) = 0
      do lm=2,numm(iam)
         lnstart(lm) = lnstart(lm-1) + nlen(locm(lm-1,iam))
      enddo
!   
      return
   end subroutine decomp_wavenumbers

!========================================================================

  subroutine spmdbuf 
!----------------------------------------------------------------------- 
! 
! Purpose: allocate spmd pack buffers used in collective communications
! 
! Author: CCM Core Group
!
! Note: Call after phys_grid_init
! 
!-----------------------------------------------------------------------
     use error_messages, only: alloc_err
     use comspe,         only: nlen, maxm
     use constituents,   only: pcnst
!-----------------------------------------------------------------------
!
! Local workspace
!
     integer :: maxcount(5),m
     integer :: length,i,lm,istat1,istat2
     integer :: bsiz, glb_bsiz       ! buffer size (in bytes)
!
! realloc4a max: 8  2 plev*numm*numlats (e.g. tdyn)
!                1  2     *numm*numlats (bpstr)
!
     maxcount(1) = (npes-1)*maxlats*(2*maxm*(plev*8 + 1))
!
! realloc4b max: 8  2 plev*numm*numlats (e.g. vort)
!                4  2     *numm*numlats (e.g. dps)
!
     maxcount(2) = (npes-1)*maxlats*(2*maxm*(plev*8 + 4))
!
! realloc5 max: 6 numlats         (e.g. tmass)
!               5 numlats  *pcnst (e.g. hw1lat)
!               2 4*numlats*pcnst (e.g. hw2al)
!
     maxcount(3) = npes*maxlats*(6 + (5 + 2*4)*pcnst)
!
! realloc7 max: 3 plev *numlats    (e.g. vmax2d)
!               5      *numlats    (e.g. psurf)
!
     maxcount(4) = npes*maxlats*(3*plev + 5)
!
! dp_coupling max:
!
     if (.not. local_dp_map) then
        maxcount(5) = (5 + pcnst)*max(block_buf_nrecs,chunk_buf_nrecs)
     else
        maxcount(5) = 0
     endif
!
     m = maxval(maxcount)
     call mpipack_size (m, mpir8, mpicom, bsiz)
     call mpiallmaxint(bsiz, glb_bsiz, 1, mpicom)
     if (masterproc) then
        write(iulog,*) 'SPMDBUF: Allocating SPMD buffers of size ',glb_bsiz
     endif
     spmdbuf_siz = glb_bsiz/8 + 1
#if (defined CAF)
     allocate(buf1(spmdbuf_siz)[*], stat=istat1)
     allocate(buf2(spmdbuf_siz)[*], stat=istat2)
#else
     allocate(buf1(spmdbuf_siz), stat=istat1)
     allocate(buf2(spmdbuf_siz), stat=istat2)
#endif
     call alloc_err( istat1, 'spmdbuf', 'buf1', spmdbuf_siz )
     call alloc_err( istat2, 'spmdbuf', 'buf2', spmdbuf_siz )
     call mpiwincreate(buf1,spmdbuf_siz*8,mpicom,buf1win)
     call mpiwincreate(buf2,spmdbuf_siz*8,mpicom,buf2win)
     buf1 = 0.0_r8
     buf2 = 0.0_r8
     return
  end subroutine spmdbuf

!========================================================================

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = numperlat*(cut(1,p)-1)
     end do
     
  end subroutine compute_gsfactors

#endif

end module spmd_dyn
