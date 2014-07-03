
module spmd_dyn
!BOP
!
! !MODULE: Subroutines to initialize SPMD implementation of CAM
!

#if (defined SPMD)

!
! !USES:
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: iam, masterproc, npes
   use pmgrid,             only: plat, plon, numbnd, &
                                 numlats, beglat, endlat, &
                                 plev, beglev, endlev, endlevp1, &
                                 endlevp, myid_y, myid_z, npr_y, npr_z, plevp, &
                                 myidxy_x, myidxy_y, nprxy_x, nprxy_y, &
                                 beglonxy, endlonxy, beglatxy, endlatxy, &
                                 twod_decomp, spmd_on, mod_transpose, mod_geopk, &
                                 mod_gatscat
   use mpishorthand,       only: mpir8, mpicom, mpiint, mpi_success
   use decompmodule,       only: decomptype, decompcreate
   use ghostmodule,        only: ghosttype
   use parutilitiesmodule, only: parpatterntype
   use fv_control_mod,     only: ct_overlap, trac_decomp
   use abortutils,         only: endrun
   use cam_logfile,        only: iulog

   implicit none

! !PUBLIC MEMBER FUNCTIONS:

   public spmdinit_dyn, decomp_wavenumbers
   public compute_gsfactors, spmdbuf, spmd_readnl

! !PUBLIC DATA MEMBERS:

   integer ::force_2d = 0                 !option to force transpose computation for 1D decomp.
   integer :: geopkblocks = 1              !number of stages to use in Z-serial non-transpose
                                                  ! geopotential method (routine geopk_d)
                                                  ! with 2D decomp.
   logical :: geopkdist  = .false.         !use a distributed method for geopotential calculation 
                                                  ! with 2D decomp.
   logical :: geopk16byte   = .false.      !use Z-parallel distributed method for geopotential 
                                                  ! calculation with 2D decomp.; otherwise, use Z-serial 
                                                  ! pipeline algorithm
   integer :: geopktrans = 0               
   integer :: npr_yz(4)                    !yz and xy decompositions
   integer :: modcomm_transpose = 0        !mod_comm transpose method
                                                  !   0 for temporary contiguous buffers
                                                  !   1 for mpi derived types
   integer :: modcomm_geopk = 0            !mod_comm geopk method
                                                  !   0 for temporary contiguous buffers
                                                  !   1 for mpi derived types
   integer :: modcomm_gatscat = 0          !mod_comm gather/scatter method
                                                  !   0 for temporary contiguous buffers
                                                  !   1 for mpi derived types
   integer :: modc_sw_dynrun = 0           !mod_comm irregular underlying communication method for dyn_run/misc
                                                  !  0 for original mp_sendirr/mp_recvirr
                                                  !  1 for mp_swapirr and point-to-point communications
                                                  !  2 for mp_swapirr and all-to-all communications
   logical :: modc_hs_dynrun = .true.      !mod_comm irreg comm handshaking for dyn_run/misc
   logical :: modc_send_dynrun = .true.    ! true for mod_comm irregular communication blocking send for
                                                  ! dyn_run/misc, false for nonblocking send
   integer :: modc_mxreq_dynrun = -1       !maximum number of nonblocking communication requests to allow
                                                  ! when using mp_swapirr and point-to-point communications for
                                                  ! dyn_run/misc
                                                  ! < 0 implies no limits
   integer :: modc_sw_cdcore = 0           !mod_comm irregular underlying communication method for cd_core/geopk
                                                  !  0 for original mp_sendirr/mp_recvirr
                                                  !  1 for mp_swapirr and point-to-point communications
                                                  !  2 for mp_swapirr and all-to-all communications
   logical :: modc_hs_cdcore = .true.      ! true for mod_comm irregular communication handshaking for cd_core/geopk
   logical :: modc_send_cdcore  = .true.   ! true for geopk_d or mod_comm irregular communication blocking send for
                                                  !  cd_core/geopk, false for nonblocking send
   integer :: modc_mxreq_cdcore = -1       ! maximum number of nonblocking communication requests to allow
                                                  !  when using mp_swapirr and point-to-point communications for
                                                  !  cd_core/geopk
                                                  !  < 0 implies no limits
   integer :: modc_sw_gather = 1           ! mod_comm irregular underlying communication method for gather
                                                  !  0 for original mp_sendirr/mp_recvirr
                                                  !  1 for mp_swapirr and point-to-point communications
                                                  !  2 for mp_swapirr and all-to-all communications
   logical :: modc_hs_gather = .true.      ! true for mod_comm irregular communication handshaking for gather
   logical :: modc_send_gather = .true.    ! true for mod_comm irregular communication blocking send for
                                                  !  gather, false for nonblocking send
   integer :: modc_mxreq_gather = 64       ! maximum number of nonblocking communication requests to allow
                                                  !  when using mp_swapirr and point-to-point communications for
                                                  !  gather
                                                  !  < 0 implies no limits
   integer :: modc_sw_scatter = 0          ! mod_comm irregular underlying communication method for scatter
                                                  !  0 for original mp_sendirr/mp_recvirr
                                                  !  1 for mp_swapirr and point-to-point communications
                                                  !  2 for mp_swapirr and all-to-all communications
   logical :: modc_hs_scatter = .false.    ! true for mod_comm irregular communication handshaking for scatter
   logical :: modc_send_scatter = .true.   ! true for mod_comm irregular communication blocking send for
                                                  !  scatter, false for nonblocking send
   integer :: modc_mxreq_scatter = -1      ! maximum number of nonblocking communication requests to allow
                                                  !  when using mp_swapirr and point-to-point communications for
                                                  !  scatter
                                                  !  < 0 implies no limits
   integer :: modc_sw_tracer = 0           ! mod_comm irregular underlying communication method for multiple tracers
                                                  !  0 for original mp_sendirr/mp_recvirr
                                                  !  1 for mp_swapirr and point-to-point communications
                                                  !  2 for mp_swapirr and all-to-all communications
   logical :: modc_hs_tracer = .true.      ! true for mod_comm irregular communication handshaking for multiple tracers
   logical :: modc_send_tracer = .true.    ! true for mod_comm irregular communication blocking send for
                                                  !  multiple tracers, false for nonblocking send
   integer :: modc_mxreq_tracer = -1       ! maximum number of nonblocking communication requests to allow
                                                  !  when using mp_swapirr and point-to-point communications for
                                                  !  multiple tracers
                                                  !  < 0 implies no limits
   integer :: modc_onetwo = 2              !one or two simultaneous mod_comm irregular communications 
!                                                  (excl. tracers)
   integer :: modc_tracers = 3             ! max number of tracers for simultaneous mod_comm irregular communications 
                                                  !  0 for original mp_sendirr/mp_recvirr communications
                                                  !  positive for special tracer routines


   logical :: local_dp_map=.false.    ! flag indicates that mapping between dynamics 
                                      !  and physics decompositions does not require 
                                      !  interprocess communication
   integer :: block_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in dynamics decomposition (including level 0)
   integer :: chunk_buf_nrecs         ! number of local grid points (lon,lat,lev)
                                      !  in physics decomposition (including level 0)

   integer :: proc(plat)              ! processor id associated with a given lat.
   integer, allocatable :: cut(:,:)   ! partition for MPI tasks
   integer, allocatable :: nlat_p(:)  ! number of latitudes per subdomain

   integer comm_y            ! communicator in latitude
   integer comm_z            ! communicator in vertical
   integer commxy_x          ! communicator in longitude (xy second. decomp.)
   integer commxy_y          ! communicator in latitude (xy second. decomp.)
   integer mpicom_yz         ! communicator for yz decomposition
   integer mpicom_nyz        ! communicator for multiple yz decomposition
   integer mpicom_xy         ! communicator for xy decomposition
   integer npes_yz           ! number of processes for yz decomposition
   integer npes_xy           ! number of processes for xy decomposition
   integer, allocatable :: lonrangexy(:,:)   ! global xy-longitude subdomain index
   integer, allocatable :: latrangexy(:,:)   ! global xy-latitude subdomain index
   type (ghosttype), save  :: ghostpe_yz, ghostpe1_yz
   type (parpatterntype)   :: ikj_xy_to_yz, ijk_yz_to_xy, ijk_xy_to_yz, &
                              pexy_to_pe, pkxy_to_pkc
!
! !DESCRIPTION: 
!   {\bf Purpose:} Subroutines to initialize SPMD implementation of CAM
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Alterations for LR SPMD mode
!   01.05.09  Mirin              2-D yz decomposition
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.12.20  Sawyer             Changed index order of Q3 decomposition
!   02.12.11  Sawyer             Use parbegin/endtransfer for transposes
!   03.05.07  Sawyer             Removed unneeded decompositions
!   06.03.01  Sawyer             Removed tracertrans-related variables
!
!EOP
!-----------------------------------------------------------------------

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

! !DESCRIPTION: Read in FV-specific namelist variables.  Must be 
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
   !----------------------------------------------------------------------
   integer color, ierror, ntemp

   namelist /spmd_fv_inparm/ npr_yz, geopktrans,        &
        geopkblocks,                                &
        force_2d, modcomm_transpose,             &
        modcomm_geopk, modcomm_gatscat,          &
        modc_sw_dynrun, modc_hs_dynrun,          &
        modc_send_dynrun, modc_mxreq_dynrun,     &
        modc_sw_cdcore, modc_hs_cdcore,          &
        modc_send_cdcore, modc_mxreq_cdcore,     &
        modc_sw_gather, modc_hs_gather,          &
        modc_send_gather, modc_mxreq_gather,     &
        modc_sw_scatter, modc_hs_scatter,        &
        modc_send_scatter, modc_mxreq_scatter,   &
        modc_sw_tracer, modc_hs_tracer,          &
        modc_send_tracer, modc_mxreq_tracer,     &
        modc_onetwo, modc_tracers

   npr_yz(1) = npes
   npr_yz(2) = 1
   npr_yz(3) = 1
   npr_yz(4) = npes

   if (masterproc) then
      write(iulog,*) 'Read in spmd_fv_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for spmd_fv_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'spmd_fv_inparm', status=ierr)
      if (ierr == 0) then  ! found spmd_fv_inparm
         read(unitn, spmd_fv_inparm, iostat=ierr)  ! read the spmd_fv_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                 ' error condition for spmd_fv_inparm' )
         end if
      end if
      close( unitn )
      call freeunit( unitn )
   endif

   call mpibcast (npr_yz            ,4,mpiint,0,mpicom)
   call mpibcast (geopktrans        ,1,mpiint,0,mpicom)
   call mpibcast (geopkblocks       ,1,mpiint,0,mpicom)
   call mpibcast (force_2d          ,1,mpiint,0,mpicom)
   call mpibcast (modcomm_transpose ,1,mpiint,0,mpicom)
   call mpibcast (modcomm_geopk     ,1,mpiint,0,mpicom)
   call mpibcast (modcomm_gatscat   ,1,mpiint,0,mpicom)
   call mpibcast (modc_sw_dynrun    ,1,mpiint,0,mpicom)
   call mpibcast (modc_hs_dynrun    ,1,mpilog,0,mpicom)
   call mpibcast (modc_send_dynrun  ,1,mpilog,0,mpicom)
   call mpibcast (modc_mxreq_dynrun ,1,mpiint,0,mpicom)
   call mpibcast (modc_sw_cdcore    ,1,mpiint,0,mpicom)
   call mpibcast (modc_hs_cdcore    ,1,mpilog,0,mpicom)
   call mpibcast (modc_send_cdcore  ,1,mpilog,0,mpicom)
   call mpibcast (modc_mxreq_cdcore ,1,mpiint,0,mpicom)
   call mpibcast (modc_sw_gather    ,1,mpiint,0,mpicom)
   call mpibcast (modc_hs_gather    ,1,mpilog,0,mpicom)
   call mpibcast (modc_send_gather  ,1,mpilog,0,mpicom)
   call mpibcast (modc_mxreq_gather ,1,mpiint,0,mpicom)
   call mpibcast (modc_sw_scatter   ,1,mpiint,0,mpicom)
   call mpibcast (modc_hs_scatter   ,1,mpilog,0,mpicom)
   call mpibcast (modc_send_scatter ,1,mpilog,0,mpicom)
   call mpibcast (modc_mxreq_scatter,1,mpiint,0,mpicom)
   call mpibcast (modc_sw_tracer    ,1,mpiint,0,mpicom)
   call mpibcast (modc_hs_tracer    ,1,mpilog,0,mpicom)
   call mpibcast (modc_send_tracer  ,1,mpilog,0,mpicom)
   call mpibcast (modc_mxreq_tracer ,1,mpiint,0,mpicom)
   call mpibcast (modc_onetwo       ,1,mpiint,0,mpicom)
   call mpibcast (modc_tracers      ,1,mpiint,0,mpicom)

   if (npr_yz(1) == npes .and. npr_yz(2) == 1 .and. npr_yz(3) == 1 .and. npr_yz(4) == npes) then
      npr_y   = npes
      npr_z   = 1
      nprxy_x = 1
      nprxy_y = npes
      if (masterproc) then
         write(iulog,*) 'WARNING : npr_yz not present - using 1-D domain decomposition'
      endif
      npes_yz = npes
      npes_xy = npes
   else
      npr_y   = npr_yz(1)
      npr_z   = npr_yz(2)
      nprxy_x = npr_yz(3)
      nprxy_y = npr_yz(4)
      npes_yz = npr_y*npr_z
      npes_xy = nprxy_x*nprxy_y
      if (masterproc) then
         write(iulog,*) 'npr_y = ', npr_y, '  npr_z = ', npr_z
         write(iulog,*) 'nprxy_x = ', nprxy_x, '  nprxy_y = ', nprxy_y
         write(iulog,*) 'npes = ', npes, '  npes_yz= ', npes_yz, '  npes_xy = ', npes_xy
      endif
      if (npes_yz > npes) then
         call endrun ('SPMD_DYN_SET : incorrect yz domain decomposition - aborting')
      endif
      if (npes_xy > npes) then
         call endrun ('SPMD_DYN_SET : incorrect xy domain decomposition - aborting')
      endif
      if (npes_xy < npes) then
         if (masterproc) then
            write(iulog,*) 'WARNING - proceeding with auxiliary dynamics processes'
         endif
      endif
      if (npes_yz < npes_xy) then
         if (masterproc) then
            write(iulog,*) 'WARNING - proceeding with smaller yz decomposition'
         endif
      endif
   endif
   if (ct_overlap .ne. 0) then
      if (npes .lt. 2*npes_yz) then
         call endrun ('SPMD_READNL: Not enough processes to overlap cd_core and trac2d')
      else
         if (masterproc) then
            write(iulog,*) 'Overlapping tracer and dynamics subcycles'
         endif
      endif
   endif
   if (trac_decomp .le. 0) then
      call endrun ('SPMDINIT_READNL: trac_decomp improperly initialized')
   endif
   if (npes .lt. trac_decomp*npes_yz) then
      call endrun ('SPMDINIT_READNL: Not enough processes to decompose tracers ')
   else
      if (masterproc) then
         write(iulog,*) 'Decomposing tracers into ', trac_decomp, ' groups'
      endif
   endif
   if (ct_overlap .gt. 0 .and. trac_decomp .gt. 1) then
      call endrun ('SPMDINIT_READNL: Cannot simultaneously overlap cd_core/trac2d and decompose tracers')
   endif
   myid_z   = iam/npr_y
   myid_y   = iam - myid_z*npr_y
   color = iam/npes_yz
   call mpi_comm_split(mpicom, color, iam, mpicom_yz, ierror)
   if (ierror /= mpi_success) then
      write(iulog,*) 'SPMD_DYN_READNL:  ERROR:  mpi_comm_split_yz failed with IER=', ierror
      call endrun
   endif
   call mpi_comm_size(mpicom_yz, ntemp, ierror)
   if (masterproc .and. ntemp .ne. npes_yz) then
      write(iulog,*) 'SPMD_DYN_READNL:  ERROR:  mpicom_yz has incorrect size of ', ntemp
   endif
   if (ct_overlap .gt. 0 .or. trac_decomp .gt. 1) then
      ! These are mutually exclusive options
      if ((ct_overlap .gt. 0 .and. iam .lt. 2*npes_yz) .or.         &
           (trac_decomp .gt. 1 .and. iam .lt. trac_decomp*npes_yz)) then
         color = 1
      else
         color = 0
      endif
      call mpi_comm_split(mpicom, color, iam, mpicom_nyz, ierror)
      if (ierror /= mpi_success) then
         write (iulog,*) 'SPMD_DYN_READNL:  ERROR:  mpi_comm_split_nyz failed with IER=', ierror
         call endrun
      endif
   else
      mpicom_nyz = mpicom_yz
   endif
   myidxy_y = iam/nprxy_x
   myidxy_x = iam - myidxy_y*nprxy_x
   color = iam/npes_xy
   call mpi_comm_split(mpicom, color, iam, mpicom_xy, ierror)
   if (ierror /= mpi_success) then
      write(iulog,*) 'SPMD_DYN_READNL:  ERROR:  mpi_comm_split_xy failed with IER=', ierror
      call endrun
   endif
   call mpi_comm_size(mpicom_xy, ntemp, ierror)
   if (ntemp .ne. npes_xy) then
      write(iulog,*) 'SPMD_DYN_READNL:  ERROR:  mpicom_xy has incorrect size of ', ntemp
   endif

   geopkdist   = .false.
   geopk16byte = .false.
   if (geopktrans .ne. 0) geopkdist   = .true.
   if (geopktrans .eq. 1) geopk16byte = .true.
#ifdef NO_CRAY_POINTERS
   if (geopk16byte) then
      call endrun ('SPMD_DYN_SET : cannot use geopk16 unless compiler supports cray pointers')
   end if
#endif
   if (masterproc) then
      write(iulog,*) 'non-transpose geopk communication method = ', geopkdist
      write(iulog,*) 'Z-parallel non-transpose geopk communication method = ', geopk16byte
   endif

   geopkblocks = max(1,geopkblocks)
   if ((masterproc) .and. (geopkdist) .and. (.not. geopk16byte)) then
      write(iulog,*) 'number of stages in Z-serial non-transpose geopk method = ', geopkblocks
   endif

   twod_decomp = 1

   if (npr_z .eq. 1 .and. nprxy_x .eq. 1 .and. force_2d .eq. 0) then
      twod_decomp = 0
      if (masterproc) then
         write(iulog,*) 'decomposition is effectively 1D - skipping transposes'
      endif
   else
      if (masterproc) then
         write(iulog,*) 'using multi-2d decomposition methodology'
      endif
   endif

   if (masterproc) then
      write(iulog,*) 'modcomm transpose method = ', mod_transpose
   endif

   if (masterproc) then
      write(iulog,*) 'modcomm geopk method = ', mod_geopk
   endif


   if (masterproc) then
      write(iulog,*) 'modcomm gatscat method = ', mod_gatscat
   endif

   if (masterproc) then
      write(iulog,*) 'modc_sw_dynrun = ', modc_sw_dynrun
   endif
   if (modc_sw_dynrun .lt. 0 .or. modc_sw_dynrun .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_dynrun')
   endif
   if (modc_sw_dynrun .gt. 0 .and. mod_transpose .gt. 0) then
      modc_sw_dynrun = 0
      if (masterproc) then
         write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_dynrun reset to 0 for consistency'
      endif
   endif


   if (masterproc) then
      write(iulog,*) 'modc_hs_dynrun = ', modc_hs_dynrun
   endif

   if (masterproc) then
      write(iulog,*) 'modc_send_dynrun = ', modc_send_dynrun
   endif
   if (masterproc) then
      write(iulog,*) 'modc_mxreq_dynrun = ', modc_mxreq_dynrun
   endif
   if (masterproc) then
      write(iulog,*) 'modc_sw_cdcore = ', modc_sw_cdcore
   endif
   if (modc_sw_cdcore .lt. 0 .or. modc_sw_cdcore .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_cdcore')
   endif
   if (modc_sw_cdcore .gt. 0 .and. (mod_transpose .gt. 0 .or. (mod_geopk .gt. 0 .and. geopk16byte))) then
      modc_sw_cdcore = 0
      if (masterproc) then
         write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_cdcore reset to 0 for consistency'
      endif
   endif
   if (masterproc) then
      write(iulog,*) 'modc_hs_cdcore = ', modc_hs_cdcore
   endif
   if (masterproc) then
      write(iulog,*) 'modc_send_cdcore = ', modc_send_cdcore
   endif
   if (masterproc) then
      write(iulog,*) 'modc_mxreq_cdcore = ', modc_mxreq_cdcore
   endif
   if (masterproc) then
      write(iulog,*) 'modc_sw_gather = ', modc_sw_gather
   endif
   if (modc_sw_gather .lt. 0 .or. modc_sw_gather .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_gather')
   endif
   if (modc_sw_gather .gt. 0 .and. mod_gatscat .gt. 0) then
      modc_sw_gather = 0
      if (masterproc) then
         write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_gather reset to 0 for consistency'
      endif
   endif
   if (masterproc) then
      write(iulog,*) 'modc_hs_gather = ', modc_hs_gather
   endif
   if (masterproc) then
      write(iulog,*) 'modc_send_gather = ', modc_send_gather
   endif

   if (masterproc) then
      write(iulog,*) 'modc_mxreq_gather = ', modc_mxreq_gather
   endif
   if (masterproc) then
      write(iulog,*) 'modc_sw_scatter = ', modc_sw_scatter
   endif
   if (modc_sw_scatter .lt. 0 .or. modc_sw_scatter .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_scatter')
   endif
   if (modc_sw_scatter .gt. 0 .and. mod_gatscat .gt. 0) then
      modc_sw_scatter = 0
      if (masterproc) then
         write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_scatter reset to 0 for consistency'
      endif
   endif
   if (masterproc) then
      write(iulog,*) 'modc_hs_scatter = ', modc_hs_scatter
   endif
   if (masterproc) then
      write(iulog,*) 'modc_send_scatter = ', modc_send_scatter
   endif
   if (masterproc) then
      write(iulog,*) 'modc_mxreq_scatter = ', modc_mxreq_scatter
   endif
   if (masterproc) then
      write(iulog,*) 'modc_sw_tracer = ', modc_sw_tracer
   endif
   if (modc_sw_tracer .lt. 0 .or. modc_sw_tracer .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_sw_tracer')
   endif
   if (modc_sw_tracer .gt. 0 .and. mod_transpose .gt. 0) then
      modc_sw_tracer = 0
      if (masterproc) then
         write (iulog,*) 'WARNING (SPMD_DYN_SET) - modc_sw_tracer reset to 0 for consistency'
      endif
   endif
   if (masterproc) then
      write(iulog,*) 'modc_hs_tracer = ', modc_hs_tracer
   endif
   if (masterproc) then
      write(iulog,*) 'modc_send_tracer = ', modc_send_tracer
   endif
   if (masterproc) then
      write(iulog,*) 'modc_mxreq_tracer = ', modc_mxreq_tracer
   endif
   if (masterproc) then
      write(iulog,*) 'modc_onetwo = ', modc_onetwo
   endif
   if (modc_onetwo .lt. 1 .or. modc_onetwo .gt. 2) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_onetwo')
   endif
   if (masterproc) then
      write(iulog,*) 'modc_tracers = ', modc_tracers
   endif
   if (modc_tracers .lt. 0) then
      call endrun ('SPMD_DYN_SET : inadmissable value of modc_tracers')
   endif

 end subroutine spmd_readnl

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: spmdinit_dyn --- SPMD initialization for dynamics
!
! !INTERFACE:

   subroutine spmdinit_dyn ()

! !USES:
      use parutilitiesmodule, only : parinit, parsplit
      use decompmodule, only : decompcreate

! !DESCRIPTION:
!
!   SPMD initialization routine: get number of cpus, processes, tids, etc
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Added LR-specific initialization
!   01.03.26  Sawyer             Added ProTeX documentation
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.10.16  Sawyer             Added Y at each Z decompositions
!   03.07.22  Sawyer             Removed decomps used by highp2
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer procid    ! processor id
      integer procids   ! processor id SH
      integer procidn   ! processor id NH
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer workleft  ! amount of work still to be parcelled out
      integer actual    ! actual amount of work parcelled out
      integer ideal     ! ideal amt of work to parcel out
      integer pesleft   ! number of procs still to be given work
      integer isum      ! running total of work parcelled out
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer xdist(1)  ! number of lons per subdomain
      integer, allocatable :: ydist(:) ! number of lats per subdomain
      integer, allocatable :: zdist(:) ! number of levels per subdomain
      integer, allocatable :: zdistq(:) ! number of levels per subdomain for Q3
      integer ier       ! error flag
      integer rank_y, size_y   !  rank and size wrt y-communicator
      integer rank_z, size_z   !  rank and size wrt z-communicator
      integer rankxy_x, sizexy_x   !  rank and size wrt xy x-communicator
      integer rankxy_y, sizexy_y   !  rank and size wrt xy y-communicator
      integer zdist1(1) ! used for misc. decomposition definitions
      integer, allocatable :: xdistxy(:) ! number of xy-longs per subdomain
      integer, allocatable :: ydistxy(:) ! number of xy-lats per subdomain
      integer, allocatable :: ydistqxy(:) ! number of xy tracer/lats per subdomain
      integer zdistxy(1)  ! number of xy-verts per subdomain
      integer j, k, vert, lonn
      integer ydistk(1)
      integer mod_maxirr

      spmd_on = 1

! Default 2D decomposition
      beglev = 1
      endlev = plev
      endlevp1 = plev + 1
      endlevp = plev + 1
      mod_maxirr = max(modc_onetwo, modc_tracers)
!
! Addition for LR dynamical core to initialize PILGRIM library
!
      call parinit(comm=mpicom, &
                   npryzxy = (/ npr_y, npr_z, nprxy_x, nprxy_y /), &
                   mod_method  = mod_transpose, &
                   mod_geopk   = mod_geopk,     &
                   mod_maxirr  = mod_maxirr,    &
                   mod_gatscat = mod_gatscat )
!
! Form separate communicators
!
      call parsplit(mpicom, myid_z, iam, comm_y, rank_y, size_y)
      call parsplit(mpicom, myid_y, iam, comm_z, rank_z, size_z)
      call parsplit(mpicom, myidxy_y, iam, commxy_x, rankxy_x, sizexy_x)
      call parsplit(mpicom, myidxy_x, iam, commxy_y, rankxy_y, sizexy_y)

!
!-----------------------------------------------------------------------
!
! Compute y decomposition
!
      allocate (ydist  (npr_y))
      allocate (nlat_p (0:npes-1))
      allocate (cut    (2,0:npes-1))

      ydist(:) = 0
      nlat_p(:) = 0
      cut(1,:) = -1
      cut(2,:) = -2

      lat = plat / npr_y
      workleft = plat - lat * npr_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 latitudes per subdomain')
      endif
!
! Be careful:  ydist is 1-based.  NCARs arrays, e.g., cut,  are 0-based
!
      do procid=1,npr_y
         ydist(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_y
            ydist(procids) = ydist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydist(procidn) = ydist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydist) /= plat ) then
         write(iulog,*)'SPMDINIT_DYN:', ydist,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(y) not zero.  Value is ',workleft
         call endrun
      end if

! Set the NCAR data structures

      lat  = 0
      do procid=0,npr_y-1
         cut(1,procid) = lat+1
         lat = lat + ydist(procid+1)
         cut(2,procid) = lat
         nlat_p(procid) = ydist(procid+1)

         if (masterproc) then
            write(iulog,*) 'nlat_p(',procid,') = ', nlat_p(procid)
         end if

         if (myid_y == procid) then
            beglat  = cut(1,myid_y)
            endlat  = cut(2,myid_y)
            numlats = ydist(procid+1)
         end if
      enddo

      do k = 1, npr_z-1
         do j = 0, npr_y-1
            procid = j + k*npr_y
            cut(1,procid) = cut(1,j)
            cut(2,procid) = cut(2,j)
            nlat_p(procid) = nlat_p(j)
         enddo
      enddo
!
! Compute z decomposition
!
      allocate (zdist ((npes-1)/npr_y+1))
      allocate (zdistq(npr_z))

      zdist(:) = 0

      vert = plev / npr_z
      workleft = plev - vert * npr_z
      if ( vert < 1 ) then
         call endrun ('SPMDINIT_DYN: less than 1 verticals per subdomain')
      endif

      do procid=1,npr_z
         zdist(procid) = vert
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_z+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_z
            zdist(procids) = zdist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               zdist(procidn) = zdist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(zdist) /= plev ) then
         write(iulog,*)'SPMDINIT_DYN:', zdist,' does not add up to ', plev
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(z) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglev = 1
      endlev = zdist(1)
      do procid = 1, myid_z
         beglev = endlev + 1
         endlev = beglev + zdist(procid+1) - 1
      enddo
      endlevp1 = endlev + 1
      endlevp = endlev
      if (myid_z == npr_z-1) endlevp = endlev + 1

      if (iam .ge. npes_yz) then
! Auxiliary processes only
        beglat = 1
        endlat = 0
        numlats = 0
        beglev = 1
        endlev = 0
        endlevp = endlev + 1
        endlevp1 = endlev + 1
      endif

!
! Compute x secondary decomposition
!
      allocate (xdistxy (nprxy_x))

      xdistxy(:) = 0

      lonn = plon / nprxy_x
      workleft = plon - lonn * nprxy_x
      if ( lonn < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-longitudes per subdomain')
      endif

      do procid=1,nprxy_x
         xdistxy(procid) = lonn
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_x+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_x
            xdistxy(procids) = xdistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               xdistxy(procidn) = xdistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(xdistxy) /= plon ) then
         write(iulog,*)'SPMDINIT_DYN:', xdistxy,' does not add up to ', plon
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(xy-x) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglonxy = 1
      endlonxy = xdistxy(1)
      do procid = 1, myidxy_x
         beglonxy = endlonxy + 1
         endlonxy = beglonxy + xdistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (lonrangexy(2,nprxy_x))
      lonrangexy(1,1) = 1
      lonrangexy(2,1) = xdistxy(1)
      do procid = 2, nprxy_x
         lonrangexy(1,procid) = lonrangexy(2,procid-1) + 1
         lonrangexy(2,procid) = lonrangexy(1,procid) + xdistxy(procid) - 1
      enddo
!
! Compute y secondary decomposition
!
      allocate (ydistxy ((npes-1)/nprxy_x+1))

      ydistxy(:) = 0

      lat = plat / nprxy_y
      workleft = plat - lat * nprxy_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-latitudes per subdomain')
      endif

      do procid=1,nprxy_y
         ydistxy(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_y
            ydistxy(procids) = ydistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydistxy(procidn) = ydistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydistxy) /= plat ) then
         write(iulog,*)'SPMDINIT_DYN:', ydistxy,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(iulog,*)'SPMDINIT_DYN: Workleft(xy-y) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglatxy = 1
      endlatxy = ydistxy(1)
      do procid = 1, myidxy_y
         beglatxy = endlatxy + 1
         endlatxy = beglatxy + ydistxy(procid+1) - 1
      enddo

      if (iam .ge. npes_xy) then
! Auxiliary processes only
        beglonxy = 1
        endlonxy = 0
        beglatxy = 1
        endlatxy = 0
      endif

! Compute global table

      allocate (latrangexy(2,nprxy_y))
      latrangexy(1,1) = 1
      latrangexy(2,1) = ydistxy(1)
      do procid = 2, nprxy_y
         latrangexy(1,procid) = latrangexy(2,procid-1) + 1
         latrangexy(2,procid) = latrangexy(1,procid) + ydistxy(procid) - 1
      enddo

!
! Do generic NCAR decomposition
!
      proc(:) = 0
      do procid=0,npr_y*npr_z-1
         if (iam == 0) then
            write(iulog,*)'procid ',procid,' assigned ', &
                 cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                 cut(1,procid),' through ',cut(2,procid)
         endif
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
      end do

      nmostlat = plat
      smostlat = 1
      if (iam .lt. npes_yz) then

! Primary processes only
!
! Number of neighbor processors needed for boundary communication.  North
! first.
!
        nmostlat = 0
        isum = 0
        do procid=myid_y+1,npr_y-1
           nmostlat = cut(2,procid)
           isum = isum + cut(2,procid) - cut(1,procid) + 1
           if (isum >= numbnd) goto 20
        end do
20      if (myid_y /= npr_y-1 .and. isum < numbnd .and. nmostlat /= plat)then
           call endrun ('SPMDINIT_DYN: Something wrong in computation of northern neighbors')
        end if

        smostlat = 0
        isum = 0
        do procid=myid_y-1,0,-1
           smostlat = cut(1,procid)
           isum = isum + cut(2,procid) - cut(1,procid) + 1
           if (isum >= numbnd) goto 30
        end do
30      if (myid_y /= 0 .and. isum < numbnd .and. smostlat /= 1) then
           call endrun ('SPMDINIT_DYN: Something wrong in computation of southern neighbors')
        end if

!        write(iulog,*)'-----------------------------------------'
!        write(iulog,*)'Number of lats passed north & south = ',numbnd
!        write(iulog,*)'Node  Partition'
!        write(iulog,*)'-----------------------------------------'
!        do procid=0,npes-1
!           write(iulog,200) procid,cut(1,procid),cut(2,procid)
!        end do
!        write(iulog,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!        write(iulog,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      endif

      deallocate (ydist)
      deallocate (zdist)

      return
!
! Formats
!
200   format(i3,4x,i3,'-',i3,7x,i3,'-',i3)

!EOC
   end subroutine spmdinit_dyn

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      call endrun ('decomp_wavenumbers() should never be called in LR dynamics')

   end subroutine decomp_wavenumbers

   subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: placeholder for buffer allocation routine 
! 
! Method: 
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      return

   end subroutine spmdbuf

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
      integer, intent(in) :: numperlat            ! number of elements per latitude
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
     
      displs(:) = 0
      do p=1,npr_y-1
         displs(p) = displs(p-1) + numperproc(p-1)
      end do

      if (npr_z > 1) then
         do p=1,npr_z-1
            displs(p*npr_y:(p+1)*npr_y-1) = displs(0:npr_y-1)
         enddo
      endif

    end subroutine compute_gsfactors

#endif

end module spmd_dyn

