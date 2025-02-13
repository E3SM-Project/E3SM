!/ ------------------------------------------------------------------- /
      MODULE WAV_COMP_MCT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     A generic cpl7 interface for WAVEWATCH III
!     using input fields from cpl7.
!
!  2. Method :
!
!     MCT component for the actual wave model (W3WAVE).
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!       NHMAX   I.P.  Maximum number of homogeneous fields.
!
!       NDSO    Int.  General output unit number (shell only).
!       NDSE    Int.  Error output unit number (shell only).
!       NDST    Int.  Test output unit number (shell only).
!       FLH     L.A.  Flags for homogeneous fields.
!       NH      I.A.  Number of times for homogeneous fields.
!       THO     I.A.  Times of homogeneous fields.
!       TIME0   I.A.  Starting time.
!       TIMEN   I.A.  Ending time.
!     ----------------------------------------------------------------
!
!       NDS, NTRACE, ..., see W3WAVE
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set nummber of data structures
!      W3SETG    Subr.   Id.    Point to data structure.
!      W3NDAT    Subr. W3WDATMD Set nummber of data structures
!      W3SETW    Subr.   Id.    Point to data structure.
!      W3NMOD    Subr. W3ADATMD Set nummber of data structures
!      W3NAUX    Subr.   Id.    Point to data structure.
!      W3NOUT    Subr. W3ODATMD Set nummber of data structures
!      W3SETO    Subr.   Id.    Point to data structure.
!      W3NINP    Subr. W3IDATMD Set nummber of data structures
!      W3SETI    Subr.   Id.    Point to data structure.
!
!      STME21    Subr. W3TIMEMD Print date and time readable.
!      DSEC21    Func.   Id.    Difference between times.
!      TICK21    Subr.   Id.    Increment time.
!
!      W3INIT    Subr. W3INITMD Wave model initialization.
!      W3WAVE    Subr. W3WAVEMD Wave model.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!     - Checks on I-O.
!     - Check on time interval.
!
!  7. Remarks :
!
!     - A rigourous input check is made in W3INIT.
!     - See W3WDAS for documentation on the set-up of the data
!       assimilation.
!
!  8. Structure :
!
!     ----------------------------------------------------------------
!
!     wav_comp_init
!
!        0.   Set up data structures.                ( W3NMOD, etc. )
!        1.   I-O setup.
!          a  For shell.
!          b  For WAVEWATCH III.
!          c  Local parameters.
!        2.   Define input fields
!        3.   Set time frame.
!        4.   Define output
!          a  Loop over types, do
!        +--------------------------------------------------------+
!        | b    Process standard line                             |
!        | c    If type 1: fields of mean wave parameters         |
!        | d    If type 2: point output                           |
!        | e    If type 3: track output                           |
!        | f    If type 4: restart files                          |
!        | g    If type 5: boundary output                        |
!        | h    If type 6: separated wave fields                  |
!        +--------------------------------------------------------+
!        5.   Initialzations
!          a  Wave model.                              ( W3INIT )
!          d  Set field times.
!
!     wav_comp_run
!
!        7.   Run model with input
!             Do until end time is reached
!        +--------------------------------------------------------+
!        | a  Determine next time interval and input fields.      |
!        |   1  Preparation                                       |
!        |      Loop over input fields                            |
!        | +------------------------------------------------------|
!        | | 2  Check if update is needed                         |
!        | | 4  Update next ending time                           |
!        | +------------------------------------------------------|
!        | b  Run wave model.                          ( W3WAVE ) |
!        | c  If requested, data assimilation.         ( W3WDAS ) |
!        | d  Final output if needed.                  ( W3WAVE ) |
!        | e  Check time                                          |
!        +--------------------------------------------------------+
!
!     wav_comp_fin
!
!     ----------------------------------------------------------------
!
!
!/ ------------------------------------------------------------------- /
      use w3gdatmd, only: dtmax, dtcfl, dtcfli, dtmin, &
                          nx, ny, nsea, nseal, mapsf, mapfs, mapsta, mapst2, x0, y0, sx, sy, xgrd, ygrd, &
                          w3nmod, w3setg, AnglD, &
                          sig, nk, zb, dmin, &
                          usspf
      use w3wdatmd, only: time, w3ndat, w3setw, wlv, va, ust, ice 
      use w3adatmd, only: ussp, w3naux, w3seta, sxx, sxy, syy, fliwnd, flcold, dw, cg, wn, hs, fp0, thp0
      use w3idatmd, only: inflags1, inflags2,w3seti, w3ninp
      USE W3IDATMD, ONLY: TC0, CX0, CY0, TCN, CXN, CYN, ICEP1, ICEP5, TI1, TI5
      USE W3IDATMD, ONLY: TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN
      USE W3IDATMD, ONLY: TIN, ICEI
      USE W3IDATMD, ONLY: TLN, WLEV
      use w3odatmd, only: w3nout, w3seto, naproc, iaproc, napout, naperr, ndse, ndso, &
                          nogrp, ngrpp, noge, idout, fnmpre, iostyp, notype, flout, &
                          fnmpre, ifile4, ofiles
      use w3servmd, only: w3xyrtn
      use wmscrpmd, only: grid_area
      use w3parall, only: init_get_isea
      use w3dispmd, only: wavnu1
      use w3triamd, only: SETUGIOBP
!/
      use w3initmd, only: w3init 
      use w3wavemd, only: w3wave
      use w3gridmd, only: w3grid
!/
      use w3iopomd, only:
      use w3iorsmd, only: w3iors
      use w3iogomd, only: w3flgrdflag
      use w3timemd, only: stme21 
      use w3cesmmd, only : casename, initfile, runtype
      use w3cesmmd, only : inst_index, inst_name, inst_suffix

      use esmf
      use mct_mod
      use seq_flds_mod

      use ww3_cpl_indices  , only : ww3_cpl_indices_set
      use ww3_cpl_indices  , only : index_x2w_Sa_u, index_x2w_Sa_v, index_x2w_Sa_tbot, index_x2w_Si_ifrac, index_x2w_si_ithick
      use ww3_cpl_indices  , only : index_x2w_So_t, index_x2w_So_u, index_x2w_So_v, index_x2w_So_bldepth, index_x2w_So_ssh
      use ww3_cpl_indices  , only : index_w2x_Sw_ustokes_wavenumber_1, index_w2x_Sw_vstokes_wavenumber_1, &
                                    index_w2x_Sw_ustokes_wavenumber_2, index_w2x_Sw_vstokes_wavenumber_2, &
                                    index_w2x_Sw_ustokes_wavenumber_3, index_w2x_Sw_vstokes_wavenumber_3, &
                                    index_w2x_Sw_ustokes_wavenumber_4, index_w2x_Sw_vstokes_wavenumber_4, &
                                    index_w2x_Sw_ustokes_wavenumber_5, index_w2x_Sw_vstokes_wavenumber_5, &
                                    index_w2x_Sw_ustokes_wavenumber_6, index_w2x_Sw_vstokes_wavenumber_6, &
                                    index_w2x_Sw_Hs, index_w2x_Sw_Fp, index_w2x_Sw_Dp


      use shr_sys_mod      , only : shr_sys_flush, shr_sys_abort
      use shr_kind_mod     , only : in=>shr_kind_in, r8=>shr_kind_r8, &
                                    cs=>shr_kind_cs, cl=>shr_kind_cl
      use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
      use seq_timemgr_mod  , only : seq_timemgr_eclockgetdata, &
                                    seq_timemgr_RestartAlarmIsOn, &
                                    seq_timemgr_historyAlarmIsOn
      use seq_infodata_mod , only : seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata, &
                                    seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                    seq_infodata_start_type_brnch
      use seq_comm_mct     , only : seq_comm_inst, seq_comm_name, seq_comm_suffix
      use shr_file_mod     , only : shr_file_setlogunit, shr_file_setloglevel, &
                                    shr_file_getlogunit, shr_file_getloglevel, &
                                    shr_file_getunit, shr_file_freeunit, shr_file_setio
      use shr_nl_mod       , only : shr_nl_find_group_name
      use shr_mpi_mod      , only : shr_mpi_bcast

!
      implicit none
!
      public :: wav_init_mct
      public :: wav_run_mct
      public :: wav_final_mct

      private

      private :: wav_setgsmap_mct
      private :: wav_domain_mct

      integer,target,save :: stdout
      integer,save :: odat(40)
      integer,save :: nds(13)
      integer,save :: mds(5)

      real,allocatable,save :: AnglDL(:)

      include "mpif.h"
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CONTAINS
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SUBROUTINE WAV_INIT_MCT( EClock, cdata, x2w_w, w2x_w, NLFilename )

      !/ ------------------------------------------------------------------- /
      !/ Parameter list
      !/

      implicit none
      type(esmf_clock), intent(inout) :: eclock
      type(seq_cdata) , intent(inout) :: cdata
      type(mct_avect) , intent(inout) :: x2w_w, w2x_w
      character(len=*), intent(in), optional :: nlfilename ! namelist filename

      !/ ------------------------------------------------------------------- /
      !/ Local PARAMETER statements
      !/

      INTEGER, PARAMETER  :: NHMAX =    200

      !/ ------------------------------------------------------------------- /
      !/ Local parameters
      !/

      integer :: n
      integer :: compid
      integer :: mpi_comm
      integer :: lsize
      integer :: shrlogunit, shrloglev
      integer :: hh,mm,ss
      integer :: dtime_sync        ! integer timestep size
      integer :: start_ymd         ! start date (yyyymmdd)
      integer :: start_tod         ! start time of day (sec)
      integer :: stop_ymd          ! stop date (yyyymmdd)
      integer :: stop_tod          ! stop time of day (sec)
      integer :: ix, iy
      integer :: is, ik
      integer :: gindex

      character(CL)            :: starttype
      type(mct_gsmap), pointer :: gsmap
      type(mct_ggrid), pointer :: dom
      type(seq_infodata_type), pointer :: infodata   ! input init object
      type(mct_aVect) :: x2w0

      integer             :: unitn            ! namelist unit number
      integer             :: ntrace(2), time0(2)
      integer             :: nu
      integer             :: timen(2), nh(4), iprt(6)
      integer             :: i,j,npts
      integer             :: ierr
      integer             :: jsea,isea
      integer             :: pnt_out_freq, grd_out_freq
      integer             :: iproc
      integer             :: imod
      integer             :: ndno
      integer, allocatable:: maptst(:,:)
      real                :: a(nhmax,4)
      real, allocatable   :: x(:), y(:)
      logical             :: flgrd(nogrp,ngrpp), flgrd2(nogrp,ngrpp)
      logical             :: flgd(nogrp), flgd2(nogrp)
      logical             :: prtfrm, flt
      logical             :: exists
      logical             :: IsMulti
      logical             :: flagstidein(4)
      character(len=*),parameter :: subname = '(wav_init_mct)'
      real                :: wlveff
      real                :: depth

      character(len=3)    :: fext
      character(len=3)    :: idtst
      character(len=10)   :: pn
      character(len=15)   :: restart_timestamp
      character(len=20)   :: strng
      character(len=23)   :: dtme21
      character(len=40), allocatable :: pnames(:)
      character(len=256) :: stafile
      character(len=1024) :: fldout, fldcou

      !/ ------------------------------------------------------------------- /

      namelist /ww3_inparm/ stafile, fldout, fldcou, pnt_out_freq, grd_out_freq


      !--------------------------------------------------------------------
      ! Initialize mpi
      !--------------------------------------------------------------------

      call seq_cdata_setptrs(cdata, id=compid, mpicom=mpi_comm, &
           gsmap=gsmap, dom=dom, infodata=infodata)
      inst_name   = seq_comm_name(compid)
      inst_index  = seq_comm_inst(compid)
      inst_suffix = seq_comm_suffix(compid)
      call ww3_cpl_indices_set()

      call mpi_barrier ( mpi_comm, ierr )
      call mpi_comm_rank(mpi_comm, iproc, ierr)
      iproc = iproc + 1

      if (iproc .eq. 1) then
         call shr_file_setio('wav_modelio.nml'//trim(inst_suffix),stdout)
      endif

      ! Redirect share output to wav log
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit(stdout)
    
      !--------------------------------------------------------------------
      ! Initialize WW3 grid
      !--------------------------------------------------------------------

      ndso   => stdout
      ndse   => stdout

      mds(1) = shr_file_getunit()
      mds(2) = shr_file_getunit()
      mds(3) = shr_file_getunit()
      mds(4) = shr_file_getunit()
      mds(5) = 10
      if ( iproc .eq. 1) then
        call w3grid(mds)
      endif

      do i=1,4
         call shr_file_freeUnit(mds(i))
      end do

      call mpi_barrier ( mpi_comm, ierr )

      !--------------------------------------------------------------------
      ! Set up data structures
      !--------------------------------------------------------------------
      
      if ( iproc .ne. 1) then
        call w3nmod ( 1, 6, 6 ) ! this is called for iproc = 1 in w3grid
      endif
      call w3ndat (    6, 6 )
      call w3naux (    6, 6 )
      if ( iproc .ne. 1) then
        call w3nout (    6, 6) ! this is called for iproc = 1 in w3grid
      end if
      call w3ninp (    6, 6 )

      call w3setg ( 1, 6, 6 )
      call w3setw ( 1, 6, 6 )
      call w3seta ( 1, 6, 6 )
      call w3seto ( 1, 6, 6 )
      call w3seti ( 1, 6, 6 )

      call mpi_comm_rank(mpi_comm, iaproc, ierr)
      call mpi_comm_size(mpi_comm, naproc, ierr)
      iaproc = iaproc + 1
      napout = 1
      naperr = 1

      call shr_mpi_bcast(usspf,mpi_comm)

      ndso   => stdout
      ndse   => stdout

      ofiles = 0

      !--------------------------------------------------------------------
      ! Initialize run type
      !--------------------------------------------------------------------

      call seq_infodata_GetData( infodata, start_type=starttype)

      if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
         runtype = "initial"
      else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
         runtype = "continue"
      else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
         runtype = "branch"
      end if

      !--------------------------------------------------------------------
      ! IO set-up
      !--------------------------------------------------------------------

      ndno = shr_file_getunit()

      ! log units
      nds( 1) = stdout               ! General output unit number ("log file")
      nds( 2) = stdout               ! Error output unit number.
      nds( 3) = stdout               ! Test output unit number.
      nds( 4) = stdout               ! "screen", i.e., direct output location

      ! input units
      nds( 5) = shr_file_getunit()   ! Model definition file unit number.
      nds( 9) = ndno                 ! Input boundary data file unit number. (FLOUT(3) flag)

      ! output units
      nds( 6) = shr_file_getunit()   ! Restart file unit number.
      nds( 7) = shr_file_getunit()   ! Grid output file unit number. (FLOUT(1) flag)
      nds( 8) = shr_file_getunit()   ! Point output file unit number. (FLOUT(2) flag)
      nds(10) = ndno                 ! Output boundary data file unit number (FLOUT(4) flag)
      nds(11) = ndno                 ! Track information file unit number. (FLOUT(5) flag)
      nds(12) = ndno                 ! Track output file unit number. (FLOUT(6) flag)
      nds(13) = ndno

      ntrace(1) =  nds(3)            ! Output unit number for trace.
      ntrace(2) =  10                ! Maximum number of trace prints. 

      if ( iaproc .eq. napout ) then
        write (ndso,900)
        call shr_sys_flush(ndso)
      endif

      if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
         if ( iaproc .eq. napout ) write(ndso,*) 'starttype: initial'
      else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
         if ( iaproc .eq. napout ) write(ndso,*) 'starttype: continue'
      else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
         if ( iaproc .eq. napout ) write(ndso,*) 'starttype: branch'
      end if
      if ( iaproc == napout) then
         write(ndso,*) trim(subname),' inst_name   = ',trim(inst_name)
         write(ndso,*) trim(subname),' inst_index  = ',inst_index
         write(ndso,*) trim(subname),' inst_suffix = ',trim(inst_suffix)
      endif

      !--------------------------------------------------------------------
      ! Define input fields
      !--------------------------------------------------------------------

      inflags1 = .false.

      inflags1(1)  = .true. ! water levels
      inflags1(2)  = .true. ! currents
      inflags1(3)  = .true. ! winds
      inflags1(4)  = .true. ! ice concentration
      inflags1(-7) = .true. ! ice thickness
      inflags1(-3) = .true. ! ice floe size

      inflags2 = inflags1

      !--------------------------------------------------------------------
      ! Set time frame
      !--------------------------------------------------------------------

      ! TIME0 = from ESMF clock
      ! NOTE - are not setting TIMEN here

      if ( iaproc .eq. napout ) then
        write (ndso,930)
        call shr_sys_flush(ndso)
      endif

      if ( runtype == "continue" .or. runtype == "branch") then      
        call seq_timemgr_EClockGetData(EClock, &
             curr_ymd=start_ymd, curr_tod=start_tod)
      else
        call seq_timemgr_EClockGetData(EClock, &
             start_ymd=start_ymd, start_tod=start_tod)
      endif

      hh = start_tod/3600
      mm = (start_tod - (hh * 3600))/60
      ss = start_tod - (hh*3600) - (mm*60)

      time0(1) = start_ymd
      time0(2) = hh*10000 + mm*100 + ss

      call seq_timemgr_EClockGetData(EClock, &
           stop_ymd=stop_ymd, stop_tod=stop_tod)

      hh = stop_tod/3600
      mm = (stop_tod - (hh * 3600))/60
      ss = stop_tod - (hh*3600) - (mm*60)

      timen(1) = stop_ymd
      timen(2) = hh*10000 + mm*100 + ss

      call stme21 ( time0 , dtme21 )
      if ( iaproc .eq. napout ) then 
        write (ndso,931) dtme21
        call shr_sys_flush(ndso)
      endif
      time = time0

      !--------------------------------------------------------------------
      ! Handle restart
      !--------------------------------------------------------------------

      if ( runtype == "continue" .or. runtype == "branch") then
         if (iaproc == napout) then

                inquire(file='rpointer.wav',exist=exists)
                if (.not.exists) then 
                   write(ndso,*) ' ERROR: rpointer file does not exist'
                   call shr_sys_abort(' wav ERROR: rpointer file missing')
                endif

                nu = shr_file_getUnit()
                open(unit=nu,file='rpointer.wav')
                read(nu,*) restart_timestamp
                close(nu)
                call shr_file_freeUnit(nu)

                inquire(file=trim(restart_timestamp)//'.restart.ww3',exist=exists)
                if (.not.exists) then 
                   write(ndso,*) ' ERROR: ww3 restart file does not exist'
                   call shr_sys_abort(' wav ERROR: restart file missing')
                else
                   write(ndso,*) 'Restart file: '//trim(restart_timestamp)//'.restart.ww3'//' found'
                endif
                call shr_sys_flush(ndso)
         endif 

         call shr_mpi_bcast(restart_timestamp,mpi_comm)
      endif

      !--------------------------------------------------------------------
      ! Read namelist
      !--------------------------------------------------------------------
      if ( iaproc .eq. napout ) then
         unitn = shr_file_getunit()
         write(ndso,*) 'Read in ww3_inparm namelist from ww3_in'//trim(inst_suffix)

         open( unitn, file='ww3_in'//trim(inst_suffix), status='old' )

         call shr_nl_find_group_name(unitn, 'ww3_inparm', status=ierr)
         if (ierr == 0) then
            read (unitn, ww3_inparm, iostat=ierr)
            if (ierr /= 0) then
               call shr_sys_abort('problem reading ww3_inparm namelist')
            end if
         end if
         close( unitn )
         call shr_file_freeUnit( unitn )
      end if
      call shr_mpi_bcast(stafile, mpi_comm)
      call shr_mpi_bcast(fldout, mpi_comm)
      call shr_mpi_bcast(fldcou, mpi_comm)
      call shr_mpi_bcast(pnt_out_freq, mpi_comm)
      call shr_mpi_bcast(grd_out_freq, mpi_comm)


      !--------------------------------------------------------------------
      ! Define output type and fields
      !--------------------------------------------------------------------

      iostyp = 2        ! gridded field
      if ( iaproc .eq. napout ) then
        write (ndso,940) 'no dedicated output process, any file system '
        call shr_sys_flush(ndso)
      endif

      ! ODAT    I.A.   I   Output data, five parameters per output type
      !                     1-5  Data for OTYPE = 1; gridded fields.
      !                          1 YYYMMDD for first output.
      !                          2 HHMMSS for first output.
      !                          3 Output interval in seconds.
      !                          4 YYYMMDD for last output.
      !                          5 HHMMSS for last output.
      !                     6-10 Id. for OTYPE = 2; point output.
      !                    11-15 Id. for OTYPE = 3; track point output.
      !                    16-20 Id. for OTYPE = 4; restart files.
      !                    21-25 Id. for OTYPE = 5; boundary data.
      !                    26-30 Id. for OTYPE = 6; partitioned data.
      !                    31-35 Id. for OTYPE = 7; coupling data.
      !                    36-40 Id. for OTYPE = 8; second restart file
      ! FLGRD   L.A.   I   Flags for gridded output.
      ! NPT     Int.   I   Number of output points
      ! X/YPT   R.A.   I   Coordinates of output points.
      ! PNAMES  C.A.   I   Output point names.

      notype = 8

      do j=1, notype
         odat(5*(j-1)+3) = 0
      end do

      ! get coupling interval
      call seq_timemgr_eclockgetdata(eclock, dtime=dtime_sync )

      ! Gridded fields
      odat(1) = time(1)     ! YYYYMMDD for first output
      odat(2) = time(2)     ! HHMMSS for first output
      odat(3) = grd_out_freq  ! output interval in sec ! changed by Adrean
      odat(4) = 99990101    ! YYYYMMDD for last output
      odat(5) = 0           ! HHMMSS for last output

      ! Point output
      odat(6) = time(1)     ! YYYYMMDD for first output
      odat(7) = time(2)     ! HHMMSS for first output
      odat(8) = pnt_out_freq    ! output interval in sec 
      odat(9) = 99990101    ! YYYYMMDD for last output
      odat(10) = 0          ! HHMMSS for last output

      ! Restart files
      odat(16) = time(1)    ! YYYYMMDD for first output
      odat(17) = time(2)    ! HHMMSS for first output
      odat(18) = dtime_sync ! output interval in sec (dummy for MPI initialization in w3init)
      odat(19) = 99990101   ! YYYYMMDD for last output
      odat(20) = 0          ! HHMMSS for last output

      ! Coupling data
      odat(31) = time(1)    ! YYYYMMDD for first output
      odat(32) = time(2)    ! HHMMSS for first output
      odat(33) = dtime_sync ! output interval in sec
      odat(34) = 99990101   ! YYYYMMDD for last output
      odat(35) = 0          ! HHMMSS for last output

      ! Output Type 1: fields of mean wave parameters gridded output
      flgrd2 = .false.
      flgd2  = .false.
      flgd = .false.

      !  G  G
      !  R  X Grp  Param Code     Output  Parameter/Group
      !  B  O Numb Numbr Name        Tag  Definition 
      !  --------------------------------------------------
      !        1                          Forcing Fields
      !   -------------------------------------------------
      !  T  T  1     1   DW         DPT   Water depth.
      !  T  T  1     2   C[X,Y]     CUR   Current velocity.
      !  T  T  1     3   UA         WND   Wind speed.
      !  T  T  1     4   AS         AST   Air-sea temperature difference.
      !  T  T  1     5   WLV        WLV   Water levels.
      !  T  T  1     6   ICE        ICE   Ice concentration.
      !  T  T  1     7   IBG        IBG   Iceberg-induced damping.
      !  T  T  1     8   D50        D50   Median sediment grain size.
      !  T  T  1     9   IC1        IC1   Ice thickness.
      !  T  T  1    10   IC5        IC5   Ice flow diameter.
      !   -------------------------------------------------
      !        2                          Standard mean wave Parameters
      !   -------------------------------------------------
      !  T  T  2     1   HS         HS    Wave height.
      !  T  T  2     2   WLM        LM    Mean wave length.
      !  T  T  2     3   T02        T02   Mean wave period (Tm0,2).
      !  T  T  2     4   T0M1       T0M1  Mean wave period (Tm0,-1).
      !  T  T  2     5   T01        T01   Mean wave period (Tm0,1).
      !  T  T  2     6   FP0        FP    Peak frequency.
      !  T  T  2     7   THM        DIR   Mean wave direction.
      !  T  T  2     8   THS        SPR   Mean directional spread.
      !  T  T  2     9   THP0       DP    Peak direction.
      !  T  T  2    10   HIG        HIG   Infragravity height
      !  T  T  2    11   STMAXE     MXE   Max surface elev (STE)
      !  T  T  2    12   STMAXD     MXES  St Dev of max surface elev (STE)
      !  T  T  2    13   HMAXE      MXH   Max wave height (STE)
      !  T  T  2    14   HCMAXE     MXHC  Max wave height from crest (STE)
      !  T  T  2    15   HMAXD      SDMH  St Dev of MXC (STE)
      !  T  T  2    16   HCMAXD     SDMHC St Dev of MXHC (STE)
      !  F  T  2    17   WBT        WBT   Dominant wave breaking probability bT
      !  F  F  2    18   FP0        TP    Peak period (from peak freq)
      !   -------------------------------------------------
      !        3                          Spectral Parameters (first 5)
      !   -------------------------------------------------
      !  F  F  3     1   EF         EF    Wave frequency spectrum
      !  F  F  3     2   TH1M       TH1M  Mean wave direction from a1,b2
      !  F  F  3     3   STH1M      STH1M Directional spreading from a1,b2
      !  F  F  3     4   TH2M       TH2M  Mean wave direction from a2,b2
      !  F  F  3     5   STH2M      STH2M Directional spreading from a2,b2
      !  F  F  3     6   WN         WN    Wavenumber array
      !   -------------------------------------------------
      !        4                          Spectral Partition Parameters 
      !   -------------------------------------------------
      !  T  T  4     1   PHS        PHS   Partitioned wave heights.
      !  T  T  4     2   PTP        PTP   Partitioned peak period.
      !  T  T  4     3   PLP        PLP   Partitioned peak wave length.
      !  T  T  4     4   PDIR       PDIR  Partitioned mean direction.
      !  T  T  4     5   PSI        PSPR  Partitioned mean directional spread.
      !  T  T  4     6   PWS        PWS   Partitioned wind sea fraction.
      !  T  T  4     7   PDP        PDP   Peak wave direction of partition.
      !  T  T  4     8   PQP        PQP   Goda peakdedness parameter of partition.
      !  T  T  4     9   PPE        PPE   JONSWAP peak enhancement factor of partition.
      !  T  T  4    10   PGW        PGW   Gaussian frequency width of partition.
      !  T  T  4    11   PSW        PSW   Spectral width of partition.
      !  T  T  4    12   PTM1       PTM10 Mean wave period (m-1,0) of partition.
      !  T  T  4    13   PT1        PT01  Mean wave period (m0,1) of partition.
      !  T  T  4    14   PT2        PT02  Mean wave period (m0,2) of partition.
      !  T  T  4    15   PEP        PEP   Peak spectral density of partition.
      !  T  T  4    16   PWST       TWS   Total wind sea fraction.
      !  T  T  4    17   PNR        PNR   Number of partitions.
      !   -------------------------------------------------
      !        5                          Atmosphere-waves layer
      !   -------------------------------------------------
      !  T  T  5     1   UST        UST   Friction velocity.
      !  F  T  5     2   CHARN      CHA   Charnock parameter
      !  F  T  5     3   CGE        CGE   Energy flux
      !  F  T  5     4   PHIAW      FAW   Air-sea energy flux
      !  F  T  5     5   TAUWI[X,Y] TAW   Net wave-supported stress
      !  F  T  5     6   TAUWN[X,Y] TWA   Negative part of the wave-supported stress
      !  F  F  5     7   WHITECAP   WCC   Whitecap coverage
      !  F  F  5     8   WHITECAP   WCF   Whitecap thickness
      !  F  F  5     9   WHITECAP   WCH   Mean breaking height
      !  F  F  5    10   WHITECAP   WCM   Whitecap moment
      !  F  F  5    11   FWS        FWS   Wind sea mean period
      !   -------------------------------------------------
      !        6                          Wave-ocean layer 
      !   -------------------------------------------------
      !  F  F  6     1   S[XX,YY,XY] SXY  Radiation stresses.
      !  F  F  6     2   TAUO[X,Y]  TWO   Wave to ocean momentum flux
      !  F  F  6     3   BHD        BHD   Bernoulli head (J term) 
      !  F  F  6     4   PHIOC      FOC   Wave to ocean energy flux
      !  F  F  6     5   TUS[X,Y]   TUS   Stokes transport
      !  F  F  6     6   USS[X,Y]   USS   Surface Stokes drift
      !  F  F  6     7   [PR,TP]MS  P2S   Second-order sum pressure 
      !  F  F  6     8   US3D       USF   Spectrum of surface Stokes drift
      !  F  F  6     9   P2SMS      P2L   Micro seism  source term
      !  F  F  6    10   TAUICE     TWI   Wave to sea ice stress
      !  F  F  6    11   PHICE      FIC   Wave to sea ice energy flux
      !  F  F  6    12   USSP       USP   Partitioned surface Stokes drift
      !   -------------------------------------------------
      !        7                          Wave-bottom layer 
      !   -------------------------------------------------
      !  F  F  7     1   ABA        ABR   Near bottom rms amplitides.
      !  F  F  7     2   UBA        UBR   Near bottom rms velocities.
      !  F  F  7     3   BEDFORMS   BED   Bedforms
      !  F  F  7     4   PHIBBL     FBB   Energy flux due to bottom friction 
      !  F  F  7     5   TAUBBL     TBB   Momentum flux due to bottom friction
      !   -------------------------------------------------
      !        8                          Spectrum parameters
      !   -------------------------------------------------
      !  F  F  8     1   MSS[X,Y]   MSS   Mean square slopes
      !  F  F  8     2   MSC[X,Y]   MSC   Spectral level at high frequency tail
      !  F  F  8     3   WL02[X,Y]  WL02  East/X North/Y mean wavelength compon
      !  F  F  8     4   ALPXT      AXT   Correl sea surface gradients (x,t)
      !  F  F  8     5   ALPYT      AYT   Correl sea surface gradients (y,t)
      !  F  F  8     6   ALPXY      AXY   Correl sea surface gradients (x,y)
      !   -------------------------------------------------
      !        9                          Numerical diagnostics  
      !   -------------------------------------------------
      !  T  T  9     1   DTDYN      DTD   Average time step in integration.
      !  T  T  9     2   FCUT       FC    Cut-off frequency.
      !  T  T  9     3   CFLXYMAX   CFX   Max. CFL number for spatial advection. 
      !  T  T  9     4   CFLTHMAX   CFD   Max. CFL number for theta-advection. 
      !  F  F  9     5   CFLKMAX    CFK   Max. CFL number for k-advection. 
      !   -------------------------------------------------
      !        10                         User defined          
      !   -------------------------------------------------
      !  F  F  10    1              U1    User defined #1. (requires coding ...)
      !  F  F  10    2              U2    User defined #1. (requires coding ...)
      !   -------------------------------------------------

      call w3flgrdflag(ndso,ndso,ndse,fldout,flgd, flgrd, iaproc,napout,ierr)
      call w3flgrdflag(ndso,ndso,ndse,fldcou,flgd2,flgrd2,iaproc,napout,ierr)

      call read_stations_file(ndso,stafile,npts,x,y,pnames)

      !--------------------------------------------------------------------
      ! Wave model initializations
      !--------------------------------------------------------------------

      ! Notes on ww3 initialization:
      ! ww3 read initialization occurs in w3iors (which is called by initmd)
      ! For a startup (including hybrid) or branch run the initial datafile is
      ! set in namelist input 'initfile'
      ! For a continue run - the initfile vluae is created from the time(1:2)
      ! array set below


      ! Set casename (in w3cesmmd)
      call seq_infodata_GetData(infodata,case_name=casename)

      call mpi_barrier ( mpi_comm, ierr )

      if ( iaproc .eq. napout ) then
        write (ndso,*) 'before w3init'
        call shr_sys_flush(ndso)
      endif

      ! Read in input data and initialize the model
      ! w3init calls w3iors which:
      ! - reads either the initfile if the run is startup or branch
      ! - constructs the filename from the casename variable and the time(:) array
      !   which is set above
      imod = 1
      IsMulti = .false.
      fext = 'ww3'
      prtfrm = .false. 
      iprt = 0
      flagstidein = .false.
      call w3init ( imod, IsMulti, fext, nds, ntrace, odat, flgrd, flgrd2, flgd, flgd2, npts, x, y,   &
           pnames, iprt, prtfrm, mpi_comm, flagstidein)

      if ( iaproc .eq. napout ) then
        write (ndso,*) 'after w3init'
        call shr_sys_flush(ndso)
      endif


      ! overwrite dt values with variables from coupler
      ! is this a problem with any things being set in w3init?
      dtmax  = real(dtime_sync)
      dtcfl  = real(dtime_sync) / 2. !checked by adrean
      dtcfli = real(dtime_sync)      !checked by adrean
      dtmin  = real(dtime_sync) / 12 !checked by adrean

      ! overwrite restart output values
      flout(4) = .false.

      ! Localize AnglDL 
      allocate(AnglDL(nseal))
      do jsea = 1,nseal
        isea = iaproc + (jsea-1)*naproc
        AnglDL(jsea) = AnglD(isea)
      enddo

      call mpi_barrier ( mpi_comm, ierr )

      !--------------------------------------------------------------------
      ! cpl7/mct initialization
      !--------------------------------------------------------------------

      ! initialize mct gsmap

      call wav_setgsmap_mct(mpi_comm, compid, gsmap)
      lsize = mct_gsmap_lsize(gsmap, mpi_comm)

      ! initialize mct domain

      call wav_domain_mct(lsize, gsmap, mpi_comm, dom)

      ! set flags in infodata
      ! wav_prognostic is set to .false. for debugging purposes only

      call seq_infodata_putdata(infodata, wav_present=.true., &
           wav_prognostic=.true., wav_nx=nx, wav_ny=ny)

      ! initialize mct attribute vectors

      call mct_avect_init(w2x_w, rlist=seq_flds_w2x_fields, lsize=lsize)
      call mct_avect_zero(w2x_w)

      call mct_avect_init(x2w_w, rlist=seq_flds_x2w_fields, lsize=lsize)
      call mct_avect_zero(x2w_w)

      ! add call to gptl timer

      if ( iaproc .eq. napout ) then
        call shr_sys_flush(ndso)
      endif
      call shr_file_setlogunit (shrlogunit)
      call shr_file_setloglevel(shrloglev)

      !--------------------------------------------------------------------
      ! Restart 
      !--------------------------------------------------------------------

      if ( runtype == "continue" .or. runtype == "branch") then
        ifile4 = 0                                   ! reset file counter
        fnmpre = './'//trim(restart_timestamp)//'.'  ! add restart timestamp to file name prefix

        ALLOCATE(MAPTST(NY,NX))
        MAPTST = MAPSTA

        CALL W3IORS ( 'READ', NDS(6), SIG(NK), 1)

        fnmpre = './'                                ! reset file name frefix to default value
        fliwnd = .false.                             ! prevent initialization to fetch-limited spectra
        flcold = .false.                             ! prevent inttialization of DTDYN and FCUT

        ! Compare MAPSTA from grid and restart
        DO IX=1, NX
          DO IY=1, NY
            IF ( ABS(MAPSTA(IY,IX)).EQ.2 .OR.                           &    
                 ABS(MAPTST(IY,IX)).EQ.2 ) THEN 
                MAPSTA(IY,IX) = SIGN ( MAPTST(IY,IX) , MAPSTA(IY,IX) )
            END IF
          END DO
        END DO
        MAPTST = MOD(MAPST2/2,2)
        MAPST2 = MAPST2 - 2*MAPTST
  
        ! Calculate depth
        DO ISEA=1, NSEA 
          IX = MAPSF(ISEA,1)
          IY = MAPSF(ISEA,2)

          WLVeff=WLV(ISEA)
          DW(ISEA) = MAX ( 0. , WLVeff-ZB(ISEA) )
          IF ( WLVeff-ZB(ISEA) .LE.0. ) THEN 
            MAPTST(IY,IX) = 1
            MAPSTA(IY,IX) = -ABS(MAPSTA(IY,IX))
          END IF
        END DO

        DO JSEA=1, NSEAL
          CALL INIT_GET_ISEA(ISEA, JSEA)
          WLVeff=WLV(ISEA)
          DW(ISEA) = MAX ( 0. , WLVeff-ZB(ISEA) )
          IF ( WLVeff-ZB(ISEA) .LE.0. ) THEN 
            VA(:,JSEA) = 0. 
          END IF
        END DO

        MAPST2 = MAPST2 + 2*MAPTST
        DEALLOCATE(MAPTST)

        ! Fill wavenumber and group velocity arrays 
        DO IS=0, NSEA 
          IF (IS.GT.0) THEN 
            DEPTH  = MAX ( DMIN , DW(IS) )
          ELSE 
            DEPTH = DMIN 
          END IF 
          DO IK=0, NK+1 
            CALL WAVNU1(SIG(IK),DEPTH,WN(IK,IS),CG(IK,IS))
          END DO
        END DO

        DW(0) = 0.

        CALL SETUGIOBP

        call mpi_barrier ( mpi_comm, ierr )
      
      endif


900   FORMAT (/15X,'      *** WAVEWATCH III Program shell ***      '/ &
               15X,'==============================================='/)
901   FORMAT ( '  Comment character is ''',A,''''/)

930   FORMAT (/'  Time interval : '/                                  &
           ' --------------------------------------------------')
931   FORMAT ( '       Starting time : ',A)

940   FORMAT (/'  Output requests : '/                                &
           ' --------------------------------------------------'/ &
           '       ',A)

950   FORMAT (/'  Initializations :'/                                 &
           ' --------------------------------------------------')
951   FORMAT ( '       ',A)
1945  FORMAT ( '            Fields   : ',A)
1946  FORMAT ( '                       ',A)

    END SUBROUTINE WAV_INIT_MCT

!=====================================================================
!=====================================================================
!=====================================================================

    SUBROUTINE WAV_RUN_MCT(EClock, cdata_w, x2w_w, w2x_w)

      ! Parameters

      implicit none
      type(ESMF_Clock)            ,intent(inout) :: EClock
      type(seq_cdata)             ,intent(inout) :: cdata_w
      type(mct_aVect)             ,intent(inout) :: x2w_w
      type(mct_aVect)             ,intent(inout) :: w2x_w


      !/ ------------------------------------------------------------------- /
      !/ Local parameters

      integer :: time0(2), timen(2), ierr, i, j, ix, iy, ik
      integer :: ymd              ! current year-month-day
      integer :: tod              ! current time of day (sec)
      integer :: hh,mm,ss
      integer :: n,jsea,isea
      integer :: mpi_comm
      integer :: gindex
      integer :: nu
      integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
      type(mct_aVect) :: x2w0
      type(mct_gsmap),pointer :: gsmap
      real :: def_value
      real, dimension(:), allocatable :: cx, cy 
      real, dimension(:), allocatable :: wx, wy 
      real :: xxx
      real :: wlveff, depth

      character(len=*),parameter :: subname = '(wav_run_mct)'
      character(15) :: restart_date

      !----------------------------------------------------------------------------
      ! Reset shr logging to my log file
      !----------------------------------------------------------------------------
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (stdout)

      call seq_timemgr_EClockGetData( EClock, curr_ymd=ymd, curr_tod=tod )

      hh = tod/3600
      mm = (tod - (hh * 3600))/60
      ss = tod - (hh*3600) - (mm*60)

      timen(1) = ymd
      timen(2) = hh*10000 + mm*100 + ss

      call seq_timemgr_EClockGetData( EClock, prev_ymd=ymd, prev_tod=tod )

      hh = tod/3600
      mm = (tod - (hh * 3600))/60
      ss = tod - (hh*3600) - (mm*60)

      time0(1) = ymd
      time0(2) = hh*10000 + mm*100 + ss

      time = time0


      !--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! 7.  Model with input
      ! 7.a Determine next time interval and input fields
      ! 7.a.1 Preparation
      ! 7.a.3 Update time and fields / data

      !--- input fields associated with W3FLDG calls in ww3_shel.ftn
      !--- these arrays are global, just fill the local cells for use later
      !--- fill both the lower (0) and upper (N) bound data with the same values
      !--- fill with special values as default, these should not be used in practice
      !--- set time for input data to time0 and timen (shouldn't matter)

      ! def_value = -1.0e36
      def_value = 0.0

      if (inflags1(1)) then
         TLN  = timen
         WLEV = def_value   ! ssh
      endif

      if (inflags1(2)) then
         TC0  = time0
         TCN  = timen
         CX0  = def_value   ! ocn u current
         CXN  = def_value
         CY0  = def_value   ! ocn v current
         CYN  = def_value
      endif

      if (inflags1(3)) then
         TW0  = time0
         TWN  = timen
         WX0  = def_value   ! atm u wind
         WXN  = def_value
         WY0  = def_value   ! atm v wind
         WYN  = def_value
         DT0  = def_value   ! air temp - ocn temp
         DTN  = def_value
      endif

      if (inflags1(4)) then
         TIN  = timen
         ICEI = def_value   ! ice frac
         TI1 = timen
         ICEP1 = def_value
         TI5 = timen
         ICEP5 = 1000.0
      endif

      ! this is the global fill
      call seq_cdata_setptrs(cdata_w,gsmap=gsmap,mpicom=mpi_comm)
      call mct_aVect_gather(x2w_w,x2w0,gsmap,0,mpi_comm)
      call mct_aVect_bcast(x2w0,0,mpi_comm)
 
      if (inflags1(2)) then 
        allocate(cx(NX*NY), cy(NX*NY))
      endif
      if (inflags1(3)) then
        allocate(wx(NX*NY), wy(NX*NY))
      endif

      gindex = 0
      do IY = 1,NY
      do IX = 1,NX
         gindex = gindex + 1

         if (inflags1(2)) then
            CX(gindex)  = x2w0%rattr(index_x2w_so_u,gindex)
            CY(gindex)  = x2w0%rattr(index_x2w_so_v,gindex)
         endif

         if (inflags1(3)) then
            WX(gindex)  = x2w0%rattr(index_x2w_sa_u,gindex)
            WY(gindex)  = x2w0%rattr(index_x2w_sa_v,gindex)
         endif

      enddo
      enddo

      if (inflags1(2)) then
        call w3xyrtn(NX*NY,CX,CY,-AnglD) 
      endif
      if (inflags1(3)) then
        call w3xyrtn(NX*NY,WX,WY,-AnglD)
      endif

      ! use these loops for global copy
      gindex = 0
      do IY = 1,NY
      do IX = 1,NX
         gindex = gindex + 1
         ! use this loop to do local copy
         !      do jsea=1, nseal
         !         isea = iaproc + (jsea-1)*naproc
         !         IX  = MAPSF(ISEA,1)
         !         IY  = MAPSF(ISEA,2)
         !         gindex = ix + (iy-1)*nx

         if (inflags1(1)) then
            WLEV(IX,IY) = x2w0%rattr(index_x2w_so_ssh,gindex)

            ! Fill wavenumber and group velocity arrays
            ISEA = MAPFS(IY,IX)
            WLV(ISEA) = WLEV(IX,IY)
            WLVeff=WLV(ISEA)
            DW(ISEA) = MAX ( 0. , WLVeff-ZB(ISEA) )
            IF (ISEA.GT.0) THEN 
              DEPTH  = MAX ( DMIN , DW(ISEA) )
            ELSE 
              DEPTH = DMIN 
            END IF 
            DO IK=0, NK+1 
              CALL WAVNU1(SIG(IK),DEPTH,WN(IK,ISEA),CG(IK,ISEA))
            END DO
         endif

         if (inflags1(2)) then
            CX0(IX,IY)  = CX(gindex) 
            CXN(IX,IY)  = CX0(IX,IY)
            CY0(IX,IY)  = CY(gindex)
            CYN(IX,IY)  = CY0(IX,IY)
         endif

         if (inflags1(3)) then
            WX0(IX,IY)  = WX(gindex) 
            WXN(IX,IY)  = WX0(IX,IY)
            WY0(IX,IY)  = WY(gindex) 
            WYN(IX,IY)  = WY0(IX,IY)
            !DT0(IX,IY)  = x2w0%rattr(index_x2w_sa_tbot,gindex) - x2w0%rattr(index_x2w_so_t,gindex)
            !DTN(IX,IY)  = DT0(IX,IY)
         endif

         if (inflags1(4)) then
            ICEI(IX,IY) = x2w0%rattr(index_x2w_si_ifrac,gindex)
            ICEP1(IX,IY) = x2w0%rattr(index_x2w_si_ithick,gindex)
            !ICEP5(IX,IY) = x2w0%rattr(index_x2w_si_ifloe,gindex)

         endif

      enddo
      enddo

      call mct_aVect_clean(x2w0)

      ! 7.b Run the wave model for the given interval

      call w3wave ( 1, odat, timen )

      ! rotate stokes drift
      do i = 1,usspf(2)
        call w3xyrtn(nseal,USSP(1:nseal,i),USSP(1:nseal,nk+i),AnglDL)
      enddo

      ! copy ww3 data to coupling datatype
      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         IX  = MAPSF(ISEA,1)
         IY  = MAPSF(ISEA,2)
         if (MAPSTA(IY,IX) .eq. 1) then

             if (wav_ocn_coup .eq. 'twoway') then
                w2x_w%rattr(index_w2x_Sw_Hs,jsea) = HS(jsea)
                w2x_w%rattr(index_w2x_Sw_Fp,jsea) = FP0(jsea)
                w2x_w%rattr(index_w2x_Sw_Dp,jsea) = THP0(jsea)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_1,jsea) = USSP(jsea,1)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_1,jsea) = USSP(jsea,nk+1)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_2,jsea) = USSP(jsea,2)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_2,jsea) = USSP(jsea,nk+2)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_3,jsea) = USSP(jsea,3)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_3,jsea) = USSP(jsea,nk+3)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_4,jsea) = USSP(jsea,4)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_4,jsea) = USSP(jsea,nk+4)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_5,jsea) = USSP(jsea,5)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_5,jsea) = USSP(jsea,nk+5)

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_6,jsea) = USSP(jsea,6)
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_6,jsea) = USSP(jsea,nk+6)
             endif
          else
             if (wav_ocn_coup .eq. 'twoway') then
                w2x_w%rattr(index_w2x_Sw_Hs,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_Fp,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_Dp,jsea) = 0.0
            
                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_1,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_1,jsea) = 0.0

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_2,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_2,jsea) = 0.0

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_3,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_3,jsea) = 0.0

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_4,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_4,jsea) = 0.0

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_5,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_5,jsea) = 0.0

                w2x_w%rattr(index_w2x_Sw_ustokes_wavenumber_6,jsea) = 0.0
                w2x_w%rattr(index_w2x_Sw_vstokes_wavenumber_6,jsea) = 0.0
          endif
        endif
      enddo

      !----------------------------------------------------------------------------
      ! Reset shr logging to original values
      !----------------------------------------------------------------------------
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)
      if ( iaproc .eq. napout ) then
        call shr_sys_flush(stdout)
      endif

      !----------------------------------------------------------------------------
      ! Restart file output 
      !----------------------------------------------------------------------------
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
        if (iaproc == 1) then
          write(restart_date,"(i8.8,'.'i6.6)") time(1),time(2)
          nu = shr_file_getUnit()
          open(nu,file='rpointer.wav',form='formatted')
          write(nu,'(a)') restart_date
          close(nu)
          call shr_file_freeUnit(nu)        
        endif
        CALL W3IORS ('HOT', NDS(6), XXX, 1, .TRUE. )
      endif


    END SUBROUTINE WAV_RUN_MCT

!=====================================================================
!=====================================================================
!=====================================================================

    SUBROUTINE WAV_FINAL_MCT(EClock, cdata_w, x2w_w, w2x_w)

      ! Parameters

      implicit none
      type(ESMF_Clock)            ,intent(inout) :: EClock
      type(seq_cdata)             ,intent(inout) :: cdata_w
      type(mct_aVect)             ,intent(inout) :: x2w_w
      type(mct_aVect)             ,intent(inout) :: w2x_w

      character(len=*),parameter :: subname = '(wav_final_mct)'

      ! do nothing now

    END SUBROUTINE WAV_FINAL_MCT

!=====================================================================
!=====================================================================
!=====================================================================

    subroutine wav_setgsmap_mct(mpi_comm, compid, gsmap)

      !/ ------------------------------------------------------------------- /
      !use w3gdatmd, only: nx, ny, nseal
      !use w3odatmd, only: naproc
      implicit none
      !/
      !/ ------------------------------------------------------------------- /
      !/ parameter list
      !/
      integer        , intent(in)  :: mpi_comm
      integer        , intent(in)  :: compid
      type(mct_gsmap), intent(out) :: gsmap
      !/
      !/ ------------------------------------------------------------------- /
      !/ local parameters
      !/
      integer, allocatable :: gindex(:)
      integer :: n,jsea,isea,ix,iy
      character(len=*),parameter :: subname = '(wav_setgsmap_mct)'
      ! -------------------------------------------------------------------- /

      allocate(gindex(nseal))
      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         gindex(jsea) = ix + (iy-1)*nx
      end do
      call mct_gsmap_init( gsmap, gindex, mpi_comm, compid, nseal, nx*ny)
      deallocate(gindex)

    end subroutine wav_setgsmap_mct

!=====================================================================
!=====================================================================
!=====================================================================

    subroutine wav_domain_mct(lsize, gsmap, mpi_comm, dom)

      implicit none
      integer        , intent(in)   :: lsize
      type(mct_gsmap), intent(in)   :: gsmap
      integer        , intent(in)  :: mpi_comm
      type(mct_ggrid), intent(inout):: dom

      integer  :: n,i,ix,iy,isea,jsea   ! indices
      real(r8) :: lon, lat, mask
      real(r8), pointer  :: data(:)     ! temporary
      integer , pointer  :: idata(:)    ! temporary
      real(r8), parameter:: rad2deg = 180.0_r8/shr_const_pi
      real(r8), parameter:: deg2rad = shr_const_pi/180.0_r8
      character(len=*),parameter :: subname = '(wav_domain_mct)'

      ! initialize mct domain type
      ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
      ! note that in addition land carries around landfrac for the purposes of domain checking

      call mct_ggrid_init( ggrid=dom, coordchars=trim(seq_flds_dom_coord), &
           otherchars=trim(seq_flds_dom_other), lsize=lsize )

      ! allocate memory

      allocate(data(lsize))

      ! determine global gridpoint number attribute, globgridnum, which is set automatically by mct

      call mct_gsMap_orderedPoints(gsmap, iaproc-1, idata)
      call mct_gGrid_importIattr(dom,'GlobGridNum',idata,lsize)

      ! determine domain (numbering scheme is: west to east and south to north to south pole)
      ! initialize attribute vector with special value

      data(:) = -9999.0_r8
      call mct_ggrid_importrattr(dom,"lat"  ,data,lsize)
      call mct_ggrid_importrattr(dom,"lon"  ,data,lsize)
      call mct_ggrid_importrattr(dom,"area" ,data,lsize)
      call mct_ggrid_importrattr(dom,"aream",data,lsize)
      data(:) = 0.0_r8
      call mct_ggrid_importrattr(dom,"mask" ,data,lsize)

      ! fill in correct values for domain components
      ! note aream will be filled in in the atm-lnd mapper
      ! sx, sy  real  i  grid increments (deg.).

      if (iaproc .ne. 1) then
        allocate(grid_area(nx))
      endif
      call shr_mpi_bcast(grid_area, mpi_comm)

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         lon = xgrd(iy,ix) 
         data(jsea) = lon
         !write(stdout,*)' jsea= ',jsea,' lon is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lon",data,lsize)

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         lat = ygrd(iy,ix)
         data(jsea) = lat
         !write(stdout,*)' jsea= ',jsea,' lat is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lat",data,lsize)

      do jsea = 1,nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         data(jsea) = grid_area(ix)*deg2rad**2
         !write(stdout,*)' jsea= ',jsea,' area is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"area",data,lsize)

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         ! QL, 150827, should be 1 for all sea point
         !mask = mapsta(iy,ix)
         if (mapsta(iy,ix) .ne. 0) then
            mask = 1.0_r8
         else
            mask = 0.0_r8
         end if
         data(jsea) = mask
         !write(stdout,*)' jsea= ',jsea,' mask is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"mask",data,lsize)
      call mct_ggrid_importrattr(dom,"frac",data,lsize)

      n = mct_aVect_lSize(dom%data)

      deallocate(data)
      deallocate(idata)

    end subroutine wav_domain_mct

!=====================================================================
!=====================================================================
!=====================================================================

    subroutine read_stations_file(ndso,fname,npts,x,y,pnames)

    implicit none

    integer, intent(in) :: ndso
    character(*), intent(in) :: fname
    integer, intent(out) :: npts
    real, dimension(:), allocatable, intent(out) :: x
    real, dimension(:), allocatable, intent(out) :: y
    character(len=40), dimension(:), allocatable, intent(out) :: pnames

    integer :: ndsl
    integer :: ipts
    integer :: iloop
    integer :: ierr
    character(len=256) :: tmpline,test
    character(len=1) :: comstr
    real :: xx, yy
    character(len=40) :: pn

    ! Adapted from ww3_shel.ftn

    NDSL = shr_file_getunit()
    COMSTR = "$"

    OPEN(NDSL, FILE=TRIM(ADJUSTL(fname)), FORM='FORMATTED', STATUS='OLD', ERR=2104, IOSTAT=IERR)
    
    ! first loop to count the number of points
    ! second loop to allocate the array and store the points
    IPTS = 0
    DO ILOOP=1,2
      REWIND (NDSL)
    
      IF ( ILOOP.EQ.2) THEN 
        NPTS = IPTS 
        IF ( NPTS.GT.0 ) THEN 
          ALLOCATE ( X(NPTS), Y(NPTS), PNAMES(NPTS) )
          IPTS = 0 ! reset counter to be reused for next do loop
        ELSE 
          ALLOCATE ( X(1), Y(1), PNAMES(1) )
          return 
        END IF
      END IF
    
      DO   
        READ (NDSL,*,ERR=2004,IOSTAT=IERR) TMPLINE
        ! if end of file or stopstring, then exit
        IF ( IERR.NE.0 .OR. INDEX(TMPLINE,"STOPSTRING").NE.0 ) EXIT 
        ! leading blanks removed and placed on the right
        TEST = ADJUSTL ( TMPLINE )
        IF ( TEST(1:1).EQ.COMSTR .OR. LEN_TRIM(TEST).EQ.0 ) THEN 
          ! if comment or blank line, then skip
          CYCLE
        ELSE 
          ! otherwise, backup to beginning of line
          BACKSPACE ( NDSL, ERR=2004, IOSTAT=IERR)
          READ (NDSL,*,ERR=2004,IOSTAT=IERR) XX, YY, PN
        END IF
        IPTS = IPTS + 1
        IF ( ILOOP .EQ. 1 ) CYCLE
        IF ( ILOOP .EQ. 2 ) THEN 
          X(IPTS)      = XX 
          Y(IPTS)      = YY 
          PNAMES(IPTS) = PN 
          IF ( IAPROC .EQ. NAPOUT ) THEN 
            IF ( IPTS .EQ. 1 ) THEN 
              WRITE (NDSO,2945) XX, YY, PN
            ELSE 
              WRITE (NDSO,2946) IPTS, XX, YY, PN
            END IF
          END IF
        END IF ! ILOOP.EQ.2
      END DO ! end of file                      
    END DO ! ILOOP
    CLOSE(NDSL)
    call shr_file_freeUnit(NDSL)
    

 2945 FORMAT ( '            Point  1 : ',2F8.2,2X,A)
 2946 FORMAT ( '              ',I6,' : ',2F8.2,2X,A)
 1104 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
               '     ERROR IN OPENING POINT FILE'/                    &
               '     IOSTAT =',I5/)
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3SHEL : *** '/           &
               '     ERROR IN READING FROM POINT FILE'/               &
               '     IOSTAT =',I5/)

 2104 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSO,1104) IERR 
 2004 CONTINUE
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSO,1004) IERR 

    end subroutine read_stations_file

  END MODULE WAV_COMP_MCT
