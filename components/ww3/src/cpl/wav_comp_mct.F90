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
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!       !/MPI   Id.
!
!       !/LLG   Spherical grid.
!       !/XYG   Cartesian grid.
!       !/MGW   Moving grid wind correction.
!       !/MGP   Moving grid propagation correction.
!
!       !/T     Enable test output.
!       !/O7    Echo input homogeneous fields.
!
!       !/NCO   NCEP NCO modifications for operational implementation.
!
!       !/F90   Timer function included for F90.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      use w3gdatmd, only: dtmax, dtcfl, dtcfli, dtmin, nx, ny, nseal, mapsf, mapsta, x0, y0, sx, sy, w3nmod, w3setg
      use w3wdatmd, only: time, w3ndat, w3setw
      use w3adatmd, only: ussx, ussy, w3naux, w3seta, sxx, sxy, syy !SB, lamult
      use w3idatmd, only: inflags1, w3seti, w3ninp
      USE W3IDATMD, ONLY: TC0, CX0, CY0, TCN, CXN, CYN
      USE W3IDATMD, ONLY: TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN
      USE W3IDATMD, ONLY: TIN, ICEI
      USE W3IDATMD, ONLY: TLN, WLEV
      USE W3IDATMD, ONLY: !SB, HML   ! QL, 150525, mixing layer depth
      use w3odatmd, only: w3nout, w3seto, naproc, iaproc, napout, naperr,             &
                          nogrp, ngrpp, noge, idout, fnmpre, iostyp, notype
!/
      use w3initmd, only: w3init 
      use w3wavemd, only: w3wave 
!/
      use w3iopomd, only:
      use w3timemd, only: stme21 
      use w3cesmmd, only : casename, initfile, rstwr, runtype, histwr, outfreq
      use w3cesmmd, only : inst_index, inst_name, inst_suffix

      use esmf
      use mct_mod
      use seq_flds_mod

      use ww3_cpl_indices  , only : ww3_cpl_indices_set
      use ww3_cpl_indices  , only : index_x2w_Sa_u, index_x2w_Sa_v, index_x2w_Sa_tbot, index_x2w_Si_ifrac
      use ww3_cpl_indices  , only : index_x2w_So_t, index_x2w_So_u, index_x2w_So_v, index_x2w_So_bldepth, index_x2w_So_ssh
      use ww3_cpl_indices  , only : index_w2x_Sw_ustokes, index_w2x_Sw_vstokes

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

      integer,save :: stdout

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
      ! QL, 150629, calculating restart interval
      integer :: stop_ymd          ! stop date (yyyymmdd)
      integer :: stop_tod          ! stop time of day (sec)
      integer :: ix, iy

      character(CL)            :: starttype
      type(mct_gsmap), pointer :: gsmap
      type(mct_ggrid), pointer :: dom
      type(seq_infodata_type), pointer :: infodata   ! input init object

      integer             :: unitn            ! namelist unit number
      integer             :: ndso, ndse, nds(13), ntrace(2), time0(2)
      integer             :: timen(2), odat(35), nh(4), iprt(6)
      integer             :: i,j,npts
      integer             :: ierr
      integer             :: jsea,isea
      real                :: a(nhmax,4)
      real, allocatable   :: x(:), y(:)
      logical             :: flogrd(nogrp,ngrpp), flogrd2(nogrp,ngrpp)
      logical             :: flogd(nogrp), flogd2(nogrp)
      logical             :: prtfrm, flt
      character(len=*),parameter :: subname = '(wav_init_mct)'

      character(len=3)    :: idstr(8), idtst
      character(len=10)   :: pn
      character(len=13)   :: idflds(8)
      character(len=20)   :: strng
      character(len=23)   :: dtme21
      character(len=30)   :: idotyp(6)
      character(len=10), allocatable :: pnames(:)

      !/ ------------------------------------------------------------------- /

      DATA IDFLDS / 'water levels ' , 'currents     ' ,               &
                    'winds        ' , 'ice fields   ' ,               &
                    'mean param.  ' , '1D spectra   ' ,               &
                    '2D spectra   ' , 'moving grid  ' /
      DATA IDOTYP / 'Fields of mean wave parameters' ,                &
                    'Point output                  ' ,                &
                    'Track point output            ' ,                &
                    'Restart files                 ' ,                &
                    'Nesting data                  ' ,                &
                    'Partitioned wave field data   ' /
      DATA IDSTR  / 'LEV', 'CUR', 'WND', 'ICE', 'DT0', 'DT1', 'DT2',  &
                    'MOV' /

      namelist /ww3_inparm/ initfile, outfreq

      !--------------------------------------------------------------------
      ! Set up data structures
      !--------------------------------------------------------------------

      call w3nmod ( 1, 6, 6 )
      call w3ndat (    6, 6 )
      call w3naux (    6, 6 )
      call w3nout (    6, 6 )
      call w3ninp (    6, 6 )

      call w3setg ( 1, 6, 6 )
      call w3setw ( 1, 6, 6 )
      call w3seta ( 1, 6, 6 )
      call w3seto ( 1, 6, 6 )
      call w3seti ( 1, 6, 6 )

      !--------------------------------------------------------------------
      ! Initialize mpi
      !--------------------------------------------------------------------

      call seq_cdata_setptrs(cdata, id=compid, mpicom=mpi_comm, &
           gsmap=gsmap, dom=dom, infodata=infodata)
      inst_name   = seq_comm_name(compid)
      inst_index  = seq_comm_inst(compid)
      inst_suffix = seq_comm_suffix(compid)
      call ww3_cpl_indices_set()

      call mpi_comm_size(mpi_comm, naproc, ierr)
      call mpi_comm_rank(mpi_comm, iaproc, ierr)
      iaproc = iaproc + 1

      !--------------------------------------------------------------------
      ! Initialize run type
      ! QL, 150525
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

      napout = 1
      naperr = 1

      ! 1.b For WAVEWATCH III (See W3INIT) ??? ask adrean if i am missing something
      !
      ! The following units are referenced in module w3initmd
      ! NDS(1) ! OUTPUT LOG: General output unit number ("log file") (NDS0)
      ! NDS(2) ! OUTPUT LOG: Error output unit number (NDSE)
      ! NDS(3) ! OUTPUT LOG: Test output unit number (NDST)
      ! NDS(4) ! OUTPUT LOG: Unit for 'direct' output (SCREEN)
      !
      ! NDS(5) ! INPUT: mod_def.ww3 file (model definition) unit number
      ! NDS(9) ! INPUT: unit for read in boundary conditions (based on FLBPI)
      !
      ! The following units are referenced in module w3wavemd for output
      ! NDS( 6) ! OUTPUT DATA: restart(N).ww3 file (model restart) unit number
      ! NDS( 7) ! unit for output for FLOUT(1) flag
      ! NDS( 8) ! unit for output for FLOUT(2) flag
      ! NDS(11) ! unit for output for FLOUT(3) flag
      ! NDS(12) ! unit for output for FLOUT(3) flag
      ! NDS(13) ! unit for output for FLOUT(6) flag
      ! NDS(10) ! unit for output for FLOUT(5) flag

      if (iaproc .eq. napout) then
         stdout = shr_file_getunit()
         call shr_file_setio('wav_modelio.nml'//trim(inst_suffix),stdout)
      else
         stdout = 6
      endif

      nds( 1) = stdout
      nds( 2) = stdout
      nds( 3) = stdout
      nds( 4) = stdout
      nds( 5) = shr_file_getunit()
      nds( 6) = shr_file_getunit()
      nds( 7) = shr_file_getunit()
      nds( 8) = shr_file_getunit()
      nds( 9) = shr_file_getunit()
      nds(10) = shr_file_getunit()
      nds(11) = shr_file_getunit()
      nds(12) = shr_file_getunit()
      nds(13) = shr_file_getunit()

      ndso      =  stdout
      ndse      =  stdout
      ntrace(1) =  nds(3)
      ntrace(2) =  10

      ! Redirect share output to wav log
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (ndso)

      if ( iaproc .eq. napout ) write (ndso,900)
      call shr_sys_flush(ndso)

      if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
         write(ndso,*) 'starttype: initial'
      else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
         write(ndso,*) 'starttype: continue'
      else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
         write(ndso,*) 'starttype: branch'
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
      ! QL, 150525, flags for passing variables from coupler to ww3,
      !             lev, curr, wind, ice and mixing layer depth on
      inflags1(1:5) = .true.
!      inflags1(1:4) = .true.   !changed by Adrean (lev,curr,wind,ice on)
!      inflags1(3:4) = .true.   !changed by Adrean (wind,ice on)

      !--------------------------------------------------------------------
      ! Set time frame
      !--------------------------------------------------------------------

      ! TIME0 = from ESMF clock
      ! NOTE - are not setting TIMEN here

      if ( iaproc .eq. napout ) write (ndso,930)
      call shr_sys_flush(ndso)

      ! QL, 150525, initial run or restart run
      if ( runtype .eq. "initial") then
         call seq_timemgr_EClockGetData(EClock, &
              start_ymd=start_ymd, start_tod=start_tod)
      else
         call seq_timemgr_EClockGetData(EClock, &
              curr_ymd=start_ymd, curr_tod=start_tod)
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
      if ( iaproc .eq. napout ) write (ndso,931) dtme21
      call shr_sys_flush(ndso)
      time = time0

      !--------------------------------------------------------------------
      ! Define output type and fields
      !--------------------------------------------------------------------

      iostyp = 1        ! gridded field
      write (ndso,940) 'no dedicated output process, any file system '
      call shr_sys_flush(ndso)

      ! TODO - need to enable restart files in run
      ! Actually will need a new restart flag - since all of the ODAT
      ! should be set to 0 - since they are initializated in w3initmd
      ! ODAT    I.A.   I   Output data, five parameters per output type
      !                          1 YYYMMDD for first output.
      !                          2 HHMMSS for first output.
      !                          3 Output interval in seconds.
      !                          4 YYYMMDD for last output.
      !                          5 HHMMSS for last output.
      !                     1-5  Data for OTYPE = 1; gridded fields.
      !                     6-10 Id.  for OTYPE = 2; point output.
      !                    11-15 Id.  for OTYPE = 3; track point output.
      !                    16-20 Id.  for OTYPE = 4; restart files.
      !                    21-25 Id.  for OTYPE = 5; boundary data.
      !                    31-35 Id.  for OTYPE = 5; coupling data.
      ! FLGRD   L.A.   I   Flags for gridded output.
      ! NPT     Int.   I   Number of output points
      ! X/YPT   R.A.   I   Coordinates of output points.
      ! PNAMES  C.A.   I   Output point names.

      npts = 0
      allocate ( x(1), y(1), pnames(1) )
      pnames(1) = ' '

      notype = 7

      do j=1, notype
         odat(5*(j-1)+3) = 0
      end do

      ! QL, 160823, initialize flag for restart
      rstwr = .false.
      ! QL, 160601, initialize flag for history file
      histwr = .false.

      ! QL, 160601, get coupling interval
      call seq_timemgr_eclockgetdata(eclock, dtime=dtime_sync )
      !DEBUG
      ! Hardwire gridded output for now
      ! first output time stamp is now read from file
      ! QL, 150525, 1-5 for history files, 16-20 for restart files
      !     150629, restart output interval is set to the total time of run
      !     150823, restart is taken over by rstwr
      !     160601, output interval is set to coupling interval, so that
      !             variables calculated in W3IOGO could be updated at
      !             every coupling interval
      odat(1) = time(1)     ! YYYYMMDD for first output
      odat(2) = time(2)     ! HHMMSS for first output
      odat(3) = dtime_sync  ! output interval in sec ! changed by Adrean
      odat(4) = 99990101    ! YYYYMMDD for last output
      odat(5) = 0           ! HHMMSS for last output
      odat(16) = time(1)    ! YYYYMMDD for first output
      odat(17) = time(2)    ! HHMMSS for first output
      odat(18) = 86400*5    ! output interval in sec
      odat(19) = 99990101   ! YYYYMMDD for last output
      odat(20) = 0          ! HHMMSS for last output
      !DEBUG

      ! Output Type 1: fields of mean wave parameters gridded output
      flogrd2 = .false.
      flogd2  = .false.
      flogd = .false.

      ! 1) Forcing fields
      flogrd(1,1)  = .false. ! Water depth         
      flogrd(1,2)  = .false. ! Current vel.        
      flogrd(1,3)  = .true.  ! Wind speed          
      flogrd(1,4)  = .false. ! Air-sea temp. dif.  
      flogrd(1,5)  = .false. ! Water level         
      flogrd(1,6)  = .false. ! Ice concentration   
      flogrd(1,7)  = .false. ! Iceberg damp coeffic

      ! 2) Standard mean wave parameters
      flogrd(2,1)  = .true.  ! Wave height         
      flogrd(2,2)  = .false. ! Mean wave length    
      flogrd(2,3)  = .false. ! Mean wave period(+2)
      flogrd(2,4)  = .false. ! Mean wave period(-1)
      flogrd(2,5)  = .false. ! Mean wave period(+1)
      flogrd(2,6)  = .true.  ! Peak frequency      
      flogrd(2,7)  = .false. ! Mean wave dir. a1b1 
      flogrd(2,8)  = .false. ! Mean dir. spr. a1b1 
      flogrd(2,9)  = .true.  ! Peak direction      
      flogrd(2,10) = .false. ! Infragravity height
      flogrd(2,11) = .false. ! Space-Time Max E   
      flogrd(2,12) = .false. ! Space-Time Max Std 
      flogrd(2,13) = .false. ! Space-Time Hmax    
      flogrd(2,14) = .false. ! Spc-Time Hmax^crest
      flogrd(2,15) = .false. ! STD Space-Time Hmax
      flogrd(2,16) = .false. ! STD ST Hmax^crest  

      ! 3) Frequency-dependent standard parameters
      flogrd(3,1)  = .false. ! 1D Freq. Spectrum   
      flogrd(3,2)  = .false. ! Mean wave dir. a1b1 
      flogrd(3,3)  = .false. ! Mean dir. spr. a1b1 
      flogrd(3,4)  = .false. ! Mean wave dir. a2b2 
      flogrd(3,5)  = .false. ! Mean dir. spr. a2b2 
      flogrd(3,6)  = .false. ! Wavenumber array    

      ! 4) Spectral Partitions parameters
      flogrd(4,1) = .false. ! Part. wave heigt    
      flogrd(4,2) = .false. ! Part. peak period   
      flogrd(4,3) = .false. ! Part. peak wave l.  
      flogrd(4,4) = .false. ! Part. mean direction
      flogrd(4,5) = .false. ! Part. dir. spread   
      flogrd(4,6) = .false. ! Part. wind sea frac.
      flogrd(4,7) = .false. ! Total wind sea frac.
      flogrd(4,8) = .false. ! Number of partitions

      ! 5) Atmosphere-waves layer
      flogrd(5,1) = .false. ! Friction velocity   
      flogrd(5,2) = .false. ! Charnock parameter  
      flogrd(5,3) = .false. ! Energy flux         
      flogrd(5,4) = .false. ! Wind-wave enrgy flux
      flogrd(5,5) = .false. ! Wind-wave net mom. f
      flogrd(5,6) = .false. ! Wind-wave neg.mom.f.
      flogrd(5,7) = .false. ! Whitecap coverage   
      flogrd(5,8) = .false. ! Whitecap mean thick.
      flogrd(5,9) = .false. ! Mean breaking height
      flogrd(5,10) = .false. ! Dominant break prob 

      ! 6) Wave-ocean layer
      flogrd(6,1) = .false. ! Radiation stresses  
      flogrd(6,2) = .false. ! Wave-ocean mom. flux
      flogrd(6,3) = .false. ! wave ind p Bern Head
      flogrd(6,4) = .false. ! Wave-ocean TKE  flux
      flogrd(6,5) = .false. ! Stokes transport    
      flogrd(6,6) = .false. ! Stokes drift at z=0 
      flogrd(6,7) = .false. ! 2nd order pressure  
      flogrd(6,8) = .false. ! Stokes drft spectrum
      flogrd(6,9) = .false. ! 2nd ord press spectr
      flogrd(6,10) = .false. ! Wave-ice mom. flux  
      flogrd(6,11) = .false. ! Wave-ice energy flux

      ! 7) Wave-bottom layer
      flogrd(7,1) = .false. ! Bottom rms ampl.    
      flogrd(7,2) = .false. ! Bottom rms velocity 
      flogrd(7,3) = .false. ! Bedform parameters  
      flogrd(7,4) = .false. ! Energy diss. in WBBL
      flogrd(7,5) = .false. ! Moment. loss in WBBL

      ! 8) Spectrum parameters
      flogrd(8,1) = .false. ! Mean square slopes  
      flogrd(8,2) = .false. ! Phillips tail const
      flogrd(8,3) = .false. ! Lx-Ly mean wvlength
      flogrd(8,4) = .false. ! Surf grad correl XT
      flogrd(8,5) = .false. ! Surf grad correl YT
      flogrd(8,6) = .false. ! Surf grad correl XY
      flogrd(8,7) = .false. ! Surface crest param

      ! 9) Numerical diagnostics
      flogrd(9,1) = .false. ! Avg. time step.     
      flogrd(9,2) = .false. ! Cut-off freq.       
      flogrd(9,3) = .false. ! Maximum spatial CFL 
      flogrd(9,4) = .false. ! Maximum angular CFL 
      flogrd(9,5) = .false. ! Maximum k advect CFL

      do i = 1, nogrp
        if(any(flogrd(i,:))) then
          flogd(i) = .true.
        end if
      end do

      if ( iaproc .eq. napout ) then
         flt = .true.
         do i=1, nogrp
            do j = 1, noge(i)
              if ( flogrd(i,j) ) then
                 if ( flt ) then
                    write (ndso,1945) idout(i,j)
                    flt    = .false.
                 else
                    write (ndso,1946) idout(i,j)
                 end if
              end if
            end do
         end do
         if ( flt ) write (ndso,1945) 'no fields defined'
      end if

      call shr_sys_flush(ndso)

      !--------------------------------------------------------------------
      ! Wave model initializations
      !--------------------------------------------------------------------

      ! Notes on ww3 initialization:
      ! ww3 read initialization occurs in w3iors (which is called by initmd)
      ! For a startup (including hybrid) or branch run the initial datafile is
      ! set in namelist input 'initfile'
      ! For a continue run - the initfile vluae is created from the time(1:2)
      ! array set below

      if ( iaproc .eq. napout ) write (ndso,950)
      if ( iaproc .eq. napout ) write (ndso,951) 'wave model ...'
      call shr_sys_flush(ndso)

      ! Read namelist (set initfile in w3cesmmd)
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
      call shr_mpi_bcast(initfile, mpi_comm)
      call shr_mpi_bcast(outfreq, mpi_comm)

      ! Set casename (in w3cesmmd)
      call seq_infodata_GetData(infodata,case_name=casename)

      ! Read in input data and initialize the model
      ! w3init calls w3iors which:
      ! - reads either the initfile if the run is startup or branch
      ! - constructs the filename from the casename variable and the time(:) array
      !   which is set above
      call w3init ( 1, .false., 'ww3', nds, ntrace, odat, flogrd, flogrd2, flogd, flogd2, npts, x, y,   &
           pnames, iprt, prtfrm, mpi_comm )
      call shr_sys_flush(ndso)

      ! overwrite dt values with variables from coupler
      ! is this a problem with any things being set in w3init?
      dtmax  = real(dtime_sync)
      dtcfl  = real(dtime_sync) / 2. !checked by adrean
      dtcfli = real(dtime_sync)      !checked by adrean
      dtmin  = real(dtime_sync) / 12 !checked by adrean

      call mpi_barrier ( mpi_comm, ierr )

      !--------------------------------------------------------------------
      ! cpl7/mct initialization
      !--------------------------------------------------------------------

      ! initialize mct gsmap

      call wav_setgsmap_mct(mpi_comm, compid, gsmap)
      lsize = mct_gsmap_lsize(gsmap, mpi_comm)
      call shr_sys_flush(ndso)

      ! initialize mct domain

      call wav_domain_mct(lsize, gsmap, dom)

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

      ! QL, 150823, send initial state to driver
      ! QL, 160611, initial values for lamult, ustokes and vstokes
      do jsea=1, nseal
          w2x_w%rattr(index_w2x_Sw_ustokes,jsea) = 0.
          w2x_w%rattr(index_w2x_Sw_vstokes,jsea) = 0.
      enddo

      ! end redirection of share output to wav log

      call shr_sys_flush(ndso)
      call shr_file_setlogunit (shrlogunit)
      call shr_file_setloglevel(shrloglev)


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

      integer :: time0(2), timen(2), ierr, i, j, ix, iy
      integer :: ymd              ! current year-month-day
      integer :: tod              ! current time of day (sec)
      integer :: hh,mm,ss
      integer :: n,jsea,isea
      integer :: mpi_comm
      integer :: gindex
      integer(IN)   :: shrlogunit, shrloglev ! original log unit and level
      type(mct_aVect) :: x2w0
      type(mct_gsmap),pointer :: gsmap
      real :: def_value

      character(len=*),parameter :: subname = '(wav_run_mct)'

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

      ! QL, 171107, output every outfreq hours
      if (outfreq .gt. 0 .and. mod(hh, outfreq) .eq. 0 ) then
          histwr = .true.
      else
          histwr = seq_timemgr_historyAlarmIsOn(EClock)
      end if


      hh = tod/3600
      mm = (tod - (hh * 3600))/60
      ss = tod - (hh*3600) - (mm*60)

      time0(1) = ymd
      time0(2) = hh*10000 + mm*100 + ss

      time = time0

      ! QL, 150823, set flag for writing restart file
      rstwr = seq_timemgr_RestartAlarmIsOn(EClock)

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
      endif

      ! this is the global fill
      call seq_cdata_setptrs(cdata_w,gsmap=gsmap,mpicom=mpi_comm)
      call mct_aVect_gather(x2w_w,x2w0,gsmap,0,mpi_comm)
      call mct_aVect_bcast(x2w0,0,mpi_comm)

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
         endif

         if (inflags1(2)) then
            CX0(IX,IY)  = x2w0%rattr(index_x2w_so_u,gindex)
            CXN(IX,IY)  = CX0(IX,IY)
            CY0(IX,IY)  = x2w0%rattr(index_x2w_so_v,gindex)
            CYN(IX,IY)  = CY0(IX,IY)
         endif

         if (inflags1(3)) then
            WX0(IX,IY)  = x2w0%rattr(index_x2w_sa_u,gindex)
            WXN(IX,IY)  = WX0(IX,IY)
            WY0(IX,IY)  = x2w0%rattr(index_x2w_sa_v,gindex)
            WYN(IX,IY)  = WY0(IX,IY)
            DT0(IX,IY)  = x2w0%rattr(index_x2w_sa_tbot,gindex) - x2w0%rattr(index_x2w_so_t,gindex)
            DTN(IX,IY)  = DT0(IX,IY)
         endif

         if (inflags1(4)) then
            ICEI(IX,IY) = x2w0%rattr(index_x2w_si_ifrac,gindex)
         endif

         ! QL, 150525, get mixing layer depth from coupler
         ! SB, if (inflags1(5)) then
         ! SB,    HML(IX,IY) = max(x2w0%rattr(index_x2w_so_bldepth,gindex), 5.)
         ! SB, endif

      enddo
      enddo

      call mct_aVect_clean(x2w0)

      ! 7.b Run the wave model for the given interval

      !      write(stdout,*) 'wrm tcx7'
      !      call shr_sys_flush(stdout)
      !      call mpi_barrier(mpi_comm,ierr)
      !      write(stdout,*) 'wrm tcx7a'
      !      call shr_sys_flush(stdout)

      call w3wave ( 1, timen )

      ! copy ww3 data to coupling datatype
      ! QL, 150612, copy enhancement factor, uStokes, vStokes to coupler
      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         IX  = MAPSF(ISEA,1)
         IY  = MAPSF(ISEA,2)
         if (MAPSTA(IY,IX) .eq. 1) then
             w2x_w%rattr(index_w2x_Sw_ustokes,jsea) = USSX(ISEA)
             w2x_w%rattr(index_w2x_Sw_vstokes,jsea) = USSY(ISEA)
          else
             w2x_w%rattr(index_w2x_Sw_ustokes,jsea) = 0.0
             w2x_w%rattr(index_w2x_Sw_vstokes,jsea) = 0.0
          endif
      enddo

      !      write(stdout,*) 'wrm tcx8'
      !      call shr_sys_flush(stdout)

      !----------------------------------------------------------------------------
      ! Reset shr logging to original values
      !----------------------------------------------------------------------------
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)
      call shr_sys_flush(stdout)

      ! write(stdout,*) 'wrm tcx9'
      ! call shr_sys_flush(stdout)

      ! TODO Put in gptl timer calls

      ! Formats

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

    subroutine wav_domain_mct(lsize, gsmap, dom)

      implicit none
      integer        , intent(in)   :: lsize
      type(mct_gsmap), intent(in)   :: gsmap
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


      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         lon = x0 + real(ix-1)*sx
         data(jsea) = lon
         !write(stdout,*)' jsea= ',jsea,' lon is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lon",data,lsize)

      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         lat = y0 + real(iy-1)*sy
         data(jsea) = lat
         !write(stdout,*)' jsea= ',jsea,' lat is ',data(jsea)
      end do
      call mct_ggrid_importrattr(dom,"lat",data,lsize)

      do jsea = 1,nseal
         isea = iaproc + (jsea-1)*naproc
         ix = mapsf(isea,1)
         iy = mapsf(isea,2)
         lat = y0 + real(iy-1)*sy
         data(jsea) = sx*deg2rad*sy*deg2rad*cos(lat*deg2rad)
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

  END MODULE WAV_COMP_MCT
