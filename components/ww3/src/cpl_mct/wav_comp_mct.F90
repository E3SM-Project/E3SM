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
      use w3gdatmd
      use w3wdatmd, only: time, w3ndat, w3dimw, w3setw
      use w3adatmd
      use w3idatmd, only: flags, w3seti, w3ninp 
      USE W3IDATMD, ONLY: TC0, CX0, CY0, TCN, CXN, CYN
      USE W3IDATMD, ONLY: TW0, WX0, WY0, DT0, TWN, WXN, WYN, DTN
      USE W3IDATMD, ONLY: TIN, ICEI
      USE W3IDATMD, ONLY: TLN, WLEV
      use w3odatmd, only: w3nout, w3seto, naproc, iaproc, napout, naperr,             &
                          nogrd, idout, fnmpre, iostyp
!/
      use w3initmd
      use w3wavemd
!/
      use w3iopomd
      use w3timemd
      use w3cesmmd, only : casename

      use esmf
      use mct_mod 
      use seq_flds_mod
      use ww3_cpl_indices  , only :  ww3_cpl_indices_set
      use ww3_cpl_indices  , only :  &
         index_x2w_Sa_u, index_x2w_Sa_v, index_x2w_Sa_tbot, &
         index_x2w_Si_ifrac, &
         index_x2w_So_t, index_x2w_So_u, index_x2w_So_v, index_x2w_So_bldepth, &
         index_w2x_Sw_lamult, index_w2x_Sw_ustokes, index_w2x_Sw_vstokes, index_w2x_Sw_hstokes

      use shr_sys_mod      , only : shr_sys_flush 
      use shr_kind_mod     , only : in=>shr_kind_in, r8=>shr_kind_r8, &
                                    cs=>shr_kind_cs, cl=>shr_kind_cl
      use seq_cdata_mod    , only : seq_cdata, seq_cdata_setptrs
      use seq_timemgr_mod  , only : seq_timemgr_eclockgetdata
      use seq_infodata_mod , only : seq_infodata_type, seq_infodata_getdata, seq_infodata_putdata, &
                                    seq_infodata_start_type_start, seq_infodata_start_type_cont,   &
                                    seq_infodata_start_type_brnch
      use seq_comm_mct     , only : seq_comm_inst, seq_comm_name, seq_comm_suffix
      use shr_file_mod     , only : shr_file_setlogunit, shr_file_setloglevel, &
                                    shr_file_getlogunit, shr_file_getloglevel, &
                                    shr_file_getunit, shr_file_setio
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
      integer,save :: inst_index            ! number of current instance (ie. 1)
      character(len=16),save :: inst_name   ! fullname of current instance (ie. "wav_0001")
      character(len=16),save :: inst_suffix ! char string associated with instance

      include "mpif.h"
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CONTAINS
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SUBROUTINE WAV_INIT_MCT( EClock, cdata, x2w, w2x, NLFilename )

      !/ ------------------------------------------------------------------- /
      !/ Parameter list
      !/

      implicit none
      TYPE(ESMF_CLOCK), INTENT(INOUT) :: ECLOCK
      TYPE(SEQ_CDATA) , INTENT(INOUT) :: CDATA
      TYPE(MCT_AVECT) , INTENT(INOUT) :: X2W, W2X
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: NLFILENAME ! NAMELIST FILENAME

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

      character(CL)            :: starttype
      type(mct_gsmap), pointer :: gsmap
      type(mct_ggrid), pointer :: dom
      type(seq_infodata_type), pointer :: infodata   ! input init object

      integer             :: ndso, ndse, nds(13), ntrace(2), time0(2), &
                             timen(2), odat(30), nh(4), iprt(6)
      integer             :: i,j,npts
      integer             :: ierr
      real                :: a(nhmax,4)
      real, allocatable   :: x(:), y(:)
      logical             :: flgrd(nogrd), prtfrm, flt 
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
      call shr_file_setLogUnit (stdout)

      if ( iaproc .eq. napout ) write (ndso,900)
      call shr_sys_flush(stdout)

      !--------------------------------------------------------------------
      ! Define input fields
      !--------------------------------------------------------------------

      flags = .false.
!      flags(1:4) = .true.   !changed by Adrean (lev,curr,wind,ice on)
      flags(3:4) = .true.   !changed by Adrean (wind,ice on)

      !--------------------------------------------------------------------
      ! Set time frame
      !--------------------------------------------------------------------

      ! TIME0 = from ESMF clock
      ! NOTE - are not setting TIMEN here

      if ( iaproc .eq. napout ) write (ndso,930)
      call shr_sys_flush(stdout)

      call seq_timemgr_EClockGetData(EClock, &
           start_ymd=start_ymd, start_tod=start_tod)

      hh = start_tod/3600
      mm = (start_tod - (hh * 3600))/60
      ss = start_tod - (hh*3600) - (mm*60) 

      time0(1) = start_ymd
      time0(2) = hh*10000 + mm*100 + ss
      call stme21 ( time0 , dtme21 )
      if ( iaproc .eq. napout ) write (ndso,931) dtme21
      call shr_sys_flush(stdout)
      time = time0

      !--------------------------------------------------------------------
      ! Define output type and fields
      !--------------------------------------------------------------------

      iostyp = 1
      write (ndso,940) 'no dedicated output process, any file system '
      call shr_sys_flush(stdout)

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
      ! FLGRD   L.A.   I   Flags for gridded output.
      ! NPT     Int.   I   Number of output points
      ! X/YPT   R.A.   I   Coordinates of output points.
      ! PNAMES  C.A.   I   Output point names.

      npts = 0
      allocate ( x(1), y(1), pnames(1) )
      pnames(1) = ' '

      do j=1, 6
         odat(5*(j-1)+3) = 0
         !DEBUG
         ! Hardwire gridded output for now 
         odat(1) = 10101
         odat(2) = 0
         odat(3) = 1800   ! changed by Adrean
         odat(4) = 99990101
         odat(5) = 0
         !DEBUG
      end do
                                                                          
      ! Output Type 1: fields of mean wave parameters gridded output

      flgrd( 1) = .false. !   1. depth (m)                              
      flgrd( 2) = .false. !   2. mean current vel (vec, m/s)            
      flgrd( 3) = .true.  !   3. mean wind vel (vec, m/s)               
      flgrd( 4) = .true.  !   4. air-sea temp diff (deg C)              
      flgrd( 5) = .true.  !   5. skin friction vel (scalar, m/s)        
      flgrd( 6) = .true.  !   6. significant wave height (m)            
      flgrd( 7) = .false. !   7. mean wave length (m)                   
      flgrd( 8) = .true.  !   8. mean wave period (Tn1, s)              
      flgrd( 9) = .true.  !   9. mean wave dir (deg: met conv)          
      flgrd(10) = .false. !  10. mean dir spread (deg: )                
      flgrd(11) = .false. !  11. peak freq (Hz)                         
      flgrd(12) = .false. !  12. peak dir (deg: )                       
      flgrd(13) = .false. !  13. peak freq of wind-sea part             
      flgrd(14) = .false. !  14. wind-sea dir (deg: met conv)           
      flgrd(15) = .false. !  15. wave height of partitions               
      flgrd(16) = .false. !  16. peak periods of partitions             
      flgrd(17) = .false. !  17. peak wave length of partitions         
      flgrd(18) = .false. !  18. mean dir of partitions                 
      flgrd(19) = .false. !  19. dir spread of partitions               
      flgrd(20) = .false. !  20. wind-sea frac of partitions             
      flgrd(21) = .false. !  21. wind-sea frac of entire spec           
      flgrd(22) = .false. !  22. number of partitions                   
      flgrd(23) = .false. !  23. average time step (s)                  
      flgrd(24) = .false. !  24. cut-off freq (Hz)                      
      flgrd(25) = .true.  !  25. ice concentration (frac)               
      flgrd(26) = .false. !  26. water level (m?)                       
      flgrd(27) = .false. !  27. near-bottom rms exclusion amp          
      flgrd(28) = .false. !  28. near-bottom rms orbital vel            
      flgrd(29) = .false. !  29. radiation stresses                     
      flgrd(30) = .false. !  30. user defined (1)                       
      flgrd(31) = .false. !  31. user defined (2)                       

      if ( iaproc .eq. napout ) then
         flt = .true.
         do i=1, nogrd
            if ( flgrd(i) ) then
               if ( flt ) then
                  write (ndso,1945) idout(i)
                  flt    = .false.
               else
                  write (ndso,1946) idout(i)
               end if
            end if
         end do
         if ( flt ) write (ndso,1945) 'no fields defined'
      end if
      call shr_sys_flush(stdout)

      !--------------------------------------------------------------------
      ! Wave model initializations
      !--------------------------------------------------------------------

      if ( iaproc .eq. napout ) write (ndso,950)
      if ( iaproc .eq. napout ) write (ndso,951) 'wave model ...'
      call shr_sys_flush(stdout)

      call seq_infodata_GetData(infodata,case_name=casename)

      call w3init ( 1, 'ww3', nds, ntrace, odat, flgrd, npts, x, y,   &
           pnames, iprt, prtfrm, mpi_comm )
      call shr_sys_flush(stdout)

      ! overwrite dt values with variables from coupler
      ! is this a problem with any things being set in w3init?
      call seq_timemgr_eclockgetdata(eclock, dtime=dtime_sync )      
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
      write(stdout,*)'lsize= ',lsize
      call shr_sys_flush(stdout)

      ! initialize mct domain

      call wav_domain_mct(lsize, gsmap, dom)
      
      ! set flags in infodata
      ! wav_prognostic is set to .false. for debugging purposes only

      call seq_infodata_putdata(infodata, wav_present=.true., &
           wav_prognostic=.true., wav_nx=nx, wav_ny=ny)

      ! initialize mct attribute vectors
      
      call mct_avect_init(w2x, rlist=seq_flds_w2x_fields, lsize=lsize)
      call mct_avect_zero(w2x)
      
      call mct_avect_init(x2w, rlist=seq_flds_x2w_fields, lsize=lsize)
      call mct_avect_zero(x2w)

      ! add call to gptl timer

      ! end redirection of share output to wav log

      call shr_sys_flush(stdout)
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

!      write(stdout,*) 'wrm tcx1'
!      call shr_sys_flush(stdout)

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

      write(stdout,*)'time0= ',time0
      write(stdout,*)'timen= ',timen

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

!      def_value = -1.0e36
      def_value = 0.0

!      write(stdout,*) 'wrm tcx5'
!      call shr_sys_flush(stdout)

      if (flags(1)) then
         TLN  = timen
         WLEV = def_value   ! ssh
      endif

      if (flags(2)) then
         TC0  = time0
         TCN  = timen
         CX0  = def_value   ! ocn u current   
         CXN  = def_value   
         CY0  = def_value   ! ocn v current
         CYN  = def_value      
      endif

      if (flags(3)) then
         TW0  = time0
         TWN  = timen
         WX0  = def_value   ! atm u wind
         WXN  = def_value      
         WY0  = def_value   ! atm v wind
         WYN  = def_value      
         DT0  = def_value   ! air temp - ocn temp
         DTN  = def_value
      endif

      if (flags(4)) then
         TIN  = timen
         ICEI = def_value   ! ice frac
      endif

!      write(stdout,*) 'wrm tcx6'
!      call shr_sys_flush(stdout)

#if (1 == 0) 
! tcraig, this is the local fill only, really only need 1d
      do jsea=1, nseal
         isea = iaproc + (jsea-1)*naproc
         IX  = MAPSF(ISEA,1)
         IY  = MAPSF(ISEA,2)

         if (flags(1)) then
            WLEV(IX,IY) = 0.0
         endif

         if (flags(2)) then
            CX0(IX,IY)  = x2w_w%rattr(index_x2w_so_u,JSEA)
            CXN(IX,IY)  = CX0(IX,IY)
            CY0(IX,IY)  = x2w_w%rattr(index_x2w_so_v,JSEA)
            CYN(IX,IY)  = CY0(IX,IY)
         endif

         if (flags(3)) then
            WX0(IX,IY)  = x2w_w%rattr(index_x2w_sa_u,JSEA)
            WXN(IX,IY)  = WX0(IX,IY)
            WY0(IX,IY)  = x2w_w%rattr(index_x2w_sa_v,JSEA)
            WYN(IX,IY)  = WY0(IX,IY)
            DT0(IX,IY)  = x2w_w%rattr(index_x2w_sa_tbot,JSEA) - x2w_w%rattr(index_x2w_so_t,JSEA)
            DTN(IX,IY)  = DT0(IX,IY)
         endif

         if (flags(4)) then
            ICEI(IX,IY) = x2w_w%rattr(index_x2w_si_ifrac,JSEA)
         endif

      enddo

#else
! tcraig, this is the global fill
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

         if (flags(1)) then
            WLEV(IX,IY) = 0.0
         endif

         if (flags(2)) then
            CX0(IX,IY)  = x2w0%rattr(index_x2w_so_u,gindex)
            CXN(IX,IY)  = CX0(IX,IY)
            CY0(IX,IY)  = x2w0%rattr(index_x2w_so_v,gindex)
            CYN(IX,IY)  = CY0(IX,IY)
         endif

         if (flags(3)) then
            WX0(IX,IY)  = x2w0%rattr(index_x2w_sa_u,gindex)
            WXN(IX,IY)  = WX0(IX,IY)
            WY0(IX,IY)  = x2w0%rattr(index_x2w_sa_v,gindex)
            WYN(IX,IY)  = WY0(IX,IY)
            DT0(IX,IY)  = x2w0%rattr(index_x2w_sa_tbot,gindex) - x2w0%rattr(index_x2w_so_t,gindex)
            DTN(IX,IY)  = DT0(IX,IY)
         endif

         if (flags(4)) then
            ICEI(IX,IY) = x2w0%rattr(index_x2w_si_ifrac,gindex)
         endif

      enddo
      enddo

      call mct_aVect_clean(x2w0)
#endif

      ! 7.b Run the wave model for the given interval

!      write(stdout,*) 'wrm tcx7'
!      call shr_sys_flush(stdout)
!      call mpi_barrier(mpi_comm,ierr)
!      write(stdout,*) 'wrm tcx7a'
!      call shr_sys_flush(stdout)

      call w3wave ( 1, timen )

!      copy ww3 data to coupling datatype
!      do jsea=1, nseal
!         isea = iaproc + (jsea-1)*naproc
!         IX  = MAPSF(ISEA,1)
!         IY  = MAPSF(ISEA,2)
!         w2x_w%rattr(index_w2x_Sw_lamult,jsea) = ??
!         w2x_w%rattr(index_w2x_Sw_ustokes,jsea) = ??
!         w2x_w%rattr(index_w2x_Sw_vstokes,jsea) = ??
!         w2x_w%rattr(index_w2x_Sw_hstokes,jsea) = ??
!      enddo

!      write(stdout,*) 'wrm tcx8'
!      call shr_sys_flush(stdout)

      !----------------------------------------------------------------------------
      ! Reset shr logging to original values
      !----------------------------------------------------------------------------
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)
      call shr_sys_flush(stdout)

!      write(stdout,*) 'wrm tcx9'
!      call shr_sys_flush(stdout)

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
         mask = mapsta(iy,ix)
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



