
      program mozart_pp
!-----------------------------------------------------------------------
!        ... Mozart chemistry pre-processor
!-----------------------------------------------------------------------

      use io
      use elements
      use mass_diags
      use var_mod
      use rxt_mod
      use lin_matrix
      use nln_matrix
      use lu_factor
      use lu_solve
      use prod_loss
      use ind_prod
      use set_rxt_rates
      use mo_spat_dims
      use mo_ver_hdr
      use mo_files_hdr
      use mo_hist_out
      use mo_chem,     only : chem
      use sp_mods,     only : sparsity
      use mo_ver_opts, only : ver_opts
      use rxt_equations_mod

      implicit none

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, target  ::  nind(200)
      integer, pointer, dimension(:) ::  nbeg, nend
      integer  ::  dimensions(6) = (/ 128, 64, 18, 1, 1, 32 /)

      integer  ::  plon, plonl, plat, plev          ! spatial dimensions of simulation
      integer  ::  jintmx, nxpt                     ! slt parameters for bounds and array "padding"
      integer  ::  sub_cnt = 0                      ! count of user subroutines

      integer  ::  class, clsndx
      integer  ::  grp_rows, rel_rows
!-----------------------------------------------------------------------
!        ... Iteration counts are as follows:
!            (1) == "hov" iteration count
!            (2) == "implicit" iteration count
!            (3) == "implicit" jacobian update count( first count iterations)
!-----------------------------------------------------------------------
      character(len=1), parameter :: on = 'y'
      character(len=1), parameter :: off = 'n'

      integer  ::  iter_counts(4) = (/ 7, 4, 2, 5 /)

      character(len=320) :: lib_src(350)
      character(len=320) :: chem_src(50)
      character(len=320) :: filename(100)
      character(len=320) :: filepath(100)
      character(len=320) :: sub_names(100)
      character(len=80)  :: iout(100)
      character(len=320) :: mod_names(100)
      character(len=320) :: mod_paths(100)
      character(len=320) :: mod_src(100)
      character(len=64)  :: histinp(4)
      character(len=64)  :: histout(6)
      character(len=16)  :: jobctl(8)
      character(len=16)  :: wrk_rxt(10)
      character(len=10)  :: clshdr(5)
      character(len=16)   :: wrk_chr(10)
      character(len=580) :: command, cpp_command
      character(len=256) :: errcom, filout, filin
      character(len=64)  :: oper_flpth
      character(len=64)  :: cpp_dir, cpp_opts
      character(len=64)  :: wrk_dir
      character(len=64)  :: tar_flnm
      character(len=64)  :: subfile
      character(len=64)  :: filenm
      character(len=64)  :: tmp_filenm
      character(len=16)  :: param
      character(len=16)  :: hostname
      character(len=16)  :: jobname
    
      character(len=16) ::  machine   = 'IBM'
      character(len=16) ::  march     = 'SCALAR'
      character(len=16) ::  model     = 'CAM'
      character(len=16) ::  arch_type = 'HYBRID'
      character(len=16) ::  char
      character(len=4) ::  ftunit = 'ft'
      character(len=1) ::  errflg

      integer, allocatable :: mask(:)
      integer ::  entry(11)
      integer ::  filelines(5)
      integer ::  dyn_hst_fld_cnt(2)                        ! multi and single level field count
      integer ::  ratind(2)
      integer ::  additions(2)
      integer ::  multiplications(2)
      integer ::  nzcnt(2)  = 0                             ! number of nonzero entries in lu
      integer ::  file_cnt, hst_file_cnt
      integer ::  nchar, k, noff, m, j, l, il, iu
      integer ::  i, indx, ntab, ios, astat
      integer ::  spcno, counter, rxno, col, retcod, length, place
      integer ::  fixrows, rxmrows, pcelrows, pceprows
      integer ::  cpucnt = 1                                ! number of cpu's
      
      logical, target      ::  options(20)
      logical, allocatable :: lin_mat_pat(:)
      logical, pointer     :: usemods
      logical ::  null_flag
      logical ::  found
      logical ::  lexist
      logical ::  vec_ftns = .false.            ! vector functions
      logical ::  radj_flag
      logical ::  ohstflag                      ! output history tape flag
      logical ::  diagprnt = .false.            ! chktrc or negtrc diagnostics printout
      logical ::  tavgprnt = .false.            ! time averaged printout
      logical ::  longnames = .false.           ! do not use long names

      integer :: rxt_tag_cnt
      character(len=32) :: rxt_rates_conv_file =  'mo_rxt_rates_conv.F90'

      type(SPARSITY) :: sparse(2)

!----------------------------------------------------------------------------------
!	... Function declarations
!----------------------------------------------------------------------------------
      integer  ::   LENOF
      integer  ::   STRLEN
      
!----------------------------------------------------------------------------------
!       ... The options array has the following mapping:
!
!       (1) Chemistry (on/off)            (2) Target machine == cray (yes/no)
!       (3) Diffusion (on/off)            (4) Convection (on/off)
!       (5) Iter norms (on/off)           (6) Conservation (on/off)
!       (7) Source code (yes/no)          (8) Submission files (yes/no)
!       (9) Execution (yes/no)           (10) SLT fixer (on/off)
!      (11) Multitasking (yes/no)        (12) Rxt rate lookup table (on/off)
!      (13) Relative humidity (yes/no)   (14) New compiler (yes/no)
!      (15) Height field (yes/no )       (16) User "hook" (yes/no)
!      (17) Use f90 modules (yes/no)     (18) Make and use f90 names module (yes/no)
! (19 - 20) Unused
!
!           Iter norms, Execution, and Rxt rate lookup default to off
!----------------------------------------------------------------------------------
      data options / 4*.true., .false., 3*.true., .false., 2*.true., 2*.false., .true., 6*.false. /
      data clshdr / 'Explicit', 'Ebi', 'Hov', 'Implicit', 'Rodas' /

!-----------------------------------------------------------------------
!        ... Initialize pointers and data
!-----------------------------------------------------------------------
      usemods => options(17)
      nbeg => nind(1:100)
      nend => nind(101:200)
      jobctl(:) = ' '
      histout(:) = ' '
      histinp(:) = ' '
      histinp(4) = 'LONG'
      wrk_dir = '$TMPDIR'
      subfile = ' '
      entry(:)  = 0
      filelines(:) = 0
      src_dir = '../bkend/'

!-----------------------------------------------------------------------
!        ... Set default filenames/paths
!-----------------------------------------------------------------------
      output_path = '../output/'
      input_path  = '../input/'
      temp_path   = '../tmp/'
      sim_dat_path = output_path
      procout_path = output_path
      procfiles_path   = '../procfiles/'
      sim_dat_filename = 'sim.dat'
      sim_dat_filespec = trim(sim_dat_path) // 'sim.dat'
      cpp_dir  = ' '
      cpp_opts = '-P -C -I.'

!-----------------------------------------------------------------------
!        ... Assign default input, output units
!-----------------------------------------------------------------------
      lin  = 5
      lout = 6

!-----------------------------------------------------------------------
!        ... Get arguments
!-----------------------------------------------------------------------
      filin  = ' '
      filout = ' '
      call getarg( 1, filin )
      call getarg( 2, filout )

!-----------------------------------------------------------------------
!        ... No input filespec on command line; request input
!-----------------------------------------------------------------------
      if( filin == ' ' ) then
         write(*,'('' Enter filespec of input file'')')
         read(*,'(a80)') filin
         if( filin == ' ' ) then
            filin = './mozart2.inp'
         end if
      end if
      open( unit = 5, &
            file = trim( filin ), &
            status = 'old', &
            iostat = ios )
      if( ios /= 0 ) then
	 write(*,*) ' Failed to open file ',trim( filin )
	 write(*,*) ' Error code = ',ios
	 stop
      end if

      call cardin( lin, buff, nchar )
      buffh = buff
      call upcase( buffh )

!-----------------------------------------------------------------------
!        ... Check for input overide and process if present
!            if input unit 5 is overriden take 
!            all simulation input from lin
!-----------------------------------------------------------------------
      if( buffh(:18) == 'INPUT_UNIT_NUMBER=' ) then
         errflg = on
         call intcon( buff(19:nchar), &
                      nchar - 18, &
                      lin, &
                      retcod )
         if( retcod /= 0 ) then
            errcom = buff(19:nchar) // ' is an invalid unit number@'
         else
            if( lin <= 0 ) then
               errcom = buff(19:nchar) // ' is an invalid unit number@'
            else if( lin <= 3 ) then
               errcom = buff(19:nchar) // ' is a reserved unit number@'
            else if( lin == 6 ) then
               errcom = buff(19:nchar) // ' is a reserved unit number@'
            else if( lin >= 100 ) then
               errcom = buff(19:nchar) // ' is an invalid unit number@'
            else
               errflg = off
            end if
         end if

         if( errflg == on ) then
            call errmes( errcom, 6, param, 8, buff )
         end if

!-----------------------------------------------------------------------
!        ... Check for input file override and process if present
!-----------------------------------------------------------------------
         call cardin ( 5, buff, nchar )
         buffh = buff
         call upcase( buffh )
         if( buffh(:15) == 'INPUT_FILESPEC=' ) then
            filin = buff(16:nchar)
            if( lin <= 10 ) then
               write (ftunit(3:4),'(''0'',i1)') lin
            else
               write (ftunit(3:4),'(i2)') lin
            end if
            close( unit = 5 )
            open( unit = lin, &
                  file = trim( filin ), &
                  status = 'old', &
                  iostat = ios )
            if( ios /= 0 ) then
	       write(*,*) ' Failed to open file ',trim( filin )
	       write(*,*) ' Error code = ',ios
	       stop
            end if
            call cardin( lin, buff, nchar )
            buffh = buff
            call upcase( buffh )
         else
            call errmes( ' ** Input reassignment requires both a unit number and filename@', &
                         lout, &
                         param, &
                         8, &
                         buff )
         end if
      end if

!-----------------------------------------------------------------------
!        ... Check for simulation start card (begsim)
!-----------------------------------------------------------------------
      if( buffh /= 'BEGSIM' ) then
         call errmes ( ' ** first card not begsim **@', &
                       lout, &
                       param, &
                       8, &
                       buff )
      end if

      call cardin ( lin, buff, nchar )
      buffh = buff
      call upcase( buffh )

!-----------------------------------------------------------------------
!        ... Check for doc file overide and process if present
!-----------------------------------------------------------------------
      if( buffh(:19) == 'OUTPUT_UNIT_NUMBER=' ) then
         errflg = on
         call intcon( buff(20:nchar), &
                      nchar - 19, &
                      lout, &
                      retcod )
         if( retcod /= 0 ) then
            errcom = buff(20:nchar) // ' is an invalid unit number@'
         else
            if( lout <= 0 ) then
               errcom = buff(20:nchar) // ' is an invalid unit number@'
            else if( lout <= 3 ) then
               errcom = buff(20:nchar) // ' is a reserved unit number@'
            else if( lout == lin .or. lout == 6 ) then
               errcom = buff(20:nchar) // ' is a reserved unit number@'
            else if( lout >= 100 ) then
               errcom = buff(20:nchar) // ' is an invalid unit number@'
            else
               errflg = off
            end if
         end if

!-----------------------------------------------------------------------
!        ... Error in assigning output unit number
!-----------------------------------------------------------------------
         if( errflg == on ) then
            call errmes ( errcom, 6, param, 8, buff )
         end if

!-----------------------------------------------------------------------
!        ... Set the output unit number
!-----------------------------------------------------------------------
         if( lout <= 10 ) then
            write (ftunit(3:4),'(''0'',i1)') lout
         else
            write (ftunit(3:4),'(i2)') lout
         end if

         call cardin ( lin, buff, nchar )
         buffh = buff
         call upcase( buffh )

!-----------------------------------------------------------------------
!        ... get output_file
!-----------------------------------------------------------------------
         if( buffh(:12) == 'OUTPUT_FILE=' ) then
            filout = buff(13:nchar)
            call cardin( lin, buff, nchar )
            buffh = buff
            call upcase( buffh )
         else
            call errmes( ' ** Output reassignment requires both a unit number and filename@', &
                         lout, &
                         param, &
                         8, &
                         buff )
         end if
      else
!-----------------------------------------------------------------------
!        ... Assign output unit
!-----------------------------------------------------------------------
         if( filout == ' ' ) then
            write(*,'('' Enter filename of output file'')')
            read(*,'(a80)') filout
            if( filout == ' ' ) then
               filout = 'mozart2.doc'
            end if
            filout = trim(output_path) // trim( filout )
         end if
         open( unit   = lout, &
               file   = trim( filout ), &
               status = 'new', &
               iostat = ios )
         if( ios /= 0 ) then
	    write(*,*) ' Failed to open file ',trim(filout)
	    write(*,*) ' Error code = ',ios
	    stop
         end if
      end if

      do
!-----------------------------------------------------------------------
!        ... Check for procout path override
!-----------------------------------------------------------------------
         if( buffh(:13) == 'PROCOUT_PATH=' ) then
            procout_path = buff(14:nchar)
!-----------------------------------------------------------------------
!        ... Check for output path override
!-----------------------------------------------------------------------
         else if( buffh(:15) == 'PROCFILES_PATH=' ) then
            procfiles_path = buff(16:nchar)
         else if( buffh(:12) == 'OUTPUT_PATH=' ) then
            output_path = buff(13:nchar)
!-----------------------------------------------------------------------
!        ... Check for tmp_path override
!-----------------------------------------------------------------------
         else if( buffh(:10) == 'TEMP_PATH=' ) then
             temp_path= buff(11:nchar)
!-----------------------------------------------------------------------
!        ... Check for sim_dat_path override
!-----------------------------------------------------------------------
         else if( buffh(:13) == 'SIM_DAT_PATH=' ) then
            sim_dat_path = buff(14:nchar)
!-----------------------------------------------------------------------
!        ... Check for src_path override
!-----------------------------------------------------------------------
         else if( buffh(:9) == 'SRC_PATH=' ) then
            src_dir = buff(10:nchar)
!-----------------------------------------------------------------------
!        ... Check for sim_dat_filename override
!-----------------------------------------------------------------------
         else if( buffh(:17) == 'SIM_DAT_FILENAME=' ) then
            sim_dat_filename = buff(18:nchar)
	    sim_dat_filespec = trim( sim_dat_path ) // trim( sim_dat_filename )
	    inquire( file = trim( sim_dat_filespec ), exist = lexist )
	    if( lexist ) then
	       call system( 'rm -f ' // trim( sim_dat_filespec ) )
	    end if
	 else
	    exit
         end if
         call cardin ( lin, buff, nchar )
         buffh = buff
         call upcase( buffh )
      end do

!-----------------------------------------------------------------------
!        ... Find cpp preprocessor
!-----------------------------------------------------------------------
      command = 'whereis cpp > ' // trim(temp_path) // 'cpp.path'
      call system( trim( command ) )
      open( unit=20, file=trim(temp_path)//'cpp.path', iostat=ios )
      if( ios /= 0 ) then
	 write(*,*) ' Failed to locate cpp path'
	 stop
      end if
      read(20,'(a)',iostat=ios) iout(1)
      if( ios /= 0 ) then
	 write(*,*) ' Failed to read cpp path'
	 stop
      end if
      if( iout(1)(1:4) == 'cpp:' ) then
	 nbeg(1) = index( trim(iout(1)), ' ' ) + 1
	 nend(1) = index( trim(iout(1)(nbeg(1):)), ' ' )
	 nend(1) = nbeg(1) + nend(1) - 1
	 cpp_dir = iout(1)(nbeg(1):nend(1))
	 close( 20 )
      else
	 write(*,*) ' Failed to locate cpp path'
	 stop
      end if

!-----------------------------------------------------------------------
!        ... Check for output file override and process if present
!-----------------------------------------------------------------------
      filout = trim(output_path) // filout 
      open( unit   = lout, &
           file   = trim( filout ), &
           status = 'replace', &
           iostat = ios )
      if( ios /= 0 ) then
         write(*,*) ' Failed to open file ',trim(filout)
         write(*,*) ' Error code = ',ios
         stop
      end if

!-----------------------------------------------------------------------
!        ... Check for comments and process if present
!-----------------------------------------------------------------------
      if( trim(buffh) == 'COMMENTS' ) then
         k = 1
         do 
            call cardin( lin, buff, nchar )
            buffh = buff
            call upcase( buffh )
            if( buffh == 'ENDCOMMENTS' ) then
               exit
            end if
	    iout(k) = buff
            k = k + 1
         end do

         k = k - 1
         noff = 100

         do m = 1,k
            buff = iout(m)
            do j = 1,80
               if( buff(j:j) /= ' ' ) then
                  exit
               end if
            end do
            l = j
            do j = 80,l,-1
               if( buff(j:j) /= ' ' ) then
                  nchar = j - l + 1
                  nchar = 40 - nchar/2
                  nbeg(m) = l
                  nend(m) = j
                  noff = MIN( nchar, noff )
                  exit
               end if
            end do
         end do

         do m = 1,k
            buff = iout(m)
            iout(m) = ' '
            il = nbeg(m)
            if( il /= 0 ) then
               iu = nend(m)
               iout(m)(noff:) = buff(il:iu)
            end if
         end do

!-----------------------------------------------------------------------
!        ... Write out the comments
!-----------------------------------------------------------------------
         write(lout,*) ' '
         write(lout,*) ' '
         write(lout,1565)
         write(lout,1571)
         write(lout,1571)
         write(lout,1567) (iout(m),m = 1,k)
         write(lout,1571)
         write(lout,1571)
         write(lout,1565)
         do m = 1,k
            iout(m) = ' '
         end do
         call cardin ( lin, &
                       buff, &
                       nchar )
         buffh = buff
         call upcase( buffh )
      end if

!-----------------------------------------------------------------------
!        ... Clean the temp work directory
!-----------------------------------------------------------------------
      call system( 'rm -f ' // trim( temp_path ) // '*' )
!-----------------------------------------------------------------------
!        ... Initialize the variables
!-----------------------------------------------------------------------
      call var_ini
!-----------------------------------------------------------------------
!        ... Initialize the reactions
!-----------------------------------------------------------------------
      call rxt_ini
!-----------------------------------------------------------------------
!        ... The species symbol list processing
!-----------------------------------------------------------------------
      call symbol( iout )
      ntab = maxval( spccnt(1:5) )

!-----------------------------------------------------------------------
!       ... Get variable mass
!-----------------------------------------------------------------------
      call iniele
      do i = 1,spccnt(1)
	 if( aliases(i) /= ' ' ) then
	    mass(i) = com_mass( aliases(i) )
	    c_mass(i) = com_mass( aliases(i),carbmass=.true. )
	 else
	    mass(i) = com_mass( solsym(i) )
	    c_mass(i) = com_mass( solsym(i),carbmass=.true. )
	 end if
      end do
      do i = spccnt(1)+1, spccnt(1)+spccnt(3)
	 if( aliases(i) /= ' ' ) then
	    mass(i) = com_mass( aliases(i) )
	    c_mass(i) = com_mass( aliases(i),carbmass=.true. )
	 else
	    mass(i) = com_mass( fixsym(i-spccnt(1)) )
	    c_mass(i) = com_mass( fixsym(i-spccnt(1)),carbmass=.true. )
	 end if
      end do

!-----------------------------------------------------------------------
!       ... Form individual group members
!-----------------------------------------------------------------------
      i = 0
      do l = 1,ngrp
         do k = 1,grpcnt(l)
            j = grpmap(k,l)
            if( j <= 1999 ) then
               j = j - 1000
               i = i + 1
	       grp_mem_sym(i) = solsym(j)
            end if
         end do
      end do
      grp_mem_cnt = i

!-----------------------------------------------------------------------
!        ... Now begin group modification process by making
!            new species numbering and group association tables
!-----------------------------------------------------------------------
      counter = 0
      do i = 1,ngrp
         do j = 1,grpcnt(i)
            indx = mod( grpmap(j,i),1000 )
            grpflg(indx) = i
            counter = counter + 1
            mem2grp_map(counter) = i
            grp_rat_ind(indx) = counter
         end do
         new_solsym(i) = grpsym(i)
      end do

      do i = 1,relcnt
         indx = relmap(i,1)
         rel_flg(indx) = i
      end do
      
      indx = ngrp
      do i = 1,nq
         if( grpflg(i) /= 0 .or. rel_flg(i) /= 0 ) then
            cycle
         else
            indx = indx + 1
            newind(i) = indx
            new_solsym(indx) = solsym(i)
         end if
      end do
      new_nq = indx
      
      write(lout,*) ' '
      write(lout,*) ' '
      write(lout,230)
      do j = 1,nq
	 if( aliases(j) == ' ' ) then
            write(lout,231) j, solsym(j)
	 else
            write(lout,'(6x,''('',i3,'')'',2x,a16,3x,''('',a,'')'')') j, solsym(j), trim( aliases(j) )
	 end if
      end do
      if( relcnt /= 0 ) then
         write(lout,235)
         do j = 1,relcnt
            length = STRLEN( solsym(relmap(j,1)) )
            buff = ' '
            write(buff,'(6x,''('',i2,'')'')') j
            buff(STRLEN(buff)+2:) = solsym(relmap(j,1))(:length) // ' ~ ' // solsym(relmap(j,2))
            write(lout,'(a)') buff(:STRLEN(buff))
         end do
      end if
      if( nfs /= 0 ) then
         write(lout,*) ' '
         write(lout,*) ' '
         write(lout,232)
         write(lout,231) (j, fixsym(j), j = 1,nfs)
      end if
      if( ncol /= 0 ) then
         write(lout,*) ' '
         write(lout,*) ' '
         write(lout,236)
         write(lout,238) (j, colsym(j), colub(j), j = 1,ncol)
      end if
      if( ngrp /= 0 ) then
         write(lout,*) ' '
         write(lout,*) ' '
         write (lout,237)
         write (lout,'(2x,''('',i2,'')'',2x,a80)') ( j, iout(j), j = 1,ngrp )
         do j = 1,ngrp
           iout(j) = ' '
         end do
      end if

!-----------------------------------------------------------------------
!        ... Write out group modified species list
!-----------------------------------------------------------------------
      if( ngrp /= 0 ) then
         write(lout,*) ' '
         write(lout,*) ' '
         write(lout,*) 'Advected species'
         write(lout,231) (j, new_solsym(j), j = 1,new_nq)
      end if

!-----------------------------------------------------------------------
!        ... Define the solution classes
!-----------------------------------------------------------------------
      call sol_cls( iout )
     
!-----------------------------------------------------------------------
!        ... Write out class lists
!-----------------------------------------------------------------------
      write(lout,*) ' '
      write(lout,'(''Class List'')')
      write(lout,'(''=========='')')
      do k = 1,5
         if( clscnt(k) /= 0 ) then
            if( k > 1 ) then
               write(lout,*) ' '
            end if
            write(lout,'(1x,a10)') clshdr(k)
            write(lout,'('' --------'')')
            write(lout,231) (j, new_solsym(clsmap(j,k,2)), j = 1,clscnt(k) )
         end if
      end do
!-----------------------------------------------------------------------
!        End of the variable list processing
!------------------------------------------------------------------------

!=======================================================================
!        ... Chemistry processing
!=======================================================================
      call cardin( lin, buff, nchar )
      call upcase( buff )
      if( buff == 'CHEMISTRY' ) then
!-----------------------------------------------------------------------
!        ... Set the reactions and rates
!-----------------------------------------------------------------------
         call chem
         options(1) = .true.
         gascnt = rxntot - phtcnt
         do i = 1,rxntot
            if( rxt_tag(i) /= ' ' ) then
               rxt_has_tag(i) = .true.
            end if
         end do
         call cardin( lin, &
                      buff, &
                      nchar )
         call upcase( buff )
      else
         options(1) = .false.
      end if

!-----------------------------------------------------------------------
!        ... Transform the "hetero" reaction map
!            The 1st column is the new species number
!            if the species is a group member then the second column
!            indicates the species number within the group ( the 1st col)
!-----------------------------------------------------------------------
      do j = 1,hetcnt
         spcno = hetmap(j,1)
         if( grpflg(spcno) /= 0 ) then
            hetmap(j,1) = grpflg(spcno)
            hetmap(j,2) = grp_rat_ind(spcno)
         else
            hetmap(j,1) = newind(spcno)
            hetmap(j,2) = 0
         end if
      end do

!-----------------------------------------------------------------------
!        ... Then the "extraneous" reaction map
!-----------------------------------------------------------------------
      do j = 1,usrcnt
         spcno = usrmap(j)
         if( grpflg(spcno) /= 0 ) then
            usrmap(j) = grpflg(spcno)
         else
            usrmap(j) = newind(spcno)
         end if
      end do
     
!=======================================================================
!        ... The run parameters processing section
!=======================================================================
      if( buff == 'ENDSIM' ) then
         go to 292
      else if( buff == 'SIMULATIONPARAMETERS' ) then
         do
            call cardin( lin, &
                         buff, &
                         nchar )
            call upcase( buff )
            if( buff == 'SPATIALDIMENSIONS' ) then
               if( entry(1) /= 0 ) then
                  call errmes( ' spatial dimensions already' &
                            // '  prescribed@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(1) = 1
                  call spat_dims( buff, dimensions )
                  plon = dimensions(1)
                  plat = dimensions(2)
                  plev = dimensions(3)
                  nxpt = dimensions(4)
                  jintmx = dimensions(5)
                  plonl = dimensions(6)
               end if
            else if( buff == 'MASSDIAGNOSTICS' ) then
               if( entry(11) /= 0 ) then
                  call errmes( ' Mass diagsnostics already prescribed@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else if( entry(1) == 0 ) then
                  call errmes( ' Spatial dimensions must be done before mass diags@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(11) = 1
                  call mass_diagnostics( spcsym, &
                                         spccnt, &
                                         plon, plev, plat )
               end if
            else if( buff == 'VERSIONOPTIONS' ) then
               if( entry(2) /= 0 ) then
                  call errmes( ' Version options already prescribed@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(2) = 1
                  call ver_opts( options(2:), model, machine, march, arch_type, &
                                 wrk_dir, cpp_dir, cpp_opts, subfile, diagprnt, &
                                 tavgprnt, cpucnt, vec_ftns )
!-----------------------------------------------------------------------
!        ... Write out the species and reaction id files
!-----------------------------------------------------------------------
		  if( options(18) ) then
                     call make_name_mod
                     call make_rxt_name_mod
                     !call make_het_name_mod
                  end if
               end if
            else if( buff == 'EXECUTIONOPTIONS' ) then
               if( entry(4) /= 0 ) then
                  call errmes( ' Exec options already prescribed@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(4) = 1
                  call exe_opts( options(8), &
                                 lin, &
                                 lout )
               end if
            else if( buff == 'USERSUBROUTINES' ) then
               if( entry(3) /= 0 ) then
                  call errmes( ' Subroutines already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(3) = 1
                  call usrsubs( sub_names, &
                                sub_cnt )
!-----------------------------------------------------------------------
!       ... Parse user file pathnames
!-----------------------------------------------------------------------
                  if( sub_cnt /= 0 ) then
                     do i = 1,sub_cnt
                        call parse_flpth( sub_names(i), &
                                          filename(i), &
                                          filepath(i) )
                     end do
                     do i = 1,sub_cnt
	                if( index( filename(i), '.mod', back = .true. ) /= 0 ) then
			   options(17) = .true.                 ! force fortran90
			   usemods = .true.
			   exit
			end if
                     end do
                  end if
               end if
            else if( buff == 'JOBCONTROL' ) then
               if( entry(5) /= 0 ) then
                  call errmes( ' Job control already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(5) = 1
                  call job_ctl(  lin, &
                                 lout, &
                                 jobctl )
               end if
            else if( buff == 'NUMERICALCONTROL' ) then
               if( entry(10) /= 0 ) then
                  call errmes( ' Numerical control already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(10) = 1
                  call num_ctl( iter_counts )
               end if
            else if( buff == 'INPUTS' ) then
               if( entry(6) /= 0 ) then
                  call errmes( ' Inputs already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(6) = 1
                  call hist_inp( lin, &
                                 lout, &
                                 histinp, &
				 dyn_hst_fld_cnt )
               end if
            else if( buff == 'OUTPUTS' ) then
               if( entry(7) /= 0 ) then
                  call errmes( ' Outputs already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(6) = 1
                  call hist_out( histout, longnames, hst_file_cnt )
               end if
            else if( buff == 'SURFACEFLUX' ) then
               if( entry(8) /= 0 ) then
                  call errmes( ' Surf flux already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(8) = 1
                  call srfflx( lin, &
                               lout, &
                               new_nq, &
                               new_solsym, &
                               srf_flx_map, &
                               srf_flx_cnt, &
			       1 )
               end if
            else if( buff == 'SURFACEDEPOSITION' ) then
               if( entry(9) /= 0 ) then
                  call errmes( ' Surf Dep already specified@', &
                               lout, &
                               char, &
                               1, &
                               buff )
               else
                  entry(9) = 1
                  call srfflx( lin, &
                               lout, &
                               new_nq, &
                               new_solsym, &
                               dvel_map, &
                               dvel_cnt, &
			       2 )
               end if
            else if( buff == 'ENDSIMULATIONPARAMETERS' ) then
               exit
            else if( buff /= 'ENDPAR' ) then
               call errmes ( ' endsim card missing@', &
                             lout, &
                             char, &
                             1, &
                             buff )
            end if
         end do
      else
         call errmes ( ' endsim card missing@', &
                       lout, &
                       char, &
                       1, &
                       buff )
      end if
292   continue

Has_chemistry: &
      if( options(1) ) then         ! do only if there is chemistry
!=======================================================================
!        ... Weed out the proportional products in all reaction maps
!=======================================================================
         allocate( mask(max(7,prd_limp1)),stat=astat )
	 if( astat /= 0 ) then
	    write(lout,*) 'Failed to allocate the mask array; error = ',astat
	    stop 'abort'
	 end if
!-----------------------------------------------------------------------
!        ... First the "independent" production map
!-----------------------------------------------------------------------
         do j = 1,prdcnt
            do k = 2,prd_limp1
               spcno = prdmap(j,k)
               if( spcno == 0 ) then
                  mask(k) = -1
                  exit
               else if( rel_flg(spcno) == 0 ) then
                  mask(k) = 1
               else
                  mask(k) = 0
               end if
            end do
            place = 1
            do k = 2,prd_limp1
               if( mask(k) == -1 ) then
                  prdmap(j,place+1:prd_limp1) = 0
                  exit
               else if( mask(k) == 1 ) then
                  place = place + 1
                  prdmap(j,place) = prdmap(j,k)
               end if
            end do
         end do

!-----------------------------------------------------------------------
!        ... Then the "regular" reaction map
!-----------------------------------------------------------------------
         do i = 1,2
            do j = 1,rxmcnt(i)
               do k = i+2,i+prd_limp1
                  spcno = rxmap(j,k,i)
                  if( spcno == 0 ) then
                     mask(k) = -1
                     exit
                  else if( rel_flg(spcno) == 0 ) then
                     mask(k) = 1
                  else
                     mask(k) = 0
                  end if
               end do
               place = i + 1
               do k = i+2,i+prd_limp1
                  if( mask(k) == -1 ) then
                     rxmap(j,place+1:i+prd_limp1,i) = 0
                     exit
                  else if( mask(k) == 1 ) then
                     place = place + 1
                     rxmap(j,place,i) = rxmap(j,k,i)
                  end if
               end do
            end do
         end do

!-----------------------------------------------------------------------
!        ... Now xform all "proportional" reactants to proportional species
!        NOTE! The proportional reactants are replaced by the NEGATIVE
!              index of the species they are proportional to
!-----------------------------------------------------------------------
         do i = 1,2
            do j = 1,rxmcnt(i)
               counter = 0
               do k = 2,i+1                                     ! only do the reactants
                  spcno = rxmap(j,k,i)
                  if( rel_flg(spcno) /= 0 ) then
                     counter = counter + 1
                     rxmap(j,k,i) = -relmap(rel_flg(spcno),2)
                     ratind(counter) = spcno
                  end if
               end do
               if( counter /= 0 ) then
                  rxno = rxmap(j,1,i)
                  rel_rxt_cnt(counter) = rel_rxt_cnt(counter) + 1
                  indx = rel_rxt_cnt(counter)
                  rel_rxt_map(indx,1,counter) = rxno   !the reaction number
                  do l = 1,counter
                     rel_rxt_map(indx,l+1,counter) = rel_flg(ratind(l))
                  end do
                  rel_rxt_pntr(rxno,1) = counter
                  rel_rxt_pntr(rxno,2) = indx
               end if
            end do
         end do

!-----------------------------------------------------------------------
!        Now do the actual reaction matrix transforms
!        The first phase just does the basic x-form.
!        The second phase scans resultant maps to "eliminate"
!        matching product and reactant species in the
!        same reaction.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        ... First the "independent" production map
!-----------------------------------------------------------------------
         do j = 1,prdcnt
            do k = 2,prd_limp1
               spcno = prdmap(j,k)
               if( spcno == 0 ) then
                  exit
               else if( grpflg(spcno) /= 0 ) then
                  prdmap(j,k) = grpflg(spcno)
               else
                  prdmap(j,k) = newind(spcno)
               end if
            end do
         end do

!-----------------------------------------------------------------------
!        ... Then the "regular" reaction map
!-----------------------------------------------------------------------
         do i = 1,2
            do j = 1,rxmcnt(i)
               counter = 0
               do k = 2,i+prd_limp1
                  spcno = abs( rxmap(j,k,i) )
                  if( spcno == 0 ) then
                     exit
                  else if( grpflg(spcno) /= 0 ) then
                     rxmap(j,k,i) = SIGN( grpflg(spcno), rxmap(j,k,i) )
                     if( i == 1 ) then
                        if( k == 2 ) then
                           counter = counter + 1
                           ratind(counter) = spcno
                        end if
                     else
                        if( k <= 3 ) then
                           counter = counter + 1
                           ratind(counter) = spcno
                        end if
                     end if
                  else
                     rxmap(j,k,i) = SIGN( newind(spcno), rxmap(j,k,i) )
                  end if
               end do
               if( counter /= 0 ) then
                  rxno = rxmap(j,1,i)
                  grp_rat_cnt(counter) = grp_rat_cnt(counter) + 1
                  indx = grp_rat_cnt(counter)
                  grp_rat_map(indx,1,counter) = rxno
                  do l = 1,counter
                     grp_rat_map(indx,l+1,counter) = grp_rat_ind(ratind(l))
                  end do
                  rxt_to_grp_map(rxno,1) = counter
                  rxt_to_grp_map(rxno,2) = indx
               end if
            end do
         end do

!-----------------------------------------------------------------------
!        Scan reaction matrix to eliminate equally weigthed reactants
!        and products by setting the index = -index
!-----------------------------------------------------------------------
         do i = 1,2
            do j = 1,rxmcnt(i)
               do k = i+2,i+prd_limp1
                  spcno = rxmap(j,k,i)
                  if( spcno == 0 ) then
                     exit
                  else
                     col = pcoeff_ind(rxmap(j,1,i))
                     if( col /= 0 ) then
                        if( pcoeff(k-(i+1),col) /= 1. ) then
                           cycle
                        end if
                     end if
                     do l = 2,i+1
                        if( spcno == rxmap(j,l,i) ) then
                           rxmap(j,l,i) = -rxmap(j,l,i)
                           rxmap(j,k,i) = -spcno
                           exit
                        end if
                     end do
                  end if
               end do
            end do
         end do

!-----------------------------------------------------------------------
!        Scan reaction matrix to detect "null" reactions
!        and eliminate such reactions from the following maps:
!           1. groups
!           2. relationships
!           3. reactions
!-----------------------------------------------------------------------
         do i = 1,2
            place = 1
            do j = 1,rxmcnt(i)
               null_flag = .true.         ! assume a null reaction
               do k = 2,i+prd_limp1
                  if( rxmap(j,k,i) > 0 ) then
                     null_flag = .false.  ! not a null reaction
                     exit
                  end if
               end do
               if( null_flag ) then       ! remove from lists if null
                  rxno = rxmap(j,1,i)
                  rxt_to_grp_map(rxno,1:2) = 0
                  rel_rxt_pntr(rxno,1:2)   = 0
               else                       ! a non-null reaction; keep it
                  rxmap(place,1:i+prd_limp1,i) = rxmap(j,1:i+prd_limp1,i)
                  place = place + 1
               end if
            end do
            rxmcnt(i) = place - 1
         end do

!-----------------------------------------------------------------------
!        ... Form the solution class reaction maps
!-----------------------------------------------------------------------
         call cls_maps
     
!-----------------------------------------------------------------------
!        ... Order class reaction map reactants for the nonlinear reactions
!-----------------------------------------------------------------------
         do i = 1,5
	    if( cls_rxt_cnt(3,i) /= 0 ) then
	       do k = sum( cls_rxt_cnt(1:2,i) )+1,sum( cls_rxt_cnt(1:3,i) )
		  if( (abs(cls_rxt_map(k,2,i)) > abs(cls_rxt_map(k,3,i)) .and. &
		       cls_rxt_map(k,3,i) > 0 ) .or. cls_rxt_map(k,2,i) <= 0 ) then
		     m = cls_rxt_map(k,3,i)
		     cls_rxt_map(k,3,i) = cls_rxt_map(k,2,i)
		     cls_rxt_map(k,2,i) = m
		  end if
	       end do
	    end if
	 end do
!=======================================================================
!        ... Call the code writing utilities
!=======================================================================
!-----------------------------------------------------------------------
!        ... Force permutation for explicit method
!-----------------------------------------------------------------------
         if( clscnt(1) /= 0 ) then
	    permute(:clscnt(1),1) = (/ (i,i=1,clscnt(1)) /)
	 end if
!-----------------------------------------------------------------------
!        ... The iterated Euler backward and "Hov" methods
!-----------------------------------------------------------------------
	 do class = 2,3
            if( clscnt(class) /= 0 ) then
	       clsndx = class - 1
	       allocate( sparse(clsndx)%mat_sp_pat(clscnt(class),clscnt(class)),stat=astat )
	       if( astat /= 0 ) then
	          write(lout,*) 'Failed to allocate the matrix sparsity pattern array'
	          write(lout,*) 'for class = ',class,' ; error = ',astat
	          stop 'abort'
	       end if
	       allocate( sparse(clsndx)%lu_sp_pat(clscnt(class),clscnt(class)),stat=astat )
	       if( astat /= 0 ) then
	          write(lout,*) 'Failed to allocate the lu sparsity pattern array'
	          write(lout,*) 'for class = ',class,' ; error = ',astat
	          stop 'abort'
	       end if
	       call sparsity_pat( clscnt(class), clsmap(1,class,2), cls_rxt_cnt(1,class), &
			          cls_rxt_map(1,1,class), sparse(clsndx)%mat_sp_pat )
	       sparse(clsndx)%lu_sp_pat(:,:) = sparse(clsndx)%mat_sp_pat(:,:)
	       call diag_mark( clscnt(class), sparse(clsndx)%lu_sp_pat, permute(1,class) )
	       permute_orig(:clscnt(class),class-1) = permute(:clscnt(class),class)
	       deallocate( sparse(clsndx)%mat_sp_pat )
	       deallocate( sparse(clsndx)%lu_sp_pat )
	    end if
	 end do
!-----------------------------------------------------------------------
!        ... The sparse matrix backward Euler method
!-----------------------------------------------------------------------
sparse_matrix_loop : &
	 do class = 4,5
	    k = max( 1,clscnt(class) )
	    clsndx = class - 3
	    allocate( sparse(clsndx)%mat_sp_pat(k,k), stat=astat )
	    if( astat /= 0 ) then
	       write(lout,*) 'Failed to allocate the matrix sparsity pattern array'
	       write(lout,*) 'for class = ',class,' ; error = ',astat
	       stop 'abort'
	    end if
	    allocate( sparse(clsndx)%lu_sp_pat(k,k), stat=astat )
	    if( astat /= 0 ) then
	       write(lout,*) 'Failed to allocate the lu sparsity pattern array'
	       write(lout,*) 'for class = ',class,' ; error = ',astat
	       stop 'abort'
	    end if
	    allocate( sparse(clsndx)%mat_sp_map(k,k), stat=astat )
	    if( astat /= 0 ) then
	       write(lout,*) 'Failed to allocate the matrix sparsity map array'
	       write(lout,*) 'for class = ',class,' ; error = ',astat
	       stop 'abort'
	    end if
	    allocate( sparse(clsndx)%diag_map(k), stat=astat )
	    if( astat /= 0 ) then
	       write(lout,*) 'Failed to allocate the matrix diagonal map array'
	       write(lout,*) 'for class = ',class,' ; error = ',astat
	       stop 'abort'
	    end if
            if( clscnt(class) /= 0 ) then
!-----------------------------------------------------------------------
!        ... Determine original jacobian sparsity
!-----------------------------------------------------------------------
	       call sparsity_pat( clscnt(class), &
			          clsmap(1,class,2), &
			          cls_rxt_cnt(1,class), &
			          cls_rxt_map(1,1,class), &
			          sparse(clsndx)%mat_sp_pat )
!              call draw_mat( clscnt(class), mat_sp_pat(1,1,class-3) )
	       sparse(clsndx)%lu_sp_pat(:,:) = sparse(clsndx)%mat_sp_pat(:,:)
!-----------------------------------------------------------------------
!        ... Reorder according to diagonal Markowitz
!-----------------------------------------------------------------------
	       call diag_mark( clscnt(class), sparse(clsndx)%lu_sp_pat, permute(1,class) )
!-----------------------------------------------------------------------
!        ... Permute the original sparsity pattern
!-----------------------------------------------------------------------
	       call perm_mat( clscnt(class), sparse(clsndx)%lu_sp_pat, permute(1,class) )
	       sparse(clsndx)%mat_sp_pat(:,:) = sparse(clsndx)%lu_sp_pat(:,:)
!              call draw_mat( clscnt(class), lu_sp_pat )
!-----------------------------------------------------------------------
!        ... Symbolic factorization; includes fillin
!-----------------------------------------------------------------------
	       call sym_fac( clscnt(class), sparse(clsndx)%lu_sp_pat, additions, multiplications )
!-----------------------------------------------------------------------
!        ... Make column oriented non-zero "map"
!-----------------------------------------------------------------------
               nzcnt(class-3) = COUNT( sparse(clsndx)%lu_sp_pat(:,:) )
	       sparse(clsndx)%mat_sp_map(:,:) = 0
	       k = 0
	       do j = 1,clscnt(class)
	          do i = 1,clscnt(class)
                     if( sparse(clsndx)%lu_sp_pat(i,j) ) then
		        k = k + 1
		        sparse(clsndx)%mat_sp_map(i,j) = k
		        if( i == j ) then
		           sparse(clsndx)%diag_map(j) = k
		        end if
		     end if
	          end do
	       end do
            end if
!-----------------------------------------------------------------------
!        ... Write the factorization code
!-----------------------------------------------------------------------
            if( class /= 5 .or. model == 'MOZART' ) then
               call make_lu_fac( clscnt(class), class, sparse(clsndx)%lu_sp_pat, &
	                         sparse(clsndx)%mat_sp_pat, sparse(clsndx)%mat_sp_map, march, model )
!-----------------------------------------------------------------------
!        ... Write the solver code
!-----------------------------------------------------------------------
               call make_lu_slv( clscnt(class), class, sparse(clsndx)%lu_sp_pat, march, model )
	    end if
	    if( associated( sparse(clsndx)%lu_sp_pat ) ) then
	       deallocate( sparse(clsndx)%lu_sp_pat )
	    end if
         end do sparse_matrix_loop
!-----------------------------------------------------------------------
!        ... Make "from-to" permutation
!-----------------------------------------------------------------------
	 do class = 2,5
	    do j = 1,clscnt(class)
	       do i = 1,clscnt(class)
	          if( permute(i,class) == j ) then
	             permutation(j) = i
		  end if
	       end do
	    end do
	    if( clscnt(class) /= 0 ) then
	       permute(:clscnt(class),class) = permutation(:clscnt(class))
	    end if
	 end do
!-----------------------------------------------------------------------
!        ... Make reaction scheme dependent prod & loss code
!-----------------------------------------------------------------------
         call pl_code( new_nq, clscnt, clsmap, cls_rxt_cnt, cls_rxt_map, &
                       pcoeff_ind, pcoeff, permute, march, model )
         cls_ind_prdcnt = sum( cls_rxt_cnt(1,1:5) )
!-----------------------------------------------------------------------
!        ... Make reaction scheme independent prod & loss code
!-----------------------------------------------------------------------
         call ipd_code( new_nq, clscnt, clsmap, cls_rxt_cnt, extcnt, &
                        cls_rxt_map, pcoeff_ind, pcoeff, permute, model, march )
!-----------------------------------------------------------------------
!        ... Make tabular reaction rates
!-----------------------------------------------------------------------
         if( options(12) ) then
            call make_rate_tab( rxparm, rxptab, rxpcnt )
         else
            call make_rate( sym_rates, rxptab, rxpcnt, machine, vec_ftns, &
                            model, march )
         end if
         if( rxmcnt(2) /= 0 .or. fixcnt(2) /= 0 ) then
            radj_flag = .true.
         else
            radj_flag = .false.
            do i = 1,fixcnt(1)
               if( abs(fixmap(i,1,1)) > phtcnt ) then
                  radj_flag = .true.
               end if
            end do
         end if
!-----------------------------------------------------------------------
!        ... Make reaction rate adjustment routines
!-----------------------------------------------------------------------
         call make_radj( fixmap, fixcnt, rxmap(1,1,2), rxmcnt(2), phtcnt, &
                         model, march )
         call make_padj( fixmap, fixcnt(1), phtcnt, model, march )
         call make_rmod( rel_rxt_pntr, rel_rxt_map, rxt_to_grp_map, grp_rat_map, hetmap(1,2), &
                         hetcnt, rxntot, model, march )
         do class = 4,5
            clsndx = class - 3
            allocate( lin_mat_pat(nzcnt(clsndx)), stat=astat )
            if( astat /= 0 ) then
               write(lout,*) 'Failed to allocate the lin matrix map array'
               write(lout,*) 'for class = ',class,' ; error = ',astat
               stop 'abort'
            end if
            if( class /= 5 .or. model == 'MOZART' ) then
               call make_lin( clscnt(class), clsmap, cls_rxt_cnt(1,class), cls_rxt_map(1,1,class), pcoeff_ind, &
                              pcoeff, permute(1,class), sparse(clsndx)%mat_sp_map, class,  &
                              lin_mat_pat, march, model )
               call make_nln( clscnt(class), clsmap, cls_rxt_cnt(1,class), cls_rxt_map(1,1,class), pcoeff_ind, &
                              pcoeff, permute(1,class), sparse(clsndx)%mat_sp_map, class, &
                              lin_mat_pat, nzcnt(clsndx), sparse(clsndx)%diag_map, march, model )
	    end if
	    if( associated( sparse(clsndx)%mat_sp_pat ) ) then
	       deallocate( sparse(clsndx)%mat_sp_pat )
	    end if
	    if( allocated( lin_mat_pat ) ) then
	       deallocate( lin_mat_pat )
	    end if
         end do
!-----------------------------------------------------------------------
!        ... Make group members vmr subroutine
!-----------------------------------------------------------------------
         call mak_grp_vmr( grp_mem_cnt, mem2grp_map, model, march )
!-----------------------------------------------------------------------
!        ... Writeout the surface flux and depos vel info
!-----------------------------------------------------------------------
         if( srf_flx_cnt /= 0 ) then
            write(lout,*) ' '
            write(lout,'('' Species with non-zero surface emission'')')
            do i = 1,srf_flx_cnt
               write(lout,'(1x,''('',i2,'')'',3x,a8)') i, new_solsym(srf_flx_map(i))
            end do
         end if
         if( dvel_cnt /= 0 ) then
            write(lout,*) ' '
            write(lout,'('' Species with non-zero dry deposition flux'')')
            do i = 1,dvel_cnt
               write(lout,'(1x,''('',i2,'')'',3x,a8)') i, new_solsym(dvel_map(i))
            end do
         end if
!-----------------------------------------------------------------------
!        ... Call the equation reporting utility
!-----------------------------------------------------------------------
         call equation_rep( new_nq, new_solsym, nfs, fixsym, prdcnt, &
                            prdmap, rxntot, rxmcnt, rxmap, pcoeff_cnt, &
                            pcoeff_ind, pcoeff, fixcnt, fixmap, phtcnt )

         call write_rxt_out_code ( &
              rxmcnt, &
              rxmap, &
              fixmap, &
              solsym, &
              fixsym, &
              prdcnt, &
              prdmap, &
              rxntot, &
              phtcnt, & 
              rxt_rates_conv_file )

      end if Has_chemistry

!=======================================================================
!     ... This is for the new CTM interface; the old driver
!         data file can still be output for potential diagnostics
!=======================================================================
      open( unit = 20, &
            file   = trim( sim_dat_filespec ), &
            status = 'replace', &
            iostat = ios )
      if( ios /= 0 ) then
	 write(*,*) ' Failed to open file '// trim( sim_dat_filespec )
	 write(*,*) ' Error code = ',ios
	 stop
      end if
      if( model == 'MOZART' ) then
         do class = 1,5
	    if( class == 2 .or. class == 3 ) then
	       cycle
	    end if
            if( clscnt(class) > 0 ) then
                  write(20,508) cls_rxt_cnt(:,class)
                  write(20,522) clsmap(:clscnt(class),class,2)
               if( class >= 4 ) then
                  write(20,522) permute(:clscnt(class),class)
                  write(20,522) sparse(class-3)%diag_map(:clscnt(class))
               end if
            end if
         end do
!-----------------------------------------------------------------------
!        ... Write the "class" maps & species symbols
!-----------------------------------------------------------------------
         temp_mass(:) = 0.
         if( nq > 0 ) then
            do i = 1,nq
	       if( newind(i) /= 0 ) then
	          temp_mass(newind(i)) = mass(i)
	       end if
            end do
            write(20,*) temp_mass(:new_nq)
         end if
         if( ngrp > 0 ) then
	    do i = 1,nq
	       if( grp_rat_ind(i) /= 0 ) then
	          temp_mass(grp_rat_ind(i)) = mass(i)
	       end if
	    end do
            write(20,*) temp_mass(:grp_mem_cnt)
         end if
         if( new_nq > 0 ) then
            write(20,'(10a16)') new_solsym(:new_nq)
         end if
         if( grp_mem_cnt > 0 ) then
            write(20,'(i4)') ngrp
            write(20,'(20i4)') grpcnt(:ngrp)
            write(20,'(10a16)') grpsym(:ngrp)
            write(20,'(10a16)') grp_mem_sym(:grp_mem_cnt)
         end if
         write(20,'(i4)') srf_flx_cnt
         if( srf_flx_cnt > 0 ) then
	    do m = 1,(srf_flx_cnt-1)/10+1
	       wrk_chr(:) = ' '
	       il = (m-1)*10 + 1
	       iu = min( 10*m,srf_flx_cnt )
	       do i = il,iu
	          wrk_chr(i-il+1) = new_solsym(srf_flx_map(i))
	       end do
               write(20,'(10a16)') wrk_chr
	    end do
         end if
         write(20,'(i4)') dvel_cnt
         if( dvel_cnt > 0 ) then
	    do m = 1,(dvel_cnt-1)/10+1
	       wrk_chr(:) = ' '
	       il = (m-1)*10 + 1
	       iu = min( 10*m,dvel_cnt )
	       do i = il,iu
	          wrk_chr(i-il+1) = new_solsym(dvel_map(i))
	       end do
               write(20,'(10a16)') wrk_chr
	    end do
         end if
!-----------------------------------------------------------------------
!        ... Write the wet removal species
!-----------------------------------------------------------------------
         if( hetcnt > 0 ) then
	    do m = 1,(hetcnt-1)/10+1
	       wrk_chr(:) = ' '
	       il = (m-1)*10 + 1
	       iu = min( 10*m,hetcnt )
	       do i = il,iu
	          wrk_chr(i-il+1) = new_solsym(hetmap(i,1))
	       end do
               write(20,'(10a16)') wrk_chr
	    end do
         end if
!-----------------------------------------------------------------------
!        ... Write the ext frcing species
!-----------------------------------------------------------------------
         if( usrcnt > 0 ) then
	    do m = 1,(usrcnt-1)/10+1
	       wrk_chr(:) = ' '
	       il = (m-1)*10 + 1
	       iu = min( 10*m,usrcnt )
	       do i = il,iu
	          wrk_chr(i-il+1) = new_solsym(usrmap(i))
	       end do
               write(20,'(10a16)') wrk_chr
	    end do
         end if
!-----------------------------------------------------------------------
!        ... Write the rxt aliases
!-----------------------------------------------------------------------
         k = count( rxt_has_alias(:rxntot) )
         write(20,'(i4)') k
         if( k > 0 ) then
	    wrk_chr(:) = ' '
	    i = 0
	    do m = 1,rxntot
	       if( rxt_has_alias(m) ) then
	          i = i + 1
	          wrk_rxt(i) = rxt_tag(m)
	          if( i == 5 ) then
                     write(20,'(5a16)') wrk_rxt(:5)
		     i = 0
	          end if
	       end if
	    end do
	    if( i /= 0 ) then
               write(20,'(5a16)') wrk_rxt(:i)
	    end if
	    i = 0
	    do m = 1,rxntot
	       if( rxt_has_alias(m) ) then
	          i = i + 1
	          nind(i) = m
	          if( i == 20 ) then
                     write(20,'(20i4)') nind(:20)
		     i = 0
	          end if
	       end if
	    end do
	    if( i /= 0 ) then
               write(20,'(20i4)') nind(:i)
	    end if
         end if
      else
         write(20,510) clscnt(:)
         write(20,508) cls_rxt_cnt(:,:)
         do class = 1,5
	    if( class == 2 .or. class == 3 ) then
	       cycle
	    end if
            if( clscnt(class) > 0 ) then
               write(20,522) clsmap(:clscnt(class),class,2)
            end if
         end do
         temp_mass(:) = 0.
         if( nq > 0 ) then
            do i = 1,nq
	       if( newind(i) /= 0 ) then
	          temp_mass(newind(i)) = mass(i)
	       end if
            end do
            write(20,*) temp_mass(:new_nq)
         end if
         if( new_nq > 0 ) then
            write(20,'(10a16)') new_solsym(:new_nq)
         end if
         do class = 2,5
            if( clscnt(class) > 0 ) then
               write(20,522) permute(:clscnt(class),class)
               if( class > 3 ) then
                  write(20,522) sparse(class-3)%diag_map(:clscnt(class))
               end if
            end if
         end do
         call make_sim_dat( model, sparse )
      end if
      close( 20)

     rxt_tag_cnt = count( rxt_has_tag )

!-----------------------------------------------------------------------
!        ... Write the chemistry header file
!-----------------------------------------------------------------------
      call chm_hdr( rxt_tag_cnt, hetcnt, usrcnt, cls_rxt_cnt, radj_flag, phtcnt, &
                    rxpcnt, rxparm, rxntot, ncol, nfs, nslvd, &
                    indexm, indexh2o, new_nq, relcnt, grp_mem_cnt, &
                    clscnt, iter_counts, nzcnt, vec_ftns, machine, options(1) )

!-----------------------------------------------------------------------
!        ... Write the resolution header file
!-----------------------------------------------------------------------
!      call res_hdr( plon, plonl, plat, plev, jintmx, &
!                    nxpt, arch_type, cpucnt )

!-----------------------------------------------------------------------
!        ... Write the version header file
!-----------------------------------------------------------------------
      ptplen = histout_cnt(1,1,1) + histout_cnt(2,1,1) + histout_cnt(5,1,1) &
             + histout_cnt(3,1,1) + histout_cnt(4,1,1) + histout_cnt(7,1,1) &
             + histout_cnt(1,2,1) + histout_cnt(2,2,1) + histout_cnt(5,2,1) &
             + histout_cnt(3,2,1) + histout_cnt(4,2,1) + histout_cnt(7,2,1)
      if( ptplen /= 0 .and. histout(2) /= ' ' ) then
         ohstflag = .true.
      else
         ohstflag = .false.
      end if
      call ver_hdr( options, plon, plonl, plev, machine, &
                    model, arch_type, ohstflag, diagprnt, tavgprnt, srf_flx_cnt, &
		    hetcnt, rxntot, clscnt, nzcnt, new_nq, dvel_cnt )

!-----------------------------------------------------------------------
!        ... Write the slt header file
!-----------------------------------------------------------------------
      call slt_hdr( options(2), options(11), cpucnt, machine )

!-----------------------------------------------------------------------
!        ... Write the history tape header file
!-----------------------------------------------------------------------
      call hist_hdr( hst_file_cnt, histout_cnt, histout_map, user_hst_names, histinp(4), &
                     dyn_hst_fld_cnt, spcsym, spccnt, hetmap, usrmap, &
		     ptplen, sim_dat_filespec, model )
      if( model == 'MOZART' ) then
!-----------------------------------------------------------------------
!       ... Special section for invariants
!-----------------------------------------------------------------------
         open( unit = 20, &
               file     = trim( sim_dat_filespec ), &
               status   = 'old', &
               position = 'append', &
               iostat   = ios )
         if( ios /= 0 ) then
	    write(*,*) ' Failed to open file '// trim( sim_dat_filespec )
	    write(*,*) ' Error code = ',ios
	    stop
         end if
!-----------------------------------------------------------------------
!        ... Write the invariants
!-----------------------------------------------------------------------
         if( nfs /= 0 ) then
	    do m = 1,(nfs-1)/10+1
	       wrk_chr(:) = ' '
	       il = (m-1)*10 + 1
	       iu = min( 10*m,nfs )
	       do i = il,iu
	          wrk_chr(i-il+1) = fixsym(i)
	       end do
               write(20,'(10a16)') wrk_chr
	    end do
         end if
         close( 20 )
      end if

!-----------------------------------------------------------------------
!       ... Write the files.h file
!-----------------------------------------------------------------------
      call files_hdr

!-----------------------------------------------------------------------
!       ... Form "include" files stub file
!-----------------------------------------------------------------------
      inquire( file = trim( temp_path ) // 'wrk.stub.F', exist = lexist )
      if( lexist ) then
	 call system( 'rm -f ' // trim( temp_path ) // 'wrk.stub.F' )
      end if
      open( unit = 3, &
            file = trim( temp_path ) // 'wrk.stub.F', &
            iostat = ios )
      if( ios /= 0 ) then
	 write(*,*) ' Failed to open wrk file; terminating'
	 write(*,*) ' Error code = ',ios
	 stop
      end if
      write(3,'(''# include <version.h>'')') 
!!$      write(3,'(''# include <res.h>'')') 
      write(3,'(''# include <slt.h>'')') 
      write(3,'(''# include <chem.h>'')')
      write(3,'(''# include <hist.h>'')') 
      close( 3 )
      call system( 'rm -f wrk.F' )

!-----------------------------------------------------------------------
!       ... Check for cpp utility
!-----------------------------------------------------------------------
      cpp_command = trim( cpp_dir ) // ' ' // trim( cpp_opts )
      inquire( file = trim( cpp_dir ), exist = lexist )
      if( .not. lexist ) then
	 buff = ' '
         call errmes( ' ** cpp not in #', &
                      lout, &
                      trim( cpp_dir ), &
                      len_trim( cpp_dir ), &
                      buff )
      end if
!-----------------------------------------------------------------------
!       ... Check for fortran 90 modules
!-----------------------------------------------------------------------
      inquire( file = trim( temp_path ) // 'mod.src.files.PP', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim( temp_path ) // 'mod.src.files.PP' )
      end if
!-----------------------------------------------------------------------
!       ... Get mozart module files
!-----------------------------------------------------------------------
      if( model == 'MOZART' ) then
         inquire( file = trim( src_dir ) // 'mozart.mod.files.PP', exist = lexist )
         if( .not. lexist ) then
            call errmes( ' ** Module source file missing@', 6, param, 8, buff )
         end if
         call system( 'cat ' // trim( src_dir ) // 'mozart.mod.files.PP  > ' &
	                     // trim( temp_path ) // 'mod.src.files.PP' )
      else
         inquire( file = trim( src_dir ) // 'cam.mod.files.PP', exist = lexist )
         if( .not. lexist ) then
            call errmes( ' ** Module source file missing@', 6, param, 8, buff )
         end if
         call system( 'cat ' // trim( src_dir ) // 'cam.mod.files.PP  > ' &
	                     // trim( temp_path ) // 'mod.src.files.PP' )
      end if
      inquire( file = trim( temp_path ) // 'mod.src.files', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim( temp_path ) // 'mod.src.files' )
      end if
      call system( trim( cpp_command ) // ' ' // trim( temp_path ) // 'mod.src.files.PP > ' &
				       // trim( temp_path ) // 'mod.src.files' )
      close(2)
      open( unit = 2, &
            file = trim( temp_path ) // 'mod.src.files', &
            status = 'old', &
            position = 'rewind', &
	    iostat = ios )
      if( ios /= 0 ) then
         write(lout,*) ' Failed to open mod.src.files file'
	 write(lout,*) ' Error code = ',ios
	 stop
      end if
      file_cnt = 1
      do k = 1,500
         read(2,'(a320)',end=1105) mod_src(file_cnt)
         if( mod_src(file_cnt) /= ' ' ) then
	    filelines(5) = filelines(5) + 1
	    file_cnt = file_cnt + 1
	 end if
      end do
1105  continue
      close( 2 )
!-----------------------------------------------------------------------
!       ... Check for species names module
!-----------------------------------------------------------------------
      if( options(18) ) then
	 k = 1
!!$	 mod_paths(k) = './'
	 mod_paths(k) = ' '
	 mod_names(k) = trim( temp_path ) // 'spc_names.mod'
         k = 2
!!$         mod_paths(k) = './'
	 mod_paths(k) = ' '
         mod_names(k) = trim( temp_path ) // 'rxt_names.mod'
!!$         k = 3
!!$         mod_paths(k) = './'
!!$	 mod_paths(k) = ' '
!!$         mod_names(k) = trim( temp_path ) // 'het_names.mod'
      else
	 k = 0
      end if
!-----------------------------------------------------------------------
!       ... Check user files for any module files
!-----------------------------------------------------------------------
      do i = 1,sub_cnt
	 if( index( filename(i), 'mod', back = .true. ) == len_trim(filename(i))-2 ) then
	    k = k + 1
	    mod_paths(k) = trim( filepath(i) )
	    mod_names(k) = trim( filename(i) )
	    nend(i) = 0
	 else
	    nend(i) = 1
	 end if
      end do
      filelines(2) = k
!-----------------------------------------------------------------------
!       ... Remove .mod files from user subroutine lists
!-----------------------------------------------------------------------
      k = 0
      do i = 1,sub_cnt
	 if( nend(i) == 1 ) then
	    k = k + 1
	    filepath(k) = filepath(i)
	    filename(k) = filename(i)
	 end if
      end do
      sub_cnt = k
!-----------------------------------------------------------------------
!       ... Check for user file "overrides"
!-----------------------------------------------------------------------
      if( filelines(2) /= 0 ) then
         call sub_scan( filelines(5), &
                        mod_src, &
                        mod_paths, &
                        mod_names, &
                        filelines(2) )
      end if
!-----------------------------------------------------------------------
!       ... Form and preprocess the modules
!-----------------------------------------------------------------------
      call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
      do i = 1,filelines(5)
	 command = 'cat ' // trim( mod_src(i) ) // ' >> wrk.F'
         call system( trim( command ) )
      end do
      if( filelines(2) > 0 ) then
         do i = 1,filelines(2)
            command = 'cat ' // trim( mod_paths(i) ) // trim( mod_names(i) ) // ' >> wrk.F'
            call system( trim( command ) )
         end do
      end if
      inquire( file = trim(procout_path)//'moz.mods.F90', exist = lexist )
      if( lexist ) then
         call system( 'rm -f '//trim(procout_path)//'moz.mods.F90' )
      end if
      call system( trim( cpp_command ) // ' wrk.F > '//trim(procout_path)//'moz.mods.F90' )
      call system( 'rm -f wrk.F' )

!-----------------------------------------------------------------------
!       ... for CAM,WRF remove tar file if it exists
!-----------------------------------------------------------------------
      if( model /= 'MOZART' ) then
         if( model == 'CAM' ) then
            tar_flnm = 'cam.subs.tar'
         else if( model == 'WRF' ) then
            tar_flnm = 'wrf.subs.tar'
         end if
         inquire( file = trim( temp_path ) // trim(tar_flnm), exist = lexist )
         if( lexist ) then
            call system( 'rm -f ' // trim(temp_path ) // trim(tar_flnm) )
         end if
!-----------------------------------------------------------------------
!       ... add module files to cam tar file
!-----------------------------------------------------------------------
         if( filelines(5) /= 0 ) then
            do i = 1,filelines(5)
               call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
               call system( 'cat '// trim( mod_src(i) ) // ' >> wrk.F' )
!              write(*,*) 'cpp file ',trim(mod_src(i))
	       il = index( mod_src(i), '/', back = .true. ) + 1
	       iu = index( mod_src(i), '.mod', back = .true. ) - 1
               select case( mod_src(i)(il:iu) )
                  case( 'mo_chem' )
                     tmp_filenm = 'chem_mods'
                  case default
                     tmp_filenm = mod_src(i)(il:iu)
               end select 
               filenm = trim(tmp_filenm) // '.F90'
!              write(*,*) 'tar file ',trim(filenm)
               call system( trim( cpp_command ) // ' wrk.F > '// trim(filenm) )
               if( i == 1 ) then
                  call system( 'tar -cf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               else
                  call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               end if
               call system( 'rm -f wrk.F' )
               call system( 'rm -f ' // trim(filenm) )
            end do
         end if
         if( filelines(2) > 0 ) then
            do i = 1,filelines(2)
               call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
               call system( 'cat '// trim( mod_paths(i) ) // trim(mod_names(i) ) // ' >> wrk.F' )
!              write(*,*) 'cpp file ',trim(mod_paths(i)) // trim(mod_names(i))
               il = index( mod_names(i), '/', back = .true. ) + 1
               iu = index( mod_names(i), '.mod', back = .true. ) - 1
               select case( mod_names(i)(il:iu) )
               case( 'spc_names' )
                  tmp_filenm = 'm_spc_id'
               case( 'rxt_names' )
                  tmp_filenm = 'm_rxt_id'
               case( 'het_names' )
                  tmp_filenm = 'm_het_id'
               case default
                  tmp_filenm = mod_src(i)(il:iu)
               end select 
               filenm = trim(tmp_filenm) // '.F90'
!              write(*,*) 'tar file ',trim(filenm)
               call system( trim( cpp_command ) // ' wrk.F > '// trim(filenm) )
               if( filelines(5) == 0 .and. i == 1 ) then
                  call system( 'tar -cf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               else
                  call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               end if
               call system( 'rm -f wrk.F' )
               call system( 'rm -f ' // trim(filenm) )
            end do
         end if
      end if

!-----------------------------------------------------------------------
!       ... Get all source files
!-----------------------------------------------------------------------
      if( options(7) ) then
!-----------------------------------------------------------------------
!       ... Get "main" library files
!-----------------------------------------------------------------------
         inquire( file = trim( temp_path ) // 'lib.src.files.PP', exist = lexist )
         if( lexist ) then
            call system( 'rm -f ' // trim(temp_path ) // 'lib.src.files.PP' )
         end if
!-----------------------------------------------------------------------
!       ... Get "chemistry" library files
!-----------------------------------------------------------------------
         if( model == 'MOZART' ) then
            inquire( file = trim( src_dir ) // 'mozart.src.files.PP', exist = lexist )
         else if( model == 'CAM' ) then
            inquire( file = trim( src_dir ) // 'cam.src.files.PP', exist = lexist )
         end if
         if( .not. lexist ) then
            call errmes( ' ** Chem source file missing@', 6, param, 8, buff )
         end if
         if( model == 'MOZART' ) then
            call system( 'cat ' // trim( src_dir ) // 'mozart.src.files.PP ' // ' > ' // trim( temp_path ) // 'lib.src.files.PP' )
         else if( model == 'CAM' ) then
            call system( 'cat ' // trim( src_dir ) // 'cam.src.files.PP ' // ' > ' // trim( temp_path ) // 'lib.src.files.PP' )
         end if
         command = trim( cpp_command ) // ' ' // trim( temp_path ) // 'lib.src.files.PP > ' // trim( temp_path ) // 'lib.src.files'
         inquire( file = trim( temp_path ) // 'lib.src.files', exist = lexist )
         if( lexist ) then
            call system( 'rm -f ' // trim(temp_path ) // 'lib.src.files' )
         end if
         call system( trim( command ) )
         close(2)
         open( unit = 2, &
               file = trim( temp_path ) // 'lib.src.files', &
               status = 'old', &
               iostat = ios )
         if( ios /= 0 ) then
            write(lout,*) ' Failed to open lib.src.files file'
            write(lout,*) ' Error code = ',ios
            stop
         end if
         file_cnt = 1
         do k = 1,500
            read(2,'(a320)',end=1005) lib_src(file_cnt)
            if( lib_src(file_cnt) /= ' ' ) then
               filelines(1) = filelines(1) + 1
               file_cnt = file_cnt + 1
            end if
         end do
1005     continue
!-----------------------------------------------------------------------
!       ... set cam implicit solvers
!-----------------------------------------------------------------------
         if( model == 'CAM' ) then
            if( march == 'SCALAR' ) then
               lib_src(file_cnt)   = trim(procfiles_path) // 'mo_imp_sol_scalar.F90'
            else if( march == 'CACHE' ) then
               lib_src(file_cnt)   = trim(procfiles_path) // 'mo_imp_sol_cache.F90'
            else if( march == 'VECTOR' ) then
               lib_src(file_cnt)   = trim(procfiles_path) // 'mo_imp_sol_vector.F90'
            end if
	    file_cnt = file_cnt + 1
            filelines(1) = filelines(1) + 1
         end if
!-----------------------------------------------------------------------
!       ... check for iterative convergence norms
!-----------------------------------------------------------------------
         if( options(5) ) then
            lib_src(file_cnt)   = 'del_norm.F'
            lib_src(file_cnt+1) = 'it_norm.F'
            filelines(1) = filelines(1) + 2
         end if
!-----------------------------------------------------------------------
!       ... Check for user file "overrides"
!-----------------------------------------------------------------------
         if( sub_cnt /= 0 ) then
            call sub_scan( filelines(1), lib_src, filepath, filename, sub_cnt )
         end if
!-----------------------------------------------------------------------
!       ... Form main lib portion of ccmpp file
!-----------------------------------------------------------------------
	 call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
         if( filelines(1) /= 0 ) then
            do i = 1,filelines(1)
	       command = 'cat '// trim( lib_src(i) ) // ' >> wrk.F'
	       call system( trim( command ) )
            end do
         end if

!-----------------------------------------------------------------------
!       ... Get user specified files
!-----------------------------------------------------------------------
         if( sub_cnt > 0 ) then
            do i = 1,sub_cnt
	       command = 'cat ' // trim( filepath(i) ) // trim( filename(i) ) // ' >> wrk.F'
	       call system( trim( command ) )
            end do
         end if
      end if
      inquire( file = trim( procout_path ) // 'moz.subs.F90', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim(procout_path ) // 'moz.subs.F90' )
      end if
      call system( trim( cpp_command ) // ' wrk.F > '// trim(procout_path) // 'moz.subs.F90' )
      call system( 'rm -f wrk.F' )

!-----------------------------------------------------------------------
!       ... add source files to cam tar file
!-----------------------------------------------------------------------
      if( model /= 'MOZART' ) then
         if( filelines(1) /= 0 ) then
            do i = 1,filelines(1)
	       il = index( lib_src(i), '/', back = .true. ) + 1
	       iu = index( lib_src(i), '.F', back = .true. ) - 1
               if( lib_src(i)(il:iu) /= 'mo_setrxt' .and. lib_src(i)(il:iu) /= 'mo_sim_dat' ) then
                  call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
!                 write(*,*) 'cpp file ',trim(lib_src(i))
               end if
               call system( 'cat '// trim( lib_src(i) ) // ' >> wrk.F' )
               if( lib_src(i)(il:iu) == 'mo_imp_sol_scalar' .or. &
                   lib_src(i)(il:iu) == 'mo_imp_sol_cache' .or. &
                   lib_src(i)(il:iu) == 'mo_imp_sol_vector' ) then
                   filenm = 'mo_imp_sol.F90'
               else
                  filenm = lib_src(i)(il:iu) // '.F90'
               end if
               if( lib_src(i)(il:iu) /= 'mo_setrxt' .and. lib_src(i)(il:iu) /= 'mo_sim_dat' ) then
                  call system( trim( cpp_command ) // ' wrk.F > '// trim(filenm) )
               else
                  call system( 'cp wrk.F ' // trim(filenm) )
               end if
!              write(*,*) 'tar file ',trim(filenm)
               call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               call system( 'rm -f wrk.F' )
               call system( 'rm -f ' // trim(filenm) )
            end do
         end if
         if( sub_cnt > 0 ) then
            do i = 1,sub_cnt
               call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
               call system( 'cat '// trim( filepath(i) ) // trim(filename(i) ) // ' >> wrk.F' )
!              write(*,*) 'cpp file ',trim(filepath(i)) // trim(filename(i))
	       il = index( filename(i), '/', back = .true. ) + 1
	       iu = index( filename(i), '.F', back = .true. ) - 1
               filenm = filename(i)(il:iu) // '.F90'
               write(*,*) 'tar file ',trim(filenm)
               call system( trim( cpp_command ) // ' wrk.F > '// trim(filenm) )
               call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               call system( 'rm -f wrk.F' )
               call system( 'rm -f ' // trim(filenm) )
            end do
         end if
      end if

!-----------------------------------------------------------------------
!       ... Get all matrix and production/loss files
!-----------------------------------------------------------------------
      inquire( file = trim( temp_path ) // 'lib.mat.files.PP', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim( temp_path ) // 'lib.mat.files.PP' )
      end if
      inquire( file = trim( src_dir ) // 'mozart.mat.files.PP', exist = lexist )
      if( .not. lexist ) then
         call errmes( ' ** Matrix source file missing@', 6, param, 8, buff )
      end if
      call system( 'cat ' // trim( src_dir ) // 'mozart.mat.files.PP ' // ' > ' &
		          // trim( temp_path) // 'lib.mat.files.PP' )
      inquire( file = trim( temp_path ) // 'lib.mat.files', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim( temp_path ) // 'lib.mat.files' )
      end if
      call system( trim( cpp_command ) // ' ' // trim( temp_path ) // 'lib.mat.files.PP > ' &
				       // trim( temp_path ) // 'lib.mat.files' )
      close(2)
      open( unit = 2, &
            file = trim( temp_path ) // 'lib.mat.files', &
            status = 'old', &
            iostat = ios )
      if( ios /= 0 ) then
	 write(lout,*) ' Failed to open ',trim(temp_path) // 'lib.mat.files file'
	 write(lout,*) ' Error code = ',ios
	 stop
      end if
      file_cnt = 1
      do k = 1,500
         read(2,'(a320)',end=1015) lib_src(file_cnt)
	 if( lib_src(file_cnt) /= ' ' ) then
	    lib_src(file_cnt) = trim( temp_path ) // trim( lib_src(file_cnt) )
	    filelines(3) = filelines(3) + 1
	    file_cnt = file_cnt + 1
	 end if
      end do
1015  continue
!-----------------------------------------------------------------------
!       ... Form main lib portion of ccmpp file
!-----------------------------------------------------------------------
      call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
      if( filelines(3) /= 0 ) then
         do i = 1,filelines(3)
	    command = 'cat '// trim( lib_src(i) ) // ' >> wrk.F'
	    call system( trim( command ) )
         end do
      end if

      inquire( file = trim( procout_path ) // 'moz.mat.F90', exist = lexist )
      if( lexist ) then
         call system( 'rm -f ' // trim(procout_path ) // 'moz.mat.F90' )
      end if
      call system( trim( cpp_command ) // ' wrk.F > '//trim(procout_path)//'moz.mat.F90' )
      call system( 'rm -f wrk.F' )

!-----------------------------------------------------------------------
!       ... add matrix source files to cam tar file
!-----------------------------------------------------------------------
      if( model /= 'MOZART' ) then
         if( filelines(3) /= 0 ) then
            do i = 1,filelines(3)
               call system( 'cat ' // trim( temp_path ) // 'wrk.stub.F > wrk.F' )
               call system( 'cat '// trim( lib_src(i) ) // ' >> wrk.F' )
!              write(*,*) 'cpp file ',trim(lib_src(i))
               il = index( lib_src(i), '/', back = .true. ) + 1
               iu = index( lib_src(i), '.F', back = .true. ) - 1
               select case( lib_src(i)(il:iu) )
                  case( 'prd_loss' )
                     mod_src(1) = 'prod_loss'
                  case( 'lu_fac' )
                     mod_src(1) = 'lu_factor'
                  case( 'lu_slv' )
                     mod_src(1) = 'lu_solve'
                  case( 'linmat' )
                     mod_src(1) = 'lin_matrix'
                  case( 'nlnmat' )
                     mod_src(1) = 'nln_matrix'
                  case default
                     mod_src(1) = lib_src(i)(il:iu)
               end select 
               filenm = 'mo_' // trim(mod_src(1)) // '.F90'
!              write(*,*) 'tar file ',trim(filenm)
               call system( trim( cpp_command ) // ' wrk.F > '// trim(filenm) )
               call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(filenm) )
               call system( 'rm -f wrk.F' )
               call system( 'rm -f ' // trim(filenm) )
            end do

            call system( 'cp '// trim(temp_path) // trim(rxt_rates_conv_file) // ' .' )
            call system( 'tar -rf ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(rxt_rates_conv_file) )
            call system( 'rm -f ' // trim(rxt_rates_conv_file) )

            call system( 'mv ' // trim(temp_path) // trim(tar_flnm) // ' ' // trim(output_path) // '.' )
         end if
      end if


!-----------------------------------------------------------------------
!       ... Write the params.h file
!-----------------------------------------------------------------------
      call params_hdr( plon, plonl, plat, plev, phtcnt, &
                       rxntot, new_nq, grp_mem_cnt, histout_cnt, options(1), &
                       options(3), options(4), arch_type, trim(procout_path)//'params.h' )
!-----------------------------------------------------------------------
!       ... Clean up this directory
!-----------------------------------------------------------------------
      call system( 'mv *.h ' // trim( temp_path ) // '.' )

      write(*,*) ' '
      write(*,*) '================================================'
      write(*,*) 'CAM-Chem preprocessor has successfully completed'
      write(*,*) '================================================'
      write(*,*) ' '

!-----------------------------------------------------------------------
!       ... Format statements
!-----------------------------------------------------------------------
101   format('0 *****Species header must be first card *****')
102   format('0 *****Solution must follow species card *****')

202   format(6x,i3,2x,i3)
204   format(6x,i3,2x,i3,2x,i3)
206   format(6x,i3,2x,i3,2x,i3,2x,i3)
208   format(6x,i3,2x,i3,2x,i3,2x,i3,2x,i3)
201   format('0     the unimolecular fixed map'/6x,'rxn',2x,'fsn')
203   format('0     the bimolecular fixed map'/6x,'rxn',2(2x,'fsn'))
205   format('0     the production map'/6x,'rxn',2(2x,'psn'))
209   format('0     the unimolecular reaction map'/6x,'rxn',2x,'rsn',2(2x,'psn'))
210   format(6x,i3,2x,i3,2x,i3,2x,i3,2x,i3,2x,i3)
211   format('0     the bimolecualr reaction map'/6x,'rxn',2(2x,'rsn'),2(2x,'psn'))
213   format('0     the pce loss map'/6x,'pcn',2x,'rxn',2(2x,'psn'))
215   format('0     the pce,sol map'/6x,'pcn',2x,'rxn',2x,'rsn',2(2x,'psn'))
217   format('0     pure prod map for pces'/6x,'pcn',2x,'rxn',2x,'ind')
219   format('0     the linear prod map for pces'/6x,'pcn',2x,'rxn',2x,'rsn',2x,'ind')
221   format('0     the quadratic prod map for pces'/6x,'pcn',2x,'rxn',2x,'rsn',2x,'rsn',2x,'ind')
230   format(5x,'Solution species')
231   format(6x,'(',i3,')',2x,a16)
232   format(5x,'Invariant species')
235   format(5x,'Relationships')
236   format(5x,'Column integrals')
237   format(5x,'Groups')
238   format(3x,'(',i2,')',2x,a16,' - ',1pe10.3)

502   format(10i4)
504   format(2i4)
506   format(3i4)
508   format(4i4)
510   format(5i4)
512   format(6i4)
514   format(i4)
516   format(1x,2e10.4)
519   format(5e16.8)
522   format(20i4)

600   format('0  upper bndy conds'/2x,'species  d  n')
602   format(1x,a16,2i3)
604   format('0  lower bndy conds'/2x,'species  d  n')
606   format('0  upper bndy flux'/2x,'species     day        night   ')
608   format(1x,a16,1p,2e12.4)
610   format('0  lower bndy flux'/2x,'species     day        night   ')
612   format('0  upper bndy dir constants'/2x,'species',5x,'day',8x,'night   ')
614   format('0  lower bndy dir constants'/2x,'species',5x,'day',8x,'night   ')
616   format('0  aust coefficients'/2x,'species',5x,'day',8x,'night   ')
618   format('0  the time increments')
620   format(1x,'(',i2,')',2x,1pe12.4)

1565  format(11x,'|--------------------------------------------------------------------------------------------------|')
1571  format(11x,'| ',96x,' |')
1567  format(11x,'| ',8x,a80,8x,' |')
1569  format(11x,'|**************************************************************************************************|')

2502  format('0  rxn      a0          b0')
2504  format(3x,i3,1p,2e12.4)
2506  format('1',14x,'boundary conditions'/'+',14x,'________ __________' &
      /'0',12x,'upper boundary',6x,'lower boundary'/'+',12x, &
      '_____ ________',6x,'_____ ________'/'  species',5x,'day',6x, &
      'night',6x,'day',6x,'night'/'+ _______',5x,'___',6x,'_____', &
      6x,'___',6x,'_____')
2508  format(1x,10a16)
2270  format('0',10a16)

4000  format('0*** group table ***')
4002  format('0   group no',i4)
4004  format(1x,f3.1,i4)
4010  format('0*** column table ***')
4012  format(1x,2i4)
4020  format('0*** printout table ***')
4022  format('0 index',i4,' type',i4,' ic flag',i4)
4024  format('0 st fine prt',1pe12.4,' st course prt',e12.4)
4026  format('0 fine prt grid'/(1x,1pe12.4))
4028  format('0 course prt grid'/(1x,1pe12.4))
4030  format('0 column count'/(1x,i4))
4032  format('0 print directory'/(1x,i4))

      end program mozart_pp
