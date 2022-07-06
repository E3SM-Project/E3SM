
      subroutine JCL( jobctl, &
                      histinp, &
                      histout, &
                      options, &
                      sub_cnt, &
                      imp_cls_cnt, &
                      hostname, &
                      jobname, &
                      machine, &
                      wrk_dir, &
                      subfile, &
                      cpucnt )

      use IO, only : lin, lout

!-------------------------------------------------------
!	... Input arguments
!-------------------------------------------------------
      integer, intent(in) ::        sub_cnt       ! user subroutine count
      integer, intent(in) ::        imp_cls_cnt   ! implicit soln class count
      integer, intent(in) ::        cpucnt        ! cpu count for distributed environs

      character(len=16), intent(in) ::   hostname      ! user machine hostname
      character(len=16), intent(in) ::   jobname       ! unique name for file
      character(len=16), intent(in) ::   jobctl(8)     ! job control parms
      character(len=64), intent(in) ::   histinp(3)    ! "history" tape input parms
      character(len=64), intent(in) ::   histout(6)    ! "history" tape output parms
      character(len=64), intent(in) ::   wrk_dir       ! working directory on target machine
      character(len=64), intent(inout) ::   subfile    ! submission filespec
      character(len=16), intent(in)  ::   machine       ! target machine "name"

      logical, intent(in), target   ::   options(*)    ! run options

!-------------------------------------------------------
!	... Local variables
!-------------------------------------------------------
      integer       :: spos
      real          :: time, days, seconds
      real          :: dtime
      character(len=320) :: buff
      character(len=64)  :: caps
      character(len=64)  :: fpth
      character(len=64)  :: fname
      character(len=6)   :: date
      character(len=1)   :: char
      logical       :: lexist
      logical, pointer :: f90
      logical, pointer :: usemods

!-------------------------------------------------------
!	... Function declarations
!-------------------------------------------------------
      integer    :: STRLEN
      logical    :: ISNUM

!-----------------------------------------------------
!	... The options array has the following mapping:
!
!       (1) Chemistry (on/off)            (2) Target machine == cray (yes/no)
!       (3) Diffusion (on/off)            (4) Convection (on/off)
!       (5) Iter norms (on/off)           (6) Conservation (on/off)
!       (7) Source code (yes/no)          (8) Submission files (yes/no)
!       (9) Execution (yes/no)           (10) SLT fixer (on/off)
!      (11) Multitasking (yes/no)        (12) Rxt rate lookup table (on/off)
!      (13) Relative humidity (yes/no)   (14) New compiler (yes/no)
!      (15) Height field (yes/no )       (16) User "hook" (yes/no)
!      (17) Use f90 modules (yes/no)     (18) - (20) Unused
!-----------------------------------------------------

      f90 => options(14)
      usemods => options(17)
!-------------------------------------------------------
!	... Open "script" file and make script
!-------------------------------------------------------
      if( machine == 'CRAYYMP' ) then
	 if( subfile == ' ' ) then
	    subfile = 'ctm.ymp.job'
	 end if
         INQUIRE( file = subfile, exist = lexist )
         if( lexist ) then
	    call SYSTEM( 'rm ' // subfile(:LEN_TRIM(subfile)) )
         end if
         CLOSE( 3 )
         OPEN( unit = 3, &
               file = subfile, &
               status = 'new' )

         buff = '#QSUB -s /bin/csh'
      else if( machine == 'CRAY3' ) then
	 if( subfile == ' ' ) then
	    subfile = 'ctm.c3.job'
	 end if
         INQUIRE( file  = subfile, exist = lexist )
         if( lexist ) then
	    call SYSTEM( 'rm ' // subfile(:LEN_TRIM(subfile)) )
         end if
         CLOSE( 3 )
         OPEN( unit = 3, &
               file = subfile, &
               status = 'new' )

         buff = '#QSUB -s /bin/csh'
      else if( machine == 'T3D' ) then
	 if( subfile == ' ' ) then
	    subfile = 'ctm.t3d.job'
	 end if
         INQUIRE( file  = subfile, exist = lexist )
         if( lexist ) then
	    call SYSTEM( 'rm ' // subfile(:LEN_TRIM(subfile)) )
         end if
         CLOSE( 3 )
         OPEN( unit = 3, &
               file = subfile, &
               status = 'new' )

         buff = '#QSUB -s /bin/csh'
      else if( machine == 'RS6000' ) then
	 if( subfile == ' ' ) then
	    subfile = 'ctm.rs6000.job'
	 end if
         INQUIRE( file  = subfile, exist = lexist )
         if( lexist ) then
	    call SYSTEM( 'rm ' // subfile(:LEN_TRIM(subfile)) )
         end if
         CLOSE( 3 )
         OPEN( unit = 3, &
               file = subfile, &
               status = 'new' )

         buff = '#! /bin/csh'
      end if
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... Check for account override
!-------------------------------------------------------
      if( jobctl(5) /= ' ' .and. machine /= 'CRAY3' ) then
         buff = '#QSUB -A ' // jobctl(5)(:STRLEN(jobctl(5)))
         write(3,100) buff(:STRLEN(buff))
      end if
!-------------------------------------------------------
!	... "Write" the que
!-------------------------------------------------------
      if( options(2) ) then
         if( jobctl(8) /= ' ' ) then
	    caps = jobctl(8)
	    call UPCASE( caps )
	    if( caps(:4) == 'PREM' ) then
               buff = '#QSUB -q prem'
	    else if( caps(:4) == 'ECON' ) then
               buff = '#QSUB -q econ'
	    else if( caps(:3) == 'REG' ) then
               buff = '#QSUB -q reg'
	    else if( caps(:2) == 'LM' ) then
               buff = '#QSUB -q lm'
	    else
               buff = '#QSUB -q reg'
	    endif
         else
            buff = '#QSUB -q prem'
         end if
         write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... "Write" the cpu limit
!-------------------------------------------------------
         if( options(2)  ) then
	    if( cpucnt > 1 ) then
	       buff = '#QSUB -la '
	       if( machine /= 'CRAY3' ) then
	          write(buff(11:12),'(i2)') MIN( cpucnt,16 )
	       end if
	       buff(13:) = 'cpus'
               write(3,100) buff(:STRLEN(buff))
               buff = ' '
            end if
         end if

!-------------------------------------------------------
!	... "Write" the time limit
!-------------------------------------------------------
         call TIMCON( jobctl(2)(:STRLEN(jobctl(2))), &
                      time, &
                      lout )
         write(buff,'(''#QSUB -lT '',i6)') INT(time)
         write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... "Write" the memory limit
!-------------------------------------------------------
         buff = '#QSUB -lM ' // jobctl(4)(:STRLEN(jobctl(4))) // 'Mw'
         write(3,100) buff(:STRLEN(buff))

         buff = 'set echo'
         write(3,100) buff(:STRLEN(buff))
         buff = 'set timestamp'
         write(3,100) buff(:STRLEN(buff))
         buff = 'date'
         write(3,100) buff(:STRLEN(buff))
      end if
!-------------------------------------------------------
!	... Change to working directory
!-------------------------------------------------------
      if( options(2) ) then
	 if( wrk_dir == ' ' ) then
            buff = 'cd $TMPDIR'
	 else
            buff = 'cd ' // wrk_dir(:LEN_TRIM(wrk_dir))
	 end if
      else if( machine == 'RS6000' ) then
         buff = 'set echo'
         write(3,100) buff(:STRLEN(buff))
         buff = 'set timestamp'
         write(3,100) buff(:STRLEN(buff))
         buff = 'date'
         write(3,100) buff(:STRLEN(buff))
	 buff = 'if( ! -e /usr/tmp/stacy ) then'
         write(3,100) buff(:STRLEN(buff))
	 buff = '   mkdir /usr/tmp/stacy'
         write(3,100) buff(:STRLEN(buff))
	 buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
         buff = 'cd /usr/tmp/stacy'
      end if
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... "Write" the namelist inputs
!-------------------------------------------------------
      buff = 'if ( -e ctm.dat ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   rm ctm.dat'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))
      buff = 'cat > ctm.dat << ''END1'''
      write(3,100) buff(:STRLEN(buff))
      buff = 'ctm off-line (ver 2.0 ) : case = ' // jobctl(6)(:STRLEN(jobctl(6)))
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... "Experiment" definition
!-------------------------------------------------------
      buff = ' &EXPDEF'
      write(3,100) buff(:STRLEN(buff))
      caps = jobctl(7)
      call UPCASE( caps )
      if( caps == 'TRUE' ) then
         buff = ' NSREST = 1,'
         write(3,100) buff(:STRLEN(buff))
      end if
      if( histout(5) /= ' ' ) then
         buff = ' WPASS = ''' // histout(5)(:STRLEN(histout(5))) // ''','
         write(3,100) buff(:STRLEN(buff))
      end if
      if( histout(1) /= ' ' ) then
         call TIMCON( histout(1)(:STRLEN(histout(1))), &
                      time, &
                      lout )
         write(buff,'('' IRT = '',i5,'','')') INT( time/8.64e4 + .01 )
         write(3,100) buff(:STRLEN(buff))
      end if
      buff = ' CASEID = ''' // jobctl(6)(:STRLEN(jobctl(6))) // ''','
      write(3,100) buff(:STRLEN(buff))
      buff = ' /'
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... "Run" definition
!-------------------------------------------------------
      buff = ' &NEWRUN'
      write(3,100) buff(:STRLEN(buff))
      call TIMCON( jobctl(1)(:STRLEN(jobctl(1))), &
                   dtime, &
                   lout )
      write(buff,'('' DTIME = '',i5,'','')') INT(dtime)
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... Simulation duration in time steps
!-------------------------------------------------------
      spos = STRLEN( jobctl(3) )
      char = jobctl(3)(spos:spos)
      if( ISNUM( char ) ) then
         buff = ' NESTEP = ' // jobctl(3)(:spos) // ','
      else
         call TIMCON( jobctl(3)(:spos), &
                      time, &
                      lout )
         write(buff,'('' NESTEP = '',i5,'','')') INT( (time + .01)/dtime )
      end if
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... History tape output frequency
!-------------------------------------------------------
      if( histout(2) /= ' ' ) then
         spos = STRLEN( histout(2) )
         char = histout(2)(spos:spos)
         if( char >= '0' .and. char <= '9' ) then
	    buff = ' NNUMWT = ' // histout(2)(:spos) // ','
	 else
            call TIMCON( histout(2)(:spos), &
                         time, &
                         lout )
            write(buff,'('' NNUMWT = '',i5,'','')') INT( (time + .01)/dtime )
	 end if
      else
	 buff = ' NNUMWT = 0,'
      end if
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... Simulation printout frequency
!-------------------------------------------------------
      if( histout(6) /= ' ' ) then
         spos = STRLEN( histout(6) )
         char = histout(6)(spos:spos)
         if( char >= '0' .and. char <= '9' ) then
	    buff = ' PRFREQ = ' // histout(6)(:spos) // ','
	 else
            call TIMCON( histout(6)(:spos), &
                         time, &
                         lout )
            write(buff,'('' PRFREQ = '',i5,'','')') INT( (time + .01)/dtime )
	 end if
         write(3,100) buff(:STRLEN(buff))
      end if
!-------------------------------------------------------
!	... Output history tape density
!-------------------------------------------------------
      if( histout(4) /= ' ' ) then
	 buff = ' NDENS = ' // histout(4)(:STRLEN(histout(4))) // ','
         write(3,100) buff(:STRLEN(buff))
      end if
      if( options(10) ) then
	 buff = ' LFIXER = .TRUE.,'
         write(3,100) buff(:STRLEN(buff))
	 buff = ' LIMFIX = .TRUE.,'
         write(3,100) buff(:STRLEN(buff))
      end if
      if( histout(3) /= ' ' ) then
         buff = ' STFNUM = ' // histout(3)(:STRLEN(histout(3)))
         write(3,100) buff(:STRLEN(buff))
      end if
      buff = ' /'
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... "Input" definition
!-------------------------------------------------------
      buff = ' &INPUT'
      write(3,100) buff(:STRLEN(buff))
      call PARSE_FLPTH( histinp(1), fname, fpth )
      spos = STRLEN(fpth)
      if( fpth(spos:spos) == '/' ) then
	 fpth(spos:spos) = ' '
      end if
      buff = ' MSPATH = ''' // fpth(:STRLEN(fpth)) // ''','
      write(3,100) buff(:STRLEN(buff))
      buff = ' MSFN = ''' // fname(:STRLEN(fname)) // ''','
      write(3,100) buff(:STRLEN(buff))
      spos = STRLEN(histinp(3))
      if( spos /= 0 ) then
         buff = ' ICFLNM = ''' // histinp(3)(:spos) // ''','
         write(3,100) buff(:STRLEN(buff))
      end if
      call TIMCON_D( histinp(2)(:STRLEN(histinp(2))), &
                     days, &
                     seconds)
      call MKDATE( days, date )
      buff = ' ICDATE = ' // date // ','
      write(3,100) buff(:STRLEN(buff))
      write(buff,'('' ICSEC = '',i5,'','')') INT( seconds )
      write(3,100) buff(:STRLEN(buff))
      buff = ' /'
      write(3,100) buff(:STRLEN(buff))

      buff = '''END1'''
      write(3,100) buff(:STRLEN(buff))

      buff = ' '
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... Write the sim.dat file
!-------------------------------------------------------
      buff = 'if ( -e sim.dat ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   rm sim.dat'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))
      buff = 'cat > sim.dat << ''END1'''
      write(3,100) buff(:STRLEN(buff))
      CLOSE(2)
      OPEN( unit   = 2, &
            file   = 'sim.dat', &
            status = 'old' )
      do
	 read(2,100,end=200) buff
         write(3,100) buff(:STRLEN(buff))
      end do

200   continue
      buff = '''END1'''
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... Write the modules file
!-------------------------------------------------------
      if( usemods ) then
         buff = 'if ( -e ctm.mods.f ) then'
         write(3,100) buff(:STRLEN(buff))
         buff = '   rm ctm.mods.f'
         write(3,100) buff(:STRLEN(buff))
         buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
         buff = 'cat > ctm.mods.f << ''END1'''
         write(3,100) buff(:STRLEN(buff))
         CLOSE(2)
         OPEN( unit   = 2, &
               file   = 'ctm.mods.f', &
               status = 'old' )
         do
	    read(2,100,end=211) buff
            write(3,100) buff(:STRLEN(buff))
         end do
211      continue
         buff = '''END1'''
         write(3,100) buff(:STRLEN(buff))
      end if

!-------------------------------------------------------
!	... Write the main file
!-------------------------------------------------------
      buff = 'if ( -e ctm.main.f ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   rm ctm.main.f'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))
      buff = 'cat > ctm.main.f << ''END1'''
      write(3,100) buff(:STRLEN(buff))
      CLOSE(2)
      OPEN( unit   = 2, &
            file   = 'ctm.main.f', &
            status = 'old' )
      do
	 read(2,100,end=210) buff
         write(3,100) buff(:STRLEN(buff))
      end do

210   continue
      buff = '''END1'''
      write(3,100) buff(:STRLEN(buff))

      if( options(7) .or. sub_cnt /= 0 ) then
!-------------------------------------------------------
!	... Write the subroutine file
!-------------------------------------------------------
         buff = 'if ( -e ctm.subs.f ) then'
         write(3,100) buff(:STRLEN(buff))
         buff = '   rm ctm.subs.f'
         write(3,100) buff(:STRLEN(buff))
         buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
         buff = 'cat > ctm.subs.f << ''END1'''
         write(3,100) buff(:STRLEN(buff))
         CLOSE(2)
         OPEN( unit   = 2, &
               file   = 'ctm.subs.f', &
               status = 'old' )
         do
	    read(2,100,end=220) buff
            write(3,100) buff(:STRLEN(buff))
         end do

220      continue
         buff = '''END1'''
         write(3,100) buff(:STRLEN(buff))
      end if

!-------------------------------------------------------
!	... Compile the sources
!-------------------------------------------------------
      buff = ' '
      write(3,100) buff(:STRLEN(buff))
      if( options(2) ) then
!-------------------------------------------------------
!	... Compile for a Cray ( not T3D )
!-------------------------------------------------------
         buff = 'ja'
         write(3,100) buff(:STRLEN(buff))
	 if( machine /= 'CRAY3' ) then
	    if( usemods ) then
               buff = 'f90 -c ctm.mods.f'
               write(3,100) buff(:STRLEN(buff))
               buff = 'if ( $status ) then'
               write(3,100) buff(:STRLEN(buff))
               buff = '   echo ctm.mods.f compile failed'
               write(3,100) buff(:STRLEN(buff))
               buff = '   goto errexit'
               write(3,100) buff(:STRLEN(buff))
               buff = 'endif'
               write(3,100) buff(:STRLEN(buff))
	    end if
            buff = 'f90 -c ctm.main.f'
	 else
	    if( options(11) ) then
               buff = '/u0/cs/bin/f77 -c -h fmp ctm.main.f'
	    else
               buff = '/u0/cs/bin/f77 -c ctm.main.f'
	    end if
	 end if
         write(3,100) buff(:STRLEN(buff))
      else if( machine == 'RS6000' ) then
         buff = 's2d ctm.main.f > Main.f'
         write(3,100) buff(:STRLEN(buff))
         buff = 'if ( $status ) then'
         write(3,100) buff(:STRLEN(buff))
         buff = '   echo ctm.main.f s2d failed'
         write(3,100) buff(:STRLEN(buff))
         buff = '   goto errexit'
         write(3,100) buff(:STRLEN(buff))
         buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
         buff = 's2d ctm.subs.f > Subs.f'
         write(3,100) buff(:STRLEN(buff))
         buff = 'if ( $status ) then'
         write(3,100) buff(:STRLEN(buff))
         buff = '   echo ctm.subs.f s2d failed'
         write(3,100) buff(:STRLEN(buff))
         buff = '   goto errexit'
         write(3,100) buff(:STRLEN(buff))
         buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
         buff = 'xlf -c -O3 Main.f'
         write(3,100) buff(:STRLEN(buff))
      end if
      buff = 'if ( $status ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   echo ctm.main.f compile failed'
      write(3,100) buff(:STRLEN(buff))
      buff = '   goto errexit'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))
      
      if( options(7) .or. sub_cnt /= 0 ) then
	 if( options(2) ) then
	    if( machine /= 'CRAY3' ) then
               buff = 'f90 -c ctm.subs.f'
	    else
               buff = '/u0/cs/bin/f77 -c ctm.subs.f'
	    end if
	 else if( machine == 'RS6000' ) then
            buff = 'xlf -c -O3 Subs.f'
	 end if
         write(3,100) buff(:STRLEN(buff))
         buff = 'if ( $status ) then'
         write(3,100) buff(:STRLEN(buff))
         buff = '   echo ctm.subs.f compile failed'
         write(3,100) buff(:STRLEN(buff))
         buff = '   goto errexit'
         write(3,100) buff(:STRLEN(buff))
         buff = 'endif'
         write(3,100) buff(:STRLEN(buff))
      end if

!-------------------------------------------------------
!	... Form the executable
!-------------------------------------------------------
      buff = ' '
      write(3,100) buff(:STRLEN(buff))
      if( options(2) ) then
         if( options(7) ) then
	    if( imp_cls_cnt /= 0 .and. machine /= 'CRAY3' ) then
	       if( f90 ) then
	          buff = 'f90 -Wl"-f indef" ctm.main.o ctm.subs.o /crestone/u2/stacy/ctm/lib/fsim.o \\'
	       else
	          buff = 'segldr -f indef ctm.main.o ctm.subs.o /crestone/u2/stacy/ctm/lib/fsim.o \\'
	       end if
	    else
	       buff = 'segldr -f indef ctm.main.o ctm.subs.o \\ '
	    end if
            write(3,100) buff(:STRLEN(buff))
	    if( f90 ) then
	       buff = '       -L /usr/local/lib -lncaro,hpf,mss'
	    else
	       buff = '       -L /lib,/usr/lib,/usr/local/lib \\'
	    end if
            write(3,100) buff(:STRLEN(buff))
	    if( machine /= 'CRAY3' ) then
	       if( .not. f90 ) then
	          buff = '       -l ncaro,hpf,mss -M ,s'
	       end if
	    else
	       if( options(11) ) then
	          buff = '       -l ncaro,hpf,mss,net,pll -M ,s'
	       else
	          buff = '       -l ncaro,hpf,mss,net -M ,s'
	       end if
	    end if
         else if( sub_cnt == 0 ) then
	    if( imp_cls_cnt /= 0 .and. machine /= 'CRAY3' ) then
	       buff = 'segldr -f indef ctm.main.o /crestone/u2/stacy/ctm/lib/fsim.o \\'
	    else
	       buff = 'segldr -f indef ctm.main.o \\'
	    end if
            write(3,100) buff(:STRLEN(buff))
	    buff = '       -L /lib,/usr/lib,/usr/local/lib,/u0/stacy/ctm/lib \\'
            write(3,100) buff(:STRLEN(buff))
	    buff = '       -l ncaro,hpf,mss,ctm -M ,s'
         else
	    buff = 'segldr -f indef ctm.main.o ctm.subs.o \\ '
            write(3,100) buff(:STRLEN(buff))
	    buff = '       -L /lib,/usr/lib,/usr/local/lib,/u2/stacy/ctm/lib \\'
            write(3,100) buff(:STRLEN(buff))
	    buff = '       -l ncaro,hpf,mss,ctm -M ,s'
         end if
      else if( machine == 'RS6000' ) then
	 if( cpucnt == 1 ) then
            buff = 'xlf Main.o Subs.o -L/usr/local/lib -lmss -lncaru -lessl'
	 else
            buff = 'xlf Main.o Subs.o -L /usr/lpp/pvm/lib -lpvm -lf2c \\'
            write(3,100) buff(:STRLEN(buff))
            buff = '                  -L/usr/local/lib -lmss -lncaru -lessl \\'
            write(3,100) buff(:STRLEN(buff))
            buff = '                  -bI:/usr/lpp/pvm/lib/pvme.exp'
	 end if
      end if
      if( .not. f90 ) then
         write(3,100) buff(:STRLEN(buff))
      end if
      buff = 'if ( $status ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   echo segldr failed'
      write(3,100) buff(:STRLEN(buff))
      buff = '   goto errexit'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... Execution
!-------------------------------------------------------
      buff = ' '
      write(3,100) buff(:STRLEN(buff))
      buff = 'if ( -e ctm.out.$$ ) then'
      write(3,100) buff(:STRLEN(buff))
      buff = '   rm ctm.out.$$'
      write(3,100) buff(:STRLEN(buff))
      buff = 'endif'
      write(3,100) buff(:STRLEN(buff))
      buff = ' '
      write(3,100) buff(:STRLEN(buff))
      if( options(2)  ) then
	 if( cpucnt > 1 ) then
	    buff = 'setenv NCPUS '
	    if( machine /= 'CRAY3' ) then
	       write(buff(14:),'(i2)') MIN( cpucnt,16 )
	    else
	       write(buff(14:),'(i1)') MIN( cpucnt,2 )
	    end if
            write(3,100) buff(:STRLEN(buff))
            buff = ' '
         end if
      end if
      buff = 'a.out < ctm.dat > ctm.out.$$'
      write(3,100) buff(:STRLEN(buff))

!-------------------------------------------------------
!	... Disperse printouts; normal termination
!-------------------------------------------------------
      buff = ' '
      write(3,100) buff(:STRLEN(buff))
      buff = 'rcp ctm.out.$$  \\'
      write(3,100) buff(:STRLEN(buff))
      if( hostname /= 'acd' ) then
         buff = '    ' // hostname(:STRLEN(hostname)) // '.acd.ucar.edu:rje/ctm.out.$$'
      else
         buff = '         acd.ucar.edu:rje/ctm.out.$$'
      end if
      write(3,100) buff(:STRLEN(buff))
      if( options(2) ) then
	 if( machine /= 'CRAY3' ) then
            buff = 'ja -schflt > accnting.$$'
	 else
            buff = 'ja -scflt > accnting.$$'
	 end if
         write(3,100) buff(:STRLEN(buff))
         buff = 'rcp  accnting.$$ \\'
         write(3,100) buff(:STRLEN(buff))
         if( hostname /= 'acd' ) then
            buff = '    ' // hostname(:STRLEN(hostname)) // '.acd.ucar.edu:rje/accnting.$$'
         else
            buff = '              acd.ucar.edu:rje/accnting.$$'
         end if
         write(3,100) buff(:STRLEN(buff))
         if( hostname /= 'acd' .and. options(9) ) then
	    buff = 'cd $home'
            write(3,100) buff(:STRLEN(buff))
	    buff = 'echo "rcp /usr/tmp/O' // jobname(:STRLEN(jobname)) // ' ' // hostname(:STRLEN(hostname))  &
                   // '.acd.ucar.edu:rje/O' // jobname(:STRLEN(jobname)) // '" | at now + 2 minute'
            write(3,100) buff(:STRLEN(buff))
         end if
      end if
      buff = 'exit( 0 )'
      write(3,100) buff(:STRLEN(buff))

      buff = ' '
      write(3,100) buff(:STRLEN(buff))
!-------------------------------------------------------
!	... Disperse printouts; error termination
!-------------------------------------------------------
      buff = 'errexit:'
      write(3,100) buff(:STRLEN(buff))
      buff = 'rcp ctm.out.$$ \\'
      write(3,100) buff(:STRLEN(buff))
      if( hostname /= 'acd' ) then
         buff = '    ' // hostname(:STRLEN(hostname)) // '.acd.ucar.edu:rje/ctm.out.$$'
      else
         buff = '         acd.ucar.edu:rje/ctm.out.$$'
      end if
      write(3,100) buff(:STRLEN(buff))
      if( options(2) ) then
	 if( machine /= 'CRAY3' ) then
            buff = 'ja -schflt > accnting.$$'
	 else
            buff = 'ja -scflt > accnting.$$'
	 end if
         write(3,100) buff(:STRLEN(buff))
         buff = 'rcp  accnting.$$ \\'
         write(3,100) buff(:STRLEN(buff))
         if( hostname /= 'acd' ) then
         buff = '    ' // hostname(:STRLEN(hostname)) // '.acd.ucar.edu:rje/accnting.$$'
         else
            buff = '              acd.ucar.edu:rje/accnting.$$'
         end if
         write(3,100) buff(:STRLEN(buff))
         if( hostname /= 'acd' .and. options(9) ) then
	    buff = 'cd $home'
            write(3,100) buff(:STRLEN(buff))
	    buff = 'echo "rcp /usr/tmp/O' // jobname(:STRLEN(jobname)) // ' ' // hostname(:STRLEN(hostname))  &
                   // '.acd.ucar.edu:rje/O' // jobname(:STRLEN(jobname)) // '" | at now + 2 minute'
            write(3,100) buff(:STRLEN(buff))
         end if
      end if
      buff = 'exit( -1 )'
      write(3,100) buff(:STRLEN(buff))
      CLOSE(3)

100   format(a)

      end
