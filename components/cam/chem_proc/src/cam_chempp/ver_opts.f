
      module mo_ver_opts

      private
      public :: ver_opts

      contains

      subroutine ver_opts( options, model, machine, march, arch_type, &
                           wrk_dir, cpp_dir, cpp_opts, subfile, diagprnt, &
                           tavgprnt, cpucnt, vec_ftns )
!-----------------------------------------------------------------------
!   	... Set the simulation options
!-----------------------------------------------------------------------

      use io

      implicit none

!-----------------------------------------------------------------------
!   	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(inout) ::      cpucnt
      character(len=16), intent(inout)  ::  machine
      character(len=16), intent(inout)  ::  model
      character(len=16), intent(inout)  ::  march
      character(len=16), intent(inout)  ::  arch_type
      character(len=64), intent(out)   ::  wrk_dir
      character(len=64), intent(inout) ::  cpp_dir, cpp_opts
      character(len=64), intent(out)   ::  subfile
      logical, intent(inout)           ::  diagprnt
      logical, intent(inout)           ::  tavgprnt
      logical, intent(inout)           ::  vec_ftns
      logical, intent(inout)           ::  options(:)

!-----------------------------------------------------------------------
!   	... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: maxparm = 26

      integer      :: kpar, nchar, k
      integer      :: err
      logical      :: entered(maxparm)

      character(len=20) :: keywrd
      character(len=20) :: parkey(maxparm)

      logical      :: found

      integer  ::  lenof

      parkey(1) = 'MACHINE'
      parkey(2) = 'DIFFUSION'
      parkey(3) = 'CONVECTION'
      parkey(4) = 'NORMS'
      parkey(5) = 'CONSERVATION'
      parkey(6) = 'SOURCECODE'
      parkey(7) = 'CPUS'
      parkey(8) = 'MULTITASK'
      parkey(9) = 'FIXER'
      parkey(10) = 'DIAGPRNT'
      parkey(11) = 'RXTNLOOKUP'
      parkey(12) = 'RELHUM'
      parkey(13) = 'F90'
      parkey(14) = 'GEOHEIGHT'
      parkey(15) = 'USERHOOK'
      parkey(16) = 'MODULES'
      parkey(17) = 'WORKDIR'
      parkey(18) = 'NAMEMOD'
      parkey(19) = 'SUBFILE'
      parkey(20) = 'TAVGPRNT'
      parkey(21) = 'VEC_FTNS'
      parkey(22) = 'ARCHITECTURE'
      parkey(23) = 'CPP_DIR'
      parkey(24) = 'CPP_OPTS'
      parkey(25) = 'MODEL'
      parkey(26) = 'MODEL_ARCHITECTURE'

      entered = .false.

!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      do
         call cardin( lin, buff, nchar )
	 buffh = buff
	 call upcase ( buffh )
         if( buffh == 'ENDVERSIONOPTIONS' ) then
	    if( .not. options(1) ) then       ! not a cray target
	       if( cpucnt > 1 ) then
	          options(10) = .true.        ! "distributed" processing
	       else
	          options(10) = .false.       ! no dist processing
	       end if
	    end if
	    if( .not. options(13) ) then      ! if not fortran90 then no modules
	       options(16) = .false.
	       options(17) = .false.
	    end if
	    if( machine /= 'IBM' ) then
	       vec_ftns = .false.
	    end if
	    if( machine == 'NEC' ) then
               march = 'VECTOR'
	    else if( machine == 'INTEL' ) then
               march = 'SCALAR'
	    end if
	    if( model == 'WRF' ) then
               march = 'SCALAR'
	    end if
            exit
	 end if
	 k = index( buffh(:nchar), '=' )
         if( k /= 0 ) then
	    keywrd = buffh(:k-1)
	    found  = .false.
            do kpar = 1,maxparm
               if( keywrd == parkey(kpar) ) then
		  found = .true.
	          exit
	       end if
	    end do
	 else
            call errmes ( ' option specification has no = operator@', lout, buff, 1, buff )
         end if
	 if( .not. found ) then
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call errmes ( ' # is an invalid options' &
                       // ' parameter keyword@', lout, keywrd, &
                          LENOF(20,keywrd), buff )
         end if

!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
         if( entered(kpar) ) then
            call errmes( '0 *** # has already been specified@', &
                          lout, parkey(kpar), k, ' ' )
         end if

!-----------------------------------------------------------------------
!     	... Set individual options
!-----------------------------------------------------------------------
         if( kpar == 1 ) then
	    machine = buffh(k+1:nchar)
	    if( machine /= 'CRAYYMP' .and. machine /= 'CRAY2' .and. machine /= 'CRAY3' &
	                             .and. machine /= 'J90' .and. machine /= 'C90' ) then
	       options(1) = .false.
	    end if
         else if( kpar == 6 ) then
	    if( buffh(k+1:nchar) /= 'FULL' ) then
	       options(6) = .false.
	    end if
         else if( kpar == 7 ) then
	    call intcon( buffh(k+1:), &
                         nchar - k, &
                         cpucnt, &
                         err )
	    if( err /= 0 ) then
	       call errmes( ' # is not a valid number@', &
                            lout, &
                            buffh(k+1:), &
                            nchar - k, &
                            buff )
	    end if
         else if( kpar == 21 ) then
	    if( buffh(k+1:nchar) == 'ON' .or. buffh(k+1:nchar) == 'YES' ) then
	       vec_ftns = .true.
	    end if
         else if( kpar == 22 ) then
	    arch_type = buffh(k+1:nchar)
	    if( arch_type /= 'OMP' .and. arch_type /= 'MPI' .and. arch_type /= 'HYBRID' ) then
               call errmes( '0 *** # is not a valid architecture type@', lout, arch_type, 8, ' ' )
	    end if
         else if( kpar == 25 ) then
	    model = buffh(k+1:nchar)
	    if( model /= 'MOZART' .and. model /= 'CAM' .and. model /= 'WRF' ) then
               call errmes( '# is not a valid model type@', lout, buff(k+1:), nchar-k , ' ' )
	    end if
         else if( kpar == 26 ) then
	    march = buffh(k+1:nchar)
	    if( march /= 'SCALAR' .and. march /= 'CACHE' .and. march /= 'VECTOR' ) then
               call errmes( '# is not a valid model architecture type@', lout, buff(k+1:), nchar-k , ' ' )
	    end if
	 else
	    if( buffh(k+1:nchar) == 'ON' .or. buffh(k+1:nchar) == 'YES' ) then
	       if( kpar == 10) then
		  diagprnt = .true.
	       else if( kpar == 20) then
		  tavgprnt = .true.
	       else if( kpar == 18 ) then
	          options(17) = .true.
	       else if( kpar /= 8 .and. kpar /= 9 ) then
	          options(kpar) = .true.
	       else if( kpar == 8 ) then
	          options(10) = .true.
	       else if( kpar == 9 ) then
	          options(9) = .true.
	       end if
	    else
	       if( kpar == 17) then
		  wrk_dir = buff(k+1:nchar)
	       else if( kpar == 23) then
		  cpp_dir = buff(k+1:nchar)
	       else if( kpar == 24) then
		  cpp_opts = buff(k+1:nchar)
	       else if( kpar == 19) then
		  subfile = buff(k+1:nchar)
	       else if( kpar == 10) then
		  diagprnt = .false.
	       else if( kpar == 20) then
		  tavgprnt = .false.
	       else if( kpar /= 8 .and. kpar /= 9 ) then
	          options(kpar) = .false.
	       else if( kpar == 8 ) then
	          options(10) = .false.
	       else if( kpar == 9 ) then
	          options(9) = .false.
	       end if
	    end if
	 end if
	 entered(kpar) = .true.
      end do

      end subroutine ver_opts

      end module mo_ver_opts
