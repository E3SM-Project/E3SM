program piocprnc
  use netcdf
  use filestruct
  use compare_vars_mod
#ifdef NAGFOR
  use f90_unix
#endif
  implicit none

  integer :: nargs, n
  character(len=1024) :: arg = ''                     ! cmd-line argument
  character(len=1024) :: fname(2) = ' '          ! input filenames
  integer :: nchars
  integer :: numcases=1
  integer :: ierr
!  integer, external :: iargc
  type(file_t) :: file(2)
  type(dim_t) :: dimoptions(12)
  integer :: dimoptioncnt
  integer :: nvars, ndiffs, nfilldiffs, numnotfound
  integer :: num_sizes_differ
  integer :: num_not_analyzed

  n = 1
  do while(n <= 12)
     dimoptions(n)%name = ''
     n = n + 1
  end do

!
! Parse arg list
!


   nargs = command_argument_count ()
   dimoptioncnt=0
   ignoretime=.false.
   n = 1
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      select case (arg)
      case ('-v')
         verbose = .true.
      case ('-d')
         call getarg(n, arg)
         n=n+1
         dimoptioncnt=dimoptioncnt+1
         call parsearg(arg, dimoptions(dimoptioncnt)%name, dimoptions(dimoptioncnt)%start, dimoptions(dimoptioncnt)%kount)

      case ('-m')
         ignoretime=.true.
      case default
         if (fname(1) == ' ') then
            fname(1) = arg(1:len_trim(arg))
            nchars = len_trim (fname(1))
            write (6,*) 'file 1=',fname(1)(1:nchars)
         else if (fname(2)==' ') then
            fname(2) = arg(1:len_trim(arg))
            nchars = len_trim (fname(2))
            write (6,*) 'file 2=',fname(2)(1:nchars)
            numcases = 2
         else
            call usage_exit (' ')
         end if
      end select
   end do
!
! Must have at least 1 file input
!
   if (fname(1) == ' ') then
      call usage_exit ('You must enter at least 1 input file')
   end if

!
!  Read the files and initialize file_t
!
   do n=1, numcases
      ierr = nf90_open(fname(n),NF90_NOWRITE, file(n)%fh)
      if(ierr /= NF90_NOERR) then
         stop 'Failed to open file '
      endif
      if(dimoptioncnt>0) then
         call init_file_struct( file(n), dimoptions(1:dimoptioncnt) )
      else
         call init_file_struct( file(n))
      end if
   end do

   if(numcases==2) then
      call compare_metadata(file(1), file(2))

      call compare_dimensions( file(1)%dim, file(2)%dim)

      call match_vars( file(1), file(2) )
   end if
   call compare_vars(numcases, file, nvars, ndiffs, nfilldiffs, &
        num_sizes_differ, num_not_analyzed, numnotfound)


!
! Summarize results
!
     write(6,806)
      write(6,*) ' '
      write(6,700) 'SUMMARY of cprnc:'
      if(numcases==1) then
         write(6,700) '  A total number of ',nvars,' fields in file 1 were analyzed (non-compare mode)'
         write(6,700) '  A total number of ',num_not_analyzed, &
              ' fields in file 1 could not be analyzed'
      else
         write(6,700) ' A total number of ',nvars,' fields were compared'
         write(6,700) '          of which ',ndiffs,' had non-zero differences'
         write(6,700) '               and ',nfilldiffs,' had differences in fill patterns'
         write(6,700) '               and ',num_sizes_differ,' had different dimension sizes'
         write(6,700) ' A total number of ',num_sizes_differ + num_not_analyzed, &
              ' fields could not be analyzed'
         write(6,700) ' A total number of ',numnotfound,' fields on file 1 were not found on file2.'
         if (nvars > 0 .and. ndiffs == 0 .and. nfilldiffs == 0 .and. &
              num_sizes_differ == 0 .and. num_not_analyzed<nvars) then
            write(6,700) '  diff_test: the two files seem to be IDENTICAL '
         else
            write(6,700) '  diff_test: the two files seem to be DIFFERENT '
         end if
      end if
      write(6,*) ' '
700   format(a,i6,a)
806   format(132('*'))



 contains
   subroutine usage_exit (arg)
     implicit none

     character(len=*), intent(in) :: arg

     if (arg /= ' ') write (6,*) arg
     write(6,*)'Usage: cprnc  [-m] [-v] [-d dimname:start[:count]] file1 [file2]'
     write(6,*)'-v: Verbose output'
     write(6,*)'-m: Ignore time variable and just match contents (default is to match the values in variable time.)'
     write(6,*)'-d dimname:start[:count]: Print variable values for the specified dimension index start and count. If not present,'
     write(6,*)'                          count will default to 1.  If count is < 0 then count will be set to dimsize-start'
     write(6,*)'                           '

     stop 999
   end subroutine usage_exit


   subroutine parsearg (arg, dimname, v1, v2)
     !-------------------------------------------------------------------------------------------
     ! Purpose: Parse cmd line args about printing.
     !
     ! Method:  Input is expected in the form: dimname:number1[:number2] where dimname is expected to
     !           be the name of a dimension in the input file(s), number1 is the starting position in that
     !           dimension to be evaluated and number2 is the number of values to read in the dimension
     !           if number2 is missing all remaining values are read.
     !
     !-------------------------------------------------------------------------------------------
     implicit none

     character(len=*), intent(in) :: arg    ! cmd line arg expected of the form 'num1:num2' or 'num1'

     character(len=*), intent(out) :: dimname
     integer, intent(out) :: v1             ! e.g. num1 from above example
     integer, intent(out) :: v2             ! e.g. num2 from above example

     integer :: i, j                        ! indices through arg
     integer :: ierr                        ! io error status

     !
     ! First get a dimension name
     !
     dimname = ' '
     i = scan(arg,':')
     dimname(1:i-1)=arg(1:i-1)
     i=i+1

     !
     ! now try to get an integer number for everything up to ":"
     !
     j=i
     do while (j < len(arg) .and. arg(j:j) >= '0' .and. arg(j:j) <= '9' .and. arg(j:j) /= ':')
        j = j + 1
     end do
     read (arg(i:j-1), '(i5)') v1
     !
     ! Next, if ":" comes after the number, look for the next number
     !
     i=j

     if (arg(i:i) == ':') then
        j = i + 1
        do while (j < len(arg) .and. scan(arg(j:j),"-0123456789")>0)
           j = j + 1
        end do
        read (arg(i+1:j-1), '(i5)', iostat=ierr) v2
        !
        ! On unexpected input set v2 = -1, e.g. "-d lon:2:blah" will mean get all lons > 1
        !
        if (ierr /= 0) then
           v2 = -1
        end if
     else
        !
        ! ":" not present. Interpret for example '-d lon:2' to mean '-d lon:2:1'
        !
        v2 = 1
     end if
     if(verbose) print *,__FILE__,__LINE__,trim(dimname),v1,v2
     return
   end subroutine parsearg

 end program piocprnc
