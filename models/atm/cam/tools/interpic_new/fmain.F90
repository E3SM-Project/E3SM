program fmain

   use control,    only: verbose, silent, prec_override, prec_out, &
                         set_user_skip, set_user_include
   use dimensions, only: is_ewdim, is_nsdim, is_zdim, is_ncoldim, ncoldimid

   implicit none

   include 'netcdf.inc'

   ! Local workspace

   character(len=256)  :: arg, file1, file2, template
   character(len=256)  :: exclude = ' ', inc_template = ' '
   character(len=1024) :: cmdline
   character(len=nf_max_name) :: name, namei
   character(len=nf_max_name) :: attname

   integer :: attlen
   integer :: ncidi, ncido, ncidt
   integer :: ndims, ndimsi
   integer :: unlimdimid, unlimdimidi
   integer :: dimlen
   integer :: dimid
   integer :: nvars, nvarsi
   integer :: ngatts, ngattsi
   integer :: ret
   integer :: n
   integer :: ntime = 1
   integer :: nargs
   integer :: nprec

   integer iargc
   external iargc

   ! Default settings before parsing argument list

   file1 = ' '
   file2 = ' '
   template = ' '
   verbose = .false.
   silent = .false.

   nargs = iargc()
   n = 1
   cmdline = 'interpic '
   do while (n .le. nargs)
      arg = ' '
      call getarg(n, arg)
      n = n + 1

      select case (arg)
      case ('-e')
         call getarg(n, arg)
         n = n + 1
         exclude = arg
         cmdline = trim(cmdline) // ' -e ' // trim(arg)
      case ('-i')
         call getarg(n, arg)
         n = n + 1
         inc_template = arg
         cmdline = trim(cmdline) // ' -i ' // trim(arg)
      case ('-p')
         prec_override = .true.
         call getarg(n, arg)
         n = n + 1
         read(arg,'(i1)') nprec
         cmdline = trim(cmdline) // ' -p ' // trim(arg)
      case ('-s')
         silent = .true.
         cmdline = trim(cmdline) // ' -s'
      case ('-t')
         call getarg(n, arg)
         n = n + 1
         template = arg
         cmdline = trim(cmdline) // ' -t ' // trim(template)
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         if (file1 .eq. ' ') then
            file1 = arg
         else if (file2 .eq. ' ') then
            file2 = arg
         else
            write (6,*) 'Argument ', arg,' is not known'
            call usage_exit(' ')
         end if
         cmdline = trim(cmdline) // ' ' // trim(arg)
      end select
   end do
  
   if (file1 .eq. ' '  .or.  file2 .eq. ' ') then
      call usage_exit('Must enter an input file and an output file')
   else if (silent .and. verbose) then
      call usage_exit('-s cannot be specified with -v')
   end if

   if (len_trim(exclude) > 0) then
      call set_user_skip(exclude)
   end if

   if (len_trim(inc_template) > 0) then
      call set_user_include(inc_template)
   end if

   if (prec_override) then
      if (nprec == 4) then
         prec_out = NF_FLOAT
      else if (nprec == 8) then
         prec_out = NF_DOUBLE
      else
         call usage_exit('-p can only be set to 4 or 8')
      end if
   end if

   ! Open input and output netcdf files
   call wrap_open(file1, NF_NOWRITE, ncidi)
   call wrap_open(template, NF_NOWRITE, ncidt)
   call wrap_create (file2, IOR(NF_CLOBBER,NF_64BIT_OFFSET), ncido)

   ! Determine space dimensions of output file from template file.
   ! Currently all spatial dimensions of the output file must exist
   ! on the template file.
   call wrap_inq(ncidt, ndims, nvars, ngatts, unlimdimid)

   do n=1,ndims

      call wrap_inq_dim(ncidt, n, name, dimlen)
    
      if (is_ncoldim(name) .or. is_ewdim(name) .or. &
          is_nsdim(name) .or. is_zdim (name) ) then

         if (nf_def_dim(ncido, name, dimlen, dimid) /= NF_NOERR) then
            print *, 'fmain: ERROR return from nf_def_dim: failed to define '//trim(name)
         end if

         if (verbose) then
            print *, 'define dim from template file: '//trim(name)//' len=',dimlen
         end if

         if (is_ncoldim(name)) ncoldimid = dimid

      end if

   end do

   ! Determine all non-spatial dimensions from the input file
   call wrap_inq(ncidi, ndimsi, nvarsi, ngattsi, unlimdimidi)

   do n=1,ndimsi

      call wrap_inq_dim(ncidi, n, name, dimlen)

      if (.not. is_ncoldim(name) .and. .not. is_ewdim(name) .and. &
          .not. is_nsdim(name) .and. .not. is_zdim (name) ) then

         ! check for time dimension
         if (trim(name) == 'time') then
            ntime = dimlen

            ! Reset dimlen for unlimited dimension
            ! make sure time is marked UNLIMITED.  netcdf3, 64bit, only allows
            ! 4GB variables.  if time is UNLIMITED, each timeslice limited to 4GB
            ! otherwise, entire variable must be < 4GB
            dimlen = NF_UNLIMITED
            if (n /= unlimdimidi) then
               print *, 'fmain: INFO: change the time dimension from fixed length to unlimited'
            end if

         end if

         if (nf_def_dim(ncido, name, dimlen, dimid) /= NF_NOERR) then
            print *, 'fmain: ERROR return from nf_def_dim: failed to define '//trim(name)
         end if

         if (verbose) then
            print *, 'define dim from input file: '//trim(name)//' len=',dimlen
         end if

      end if
     
   end do

   ! Copy global attributes from input to output file
   do n=1,ngattsi

      call wrap_inq_attname(ncidi, NF_GLOBAL, n, attname)

      if (verbose) then
         write(6,*) 'copy '//trim(attname)
      end if

      call wrap_copy_att(ncidi, NF_GLOBAL, attname, ncido, NF_GLOBAL)

   end do

   ! Add global attributes for interpic
   call addglobal(ncido, cmdline)

   ! Special cases: spatial coordinate variables and variables that depend on spatial
   ! dimensions will be copied from the template file to output file
   if (verbose) then
      write(6,*) 'handle_special_cases from template file:'
   end if
   call handle_special_cases(ncidt, ncido)

   ! Copying spatial coordinate variables and variables that depend on spatial
   ! dimensions from the input file to output is currently disabled.
   !
   !if (verbose) then
   !   write(6,*) 'handle_special_cases from input file:'
   !end if
   !call handle_special_cases(ncidi, ncido)

   ! Call driver code to do the interpolations and/or copies
   call driver(ncidi, ncido, ncidt, ntime)

   if (nf_close(ncido) /= nf_noerr) then
      call err_exit('error from nf_close')
   end if
   if (verbose) then
      write(6,*) 'Finished.'
   end if

end program fmain

subroutine usage_exit (arg)
   implicit none
   character*(*) arg

   if (arg.ne.' ') write (6,*) arg
   write (6,*) 'Usage: interpic [-e "var1,...,varn"] [-i "var1,...,varn"] &
        &[-p (4 or 8)] [-s] [-v]'
   write (6,*) '                -t template infile outfile'
   write (6,*) 'OPTIONS:'
   write (6,*) '-e   comma separated list of variables to exclude when &
        &copying from input file to output file'
   write (6,*) '-i   comma separated list of variables to include when &
        &copying from template file to output file'
   write (6,*) '-p   set precision of data variables in output file'
   write (6,*) '-s   set silent mode'
   write (6,*) '-v   set verbose mode'
   stop 999
end subroutine usage_exit
