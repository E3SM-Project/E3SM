program fmain
  use control
  use dimensions
!
! $Id$
!
  implicit none

  include 'netcdf.inc'
!
! Local workspace
!
  character*80 :: arg, file1, file2, template
  character*1500 :: text
  character*256 :: cmdline
  character*(nf_max_name) :: name, namei
  character*(nf_max_name) :: attname

  integer, parameter :: maxny = 1000

  integer attlen
  integer ncidi, ncido, ncidt
  integer ndims, ndimsi
  integer dimlen
  integer dimid
  integer nvars, nvarsi
  integer ngatts, ngattsi
  integer ret
  integer n
  integer ntime
  integer nargs
  integer nx, ny
  data nx, ny /0, 0/

  integer iargc
  external iargc
!
! Default settings before parsing argument list
!
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
    call getarg (n, arg)
    n = n + 1

    select case (arg)
    case ('-s')
      silent = .true.
      cmdline = trim(cmdline) // ' -s'
    case ('-t')
      call getarg (n, arg)
      n = n + 1
      template = arg
      cmdline = trim(cmdline) // ' -t ' // trim(template)
    case ('-v')
      verbose = .true.
      cmdline = trim(cmdline) // ' -v'
    case ('-x')
      call getarg (n, arg)
      n = n + 1
      read(nx,'(i4)') arg
      cmdline = trim(cmdline) // ' -x ' // trim(arg)
    case ('-y')
      call getarg (n, arg)
      n = n + 1
      read(ny,'(i4)') arg
      cmdline = trim(cmdline) // ' -y ' // trim(arg)
    case default
      if (file1 .eq. ' ') then
        file1 = arg
      else if (file2 .eq. ' ') then
        file2 = arg
      else
        write (6,*) 'Argument ', arg,' is not known'
        call usage_exit (' ')
      end if
      cmdline = trim(cmdline) // ' ' // trim(arg)
    end select
  end do
  
  if ((nx > 0 .or. ny > 0) .and. template /= ' ') then
    call usage_exit ('Cannot specify template file in addition to nx or ny')
  else if (file1.eq.' ' .or. file2.eq.' ') then
    call usage_exit ('Must enter an input file and an output file')
  else if (silent .and. verbose) then
    call usage_exit ('-s cannot be specified with -v')
  end if
!
! Open input and output netcdf files
!
  call wrap_open (file1, NF_NOWRITE, ncidi)
  call wrap_open (template, NF_NOWRITE, ncidt)
  call wrap_create (file2, NF_CLOBBER, ncido)
!
! Copy dimension and attribute information from template file to output file
!
  call wrap_inq (ncidt, ndims, nvars, ngatts, unlimdimid)
!
! Determine space and time dimensions of output file from template file
!
  do n=1,ndims
    call wrap_inq_dim (ncidt, n, name, dimlen)

    if (n == unlimdimid) then
      call wrap_def_dim (ncido, name, NF_UNLIMITED, dimid)
    else
      call wrap_def_dim (ncido, name, dimlen, dimid)
      if (is_ewdim (name)) then
        call add_dim (ewdim, name, dimlen)
      else if (is_nsdim (name)) then
        call add_dim (nsdim, name, dimlen)
      else if (is_zdim (name)) then
        call add_dim (zdim, name, dimlen)
      end if
    end if

    if (dimid /= n) then
      call err_exit ('Input dimid not equal to output dimid')
    end if
  end do
!
! Determine space and time dimensions of input file
!
  call wrap_inq (ncidi, ndimsi, nvarsi, ngattsi, unlimdimidi)

  ret = nf_inq_dimlen (ncidi, unlimdimidi, ntime)
  if (ret /= NF_NOERR) then
    ntime = 1
    write(6,*)'INFO: Input file has no unlimited dimension'
  end if

  do n=1,ndimsi
    call wrap_inq_dim (ncidi, n, namei, dimlen)

    if (n /= unlimdimidi) then
      if (is_ewdim (namei)) then
        call add_dim (ewdimi, namei, dimlen)
      else if (is_nsdim (namei)) then
        call add_dim (nsdimi, namei, dimlen)
      else if (is_zdim (namei)) then
        call add_dim (zdimi, namei, dimlen)
      end if
    end if
  end do
!
! Copy global attributes from template file
!
  do n=1,ngatts
    call wrap_inq_attname (ncidt, NF_GLOBAL, n, attname)
    ret = nf_inq_attlen (ncidt, NF_GLOBAL, attname, attlen)

    if (attlen > len(text)) then
      write(6,*) 'attribute ',trim(attname),' too long'
      stop 999
    end if

    if (attname == 'case') then
      text = ' '
      call wrap_get_att_text (ncidt, NF_GLOBAL, attname, text)
      write(6,*)'case =',trim(text)
    else if (attname == 'title') then
      text = ' '
      call wrap_get_att_text (ncidt, NF_GLOBAL, attname, text)
      write(6,*)'title =',trim(text)
    end if
    call wrap_copy_att (ncidt, NF_GLOBAL, attname, ncido, NF_GLOBAL)
  end do

  if (ny > maxny) then
    write(6,*)'maxny too small: recompile with this parameter > ',ny
    stop 999
  end if
!
! Add global attributes for interpic
!
  call addglobal (ncido, cmdline)
!
! Special cases: coordinate variables and offshoots (e.g. nlon, rlon, hyai) will be 
! copied from template file to output file.
!
  call handle_special_cases (ncidt, ncido, nvars, ntime, unlimdimid)
!
! Call driver code to do the interpolations and/or copies
!
  call driver (ncidi, ncido, ncidt, nvars, ntime)

  if (nf_close (ncido) /= nf_noerr) then
    call err_exit ('error from nf_close')
  end if

  stop
end program

subroutine usage_exit (arg)
  implicit none
  character*(*) arg

  if (arg.ne.' ') write (6,*) arg
  write (6,*) 'Usage: interpic [-s] [-v] ', &
              '-t template infile outfile'
  stop 999
end subroutine

