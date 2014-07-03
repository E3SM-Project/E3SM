program fmain

!-------------------------------------------------------------------------------
! PURPOSE:
! o given a SCRIP map matrix data file, map a field
!
! NOTES:
! o all output data is base on the "_a" grid, the "_b" grid is ignored
!-------------------------------------------------------------------------------

  implicit none
  include 'netcdf.inc'

  integer,parameter :: R8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: R4 = selected_real_kind( 6) ! 4 byte real
  integer,parameter :: RN = kind(1.0)              ! native real
  integer,parameter :: I8 = selected_int_kind (13) ! 8 byte integer
  integer,parameter :: I4 = selected_int_kind ( 6) ! 4 byte integer
  integer,parameter :: IN = kind(1)                ! native integer
  integer,parameter :: CS = 80                     ! short char
  integer,parameter :: CL = 256                    ! long char
  integer,parameter :: CX = 512                    ! extra-long char
  integer :: n                      ! index 
  integer :: nargs                  ! number of arguments  
  integer, external  :: iargc       ! number of arguments function
  character(LEN=512) :: arg         ! input argument
  character(LEN=512) :: cmdline     ! input command line
  character(LEN=512) :: fmap        ! file name ( input nc file)
  character(LEN=512) :: fn_out,fn_in     ! temporary 
  character(LEN=512) :: var_out,var_in     ! temporary
  character(LEN=512) :: usercomment ! user comment 
  character(LEN= 8)  :: cdate       ! wall clock date
  character(LEN=10)  :: ctime       ! wall clock time
  !----------------------------------------------------

  fmap        = 'null'
  fn_in       = 'null'
  var_in      = 'null'
  fn_out      = 'null'
  var_out     = 'null'
  usercomment = 'null'

  nargs = iargc()
  if (nargs == 0) then
     write(6,*)'invoke map_field -h for usage'
     stop
  end if

  cmdline = 'map_field '
  n = 1
  do while (n <= nargs)
    arg = ' '
    call getarg (n, arg)
    n = n + 1

    select case (arg)
    case ('-m')
       ! input mapping file
       call getarg (n, arg)
       n = n + 1
       fmap = trim(arg)
       cmdline = trim(cmdline) // ' -m ' // trim(arg)
    case ('-of')
       ! output filename
       call getarg (n, arg)
       n = n + 1
       fn_out = trim(arg)
       cmdline = trim(cmdline) // ' -of ' // trim(arg)
    case ('-ov')
       ! output variable name
       call getarg (n, arg)
       n = n + 1
       var_out = trim(arg)
       cmdline = trim(cmdline) // ' -ov ' // trim(arg)
    case ('-if')
       ! input filename
       call getarg (n, arg)
       n = n + 1
       fn_in = trim(arg)
       cmdline = trim(cmdline) // ' -if ' // trim(arg)
    case ('-iv')
       ! input variable name
       call getarg (n, arg)
       n = n + 1
       var_in = trim(arg)
       cmdline = trim(cmdline) // ' -iv ' // trim(arg)
    case ('-c')
       ! user comment
       call getarg (n, arg)
       n = n + 1
       usercomment = trim(arg)
    case ('-h') 
       call usage_exit (' ')
    case default
       write (6,*) 'Argument ', arg,' is not known'
       call usage_exit (' ')
    end select
  end do
  
  if (fmap == 'null' .or. fn_out == 'null' .or. fn_in== 'null') then
    call usage_exit ('Must specify all the following arguments')
  end if

  call date_and_time(cdate,ctime)

  call map_field (fmap, fn_in, var_in, fn_out, var_out, usercomment)

!--------------
contains
!--------------

  subroutine map_field(fmap, fn_in, var_in, fn_out, var_out, usercomment)

   implicit none

   !--- includes ---
   include 'netcdf.inc'       ! netCDF defs

   character(LEN=*), intent(in) :: fmap        ! file name ( input nc file)
   character(LEN=*), intent(in) :: fn_in       ! file name
   character(LEN=*), intent(in) :: var_in      ! var name 
   character(LEN=*), intent(in) :: fn_out      ! file name
   character(LEN=*), intent(in) :: var_out     ! var name 
   character(LEN=*), intent(in) :: usercomment ! user comment from namelist

   !--- domain data ---
   integer         ::   nb        ! size of 1d domain
   integer         ::   nbi       ! size of i-axis of 2d domain
   integer         ::   nbj       ! size of j-axis of 2d domain
   integer         ::   na        ! size of source array
   integer         ::   nai       ! size of i-axis of 2d domain
   integer         ::   naj       ! size of j-axis of 2d domain
   integer         ::   nd        ! size of fld1
   real(r8) ,pointer :: xc(  :)   ! x-coords of center    
   real(r8) ,pointer :: yc(  :)   ! y-coords of center   
   real(r8) ,pointer :: xv(:,:)   ! x-coords of verticies 
   real(r8) ,pointer :: yv(:,:)   ! y-coords of verticies 
   real(r8) ,pointer :: area(:)   ! cell area
   real(r8) ,pointer :: fld1(:)   ! fld1
   real(r8) ,pointer :: fld12(:,:)  ! fld1 2d
   real(r8) ,pointer :: fld2(:)   ! fld2
   integer  ,pointer :: src_grid_dims(:)
   integer  ,pointer :: dst_grid_dims(:)


   !--- for mapping ---
   integer          :: ns        ! size of wgts list
   integer ,pointer :: col (  :) ! column index
   integer ,pointer :: row (  :) ! row index
   real(r8),pointer :: S   (  :) ! wgts

    !--- local ---
   character(LEN=CL)     :: flong       ! long name
   character(LEN=CL)     :: str_grido   ! global attribute str - grid_file_ocn
   character(LEN=CL)     :: str_grida   ! global attribute str - grid_file_atm
   integer               :: fid         ! nc file     ID
   integer               :: i,j,k       ! generic indicies
   integer               :: attnum      ! attribute number  
   character(LEN=CL)     :: units       ! netCDF attribute name string
   integer               :: vid         ! nc variable ID
   integer               :: did         ! nc dimension ID
   integer               :: dimids(2)   ! nc dimension ID
   integer               :: rcode       ! routine return error code
   integer               :: dst_grid_rank, src_grid_rank
   real(r8),parameter    :: pi  = 3.14159265358979323846
   real(r8),parameter    :: c0  = 0.00000000000000000000
   real(r8),parameter    :: c1  = 1.00000000000000000000
   real(r8),parameter    :: c90 = 90.0000000000000000000

   !--- formats ---
   !   character(LEN=*),parameter :: F00 = "(120a )"
   !   character(LEN=*),parameter :: F02 = "(a,5i6,i12)"
   !   character(LEN=*),parameter :: F10=&
   !   & "('Data created: 'i4,'-',i2,2('-',i2),' ',i2,2('-',i2),' ')"

   !----------------------------------------------------------------------------

   write(6,*) 'fmap   = ',trim(fmap)
   write(6,*) 'fn_in  = ',trim(fn_in)
   write(6,*) 'var_in = ',trim(var_in)
   write(6,*) 'fn_out = ',trim(fn_out)
   write(6,*) 'var_out= ',trim(var_out)
   write(6,*) 'usercomment= ',trim(usercomment)

   !----------------------------------------------------------------------------
   write(6,*) ' '
   write(6,*) 'input SCRIP data...'
   !----------------------------------------------------------------------------

      write(6,*) ' '
      write(6,*) 'input file  = ',fmap(1:len_trim(fmap))
      call check_ret(nf_open(fmap(1:len_trim(fmap)),NF_NOWRITE,fid))
      write(6,*) 'open ',trim(fmap)

      str_grido = 'unknown'
      str_grida = 'unknown'
      
      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_ocn', str_grido)
      if ( rcode == nf_enotatt ) then
         rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_src', str_grido)
      end if
      
      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_atm', str_grida)
      if ( rcode == nf_enotatt ) then
         rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_dst', str_grida)
      end if

      write(6,*) 'grid_file_ocn= ',trim(str_grido)
      write(6,*) 'grid_file_atm= ',trim(str_grida)
      
      !----------------------------------------------
      ! get domain info 
      !----------------------------------------------

      call check_ret(nf_inq_dimid (fid, 'n_b', did))
      call check_ret(nf_inq_dimlen(fid, did   , nb))
      rcode = nf_inq_dimid (fid, 'ni_b', did)
      if (rcode == 0) then
         call check_ret(nf_inq_dimlen(fid, did, nbi))
      else
         nbi = n
      end if
      rcode = nf_inq_dimid (fid, 'nj_b', did)
      if (rcode == 0) then
         call check_ret(nf_inq_dimlen(fid, did, nbj))
      else
         nbj = 1
      end if
      if (nbi == 1 .and. nbj == 0) then
         nbi = n
         nbj = 1
      end if
      
      call check_ret(nf_inq_dimid (fid, 'n_s', did))
      call check_ret(nf_inq_dimlen(fid, did   , ns))
      call check_ret(nf_inq_dimid (fid, 'n_a', did))
      call check_ret(nf_inq_dimlen(fid, did   , na))
      call check_ret(nf_inq_dimid (fid, 'dst_grid_rank', did))
      call check_ret(nf_inq_dimlen(fid, did, dst_grid_rank))
      call check_ret(nf_inq_dimid (fid, 'src_grid_rank', did))
      call check_ret(nf_inq_dimlen(fid, did, src_grid_rank)) 

      allocate(src_grid_dims(src_grid_rank), dst_grid_dims(dst_grid_rank))
      call check_ret(nf_inq_varid(fid,'src_grid_dims', vid))
      call check_ret(nf_get_var_int(fid,vid,src_grid_dims)) 
      call check_ret(nf_inq_varid(fid,'dst_grid_dims', vid))
      call check_ret(nf_get_var_int(fid,vid,dst_grid_dims)) 

      if (dst_grid_rank == 2) then
         nbi = dst_grid_dims(1)
         nbj = dst_grid_dims(2)
      end if
      nb = nbi * nbj

      if (src_grid_rank == 2) then
         nai = src_grid_dims(1)
         naj = src_grid_dims(2)
         na = nai * naj
      end if

      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_atm', str_grida)
      if ( rcode == nf_enotatt ) then
         call check_ret(nf_get_att_text(fid, NF_GLOBAL, 'grid_file_dst', str_grida))
      end if
      
      write(6,*)'na,nai,naj,nb,nbi,nbj,ns= ',na,nai,naj,nb,nbi,nbj,ns
      
      allocate(fld1(na)) 
      allocate(fld2(nb))

      allocate(col(ns))
      allocate(row(ns))
      allocate(S(ns))

      !--- read weights ---

      call check_ret(nf_inq_varid(fid,'col', vid ))   
      call check_ret(nf_get_var_int(fid,vid, col ))   
      call check_ret(nf_inq_varid(fid,'row', vid ))   
      call check_ret(nf_get_var_int(fid,vid, row ))   
      call check_ret(nf_inq_varid(fid,'S', vid ))     
      call check_ret(nf_get_var_double(fid,vid,S))    

      call check_ret(nf_close(fid))

      !--- read fld1 ---

      write(6,*) ' '
      write(6,*) 'input file  = ',fn_in(1:len_trim(fn_in))
      call check_ret(nf_open(fn_in(1:len_trim(fn_in)),NF_NOWRITE,fid))
      write(6,*) 'open ',trim(fn_in)
      call check_ret(nf_inq_varid(fid,trim(var_in), vid ))
      call check_ret(nf_inq_varndims(fid,vid,nd))
      if (nd == 1) then       
        call check_ret(nf_inq_vardimid(fid,vid,dimids))
        call check_ret(nf_inq_dimlen(fid,dimids(1),nai))
        if (nai /= na) then
           write(6,*) 'error nai size ',nai,na
           stop
        endif
        call check_ret(nf_get_var_double(fid,vid,fld1 ))  
      elseif (nd == 2) then
        call check_ret(nf_inq_vardimid(fid,vid,dimids))
        call check_ret(nf_inq_dimlen(fid,dimids(1),nai))
        call check_ret(nf_inq_dimlen(fid,dimids(2),naj))
        if (nai*naj /= na) then
           write(6,*) 'error nai naj size ',nai,naj,nai*naj,na
           stop
        endif
        allocate(fld12(nai,naj))
        call check_ret(nf_get_var_double(fid,vid,fld12 ))  
        n = 0
        do j = 1,naj
        do i = 1,nai
           n = n + 1
           fld1(n) = fld12(i,j)
        enddo
        enddo
        deallocate(fld12)
      else
        write(6,*) 'error nd ',nd
        stop
      endif
      call check_ret(nf_close(fid))

      !----------------------------------------------------------------------------
      write(6,*) 'compute fld2'
      !----------------------------------------------------------------------------

      fld2 = c0
      do k = 1,ns
         fld2(row(k)) = fld2(row(k)) + fld1(col(k))*S(k)
      enddo


      !-----------------------------------------------------------------
      ! create a new nc files
      !-----------------------------------------------------------------
      
      write(6,*) ' '
      write(6,*) 'output output file...'

      units = 'unitless'
      flong = 'mapped '//trim(var_in)

      call check_ret(nf_create(fn_out(1:len_trim(fn_out)),NF_CLOBBER,fid))
      write(6,*) 'write ',trim(fn_out)
      call write_file(fid, fn_out, units, nb, nbi, nbj, &
              fld2, flong, fmap, str_grido, str_grida)
      call check_ret(nf_close(fid))
      write(6,*) 'successfully created domain file ', trim(fn_out)

  end subroutine map_field

!===========================================================================

  subroutine check_ret(ret)
    ! Check return status from netcdf call
    implicit none
    integer, intent(in) :: ret
    
    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error with rcode = ', ret,' error = ', nf_strerror(ret)
       call abort()
    end if
    
  end subroutine check_ret

  subroutine usage_exit (arg)
    implicit none
    character(len=*) :: arg
    if (arg /= ' ') then
       write (6,*) arg
    end if
    write(6,*) ' Purpose:'
    write(6,*) '    Map data from one grid to another given a mapping file '
    write(6,*) '    and a source field in another netcdf file '
    write(6,*) ' '
    write(6,*) ' Usage: '
    write(6,*) '    map_field  -m <filemap>'
    write(6,*) '                -if <input_file>'
    write(6,*) '                -iv <input_varname>'
    write(6,*) '                -of <output_file>'
    write(6,*) '                -ov <output_varname>'
    write(6,*) '                [-c <usercomment>]'
    write(6,*) ' '
    write(6,*) ' Where: '
    write(6,*) '    filemap = input mapping filename'
    write(6,*) '    input_file = input file where input field is read'
    write(6,*) '    input_varname = name of the variable in the input_file to read'
    write(6,*) '    output_file = where mapped field is written'
    write(6,*) '    input_varname = name of the mapped field variable in the output_file'
    write(6,*) '    usercomment = optional, netcdf global attribute (character string)'
    write(6,*) ' '
    write(6,*) '  NOTE that the output_file is always clobbered when using this tool'
    write(6,*) ' '
    stop 
  end subroutine usage_exit

!===========================================================================

  subroutine write_file(fid, fout, units, n, ni, nj, &
              fld2, fld2long, fmap, str_grido, str_grida)
       
    implicit none
    
    !--- includes ---
    include 'netcdf.inc'       ! netCDF defs
    
    integer         , intent(in) :: fid          ! nc file  ID
    character(LEN=*), intent(in) :: fout       ! file name ( output nc file)
    character(LEN=*), intent(in) :: units        ! netCDF attribute name string
    integer         , intent(in) :: n            ! size of 1d domain
    integer         , intent(in) :: ni           ! size of i-axis of 2d domain
    integer         , intent(in) :: nj           ! size of j-axis of 2d domain
    real(r8)        , pointer    :: fld2(:)      ! output field
    character(LEN=*), intent(in) :: fld2long     ! global attribute str - grid_file_ocn
    character(LEN=*), intent(in) :: fmap         ! global attribute str - grid_file_ocn
    character(LEN=*), intent(in) :: str_grido    ! global attribute str - grid_file_ocn
    character(LEN=*), intent(in) :: str_grida    ! global attribute str - grid_file_atm
    
    !--- local ---
    character(LEN=CL)     :: host        ! hostname of machine running on
    character(LEN=CL)     :: str         ! fixed    length char string
    character(LEN=CL)     :: str_title   ! global attribute str - title
    character(LEN=CL)     :: str_source  ! global attribute str - source
    character(LEN=CL)     :: user        ! user name
    integer               :: strlen      ! (trimmed) length of string
    integer               :: vid         ! nc variable ID
    integer               :: did         ! nc dimension ID
    integer               :: vdid(3)     ! vector of nc dimension ID
    integer               :: rcode       ! routine return error code
    
    character(*),parameter:: version = 'SVN $Id: map_field.F90 39483 2012-08-16 17:35:52Z mlevy@ucar.edu $'

    ! global attributes
    str   = 'map_field data: '
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str))
    
    str   = 'CF-1.0'
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'Conventions',len_trim(str),str))
    
    str = trim(version)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'source_code',len_trim(str),str))
    
    str = ' $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/map_field/trunk/src/map_field.F90 $'
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'SVN_url',len_trim(str),str))
    
#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif
    
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'Compiler_Optimized',len_trim(str),str))
   
    call sys_getenv('HOST',host,rcode)
    if (rcode == 0) then
       call check_ret(nf_put_att_text(fid,NF_GLOBAL,'hostname' ,len_trim(host),host))
    else
       call sys_getenv('HOSTNAME',host,rcode)
       if (rcode == 0) then
          call check_ret(nf_put_att_text(fid,NF_GLOBAL,'hostname' ,len_trim(host),host))
       else
          write(6,*) 'WARNING: could not determine hostname, so that information will not be stored in netCDF attribute. To avoid this warning in the future, set environment variable HOST or HOSTNAME.'
       end if
    end if

    call sys_getenv('LOGNAME',user,rcode)
    if (rcode /= 0) then
       write(6,*) ' ERROR: getting LOGNAME'
       stop
    end if
    str = 'created by '//trim(user)//', '//cdate(1:4)//'-'//cdate(5:6)//'-'//cdate(7:8) &
         &                //' '//ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'history' ,len_trim(str),str))
    
    str = trim(fout)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'source' ,len_trim(str),str))
    str = trim(fmap)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'map_file',len_trim(str),str))
    str = trim(str_grido)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'map_grid_file_ocn',len_trim(str),str))
    str = trim(str_grida)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'map_grid_file_atm',len_trim(str),str))
    str = trim(usercomment)
    if ( str(1:4) /= 'null' ) then
       call check_ret(nf_put_att_text(fid,NF_GLOBAL,'user_comment',len_trim(str),str))
    end if

    !-----------------------------------------------------------------
    ! dimension data
    !-----------------------------------------------------------------
    call check_ret(nf_def_dim(fid, 'n' , n , did)) ! # of points total
    call check_ret(nf_def_dim(fid, 'ni', ni, did)) ! # of points wrt i
    vdid(1) = did
    call check_ret(nf_def_dim(fid, 'nj', nj, did)) ! # of points wrt j
    vdid(2) = did

    !-----------------------------------------------------------------
    ! define data -- coordinates, input grid
    !-----------------------------------------------------------------

    call check_ret(nf_def_var  (fid,trim(var_out),NF_DOUBLE ,2,vdid,vid))
    str   = trim(fld2long)
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = trim(units)
    call check_ret(nf_put_att_text(fid,vid,"units"      ,len_trim(str),str))

    call check_ret(nf_enddef(fid))

    call check_ret(nf_inq_varid(fid,trim(var_out),vid)) 
    call check_ret(nf_put_var_double(fid,  vid ,fld2)) 

  end subroutine write_file

SUBROUTINE sys_getenv(name, val, rcode)

   IMPLICIT none

   !----- arguments -----
   character(*)        ,intent(in)  :: name    ! env var name
   character(*)        ,intent(out) :: val     ! env var value
   integer,intent(out) :: rcode   ! return code

   !----- local -----
   integer             :: lenname ! length of env var name
   integer             :: lenval  ! length of env var value
   character(len=80)           :: tmpval  ! temporary env var value

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_getenv) '
   character(*),parameter :: F00     = "('(shr_sys_getenv) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

   lenname=len_trim(name)

#if (defined IRIX64 || defined CRAY || defined UNICOSMP)

   call pxfgetenv(name, lenname, val, lenval, rcode)

#elif (defined AIX || defined OSF1 || defined SUNOS || defined LINUX || defined NEC_SX)
#ifdef F2003
    call GET_ENVIRONMENT_VARIABLE(trim(name),value=val, length=lenval, status=rcode)
#else
   call getenv(trim(name),tmpval)
   val=trim(tmpval)
   rcode = 0
   if (len_trim(val) ==  0         ) rcode = 1
   if (len_trim(val) >  CL) rcode = 2
#endif
#else

   write(*,F00) 'ERROR: no implementation of getenv for this architecture'
   stop subname//'no implementation of getenv for this machine'

#endif

END SUBROUTINE sys_getenv


end program fmain
