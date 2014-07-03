program fmain

!-------------------------------------------------------------------------------
! PURPOSE:
! o given a SCRIP map matrix data file, create datm/dlnd/docn/dice domain data files
!
! NOTES:
! o all output data is base on the "_a" grid, the "_b" grid is ignored
! o to compile on an NCAR's SGI, tempest (Dec 2004): 
!   unix> f90 -64 -mips4 -r8 -i4 -lfpe -I/usr/local/include Make_domain.F90 \
!         -L/usr/local/lib64/r4i4 -lnetcdf
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
  integer :: set_fv_pole_yc         ! fix pole ycs on this grid [0,1,2]
  integer :: nargs                  ! number of arguments  
  integer, external  :: iargc       ! number of arguments function
  character(LEN=512) :: arg         ! input argument
  character(LEN=512) :: cmdline     ! input command line
  character(LEN=512) :: fmap        ! file name ( input nc file)
  character(LEN=512) :: fn1_out     ! temporary 
  character(LEN=512) :: fn2_out     ! temporary
  character(LEN=512) :: fn1_out_ocn ! file name (output nc file) for grid _a 
  character(LEN=512) :: fn2_out_lnd ! file name (output nc file) for grid _b (lnd fraction)
  character(LEN=512) :: fn2_out_ocn ! file name (output nc file) for grid _b (ocn fraction)
  character(LEN=512) :: usercomment ! user comment 
  character(LEN= 8)  :: cdate       ! wall clock date
  character(LEN=10)  :: ctime       ! wall clock time
  !----------------------------------------------------

  set_fv_pole_yc = 0
  fmap        = 'null'
  fn1_out     = 'null'
  fn2_out     = 'null'
  usercomment = 'null'

  nargs = iargc()
  if (nargs == 0) then
     write(6,*)'invoke gen_domain -h for usage'
     stop
  end if

  cmdline = 'gen_domain '
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
    case ('-o')
       ! output ocean grid name
       call getarg (n, arg)
       n = n + 1
       fn1_out = trim(arg)
       cmdline = trim(cmdline) // ' -o ' // trim(arg)
    case ('-l')
       ! output land grid name
       call getarg (n, arg)
       n = n + 1
       fn2_out = trim(arg)
       cmdline = trim(cmdline) // ' -l ' // trim(arg)
    case ('-p')
       ! set pole on this grid [0,1,2]
       call getarg (n, arg)
       n = n + 1
       set_fv_pole_yc = ichar(trim(arg))
       write(6,*)'set_fv_pole_yc is ',set_fv_pole_yc
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
  
  if (fmap == 'null' .or. fn1_out == 'null' .or. fn2_out== 'null') then
    call usage_exit ('Must specify all the following arguments')
  end if

  call date_and_time(cdate,ctime)

  fn1_out_ocn = 'domain.ocn.' // trim(fn1_out)                         // '.' // cdate(3:8) // '.nc'
  fn2_out_lnd = 'domain.lnd.' // trim(fn2_out) // '_' // trim(fn1_out) // '.' // cdate(3:8) // '.nc'
  fn2_out_ocn = 'domain.ocn.' // trim(fn2_out) // '_' // trim(fn1_out) // '.' // cdate(3:8) // '.nc'

  call gen_domain (fmap, fn1_out_ocn, fn2_out_lnd, fn2_out_ocn, set_fv_pole_yc, usercomment)

contains

  subroutine gen_domain(fmap, fn1_out_ocn, fn2_out_lnd, fn2_out_ocn, set_fv_pole_yc, usercomment)

   implicit none

   !--- includes ---
   include 'netcdf.inc'       ! netCDF defs

   character(LEN=*), intent(in) :: fmap        ! file name ( input nc file)
   character(LEN=*), intent(in) :: fn1_out_ocn ! file name (output nc file) for grid _a
   character(LEN=*), intent(in) :: fn2_out_lnd ! file name (output nc file) for grid _b (lnd frac)
   character(LEN=*), intent(in) :: fn2_out_ocn ! file name (output nc file) for grid _b (ocn frac)
   integer         , intent(in) :: set_fv_pole_yc
   character(LEN=*), intent(in) :: usercomment ! user comment from namelist

   !--- domain data ---
   integer         ::   n         ! size of 1d domain
   integer         ::   ni        ! size of i-axis of 2d domain
   integer         ::   nj        ! size of j-axis of 2d domain
   integer         ::   nv = 4    ! assume retalinear grid
   real(r8) ,pointer :: xc(  :)   ! x-coords of center    
   real(r8) ,pointer :: yc(  :)   ! y-coords of center   
   real(r8) ,pointer :: xv(:,:)   ! x-coords of verticies 
   real(r8) ,pointer :: yv(:,:)   ! y-coords of verticies 
   real(r8) ,pointer :: area(:)   ! cell area
   integer  ,pointer :: lmask(:)  ! domain mask           
   real(r8) ,pointer :: lfrac(:)  ! cell fraction
   integer  ,pointer :: omask(:)  ! domain mask           
   real(r8) ,pointer :: ofrac(:)  ! cell fraction
   integer  ,pointer :: src_grid_dims(:)
   integer  ,pointer :: dst_grid_dims(:)


   !--- for mapping ---
   logical          :: complf    ! flag for computing landfrac
   integer          :: ns        ! size of wgts list
   integer ,pointer :: col (  :) ! column index
   integer ,pointer :: row (  :) ! row index
   real(r8),pointer :: S   (  :) ! wgts
   integer          :: na        ! size of source array
   integer ,pointer :: mask_a(:) ! mask of source array, integer
   real(r8),pointer :: frac_a(:) ! mask of source array, real
   real(r8)         :: eps       ! allowable frac error
   real(r8)         :: lfrac_min ! min frac value before being set to fminval
   real(r8)         :: lfrac_max ! max frac value before being set to fmaxval
   real(r8)         :: fminval   ! set frac to zero if frac < fminval
   real(r8)         :: fmaxval   ! set frac to one  if frac > fmaxval

    !--- local ---
   character(LEN=CL)     :: str_da      ! global attribute str - domain_a
   character(LEN=CL)     :: str_db      ! global attribute str - domain_b
   character(LEN=CL)     :: str_grido   ! global attribute str - grid_file_ocn
   character(LEN=CL)     :: str_grida   ! global attribute str - grid_file_atm
   integer               :: fid         ! nc file     ID
   integer               :: i,j,k       ! generic indicies
   integer               :: attnum      ! attribute number  
   character(LEN=CL)     :: fn_out      ! current output file name
   character(LEN=CL)     :: fn_out_lnd  ! current output file name
   character(LEN=CL)     :: fn_out_ocn  ! current output file name
   integer               :: nf          ! number of files counter
   character(LEN=CL)     :: units_xc, units_yc  ! netCDF attribute name string
   character(LEN=2)      :: suffix      ! suffix _a or _b sets input grid
   logical               :: pole_fix    ! fix pole ycs
   integer               :: vid         ! nc variable ID
   integer               :: did         ! nc dimension ID
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

   eps = 1.0e-12
   fminval = 0.001
   fmaxval = c1
   write(6,*) 'fmap   = ',trim(fmap)
   write(6,*) 'fn1_out_ocn= ',trim(fn1_out_ocn)
   write(6,*) 'fn2_out_lnd= ',trim(fn2_out_lnd)
   write(6,*) 'fn2_out_ocn= ',trim(fn2_out_ocn)
   write(6,*) 'usercomment= ',trim(usercomment)
   write(6,*) 'eps    = ',eps
   write(6,*) 'fminval= ',fminval
   write(6,*) 'fmaxval= ',fmaxval
   write(6,*) 'set_fv_pole_yc = ',set_fv_pole_yc

   !----------------------------------------------------------------------------
   write(6,*) ' '
   write(6,*) 'input SCRIP data...'
   !----------------------------------------------------------------------------

   do nf = 1,2

      if (nf == 1) then
         suffix = '_a'
         fn_out = fn1_out_ocn
         complf = .false.
      elseif (nf == 2) then
         suffix = '_b'
         fn_out_lnd = fn2_out_lnd
         fn_out_ocn = fn2_out_ocn
         complf = .true.
         pole_fix = .false.
      else
         write(6,*) ' ERROR: nf loop error '
         stop
      endif

      pole_fix = .false.
      if (nf == set_fv_pole_yc) pole_fix = .true.
      write(6,*)'pole_fix = ',pole_fix
      
      write(6,*) ' '
      write(6,*) 'input file  = ',fmap(1:len_trim(fmap))
      call check_ret(nf_open(fmap(1:len_trim(fmap)),NF_NOWRITE,fid))
      write(6,*) 'open ',trim(fmap)

      str_da = 'unknown'
      str_db = 'unknown'
      str_grido = 'unknown'
      str_grida = 'unknown'
      
      call check_ret(nf_get_att_text(fid, NF_GLOBAL, 'domain_a', str_da))
      call check_ret(nf_get_att_text(fid, NF_GLOBAL, 'domain_b', str_db))
      
      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_ocn', str_grido)
      if ( rcode == nf_enotatt ) then
         rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_src', str_grido)
      end if
      
      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_atm', str_grida)
      if ( rcode == nf_enotatt ) then
         rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_dst', str_grida)
      end if

      write(6,*) 'domain_a     = ',trim(str_da)
      write(6,*) 'domain_b     = ',trim(str_db)
      write(6,*) 'grid_file_ocn= ',trim(str_grido)
      write(6,*) 'grid_file_atm= ',trim(str_grida)
      
      !----------------------------------------------
      ! get domain info 
      !----------------------------------------------
      if (trim(suffix) == '_b') then
         call check_ret(nf_inq_varid(fid,'dst_grid_dims', vid)) 
         call check_ret(nf_get_var_int(fid,vid,dst_grid_dims ))
      end if

      call check_ret(nf_inq_dimid (fid, 'n'//trim(suffix) , did))
      call check_ret(nf_inq_dimlen(fid, did   , n))
      call check_ret(nf_inq_dimid (fid, 'nv'//trim(suffix), did))
      call check_ret(nf_inq_dimlen(fid, did   , nv))
      rcode = nf_inq_dimid (fid, 'ni'//trim(suffix), did)
      if (rcode == 0) then
         call check_ret(nf_inq_dimlen(fid, did, ni))
      else
         ni = n
      end if
      rcode = nf_inq_dimid (fid, 'nj'//trim(suffix), did)
      if (rcode == 0) then
         call check_ret(nf_inq_dimlen(fid, did, nj))
      else
         nj = 1
      end if
      if (ni == 1 .and. nj == 0) then
         ni = n
         nj = 1
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

      rcode = nf_get_att_text(fid, NF_GLOBAL, 'grid_file_atm', str_grida)
      if ( rcode == nf_enotatt ) then
         call check_ret(nf_get_att_text(fid, NF_GLOBAL, 'grid_file_dst', str_grida))
      end if

      
      write(6,*)'n,nv,ni,nj,na,ns= ',n,nv,ni,nj,na,ns
      
      allocate(xc(n))     ! x-coordinates of center for either _a or _b grid
      allocate(yc(n))     ! y-coordinates of center for either _a or _b grid
      allocate(xv(nv,n))  ! x-coordinates of verticies for either _a or _b grid
      allocate(yv(nv,n))  ! y-coordinates of verticies for either _a or _b grid
      allocate(area(n))   ! grid cell area

      call check_ret(nf_inq_varid(fid,'xc'//trim(suffix), vid))
      call check_ret(nf_get_att_text(fid, vid, 'units', units_xc))
      call check_ret(nf_get_var_double(fid,vid,  xc ))
      call check_ret(nf_inq_varid(fid,'yc'//trim(suffix), vid))
      call check_ret(nf_get_att_text(fid, vid, 'units', units_yc))
      call check_ret(nf_get_var_double(fid,vid,  yc ))
      call check_ret(nf_inq_varid(fid,'xv'//trim(suffix), vid))
      call check_ret(nf_get_var_double(fid,vid,  xv ))
      call check_ret(nf_inq_varid(fid,'yv'//trim(suffix), vid))
      call check_ret(nf_get_var_double(fid,vid,  yv ))
      call check_ret(nf_inq_varid(fid,'area'//trim(suffix), vid ))
      call check_ret(nf_get_var_double(fid,vid,area ))
      
      !--- set default ocean frac ---

      allocate(omask(n))  ! domain mask ocean
      allocate(ofrac(n))  ! area frac of mask "_a" on grid "_b" or float(mask)

      if (.not. complf) then

         ! Determine ocn mask on ocn grid
         call check_ret(nf_inq_varid(fid,'mask'//trim(suffix), vid ))
         call check_ret(nf_get_var_int   (fid,vid,omask ))
         ofrac(:) = c0
         where (omask /= 0) ofrac = c1

      else
         !----------------------------------------------------------------------------
         write(6,*) 'compute frac'
         !----------------------------------------------------------------------------
         allocate(col(ns))
         allocate(row(ns))
         allocate(S(ns))
         allocate(mask_a(na))
         allocate(frac_a(na))

         allocate(lmask(n))  ! domain mask land
         allocate(lfrac(n))  ! area frac of mask "_a" on grid "_b" or float(mask)

         call check_ret(nf_inq_varid(fid,'col', vid ))   
         call check_ret(nf_get_var_int(fid,vid, col ))   
         call check_ret(nf_inq_varid(fid,'row', vid ))   
         call check_ret(nf_get_var_int(fid,vid, row ))   
         call check_ret(nf_inq_varid(fid,'S', vid ))     
         call check_ret(nf_get_var_double(fid,vid,S))    
         call check_ret(nf_inq_varid(fid,'mask_a', vid ))
         call check_ret(nf_get_var_int(fid,vid,mask_a ))  
         frac_a = c0
         where (mask_a /= 0) frac_a = c1
         !--- compute ocean fraction on atm grid ---
         ofrac = c0
         do k = 1,ns
            ofrac(row(k)) = ofrac(row(k)) + frac_a(col(k))*S(k)
         enddo
         !--- convert to land fraction, 1.0-frac and ---
         !--- trap errors and modify computed frac ---
         lmask(:) = 0
         omask(:) = 1
         lfrac_min = fmaxval
         lfrac_max = fminval
         do k = 1,n
            lfrac(k) = c1 - ofrac(k)
            lfrac_min = min(lfrac_min,lfrac(k))
            lfrac_max = max(lfrac_max,lfrac(k))
            if (lfrac(k) > fmaxval) lfrac(k) = c1
            if (lfrac(k) < fminval) lfrac(k) = c0   ! extra requirement for landfrac
            ofrac(k) = c1 - lfrac(k)
            if (lfrac(k) /= c0) then
               lmask(k) = 1
            end if
            if (ofrac(k) == c0) then
               omask(k) = 0
            end if
         enddo
         write(6,*) '----------------------------------------------------------------------'
         write(6,*) 'IMPORTANT: note original min/max frac and decide if that''s acceptable'
         write(6,*) 'original lfrac clipped above by       : ',fmaxval
         write(6,*) 'original reset to zero when less than : ',fminval
         write(6,*) 'original min, max lfrac : ',lfrac_min,lfrac_max
         write(6,*) 'final min, max llfrac   : ',minval(lfrac),maxval(lfrac)
         write(6,*) '----------------------------------------------------------------------'
      endif
      
      call check_ret(nf_close(fid))

      !-----------------------------------------------------------------
      ! adjust j = 1 and j = nj lats to -+ 90 degrees
      !-----------------------------------------------------------------
      
      if (pole_fix) then
         write(6,*)'ni,nj= ',ni,nj
         if (ni > 1 .and. nj == 1) then
            if (dst_grid_rank /= 2) then 
               write(6,*)'pole_fix not appropriate for unstructured grid'
               stop
            end if
         end if
         do i = 1,dst_grid_dims(1)
            yc(i)      = -c90
            yc(n-dst_grid_dims(1)+i) =  c90
         enddo
      endif
      
      !-----------------------------------------------------------------
      ! create a new nc files
      !-----------------------------------------------------------------
      
      write(6,*) ' '
      write(6,*) 'output domain data...'

      if (n /= ni*nj) then
         STOP 'n'
      end if
      write(6,*) 'nf = ', nf
      if (nf == 1) then
         if (src_grid_rank == 2) then
            ni = src_grid_dims(1)
            nj = src_grid_dims(2)
         end if
         write(6,*) 'create ',trim(fn_out)
         call check_ret(nf_create(fn_out(1:len_trim(fn_out)),NF_CLOBBER,fid))
         write(6,*) 'write ',trim(fn_out)
         call write_file(fid, fmap, units_xc, units_yc, n, ni, nj, &
              xc, yc, xv, yv, area, omask, ofrac, suffix, eps, pole_fix, &
              fmaxval, fminval, str_da, str_db, str_grido, str_grida)
         call check_ret(nf_close(fid))
         write(6,*) 'successfully created domain file ', trim(fn_out)
      else if (nf == 2) then
         if (dst_grid_rank == 2) then
            ni = dst_grid_dims(1)
            nj = dst_grid_dims(2)
         end if
         call check_ret(nf_create(fn_out_lnd(1:len_trim(fn_out_lnd)),NF_CLOBBER,fid))
         write(6,*) 'write ',trim(fn_out_lnd)
         call write_file(fid, fmap, units_xc, units_yc, n, ni, nj, &
              xc, yc, xv, yv, area, lmask, lfrac, suffix, eps, pole_fix, &
              fmaxval, fminval, str_da, str_db, str_grido, str_grida)
         call check_ret(nf_close(fid))
         write(6,*) 'successfully created domain file ', trim(fn_out_lnd)

         call check_ret(nf_create(fn_out_ocn(1:len_trim(fn_out_ocn)),NF_CLOBBER,fid))
         write(6,*) 'write ',trim(fn_out_ocn)
         call write_file(fid, fmap, units_xc, units_yc, n, ni, nj, &
              xc, yc, xv, yv, area, omask, ofrac, suffix, eps, pole_fix, &
              fmaxval, fminval, str_da, str_db, str_grido, str_grida)
         call check_ret(nf_close(fid))
         write(6,*) 'successfully created domain file ', trim(fn_out_ocn)
      end if
      
   enddo

  end subroutine gen_domain

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
    write(6,*) '    Given a SCRIP map matrix data file from the ocean grid '
    write(6,*) '    (where the mask is defined) to the land grid, gen_domain '
    write(6,*) '    creates land and ocean domain files'
    write(6,*) '    These files are currently used by '
    write(6,*) '       datm, dlnd, dice, docn, clm, cice(prescribed mode)'
    write(6,*) ' '
    write(6,*) ' Usage: '
    write(6,*) '    gen_domain  -m <filemap>'
    write(6,*) '                -o <gridocn>'
    write(6,*) '                -l <gridlnd>'
    write(6,*) '                [-p set_fv_pole_yc]'
    write(6,*) '                [-c <usercomment>]'
    write(6,*) ' '
    write(6,*) ' Where: '
    write(6,*) '    filemap = input conservative mapping file name (from ocn->atm)'
    write(6,*) '    gridocn = output ocean grid name'
    write(6,*) '    gridlnd = output land  grid name'
    write(6,*) '    set_fv_pole_yc = [0,1,2] ~ optional, default = 0'
    write(6,*) '    usercomment = optional, netcdf global attribute (character string)'
    write(6,*) ' '
    write(6,*) ' The following output domain files are created:'
    write(6,*) '    domain.lnd.gridlnd_gridocn.nc'
    write(6,*) '      land domain file on the land grid with a '
    write(6,*) '      land fraction corresponding to '
    write(6,*) '      (1-gridocn) mask mapped to the land grid'
    write(6,*) '    domain.ocn.gridlnd_gridocn.nc'
    write(6,*) '      ocean domain on the land grid with an '
    write(6,*) '      ocean fraction corresponding to the '
    write(6,*) '      gridocn mask mapped to the land grid'
    write(6,*) '      this is used when both atm,lnd,ice,ocn are all on the'
    write(6,*) '      same grid (F compset)'
    write(6,*) '    domain.ocn.gridocn.nc'
    write(6,*) '      ocean domain on the ocean grid '
    write(6,*) ' '
    stop 
  end subroutine usage_exit

!===========================================================================

  subroutine write_file(fid, fmap, units_xc, units_yc, n, ni, nj, &
       xc, yc, xv, yv, area, mask, frac, suffix, eps, pole_fix, &
       fmaxval, fminval, str_da, str_db, str_grido, str_grida)
       
    implicit none
    
    !--- includes ---
    include 'netcdf.inc'       ! netCDF defs
    
    integer         , intent(in) :: fid          ! nc file  ID
    character(LEN=*), intent(in) :: fmap         ! file name ( input nc file)
    character(LEN=*), intent(in) :: units_xc, units_yc  ! netCDF attribute name string
    integer         , intent(in) :: n            ! size of 1d domain
    integer         , intent(in) :: ni           ! size of i-axis of 2d domain
    integer         , intent(in) :: nj           ! size of j-axis of 2d domain
    real(r8)        , pointer    :: xc(:)        ! x-coords of center    
    real(r8)        , pointer    :: yc(:)        ! y-coords of center   
    real(r8)        , pointer    :: xv(:,:)      ! x-coords of verticies 
    real(r8)        , pointer    :: yv(:,:)      ! y-coords of verticies 
    real(r8)        , pointer    :: area(:)      ! cell area
    integer         , pointer    :: mask(:)      ! domain mask           
    real(r8)        , pointer    :: frac(:)      ! cell fraction
    character(LEN=*), intent(in) :: suffix       ! suffix _a or _b sets input grid
    real(r8)        , intent(in) :: eps          ! allowable frac error
    logical         , intent(in) :: pole_fix     ! fix pole ycs
    real(r8)        , intent(in) :: fminval      ! set frac to zero if frac < fminval
    real(r8)        , intent(in) :: fmaxval      ! set frac to one  if frac > fmaxval
    character(LEN=*), intent(in) :: str_da       ! global attribute str - domain_a
    character(LEN=*), intent(in) :: str_db       ! global attribute str - domain_b
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
    
    real(r8),parameter    :: pi  = 3.14159265358979323846
    character(*),parameter:: version = 'SVN $Id: gen_domain.F90 41914 2012-11-13 21:58:37Z mlevy@ucar.edu $'

    ! global attributes
    str   = 'CESM domain data: '
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str))
    
    str   = 'CF-1.0'
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'Conventions',len_trim(str),str))
    
    str = trim(version)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'source_code',len_trim(str),str))
    
    str = ' $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/gen_domain/trunk/src/gen_domain.F90 $'
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
    
    str = trim(fmap)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'source' ,len_trim(str),str))
    str = trim(str_da)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'map_domain_a',len_trim(str),str))
    str = trim(str_db)
    call check_ret(nf_put_att_text(fid,NF_GLOBAL,'map_domain_b',len_trim(str),str))
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
    call check_ret(nf_def_dim(fid, 'nj', nj, did)) ! # of points wrt j
    call check_ret(nf_def_dim(fid, 'nv',  4, did)) ! # of verticies per cell

    !-----------------------------------------------------------------
    ! define data -- coordinates, input grid
    !-----------------------------------------------------------------

    call check_ret(nf_inq_dimid(fid,'n' , did   ))
    call check_ret(nf_inq_dimid(fid,'ni',vdid(1)))
    call check_ret(nf_inq_dimid(fid,'nj',vdid(2)))

    call check_ret(nf_def_var  (fid,'xc',NF_DOUBLE,2,vdid,vid))
    str   = 'longitude of grid cell center'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'degrees_east'
    call check_ret(nf_put_att_text(fid,vid,"units"    ,len_trim(str),str))
    str   = 'xv'
    call check_ret(nf_put_att_text(fid,vid,"bounds"   ,len_trim(str),str))

    call check_ret(nf_def_var  (fid,'yc',NF_DOUBLE,2,vdid,vid))
    str   = 'latitude of grid cell center'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'degrees_north'
    call check_ret(nf_put_att_text(fid,vid,"units"    ,len_trim(str),str))
    str   = 'yv'
    call check_ret(nf_put_att_text(fid,vid,"bounds"   ,len_trim(str),str))
    if (pole_fix) then
       write(str,*) 'set_fv_pole_yc ON, yc = -+90 at j=1,j=nj'
       call check_ret(nf_put_att_text(fid,vid,'filter1' ,len_trim(str),str))
    endif

    call check_ret(nf_inq_dimid(fid,'nv',vdid(1)))
    call check_ret(nf_inq_dimid(fid,'ni',vdid(2)))
    call check_ret(nf_inq_dimid(fid,'nj',vdid(3)))

    call check_ret(nf_def_var  (fid,'xv',NF_DOUBLE,3,vdid,vid))
    str   = 'longitude of grid cell verticies'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'degrees_east'
    call check_ret(nf_put_att_text(fid,vid,"units"    ,len_trim(str),str))

    call check_ret(nf_def_var  (fid,'yv',NF_DOUBLE,3,vdid,vid))
    str   = 'latitude of grid cell verticies'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'degrees_north'
    call check_ret(nf_put_att_text(fid,vid,"units"    ,len_trim(str),str))

    call check_ret(nf_inq_dimid(fid,'ni',vdid(1)))
    call check_ret(nf_inq_dimid(fid,'nj',vdid(2)))

    call check_ret(nf_def_var  (fid,'mask',NF_INT ,2,vdid,vid))
    str   = 'domain mask'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'unitless'
    call check_ret(nf_put_att_text(fid,vid,"note"    ,len_trim(str),str))
    str   = 'xc yc'
    call check_ret(nf_put_att_text(fid,vid,"coordinates",len_trim(str),str))
    str   = '0 value indicates cell is not active'
    call check_ret(nf_put_att_text(fid,vid,"comment",len_trim(str),str))

    call check_ret(nf_def_var  (fid,'area',NF_DOUBLE,2,vdid,vid))
    str   = 'area of grid cell in radians squared'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'xc yc'
    call check_ret(nf_put_att_text(fid,vid,"coordinates",len_trim(str),str))
    str   = 'radian2'
    call check_ret(nf_put_att_text(fid,vid,"units"    ,len_trim(str),str))

    call check_ret(nf_def_var  (fid,'frac',NF_DOUBLE ,2,vdid,vid))
    str   = 'fraction of grid cell that is active'
    call check_ret(nf_put_att_text(fid,vid,"long_name",len_trim(str),str))
    str   = 'xc yc'
    call check_ret(nf_put_att_text(fid,vid,"coordinates",len_trim(str),str))
    str   = 'unitless'
    call check_ret(nf_put_att_text(fid,vid,"note"     ,len_trim(str),str))
    write(str,'(a,g14.7)') 'error if frac> 1.0+eps or frac < 0.0-eps; eps =',eps
    call check_ret(nf_put_att_text(fid,vid,'filter1' ,len_trim(str),str))
    write(str,'(a,g14.7,a,g14.7)') 'limit frac to [fminval,fmaxval]; fminval=',fminval,' fmaxval=',fmaxval
    call check_ret(nf_put_att_text(fid,vid,'filter2' ,len_trim(str),str))

    call check_ret(nf_enddef(fid))

    if (units_xc(1:7) == 'radians') then
       xc = xc * 180._r8 / pi
       xv = xv * 180._r8 / pi
    end if
    if (units_yc(1:7) == 'radians') then
       yc = yc * 180._r8 / pi
       yv = yv * 180._r8 / pi
    end if

    call check_ret(nf_inq_varid(fid, 'xc', vid)) 
    call check_ret(nf_put_var_double(fid,  vid , xc))

    call check_ret(nf_inq_varid(fid,  'yc',vid)) 
    call check_ret(nf_put_var_double(fid,  vid , yc)) 

    call check_ret(nf_inq_varid(fid,  'xv',vid)) 
    call check_ret(nf_put_var_double(fid,  vid , xv)) 

    call check_ret(nf_inq_varid(fid,  'yv',vid)) 
    call check_ret(nf_put_var_double(fid,  vid , yv)) 

    call check_ret(nf_inq_varid(fid,'mask',vid)) 
    call check_ret(nf_put_var_int   (fid,  vid ,mask)) 

    call check_ret(nf_inq_varid(fid,'area',vid)) 
    call check_ret(nf_put_var_double(fid,  vid ,area)) 

    call check_ret(nf_inq_varid(fid,'frac',vid)) 
    call check_ret(nf_put_var_double(fid,  vid ,frac)) 

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
