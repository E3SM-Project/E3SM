! Utility functions in support of PIO io interface
module cam_pio_utils

  use pio, only : io_desc_t, iosystem_desc_t, file_desc_t, pio_double, pio_real, pio_freedecomp, &
       pio_offset
  use shr_kind_mod, only : r8=>shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use cam_logfile,      only: iulog
  use perf_mod,         only: t_startf, t_stopf
  use spmd_utils,       only: masterproc
  use cam_history_support, only : fillvalue, hist_coords, column_info, max_chars, field_info

  implicit none
  private
  save

  public :: column_info, cam_pio_openfile, cam_pio_createfile
  public :: get_decomp
  public :: get_phys_decomp
  public :: get_dyn_decomp
  public :: init_pio_subsystem  ! called from cam_comp
  public :: clean_iodesc_list


  integer, parameter, public :: dyn_stagger_decomp=102,dyn_decomp=101,phys_decomp=100

  integer :: pio_iotype

  ! This variable should be private ?
  type(iosystem_desc_t), pointer, public :: pio_subsystem => null()

  type iodesc_list
     integer(i8) :: tag
     type(io_desc_t), pointer :: iodesc => null()
     type(iodesc_list), pointer :: next => null()
  end type iodesc_list
  type(iodesc_list), target :: iodesc_list_top

!  logical :: dumpit=.false. ! for debugging only

!  Forward these cam_history_support variables
  public :: max_chars, fillvalue

  interface get_phys_decomp
     module procedure get_phys_decomp_md1d
     module procedure get_phys_decomp_mdnd
  end interface

contains

  subroutine init_pio_subsystem(nlfilename)
    use shr_pio_mod,   only: shr_pio_getiosys, shr_pio_getiotype
    use pio,          only: pio_rearr_box, pio_init
    use spmd_utils,   only: iam, mpicom
    use cam_instance, only: atm_id
    character(len=*) nlfilename

    pio_subsystem => shr_pio_getiosys(atm_id)
    pio_iotype =  shr_pio_getiotype(atm_id)

  end subroutine init_pio_subsystem

  subroutine get_decomp (iodesc, field, dtype,  numlev_in, column) 
    use dyn_grid,   only : get_horiz_grid_dim_d
    use dycore, only : dycore_is

    type(field_info), intent(in) :: field
    integer, intent(in) :: dtype
    
    type(column_info), intent(in), optional :: column
    integer, optional, intent(in) :: numlev_in

    character(len=3) :: memorder
    type(io_desc_t), pointer :: iodesc

    integer :: hdim1_d, hdim2_d
    integer :: mdims(8), mdimcnt, i


    if(present(numlev_in)) then
       mdimcnt=1
       mdims(1)=numlev_in
    else if(associated(field%mdims)) then
       mdimcnt=size(field%mdims)
       do i=1,mdimcnt
          mdims(i)= hist_coords(field%mdims(i))%dimsize
       end do
    else
       mdimcnt=1
       mdims(1)=1
    end if


    call t_startf('get_decomp')

    memorder = 'xzy'
    if(present(column)) then
       if(field%decomp_type==phys_decomp) then
          call get_phys_decomp(iodesc, 1,mdims(1:mdimcnt),1,dtype,column=column)
       else  if(associated(column%hmask)) then
          call get_dyn_decomp( iodesc, column%column_cnt, 1, mdims(1), 0, dtype, memorder_in='xzy',column_in=column)
       else
          call get_dyn_decomp( iodesc, column%num_lons, column%num_lons, mdims(1), 0, dtype, memorder_in='xzy',column_in=column)
       end if

    else
       if(field%decomp_type==phys_decomp) then
          call get_phys_decomp(iodesc, 1,mdims(1:mdimcnt),1,dtype)
       else if(field%decomp_type==dyn_decomp) then
          call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
          call get_dyn_decomp( iodesc, hdim1_d, hdim2_d, mdims(1), 0, dtype, memorder_in='xzy')
       else if(field%decomp_type==dyn_stagger_decomp) then
          call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
          call get_dyn_decomp( iodesc, hdim1_d, hdim2_d-1, mdims(1), 0, dtype, memorder_in='xzy')
       end if
    end if

     call t_stopf('get_decomp')

   end subroutine get_decomp



  subroutine get_phys_decomp_md1d(iodesc, fdim, mdim, ldim, dtype, fileorder_in, column)

    type(io_desc_t), pointer :: iodesc
    integer, intent(in) :: fdim, mdim, ldim
    type(column_info), optional :: column
    character(len=3),optional, intent(in) :: fileorder_in
    integer, intent(in) :: dtype
    character(len=3) :: fileorder


    if(present(fileorder_in)) then
       fileorder=fileorder_in
    else
       fileorder='xyz'
    end if

    if(present(column)) then
       call get_phys_decomp_mdnd(iodesc,fdim,(/mdim/),ldim,dtype,fileorder_in=fileorder,column=column)
    else
       call get_phys_decomp_mdnd(iodesc,fdim,(/mdim/),ldim,dtype,fileorder_in=fileorder)
    end if
  end subroutine get_phys_decomp_md1d

 subroutine get_phys_decomp_mdnd(iodesc, fdim, mdim, ldim, dtype, fileorder_in, column)
   use pio, only : pio_initdecomp,  pio_offset, pio_setdebuglevel
   use pio_support, only : pio_writedof
   use dyn_grid,   only : get_horiz_grid_dim_d
   use dycore, only : dycore_is
   use spmd_utils,       only : iam, mpicom


    type(io_desc_t), pointer :: iodesc
    integer, intent(in) :: fdim, mdim(:), ldim
    type(column_info), optional :: column
    integer, intent(in) :: dtype
    character(len=3),optional, intent(in) :: fileorder_in
    character(len=3) :: fileorder

    integer :: dimlens(8)    
    integer :: dimcnt
    integer :: hdim1_d, hdim2_d
    integer(kind=pio_offset), pointer :: ldof(:)
    integer :: ierr, ldecomp
    logical :: twodhorizontal, found, oldcolumns
    integer :: mdimsize, mdimprod, i

    call t_startf('get_phys_decomp')
    
    oldcolumns=.false.
    twodhorizontal = .true.
    if(dycore_is('UNSTRUCTURED')) twodhorizontal = .false.

    if(present(fileorder_in)) then
       fileorder=fileorder_in
    else
       fileorder='xyz'
    end if

    mdimsize = size(mdim)
    mdimprod=product(mdim)


    if(present(column)) then
       if(associated(column%hmask)) then
          hdim1_d = column%column_cnt
          hdim2_d = 1
          twodhorizontal = .false.
       else
          oldcolumns=.true.
          hdim1_d=column%num_lons
          hdim2_d=column%num_lats
       end if
    else 
       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
    end if
    dimlens(1)=hdim1_d
    dimcnt=1

    if(fdim>1) then
       dimcnt=dimcnt+1
       dimlens(dimcnt)=fdim
    end if

    if(fileorder=='xyz') then
       ldecomp=5
       if(twodhorizontal) then
          dimcnt=dimcnt+1
          dimlens(dimcnt)=hdim2_d
       end if

       if(sum(mdim)> mdimsize) then
          dimcnt=dimcnt+1
          dimlens(dimcnt:dimcnt+mdimsize-1)=mdim(1:mdimsize)
          dimcnt=dimcnt+mdimsize-1
       end if
    else
       ldecomp=7
       if(sum(mdim)> mdimsize) then
          dimcnt=dimcnt+1
          dimlens(dimcnt:dimcnt+mdimsize-1)=mdim(1:mdimsize)
          dimcnt=dimcnt+mdimsize-1
       end if
       
       if(twodhorizontal) then
          dimcnt=dimcnt+1
          dimlens(dimcnt)=hdim2_d
       end if
    end if
    if(ldim>1) then
       dimcnt=dimcnt+1
       dimlens(dimcnt)=ldim
    end if
    found=.false.
    if(.not. oldcolumns) then
       if(present(column)) then
          ldecomp = ldecomp+10*column%htape
       end if
       call find_iodesc(dimcnt, dimlens(1:dimcnt), dtype, ldecomp, iodesc, found)
    end if
    if(.not.found) then
       if(oldcolumns) then
          ldof => get_column_ldof(phys_decomp, column, mdimprod, fileorder_in=fileorder)
       else if(present(column)) then
          ldof => get_phys_ldof(hdim1_d, hdim2_d, fdim, mdimprod,ldim, fileorder, column%hmask)           
       else
          ldof => get_phys_ldof(hdim1_d, hdim2_d, fdim, mdimprod,ldim, fileorder)           
       end if
       call pio_initdecomp(pio_subsystem, dtype,dimlens(1:dimcnt), ldof, iodesc)
       deallocate(ldof)
    end if

    call t_stopf('get_phys_decomp')


  end subroutine get_phys_decomp_mdnd


  subroutine find_iodesc(dimcnt, dimlens, dtype, decomptype, iodesc, found)
    use dycore, only: dycore_is
    integer, intent(in) :: dimcnt
    integer, intent(in) :: dimlens(dimcnt)
    integer, intent(in) :: dtype
    integer, intent(in) :: decomptype
    type(io_desc_t), pointer :: iodesc
    logical,intent(out) :: found
    type(iodesc_list), pointer :: this, prev
    integer :: i
    integer(i8) :: tag, j
    integer :: multiplier = 1000

    found = .false.
    this => iodesc_list_top

    if(dycore_is('UNSTRUCTURED')) multiplier = 10000

    j=1
    tag = 0
    do i=1,dimcnt
       tag = tag+dimlens(i)*j
       j=j*multiplier
    end do
    tag = tag+j*int((dimcnt*1000+dtype*100+decomptype),i8)


    do while(associated(this) .and. .not. found)
       if(tag==this%tag) then
          found=.true.
          iodesc => this%iodesc
       else
          prev=>this
          this=>this%next
       end if
    end do
    if(.not.found) then
       this=>prev
       if(associated(this%iodesc)) then
          allocate(this%next)
          this=>this%next
       end if
       allocate(this%iodesc)
       this%tag = tag
       iodesc=>this%iodesc
       if(masterproc) write(iulog,*) 'Creating new decomp: ',this%tag
       nullify(this%next)
    end if
!    if(masterproc) write(iulog,*) 'Using decomp: ',this%tag
    
  end subroutine find_iodesc


! Deallocate all entries in the iodesc list

  subroutine clean_iodesc_list()
    type(iodesc_list), pointer :: this, prev


    if(associated(iodesc_list_top%iodesc)) then

       this => iodesc_list_top
       iodesc_list_top%tag = -1
       call pio_freedecomp(pio_subsystem, this%iodesc)
       deallocate(this%iodesc)
       nullify(this%iodesc)
       this => this%next
       nullify(iodesc_list_top%next)
       
       do while(associated(this))
          call pio_freedecomp(pio_subsystem, this%iodesc)
          deallocate(this%iodesc)
          prev=>this
          this=>this%next
          deallocate(prev)

       end do
    end if
  end subroutine clean_iodesc_list





  subroutine get_dyn_decomp(iodesc, hdim1, hdim2, nlev, ncnst, dtype, memorder_in, fileorder_in,column_in)
    use dycore, only : dycore_is
    use pio, only : pio_initdecomp, pio_setdebuglevel
    use abortutils, only : endrun
    type(io_desc_t), pointer :: iodesc
    integer, intent(in) :: hdim1, hdim2, nlev, ncnst
    integer, intent(in) :: dtype
    character(len=*), optional :: memorder_in
    character(len=*), optional :: fileorder_in
    type(column_info), optional :: column_in

    character(len=3) :: memorder, fileorder
    logical :: twodhorizontal, found
    integer :: dimcnt, dimlens(4)
    integer(kind=pio_offset), pointer :: ldof(:)
    integer :: tdecomp
    logical :: oldcolumns

    call t_startf('get_dyn_decomp')


    twodhorizontal = .true.
    oldcolumns=.false.

    dimlens(1)=hdim1
    if(dycore_is('UNSTRUCTURED')) twodhorizontal = .false.

    if(present(column_in)) then
       if(associated(column_in%hmask)) then
          oldcolumns=.false.
          twodhorizontal = .false.
       else
          oldcolumns=.true.
          dimlens(1)=column_in%num_lons
       end if
    end if

    if(present(memorder_in)) then
       memorder=memorder_in
    else if(dycore_is('LR')) then
       memorder='xyz'
    else
       memorder='xzy'
    end if
    if(memorder.eq.'xyz') then
       tdecomp=1
    else
       tdecomp=3
    end if

    if(present(fileorder_in)) then
       fileorder = fileorder_in
    else
       fileorder = 'xyz'
    end if


    dimcnt=1
    if(fileorder=='xyz') then
       tdecomp=tdecomp+1
       if(twodhorizontal) then
          dimcnt=dimcnt+1
          if(oldcolumns) then
             dimlens(dimcnt)=column_in%num_lats
          else
             dimlens(dimcnt)=hdim2
          end if
       end if
       if(nlev>1) then
          dimcnt=dimcnt+1
          dimlens(dimcnt)=nlev
       end if
    else
       if(nlev>1) then
          dimcnt=dimcnt+1
          dimlens(dimcnt)=nlev
       end if
       if(twodhorizontal) then
          dimcnt=dimcnt+1
          if(oldcolumns) then
             dimlens(dimcnt)=column_in%num_lats
          else
             dimlens(dimcnt)=hdim2
          end if
       end if
    end if


    if(ncnst>0) then
       dimcnt=dimcnt+1
       dimlens(dimcnt)=ncnst
    end if
    found=.false.
    if(.not.oldcolumns) then
       call find_iodesc(dimcnt, dimlens(1:dimcnt), dtype, tdecomp, iodesc, found)
    end if
    if(.not. found) then       
       if(oldcolumns) then
          ldof => get_column_ldof(dyn_decomp, column_in, nlev, fileorder, memorder)
       else if(present(column_in)) then
          ldof => get_dyn_ldof(column_in%column_cnt, 1, nlev*max(1,ncnst),memorder_in=memorder, &
               fileorder_in=fileorder,column_in=column_in)
       else
          ldof => get_dyn_ldof(hdim1, hdim2, nlev*max(1,ncnst),memorder_in=memorder,fileorder_in=fileorder)
       end if
       call pio_initdecomp(pio_subsystem, dtype, dimlens(1:dimcnt), ldof, iodesc)

       deallocate(ldof)
    end if


    call t_stopf('get_dyn_decomp')

  end subroutine get_dyn_decomp

  !
  ! Get the integer mapping of a variable in the dynamics decomp in memory.  
  ! The canonical ordering is as on the file. A 0 value indicates that the
  ! variable is not on the file (eg halo or boundary values)
  !
  function get_dyn_ldof(hdim1_d, hdim2_d, nlev, fileorder_in, memorder_in, column_in) result(ldof)
    use dyn_grid, only : get_gcol_block_d, get_block_owner_d, get_dyn_grid_parm, &
         get_gcol_block_cnt_d, get_block_gcol_cnt_d, get_block_gcol_d, get_block_bounds_d, get_block_ldof_d
    use spmd_utils, only : iam
    use dycore, only : dycore_is
    use abortutils, only : endrun

    integer, intent(in) :: hdim1_d, hdim2_d, nlev
    integer(kind=pio_offset), pointer :: ldof(:)
    integer :: i, block_cnt, max_block_cnt, b, k, ii, j
    integer :: lcnt, ngcols
    integer, allocatable :: gcols(:)
    character(len=3),optional, intent(in) :: fileorder_in, memorder_in
    type(column_info), optional, intent(in) :: column_in

    character(len=3) :: fileorder, memorder
    integer :: bfirst, blast, ncols

    integer :: beglatxy, beglonxy, endlatxy, endlonxy, plat
    integer :: nelemd, npsq
     
    logical, allocatable :: myblock(:)
    integer, pointer :: hmask(:,:)
    logical :: newcolumns


    beglonxy = get_dyn_grid_parm('beglonxy')
    endlonxy = get_dyn_grid_parm('endlonxy')
    beglatxy = get_dyn_grid_parm('beglatxy')
    endlatxy = get_dyn_grid_parm('endlatxy')
    newcolumns=.false.

    if(present(column_in) .and. hdim2_d==1) then
       plat = 1
       if(associated(column_in%hmask_dyn)) then
          hmask => column_in%hmask_dyn
          newcolumns=.true.
       end if
    else
       plat = get_dyn_grid_parm('plat')
    end if

    if(present(fileorder_in)) then
       fileorder=fileorder_in
    else		
       fileorder='xyz'
    end if

    if(present(memorder_in)) then
       memorder=memorder_in
    else if(dycore_is('LR')) then
       memorder='xyz'
    else
       memorder='xzy'
    end if

    ngcols = hdim1_d*hdim2_d

    if(dycore_is('SE')) then
       !
       ! dyn_decomp variables in spectral element are stored as (1:npsq,nlev,nelemd)
       !
       call get_block_ldof_d(nlev,ldof)


    else
       lcnt=(endlatxy-beglatxy+1)*nlev*(endlonxy-beglonxy+1)
       allocate(ldof(lcnt))
       lcnt=0
       ldof(:)=0	
       if(newcolumns) then
          if(memorder.eq.'xzy') then
             do j=1,endlatxy-beglatxy+1
                do k=0,nlev-1
                   do i=1,endlonxy-beglonxy+1
                      lcnt=lcnt+1
                      if(hmask(i,j)>0) then
                         ldof(lcnt) = hmask(i,j)+k*hdim1_d
                      end if
                   end do
                end do
             end do
          else
             do k=0,nlev-1
                do j=1,endlatxy-beglatxy+1
                   do i=1,endlonxy-beglonxy+1
                      lcnt=lcnt+1
                      if(hmask(i,j)>0) then
                         ldof(lcnt) = hmask(i,j)+k*hdim1_d
                      end if
                   end do
                end do
             end do
          end if
       else
          if(memorder.eq.'xzy') then
             do j=beglatxy,endlatxy
                do k=1,nlev
                   do i=beglonxy, endlonxy
                      lcnt=lcnt+1
                      if(j.eq.1.and.hdim2_d==plat-1) then
                         ldof(lcnt)=0
                      else
                         if(fileorder.eq.'xyz') then
                            ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d+(k-1)*hdim1_d*hdim2_d
                         else if(fileorder.eq.'xzy') then
                            ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d*nlev+(k-1)*hdim1_d
                         end if
                      end if
                   end do
                end do
             end do
          else  ! if(memorder.eq.'xyz') then
             do k=1,nlev
                do j=beglatxy,endlatxy
                   do i=beglonxy, endlonxy
                      lcnt=lcnt+1
                      if(j.eq.1.and.hdim2_d==plat-1) then
                         ldof(lcnt)=0
                      else if(hdim2_d>1) then
                         if(fileorder.eq.'xyz') then
                            ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d+(k-1)*hdim1_d*hdim2_d
                         else 
                            ldof(lcnt)=i+(j-(plat-hdim2_d+1))*hdim1_d*nlev+(k-1)*hdim1_d
                         end if
                      else  ! lon, lev decomp used for history nacs
                         ldof(lcnt)=i+(k-1)*hdim1_d
                      end if
                   end do
                end do
             end do
          end if
       end if
    end if



  end function get_dyn_ldof




  subroutine cam_pio_openfile(file, fname, mode, is_init)
    use pio, only : pio_openfile, file_desc_t, pio_noerr, pio_noclobber
    use abortutils, only : endrun
    type(file_desc_t), intent(inout), target :: file
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode
    logical, optional, intent(in) :: is_init

    integer :: ierr

    ierr = pio_openfile(pio_subsystem, file, pio_iotype, fname, mode)

    if(ierr/= PIO_NOERR) then
       call endrun('Failed to open restart file to read')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened existing file ', trim(fname), file%fh
    end if

  end subroutine cam_pio_openfile

  subroutine cam_pio_createfile(file, fname, mode)
    use pio, only : pio_createfile, file_desc_t, pio_noerr, pio_clobber, pio_64bit_offset
    use cam_control_mod, only : use_64bit_nc
    use abortutils, only : endrun
    type(file_desc_t), intent(inout) :: file
    character(len=*), intent(in) :: fname
    integer, intent(in) :: mode
    integer :: ierr


    if(use_64bit_nc) then
       ierr = pio_createfile(pio_subsystem, file, pio_iotype, fname, ior(PIO_CLOBBER,PIO_64BIT_OFFSET))
    else
       ierr = pio_createfile(pio_subsystem, file, pio_iotype, fname, PIO_CLOBBER)
    end if

    if(ierr/= PIO_NOERR) then
       call endrun('Failed to open restart file to write')
    else if(pio_subsystem%io_rank==0) then
       write(iulog,*) 'Opened file ', trim(fname),  ' to write', file%fh
    end if


  end subroutine cam_pio_createfile

#ifdef OLD
  subroutine set_pio_parameters(io_type_name)
    use pio, only : iotype_netcdf, iotype_pnetcdf, iotype_binary, iotype_pbinary, &
         pio_iotype_netcdf4p, pio_iotype_netcdf4c
    use spmd_utils, only : nsmps, npes
    use cam_logfile,      only: iulog
    use shr_string_mod,   only: shr_string_toUpper

    character(len=*), intent(in) :: io_type_name
    character(len=16) :: tmpname

    tmpname = shr_string_toupper(io_type_name)


    if(tmpname .eq. 'NETCDF') then
       io_type = iotype_netcdf
    else if(tmpname .eq. 'PNETCDF') then
       io_type = iotype_pnetcdf
    else if(tmpname .eq. 'NETCDF4P') then
       io_type = pio_iotype_netcdf4p
    else if(tmpname .eq. 'NETCDF4C') then
       io_type = pio_iotype_netcdf4c
    else
       write(iulog,*) 'Bad io_type argument - using iotype_netcdf'
       io_type=iotype_netcdf
    end if

    if(io_stride>0.and.num_iotasks<0) then
       num_iotasks = npes/io_stride
    else if(num_iotasks>0 .and. io_stride<0) then
       io_stride = npes/num_iotasks
    else if(num_iotasks<0 .and. io_stride<0) then
       io_stride = max(min(npes,4),npes/8)
       num_iotasks = npes/io_stride
    end if


    if(io_stride*num_iotasks> npes .or. io_stride<=0 .or. num_iotasks<=0) then
       write(iulog,*) 'io_stride or num_iotasks out of bounds - resetting to defaults ',io_stride, num_iotasks 
       io_stride = max(1,npes/(nsmps*4))
       num_iotasks = npes/io_stride
    end if
    write(iulog,*) 'Using io_type ',tmpname,' with stride ',io_stride,' and ',num_iotasks,' io tasks'

  end subroutine set_pio_parameters
#endif
  !
  ! Get the integer mapping of a variable in the physics decomp in memory.  
  ! The canonical ordering is as on the file. A 0 value indicates that the
  ! variable is not on the file (eg halo or boundary values) hmask can be used 
  ! to define a specific subset of global size hsize 
  !

  function get_phys_ldof(hdim1_d,hdim2_d,fdim,mdim,ldim , fileorder_in, hmask) result(ldof)
    use phys_grid, only : get_gcol_all_p, get_ncols_p, get_lon_all_p, get_lat_all_p
    use ppgrid, only : pcols, begchunk, endchunk
    use dycore, only : dycore_is
    use spmd_utils, only : mpicom, iam
    use abortutils, only : endrun

    integer, intent(in) :: hdim1_d, hdim2_d, fdim,mdim,ldim
    character(len=3),optional, intent(in) :: fileorder_in
    integer, optional    :: hmask(pcols,begchunk:endchunk)
    integer(kind=pio_offset), pointer :: ldof(:)

    integer :: hsize

    integer :: gcols(pcols), ilat(pcols), ilon(pcols)
    integer :: ncols, lchnk, i, k, lcnt, f, ii
    character(len=3) :: fileorder
    integer :: ierr
    logical, pointer :: localhmask(:,:)

    allocate(localhmask(pcols,begchunk:endchunk))

    if(present(hmask)) then
       localhmask=(hmask>0)
    else
       localhmask=.false.
       do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          localhmask(1:ncols,lchnk)=.true.
       end do
    end if
    hsize= hdim1_d*hdim2_d
    


    if(present(fileorder_in)) then
       fileorder=fileorder_in
    else		
       fileorder='xyz'
    end if


    allocate(ldof(pcols*(endchunk-begchunk+1)*fdim*mdim*ldim))

    if(fileorder.eq.'xzy') then
       lcnt=0
       do lchnk=begchunk,endchunk
!          ncols = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, pcols, ilon)
          call get_lat_all_p(lchnk, pcols, ilat)
          do f=1,fdim
             do k=1,mdim
                do i=1,pcols
                   do ii=1,ldim
                      lcnt=lcnt+1
                      if(localhmask(i,lchnk)) then
!                      if(i<=ncols) then
                         ldof(lcnt) = (ilon(i)-1)*ldim+(k-1)*hdim1_d*ldim+ &
                              (ilat(i)-1)*hdim1_d*mdim*ldim+(f-1)*hsize*mdim*ldim+ii
                      else
                         ldof(lcnt)=0
                      end if

                   end do
                end do
             end do
          end do
       end do
    else
       lcnt=0
       do ii=1,ldim
          do lchnk=begchunk,endchunk
             if(present(hmask)) then
                gcols = hmask(:,lchnk)
             else
                call get_gcol_all_p(lchnk, pcols, gcols)
             end if


             do k=1,mdim
                do i=1,pcols
                   do f=1,fdim
                      lcnt=lcnt+1
                      if(localhmask(i,lchnk)) then
                         ldof(lcnt) = (f-1)+fdim*(gcols(i)+(k-1)*hsize+(ii-1)*hsize*mdim)
                      else
                         ldof(lcnt)=0
                      end if
                   end do
                end do
             end do
          end do
       end do
    end if

  end function get_phys_ldof

  function get_column_ldof(decomp_type, column, nlev, fileorder_in, memorder_in) result(ldof)
    use dyn_grid, only : get_dyn_grid_parm
    use phys_grid, only : get_gcol_all_p, get_ncols_p, get_lon_all_p, get_lat_all_p
    use ppgrid, only : pcols, begchunk, endchunk

    integer, intent(in) :: decomp_type
    type(column_info), intent(in) :: column
    integer, intent(in) :: nlev
    integer :: hsize
    integer(kind=pio_offset), pointer :: ldof(:)
    integer :: gcols(pcols), ilat(pcols), ilon(pcols)
    integer :: ncols, lchnk, i, j, k, lcnt, num_lats, num_lons
    integer :: lon1, lon2, lat1, lat2, beglatxy, endlatxy, beglonxy, endlonxy
    character(len=3),optional, intent(in) :: fileorder_in, memorder_in
    character(len=3) :: fileorder, memorder


    if(present(fileorder_in)) then
       fileorder=fileorder_in
    else		
       fileorder='xyz'
    end if
    if(present(memorder_in)) then
       memorder=memorder_in
    else		
       memorder='xzy'
    end if

    num_lons= column%num_lons
    num_lats= column%num_lats
    lat1=column%columnlat(1)
    lat2=column%columnlat(2)
    lon1=column%columnlon(1)
    lon2=column%columnlon(2)
    hsize= num_lats*num_lons

    
    lcnt=0

    if(decomp_type==phys_decomp) then
       allocate(ldof(pcols*(endchunk-begchunk+1)*nlev))

       do lchnk=begchunk,endchunk
          ncols = get_ncols_p(lchnk)
          call get_lon_all_p(lchnk, pcols, ilon)
          call get_lat_all_p(lchnk, pcols, ilat)

          if(fileorder.eq.'xzy') then
             do k=1,nlev
                do i=1,pcols
                   lcnt=lcnt+1
                   if(i<=ncols.and.&
                        ilat(i)>=lat1.and. &
                        ilat(i)<=lat2.and. &
                        ilon(i)>=lon1.and. &
                        ilon(i)<=lon2) then
                      ldof(lcnt) = 1+(ilon(i)-lon1)+(k-1)*num_lons+&
                           (ilat(i)-lat1)*num_lons*nlev
                   else
                      ldof(lcnt)=0
                   end if
                end do
             end do
          else
             do k=1,nlev
                do i=1,pcols
                   lcnt=lcnt+1
                   if(i<=ncols.and.&
                        ilat(i)>=lat1.and. &
                        ilat(i)<=lat2.and. &
                        ilon(i)>=lon1.and. &
                        ilon(i)<=lon2) then
                      ldof(lcnt) = 1+(ilon(i)-lon1)+(k-1)*hsize+&
                           (ilat(i)-lat1)*num_lons
                   else
                      ldof(lcnt)=0
                   end if
                end do
             end do
          end if
       end do
    else
       beglonxy = get_dyn_grid_parm('beglonxy')
       endlonxy = get_dyn_grid_parm('endlonxy')
       beglatxy = get_dyn_grid_parm('beglatxy')
       endlatxy = get_dyn_grid_parm('endlatxy')
       allocate(ldof((endlonxy-beglonxy+1)*(endlatxy-beglatxy+1)*nlev))
       ldof = 0
       if(memorder.eq.'xzy') then
          do j=beglatxy,endlatxy
             do k=0,nlev-1
                do i=beglonxy, endlonxy
                   lcnt=lcnt+1
                   if(  j>=lat1.and. &
                        j<=lat2.and. &
                        i>=lon1.and. &
                        i<=lon2) then

                      if(fileorder.eq.'xyz') then
                         ldof(lcnt)=1+(i-lon1)+(j-lat1)*num_lons+k*hsize
                      else if(fileorder.eq.'xzy') then
                         ldof(lcnt)=1+(i-lon1)+(j-lat1)*num_lons*nlev+k*num_lons
                      end if
                   end if
                end do
             end do
          end do
       else
       
         do k=0,nlev-1
             do j=beglatxy,endlatxy
                do i=beglonxy, endlonxy
                   lcnt=lcnt+1
                   if(  j>=lat1.and. &
                        j<=lat2.and. &
                        i>=lon1.and. &
                        i<=lon2) then
                      if(fileorder.eq.'xyz') then
                         ldof(lcnt)=1+(i-lon1)+(j-lat1)*num_lons+k*hsize
                      else 
                         ldof(lcnt)=1+(i-lon1)+(j-lat1)*num_lons*nlev+k*num_lons
                      end if
                   end if
                end do
             end do
          end do
       end if
    end if
    
!    print *,__FILE__,__LINE__,hsize,nlev,sum(ldof),maxval(ldof)

  end function get_column_ldof


end module cam_pio_utils
