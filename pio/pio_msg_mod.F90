#define __PIO_FILE__ "pio_msg_mod.F90"
module pio_msg_mod
  use pio_kinds
  use pio_types
  use pio_support, only : piodie, DebugAsync

  implicit none
  private
  public :: pio_msg_handler_init, pio_msg_handler


  public :: add_to_file_list, lookupfile, delete_from_file_list, lookupiodesc, add_to_iodesc_list

! PIO ASYNC MESSAGE TAGS
   integer, parameter, public :: pio_msg_create_file = 300
   integer, parameter, public :: pio_msg_open_file = 301
   integer, parameter, public :: pio_msg_close_file = 302
   integer, parameter, public :: pio_msg_def_dim = 310
   integer, parameter, public :: pio_msg_def_var = 312
   integer, parameter, public :: pio_msg_enddef = 313
   integer, parameter, public :: pio_msg_redef = 314
   integer, parameter, public :: pio_msg_initdecomp_dof = 315
   
   integer, parameter, public :: pio_msg_writedarray = 320
   integer, parameter, public :: pio_msg_readdarray = 325
   integer, parameter, public :: pio_msg_inquire = 330
   integer, parameter, public :: pio_msg_inq_att = 331
   integer, parameter, public :: pio_msg_inq_attname = 332
   integer, parameter, public :: pio_msg_inq_varid = 333
   integer, parameter, public :: pio_msg_inq_varname = 334
   integer, parameter, public :: pio_msg_inq_vardimid = 335
   integer, parameter, public :: pio_msg_inq_varnatts = 336
   integer, parameter, public :: pio_msg_inq_varndims = 337
   integer, parameter, public :: pio_msg_inq_vartype = 338
   integer, parameter, public :: pio_msg_inq_dimid = 339
   integer, parameter, public :: pio_msg_inq_dimlen = 340
   integer, parameter, public :: pio_msg_inq_dimname = 341
   integer, parameter, public :: pio_msg_inq_attlen = 342
   integer, parameter, public :: pio_msg_seterrorhandling = 350

   integer, parameter, public :: pio_msg_getvar1 = 360
   integer, parameter, public :: pio_msg_getvar_0d = 361
   integer, parameter, public :: pio_msg_getvar_1d = 362
   integer, parameter, public :: pio_msg_getvar_2d = 363
   integer, parameter, public :: pio_msg_getvar_3d = 364
   integer, parameter, public :: pio_msg_getvar_4d = 365
   integer, parameter, public :: pio_msg_getvar_5d = 366

   integer, parameter, public :: pio_msg_getvara_1d = 367
   integer, parameter, public :: pio_msg_getvara_2d = 368
   integer, parameter, public :: pio_msg_getvara_3d = 369
   integer, parameter, public :: pio_msg_getvara_4d = 370
   integer, parameter, public :: pio_msg_getvara_5d = 371

   integer, parameter, public :: pio_msg_putvar1 = 380
   integer, parameter, public :: pio_msg_putvar_0d = 381
   integer, parameter, public :: pio_msg_putvar_1d = 382
   integer, parameter, public :: pio_msg_putvar_2d = 383
   integer, parameter, public :: pio_msg_putvar_3d = 384
   integer, parameter, public :: pio_msg_putvar_4d = 385
   integer, parameter, public :: pio_msg_putvar_5d = 386

   integer, parameter, public :: pio_msg_putvara_1d = 387
   integer, parameter, public :: pio_msg_putvara_2d = 388
   integer, parameter, public :: pio_msg_putvara_3d = 389
   integer, parameter, public :: pio_msg_putvara_4d = 390
   integer, parameter, public :: pio_msg_putvara_5d = 391

   integer, parameter, public :: pio_msg_getatt = 400
   integer, parameter, public :: pio_msg_getatt_1d = 401
   integer, parameter, public :: pio_msg_putatt = 402
   integer, parameter, public :: pio_msg_putatt_1d = 403

   integer, parameter, public :: PIO_MSG_SYNC_FILE = 500
   integer, parameter, public :: pio_msg_exit = 999   

   
   type :: file_desc_list
      type(file_desc_t), pointer :: file
      type(file_desc_list), pointer :: next
   end type file_desc_list

   type(file_desc_list), target :: top_file

   type :: io_desc_list
      type(io_desc_t), pointer :: iodesc
      type(io_desc_list), pointer :: next
   end type io_desc_list

   type(io_desc_list), target :: top_iodesc

   integer :: io_comm, iorank


contains

  subroutine pio_msg_handler_init(io_comm_in, io_rank_in)
    integer, intent(in) :: io_comm_in, io_rank_in

    io_comm = io_comm_in
    iorank = io_rank_in
    
    nullify(top_file%file)
    nullify(top_file%next)
    nullify(top_iodesc%iodesc)
    nullify(top_iodesc%next)

  end subroutine pio_msg_handler_init


  subroutine pio_msg_handler(numcomps, iosystem)
!    use pio_types, only : 
#ifdef TIMING
    use perf_mod        ! _EXTERNAL
#endif

    implicit none
    integer, intent(in) :: numcomps
    type(iosystem_desc_t), target :: iosystem(numcomps)
    type(iosystem_desc_t), pointer :: ios
    integer :: msg = 0, ierr
    include 'mpif.h' ! _EXTERNAL
 
    integer :: status(MPI_STATUS_SIZE)
    integer :: req(numcomps)
    integer :: index

#ifdef TIMING    
    call t_startf('pio_msg_mod')
#endif
    if(iorank==0) then
       do index=1,numcomps
          ios=>iosystem(index)
          if(ios%io_rank==0) then
             call mpi_irecv(msg, 1, mpi_integer, ios%comproot, 1, ios%union_comm, req(index), ierr)       
          end if
       enddo
    end if
    do while(msg /= pio_msg_exit)
       if(iorank==0) then
          if(Debugasync) print *,__PIO_FILE__,__LINE__, ' waiting'
          call mpi_waitany(numcomps, req, index, status, ierr)
          if(Debugasync) print *,__PIO_FILE__,__LINE__, ' recieved on ', index
       end if

       call mpi_bcast(index, 1, mpi_integer, 0, io_comm, ierr)
       ios => iosystem(index)

       if(Debugasync) print *,__PIO_FILE__,__LINE__, index, ios%intercomm
       call mpi_bcast(msg, 1, mpi_integer, 0, io_comm, ierr)
       if(Debugasync) print *,__PIO_FILE__,__LINE__,msg


       select case(msg) 
       case (PIO_MSG_CREATE_FILE)
          call create_file_handler(ios)
       case (PIO_MSG_OPEN_FILE)
          call open_file_handler(ios)
       case (PIO_MSG_INITDECOMP_DOF)
          call initdecomp_dof_handler(ios)
       case (PIO_MSG_WRITEDARRAY)
          call writedarray_handler(ios)
       case (PIO_MSG_READDARRAY)
          call readdarray_handler(ios)
       case (PIO_MSG_SETERRORHANDLING)
          call seterrorhandling_handler(ios)
       case (PIO_MSG_GETVAR1)
          call var1_handler(ios, msg)
       case (PIO_MSG_GETVAR_0d)
          call var_0d_handler(ios, msg)
       case (PIO_MSG_GETVAR_1d)
          call var_1d_handler(ios, msg)
       case (PIO_MSG_GETVAR_2d)
          call var_2d_handler(ios, msg)
       case (PIO_MSG_GETVAR_3d)
          call var_3d_handler(ios, msg)
       case (PIO_MSG_GETVAR_4d)
          call var_4d_handler(ios, msg)
       case (PIO_MSG_GETVAR_5d)
          call var_5d_handler(ios, msg)
       case (PIO_MSG_GETVARA_1d)
          call vara_1d_handler(ios, msg)
       case (PIO_MSG_GETVARA_2d)
          call vara_2d_handler(ios, msg)
       case (PIO_MSG_GETVARA_3d)
          call vara_3d_handler(ios, msg)
       case (PIO_MSG_GETVARA_4d)
          call vara_4d_handler(ios, msg)
       case (PIO_MSG_GETVARA_5d)
          call vara_5d_handler(ios, msg)

       case (PIO_MSG_PUTVAR1)
          call var1_handler(ios, msg)
       case (PIO_MSG_PUTVAR_0d)
          call var_0d_handler(ios, msg)
       case (PIO_MSG_PUTVAR_1d)
          call var_1d_handler(ios, msg)
       case (PIO_MSG_PUTVAR_2d)
          call var_2d_handler(ios, msg)
       case (PIO_MSG_PUTVAR_3d)
          call var_3d_handler(ios, msg)
       case (PIO_MSG_PUTVAR_4d)
          call var_4d_handler(ios, msg)
       case (PIO_MSG_PUTVAR_5d)
          call var_5d_handler(ios, msg)

       case (PIO_MSG_PUTVARA_1d)
          call vara_1d_handler(ios, msg)
       case (PIO_MSG_PUTVARA_2d)
          call vara_2d_handler(ios, msg)
       case (PIO_MSG_PUTVARA_3d)
          call vara_3d_handler(ios, msg)
       case (PIO_MSG_PUTVARA_4d)
          call vara_4d_handler(ios, msg)
       case (PIO_MSG_PUTVARA_5d)
          call vara_5d_handler(ios, msg)
       case (PIO_MSG_GETATT)
          call att_handler(ios, msg)
       case (PIO_MSG_GETATT_1D)
          call att_1d_handler(ios, msg)
       case (PIO_MSG_PUTATT)
          call att_handler(ios, msg)
       case (PIO_MSG_PUTATT_1D)
          call att_1d_handler(ios, msg)          
       case (PIO_MSG_EXIT)
          call finalize_handler(ios)
          print *,'PIO Exiting'
       case default
          call pio_callback_handler(ios,msg)
       end select   
       if(iorank==0) then
          call mpi_irecv(msg, 1, mpi_integer, ios%comproot, 1, ios%union_comm, req(index), ierr)
       end if

    end do

#ifdef TIMING
    call t_stopf('pio_msg_mod')
    call t_finalizef()
#endif


    if(Debugasync) print *,__PIO_FILE__,__LINE__
    call mpi_finalize(ierr)
    stop

  end subroutine pio_msg_handler


  subroutine add_to_file_list(file)
    type(file_desc_t), pointer :: file
    type(file_desc_list), pointer :: list_item

    list_item=> top_file

    if(associated(list_item%file)) then
       do while(associated(list_item%file) .and. associated(list_item%next))
       if(Debugasync) print *,__PIO_FILE__,__LINE__,list_item%file%fh
          list_item => list_item%next
       end do
       if(associated(list_item%file)) then
          allocate(list_item%next)
          list_item=>list_item%next
          nullify(list_item%next)
       end if
    end if
    if(Debugasync) print *,__PIO_FILE__,__LINE__,file%fh
    list_item%file => file

  end subroutine add_to_file_list


  subroutine add_to_iodesc_list(iodesc)
    type(io_desc_t), pointer :: iodesc
    type(io_desc_list), pointer :: list_item
    integer :: id

    list_item=> top_iodesc

    id = 0
    if(associated(list_item%iodesc)) then
       do while(associated(list_item%iodesc) .and. associated(list_item%next))
          list_item => list_item%next
       end do
       if(associated(list_item%iodesc)) then
          id = max(id, list_item%iodesc%async_id+1)
          allocate(list_item%next)
          list_item=>list_item%next
          nullify(list_item%next)
       end if
    end if
    iodesc%async_id=id
    list_item%iodesc => iodesc

  end subroutine add_to_iodesc_list


  subroutine delete_from_file_list(fh)
    integer, intent(in) :: fh
    type(file_desc_list), pointer :: list_item
    integer :: fh1

    fh1 = abs(fh)

    list_item=> top_file

    do while(associated(list_item%file) )
       if(abs(list_item%file%fh) == fh1) then
          nullify(list_item%file)
          exit
       end if
       if(associated(list_item%next)) then
          list_item=>list_item%next
       else
          call piodie(__PIO_FILE__,__LINE__)
       end if
    end do




  end subroutine delete_from_file_list



  function lookupfile(fh) result(file)
    type(file_desc_t), pointer :: file
    integer, intent(in) :: fh
    type(file_desc_list), pointer :: list_item

    integer :: fh1
    
    fh1 = abs(fh)

    list_item=> top_file
    
    do while(associated(list_item%file) )
       if(abs(list_item%file%fh) == fh1) then
          file => list_item%file
          exit
       end if
       list_item=>list_item%next
    end do


  end function lookupfile

  function lookupiodesc(async_id) result(iodesc)
    type(io_desc_t), pointer :: iodesc
    integer, intent(in) :: async_id
    type(io_desc_list), pointer :: list_item


    list_item=> top_iodesc
    
    do while(associated(list_item%iodesc) )
       if(abs(list_item%iodesc%async_id) == async_id) then
          iodesc => list_item%iodesc
          exit
       end if
       list_item=>list_item%next
    end do


  end function lookupiodesc

end module pio_msg_mod

