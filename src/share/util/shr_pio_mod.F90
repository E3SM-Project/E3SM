module shr_pio_mod
  use pio
  use shr_kind_mod, only : shr_kind_CS, shr_kind_cl, shr_kind_in
  use shr_file_mod, only : shr_file_getunit, shr_file_freeunit
  use shr_log_mod,  only : shr_log_unit
  use shr_mpi_mod,  only : shr_mpi_bcast, shr_mpi_chkerr
  use shr_sys_mod,  only : shr_sys_abort
#ifndef NO_MPIMOD
  use mpi, only : mpi_comm_null, mpi_comm_world
#endif
  implicit none
#ifdef NO_MPIMOD
#include <mpif.h>
#endif
  private
  public :: shr_pio_init1
  public :: shr_pio_init2
  public :: shr_pio_getiosys
  public :: shr_pio_getiotype
  public :: shr_pio_getioroot
  public :: shr_pio_finalize

  interface shr_pio_getiotype
     module procedure shr_pio_getiotype_fromid, shr_pio_getiotype_fromname
  end interface
  interface shr_pio_getiosys
     module procedure shr_pio_getiosys_fromid, shr_pio_getiosys_fromname
  end interface
  interface shr_pio_getioroot
     module procedure shr_pio_getioroot_fromid, shr_pio_getioroot_fromname
  end interface
  interface shr_pio_getindex
     module procedure shr_pio_getindex_fromid, shr_pio_getindex_fromname
  end interface



  type pio_comp_t
     integer :: compid
     integer :: pio_root
     integer :: pio_stride
     integer :: pio_numiotasks
     integer :: pio_iotype
     integer :: pio_rearranger
  end type pio_comp_t

  character(len=16), allocatable :: io_compname(:)
  type(pio_comp_t), allocatable :: pio_comp_settings(:)
  type (iosystem_desc_t), allocatable, target :: iosystems(:)
  integer :: io_comm
  logical :: pio_async_interface
  integer, allocatable :: io_compid(:)
  integer :: pio_debug_level=0, pio_blocksize=0
  integer(kind=pio_offset_kind) :: pio_buffer_size_limit=-1
  type(pio_rearr_opt_t)  :: pio_rearr_opts
  integer :: total_comps=0

#define DEBUGI 1

#ifdef DEBUGI
  integer :: drank
#endif


contains
!>
!! @public
!! @brief should be the first routine called after mpi_init.
!! It reads the pio default settings from file drv_in, namelist pio_default_inparm
!! and, if pio_async_interface is true, splits the IO tasks away from the
!! Compute tasks.  It then returns the new compute comm in
!! Global_Comm and sets module variable io_comm.
!!
!<
  subroutine shr_pio_init1(ncomps, nlfilename, Global_Comm)
    integer, intent(in) :: ncomps
    character(len=*) :: nlfilename
    integer, intent(inout) :: Global_Comm


    integer :: i, pio_root, pio_stride, pio_numiotasks, pio_iotype, pio_rearranger
    integer :: mpigrp_world, mpigrp, ierr, mpicom
    character(*),parameter :: subName =   '(shr_pio_init1) '
    integer :: pelist(3,1)

    call shr_pio_read_default_namelist(nlfilename, Global_Comm, pio_stride, pio_root, pio_numiotasks, &
         pio_iotype, pio_async_interface, pio_rearranger)


    call MPI_comm_rank(Global_Comm, drank, ierr)

    io_comm = MPI_COMM_NULL
    allocate(pio_comp_settings(ncomps))
    do i=1,ncomps
       pio_comp_settings(i)%pio_root = pio_root
       pio_comp_settings(i)%pio_stride = pio_stride
       pio_comp_settings(i)%pio_numiotasks = pio_numiotasks
       pio_comp_settings(i)%pio_iotype = pio_iotype
       pio_comp_settings(i)%pio_rearranger = pio_rearranger
    end do
    if(pio_async_interface) then
#ifdef NO_MPI2
       call shr_sys_abort(subname//':: async IO requires an MPI2 compliant MPI library')
#else

       pelist(1,1) = pio_root
       pelist(2,1) = pio_root + (pio_numiotasks-1)*pio_stride
       pelist(3,1) = pio_stride

       call mpi_comm_group(GLOBAL_COMM, mpigrp_world, ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
       call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')
       call mpi_comm_create(GLOBAL_COMM, mpigrp, io_comm, ierr)

       call mpi_group_range_excl(mpigrp_world, 1, pelist, mpigrp,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')
       call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
       Global_COMM=mpicom

       print *,__FILE__,__LINE__,subname, ' complete'
#endif
    end if
    total_comps = ncomps
  end subroutine shr_pio_init1
!>
!! @public
!! @brief if pio_async_interface is true, tasks in io_comm do not return from this subroutine.
!!
!! if pio_async_interface is false each component namelist pio_inparm is read from compname_modelio.nml
!! Then a subset of each components compute tasks are Identified as IO tasks using the root, stride and count
!! variables to select the tasks.
!!
!<


  subroutine shr_pio_init2(comp_id, comp_name, comp_iamin, comp_comm, comp_comm_iam)
    use shr_string_mod, only : shr_string_toLower
    integer, intent(in) :: comp_id(:)
    logical, intent(in) :: comp_iamin(:)
    character(len=*), intent(in) :: comp_name(:)
    integer, intent(in) ::  comp_comm(:), comp_comm_iam(:)
    integer :: i
    integer :: ncomps
    character(len=shr_kind_cl) :: nlfilename, cname
    type(iosystem_desc_t) :: iosys
    character(*), parameter :: subName = '(shr_pio_init2) '

    if(pio_debug_level>0) then
       if(comp_comm_iam(1)==0) then
          write(shr_log_unit,*) 'Setting pio_debuglevel : ',pio_debug_level
       end if
       call pio_setdebuglevel(pio_debug_level)
    endif
    ! 0 is a valid value of pio_buffer_size_limit
    if(pio_buffer_size_limit>=0) then
       if(comp_comm_iam(1)==0) then
          write(shr_log_unit,*) 'Setting pio_buffer_size_limit : ',pio_buffer_size_limit
       end if
       call pio_set_buffer_size_limit(pio_buffer_size_limit)
    endif
    if(pio_blocksize>0) then
       if(comp_comm_iam(1)==0) then
          write(shr_log_unit,*) 'Setting pio_blocksize : ',pio_blocksize
       end if
       call pio_set_blocksize(pio_blocksize)
    endif




    allocate(io_compid(total_comps), io_compname(total_comps))

    io_compid = comp_id
    io_compname = comp_name
    allocate(iosystems(total_comps))

    if(pio_async_interface) then
       call pio_init(total_comps,mpi_comm_world, comp_comm, io_comm, iosystems, rearr_opts=pio_rearr_opts)
       i=1
    else
       do i=1,total_comps
          if(comp_iamin(i)) then
             cname = comp_name(i)
	     if(len_trim(cname) <= 3) then
                nlfilename=trim(shr_string_toLower(cname))//'_modelio.nml'
             else
                nlfilename=trim(shr_string_toLower(cname(1:3)))//'_modelio.nml_'//cname(4:8)
             endif

             call shr_pio_read_component_namelist(nlfilename , comp_comm(i), pio_comp_settings(i)%pio_stride, &
                  pio_comp_settings(i)%pio_root, pio_comp_settings(i)%pio_numiotasks, &
                  pio_comp_settings(i)%pio_iotype, pio_comp_settings(i)%pio_rearranger)
             call pio_init(comp_comm_iam(i), comp_comm(i), pio_comp_settings(i)%pio_numiotasks, 0, &
                  pio_comp_settings(i)%pio_stride, &
                  pio_comp_settings(i)%pio_rearranger, iosystems(i), &
                  base=pio_comp_settings(i)%pio_root, rearr_opts=pio_rearr_opts)
             if(comp_comm_iam(i)==0) then
                write(shr_log_unit,*) io_compname(i),' : pio_numiotasks = ',pio_comp_settings(i)%pio_numiotasks
                write(shr_log_unit,*) io_compname(i),' : pio_stride = ',pio_comp_settings(i)%pio_stride
                write(shr_log_unit,*) io_compname(i),' : pio_root = ',pio_comp_settings(i)%pio_root
                write(shr_log_unit,*) io_compname(i),' : pio_iotype = ',pio_comp_settings(i)%pio_iotype
             end if


          end if
       end do
    end if
    do i=1,total_comps
       if(comp_comm_iam(i)==0) then
          write(shr_log_unit,*) io_compname(i),' : pio_numiotasks = ',pio_comp_settings(i)%pio_numiotasks
          write(shr_log_unit,*) io_compname(i),' : pio_stride = ',pio_comp_settings(i)%pio_stride
          write(shr_log_unit,*) io_compname(i),' : pio_rearranger = ',pio_comp_settings(i)%pio_rearranger
          write(shr_log_unit,*) io_compname(i),' : pio_root = ',pio_comp_settings(i)%pio_root
          write(shr_log_unit,*) io_compname(i),' : pio_iotype = ',pio_comp_settings(i)%pio_iotype
       end if
    enddo


  end subroutine shr_pio_init2



!===============================================================================
  subroutine shr_pio_finalize(  )
    integer :: ierr
    integer :: i
    logical :: active
    do i=1,total_comps
!       print *,__FILE__,__LINE__,drank,i,iosystems(i)%iosysid
       call pio_finalize(iosystems(i), ierr)
    end do

  end subroutine shr_pio_finalize

!===============================================================================
  function shr_pio_getiotype_fromid(compid) result(io_type)
    integer, intent(in) :: compid
    integer :: io_type

    io_type = pio_comp_settings(shr_pio_getindex(compid))%pio_iotype

  end function shr_pio_getiotype_fromid


  function shr_pio_getiotype_fromname(component) result(io_type)
    ! 'component' must be equal to some element of io_compname(:)
    ! (but it is case-insensitive)
    character(len=*), intent(in) :: component
    integer :: io_type

    io_type = pio_comp_settings(shr_pio_getindex(component))%pio_iotype

  end function shr_pio_getiotype_fromname

!===============================================================================
  function shr_pio_getioroot_fromid(compid) result(io_root)
    ! 'component' must be equal to some element of io_compname(:)
    ! (but it is case-insensitive)
    integer, intent(in) :: compid
    integer :: io_root

    io_root = pio_comp_settings(shr_pio_getindex(compid))%pio_root

  end function shr_pio_getioroot_fromid

  function shr_pio_getioroot_fromname(component) result(io_root)
    ! 'component' must be equal to some element of io_compname(:)
    ! (but it is case-insensitive)
    character(len=*), intent(in) :: component
    integer :: io_root

    io_root = pio_comp_settings(shr_pio_getindex(component))%pio_root


  end function shr_pio_getioroot_fromname


!===============================================================================

  !! Given a component name, return the index of that component.
  !! This is the index into io_compid, io_compname, comp_pio_iotype, etc.
  !! If the given component is not found, return -1

  integer function shr_pio_getindex_fromid(compid) result(index)
     implicit none
     integer, intent(in) :: compid
     integer :: i

     index = -1
     do i=1,total_comps
        if(io_compid(i)==compid) then
          index = i
          exit
       end if
    end do

    if(index<0) then
       call shr_sys_abort('shr_pio_getindex :: compid out of allowed range')
    end if
  end function shr_pio_getindex_fromid


  integer function shr_pio_getindex_fromname(component) result(index)
     use shr_string_mod, only : shr_string_toupper

     implicit none

     ! 'component' must be equal to some element of io_compname(:)
     ! (but it is case-insensitive)
     character(len=*), intent(in) :: component

     character(len=len(component)) :: component_ucase
     integer :: i

     ! convert component name to upper case in order to match case in io_compname
     component_ucase = shr_string_toUpper(component)

     index = -1  ! flag for not found
     do i=1,size(io_compname)
        if (trim(component_ucase) == trim(io_compname(i))) then
           index = i
           exit
        end if
     end do
    if(index<0) then
       call shr_sys_abort(' shr_pio_getindex:: compid out of allowed range')
    end if
   end function shr_pio_getindex_fromname

  function shr_pio_getiosys_fromid(compid) result(iosystem)
    ! 'component' must be equal to some element of io_compname(:)
    ! (but it is case-insensitive)
    integer, intent(in) :: compid
    type(iosystem_desc_t), pointer :: iosystem


    iosystem => iosystems(shr_pio_getindex(compid))

  end function shr_pio_getiosys_fromid

  function shr_pio_getiosys_fromname(component) result(iosystem)
    ! 'component' must be equal to some element of io_compname(:)
    ! (but it is case-insensitive)
    character(len=*), intent(in) :: component
    type(iosystem_desc_t), pointer :: iosystem

    iosystem => iosystems(shr_pio_getindex(component))

  end function shr_pio_getiosys_fromname

!===============================================================================



  subroutine shr_pio_read_default_namelist(nlfilename, Comm, pio_stride, pio_root, pio_numiotasks, &
       pio_iotype, pio_async_interface, pio_rearranger)

    character(len=*), intent(in) :: nlfilename
    integer, intent(in) :: Comm
    logical, intent(out) :: pio_async_interface
    integer, intent(out) :: pio_stride, pio_root, pio_numiotasks, pio_iotype, pio_rearranger

    character(len=shr_kind_cs) :: pio_typename
    character(len=shr_kind_cs) :: pio_rearr_comm_type, pio_rearr_comm_fcd
    integer :: pio_rearr_comm_max_pend_req_comp2io
    logical :: pio_rearr_comm_enable_hs_comp2io, pio_rearr_comm_enable_isend_comp2io
    integer :: pio_rearr_comm_max_pend_req_io2comp
    logical :: pio_rearr_comm_enable_hs_io2comp, pio_rearr_comm_enable_isend_io2comp
    character(*),parameter :: subName =   '(shr_pio_read_default_namelist) '

    integer :: iam, ierr, npes, unitn
    logical :: iamroot
    namelist /pio_default_inparm/ pio_stride, pio_root, pio_numiotasks, &
          pio_typename, pio_async_interface, pio_debug_level, pio_blocksize, &
          pio_buffer_size_limit, pio_rearranger, &
          pio_rearr_comm_type, pio_rearr_comm_fcd, &
          pio_rearr_comm_max_pend_req_comp2io, pio_rearr_comm_enable_hs_comp2io, &
          pio_rearr_comm_enable_isend_comp2io, &
          pio_rearr_comm_max_pend_req_io2comp, pio_rearr_comm_enable_hs_io2comp, &
          pio_rearr_comm_enable_isend_io2comp


    call mpi_comm_rank(Comm, iam  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')
    call mpi_comm_size(Comm, npes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    if(iam==0) then
       iamroot=.true.
    else
       iamroot=.false.
    end if

    !--------------------------------------------------------------------------
    ! read io nml parameters
    !--------------------------------------------------------------------------
    pio_stride   = -99 ! set based on pio_numiotasks value when initialized < 0
    pio_numiotasks = -99 ! set based on pio_stride   value when initialized < 0
    pio_root     = -99
    pio_typename = 'nothing'
    pio_blocksize= -99  ! io blocking size set internally in pio when < 0
    pio_buffer_size_limit = -99 ! io task memory buffer maximum set internally in pio when < 0
    pio_debug_level = 0 ! no debug info by default
    pio_async_interface = .false.   ! pio tasks are a subset of component tasks
    pio_rearranger = PIO_REARR_SUBSET

    if(iamroot) then
       unitn=shr_file_getunit()
       open( unitn, file=trim(nlfilename), status='old' , iostat=ierr)
       if(ierr/=0) then
          write(shr_log_unit,*) 'File ',trim(nlfilename),' not found, setting default values.'
       else
          ierr = 1
          do while( ierr /= 0 )
             read(unitn,nml=pio_default_inparm,iostat=ierr)
             if (ierr < 0) then
                call shr_sys_abort( subname//':: namelist read returns an'// &
                     ' end of file or end of record condition '//trim(nlfilename) )
             end if
          end do
          close(unitn)
          call shr_file_freeUnit( unitn )

          call shr_pio_getiotypefromname(pio_typename, pio_iotype, pio_iotype_netcdf)
       end if
    end if

     call shr_pio_namelist_set(npes, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype, &
          iamroot, pio_rearranger)

    call shr_mpi_bcast(pio_debug_level, Comm)
    call shr_mpi_bcast(pio_blocksize, Comm)
    call shr_mpi_bcast(pio_buffer_size_limit, Comm)
    call shr_mpi_bcast(pio_async_interface, Comm)
    call shr_mpi_bcast(pio_rearranger, Comm)

     call shr_pio_rearr_opts_set(Comm, pio_rearr_comm_type, pio_rearr_comm_fcd, &
           pio_rearr_comm_max_pend_req_comp2io, pio_rearr_comm_enable_hs_comp2io, &
           pio_rearr_comm_enable_isend_comp2io, &
           pio_rearr_comm_max_pend_req_io2comp, pio_rearr_comm_enable_hs_io2comp, &
           pio_rearr_comm_enable_isend_io2comp, pio_numiotasks)

  end subroutine shr_pio_read_default_namelist

  subroutine shr_pio_read_component_namelist(nlfilename, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype, pio_rearranger)
    character(len=*), intent(in) :: nlfilename
    integer, intent(in) :: Comm

    integer, intent(inout) :: pio_stride, pio_root, pio_numiotasks, pio_iotype, pio_rearranger
    character(len=SHR_KIND_CS) ::  pio_typename
    integer :: unitn

    integer :: iam, ierr, npes
    logical :: iamroot
    character(*),parameter :: subName =   '(shr_pio_read_component_namelist) '
    integer :: pio_default_stride, pio_default_root, pio_default_numiotasks, pio_default_iotype
    integer :: pio_default_rearranger

    namelist /pio_inparm/ pio_stride, pio_root, pio_numiotasks, &
         pio_typename, pio_rearranger



    call mpi_comm_rank(Comm, iam  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')
    call mpi_comm_size(Comm, npes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    if(iam==0) then
       iamroot=.true.
    else
       iamroot=.false.
    end if

    pio_default_stride = pio_stride
    pio_default_root = pio_root
    pio_default_numiotasks = pio_numiotasks
    pio_default_iotype = pio_iotype
    pio_default_rearranger = pio_rearranger


    !--------------------------------------------------------------------------
    ! read io nml parameters
    !--------------------------------------------------------------------------
    pio_stride   = -99 ! set based on pio_numiotasks value when initialized < 0
    pio_numiotasks = -99 ! set based on pio_stride   value when initialized < 0
    pio_root     = -99
    pio_typename = 'nothing'
    pio_rearranger = -99

    if(iamroot) then
       unitn=shr_file_getunit()
       open( unitn, file=trim(nlfilename), status='old' , iostat=ierr)
       if( ierr /= 0) then
          write(shr_log_unit,*) 'No ',trim(nlfilename),' found, using defaults for pio settings'
           pio_stride     = pio_default_stride
           pio_root       = pio_default_root
           pio_numiotasks = pio_default_numiotasks
           pio_iotype     = pio_default_iotype
           pio_rearranger = pio_default_rearranger
       else
          ierr = 1
          do while( ierr /= 0 )
             read(unitn,nml=pio_inparm,iostat=ierr)
             if (ierr < 0) then
                call shr_sys_abort( subname//':: namelist read returns an'// &
                     ' end of file or end of record condition' )
             end if
          end do
          close(unitn)
          call shr_file_freeUnit( unitn )

          call shr_pio_getiotypefromname(pio_typename, pio_iotype, pio_default_iotype)
       end if
       if(pio_stride== -99) pio_stride = pio_default_stride
       if(pio_root == -99) pio_root = pio_default_root
       if(pio_rearranger == -99) pio_rearranger = pio_default_rearranger
       if(pio_numiotasks == -99) then
          pio_numiotasks = npes/pio_stride
       endif
    endif



    call shr_pio_namelist_set(npes, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype, &
         iamroot, pio_rearranger)


  end subroutine shr_pio_read_component_namelist

  subroutine shr_pio_getiotypefromname(typename, iotype, defaulttype)
    use shr_string_mod, only : shr_string_toupper
    character(len=*), intent(inout) :: typename
    integer, intent(out) :: iotype
    integer, intent(in) :: defaulttype

    typename = shr_string_toupper(typename)
    if      ( typename .eq. 'NETCDF' ) then
       iotype = pio_iotype_netcdf
    else if ( typename .eq. 'PNETCDF') then
       iotype = pio_iotype_pnetcdf
    else if ( typename .eq. 'NETCDF4P') then
       iotype = pio_iotype_netcdf4p
    else if ( typename .eq. 'NETCDF4C') then
       iotype = pio_iotype_netcdf4c
!  Not yet supported
!    else if ( typename .eq. 'VDC') then
!       iotype = pio_iotype_vdc
    else if ( typename .eq. 'NOTHING') then
       iotype = defaulttype
    else if ( typename .eq. 'DEFAULT') then
       iotype = defaulttype
    else
       write(shr_log_unit,*) 'shr_pio_mod: WARNING Bad io_type argument - using iotype_netcdf'
       iotype=pio_iotype_netcdf
    end if

  end subroutine shr_pio_getiotypefromname

!===============================================================================
  subroutine shr_pio_namelist_set(npes,mycomm, pio_stride, pio_root, pio_numiotasks, &
       pio_iotype, iamroot, pio_rearranger)
    integer, intent(in) :: npes, mycomm
    integer, intent(inout) :: pio_stride, pio_root, pio_numiotasks
    integer, intent(inout) :: pio_iotype, pio_rearranger
    logical, intent(in) :: iamroot
    character(*),parameter :: subName =   '(shr_pio_namelist_set) '

    call shr_mpi_bcast(pio_iotype  , mycomm)
    call shr_mpi_bcast(pio_stride  , mycomm)
    call shr_mpi_bcast(pio_root    , mycomm)
    call shr_mpi_bcast(pio_numiotasks, mycomm)
    call shr_mpi_bcast(pio_rearranger, mycomm)

    if (pio_root<0) then
       pio_root = 1
    endif
    pio_root = min(pio_root,npes-1)

! If you are asking for parallel IO then you should use at least two io pes
    if(npes > 1 .and. pio_numiotasks == 1 .and. &
         (pio_iotype .eq. PIO_IOTYPE_PNETCDF .or. &
         pio_iotype .eq. PIO_IOTYPE_NETCDF4P)) then
       pio_numiotasks = 2
    endif

    !--------------------------------------------------------------------------
    ! check/set/correct io pio parameters
    !--------------------------------------------------------------------------
    if (pio_stride>0.and.pio_numiotasks<0) then
       pio_numiotasks = max(1,npes/pio_stride)
    else if(pio_numiotasks>0 .and. pio_stride<0) then
       pio_stride = max(1,npes/pio_numiotasks)
    else if(pio_numiotasks<0 .and. pio_stride<0) then
       pio_stride = max(1,npes/4)
       pio_numiotasks = max(1,npes/pio_stride)
    end if
    if(pio_stride == 1) then
       pio_root = 0
    endif
    if(pio_rearranger .ne. PIO_REARR_SUBSET .and. pio_rearranger .ne. PIO_REARR_BOX) then
       write(shr_log_unit,*) 'pio_rearranger value, ',pio_rearranger,&
            ', not supported - using PIO_REARR_BOX'
       pio_rearranger = PIO_REARR_BOX

    endif


    if (pio_root + (pio_stride)*(pio_numiotasks-1) >= npes .or. &
         pio_stride<=0 .or. pio_numiotasks<=0 .or. pio_root < 0 .or. &
         pio_root > npes-1) then
       if(npes<100) then
          pio_stride = max(1,npes/4)
       else if(npes<1000) then
          pio_stride = max(1,npes/8)
       else
          pio_stride = max(1,npes/16)
       end if
       if(pio_stride>1) then
          pio_numiotasks = npes/pio_stride
          pio_root = min(1,npes-1)
       else
          pio_numiotasks = npes
          pio_root = 0
       end if
       if( iamroot) then
          write(shr_log_unit,*) 'pio_stride, iotasks or root out of bounds - resetting to defaults: ',&
               pio_stride,pio_numiotasks, pio_root
       end if
    end if

  end subroutine shr_pio_namelist_set

  ! This subroutine sets the global PIO rearranger options
  ! The input args that represent the rearranger options are valid only
  ! on the root proc of comm
  ! The rearranger options are passed to PIO_Init() in shr_pio_init2()
  subroutine shr_pio_rearr_opts_set(comm, pio_rearr_comm_type, pio_rearr_comm_fcd, &
          pio_rearr_comm_max_pend_req_comp2io, pio_rearr_comm_enable_hs_comp2io, &
          pio_rearr_comm_enable_isend_comp2io, &
          pio_rearr_comm_max_pend_req_io2comp, pio_rearr_comm_enable_hs_io2comp, &
          pio_rearr_comm_enable_isend_io2comp, &
          pio_numiotasks)
    integer(SHR_KIND_IN), intent(in) :: comm
    character(len=shr_kind_cs), intent(in) :: pio_rearr_comm_type, pio_rearr_comm_fcd
    integer, intent(in) :: pio_rearr_comm_max_pend_req_comp2io
    logical, intent(in) :: pio_rearr_comm_enable_hs_comp2io
    logical, intent(in) :: pio_rearr_comm_enable_isend_comp2io
    integer, intent(in) :: pio_rearr_comm_max_pend_req_io2comp
    logical, intent(in) :: pio_rearr_comm_enable_hs_io2comp
    logical, intent(in) :: pio_rearr_comm_enable_isend_io2comp
    integer, intent(in) :: pio_numiotasks

    character(*), parameter :: subname = '(shr_pio_rearr_opts_set) '
    integer, parameter :: NUM_REARR_COMM_OPTS = 8
    integer, parameter :: PIO_REARR_COMM_DEF_MAX_PEND_REQ = 64
    integer(SHR_KIND_IN), dimension(NUM_REARR_COMM_OPTS) :: buf
    integer :: rank, ierr

    call mpi_comm_rank(comm, rank, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    buf = 0
    ! buf(1) = comm_type
    ! buf(2) = comm_fcd
    ! buf(3) = max_pend_req_comp2io
    ! buf(4) = enable_hs_comp2io
    ! buf(5) = enable_isend_comp2io
    ! buf(6) = max_pend_req_io2comp
    ! buf(7) = enable_hs_io2comp
    ! buf(8) = enable_isend_io2comp
    if(rank == 0) then
      ! buf(1) = comm_type
      select case(pio_rearr_comm_type)
        case ("p2p")
        case ("default")
          buf(1) = pio_rearr_comm_p2p
        case ("coll")
          buf(1) = pio_rearr_comm_coll
        case default
          write(shr_log_unit,*) "Invalid PIO rearranger comm type, ", pio_rearr_comm_type
          write(shr_log_unit,*) "Resetting PIO rearrange comm type to p2p"
          buf(1) = pio_rearr_comm_p2p
      end select

      ! buf(2) = comm_fcd
      select case(pio_rearr_comm_fcd)
        case ("2denable")
        case ("default")
          buf(2) = pio_rearr_comm_fc_2d_enable
        case ("io2comp")
          buf(2) = pio_rearr_comm_fc_1d_io2comp
        case ("comp2io")
          buf(2) = pio_rearr_comm_fc_1d_comp2io
        case ("disable")
          buf(2) = pio_rearr_comm_fc_2d_disable
        case default
          write(shr_log_unit,*) "Invalid PIO rearranger comm flow control direction, ", pio_rearr_comm_fcd
          write(shr_log_unit,*) "Resetting PIO rearrange comm flow control direction to 2denable"
          buf(2) = pio_rearr_comm_fc_2d_enable
      end select

      ! buf(3) = max_pend_req_comp2io
      if((pio_rearr_comm_max_pend_req_comp2io <= 0) .and. &
          (pio_rearr_comm_max_pend_req_comp2io /= PIO_REARR_COMM_UNLIMITED_PEND_REQ)) then

        ! Small multiple of pio_numiotasks has proven to perform
        ! well empirically, and we do not want to allow maximum for
        ! very large process count runs. Can improve this by
        ! communicating between iotasks first, and then non-iotasks
        ! to iotasks (TO DO)
        write(shr_log_unit, *) "Invalid PIO rearranger comm max pend req (comp2io), ",&
             pio_rearr_comm_max_pend_req_comp2io
        write(shr_log_unit, *) "Resetting PIO rearranger comm max pend req (comp2io) to ", &
             max(PIO_REARR_COMM_DEF_MAX_PEND_REQ, 2 * pio_numiotasks)
        buf(3) = max(PIO_REARR_COMM_DEF_MAX_PEND_REQ, 2 * pio_numiotasks)
      else
        buf(3) = pio_rearr_comm_max_pend_req_comp2io
      end if

      ! buf(4) = enable_hs_comp2io
      if(pio_rearr_comm_enable_hs_comp2io) then
        buf(4) = 1
      else
        buf(4) = 0
      end if

      ! buf(5) = enable_isend_comp2io
      if(pio_rearr_comm_enable_isend_comp2io) then
        buf(5) = 1
      else
        buf(5) = 0
      end if

      ! buf(6) = max_pend_req_io2comp
      if((pio_rearr_comm_max_pend_req_io2comp <= 0) .and. &
          (pio_rearr_comm_max_pend_req_io2comp /= PIO_REARR_COMM_UNLIMITED_PEND_REQ)) then
        write(shr_log_unit, *) "Invalid PIO rearranger comm max pend req (io2comp), ", pio_rearr_comm_max_pend_req_io2comp
        write(shr_log_unit, *) "Resetting PIO rearranger comm max pend req (io2comp) to ", PIO_REARR_COMM_DEF_MAX_PEND_REQ
        buf(6) = PIO_REARR_COMM_DEF_MAX_PEND_REQ
      else
        buf(6) = pio_rearr_comm_max_pend_req_io2comp
      end if

      ! buf(7) = enable_hs_io2comp
      if(pio_rearr_comm_enable_hs_io2comp) then
        buf(7) = 1
      else
        buf(7) = 0
      end if

      ! buf(8) = enable_isend_io2comp
      if(pio_rearr_comm_enable_isend_io2comp) then
        buf(8) = 1
      else
        buf(8) = 0
      end if

      ! Log the rearranger options
      write(shr_log_unit, *) "PIO rearranger options:"
      write(shr_log_unit, *) "  comm type     =", pio_rearr_comm_type
      write(shr_log_unit, *) "  comm fcd      =", pio_rearr_comm_fcd
      write(shr_log_unit, *) "  max pend req (comp2io)  =", pio_rearr_comm_max_pend_req_comp2io
      write(shr_log_unit, *) "  enable_hs (comp2io)     =", pio_rearr_comm_enable_hs_comp2io
      write(shr_log_unit, *) "  enable_isend (comp2io)  =", pio_rearr_comm_enable_isend_comp2io
      write(shr_log_unit, *) "  max pend req (io2comp)  =", pio_rearr_comm_max_pend_req_io2comp
      write(shr_log_unit, *) "  enable_hs (io2comp)    =", pio_rearr_comm_enable_hs_io2comp
      write(shr_log_unit, *) "  enable_isend (io2comp)  =", pio_rearr_comm_enable_isend_io2comp
    end if

    call shr_mpi_bcast(buf, comm)

    ! buf(1) = comm_type
    ! buf(2) = comm_fcd
    ! buf(3) = max_pend_req_comp2io
    ! buf(4) = enable_hs_comp2io
    ! buf(5) = enable_isend_comp2io
    ! buf(6) = max_pend_req_io2comp
    ! buf(7) = enable_hs_io2comp
    ! buf(8) = enable_isend_io2comp
    pio_rearr_opts%comm_type = buf(1)
    pio_rearr_opts%fcd = buf(2)
    pio_rearr_opts%comm_fc_opts_comp2io%max_pend_req = buf(3)
    if(buf(4) == 0) then
      pio_rearr_opts%comm_fc_opts_comp2io%enable_hs = .false.
    else
      pio_rearr_opts%comm_fc_opts_comp2io%enable_hs = .true.
    end if
    if(buf(5) == 0) then
      pio_rearr_opts%comm_fc_opts_comp2io%enable_isend = .false.
    else
      pio_rearr_opts%comm_fc_opts_comp2io%enable_isend = .true.
    end if
    pio_rearr_opts%comm_fc_opts_io2comp%max_pend_req = buf(6)
    if(buf(7) == 0) then
      pio_rearr_opts%comm_fc_opts_io2comp%enable_hs = .false.
    else
      pio_rearr_opts%comm_fc_opts_io2comp%enable_hs = .true.
    end if
    if(buf(8) == 0) then
      pio_rearr_opts%comm_fc_opts_io2comp%enable_isend = .false.
    else
      pio_rearr_opts%comm_fc_opts_io2comp%enable_isend = .true.
    end if
  end subroutine
!===============================================================================

end module shr_pio_mod
