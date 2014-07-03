module shr_pio_mod
  use pio
  use shr_kind_mod, only : shr_kind_CS, shr_kind_cl
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
  end type pio_comp_t

  character(len=16), allocatable :: io_compname(:)
  type(pio_comp_t), allocatable :: pio_comp_settings(:)
  type (iosystem_desc_t), allocatable, target :: iosystems(:)
  integer :: io_comm
  logical :: pio_async_interface
  integer, allocatable :: io_compid(:)
  integer :: pio_debug_level=0, pio_blocksize=0, pio_buffer_size_limit=-1
  integer :: total_comps=0
contains
!> 
!! @public
!! @brief should be the first routine called after mpi_init. It reads the pio default settings from file drv_in, namelist pio_default_inparm
!! and, if pio_async_interface is true, splits the IO tasks away from the Compute tasks.  It then returns the new compute comm in  
!! Global_Comm and sets module variable io_comm.   
!!
!<
  subroutine shr_pio_init1(ncomps, nlfilename, Global_Comm)
    integer, intent(in) :: ncomps
    character(len=*) :: nlfilename
    integer, intent(inout) :: Global_Comm


    integer :: i, pio_root, pio_stride, pio_numiotasks, pio_iotype
    integer :: mpigrp_world, mpigrp, ierr, mpicom
    character(*),parameter :: subName =   '(shr_pio_init1) '
    integer :: pelist(3,1)

    call shr_pio_read_default_namelist(nlfilename, Global_Comm, pio_stride, pio_root, pio_numiotasks, &
         pio_iotype, pio_async_interface)

    io_comm = MPI_COMM_NULL
    allocate(pio_comp_settings(ncomps))
    do i=1,ncomps
       pio_comp_settings(i)%pio_root = pio_root
       pio_comp_settings(i)%pio_stride = pio_stride
       pio_comp_settings(i)%pio_numiotasks = pio_numiotasks
       pio_comp_settings(i)%pio_iotype = pio_iotype
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
       call pio_init(total_comps,mpi_comm_world, comp_comm, io_comm, iosystems)
       i=1
       if(comp_comm_iam(i)==0) then
          write(shr_log_unit,*) io_compname(i),' : pio_numiotasks = ',pio_comp_settings(i)%pio_numiotasks
          write(shr_log_unit,*) io_compname(i),' : pio_stride = ',pio_comp_settings(i)%pio_stride
          write(shr_log_unit,*) io_compname(i),' : pio_root = ',pio_comp_settings(i)%pio_root
          write(shr_log_unit,*) io_compname(i),' : pio_iotype = ',pio_comp_settings(i)%pio_iotype
       end if

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
                  pio_comp_settings(i)%pio_iotype)
             call pio_init(comp_comm_iam(i), comp_comm(i), pio_comp_settings(i)%pio_numiotasks, 0, &
                  pio_comp_settings(i)%pio_stride, &
                  pio_rearr_box, iosystems(i), &
                  base=pio_comp_settings(i)%pio_root)

             
!                write(shr_log_unit,*) io_compname(i),' : iosys ', iosystems(i)%comp_comm,iosystems(i)%io_comm, comp_comm_iam(i), comp_comm(i), iosystems(i)%iomaster, iosystems(i)%async_interface

                if(iosystems(i)%comp_comm==MPI_COMM_NULL) then
                   call shr_sys_abort(subname//':: '//trim(io_compname(i))//' pio_init failed ')
                end if


             if(comp_comm_iam(i)==0) then
                write(shr_log_unit,*) io_compname(i),' : pio_numiotasks = ',pio_comp_settings(i)%pio_numiotasks
                write(shr_log_unit,*) io_compname(i),' : pio_stride = ',pio_comp_settings(i)%pio_stride
                write(shr_log_unit,*) io_compname(i),' : pio_root = ',pio_comp_settings(i)%pio_root
                write(shr_log_unit,*) io_compname(i),' : pio_iotype = ',pio_comp_settings(i)%pio_iotype
                write(shr_log_unit,*) io_compname(i),' : iosys ', iosystems(i)%comp_comm,iosystems(i)%io_comm
             end if



          end if
       end do
    end if




  end subroutine shr_pio_init2



!===============================================================================
  subroutine shr_pio_finalize(  )
    integer :: ierr
    integer :: i
    
    do i=1,total_comps
       if(iosystems(i)%comp_comm /= MPI_COMM_NULL) then
       ! write(shr_log_unit,*) iosystems(i)%comp_rank, iosystems(i)%io_rank,io_compname(i),' Calling PIO finalize'
          call pio_finalize(iosystems(i), ierr)
       end if
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
       pio_iotype, pio_async_interface)
    
    character(len=*), intent(in) :: nlfilename
    integer, intent(in) :: Comm
    logical, intent(out) :: pio_async_interface
    integer, intent(out) :: pio_stride, pio_root, pio_numiotasks, pio_iotype

    character(len=shr_kind_cs) :: pio_typename
    character(*),parameter :: subName =   '(shr_pio_read_default_namelist) '

    integer :: iam, ierr, npes, unitn
    logical :: iamroot

    namelist /pio_default_inparm/ pio_stride, pio_root, pio_numiotasks, &
         pio_typename, pio_async_interface, pio_debug_level, pio_blocksize, pio_buffer_size_limit



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



    call shr_pio_namelist_set(npes, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype, iamroot)
    call shr_mpi_bcast(pio_debug_level, Comm)
    call shr_mpi_bcast(pio_blocksize, Comm)
    call shr_mpi_bcast(pio_buffer_size_limit, Comm)
    call shr_mpi_bcast(pio_async_interface, Comm)
    

  end subroutine shr_pio_read_default_namelist

  subroutine shr_pio_read_component_namelist(nlfilename, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype)
    character(len=*), intent(in) :: nlfilename
    integer, intent(in) :: Comm

    integer, intent(inout) :: pio_stride, pio_root, pio_numiotasks, pio_iotype
    character(len=SHR_KIND_CS) ::  pio_typename
    integer :: unitn
    
    integer :: iam, ierr, npes
    logical :: iamroot
    character(*),parameter :: subName =   '(shr_pio_read_component_namelist) '
    integer :: pio_default_stride, pio_default_root, pio_default_numiotasks, pio_default_iotype

    namelist /pio_inparm/ pio_stride, pio_root, pio_numiotasks, &
         pio_typename



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


    !--------------------------------------------------------------------------
    ! read io nml parameters 
    !--------------------------------------------------------------------------
    pio_stride   = -99 ! set based on pio_numiotasks value when initialized < 0
    pio_numiotasks = -99 ! set based on pio_stride   value when initialized < 0
    pio_root     = -99
    pio_typename = 'nothing'

    if(iamroot) then
       unitn=shr_file_getunit()
       open( unitn, file=trim(nlfilename), status='old' , iostat=ierr)
       if( ierr /= 0) then
          write(shr_log_unit,*) 'No ',trim(nlfilename),' found, using defaults for pio settings'
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

          if(pio_stride== -99) pio_stride = pio_default_stride
          if(pio_root == -99) pio_root = pio_default_root
          if(pio_numiotasks == -99) then
#if defined(BGP) || defined(BGL)
             if(pio_default_numiotasks < 0 ) then
                pio_numiotasks = pio_default_numiotasks
             else
                pio_numiotasks = min(pio_default_numiotasks,npes/pio_stride)
             end if
#else
             pio_numiotasks = min(pio_default_numiotasks,npes/pio_stride)
#endif
          endif
       endif
    end if

    call shr_pio_namelist_set(npes, Comm, pio_stride, pio_root, pio_numiotasks, pio_iotype, iamroot)


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
  subroutine shr_pio_namelist_set(npes,mycomm, pio_stride, pio_root, pio_numiotasks, pio_iotype, iamroot)
    integer, intent(in) :: npes, mycomm
    integer, intent(inout) :: pio_stride, pio_root, pio_numiotasks
    integer, intent(inout) :: pio_iotype
    logical, intent(in) :: iamroot
    character(*),parameter :: subName =   '(shr_pio_namelist_set) '

    call shr_mpi_bcast(pio_iotype  , mycomm)
    call shr_mpi_bcast(pio_stride  , mycomm)
    call shr_mpi_bcast(pio_root    , mycomm)
    call shr_mpi_bcast(pio_numiotasks, mycomm)

    if (pio_root<0) then
       pio_root = 1
    endif
    pio_root = min(pio_root,npes-1)


#if defined(BGP) || defined(BGL)
! On Bluegene machines a negative numiotasks refers to the number of iotasks per ionode, this code
! allows for that special case.
    if(pio_numiotasks<0) then
       pio_stride=0       
    else
#endif





    !--------------------------------------------------------------------------
    ! check/set/correct io pio parameters
    !--------------------------------------------------------------------------
    if (pio_stride>0.and.pio_numiotasks<0) then
       pio_numiotasks = npes/pio_stride
    else if(pio_numiotasks>0 .and. pio_stride<0) then
       pio_stride = npes/pio_numiotasks
    else if(pio_numiotasks<0 .and. pio_stride<0) then
       pio_stride = 4
       pio_numiotasks = npes/pio_stride
       pio_numiotasks = max(1, pio_numiotasks)
    end if


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
#if defined(BGP) || defined(BGL)
    end if
#endif

  end subroutine shr_pio_namelist_set

!===============================================================================

end module shr_pio_mod
