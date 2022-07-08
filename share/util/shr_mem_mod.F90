MODULE shr_mem_mod

  use shr_kind_mod, only : shr_kind_r8
  use shr_log_mod, only: s_logunit => shr_log_Unit
  use shr_sys_mod, only: shr_sys_abort

  implicit none
  private

  ! PUBLIC: Public interfaces

  public ::  shr_mem_getusage, &
       shr_mem_init

  ! PUBLIC: Public interfaces

  real(shr_kind_r8) :: mb_blk = 0.0_shr_kind_r8

  !===============================================================================
CONTAINS
  !===============================================================================

  subroutine shr_mem_init(prt, strbuf)

    implicit none

    !----- arguments -----

    logical, optional :: prt
    character(len=*), optional :: strbuf
    !----- local -----

    ! --- Memory stats ---
    integer :: msize                   ! memory size (high water)
    integer :: mrss0,mrss1,mrss2       ! temporary rss
    integer :: mshare,mtext,mdatastack
    logical :: lprt
    integer :: ierr

    integer :: GPTLget_memusage

    real(shr_kind_r8),allocatable :: mem_tmp(:)

    character(*),parameter :: subname = "(shr_mem_init)"
    !---------------------------------------------------

    lprt = .false.
    if (present(prt)) then
       lprt = prt
    endif

    ierr = GPTLget_memusage (msize, mrss0, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss0 failed')

    allocate(mem_tmp(1024*1024), stat=ierr)    ! 1 MWord, 8 MB
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': allocate failed')

    mem_tmp = -1.0
    ierr = GPTLget_memusage (msize, mrss1, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss1 failed')

    deallocate(mem_tmp, stat=ierr)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': deallocate failed')

    ierr = GPTLget_memusage (msize, mrss2, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss2 failed')

    mb_blk = 0.0_shr_kind_r8
    if (mrss1 - mrss0 > 0) then
       mb_blk = (8.0_shr_kind_r8)/((mrss1-mrss0)*1.0_shr_kind_r8)
    endif

    if (lprt) then
       write(s_logunit,'(A,f16.2)') '8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
       write(s_logunit,'(A,f16.2)') '8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
       write(s_logunit,'(A,f16.2)') 'Memory block size conversion in bytes is ',mb_blk*1024.0_shr_kind_r8*1024.0_shr_kind_r8
    endif
    if (present(strbuf)) then
       write(strbuf,'(3(A,f16.2))') '8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk, &
            '\n8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk, &
            '\nMemory block size conversion in bytes is ',mb_blk*1024.0_shr_kind_r8*1024.0_shr_kind_r8
    endif


  end subroutine shr_mem_init

  !===============================================================================

  subroutine shr_mem_getusage(r_msize,r_mrss,prt)

    implicit none

    !----- arguments ---
    real(shr_kind_r8) :: r_msize,r_mrss
    logical, optional :: prt

    !----- local ---
    integer :: msize,mrss
    integer :: mshare,mtext,mdatastack
    integer :: ierr
    integer :: GPTLget_memusage, GPTLprint_memusage

    !---------------------------------------------------

    ierr = GPTLget_memusage (msize, mrss, mshare, mtext, mdatastack)
    r_msize = msize / 1024.0_shr_kind_r8
    r_mrss  = mrss  / 1024.0_shr_kind_r8

    if (present(prt)) then
       if (prt) then
          ierr = GPTLprint_memusage(' ')
       endif
    endif


  end subroutine shr_mem_getusage

  !===============================================================================

END MODULE shr_mem_mod
