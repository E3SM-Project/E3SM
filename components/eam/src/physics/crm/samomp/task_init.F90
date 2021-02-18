module task_init_mod
  implicit none

contains

  subroutine task_init

    !   Check things, initialize multitasking:

    use grid
    implicit none

    integer itasks,ntasks

    if(YES3D .ne. 1 .and. YES3D .ne. 0) then
      print*,'YES3D is not 1 or 0. STOP'
      stop
    endif

    if(YES3D .eq. 1 .and. ny_gl .lt. 4) then
      print*,'ny_gl is too small for a 3D case.STOP'
      stop
    endif

    if(YES3D .eq. 0 .and. ny_gl .ne. 1) then
      print*,'ny_gl should be 1 for a 2D case. STOP'
      stop
    endif

    if(nsubdomains.eq.1) then

      rank =0
      ntasks = 1
      dompi = .false.

    else

      !  call task_start(rank, ntasks)

      !  dompi = .true.

      !  call systemf('hostname')

      !  if(ntasks.ne.nsubdomains) then
      !    if(masterproc) print *,'number of processors is not equal to nsubdomains!',&
      !             '  ntasks=',ntasks,'   nsubdomains=',nsubdomains
      !    call task_abort()
      !  endif

      !  call task_barrier()

      !  call task_ranks()

    end if ! nsubdomains.eq.1

    masterproc = rank.eq.0

  end

end module task_init_mod
