module task_util_mod
  use params, only: crm_rknd
  implicit none

contains

  subroutine task_start(rank,numtasks)
    integer rank,numtasks
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_abort()
    print*,'Aborting the program...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_finish()
    print*,'program is finished...'
    stop
  end

  !----------------------------------------------------------------------
  subroutine task_barrier()
    return
  end

  !----------------------------------------------------------------------

  subroutine task_bcast_float(rank_from,buffer,length)
    implicit none
    integer rank_from       ! broadcasting task's rank
    real(crm_rknd) buffer(*)          ! buffer of data
    integer length          ! buffers' length
    print*, 'MPIsndf call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_send_float(rank_to,buffer,length,tag,request)
    implicit none
    integer rank_to     ! receiving task's rank
    real(crm_rknd) buffer(*)      ! buffer of data
    integer length      ! buffers' length
    integer tag     ! tag of the message
    integer request     ! request id
    print*, 'MPIsndf call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_send_integer(rank_to,buffer,length,tag,request)
    implicit none
    integer rank_to     ! receiving task's rank
    integer buffer(*)   ! buffer of data
    integer length      ! buffers' length
    integer tag     ! tag of the message
    integer request
    print*, 'MPIsndi call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_send_character(rank_to,buffer,length,tag,request)
    implicit none
    integer rank_to     ! receiving task's rank
    character*1 buffer(*)   ! buffer of data
    integer length      ! buffers' length
    integer tag     ! tag of the message
    integer request
    print*, 'MPIsndi call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_receive_float(buffer,length,request)
    real(crm_rknd) buffer(*)      ! buffer of data
    integer length      ! buffers' length
    integer request
    print*, 'MPIrcvf call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_receive_charcater(buffer,length,request)
    character*1 buffer(*)   ! buffer of data
    integer length      ! buffers' length
    integer request
    print*, 'MPIrcvi call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_receive_integer(buffer,length,request)
    integer buffer(*)   ! buffer of data
    integer length      ! buffers' length
    integer request
    print*, 'MPIrcvi call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_bsend_float(rank_to,buffer,length,tag)
    integer rank_to         ! receiving task's rank
    real(crm_rknd) buffer(*)          ! buffer of data
    integer length          ! buffers' length
    integer tag             ! tag of the message
    print*, 'MPI call from a single task program! Exiting...'
    stop
    return
  end

  !----------------------------------------------------------------------

  subroutine task_wait(request,rank,tag)
    integer request
    integer rank, tag
    return
  end

  !----------------------------------------------------------------------

  subroutine task_waitall(count,reqs,ranks,tags)
    integer count,reqs(count)
    integer ranks(count),tags(count)
    return
  end

  !----------------------------------------------------------------------

  subroutine task_test(request,flag,rank,tag)
    integer request
    integer rank, tag
    logical flag
    print*, 'MPItst call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_sum_real(buffer1,buffer2,length)
    real(crm_rknd) buffer1(*) ! buffer of data
    real(crm_rknd) buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end

  !----------------------------------------------------------------------

  subroutine task_sum_real8(buffer1,buffer2,length)
    real(8) buffer1(*) ! buffer of data
    real(8) buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_sum_integer(buffer1,buffer2,length)
    integer buffer1(*) ! buffer of data
    integer buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_max_real(buffer1,buffer2,length)
    real(crm_rknd) buffer1(*) ! buffer of data
    real(crm_rknd) buffer2(*) ! buffer of data
    integer length      ! buffers' length
    return
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_max_integer(buffer1,buffer2,length)
    integer buffer1(*) ! buffer of data
    integer buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_min_real(buffer1,buffer2,length)
    real(crm_rknd) buffer1(*) ! buffer of data
    real(crm_rknd) buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_min_integer(buffer1,buffer2,length)
    integer buffer1(*) ! buffer of data
    integer buffer2(*) ! buffer of data
    integer length      ! buffers' length
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------

  subroutine task_receive_character(buffer,length,request)
    character*1 buffer(*)   ! buffer of data
    integer length          ! buffers' length
    integer request
    print*, 'MPI call from a single task program! Exiting...'
    stop
  end
  !----------------------------------------------------------------------
  subroutine task_rank_to_index (rank,i,j)
    integer rank, i, j
    i=0
    j=0
  end
  !----------------------------------------------------------------------
  subroutine task_bound_duvdt ()
    return
  end
  !----------------------------------------------------------------------
  subroutine task_boundaries(flag)
    integer flag
  end


end module task_util_mod
