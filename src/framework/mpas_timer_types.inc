        type timer_node
          character (len=StrKIND) :: timer_name
          logical :: running, printable
          integer :: levels, calls, nlen
          real (kind=R8KIND) :: start_time, end_time, total_time
          real (kind=RKIND) :: max_time, min_time, avg_time
          real (kind=RKIND) :: efficiency
          type (timer_node), pointer :: next
        end type timer_node

