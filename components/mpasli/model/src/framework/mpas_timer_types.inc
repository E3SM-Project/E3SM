        type mpas_timer_root
           type (mpas_timer_node), pointer :: root_timer => null()
           type (mpas_timer_node), pointer :: current_timer => null()
           type (dm_info), pointer :: dminfo => null()
        end type mpas_timer_root

        type mpas_timer_node
          character (len=StrKIND) :: timer_name
          logical :: printed
          integer :: nlen
          logical, pointer, dimension(:) :: running
          integer :: calls
          real (kind=R8KIND), pointer, dimension(:) :: start_time, end_time, total_time
          real (kind=RKIND), pointer, dimension(:) :: max_time, min_time, avg_time
          type (mpas_timer_node), pointer :: next => null()
          type (mpas_timer_node), pointer :: child => null()
          type (mpas_timer_node), pointer :: parent => null()
        end type mpas_timer_node
