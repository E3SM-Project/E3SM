   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: using_external_comm

     ! Add variables specific to block decomposition. {{{
     ! These are used in mpas_block_decomp.F
     integer :: total_blocks
     logical :: explicitDecomp
     integer, dimension(:), allocatable :: block_proc_list
     integer, dimension(:), allocatable :: block_local_id_list
     !}}}
   end type dm_info


   type mpas_exchange_list
     integer :: endPointID
     integer :: nlist
     integer, dimension(:), pointer :: srcList
     integer, dimension(:), pointer :: destList
     type (mpas_exchange_list), pointer :: next
   end type mpas_exchange_list


   type mpas_exchange_list_pointer
     type (mpas_exchange_list), pointer :: exchList
   end type mpas_exchange_list_pointer


   type mpas_multihalo_exchange_list
     type (mpas_exchange_list_pointer), dimension(:), pointer :: halos
     ! Pointers to the mulithalo exchange lists for this variable on the prev and next blocks on this processor
     type (mpas_multihalo_exchange_list), pointer :: prev, next
   end type mpas_multihalo_exchange_list


   type mpas_communication_list
     integer :: procID
     integer :: nlist
     real (kind=RKIND), dimension(:), pointer :: rbuffer
     integer, dimension(:), pointer :: ibuffer
     integer :: reqID
     type (mpas_communication_list), pointer :: next
   end type mpas_communication_list
