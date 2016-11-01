   integer, parameter :: MPAS_DMPAR_NOERR = 0
   integer, parameter :: MPAS_DMPAR_MISSING_GROUP = 1
   integer, parameter :: MPAS_DMPAR_EXISTING_GROUP = 2
   integer, parameter :: MPAS_DMPAR_MISSING_FIELD = 3
   integer, parameter :: MPAS_DMPAR_FIELD_TIMELEVEL_ERR = 4
   integer, parameter :: MPAS_DMPAR_FIELD_HALO_ERR = 5
   integer, parameter :: MPAS_DMPAR_BUFFER_EXISTS = 6

   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: initialized_mpi

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
     integer :: bufferOffset
     real (kind=RKIND), dimension(:), pointer :: rbuffer => null()
     integer, dimension(:), pointer :: ibuffer => null()
     integer :: reqID
     type (mpas_communication_list), pointer :: next => null()
     logical :: received
   end type mpas_communication_list

   type mpas_exchange_field_list
      integer :: nLen
      character (len=StrKIND) :: fieldName
      integer :: fieldType, nDims
      logical, dimension(:), pointer :: haloLayers => null()
      logical, dimension(:), pointer :: timeLevels => null()
      type (mpas_exchange_field_list), pointer :: next => null()
   end type mpas_exchange_field_list

   type mpas_exchange_group
      integer :: nLen
      character (len=StrKIND) :: groupName
      type (mpas_pool_type), pointer :: fieldPool => null()
      type (mpas_exchange_field_list), pointer :: fieldList => null()

      ! Each group should have it's own communication lists, since they have the buffers for the communications
      type (mpas_communication_list), pointer :: sendList => null()
      type (mpas_communication_list), pointer :: recvList => null()
      type (mpas_exchange_group), pointer :: next => null()
   end type mpas_exchange_group
