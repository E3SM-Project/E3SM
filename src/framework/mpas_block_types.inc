
   ! Type for storing (possibly architecture specific) information concerning to parallelism
   type parallel_info
      type (mpas_multihalo_exchange_list), pointer :: cellsToSend            ! List of types describing which cells to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToRecv            ! List of types describing which cells to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: cellsToCopy            ! List of types describing which cells to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: edgesToSend            ! List of types describing which edges to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToRecv            ! List of types describing which edges to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: edgesToCopy            ! List of types describing which edges to copy from other blocks

      type (mpas_multihalo_exchange_list), pointer :: verticesToSend         ! List of types describing which vertices to send to other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToRecv         ! List of types describing which vertices to receive from other blocks
      type (mpas_multihalo_exchange_list), pointer :: verticesToCopy         ! List of types describing which vertices to copy from other blocks
   end type parallel_info


   ! Derived type for storing part of a domain; used as a basic unit of work for a process
   type block_type

      integer :: blockID   ! Unique global ID number for this block
      integer :: localBlockID  ! Unique local ID number for this block

      ! Pointer to domain that owns the block
      type (domain_type), pointer :: domain

      ! Data structures for exchange lists for a block (exchanging standard fields)
      type (parallel_info), pointer :: parinfo

      ! Loop bounds for threaded loops
      integer :: nThreads = 0
      integer, dimension(:), pointer :: cellThreadStart => null()
      integer, dimension(:), pointer :: cellThreadEnd => null()
      integer, dimension(:), pointer :: cellSolveThreadStart => null()
      integer, dimension(:), pointer :: cellSolveThreadEnd => null()
      integer, dimension(:), pointer :: edgeThreadStart => null()
      integer, dimension(:), pointer :: edgeThreadEnd => null()
      integer, dimension(:), pointer :: edgeSolveThreadStart => null()
      integer, dimension(:), pointer :: edgeSolveThreadEnd => null()
      integer, dimension(:), pointer :: vertexThreadStart => null()
      integer, dimension(:), pointer :: vertexThreadEnd => null()
      integer, dimension(:), pointer :: vertexSolveThreadStart => null()
      integer, dimension(:), pointer :: vertexSolveThreadEnd => null()

      ! Linked list pointers
      type (block_type), pointer :: prev => null()
      type (block_type), pointer :: next => null()

      ! Data structures for core specific data
      type (mpas_pool_type), pointer :: structs, dimensions, configs, packages

      ! Data structures for IO infrastructure
      type (mpas_pool_type), pointer :: allFields, allStructs

      ! Data types for particle lists {{{

      ! Define the particle list on a block
      type (mpas_particle_list_type), pointer :: particlelist => null()

      ! Define the list of neighboring blocks on a block, for exchanging particles
      integer, dimension(:), pointer :: blockNeighs => null()
      ! Define the list of neighboring processors on a block, for exchanging particles
      integer, dimension(:), pointer :: procNeighs => null()

      !}}}

   end type block_type
