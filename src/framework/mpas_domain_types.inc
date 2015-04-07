
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

      type (domain_type), pointer :: domain

      type (parallel_info), pointer :: parinfo

      type (block_type), pointer :: prev => null()
      type (block_type), pointer :: next => null()

      type (mpas_pool_type), pointer :: structs, dimensions, configs, packages
      type (mpas_pool_type), pointer :: allFields, allStructs
   end type block_type


   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
      type (mpas_pool_type), pointer :: configs, packages

      type (MPAS_Clock_type), pointer :: clock
	  type (MPAS_streamManager_type), pointer :: streamManager
	  type (mpas_decomp_list), pointer :: decompositions
   
      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo
#include "model_variables.inc"
      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius
      character (len=StrKIND*2) :: history !< History attribute, read in from input file.
      character (len=StrKIND) :: mesh_id !< mesh_id attribute, randomly generated
      character (len=StrKIND) :: Conventions !< Conventions attribute, read in from input file.
      character (len=StrKIND) :: source !< source attribute, read in from input file.
      character (len=StrKIND) :: mesh_spec !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id !< parent_id attribute, read in from input file.
   end type domain_type
