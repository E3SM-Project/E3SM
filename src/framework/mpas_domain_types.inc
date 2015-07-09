   ! Derived type for storing list of blocks from a domain to be handled by a process
   type domain_type
      type (block_type), pointer :: blocklist
      type (mpas_pool_type), pointer :: configs, packages

      type (MPAS_Clock_type), pointer :: clock
      type (MPAS_streamManager_type), pointer :: streamManager
      type (mpas_decomp_list), pointer :: decompositions => null()

      ! Also store parallelization info here
      type (dm_info), pointer :: dminfo

      ! Domain specific constants
      logical :: on_a_sphere
      real (kind=RKIND) :: sphere_radius
      character (len=StrKIND) :: namelist_filename !< Constant: Name of namelist file
      character (len=StrKIND) :: streams_filename !< Constant: Name of stream configuration file
      character (len=StrKIND) :: mesh_spec !< mesh_spec attribute, read in from input file.
      character (len=StrKIND) :: parent_id !< parent_id attribute, read in from input file.

      ! Back pointer to core
      type (core_type), pointer :: core => null()

      ! Domain_type is a linked list
      type (domain_type), pointer :: next => null()
   end type domain_type
