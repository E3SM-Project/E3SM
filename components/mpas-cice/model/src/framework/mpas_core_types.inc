   abstract interface
      function mpas_setup_namelist_function(configs, namelistFilename, dminfo) result(iErr)
         use mpas_kind_types
         import mpas_pool_type
         import dm_info

         type (mpas_pool_type), intent(inout) :: configs
         character (len=*), intent(in) :: namelistFilename
         type (dm_info), intent(in) :: dminfo
         integer :: iErr
      end function mpas_setup_namelist_function
   end interface

   abstract interface
      function mpas_define_packages_function(packages) result(iErr)
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: packages
         integer :: iErr
      end function mpas_define_packages_function
   end interface

   abstract interface
      function mpas_setup_packages_function(configs, packages, iocontext) result(iErr)
         import mpas_pool_type
		 import mpas_io_context_type

         type (mpas_pool_type), intent(inout) :: configs
         type (mpas_pool_type), intent(inout) :: packages
         type (mpas_io_context_type), intent(inout) :: iocontext
         integer :: iErr
      end function mpas_setup_packages_function
   end interface

   abstract interface
      function mpas_setup_decompositions_function(decompList) result(iErr)
         import mpas_decomp_list

         type (mpas_decomp_list), pointer :: decompList
         integer :: iErr
      end function mpas_setup_decompositions_function
   end interface

   abstract interface
      function mpas_get_mesh_stream_function(configs, stream) result(iErr)
         use mpas_kind_types
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: configs
         character (len=StrKIND), intent(out) :: stream
         integer :: iErr
      end function mpas_get_mesh_stream_function
   end interface

   abstract interface
      function mpas_setup_clock_function(clock, configs) result(iErr)
         import mpas_clock_type
         import mpas_pool_type

         type (mpas_clock_type), intent(inout) :: clock
         type (mpas_pool_type), intent(inout) :: configs
         integer :: iErr
      end function mpas_setup_clock_function
   end interface

   abstract interface
      function mpas_setup_log_function(logInfo, domain) result(iErr)
         import mpas_log_type
         import domain_type

         type (mpas_log_type), pointer, intent(inout) :: logInfo
         type (domain_type), pointer, intent(in) :: domain
         integer :: iErr
      end function mpas_setup_log_function
   end interface

   abstract interface
      function mpas_setup_immutable_streams_function(manager) result(iErr)
         import mpas_streamManager_type

         type (mpas_streamManager_type), pointer :: manager
         integer :: iErr
      end function mpas_setup_immutable_streams_function
   end interface

   abstract interface
      function mpas_setup_block_function(block) result(iErr)
         import block_type

         type (block_type), pointer :: block
         integer :: iErr
      end function mpas_setup_block_function
   end interface

   abstract interface
      function mpas_setup_derived_dimensions_function(readDimensions, dimensionPool, configPool) result(iErr)
         import mpas_pool_type

         type (mpas_pool_type), intent(inout) :: readDimensions
         type (mpas_pool_type), intent(inout) :: dimensionPool
         type (mpas_pool_type), intent(inout) :: configPool
         integer :: iErr
      end function mpas_setup_derived_dimensions_function
   end interface

   abstract interface
      function mpas_setup_decomposed_dimensions_function(block, streamManager, readDimensions, dimensionPool, totalBlocks) result(iErr)
		 import block_type
		 import mpas_streamManager_type
         import mpas_pool_type

		 type (block_type), intent(inout) :: block
		 type (mpas_streamManager_type), intent(inout) :: streamManager
         type (mpas_pool_type), intent(inout) :: readDimensions
         type (mpas_pool_type), intent(inout) :: dimensionPool
		 integer, intent(in) :: totalBlocks
         integer :: iErr
      end function mpas_setup_decomposed_dimensions_function
   end interface

   abstract interface
      function mpas_core_init_function(domain, timeStamp) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         character (len=*), intent(out) :: timeStamp
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_core_run_function(domain) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         integer :: iErr
      end function
   end interface

   abstract interface
      function mpas_core_finalize_function(domain) result(iErr)
         import domain_type
         type (domain_type), intent(inout) :: domain
         integer :: iErr
      end function
   end interface

   type core_type
      type (domain_type), pointer :: domainlist => null()

      character (len=StrKIND) :: modelName !< Constant: Name of model
      character (len=StrKIND) :: coreName !< Constant: Name of core
      character (len=StrKIND) :: modelVersion !< Constant: Version number
      character (len=StrKIND) :: executableName !< Constant: Name of executable generated at build time.
      character (len=StrKIND) :: git_version !< Constant: Version string from git-describe.
      character (len=StrKIND*2) :: history !< History attribute, read in from input file.
      character (len=StrKIND) :: Conventions !< Conventions attribute, read in from input file.
      character (len=StrKIND) :: source !< source attribute, read in from input file.

      ! Core init, run, and finalize function pointers
      procedure (mpas_core_init_function), pointer, nopass :: core_init => null()
      procedure (mpas_core_run_function), pointer, nopass :: core_run => null()
      procedure (mpas_core_finalize_function), pointer, nopass :: core_finalize => null()

      ! Core framework function pointers
      procedure (mpas_setup_namelist_function), pointer, nopass :: setup_namelist => null()
      procedure (mpas_define_packages_function), pointer, nopass :: define_packages => null()
      procedure (mpas_setup_packages_function), pointer, nopass :: setup_packages => null()
      procedure (mpas_setup_decompositions_function), pointer, nopass :: setup_decompositions => null()
      procedure (mpas_get_mesh_stream_function), pointer, nopass :: get_mesh_stream => null()
      procedure (mpas_setup_clock_function), pointer, nopass :: setup_clock => null()
      procedure (mpas_setup_log_function), pointer, nopass :: setup_log => null()
      procedure (mpas_setup_block_function), pointer, nopass :: setup_block => null()
      procedure (mpas_setup_immutable_streams_function), pointer, nopass :: setup_immutable_streams => null()
      procedure (mpas_setup_derived_dimensions_function), pointer, nopass :: setup_derived_dimensions => null()
      procedure (mpas_setup_decomposed_dimensions_function), pointer, nopass :: setup_decomposed_dimensions => null()

      ! core_type is a linked list
      type (core_type), pointer :: next => null()
   end type core_type
