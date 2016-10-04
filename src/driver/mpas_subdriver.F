! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_subdriver

   use mpas_framework
   use mpas_kind_types

   use mpas_derived_types, only: dm_info, domain_type

#ifdef CORE_ATMOSPHERE
   use atm_core_interface
#endif
#ifdef CORE_CICE
   use cice_core_interface
#endif
#ifdef CORE_INIT_ATMOSPHERE
   use init_atm_core_interface
#endif
#ifdef CORE_LANDICE
   use li_core_interface
#endif
#ifdef CORE_OCEAN
   use ocn_core_interface
#endif
#ifdef CORE_SW
   use sw_core_interface
#endif
#ifdef CORE_TEST
   use test_core_interface
#endif

   type (core_type), pointer :: corelist => null()
   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain_ptr

   contains


   subroutine mpas_init()

      use mpas_stream_manager, only : MPAS_stream_mgr_init, MPAS_build_stream_filename
      use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
      use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string
      use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time
      use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
 
      implicit none

      integer :: iArg, nArgs
      logical :: readNamelistArg, readStreamsArg
      character(len=StrKIND) :: argument, namelistFile, streamsFile
      character(len=StrKIND) :: timeStamp
      integer :: ierr
      integer :: blockID

      character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
      integer(kind=c_int) :: c_comm
      integer(kind=c_int) :: c_ierr
      type (c_ptr) :: mgr_p
      character(len=StrKIND) :: mesh_stream
      character(len=StrKIND) :: mesh_filename
      character(len=StrKIND) :: mesh_filename_temp
      character(len=StrKIND) :: ref_time_temp
      character(len=StrKIND) :: filename_interval_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_iotype
      type (MPAS_Time_type) :: start_time
      type (MPAS_Time_type) :: ref_time
      type (MPAS_TimeInterval_type) :: filename_interval
      character(len=StrKIND) :: start_timestamp
      character(len=StrKIND) :: iotype
      logical :: streamsExists
      integer :: mesh_iotype

      interface
         subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
            use iso_c_binding, only : c_char, c_ptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            type (c_ptr), intent(inout) :: mgr_p
            integer(kind=c_int), intent(inout) :: comm
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_parser

         subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
            use iso_c_binding, only : c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            character(kind=c_char), dimension(*), intent(in) :: streamname
            integer(kind=c_int), intent(inout) :: comm
            character(kind=c_char), dimension(*), intent(out) :: filename
            character(kind=c_char), dimension(*), intent(out) :: ref_time
            character(kind=c_char), dimension(*), intent(out) :: filename_interval
            character(kind=c_char), dimension(*), intent(out) :: io_type
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_get_attributes
      end interface

      readNamelistArg = .false.
      readStreamsArg = .false.

      nArgs = command_argument_count()
      iArg = 1
      do while (iArg < nArgs)
         call get_command_argument(iArg, argument)
         if (len_trim(argument) == 0) exit

         if ( trim(argument) == '-n' ) then
            iArg = iArg + 1
            readNamelistArg = .true.
            call get_command_argument(iArg, namelistFile)
            if ( len_trim(namelistFile) == 0 ) then
                write(0,*) 'ERROR: The -n argument requires a namelist file argument.'
                stop
            else if ( trim(namelistFile) == '-s' ) then
                write(0,*) 'ERROR: The -n argument requires a namelist file argument.'
                stop
            end if
         else if ( trim(argument) == '-s' ) then
            iArg = iArg + 1
            readStreamsArg = .true.
            call get_command_argument(iArg, streamsFile)
            if ( len_trim(streamsFile) == 0 ) then
                write(0,*) 'ERROR: The -s argument requires a streams file argument.'
                stop
            else if ( trim(streamsFile) == '-n' ) then
                write(0,*) 'ERROR: The -s argument requires a streams file argument.'
                stop
            end if
         end if

         iArg = iArg + 1
      end do

      allocate(corelist)
      nullify(corelist % next)

      allocate(corelist % domainlist)
      nullify(corelist % domainlist % next)

      domain_ptr => corelist % domainlist
      domain_ptr % core => corelist

      call mpas_allocate_domain(domain_ptr)

      !
      ! Initialize infrastructure
      !
      call mpas_framework_init_phase1(domain_ptr % dminfo)

#ifdef CORE_ATMOSPHERE
      call atm_setup_core(corelist)
      call atm_setup_domain(domain_ptr)
#endif
#ifdef CORE_CICE
      call cice_setup_core(corelist)
      call cice_setup_domain(domain_ptr)
#endif
#ifdef CORE_INIT_ATMOSPHERE
      call init_atm_setup_core(corelist)
      call init_atm_setup_domain(domain_ptr)
#endif
#ifdef CORE_LANDICE
      call li_setup_core(corelist)
      call li_setup_domain(domain_ptr)
#endif
#ifdef CORE_OCEAN
      call ocn_setup_core(corelist)
      call ocn_setup_domain(domain_ptr)
#endif
#ifdef CORE_SW
      call sw_setup_core(corelist)
      call sw_setup_domain(domain_ptr)
#endif
#ifdef CORE_TEST
      call test_setup_core(corelist)
      call test_setup_domain(domain_ptr)
#endif

      if ( readNamelistArg ) then
         domain_ptr % namelist_filename = namelistFile
      end if

      if ( readStreamsArg ) then
         domain_ptr % streams_filename = streamsFile
      end if

      ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Namelist setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      call mpas_framework_init_phase2(domain_ptr)

      ierr = domain_ptr % core % define_packages(domain_ptr % packages)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Package definition failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % iocontext)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Package setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Decomposition setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Clock setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      write(stderrUnit,*) 'Reading streams configuration from file '//trim(domain_ptr % streams_filename)
      inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

      if ( .not. streamsExists ) then
         call mpas_dmpar_global_abort('ERROR: Streams file '//trim(domain_ptr % streams_filename)//' does not exist.')
      end if

      call mpas_timer_start('total time')
      call mpas_timer_start('initialize')

      !
      ! Using information from the namelist, a graph.info file, and a file containing
      !    mesh fields, build halos and allocate blocks in the domain
      !
      ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Failed to find mesh stream for core ' // trim(domain_ptr % core % coreName))
      end if

      call mpas_f_to_c_string(domain_ptr % streams_filename, c_filename)
      call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
      c_comm = domain_ptr % dminfo % comm
      call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                     c_mesh_filename_temp, c_ref_time_temp, &
                                     c_filename_interval_temp, c_iotype, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(domain_ptr % dminfo)
      end if
      call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
      call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
      call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
      call mpas_c_to_f_string(c_iotype, iotype)

      if (trim(iotype) == 'pnetcdf') then
         mesh_iotype = MPAS_IO_PNETCDF
      else if (trim(iotype) == 'pnetcdf,cdf5') then
         mesh_iotype = MPAS_IO_PNETCDF5
      else if (trim(iotype) == 'netcdf') then
         mesh_iotype = MPAS_IO_NETCDF
      else if (trim(iotype) == 'netcdf4') then
         mesh_iotype = MPAS_IO_NETCDF4
      else
         mesh_iotype = MPAS_IO_PNETCDF
      end if

      start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
      if ( trim(ref_time_temp) == 'initial_time' ) then
          call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
      end if

      blockID = -1
      if ( trim(filename_interval_temp) == 'none' ) then
          call mpas_expand_string(ref_time_temp, blockID, mesh_filename_temp, mesh_filename)
      else
          call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
          call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
          call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, blockID, mesh_filename, ierr)
      end if
      write(stderrUnit, *) ' ** Attempting to bootstrap MPAS framework using stream: ', trim(mesh_stream)
      call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

      !
      ! Set up run-time streams
      !
      call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % ioContext, domain_ptr % clock, &
                                domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

      call add_stream_attributes(domain_ptr)

      ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Immutable streams setup failed for core ' // trim(domain_ptr % core % coreName))
      end if

      mgr_p = c_loc(domain_ptr % streamManager)
      call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
      if (c_ierr /= 0) then
         call mpas_dmpar_abort(domain_ptr % dminfo)
      end if

      !
      ! Finalize the setup of blocks and fields
      !
      call mpas_bootstrap_framework_phase2(domain_ptr)

      !
      ! Initialize core
      !
      iErr = domain_ptr % core % core_init(domain_ptr, timeStamp)
      if ( ierr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Core init failed for core ' // trim(domain_ptr % core % coreName))
      end if

      call mpas_timer_stop('initialize')

   end subroutine mpas_init


   subroutine mpas_run()

      implicit none

      integer :: iErr

      iErr = domain_ptr % core % core_run(domain_ptr)
      if ( iErr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Core run failed for core ' // trim(domain_ptr % core % coreName))
      end if

   end subroutine mpas_run


   subroutine mpas_finalize()

      use mpas_stream_manager, only : MPAS_stream_mgr_finalize
   
      implicit none

      integer :: iErr


      !
      ! Finalize core
      !
      iErr = domain_ptr % core % core_finalize(domain_ptr)
      if ( iErr /= 0 ) then
         call mpas_dmpar_global_abort('ERROR: Core finalize failed for core ' // trim(domain_ptr % core % coreName))
      end if

      call mpas_timer_stop('total time')
      call mpas_timer_write_header()
      call mpas_timer_write()
      call mpas_timer_finalize(domain_ptr)


      !
      ! Finalize infrastructure
      !
      call MPAS_stream_mgr_finalize(domain_ptr % streamManager)

      call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr)

      deallocate(corelist % domainlist)
      deallocate(corelist)

   end subroutine mpas_finalize


   subroutine add_stream_attributes(domain)

      use mpas_stream_manager, only : MPAS_stream_mgr_add_att

      implicit none

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt

      integer :: local_ierr

      if (domain % dminfo % nProcs < 10) then
          write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100) then
          write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 1000) then
          write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 10000) then
          write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100000) then
          write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else
          write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      end if
     
      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source', domain % core % source)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'is_periodic', domain % is_periodic)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'x_period', domain % x_period)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'y_period', domain % y_period)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain %  parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes

end module mpas_subdriver
