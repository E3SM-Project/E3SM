!=======================================================================
!
!BOP
!
! !MODULE: ice_restart - ice model restart files
!
! !DESCRIPTION:
!
! Read and write ice model restart files
!
! !REVISION HISTORY:
!  SVN:$Id: ice_restart.F90 53 2007-02-08 00:02:16Z dbailey $
!
! authors Elizabeth C. Hunke, LANL
!         William H. Lipscomb LANL
!
! 2004-05: Block structure added by William Lipscomb
!          Restart module separated from history module
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Converted to free source form (F90) 
! 2008 ECH: Rearranged order in which internal stresses are written and read
! 
! !INTERFACE:
!
      module ice_restart
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task
      use ice_blocks
      use ice_read_write
      use ice_fileunits
      use ice_timers
      use ice_exit, only: abort_ice	
!
!EOP
!
      implicit none
      save

      character(len=char_len_long) :: &
         ice_ic      ! method of ice cover initialization
                     ! 'default'  => latitude and sst dependent
                     ! 'none'     => no ice
                     ! note:  restart = .true. overwrites

      logical (kind=log_kind) :: &
         restart ! ONLY USED if CCSMCOUPLED is  not defined
                 ! if true, initialize using restart file instead of defaults
                 ! ice_forcing uses this variable, so cannot use a ccp if-def here
                 
      character (len=char_len) :: &       	 
         runtype           ! ONLY USED if CCSMCOUPLED is defined
                           ! initial, continue, branch or hybrid 
                           ! branch/hybrid applies to ccsm concurrent mode
                           ! branch applies to ccsm sequential model
                           ! branch or hybrid do not apply to stand-alone cice

      character (len=char_len_long) :: &
         restart_file  , & ! output file prefix for restart dump
         restart_dir   , & ! directory name for restart dump
         runid             ! identifier for CCSM coupled run

      character (len=char_len) :: &
         restart_format, & ! format of restart files 'nc' or 'bin'
         restart_format_in ! format of restart files 'nc' or 'bin'

      character (len=char_len_long) :: &
         pointer_file      ! input pointer file for restarts

      character (len=char_len) :: &       	 
         resttype           ! type of restart format ('new' or 'old')

      real (kind=dbl_kind), private, &
         dimension(nx_block,ny_block,max_blocks) :: &
         work1

      integer (kind=int_kind) :: ncid ! netcdf restart file id

      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine

      character (len=1) :: nchar

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine dumpfile(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for a restart
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_flux
      use ice_grid
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_dyn_evp
      use ice_ocean, only: oceanmixed_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: dimid_ni, dimid_nj, dimid_ncat, &
                                 dimid_ntilyr, dimid_ntslyr
      integer (kind=int_kind), allocatable :: dims(:)

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         if (restart_format == 'nc') then
            write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a)') &
                 restart_dir(1:lenstr(restart_dir)), &
                 restart_file(1:lenstr(restart_file)),'.', &
                 iyear,'-',month,'-',mday,'-',sec,'.nc'
         else
            write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
                 restart_dir(1:lenstr(restart_dir)), &
                 restart_file(1:lenstr(restart_file)),'.', &
                 iyear,'-',month,'-',mday,'-',sec
         endif
      end if
         
      if (restart_format == 'nc') filename = trim(filename) // '.nc'

      ! write pointer (path/file)
      if (my_task == master_task) then
        open(nu_rst_pointer,file=pointer_file)
        write(nu_rst_pointer,'(a)') filename
        close(nu_rst_pointer)
      endif

      ! begin writing restart data
      if (restart_format == "nc") then

         if (my_task == master_task) then
#ifdef _HIRES
            status = nf90_create(filename, NF90_64BIT_OFFSET, ncid)
#else
            status = nf90_create(filename, nf90_clobber, ncid)
#endif
            if (status /= nf90_noerr) call abort_ice( &
               'ice: Error creating restart file '//filename)

            status = nf90_put_att(ncid,nf90_global,'istep1',istep1)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute istep1')
            status = nf90_put_att(ncid,nf90_global,'time',time)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute time')
            status = nf90_put_att(ncid,nf90_global,'time_forc',time_forc)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute time_forc')

            status = nf90_def_dim(ncid,'ni',nx_global,dimid_ni)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining dim ni')

            status = nf90_def_dim(ncid,'nj',ny_global,dimid_nj)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining dim nj')

            status = nf90_def_dim(ncid,'ncat',ncat,dimid_ncat)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining dim ncat')

            status = nf90_def_dim(ncid,'ntilyr',ntilyr,dimid_ntilyr)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining dim ntilyr')

            status = nf90_def_dim(ncid,'ntslyr',ntslyr,dimid_ntslyr)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining dim ntslyr')

            write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         endif

         call broadcast_scalar(ncid, master_task)
      else
         filename = trim(filename_spec)
         call ice_open(nu_dump,filename,0)

         if (my_task == master_task) then
           write(nu_dump) istep1,time,time_forc
           write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
           write(nu_diag,*) 'Restart written ',istep1,time,time_forc
         endif

      endif

      diag = .true.

      if (restart_format == 'nc') then

      if (my_task == master_task) then

      allocate(dims(3))

      dims(1) = dimid_ni
      dims(2) = dimid_nj
      dims(3) = dimid_ncat

      call define_rest_field(ncid,'aicen',dims)
      call define_rest_field(ncid,'vicen',dims)
      call define_rest_field(ncid,'vsnon',dims)
      call define_rest_field(ncid,'Tsfcn',dims)

      if (tr_aero) then
         do k=1,n_aero
            write(nchar,'(i1.1)') k
            call define_rest_field(ncid,'aerosnossl'//nchar, dims)
            call define_rest_field(ncid,'aerosnoint'//nchar, dims)
            call define_rest_field(ncid,'aeroicessl'//nchar, dims)
            call define_rest_field(ncid,'aeroiceint'//nchar, dims)
         enddo
      endif

      if (tr_iage) call define_rest_field(ncid,'iage',dims)

      if (tr_FY)   call define_rest_field(ncid,'FY',dims)

      if (tr_pond) then
         call define_rest_field(ncid,'volpn',dims)
         call define_rest_field(ncid,'apondn',dims)
         call define_rest_field(ncid,'hpondn',dims)
      endif

      dims(3) = dimid_ntilyr
      call define_rest_field(ncid,'eicen',dims)
      dims(3) = dimid_ntslyr
      call define_rest_field(ncid,'esnon',dims)

      deallocate(dims)

      allocate(dims(2))
      dims(1) = dimid_ni
      dims(2) = dimid_nj

      call define_rest_field(ncid,'uvel',dims)
      call define_rest_field(ncid,'vvel',dims)

      call define_rest_field(ncid,'coszen',dims)
      call define_rest_field(ncid,'scale_factor',dims)
      call define_rest_field(ncid,'swvdr',dims)
      call define_rest_field(ncid,'swvdf',dims)
      call define_rest_field(ncid,'swidr',dims)
      call define_rest_field(ncid,'swidf',dims)

      call define_rest_field(ncid,'strocnxT',dims)
      call define_rest_field(ncid,'strocnyT',dims)

      call define_rest_field(ncid,'stressp_1',dims)
      call define_rest_field(ncid,'stressp_2',dims)
      call define_rest_field(ncid,'stressp_3',dims)
      call define_rest_field(ncid,'stressp_4',dims)

      call define_rest_field(ncid,'stressm_1',dims)
      call define_rest_field(ncid,'stressm_2',dims)
      call define_rest_field(ncid,'stressm_3',dims)
      call define_rest_field(ncid,'stressm_4',dims)

      call define_rest_field(ncid,'stress12_1',dims)
      call define_rest_field(ncid,'stress12_2',dims)
      call define_rest_field(ncid,'stress12_3',dims)
      call define_rest_field(ncid,'stress12_4',dims)

      call define_rest_field(ncid,'iceumask',dims)

      if (oceanmixed_ice) then
         call define_rest_field(ncid,'sst',dims)
         call define_rest_field(ncid,'frzmlt',dims)
      endif

      deallocate(dims)

      status = nf90_enddef(ncid)

      endif ! master_task

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------

      do n=1,ncat
         call ice_write_nc(ncid,n,'aicen',aicen(:,:,n,:),'rda8',diag)
         call ice_write_nc(ncid,n,'vicen',vicen(:,:,n,:),'rda8',diag)
         call ice_write_nc(ncid,n,'vsnon',vsnon(:,:,n,:),'rda8',diag)
         call ice_write_nc(ncid,n,'Tsfcn',trcrn(:,:,nt_Tsfc,n,:),'rda8',diag)
      enddo

      do k=1,ntilyr
         call ice_write_nc(ncid,k,'eicen',eicen(:,:,k,:),'rda8',diag)
      enddo

      do k=1,ntslyr
         call ice_write_nc(ncid,k,'esnon',esnon(:,:,k,:),'rda8',diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      call ice_write_nc(ncid,1,'uvel',uvel,'rda8',diag)
      call ice_write_nc(ncid,1,'vvel',vvel,'rda8',diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      call ice_write_nc(ncid,1,'coszen',coszen,'rda8',diag)
      call ice_write_nc(ncid,1,'scale_factor',scale_factor,'rda8',diag)
      call ice_write_nc(ncid,1,'swvdr',swvdr,'rda8',diag)
      call ice_write_nc(ncid,1,'swvdf',swvdf,'rda8',diag)
      call ice_write_nc(ncid,1,'swidr',swidr,'rda8',diag)
      call ice_write_nc(ncid,1,'swidf',swidf,'rda8',diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call ice_write_nc(ncid,1,'strocnxT',strocnxT,'rda8',diag)
      call ice_write_nc(ncid,1,'strocnyT',strocnyT,'rda8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call ice_write_nc(ncid,1,'stressp_1',stressp_1,'rda8',diag)
      call ice_write_nc(ncid,1,'stressp_2',stressp_2,'rda8',diag)
      call ice_write_nc(ncid,1,'stressp_3',stressp_3,'rda8',diag)
      call ice_write_nc(ncid,1,'stressp_4',stressp_4,'rda8',diag)

      call ice_write_nc(ncid,1,'stressm_1',stressm_1,'rda8',diag)
      call ice_write_nc(ncid,1,'stressm_2',stressm_2,'rda8',diag)
      call ice_write_nc(ncid,1,'stressm_3',stressm_3,'rda8',diag)
      call ice_write_nc(ncid,1,'stressm_4',stressm_4,'rda8',diag)

      call ice_write_nc(ncid,1,'stress12_1',stress12_1,'rda8',diag)
      call ice_write_nc(ncid,1,'stress12_2',stress12_2,'rda8',diag)
      call ice_write_nc(ncid,1,'stress12_3',stress12_3,'rda8',diag)
      call ice_write_nc(ncid,1,'stress12_4',stress12_4,'rda8',diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            if (iceumask(i,j,iblk)) work1(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_write_nc(ncid,1,'iceumask',work1,'rda8',diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
         call ice_write_nc(ncid,1,'sst',sst,'rda8',diag)
         call ice_write_nc(ncid,1,'frzmlt',frzmlt,'rda8',diag)
      endif

      if (tr_aero) then
         do k=1,n_aero
         do n=1,ncat
            write(nchar,'(i1.1)') k
            call ice_write_nc(ncid,n,'aerosnossl'//nchar, &
                trcrn(:,:,nt_aero  +(k-1)*4,n,:),'rda8',diag)
            call ice_write_nc(ncid,n,'aerosnoint'//nchar, &
                trcrn(:,:,nt_aero+1+(k-1)*4,n,:),'rda8',diag)
            call ice_write_nc(ncid,n,'aeroicessl'//nchar, &
                trcrn(:,:,nt_aero+2+(k-1)*4,n,:),'rda8',diag)
            call ice_write_nc(ncid,n,'aeroiceint'//nchar, &
                trcrn(:,:,nt_aero+3+(k-1)*4,n,:),'rda8',diag)
         enddo
         enddo
      endif

      if (tr_iage) then
         do n=1,ncat
            call ice_write_nc(ncid,n,'iage',trcrn(:,:,nt_iage,n,:),'rda8',diag)
         enddo
      endif

      if (tr_FY) then
         do n=1,ncat
            call ice_write_nc(ncid,n,'FY',trcrn(:,:,nt_FY,n,:),'rda8',diag)
         enddo
      endif

      if (tr_pond) then
         do n=1,ncat
            call ice_write_nc(ncid,n,'volpn',trcrn(:,:,nt_volpn,n,:),'rda8', &
                              diag)
            call ice_write_nc(ncid,n,'apondn',apondn(:,:,n,:),'rda8',diag)
            call ice_write_nc(ncid,n,'hpondn',hpondn(:,:,n,:),'rda8',diag)
         enddo
      endif

      if (my_task == master_task) then
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call abort_ice( &
                       'ice: Error closing netCDF restart file')
         write(nu_diag,*) 'Restart written ',istep1,time,time_forc
      endif

      else ! binary restart files

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------

      do n=1,ncat
         call ice_write(nu_dump,0,aicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vicen(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,vsnon(:,:,n,:),'ruf8',diag)
         call ice_write(nu_dump,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag)
      enddo

      do k=1,ntilyr
         call ice_write(nu_dump,0,eicen(:,:,k,:),'ruf8',diag)
      enddo

      do k=1,ntslyr
         call ice_write(nu_dump,0,esnon(:,:,k,:),'ruf8',diag)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,uvel,'ruf8',diag)
      call ice_write(nu_dump,0,vvel,'ruf8',diag)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,coszen,'ruf8',diag)
      call ice_write(nu_dump,0,scale_factor,'ruf8',diag)
      call ice_write(nu_dump,0,swvdr,'ruf8',diag)
      call ice_write(nu_dump,0,swvdf,'ruf8',diag)
      call ice_write(nu_dump,0,swidr,'ruf8',diag)
      call ice_write(nu_dump,0,swidf,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,strocnxT,'ruf8',diag)
      call ice_write(nu_dump,0,strocnyT,'ruf8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------
      call ice_write(nu_dump,0,stressp_1,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_3,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_2,'ruf8',diag)
      call ice_write(nu_dump,0,stressp_4,'ruf8',diag)

      call ice_write(nu_dump,0,stressm_1,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_3,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_2,'ruf8',diag)
      call ice_write(nu_dump,0,stressm_4,'ruf8',diag)

      call ice_write(nu_dump,0,stress12_1,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_3,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_2,'ruf8',diag)
      call ice_write(nu_dump,0,stress12_4,'ruf8',diag)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            work1(i,j,iblk) = c0
            if (iceumask(i,j,iblk)) work1(i,j,iblk) = c1
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      call ice_write(nu_dump,0,work1,'ruf8',diag)

      ! for mixed layer model
      if (oceanmixed_ice) then
         call ice_write(nu_dump,0,sst,'ruf8',diag)
         call ice_write(nu_dump,0,frzmlt,'ruf8',diag)
      endif

      if (my_task == master_task) then

         write(nu_dump) filename_volpn
         write(nu_dump) filename_aero
         write(nu_dump) filename_iage
         write(nu_dump) filename_FY

         close(nu_dump)

      endif

      endif ! netcdf or binary restart files

      end subroutine dumpfile

!=======================================================================
!BOP
!
! !IROUTINE: define_rest_field
!
! !INTERFACE:
!
      subroutine define_rest_field(ncid, vname, dims)
!
! !DESCRIPTION:
!
! Defines a restart field
!
! !REVISION HISTORY:
!
! author David A Bailey, NCAR
!
! !USES:

      character (len=*), intent(in) :: vname

      integer (kind=int_kind), intent(in) :: ncid

      integer (kind=int_kind), intent(in) :: dims(:)

      integer (kind=int_kind) :: status, varid

      status = nf90_def_var(ncid,trim(vname),nf90_double,dims,varid)
      if (status /= nf90_noerr) call abort_ice( &
                    'ice: Error defining variable '//trim(vname))

      end subroutine define_rest_field

!=======================================================================
!BOP
!
! !IROUTINE: restartfile  - restarts from a dumpfile
!
! !INTERFACE:
!
      subroutine restartfile(ice_ic)
!
! !DESCRIPTION:
!
! Restarts from a dump
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_boundary
      use ice_domain_size
      use ice_domain
      use ice_calendar, only: istep0, istep1, time, time_forc, calendar
      use ice_flux
      use ice_state
      use ice_grid, only: tmask, umask
      use ice_itd
      use ice_ocean, only: oceanmixed_ice
      use ice_work, only: work_g1, work_g2
      use ice_gather_scatter, only: scatter_global_stress
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=*), optional :: ice_ic
!EOP
!
      integer (kind=int_kind) :: &
         i, j, k, n, it, iblk ! counting indices

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag, hit_eof

      integer (kind=int_kind) :: &
         nrec

      if (present(ice_ic)) then 
         filename = ice_ic
      else
         if (my_task == master_task) then
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)
            write(nu_diag,*) 'Read ',pointer_file(1:lenstr(pointer_file))
         endif
         call broadcast_scalar(filename, master_task)
      endif

      ! read restart file

      restart_format_in = restart_format
      if (index(filename,'nc') == 0) restart_format = 'bin'

      ! Initialize all tracer fields to zero and read in from
      ! restart when available.
      if (tr_iage) trcrn(:,:,nt_iage, :,:) = c0
      if (tr_FY)   trcrn(:,:,nt_FY,   :,:) = c0
      if (tr_aero) trcrn(:,:,nt_aero:nt_aero+n_aero*4-1,:,:) = c0

      ! Need to initialize ponds in all cases.
      trcrn(:,:,nt_volpn,:,:) = c0
      apondn(:,:,:,:) = c0
      hpondn(:,:,:,:) = c0

      if (restart_format == 'nc') then
         call ice_open_nc(filename, ncid)

         if (my_task == master_task) then
            write(nu_diag,*) 'Using restart dump=', trim(filename)
            ! Need to read istep1 into istep0 here.
            status = nf90_get_att(ncid,nf90_global,'istep1',istep0)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute istep1')
            status = nf90_get_att(ncid,nf90_global,'time',time)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute time')
            status = nf90_get_att(ncid,nf90_global,'time_forc',time_forc)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error in global attribute time_forc')
            write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
         endif

         resttype = 'new'

      else
         ! determine format of restart file

         resttype = restformat(nu_restart,filename) 

         call ice_open(nu_restart,filename,0)

         if (my_task == master_task) then
            write(nu_diag,*) 'Using restart dump=', trim(filename)
            read (nu_restart) istep0,time,time_forc
            write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
         endif
      endif

      call calendar(time)

      call broadcast_scalar(istep0,master_task)

      istep1 = istep0

      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)

      diag = .true.     ! write min/max diagnostics for field
	

      if (restart_format == 'nc') then

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'
         call ice_read_nc(ncid,n,'aicen',aicen(:,:,n,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read_nc(ncid,n,'vicen',vicen(:,:,n,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read_nc(ncid,n,'vsnon',vsnon(:,:,n,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read_nc(ncid,n,'Tsfcn',trcrn(:,:,nt_Tsfc,n,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max eicen for each layer'

      do k=1,ntilyr
         call ice_read_nc(ncid,k,'eicen',eicen(:,:,k,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max esnon for each layer'

      do k=1,ntslyr
         call ice_read_nc(ncid,k,'esnon',esnon(:,:,k,:),diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read_nc(ncid,1,'uvel',uvel,diag, &
         field_type=field_type_vector,field_loc=field_loc_NEcorner)
      call ice_read_nc(ncid,1,'vvel',vvel,diag, &
         field_type=field_type_vector,field_loc=field_loc_NEcorner)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
           write(nu_diag,*) 'radiation fields'

      call ice_read_nc(ncid,1,'coszen',coszen,diag, &
                    field_loc_center, field_type_scalar)
      call ice_read_nc(ncid,1,'scale_factor',scale_factor,diag, &
                    field_loc_center, field_type_scalar)
      call ice_read_nc(ncid,1,'swvdr',swvdr,diag, &
                    field_loc_center, field_type_scalar)
      call ice_read_nc(ncid,1,'swvdf',swvdf,diag, &
                    field_loc_center, field_type_scalar)
      call ice_read_nc(ncid,1,'swidr',swidr,diag, &
                    field_loc_center, field_type_scalar)
      call ice_read_nc(ncid,1,'swidf',swidf,diag, &
                    field_loc_center, field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read_nc(ncid,1,'strocnxT',strocnxT,diag)
      call ice_read_nc(ncid,1,'strocnyT',strocnyT,diag)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
         allocate(work_g2(1,1))   ! to save memory
      endif

      call ice_read_global_nc(ncid,1,'stressp_1',work_g1,diag) ! stressp_1
      call ice_read_global_nc(ncid,1,'stressp_3',work_g2,diag) ! stressp_3
      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global_nc(ncid,1,'stressp_2',work_g1,diag) ! stressp_2
      call ice_read_global_nc(ncid,1,'stressp_4',work_g2,diag) ! stressp_4
      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global_nc(ncid,1,'stressm_1',work_g1,diag) ! stressm_1
      call ice_read_global_nc(ncid,1,'stressm_3',work_g2,diag) ! stressm_3
      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global_nc(ncid,1,'stressm_2',work_g1,diag) ! stressm_2
      call ice_read_global_nc(ncid,1,'stressm_4',work_g2,diag) ! stressm_4
      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global_nc(ncid,1,'stress12_1',work_g1,diag) ! stress12_1
      call ice_read_global_nc(ncid,1,'stress12_3',work_g2,diag) ! stress12_3
      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global_nc(ncid,1,'stress12_2',work_g1,diag) ! stress12_2
      call ice_read_global_nc(ncid,1,'stress12_4',work_g2,diag) ! stress12_4
      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read_nc(ncid,1,'iceumask',work1,diag)

      iceumask(:,:,:) = .false.
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call ice_read_nc(ncid,1,'sst',sst,diag)
         call ice_read_nc(ncid,1,'frzmlt',frzmlt,diag)
      endif

      if (tr_aero) then
         if (my_task == master_task) &
              write(nu_diag,*) 'min/max aerosols'

         do k=1,n_aero
         do n=1,ncat
            write(nchar,'(i1.1)') k
            call ice_read_nc(ncid,n,'aerosnossl'//nchar, &
               trcrn(:,:,nt_aero  +(k-1)*4,n,:),diag,    &
               field_loc_center, field_type_scalar)
            call ice_read_nc(ncid,n,'aerosnoint'//nchar, &
               trcrn(:,:,nt_aero+1+(k-1)*4,n,:),diag,    &
               field_loc_center, field_type_scalar)
            call ice_read_nc(ncid,n,'aeroicessl'//nchar, &
               trcrn(:,:,nt_aero+2+(k-1)*4,n,:),diag,    &
               field_loc_center, field_type_scalar)
            call ice_read_nc(ncid,n,'aeroiceint'//nchar, &
               trcrn(:,:,nt_aero+3+(k-1)*4,n,:),diag,    &
               field_loc_center, field_type_scalar)
         enddo
         enddo
      endif

      if (tr_iage) then
         if (my_task == master_task) &
              write(nu_diag,*) 'min/max ice age'
         do n=1,ncat
            call ice_read_nc(ncid,n,'iage',trcrn(:,:,nt_iage,n,:),diag, &
               field_loc_center, field_type_scalar)
         enddo
      endif

      if (tr_FY) then
         if (my_task == master_task) &
              write(nu_diag,*) 'min/max FY area'
         do n=1,ncat
            call ice_read_nc(ncid,n,'FY',trcrn(:,:,nt_FY,n,:),diag, &
               field_loc_center, field_type_scalar)
         enddo
      endif

      if (tr_pond) then
         if (my_task == master_task) &
              write(nu_diag,*) 'min/max melt ponds'
         do n=1,ncat
            call ice_read_nc(ncid,n,'volpn',trcrn(:,:,nt_volpn,n,:),diag, &
               field_loc_center, field_type_scalar)
            call ice_read_nc(ncid,n,'apondn',apondn(:,:,n,:),diag)
            call ice_read_nc(ncid,n,'hpondn',hpondn(:,:,n,:),diag)
         enddo
      endif

      else ! binary restart

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      do n=1,ncat
         if (my_task == master_task) &
              write(nu_diag,*) 'cat ',n, &
                               ' min/max area, vol ice, vol snow, Tsfc'

         call ice_read(nu_restart,0,aicen(:,:,n,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read(nu_restart,0,vicen(:,:,n,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read(nu_restart,0,vsnon(:,:,n,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
         call ice_read(nu_restart,0,trcrn(:,:,nt_Tsfc,n,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max eicen for each layer'
      do k=1,ntilyr
         call ice_read(nu_restart,0,eicen(:,:,k,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max esnon for each layer'
      do k=1,ntslyr
         call ice_read(nu_restart,0,esnon(:,:,k,:),'ruf8',diag, &
            field_type=field_type_scalar,field_loc=field_loc_center)
      enddo

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      call ice_read(nu_restart,0,uvel,'ruf8',diag, &
         field_type=field_type_vector,field_loc=field_loc_NEcorner)
      call ice_read(nu_restart,0,vvel,'ruf8',diag, &
         field_type=field_type_vector,field_loc=field_loc_NEcorner)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      if (trim(resttype) == 'new') then
         if (my_task == master_task) &
              write(nu_diag,*) 'radiation fields'

         call ice_read(nu_restart,0,coszen,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,scale_factor,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,swvdr,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,swvdf,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,swidr,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
         call ice_read(nu_restart,0,swidf,'ruf8',diag, &
                       field_loc_center, field_type_scalar)
      end if

      if (trim(resttype) == 'old') then
         if (my_task == master_task) &
              write(nu_diag,*) 'min/max fresh water and heat flux components'

         call ice_read(nu_restart,0,fresh,'ruf8',diag)
         call ice_read(nu_restart,0,fsalt,'ruf8',diag)
         call ice_read(nu_restart,0,fhocn,'ruf8',diag)
      endif

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      call ice_read(nu_restart,0,strocnxT,'ruf8',diag)
      call ice_read(nu_restart,0,strocnyT,'ruf8',diag)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         allocate(work_g2(nx_global,ny_global))
      else
         allocate(work_g1(1,1))
         allocate(work_g2(1,1))   ! to save memory
      endif

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_3
      call scatter_global_stress(stressp_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressp_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressp_4
      call scatter_global_stress(stressp_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressp_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_3
      call scatter_global_stress(stressm_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stressm_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stressm_4
      call scatter_global_stress(stressm_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stressm_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_1
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_3
      call scatter_global_stress(stress12_1, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_3, work_g2, work_g1, &
                                 master_task, distrb_info)

      call ice_read_global(nu_restart,0,work_g1,'ruf8',diag) ! stress12_2
      call ice_read_global(nu_restart,0,work_g2,'ruf8',diag) ! stress12_4
      call scatter_global_stress(stress12_2, work_g1, work_g2, &
                                 master_task, distrb_info)
      call scatter_global_stress(stress12_4, work_g2, work_g1, &
                                 master_task, distrb_info)

      deallocate (work_g1, work_g2)

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      call ice_read(nu_restart,0,work1,'ruf8',diag)

      iceumask(:,:,:) = .false.
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      ! for mixed layer model
      if (oceanmixed_ice) then

         if (my_task == master_task) &
              write(nu_diag,*) 'min/max sst, frzmlt'

         call ice_read(nu_restart,0,sst,'ruf8',diag)
         call ice_read(nu_restart,0,frzmlt,'ruf8',diag)
      endif

      if (my_task == master_task) then

         read(nu_restart, end=99) filename_volpn
         read(nu_restart, end=99) filename_aero
         read(nu_restart, end=99) filename_iage
         read(nu_restart, end=99) filename_FY

   99    continue
      endif

      if (my_task == master_task) then
         write(nu_diag,'(a,a)') 'filename_volpn: ',filename_volpn
         write(nu_diag,'(a,a)') 'filename_aero : ',filename_aero
         write(nu_diag,'(a,a)') 'filename_iage : ',filename_iage
         write(nu_diag,'(a,a)') 'filename_FY   : ',filename_FY
      endif
       
      if (my_task == master_task) close(nu_restart)

      call broadcast_scalar(filename_volpn, master_task)
      call broadcast_scalar(filename_aero,  master_task)
      call broadcast_scalar(filename_iage,  master_task)
      call broadcast_scalar(filename_FY,    master_task)

      restart_format = restart_format_in

      endif ! netcdf or binary restart

      !-----------------------------------------------------------------
      ! Ensure unused stress values in west and south ghost cells are 0
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, nghost
         do i = 1, nx_block
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
         do j = 1, ny_block
         do i = 1, nghost
            stressp_1 (i,j,iblk) = c0
            stressp_2 (i,j,iblk) = c0
            stressp_3 (i,j,iblk) = c0
            stressp_4 (i,j,iblk) = c0
            stressm_1 (i,j,iblk) = c0
            stressm_2 (i,j,iblk) = c0
            stressm_3 (i,j,iblk) = c0
            stressm_4 (i,j,iblk) = c0
            stress12_1(i,j,iblk) = c0
            stress12_2(i,j,iblk) = c0
            stress12_3(i,j,iblk) = c0
            stress12_4(i,j,iblk) = c0
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

      ! zero out prognostic fields at land points
      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (.not. tmask(i,j,iblk)) then
               aicen(i,j,:,iblk) = c0
               vicen(i,j,:,iblk) = c0
               vsnon(i,j,:,iblk) = c0
               trcrn(i,j,nt_Tsfc,:,iblk) = c0
               eicen(i,j,:,iblk) = c0
               esnon(i,j,:,iblk) = c0
            endif
            if (.not. umask(i,j,iblk)) then
               uvel(i,j,iblk) = c0
               vvel(i,j,iblk) = c0
            endif
         enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblk,j,i)
      do iblk = 1, nblocks

         call aggregate (nx_block, ny_block, &
                         aicen(:,:,:,iblk),  &
                         trcrn(:,:,:,:,iblk),&
                         vicen(:,:,:,iblk),  &
                         vsnon(:,:,:,iblk),  &
                         eicen(:,:,:,iblk),  &
                         esnon(:,:,:,iblk),  &
                         aice (:,:,  iblk),  &
                         trcr (:,:,:,iblk),  &
                         vice (:,:,  iblk),  &
                         vsno (:,:,  iblk),  &
                         eice (:,:,  iblk),  &
                         esno (:,:,  iblk),  &
                         aice0(:,:,  iblk),  &
                         tmask(:,:,  iblk),  &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo
      !$OMP END PARALLEL DO

      end subroutine restartfile

!=======================================================================
!BOP
!
! !IROUTINE: integer function lenstr(label) - compute length string
!
! !INTERFACE:
!
      integer function lenstr(label)
!
! !DESCRIPTION:
!
! Compute length of string by finding first non-blank
! character from the right.
!
! !REVISION HISTORY:
!
! author:   ?
!
! !INPUT/OUTPUT PARAMETERS:
!
      character*(*) label
!
!EOP
!
      integer (kind=int_kind) :: &
         length, & ! length of character string
         n         ! loop index

      length = len(label)
      do n=length,1,-1
        if( label(n:n) /= ' ' ) exit
      enddo
      lenstr = n

      end function lenstr

!=======================================================================
!BOP
!
! !IROUTINE: character function restformat - determine format of restart file
!
! !INTERFACE:
!
      character(len=char_len) function restformat(nu, filename)
!
! !DESCRIPTION:
!
! Determine number of records in restart file
!
! !REVISION HISTORY:
!
! author: Mariana Vertenstein 12/2008
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent(in) :: nu
      character(len=char_len), intent(in) :: filename
!
!EOP
!
      logical :: readdata
      integer :: nrec
      integer :: idummy, ios

      integer, parameter :: nrecold = ncat*4+ntilyr+ntslyr+21
      integer, parameter :: nrecnew = ncat*4+ntilyr+ntslyr+27

      if (my_task == master_task) then
         call ice_open(nu_restart,filename,0)

         readdata = .true.
         nrec = 0
         do while (readdata)
            read(nu, iostat=ios) idummy 
            if (ios < 0) then
               readdata = .false.
            else
               nrec = nrec + 1
            end if
         end do

         if (nrec == nrecold) then
            restformat = 'old'
         else if (nrec >= nrecnew) then
            restformat = 'new'
	 else
            call abort_ice ( & 
                 'restformat: number of records on restart file not supported')
         end if
         write(nu_diag,*)'Number of records restart file = ',nrec,' restart format is = ',trim(restformat)
         close(nu)
      end if
      call broadcast_scalar(restformat,master_task)

    end function restformat

!=======================================================================

      end module ice_restart

!=======================================================================
