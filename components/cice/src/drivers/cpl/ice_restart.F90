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
      use ice_boundary
      use ice_calendar, only: sec, month, mday, nyr, istep0, istep1, &
                              time, time_forc, idate, year_init, calendar
      use ice_read_write
      use ice_fileunits
      use ice_timers
      use ice_exit, only: abort_ice	
!
!EOP
!
      implicit none
      save

      interface dumpfile
         module procedure dumpfile_bin
         module procedure dumpfile_pio
      end interface

      interface restartfile
         module procedure restartfile_bin
         module procedure restartfile_pio
      end interface

      character (len=char_len) :: &       	 
         resttype           ! type of restart format ('new' or 'old')

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
         restart_format    ! format of restart files 'nc' or 'bin'

      character (len=char_len_long) :: &
         pointer_file      ! input pointer file for restarts

      real (kind=dbl_kind), private, &
         dimension(nx_block,ny_block,max_blocks) :: &
         work1

      logical (kind=log_kind) :: lcdf64

      integer (kind=int_kind) :: ncid ! netcdf restart file id

      integer (kind=int_kind) :: &
        status        ! status variable from netCDF routine

      character (len=1) :: nchar

!=======================================================================

      contains

!=======================================================================

!=======================================================================
!---subroutines write/read Fortran unformatted/netcdf data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile - dumps all fields required for restart
!
! !INTERFACE:
!
    subroutine dumpfile_bin()
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
      use ice_state
      use ice_dyn_evp
      use ice_blocks, only : block, get_block, nx_block, ny_block
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      logical (kind=log_kind) :: diag

      character(len=char_len_long) :: filename

      iyear = nyr + year_init - 1
      imonth = month
      iday = mday
      
      write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
            restart_dir(1:lenstr(restart_dir)), &
            restart_file(1:lenstr(restart_file)),'.', &
            iyear,'-',month,'-',mday,'-',sec

      call ice_open(nu_dump,filename,0)

      if (my_task == master_task) then
         write(nu_dump) istep1,time,time_forc
         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
         write(nu_diag,*) 'Restart written ',istep1,time,time_forc
      endif

      diag = .true.

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

      if (my_task == master_task) then
         write(nu_dump) filename_volpn
         write(nu_dump) filename_aero
         write(nu_dump) filename_iage
         write(nu_dump) filename_FY
         write(nu_dump) filename_lvl

         close(nu_dump)
      endif

    end subroutine dumpfile_bin

!=======================================================================
!
!BOP
!
! !IROUTINE: dumpfile_pio - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine dumpfile_pio(filename_spec)
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
      use ice_state
      use ice_dyn_evp
      use ice_blocks, only : block, get_block, nx_block, ny_block
      use ice_pio	
      use pio
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in) :: filename_spec
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          ilo, ihi, jlo, jhi,   & ! counting indices
          lon, lat                ! global indices

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      integer (kind=int_kind) :: dimid_ni, dimid_nj, dimid_ncat, &
                                 dimid_ntilyr, dimid_ntslyr

      integer (kind=int_kind), allocatable :: dims(:)

      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d
      type(io_desc_t)       :: iodesc3d_ncat
      type(io_desc_t)       :: iodesc3d_ntilyr
      type(io_desc_t)       :: iodesc3d_ntslyr
      type(var_desc_t)      :: varid

      type(block) :: this_block 

      ! construct path/file
      filename = trim(filename_spec) // ".nc"
         
      ! write pointer (path/file)
      if (my_task == master_task) then
        open(nu_rst_pointer,file=pointer_file)
        write(nu_rst_pointer,'(a)') filename
        close(nu_rst_pointer)
      endif

      ! begin writing restart data

      File%fh=-1
      call ice_pio_init(mode='write',filename=trim(filename), File=File, &
           clobber=.true., cdf64=lcdf64 )


      status = pio_put_att(File,pio_global,'istep1',istep1)
      status = pio_put_att(File,pio_global,'time',time)
      status = pio_put_att(File,pio_global,'time_forc',time_forc)
      status = pio_put_att(File,pio_global,'nyr',nyr)
      status = pio_put_att(File,pio_global,'month',month)
      status = pio_put_att(File,pio_global,'mday',mday)
      status = pio_put_att(File,pio_global,'sec',sec)

      status = pio_def_dim(File,'ni',nx_global,dimid_ni)
      status = pio_def_dim(File,'nj',ny_global,dimid_nj)
      status = pio_def_dim(File,'ncat',ncat,dimid_ncat)
      status = pio_def_dim(File,'ntilyr',ntilyr,dimid_ntilyr)
      status = pio_def_dim(File,'ntslyr',ntslyr,dimid_ntslyr)
      if (my_task == master_task) then
         write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif
      diag = .true.

      allocate(dims(3))

      dims(1) = dimid_ni
      dims(2) = dimid_nj
      dims(3) = dimid_ncat

      call define_rest_field(File,'aicen',dims)
      call define_rest_field(File,'vicen',dims)
      call define_rest_field(File,'vsnon',dims)
      call define_rest_field(File,'Tsfcn',dims)

      if (tr_aero) then
         do k=1,n_aero
            write(nchar,'(i1.1)') k
            call define_rest_field(File,'aerosnossl'//nchar, dims)
            call define_rest_field(File,'aerosnoint'//nchar, dims)
            call define_rest_field(File,'aeroicessl'//nchar, dims)
            call define_rest_field(File,'aeroiceint'//nchar, dims)
         enddo
      endif

      if (tr_iage) then
         call define_rest_field(File,'iage',dims)
      end if

      if (tr_FY) then
         call define_rest_field(File,'FY',dims)
      end if

      if (tr_lvl) then
         call define_rest_field(File,'alvl',dims)
         call define_rest_field(File,'vlvl',dims)
      end if

      if (tr_pond) then
         call define_rest_field(File,'volpn' ,dims)
         call define_rest_field(File,'apondn',dims)
         call define_rest_field(File,'hpondn',dims)
      end if

      dims(3) = dimid_ntilyr
      call define_rest_field(File,'eicen',dims)
      dims(3) = dimid_ntslyr
      call define_rest_field(File,'esnon',dims)

      deallocate(dims)

      allocate(dims(2))
      dims(1) = dimid_ni
      dims(2) = dimid_nj

      call define_rest_field(File,'uvel',dims)
      call define_rest_field(File,'vvel',dims)

      call define_rest_field(File,'coszen',dims)
      call define_rest_field(File,'scale_factor',dims)
      call define_rest_field(File,'swvdr',dims)
      call define_rest_field(File,'swvdf',dims)
      call define_rest_field(File,'swidr',dims)
      call define_rest_field(File,'swidf',dims)

      call define_rest_field(File,'strocnxT',dims)
      call define_rest_field(File,'strocnyT',dims)

      call define_rest_field(File,'stressp_1',dims)
      call define_rest_field(File,'stressp_2',dims)
      call define_rest_field(File,'stressp_3',dims)
      call define_rest_field(File,'stressp_4',dims)

      call define_rest_field(File,'stressm_1',dims)
      call define_rest_field(File,'stressm_2',dims)
      call define_rest_field(File,'stressm_3',dims)
      call define_rest_field(File,'stressm_4',dims)

      call define_rest_field(File,'stress12_1',dims)
      call define_rest_field(File,'stress12_2',dims)
      call define_rest_field(File,'stress12_3',dims)
      call define_rest_field(File,'stress12_4',dims)

      call define_rest_field(File,'iceumask',dims)

      status = pio_enddef(File)

      deallocate(dims)

      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      call ice_pio_initdecomp(ndim3=ncat  , iodesc=iodesc3d_ncat)

      status = pio_inq_varid(File,'aicen',varid)
      call pio_write_darray(File, varid, iodesc3d_ncat, aicen(:,:,:,:), status, fillval=c0)
      
      status = pio_inq_varid(File,'vicen',varid)
      call pio_write_darray(File, varid, iodesc3d_ncat, vicen, status, fillval=c0)
      
      status = pio_inq_varid(File,'vsnon',varid)
      call pio_write_darray(File, varid, iodesc3d_ncat, vsnon, status, fillval=c0)

      status = pio_inq_varid(File,'Tsfcn',varid)
      call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_Tsfc,:,:), status, fillval=c0)

      call ice_pio_initdecomp(ndim3=ntilyr, iodesc=iodesc3d_ntilyr)
      status = pio_inq_varid(File,'eicen',varid)
      call pio_write_darray(File, varid, iodesc3d_ntilyr, eicen, status, fillval=c0)
      call PIO_freeDecomp(File,iodesc3d_ntilyr)

      call ice_pio_initdecomp(ndim3=ntslyr, iodesc=iodesc3d_ntslyr)
      status = pio_inq_varid(File,'esnon',varid)
      call pio_write_darray(File, varid, iodesc3d_ntslyr, esnon, status, fillval=c0)

      call PIO_freeDecomp(File,iodesc3d_ntslyr)

      call ice_pio_initdecomp(iodesc=iodesc2d)
      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      status = pio_inq_varid(File,'uvel',varid)
      call pio_write_darray(File, varid, iodesc2d, uvel, status, fillval=c0)

      status = pio_inq_varid(File,'vvel',varid)
      call pio_write_darray(File, varid, iodesc2d, vvel, status, fillval=c0)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------
      status = pio_inq_varid(File,'coszen',varid)
      call pio_write_darray(File, varid, iodesc2d, coszen, status, fillval=c0)

      status = pio_inq_varid(File,'scale_factor',varid)
      call pio_write_darray(File, varid, iodesc2d, scale_factor, status, fillval=c0)

      status = pio_inq_varid(File,'swvdr',varid)
      call pio_write_darray(File, varid, iodesc2d, swvdr, status, fillval=c0)

      status = pio_inq_varid(File,'swvdf',varid)
      call pio_write_darray(File, varid, iodesc2d, swvdf, status, fillval=c0)

      status = pio_inq_varid(File,'swidr',varid)
      call pio_write_darray(File, varid, iodesc2d, swidr, status, fillval=c0)

      status = pio_inq_varid(File,'swidf',varid)
      call pio_write_darray(File, varid, iodesc2d, swidf, status, fillval=c0)

      !-----------------------------------------------------------------
      ! ocean stress (for bottom heat flux in thermo)
      !-----------------------------------------------------------------
      status = pio_inq_varid(File,'strocnxT',varid)
      call pio_write_darray(File, varid, iodesc2d, strocnxT, status, fillval=c0)

      status = pio_inq_varid(File,'strocnyT',varid)
      call pio_write_darray(File, varid, iodesc2d, strocnyT, status, fillval=c0)

      !-----------------------------------------------------------------
      ! internal stress
      !-----------------------------------------------------------------

      status = pio_inq_varid(File,'stressp_1',varid)
      call pio_write_darray(File, varid, iodesc2d, stressp_1, status, fillval=c0)

      status = pio_inq_varid(File,'stressp_2',varid)
      call pio_write_darray(File, varid, iodesc2d, stressp_2, status, fillval=c0)

      status = pio_inq_varid(File,'stressp_3',varid)
      call pio_write_darray(File, varid, iodesc2d, stressp_3, status, fillval=c0)

      status = pio_inq_varid(File,'stressp_4',varid)
      call pio_write_darray(File, varid, iodesc2d, stressp_4, status, fillval=c0)

      status = pio_inq_varid(File,'stressm_1',varid)
      call pio_write_darray(File, varid, iodesc2d, stressm_1, status, fillval=c0)

      status = pio_inq_varid(File,'stressm_2',varid)
      call pio_write_darray(File, varid, iodesc2d, stressm_2, status, fillval=c0)
      
      status = pio_inq_varid(File,'stressm_3',varid)
      call pio_write_darray(File, varid, iodesc2d, stressm_3, status, fillval=c0)

      status = pio_inq_varid(File,'stressm_4',varid)
      call pio_write_darray(File, varid, iodesc2d, stressm_4, status, fillval=c0)

      status = pio_inq_varid(File,'stress12_1',varid)
      call pio_write_darray(File, varid, iodesc2d, stress12_1, status, fillval=c0)

      status = pio_inq_varid(File,'stress12_2',varid)
      call pio_write_darray(File, varid, iodesc2d, stress12_2, status, fillval=c0)

      status = pio_inq_varid(File,'stress12_3',varid)
      call pio_write_darray(File, varid, iodesc2d, stress12_3, status, fillval=c0)

      status = pio_inq_varid(File,'stress12_4',varid)
      call pio_write_darray(File, varid, iodesc2d, stress12_4, status, fillval=c0)

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
      status = pio_inq_varid(File,'iceumask',varid)
      call pio_write_darray(File, varid, iodesc2d, work1, status, fillval=c0)
      call PIO_freeDecomp(File,iodesc2d)

      if (tr_aero) then
         do k=1,n_aero
            write(nchar,'(i1.1)') k

            status = pio_inq_varid(File,'aerosnossl'//nchar,varid)
            call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_aero+  (k-1)*4,:,:), status, fillval=c0)

            status = pio_inq_varid(File,'aerosnoint'//nchar,varid)
            call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_aero+1+(k-1)*4,:,:), status, fillval=c0)

            status = pio_inq_varid(File,'aeroicessl'//nchar,varid)
            call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_aero+2+(k-1)*4,:,:), status, fillval=c0)

            status = pio_inq_varid(File,'aeroiceint'//nchar,varid)
            call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_aero+3+(k-1)*4,:,:), status, fillval=c0)
         enddo
      endif

      if (tr_iage) then
         status = pio_inq_varid(File,'iage',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_iage,:,:), status, fillval=c0)
      endif

      if (tr_FY) then
         status = pio_inq_varid(File,'FY',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_FY,:,:), status, fillval=c0)
      endif

      if (tr_lvl) then
         status = pio_inq_varid(File,'alvl',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_alvl,:,:), status, fillval=c0)
         status = pio_inq_varid(File,'vlvl',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_vlvl,:,:), status, fillval=c0)
      endif

      if (tr_pond) then
         status = pio_inq_varid(File,'volpn',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_volpn,:,:), status, fillval=c0)

         status = pio_inq_varid(File,'apondn',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, apondn, status, fillval=c0)

         status = pio_inq_varid(File,'hpondn',varid)
         call pio_write_darray(File, varid, iodesc3d_ncat, hpondn, status, fillval=c0)
      endif

      call PIO_freeDecomp(File,iodesc3d_ncat)

      call pio_closefile(File)

      if (my_task == master_task) then
         write(nu_diag,*) 'Restart written ',istep1,time,time_forc
      endif

    end subroutine dumpfile_pio

!=======================================================================
!BOP
!
! !IROUTINE: define_rest_field
!
! !INTERFACE:
!
      subroutine define_rest_field(File, vname, dims)
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

        use pio
!
! !INPUT/OUTPUT PARAMETERS:
!
        type(file_desc_t)      , intent(in)  :: File
        character (len=*)      , intent(in)  :: vname
        integer (kind=int_kind), intent(in)  :: dims(:)
        !
!EOP
!
        type(var_desc_t) :: varid

        status = pio_def_var(File,trim(vname),pio_double,dims,varid)
        
      end subroutine define_rest_field

!=======================================================================
!BOP
!
! !IROUTINE: restartfile  - restarts from a dumpfile
!
! !INTERFACE:
!
      subroutine restartfile_bin(ice_ic)
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
      use ice_flux
      use ice_state
      use ice_grid, only: tmask, umask, grid_type
      use ice_itd
      use ice_work, only: work_g1, work_g2
      use ice_gather_scatter, only: scatter_global_stress, gather_global
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=*), optional :: ice_ic
!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          ilo, ihi, jlo, jhi,   & ! counting indices
          lon, lat,             & ! global indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

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

      ! Initialize all tracer fields to zero and read in from
      ! restart when available.
      if (tr_iage) trcrn(:,:,nt_iage, :,:) = c0
      if (tr_FY)   trcrn(:,:,nt_FY,   :,:) = c0
      if (tr_lvl)  trcrn(:,:,nt_alvl, :,:) = c1
      if (tr_lvl)  trcrn(:,:,nt_vlvl, :,:) = c1
      if (tr_aero) trcrn(:,:,nt_aero:nt_aero+n_aero*4-1,:,:) = c0

      ! Need to initialize ponds in all cases.
      trcrn(:,:,nt_volpn,:,:) = c0
      apondn(:,:,:,:) = c0
      hpondn(:,:,:,:) = c0

      ! determine format of binary restart file
      resttype = restformat(nu_restart,filename) 
      
      call ice_open(nu_restart,filename,0)

      if (my_task == master_task) then
         write(nu_diag,*) 'Using restart dump=', trim(filename)
         read (nu_restart) istep0,time,time_forc
         write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
      endif
      
      call broadcast_scalar(istep0,master_task)
      
      istep1 = istep0
      
      call broadcast_scalar(time,master_task)
      call broadcast_scalar(time_forc,master_task)
      
      diag = .true.     ! write min/max diagnostics for field
	
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

      if (my_task == master_task) then

         read(nu_restart, end=99) filename_volpn
         read(nu_restart, end=99) filename_aero
         read(nu_restart, end=99) filename_iage
         read(nu_restart, end=99) filename_FY
         read(nu_restart, end=99) filename_lvl

   99    continue
      endif

      if (my_task == master_task) then
         write(nu_diag,'(a,a)') 'filename_volpn: ',filename_volpn
         write(nu_diag,'(a,a)') 'filename_aero : ',filename_aero
         write(nu_diag,'(a,a)') 'filename_iage : ',filename_iage
         write(nu_diag,'(a,a)') 'filename_FY   : ',filename_FY
         write(nu_diag,'(a,a)') 'filename_lvl  : ',filename_lvl
      endif
       
      if (my_task == master_task) close(nu_restart)

      call broadcast_scalar(filename_volpn, master_task)
      call broadcast_scalar(filename_aero,  master_task)
      call broadcast_scalar(filename_iage,  master_task)
      call broadcast_scalar(filename_FY,    master_task)
      call broadcast_scalar(filename_lvl,   master_task)

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
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

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
                         ntrcr, trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo
      !$OMP END PARALLEL DO

    end subroutine restartfile_bin

!=======================================================================
!BOP
!
! !IROUTINE: restartfile  - restarts from a dumpfile
!
! !INTERFACE:
!
      subroutine restartfile_pio(usepio, ice_ic)
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
      use ice_flux
      use ice_state
      use ice_grid, only: tmask, umask, grid_type
      use ice_itd
      use ice_pio	
      use pio
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical, intent(in) :: usepio 
      character(len=*), intent(in), optional :: ice_ic
!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          ilo, ihi, jlo, jhi,   & ! counting indices
          lon, lat                ! global indices

      character(len=char_len_long) :: &
         filename, filename0

      logical (kind=log_kind) :: &
         diag

      type(block) :: this_block 

      type(iosystem_desc_t) :: pio_subsystem
      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d
      type(io_desc_t)       :: iodesc3d_ncat
      type(io_desc_t)       :: iodesc3d_ntilyr
      type(io_desc_t)       :: iodesc3d_ntslyr
      type(var_desc_t)      :: varid

      resttype = 'new'

      ! Initialize all tracer fields to zero and read in from
      ! restart when available.
      if (tr_iage) trcrn(:,:,nt_iage, :,:) = c0
      if (tr_FY)   trcrn(:,:,nt_FY,   :,:) = c0
      if (tr_lvl)  trcrn(:,:,nt_alvl, :,:) = c1
      if (tr_lvl)  trcrn(:,:,nt_vlvl, :,:) = c1
      if (tr_aero) trcrn(:,:,nt_aero:nt_aero+n_aero*4-1,:,:) = c0

      ! Need to initialize ponds in all cases.
      trcrn(:,:,nt_volpn,:,:) = c0
      apondn(:,:,:,:) = c0
      hpondn(:,:,:,:) = c0

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

      ! Read restart file

      File%fh=-1
      call ice_pio_init(mode='read', filename=trim(filename), File=File)
      
      call ice_pio_initdecomp(iodesc=iodesc2d)
      call ice_pio_initdecomp(ndim3=ncat  , iodesc=iodesc3d_ncat)
      call ice_pio_initdecomp(ndim3=ntilyr, iodesc=iodesc3d_ntilyr)
      call ice_pio_initdecomp(ndim3=ntslyr, iodesc=iodesc3d_ntslyr)

      if (my_task == master_task) then
         write(nu_diag,*) 'Using restart dump=', trim(filename)
      end if
      
      status = pio_get_att(File, pio_global, 'istep1', istep1)
      status = pio_get_att(File, pio_global, 'time', time)
      status = pio_get_att(File, pio_global, 'time_forc', time_forc)
      call pio_seterrorhandling(File, PIO_BCAST_ERROR)
      status = pio_get_att(File, pio_global, 'nyr', nyr)
      if (status == PIO_noerr) then
         status = pio_get_att(File, pio_global, 'month', month)
         status = pio_get_att(File, pio_global, 'mday', mday)
         status = pio_get_att(File, pio_global, 'sec', sec)
      endif
      call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      
      if (my_task == master_task) then
         write(nu_diag,*) 'Restart read at istep=',istep0,time,time_forc
      endif
      
      istep1 = istep0

      diag = .true.     ! write min/max diagnostics for field
	
      !-----------------------------------------------------------------
      ! state variables
      !-----------------------------------------------------------------
      if (my_task == master_task) &
         write(nu_diag,*) ' min/max area, vol ice, vol snow, Tsfc'

      status = pio_inq_varid(File,'aicen',varid)
      call pio_read_darray(File, varid, iodesc3d_ncat, aicen, status)

      status = pio_inq_varid(File,'vicen',varid)
      call pio_read_darray(File, varid, iodesc3d_ncat, vicen, status)

      status = pio_inq_varid(File,'vsnon',varid)
      call pio_read_darray(File, varid, iodesc3d_ncat, vsnon, status)

      status = pio_inq_varid(File,'Tsfcn',varid)
      call pio_read_darray(File, varid, iodesc3d_ncat, trcrn(:,:,nt_Tsfc,:,:), status)

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max eicen for each layer'

      status = pio_inq_varid(File,'eicen',varid)
      call pio_read_darray(File, varid, iodesc3d_ntilyr, eicen, status)

      if (my_task == master_task) &
           write(nu_diag,*) 'min/max esnon for each layer'

      status = pio_inq_varid(File,'esnon',varid)
      call pio_read_darray(File, varid, iodesc3d_ntslyr, esnon, status)

      !-----------------------------------------------------------------
      ! velocity
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max velocity components'

      status = pio_inq_varid(File,'uvel',varid)
      call pio_read_darray(File, varid, iodesc2d, uvel, status)
      call ice_HaloUpdate (uvel,               halo_info,     &
                           field_loc_NEcorner, field_type_vector)

      status = pio_inq_varid(File,'vvel',varid)
      call pio_read_darray(File, varid, iodesc2d, vvel, status)
      call ice_HaloUpdate (vvel,               halo_info,     &
                           field_loc_NEcorner, field_type_vector)

      !-----------------------------------------------------------------
      ! radiation fields
      !-----------------------------------------------------------------

      if (my_task == master_task) &
           write(nu_diag,*) 'radiation fields'

      status = pio_inq_varid(File,'coszen',varid)
      call pio_read_darray(File, varid, iodesc2d, coszen, status)
      call ice_HaloUpdate (coszen,               halo_info,     &
                           field_loc_center,    field_type_scalar)

      status = pio_inq_varid(File,'scale_factor',varid)
      call pio_read_darray(File, varid, iodesc2d, scale_factor, status)
      call ice_HaloUpdate (scale_factor,        halo_info,     &
                           field_loc_center,    field_type_scalar)

      status = pio_inq_varid(File,'swvdr',varid)
      call pio_read_darray(File, varid, iodesc2d, swvdr, status)
      call ice_HaloUpdate (swvdr,               halo_info,     &
                           field_loc_center,    field_type_scalar)

      status = pio_inq_varid(File,'swvdf',varid)
      call pio_read_darray(File, varid, iodesc2d, swvdf, status)
      call ice_HaloUpdate (swvdf,               halo_info,     &
                           field_loc_center,    field_type_scalar)

      status = pio_inq_varid(File,'swidr',varid)
      call pio_read_darray(File, varid, iodesc2d, swidr, status)
      call ice_HaloUpdate (swidr,               halo_info,     &
                           field_loc_center,    field_type_scalar)

      status = pio_inq_varid(File,'swidf',varid)
      call pio_read_darray(File, varid, iodesc2d, swidf, status)
      call ice_HaloUpdate (swidf,               halo_info,     &
                           field_loc_center,    field_type_scalar)

      !-----------------------------------------------------------------
      ! ocean stress
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'min/max ocean stress components'

      status = pio_inq_varid(File,'strocnxT',varid)
      call pio_read_darray(File, varid, iodesc2d, strocnxT, status)

      status = pio_inq_varid(File,'strocnyT',varid)
      call pio_read_darray(File, varid, iodesc2d, strocnyT, status)

      !-----------------------------------------------------------------
      ! internal stress
      ! The stress tensor must be read and scattered in pairs in order
      ! to properly match corner values across a tripole grid cut.
      !-----------------------------------------------------------------
      if (my_task == master_task) write(nu_diag,*) &
           'internal stress components'
      
      status = pio_inq_varid(File,'stressp_1',varid)
      call pio_read_darray(File, varid, iodesc2d, stressp_1, status)

      status = pio_inq_varid(File,'stressp_3',varid)
      call pio_read_darray(File, varid, iodesc2d, stressp_3, status)

      status = pio_inq_varid(File,'stressp_2',varid)
      call pio_read_darray(File, varid, iodesc2d, stressp_2, status)

      status = pio_inq_varid(File,'stressp_4',varid)
      call pio_read_darray(File, varid, iodesc2d, stressp_4, status)

      status = pio_inq_varid(File,'stressm_1',varid)
      call pio_read_darray(File, varid, iodesc2d, stressm_1, status)

      status = pio_inq_varid(File,'stressm_3',varid)
      call pio_read_darray(File, varid, iodesc2d, stressm_3, status)

      status = pio_inq_varid(File,'stressm_2',varid)
      call pio_read_darray(File, varid, iodesc2d, stressm_2, status)

      status = pio_inq_varid(File,'stressm_4',varid)
      call pio_read_darray(File, varid, iodesc2d, stressm_4, status)

      status = pio_inq_varid(File,'stress12_1',varid)
      call pio_read_darray(File, varid, iodesc2d, stress12_1, status)

      status = pio_inq_varid(File,'stress12_3',varid)
      call pio_read_darray(File, varid, iodesc2d, stress12_3, status)

      status = pio_inq_varid(File,'stress12_2',varid)
      call pio_read_darray(File, varid, iodesc2d, stress12_2, status)

      status = pio_inq_varid(File,'stress12_4',varid)
      call pio_read_darray(File, varid, iodesc2d, stress12_4, status)

      call ice_HaloUpdate(stressp_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressp_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressp_2, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressp_4, halo_info, &
                           field_loc_center,  field_type_scalar)

      call ice_HaloUpdate(stressm_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressm_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressm_2, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stressm_4, halo_info, &
                           field_loc_center,  field_type_scalar)

      call ice_HaloUpdate(stress12_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stress12_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stress12_2, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate(stress12_4, halo_info, &
                           field_loc_center,  field_type_scalar)

      ! Special halo updates for tripole grid

      if (trim(grid_type) == 'tripole') then

      call ice_HaloUpdate_stress(stressp_1, stressp_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressp_3, stressp_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressp_2, stressp_4, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressp_4, stressp_2, halo_info, &
                           field_loc_center,  field_type_scalar)

      call ice_HaloUpdate_stress(stressm_1, stressm_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressm_3, stressm_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressm_2, stressm_4, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stressm_4, stressm_2, halo_info, &
                           field_loc_center,  field_type_scalar)

      call ice_HaloUpdate_stress(stress12_1, stress12_3, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stress12_3, stress12_1, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stress12_2, stress12_4, halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate_stress(stress12_4, stress12_2, halo_info, &
                           field_loc_center,  field_type_scalar)

      endif

      !-----------------------------------------------------------------
      ! ice mask for dynamics
      !-----------------------------------------------------------------
      if (my_task == master_task) &
           write(nu_diag,*) 'ice mask for dynamics'

      status = pio_inq_varid(File,'iceumask',varid)
      call pio_read_darray(File, varid, iodesc2d, work1, status)
      call ice_HaloUpdate (work1,             halo_info, &
                           field_loc_center,  field_type_scalar)

      iceumask(:,:,:) = .false.
      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            if (work1(i,j,iblk) > p5) iceumask(i,j,iblk) = .true.
         enddo
      enddo
      enddo

      if (tr_aero) then
         call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         status = pio_inq_varid(File,'aerosnossl1',varid)
         if (status == PIO_noerr) then
            do k=1,n_aero
               write(nchar,'(i1.1)') k

               status = pio_inq_varid(File,'aerosnossl'//nchar,varid)
               call pio_read_darray(File, varid, iodesc3d_ncat, &
                                    trcrn(:,:,nt_aero+(k-1)*4,:,:), status)

               status = pio_inq_varid(File,'aerosnoint'//nchar,varid)
               call pio_read_darray(File, varid, iodesc3d_ncat, &
                                    trcrn(:,:,nt_aero+1+(k-1)*4,:,:), status)
   
               status = pio_inq_varid(File,'aeroicessl'//nchar,varid)
               call pio_read_darray(File, varid, iodesc3d_ncat, &
                                    trcrn(:,:,nt_aero+2+(k-1)*4,:,:), status)

               status = pio_inq_varid(File,'aeroiceint'//nchar,varid)
               call pio_read_darray(File, varid, iodesc3d_ncat, &
                                    trcrn(:,:,nt_aero+3+(k-1)*4,:,:), status)
            enddo
         endif
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      endif

      if (tr_iage) then
         call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         status = pio_inq_varid(File,'iage',varid)
         if (status == PIO_noerr) then
            call pio_read_darray(File, varid, iodesc3d_ncat, &
	                         trcrn(:,:,nt_iage,:,:), status)
         endif
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      endif

      if (tr_FY) then
         call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         status = pio_inq_varid(File,'FY',varid)
         if (status == PIO_noerr) then
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 trcrn(:,:,nt_FY,:,:), status)
         endif
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      endif

      if (tr_lvl) then
         call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         status = pio_inq_varid(File,'alvl',varid)
         if (status == PIO_noerr) then
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 trcrn(:,:,nt_alvl,:,:), status)
         endif
         status = pio_inq_varid(File,'vlvl',varid)
         if (status == PIO_noerr) then
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 trcrn(:,:,nt_vlvl,:,:), status)
         endif
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      endif

      if (tr_pond) then
         call pio_seterrorhandling(File, PIO_BCAST_ERROR)
         status = pio_inq_varid(File,'volpn',varid)
         if (status == PIO_noerr) then
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 trcrn(:,:,nt_volpn,:,:), status)
         
            status = pio_inq_varid(File,'apondn',varid)
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 apondn(:,:,:,:), status)
         
            status = pio_inq_varid(File,'hpondn',varid)
            call pio_read_darray(File, varid, iodesc3d_ncat, &
                                 hpondn(:,:,:,:), status)
         endif
         call pio_seterrorhandling(File, PIO_INTERNAL_ERROR)
      endif

      call pio_closefile(File)

      call PIO_freeDecomp(File,iodesc2d)
      call PIO_freeDecomp(File,iodesc3d_ncat)
      call PIO_freeDecomp(File,iodesc3d_ntilyr)
      call PIO_freeDecomp(File,iodesc3d_ntslyr)

      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)

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
      ! Ensure ice is binned in correct categories
      ! (should not be necessary unless restarting from a run with
      !  different category boundaries).
      !
      ! If called, this subroutine does not give exact restart.
      !-----------------------------------------------------------------
!!!      call cleanup_itd

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
                         ntrcr, trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo
      !$OMP END PARALLEL DO

    end subroutine restartfile_pio

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
      character(len=char_len_long), intent(in) :: filename
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
