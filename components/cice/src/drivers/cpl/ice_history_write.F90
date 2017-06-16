!=======================================================================
!
!BOP
!
! !MODULE: ice_history_write - ice model history files writes using pio
!
! Output files: netCDF 
!
! !INTERFACE:
!
      module ice_history_write
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
      use ice_blocks
      use ice_grid
      use ice_fileunits
      use ice_history_fields
!
!EOP
!
      implicit none
      save

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: icecdf - write netCDF history file
!
! !INTERFACE:
!
      subroutine icecdf(ns)
!
! !DESCRIPTION:
!
! write netCDF history file
!
! !REVISION HISTORY:
!
! authors:   Mariana Vertenstein, NCAR
!
! !USES:
!
      use shr_mpi_mod
!     use shr_mem_mod, only : shr_get_memusage
      use ice_gather_scatter
      use ice_domain_size
      use ice_constants
      use ice_calendar, only: time, sec, idate, idate0, nyr, month, &
                              mday, write_ic, histfreq, histfreq_n, &
                              year_init, new_year, new_month, new_day, &
                              dayyr, daymo, days_per_year
      use ice_work, only: work_g1, work_gr, work_gr3
      use ice_restart, only: lenstr, runid, lcdf64
      use ice_domain, only: distrb_info
      use ice_itd, only: c_hi_range
      use ice_exit
      use ice_pio	
      use pio
      use shr_sys_mod, only : shr_sys_flush
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind), intent(in) :: ns

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat
      integer (kind=int_kind) :: i,j,n,k, &
         status,imtid,jmtid,timid, &
         length,nvertexid,ivertex
      integer (kind=int_kind), dimension(2) :: dimid2
      integer (kind=int_kind), dimension(3) :: dimid3
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      real (kind=real_kind) :: ltime
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm)

      integer (kind=int_kind) :: iyear, imonth, iday
      integer (kind=int_kind) :: icategory,ind,i_aice,boundid

      character (char_len) :: start_time,current_date,current_time
      character (len=16) :: c_aice
      character (len=8) :: cdate

      type(file_desc_t)     :: File
      type(io_desc_t)       :: iodesc2d, iodesc3d
      type(var_desc_t)      :: varid

! Info for lat, lon and time invariant variables
      ! 4 coordinate variables: TLON, TLAT, ULON, ULAT
      INTEGER (kind=int_kind), PARAMETER :: ncoord = 4

      ! 4 vertices in each grid cell
      INTEGER (kind=int_kind), PARAMETER :: nverts = 4

      ! 4 variables describe T, U grid boundaries:
      ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
      INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 4

      TYPE coord_attributes         ! netcdf coordinate attributes
        character (len=11)   :: short_name
        character (len=45)   :: long_name
        character (len=20)   :: units
      END TYPE coord_attributes

      TYPE req_attributes         ! req'd netcdf attributes
        type (coord_attributes) :: req
        character (len=20)   :: coordinates
      END TYPE req_attributes

      TYPE(req_attributes), dimension(nvar) :: var
      TYPE(coord_attributes), dimension(ncoord) :: coord_var
      TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
      CHARACTER (char_len), dimension(ncoord) :: coord_bounds

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
           work1
      real (kind=dbl_kind), dimension(nverts,nx_block,ny_block,max_blocks) :: &
           work3
      character(len=char_len_long) :: &
           filename

      real(kind=dbl_kind) :: mrss, mrss0,msize,msize0

      integer (kind=int_kind), dimension(2) ::  &
         bnd_start,bnd_length          ! dimension quantities for netCDF


      !-------------------------
      !  Test memory usage
      !-------------------------
!     call shr_get_memusage(msize,mrss)
!     call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,'ice_init_mct mrss0')
!     call shr_mpi_max(msize,msize0,MPI_COMM_ICE,'ice_init_mct msize0')
!     if(my_task == master_task) then
!       write(nu_diag,105) 'icecdf: [start of subroutine] memory_write: memory:= ',  &
!    msize0,' MB (highwater) ',mrss0,' MB (usage)'
!     endif

      if (my_task == master_task) then
         call construct_filename(ncfile(ns),'nc',ns)
         
         ! add local directory path name to ncfile
         if (write_ic) then
            ncfile(ns) = trim(incond_dir)//ncfile(ns)
         else
            ncfile(ns) = trim(history_dir)//ncfile(ns)
         endif
         filename = ncfile(ns)
      end if
      call broadcast_scalar(filename, master_task)
      
      ! create file
      File%fh=-1
      call ice_pio_init(mode='write', filename=trim(filename), File=File, &
	clobber=.true., cdf64=lcdf64)

      call ice_pio_initdecomp(iodesc=iodesc2d)
      call ice_pio_initdecomp(ndim3=nverts, inner_dim=.true., iodesc=iodesc3d)

      !-------------------------
      !  Test memory usage
      !-------------------------
!     call shr_get_memusage(msize,mrss)
!     call shr_mpi_max(mrss, mrss0, MPI_COMM_ICE,'ice_init_mct mrss0')
!     call shr_mpi_max(msize,msize0,MPI_COMM_ICE,'ice_init_mct msize0')
!     if(my_task == master_task) then
!       write(nu_diag,105) 'icecdf: [after initdecomp] memory_write: memory:= ',  &
!    msize0,' MB (highwater) ',mrss,' MB (usage)'
!     endif

      ltime = time/int(secday)

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

      if (hist_avg .and. histfreq(ns) /= '1') then
         status = pio_def_dim(File,'d2',2, boundid)
      endif
      
      status = pio_def_dim(File,'ni',nx_global,imtid)
      status = pio_def_dim(File,'nj',ny_global,jmtid)
      status = pio_def_dim(File,'time',pio_unlimited,timid)
      status = pio_def_dim(File,'nvertices',nverts,nvertexid)
     
      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

      status = pio_def_var(File,'time',pio_real,(/timid/),varid)
      status = pio_put_att(File,varid,'long_name','model time')
      
      write(cdate,'(i8.8)') idate0
      write(title,'(a,a,a,a,a,a,a)') 'days since ', &
           cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
      status = pio_put_att(File,varid,'units',title)
      
      if (days_per_year == 360) then
         status = pio_put_att(File,varid,'calendar','360_day')
      else
         status = pio_put_att(File,varid,'calendar','noleap')
      endif
      
      if (hist_avg .and. histfreq(ns) /= '1') then
         status = pio_put_att(File,varid,'bounds','time_bounds')
      endif

      !-----------------------------------------------------------------
      ! Define attributes for time bounds if hist_avg is true
      !-----------------------------------------------------------------

      if (hist_avg .and. histfreq(ns) /= '1') then
         dimid2(1) = boundid
         dimid2(2) = timid
         status = pio_def_var(File,'time_bounds',pio_real,dimid2,varid)
         status = pio_put_att(File,varid,'long_name', &
              'boundaries for time-averaging interval')
         write(cdate,'(i8.8)') idate0
         write(title,'(a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
         status = pio_put_att(File,varid,'units',title)
      endif

      !-----------------------------------------------------------------
      ! define information for required time-invariant variables
      !-----------------------------------------------------------------

      ind = 0
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLON', &
                       'T grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lont_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('TLAT', &
                       'T grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latt_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULON', &
                       'U grid center longitude', 'degrees_east')
      coord_bounds(ind) = 'lonu_bounds'
      ind = ind + 1
      coord_var(ind) = coord_attributes('ULAT', &
                       'U grid center latitude',  'degrees_north')
      coord_bounds(ind) = 'latu_bounds'

      !-----------------------------------------------------------------
      ! define information for optional time-invariant variables
      !-----------------------------------------------------------------

      var(n_tmask)%req = coord_attributes('tmask', &
                  'mask of T grid cells', 'unitless')
      var(n_tmask)%coordinates = 'TLON TLAT'
      var(n_blkmask)%req = coord_attributes('blkmask', &
                  'blockid of T grid cells', 'unitless')
      var(n_blkmask)%coordinates = 'TLON TLAT'
      var(n_tarea)%req = coord_attributes('tarea', &
                  'area of T grid cells', 'm^2')
      var(n_tarea)%coordinates = 'TLON TLAT'
      var(n_uarea)%req = coord_attributes('uarea', &
                  'area of U grid cells', 'm^2')
      var(n_uarea)%coordinates = 'ULON ULAT'
      var(n_dxt)%req = coord_attributes('dxt', &
                  'T cell width through middle', 'm')
      var(n_dxt)%coordinates = 'TLON TLAT'
      var(n_dyt)%req = coord_attributes('dyt', &
                  'T cell height through middle', 'm')
      var(n_dyt)%coordinates = 'TLON TLAT'
      var(n_dxu)%req = coord_attributes('dxu', &
                  'U cell width through middle', 'm')
      var(n_dxu)%coordinates = 'ULON ULAT'
      var(n_dyu)%req = coord_attributes('dyu', &
                  'U cell height through middle', 'm')
      var(n_dyu)%coordinates = 'ULON ULAT'
      var(n_HTN)%req = coord_attributes('HTN', &
                  'T cell width on North side','m')
      var(n_HTN)%coordinates = 'TLON TLAT'
      var(n_HTE)%req = coord_attributes('HTE', &
                  'T cell width on East side', 'm')
      var(n_HTE)%coordinates = 'TLON TLAT'
      var(n_ANGLE)%req = coord_attributes('ANGLE', &
                  'angle grid makes with latitude line on U grid', &
                  'radians')
      var(n_ANGLE)%coordinates = 'ULON ULAT'
      var(n_ANGLET)%req = coord_attributes('ANGLET', &
                  'angle grid makes with latitude line on T grid', &
                  'radians')
      var(n_ANGLET)%coordinates = 'TLON TLAT'

      ! These fields are required for CF compliance
      ! dimensions (nx,ny,nverts)
      var_nverts(n_lont_bnds) = coord_attributes('lont_bounds', &
                  'longitude boundaries of T cells', 'degrees_east')
      var_nverts(n_latt_bnds) = coord_attributes('latt_bounds', &
                  'latitude boundaries of T cells', 'degrees_north')
      var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds', &
                  'longitude boundaries of U cells', 'degrees_east')
      var_nverts(n_latu_bnds) = coord_attributes('latu_bounds', &
                  'latitude boundaries of U cells', 'degrees_north')

      !-----------------------------------------------------------------
      ! define attributes for time-invariant variables
      !-----------------------------------------------------------------

      dimid3(1) = imtid
      dimid3(2) = jmtid
      dimid3(3) = timid
      
      dimid2(1) = imtid
      dimid2(2) = jmtid

      do i = 1, ncoord
         status = pio_def_var(File, coord_var(i)%short_name, pio_real, &
              dimid2, varid)

         status = pio_put_att(File, varid,'long_name',coord_var(i)%long_name)
         status = pio_put_att(File, varid, 'units', coord_var(i)%units)
         status = pio_put_att(File, varid, 'missing_value', spval)
         status = pio_put_att(File, varid,'_FillValue',spval)
         if (coord_var(i)%short_name == 'ULAT') then
            status = pio_put_att(File,varid,'comment', &
                 'Latitude of NE corner of T grid cell')
         endif

         if (f_bounds) then
            status = pio_put_att(File, varid, 'bounds', coord_bounds(i))
         endif
      enddo

!      ! Attributes for tmask defined separately, since it has no units
!      if (igrd(n_tmask)) then
!         status = pio_def_var(File,'tmask', pio_real, dimid2, varid)
!         status = pio_put_att(File, varid, 'long_name', 'ocean grid mask')
!         status = pio_put_att(File, varid, 'coordinates', 'TLON TLAT')
!         status = pio_put_att(File, varid, 'comment', '0 = land, 1 = ocean')
!      endif

      do i = 1, nvar       ! note: n_tmask=1
         if (igrd(i)) then
            status = pio_def_var(File, var(i)%req%short_name, &
                 pio_real, dimid2, varid)
            status = pio_put_att(File, varid, 'long_name', var(i)%req%long_name)
            status = pio_put_att(File, varid, 'units', var(i)%req%units)
            status = pio_put_att(File, varid, 'coordinates', var(i)%coordinates)
            status = pio_put_att(File, varid, 'missing_value', spval)
            status = pio_put_att(File, varid,'_FillValue',spval)
         endif
      enddo
        
      ! Fields with dimensions (nverts,nx,ny)
      dimid_nverts(1) = nvertexid
      dimid_nverts(2) = imtid
      dimid_nverts(3) = jmtid
      do i = 1, nvar_verts
         if (f_bounds) then
            status = pio_def_var(File, var_nverts(i)%short_name, &
                 pio_real,dimid_nverts, varid)
            status = &
                 pio_put_att(File,varid, 'long_name', var_nverts(i)%long_name)
            status = &
                 pio_put_att(File, varid, 'units', var_nverts(i)%units)
         endif
      enddo
      
      do n=1,num_avail_hist_fields

         if (avail_hist_fields(n)%vhistfreq == histfreq(ns).or.write_ic) then

            status  = pio_def_var(File, avail_hist_fields(n)%vname, &
                 pio_real, dimid3, varid)
            status = pio_put_att(File,varid,'units', &
                 avail_hist_fields(n)%vunit)
            status = pio_put_att(File,varid, 'long_name', &
                 avail_hist_fields(n)%vdesc)

            status = pio_put_att(File,varid,'coordinates', &
                 avail_hist_fields(n)%vcoord)
            status = pio_put_att(File,varid,'cell_measures', &
                 avail_hist_fields(n)%vcellmeas)
            status = pio_put_att(File,varid,'missing_value',spval)
            status = pio_put_att(File,varid,'_FillValue',spval)

            !-----------------------------------------------------------------
            ! Append ice thickness range to aicen comments
            !-----------------------------------------------------------------

            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (i_aice > 4 .and. c_aice(1:5) == 'aicen') then
              read(c_aice(6:9), '(i3)') icategory
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            status = pio_put_att(File,varid,'comment', &
                 avail_hist_fields(n)%vcomment)

            !-----------------------------------------------------------------
            ! Add cell_methods attribute to variables if averaged
            !-----------------------------------------------------------------
            if (hist_avg .and. histfreq(ns) /= '1') then
               if (TRIM(avail_hist_fields(n)%vname)/='sig1' &
               .or.TRIM(avail_hist_fields(n)%vname)/='sig2') then
                  status = pio_put_att(File,varid,'cell_methods','time: mean')
               endif
            endif

            ! Need divu and shear as monthly means for CMIP/IPCC.
            if (histfreq(ns) == '1'     .or. .not. hist_avg      &
!                .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                 .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      & 
                 .or. n==n_trsig(ns)                             &
                 .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                 .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)    &
                 .or. n==n_FY(ns)) then
               status = pio_put_att(File,varid,'time_rep','instantaneous')
            else
               status = pio_put_att(File,varid,'time_rep','averaged')
            endif

         endif
      enddo  ! num_avail_hist_fields

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
      status = pio_put_att(File,pio_global,'title',runid)
#else
      title  = 'sea ice model output for CICE'
      status = pio_put_att(File,pio_global,'title',title)
#endif
      title = 'Diagnostic and Prognostic Variables'
      status = pio_put_att(File,pio_global,'contents',title)

      title  = 'sea ice model: Community Ice Code (CICE)'
      status = pio_put_att(File,pio_global,'source',title)

      write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
      status = pio_put_att(File,pio_global,'comment',title)

      write(title,'(a,i8)') 'File written on model date ',idate
      status = pio_put_att(File,pio_global,'comment2',title)

      write(title,'(a,i6)') 'seconds elapsed into model date: ',sec
      status = pio_put_att(File,pio_global,'comment3',title)

      title = 'CF-1.0'
      status =  &
           pio_put_att(File,pio_global,'conventions',title)
      if (my_task == master_task) then
      call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
           current_date(7:8), current_time(1:2), &
           current_time(3:4)
1000    format('This dataset was created on ', &
           a,'-',a,'-',a,' at ',a,':',a)
      end if
      call broadcast_scalar(start_time, master_task)
      status = pio_put_att(File,pio_global,'history',start_time)
      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

      status = pio_enddef(File)

      !-----------------------------------------------------------------
      ! write time variable
      !-----------------------------------------------------------------

      status = pio_inq_varid(File,'time',varid)
      status = pio_put_var(File,varid,ltime)
      
      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

      if (hist_avg .and. histfreq(ns) /= '1') then
         status = pio_inq_varid(File,'time_bounds',varid)
         time_bounds=(/time_beg(ns),time_end(ns)/)
         bnd_start=(/1,1/)
         bnd_length=(/2,1/)
         status = pio_put_var(File,varid,start=bnd_start(:),count=bnd_length(:),ival=time_bounds) 
      endif
      
      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

      do i = 1,ncoord
         call broadcast_scalar(coord_var(i)%short_name,master_task)
         
         SELECT CASE (coord_var(i)%short_name)
         CASE ('TLON')
            work1(:,:,:) = tlon(:,:,:)
            ! Convert T grid longitude from -180 -> 180 to 0 to 360
            work1 = work1*rad_to_deg + c360    ! single precision
            where (work1 > c360) work1 = work1 - c360
            where (work1 < c0  ) work1 = work1 + c360
         CASE ('TLAT')
            work1(:,:,:) = tlat(:,:,:)
            work1 = work1*rad_to_deg
         CASE ('ULON')
            work1(:,:,:) = ulon(:,:,:)
            work1 = work1*rad_to_deg
         CASE ('ULAT')
            work1(:,:,:) = ulat(:,:,:)
            work1 = work1*rad_to_deg
         END SELECT
          
         status = pio_inq_varid(File, coord_var(i)%short_name, varid)
         call pio_write_darray(File, varid, iodesc2d, &
                               work1, status, fillval=spval_dbl)
      enddo

      !-----------------------------------------------------------------
      ! write grid mask, area and rotation angle
      !-----------------------------------------------------------------

!      if (igrd(n_tmask)) then
!         work1(:,:,:) = hm(:,:,:)
!         status = pio_inq_varid(File, 'tmask', varid)
!         call pio_write_darray(File, varid, iodesc2d, &
!                               work1, status, fillval=spval_dbl)
!      endif

      do i = 1,nvar
         if (igrd(i)) then
            call broadcast_scalar(var(i)%req%short_name,master_task)

            SELECT CASE (var(i)%req%short_name)
            CASE ('tmask')
               work1(:,:,:) = hm(:,:,:)
            CASE ('blkmask')
               work1(:,:,:) = bm(:,:,:)
            CASE ('tarea')
               work1(:,:,:) = tarea(:,:,:)
            CASE ('uarea')
               work1(:,:,:) = uarea(:,:,:)
            CASE ('dxu')
               work1(:,:,:) = dxu(:,:,:)
            CASE ('dyu')
               work1(:,:,:) = dyu(:,:,:)
            CASE ('dxt')
               work1(:,:,:) = dxt(:,:,:)
            CASE ('dyt')
               work1(:,:,:) = dyt(:,:,:)
            CASE ('HTN')
               work1(:,:,:) = HTN(:,:,:)
            CASE ('HTE')
               work1(:,:,:) = HTE(:,:,:)
            CASE ('ANGLE')
               work1(:,:,:) = ANGLE(:,:,:)
            CASE ('ANGLET')
               work1(:,:,:) = ANGLET(:,:,:)
            END SELECT

            status = pio_inq_varid(File, var(i)%req%short_name, varid)
            call pio_write_darray(File, varid, iodesc2d, &
                                  work1, status, fillval=spval_dbl)
         endif
      enddo

      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
         work3(:,:,:,:) = c0
         do i = 1, nvar_verts
            call broadcast_scalar(var_nverts(i)%short_name,master_task)
            SELECT CASE (var_nverts(i)%short_name)
            CASE ('lont_bounds')
               do ivertex = 1, nverts
                  work3(ivertex,:,:,:) = lont_bounds(ivertex,:,:,:)
               enddo
            CASE ('latt_bounds')
               do ivertex = 1, nverts
                  work3(ivertex,:,:,:) = latt_bounds(ivertex,:,:,:)
               enddo
            CASE ('lonu_bounds')
               do ivertex = 1, nverts
                  work3(ivertex,:,:,:) = lonu_bounds(ivertex,:,:,:)
               enddo
            CASE ('latu_bounds')
               do ivertex = 1, nverts
                  work3(ivertex,:,:,:) = latu_bounds(ivertex,:,:,:)
               enddo
            END SELECT
            
            status = pio_inq_varid(File, var_nverts(i)%short_name, varid) 
            call pio_write_darray(File, varid, iodesc3d, &
                                  work3, status, fillval=spval_dbl)
         enddo
      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      do n=1,num_avail_hist_fields
         if (avail_hist_fields(n)%vhistfreq == histfreq(ns).or.write_ic) then
            status  = pio_inq_varid(File,avail_hist_fields(n)%vname,varid)
            if (status /= PIO_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            work1(:,:,:) = aa(:,:,n,:)
            call pio_setframe(File,varid, int(1,kind=PIO_OFFSET_KIND))
            call pio_write_darray(File, varid, iodesc2d,&
                                  work1, status, fillval=spval_dbl)
         endif
      enddo ! num_avail_hist_fields

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      call  pio_closefile(File)
      if (my_task == master_task) then
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif
      !-------------------------
      !  Test memory usage
      !-------------------------
!     call shr_get_memusage(msize,mrss)
!     call shr_mpi_max(mrss, mrss0,  MPI_COMM_ICE, 'ice_init_mct mrss0')
!     call shr_mpi_max(msize,msize0, MPI_COMM_ICE, 'ice_init_mct msize0')
!     if(my_task == master_task) then
!       write(nu_diag,105) 'icecdf: [before freedecomp] memory_write: memory:= ', &
!    msize0,' MB (highwater) ',mrss0,' MB (usage)'
!     endif
      
      ! -------------------------
      ! clean-up PIO descriptors
      ! -------------------------
      call pio_freedecomp(File,iodesc2d)
      call pio_freedecomp(File,iodesc3d)

      !-------------------------
      !  Test memory usage
      !-------------------------
!     call shr_get_memusage(msize,mrss)
!     call shr_mpi_max(mrss, mrss0,  MPI_COMM_ICE, 'ice_init_mct mrss0')
!     call shr_mpi_max(msize,msize0, MPI_COMM_ICE, 'ice_init_mct msize0')
!     if(my_task == master_task) then
!       write(nu_diag,105) 'icecdf: [end of subroutine] memory_write: memory:= ', &
!    msize0,' MB (highwater) ',mrss0,' MB (usage)'
!     endif

  105  format( A, f10.2, A, f10.2, A)

      end subroutine icecdf

!=======================================================================
!
!BOP
!
! !IROUTINE: icebin - write binary history file
! This routine writes fewer grid variables compared with the netcdf
! version, to reduce file size.  Grid variables can be obtained from
! the original grid input files.
!
! !INTERFACE:
!
      subroutine icebin(ns)
!
! !DESCRIPTION:
!
! write binary history file
!
! !REVISION HISTORY:
!
! authors:   E.C.Hunke, LANL
!
! !USES:
!
      use ice_gather_scatter
      use ice_domain_size
      use ice_constants
      use ice_restart, only: lenstr, runid
      use ice_itd, only: c_hi_range
      use ice_calendar, only: write_ic, dayyr, histfreq
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind), intent(in) :: ns

      integer (kind=int_kind) :: i,j,n,nrec,nbits
      character (char_len) :: title
      character (char_len_long) :: ncfile(max_nstrm), hdrfile

      integer (kind=int_kind) :: icategory,i_aice

      character (char_len) :: current_date,current_time
      character (len=16) :: c_aice
      logical (kind=log_kind) :: dbug

      dbug = .false.

      if (my_task == master_task) then

        call construct_filename(ncfile(ns),'da',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif
        hdrfile = trim(ncfile(ns))//'.hdr'

        !-----------------------------------------------------------------
        ! create history files
        !-----------------------------------------------------------------
        nbits = 32 ! single precision
        call ice_open(nu_history, ncfile(ns), nbits) ! direct access
        open(nu_hdr,file=hdrfile,form='formatted',status='unknown') ! ascii

!echmod call ice_write(nu_history, nrec, work, rda8 or ida4, dbug)

        title  = 'sea ice model: Community Ice Code (CICE)'
        write (nu_hdr, 999) 'source',title,' '

        write (nu_hdr, 999) 'file name contains model date',trim(ncfile(ns)),' '
#ifdef CCSMCOUPLED
        write (nu_hdr, 999) 'runid',runid,' '
#endif
        write (nu_hdr, 999) 'calendar','noleap',' '
        write (title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        write (nu_hdr, 999) 'comment',title,' '
        write (nu_hdr, 999) 'conventions','CICE',' '
        write (nu_hdr, 997) 'missing_value',spval
        write (nu_hdr, 997) '_FillValue',spval

        call date_and_time(date=current_date, time=current_time)
        write (nu_hdr,1000) current_date(1:4), current_date(5:6), &
                            current_date(7:8), current_time(1:2), &
                            current_time(3:4), current_time(5:8)
        write (nu_hdr, *  ) ' '
        write (nu_hdr, *  ) 'Grid size:'
        write (nu_hdr, 998) '  ni',nx_global
        write (nu_hdr, 998) '  nj',ny_global
   
        write (nu_hdr, *  ) 'Grid variables: (left column = nrec)'
        nrec = 1
        write (nu_hdr, 996) nrec,'tarea','area of T grid cells','m^2'
        write (nu_hdr, *  ) 'History variables: (left column = nrec)'
      endif  ! my_task = master_task

      call ice_write(nu_history, nrec, tarea, 'rda4', dbug)

      do n=1,num_avail_hist_fields

          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          nrec = nrec + 1
          if (my_task == master_task) then
            write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vdesc),trim(avail_hist_fields(n)%vunit)

            ! Append ice thickness range to aicen comments
            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (i_aice > 4 .and. c_aice(1:5) == 'aicen') then
              read(c_aice(6:9), '(i3)') icategory
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            write (nu_hdr, 995) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vcomment)

            ! Need divu and shear as monthly means for CMIP/IPCC.
            if (histfreq(ns) == '1'     .or. .not. hist_avg     &
!               .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                .or. n==n_trsig(ns)                             &
                .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)    &
                .or. n==n_FY(ns)) then
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','instantaneous'
            else
               write (nu_hdr, 996) nrec,trim(avail_hist_fields(n)%vname), &
                  'time_rep','averaged'
            endif
          endif

          call ice_write(nu_history, nrec, aa(:,:,n,:), 'rda4', dbug)

          endif
      enddo ! num_avail_hist_fields

995     format(i3,2x,a,' comment: ',a)
996     format(i3,2x,a,': ',a,',',2x,a)
997     format(a,': ',es13.6)
998     format(a,': ',i6)
999     format(a,': ',a,2x,a)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

      if (my_task == master_task) then
        close (nu_hdr)     ! header file
        close (nu_history) ! data file
        write (nu_diag,*) ' '
        write (nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif

      end subroutine icebin

      end module ice_history_write
