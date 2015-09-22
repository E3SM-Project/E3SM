!=======================================================================
!
!BOP
!
! !MODULE: ice_history_ncfbin - ice model history files writes using pio
!
! Output files: netCDF /binary
!
! !INTERFACE:
!
      module ice_history_write
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task
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
! authors:   E.C.Hunke, LANL
!            Bruce P. Briegleb, NCAR
!
! !USES:
!
#ifdef ncdf

      use ice_gather_scatter
      use ice_domain_size
      use ice_constants
      use ice_calendar, only: time, sec, idate, idate0, nyr, month, &
                              mday, write_ic, histfreq, histfreq_n, &
                              year_init, new_year, new_month, new_day, &
                              dayyr, daymo, days_per_year
      use ice_work, only: work_g1, work_gr, work_gr3
      use ice_restart, only: lenstr, runid
      use ice_domain, only: distrb_info
      use ice_itd, only: c_hi_range
      use ice_exit
      use netcdf
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind), intent(in) :: ns

      integer (kind=int_kind) :: i,j,n, &
         ncid,status,imtid,jmtid,timid,varid, &
         length,nvertexid,ivertex
      integer (kind=int_kind), dimension(3) :: dimid
      integer (kind=int_kind), dimension(3) :: dimid_nverts
      real (kind=real_kind) :: ltime
      character (char_len) :: title
      character (char_len_long) :: ncfile(nstreams)

      integer (kind=int_kind) :: iyear, imonth, iday
      integer (kind=int_kind) :: icategory,ind,i_aice,boundid

      character (char_len) :: start_time,current_date,current_time
      character (len=16) :: c_aice
      character (len=8) :: cdate

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

      if (my_task == master_task) then

        ltime=time/int(secday)

        call construct_filename(ncfile(ns),'nc',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
          ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
          ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif

        ! create file
#ifdef _HIRES
        status = nf90_create(ncfile(ns), NF90_64BIT_OFFSET, ncid)
#else
        status = nf90_create(ncfile(ns), nf90_clobber, ncid)
#endif
        if (status /= nf90_noerr) call abort_ice( &
           'ice: Error creating history ncfile '//ncfile(ns))

      !-----------------------------------------------------------------
      ! define dimensions
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_def_dim(ncid,'d2',2,boundid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error defining dim d2')
        endif

        status = nf90_def_dim(ncid,'ni',nx_global,imtid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim ni')

        status = nf90_def_dim(ncid,'nj',ny_global,jmtid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim nj')

        status = nf90_def_dim(ncid,'time',NF90_UNLIMITED,timid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim time')

        status = nf90_def_dim(ncid,'nvertices',nverts,nvertexid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining dim nverts')

      !-----------------------------------------------------------------
      ! define coordinate variables
      !-----------------------------------------------------------------

        status = nf90_def_var(ncid,'time',nf90_float,timid,varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error defining var time')

        status = nf90_put_att(ncid,varid,'long_name','model time')
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time long_name')

        write(cdate,'(i8.8)') idate0
        write(title,'(a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
        status = nf90_put_att(ncid,varid,'units',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time units')

        if (days_per_year == 360) then
           status = nf90_put_att(ncid,varid,'calendar','360_day')
           if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        else
           status = nf90_put_att(ncid,varid,'calendar','noleap')
           if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        endif

        if (hist_avg) then
          status = nf90_put_att(ncid,varid,'bounds','time_bounds')
          if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time bounds')
        endif

      !-----------------------------------------------------------------
      ! Define attributes for time bounds if hist_avg is true
      !-----------------------------------------------------------------

        if (hist_avg) then
          dimid(1) = boundid
          dimid(2) = timid
          status = nf90_def_var(ncid,'time_bounds',nf90_float,dimid(1:2),varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error defining var time_bounds')
          status = nf90_put_att(ncid,varid,'long_name', &
                                'boundaries for time-averaging interval')
          if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds long_name')
          write(cdate,'(i8.8)') idate0
          write(title,'(a,a,a,a,a,a,a)') 'days since ', &
                cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
          status = nf90_put_att(ncid,varid,'units',title)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds units')
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

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid

        do i = 1, ncoord
          status = nf90_def_var(ncid, coord_var(i)%short_name, nf90_float, &
                                dimid(1:2), varid)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining short_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'long_name',coord_var(i)%long_name)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining long_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid, varid, 'units', coord_var(i)%units)
          if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//coord_var(i)%short_name)
          if (coord_var(i)%short_name == 'ULAT') then
             status = nf90_put_att(ncid,varid,'comment', &
                  'Latitude of NE corner of T grid cell')
             if (status /= nf90_noerr) call abort_ice( &
                   'Error defining comment for '//coord_var(i)%short_name)
          endif
          if (f_bounds) then
             status = nf90_put_att(ncid, varid, 'bounds', coord_bounds(i))
             if (status /= nf90_noerr) call abort_ice( &
                 'Error defining bounds for '//coord_var(i)%short_name)
          endif
        enddo

        ! Attributes for tmask defined separately, since it has no units
        if (igrd(n_tmask)) then
           status = nf90_def_var(ncid, 'tmask', nf90_float, dimid(1:2), varid)
           if (status /= nf90_noerr) call abort_ice( &
                         'ice: Error defining var tmask')
           status = nf90_put_att(ncid,varid, 'long_name', 'ocean grid mask')
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask long_name')
           status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask units')
           status = nf90_put_att(ncid,varid,'comment', '0 = land, 1 = ocean')
           if (status /= nf90_noerr) call abort_ice('ice Error: tmask comment')
        endif

        do i = 2, nvar       ! note: n_tmask=1
          if (igrd(i)) then
             status = nf90_def_var(ncid, var(i)%req%short_name, &
                                   nf90_float, dimid(1:2), varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//var(i)%req%short_name)
             status = nf90_put_att(ncid,varid, 'long_name', var(i)%req%long_name)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'units', var(i)%req%units)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//var(i)%req%short_name)
             status = nf90_put_att(ncid, varid, 'coordinates', var(i)%coordinates)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining coordinates for '//var(i)%req%short_name)
          endif
        enddo

        ! Fields with dimensions (nverts,nx,ny)
        dimid_nverts(1) = nvertexid
        dimid_nverts(2) = imtid
        dimid_nverts(3) = jmtid
        do i = 1, nvar_verts
          if (f_bounds) then
             status = nf90_def_var(ncid, var_nverts(i)%short_name, &
                                   nf90_float,dimid_nverts, varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining variable '//var_nverts(i)%short_name)
             status = &
             nf90_put_att(ncid,varid, 'long_name', var_nverts(i)%long_name)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining long_name for '//var_nverts(i)%short_name)
             status = &
             nf90_put_att(ncid, varid, 'units', var_nverts(i)%units)
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//var_nverts(i)%short_name)
          endif
        enddo

        do n=1,num_avail_hist_fields

          if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

            status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                         nf90_float, dimid, varid)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining variable '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'units', &
                        avail_hist_fields(n)%vunit)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining units for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid, 'long_name', &
                        avail_hist_fields(n)%vdesc)
                   
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining long_name for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'coordinates', &
                        avail_hist_fields(n)%vcoord)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining coordinates for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'cell_measures', &
                        avail_hist_fields(n)%vcellmeas)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining cell measures for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining mising_value for '//avail_hist_fields(n)%vname)
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining _FillValue for '//avail_hist_fields(n)%vname)

      !-----------------------------------------------------------------
      ! Append ice thickness range to aicen comments
      !-----------------------------------------------------------------
            c_aice = TRIM(avail_hist_fields(n)%vname)
            i_aice = lenstr(c_aice)
            if (c_aice(1:4) == 'aice' .and. i_aice > 4 ) then
              if (i_aice == 5) then         ! categories 1-9
                read(c_aice(i_aice:i_aice), '(i1)') icategory
              else                          ! categories > 9
                read(c_aice(i_aice-1:i_aice), '(i2)') icategory
              endif
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            status = nf90_put_att(ncid,varid,'comment', &
                        avail_hist_fields(n)%vcomment)
            if (status /= nf90_noerr) call abort_ice( &
               'Error defining comment for '//avail_hist_fields(n)%vname)
      !-----------------------------------------------------------------
      ! Add cell_methods attribute to variables if averaged
      !-----------------------------------------------------------------
            if (hist_avg) then
              if (TRIM(avail_hist_fields(n)%vname)/='sig1' &
              .or.TRIM(avail_hist_fields(n)%vname)/='sig2') then
                status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                if (status /= nf90_noerr) call abort_ice( &
                 'Error defining cell methods for '//avail_hist_fields(n)%vname)
              endif
            endif

            if (histfreq(ns) == '1'     .or. .not. hist_avg &
                .or. n==n_divu      .or. n==n_shear     &  ! snapshots
                .or. n==n_sig1      .or. n==n_sig2 .or. n==n_trsig &
                .or. n==n_mlt_onset .or. n==n_frz_onset &
                .or. n==n_hisnap    .or. n==n_aisnap .or. n==n_FY) then
               status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
            else
               status = nf90_put_att(ncid,varid,'time_rep','averaged')
            endif

           endif 
        enddo  ! num_avail_hist_fields

      !-----------------------------------------------------------------
      ! global attributes
      !-----------------------------------------------------------------
      ! ... the user should change these to something useful ...
      !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
        status = nf90_put_att(ncid,nf90_global,'title',runid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#else
        title  = 'sea ice model output for CICE'
        status = nf90_put_att(ncid,nf90_global,'title',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#endif
        title = 'Diagnostic and Prognostic Variables'
        status = nf90_put_att(ncid,nf90_global,'contents',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute contents')

        title  = 'sea ice model: Community Ice Code (CICE)'
        status = nf90_put_att(ncid,nf90_global,'source',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute source')

        write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        status = nf90_put_att(ncid,nf90_global,'comment',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute comment')

        write(title,'(a,i8)') 'File written on model date ',idate
        status = nf90_put_att(ncid,nf90_global,'comment2',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date1')

        write(title,'(a,i6)') 'seconds elapsed into model date: ',sec
        status = nf90_put_att(ncid,nf90_global,'comment3',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date2')

        title = 'CF-1.0'
        status =  &
             nf90_put_att(ncid,nf90_global,'conventions',title)
        if (status /= nf90_noerr) call abort_ice( &
             'Error in global attribute conventions')

        call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
                               current_date(7:8), current_time(1:2), &
                               current_time(3:4), current_time(5:8)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

        status = nf90_put_att(ncid,nf90_global,'history',start_time)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute history')

      !-----------------------------------------------------------------
      ! end define mode
      !-----------------------------------------------------------------

        status = nf90_enddef(ncid)
        if (status /= nf90_noerr) call abort_ice('ice: Error in nf90_enddef')

      !-----------------------------------------------------------------
      ! write time variable
      !-----------------------------------------------------------------

        status = nf90_inq_varid(ncid,'time',varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error getting time varid')
        status = nf90_put_var(ncid,varid,ltime)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error writing time variable')

      !-----------------------------------------------------------------
      ! write time_bounds info
      !-----------------------------------------------------------------

        if (hist_avg) then
          status = nf90_inq_varid(ncid,'time_bounds',varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error getting time_bounds id')
          status = nf90_put_var(ncid,varid,time_beg(ns),start=(/1/))
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing time_beg')
          status = nf90_put_var(ncid,varid,time_end(ns),start=(/2/))
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing time_end')
        endif

      endif                     ! master_task

      if (my_task==master_task) then
         allocate(work_g1(nx_global,ny_global))
         allocate(work_gr(nx_global,ny_global))
      else
         allocate(work_gr(1,1))   ! to save memory
         allocate(work_g1(1,1))
      endif

      work_g1(:,:) = c0

      !-----------------------------------------------------------------
      ! write coordinate variables
      !-----------------------------------------------------------------

        do i = 1,ncoord
          call broadcast_scalar(coord_var(i)%short_name,master_task)
          SELECT CASE (coord_var(i)%short_name)
            CASE ('TLON')
              call gather_global(work_g1,TLON,master_task,distrb_info)
              if (my_task == master_task) then
              ! Convert T grid longitude from -180 -> 180 to 0 to 360
                 work_gr = work_g1*rad_to_deg + c360    ! single precision
                 where (work_gr > c360) work_gr = work_gr - c360
                 where (work_gr < c0 )  work_gr = work_gr + c360
              endif
            CASE ('TLAT')
              call gather_global(work_g1,TLAT,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
            CASE ('ULON')
              call gather_global(work_g1,ULON,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
            CASE ('ULAT')
              call gather_global(work_g1,ULAT,master_task,distrb_info)
              if (my_task == master_task) work_gr = work_g1*rad_to_deg
          END SELECT
          
          if (my_task == master_task) then
             status = nf90_inq_varid(ncid, coord_var(i)%short_name, varid)
             if (status /= nf90_noerr) call abort_ice( &
                  'ice: Error getting varid for '//coord_var(i)%short_name)
             status = nf90_put_var(ncid,varid,work_gr)
             if (status /= nf90_noerr) call abort_ice( &
                           'ice: Error writing'//coord_var(i)%short_name)
          endif
        enddo

      !-----------------------------------------------------------------
      ! write grid mask, area and rotation angle
      !-----------------------------------------------------------------

      if (igrd(n_tmask)) then
      call gather_global(work_g1, hm, master_task, distrb_info)
      if (my_task == master_task) then
        work_gr=work_g1
        status = nf90_inq_varid(ncid, 'tmask', varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error getting varid for tmask')
        status = nf90_put_var(ncid,varid,work_gr)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error writing variable tmask')
      endif
      endif

      do i = 2,nvar
        if (igrd(i)) then
        call broadcast_scalar(var(i)%req%short_name,master_task)
        SELECT CASE (var(i)%req%short_name)
          CASE ('tarea')
            call gather_global(work_g1, tarea, master_task, distrb_info)
          CASE ('uarea')
            call gather_global(work_g1, uarea, master_task, distrb_info)
          CASE ('dxu')
            call gather_global(work_g1,   dxu, master_task, distrb_info)
          CASE ('dyu')
            call gather_global(work_g1,   dyu, master_task, distrb_info)
          CASE ('dxt')
            call gather_global(work_g1,   dxt, master_task, distrb_info)
          CASE ('dyt')
            call gather_global(work_g1,   dyt, master_task, distrb_info)
          CASE ('HTN')
            call gather_global(work_g1,   HTN, master_task, distrb_info)
          CASE ('HTE')
            call gather_global(work_g1,   HTE, master_task, distrb_info)
          CASE ('ANGLE')
            call gather_global(work_g1, ANGLE, master_task, distrb_info)
          CASE ('ANGLET')
            call gather_global(work_g1, ANGLET,master_task, distrb_info)
        END SELECT

        if (my_task == master_task) then
          work_gr=work_g1
          status = nf90_inq_varid(ncid, var(i)%req%short_name, varid)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error getting varid for '//var(i)%req%short_name)
          status = nf90_put_var(ncid,varid,work_gr)
          if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing variable '//var(i)%req%short_name)
        endif
        endif
      enddo

      deallocate(work_gr)
      !----------------------------------------------------------------
      ! Write coordinates of grid box vertices
      !----------------------------------------------------------------

      if (f_bounds) then
      if (my_task==master_task) then
         allocate(work_gr3(nverts,nx_global,ny_global))
      else
         allocate(work_gr3(1,1,1))   ! to save memory
      endif

      work_gr3(:,:,:) = c0
      work1   (:,:,:) = c0

      do i = 1, nvar_verts
        call broadcast_scalar(var_nverts(i)%short_name,master_task)
        SELECT CASE (var_nverts(i)%short_name)
        CASE ('lont_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = lont_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latt_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latt_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('lonu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = lonu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        CASE ('latu_bounds')
        do ivertex = 1, nverts
           work1(:,:,:) = latu_bounds(ivertex,:,:,:)
           call gather_global(work_g1, work1, master_task, distrb_info)
           if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
        enddo
        END SELECT

        if (my_task == master_task) then
          status = nf90_inq_varid(ncid, var_nverts(i)%short_name, varid)
          if (status /= nf90_noerr) call abort_ice( &
             'ice: Error getting varid for '//var_nverts(i)%short_name)
          status = nf90_put_var(ncid,varid,work_gr3)
          if (status /= nf90_noerr) call abort_ice( &
             'ice: Error writing variable '//var_nverts(i)%short_name)
        endif
      enddo
      deallocate(work_gr3)
      endif

      !-----------------------------------------------------------------
      ! write variable data
      !-----------------------------------------------------------------

      if (my_task==master_task) then
         allocate(work_gr3(nx_global,ny_global,1))
      else
         allocate(work_gr3(1,1,1))   ! to save memory
      endif

      do n=1,num_avail_hist_fields

         if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) then

          call gather_global(work_g1, aa(:,:,n,:), &
                             master_task, distrb_info)
          if (my_task == master_task) then
            work_gr3(:,:,1) = work_g1(:,:)
            status  = nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid)
            if (status /= nf90_noerr) call abort_ice( &
               'ice: Error getting varid for '//avail_hist_fields(n)%vname)
            status  = nf90_put_var(ncid,varid,work_gr3, &
                                   count=(/nx_global,ny_global,1/))
            if (status /= nf90_noerr) call abort_ice( &
               'ice: Error writing variable '//avail_hist_fields(n)%vname)
          endif

         endif
      enddo ! num_avail_hist_fields

      deallocate(work_gr3)
      deallocate(work_g1)

      !-----------------------------------------------------------------
      ! close output dataset
      !-----------------------------------------------------------------

      if (my_task == master_task) then
         status = nf90_close(ncid)
         if (status /= nf90_noerr) call abort_ice( &
                       'ice: Error closing netCDF history file')
         write(nu_diag,*) ' '
         write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
      endif
#endif

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
      character (char_len_long) :: ncfile(nstreams), hdrfile

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
            if (c_aice(1:4) == 'aice' .and. i_aice > 4 ) then
              if (i_aice == 5) then         ! categories 1-9
                read(c_aice(i_aice:i_aice), '(i1)') icategory
              else                          ! categories > 9
                read(c_aice(i_aice-1:i_aice), '(i2)') icategory
              endif
              avail_hist_fields(n)%vcomment = &
                 'Ice range: '//c_hi_range(icategory)
            endif
            write (nu_hdr, 995) nrec,trim(avail_hist_fields(n)%vname), &
               trim(avail_hist_fields(n)%vcomment)

            if (histfreq(ns) == '1'     .or. .not. hist_avg &
                .or. n==n_divu      .or. n==n_shear     &  ! snapshots
                .or. n==n_sig1      .or. n==n_sig2 .or. n==n_trsig &
                .or. n==n_mlt_onset .or. n==n_frz_onset &
                .or. n==n_hisnap    .or. n==n_aisnap) then
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
