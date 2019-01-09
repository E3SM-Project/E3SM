module native_mapping
!
!  Create mapping files using the SE basis functions.  This module looks for the namelist 'native_mapping' in
!  file NLFileName (usually atm_in) and reads from it a list of up to maxoutgrids grid description files
!  It then creates a grid mapping file from the currently defined SE grid to the grid described in each file
!  using the SE basis functions.   The output mapping file name is generated based on the SE model resolution
!  and the input grid file name and ends in '_date_native.nc'
!
  use cam_logfile,       only : iulog
  use shr_kind_mod,      only : r8 => shr_kind_r8, shr_kind_cl
  use shr_const_mod,     only : pi=>shr_const_pi
  use cam_abortutils,    only : endrun
  use spmd_utils,        only : iam, masterproc, mpi_character, mpi_logical, mpi_integer, mpi_max, &
                                mpicom, mstrid=>masterprocid

  implicit none
  private
  public :: native_mapping_readnl, create_native_mapping_files, do_native_mapping

  integer, parameter :: maxoutgrids=5
  character(len=shr_kind_cl) :: native_mapping_outgrids(maxoutgrids)
  logical, protected :: do_native_mapping

!=============================================================================================
contains
!=============================================================================================

subroutine native_mapping_readnl(NLFileName)

   use units,          only : getunit, freeunit
   use namelist_utils, only : find_group_name

   character(len=*), intent(in) :: NLFileName

   character(len=shr_kind_cl) :: mappingfile, fname

   namelist /native_mapping_nl/  native_mapping_outgrids
   integer :: nf, unitn, ierr
   logical :: exist
   character(len=*), parameter ::  sub="native_mapping_readnl"
   !-----------------------------------------------------------------------------

   do_native_mapping=.false.

   do nf=1,maxoutgrids
      native_mapping_outgrids(nf)=''
   enddo

   if(masterproc) then
      exist=.true.
      write(iulog,*) sub//': Check for native_mapping_nl namelist in ',trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      call find_group_name(unitn, 'native_mapping_nl', status=ierr)
      if(ierr/=0) then
         write(iulog,*) sub//': No native_mapping_nl namelist found'
         exist=.false.
      end if
      if(exist) then
         read(unitn, native_mapping_nl, iostat=ierr)
         if(ierr/=0) then
            call endrun(sub//': namelist read returns an error condition for native_mapping_nl')
         end if
         if(len_trim(native_mapping_outgrids(1))==0) exist=.false.
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(exist, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: exist")

   if(.not. exist) return

   call mpi_bcast(native_mapping_outgrids, maxoutgrids*shr_kind_cl, mpi_character, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: native_mapping_outgrids")

   do_native_mapping=.true.

end subroutine native_mapping_readnl

!=============================================================================================

subroutine create_native_mapping_files(par, elem, maptype, ncol, clat, clon, areaa)

    use parallel_mod, only : parallel_t, global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use cam_pio_utils, only : cam_pio_openfile, cam_pio_createfile
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t, config_thread_region
    use pio, only : pio_noerr, pio_openfile, pio_createfile, pio_closefile, &
         pio_get_var, pio_put_var, pio_write_darray,pio_int, pio_double, &
         pio_def_var, pio_put_att, pio_global, file_desc_t, var_desc_t, &
         io_desc_t, pio_internal_error,pio_inq_dimlen, pio_inq_varid, &
         pio_get_att, pio_enddef, pio_bcast_error,pio_internal_error, &
         pio_def_dim, pio_inq_dimid, pio_seterrorhandling, pio_initdecomp
    use quadrature_mod, only : quadrature_t, gauss, gausslobatto
    use interpolate_mod, only : interpdata_t, cube_facepoint_ne, interpolate_scalar, set_interp_parameter, interp_init, &
         get_interp_parameter
    use coordinate_systems_mod, only : spherical_polar_t, cartesian2d_t
    use dimensions_mod, only : nelemd, ne, np, npsq, nelem
    use reduction_mod, only : ParallelMin,ParallelMax
    use cube_mod, only : convert_gbl_index
    use infnan, only : isnan
    use dof_mod, only : CreateMetaData
    use thread_mod,     only: omp_get_thread_num
    use datetime_mod, only: datetime


    use cam_history_support, only : fillvalue


    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in) :: elem(:)
    character(len=*), intent(in) :: maptype
    integer,          intent(in) :: ncol
    real(r8),         intent(in) :: clat(ncol)
    real(r8),         intent(in) :: clon(ncol)
    real(r8),         intent(in) :: areaa(ncol)

    character(len=shr_kind_cl) :: mappingfile, fname


    type(hybrid_t) :: hybrid
    logical :: exist

    type (spherical_polar_t) :: sphere
    type(file_desc_t) :: ogfile, agfile
    type (interpdata_t)  :: interpdata(nelemd)
    integer :: ierr, dimid, npts, vid
    real(r8), allocatable :: lat(:), lon(:)
    integer :: i, ii, ie2, je2, ie, je, face_no, face_no2, k, j, n, ngrid, tpts, nf, number
    real(r8) :: countx, count_max, count_total
    integer :: fdofp(np,np,nelemd)
    type (cartesian2D_t) :: cart
    real(r8) :: f(np,np)
    real(r8), allocatable :: h(:), h1d(:)
    integer, allocatable :: grid_imask(:), row(:), col(:), ldof(:), dg_dims(:)
    integer :: ns_dim, cnt, na_dim, nb_dim, sg_dim, dg_dim
    type(var_desc_t) :: rowid, colid, sid, xca_id, yca_id, xcb_id, ycb_id, maskb_id, maska_id
    type(var_desc_t) :: areaA_id, areaB_id, dg_id, sg_id
    type(io_desc_t) :: iodesci, iodescd
    character(len=12) :: unit_str
    real(r8), allocatable :: areaB(:)
    integer :: cntperelem_in(nelem),  cntperelem_out(nelem)
    integer :: ithr, dg_rank, substr1, substr2

    type(interpdata_t), pointer :: mapping_interpolate(:)
    character(len=8) :: cdate, ctime
    integer :: olditype, oldnlat, oldnlon, itype



    if(.not. do_native_mapping) return

    if (maptype=='native') then
       itype=0
    else if (maptype=='bilin') then
       itype=1
    else
       call endrun('bad interp_type')
    endif




    if(iam > par%nprocs) then
       ! The special case of npes_se < npes_cam is not worth dealing with here
       call endrun('Native mapping code requires npes_se==npes_cam')
    end if


    call interp_init()


    oldnlon = get_interp_parameter('nlon')
    oldnlat = get_interp_parameter('nlat')
    olditype = get_interp_parameter('itype')

    call datetime(cdate, ctime)

    do nf=1,maxoutgrids
       fname = native_mapping_outgrids(nf)
       if(masterproc) then
          write(iulog,*) 'looking for target grid = ',trim(fname)
       endif
       if(len_trim(fname)==0) cycle
       inquire(file=fname,exist=exist)
       if(.not. exist) then
          write(iulog,*) 'WARNING: Could not find or open grid file ',fname
          cycle
       end if
       if(masterproc) then
          write(iulog,*) 'Creating ',trim(maptype),' mapping to grid ',fname
       endif
       call cam_pio_openfile( ogfile, fname, 0)

       ierr = pio_inq_dimid( ogfile, 'grid_size', dimid)
       ierr = pio_inq_dimlen( ogfile, dimid, npts)
       allocate(lat(npts), lon(npts), grid_imask(npts), areab(npts))

       ierr = pio_inq_dimid( ogfile, 'grid_rank', dimid)
       ierr = pio_inq_dimlen(ogfile, dimid, dg_rank)
       allocate(dg_dims(dg_rank))
       ierr = pio_inq_varid( ogfile, 'grid_dims', vid)
       ierr = pio_get_var( ogfile, vid, dg_dims)


       ierr = pio_inq_varid( ogfile, 'grid_center_lat', vid)
       ierr = pio_get_var(ogfile, vid, lat)
       ierr = pio_get_att(ogfile, vid, 'units', unit_str)

       ierr = pio_inq_varid( ogfile, 'grid_center_lon', vid)
       ierr = pio_get_var(ogfile, vid, lon)

       call pio_seterrorhandling(ogfile, PIO_BCAST_ERROR)
       ierr = pio_inq_varid( ogfile, 'grid_area', vid)
       call pio_seterrorhandling(ogfile, PIO_INTERNAL_ERROR)
       if(ierr == PIO_NOERR) then
          ierr = pio_get_var(ogfile, vid, areaB)
       else
          areaB=fillvalue
       end if

       if(unit_str .eq. 'degrees') then
          lat = lat * pi/180_r8
          lon = lon * pi/180_r8
       end if

       ierr = pio_inq_varid( ogfile, 'grid_imask', vid)
       ierr = pio_get_var(ogfile, vid, grid_imask)
       call pio_closefile(ogfile)

       do ie=1,nelemd
          interpdata(ie)%n_interp=0
       end do

       call set_interp_parameter('itype',itype)   ! itype=0 native, 1 for bilinear
       if(lon(1)==lon(2)) then
          call set_interp_parameter('nlon',dg_dims(1))
          call set_interp_parameter('nlat',dg_dims(2))
       else
          call set_interp_parameter('nlon',dg_dims(2))
          call set_interp_parameter('nlat',dg_dims(1))
       end if





!       call setup_latlon_interp(elem, cam_interpolate, hybrid, 1, nelemd)
       ! go through once, counting the number of points on each element

       sphere%r=1
       do i=1,npts
          if(grid_imask(i)==1) then
             sphere%lat=lat(i)
             sphere%lon=lon(i)
             call cube_facepoint_ne(sphere, ne, cart, number)  ! new interface
             if (number /= -1) then
                do ii=1,nelemd
                   if (number == elem(ii)%vertex%number) then
                      interpdata(ii)%n_interp = interpdata(ii)%n_interp + 1
                      exit
                   endif
                enddo
             endif


             if(masterproc) then
                if(mod(i,npts/10).eq.1) then
                   print *,'finished point ',i,' of ',npts
                endif
             end if
          end if
       enddo

       hybrid = config_thread_region(par,'serial')
!       ithr=omp_get_thread_num()
!       hybrid = hybrid_create(par,ithr,1)



       ! check if every point in interpolation grid was claimed by an element:
       countx=sum(interpdata(1:nelemd)%n_interp)
       global_shared_buf(1,1) = countx
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm, nsize=1)
       count_total = global_shared_sum(1)
       tpts = sum(grid_imask)
       if (count_total /= tpts ) then
          write(iulog,*)__FILE__,__LINE__,iam, count_total, tpts, npts
          call endrun('Error setting up interpolation grid count_total<>npts')
       endif

       countx=maxval(interpdata(1:nelemd)%n_interp)
       count_max = ParallelMax(countx,hybrid)

       if (masterproc) then
          write(iulog,'(a,f8.1)') 'Average number of interpolation points per element: ',count_total/real(6*ne*ne)
          write(iulog,'(a,f8.0)') 'Maximum number of interpolation points on any element: ',count_max
       endif


       ! allocate storage
       do ii=1,nelemd
          ngrid = interpdata(ii)%n_interp
          allocate(interpdata(ii)%interp_xy( ngrid ) )
          allocate(interpdata(ii)%ilat( ngrid ) )
          allocate(interpdata(ii)%ilon( ngrid ) )
          interpdata(ii)%n_interp=0  ! reset counter
       enddo

       ! now go through the list again, adding the coordinates
       ! if this turns out to be slow, then it can be done in the loop above
       ! but we have to allocate and possibly resize the interp_xy() array.
       do i=1,npts
          if(grid_imask(i)==1) then
             sphere%lat=lat(i)
             sphere%lon=lon(i)
             call cube_facepoint_ne(sphere, ne, cart, number)  ! new interface
             if (number /= -1) then
                do ii=1,nelemd
                   if (number == elem(ii)%vertex%number) then
                      ngrid = interpdata(ii)%n_interp + 1
                      interpdata(ii)%n_interp = ngrid
                      interpdata(ii)%interp_xy( ngrid ) = cart
                      interpdata(ii)%ilon( ngrid ) = i
                      interpdata(ii)%ilat( ngrid ) = i
                   endif
                enddo
             endif
          end if
       end do


       allocate(h(int(countx)))
       allocate(h1d(int(countx)*npsq*nelemd))
       allocate(row(int(countx)*npsq*nelemd))
       allocate(col(int(countx)*npsq*nelemd))

       row = 0
       col = 0

       ngrid=0
       cntperelem_in=0
       call CreateMetaData(hybrid%par, elem, fdofp=fdofp)

       do ie=1,nelemd
          ii=0
          do j=1,np
             do i=1,np
                ii=ii+1
                f = 0.0_R8
                f(i,j) = 1.0_R8
                h = 0
                call interpolate_scalar(interpdata(ie), f, np, 0, h(:))

                do n=1,interpdata(ie)%n_interp
                   if(any(isnan(h ))) then

                      call endrun('nan generated')
                   end if
                   if(h(n)/=0) then
                      ngrid=ngrid+1
                      h1d(ngrid) = h(n)
                      row(ngrid) = interpdata(ie)%ilon(n)
                      col(ngrid) =  fdofp(i,j,ie)
                      cntperelem_in(elem(ie)%Globalid)=cntperelem_in(elem(ie)%Globalid)+1
                   end if
                enddo

             enddo
          end do
       end do

       countx=ngrid
       global_shared_buf(1,1) = countx
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm, nsize=1)
       count_total = global_shared_sum(1)


       call mpi_allreduce(cntperelem_in, cntperelem_out, nelem, MPI_INTEGER, MPI_MAX, par%comm, ierr)


       allocate(ldof(ngrid))
       ldof = 0
       ii=1
       do ie=1,nelemd
          if(elem(ie)%GlobalID==1) then
             cnt = 0
          else
             cnt = sum(cntperelem_out(1:elem(ie)%globalid-1))
          endif
          do i=1,cntperelem_out(elem(ie)%globalid)
             ldof(ii) = cnt+i
             ii=ii+1
          end do
       end do

       deallocate(h)

       ngrid = int(count_total)

       substr1 = index(fname,'/',BACK=.true.)
       substr2 = index(fname,'.nc',BACK=.true.)

       if(ne<100) then
          write(mappingfile,113) ne,np,fname(substr1+1:substr2-1),trim(maptype),cdate(7:8),cdate(1:2),cdate(4:5)
       else if(ne<1000) then
          write(mappingfile,114) ne,np,fname(substr1+1:substr2-1),trim(maptype),cdate(7:8),cdate(1:2),cdate(4:5)
       else
          write(mappingfile,115) ne,np,fname(substr1+1:substr2-1),trim(maptype),cdate(7:8),cdate(1:2),cdate(4:5)
       end if

113    format('map_ne',i2.2,'np',i1,'_to_',a,'_',a,'_',3a2,'.nc')
114    format('map_ne',i3.3,'np',i1,'_to_',a,'_',a,'_',3a2,'.nc')
115    format('map_ne',i4.4,'np',i1,'_to_',a,'_',a,'_',3a2,'.nc')

       call cam_pio_createfile( ogfile,mappingfile , 0)

       ierr = pio_def_dim( ogfile, 'n_a', ncol, na_dim)
       ierr = pio_def_dim( ogfile, 'n_b', npts, nb_dim)
       ierr = pio_def_dim( ogfile, 'n_s', ngrid, ns_dim)

       ierr = pio_def_dim( ogfile, 'src_grid_rank', 1, sg_dim)
       ierr = pio_def_var( ogfile, 'src_grid_dims',pio_int, (/sg_dim/),sg_id)

       ierr = pio_def_dim( ogfile, 'dst_grid_rank',dg_rank, dg_dim)
       ierr = pio_def_var( ogfile, 'dst_grid_dims',pio_int, (/dg_dim/),dg_id)





       ierr = pio_def_var( ogfile, 'col', pio_int, (/ns_dim/), colid)
       ierr = pio_def_var( ogfile, 'row', pio_int, (/ns_dim/), rowid)
       ierr = pio_def_var( ogfile, 'S', pio_double, (/ns_dim/), sid)

       ierr = pio_def_var( ogfile, 'xc_a', pio_double, (/na_dim/), xca_id)
       ierr = pio_def_var( ogfile, 'yc_a', pio_double, (/na_dim/), yca_id)

       ierr = pio_def_var( ogfile, 'xc_b', pio_double, (/nb_dim/), xcb_id)
       ierr = pio_def_var( ogfile, 'yc_b', pio_double, (/nb_dim/), ycb_id)

       ierr = pio_def_var( ogfile, 'area_a', pio_double, (/na_dim/), areaA_id)
       ierr = pio_def_var( ogfile, 'area_b', pio_double, (/nb_dim/), areaB_id)
       ierr = pio_put_att( ogfile, areaB_id, '_FillValue',fillvalue)

       ierr = pio_def_var( ogfile, 'mask_a', pio_int, (/na_dim/), maska_id)
       ierr = pio_def_var( ogfile, 'mask_b', pio_int, (/nb_dim/), maskb_id)



       ierr = pio_put_att( ogfile, xca_id, 'units','radians')
       ierr = pio_put_att( ogfile, yca_id, 'units','radians')
       ierr = pio_put_att( ogfile, xcb_id, 'units','radians')
       ierr = pio_put_att( ogfile, ycb_id, 'units','radians')

       ierr = pio_put_att( ogfile, PIO_GLOBAL, 'title', 'SE NATIVE Regridding Weights')
       ierr = pio_put_att( ogfile, PIO_GLOBAL, 'normalization', 'none')
       if (itype==0 ) then
          ierr = pio_put_att( ogfile, PIO_GLOBAL, 'map_method', 'Spectral-Element remapping')
       else if (itype==1) then
          ierr = pio_put_att( ogfile, PIO_GLOBAL, 'map_method', 'Bilinear remapping')
       endif
       ierr = pio_put_att( ogfile, PIO_GLOBAL, 'conventions', 'NCAR-CSM')

       ierr = pio_put_att( ogfile, PIO_GLOBAL, 'grid_file_out', fname  )
       ierr = pio_put_att( ogfile, PIO_GLOBAL, 'grid_file_atm', 'none - model generated')


       ierr = pio_enddef ( ogfile )

       ierr = pio_put_var(ogfile, sg_id, ncol)
       ierr = pio_put_var(ogfile, dg_id, dg_dims(1:dg_rank))


       call pio_initdecomp( ogfile%iosystem, pio_int, (/ngrid/), ldof, iodesci)
       call pio_initdecomp( ogfile%iosystem, pio_double, (/ngrid/), ldof, iodescd)

       call pio_write_darray(ogfile, colid, iodesci, col, ierr)
       call pio_write_darray(ogfile, rowid, iodesci, row, ierr)
       call pio_write_darray(ogfile, sid, iodescd, h1d, ierr)


       ierr = pio_put_var(ogfile, xcb_id, lon)
       ierr = pio_put_var(ogfile, ycb_id, lat)

       ierr = pio_put_var(ogfile, xca_id, clon)
       ierr = pio_put_var(ogfile, yca_id, clat)

       ierr = pio_put_var(ogfile, maskb_id, grid_imask)
       deallocate(grid_imask)

       ierr = pio_put_var(ogfile, areaA_id, areaA)
       ierr = pio_put_var(ogfile, areaB_id, areaB)
       deallocate(areaB)

       allocate(grid_imask(ncol))
       grid_imask=1

       ierr = pio_put_var(ogfile, maska_id, grid_imask)

       call pio_closefile(ogfile)

       deallocate(grid_imask, lat,lon, h1d, col, row, dg_dims, ldof)
       do ii=1,nelemd
          if(associated(interpdata(ii)%interp_xy))then
             deallocate(interpdata(ii)%interp_xy)
          endif
          if(associated(interpdata(ii)%ilat))then
             deallocate(interpdata(ii)%ilat)
          endif
          if (associated(interpdata(ii)%ilon))then
             deallocate(interpdata(ii)%ilon)
          endif
       end do


    end do

    call set_interp_parameter('itype',olditype)
    call set_interp_parameter('nlon',oldnlon)
    call set_interp_parameter('nlat',oldnlat)


  end subroutine create_native_mapping_files





end module native_mapping
