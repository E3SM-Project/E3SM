module native_mapping
!
!  Create mapping files using the SE basis functions.  This module looks for the namelist 'native_mapping' in 
!  file NLFileName (usually atm_in) and reads from it a list of up to maxoutgrids grid description files
!  It then creates a grid mapping file from the currently defined SE grid to the grid described in each file
!  using the SE basis functions.   The output mapping file name is generated based on the SE model resolution
!  and the input grid file name and ends in '_date_native.nc'
! 
  use kinds,            only: iulog, r8=>real_kind
  use parallel_mod,     only: haltmp, iam, mpi_character, mpi_logical, mpi_integer, mpi_max
  use control_mod,      only: cubed_sphere_map,dd_pi
  use common_io_mod,    only: infilenames, PIOFS,io_stride, num_io_procs, num_agg, output_dir

  implicit none
  private
  public :: create_native_mapping_files

  integer, parameter :: maxoutgrids=5
  character(len=240) :: native_mapping_outgrids(maxoutgrids)
  
contains

  subroutine create_native_mapping_files(hybrid, elem, maptype)
    use parallel_mod, only : parallel_t, global_shared_buf, global_shared_sum
    use global_norms_mod, only: wrap_repro_sum
    use element_mod, only : element_t
    use hybrid_mod, only : hybrid_t, hybrid_create
    use pio, only : pio_noerr, pio_openfile, pio_createfile, pio_closefile, &
         pio_get_var, pio_put_var, pio_write_darray,pio_int, pio_double, &
         pio_def_var, pio_put_att, pio_global, file_desc_t, var_desc_t, &
         io_desc_t, pio_internal_error,pio_inq_dimlen, pio_inq_varid, &
         pio_get_att, pio_enddef, pio_bcast_error,pio_internal_error, &
         pio_def_dim, pio_inq_dimid, pio_seterrorhandling, pio_initdecomp, &
         iotype_netcdf, pio_clobber, pio_setdebuglevel, pio_rearr_box, pio_init
    use quadrature_mod, only : quadrature_t, gauss, gausslobatto
    use interpolate_mod, only : interpdata_t, find_ref_coordinates, interpolate_scalar, set_interp_parameter, &
         interp_init
    use coordinate_systems_mod, only : spherical_polar_t, cartesian2d_t
    use dimensions_mod, only : nelemd, ne, np, npsq, nelem
    use reduction_mod, only : ParallelMin,ParallelMax
    use cube_mod, only : convert_gbl_index
    !use dyn_grid, only : get_horiz_grid_d
    use dof_mod, only : CreateMetaData
    use shr_const_mod, only : shr_const_pi
    use thread_mod,     only: omp_get_thread_num

    character(len=*), intent(in) :: maptype
    type(element_t), intent(in) :: elem(:)
    type(hybrid_t) :: hybrid

    character(len=240) :: mappingfile, fname	
    logical :: exist
    type (spherical_polar_t) :: sphere
    type(file_desc_t) :: ogfile
    type (interpdata_t)  :: interpdata(nelemd)
    integer :: ierr, dimid, npts, vid, ncol
    real(r8), allocatable :: lat(:), lon(:), clat(:), clon(:)
    integer :: i, ii, gid, ie2, je2, ie, je, face_no, face_no2, k, j, n, ngrid, tpts, nf
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
    character(len=12) :: unit_str,nestr,npstr
    real(r8), allocatable :: areaA(:), areaB(:)
    integer :: cntperelem_in(nelem),  cntperelem_out(nelem)
    integer :: ithr, dg_rank, substr1, substr2

    type(interpdata_t), pointer :: mapping_interpolate(:)
    integer, allocatable :: gid_list(:),gid_in(:)
    integer :: itype
    real(r8) :: fill_double=1.e36_r8
    real(r8) :: rad2deg 
    rad2deg = 180/dd_pi


    if (maptype=='native') then
       itype=0
    else if (maptype=='bilin') then
       itype=1
    else
       call haltmp('bad interp_type')
    endif

    if(.not.hybrid%par%dynproc) then
       ! The special case of npes_se < npes_cam is not worth dealing with here
       call haltmp('Native mapping code requires npes_se==npes_cam')
    end if
    
    call interp_init()

    call PIO_Init(hybrid%par%rank, hybrid%par%comm, num_io_procs, num_agg, &
         io_stride, pio_rearr_box, PIOFS)



    do nf=1,1  ! maxoutgrids
       !fname = native_mapping_outgrids(nf) 
       fname = infilenames(1)
       if(hybrid%masterthread) write(iulog,*) ""
       if(hybrid%masterthread) write(iulog,*) "map type =",trim(maptype)
       if(hybrid%masterthread) write(iulog,*) "detination grid=",trim(fname)
       if(len_trim(fname)==0) cycle
       inquire(file=fname,exist=exist)
       if(.not. exist) then
          write(iulog,*) 'WARNING: Could not find or open grid file ',trim(fname)
          cycle
       end if
       ierr = pio_openfile( PIOFS, ogfile, iotype_netcdf, trim(fname))

       ierr = pio_inq_dimid( ogfile, 'grid_size', dimid)
       ierr = pio_inq_dimlen( ogfile, dimid, npts)
       allocate(lat(npts), lon(npts), grid_imask(npts), areab(npts))
       allocate(gid_list(npts))

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
          areaB=fill_double
       end if

       if(unit_str .eq. 'degrees') then
          lat = lat * shr_const_pi/180_r8
          lon = lon * shr_const_pi/180_r8
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
       if (hybrid%masterthread) print *,'locating all interpolation points...'
       sphere%r=1    
       do i=1,npts
          gid_list(i)=-1
          if(grid_imask(i)==1) then
             sphere%lat=lat(i)
             sphere%lon=lon(i)
             call find_ref_coordinates(sphere,elem,ii,gid,cart)
             if (ii /= -1) then 
                interpdata(ii)%n_interp = interpdata(ii)%n_interp + 1
                gid_list(i)=gid
             endif

             if(hybrid%masterthread) then
                if(mod(i,npts/10).eq.1) then
                   print *,'finished point ',i,' of ',npts
                endif
             end if
          end if
       enddo

       ! remove duplicates.  a point claimed by two elements will be owned by the
       ! element with the largest gid:
       allocate(gid_in(npts))
       gid_in=gid_list
       call mpi_allreduce(gid_in,gid_list,npts, MPI_INTEGER, MPI_MAX, hybrid%par%comm, ierr)
       deallocate(gid_in)


       countx=maxval(interpdata(1:nelemd)%n_interp)
       count_max = ParallelMax(countx,hybrid)

       ! allocate storage, including duplicates
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
             call find_ref_coordinates(sphere,elem,ii,gid,cart)
             if (ii /= -1 .and. gid==gid_list(i)) then 
                ngrid = interpdata(ii)%n_interp + 1
                interpdata(ii)%n_interp = ngrid
                interpdata(ii)%interp_xy( ngrid ) = cart
                interpdata(ii)%ilon( ngrid ) = i
                interpdata(ii)%ilat( ngrid ) = i
             endif
          end if
       end do


       ! check if every point in interpolation grid was claimed by an element:
       countx=sum(interpdata(1:nelemd)%n_interp)
       global_shared_buf(1,1) = countx
       call wrap_repro_sum(nvars=1, comm=hybrid%par%comm, nsize=1)
       count_total = global_shared_sum(1)
       tpts = sum(grid_imask)
       if (hybrid%masterthread) then
          write(iulog,'(a,f8.2)') 'Average number of interpolation points per element: ',count_total/real(nelem)
          write(iulog,'(a,f8.0)') 'Maximum number of interpolation points on any element: ',count_max
       endif
       if (count_total /= tpts ) then
          write(iulog,*)__FILE__,__LINE__,iam, count_total, tpts, npts
          call haltmp('Error setting up interpolation grid count_total<>npts') 
       endif



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
                call interpolate_scalar(interpdata(ie), f, np, h(:))

                do n=1,interpdata(ie)%n_interp
                   if(any(h/=h)) then
                      call haltmp('nan generated')
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


       call mpi_allreduce(cntperelem_in, cntperelem_out, nelem, MPI_INTEGER, MPI_MAX, hybrid%par%comm, ierr)


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
       ncol =  2+nelem*(np-1)**2

       allocate(areaA(ncol))
       allocate(clat(ncol),clon(ncol))
       !call get_horiz_grid_d(ncol, clat_d_out=clat, clon_d_out=clon, area_d_out=areaA)
       ! todo: compute these in HOMME
       clat=0
       clon=0
       areaA=0

       substr1 = index(fname,'/',BACK=.true.)
       substr2 = index(fname,'.nc',BACK=.true.)

       write(nestr,*) ne
       write(npstr,*) np
       mappingfile = trim(output_dir) // &
            'map_ne' // trim(adjustl(nestr)) // 'np' // trim(adjustl(npstr)) // &
            '_to_' // fname(substr1+1:substr2-1) // '_' // trim(maptype) // ".nc"

       ierr = pio_createfile(PIOFS, ogfile, iotype_netcdf, mappingfile,PIO_CLOBBER)

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
       ierr = pio_put_att( ogfile, areaB_id, '_FillValue',fill_double)

       ierr = pio_def_var( ogfile, 'mask_a', pio_int, (/na_dim/), maska_id)
       ierr = pio_def_var( ogfile, 'mask_b', pio_int, (/nb_dim/), maskb_id)



       ierr = pio_put_att( ogfile, xca_id, 'units','degrees')
       ierr = pio_put_att( ogfile, yca_id, 'units','degrees')
       ierr = pio_put_att( ogfile, xcb_id, 'units','degrees')
       ierr = pio_put_att( ogfile, ycb_id, 'units','degrees')

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

       lon=lon*rad2deg
       lat=lat*rad2deg
       clon=clon*rad2deg
       clat=clat*rad2deg

       ierr = pio_put_var(ogfile, xcb_id, lon)
       ierr = pio_put_var(ogfile, ycb_id, lat)   

       ierr = pio_put_var(ogfile, xca_id, clon)
       ierr = pio_put_var(ogfile, yca_id, clat)
       ierr = pio_put_var(ogfile, maskb_id, grid_imask)
       deallocate(grid_imask)

       ierr = pio_put_var(ogfile, areaA_id, areaA)
       ierr = pio_put_var(ogfile, areaB_id, areaB)
       deallocate(areaA,areaB)

       allocate(grid_imask(ncol))
       grid_imask=1

       ierr = pio_put_var(ogfile, maska_id, grid_imask)

       call pio_closefile(ogfile)

       deallocate(grid_imask, lat,lon, clat, clon, h1d, col, row, dg_dims, ldof)
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



  end subroutine create_native_mapping_files





end module native_mapping


