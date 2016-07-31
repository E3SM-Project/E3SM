module interp_mod
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use dimensions_mod,      only: nelemd, np, ne
  use interpolate_mod,     only: interpdata_t
  use interpolate_mod,     only: interp_lat => lat, interp_lon => lon
  use interpolate_mod,     only: interp_gweight => gweight
  use dyn_grid,            only: elem
  use spmd_utils,          only: masterproc, iam
  use cam_history_support, only: fillvalue
  use hybrid_mod,          only: hybrid_t, hybrid_create
  use cam_abortutils,      only: endrun

  implicit none
  private
  save

  public :: setup_history_interpolation
  public :: set_interp_hfile
  public :: write_interpolated

  interface write_interpolated
     module procedure write_interpolated_scalar
     module procedure write_interpolated_vector
  end interface

  ! hybrid is created in setup_history_interpolation
  type(hybrid_t) :: hybrid

! A type to hold interpdata info for each interpolated history file
type cam_interpolate_t
  type(interpdata_t), pointer :: interpdata(:) => NULL()
end type cam_interpolate_t

type(cam_interpolate_t), pointer :: interpdata_set(:) => NULL()  ! all files
type(interpdata_t),      pointer :: cam_interpolate(:) => NULL() ! curr. file

CONTAINS

  subroutine setup_history_interpolation(interp_ok, mtapes, interp_output,    &
       interp_info)

    use cam_history_support, only: interp_info_t
    use cam_history_support, only: interp_type_native
    use cam_history_support, only: interp_type_bilinear
    use cam_history_support, only: interp_gridtype_equal_poles
    use cam_history_support, only: interp_gridtype_gauss
    use cam_history_support, only: interp_gridtype_equal_nopoles
    use cam_grid_support,    only: horiz_coord_t, horiz_coord_create, iMap
    use cam_grid_support,    only: cam_grid_register, cam_grid_attribute_register
    use cam_grid_support,    only: max_hcoordname_len
    use interpolate_mod,     only: get_interp_lat, get_interp_lon
    use interpolate_mod,     only: get_interp_parameter, set_interp_parameter
    use interpolate_mod,     only: get_interp_gweight, setup_latlon_interp
    use parallel_mod,        only: par
    use thread_mod,          only: omp_get_thread_num
  
    ! Dummy arguments
    logical,             intent(inout) :: interp_ok
    integer,             intent(in)    :: mtapes
    logical,             intent(in)    :: interp_output(:)
    type(interp_info_t), intent(inout) :: interp_info(:)

    ! Local variables
    integer                            :: ithr, i, j
    real(r8),            pointer       :: w(:)
    integer(iMap),       pointer       :: grid_map(:,:)
    type(horiz_coord_t), pointer       :: lat_coord
    type(horiz_coord_t), pointer       :: lon_coord
    character(len=max_hcoordname_len)  :: gridname

    if (associated(cam_interpolate)) then
      do i = 1, size(cam_interpolate)
        if (associated(interpdata_set(ithr)%interpdata)) then
          deallocate(interpdata_set(ithr)%interpdata)
          nullify(interpdata_set(ithr)%interpdata)
        end if
      end do
      deallocate(cam_interpolate)
      nullify(cam_interpolate)
    end if
    nullify(grid_map)

    ! For this dycore, interpolated output should be OK
    interp_ok = (iam < par%nprocs)

    if (interp_ok) then
      ithr = omp_get_thread_num()
      hybrid = hybrid_create(par,ithr,1)
       
      if(any(interp_output)) then
        allocate(interpdata_set(mtapes))
        do i = 1, mtapes
          if (interp_output(i)) then
            if ( (interp_info(i)%interp_nlon == 0) .or.                       &
                 (interp_info(i)%interp_nlat == 0)) then
              ! compute interpolation grid based on number of points around equator
              call set_interp_parameter('auto', (4 * ne * (np-1)))
              interp_info(i)%interp_nlat = get_interp_parameter('nlat')
              interp_info(i)%interp_nlon = get_interp_parameter('nlon')
            else
              call set_interp_parameter('nlat', interp_info(i)%interp_nlat)
              call set_interp_parameter('nlon', interp_info(i)%interp_nlon)
            end if
            call set_interp_parameter('itype', interp_info(i)%interp_type)
            call set_interp_parameter('gridtype', interp_info(i)%interp_gridtype)

            allocate(interpdata_set(i)%interpdata(nelemd))
            ! Reset pointers in the interpolate module so they are not
            ! overwritten
            nullify(interp_lat)
            nullify(interp_lon)
            nullify(interp_gweight)
            call setup_latlon_interp(elem, interpdata_set(i)%interpdata, par)
            ! Create the grid coordinates
            lat_coord => horiz_coord_create('lat', '',                        &
                 interp_info(i)%interp_nlat, 'latitude', 'degrees_north',     &
                 1, interp_info(i)%interp_nlat, get_interp_lat())
            lon_coord => horiz_coord_create('lon', '',                        &
                 interp_info(i)%interp_nlon, 'longitude', 'degrees_east',     &
                 1, interp_info(i)%interp_nlon, get_interp_lon())
            ! Create a grid for this history file
            write(gridname, '(a,i0)') 'interp_out_', i
            interp_info(i)%grid_id = 200 + i
            call cam_grid_register(trim(gridname), interp_info(i)%grid_id,    &
                 lat_coord, lon_coord, grid_map, unstruct=.false.)
            interp_info(i)%gridname = trim(gridname)
            ! Add grid attributes
            allocate(w(get_interp_parameter('nlat')))
            w = get_interp_gweight()
            select case(interp_info(i)%interp_gridtype)
            case(interp_gridtype_equal_poles)
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_outputgridtype', 'equally spaced with poles')
              call cam_grid_attribute_register(trim(gridname), 'w',          &
                   'latitude weights', 'lat', w)
            case(interp_gridtype_equal_nopoles)
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_outputgridtype', 'equally spaced no poles')
              call cam_grid_attribute_register(trim(gridname), 'gw',          &
                   'latitude weights', 'lat', w)
            case(interp_gridtype_gauss)
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_outputgridtype', 'Gauss')
              call cam_grid_attribute_register(trim(gridname), 'gw',          &
                   'gauss weights', 'lat', w)
            case default
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_outputgridtype',                                   &
                   'Unknown interpolation output grid type',                  &
                   interp_info(i)%interp_gridtype)
            end select
            nullify(w) ! belongs to attribute
            if(interp_info(i)%interp_type == interp_type_native) then
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_type', 'se basis functions')
            else if(interp_info(i)%interp_type == interp_type_bilinear) then
              call cam_grid_attribute_register(trim(gridname),                &
                   'interp_type', 'bilinear')
            else
              call cam_grid_attribute_register(trim(gridname), 'interp_type', &
                   'Unknown interpolation type', interp_info(i)%interp_type)
            end if
            ! Store the data pointers for reuse later
            interp_info(i)%interp_lat => interp_lat
            interp_info(i)%interp_lon => interp_lon
            interp_info(i)%interp_gweight => interp_gweight
          end if
        end do
      end if
    end if

  end subroutine setup_history_interpolation

  subroutine set_interp_hfile(hfilenum, interp_info)
    use cam_history_support, only: interp_info_t
    use interpolate_mod,     only: set_interp_parameter

    ! Dummy arguments
    integer,             intent(in)    :: hfilenum
    type(interp_info_t), intent(inout) :: interp_info(:)

    if (.not. associated(interpdata_set)) then
      call endrun('SET_INTERP_HFILE: interpdata_set not allocated')
    else if ((hfilenum < 1) .or. (hfilenum > size(interpdata_set))) then
      call endrun('SET_INTERP_HFILE: hfilenum out of range')
    else if (hfilenum > size(interp_info)) then
      call endrun('SET_INTERP_HFILE: hfilenum out of range')
    else
     cam_interpolate => interpdata_set(hfilenum)%interpdata
      interp_lat => interp_info(hfilenum)%interp_lat
      interp_lon => interp_info(hfilenum)%interp_lon
      interp_gweight => interp_info(hfilenum)%interp_gweight
      call set_interp_parameter('nlat', interp_info(hfilenum)%interp_nlat)
      call set_interp_parameter('nlon', interp_info(hfilenum)%interp_nlon)
      call set_interp_parameter('itype', interp_info(hfilenum)%interp_type)
      call set_interp_parameter('gridtype', interp_info(hfilenum)%interp_gridtype)
    end if
  end subroutine set_interp_hfile

  subroutine write_interpolated_scalar(File, varid, fld, numlev, data_type, decomp_type)
    use pio,              only: file_desc_t, var_desc_t
    use pio,              only: iosystem_desc_t
    use pio,              only: pio_initdecomp, pio_freedecomp
    use pio,              only: io_desc_t, pio_write_darray
    use interpolate_mod,  only: interpolate_scalar
    use cam_instance,     only: atm_id
    use spmd_dyn,         only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
    use ppgrid,           only: begchunk, endchunk, pcols, pverp, pver
    use phys_grid,        only: get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
                               transpose_chunk_to_block
    use dyn_grid,         only: get_gcol_block_d
    use dimensions_mod,   only: npsq
    use dof_mod,          only: PutUniquePoints
    use interpolate_mod,  only: get_interp_parameter
    use shr_pio_mod,      only: shr_pio_getiosys
    use edgetype_mod,     only : edgebuffer_t
    use edge_mod,         only: edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use bndry_mod,        only: bndry_exchangeV
    use parallel_mod,     only: par
    use cam_grid_support, only: cam_grid_id

    
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varid
    real(r8), intent(in) :: fld(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type

    type(io_desc_t) :: iodesc

    integer :: lchnk, i, j, icol, ncols, pgcols(pcols), ierr
    integer :: idmb1(1), idmb2(1), idmb3(1)
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    integer :: phys_decomp

    real(r8), pointer :: dest(:,:,:,:) 
    real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:)
    real(r8) :: fld_dyn(npsq,numlev,nelemd)
    integer :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer :: nlon, nlat, ncol
    logical :: usefillvalues
    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t) :: edgebuf              ! edge buffer

    usefillvalues=.false.
    phys_decomp = cam_grid_id('physgrid')

    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys(atm_id)

    if(decomp_type==phys_decomp) then
       fld_dyn = -999_R8
       if(local_dp_map) then
          !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff)
          do lchnk=begchunk,endchunk
             ncols=get_ncols_p(lchnk)
             call get_gcol_all_p(lchnk,pcols,pgcols)
             do icol=1,ncols
                call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
                ie = idmb3(1)
                ioff=idmb2(1)
                do k=1,numlev
                   fld_dyn(ioff,k,ie) = fld(icol, k, lchnk-begchunk+1)
                end do
             end do

          end do
       else

          allocate( bbuffer(block_buf_nrecs*numlev) )
          allocate( cbuffer(chunk_buf_nrecs*numlev) )

          !$omp parallel do private (lchnk, ncols, cpter, i, icol)
          do lchnk = begchunk,endchunk
             ncols = get_ncols_p(lchnk)

             call chunk_to_block_send_pters(lchnk,pcols,pverp,1,cpter)

             do i=1,ncols
                cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
             end do

             do k=1,numlev
                do icol=1,ncols
                   cbuffer(cpter(icol,k-1)) = fld(icol,k,lchnk-begchunk+1)
                end do
             end do

          end do

          call transpose_chunk_to_block(1, cbuffer, bbuffer)
          if(iam < par%nprocs) then
             !$omp parallel do private (ie, bpter, icol)
             do ie=1,nelemd

                call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pverp,1,bpter)
                ncols = elem(ie)%idxp%NumUniquePts
                do k = 1, numlev
                  do icol=1,ncols
                    fld_dyn(icol,k,ie) = bbuffer(bpter(icol,k-1))
                  end do
                end do

             end do
          end if
          deallocate( bbuffer )
          deallocate( cbuffer )

       end if
       allocate(dest(np,np,numlev,nelemd))
       call initEdgeBuffer(hybrid%par,edgebuf, elem,numlev, numthreads_in=1)

       do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, numlev, fld_dyn(1:ncols,:,ie), dest(:,:,:,ie))
          call edgeVpack(edgebuf, dest(:,:,:,ie), numlev, 0, ie)
       enddo
       if(iam < par%nprocs) then
          call bndry_exchangeV(par, edgebuf)
       end if
       do ie=1,nelemd
          call edgeVunpack(edgebuf, dest(:,:,:,ie), numlev, 0, ie)
       end do
       call freeEdgeBuffer(edgebuf)
       usefillvalues = any(dest == fillvalue)
    else
      usefillvalues=any(fld==fillvalue)
      allocate(dest(np,np,numlev,1))
    end if

     ncnt_out = sum(cam_interpolate(1:nelemd)%n_interp)
     allocate(fldout(ncnt_out,numlev))
    allocate(idof(ncnt_out*numlev))
    fldout = -999_r8
    idof = 0
    st = 1
    
    

    do ie=1,nelemd
       ncol = cam_interpolate(ie)%n_interp
       do k=0,numlev-1
          do i=1,ncol
             idof(st+i-1+k*ncnt_out)=cam_interpolate(ie)%ilon(i)+nlon*(cam_interpolate(ie)%ilat(i)-1)+nlon*nlat*k
          enddo
       enddo
       
       ! Now that we have the field on the dyn grid we need to interpolate
       en = st+cam_interpolate(ie)%n_interp-1
       if(decomp_type==phys_decomp) then
          if(usefillvalues) then
             call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), np, numlev, fldout(st:en,:), fillvalue) 
          else
             call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), np, numlev, fldout(st:en,:)) 
          end if
       else
          do j=1,np
             do i=1,np
                dest(i,j,:,1) = fld(i+(j-1)*np,:,ie)
             end do
          end do
          if(usefillvalues) then
             call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,1), &
                  np, numlev, fldout(st:en,:), fillvalue) 
          else
             call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,1), &
                  np, numlev, fldout(st:en,:)) 
          end if
       end if


       st = en+1
    end do


    if(numlev==1) then
       call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat/), idof, iodesc)
    else
       call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat,numlev/), idof, iodesc)
    end if
    call pio_write_darray(File, varid, iodesc, fldout, ierr)

    deallocate(dest)

    deallocate(fldout)
    deallocate(idof)
    call pio_freedecomp(file,iodesc)

  end subroutine write_interpolated_scalar

  subroutine write_interpolated_vector(File, varidu, varidv, fldu, fldv, numlev, data_type, decomp_type)
    use pio,              only: file_desc_t, var_desc_t
    use pio,              only: iosystem_desc_t
    use pio,              only: pio_initdecomp, pio_freedecomp
    use pio,              only: io_desc_t, pio_write_darray
    use interpolate_mod,  only: interpolate_vector
    use spmd_dyn,         only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
    use ppgrid,           only : begchunk, endchunk, pcols, pverp, pver
    use phys_grid,        only : get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
         transpose_chunk_to_block

    use dyn_grid,         only: get_gcol_block_d
    use dimensions_mod,   only: npsq
    use dof_mod,          only : PutUniquePoints
    use interpolate_mod,  only : get_interp_parameter
    use shr_pio_mod,      only : shr_pio_getiosys
    use edgetype_mod,     only : edgebuffer_t
    use edge_mod,         only : edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use bndry_mod,        only : bndry_exchangeV
    use parallel_mod,     only: par
    use cam_grid_support, only: cam_grid_id
    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varidu, varidv
    real(r8), intent(in) :: fldu(:,:,:), fldv(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type

    type(io_desc_t) :: iodesc

    integer :: lchnk, i, j, icol, ncols, pgcols(pcols), ierr
    integer :: idmb1(1), idmb2(1), idmb3(1)
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data
    integer :: phys_decomp

    real(r8), allocatable :: dest(:,:,:,:,:)
    real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:,:)
    real(r8) :: fld_dyn(npsq,2,numlev,nelemd)
    integer :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer :: nlon, nlat, ncol
    logical :: usefillvalues

    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t) :: edgebuf              ! edge buffer

    usefillvalues=.false.
    phys_decomp = cam_grid_id('physgrid')

    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys('ATM')
    fld_dyn = -999_R8
    if(decomp_type==phys_decomp) then
       allocate(dest(np,np,2,numlev,nelemd))
       if(local_dp_map) then
          !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff)
          do lchnk=begchunk,endchunk
             ncols=get_ncols_p(lchnk)
             call get_gcol_all_p(lchnk,pcols,pgcols)
             
             do icol=1,ncols
                call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
                ie = idmb3(1)
                ioff=idmb2(1)
                do k=1,numlev
                   fld_dyn(ioff,1,k,ie)      = fldu(icol, k, lchnk-begchunk+1)
                   fld_dyn(ioff,2,k,ie)      = fldv(icol, k, lchnk-begchunk+1)
                end do
             end do
             
          end do
       else

          allocate( bbuffer(2*block_buf_nrecs*numlev) )
          allocate( cbuffer(2*chunk_buf_nrecs*numlev) )

          !$omp parallel do private (lchnk, ncols, cpter, i, icol)
          do lchnk = begchunk,endchunk
             ncols = get_ncols_p(lchnk)
             
             call chunk_to_block_send_pters(lchnk,pcols,pverp,2,cpter)
             
             do i=1,ncols
                cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
             end do

             do icol=1,ncols
                do k=1,numlev
                   cbuffer(cpter(icol,k-1))   = fldu(icol,k,lchnk-begchunk+1)
                   cbuffer(cpter(icol,k-1)+1) = fldv(icol,k,lchnk-begchunk+1)
                end do
             end do

          end do

          call transpose_chunk_to_block(2, cbuffer, bbuffer)
          if(iam < par%nprocs) then
             !$omp parallel do private (ie, bpter, icol)
             do ie=1,nelemd
                
                call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pverp,2,bpter)
                ncols = elem(ie)%idxp%NumUniquePts
                do icol=1,ncols
                   do k=1,numlev
                      fld_dyn(icol,1,k,ie) = bbuffer(bpter(icol,k-1))
                      fld_dyn(icol,2,k,ie) = bbuffer(bpter(icol,k-1)+1)
                   enddo
                end do
                
             end do
          end if
          deallocate( bbuffer )
          deallocate( cbuffer )

       end if
       call initEdgeBuffer(hybrid%par,edgebuf,elem, 2*numlev, numthreads_in=1)

       do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, 2, numlev, fld_dyn(1:ncols,:,1:numlev,ie), dest(:,:,:,1:numlev,ie))
          
          call edgeVpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, ie)
       enddo
       if(iam < par%nprocs) then
          call bndry_exchangeV(par, edgebuf)
       end if

       do ie=1,nelemd
          call edgeVunpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, ie)
       enddo
       call freeEdgeBuffer(edgebuf)
       usefillvalues = any(dest==fillvalue)
    else
       usefillvalues = (any(fldu==fillvalue) .or. any(fldv==fillvalue))
       allocate(dest(np,np,2,numlev,1))
    endif
    ncnt_out = sum(cam_interpolate(1:nelemd)%n_interp)
    allocate(fldout(ncnt_out,numlev,2))
    allocate(idof(ncnt_out*numlev))
    
    fldout = -999_r8
    idof = 0
    st = 1
    do ie=1,nelemd
       ncol = cam_interpolate(ie)%n_interp
       do k=0,numlev-1
          do i=1,ncol
             idof(st+i-1+k*ncnt_out)=cam_interpolate(ie)%ilon(i)+nlon*(cam_interpolate(ie)%ilat(i)-1)+nlon*nlat*k
          enddo
       enddo
       

    ! Now that we have the field on the dyn grid we need to interpolate
       en = st+cam_interpolate(ie)%n_interp-1
       if(decomp_type==phys_decomp) then
          if(usefillvalues) then
             call interpolate_vector(cam_interpolate(ie),elem(ie), &
                  dest(:,:,:,:,ie), numlev, fldout(st:en,:,:), 0, fillvalue) 
          else
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,ie), numlev, fldout(st:en,:,:), 0) 
          endif
       else
          do k=1,numlev
             do j=1,np
                do i=1,np
                   dest(i,j,1,k,1) = fldu(i+(j-1)*np,k,ie)
                   dest(i,j,2,k,1) = fldv(i+(j-1)*np,k,ie)
                end do
             end do
          end do
          if(usefillvalues) then
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,1), numlev, fldout(st:en,:,:), 0, fillvalue) 
          else
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,1), numlev, fldout(st:en,:,:), 0) 
          end if
       end if

       st = en+1
    end do

    if(numlev==1) then
       call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat/), idof, iodesc)
    else
       call pio_initdecomp(pio_subsystem, data_type, (/nlon,nlat,numlev/), idof, iodesc)
    end if

    call pio_write_darray(File, varidu, iodesc, fldout(:,:,1), ierr)

    call pio_write_darray(File, varidv, iodesc, fldout(:,:,2), ierr)


    deallocate(fldout)
    deallocate(idof)
    deallocate(dest)
    call pio_freedecomp(file,iodesc)

  end subroutine write_interpolated_vector

end module interp_mod

