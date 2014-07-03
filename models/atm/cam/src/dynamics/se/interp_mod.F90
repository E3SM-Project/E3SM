module interp_mod
  use cam_logfile, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  use dimensions_mod, only : nelemd, np
  use interpolate_mod, only : interpolate_scalar, setup_latlon_interp, set_interp_parameter, get_interp_lat, get_interp_lon, &
       var_is_vector_uvar, var_is_vector_vvar, interpolate_vector, interpdata_t, get_interp_gweight
  use dyn_grid,       only : elem, w
  use spmd_utils,       only : masterproc, iam
  use cam_pio_utils,  only: phys_decomp, fillvalue
  use hybrid_mod,     only : hybrid_t, hybrid_create
  use abortutils, only: endrun

  implicit none
  private
  type(interpdata_t), pointer :: cam_interpolate(:)

  public get_interp_lat, get_interp_lon, setup_history_interpolation, write_interpolated
  public var_is_vector_uvar, var_is_vector_vvar, latlon_interpolation, add_interp_attributes

  interface write_interpolated
     module procedure write_interpolated_scalar
     module procedure write_interpolated_vector
  end interface
  type(hybrid_t) :: hybrid

contains

  subroutine add_interp_attributes(file)
    use pio, only : file_desc_t, pio_put_att, pio_global
    use interpolate_mod, only : get_interp_parameter
    type(file_desc_t) :: file

    integer :: ierr
    integer :: itmp
    
    itmp = get_interp_parameter('itype')
    if(itmp == 0) then
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', 'se basis functions')
    else if(itmp == 1) then
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', 'bilinear')
    else
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_type', itmp)
    end if
    
    itmp = get_interp_parameter('gridtype')
    select case(itmp)
    case(1)
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'equally spaced with poles')
    case(2)        
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'Gauss')
    case(3)
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', 'equally spaced no poles')
    case default
       ierr = pio_put_att(file, PIO_GLOBAL, 'interp_outputgridtype', itmp)
    end select


  end subroutine add_interp_attributes

  subroutine setup_history_interpolation(mtapes)

    use dyn_comp, only : dom_mt
    use parallel_mod,   only: par
    use thread_mod,     only: omp_get_thread_num
    use interpolate_mod, only : interpolate_analysis, get_interp_parameter
    implicit none
  
    integer, intent(in) :: mtapes
    integer :: ithr, nthreads

    if(iam>= par%nprocs) return

    ithr=omp_get_thread_num()
    hybrid = hybrid_create(par,ithr,1)
       
    if(any(interpolate_analysis(1:mtapes))) then
       allocate(cam_interpolate(nelemd))
       call setup_latlon_interp(elem, cam_interpolate, par)
       allocate(w(get_interp_parameter('nlat')))
       w = get_interp_gweight()
    end if

  end subroutine setup_history_interpolation

  function latlon_interpolation(t)
    use interpolate_mod, only : interpolate_analysis
    integer, intent(in) :: t

    logical :: latlon_interpolation

    latlon_interpolation = interpolate_analysis(t)
  end function latlon_interpolation



  subroutine write_interpolated_scalar(File, varid, fld, numlev, data_type, decomp_type) 
    use pio, only : file_desc_t, io_desc_t, var_desc_t, pio_write_darray, iosystem_desc_t, &
         pio_initdecomp, pio_freedecomp, pio_setdebuglevel
    use cam_instance, only: atm_id
    use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
    use ppgrid, only : begchunk, endchunk, pcols, pver
    use phys_grid, only : get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
         transpose_chunk_to_block
    use dyn_grid,       only: get_gcol_block_d
    use dimensions_mod, only: npsq
    use element_mod, only : element_t
    use dof_mod, only : PutUniquePoints
    use interpolate_mod, only : get_interp_parameter
    use shr_pio_mod, only : shr_pio_getiosys
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use bndry_mod, only : bndry_exchangeV
    use parallel_mod,   only: par
    use abortutils, only : endrun
    
    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varid
    real(r8), intent(in) :: fld(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type

    type(io_desc_t) :: iodesc

    integer :: lchnk, i, j, m, icol, ncols, pgcols(pcols), ierr
    integer :: idmb1(1), idmb2(1), idmb3(1)
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8), pointer :: dest(:,:,:,:) 
    real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:)
    real(r8) :: fld_dyn(npsq,numlev,nelemd)
    integer :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer :: nlon, nlat, ncol
    logical :: usefillvalues=.false.
    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t) :: edgebuf              ! edge buffer



    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys(atm_id)

    if(decomp_type==phys_decomp) then
       fld_dyn = -999_R8
       if(local_dp_map) then
          !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, m)
          do lchnk=begchunk,endchunk
             ncols=get_ncols_p(lchnk)
             call get_gcol_all_p(lchnk,pcols,pgcols)
             do icol=1,ncols
                call get_gcol_block_d(pgcols(icol),1,idmb1,idmb2,idmb3)
                ie = idmb3(1)
                ioff=idmb2(1)
                do k=1,numlev
                   fld_dyn(ioff,k,ie)      = fld(icol, k, lchnk-begchunk+1)
                end do
             end do
             
          end do
       else

          allocate( bbuffer(block_buf_nrecs*numlev) )
          allocate( cbuffer(chunk_buf_nrecs*numlev) )

          !$omp parallel do private (lchnk, ncols, cpter, i, icol)
          do lchnk = begchunk,endchunk
             ncols = get_ncols_p(lchnk)
             
             call chunk_to_block_send_pters(lchnk,pcols,pver+1,1,cpter)
             
             do i=1,ncols
                cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
             end do

             do icol=1,ncols
                
                cbuffer   (cpter(icol,:))     = fld(icol,:,lchnk-begchunk+1)
             end do

          end do

          call transpose_chunk_to_block(1, cbuffer, bbuffer)
          if(iam < par%nprocs) then
!$omp parallel do private (ie, bpter, icol)
             do ie=1,nelemd
          
                call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,1,bpter)
                ncols = elem(ie)%idxp%NumUniquePts
                do icol=1,ncols
                   fld_dyn   (icol,:,ie)   = bbuffer(bpter(icol,:))
                end do

             end do
          end if
          deallocate( bbuffer )
          deallocate( cbuffer )

       end if
       allocate(dest(np,np,numlev,nelemd))
       call initEdgeBuffer(edgebuf, numlev)

       do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, numlev, fld_dyn(1:ncols,:,ie), dest(:,:,:,ie))
          call edgeVpack(edgebuf, dest(:,:,:,ie), numlev, 0, elem(ie)%desc)
       enddo
       if(iam < par%nprocs) then
          call bndry_exchangeV(par, edgebuf)
       end if
       do ie=1,nelemd
          call edgeVunpack(edgebuf, dest(:,:,:,ie), numlev, 0, elem(ie)%desc)
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
    use pio, only : file_desc_t, io_desc_t, var_desc_t, pio_write_darray, iosystem_desc_t, &
         pio_initdecomp, pio_freedecomp, pio_setdebuglevel
    use spmd_dyn,       only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
    use ppgrid, only : begchunk, endchunk, pcols, pver
    use phys_grid, only : get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
         transpose_chunk_to_block
    use dyn_grid,       only: get_gcol_block_d
    use dimensions_mod, only: npsq
    use element_mod, only : element_t
    use dof_mod, only : PutUniquePoints
    use interpolate_mod, only : get_interp_parameter
    use shr_pio_mod, only : shr_pio_getiosys
    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use bndry_mod, only : bndry_exchangeV
    use parallel_mod,   only: par
    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t), intent(inout) :: varidu, varidv
    real(r8), intent(in) :: fldu(:,:,:), fldv(:,:,:)
    integer, intent(in) :: numlev, data_type, decomp_type

    type(io_desc_t) :: iodesc

    integer :: lchnk, i, j, m, icol, ncols, pgcols(pcols), ierr
    integer :: idmb1(1), idmb2(1), idmb3(1)
    integer :: bpter(npsq,0:pver)    ! offsets into block buffer for packing data
    integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8), allocatable :: dest(:,:,:,:,:)
    real(r8), pointer :: bbuffer(:), cbuffer(:), fldout(:,:,:)
    real(r8) :: fld_dyn(npsq,2,numlev,nelemd)
    integer :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer :: nlon, nlat, ncol
    logical :: usefillvalues=.false.

    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t) :: edgebuf              ! edge buffer



    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys('ATM')
    fld_dyn = -999_R8
    if(decomp_type==phys_decomp) then
       allocate(dest(np,np,2,numlev,nelemd))
       if(local_dp_map) then
          !$omp parallel do private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff, m)
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
             
             call chunk_to_block_send_pters(lchnk,pcols,pver+1,2,cpter)
             
             do i=1,ncols
                cbuffer(cpter(i,1):cpter(i,1)) = 0.0_r8
             end do

             do icol=1,ncols
                do k=1,numlev
                   cbuffer   (cpter(icol,k))     = fldu(icol,k,lchnk-begchunk+1)
                   cbuffer   (cpter(icol,k)+1)   = fldv(icol,k,lchnk-begchunk+1)
                end do
             end do

          end do

          call transpose_chunk_to_block(2, cbuffer, bbuffer)
          if(iam < par%nprocs) then
             !$omp parallel do private (ie, bpter, icol)
             do ie=1,nelemd
                
                call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pver+1,2,bpter)
                ncols = elem(ie)%idxp%NumUniquePts
                do icol=1,ncols
                   do k=1,numlev
                      fld_dyn   (icol,1,k,ie)   = bbuffer(bpter(icol,k))
                      fld_dyn   (icol,2,k,ie)   = bbuffer(bpter(icol,k)+1)
                   enddo
                end do
                
             end do
          end if
          deallocate( bbuffer )
          deallocate( cbuffer )

       end if
       call initEdgeBuffer(edgebuf, 2*numlev)

       do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, 2, numlev, fld_dyn(1:ncols,:,:,ie), dest(:,:,:,:,ie))
          
          call edgeVpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, elem(ie)%desc)
       enddo
       if(iam < par%nprocs) then
          call bndry_exchangeV(par, edgebuf)
       end if

       do ie=1,nelemd
          call edgeVunpack(edgebuf, dest(:,:,:,:,ie), 2*numlev, 0, elem(ie)%desc)
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
                  dest(:,:,:,:,ie), np, numlev, fldout(st:en,:,:), 0, fillvalue) 
          else
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,ie), np, numlev, fldout(st:en,:,:), 0) 
          endif
       else
          do k=1,numlev
             do j=1,np
                do i=1,np
                   dest(i,j,1,k,1) = fldu(i+(j-1)*np,k,ie)
                   dest(i,j,1,k,1) = fldv(i+(j-1)*np,k,ie)
                end do
             end do
          end do
          if(usefillvalues) then
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,1), np, numlev, fldout(st:en,:,:), 0, fillvalue) 
          else
             call interpolate_vector(cam_interpolate(ie),elem(ie),&
                  dest(:,:,:,:,1), np, numlev, fldout(st:en,:,:), 0) 
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

