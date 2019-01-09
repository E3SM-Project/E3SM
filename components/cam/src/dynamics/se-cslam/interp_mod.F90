module interp_mod
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use dimensions_mod,      only: nelemd, np, ne
  use interpolate_mod,     only: interpdata_t
  use interpolate_mod,     only: interp_lat => lat, interp_lon => lon
  use interpolate_mod,     only: interp_gweight => gweight
  use dyn_grid,            only: elem,fvm
  use spmd_utils,          only: masterproc, iam
  use cam_history_support, only: fillvalue
  use hybrid_mod,          only: hybrid_t, config_thread_region
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
!JMD  This is strange ithr used before it is set
      ithr = 0
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
      hybrid = config_thread_region(par,'serial')
!      ithr = omp_get_thread_num()
!      hybrid = hybrid_create(par,ithr,1)

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
              call cam_grid_attribute_register(trim(gridname), 'w',           &
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
    use dimensions_mod,   only: npsq, fv_nphys,nc,nhc,nhc_phys
    use dof_mod,          only: PutUniquePoints
    use interpolate_mod,  only: get_interp_parameter
    use shr_pio_mod,      only: shr_pio_getiosys
    use edge_mod,         only: edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use edgetype_mod,     only: EdgeBuffer_t
    use bndry_mod,        only: bndry_exchange
    use parallel_mod,     only: par
    use thread_mod,       only: horz_num_threads
    use cam_grid_support, only: cam_grid_id
    use hybrid_mod,       only: hybrid_t,config_thread_region, get_loop_ranges
    use fvm_mapping,      only: fvm2dyn,phys2dyn
    use fvm_mod,          only: fill_halo_and_extend_panel

    type(file_desc_t), intent(inout) :: File
    type(var_desc_t) , intent(inout) :: varid
    real(r8),          intent(in)    :: fld(:,:,:)
    integer,           intent(in)    :: numlev, data_type, decomp_type
    !
    ! local variables
    !
    type(io_desc_t)                :: iodesc
    type(hybrid_t)                 :: hybrid
    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t)            :: edgebuf              ! edge buffer


    integer              :: lchnk, i, j, icol, ncols, pgcols(pcols), ierr
    integer              :: idmb1(1), idmb2(1), idmb3(1), nets, nete
    integer, allocatable :: bpter(:,:)! offsets into block buffer for packing data
    integer              :: cpter(pcols,0:pver)    ! offsets into chunk buffer for unpacking data
    integer              :: phys_decomp, fvm_decomp,gll_decomp

    real(r8), pointer     :: dest(:,:,:,:)
    real(r8), pointer     :: bbuffer(:), cbuffer(:), fldout(:,:)
    real(r8), allocatable :: fld_dyn(:,:,:), fld_tmp(:,:,:,:,:)

    integer          :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer          :: nlon, nlat, ncol,nsize,nhalo,nhcc
    logical          :: usefillvalues

    usefillvalues=.false.

    phys_decomp = cam_grid_id('physgrid')
    fvm_decomp  = cam_grid_id('FVM')
    gll_decomp  = cam_grid_id('GLL')
    !
    ! There are 2 main scenarios regarding decomposition:
    !
    !   decomp_type==phys_decomp: we need to move data from physics decomposition to dynamics decomposition
    !   else                    : data is on dynamics decomposition
    !
    if (decomp_type==phys_decomp) then
      if (fv_nphys>0) then
        !
        ! note that even if fv_nphys<4 then SIZE(fld,DIM=1)=PCOLS
        !
        nsize = fv_nphys
        nhalo = 1!for bilinear only a halo of 1 is needed
        nhcc  = nhc_phys
      else
        nsize = np
        nhalo = 0!no halo needed (lat-lon point always surrounded by GLL points)
        nhcc  = 0
      end if
    else if (decomp_type==fvm_decomp) then
      !
      ! CSLAM grid output
      !
      nsize = nc
      nhalo = 1!for bilinear only a halo of 1 is needed
      nhcc  = nhc
    else if (decomp_type==gll_decomp) then
      nsize = np
      nhalo = 0!no halo needed (lat-lon point always surrounded by GLL points)
      nhcc  = 0
    else
      call endrun('write_interpolated_scalar: unknown decomp_type')
    end if
    allocate(fld_dyn(nsize*nsize,numlev,nelemd))
    allocate(fld_tmp(1-nhcc:nsize+nhcc,1-nhcc:nsize+nhcc,numlev,1,nelemd))
    allocate(dest(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo,numlev,nelemd))

    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys(atm_id)

    if(decomp_type==phys_decomp) then
      fld_dyn = -999_R8
      if(local_dp_map) then
        !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, ioff,k)
        do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)
          call get_gcol_all_p(lchnk,pcols,pgcols)
          if (fv_nphys>0) ncols = fv_nphys*fv_nphys
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
        allocate( bbuffer(block_buf_nrecs*numlev) )!xxx Steve: this is different that dp_coupling? (no numlev in dp_coupling)
        allocate( cbuffer(chunk_buf_nrecs*numlev) )

        !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, i, k, icol)
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
          if (fv_nphys>0) then
             allocate(bpter(fv_nphys*fv_nphys,0:pver))
          else
             allocate(bpter(npsq,0:pver))
          end if
          !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, k, ncols, icol)
          do ie=1,nelemd
            if (fv_nphys>0) then
              call chunk_to_block_recv_pters(elem(ie)%GlobalID,fv_nphys*fv_nphys,pverp,1,bpter)
              ncols = fv_nphys*fv_nphys
            else
              call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pverp,1,bpter)
              ncols = elem(ie)%idxp%NumUniquePts
            end if
            do k = 1, numlev
              do icol=1,ncols
                fld_dyn(icol,k,ie) = bbuffer(bpter(icol,k-1))
              end do
            end do
          end do
        end if
        deallocate( bbuffer )
        deallocate( cbuffer )
        deallocate( bpter   )

      end if!local_dp_map
      if (fv_nphys>0) then
        do ie = 1, nelemd
          fld_tmp(1:nsize,1:nsize,:,1,ie) = RESHAPE(fld_dyn(:,:,ie),(/nsize,nsize,numlev/))
        end do
      else
        call initEdgeBuffer(par, edgebuf, elem, numlev,nthreads=1)

        do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, numlev, fld_dyn(1:ncols,1:numlev,ie), fld_tmp(:,:,1:numlev,1,ie))
          call edgeVpack(edgebuf, fld_tmp(:,:,1:numlev,1,ie), numlev, 0, ie)
        end do
        if(iam < par%nprocs) then
          call bndry_exchange(par, edgebuf,location='write_interpolated_scalar')
        end if
        do ie=1,nelemd
          call edgeVunpack(edgebuf, fld_tmp(:,:,1:numlev,1,ie), numlev, 0, ie)
        end do
        call freeEdgeBuffer(edgebuf)
        usefillvalues = any(fld_tmp == fillvalue)
      end if
    else
      !
      ! not physics decomposition
      !
      do ie = 1, nelemd
        fld_tmp(1:nsize,1:nsize,1:numlev,1,ie) = RESHAPE(fld(1:nsize*nsize,1:numlev,ie),(/nsize,nsize,numlev/))
      end do
    end if
    deallocate(fld_dyn)
    !
    ! code for non-GLL grids: need to fill halo and interpolate (if on panel edge/corner) for bilinear interpolation
    !
    if (decomp_type==fvm_decomp.or.(fv_nphys>0.and.decomp_type==phys_decomp)) then
      !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
      !JMD        hybrid = config_thread_region(par,'horizontal')
      hybrid = config_thread_region(par,'serial')
      call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
      call fill_halo_and_extend_panel(elem(nets:nete),fvm(nets:nete),&
           fld_tmp(:,:,:,:,nets:nete),hybrid,nets,nete,nsize,nhcc,nhalo,numlev,1,.true.,.true.)
    end if
    !
    ! WARNING - 1:nelemd and nets:nete
    !
    !$OMP MASTER   !JMD
    dest(:,:,:,1:nelemd) = fld_tmp(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo,:,1,1:nelemd)
    !$OMP END MASTER
    deallocate(fld_tmp)
    !
    !***************************************************************************
    !
    ! now data is on dynamics decomposition
    !
    !***************************************************************************
    !
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
      if(usefillvalues) then
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), nsize, nhalo, numlev, fldout(st:en,:), fillvalue)
      else
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,:,ie), nsize, nhalo, numlev, fldout(st:en,:))
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
    use cam_instance,     only: atm_id
    use interpolate_mod,  only: interpolate_scalar, vec_latlon_to_contra,get_interp_parameter
    use spmd_dyn,         only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
    use ppgrid,           only: begchunk, endchunk, pcols, pverp, pver
    use phys_grid,        only: get_gcol_all_p, get_ncols_p, chunk_to_block_send_pters, chunk_to_block_recv_pters, &
         transpose_chunk_to_block
    use dyn_grid,         only: get_gcol_block_d
    use hybrid_mod,       only: hybrid_t,config_thread_region, get_loop_ranges
    use dimensions_mod,   only: npsq, fv_nphys,nc,nhc,nhc_phys
    use dof_mod,          only: PutUniquePoints
    use shr_pio_mod,      only: shr_pio_getiosys
    use edge_mod,         only: edgevpack, edgevunpack, initedgebuffer, freeedgebuffer
    use edgetype_mod,     only: EdgeBuffer_t
    use bndry_mod,        only: bndry_exchange
    use parallel_mod,     only: par
    use thread_mod,       only: horz_num_threads
    use cam_grid_support, only: cam_grid_id
    use fvm_mod,          only: fill_halo_and_extend_panel
    use control_mod,      only: cubed_sphere_map
    use cube_mod,         only: dmap

    implicit none
    type(file_desc_t), intent(inout) :: File
    type(var_desc_t),  intent(inout) :: varidu, varidv
    real(r8),          intent(in)    :: fldu(:,:,:), fldv(:,:,:)
    integer,           intent(in)    :: numlev, data_type, decomp_type

    type(hybrid_t)                 :: hybrid
    type(io_desc_t)                :: iodesc
    type(iosystem_desc_t), pointer :: pio_subsystem
    type (EdgeBuffer_t)            :: edgebuf              ! edge buffer

    integer              :: lchnk, i, j, icol, ncols, pgcols(pcols), ierr, nets, nete
    integer              :: idmb1(1), idmb2(1), idmb3(1)
    integer, allocatable :: bpter(:,:)    ! offsets into block buffer for packing data
    integer              :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

    real(r8), allocatable :: dest(:,:,:,:,:)
    real(r8), pointer     :: bbuffer(:), cbuffer(:), fldout(:,:,:)
    real(r8), allocatable :: fld_dyn(:,:,:,:),fld_tmp(:,:,:,:,:)

    integer          :: st, en, ie, ioff, ncnt_out, k
    integer, pointer :: idof(:)
    integer          :: nlon, nlat, ncol,nsize,nhalo,nhcc
    logical          :: usefillvalues
    integer          :: phys_decomp, fvm_decomp,gll_decomp
    real (r8)        :: D(2,2)   ! derivative of gnomonic mapping
    real (r8)        :: v1,v2

    usefillvalues=.false.

    phys_decomp = cam_grid_id('physgrid')
    fvm_decomp  = cam_grid_id('FVM')
    gll_decomp  = cam_grid_id('GLL')
    !
    ! There are 2 main scenarios regarding decomposition:
    !
    !   decomp_type==phys_decomp: we need to move data from physics decomposition to dynamics decomposition
    !   else                    : data is on dynamics decomposition
    !
    if (decomp_type==phys_decomp) then
      if (fv_nphys>0) then
        !
        ! note that even if fv_nphys<4 then SIZE(fld,DIM=1)=npsq
        !
        nsize = fv_nphys
        nhalo = 1!for bilinear only a halo of 1 is needed
        nhcc  = nhc_phys
      else
        nsize = np
        nhalo = 0!no halo needed (lat-lon point always surrounded by GLL points)
        nhcc  = 0
      end if
    else if (decomp_type==fvm_decomp) then
      !
      ! CSLAM grid output
      !
      nsize = nc
      nhalo = 1!for bilinear only a halo of 1 is needed
      nhcc  = nhc
    else if (decomp_type==gll_decomp) then
      nsize = np
      nhalo = 0!no halo needed (lat-lon point always surrounded by GLL points)
      nhcc  = 0
    else
      call endrun('write_interpolated_scalar: unknown decomp_type')
    end if
    allocate(fld_dyn(nsize*nsize,2,numlev,nelemd))
    allocate(fld_tmp(1-nhcc:nsize+nhcc,1-nhcc:nsize+nhcc,2,numlev,nelemd))
    allocate(dest(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo,2,numlev,nelemd))

    nlon=get_interp_parameter('nlon')
    nlat=get_interp_parameter('nlat')
    pio_subsystem => shr_pio_getiosys(atm_id)
    fld_dyn = -999_R8
    if(decomp_type==phys_decomp) then
      if(local_dp_map) then
        !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, pgcols, icol, idmb1, idmb2, idmb3, ie, k, ioff)
        do lchnk=begchunk,endchunk
          ncols=get_ncols_p(lchnk)
          call get_gcol_all_p(lchnk,pcols,pgcols)
          if (fv_nphys>0) ncols = fv_nphys*fv_nphys
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
        !$omp parallel do num_threads(horz_num_threads) private (lchnk, ncols, cpter, i, k, icol)
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
          if (fv_nphys>0) then
            allocate(bpter(fv_nphys*fv_nphys,0:pver))
          else
            allocate(bpter(npsq,0:pver))
          end if
          !$omp parallel do num_threads(horz_num_threads) private (ie, bpter, k, icol)
          do ie=1,nelemd
            if (fv_nphys>0) then
              call chunk_to_block_recv_pters(elem(ie)%GlobalID,fv_nphys*fv_nphys,pverp,2,bpter)
              ncols = fv_nphys*fv_nphys
            else
              call chunk_to_block_recv_pters(elem(ie)%GlobalID,npsq,pverp,2,bpter)
              ncols = elem(ie)%idxp%NumUniquePts
            end if
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
        deallocate( bpter   )
      end if!local_dp_map
      if (fv_nphys>0) then
        do ie = 1, nelemd
          fld_tmp(1:nsize,1:nsize,:,:,ie) = RESHAPE(fld_dyn(:,:,:,ie),(/nsize,nsize,2,numlev/))
        end do
      else
        call initEdgeBuffer(par, edgebuf, elem, 2*numlev,nthreads=1)

        do ie=1,nelemd
          ncols = elem(ie)%idxp%NumUniquePts
          call putUniquePoints(elem(ie)%idxP, 2, numlev, fld_dyn(1:ncols,:,1:numlev,ie), fld_tmp(:,:,:,1:numlev,ie))
          call edgeVpack(edgebuf, fld_tmp(:,:,:,:,ie), 2*numlev, 0, ie)
        enddo
        if(iam < par%nprocs) then
          call bndry_exchange(par, edgebuf,location='write_interpolated_vector')
        end if

        do ie=1,nelemd
          call edgeVunpack(edgebuf, fld_tmp(:,:,:,:,ie), 2*numlev, 0, ie)
        enddo
        call freeEdgeBuffer(edgebuf)
        usefillvalues = any(fld_tmp==fillvalue)
      end if
    else
      !
      ! not physics decomposition
      !
      usefillvalues = (any(fldu(1:nsize:1,nsize,:)==fillvalue) .or. any(fldv(1:nsize:1,nsize,:)==fillvalue))
      do ie = 1, nelemd
        fld_tmp(1:nsize,1:nsize,1,1:numlev,ie) = RESHAPE(fldu(1:nsize*nsize,1:numlev,ie),(/nsize,nsize,numlev/))
        fld_tmp(1:nsize,1:nsize,2,1:numlev,ie) = RESHAPE(fldv(1:nsize*nsize,1:numlev,ie),(/nsize,nsize,numlev/))
      end do
    endif
    deallocate(fld_dyn)
    !
    !***************************************************************************
    !
    ! now data is on dynamics decomposition
    !
    !***************************************************************************
    !
    if (decomp_type==fvm_decomp.or.(fv_nphys>0.and.decomp_type==phys_decomp)) then
      !
      !***************************************************************************
      !
      ! code for non-GLL grids: need to fill halo and interpolate
      ! (if on panel edge/corner) for bilinear interpolation
      !
      !***************************************************************************
      !

      !JMD $OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete,n)
      !JMD        hybrid = config_thread_region(par,'horizontal')
      hybrid = config_thread_region(par,'serial')
      call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
      call fill_halo_and_extend_panel(elem(nets:nete),fvm(nets:nete),&
           fld_tmp(:,:,:,:,nets:nete),hybrid,nets,nete,nsize,nhcc,nhalo,2,numlev,.true.,.false.)
      do ie=nets,nete
        call vec_latlon_to_contra(elem(ie),nsize,nhcc,numlev,fld_tmp(:,:,:,:,ie),fvm(ie))
      end do
      call fill_halo_and_extend_panel(elem(nets:nete),fvm(nets:nete),&
           fld_tmp(:,:,:,:,nets:nete),hybrid,nets,nete,nsize,nhcc,nhalo,2,numlev,.false.,.true.)
    else
      do ie=1,nelemd
        call vec_latlon_to_contra(elem(ie),nsize,nhcc,numlev,fld_tmp(:,:,:,:,ie))
      end do
    end if
    !
    ! WARNING - 1:nelemd and nets:nete
    !
    !$OMP MASTER   !JMD
    dest(:,:,:,:,1:nelemd) = fld_tmp(1-nhalo:nsize+nhalo,1-nhalo:nsize+nhalo,:,:,1:nelemd)
    !$OMP END MASTER
    deallocate(fld_tmp)
    !
    !***************************************************************************
    !
    ! do mapping from source grid to latlon grid
    !
    !***************************************************************************
    !
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
      if(usefillvalues) then
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,1,:,ie), nsize, nhalo, numlev, fldout(st:en,:,1), fillvalue)
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,2,:,ie), nsize, nhalo, numlev, fldout(st:en,:,2), fillvalue)
      else
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,1,:,ie), nsize, nhalo, numlev, fldout(st:en,:,1))
        call interpolate_scalar(cam_interpolate(ie),dest(:,:,2,:,ie), nsize, nhalo, numlev, fldout(st:en,:,2))
      end if
      !
      ! convert from contravariant components to lat-lon
      !
      do i=1,cam_interpolate(ie)%n_interp
        ! convert fld from contra->latlon
        call dmap(D,cam_interpolate(ie)%interp_xy(i)%x,cam_interpolate(ie)%interp_xy(i)%y,&
             elem(ie)%corners3D,cubed_sphere_map,elem(ie)%corners,elem(ie)%u2qmap,elem(ie)%facenum)
        ! convert fld from contra->latlon
        do k=1,numlev
          v1 = fldout(st+i-1,k,1)
          v2 = fldout(st+i-1,k,2)
          fldout(st+i-1,k,1)=D(1,1)*v1 + D(1,2)*v2
          fldout(st+i-1,k,2)=D(2,1)*v1 + D(2,2)*v2
        end do
      end do
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
