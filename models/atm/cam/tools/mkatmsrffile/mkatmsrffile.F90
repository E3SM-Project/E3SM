!===============================================================================
! SVN $Id:  $
! SVN $URL: $
! 12/06/2010  Jim Edwards jedwards@ucar.edu
! Interpolate files needed for cam atmosphere dry deposition to model grid
!===============================================================================

program mkatmsrffile
  use mpi
  use pio
  use mct_mod
  use shr_mct_mod
  use shr_kind_mod, only : r8=>shr_kind_r8, shr_kind_cl
  implicit none
  integer :: ierr, npes, iam, npft
  integer, pointer :: sfcgindex(:), atmgindex(:)

  type rptr
     real(r8), pointer :: fld(:)
  end type rptr

  type(rptr), pointer :: soilw(:), pft(:), apft(:), asoilw(:)


  type(iosystem_desc_t) :: iosystem

  type(mct_gsmap) :: gsMap_srf, gsMap_atm
  type(mct_SMatP) :: sMatP
  type(mct_aVect), target :: srf_av, atm_av



  integer, pointer :: comps(:) ! array with component ids
  integer, pointer :: comms(:) ! array with mpicoms
  type(io_desc_t) :: atm_iodesc, srf_iodesc

  character(len=shr_kind_cl) :: srffilename
  character(len=shr_kind_cl) :: atmfilename

  character(len=shr_kind_cl) :: soilwfilename
  character(len=shr_kind_cl) :: landfilename

  character(len=shr_kind_cl) :: outputfilename
  

  type(file_desc_t) :: landfile, newfile
  integer, pointer :: dof(:), dof2(:), dof3(:)
  real(r8), pointer :: landmask(:),lake(:), wetland(:), urban(:)
  real(r8), pointer :: alake(:), awetland(:), aurban(:), fraction_landuse(:,:)
  integer :: srfnx, atmnx, srfnxg, atmnxg, dimid, nlat, nlon, i, j, clen, index, dim1, dim2
  type(var_desc_t) :: vid, vid1, vid2
  
  character(len=*), parameter :: srffields(5) =(/"PCT_LAKE   ",&
                                                 "PCT_WETLAND",&
                                                 "PCT_URBAN  ", &
                                                 "SOILW      ", & 
                                                 "PCT_PFT    "/)
  
  
  character(len=220) :: rList
  character(len=6) :: str
  real(r8) :: total_land, fraction_soilw
  real(r8), pointer :: total_soilw(:,:)
  Character(len=*), parameter :: ConfigFileName="mkatmsrffile.rc"

  call mpi_init(ierr)

  call mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, iam, ierr)
  call pio_init(iam, MPI_COMM_WORLD, npes, 0, 1, pio_rearr_none, iosystem,base=0)
  allocate(comps(2), comms(2))
  comps(1)=1
  comps(2)=2
  call mpi_comm_dup(MPI_COMM_WORLD,comms(1),ierr)
  call mpi_comm_dup(MPI_COMM_WORLD,comms(2),ierr)
  call mct_world_init(2, MPI_COMM_WORLD, comms, comps)
  

  call I90_allLoadF(ConfigFileName,0,MPI_COMM_WORLD,ierr)

  call I90_label('srfFileName:', ierr)
  call i90_gtoken(srffilename, ierr)
  call I90_label('atmFileName:', ierr)
  call i90_gtoken(atmfilename, ierr)
  call I90_label('landFileName:', ierr)
  call i90_gtoken(landfilename, ierr)
  call I90_label('soilwFileName:', ierr)
  call i90_gtoken(soilwfilename, ierr)
  call I90_label('outputFileName:', ierr)
  call i90_gtoken(outputfilename, ierr)
  call i90_release(ierr)


  call openfile_and_initdecomp(iosystem, srffilename, npes, iam, gsmap_srf, srfnx, srfnxg)
  call openfile_and_initdecomp(iosystem, atmfilename, npes, iam, gsmap_atm, atmnx, atmnxg)

 
  call shr_mct_sMatPInitnc(sMatP,gsmap_srf, gsmap_atm, "mkatmsrffile.rc", &
       "srf2atmFmapname:","srf2atmFmaptype:", MPI_COMM_WORLD)






  ierr = pio_openfile(iosystem, landFile, pio_iotype_netcdf, landfilename, pio_noclobber)

  ierr = pio_inq_dimid(landFile, 'lon', dimid)
  ierr = pio_inq_dimlen(landfile, dimid, nlon)
  ierr = pio_inq_dimid(landFile, 'lat', dimid)
  ierr = pio_inq_dimlen(landfile, dimid, nlat)
  ierr = pio_inq_dimid(landFile, 'pft', dimid)
  ierr = pio_inq_dimlen(landfile, dimid, npft)

  call mct_gsmap_OrderedPoints(gsMap_srf, iam, Dof)

  call pio_initdecomp(iosystem, pio_double, (/nlon,nlat/), dof, srf_iodesc)

  deallocate(dof)

  rlist = ' '
  clen=1
  do i=1,npft+15
     if(i<=npft) then
        write(str,'(A,i2.2,A)') 'pft',i,':'
        rlist(clen:clen+5) = str
        clen=clen+6
     else if(i<=npft+12) then
        write(str,'(A,i2.2,A)') 'slw',i-npft,':'
        rlist(clen:clen+5) = str
        clen=clen+6
     else
        if(i==npft+13) clen=clen-1
        rlist(clen:clen+len_trim(srffields(i-12-npft))) = ':'//trim(srffields(i-12-npft))
        clen = clen+len_trim(srffields(i))+1
     end if
  end do


  call mct_aVect_init(srf_av, rlist=trim(rlist), lsize=srfnx)
  call mct_aVect_zero(srf_av)
  call mct_aVect_init(atm_av, rlist=rlist, lsize=atmnx)
  call mct_aVect_zero(atm_av)


  index = mct_avect_indexra(srf_av,'PCT_LAKE')
  lake => srf_av%rattr(index,:)
  ierr = pio_inq_varid( landFile, 'PCT_LAKE', vid) 
  call pio_read_darray(landFile, vid, srf_iodesc, lake, ierr)
  lake = lake * 0.01_r8

  ierr = pio_inq_varid( landFile, 'PCT_PFT', vid) 
  allocate(pft(npft),apft(npft))
  do i=1,npft
     write(str,'(A,i2.2)') 'pft',i

     pft(i)%fld => srf_av%rattr(mct_avect_indexra(srf_av,str(1:5)),:)
     apft(i)%fld => atm_av%rattr(mct_avect_indexra(atm_av,str(1:5)),:)

     call pio_setframe(vid,int(i,kind=PIO_OFFSET))
     call pio_read_darray(landFile, vid, srf_iodesc, pft(i)%fld, ierr)
     pft(i)%fld = pft(i)%fld * 0.01_r8
  end do


  index = mct_avect_indexra(srf_av,'PCT_WETLAND')
  wetland => srf_av%rattr(index,:)
  ierr = pio_inq_varid( landFile, 'PCT_WETLAND', vid) 
  call pio_read_darray(landFile, vid, srf_iodesc, wetland, ierr)
  wetland = wetland * 0.01_r8

  index = mct_avect_indexra(srf_av,'PCT_URBAN')
  urban => srf_av%rattr(index,:)
  ierr = pio_inq_varid( landFile, 'PCT_URBAN', vid) 
  call pio_read_darray(landFile, vid, srf_iodesc, urban, ierr)
  urban = urban * 0.01_r8

  allocate(landmask(srfnx))
  ierr = pio_inq_varid( landFile, 'LANDMASK', vid) 
  call pio_read_darray(landFile, vid, srf_iodesc, landmask, ierr)

  call pio_closefile(landfile)

!  call pio_freedecomp(iosystem, srf_iodesc)

  ierr = pio_openfile(iosystem, landFile, pio_iotype_netcdf, soilwfilename, pio_noclobber)

!  call pio_initdecomp(iosystem, pio_double, (/nlon,nlat,12/), vdof, srf_iodesc)
  ierr = pio_inq_varid( landFile, 'SOILW', vid) 
  allocate(soilw(12),asoilw(12))
  do i=1,12
     str = ' '
     write(str,'(A,i2.2)') 'slw',i     
     soilw(i)%fld => srf_av%rattr(mct_avect_indexra(srf_av,str(1:5)),:)
     asoilw(i)%fld => atm_av%rattr(mct_avect_indexra(atm_av,str(1:5)),:)

     call pio_setframe(vid,int(i,kind=PIO_OFFSET))
     call pio_read_darray(landFile, vid, srf_iodesc, soilw(i)%fld, ierr)
  end do
  call pio_closefile(landfile)
  call pio_freedecomp(iosystem, srf_iodesc)

  
  do i=1,srfnx
     if(nint(landmask(i)) == 0) then
        lake(i) = 1.0
        wetland(i) = 0.0
        urban(i) = 0.0
        do j=1,12
           soilw(j)%fld(i) = 0.0
        end do
     end if
  end do
  deallocate(landmask)

  index = mct_avect_indexra(atm_av,'PCT_LAKE')
  alake => atm_av%rattr(index,:)


  call mct_sMat_avMult( srf_av, smatP, atm_av)

  index = mct_avect_indexra(atm_av,'PCT_LAKE')
  alake => atm_av%rattr(index,:)


  index = mct_avect_indexra(atm_av,'PCT_WETLAND')
  awetland => atm_av%rattr(index,:)

  index = mct_avect_indexra(atm_av,'PCT_URBAN')
  aurban => atm_av%rattr(index,:)



  fraction_soilw=0.0

  allocate(fraction_landuse(atmnx,11))
  allocate(total_soilw(atmnx,12))

  fraction_landuse = 0.0_r8
  do i=1,atmnx
     total_soilw(i,:)=0.0
     total_land = (alake(i)+awetland(i)+aurban(i))
     do j=1,npft
        total_land=total_land+apft(j)%fld(i)
     end do
     fraction_soilw = total_land - (alake(i)+wetland(i))
     if(total_land < 1.0_r8) then
        alake(i) = alake(i) + (1.0_r8 - total_land)
     end if
     
!     print *,i,fraction, fraction_soilw
!     if(abs(fraction-1.0_r8) > 0.1_r8) then
!        print *, i, fraction, alake(i), awetland(i), aurban(i), (apft(j)%fld(i), j=1,npft)
!     end if

     do j=1,12
        total_soilw(i,j) = total_soilw(i,j) + asoilw(j)%fld(i) * fraction_soilw
     end do

     fraction_landuse(i,1) = aurban(i)
     fraction_landuse(i,2) = apft(16)%fld(i) + apft(17)%fld(i)
     fraction_landuse(i,3) = apft(13)%fld(i) + apft(14)%fld(i) + apft(15)%fld(i)
     fraction_landuse(i,4) = apft(5)%fld(i) + apft(6)%fld(i) + apft(7)%fld(i)+ apft(8)%fld(i) + apft(9)%fld(i)
     fraction_landuse(i,5) = apft(2)%fld(i) + apft(3)%fld(i) + apft(4)%fld(i)
     fraction_landuse(i,6) = awetland(i)
     fraction_landuse(i,7) = alake(i)
     fraction_landuse(i,8) = apft(1)%fld(i)
     fraction_landuse(i,11) = apft(10)%fld(i) + apft(11)%fld(i) + apft(12)%fld(i)

     if(abs(sum(fraction_landuse(i,:)-1._r8)) > 0.001_r8) then
        fraction_landuse(i,:) = fraction_landuse(i,:)/sum(fraction_landuse(i,:))
     end if
     

  end do

  ierr = pio_createfile(iosystem, newFile, pio_iotype_netcdf, trim(outputfilename), pio_clobber)

  ierr = pio_def_dim(newFile, 'ncol', atmnxg, dim1)
  ierr = pio_def_dim(newFile, 'class',11, dim2)

  ierr = pio_def_var(newFile, 'fraction_landuse', pio_double, (/dim1,dim2/), vid1) 
  ierr = pio_def_dim(newFile, 'month',12, dim2)

  ierr = pio_def_var(newFile, 'soilw', pio_double, (/dim1,dim2/), vid2) 

  ierr = pio_enddef(newFile)
  
  call mct_gsmap_OrderedPoints(gsMap_atm, iam, Dof)



  allocate(dof2(atmnx*12))
  do j=1,12
     do i=1,atmnx
        dof2(i+(j-1)*atmnx) = dof(i)+(j-1)*atmnxg
     end do
  end do

  call pio_initdecomp(iosystem, pio_double, (/atmnxg,11/), dof2(1:11*atmnx-1), atm_iodesc)

  call pio_write_darray(newFile, vid1, atm_iodesc, fraction_landuse ,ierr)
  call pio_freedecomp(newfile, atm_iodesc)

  call pio_initdecomp(iosystem, pio_double, (/atmnxg,12/), dof2, atm_iodesc)

  call pio_write_darray(newFile, vid2, atm_iodesc, total_soilw ,ierr)
  call pio_freedecomp(newfile, atm_iodesc)


  call pio_closefile(newFile)
  call pio_finalize(iosystem, ierr)

  deallocate(comps, comms)
  call mpi_finalize(ierr)


contains

  subroutine openfile_and_initdecomp(iosystem, filename, npes, iam, gsmap, nx, nxg)
    type(iosystem_desc_t) :: iosystem
    character(len=*), intent(in) :: filename
    integer, intent(in) :: npes, iam
    type(mct_gsmap), intent(out) :: gsmap
    integer, intent(out) :: nx, nxg
    integer, pointer :: gindex(:)



    gindex => get_grid_index(iosystem, filename, npes, iam, nx, nxg)
    call mct_gsMap_init( gsMap, gindex, MPI_COMM_WORLD,1 , nx, nxg)
    deallocate(gindex)



  end subroutine openfile_and_initdecomp


  function get_grid_index(iosystem, filename, npes, iam, nx, nxg) result(gindex)
    use pio
    implicit none
    character(len=*), intent(in) :: filename
    type(iosystem_desc_t) :: iosystem
    integer, intent(in) :: npes, iam
    integer, intent(out) :: nx, nxg

    type(file_desc_t) :: file
    integer, pointer :: gindex(:)


    integer :: dimid, ierr, add1=0, i, start_offset

    ierr = pio_openfile(iosystem, File, PIO_IOTYPE_NETCDF, filename, PIO_NOCLOBBER)

    ierr = pio_inq_dimid(File, 'grid_size', dimid)
    ierr = pio_inq_dimlen(File, dimid, nxg)

    nx = nxg/npes

    if(nx*npes < nxg-(npes-iam-1)) then
       start_offset = nxg-(npes-iam-1)-(nx*npes)-1
       add1 = 1
    else
       add1 = 0
       start_offset=0
    end if
    allocate(gindex(nx+add1))
    do i=1,nx+add1
       gindex(i)=i+iam*nx+start_offset
    end do



  call pio_closefile(FILE)

  end function get_grid_index

end program mkatmsrffile




