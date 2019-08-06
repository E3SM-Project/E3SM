!===============================================================================
! 12/06/2010  Jim Edwards jedwards@ucar.edu
! 06/24/2018  Balwinder Singh (balwinder.singh -at- pnnl (dot) gov) removed pio, 
!             mct, mpi and csm_share dependencies
!
! Interpolate files needed for cam atmosphere dry deposition to model grid
!===============================================================================

program mkatmsrffile

  use netcdf, only: nf90_open, nf90_inq_dimid, nf90_inquire_dimension,      &
       nf90_inq_varid, nf90_get_var, nf90_close, NF90_NOWRITE, nf90_create, &
       nf90_def_dim, nf90_def_var, nf90_enddef, nf90_put_var, NF90_NOERR,   &
       nf90_strerror, NF90_DOUBLE, nf90_netcdf4

  implicit none

  !Define datatypes (from csm share library)
  integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real
  integer,parameter :: shr_kind_cl = 256           ! long char
  integer,parameter :: shr_kind_cx = 512           ! extra-long char

  !Other parameters
  real(r8),parameter :: huge_real = huge(1.0_r8)   ! used to initialize variables
  
  !Other variables and pointers
  type rptr
     real(r8), allocatable :: fld(:)
  end type rptr

  type(rptr), pointer :: soilw(:), pft(:), apft(:), asoilw(:)

  character(len=shr_kind_cl) :: srffilename, atmfilename
  character(len=shr_kind_cl) :: soilwfilename, landfilename, outputfilename

  character(len=shr_kind_cx) :: srf2atmFmapname

  integer :: i_lp, j_lp, k_lp                       !loop indices
  integer :: srfnx, atmnx, dimid, nlat, nlon, dim1, dim2, dim3,npft, ntime
  integer :: ncid_map, ncid_land, ncid_soil, ncid_out
  integer :: varid, varid1, varid2, n_a, n_b, num_elements
  integer :: total_grd_pts, irow, icol, nclass

  integer, pointer :: col(:),row(:)

  real(r8) :: total_land, fraction_soilw, rwgt

  real(r8), pointer :: landmask(:),lake(:), wetland(:), urban(:)
  real(r8), pointer :: alake(:), awetland(:), aurban(:), fraction_landuse(:,:)
  real(r8), pointer :: total_soilw(:,:)
  real(r8), pointer :: tmp2d(:,:), tmp3d(:,:,:)
  real(r8), pointer :: wgt(:)


  !--------------------------------------------------
  ! READ NAMELIST
  !--------------------------------------------------

  namelist /input/srfFileName,atmFileName,landFileName,soilwFileName,srf2atmFmapname,outputFileName 
  
  open(101, file = 'nml_atmsrf', status = 'old')
  read(101, nml=input)
  close(101)

  !--------------------------------------------------
  ! READ GRID SIZES
  !--------------------------------------------------

  !Read grid sizes of srf and atm files
  call openfile_and_initdecomp(srffilename, srfnx)
  call openfile_and_initdecomp(atmfilename, atmnx)

  !--------------------------------------------------
  ! READ MAP FILE
  !--------------------------------------------------

  !Read map file for weights and other parameters required for remapping
  call nc_check(nf90_open(trim(srf2atmFmapname), NF90_NOWRITE, ncid_map), __LINE__)
  call nc_check(nf90_inq_dimid(ncid_map,"n_a", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_map, dimid, len = n_a), __LINE__ )

  !Sanity check
  if(n_a .ne. srfnx) then
     print*,'ERROR: Map dimension (n_a) is not equal to srf grid size dimension (grid_size)'
     print*,'n_a=',n_a,' and srf grid_size=',srfnx
     print*,'Exiting'
     call exit(1)
  endif

  call nc_check(nf90_inq_dimid(ncid_map,"n_b", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_map, dimid, len = n_b), __LINE__ )

  !Sanity check
  if(n_b .ne. atmnx) then
     print*,'ERROR: Map dimension (n_b) is not equal to atm grid size dimension (grid_size)'
     print*,'n_b=',n_b,' and atm grid_size=',atmnx
     print*,'Exiting'
     call exit(1)
  endif

  call nc_check(nf90_inq_dimid(ncid_map,"n_s", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_map, dimid, len = num_elements), __LINE__ )

  !allocate and read map file variables
  allocate(col(num_elements))
  call nc_check(nf90_inq_varid(ncid_map,'col',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_map,varid,col), __LINE__)

  allocate(row(num_elements))
  call nc_check(nf90_inq_varid(ncid_map,'row',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_map,varid,row), __LINE__)

  allocate(wgt(num_elements))
  call nc_check(nf90_inq_varid(ncid_map,'S',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_map,varid,wgt), __LINE__)

  !Close map file
  call nc_check(nf90_close ( ncid_map ), __LINE__)

  !--------------------------------------------------
  ! READ LAND INPUT FILE
  !--------------------------------------------------

  !Read Land file
  call nc_check(nf90_open(landfilename, NF90_NOWRITE, ncid_land), __LINE__)
  
  call nc_check(nf90_inq_dimid(ncid_land, "lon", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_land, dimid, len = nlon), __LINE__ )
  call nc_check(nf90_inq_dimid(ncid_land, "lat", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_land, dimid, len = nlat), __LINE__ )

  !For reshaping arrays to 1d compute total # of grid points
  total_grd_pts = nlon*nlat

  !Sanity check
  if(total_grd_pts .ne. srfnx) then
     print*,'ERROR: Land file nlat*nlon is not equal to srf grid size dimension (grid_size)'
     print*,'nlon*nlat=',total_grd_pts,' and srf grid_size=',srfnx
     print*,'Exiting'
     call exit(1)
  endif

  call nc_check(nf90_inq_dimid(ncid_land, "pft", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_land, dimid, len = npft), __LINE__ )

  !Allocate temporary variables to read data from the netcdf file
  allocate(tmp2d(nlon,nlat),tmp3d(nlon,nlat,npft))
  tmp2d(:,:)   = huge_real
  tmp3d(:,:,:) = huge_real

  allocate(lake(total_grd_pts))
  call nc_check(nf90_inq_varid(ncid_land,'PCT_LAKE',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_land,varid,tmp2d), __LINE__)
  lake(:) = reshape(tmp2d,(/total_grd_pts/)) !reshape to 1d array
  lake = lake * 0.01_r8


  call nc_check(nf90_inq_varid(ncid_land,'PCT_PFT',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_land,varid,tmp3d), __LINE__)

  allocate(pft(npft),apft(npft)) ! apft is allocated here for atm

  !Storing in a data structure
  do i_lp = 1, npft
     !Allocate land data structure and store read in values
     allocate(pft(i_lp)%fld(total_grd_pts))     
     pft(i_lp)%fld = reshape(tmp3d(:,:,i_lp),(/total_grd_pts/)) * 0.01_r8

     !Allocate and initialize corresponding atm data structure
     allocate(apft(i_lp)%fld(atmnx))
     apft(i_lp)%fld(:) = 0.0_r8
  end do

  tmp2d(:,:)   = huge_real  !Reinitialize tmp2d to inf
  allocate(wetland(total_grd_pts))
  call nc_check(nf90_inq_varid(ncid_land,'PCT_WETLAND',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_land,varid,tmp2d), __LINE__)
  wetland = reshape(tmp2d,(/total_grd_pts/))
  wetland = wetland * 0.01_r8

  tmp2d(:,:)   = huge_real  !Reinitialize tmp2d to inf
  allocate(urban(total_grd_pts))
  call nc_check(nf90_inq_varid(ncid_land,'PCT_URBAN',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_land,varid,tmp2d), __LINE__)

  urban = reshape(tmp2d,(/total_grd_pts/))
  urban = urban * 0.01_r8

  tmp2d(:,:)   = huge_real    !Reinitialize tmp2d to inf
  allocate(landmask(srfnx))
  call nc_check(nf90_inq_varid(ncid_land,'LANDMASK',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_land,varid,tmp2d), __LINE__)
  landmask = reshape(tmp2d,(/total_grd_pts/))

  deallocate(tmp2d)

  !close land file
  call nc_check(nf90_close ( ncid_land ), __LINE__)

  !--------------------------------------------------
  ! READ SOIL FILE
  !--------------------------------------------------

  !open soil file
  call nc_check(nf90_open(soilwfilename, NF90_NOWRITE, ncid_soil), __LINE__)
  call nc_check(nf90_inq_dimid(ncid_soil, "time", dimid), __LINE__ )
  call nc_check(nf90_inquire_dimension(ncid_land, dimid, len = ntime), __LINE__ )

  deallocate(tmp3d) !deallocate as shape of tmp3d will change for next read

  allocate(tmp3d(nlon,nlat,ntime))
  call nc_check(nf90_inq_varid(ncid_soil,'SOILW',varid), __LINE__)
  call nc_check(nf90_get_var(ncid_soil,varid,tmp3d), __LINE__)

  !close soil file
  call nc_check(nf90_close ( ncid_soil ), __LINE__)

  allocate(soilw(ntime),asoilw(ntime)) ! allocate corresponding atm var asoilw as well
  do i_lp = 1, ntime
     !Allocate land data structure and store read in values
     allocate(soilw(i_lp)%fld(total_grd_pts))
     soilw(i_lp)%fld = reshape(tmp3d(:,:,i_lp),(/total_grd_pts/))

     !Allocate and initialize corresponding atm data structure
     allocate(asoilw(i_lp)%fld(atmnx))
     asoilw(i_lp)%fld(:) = 0.0_r8
  end do
  deallocate(tmp3d)

  do i_lp = 1, srfnx
     if(nint(landmask(i_lp)) == 0) then
        lake(i_lp) = 1.0
        wetland(i_lp) = 0.0
        urban(i_lp) = 0.0
        do j_lp=1,ntime
           soilw(j_lp)%fld(i_lp) = 0.0
        end do
     end if
  end do

  !Deallocate memory
  deallocate(landmask)

  !--------------------------------------------------
  ! Remap to atm grid
  !--------------------------------------------------

  allocate(alake(atmnx), awetland(atmnx), aurban(atmnx))
  alake(:)    = 0.0_r8
  awetland(:) = 0.0_r8
  aurban(:)   = 0.0_r8 
  do i_lp = 1, num_elements 
     irow = row(i_lp)
     icol = col(i_lp)
     rwgt = wgt(i_lp)

     !Following equation is obtained from : cime/src/externals/mct/mct/m_MatAttrVectMul.F90 (line 254, master hash:8f364f4b926)
     alake(irow) =  alake(irow) + rwgt*lake(icol)
     awetland(irow) =  awetland(irow) + rwgt*wetland(icol)
     aurban(irow) =  aurban(irow) + rwgt*urban(icol)

     do j_lp = 1, npft
        apft(j_lp)%fld(irow) = apft(j_lp)%fld(irow) + rwgt*pft(j_lp)%fld(icol)
     enddo

     do j_lp = 1, ntime 
        asoilw(j_lp)%fld(irow) = asoilw(j_lp)%fld(irow) + rwgt*soilw(j_lp)%fld(icol)
     enddo
  enddo

  !Deallocate memory
  deallocate(col, row, wgt, lake, pft, wetland, urban, soilw)

  !--------------------------------------------------
  ! Compute fields to output
  !--------------------------------------------------
  fraction_soilw = 0.0_r8

  nclass = 11
  allocate(fraction_landuse(atmnx,nclass)) 
  allocate(total_soilw(atmnx,ntime))

  fraction_landuse = 0.0_r8
  do i_lp=1,atmnx
     total_soilw(i_lp,:)=0.0
     total_land = (alake(i_lp)+awetland(i_lp)+aurban(i_lp))
     do j_lp=1,npft
        total_land=total_land+apft(j_lp)%fld(i_lp)
     end do
     fraction_soilw = total_land - (alake(i_lp)+awetland(i_lp))
     if(total_land < 1.0_r8) then
        alake(i_lp) = alake(i_lp) + (1.0_r8 - total_land)
     end if
     

     do j_lp=1,ntime
        total_soilw(i_lp,j_lp) = total_soilw(i_lp,j_lp) + asoilw(j_lp)%fld(i_lp) * fraction_soilw
     end do

     fraction_landuse(i_lp,1 ) = aurban(i_lp)
     fraction_landuse(i_lp,2 ) = apft(16)%fld(i_lp) + apft(17)%fld(i_lp)
     fraction_landuse(i_lp,3 ) = apft(13)%fld(i_lp) + apft(14)%fld(i_lp) + apft(15)%fld(i_lp)
     fraction_landuse(i_lp,4 ) = apft(5 )%fld(i_lp) + apft(6 )%fld(i_lp) + apft(7)%fld(i_lp )+ apft(8)%fld(i_lp) + apft(9)%fld(i_lp)
     fraction_landuse(i_lp,5 ) = apft(2 )%fld(i_lp) + apft(3 )%fld(i_lp) + apft(4)%fld(i_lp )
     fraction_landuse(i_lp,6 ) = awetland(i_lp)
     fraction_landuse(i_lp,7 ) = alake(i_lp)
     fraction_landuse(i_lp,8 ) = apft(1 )%fld(i_lp)
     fraction_landuse(i_lp,11) = apft(10)%fld(i_lp) + apft(11)%fld(i_lp) + apft(12)%fld(i_lp)

     if(abs(sum(fraction_landuse(i_lp,:)-1._r8)) > 0.001_r8) then
        fraction_landuse(i_lp,:) = fraction_landuse(i_lp,:)/sum(fraction_landuse(i_lp,:))
     end if     
  end do

  !Deallocate memory
  deallocate(apft, awetland, aurban, asoilw, alake)

  !--------------------------------------------------
  ! Create output file and add fields
  !--------------------------------------------------

  call nc_check(nf90_create(trim(outputfilename), nf90_netcdf4, ncid = ncid_out), __LINE__ )   

  call nc_check(nf90_def_dim(ncid_out, 'ncol', atmnx, dim1), __LINE__ )
  call nc_check(nf90_def_dim(ncid_out, 'class', nclass, dim2), __LINE__ )
  call nc_check(nf90_def_dim(ncid_out, 'month', ntime, dim3), __LINE__)
  call nc_check(nf90_def_var(ncid_out, 'soilw', NF90_DOUBLE, (/dim1,dim3/), varid2), __LINE__)
  call nc_check(nf90_def_var(ncid_out, 'fraction_landuse', NF90_DOUBLE, (/dim1,dim2/), varid1),__LINE__ )
  
  call nc_check(nf90_enddef(ncid_out), __LINE__ )

  call nc_check(nf90_put_var(ncid_out, varid1, fraction_landuse), __LINE__)
  call nc_check(nf90_put_var(ncid_out, varid2, total_soilw ), __LINE__)

  call nc_check(nf90_close(ncid_out), __LINE__ )

  !Deallocate memory
  deallocate(fraction_landuse,total_soilw)

contains

  subroutine openfile_and_initdecomp(filename, nx)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: nx

    call nc_check(nf90_open(filename, NF90_NOWRITE, ncid_out), __LINE__)

    call nc_check(nf90_inq_dimid(ncid_out, "grid_size", dimid), __LINE__ )
    call nc_check(nf90_inquire_dimension(ncid_out, dimid, len = nx), __LINE__ )
    call nc_check(nf90_close ( ncid_out ), __LINE__)
  end subroutine openfile_and_initdecomp

  subroutine nc_check(status, line_num)
    integer, intent(in) :: status
    integer, intent(in) :: line_num !line number on which error occured
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print*,'Error at line ',line_num, ' of file ', __FILE__
      print*,'Exiting'
      call exit(1)
    end if
  end subroutine nc_check  

end program mkatmsrffile

!Q: Why total_soilw output field has huge values?
!A: 12th month of soilw var in the input file clim_soilw.nc has huge values, which
!   are then remapped to asoilw and then to total_soilw var
