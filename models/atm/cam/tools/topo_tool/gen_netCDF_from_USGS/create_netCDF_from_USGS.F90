!
!  DATE CODED:   Nov 4, 2011
!  DESCRIPTION:  This program reads USGS 30-sec terrain dataset in 33 tiles and converts
!                them one single NetCDF file.
!
!  Author: Peter Hjort Lauritzen (pel@ucar.edu) 
!                USGS read in code is from Jiundar Chern (jchern@dao.gsfc.nasa.gov)
!
!  ROUTINES CALLED:
!       netcdf routines
!
!  COMPILING:
!
      program convterr
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none
!
        integer,  parameter :: ntile = 33         ! number of tiles in USGS GTOPO30 dataset
        integer,  parameter :: im = 43200         ! total grids in x direction of 30-sec global dataset
        integer,  parameter :: jm = 21600         ! total grids in y direction of 30-sec global dataset
        real(r8), parameter :: dx = 1.0/120.0     ! space interval for 30-sec data (in degree)

        character (len=7)  :: nmtile(ntile)     ! name of each tile
        integer :: ncols,nrows                  ! number of columns and rows for 30-sec tile
        integer :: nodata                       ! integer for ocean point
        real(r8):: ulxmap       ! longitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8):: ulymap       ! latitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8):: lon_start    ! longitude at the center of grid (1,1) in the 30-sec netCDF global data
        real(r8):: lat_start    ! latitude at the center of grid (1,1) in the 30-sec netCDF global data
        real(r8):: lonsw        ! longitude at the center of southwest corner cell in the 30-sec tile
        real(r8):: latsw        ! latitude at the center of southwest corner cell in the 30-sec tile
        integer :: i1,j1        ! the (i,j) point of the southwest corner of the 30-sec tile in the global grid

        integer*2,  allocatable, dimension(:,:) :: terr          ! global 30-sec terrain data
        integer*1,  allocatable, dimension(:,:) :: land_fraction ! global 30-sec land fraction

        integer :: alloc_error,dealloc_error
        integer :: i,j,n                           ! index
        integer*2, allocatable, dimension(:,:)  :: iterr        ! terrain data for 30-sec tile
        integer*2,  allocatable, dimension(:,:) :: terr_tile    ! terrain data for 30-sec tile
        integer*1,  allocatable, dimension(:,:) :: land_fraction_tile
! 
        lat_start=-90.0 + 0.5 * dx
        lon_start=0.5*dx
        !
        ! Initialize each tile name
        !
        nmtile(1) = 'W180N90'
        nmtile(2) = 'W140N90'
        nmtile(3) = 'W100N90'
        nmtile(4) = 'W060N90'
        nmtile(5) = 'W020N90'
        nmtile(6) = 'E020N90'
        nmtile(7) = 'E060N90'
        nmtile(8) = 'E100N90'
        nmtile(9) = 'E140N90'

        nmtile(10) = 'W180N40'
        nmtile(11) = 'W140N40'
        nmtile(12) = 'W100N40'
        nmtile(13) = 'W060N40'
        nmtile(14) = 'W020N40'
        nmtile(15) = 'E020N40'
        nmtile(16) = 'E060N40'
        nmtile(17) = 'E100N40'
        nmtile(18) = 'E140N40'

        nmtile(19) = 'W180S10'
        nmtile(20) = 'W140S10'
        nmtile(21) = 'W100S10'
        nmtile(22) = 'W060S10'
        nmtile(23) = 'W020S10'
        nmtile(24) = 'E020S10'
        nmtile(25) = 'E060S10'
        nmtile(26) = 'E100S10'
        nmtile(27) = 'E140S10'
 
        nmtile(28) = 'W180S60'
        nmtile(29) = 'W120S60'
        nmtile(30) = 'W060S60'
        nmtile(31) = 'W000S60'
        nmtile(32) = 'E060S60'
        nmtile(33) = 'E120S60'


        allocate ( land_fraction(im,jm),stat=alloc_error )
        if( alloc_error /= 0 ) then
          print*,'Program could not allocate space for land_fraction'
          stop
        end if

        allocate ( terr(im,jm),stat=alloc_error )
        if( alloc_error /= 0 ) then
          print*,'Program could not allocate space for terr'
          stop
        end if
  
        do j = 1, jm
          do i = 1, im
            terr(i,j) = -999999.0
            land_fraction(i,j) = -99.0
          end do
        end do

        do n = 1,ntile
!
! Read header for each tile
!
           call rdheader(nmtile(n),nrows,ncols,nodata,ulxmap,ulymap)
 
!
! Allocate space for array iterr
!
           allocate ( iterr(ncols,nrows),stat=alloc_error )
           if( alloc_error /= 0 ) then
              print*,'Program could not allocate space for iterr'
              stop
           end if
!
! Read terr data for each tile
!
           call rdterr(nmtile(n),nrows,ncols,iterr)
!
! Allocate space for arrays terr_tile and psea10m
!
           allocate ( terr_tile(ncols,nrows),stat=alloc_error )
           if( alloc_error /= 0 ) then
              print*,'Program could not allocate space for terr_tile'
              stop
           end if
           allocate ( land_fraction_tile(ncols,nrows),stat=alloc_error )
           if( alloc_error /= 0 ) then
              print*,'Program could not allocate space for land_fraction_tile'
              stop
           end if
!
! Expand Caspian Sea for tiles 6 and 15
!
           if(nmtile(n).eq.'E020N90')call expand_sea(ncols,nrows,iterr,nodata,3600,5300)
           if(nmtile(n).eq.'E020N90')call expand_sea(ncols,nrows,iterr,nodata,4088,5874)
           if(nmtile(n).eq.'E020N40')call expand_sea(ncols,nrows,iterr,nodata,3600,1)
           print *, "min and maxiterr: ", minval(iterr), maxval(iterr)
!
! area average of 30-sec tile to 30-sec tile
!
           call avg(ncols,nrows,iterr,nodata,ulymap,dx,terr_tile,land_fraction_tile)

!
! Print some info on the fields
           print *, "min and max elevations: ", minval(terr_tile), maxval(terr_tile)
           print *, "min and max land_fraction: ", minval(land_fraction_tile), maxval(land_fraction_tile)
!
! fit the 30-sec tile into global 30-sec dataset
!

          latsw= ulymap - (nrows-1) * dx
          lonsw = ulxmap
          if( lonsw < 0.0 ) lonsw=360.0+lonsw
          i1 = nint( (lonsw - lon_start) / dx )+1
          if( i1 <= 0 ) i1 = i1 + im
          if( i1 > im ) i1 = i1 - im
          j1 = nint( (latsw- lat_start) / dx )+1
   
!          print*,'ulymap,ulxmap,latsw10,lonsw = ',ulymap,ulxmap,latsw10,lonsw
!          print*,'i1,j1 = ', i1,j1
 
          call fitin(ncols,nrows,terr_tile,land_fraction_tile,i1,j1,im,jm,terr,land_fraction)
!
! Deallocate working space for arrays iterr, terr_tile and psea10m
!
          deallocate ( iterr,terr_tile,land_fraction_tile,stat=dealloc_error )
          if( dealloc_error /= 0 ) then
             print*,'Unexpected deallocation error for arrays iterr,terr_tile'
             stop
          end if

        end do
        WRITE(*,*) 'done reading in USGS data'
!
! Print some info on the fields
        print *, "min and max elevations: ", minval(terr), maxval(terr)
        print *, "min and max land frac:  ", minval(land_fraction), maxval(land_fraction)
!
!  Write 30-sec terrain dataset, and land_fraction to NetCDF file
!
!        call wrtncdf(im,jm,terr,land_fraction,dx)
        call wrtncdf(im,jm,terr,land_fraction,dx,100)
      end program convterr

      subroutine rdheader(nmtile,nrows,ncols,nodata,ulxmap,ulymap)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine read the header of USGA Global30 sec TOPO data set.
!
        implicit none
!
! Dummy arguments
!
        character (len=7), intent(in) :: nmtile   ! name of the tile
        integer,  intent(out) :: nrows    ! number of rows
        integer,  intent(out) :: ncols    ! number of column
        integer,  intent(out) :: nodata   ! integer for ocean data point
        real(r8), intent(out) :: ulxmap
        real(r8), intent(out) :: ulymap
!
! Local variables
!
        character (len=11) :: flheader   ! file name of the header
        character (len=13) :: chars      ! dummy character
        
        flheader=nmtile//'.HDR'
 
        print*,'flheader = ', flheader
!
! Open GTOPO30 Header File
!
        open(unit=10,file=flheader,status='old',form='formatted')
!
! Read GTOPO30 Header file
!
        read (10, *)
        read (10, *)
        read (10, *) chars,nrows
        print*,chars,' = ',nrows
        read (10, *) chars,ncols
        print*,chars,' = ',ncols
        read (10, *)
        read (10, *)
        read (10, *)
        read (10, *)
        read (10, *)
        read (10, *) chars,nodata
        print*,chars,' = ',nodata
        read (10, *) chars,ulxmap
        print*,chars,' = ',ulxmap
        read (10, *) chars,ulymap
        print*,chars,' = ',ulymap
        close(10)

      end subroutine rdheader

      subroutine rdterr(nmtile,nrows,ncols,iterr)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine read the USGS Global 30-sec terrain data  for each tile.
!
        implicit none
!
! Dummy arguments
!
        character (len=7), intent(in) :: nmtile   ! name of the tile
        integer, intent(in) :: nrows    ! number of rows
        integer, intent(in) :: ncols    ! number of column
        integer*2, dimension(ncols,nrows), intent(out) :: iterr  ! terrain data
!
! Local variables
!
        character (len=11) :: flterr   ! file name for each terr dataset
        integer :: io_error            ! I/O status
        integer :: i,j                 ! Index
        integer :: length              ! record length
        
        flterr=nmtile//'.DEM'
 
!        print*,'flterr = ', flterr
!        print*,'nrows,ncols = ',nrows,ncols
!
! Open GTOPO30 Terrain dataset File
!

        length = 2 * ncols * nrows
        io_error=0
        open(unit=11,file=flterr,access='direct',recl=length,iostat=io_error)
        if( io_error /= 0 ) then
           print*,'Open file error in subroutine rdterr'
           print*,'iostat = ', io_error
           stop
        end if
!
! Read GTOPO30 Terrain data file
!
        read (11,rec=1,iostat=io_error) ((iterr(i,j),i=1,ncols),j=1,nrows)
!
        if( io_error /= 0 ) then
           print*,'Data file error in subroutine rdterr'
           print*,'iostat = ', io_error
           stop
        end if
!
! Print some info on the fields
           print *, "min and max elevations: ", minval(iterr), maxval(iterr)
!
! Correct missing data in source files
!
! Missing data near dateline

        if( nmtile == 'W180S60' ) then
           do j = 1, nrows
              iterr(1,j) = iterr(2,j)
           end do
        else if (nmtile == 'E120S60') then
           do j = 1, nrows
              iterr(ncols-1,j) = iterr(ncols-2,j)
              iterr(ncols,j)   = iterr(ncols-2,j)
           end do
        end if
!
! Missing data at the southermost row near South pole 
!
        if( nmtile == 'E060S60' .or.  nmtile == 'E120S60' .or. nmtile == 'W000S60' .or. &
            nmtile == 'W060S60' .or.  nmtile == 'W120S60' .or. nmtile == 'W180S60' ) then
           do i=1,ncols
              iterr(i,nrows) = iterr(i,nrows-1)
           end do
        end if
!
!        print*,'iterr(1,1),iterr(ncols,nrows) = ', &
!                iterr(1,1),iterr(ncols,nrows)
 
        close (11)
      end subroutine rdterr

      subroutine avg(ncols,nrows,iterr,nodata,ulymap,dx,terr_tile,land_fraction_tile)
        use shr_kind_mod, only: r8 => shr_kind_r8
        implicit none
!
! Dummy arguments
!
        integer, intent(in) :: ncols    ! number of column for 30-sec tile
        integer, intent(in) :: nrows    ! number of rows for 30-sec tile
        integer*2, dimension(ncols,nrows), intent(inout) :: iterr    ! terrain data for 30-sec tile
        integer, intent(in) :: nodata   ! integer for ocean data point
        real(r8),intent(in) :: ulymap   ! latitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8),intent(in) :: dx    ! spacing interval for 30-sec data (in degree)
        integer*2, dimension(ncols,nrows), intent(out) :: terr_tile    !  terrain data for 30-sec tile
        integer*1, dimension(ncols,nrows), intent(out) :: land_fraction_tile  
!
! Local variables
!
        real(r8) :: lats,latn    ! latitudes (in rad) for ths south and north edges of each 30-sec cell
        real(r8) :: wt           ! area weighting of each 30-sec cell
        real(r8) :: wt_tot       ! total weighting of each 30-sec cell
        real(r8) :: sumterr      ! summation of terrain height of each 30-sec cell
        real(r8) :: sumsea       ! summation of sea coverage of each 30-sec cell
        real(r8) :: pi           ! pi=3.1415
        real(r8) :: latul        ! latitude of the upper-left coner of 30-sec tile
        integer :: n1,itmp,i1,i2,j1,j2    ! temporary working spaces
        integer :: i,j,ii,jj              ! index
        logical, dimension(ncols,nrows) :: oflag
 
        pi = 4.0 * atan(1.0)
!
        n1 = ncols / ncols
        print*,'ncols,ncols,n1 = ',ncols,ncols,n1

        itmp = nint( ulymap + 0.5 * dx )
        latul = itmp
        print*,'ulymap,latul = ', ulymap,latul
        oflag = .false.

        do j = 1, nrows
           j1 = j
           j2 = j
           do i = 1, ncols
              i1 = i
              i2 = i   
              terr_tile(i,j) = 0
              land_fraction_tile(i,j) = 1
              if ( iterr(i,j) == nodata ) then
                land_fraction_tile(i,j) = 0
              else
                if ( iterr(i,j) .lt.nodata ) then
                  ! this can only happen in the expand_sea routine
                  land_fraction_tile(i,j) = 0
                  iterr(i,j) = iterr(i,j) - nodata - nodata
                endif
                terr_tile(i,j) = iterr(i,j)
              end if
            end do
          end do

      end subroutine avg

      subroutine expand_sea(ncols,nrows,iterr,nodata,startx,starty)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine reduces the resolution of the terrain data from 30-sec to 30-sec and
! compute the percentage of ocean cover (psea10m)
!
        implicit none
!
! Dummy arguments
!
        integer, intent(in) :: ncols    ! number of column for 30-sec tile
        integer, intent(in) :: nrows    ! number of rows for 30-sec tile
        integer*2, dimension(ncols,nrows), intent(inout) :: iterr    ! terrain data for 30-sec tile
        integer, intent(in) :: nodata   ! integer for ocean data point
        integer, intent(in) :: startx, starty ! where to begin the sea
!
! Local variables
!
        real(r8):: maxh
        integer :: i,j,per,ii,jj       ! index
        logical, dimension(0:ncols+1,0:nrows+1) :: flag  ! terrain data for 30-sec tile
        logical :: found

        flag = .false.

        maxh = iterr(startx,starty)
        
        iterr(startx,starty) =  iterr(startx,starty) + nodata + nodata
        flag(startx-1:startx+1,starty-1:starty+1) = .true.

        per = 0
        print *, 'expanding sea at ',maxh,' m '

2112    per = per + 1
        found = .false.
        do j = starty - per, starty + per, per*2
           do i = startx - per, startx + per
              if(i.ge.1.and.i.le.ncols.and.j.ge.1.and.j.le.nrows)then
                 if( iterr(i,j).eq.maxh .and. flag(i,j) ) then
                    iterr(i,j) =  iterr(i,j) + nodata + nodata
                    flag(i-1:i+1,j-1:j+1) = .true.
                    found = .true.
                 endif
              endif
           end do
        end do

        do i = startx - per, startx + per, per*2
           do j = starty - per + 1, starty + per - 1
              if(i.ge.1.and.i.le.ncols.and.j.ge.1.and.j.le.nrows)then
                 if( iterr(i,j).eq.maxh .and. flag(i,j) ) then
                    iterr(i,j) =  iterr(i,j) + nodata + nodata
                    flag(i-1:i+1,j-1:j+1) = .true.
                    found = .true.
                 endif
              endif
           end do
        end do
        if (found)goto 2112
        print *, 'done with expand_sea'
        return

      end subroutine expand_sea

      subroutine fitin(ncols,nrows,terr_tile,land_fraction_tile,i1,j1,im,jm,terr,land_fraction)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine put 30-sec tile into the global dataset
!
        implicit none
!
! Dummy arguments
!
        integer, intent(in) :: ncols       ! number of columns for 30-sec tile
        integer, intent(in) :: nrows       ! number of rows for 30-sec tile
        integer*2, dimension(ncols,nrows), intent(in) :: terr_tile    !  terrain data for 30-sec tile
        integer*1, dimension(ncols,nrows), intent(in) :: land_fraction_tile    
        integer, intent(in) :: i1,j1        ! the (i,j) point of the southwest corner of the 30-sec tile 
                                            ! in the global grid
        integer, intent(in) :: im,jm    ! the dimensions of the 30-sec global dataset
        integer*2,dimension(im,jm), intent(out) :: terr           ! global 30-sec terrain data
        integer*1,dimension(im,jm), intent(out) :: land_fraction  ! global 30-sec land fraction
!
! Local variables
!
        integer :: i,j,ii,jj              ! index

        do j = 1, nrows
           jj = j1 + (nrows - j)
           do i = 1, ncols
              ii = i1 + (i-1)

              if( i == 1 .and. j == 1 ) &
                 print*,'i,j,ii,jj = ',i,j,ii,jj
              if( i == ncols .and. j == nrows ) &
                 print*,'i,j,ii,jj = ',i,j,ii,jj

              if( ii > im ) ii = ii - im
              terr(ii,jj) = terr_tile(i,j)
              land_fraction(ii,jj) = land_fraction_tile(i,j)
           end do
        end do
      end subroutine fitin

      subroutine wrtncdf(im,jm,terr,land_fraction,dx)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine save 30-sec terrain data, land fraction to NetCDF file
!
        implicit none

#     include         <netcdf.inc>

!
! Dummy arguments
!
        integer, intent(in) :: im,jm    ! the dimensions of the 30-sec global dataset
        integer*2,dimension(im,jm), intent(in) :: terr       ! global 30-sec terrain data
        integer*1,dimension(im,jm), intent(in) :: land_fraction !global 30-sec land fraction
        real(r8),                      intent(in) :: dx
!
! Local variables
!
        real(r8),dimension(im) :: lonar       ! longitude array
        real(r8),dimension(im) :: latar       ! latitude array
        character (len=32) :: fout       ! NetCDF output file
        integer            :: foutid     ! Output file id
        integer            :: lonid, lonvid
        integer            :: latid, latvid
        integer            :: htopoid
        integer            :: landfid
        integer, dimension(2) :: htopodim,landfdim
        integer            :: status    ! return value for error control of netcdf routin
        integer            :: i,j
        character (len=8)  :: datestring

        integer*2,dimension(im,jm) :: h       ! global 30-sec terrain data
        integer*1,dimension(im,jm) :: lnd


!
! Fill lat and lon arrays
!  
        do i = 1,im
           lonar(i)=  dx * (i-0.5)
        enddo
        do j = 1,jm
           latar(j)= -90.0 + dx * (j-0.5)
        enddo

        do j=1,jm
          do i=1,im
            h(i,j) = terr(i,j)
            lnd(i,j) = land_fraction(i,j)
          end do
        end do

        fout='usgs-rawdata.nc'
!
!  Create NetCDF file for output
!
        print *,"Create NetCDF file for output"
        status = nf_create (fout, NF_64BIT_OFFSET , foutid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create dimensions for output
!
        print *,"Create dimensions for output"
        status = nf_def_dim (foutid, 'lon', im, lonid)
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_def_dim (foutid, 'lat', jm, latid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create variable for output
!
        print *,"Create variable for output"
        htopodim(1)=lonid
        htopodim(2)=latid
        status = nf_def_var (foutid,'htopo', NF_INT, 2, htopodim, htopoid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
        landfdim(1)=lonid
        landfdim(2)=latid
        status = nf_def_var (foutid,'landfract', NF_INT, 2, landfdim, landfid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lat', NF_DOUBLE, 1, latid, latvid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lon', NF_DOUBLE, 1, lonid, lonvid)
        if (status .ne. NF_NOERR) call handle_err(status)

!
! Create attributes for output variables
!
        status = nf_put_att_text (foutid,htopoid,'long_name', 41, '30-sec elevation from USGS 30-sec dataset')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,htopoid,'units', 5, 'meter')
        if (status .ne. NF_NOERR) call handle_err(status)
!
        status = nf_put_att_text (foutid,landfid,'long_name', 23, '30-second land fraction')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,landfid,'units', 14, 'fraction (0-1)')
        if (status .ne. NF_NOERR) call handle_err(status)
!
        status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,latvid,'units', 21, 'cell center locations')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,lonvid,'units' , 21, 'cell center locations')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,NF_GLOBAL,'source', 27, 'USGS 30-sec dataset GTOPO30')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,NF_GLOBAL,'title',  24, '30-second USGS topo data')
        if (status .ne. NF_NOERR) call handle_err(status)
        call DATE_AND_TIME(DATE=datestring)
        status = nf_put_att_text (foutid,NF_GLOBAL,'history',25, 'Written on date: ' // datestring )
        if (status .ne. NF_NOERR) call handle_err(status)
        
!
! End define mode for output file
!
        status = nf_enddef (foutid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Write variable for output
!
        print*,"writing terrain data"
        status = nf_put_var_int2 (foutid, htopoid, h)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing terrain data"
!
        status = nf_put_var_int1 (foutid, landfid, lnd)
        if (status .ne. NF_NOERR) call handle_err(status)
!
        print*,"writing lat data"
        status = nf_put_var_double (foutid, latvid, latar)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing lat data"

        print*,"writing lon data"
        status = nf_put_var_double (foutid, lonvid, lonar)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing lon data"
!
! Close output file
!
        print *,"close file"
        status = nf_close (foutid)
        if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine wrtncdf


      !
      ! same as wrtncdf but the output is coarsened
      !
      subroutine wrtncdf_coarse(im,jm,terr,land_fraction,dx,ic)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine save 30-sec terrain data, land fraction to NetCDF file
!
        implicit none

#     include         <netcdf.inc>

!
! Dummy arguments
!
        integer, intent(in) :: im,jm    ! the dimensions of the 30-sec global dataset
        integer, intent(in) :: ic       ! coarsening factor
        integer*2,dimension(im,jm), intent(in) :: terr       ! global 30-sec terrain data
        integer*1,dimension(im,jm), intent(in) :: land_fraction !global 30-sec land fraction
        real(r8),                      intent(in) :: dx
!
! Local variables
!
        real(r8),dimension(im/ic) :: lonar       ! longitude array
        real(r8),dimension(im/ic) :: latar       ! latitude array
        character (len=32) :: fout       ! NetCDF output file
        integer            :: foutid     ! Output file id
        integer            :: lonid, lonvid
        integer            :: latid, latvid
        integer            :: htopoid
        integer            :: landfid
        integer, dimension(2) :: htopodim,landfdim
        integer            :: status    ! return value for error control of netcdf routin
        integer            :: i,j
        character (len=8)  :: datestring

        integer*2,dimension(im/ic,jm/ic) :: h       ! global 30-sec terrain data
        integer*1,dimension(im/ic,jm/ic) :: lnd


!
! Fill lat and lon arrays
!  
        do i = 1,im/ic
           lonar(i)=  real(ic)*dx * (i-0.5)
        enddo
        do j = 1,jm/ic
           latar(j)= -90.0 + real(ic)*dx * (j-0.5)
        enddo

        do j=1,jm/ic
          do i=1,im/ic
            h(i,j) = terr(i*ic,j*ic)
            lnd(i,j) = land_fraction(i*ic,j*ic)
          end do
        end do

        fout='usgs-lowres.nc'
!
!  Create NetCDF file for output
!
        print *,"Create NetCDF file for output"
        status = nf_create (fout, NF_64BIT_OFFSET , foutid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create dimensions for output
!
        print *,"Create dimensions for output"
        status = nf_def_dim (foutid, 'lon', im/ic, lonid)
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_def_dim (foutid, 'lat', jm/ic, latid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create variable for output
!
        print *,"Create variable for output"
        htopodim(1)=lonid
        htopodim(2)=latid
        status = nf_def_var (foutid,'htopo', NF_INT, 2, htopodim, htopoid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
        landfdim(1)=lonid
        landfdim(2)=latid
        status = nf_def_var (foutid,'landfract', NF_INT, 2, landfdim, landfid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lat', NF_DOUBLE, 1, latid, latvid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lon', NF_DOUBLE, 1, lonid, lonvid)
        if (status .ne. NF_NOERR) call handle_err(status)

!
! Create attributes for output variables
!
        status = nf_put_att_text (foutid,htopoid,'long_name', 41, '30-sec elevation from USGS 30-sec dataset')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,htopoid,'units', 5, 'meter')
        if (status .ne. NF_NOERR) call handle_err(status)
!
        status = nf_put_att_text (foutid,landfid,'long_name', 23, '30-second land fraction')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,landfid,'units', 14, 'fraction (0-1)')
        if (status .ne. NF_NOERR) call handle_err(status)
!
        status = nf_put_att_text (foutid,latvid,'long_name', 8, 'latitude')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,latvid,'units', 13, 'degrees_north')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,latvid,'units', 21, 'cell center locations')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,lonvid,'long_name', 9, 'longitude')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,lonvid,'units', 12, 'degrees_east')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,lonvid,'units' , 21, 'cell center locations')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,NF_GLOBAL,'source', 27, 'USGS 30-sec dataset GTOPO30')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,NF_GLOBAL,'title',  24, '30-second USGS topo data')
        if (status .ne. NF_NOERR) call handle_err(status)
        call DATE_AND_TIME(DATE=datestring)
        status = nf_put_att_text (foutid,NF_GLOBAL,'history',25, 'Written on date: ' // datestring )
        if (status .ne. NF_NOERR) call handle_err(status)
        
!
! End define mode for output file
!
        status = nf_enddef (foutid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Write variable for output
!
        print*,"writing terrain data"
        status = nf_put_var_int2 (foutid, htopoid, h)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing terrain data"
!
        status = nf_put_var_int1 (foutid, landfid, lnd)
        if (status .ne. NF_NOERR) call handle_err(status)
!
        print*,"writing lat data"
        status = nf_put_var_double (foutid, latvid, latar)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing lat data"

        print*,"writing lon data"
        status = nf_put_var_double (foutid, lonvid, lonar)
        if (status .ne. NF_NOERR) call handle_err(status)
        print*,"done writing lon data"
!
! Close output file
!
        print *,"close file"
        status = nf_close (foutid)
        if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine wrtncdf_coarse
!************************************************************************
!!handle_err
!************************************************************************
!
!!ROUTINE:      handle_err
!!DESCRIPTION:  error handler
!--------------------------------------------------------------------------

      subroutine handle_err(status)

        implicit         none

#     include          <netcdf.inc>

        integer          status

        if (status .ne. nf_noerr) then
          print *, nf_strerror(status)
          stop 'Stopped'
        endif

      end subroutine handle_err



