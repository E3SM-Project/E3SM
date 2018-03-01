!
!  DATE CODED:   Oct 17, 2000
!  DESCRIPTION:  This program reads USGS 30-sec terrain dataset in 33 tiles and converts
!                them to 10-min resolution global dataset in one single NetCDF file.
!
!  Author: Jiundar Chern (jchern@dao.gsfc.nasa.gov)
!
!  ** Modified November, 2003 ***
!  This code has been modified by Jim McCaa (jmccaa@ucar.edu) for use at NCAR.  
!  In particular:
!  1) Paths and compiler options have been changed.
!  2) The code now generates a Caspian Sea based on elevation, and reports these points
!     as ocean.  This is done through three calls to the new routine expand_sea.
!
!  ** Modified February 4, 2005 B.A. Boville ***
!
!  ROUTINES CALLED:
!       netcdf routines
!
!  COMPILING:
!
! NCAR SGI (chinookfe) f90 -I/usr/local/include -O -64 -mips4 -bytereclen -s
!          -o gtopo30_to_10min gtopo30_to_10min.F90 -L/usr/local/lib64/r4i4 -lnetcdf -r8

! NASA DAO SGI: f90 -I/ford1/local/IRIX64/netcdf/include  -O -64 -mips4 -bytereclen -s
!          -o gtopo30_to_10min  gtopo30_to_10min.F90 -L/ford1/local/IRIX64/netcdf/lib -lnetcdf -r8

      program convterr
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This program converts USGS 30-sec terrain data set to 10-min resolution
! terrain data set.
!
        implicit none
!
        integer,  parameter :: ntile = 33        ! number of tiles in USGS GTOPO30 dataset
        integer,  parameter :: im10 = 2160       ! total grids in x direction of 10-min global dataset
        integer,  parameter :: jm10 = 1080       ! total grids in y direction of 10-min global dataset
        real(r8), parameter :: dx30s = 1.0/120.0 ! space interval for 30-sec data (in degree)
        real(r8), parameter :: dx10m = 1.0/6.0   ! space interval for 10-min data (in degree)
 
        character (len=7)  :: nmtile(ntile)     ! name of each tile
        integer :: ncols,nrows                  ! number of columns and rows for 30-sec tile
        integer :: nodata                       ! integer for ocean point
        integer :: ncol10,nrow10                ! number of columns and rows for 10-min tile
        real(r8):: ulxmap       ! longitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8):: ulymap       ! latitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8):: lon1_10m     ! longitude at the center of grid (1,1) in the 10-min global data
        real(r8):: lat1_10m     ! latitude at the center of grid (1,1) in the 10-min global data
        real(r8):: lonsw10      ! longitude at the center of southwest corner cell in the 10-min tile
        real(r8):: latsw10      ! latitude at the center of southwest corner cell in the 10-min tile
        integer :: i1,j1        ! the (i,j) point of the southwest corner of the 10-min tile in the global grid
        real(r8), dimension(im10,jm10)    :: terr          ! global 10-min terrain data
        real(r8), dimension(im10,jm10)    :: variance      ! global 10-min variance of elevation
        real(r8), dimension(im10,jm10)    :: land_fraction !global 10-min land fraction

        integer :: alloc_error,dealloc_error
        integer :: i,j,n                           ! index
        integer*2, allocatable, dimension(:,:) :: iterr     ! terrain data for 30-sec tile
        real(r8),  allocatable, dimension(:,:) :: terr10m   ! terrain data for 10-min tile
        real(r8),  allocatable, dimension(:,:) :: psea10m   ! percentage of ocaen for 10-min tile
        real(r8),  allocatable, dimension(:,:) :: var10m    ! variance of 30-sec elevations for 10-min tile
! 
        lat1_10m=-90.0 + 0.5 * dx10m
        lon1_10m=0.5*dx10m
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
  
        do j = 1, jm10
          do i = 1, im10
            terr(i,j) = -9999.0
            variance(i,j) = -9999.0
            land_fraction(i,j) = -9999.0
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
! Allocate space for arrays terr10m and psea10m
!
           nrow10 =nrows*dx30s/dx10m
           ncol10 =ncols*dx30s/dx10m
           allocate ( terr10m(ncol10,nrow10),psea10m(ncol10,nrow10),var10m(ncol10,nrow10),stat=alloc_error )
           if( alloc_error /= 0 ) then
              print*,'Program could not allocate space for terr10m, psea10m, and var10m'
              stop
           end if
!
! Expand Caspian Sea for tiles 6 and 15
!
           if(nmtile(n).eq.'E020N90')call expand_sea(ncols,nrows,iterr,nodata,3600,5300)
           if(nmtile(n).eq.'E020N90')call expand_sea(ncols,nrows,iterr,nodata,4088,5874)
           if(nmtile(n).eq.'E020N40')call expand_sea(ncols,nrows,iterr,nodata,3600,1)
!
! area average of 30-sec tile to 10-min tile
!
           call avg(ncols,nrows,iterr,nodata,ulymap,dx30s,ncol10,nrow10,terr10m,psea10m,var10m)

!
! Print some info on the fields
           print *, "min and max elevations: ", minval(terr10m), maxval(terr10m)
           print *, "min and max variacnes:  ", minval(var10m) , maxval(var10m)
           print *, "min and max land frac:  ", minval(psea10m), maxval(psea10m)
!
! fit the 10-min tile into global 10-min dataset
! Note: the 30-sec and 10-min tiles are scaned from north to south, the global 10-min dataset are 
!       scaned from south to north (90S to 90N) and east to west (0E to -0.1666667W) 
!
          latsw10 = nint(ulymap + 0.5 * dx30s) - nrow10 * dx10m + 0.5 * dx10m
          lonsw10 = nint(ulxmap - 0.5 * dx30s) + 0.5 * dx10m
          if( lonsw10 < 0.0 ) lonsw10=360.0+lonsw10
          i1 = nint( (lonsw10 - lon1_10m) / dx10m )+1
          if( i1 <= 0 ) i1 = i1 + im10
          if( i1 > im10 ) i1 = i1 - im10
          j1 = nint( (latsw10 - lat1_10m) / dx10m )+1
   
!          print*,'ulymap,ulxmap,latsw10,lonsw10 = ',ulymap,ulxmap,latsw10,lonsw10
!          print*,'i1,j1 = ', i1,j1
 
          call fitin(ncol10,nrow10,terr10m,psea10m,var10m,i1,j1,im10,jm10,terr,variance,land_fraction)
!
! Deallocate working space for arrays iterr, terr10m and psea10m
!
          deallocate ( iterr,terr10m,psea10m,var10m,stat=dealloc_error )
          if( dealloc_error /= 0 ) then
             print*,'Unexpected deallocation error for arrays iterr,terr10m,psea10m,var10m'
             stop
          end if

        end do

!
! Print some info on the fields
        print *, "min and max elevations: ", minval(terr), maxval(terr)
        print *, "min and max variances:  ", minval(variance), maxval(variance)
        print *, "min and max land frac:  ", minval(land_fraction), maxval(land_fraction)
!
!  Write 10-min terrain dataset, variance  and land_fraction to NetCDF file
!
        call wrtncdf(im10,jm10,terr,variance, land_fraction,dx10m)

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

      subroutine avg(ncols,nrows,iterr,nodata,ulymap,dx30s,ncol10,nrow10,terr10m,psea10m,var10m)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine reduces the resolution of the terrain data from 30-sec to 10-min and
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
        real(r8),intent(in) :: ulymap   ! latitude at the center of the upper-left corner cell in the 30-sec tile
        real(r8),intent(in) :: dx30s    ! spacing interval for 30-sec data (in degree)
        integer, intent(in) :: nrow10   ! number of rows for 10-min tile
        integer, intent(in) :: ncol10   ! number of columns for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(out) :: terr10m    !  terrain data for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(out) :: psea10m    !  percentage ocean coverage for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(out) :: var10m     !  variance of 30-sec elevations
!
! Local variables
!
        real(r8) :: lats,latn    ! latitudes (in rad) for ths south and north edges of each 30-sec cell
        real(r8) :: wt           ! area weighting of each 30-sec cell
        real(r8) :: wt_tot       ! total weighting of each 10-min cell
        real(r8) :: sumterr      ! summation of terrain height of each 10-min cell
        real(r8) :: sumsea       ! summation of sea coverage of each 10-min cell
        real(r8) :: pi           ! pi=3.1415
        real(r8) :: latul        ! latitude of the upper-left coner of 30-sec tile
        integer :: n1,itmp,i1,i2,j1,j2    ! temporary working spaces
        integer :: i,j,ii,jj              ! index
        logical, dimension(ncols,nrows) :: oflag
 
        pi = 4.0 * atan(1.0)
!
        n1 = ncols / ncol10
        print*,'ncols,ncol10,n1 = ',ncols,ncol10,n1

        itmp = nint( ulymap + 0.5 * dx30s )
        latul = itmp
        print*,'ulymap,latul = ', ulymap,latul
        oflag = .false.

        do j = 1, nrow10
           j1 = (j-1) * n1 + 1
           j2 = j     * n1
           do i = 1, ncol10
              i1 = (i-1) * n1 + 1
              i2 = i     * n1
              wt_tot = 0.0
              sumterr = 0.0
              sumsea = 0.0

              do jj = j1, j2
                 latn = ( latul - (jj -1) * dx30s ) * pi / 180.0
                 lats = ( latul - jj      * dx30s ) * pi / 180.0
                 wt = sin( latn ) - sin( lats )
      
                 do ii = i1, i2
                    wt_tot=wt_tot+wt
                    if ( iterr(ii,jj) == nodata ) then
                       sumsea = sumsea + wt
                       oflag(ii,jj) = .true.
                    else
                       if ( iterr(ii,jj) .lt.nodata ) then
                          ! this can only happen in the expand_sea routine
                          sumsea = sumsea + wt
                          oflag(ii,jj) = .true.
                          iterr(ii,jj) = iterr(ii,jj) - nodata - nodata
                       endif
                       sumterr = sumterr + iterr(ii,jj) * wt
                    end if
                 end do
              end do

              terr10m(i,j) = sumterr / wt_tot
              psea10m(i,j) = sumsea / wt_tot

           end do
        end do

        ! Now compute variance of 30-second points
        
        do j = 1, nrow10
           j1 = (j-1) * n1 + 1
           j2 = j     * n1

           do i = 1, ncol10
              i1 = (i-1) * n1 + 1
              i2 = i     * n1

              wt_tot = 0.0
              var10m(i,j)  = 0.0
              wt = 1.0
              do jj = j1, j2
                 do ii = i1, i2
                    wt_tot = wt_tot + wt
                    if ( .not. oflag(ii,jj) ) then
                       var10m(i,j) = var10m(i,j) + wt * (iterr(ii,jj)-terr10m(i,j))**2
                    end if
                 end do
              end do
              var10m(i,j) = var10m(i,j) / wt_tot

           end do
        end do

      end subroutine avg

      subroutine expand_sea(ncols,nrows,iterr,nodata,startx,starty)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine reduces the resolution of the terrain data from 30-sec to 10-min and
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

      subroutine fitin(ncol10,nrow10,terr10m,psea10m,var10m,i1,j1,im10,jm10,terr,variance,land_fraction)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine put 10-min tile into the global dataset
!
        implicit none
!
! Dummy arguments
!
        integer, intent(in) :: ncol10       ! number of columns for 10-min tile
        integer, intent(in) :: nrow10       ! number of rows for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(in) :: terr10m    !  terrain data for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(in) :: psea10m    !  percentage ocean coverage for 10-min tile
        real(r8), dimension(ncol10,nrow10), intent(in) :: var10m     !  variance of 30-sec elev for 10-min tile
        integer, intent(in) :: i1,j1        ! the (i,j) point of the southwest corner of the 10-min tile 
                                            ! in the global grid
        integer, intent(in) :: im10,jm10    ! the dimensions of the 10-min global dataset
        real(r8),dimension(im10,jm10), intent(out) :: terr           ! global 10-min terrain data
        real(r8),dimension(im10,jm10), intent(out) :: variance       ! global 10-min variance of elev
        real(r8),dimension(im10,jm10), intent(out) :: land_fraction  ! global 10-min land fraction
!
! Local variables
!
        integer :: i,j,ii,jj              ! index

        do j = 1, nrow10
           jj = j1 + (nrow10 - j)
           do i = 1, ncol10
              ii = i1 + (i-1)
              if( ii > im10 ) ii = ii - im10
              terr(ii,jj) = terr10m(i,j)
              land_fraction(ii,jj) = 1.0 - psea10m(i,j)
              variance(ii,jj) = var10m(i,j)
              if( i == 1 .and. j == 1 ) &
                 print*,'i,j,ii,jj = ',i,j,ii,jj
              if( i == ncol10 .and. j == nrow10 ) &
                 print*,'i,j,ii,jj = ',i,j,ii,jj
           end do
        end do
      end subroutine fitin

      subroutine wrtncdf(im10,jm10,terr,variance,land_fraction,dx10m)
        use shr_kind_mod, only: r8 => shr_kind_r8
!
! This subroutine save 10-min terrain data, variance, land fraction to NetCDF file
!
        implicit none

#     include         <netcdf.inc>

!
! Dummy arguments
!
        integer, intent(in) :: im10,jm10    ! the dimensions of the 10-min global dataset
        real(r8),dimension(im10,jm10), intent(in) :: terr       ! global 10-min terrain data
        real(r8),dimension(im10,jm10), intent(in) :: variance       ! global 10-min variance data
        real(r8),dimension(im10,jm10), intent(in) :: land_fraction !global 10-min land fraction
        real(r8),                      intent(in) :: dx10m
!
! Local variables
!
        real(r8),dimension(im10) :: lonar       ! longitude array
        real(r8),dimension(im10) :: latar       ! latitude array
        character (len=32) :: fout       ! NetCDF output file
        integer            :: foutid     ! Output file id
        integer            :: lonid, lonvid
        integer            :: latid, latvid
        integer            :: varianceid
        integer            :: htopoid
        integer            :: landfid
        integer, dimension(2) :: variancedim,htopodim,landfdim
        integer            :: status    ! return value for error control of netcdf routin
        integer            :: i,j
        character (len=8)  :: datestring

!
! Fill lat and lon arrays
!  
        do i = 1,im10
           lonar(i)=  dx10m * (i-0.5)
        enddo
        do j = 1,jm10
           latar(j)= -90.0 + dx10m * (j-0.5)
        enddo

        fout='topo_gtopo30_10min.nc'
!
!  Create NetCDF file for output
!
        status = nf_create (fout, NF_WRITE, foutid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create dimensions for output
!
        status = nf_def_dim (foutid, 'lon', im10, lonid)
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_def_dim (foutid, 'lat', jm10, latid)
        if (status .ne. NF_NOERR) call handle_err(status)
!
! Create variable for output
!
        variancedim(1)=lonid
        variancedim(2)=latid
        status = nf_def_var (foutid,'variance', NF_FLOAT, 2, variancedim, varianceid)
        if (status .ne. NF_NOERR) call handle_err(status)

        htopodim(1)=lonid
        htopodim(2)=latid
        status = nf_def_var (foutid,'htopo', NF_FLOAT, 2, htopodim, htopoid)
        if (status .ne. NF_NOERR) call handle_err(status)

        landfdim(1)=lonid
        landfdim(2)=latid
        status = nf_def_var (foutid,'landfract', NF_FLOAT, 2, landfdim, landfid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lat', NF_DOUBLE, 1, latid, latvid)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_def_var (foutid,'lon', NF_DOUBLE, 1, lonid, lonvid)
        if (status .ne. NF_NOERR) call handle_err(status)

!
! Create attributes for output variables
!
        status = nf_put_att_text (foutid,varianceid,'long_name', 29, 'variance of 30-sec elevations')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,varianceid,'units', 8, 'meter**2')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,htopoid,'long_name', 41, '10-min elevation from USGS 30-sec dataset')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,htopoid,'units', 5, 'meter')
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_att_text (foutid,landfid,'long_name', 23, '10-minute land fraction')
        if (status .ne. NF_NOERR) call handle_err(status)
        status = nf_put_att_text (foutid,landfid,'units', 14, 'fraction (0-1)')
        if (status .ne. NF_NOERR) call handle_err(status)

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
        status = nf_put_att_text (foutid,NF_GLOBAL,'title',  24, '10-minute USGS topo data')
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
        status = nf_put_var_double (foutid, varianceid, variance)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_var_double (foutid, htopoid, terr)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_var_double (foutid, landfid, land_fraction)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_var_double (foutid, latvid, latar)
        if (status .ne. NF_NOERR) call handle_err(status)

        status = nf_put_var_double (foutid, lonvid, lonar)
        if (status .ne. NF_NOERR) call handle_err(status)
 
!
! Close output file
!
        status = nf_close (foutid)
        if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine wrtncdf
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


