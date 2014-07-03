program regrid
!----------------------------------------------------------------------------------
!
! Purpose: Area average or linearly interpolate 1x1 degree SST and sea ice concentration
!          data to grid defined by output file.  
!
! Usage: regrid [-a] [-b begin] [-e end] -g gridfile -i icefile [-l] -o outfile -s sstfile [-v]
!
!----------------------------------------------------------------------------------
   use precision

   implicit none

   include 'netcdf.inc'
!
! Local workspace
!
   character(len=256) :: sstfilein =  ' '    ! filename for 1x1 input SST data
   character(len=256) :: icefilein =  ' '    ! filename for 1x1 input ICE data
   character(len=256) :: gridout = ' '       ! filename for output grid coordinates
   character(len=256) :: fileout = ' '       ! filename for regridded output data
   character(len=256) :: arg                 ! cmd line argument
   character(len=512) :: cmdline             ! input command line
   character(len=23)  :: scattername         ! file name for scatter plot
   character(len=19)  :: cur_timestamp
   character(len=600) :: history             ! history attribute
   character(len=256) :: units_att           ! units attribute
   character(len=3)   :: time_coord = 'old'  ! 'new' or 'old'

   logical :: uselandmask = .false.     ! Flag says to use input land mask (no longer used)
   logical :: verbose = .false.         ! Add print statements
   logical :: linear_interp             ! Flag to do linear interpolation vs. area averaging
   logical :: usr_specified_interp = .false. ! Whether cmd-line interp type specified
   logical :: scatterplots = .false.    ! Turn on for scatter plots

   integer :: cmdlen, hislen, totlen    ! character array lengths
   integer :: date(1)                   ! array to please lf95
   integer :: datesec(1)                ! array to please lf95
   integer :: mo                        ! month index

   real(r8), parameter :: iceflag = -1.8  ! sea ice flag
   real(r8) :: pertval                    ! randomizer for scatter plots
   real(r8) :: fillvaluesst               ! in case sst has fillvalue over land
   real(r8) :: fillvalueice               ! in case ice has fillvalue over land
   real(r4), parameter :: fillvalueout = 1.e36 ! placeholder in case land mask desired
!
! Netcdf file, dimension, and variable id's for input file
!
   integer :: ncidsstin = -1            ! input sst file
   integer :: ncidicein = -1            ! input sst file
   integer :: ncidout = -1              ! output file
   integer :: londimidin = -1           ! longitude dimension
   integer :: londimidicein = -1        ! longitude dimension ice file (length compared vs. sst)
   integer :: latdimidin = -1           ! latitude dimension
   integer :: latdimidicein = -1        ! latitude dimension ice file (length compared vs. sst)
   integer :: timedimidin = -1          ! time dimension
   integer :: timedimidicein = -1       ! time dimension ice file (compared vs. sst)
   integer :: lonidin = -1              ! longitude variable
   integer :: latidin = -1              ! latitude variable
   integer :: timeidin = -1             ! time variable
   integer :: timeidicein = -1          ! time variable from ice, checked vs. sst
   integer :: dateidin = -1             ! date variable
   integer :: ntimein = -1              ! number of time samples
   integer :: ntimeicein = -1           ! number of time samples ice file, checked vs. sst
   integer :: sstidin = -1              ! SST variable
   integer :: iceidin = -1              ! ICE variable
!
! Netcdf file, dimension, and variable id's for grid file
!
   integer :: ncidgrd = -1              ! grid file
   integer :: londimidgrd = -1          ! longitude dimension
   integer :: latdimidgrd = -1          ! latitude dimension
   integer :: lonidgrd = -1             ! longitude variable
   integer :: latidgrd = -1             ! latitude variable
   logical :: is_reduced_grid = .false.
   integer :: nlonidgrd = -1            ! number of longitudes (reduced grid)
   integer :: rlonidgrd = -1            ! longitude variable (reduced grid)
!
! Netcdf file, dimension, and variable id's for output file
!
   integer :: londimidout = -1          ! longitude dimension
   integer :: latdimidout = -1          ! latitude dimension
   integer :: timedimidout = -1         ! time dimension
   integer :: lonidout = -1             ! longitude variable
   integer :: latidout = -1             ! latitude variable
   integer :: nlonidout = -1            ! number of longitudes (reduced grid)
   integer :: rlonidout = -1            ! longitude variable (reduced grid)
   integer :: timeidout = -1            ! time variable
   integer :: dateidout = -1            ! date variable for model
   integer :: datesecidout = -1         ! seconds of date variable for model
   integer :: sstidout = -1             ! SST variable
   integer :: iceidout = -1             ! ICE fraction variable
   integer :: oroidsst = -1             ! ORO variable based on SST
   integer :: oroidconc = -1            ! ORO variable based on sea ice concentration
   integer :: landmaskidout = -1        ! LANDMASK variable on output grid
   integer :: dimid(3)                  ! dimension id's for variable to be defined
   integer :: startin(3) = (/1, 1, 1/)  ! starting position for variable read
   integer :: startout(3) = (/1, 1, 1/) ! starting position for variable written
   integer :: kountin(3) = (/-1, -1, 1/)! number of values for variable read
   integer :: kountout(3) = (/-1, -1, 1/)! number of values for variable written
   integer :: plonin = -1               ! longitude dimension (input grid)
   integer :: plonicein = -1            ! longitude dimension (input ice grid).  Compared to sst
   integer :: nlatin = -1               ! latitude dimension (input grid)
   integer :: nlaticein = -1            ! latitude dimension (input ice grid).  Compared to sst
   integer :: plonout = -1              ! longitude dimension (output grid)
   integer :: nlatout = -1              ! latitude dimension (output grid)
   integer :: i, j                      ! longitude, latitude indices
   integer :: ret                       ! return code
   integer :: nargs                     ! input arg
   integer :: n                         ! index loops thru input args
   integer :: itime                     ! integer equivalent of time variable
   integer :: year                      ! year from binary file
   integer :: bdate = -1                ! integer equivalent of begindate
   integer :: edate = huge(1)           ! integer equivalent of enddate

   character(len=8) :: begindate = ' '  ! date to start data values
   character(len=8) :: enddate   = ' '  ! date to start data values
!
! Data on input SST file
!
   integer, allocatable :: nlonin(:)      ! number of lons per lat on input grid

   real(r8) :: time                       ! time variable, copied from input to output
   real(r8) :: timeice                    ! time variable from ice file, checked vs. sst
   real(r8), allocatable :: lonin(:)      ! longitude on input grid (degrees)
   real(r8), allocatable :: lonin2d(:,:)  ! longitude on input grid (degrees)
   real(r8), allocatable :: latin(:)      ! latitude on input grid (degrees)
   real(r8), allocatable :: sstin(:,:)    ! mixed layer depths on input grid
   real(r8), allocatable :: icefracin(:,:)! ice fraction on input grid
!
! Data on output grid
!
   integer, allocatable :: nlonout(:)     ! number of lons per lat (reduced output grid)

   real(r8), allocatable :: latout(:)     ! latitude on output grid
   real(r8), allocatable :: lonout(:)     ! longitude on output grid
   real(r8), allocatable :: rlon(:,:)     ! longitude on (reduced) output grid
   real(r8), allocatable :: sstout(:,:)   ! mixed layer depths on output grid
   real(r8), allocatable :: orosst(:,:)   ! ocean/ice flag (1=ice 0=non-ice) based on sst
   real(r8), allocatable :: oroconc(:,:)  ! ocean/ice flag (1=ice 0=non-ice) based on ice conc.
   real(r8), allocatable :: landmasksst(:,:)  ! land flag (1=land 0=non-land) input grid
   real(r8), allocatable :: landmaskice(:,:)  ! land flag (1=land 0=non-land) input grid
   real(r8), allocatable :: landmaskout(:,:)  ! land flag (1=land 0=non-land) output grid
   real(r8), allocatable :: icefracout(:,:)   ! sea ice concentration on output grid
!
! Required by area averaging procedure
!
   integer :: mxovr = -1                  ! parameter needed by grid overlap routines
   integer :: ii,jj                       ! indices from input grid
   integer :: m                           ! index over number of overlapping points
   integer, allocatable :: iovr(:,:,:)    ! lon index of overlap input cell
   integer, allocatable :: jovr(:,:,:)    ! lat index of overlap input cell

   real(r8) :: wt                         ! weight
   real(r8) :: icewt                      ! weight for sea ice
   real(r8) :: waterwt                    ! weight for open ocean
   real(r8) :: totwt                      ! icewt + waterwt
   real(r8) :: dxin                       ! delta-x on input grid (degrees)
   real(r8) :: dxout                      ! delta-x on output grid (degrees)
   real(r8) :: dyin                       ! delta-y on input grid (degrees)
   real(r8) :: dyout                      ! delta-y on output grid (degrees)
   real(r8), allocatable :: wovr(:,:,:)   ! weight of overlap input cell
   real(r8), allocatable :: latinse(:)    ! latitude south edge on input grid (degrees)
   real(r8), allocatable :: latoutse(:)   ! latitude south edge on output grid (degrees)
   real(r8), allocatable :: loninwe(:,:)  ! longitude west edge on 2-d input grid (degrees)
   real(r8), allocatable :: lonoutwe(:,:) ! longitude west edge on 2-d output grid (degrees)
!
! For computing date info
!
   integer :: monlen(12)                  ! length of months
   data monlen/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
! Cmd line
!
   integer iargc
   external iargc

! parse command line arguments, saving them to be written to history attribute

   nargs = iargc ()
   n = 1
   cmdline = 'regrid '
   do while (n <= nargs)
      arg = ' '
      call getarg (n, arg)
      n = n + 1
      select case (arg)
      case ('-a')
         linear_interp = .false.
         usr_specified_interp = .true.
         cmdline = trim(cmdline) // ' -a'
      case ('-b')
         call getarg (n, arg)
         n = n + 1
         begindate = arg
         read (begindate,'(i8)') bdate
         cmdline = trim(cmdline) // ' -b ' // trim(begindate)
      case ('-e')
         call getarg (n, arg)
         n = n + 1
         enddate = arg
         read (enddate,'(i8)') edate
         cmdline = trim(cmdline) // ' -e ' // trim(enddate)
      case ('-g')
         call getarg (n, arg)
         n = n + 1
         gridout = arg
         cmdline = trim(cmdline) // ' -g ' // trim(gridout)
      case ('-i')
         call getarg (n, arg)
         n = n + 1
         icefilein = arg
         cmdline = trim(cmdline) // ' -i ' // trim(icefilein)
      case ('-l')
         linear_interp = .true.
         usr_specified_interp = .true.
         cmdline = trim(cmdline) // ' -l'
      case ('-o')
         call getarg (n, arg)
         n = n + 1
         fileout = arg
         cmdline = trim(cmdline) // ' -o ' // trim(fileout)
      case ('-s')
         call getarg (n, arg)
         n = n + 1
         sstfilein = arg
         cmdline = trim(cmdline) // ' -s ' // trim(sstfilein)
      case ('-v')
         verbose = .true.
         cmdline = trim(cmdline) // ' -v'
      case default
         write (6,*) 'Argument ', arg,' is not known'
         call usage_exit (' ')
      end select
   end do

   if (sstfilein == ' ') then
      call usage_exit ('No input SST file specified')
   end if

   if (icefilein == ' ') then
      call usage_exit ('No input ICE file specified')
   end if

   if (gridout == ' ') then
      call usage_exit ('No grid file specified')
   end if

   if (fileout == ' ') then
      call usage_exit ('No output file specified')
   end if
!
! Open 1x1 netcdf files for reading and get dimension info
!
   call wrap_nf_open (trim(sstfilein), nf_nowrite, ncidsstin)
   call wrap_nf_open (trim(icefilein), nf_nowrite, ncidicein)

   call wrap_nf_inq_dimid (ncidsstin, 'lon', londimidin)
   call wrap_nf_inq_dimlen (ncidsstin, londimidin, plonin)

   call wrap_nf_inq_dimid (ncidsstin, 'lat', latdimidin)
   call wrap_nf_inq_dimlen (ncidsstin, latdimidin, nlatin)

   call wrap_nf_inq_dimid (ncidsstin, 'time', timedimidin)
   call wrap_nf_inq_dimlen (ncidsstin, timedimidin, ntimein)
!
! Ensure ice file dims match those of sst
!
   call wrap_nf_inq_dimid (ncidicein, 'lon', londimidicein)
   call wrap_nf_inq_dimlen (ncidicein, londimidicein, plonicein)
   if (plonin /= plonicein) then
      write(6,*)'lon dim mismatch sst vs ice:',plonin,plonicein
      call err_exit ('regrid')
   end if

   call wrap_nf_inq_dimid (ncidicein, 'lat', latdimidicein)
   call wrap_nf_inq_dimlen (ncidicein, latdimidicein, nlaticein)
   if (nlatin /= nlaticein) then
      write(6,*)'lat dim mismatch sst vs ice:',nlatin,nlaticein
      call err_exit ('regrid')
   end if

   call wrap_nf_inq_dimid (ncidicein, 'time', timedimidicein)
   call wrap_nf_inq_dimlen (ncidicein, timedimidicein, ntimeicein)
   if (ntimein /= ntimeicein) then
      write(6,*)'time dim mismatch sst vs ice:',ntimein,ntimeicein
      call err_exit ('regrid')
   end if
!
! Allocate space for needed data from 1x1 files
!
   allocate (nlonin(nlatin))
   allocate (lonin(plonin))
   allocate (latin(nlatin))
   allocate (sstin(plonin,nlatin))
   allocate (icefracin(plonin,nlatin))
   allocate (landmasksst(plonin,nlatin))
   allocate (landmaskice(plonin,nlatin))

   nlonin(:) = plonin                ! Input grid is rectangular
!
! Read needed grid info and variables from 1x1 files
!
   call wrap_nf_inq_varid (ncidsstin, 'lon', lonidin)
   call wrap_nf_inq_varid (ncidsstin, 'lat', latidin)
   call wrap_nf_inq_varid (ncidsstin, 'time', timeidin)
   call wrap_nf_inq_varid (ncidicein, 'time', timeidicein)  ! to compare to sst
   call wrap_nf_inq_varid (ncidsstin, 'SST', sstidin)
   call wrap_nf_inq_varid (ncidicein, 'SEAICE', iceidin)

   ! Determine format of time variable in 1x1 files.
   ! Old format was 'YYYYMMDD', new format is "days since ..."
   call wrap_nf_get_att_text(ncidsstin, timeidin, 'units', units_att)
   if (units_att(1:10) == 'days since') time_coord = 'new'

   ! If 1x1 files use new time units, then get the date info from the date variable
   if (time_coord == 'new') then
      call wrap_nf_inq_varid (ncidsstin, 'date', dateidin)
   end if

   call wrap_nf_get_var_double (ncidsstin, lonidin, lonin)
   call wrap_nf_get_var_double (ncidsstin, latidin, latin)
!
! Open grid file and get dimension info
!
   call wrap_nf_open (trim(gridout), NF_NOWRITE, ncidgrd)

   call wrap_nf_inq_dimid (ncidgrd, 'lon', londimidgrd)
   call wrap_nf_inq_dimlen (ncidgrd, londimidgrd, plonout)

   call wrap_nf_inq_dimid (ncidgrd, 'lat', latdimidgrd)
   call wrap_nf_inq_dimlen (ncidgrd, latdimidgrd, nlatout)
!
! Allocate space for needed input and output data on output file
!
   allocate (latout(nlatout))
   allocate (lonout(plonout))
   allocate (rlon(plonout,nlatout))
   allocate (nlonout(nlatout))
   allocate (sstout(plonout,nlatout))
   allocate (landmaskout(plonout,nlatout))
   allocate (orosst(plonout,nlatout))
   allocate (oroconc(plonout,nlatout))
   allocate (icefracout(plonout,nlatout))
!
! Get coordinate variable IDs from grid file
!
   call wrap_nf_inq_varid (ncidgrd, 'lon', lonidgrd)
   call wrap_nf_inq_varid (ncidgrd, 'lat', latidgrd)
!
! Create output file
!
   if (nf_create(trim(fileout), NF_CLOBBER, ncidout) /= NF_NOERR) then
      write(6,*)'FAILURE: creating output file: ', fileout
      stop 999
   end if
!
! Define dimensions on output file
!
   call wrap_nf_def_dim (ncidout, 'lat',  nlatout, latdimidout)
   call wrap_nf_def_var (ncidout, 'lat',  NF_DOUBLE, 1, latdimidout, latidout)
   call wrap_nf_def_dim (ncidout, 'lon',  plonout, londimidout)
   call wrap_nf_def_var (ncidout, 'lon',  NF_DOUBLE, 1, londimidout, lonidout)
   call wrap_nf_def_dim (ncidout, 'time', NF_UNLIMITED, timedimidout)
   call wrap_nf_def_var (ncidout, 'time', NF_DOUBLE, 1, timedimidout, timeidout)
   call wrap_nf_def_var (ncidout, 'date', NF_INT, 1, timedimidout, dateidout)
   call wrap_nf_def_var (ncidout, 'datesec', NF_INT, 1, timedimidout, datesecidout)
!
! Get fillvalues for ice and sst
!
   call wrap_nf_get_att_double (ncidsstin, sstidin, '_FillValue', fillvaluesst)
   call wrap_nf_get_att_double (ncidicein, iceidin, '_FillValue', fillvalueice)
!
! Define history attribute.
!
   call get_curr_timestamp(cur_timestamp)
   history = trim(cur_timestamp) // ' ' // trim(cmdline)
   cmdlen = len_trim (history)
   call wrap_nf_put_att_text (ncidout, NF_GLOBAL, 'history', cmdlen, trim(history))
!
! Define SST on the output file
!
   dimid(1) = londimidout
   dimid(2) = latdimidout
   dimid(3) = timedimidout

   call wrap_nf_def_var (ncidout, 'SST', NF_FLOAT, 3, dimid, sstidout) 
   call wrap_nf_def_var (ncidout, 'ICEFRAC', NF_FLOAT, 3, dimid, iceidout) 
   call wrap_nf_def_var (ncidout, 'OROSST', NF_FLOAT, 3, dimid, oroidsst) 
   call wrap_nf_def_var (ncidout, 'OROCONC', NF_FLOAT, 3, dimid, oroidconc) 
!
! Copy and define attributes for SST and ICE
!
   call wrap_nf_copy_att (ncidgrd, latidgrd, 'long_name', ncidout, latidout)
   call wrap_nf_copy_att (ncidgrd, latidgrd, 'units',     ncidout, latidout)

   call wrap_nf_copy_att (ncidgrd, lonidgrd, 'long_name', ncidout, lonidout)
   call wrap_nf_copy_att (ncidgrd, lonidgrd, 'units',     ncidout, lonidout)

   if (time_coord == 'old') then
      call wrap_nf_copy_att (ncidsstin, timeidin, 'long_name', ncidout, timeidout)
   else
      call wrap_nf_copy_att (ncidsstin, timeidin, 'calendar', ncidout, timeidout)
   end if
   call wrap_nf_copy_att (ncidsstin, timeidin, 'units',     ncidout, timeidout)

   call wrap_nf_copy_att (ncidsstin, sstidin, 'long_name', ncidout, sstidout)
   call wrap_nf_copy_att (ncidsstin, sstidin, 'units',     ncidout, sstidout)
!   call wrap_nf_put_att_real (ncidout, sstidout, '_FillValue', NF_FLOAT, 1, fillvalueout)

   call wrap_nf_copy_att (ncidicein, iceidin, 'long_name', ncidout, iceidout)
   call wrap_nf_put_att_text (ncidout, iceidout, 'units', len('1'), '1')
!   call wrap_nf_put_att_real (ncidout, iceidout, '_FillValue', NF_FLOAT, 1, fillvalueout)

!
! Set longitude values for output grid (check for reduced grid)
!
   ret = nf_inq_varid (ncidgrd, 'nlon', nlonidgrd)
   if (ret == nf_noerr) then                            ! probably reduced grid
      ret = nf_inq_varid (ncidgrd, 'rlon', rlonidgrd)
      if (ret == nf_noerr) then

         is_reduced_grid = .true.

         call wrap_nf_get_var_int (ncidgrd, nlonidgrd, nlonout)
         call wrap_nf_get_var_double (ncidgrd, rlonidgrd, rlon)
         
         ! Put reduced grid variables in output file
         call wrap_nf_def_var (ncidout, 'nlon', NF_INT, 1, dimid, nlonidout) 
         call wrap_nf_def_var (ncidout, 'rlon', NF_DOUBLE, 2, dimid, rlonidout) 
         call wrap_nf_copy_att (ncidgrd, lonidgrd, 'long_name', ncidout, rlonidout)
         call wrap_nf_copy_att (ncidgrd, lonidgrd, 'units',     ncidout, rlonidout)

      else
         write(6,*)gridout, ' has nlon but not rlon => is screwed up'
         stop 999
      end if
   else                                                 ! full grid: define nlon and rlon
      do j=1,nlatout
         nlonout(j) = plonout
         do i=1,plonout
            rlon(i,j) = (i-1)*360./plonout
         end do
      end do
   end if

   ! Read latitude and longitude values for output grid
   call wrap_nf_get_var_double (ncidgrd, latidgrd, latout)
   call wrap_nf_get_var_double (ncidgrd, lonidgrd, lonout)

!
! Determine interpolation type if not specified on cmd line
! Use linear interpolation if output grid is finer than input grid.
! Otherwise use area averaging
!
   if (.not. usr_specified_interp) then
      if (nlonout(plonout/2) > nlonin(plonin/2)) then
         linear_interp = .true.
      else
         linear_interp = .false.
      end if
   end if
!
   if (linear_interp) then
      write(6,*)'Using LINEAR interpolation to output grid'
      call wrap_nf_put_att_text (ncidout, NF_GLOBAL, 'interp_type', 6, 'LINEAR')
      allocate (lonin2d(plonin,nlatin))
!
! ASSUME input grid is regular (not reduced)
!
      do j=1,nlatin
         lonin2d(:,j) = lonin(:)
      end do
   else
!
! For area averaging, compute overlap quantities.  Only output lats need to be computed
! from lat vals. on input .nc file because they might not be uniformly spaced (e.g. Gaussian)
!
      write(6,*)'Using AREA AVERAGING to output grid'
      call wrap_nf_put_att_text (ncidout, NF_GLOBAL, 'interp_type', 14, 'AREA AVERAGING')

      allocate (loninwe(plonin+1,nlatin))
      allocate (lonoutwe(plonout+1,nlatout))
      allocate (latinse(nlatin+1))
      allocate (latoutse(nlatout+1))

      dxin = 360._r8 / plonin
      dyin = 180._r8 / nlatin
!
! West edge of input grid: starts at Greenwich
! South edge of input grid: starts at south pole
!
      do i=1,plonin+1
         loninwe(i,:) = (i-1)*dxin
      end do

      do j=1,nlatin+1
         latinse(j) = -90. + (j-1)*dyin
      end do
!
! West edge of output grid: starts west of Greenwich
!
      do j=1,nlatout
         dxout = 360._r8 / nlonout(j)
         do i=1,nlonout(j)+1
            lonoutwe(i,j) = -0.5*dxout + (i-1)*dxout
         end do
      end do
      
      latoutse(1) = -90.
      do j=2,nlatout
         dyout = latout(j) - latout(j-1)
         latoutse(j) = latout(j) - 0.5*dyout
      end do
      latoutse(nlatout+1) = +90.
   
      call max_ovr (plonin, nlatin, nlonin, plonout, nlatout, nlonout, &
                    loninwe, latinse, lonoutwe, latoutse, mxovr)

      allocate (iovr(plonout,nlatout,mxovr))
      allocate (jovr(plonout,nlatout,mxovr))
      allocate (wovr(plonout,nlatout,mxovr))
      
      call map_i (plonin, nlatin, nlonin, loninwe, latinse, &
                  plonout, nlatout, nlonout, lonoutwe, latoutse, &
                  mxovr, iovr, jovr, wovr)
   end if
!
! Done with data definition
!
   if (nf_enddef (ncidout) /= nf_noerr) then
      write(6,*)'Ending define mode failed on ', fileout
      stop 999
   end if

   ! Add coordinate data to output file
   call wrap_nf_put_var_double (ncidout, latidout, latout)
   call wrap_nf_put_var_double (ncidout, lonidout, lonout)
   if (is_reduced_grid) then
      call wrap_nf_put_var_int (ncidout, nlonidout, nlonout)
      call wrap_nf_put_var_double (ncidout, rlonidout, rlon)
   end if

!
! Define lon and lat count variables for netcdf calls 
!
   kountin(1) = plonin
   kountin(2) = nlatin

   kountout(1) = plonout
   kountout(2) = nlatout
!
! Get LANDMASK variables for output grid if requested
!
   if (uselandmask) then
      call wrap_nf_inq_varid (ncidgrd, 'LANDMASK', landmaskidout)
      call wrap_nf_get_var_double (ncidgrd, landmaskidout, landmaskout)
   else
      landmaskout(:,:) = 0.
   end if
!
! Loop over time, binning 1x1 SST data to output grid
!
   do n=1,ntimein
      call wrap_nf_get_vara_double (ncidsstin, timeidin, n, 1, time)
      call wrap_nf_get_vara_double (ncidicein, timeidicein, n, 1, timeice)
!
! Ensure times must match sst vs ice
!
      if (time /= timeice) then
         write(6,*)'time mismatch sst vs ice:',time,timeice
         call err_exit ('regrid')
      end if
      if (time_coord == 'old') then
         ! date info from time variable
         itime = nint (time)
      else
         ! date info from date variable
         call wrap_nf_get_vara_int (ncidsstin, dateidin, n, 1, itime)
      end if
      if (bdate <= itime .and. edate >= itime) then
         if (verbose) then
            write(6,*)'Processing date ', itime
         end if
         startin(3) = n
         call wrap_nf_get_vara_double (ncidsstin, sstidin, startin, kountin, sstin)
         call wrap_nf_get_vara_double (ncidicein, iceidin, startin, kountin, icefracin)
!
! Turn percentage into fraction
!
         where (icefracin(:,:) /= fillvalueice)
            icefracin(:,:) = icefracin(:,:) * 0.01
         end where

         if (linear_interp) then
            call interp_driver (plonin, nlatin, 1, nlonin, lonin2d, latin, 0._r8, icefracin, &
                                plonout, nlatout, 1, nlonout, rlon, latout, 0._r8, icefracout)
            call interp_driver (plonin, nlatin, 1, nlonin, lonin2d, latin, 0._r8, sstin, &
                                plonout, nlatout, 1, nlonout, rlon, latout, 0._r8, sstout)
            do j=1,nlatout
               do i=1,nlonout(j)
                  if (nint (landmaskout(i,j)) == 1) then
                     sstout(i,j) = fillvalueout
                     icefracout(i,j) = fillvalueout
                  end if
               end do
            end do
         else
            do j=1,nlatout
               do i=1,nlonout(j)
                  if (nint (landmaskout(i,j)) == 1) then
                     sstout(i,j) = fillvalueout
                     icefracout(i,j) = fillvalueout
                  else
                     sstout(i,j) = 0.
                     icefracout(i,j) = 0.
                     icewt = 0.
                     waterwt = 0.
                     
                     do m=1,mxovr               ! overlap cell index
                        ii = iovr(i,j,m)        ! lon index (input grid) of overlap cell
                        jj = jovr(i,j,m)        ! lat index (input grid) of overlap cell
                        wt = wovr(i,j,m)        ! overlap weight
                        sstout(i,j) = sstout(i,j) + sstin(ii,jj)*wt
                        icefracout(i,j) = icefracout(i,j) + icefracin(ii,jj)*wt
!
! The next bit of code is only useful when using SST to define where ice is
!
                        if (sstin(ii,jj) > iceflag + 0.001) then
                           waterwt = waterwt + wovr(i,j,m)
                        else
                           icewt = icewt + wovr(i,j,m)
                        end if
                     end do
                  
                     totwt = waterwt + icewt
                     if (totwt < 0.9999 .or. totwt > 1.0001) then
                        write(6,*) 'regrid: bad totwt=', totwt
                        stop 999
                     end if
!
! When computing ice flag using SST's, ensure ice point if fraction is > 50%
!
!                  if (icewt >= 0.5 .and. sstout(i,j) > iceflag) then
!                     sstout(i,j) = iceflag
!                  end if
                  end if
               end do      ! i=1,nlonout(j)
            end do         ! j=1,nlatout
         end if
!
! Open output file for scatter plot at start of every year
!
         if (scatterplots .and. mod(itime/100,100) == 1) then
            year = itime/10000
            write (scattername,200) year
200         format ('scatter.data.',i4)
            open (unit=1,file=scattername,form='formatted',action='WRITE')
         end if
!
! Produce data for scatter plot, randomized by rectangular distributin bounded by
! +-0.5x1.E-02
!
         do j=1,nlatin
            do i=1,nlonin(j)
               if (icefracin(i,j) /= fillvalueice .and. sstin(i,j) /= fillvaluesst) then
                  if (verbose .and. icefracin(i,j) > 0.5 .and. sstin(i,j) > 6.) then
                     write(6,*)'i=',i,' j=',j,' icefrac,sst=',icefracin(i,j),sstin(i,j)
                  end if
                  if (scatterplots .and. icefracin(i,j) > 0.) then
                     call random_number (pertval)
                     write(1,'(2f8.4)') icefracin(i,j) + (0.5-pertval)*.01, sstin(i,j)
                  end if
               end if
            end do
         end do
         
         if (scatterplots .and. mod(itime,100) == 12) then
            close (unit=1)
         end if
!
! Define sea ice flags from 1) ice concentration, and 2) SST dataset
!
         do j=1,nlatout
            do i=1,nlonout(j)
               orosst(i,j) = 0.
               oroconc(i,j) = 0.
               if (sstout(i,j) <= iceflag) then
                  orosst(i,j) = 1.
               end if
               if (icefracout(i,j) >= 0.5) then  ! sea ice
                  oroconc(i,j) = 1.
               end if
            end do
         end do
!
! Write out SST data
!
         if (verbose) then
            write(6,*)'regrid: writing SST data for time slice ', n
         end if
!
! Build date info: use only year/month from input file because day info is incompatible 
! with CAM
!
         date(1) = (itime / 100) * 100
         mo = mod ((itime / 100), 100)
         if (mo < 1 .or. mo > 12) then
            write(6,*)'regrid: Bad mo index=', mo
            stop 999
         end if
         date(1) = date(1) + monlen(mo)/2 + 1
         datesec(1) = 0
         if (mod(monlen(mo), 2) /= 0) datesec(1) = 43200

         startout(3) = n
         call wrap_nf_put_vara_double (ncidout, timeidout, n, 1, time)
         call wrap_nf_put_vara_int (ncidout, dateidout, n, 1, date)
         call wrap_nf_put_vara_int (ncidout, datesecidout, n, 1, datesec)
         call wrap_nf_put_vara_double (ncidout, sstidout, startout, kountout, sstout)
         call wrap_nf_put_vara_double (ncidout, iceidout, startout, kountout, icefracout)
         call wrap_nf_put_vara_double (ncidout, oroidsst, startout, kountout, orosst)
         call wrap_nf_put_vara_double (ncidout, oroidconc, startout, kountout, oroconc)
      else
         if (verbose) then
            write(6,*)'Skipping date ', itime
         end if
      end if
   end do
!
! Clean up
!   
   call wrap_nf_close (ncidsstin)
   call wrap_nf_close (ncidicein)
   call wrap_nf_close (ncidgrd)
   call wrap_nf_close (ncidout)

   if (scatterplots) then
      close (unit=1)
   end if

   stop
end program regrid

subroutine usage_exit (arg)
   implicit none
   character*(*) arg
   
   if (arg /= ' ') write (6,*) arg
   write (6,*) 'Usage: regrid [-a] [-b begin] [-e end] -g gridfile -i icefile [-l] -o outfile -s sstfile [-v]'
   write (6,*) '  -a: Use area averaging (not linear interpolation). default if fine->coarse'
   write (6,*) '  -b date: (yyyymmdd) skip input dates earlier than this. default: start at the beginning'
   write (6,*) '  -e date: (yyyymmdd) skip input dates later than this. default: end at the end'
   write (6,*) '  -g file: input netcdf file containing grid for output fields'
   write (6,*) '  -i file: 1x1 input netcdf ICE concentration file'
   write (6,*) '  -l: Use linear interpolation (not area averaging). default if coarse->fine'
   write (6,*) '  -o file: output netcdf file with regridded SST and ICE'
   write (6,*) '  -s file: 1x1 input netcdf SST file'
   write (6,*) '  -v: verbose mode'
   stop 999
end subroutine usage_exit

subroutine get_curr_timestamp(time)
! return timestamp formatted as "YYYY-MM-DD HH:MM:SS"   
   character(len=19), intent(out) :: time
   integer :: t(8)
   call date_and_time(values=t)
   write(time,'(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') t(1),'-',t(2),'-',t(3),' ',&
                                                         t(5),':',t(6),':',t(7)
end subroutine get_curr_timestamp
