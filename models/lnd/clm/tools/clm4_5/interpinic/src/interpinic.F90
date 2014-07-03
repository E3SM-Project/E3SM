module interpinic

  !----------------------------------------------------------------------- 
  ! Interpolate initial conditions file from one resolution and/or landmask
  ! to another resolution and/or landmask
  !----------------------------------------------------------------------- 

  use shr_kind_mod   , only: r8 => shr_kind_r8
  use shr_const_mod  , only: SHR_CONST_PI, SHR_CONST_REARTH
  use shr_sys_mod    , only: shr_sys_flush
  implicit none

  private

  ! Public methods

  public :: interp_filei

  ! Public data

  logical, public  :: override_missing = .true. ! if you want to override missing
                                                ! types with closest bare-soil
                                                ! Otherwise it will abort.

  ! Private methods

  private :: interp_ml_real
  private :: interp_sl_real
  private :: interp_sl_int
  private :: findMinDistPFTs
  private :: findMinDistCols
  private :: findMinDistLDUs

  ! Private data
 
  ! Variables read in from input and output initial files

  integer :: nlevsno        ! maximum number of snow levels
  integer :: nlevsno1       ! maximum number of snow levels plus one
  integer :: nlevlak        ! number of lake levels
  integer :: nlevtot        ! number of soil and snow levels
  integer :: nlevtot_o       ! number of soil and snow levels on output file
  integer :: nlevgrnd        ! number of soil levels
  integer :: nlevgrnd_o      ! number of soil levels on output file
  integer :: nlevcan        ! number of canopy layers

  integer :: numpfts        ! input file number of pfts 
  integer :: numpftso       ! output file number of pfts 
  integer :: numcols        ! input file number of columns 
  integer :: numcolso       ! output file number of columns 
  integer :: numldus        ! input file number of landunits
  integer :: numlduso       ! output file number of landunits

  ! RTM river routing model
  integer, parameter :: rtmlon = 720  ! # of rtm longitudes
  integer, parameter :: rtmlat = 360  ! # of rtm latitudes
  real(r8) :: volr(rtmlon,rtmlat)     ! water volume in cell (m^3)

  ! Other parameter sizes
  integer, parameter :: numrad = 2    ! # of radiation bands

  ! Parameters

  real(r8), parameter :: spval = 1.e36_r8 ! special value for missing data (ocean)
  real(r8), parameter :: deg2rad = SHR_CONST_PI/180._r8
  real(r8), parameter :: re = SHR_CONST_REARTH
  ! These types need to agree with the types in clm_varcon.F90 in the main CLM model code
  integer,  parameter :: croptype     = 15
  integer,  parameter :: istcrop      = 2
  integer,  parameter :: istsoil      = 1
  integer,  parameter :: baresoil     = 0
  integer,  parameter :: nonurbcol    = 1

  logical , save :: allPFTSfromSameGC = .false. ! Get all PFTS from the same gridcells
  logical , save :: noAbortIfDNE      = .false. ! Do NOT abort if some input data does not exist
  integer , allocatable, save :: colindx(:)     ! Column mapping indices
  integer , allocatable, save :: pftindx(:)     ! PFT mapping indices
  integer , allocatable, save :: lduindx(:)     ! land-unit mapping indices

  logical, save :: vertinterp = .false.

  SAVE

contains

  !=======================================================================

  subroutine interp_filei (fin, fout, cmdline)

    !----------------------------------------------------------------------- 
    ! Read initial data from netCDF instantaneous initial data history file 
    !-----------------------------------------------------------------------

    use netcdf
    use shr_infnan_mod , only: shr_infnan_isnan

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments------ -----------------------------
    character(len=256), intent(in) :: fin   !input  initial dataset
    character(len=256), intent(in) :: fout   !output initial dataset
    character(len=256), intent(in) :: cmdline    !command line arguments
    ! --------------------------------------------------------------------

    ! ------------------------ local variables -----------------------------
    integer :: i,j,k,l,m,n         ! loop indices    
    integer :: ncidi               ! netCDF dataset id
    integer :: nvecin              ! input vector length
    integer :: nvars               ! number of variables
    integer :: nvecout             ! output vector length
    integer :: ncido               ! output net
    integer :: dimid               ! netCDF dimension id 
    integer :: dimidpft            ! netCDF dimension id PFT
    integer :: dimidldu            ! netCDF dimension id PFT
    integer :: dimidcols           ! netCDF dimension id columns
    integer :: dimidlak            ! netCDF dimension id lake
    integer :: dimidsno            ! netCDF dimension id snow depth
    integer :: dimidsno1           ! netCDF dimension id snow depth plus one
    integer :: dimidtot            ! netCDF dimension id total
    integer :: dimidgrnd            ! netCDF dimension id ground
    integer :: dimidrad            ! netCDF dimension id numrad
    integer :: dimidcan            ! netCDF dimension id levcan
    integer :: dimidrtmlat         ! netCDF dimension id rtmlat
    integer :: dimidrtmlon         ! netCDF dimension id rtmlon
    integer :: varid               ! netCDF variable id
    integer :: varido              ! netCDF variable id
    integer :: xtype               ! netCDF variable type
    integer :: period              ! period accumulator scalar to copy over
    integer :: ndims               ! netCDF number of dimensions
    integer :: dimids(3) = -1      ! netCDF dimension ids
    integer :: dimlen              ! input dimension length       
    integer :: ret                 ! netcdf return code
    integer :: ncformat            ! netcdf file format
    character(len=256) :: varname  !variable name
    real(r8), allocatable :: rbufmlo (:,:) !output array
    real(r8), allocatable :: rbufmco (:,:) !output array
    !--------------------------------------------------------------------

    write (6,*) 'Mapping clm initial data from input to output initial files'

    ! Open input and output initial conditions files (both just for reading now)

    call check_ret (nf90_open(fin,  NF90_NOWRITE, ncidi ))
    call check_ret (nf90_open(fout, NF90_NOWRITE, ncido ))
    call check_ret (nf_inq_format( ncido, ncformat ))
    if ( ncformat /= NF_FORMAT_64BIT )then
       write (6,*) 'error: output file is NOT in NetCDF large-file format!'
       stop
    end if

    call check_ret (nf90_inq_dimid(ncidi, "column", dimidcols ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidcols, len=numcols))
    call check_ret (nf90_inq_dimid(ncido, "column", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=numcolso))
    write (6,*) 'input numcols = ',numcols,' output numcols = ',numcolso

    call check_ret (nf90_inq_dimid(ncidi, "pft", dimidpft ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidpft, len=numpfts))
    call check_ret (nf90_inq_dimid(ncido, "pft", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=numpftso))
    write (6,*) 'input numpfts = ',numpfts,' output numpfts = ',numpftso


    call check_ret (nf90_inq_dimid(ncidi, "landunit", dimidldu ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidldu, len=numldus))
    call check_ret (nf90_inq_dimid(ncido, "landunit", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=numlduso))
    write (6,*) 'input numldus = ',numldus,' output numldus = ',numlduso

    call check_ret (nf90_inq_dimid(ncidi, "levsno", dimidsno ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidsno, len=nlevsno))
    call check_ret (nf90_inq_dimid(ncido, "levsno", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=dimlen))
    if (dimlen/=nlevsno) then
       write (6,*) 'error: input and output nlevsno values disagree'
       write (6,*) 'input nlevsno = ',nlevsno,' output nlevsno = ',dimlen
       stop
    end if

    ret = nf90_inq_dimid(ncidi, "levsno1", dimidsno1)
    if (ret == NF90_NOERR) then
       call check_ret (nf90_inquire_dimension(ncidi, dimidsno1, len=nlevsno1))
       call check_ret (nf90_inq_dimid(ncido, "levsno1", dimid ))
       call check_ret (nf90_inquire_dimension(ncido, dimid, len=dimlen))
       if (dimlen/=nlevsno1) then
          write (6,*) 'error: input and output nlevsno1 values disagree'
          write (6,*) 'input nlevsno1 = ',nlevsno1,' output nlevsno1 = ',dimlen
          stop
       end if
    else
       write (6,*) 'levsno1 dimension does NOT exist on the input dataset'
       dimidsno1 = -9999
       nlevsno1  = -9999
    end if

    ret = nf90_inq_dimid(ncidi, "levcan", dimidcan)
    if (ret == NF90_NOERR) then
       call check_ret (nf90_inquire_dimension(ncidi, dimidcan, len=nlevcan))
       call check_ret (nf90_inq_dimid(ncido, "levcan", dimid ))
       call check_ret (nf90_inquire_dimension(ncido, dimid, len=dimlen))
       if (dimlen/=nlevcan) then
          write (6,*) 'error: input and output nlevcan values disagree'
          write (6,*) 'input nlevcan = ',nlevcan,' output nlevcan = ',dimlen
          stop
       end if
    else
       write (6,*) 'levcan dimension does NOT exist on the input dataset'
       dimidcan = -9999
       nlevcan  = -9999
    end if

    call check_ret (nf90_inq_dimid(ncidi, "levlak", dimidlak ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidlak, len=nlevlak))
    call check_ret (nf90_inq_dimid(ncido, "levlak", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=dimlen))
    if (dimlen/=nlevlak) then
       write (6,*) 'error: input and output nlevlak values disagree'
       write (6,*) 'input nlevlak = ',nlevlak,' output nlevlak = ',dimlen
       stop
    end if

    call check_ret (nf90_inq_dimid(ncidi, "levtot", dimidtot ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidtot, len=nlevtot))
    call check_ret (nf90_inq_dimid(ncido, "levtot", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=nlevtot_o))
    if (nlevtot_o /= nlevtot) then
       write (6,*) 'input and output nlevtot values disagree'
       write (6,*) 'input nlevtot = ',nlevtot,' output nlevtot = ',nlevtot_o
       vertinterp = .true.
    end if

    call check_ret (nf90_inq_dimid(ncidi, "levgrnd", dimidgrnd ))
    call check_ret (nf90_inquire_dimension(ncidi, dimidgrnd, len=nlevgrnd))
    call check_ret (nf90_inq_dimid(ncido, "levgrnd", dimid ))
    call check_ret (nf90_inquire_dimension(ncido, dimid, len=nlevgrnd_o))
    if (nlevgrnd_o /= nlevgrnd) then
       write (6,*) 'input and output nlevgrnd values disagree'
       write (6,*) 'input nlevgrnd = ',nlevgrnd,' output nlevgrnd = ',nlevgrnd_o
       vertinterp = .true.
    end if
 
    ! numrad dimension
    ret = nf90_inq_dimid(ncidi, "numrad", dimidrad)
    if (ret/=NF90_NOERR) call handle_error (ret)
    ret = nf90_inquire_dimension(ncidi, dimidrad, len=dimlen)
    if (dimlen/=numrad) then
       write (6,*) 'error: input numrad dimension size does not equal ',numrad; stop
    end if
    allocate( rbufmlo(numrad,numpftso) )
    allocate( rbufmco(nlevcan,numpftso) )

    ! If RTM data exists on input file
    ret = nf_inq_varid (ncidi, 'RTM_VOLR_LIQ', varid)
    if (ret == NF_NOERR) then
       call check_ret (nf90_inq_dimid(ncidi, "rtmlon", dimidrtmlon))
       call check_ret (nf90_inquire_dimension(ncidi, dimidrtmlon, len=dimlen))
       if (dimlen/=rtmlon) then
          write (6,*) 'error: input rtmlon does not equal ',rtmlon; stop
       end if
       call check_ret (nf90_inq_dimid(ncidi, "rtmlat", dimidrtmlat))
       call check_ret (nf90_inquire_dimension(ncidi, dimidrtmlat, len=dimlen))
       if (dimlen/=rtmlat) then
          write (6,*) 'error: input rtmlat does not equal ',rtmlat; stop
       end if
    else
       dimidrtmlat = -1
       dimidrtmlon = -1
    end if
    !
    ! Check if DGVM data exists on the dataset and if so -- use all PFT data from same grid cell
    ! Otherwise use data for the same veg type from potentially different grid-cells
    !
    ret = nf90_inq_varid(ncidi, 'present', varid)
    if (ret==NF90_ENOTVAR)then
       allPFTSfromSameGC = .false.
    else if (ret/=NF90_NOERR)then
       call handle_error (ret)
    else
       allPFTSfromSameGC = .true.
    end if

    ! For each output pft, find the input pft, pftindx, that is closest

    write(6,*)'finding minimum distance for pfts'
    allocate(pftindx(numpftso))
    call findMinDistPFTs( ncidi, ncido, pftindx )

    ! For each output column, find the input column, colindx, that is closest

    write(6,*)'finding minimum distance for columns'
    allocate(colindx(numcolso))
    call findMinDistCols( ncidi, ncido, colindx )

    ! For each output landunit, find the input landunit, lduindx, that is closest

    write(6,*)'finding minimum distance for landunits'
    allocate(lduindx(numlduso))
    call findMinDistLDUs( ncidi, ncido, lduindx )
    
    ! Get list of variables
    call check_ret (nf90_inquire(ncidi, nVariables=nvars ))
    !
    ! OK now, open the output file for writing
    !
    call check_ret(nf90_close( ncido))
    call check_ret (nf90_open(fout, ior(NF90_WRITE,  NF_64BIT_OFFSET), ncido ))

    call addglobal (ncido, cmdline)

    ! Read input initial data and write output initial data
    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    write(6,*)'reading in initial dataset'
    do i = 1, nvars
       varid = i
       call check_ret (nf_inq_varname(ncidi, varid, varname ))

       ! Skip names that do NOT end in _PERIOD or _VALUE and match
       ! specific list of names
       if(((index(varname,"_PERIOD"            ) == 0 ) .and.  &
           (index(varname,"_VALUE"             ) == 0 )).and. (&
           (index(varname,"timemgr_"           ) == 1 ) .or. &
           (index(varname,"PFT_"               ) == 1 ) .or. &
           (index(varname,"grid1d_"            ) == 1 ) .or. &
           (index(varname,"cols1d_"            ) == 1 ) .or. &
           (index(varname,"pfts1d_"            ) == 1 ) .or. &
           (index(varname,"land1d_"            ) == 1 ) .or. &
           (index(varname,"type1d_"            ) == 1 ) .or. &
           (index(varname,"EFLX_LWRAD_OUT"     ) == 1 ) .or. &
           (index(varname,"FRAC_VEG_NOSNO_ALB" ) == 1 ) .or. &
           (index(varname,"tlai"               ) == 1 ) .or. &
           (index(varname,"tsai"               ) == 1 ) .or. &
           (index(varname,"elai"               ) == 1 ) .or. &
           (index(varname,"esai"               ) == 1 ) .or. &
           (index(varname,"T_REF"              ) == 1 ) .or. &
           (index(varname,"TREF"               ) == 1 ) .or. &
           (index(varname,"RTM_INPUT_LIQ"      ) == 1 ) .or. &
           (index(varname,"RTM_INPUT_ICE"      ) == 1 ) .or. &
           (index(varname,"t_ref2m"            ) == 1 ) .or. &
           (index(varname,"vf_sr"              ) == 1 ) .or. &
           (index(varname,"vf_wr"              ) == 1 ) .or. &
           (index(varname,"vf_sw"              ) == 1 ) .or. &
           (index(varname,"vf_rw"              ) == 1 ) .or. &
           (index(varname,"vf_ww"              ) == 1 ) .or. &
           (index(varname,"sabs_roof_dir"      ) == 1 ) .or. &
           (index(varname,"sabs_roof_dif"      ) == 1 ) .or. &
           (index(varname,"sabs_sunwall_dir"   ) == 1 ) .or. &
           (index(varname,"sabs_sunwall_dif"   ) == 1 ) .or. &
           (index(varname,"sabs_shadewall_dir" ) == 1 ) .or. &
           (index(varname,"sabs_shadewall_dif" ) == 1 ) .or. &
           (index(varname,"sabs_improad_dir"   ) == 1 ) .or. &
           (index(varname,"sabs_improad_dif"   ) == 1 ) .or. &
           (index(varname,"sabs_perroad_dir"   ) == 1 ) .or. &
           (index(varname,"sabs_perroad_dif"   ) == 1 ) .or. &
           (index(varname,"fpcgridold"         ) == 1 ) .or. &
           (index(varname,"htop"               ) == 1 ) .or. &
           (index(varname,"hbot"               ) == 1 ) .or. &
           (index(varname,"locfnh"             ) == 1 ) .or. &
           (index(varname,"locfnhr"            ) == 1 ) ) )then
          write (6,*) 'Skipping variable: ', trim(varname)
          cycle
       end if

       ret = nf90_inq_varid(ncido, varname, varido)
       if (ret==NF90_ENOTVAR)then
          write (6,*) 'Variable NOT found on output file: ', trim(varname)
          cycle
       else if (ret/=NF90_NOERR)then
          call handle_error (ret)
       end if
       call check_ret (nf90_inquire_variable( ncidi, varid, xtype=xtype, ndims=ndims, &
            dimids=dimids ))

       ! For scalar variables
       if ( ndims == 0 ) then
          if ( index(varname,"_PERIOD" ) /= 0 .or. trim(varname) == "RTM_NCOUNT" )then
             call check_ret(nf90_inq_varid(ncidi, varname, varid))
             call check_ret(nf90_get_var(ncidi, varid, period))

             ret = nf90_inq_varid(ncido, varname, varid)
             if (ret==NF90_ENOTVAR)then
                write (6,*) 'Variable NOT found on output file: ', trim(varname)
                cycle
             else if (ret/=NF90_NOERR)then
                call handle_error (ret)
             end if
             write (6,*) 'PERIOD variable copied over: ', trim(varname)
             call check_ret(nf90_put_var(ncido, varid, period))
          else
             write (6,*) 'Skipping scalar variable: ', trim(varname)
          end if
       ! For 1D variables
       else if ( ndims == 1 ) then

          if ( dimids(1) == dimidcols )then
             nvecin  = numcols
             nvecout = numcolso
          else if ( dimids(1) == dimidpft )then
             nvecin  = numpfts
             nvecout = numpftso
          else if ( dimids(1) == dimidldu )then
             nvecin  = numldus
             nvecout = numlduso
          else
             write (6,*) 'Skip 1D variable with unknown dimension: ', trim(varname)
             cycle
          end if

          if ( xtype == NF90_INT )then
             call interp_sl_int( varname, ncidi, ncido, nvec=nvecin, nveco=nvecout )
          else if ( xtype == NF90_DOUBLE )then
             call interp_sl_real( varname, ncidi, ncido, nvec=nvecin, nveco=nvecout )
          else
             write (6,*) 'error: variable is not of type double or integer'; stop
          end if

       ! For RTM variables
       else if ( (dimids(2) == dimidrtmlat) .and. (dimids(1) == dimidrtmlon) )then

          ! Only copy the liquid water in, leave ice alone. This seems to solve problems
          ! where the ocean blows up for negative ice flow.
          ! An alternative solution that has not been tested, yet, eliminates both "if ( index..." lines
          ! so that all RTM variables get mapped from the input to the output file. (slevis)

          ! If anything BUT RTM_VOLR_LIQ -- go to next variable
          if ( index(varname,"RTM_VOLR_LIQ" ) /= 1 )then
             write (6,*) 'Do NOT copy anything but RTM_VOLR_LIQ: ', trim(varname)
             cycle
          end if
          call check_ret(nf90_inq_varid(ncidi, varname, varid))
          call check_ret(nf90_get_var(ncidi, varid, volr))

          ret = nf90_inq_varid(ncido, varname, varid)
          if (ret==NF90_ENOTVAR)then
             write (6,*) 'Variable NOT found on output file: ', trim(varname)
             cycle
          else if (ret/=NF90_NOERR)then
             call handle_error (ret)
          end if

          !$OMP PARALLEL DO PRIVATE (i,j)
          do j  = 1, rtmlat
             do l  = 1, rtmlon
                if ( shr_infnan_isnan(volr(l,j)) ) volr(l,j) = spval
             end do
          end do
          !$OMP END PARALLEL DO

          write (6,*) 'RTM variable copied over: ', trim(varname)
          call check_ret(nf90_put_var(ncido, varid, volr))
          

       ! For 2D variables
       else if ( ndims == 2 )then

          if ( xtype /= NF90_DOUBLE )then
             write (6,*) 'error: 2D variable is not of double type:', trim(varname)
             stop
          end if
          if ( dimids(1) == dimidrad )then
             if ( dimids(2) == dimidpft )then
                ret = nf90_inq_varid(ncido, varname, varid)
                if (ret==NF90_ENOTVAR)then
                   write (6,*) 'Variable NOT found on output file: ', trim(varname)
                   cycle
                else if (ret/=NF90_NOERR)then
                   call handle_error (ret)
                end if
                call check_ret(nf90_get_var(ncido, varid, rbufmlo))
                !$OMP PARALLEL DO PRIVATE (n,k)
                do n = 1, numpftso
                   do k = 1, numrad
                      if ( shr_infnan_isnan(rbufmlo(k,n)) ) rbufmlo(k,n) = spval
                   end do
                end do
                !$OMP END PARALLEL DO
                call check_ret(nf90_put_var(ncido, varid, rbufmlo))
                write (6,*) 'copied and cleaned variable with numrad dimension: ', trim(varname)
             else
                write (6,*) 'Skipping variable with numrad dimension: ', trim(varname)
             end if
             cycle
          end if
          if ( dimids(1) == dimidcan )then
             if ( dimids(2) == dimidpft )then
                ret = nf90_inq_varid(ncido, varname, varid)
                if (ret==NF90_ENOTVAR)then
                   write (6,*) 'Variable NOT found on output file: ', trim(varname)
                   cycle
                else if (ret/=NF90_NOERR)then
                   call handle_error (ret)
                end if
                call check_ret(nf90_get_var(ncido, varid, rbufmco))
                !$OMP PARALLEL DO PRIVATE (n,k)
                do n = 1, numpftso
                   do k = 1, nlevcan
                      if ( shr_infnan_isnan(rbufmco(k,n)) ) rbufmco(k,n) = spval
                   end do
                end do
                !$OMP END PARALLEL DO
                call check_ret(nf90_put_var(ncido, varid, rbufmco))
                write (6,*) 'copied and cleaned variable with levcan dimension: ', trim(varname)
             else
                write (6,*) 'Skipping variable with levcan dimension: ', trim(varname)
             end if
             cycle
          end if
          if ( dimids(2) /= dimidcols .and. dimids(2) /= dimidpft )then
             write (6,*) 'error: variable = ', varname
             write (6,*) 'error: variables second dimension is not recognized'; stop
          end if
          if ( dimids(1) == dimidlak )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevlak, nlev_o=nlevlak, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidtot )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevtot, nlev_o=nlevtot_o, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidgrnd )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevgrnd, nlev_o=nlevgrnd_o, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidsno )then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevsno, nlev_o=nlevsno, nvec=numcols, nveco=numcolso)
          else if ( dimids(1) == dimidsno1)then
             call interp_ml_real(varname, ncidi, ncido, &
                                 nlev=nlevsno1, nlev_o=nlevsno1, nvec=numcols, nveco=numcolso)
          else
             write (6,*) 'error: variable = ', varname
             write (6,*) 'error: variables first dimension is not recognized'; stop
          end if
       else
          write (6,*) 'Skipping variable NOT 1 or 2D: ', trim(varname)
       end if
       call shr_sys_flush(6)
    end do

    ! Close input and output files

    call check_ret(nf90_close( ncidi))
    call check_ret(nf90_close( ncido))

    write (6,*) ' Successfully created initial condition file mapped from input IC file'

  end subroutine interp_filei

  !=======================================================================

  subroutine findMinDistPFTs( ncidi, ncido, pftindx )

    ! Find the PFT distances based on the column distances already calculated

    use netcdf
    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(out) :: pftindx(:)         ! vector number
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    real(r8), allocatable :: lati(:)         
    real(r8), allocatable :: loni(:)         
    real(r8), allocatable :: cos_lati(:)     
    real(r8), allocatable :: lato(:)         
    real(r8), allocatable :: lono(:)         
    real(r8), allocatable :: cos_lato(:)     
    integer , allocatable :: ltypei(:)
    integer , allocatable :: ltypeo(:)
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    real(r8) :: dx,dy,distmin,dist    
    integer  :: n,no,nmin,ier
    integer  :: ret                !NetCDF return code
    integer  :: varid              !netCDF variable id
    ! --------------------------------------------------------------------
    !
    ! Distances for PFT's to output index no
    !
    ier = 0
    allocate (lati(numpfts), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: lati'
       call shr_sys_flush(6)
       stop
    end if
    allocate (loni(numpfts), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: loni'
       call shr_sys_flush(6)
       stop
    end if
    allocate (cos_lati(numpfts), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: cos_lati'
       call shr_sys_flush(6)
       stop
    end if
    allocate (lato(numpftso), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: lato'
       call shr_sys_flush(6)
       stop
    end if
    allocate (lono(numpftso), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: lono'
       call shr_sys_flush(6)
       stop
    end if
    allocate (cos_lato(numpftso), stat=ier)
    if (ier /= 0) then
       write(6,*) 'allocation error: cos_lato'
       call shr_sys_flush(6)
       stop
    end if

    allocate (ltypei(numpfts))
    allocate (vtypei(numpfts))
    allocate (wti   (numpfts))
     
    allocate (ltypeo(numpftso))
    allocate (vtypeo(numpftso))
    allocate (wto   (numpftso))
    
    ! input 

    call check_ret(nf90_inq_varid (ncidi, 'pfts1d_lon', varid))
    call check_ret(nf90_get_var(ncidi, varid, loni))

    call check_ret(nf90_inq_varid (ncidi, 'pfts1d_lat', varid))
    call check_ret(nf90_get_var(ncidi, varid, lati))

    call check_ret(nf90_inq_varid(ncidi, 'pfts1d_ityplun', varid))
    call check_ret(nf90_get_var(ncidi, varid, ltypei))
    
    call check_ret(nf90_inq_varid( ncidi, 'pfts1d_itypveg', varid))
    call check_ret(nf90_get_var( ncidi, varid, vtypei))
    
    call check_ret(nf90_inq_varid (ncidi, 'pfts1d_wtxy', varid))
    call check_ret(nf90_get_var(ncidi, varid, wti))
    
    ! output

    call check_ret(nf90_inq_varid (ncido, 'pfts1d_lon', varid))
    call check_ret(nf90_get_var(ncido, varid, lono))

    call check_ret(nf90_inq_varid (ncido, 'pfts1d_lat', varid))
    call check_ret(nf90_get_var(ncido, varid, lato))

    call check_ret(nf90_inq_varid(ncido, 'pfts1d_ityplun', varid))
    call check_ret(nf90_get_var(ncido, varid, ltypeo))
    
    call check_ret(nf90_inq_varid( ncido, 'pfts1d_itypveg', varid))
    call check_ret(nf90_get_var( ncido, varid, vtypeo))
    
    call check_ret(nf90_inq_varid (ncido, 'pfts1d_wtxy', varid))
    call check_ret(nf90_get_var(ncido, varid, wto))

    do n = 1, numpfts
       lati(n) = lati(n)*deg2rad
       loni(n) = loni(n)*deg2rad
       cos_lati(n) = cos(lati(n))
    end do

    do n = 1, numpftso
       lato(n) = lato(n)*deg2rad
       lono(n) = lono(n)*deg2rad
       cos_lato(n) = cos(lato(n))
    end do

    write(6,*)'numpftso = ',numpftso,' numpfts= ',numpfts
    pftindx(:) = 0
    !$OMP PARALLEL DO PRIVATE (no,n,nmin,distmin,dx,dy,dist)
    do no = 1,numpftso
       if (wto(no)>0.) then 

          nmin    = 0
          distmin = spval

          do n = 1, numpfts
             if (wti(n)>0. .and. (ltypei(n) == ltypeo(no)) ) then
                if ( allPFTSfromSameGC .or. &
                     (ltypeo(no) > istsoil) .or. &
                     (vtypei(n) == vtypeo(no)) )then
                   dy   = abs(lato(no)-lati(n))*re
                   dx   = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                   dist = dx*dx + dy*dy
                   if ( dist < distmin ) then
                      distmin = dist
                      nmin    = n
                   end if
                end if
             end if
          end do

          ! If output pft type is not contained in input dataset, then use closest bare soil pft
          if ( override_missing ) then
             if (distmin == spval) then
                do n = 1, numpfts
                   if (wti(n) > 0._r8 .and. ltypei(n) == istsoil .and. vtypei(n)==baresoil) then
                      dy   = abs(lato(no)-lati(n))*re
                      dx   = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                      dist = dx*dx + dy*dy
                      if ( dist < distmin )then
                         distmin = dist
                         nmin    = n
                      end if
                   end if
                end do
                if ( distmin == spval )then
                   write(*,*) 'findMinDistPFTs: Can not find the closest pft: ',&
                        ' no,ltypeo,vtypeo=', no,ltypeo(no),vtypeo(no)
                   stop
                end if
             end if
          end if
           
          pftindx(no) = nmin
       end if  ! end if wto>0 block
    end do
    !$OMP END PARALLEL DO

    deallocate (loni)
    deallocate (lono)
    deallocate (lati)
    deallocate (lato)
    deallocate (cos_lati)
    deallocate (cos_lato)
    deallocate (ltypei) 
    deallocate (vtypei)
    deallocate (wti)   
    deallocate (ltypeo) 
    deallocate (vtypeo)
    deallocate (wto)   
    
  end subroutine findMinDistPFTs

  !=======================================================================

  subroutine findMinDistCols( ncidi, ncido, colindx )

    ! Find the minimun column distances excluding columns of different type

    use netcdf
    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(out) :: colindx(:)         ! n = colindx(no) 
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    real(r8), allocatable :: lati(:)         
    real(r8), allocatable :: loni(:)         
    real(r8), allocatable :: cos_lati(:)     
    real(r8), allocatable :: lato(:)         
    real(r8), allocatable :: lono(:)         
    real(r8), allocatable :: cos_lato(:)     
    integer , allocatable :: typei(:)
    integer , allocatable :: typeo(:)
    integer , allocatable :: typei_urb(:)
    integer , allocatable :: typeo_urb(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    real(r8) :: dx,dy,distmin,dist
    integer  :: n,no,nmin
    integer  :: varid   
    logical  :: calcmin
    integer  :: ret     
    ! --------------------------------------------------------------------

    allocate (lati(numcols))
    allocate (lato(numcolso))

    allocate (loni(numcols))
    allocate (lono(numcolso))

    allocate (typei_urb(numcols))
    allocate (typeo_urb(numcolso))

    allocate (cos_lati(numcols))
    allocate (cos_lato(numcolso))

    allocate (typei(numcols))
    allocate (typeo(numcolso))

    allocate (wti(numcols))
    allocate (wto(numcolso))

    ! input

    call check_ret(nf90_inq_varid (ncidi, 'cols1d_lon', varid))
    call check_ret(nf90_get_var(ncidi, varid, loni))

    call check_ret(nf90_inq_varid (ncidi, 'cols1d_lat', varid))
    call check_ret(nf90_get_var(ncidi, varid, lati))

    call check_ret(nf90_inq_varid (ncidi, 'cols1d_ityplun', varid))
    call check_ret(nf90_get_var(ncidi, varid, typei))

    call check_ret(nf90_inq_varid (ncidi, 'cols1d_wtxy', varid))
    call check_ret(nf90_get_var(ncidi, varid, wti))

    call check_ret(nf90_inq_varid( ncidi, 'cols1d_ityp', varid ) )
    call check_ret(nf90_get_var(ncidi, varid, typei_urb))

    ! output

    call check_ret(nf90_inq_varid (ncido, 'cols1d_lon', varid))
    call check_ret(nf90_get_var(ncido, varid, lono))

    call check_ret(nf90_inq_varid (ncido, 'cols1d_lat', varid))
    call check_ret(nf90_get_var(ncido, varid, lato))

    call check_ret(nf90_inq_varid (ncido, 'cols1d_ityplun', varid))
    call check_ret(nf90_get_var(ncido, varid, typeo))

    call check_ret(nf90_inq_varid (ncido, 'cols1d_wtxy', varid))
    call check_ret(nf90_get_var(ncido, varid, wto))

    call check_ret(nf90_inq_varid( ncido, 'cols1d_ityp', varid ))
    call check_ret(nf90_get_var(ncido, varid, typeo_urb))

    do n = 1, numcols
       lati(n) = lati(n)*deg2rad
       loni(n) = loni(n)*deg2rad
       cos_lati(n) = cos(lati(n))
    end do

    do n = 1, numcolso
       lato(n) = lato(n)*deg2rad
       lono(n) = lono(n)*deg2rad
       cos_lato(n) = cos(lato(n))
    end do

    write(6,*)'numcolso = ',numcolso
    colindx(:) = 0
    !$OMP PARALLEL DO PRIVATE (no,n,nmin,distmin,dx,dy,dist,calcmin)
    do no = 1,numcolso

       if (wto(no) > 0.) then

          distmin = spval
          nmin    = 0

          do n = 1, numcols
             calcmin = .false.
             if (wti(n) > 0.0_r8) then
                if (typei_urb(n) == nonurbcol) then
                   if (typei(n) == typeo(no)) calcmin = .true.
                else
                   if (typei(n) == typeo(no) .and. typei_urb(n) == typeo_urb(no)) calcmin = .true.
                end if
             end if
             if (calcmin) then
                dy = abs(lato(no)-lati(n))*re
                dx = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                dist = dx*dx + dy*dy
                if ( dist < distmin )then
                   distmin = dist
                   nmin = n
                end if
             end if
          end do
             
          ! If input does not have output column type than use closest soil column if override is set 
          if ( override_missing ) then
             if ( distmin == spval )then
                do n = 1, numcols
                   if (wti(n) > 0._r8 .and. typei(n)==istsoil) then
                      dy = abs(lato(no)-lati(n))*re
                      dx = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                      dist = dx*dx + dy*dy
                      if ( dist < distmin )then
                         distmin = dist
                         nmin = n
                      end if
                   end if
                end do
                if (distmin == spval) then
                   write(*,*) 'findMinDistCols: Can not find the closest column: no,typeo=',&
                        no,typeo(no)
                   stop
                end if
             end if
          end if
             
          ! Determine input column index (nmin) for the given output no value
          colindx(no) = nmin
       end if
    end do
    !$OMP END PARALLEL DO

    deallocate (lati)
    deallocate (lato)
    deallocate (loni)
    deallocate (lono)
    deallocate (cos_lati)
    deallocate (cos_lato)
    deallocate (typei)
    deallocate (typeo)
    deallocate (wti)
    deallocate (wto)
    deallocate(typei_urb)
    deallocate(typeo_urb)

  end subroutine findMinDistCols

  !=======================================================================

  subroutine findMinDistLDUs( ncidi, ncido, lduindx)

    ! Find the minimun column distances excluding columns of different type

    use netcdf
    implicit none

    ! ------------------------ arguments ---------------------------------
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(out) :: lduindx(:)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    real(r8), allocatable :: lati(:)         
    real(r8), allocatable :: loni(:)         
    real(r8), allocatable :: cos_lati(:)     
    real(r8), allocatable :: lato(:)         
    real(r8), allocatable :: lono(:)         
    real(r8), allocatable :: cos_lato(:)     
    integer , allocatable :: typei(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    real(r8) :: dx,dy,distmin,dist
    integer  :: n,no,nmin
    integer  :: varid              
    integer  :: ret                
    ! --------------------------------------------------------------------

    allocate (loni(numldus))
    allocate (lono(numlduso))

    allocate (lati(numldus))
    allocate (lato(numlduso))

    allocate (cos_lati(numldus))
    allocate (cos_lato(numlduso))

    allocate (typei(numldus))
    allocate (typeo(numlduso))

    allocate (wti(numldus))
    allocate (wto(numlduso))
    
    ! input

    call check_ret(nf90_inq_varid (ncidi, 'land1d_lon', varid))
    call check_ret(nf90_get_var (ncidi, varid, loni))

    call check_ret(nf90_inq_varid (ncidi, 'land1d_lat', varid))
    call check_ret(nf90_get_var (ncidi, varid, lati))

    call check_ret(nf90_inq_varid( ncidi, 'land1d_ityplun', varid))
    call check_ret(nf90_get_var( ncidi, varid, typei))

    call check_ret(nf90_inq_varid( ncidi, 'land1d_wtxy', varid))
    call check_ret(nf90_get_var( ncidi, varid, wti))

    ! output

    call check_ret(nf90_inq_varid (ncido, 'land1d_lon', varid))
    call check_ret(nf90_get_var (ncido, varid, lono))

    call check_ret(nf90_inq_varid (ncido, 'land1d_lat', varid))
    call check_ret(nf90_get_var (ncido, varid, lato))

    call check_ret(nf90_inq_varid( ncido, 'land1d_ityplun', varid))
    call check_ret(nf90_get_var( ncido, varid, typeo))

    call check_ret(nf90_inq_varid( ncido, 'land1d_wtxy', varid))
    call check_ret(nf90_get_var( ncido, varid, wto))
    
    do n = 1, numldus
       lati(n) = lati(n)*deg2rad
       loni(n) = loni(n)*deg2rad
       cos_lati(n) = cos(lati(n))
    end do

    do n = 1, numlduso
       lato(n) = lato(n)*deg2rad
       lono(n) = lono(n)*deg2rad
       cos_lato(n) = cos(lato(n))
    end do

    lduindx(:) = 0
    !$OMP PARALLEL DO PRIVATE (no,n,nmin,distmin,dx,dy,dist)
    do no = 1,numlduso

       if (wto(no) > 0.) then
          distmin = spval
          nmin    = 0
          
          do n = 1, numldus
             if ( (wti(n) > 0.0_r8) .and. (typei(n) == typeo(no)) ) then
                dy = abs(lato(no)-lati(n))*re
                dx = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                dist = dx*dx + dy*dy
                if ( dist < distmin ) then
                   distmin = dist
                   nmin    = n
                end if
             end if
          end do
          
          ! If input does not output landunit, then use closest soil landunit
          if ( override_missing ) then
             if ( distmin == spval )then
                do n = 1, numldus
                   if (wti(n) > 0._r8 .and. typei(n) == istsoil) then
                      dy = abs(lato(no)-lati(n))*re
                      dx = abs(lono(no)-loni(n))*re * 0.5_r8*(cos_lato(no)+cos_lati(n))
                      dist = dx*dx + dy*dy
                      if ( dist < distmin )then
                         distmin = dist
                         nmin    = n
                      end if
                   end if
                end do
             end if ! end temporary code
             if ( distmin == spval )then
                write(*,*) 'findMinDistLDUs: Can not find the closest landunit: ',&
                     'no,typeo=', no,typeo(no)
                stop
             end if
          end if
          
          lduindx(no) = nmin
       end if

    end do
    !$OMP END PARALLEL DO

    deallocate (loni)
    deallocate (lono)
    deallocate (lati)
    deallocate (lato)
    deallocate (cos_lati)
    deallocate (cos_lato)
    deallocate(typei)
    deallocate(typeo)
    deallocate(wti)
    deallocate(wto)

  end subroutine findMinDistLDUs

  !=======================================================================

  subroutine interp_ml_real (varname, ncidi, ncido, nlev, nlev_o, nvec, nveco)

    use netcdf
    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=*), intent(in) :: varname     ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nlev               ! number of levels
    integer , intent(in)  :: nlev_o             ! number of levels
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: n,no                             ! indices
    integer :: varid                            ! variable id
    real(r8), allocatable :: rbufmli (:,:)      ! input array
    real(r8), allocatable :: rbufmlo (:,:)      ! output array
    real(r8), allocatable :: wto(:)             ! weight  
    integer :: ret                              ! netCDF return code
    integer :: nlev_diff
    integer, parameter :: nlev_inactive = 5
    ! --------------------------------------------------------------------

    allocate (rbufmli(nlev,nvec))
    allocate (rbufmlo(nlev_o,nveco))
    allocate (wto(nveco))

    if (nveco == numcolso) then
       call check_ret(nf90_inq_varid(ncido, 'cols1d_wtxy', varid))
       call check_ret(nf90_get_var(ncido, varid, wto))
    else if (nveco == numpftso) then
       call check_ret(nf90_inq_varid( ncido, 'pfts1d_wtxy', varid))
       call check_ret(nf90_get_var( ncido, varid, wto))
    end if

    call check_ret(nf90_inq_varid (ncidi, trim(varname), varid))
    call check_ret(nf90_get_var(ncidi, varid, rbufmli))
    call check_ret(nf90_inq_varid (ncido, trim(varname), varid))
    call check_ret(nf90_get_var( ncido, varid, rbufmlo))

    if (nlev == nlev_o) then
       if (nvec == numcols) then
          do no = 1, nveco
             if (wto(no)>0._r8) then
                n = colindx(no)
                if (n > 0) rbufmlo(:,no) = rbufmli(:,n)
             end if
          end do
       else if (nvec == numpfts) then
          do no = 1, nveco
             if (wto(no)>0._r8) then
                n = pftindx(no)
                if (n > 0) rbufmlo(:,no) = rbufmli(:,n)
             end if
          end do
       else
          write(*,*) 'no data was written: subroutine interp_ml_real'
          stop
       end if
    else
       !!! here we repeat variables at the depth of nsoil for each of the new levels
       nlev_diff = nlev_o - nlev
       if (nlev_diff .lt. 0 ) then
          write(*,*) 'error: new grid must be longer than old grid'
          stop
       end if
       if (nvec == numcols) then
          do no = 1, nveco
             if (wto(no)>0._r8) then
                n = colindx(no)
                if (n > 0) rbufmlo(1:nlev-nlev_inactive,no) = rbufmli(1:nlev-nlev_inactive,n)
                if (n > 0) rbufmlo(1+nlev-nlev_inactive:nlev_o-nlev_inactive,no) = rbufmli(nlev-nlev_inactive,n)
                if (n > 0) rbufmlo(1+nlev_o-nlev_inactive:nlev_o,no) = rbufmli(1+nlev-nlev_inactive:nlev,n)
             end if
          end do
       else if (nvec == numpfts) then
          do no = 1, nveco
             if (wto(no)>0._r8) then
                n = pftindx(no)
                if (n > 0) rbufmlo(1:nlev-nlev_inactive,no) = rbufmli(1:nlev-nlev_inactive,n)
                if (n > 0) rbufmlo(1+nlev-nlev_inactive:nlev_o-nlev_inactive,no) = rbufmli(nlev-nlev_inactive,n)
                if (n > 0) rbufmlo(1+nlev_o-nlev_inactive:nlev_o,no) = rbufmli(1+nlev-nlev_inactive:nlev,n)
             end if
          end do
       else
          write(*,*) 'no data was written: subroutine interp_ml_real'
          stop
       end if
    endif

    call check_ret(nf90_inq_varid (ncido, varname, varid))
    call check_ret(nf90_put_var(ncido, varid, rbufmlo))

    write(*,*) 'wrote variable ',trim(varname),' to output file'

    deallocate(rbufmli)
    deallocate(rbufmlo)
    deallocate(wto)

  end subroutine interp_ml_real

  !=======================================================================

  subroutine interp_sl_real (varname, ncidi, ncido, nvec, nveco)

    use netcdf
    use shr_infnan_mod , only: shr_infnan_isnan

    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=256), intent(in) :: varname   ! input variable name 
    integer , intent(in)  :: ncidi              ! input netCdf id
    integer , intent(in)  :: ncido              ! output netCDF id  
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo,l,g,gg    ! indices
    integer :: varid                       ! variable id
    integer :: dimid                       ! dimension id
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: wto(:)
    real(r8), allocatable :: rbufsli (:)   ! input array
    real(r8), allocatable :: rbufslo (:)   ! output array
    integer :: ret                         ! NetCDF return code
    integer :: num                         ! number of gridcells NOT normalized
    integer :: numgrdso                    ! number of gridcells on output grid
    logical :: htop_var                    ! If variable name is == htop/hbot
    logical :: fpcgrid_var                 ! If variable name is == fpcgrid
    ! --------------------------------------------------------------------

    allocate (rbufsli(nvec))
    allocate (rbufslo(nveco))
    allocate (wto(nveco))

    htop_var    = .false.
    fpcgrid_var = .false.
    if (nvec == numpfts) then
       if ( trim(varname) == 'htop'    ) htop_var    = .true.
       if ( trim(varname) == 'hbot'    ) htop_var    = .true.
       if ( trim(varname) == 'fpcgrid' ) fpcgrid_var = .true.
    end if
    if ( htop_var .or. fpcgrid_var )then
       allocate (vtypeo(nveco))
       allocate (typeo(nveco))
    end if

    if (nveco == numcolso) then
       call check_ret(nf90_inq_varid (ncido, 'cols1d_wtxy', varid))
       call check_ret(nf90_get_var(ncido, varid, wto))
    else if (nveco == numpftso) then
       call check_ret(nf90_inq_varid( ncido, 'pfts1d_wtxy', varid ))
       call check_ret(nf90_get_var( ncido, varid, wto))
       if ( htop_var .or. fpcgrid_var )then
          call check_ret(nf90_inq_varid (ncido, 'pfts1d_itypveg', varid))
          call check_ret(nf90_get_var(ncido, varid, vtypeo))
          call check_ret(nf90_inq_varid (ncido, 'pfts1d_ityplun', varid))
          call check_ret(nf90_get_var( ncido, varid, typeo))
       end if
    else if (nveco == numlduso) then
       call check_ret(nf90_inq_varid (ncido, 'land1d_wtxy', varid))
       call check_ret(nf90_get_var(ncido, varid, wto))
    end if

    call check_ret(nf90_inq_varid (ncidi, varname, varid))
    call check_ret(nf90_get_var(ncidi, varid, rbufsli))
    call check_ret(nf90_inq_varid( ncido, varname, varid))
    call check_ret(nf90_get_var( ncido, varid, rbufslo))

    if ( nvec == numcols )then

       do no = 1, nveco
          if (wto(no)>0._r8) then
             n = colindx(no)
             if (n > 0) rbufslo(no) = rbufsli(n)
          end if  
       end do

    else if ( nvec == numldus )then

       do no = 1, nveco
          if (wto(no)>0._r8) then
             n = lduindx(no)
             if ( shr_infnan_isnan(rbufsli(n)) ) then
                if (n > 0) rbufsli(n) = spval
             end if
             if (n > 0) rbufslo(no) = rbufsli(n)
          end if         
       end do

    else if ( nvec == numpfts )then

       do no = 1, nveco
          if (wto(no)>0._r8) then
             !
             ! If variable-name is htop or fpcgrid
             !
             ! NB: fpcgrid and fpcgridold not needed on output file if that
             !     file will be used for a startup (not restart) simulation.
             !     I think that interpinic is intended for startup runs.
             !
             !     However,
             !     the fpcgrid interpolation/mapping needs to go on the output
             !     file as variable pfts1d_wtcol and/or PFT_WTCOL, which
             !     updates fpcgrid and fpcgridold in clm subr pftwt_init.
             !
             !     I removed fpcgrid normalization code from here because,
             !     if fpcgrid for an output gridcell all comes from a single
             !     input gridcell, which it should (check it!),
             !     then normalization is unnecessary (slevis).
             !
             if ( htop_var )then
                ! AND this is non-vegetated land-unit -- set to zero
                if( typeo(no) /= istsoil .and. typeo(no) /= istcrop )then
                   rbufslo(no) = 0.0_r8
                else
                   ! Otherwise calculate it from the nearest neighbor
                   n = pftindx(no)
                   if (n > 0) rbufslo(no) = rbufsli(n)
                end if
             else if ( fpcgrid_var )then
                ! AND this is not naturally vegetated land-unit -- set to zero
                if( typeo(no) > istsoil )then
                   if (n > 0) rbufslo(no) = 0.0_r8
                else
                   ! Otherwise calculate it from the nearest neighbor
                   n = pftindx(no)
                   if (n > 0) rbufslo(no) = rbufsli(n)
                end if
             else
                ! Otherwise calculate it from the nearest neighbor
                n = pftindx(no)
                if (n > 0) then
                   if ( shr_infnan_isnan(rbufsli(n)) ) then
                      rbufsli(n) = spval
                   end if
                   rbufslo(no) = rbufsli(n)
                else
                   rbufslo(no) = spval
                end if
             end if          !data type
          end if             !output data with positive weight
       end do                !output data land loop
    else
       write(*,*) 'subroutine interp_sl_real: no data written to variable ',varname       
       stop
    end if

    call check_ret(nf90_inq_varid (ncido, varname, varid))
    call check_ret(nf90_put_var(ncido, varid, rbufslo))

    write(*,*) 'wrote variable ', trim(varname),' to output initial file'

    deallocate(rbufsli)
    deallocate(rbufslo)

    if ( allocated(vtypei) ) deallocate(vtypei)
    if ( allocated(vtypeo) ) deallocate(vtypeo)
    if ( allocated(typeo) )  deallocate(typeo)
    deallocate(wto)

  end subroutine interp_sl_real

  !=======================================================================

  subroutine interp_sl_int (varname, ncidi, ncido, nvec, nveco)

    use netcdf
    implicit none
    include 'netcdf.inc'

    ! ------------------------ arguments ---------------------------------
    character(len=256), intent(in) :: varname
    integer , intent(in)  :: ncidi
    integer , intent(in)  :: ncido
    integer , intent(in)  :: nvec               ! number of points
    integer , intent(in)  :: nveco              ! number of points
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer :: i,j,k,n,no,ni,noo                !indices
    integer :: varid                            !variable id
    integer , allocatable :: vtypei(:)
    integer , allocatable :: vtypeo(:)
    integer , allocatable :: typeo(:)
    real(r8), allocatable :: wti(:)
    real(r8), allocatable :: wto(:)
    integer , allocatable :: ibufsli (:)        !input array
    integer , allocatable :: ibufslo (:)        !output array
    integer :: count
    integer :: ret                              !NetCDF return code
    logical :: present_var                      !If variable name is == present
    logical :: itypveg_var                      !If variable name is == itypveg
    ! --------------------------------------------------------------------

    allocate (ibufsli(nvec))
    allocate (ibufslo(nveco))
    allocate (wto(nveco))

    present_var = .false.
    itypveg_var = .false.
    if (nvec == numpfts) then
       if ( trim(varname) == 'present' ) present_var = .true.
       if ( trim(varname) == 'pfts1d_itypveg' ) itypveg_var = .true.
    end if

    if (nvec == numcols) then
       allocate (wti(nvec))
       call check_ret(nf90_inq_varid (ncidi, 'cols1d_wtxy', varid))
       call check_ret(nf90_get_var(ncidi, varid, wti))
    else if (nvec == numpfts) then
       if ( present_var .or. itypveg_var )then
          allocate (vtypei(nvec))
          call check_ret(nf90_inq_varid (ncidi, 'pfts1d_itypveg', varid))
          call check_ret(nf90_get_var(ncidi, varid, vtypei))
       end if
    end if

    if (nveco == numcolso) then
       allocate (typeo(nveco))
       call check_ret(nf90_inq_varid (ncido, 'cols1d_ityplun', varid))
       call check_ret(nf90_get_var(ncido, varid, typeo))
       call check_ret(nf90_inq_varid (ncido, 'cols1d_wtxy', varid))
       call check_ret(nf90_get_var(ncido, varid, wto))
    else if (nveco == numpftso) then
       allocate (typeo(nveco))
       allocate (vtypeo(nveco))
       call check_ret(nf90_inq_varid (ncido, 'pfts1d_itypveg', varid))
       call check_ret(nf90_get_var(ncido, varid, vtypeo))
       call check_ret(nf90_inq_varid (ncido, 'pfts1d_ityplun', varid))
       call check_ret(nf90_get_var(ncido, varid, typeo))
       call check_ret(nf90_inq_varid (ncido, 'pfts1d_wtxy', varid))
       call check_ret(nf90_get_var(ncido, varid, wto))
    end if

    call check_ret(nf90_inq_varid (ncidi, varname, varid))
    call check_ret(nf90_get_var(ncidi, varid, ibufsli))

    call check_ret(nf90_inq_varid( ncido, varname, varid))
    call check_ret(nf90_get_var( ncido, varid, ibufslo ))

    if ( nvec == numcols )then

       do no = 1, nveco
          if (wto(no)>0._r8) then
             n = colindx(no)
             if (n > 0) ibufslo(no) = ibufsli(n)
          end if 
       end do

    else if ( nvec == numpfts )then

       do no = 1, nveco
          if (wto(no)>0._r8) then
             ! If variable-name is present or itypveg 
             ! AND this is crop or non-vegetated land-unit
             if ( (present_var .or. itypveg_var) .and. &
                  ((vtypeo(no) >= croptype) .or. (typeo(no) > istsoil)) )then
                if ( present_var ) ibufslo(no) = 0
                if ( itypveg_var ) ibufslo(no) = vtypeo(no)
             ! Otherwise calculate it from the nearest neighbor
             else
                n = pftindx(no)
                if (n > 0) ibufslo(no) = ibufsli(n)
             end if          !variable type
          end if             !output data with positive weight
       end do                !output data land loop

    else

       write(*,*) 'subroutine interp_sl_int: no data written to typeo,vtypeo,no=', &
                   typeo(no),vtypeo(no),no
       stop

    end if

    call check_ret(nf90_inq_varid (ncido, varname, varid))
    call check_ret(nf90_put_var (ncido, varid, ibufslo))

    write(*,*) 'wrote variable ', trim(varname),' to output file'

    deallocate (ibufsli)
    deallocate (ibufslo)
    deallocate (wto)

    if ( allocated(vtypei) ) deallocate (vtypei)
    if ( allocated(vtypeo) ) deallocate (vtypeo)
    if ( allocated(typeo) )  deallocate (typeo)
    if ( allocated(wti) )    deallocate (wti)

  end subroutine interp_sl_int

  !=======================================================================

  subroutine addglobal (ncid, cmdline)

    implicit none
    include 'netcdf.inc'

    integer, intent(in) :: ncid
    character(len=*), intent(in) ::  cmdline

    ! ------------------------ local variables -----------------------------
    integer :: ret
    integer :: numchars
    integer :: values(8)
    integer :: hnum
    integer :: hlen
    character(len= 8) :: date
    character(len=10) :: time
    character(len= 5) :: zone
    character(len=18) :: datetime
    character(len=256):: version = &
         "$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_71/models/lnd/clm/tools/clm4_5/interpinic/src/interpinic.F90 $"
    character(len=256)  :: revision_id = "$Id: interpinic.F90 54953 2013-11-06 16:29:45Z sacks $"
    character(len=16)   :: logname
    character(len=16)   :: hostname
    character(len=256)  :: str
    character(len=16000) :: hist
#ifdef _OPENMP
    external :: OMP_GET_MAX_THREADS
    integer :: OMP_GET_MAX_THREADS
#endif
    !-----------------------------------------------------------------------

    call date_and_time (date, time, zone, values)

    datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '

    call getenv ('LOGNAME', logname)
    call getenv ('HOST', hostname)

    ret = nf_redef( ncid )
    if (ret/=NF_NOERR) call handle_error (ret)
    hlen = 0
    hist = ' '
    if (nf_inq_attid (ncid, nf_global, 'history', hnum) == nf_noerr) then
       ret = nf_inq_attlen (ncid, nf_global, 'history', hlen)
       if (ret/=NF_NOERR) call handle_error (ret)
       ret = nf_get_att_text (ncid, nf_global, 'history', hist)
       if (ret/=NF_NOERR) call handle_error (ret)
    end if

    hist = trim (hist) // char(10) // datetime // trim (logname) // ':' // &
         trim (hostname) // ':' // trim (cmdline)

    ! Add "3" to account for first newline and colons between each of 2 trimmed strings

    hlen = hlen + len(datetime) + len_trim(logname) + len_trim(hostname) + len_trim(cmdline) + 3

    if (hlen > len(hist)) then
       write(6,*)'Warning: history attribute too long: truncating'
       hlen = len(hist)
    end if

    numchars = len_trim (hist)
    ret = nf_put_att_text (ncid, nf_global, 'history', numchars, hist)
    if (ret/=NF_NOERR) call handle_error (ret)

    write(6,*) "Add SVN_version and Id to global file attributes"
    numchars = len_trim (version)
    ret = nf_put_att_text (ncid, nf_global, 'interpinic_version', numchars, version)
    if (ret/=NF_NOERR) call handle_error (ret)
    numchars = len_trim (revision_id)
    ret = nf_put_att_text (ncid, nf_global, 'interpinic_version_Id', numchars, revision_id)
    if (ret/=NF_NOERR) call handle_error (ret)

#ifdef _OPENMP
    str = 'OMP_NUM_THREADS'
    ret = nf_put_att_int (ncid, nf_global, str, NF_INT, 1, &
                   OMP_GET_MAX_THREADS() )
    if (ret/=NF_NOERR) call handle_error (ret)
    str = 'TRUE'
#else
    str = 'FALSE'
#endif
    numchars = len_trim (str)
    ret = nf_put_att_text (ncid, nf_global, 'OpenMP', numchars, str)
    if (ret/=NF_NOERR) call handle_error (ret)

#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif

    numchars = len_trim (str)
    ret = nf_put_att_text (ncid, nf_global, 'Compiler_Optimized', numchars, str)
    if (ret/=NF_NOERR) call handle_error (ret)

    ret = nf_enddef( ncid )
    if (ret/=NF_NOERR) call handle_error (ret)

  end subroutine addglobal

  !=======================================================================

  subroutine check_ret(ret)
    implicit none
    include 'netcdf.inc'
    integer, intent(in) :: ret
    if (ret /= NF_NOERR) then
       write(6,*)'netcdf error rcode = ', ret,' error = ', NF_STRERROR(ret)
       call abort()
    end if
  end subroutine check_ret

  !=======================================================================
 
  subroutine handle_error (ret)
    implicit none
    include 'netcdf.inc'
    integer ret
    write(6,*) "NetCDF error code = ", ret
    write(6,*) nf_strerror (ret)
    call abort
  end subroutine handle_error

end module interpinic
