module SatellitePhenologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! CLM Satelitte Phenology model (SP) ecosystem dynamics (phenology, vegetation). 
  ! Allow some subroutines to be used by the CLM Carbon Nitrogen model (CLMCN) 
  ! so that DryDeposition code can get estimates of LAI differences between months.
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use clm_varctl   , only : scmlat,scmlon,single_column
  use clm_varctl   , only : iulog
  use perf_mod     , only : t_startf, t_stopf
  use spmdMod      , only : masterproc
  use ncdio_pio   
  use decompMod    , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SatellitePhenology     ! CLMSP Ecosystem dynamics: phenology, vegetation
  public :: SatellitePhenologyInit ! Dynamically allocate memory
  public :: interpMonthlyVeg       ! interpolate monthly vegetation data
  public :: readAnnualVegetation   ! Read in annual vegetation (needed for Dry-deposition)
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: readMonthlyVegetation   ! read monthly vegetation data for two months
  !
  ! !PRIVATE TYPES:
  integer , private :: InterpMonths1            ! saved month index
  real(r8), private :: timwt(2)                 ! time weights for month 1 and month 2
  real(r8), private, allocatable :: mlai2t(:,:) ! lai for interpolation (2 months)
  real(r8), private, allocatable :: msai2t(:,:) ! sai for interpolation (2 months)
  real(r8), private, allocatable :: mhvt2t(:,:) ! top vegetation height for interpolation (2 months)
  real(r8), private, allocatable :: mhvb2t(:,:) ! bottom vegetation height for interpolation(2 months)
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SatellitePhenologyInit (bounds)
    !
    ! !DESCRIPTION:
    ! Dynamically allocate memory and set to signaling NaN.
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: ier    ! error code
    !-----------------------------------------------------------------------

    InterpMonths1 = -999  ! saved month index

    ier = 0
    if(.not.allocated(mlai2t)) then
       allocate (mlai2t(bounds%begp:bounds%endp,2), &
            msai2t(bounds%begp:bounds%endp,2), &
            mhvt2t(bounds%begp:bounds%endp,2), &
            mhvb2t(bounds%begp:bounds%endp,2), stat=ier)
    end if
    if (ier /= 0) then
       write(iulog,*) 'EcosystemDynini allocation error'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    mlai2t(bounds%begp : bounds%endp, :) = nan
    msai2t(bounds%begp : bounds%endp, :) = nan
    mhvt2t(bounds%begp : bounds%endp, :) = nan
    mhvb2t(bounds%begp : bounds%endp, :) = nan

  end subroutine SatellitePhenologyInit

  !-----------------------------------------------------------------------
  subroutine SatellitePhenology(bounds, num_nolakep, filter_nolakep, doalb)
    !
    ! !DESCRIPTION:
    ! Ecosystem dynamics: phenology, vegetation
    ! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
    ! height (htop).
    !
    ! !USES:
    use clmtype
    use pftvarcon, only : noveg, nc3crop, nbrdlf_dcd_brl_shrub
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                          ! bounds
    integer, intent(in) :: num_nolakep                               ! number of column non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(bounds%endp-bounds%begp+1) ! pft filter for non-lake points
    logical, intent(in) :: doalb                                     ! true = surface albedo calculation time step
    !
    ! !LOCAL VARIABLES:
                                               ! local pointers to implicit in arguments
    real(r8), pointer :: frac_sno(:)           ! fraction of ground covered by snow (0 to 1)
    integer , pointer :: pcolumn(:)            ! column index associated with each pft
    real(r8), pointer :: snow_depth(:)         ! snow height (m)
    integer , pointer :: ivt(:)                ! pft vegetation type
                                               ! local pointers to implicit out arguments
    real(r8), pointer :: tlai(:)               ! one-sided leaf area index, no burying by snow
    real(r8), pointer :: tsai(:)               ! one-sided stem area index, no burying by snow
    real(r8), pointer :: htop(:)               ! canopy top (m)
    real(r8), pointer :: hbot(:)               ! canopy bottom (m)
    real(r8), pointer :: elai(:)               ! one-sided leaf area index with burying by snow
    real(r8), pointer :: esai(:)               ! one-sided stem area index with burying by snow
    integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
                                               !    
    integer  :: fp,p,c                         ! indices
    real(r8) :: ol                             ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb                             ! fraction of canopy layer covered by snow
    !-----------------------------------------------------------------------

    if (doalb) then

       ! Assign local pointers to derived type scalar members (column-level)

       frac_sno => cps%frac_sno 
       snow_depth  => cps%snow_depth

       ! Assign local pointers to derived type scalar members (pftlevel)

       pcolumn =>pft%column
       tlai    => pps%tlai
       tsai    => pps%tsai
       elai    => pps%elai
       esai    => pps%esai
       htop    => pps%htop
       hbot    => pps%hbot
       frac_veg_nosno_alb => pps%frac_veg_nosno_alb
       ivt     =>pft%itype

       do fp = 1, num_nolakep
          p = filter_nolakep(fp)
          c = pcolumn(p)

          ! need to update elai and esai only every albedo time step so do not
          ! have any inconsistency in lai and sai between SurfaceAlbedo calls (i.e.,
          ! if albedos are not done every time step).
          ! leaf phenology
          ! Set leaf and stem areas based on day of year
          ! Interpolate leaf area index, stem area index, and vegetation heights
          ! between two monthly
          ! The weights below (timwt(1) and timwt(2)) were obtained by a call to
          ! routine InterpMonthlyVeg in subroutine NCARlsm.
          !                 Field   Monthly Values
          !                -------------------------
          ! leaf area index LAI  <- mlai1 and mlai2
          ! leaf area index SAI  <- msai1 and msai2
          ! top height      HTOP <- mhvt1 and mhvt2
          ! bottom height   HBOT <- mhvb1 and mhvb2

          tlai(p) = timwt(1)*mlai2t(p,1) + timwt(2)*mlai2t(p,2)
          tsai(p) = timwt(1)*msai2t(p,1) + timwt(2)*msai2t(p,2)
          htop(p) = timwt(1)*mhvt2t(p,1) + timwt(2)*mhvt2t(p,2)
          hbot(p) = timwt(1)*mhvb2t(p,1) + timwt(2)*mhvb2t(p,2)

          ! adjust lai and sai for burying by snow. if exposed lai and sai
          ! are less than 0.05, set equal to zero to prevent numerical
          ! problems associated with very small lai and sai.

          ! snow burial fraction for short vegetation (e.g. grasses) as in
          ! Wang and Zeng, 2007. 

          if (ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
             ol = min( max(snow_depth(c)-hbot(p), 0._r8), htop(p)-hbot(p))
             fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
          else
             fb = 1._r8 - max(min(snow_depth(c),0.2_r8),0._r8)/0.2_r8   ! 0.2m is assumed
                  !depth of snow required for complete burial of grasses
          endif

          ! area weight by snow covered fraction
          elai(p) = max(tlai(p)*(1.0_r8 - frac_sno(c)) &
                       +tlai(p)*fb*frac_sno(c), 0.0_r8)
          esai(p) = max(tsai(p)*(1.0_r8 - frac_sno(c)) &
                       +tsai(p)*fb*frac_sno(c), 0.0_r8)
          if (elai(p) < 0.05_r8) elai(p) = 0._r8
          if (esai(p) < 0.05_r8) esai(p) = 0._r8

          ! Fraction of vegetation free of snow

          if ((elai(p) + esai(p)) >= 0.05_r8) then
             frac_veg_nosno_alb(p) = 1
          else
             frac_veg_nosno_alb(p) = 0
          end if

       end do ! end of pft loop

    end if  !end of if-doalb block

  end subroutine SatellitePhenology

  !-----------------------------------------------------------------------
  subroutine interpMonthlyVeg (bounds)
    !
    ! !DESCRIPTION:
    ! Determine if 2 new months of data are to be read.
    !
    ! !USES:
    use clm_varctl      , only : fsurdat
    use clm_time_manager, only : get_curr_date, get_step_size, get_nstep
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: kyr         ! year (0, ...) for nstep+1
    integer :: kmo         ! month (1, ..., 12)
    integer :: kda         ! day of month (1, ..., 31)
    integer :: ksec        ! seconds into current date for nstep+1
    real(r8):: dtime       ! land model time step (sec)
    real(r8):: t           ! a fraction: kda/ndaypm
    integer :: it(2)       ! month 1 and month 2 (step 1)
    integer :: months(2)   ! months to be interpolated (1 to 12)
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
    !-----------------------------------------------------------------------

    dtime = get_step_size()

    call get_curr_date(kyr, kmo, kda, ksec, offset=int(dtime))

    t = (kda-0.5_r8) / ndaypm(kmo)
    it(1) = t + 0.5_r8
    it(2) = it(1) + 1
    months(1) = kmo + it(1) - 1
    months(2) = kmo + it(2) - 1
    if (months(1) <  1) months(1) = 12
    if (months(2) > 12) months(2) = 1
    timwt(1) = (it(1)+0.5_r8) - t
    timwt(2) = 1._r8-timwt(1)

    if (InterpMonths1 /= months(1)) then
       if (masterproc) then
          write(iulog,*) 'Attempting to read monthly vegetation data .....'
          write(iulog,*) 'nstep = ',get_nstep(),' month = ',kmo,' day = ',kda
       end if
       call t_startf('readMonthlyVeg')
       call readMonthlyVegetation (bounds, fsurdat, months)
       InterpMonths1 = months(1)
       call t_stopf('readMonthlyVeg')
    end if

  end subroutine interpMonthlyVeg

  !-----------------------------------------------------------------------
  ! read 12 months of veg data for dry deposition
  !-----------------------------------------------------------------------

  subroutine readAnnualVegetation (bounds)

    use clmtype
    use clm_varpar  , only : numpft
    use pftvarcon   , only : noveg
    use domainMod   , only : ldomain
    use fileutils   , only : getfil
    use clm_varctl  , only : fsurdat
    use shr_scam_mod, only : shr_scam_getCloseLatLon

    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds

    ! local vars

    type(file_desc_t) :: ncid             ! netcdf id
    real(r8), pointer :: annlai(:,:)      ! 12 months of monthly lai from input data set 
    real(r8), pointer :: mlai(:,:)        ! lai read from input files
    integer :: ier                        ! error code
    character(len=256) :: locfn           ! local file name
    integer :: g,k,l,m,n,p,ivt            ! indices
    integer :: ni,nj,ns                   ! indices
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: closelatidx,closelonidx    ! single column vars
    real(r8):: closelat,closelon          ! single column vars
    logical :: isgrid2d                   ! true => file is 2d
    character(len=32) :: subname = 'readAnnualVegetation'

    annlai    => pps%annlai

    ! Determine necessary indices

    allocate(mlai(bounds%begg:bounds%endg,0:numpft), stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname, 'allocation error ' 
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (masterproc) then
       write (iulog,*) 'Attempting to read annual vegetation data .....'
    end if

    call getfil(fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqfdims (ncid, isgrid2d, ni, nj, ns)

    if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
       write(iulog,*)trim(subname), 'ldomain and input file do not match dims '
       write(iulog,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
       write(iulog,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
       write(iulog,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    call check_dim(ncid, 'lsmpft', numpft+1)

    if (single_column) then
       call shr_scam_getCloseLatLon(locfn, scmlat, scmlon, &
            closelat, closelon, closelatidx, closelonidx)
    endif

    do k=1,12   !! loop over months and read vegetated data

       call ncd_io(ncid=ncid, varname='MONTHLY_LAI', flag='read', data=mlai, &
            dim1name=grlnd, nt=k)

       !! store data directly in clmtype structure
       !! only vegetated pfts have nonzero values
       !! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
       !! as determined in subroutine surfrd

       do p = bounds%begp,bounds%endp
          g =pft%gridcell(p)
          ivt =pft%itype(p)
          if (ivt /= noveg) then     !! vegetated pft
             do l = 0, numpft
                if (l == ivt) then
                   annlai(k,p) = mlai(g,l)
                end if
             end do
          else                       !! non-vegetated pft
             annlai(k,p) = 0._r8
          end if
       end do   ! end of loop over pfts  

    enddo ! months loop

    call ncd_pio_closefile(ncid)

    deallocate(mlai)

  endsubroutine readAnnualVegetation

  !-----------------------------------------------------------------------
  subroutine readMonthlyVegetation (bounds, fveg, months)
    !
    ! !DESCRIPTION:
    ! Read monthly vegetation data for two consec. months.
    !
    ! !USES:
    use clmtype
    use clm_varpar  , only : numpft
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
    use spmdMod     , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use shr_scam_mod, only : shr_scam_getCloseLatLon
    use clm_time_manager, only : get_nstep
    use netcdf
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    character(len=*), intent(in) :: fveg  ! file with monthly vegetation data
    integer, intent(in) :: months(2)      ! months to be interpolated (1 to 12)
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: g,n,k,l,m,p,ivt,ni,nj,ns   ! indices
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: ier                        ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon
    logical :: readvar
    real(r8), pointer :: mlai(:,:)        ! lai read from input files
    real(r8), pointer :: msai(:,:)        ! sai read from input files
    real(r8), pointer :: mhgtt(:,:)       ! top vegetation height
    real(r8), pointer :: mhgtb(:,:)       ! bottom vegetation height
    real(r8), pointer :: mlaidiff(:)      ! difference between lai month one and month two
    character(len=32) :: subname = 'readMonthlyVegetation'
    !-----------------------------------------------------------------------

    ! Determine necessary indices

    allocate(&
         mlai(bounds%begg:bounds%endg,0:numpft), &
         msai(bounds%begg:bounds%endg,0:numpft), &
         mhgtt(bounds%begg:bounds%endg,0:numpft), &
         mhgtb(bounds%begg:bounds%endg,0:numpft), &
         stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from gridcell to pft data
    ! ----------------------------------------------------------------------
   
    call getfil(fveg, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    
    if (single_column) then
       call shr_scam_getCloseLatLon (ncid, scmlat, scmlon, closelat, closelon,&
            closelatidx, closelonidx)
    endif
    
    do k=1,2   !loop over months and read vegetated data

       call ncd_io(ncid=ncid, varname='MONTHLY_LAI', flag='read', data=mlai, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_LAI NOT on fveg file'//errMsg(__FILE__, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_SAI', flag='read', data=msai, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_SAI NOT on fveg file'//errMsg(__FILE__, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', flag='read', data=mhgtt, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_HEIGHT_TOP NOT on fveg file'//errMsg(__FILE__, __LINE__))

       call ncd_io(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', flag='read', data=mhgtb, dim1name=grlnd, &
            nt=months(k), readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: MONTHLY_HEIGHT_TOP NOT on fveg file'//errMsg(__FILE__, __LINE__))

       ! Store data directly in clmtype structure
       ! only vegetated pfts have nonzero values
       ! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
       ! as determined in subroutine surfrd

       do p = bounds%begp,bounds%endp
          g =pft%gridcell(p)
          ivt =pft%itype(p)
          if (ivt /= noveg) then     ! vegetated pft
             do l = 0, numpft
                if (l == ivt) then
                   mlai2t(p,k) = mlai(g,l)
                   msai2t(p,k) = msai(g,l)
                   mhvt2t(p,k) = mhgtt(g,l)
                   mhvb2t(p,k) = mhgtb(g,l)
                end if
             end do
          else                        ! non-vegetated pft
             mlai2t(p,k) = 0._r8
             msai2t(p,k) = 0._r8
             mhvt2t(p,k) = 0._r8
             mhvb2t(p,k) = 0._r8
          end if
       end do   ! end of loop over pfts

    end do   ! end of loop over months

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       k = 2
       write(iulog,*) 'Successfully read monthly vegetation data for'
       write(iulog,*) 'month ', months(k)
       write(iulog,*)
    end if

    deallocate(mlai, msai, mhgtt, mhgtb)

    mlaidiff => pps%mlaidiff
    do p = bounds%begp,bounds%endp
       mlaidiff(p) = mlai2t(p,1)-mlai2t(p,2)
    enddo

  end subroutine readMonthlyVegetation

end module SatellitePhenologyMod
