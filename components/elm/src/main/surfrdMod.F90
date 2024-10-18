module surfrdMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains methods for reading in surface data file and determining
  ! subgrid weights
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use elm_varpar      , only : nlevsoifl, numpft, numcft
  use landunit_varcon , only : numurbl
  use elm_varcon      , only : grlnd
  use elm_varctl      , only : iulog, scmlat, scmlon, single_column, firrig_data
  use elm_varctl      , only : create_glacier_mec_landunit
  use surfrdUtilsMod  , only : check_sums_equal_1_2d, check_sums_equal_1_3d
  use surfrdUtilsMod  , only : collapse_crop_types, collapse_crop_var
  use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
  use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, check_dim, ncd_inqdid, ncd_inqdlen
  use pio

  use spmdMod                         
  use topounit_varcon , only : max_topounits, has_topounit  

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd_get_globmask  ! Reads global land mask (needed for setting domain decomp)
  public :: surfrd_get_grid      ! Read grid/ladnfrac data into domain (after domain decomp)
  public :: surfrd_get_topo      ! Read grid topography into domain (after domain decomp)
  public :: surfrd_get_data      ! Read surface dataset and determine subgrid weights
  public :: surfrd_get_grid_conn ! Reads grid connectivity information from domain file
  public :: surfrd_topounit_data ! Read topounit physical properties
  public :: surfrd_get_topo_for_solar_rad    ! Read topography dataset for TOP solar radiation parameterization
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_special             ! Read the special landunits
  private :: surfrd_veg_all             ! Read all of the vegetated landunits
  private :: surfrd_pftformat           ! Read crop pfts in file format where they are part of the vegetated land unit
  private :: surfrd_cftformat           ! Read crop pfts in file format where they are on their own landunit
  private :: surfrd_fates_nocropmod     ! Read in crop and pfts and compress them both into wt_nat_patch
                                        ! Used with FATES when no crop model active
  !
  ! !PRIVATE DATA MEMBERS:
  ! default multiplication factor for epsilon for error checks
  real(r8), private, parameter :: eps_fact = 2._r8
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine surfrd_get_globmask(filename, mask, ni, nj)
    !
    ! !DESCRIPTION:
    ! Read the surface dataset grid related information:
    ! This is the first routine called by clm_initialize 
    ! NO DOMAIN DECOMPOSITION  HAS BEEN SET YET
    !
    ! !USES:
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)    :: filename  ! grid filename
    integer         , pointer       :: mask(:)   ! grid mask 
    integer         , intent(out)   :: ni, nj    ! global grid sizes
    !
    ! !LOCAL VARIABLES:
    logical :: isgrid2d
    integer :: dimid,varid         ! netCDF id's
    integer :: ns                  ! size of grid on file
    integer :: n,i,j               ! index 
    integer :: ier                 ! error status
    type(file_desc_t)  :: ncid     ! netcdf id
    type(var_desc_t)   :: vardesc  ! variable descriptor
    character(len=256) :: varname  ! variable name
    character(len=256) :: locfn    ! local file name
    logical :: readvar             ! read variable in or not
    integer , allocatable :: idata2d(:,:)
    character(len=32) :: subname = 'surfrd_get_globmask' ! subroutine name
    !-----------------------------------------------------------------------

    if (filename == ' ') then
       mask(:) = 1
       RETURN
    end if

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions and if grid file is 2d or 1d

    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    if (masterproc) then
       write(iulog,*)'lat/lon grid flag (isgrid2d) is ',isgrid2d
    end if

    allocate(mask(ns))
    mask(:) = 1

    if (isgrid2d) then
       allocate(idata2d(ni,nj))
       idata2d(:,:) = 1
       call ncd_io(ncid=ncid, varname='LANDMASK', data=idata2d, flag='read', readvar=readvar)
       if (.not. readvar) then
          call ncd_io(ncid=ncid, varname='mask', data=idata2d, flag='read', readvar=readvar)
       end if
       if (readvar) then
          do j = 1,nj
          do i = 1,ni
             n = (j-1)*ni + i
             mask(n) = idata2d(i,j)
          enddo
          enddo
       end if
       deallocate(idata2d)
    else
       call ncd_io(ncid=ncid, varname='LANDMASK', data=mask, flag='read', readvar=readvar)
       if (.not. readvar) then
          call ncd_io(ncid=ncid, varname='mask', data=mask, flag='read', readvar=readvar)
       end if
    end if
    if (.not. readvar) call endrun( msg=' ERROR: landmask not on fatmlndfrc file'//errMsg(__FILE__, __LINE__))

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_globmask

  !-----------------------------------------------------------------------
  subroutine surfrd_get_grid(begg, endg, ldomain, filename, glcfilename)
    !
    ! !DESCRIPTION:
    ! THIS IS CALLED AFTER THE DOMAIN DECOMPOSITION HAS BEEN CREATED
    ! Read the surface dataset grid related information:
    ! o real latitude  of grid cell (degrees)
    ! o real longitude of grid cell (degrees)
    !
    ! !USES:
    use elm_varcon, only : spval, re
    use domainMod , only : domain_type, domain_init, domain_clean, lon1d, lat1d
    use fileutils , only : getfil
    use elm_varctl, only : use_pflotran
    !
    ! !ARGUMENTS:
    integer          ,intent(in)    :: begg, endg 
    type(domain_type),intent(inout) :: ldomain   ! domain to init
    character(len=*) ,intent(in)    :: filename  ! grid filename
    character(len=*) ,optional, intent(in) :: glcfilename ! glc mask filename
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid               ! netcdf id
    type(file_desc_t) :: ncidg              ! netCDF id for glcmask
    type(var_desc_t)  :: vardesc            ! variable descriptor
    integer :: beg                          ! local beg index
    integer :: end                          ! local end index
    integer :: ni,nj,ns                     ! size of grid on file
    integer :: dimid,varid                  ! netCDF id's
    integer :: start(1), count(1)           ! 1d lat/lon array sections
    integer :: ier,ret                      ! error status
    logical :: readvar                      ! true => variable is on input file 
    logical :: isgrid2d                     ! true => file is 2d lat/lon
    logical :: istype_domain                ! true => input file is of type domain
    real(r8), allocatable :: rdata2d(:,:)   ! temporary
    real(r8), allocatable :: rdata3d(:,:,:) ! temporary  ! pflotran
    character(len=16) :: vname              ! temporary
    character(len=256):: locfn              ! local file name
    integer :: n                            ! indices
    real(r8):: eps = 1.0e-12_r8             ! lat/lon error tolerance

    ! pflotran:beg-----------------------------
    integer :: j, np, nv

    ! pflotran:end-----------------------------
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    
    ! pflotran:beg-----------------------------------------------
    call ncd_inqdlen(ncid, dimid, nv, 'nv')
    if (nv>0) then
       ldomain%nv = nv
    else
       ldomain%nv = 0
    endif
    ! pflotran:end-----------------------------------------------

    ! Determine isgrid2d flag for domain
    ldomain%set = .false.
    call domain_init(ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=begg, nend=endg)

    ! Determine type of file - old style grid file or new style domain file
    call check_var(ncid=ncid, varname='LONGXY', vardesc=vardesc, readvar=readvar) 
    if (readvar) istype_domain = .false.

    call check_var(ncid=ncid, varname='xc', vardesc=vardesc, readvar=readvar) 
    if (readvar) istype_domain = .true.

    ! Read in area, lon, lat

    if (istype_domain) then
       call ncd_io(ncid=ncid, varname= 'area', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       ! convert from radians**2 to km**2
       ldomain%area = ldomain%area * (re**2)
       if (.not. readvar) call endrun( msg=' ERROR: area NOT on file'//errMsg(__FILE__, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'xc', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: xc NOT on file'//errMsg(__FILE__, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'yc', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: yc NOT on file'//errMsg(__FILE__, __LINE__))

       ! pflotran:beg-----------------------------------------------
       ! user-defined grid-cell vertices (ususally 'nv' is 4,
       ! but for future use, we set the following if condition of 'nv>=3' so that possible to use TIN grids
       if (ldomain%nv>=3 .and. use_pflotran) then
          call ncd_io(ncid=ncid, varname='xv', flag='read', data=ldomain%lonv, &
            dim1name=grlnd, readvar=readvar)
          if (.not. readvar) call endrun( msg=trim(subname)//' ERROR: xv  NOT on file'//errMsg(__FILE__, __LINE__))

          call ncd_io(ncid=ncid, varname='yv', flag='read', data=ldomain%latv, &
            dim1name=grlnd, readvar=readvar)
          if (.not. readvar) call endrun( msg=trim(subname)//' ERROR: yv  NOT on file'//errMsg(__FILE__, __LINE__))

       end if
       ! pflotran:end-----------------------------------------------
    else
       call ncd_io(ncid=ncid, varname= 'AREA', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: AREA NOT on file'//errMsg(__FILE__, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'LONGXY', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: LONGXY NOT on file'//errMsg(__FILE__, __LINE__))
       
       call ncd_io(ncid=ncid, varname= 'LATIXY', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: LATIXY NOT on file'//errMsg(__FILE__, __LINE__))

       ! pflotran:beg-----------------------------------------------
       ! user-defined grid-cell vertices (ususally 'nv' is 4,
       ! but for future use, we set the following if condition of 'nv>=3' so that possible to use TIN grids
       if (ldomain%nv>=3 .and. use_pflotran) then

          call ncd_io(ncid=ncid, varname='LONGV', flag='read', data=ldomain%lonv, &
            dim1name=grlnd, readvar=readvar)
          if (.not. readvar) call endrun( msg=trim(subname)//' ERROR: LONGV  NOT on file'//errMsg(__FILE__, __LINE__))

          call ncd_io(ncid=ncid, varname='LATIV', flag='read', data=ldomain%latv, &
            dim1name=grlnd, readvar=readvar)
          if (.not. readvar) call endrun( msg=trim(subname)//' ERROR: LATIV  NOT on file'//errMsg(__FILE__, __LINE__))

       end if
       ! pflotran:end-----------------------------------------------
    end if

    
    ! let lat1d/lon1d data available for all grid-types, if coupled with PFLOTRAN.
    if (isgrid2d .or. use_pflotran) then
       allocate(rdata2d(ni,nj), lon1d(ni), lat1d(nj))
       if (istype_domain) then
          vname = 'xc'
       else
          vname = 'LONGXY'
       end if
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lon1d(:) = rdata2d(:,1)
       if (istype_domain) then
          vname = 'yc'
       else
          vname = 'LATIXY'
       end if
       call ncd_io(ncid=ncid, varname=trim(vname), data=rdata2d, flag='read', readvar=readvar)
       lat1d(:) = rdata2d(1,:)
       deallocate(rdata2d)

       ! pflotran:beg-----------------------------------------------
       ! find the origin of ldomain, if vertices of first grid known
       if (use_pflotran) then
         ldomain%lon0 = -9999._r8
         ldomain%lat0 = -9999._r8
         if (ldomain%nv==4 .and. ldomain%nv /= huge(1)) then
          allocate(rdata3d(ni,nj,nv))
          if (istype_domain) then
             vname = 'xv'
          else
             vname = 'LONGV'
          end if

          call ncd_io(ncid=ncid, varname=trim(vname), data=rdata3d, flag='read', readvar=readvar)

          if (readvar) then
            ldomain%lon0 = 0._r8
            np=0
            do j=1,nv
               ! may have issue if mixed longitude values (i.e. 0~360 or -180~180)
               if ( ni>1 .and. &
                    ( (rdata3d(1,1,j) < lon1d(1) .and. rdata3d(1,1,j) < lon1d(2)) .or. &
                      (rdata3d(1,1,j) > lon1d(1) .and. rdata3d(1,1,j) > lon1d(2)) ) ) then
                 np = np + 1
                 ldomain%lon0 = ldomain%lon0+rdata3d(1,1,j)

               else if (ni==1 .and. rdata3d(1,1,j)<lon1d(1)) then  !either side should be OK
                 np = np + 1
                 ldomain%lon0 = ldomain%lon0+rdata3d(1,1,j)
               end if
            end do
            if (np>0) then
              ldomain%lon0 = ldomain%lon0/np
            else
              ldomain%lon0 = -9999._r8
            end if
          end if

          !
          if (istype_domain) then
             vname = 'yv'
          else
             vname = 'LATIV'
          end if
          call ncd_io(ncid=ncid, varname=trim(vname), data=rdata3d, flag='read', readvar=readvar)
          if (readvar) then
            ldomain%lat0 = 0._r8
            np=0
            do j=1,nv
               if ( nj>1 .and. &
                    ( (rdata3d(1,1,j) < lat1d(1) .and. rdata3d(1,1,j) < lat1d(2)) .or. &
                      (rdata3d(1,1,j) > lat1d(1) .and. rdata3d(1,1,j) > lat1d(2)) ) ) then
                 np = np + 1
                 ldomain%lat0 = ldomain%lat0+rdata3d(1,1,j)

               else if (nj==1 .and. rdata3d(1,1,j)<lat1d(1)) then  !either side should be OK
                 np = np + 1
                 ldomain%lat0 = ldomain%lat0+rdata3d(1,1,j)
               end if
            end do
            if (np>0) then
              ldomain%lat0 = ldomain%lat0/np
            else
              ldomain%lat0 = -9999._r8
            end if
          end if
          !
          deallocate(rdata3d)
         end if
       end if
       ! pflotran:end-----------------------------------------------
    end if  ! if (isgrid2d .or. use_pflotran)



    ! Check lat limited to -90,90

    if (minval(ldomain%latc) < -90.0_r8 .or. &
        maxval(ldomain%latc) >  90.0_r8) then
       write(iulog,*) trim(subname),' WARNING: lat/lon min/max is ', &
            minval(ldomain%latc),maxval(ldomain%latc)
       ! call endrun( msg=' ERROR: lat is outside [-90,90]'//errMsg(__FILE__, __LINE__))
       ! write(iulog,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
       !     minval(domain%latc),maxval(domain%latc)
       ! where (ldomain%latc < -90.0_r8) ldomain%latc = -90.0_r8
       ! where (ldomain%latc >  90.0_r8) ldomain%latc =  90.0_r8
    endif

    call ncd_io(ncid=ncid, varname='LANDMASK', flag='read', data=ldomain%mask, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='mask', flag='read', data=ldomain%mask, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: LANDMASK NOT on fracdata file'//errMsg(__FILE__, __LINE__))
       end if
    end if

    call ncd_io(ncid=ncid, varname='LANDFRAC', flag='read', data=ldomain%frac, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='frac', flag='read', data=ldomain%frac, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( msg=' ERROR: LANDFRAC NOT on fracdata file'//errMsg(__FILE__, __LINE__))
       end if
    end if

    ! Read xCell
    call check_var(ncid=ncid, varname='xCell', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname= 'xCell', flag='read', data=ldomain%xCell, &
            dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='yCell', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname= 'yCell', flag='read', data=ldomain%yCell, &
            dim1name=grlnd, readvar=readvar)
    endif

    call ncd_pio_closefile(ncid)

    if (present(glcfilename)) then
       if (masterproc) then
          if (glcfilename == ' ') then
             write(iulog,*) trim(subname),' ERROR: glc filename must be specified '
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
       end if
       call getfil( glcfilename, locfn, 0 )
       call ncd_pio_openfile (ncidg, trim(locfn), 0)

       ldomain%glcmask(:) = 0
       call ncd_io(ncid=ncidg, varname='GLCMASK', flag='read', data=ldomain%glcmask, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: GLCMASK NOT in file'//errMsg(__FILE__, __LINE__))

       ! Make sure the glc mask is a subset of the land mask
       do n = begg,endg
          if (ldomain%glcmask(n)==1 .and. ldomain%mask(n)==0) then
             write(iulog,*)trim(subname),&
                  'initialize1: landmask/glcmask mismatch'
             write(iulog,*)trim(subname),&
                  'glc requires input where landmask = 0, gridcell index', n
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
       enddo
       call ncd_pio_closefile(ncidg)
    endif   ! present(glcfilename)

  end subroutine surfrd_get_grid

  !-----------------------------------------------------------------------
  subroutine surfrd_get_topo(domain,filename)
    !
    ! !DESCRIPTION:
    ! Read the topo dataset grid related information:
    ! Assume domain has already been initialized and read
    !
    ! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
    use GridcellType, only : grc_pp
    
    !
    ! !ARGUMENTS:
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid      ! netcdf file id
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: dimid,varid         ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8             ! lat/lon error tolerance
    integer :: beg,end                      ! local beg,end indices
    logical             :: isgrid2d         ! true => file is 2d lat/lon
    real(r8),pointer    :: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn            ! local file name
    logical :: readvar                      ! is variable on file
    character(len=32) :: subname = 'surfrd_get_topo'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun(msg=errMsg(__FILE__, __LINE__))
       else
          write(iulog,*) 'Attempting to read lnd topo from flndtopo ',trim(filename)
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

    if (domain%ns /= ns) then
       write(iulog,*) trim(subname),' ERROR: topo file mismatch ns',&
            domain%ns,ns
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    
    beg = domain%nbeg
    end = domain%nend

    allocate(latc(beg:end),lonc(beg:end))

    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: LONGXY  NOT on topodata file'//errMsg(__FILE__, __LINE__))

    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: LONGXY  NOT on topodata file'//errMsg(__FILE__, __LINE__))

    do n = beg,end
       if (abs(latc(n)-domain%latc(n)) > eps .or. &
           abs(lonc(n)-domain%lonc(n)) > eps) then
          write(iulog,*) trim(subname),' ERROR: topo file mismatch lat,lon',latc(n),&
               domain%latc(n),lonc(n),domain%lonc(n),eps
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    enddo

    call ncd_io(ncid=ncid, varname='TOPO', flag='read', data=domain%topo, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: LONGXY  NOT on topodata file'//errMsg(__FILE__, __LINE__))

    deallocate(latc,lonc)

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_topo

  !-----------------------------------------------------------------------
  subroutine surfrd_get_data (begg, endg, ldomain, lfsurdat)
    !
    ! !DESCRIPTION:
    ! Read the surface dataset and create subgrid weights.
    ! The model's surface dataset recognizes 6 basic land cover types within a grid
    ! cell: lake, wetland, urban, glacier, glacier_mec and vegetated. The vegetated
    ! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs. These
    ! subgrid patches are read in explicitly for each grid cell. This is in
    ! contrast to LSMv1, where the PFTs were built implicitly from biome types.
    !    o real latitude  of grid cell (degrees)
    !    o real longitude of grid cell (degrees)
    !    o integer surface type: 0 = ocean or 1 = land
    !    o integer soil color (1 to 20) for use with soil albedos
    !    o real soil texture, %sand, for thermal and hydraulic properties
    !    o real soil texture, %clay, for thermal and hydraulic properties
    !    o real % of cell covered by lake    for use as subgrid patch
    !    o real % of cell covered by wetland for use as subgrid patch
    !    o real % of cell that is urban      for use as subgrid patch
    !    o real % of cell that is glacier    for use as subgrid patch
    !    o real % of cell that is glacier_mec for use as subgrid patch
    !    o integer PFTs
    !    o real % abundance PFTs (as a percent of vegetated area)
    !
    ! !USES:
    use elm_varctl  , only : create_crop_landunit, firrig_data
    use fileutils   , only : getfil
    use domainMod   , only : domain_type, domain_init, domain_clean
    use elm_varsur  , only : wt_lunit, topo_glc_mec
    use GridcellType, only : grc_pp
    use topounit_varcon, only : max_topounits, has_topounit
    !
    ! !ARGUMENTS:
    integer,          intent(in) :: begg, endg      
    type(domain_type),intent(in) :: ldomain     ! land domain
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
    !
    ! !LOCAL VARIABLES:
    type(var_desc_t)  :: vardesc              ! pio variable descriptor
    type(domain_type) :: surfdata_domain      ! local domain associated with surface dataset
    character(len=256):: locfn                ! local file name
    integer           :: n                    ! loop indices
    integer           :: ni,nj,ns             ! domain sizes
    character(len=16) :: lon_var, lat_var     ! names of lat/lon on dataset
    logical           :: readvar              ! true => variable is on dataset
    real(r8)          :: rmaxlon,rmaxlat      ! local min/max vars
    type(file_desc_t) :: ncid                 ! netcdf id
    logical           :: istype_domain        ! true => input file is of type domain
    logical           :: isgrid2d             ! true => intut grid is 2d 
    character(len=32) :: subname = 'surfrd_get_data'    ! subroutine name
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
    endif

    topo_glc_mec(:,:,:) = 0._r8

    ! Read surface data

    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Read in pft mask - this variable is only on the surface dataset - but not
    ! on the domain dataset

    call ncd_io(ncid=ncid, varname= 'PFTDATA_MASK', flag='read', data=ldomain%pftm, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: pftm NOT on surface dataset'//errMsg(__FILE__, __LINE__))
	        
    !! Read the actual number of topounits per grid    
	!call check_var(ncid=ncid, varname='topoPerGrid', vardesc=vardesc, readvar=readvar)
    !if (readvar) then
    !   call ncd_io(ncid=ncid, varname= 'topoPerGrid', flag='read', data=ldomain%num_tunits_per_grd, &
    !     dim1name=grlnd, readvar=readvar)
    !endif    
    
	! Check if fsurdat grid is "close" to fatmlndfrc grid, exit if lats/lon > 0.001

    call check_var(ncid=ncid, varname='xc', vardesc=vardesc, readvar=readvar) 
    if (readvar) then
       istype_domain = .true.
    else
       call check_var(ncid=ncid, varname='LONGXY', vardesc=vardesc, readvar=readvar) 
       if (readvar) then
          istype_domain = .false.
       else
          call endrun( msg=' ERROR: unknown domain type'//errMsg(__FILE__, __LINE__))
       end if
    end if
    if (istype_domain) then
       lon_var  = 'xc'
       lat_var  = 'yc'
    else
       lon_var  = 'LONGXY'
       lat_var  = 'LATIXY'
    end if
    if ( masterproc )then
       write(iulog,*) trim(subname),' lon_var = ',trim(lon_var),' lat_var =',trim(lat_var)
    end if

    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)
    surfdata_domain%nv = 0   ! must be initialized to 0 here prior to call 'domain_init'
    surfdata_domain%set = .false.
    call domain_init(surfdata_domain, isgrid2d, ni, nj, begg, endg, elmlevel=grlnd)

    call ncd_io(ncid=ncid, varname=lon_var, flag='read', data=surfdata_domain%lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: lon var NOT on surface dataset'//errMsg(__FILE__, __LINE__))

    call ncd_io(ncid=ncid, varname=lat_var, flag='read', data=surfdata_domain%latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: lat var NOT on surface dataset'//errMsg(__FILE__, __LINE__))

    rmaxlon = 0.0_r8
    rmaxlat = 0.0_r8
    do n = begg,endg
       if (ldomain%lonc(n)-surfdata_domain%lonc(n) > 300.) then
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-surfdata_domain%lonc(n)-360._r8))
       elseif (ldomain%lonc(n)-surfdata_domain%lonc(n) < -300.) then
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-surfdata_domain%lonc(n)+360._r8))
       else
          rmaxlon = max(rmaxlon,abs(ldomain%lonc(n)-surfdata_domain%lonc(n)))
       endif
       rmaxlat = max(rmaxlat,abs(ldomain%latc(n)-surfdata_domain%latc(n)))
    enddo
    if (rmaxlon > 0.001_r8 .or. rmaxlat > 0.001_r8) then
       write(iulog,*)' ERROR: surfdata/fatmgrid lon/lat mismatch error', rmaxlon,rmaxlat
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if
    call domain_clean(surfdata_domain)

    ! Obtain special landunit info

    call surfrd_special(begg, endg, ncid, ldomain%ns,ldomain%num_tunits_per_grd)
    
    ! Obtain vegetated landunit info

    call surfrd_veg_all(begg, endg, ncid, ldomain%ns,ldomain%num_tunits_per_grd)

    call ncd_pio_closefile(ncid)

    !call check_sums_equal_1_3d(wt_lunit, begg, 'wt_lunit', subname,ldomain%num_tunits_per_grd)
    call check_sums_equal_1_3d(wt_lunit, begg, 'wt_lunit', subname)

    if ( masterproc )then
       write(iulog,*) 'Successfully read surface boundary data'
       write(iulog,*)
    end if

  end subroutine surfrd_get_data

!-----------------------------------------------------------------------
  subroutine surfrd_special(begg, endg, ncid, ns,ntpu)
    !
    ! !DESCRIPTION:
    ! Determine weight with respect to gridcell of all special "pfts" as well
    ! as soil color and percent sand and clay
    !
    ! !USES:
    use elm_varpar      , only : maxpatch_glcmec, nlevurb
    use landunit_varcon , only : isturb_MIN, isturb_MAX, istdlak, istwet, istice, istice_mec
    use elm_varsur      , only : wt_lunit, urban_valid, wt_glc_mec, topo_glc_mec, firrig, f_surf, f_grd
    use UrbanParamsType , only : CheckUrban
    use topounit_varcon , only : max_topounits, has_topounit
    !
    ! !ARGUMENTS:
    integer          , intent(in)    :: begg, endg 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    integer          , intent(in)    :: ns     ! domain size 
    integer          , intent(in)    :: ntpu(:)     ! Number of topounits per grid
    !
    ! !LOCAL VARIABLES:
    integer  :: n,nl,nurb,g, t,tm,ti                ! indices
    integer  :: dimid,varid                ! netCDF id's
    real(r8) :: nlevsoidata(nlevsoifl)
    logical  :: found                      ! temporary for error check
    integer  :: nindx                      ! temporary for error check
    integer  :: ier                        ! error status
    logical  :: readvar
  
    real(r8),pointer :: pctgla_old(:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak_old(:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet_old(:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb_old(:,:)    ! percent of grid cell is urbanized
    integer ,pointer :: urban_region_id_old(:) 
    real(r8),pointer :: pctglc_mec_tot_old(:) ! percent of grid cell is glacier (sum over classes)
    real(r8),pointer :: pcturb_tot_old(:)  ! percent of grid cell is urban (sum over density classes)
    real(r8),pointer :: pctspec_old(:)     ! percent of spec lunits wrt gcell
    
    real(r8),pointer :: pctgla(:,:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:,:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:,:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb(:,:,:)    ! percent of grid cell is urbanized
    integer ,pointer :: urban_region_id(:,:) 
    real(r8),pointer :: pctglc_mec_tot(:,:) ! percent of grid cell is glacier (sum over classes)
    real(r8),pointer :: pcturb_tot(:,:)  ! percent of grid cell is urban (sum over density classes)
    real(r8),pointer :: pctspec(:,:)     ! percent of spec lunits wrt gcell    
    integer  :: dens_index             ! urban density index
    character(len=32) :: subname = 'surfrd_special'  ! subroutine name
    real(r8) closelat,closelon
    integer, parameter :: urban_invalid_region = 0   ! urban_region_id indicating invalid point
!-----------------------------------------------------------------------

    allocate(pctgla(begg:endg,1:max_topounits))
    allocate(pctlak(begg:endg,1:max_topounits))
    allocate(pctwet(begg:endg,1:max_topounits))
    allocate(pcturb(begg:endg,1:max_topounits,numurbl))
    allocate(pcturb_tot(begg:endg,1:max_topounits))
    allocate(urban_region_id(begg:endg,1:max_topounits))
    allocate(pctglc_mec_tot(begg:endg,1:max_topounits))
    allocate(pctspec(begg:endg,1:max_topounits))
    
    call check_dim(ncid, 'nlevsoi', nlevsoifl)

       ! Obtain non-grid surface properties of surface dataset other than percent pft

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_WETLAND  NOT on surfdata file'//errMsg(__FILE__, __LINE__))

    call ncd_io(ncid=ncid, varname='PCT_LAKE'   , flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_LAKE NOT on surfdata file'//errMsg(__FILE__, __LINE__))

    call ncd_io(ncid=ncid, varname='PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_GLACIER NOT on surfdata file'//errMsg(__FILE__, __LINE__))

    ! Read urban info
    if (nlevurb == 0) then
      ! If PCT_URBAN is not multi-density then set pcturb to zero 
      pcturb = 0._r8
      urban_valid(begg:endg,:) = .false.
      write(iulog,*)'PCT_URBAN is not multi-density, pcturb set to 0'
    else
      call ncd_io(ncid=ncid, varname='PCT_URBAN'  , flag='read', data=pcturb, &
           dim1name=grlnd, readvar=readvar)
      if (.not. readvar) call endrun( msg=' ERROR: PCT_URBAN NOT on surfdata file'//errMsg(__FILE__, __LINE__))

      call ncd_io(ncid=ncid, varname='URBAN_REGION_ID', flag='read', data=urban_region_id, &
           dim1name=grlnd, readvar=readvar)
      if (.not. readvar) call endrun( msg= ' ERROR: URBAN_REGION_ID NOT on surfdata file'//errMsg(__FILE__, __LINE__))
      where (urban_region_id == urban_invalid_region)
         urban_valid = .false.
      elsewhere
         urban_valid = .true.
      end where
    end if
    if ( nlevurb == 0 )then
       if ( any(pcturb > 0.0_r8) ) then
          call endrun( msg=' ERROR: PCT_URBAN MUST be zero when nlevurb=0'//errMsg(__FILE__, __LINE__))
       end if
    end if

    pcturb_tot(:,:) = 0._r8
    do n = 1, numurbl
       do nl = begg,endg
         do t = 1, max_topounits
           pcturb_tot(nl,t) = pcturb_tot(nl,t) + pcturb(nl,t,n)
         end do
       enddo
    enddo

    if (create_glacier_mec_landunit) then          ! call ncd_io_gs_int2d

       call check_dim(ncid, 'nglcec',   maxpatch_glcmec   )
       call check_dim(ncid, 'nglcecp1', maxpatch_glcmec+1 )

       call ncd_io(ncid=ncid, varname='PCT_GLC_MEC', flag='read', data=wt_glc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: PCT_GLC_MEC NOT on surfdata file'//errMsg(__FILE__, __LINE__))

       wt_glc_mec(:,:,:) = wt_glc_mec(:,:,:) / 100._r8
       !call check_sums_equal_1_3d(wt_glc_mec, begg, 'wt_glc_mec', subname,ntpu)
       call check_sums_equal_1_3d(wt_glc_mec, begg, 'wt_glc_mec', subname)

       call ncd_io(ncid=ncid, varname='TOPO_GLC_MEC',  flag='read', data=topo_glc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: TOPO_GLC_MEC NOT on surfdata file'//errMsg(__FILE__, __LINE__))

       topo_glc_mec(:,:,:) = max(topo_glc_mec(:,:,:), 0._r8)


       ! Put glacier area into the GLC_MEC landunit rather than the simple glacier landunit
       pctglc_mec_tot(:,:) = pctgla(:,:)
       pctgla(:,:) = 0._r8

       pctspec = pctwet + pctlak + pcturb_tot + pctglc_mec_tot

    else

       pctglc_mec_tot(:,:) = 0._r8
       pctspec = pctwet + pctlak + pcturb_tot + pctgla
 
    endif

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg,endg
       ti = (nl - begg) + 1
       if (.not. has_topounit) then
          tm = max_topounits          
       else
          tm = ntpu(ti)
       end if
       do t = 1, tm            
         if (pctspec(nl,t) > 100._r8+1.e-04_r8) then
            found = .true.
            nindx = nl
            exit
         end if
         if (found) exit
       end do
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: PFT cover>100 for nl=',nindx
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Determine wt_lunit for special landunits

    do nl = begg,endg

         wt_lunit(nl,:,istdlak)     = pctlak(nl,:)/100._r8

         wt_lunit(nl,:,istwet)      = pctwet(nl,:)/100._r8

         wt_lunit(nl,:,istice)      = pctgla(nl,:)/100._r8

         wt_lunit(nl,:,istice_mec)  = pctglc_mec_tot(nl,:)/100._r8

         do n = isturb_MIN, isturb_MAX
            dens_index = n - isturb_MIN + 1
            wt_lunit(nl,:,n)        = pcturb(nl,:,dens_index) / 100._r8
         end do
      
    end do
    
    ! Obtain firrig and surface/grnd irrigation fraction
    if (firrig_data) then
     call ncd_io(ncid=ncid, varname='FIRRIG', flag='read', data=firrig, &
          dim1name=grlnd, readvar=readvar)
     if (.not. readvar) call endrun( trim(subname)//' ERROR: FIRRIG NOT on surfdata file' )!

     call ncd_io(ncid=ncid, varname='FSURF', flag='read', data=f_surf, &
          dim1name=grlnd, readvar=readvar)
     if (.not. readvar) call endrun( trim(subname)//' ERROR: FSURF NOT on surfdata file' )!

     call ncd_io(ncid=ncid, varname='FGRD', flag='read', data=f_grd, &
          dim1name=grlnd, readvar=readvar)
     if (.not. readvar) call endrun( trim(subname)//' ERROR: FGRD NOT on surfdata file' )
    
    else
      firrig(:,:) = 0.7_r8
      f_surf(:,:) = 1.0_r8
      f_grd(:,:) = 0.0_r8
    end if
    
    call CheckUrban(begg, endg, pcturb(begg:endg,:,:), subname,ntpu)

    deallocate(pctgla,pctlak,pctwet,pcturb,pcturb_tot,urban_region_id,pctglc_mec_tot,pctspec)

  end subroutine surfrd_special

!-----------------------------------------------------------------------
  subroutine surfrd_cftformat( ncid, begg, endg, wt_cft, fert_cft, fert_p_cft, cftsize, surfpft_size )
    !
    ! !DESCRIPTION:
    !     Handle generic crop types for file format where they are on their own
    !     crop landunit and read in as Crop Function Types.
    ! !USES:
    use elm_varsur      , only : wt_nat_patch
    use elm_varpar      , only : cft_size, cft_lb, surfpft_lb
    use topounit_varcon,  only : max_topounits
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid         ! netcdf id
    integer          , intent(in)    :: begg, endg
    integer          , intent(in)    :: cftsize      ! CFT size
    real(r8), pointer, intent(inout) :: wt_cft(:,:,:)  ! CFT weights
    integer          , intent(in)    :: surfpft_size  ! natural PFT size
    real(r8), pointer, intent(inout) :: fert_cft(:,:,:)   ! Fertilizer
    real(r8), pointer, intent(inout) :: fert_p_cft(:,:,:) ! Fertilizer
    !
    ! !LOCAL VARIABLES:
    logical  :: readvar                        ! is variable on dataset
    real(r8),pointer :: array2D(:,:,:)              ! local array
    character(len=32) :: subname = 'surfrd_cftformat'! subroutine name
!-----------------------------------------------------------------------
    SHR_ASSERT_ALL((lbound(wt_cft) == (/begg,1, cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wt_cft, dim=1) == (/endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wt_cft, dim=3) >= (/cftsize+1-cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((lbound(fert_cft) == (/begg,1, cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fert_cft, dim=1) == (/endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fert_cft, dim=3) >= (/cftsize+1-cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((lbound(fert_p_cft) == (/begg,1, cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fert_p_cft, dim=1) == (/endg/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(fert_p_cft, dim=3) >= (/cftsize+1-cft_lb/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(wt_nat_patch) >= (/endg,max_topounits,surfpft_size-1+surfpft_lb/)), errMsg(__FILE__, __LINE__))

    call check_dim(ncid, 'cft', cftsize)
    call check_dim(ncid, 'natpft', surfpft_size)

    call ncd_io(ncid=ncid, varname='PCT_CFT', flag='read', data=wt_cft, &
            dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CFT NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 

    if ( cft_size > 0 )then
       call ncd_io(ncid=ncid, varname='NFERT', flag='read', data=fert_cft, &
               dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if ( masterproc ) &
                write(iulog,*) ' WARNING: NFERT NOT on surfdata file zero out'
          fert_cft = 0.0_r8
       end if
    else
       fert_cft = 0.0_r8
    end if

    if ( cft_size > 0 )then
       call ncd_io(ncid=ncid, varname='PFERT', flag='read', data=fert_p_cft, &
               dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          if ( masterproc ) &
                write(iulog,*) ' WARNING: PFERT NOT on surfdata file zero out'
          fert_p_cft = 0.0_r8
       end if
    else
       fert_p_cft = 0.0_r8
    end if

    allocate( array2D(begg:endg,1:max_topounits, 1:surfpft_size) )
    call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=array2D, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    wt_nat_patch(begg:,:, surfpft_lb:surfpft_size-1+surfpft_lb) = array2D(begg:,:,:)
    deallocate( array2D )

  end subroutine surfrd_cftformat

!-----------------------------------------------------------------------
  subroutine surfrd_pftformat( begg, endg, ncid )
    !
    ! !DESCRIPTION:
    !     Handle generic crop types for file format where they are part of the
    !     natural vegetation landunit.
    ! !USES:
    use elm_varsur      , only : fert_cft, fert_p_cft, wt_nat_patch, wt_cft
    use elm_varpar      , only : surfpft_size, cft_size, surfpft_lb, surfpft_ub
    use elm_varpar      , only : cft_lb, cft_ub
    use elm_varctl      , only : create_crop_landunit
    use topounit_varcon,  only : max_topounits
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: begg, endg
    type(file_desc_t), intent(inout) :: ncid                    ! netcdf id
    !
    ! !LOCAL VARIABLES:
    logical  :: cft_dim_exists                 ! does the dimension 'cft' exist on the dataset?
    integer  :: dimid                          ! netCDF id's
    logical  :: readvar                        ! is variable on dataset
    real(r8),pointer :: array2D(:,:,:)                 ! local 2D array
    character(len=32) :: subname = 'surfrd_pftformat'! subroutine name
!-----------------------------------------------------------------------
    SHR_ASSERT_ALL((ubound(wt_nat_patch) == (/endg,max_topounits, surfpft_size-1+surfpft_lb/)), errMsg(__FILE__, __LINE__))


    if (.not. create_crop_landunit) then
       call check_dim(ncid, 'natpft', surfpft_size)
    else
       call check_dim(ncid, 'natpft', surfpft_size + cft_size)
    endif
    ! If cft_size == 0, then we expect to be running with a surface dataset that does
    ! NOT have a PCT_CFT array (or NFERT array), and thus does not have a 'cft' dimension.
    ! Make sure that's the case.
    call ncd_inqdid(ncid, 'cft', dimid, cft_dim_exists)
    if (cft_dim_exists) then
       call endrun( msg= ' ERROR: unexpectedly found cft dimension on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(__FILE__, __LINE__))
    end if
    call ncd_io(ncid=ncid, varname='NFERT', flag='read', data=fert_cft, &
            dim1name=grlnd, readvar=readvar)
    if (readvar) then
       call endrun( msg= ' ERROR: unexpectedly found NFERT on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(__FILE__, __LINE__))
    end if
    fert_cft = 0.0_r8

   call ncd_io(ncid=ncid, varname='PFERT', flag='read', data=fert_p_cft, &
            dim1name=grlnd, readvar=readvar)
    if (readvar) then
       call endrun( msg= ' ERROR: unexpectedly found PFERT on dataset when cft_size=0'// &
               ' (if the surface dataset has a separate crop landunit, then the code'// &
               ' must also have a separate crop landunit, and vice versa)'//&
               errMsg(__FILE__, __LINE__))
    end if
    fert_p_cft = 0.0_r8

    wt_nat_patch(begg:endg, :, :) = 0.0_r8
    wt_cft(begg:endg, :, :) = 0.0_r8
    if (.not. create_crop_landunit) then
       call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=wt_nat_patch, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    else
       allocate(array2D(begg:endg,1:max_topounits,1:cft_ub+1))
       call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=array2D, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(__FILE__, __LINE__))

       wt_nat_patch(begg:endg,1:max_topounits, surfpft_lb:surfpft_ub) = array2D(begg:endg,1:max_topounits, surfpft_lb+1:surfpft_ub+1)
       wt_cft      (begg:endg,1:max_topounits, cft_lb   :cft_ub   ) = array2D(begg:endg,1:max_topounits, cft_lb+1   :cft_ub+1   )
       deallocate(array2D)
    endif

  end subroutine surfrd_pftformat

  !-----------------------------------------------------------------------
  subroutine surfrd_veg_all(begg, endg, ncid, ns,ntpu)
    !
    ! !DESCRIPTION:
    ! Determine weight arrays for non-dynamic landuse mode
    !
    ! !USES:
    use elm_varctl      , only : create_crop_landunit, use_fates
    use elm_varctl      , only : irrigate
    use elm_varpar      , only : surfpft_lb, surfpft_ub, surfpft_size, cft_lb, cft_ub, cft_size
    use elm_varpar      , only : crop_prog
    use elm_varsur      , only : wt_lunit, wt_nat_patch, wt_cft, fert_cft, fert_p_cft
    use landunit_varcon , only : istsoil, istcrop
    use pftvarcon       , only : nc3crop, nc3irrig, npcropmin
    use pftvarcon       , only : ncorn, ncornirrig, nsoybean, nsoybeanirrig
    use pftvarcon       , only : nscereal, nscerealirrig, nwcereal, nwcerealirrig
    use pftvarcon       , only : ncassava, ncotton, nfoddergrass, noilpalm, nograins, nrapeseed 
    use pftvarcon       , only : nrice, nrtubers, nsugarcane, nmiscanthus, nswitchgrass, nwillow, npoplar
    use pftvarcon       , only : ncassavairrig, ncottonirrig, nfoddergrassirrig, noilpalmirrig, nograinsirrig, nrapeseedirrig
    use pftvarcon       , only : nriceirrig, nrtubersirrig, nsugarcaneirrig, nmiscanthusirrig, nswitchgrassirrig, nwillowirrig, npoplarirrig
    use surfrdUtilsMod  , only : convert_cft_to_pft, convert_pft_to_cft
    !
    ! !ARGUMENTS:
    integer, intent(in) :: begg, endg
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id
    integer          ,intent(in)    :: ns     ! domain size
    integer          ,intent(in)    :: ntpu(:)
    !
    ! !LOCAL VARIABLES:
    integer  :: nl, t                             ! index
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: ier                            ! error status	
    integer  :: cftsize                        ! size of CFT's
    logical  :: readvar                        ! is variable on dataset
    logical  :: cft_dim_exists                 ! does the dimension 'cft' exist on the dataset?
    real(r8),pointer :: arrayl(:,:)              ! local array
    real(r8),pointer :: array2D(:,:,:)                 ! local 2D array
    real(r8),pointer :: arrayNF(:,:,:)
    real(r8),pointer :: arrayPF(:,:,:)
    character(len=32) :: subname = 'surfrd_veg_all'  ! subroutine name
!-----------------------------------------------------------------------

    call check_dim(ncid, 'lsmpft', numpft+1)

    ! This temporary array is needed because ncd_io expects a pointer, so we can't
    ! directly pass wt_lunit(begg:endg,istsoil)
    allocate(arrayl(begg:endg,max_topounits))

    call ncd_io(ncid=ncid, varname='PCT_NATVEG', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NATVEG NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    wt_lunit(begg:endg,1:max_topounits,istsoil) = arrayl(begg:endg,1:max_topounits) 

    call ncd_io(ncid=ncid, varname='PCT_CROP', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CROP NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    wt_lunit(begg:endg,1:max_topounits,istcrop) = arrayl(begg:endg,1:max_topounits) 

    deallocate(arrayl)
    
    ! Check the file format for CFT's and handle accordingly
    call ncd_inqdid(ncid, 'cft', dimid, cft_dim_exists)
    if ( cft_dim_exists .and. create_crop_landunit ) then

       ! Format where CFT's is read in a seperate landunit
       allocate(arrayNF(begg:endg, max_topounits, cft_lb:cft_size-1+cft_lb))
       allocate(arrayPF(begg:endg, max_topounits, cft_lb:cft_size-1+cft_lb))
       call surfrd_cftformat( ncid, begg, endg, wt_cft, arrayNF, arrayPF, cft_size, surfpft_size )
       fert_cft(begg:, 1:, cft_lb:)   = arrayNF(begg:,1:, cft_lb:cft_ub)
       fert_p_cft(begg:, 1:, cft_lb:) = arrayPF(begg:,1:, cft_lb:cft_ub)
       deallocate(arrayNF)
       deallocate(arrayPF)

    else if ( (.not. cft_dim_exists) .and. (.not. create_crop_landunit) )then

       ! Format where crop is part of the natural veg. landunit
       if ( masterproc ) write(iulog,*) "WARNING: The PFT format is an unsupported format that will be removed in th future!"
       call surfrd_pftformat( begg, endg, ncid )

    else if ( cft_dim_exists .and. .not. create_crop_landunit )then

       if ( use_fates ) then
          if ( masterproc ) write(iulog,*) "WARNING: When fates is on we allow new CFT based surface datasets ", &
               "to be used with create_crop_land FALSE"
          call surfrd_fates_nocropmod( ncid, begg, endg )
          ! Set the weighting on the crop patches to zero
          fert_cft(begg:endg,:,cft_lb:cft_ub) = 0.0_r8  ! cft_lb:cft_ub has a size of zero anyway...
          fert_p_cft(begg:endg,:,cft_lb:cft_ub) = 0.0_r8
          wt_cft(begg:endg,:,cft_lb:cft_ub)   = 0.0_r8  ! cft_lb:cft_ub has a size of zero anyway...
       else
          call endrun( msg=' ERROR: New format surface datasets require create_crop_landunit TRUE'//errMsg(__FILE__, __LINE__))
       end if
       
    else
       ! PFTs contain the crops but create_crop_landunit = .true.

       ! Read the data and put them in PFT and CFT
       call surfrd_pftformat(begg, endg, ncid)

       ! Update the PFT and CFT fractions
       call convert_pft_to_cft( begg, endg )
    end if

    ! Do some checking

    if ( (cft_size) == 0 .and. any(wt_lunit(begg:endg,:,istcrop) > 0._r8)) then
       ! If cft_size == 0, and thus we aren't reading PCT_CFT, then make sure PCT_CROP is
       ! 0 everywhere (PCT_CROP > 0 anywhere requires that we have a PCT_CFT array)
       call endrun( msg=' ERROR: if PCT_CROP > 0 anywhere, then cft_size must be > 0'// &
            ' (if the surface dataset has a separate crop landunit, then the code'// &
            ' must also have a separate crop landunit, and vice versa)'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! Convert from percent to fraction, check sums of nat vegetation add to 1

    if (cft_size > 0) then
       wt_cft(begg:endg,:,:) = wt_cft(begg:endg,:,:) / 100._r8
       !call check_sums_equal_1_3d(wt_cft, begg, 'wt_cft', subname,ntpu)
       call check_sums_equal_1_3d(wt_cft, begg, 'wt_cft', subname)
    end if
    wt_lunit(begg:endg,:,istsoil) = wt_lunit(begg:endg,:,istsoil) / 100._r8
    wt_lunit(begg:endg,:,istcrop) = wt_lunit(begg:endg,:,istcrop) / 100._r8
    wt_nat_patch(begg:endg,:,:)   = wt_nat_patch(begg:endg,:,:) / 100._r8
    !call check_sums_equal_1_3d(wt_nat_patch, begg, 'wt_nat_patch', subname,ntpu)
    call check_sums_equal_1_3d(wt_nat_patch, begg, 'wt_nat_patch', subname)

    ! If no irrigation, merge irrigated CFTs with rainfed
    
    if (crop_prog .and. .not. irrigate) then
       if (masterproc) then
          write(iulog,*) trim(subname)//' crop=.T. and irrigate=.F., so merging irrigated pfts with rainfed'
       end if

       if (cft_size <= 0) then
          call endrun( msg='ERROR: Trying to merge irrigated CFTs with rainfed, but cft_size <= 0'//&
               errMsg(__FILE__, __LINE__))
       end if

       do nl = begg,endg
          do t = 1, max_topounits
            ! (TODO) the following assumes that rainfed/irrigated crop are ordered side by side
            ! indexing is fixed
            wt_cft(nl,t,nc3crop)       = wt_cft(nl,t,nc3crop)  + wt_cft(nl,t,nc3irrig)
            wt_cft(nl,t,nc3irrig)      = 0._r8
            wt_cft(nl,t,ncorn)         = wt_cft(nl,t,ncorn)    + wt_cft(nl,t,ncornirrig)
            wt_cft(nl,t,ncornirrig)    = 0._r8
            wt_cft(nl,t,nscereal)      = wt_cft(nl,t,nscereal) + wt_cft(nl,t,nscerealirrig)
            wt_cft(nl,t,nscerealirrig) = 0._r8
            wt_cft(nl,t,nwcereal)      = wt_cft(nl,t,nwcereal) + wt_cft(nl,t,nwcerealirrig)
            wt_cft(nl,t,nwcerealirrig) = 0._r8
            wt_cft(nl,t,nsoybean)      = wt_cft(nl,t,nsoybean) + wt_cft(nl,t,nsoybeanirrig)
            wt_cft(nl,t,nsoybeanirrig) = 0._r8
            wt_cft(nl,t,ncassava)      = wt_cft(nl,t,ncassava) + wt_cft(nl,t,ncassavairrig)
            wt_cft(nl,t,ncassavairrig) = 0._r8
            wt_cft(nl,t,ncotton)       = wt_cft(nl,t,ncotton)  + wt_cft(nl,t,ncottonirrig)
            wt_cft(nl,t,ncottonirrig)  = 0._r8
            wt_cft(nl,t,nfoddergrass)  = wt_cft(nl,t,nfoddergrass) + wt_cft(nl,t,nfoddergrassirrig)
            wt_cft(nl,t,nfoddergrassirrig) = 0._r8
            wt_cft(nl,t,noilpalm)      = wt_cft(nl,t,noilpalm) + wt_cft(nl,t,noilpalmirrig)
            wt_cft(nl,t,noilpalmirrig) = 0._r8
            wt_cft(nl,t,nograins)      = wt_cft(nl,t,nograins) + wt_cft(nl,t,nograinsirrig)
            wt_cft(nl,t,nograinsirrig) = 0._r8
            wt_cft(nl,t,nrapeseed)       = wt_cft(nl,t,nrapeseed)  + wt_cft(nl,t,nrapeseedirrig)
            wt_cft(nl,t,nrapeseedirrig)  = 0._r8
            wt_cft(nl,t,nrice)           = wt_cft(nl,t,nrice)  + wt_cft(nl,t,nriceirrig)
            wt_cft(nl,t,nriceirrig)      = 0._r8
            wt_cft(nl,t,nrtubers)        = wt_cft(nl,t,nrtubers) + wt_cft(nl,t,nrtubersirrig)
            wt_cft(nl,t,nrtubersirrig)   = 0._r8
            wt_cft(nl,t,nsugarcane)      = wt_cft(nl,t,nsugarcane) + wt_cft(nl,t,nsugarcaneirrig)
            wt_cft(nl,t,nsugarcaneirrig) = 0._r8
            wt_cft(nl,t,nmiscanthus)     = wt_cft(nl,t,nmiscanthus) + wt_cft(nl,t,nmiscanthusirrig)
            wt_cft(nl,t,nmiscanthusirrig) = 0._r8
            wt_cft(nl,t,nswitchgrass)     = wt_cft(nl,t,nswitchgrass)  + wt_cft(nl,t,nswitchgrassirrig)
            wt_cft(nl,t,nswitchgrassirrig) = 0._r8
            wt_cft(nl,t,npoplar)         = wt_cft(nl,t,npoplar) + wt_cft(nl,t,npoplar)
            wt_cft(nl,t,npoplarirrig)    = 0._r8
            wt_cft(nl,t,nwillow)         = wt_cft(nl,t,nwillow) + wt_cft(nl,t,nwillowirrig)
            wt_cft(nl,t,nwillowirrig)    = 0._r8

          end do
       end do

       call check_sums_equal_1_3d(wt_cft, begg, 'wt_cft', subname)
    end if
    if (crop_prog) then
       ! Call collapse_crop_types: allows need to maintain only 50-pft input data
       ! For use_crop = .false. collapsing 50->16 pfts or 16->16 or some new
       !    configuration
       ! For use_crop = .true. most likely collapsing 78 to the list of crops for
       !    which the CLM includes parameterizations
       ! The call collapse_crop_types also appears in subroutine dyncrop_interp
       call collapse_crop_types(wt_cft(begg:endg,:,:), fert_cft(begg:endg,:,:), fert_p_cft(begg:endg,:,:), begg, endg, verbose=.true.)

       ! Collapse crop variables as needed
       ! The call to collapse_crop_var also appears in subroutine dyncrop_interp
       ! - fert_cft TODO Is this call redundant because it simply sets the crop
       !                 variable to 0 where is_pft_known_to_model = .false.?
       call collapse_crop_var(fert_cft(begg:endg,:,:), begg, endg)
       call collapse_crop_var(fert_p_cft(begg:endg,:,:), begg, endg)
    end if

  end subroutine surfrd_veg_all

  !-----------------------------------------------------------------------
  subroutine surfrd_get_grid_conn(filename, cellsOnCell, edgesOnCell, &
       nEdgesOnCell, areaCell, dcEdge, dvEdge, &
       nCells_loc, nEdges_loc, maxEdges)
    !
    ! !DESCRIPTION:
    ! Read grid connectivity information.
    ! NO DOMAIN DECOMPOSITION  HAS BEEN SET YET
    !
    ! !USES:
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in)  :: filename                        ! filename
    integer         , pointer     :: cellsOnCell(:,:)                ! cells-to-cell connection
    integer         , pointer     :: edgesOnCell(:,:)                ! index to determine distance between neighbors from dcEdge
    integer         , pointer     :: nEdgesOnCell(:)                 ! number of edges
    real(r8)        , pointer     :: dcEdge(:)                       ! distance between centroids of grid cells
    real(r8)        , pointer     :: dvEdge(:)                       ! distance between vertices
    real(r8)        , pointer     :: areaCell(:)                     ! area of grid cells [m^2]
    integer         , intent(out) :: nCells_loc                      ! number of local cell-to-cell connections
    integer         , intent(out) :: maxEdges                        ! max number of edges/neighbors
    integer         , intent(out) :: nEdges_loc                      ! number of edge length saved locally
    !
    ! !LOCAL VARIABLES:
    integer                      :: dimid,varid                      ! netCDF id's
    integer                      :: i                                ! index
    integer                      :: ier                              ! error status
    integer                      :: nCells                           ! global number of cell-to-cell connections
    integer                      :: nEdges                           ! global number of edges
    integer                      :: ibeg_c, iend_c                   ! beginning/ending index of data
    integer                      :: ibeg_e, iend_e                   ! beginning/ending index of data
    integer                      :: remainder                        ! temporary variable
    type(file_desc_t)            :: ncid                             ! netcdf id
    character(len=256)           :: varname                          ! variable name
    character(len=256)           :: locfn                            ! local file name
    logical                      :: readvar                          ! read variable in or not
    logical                      :: readdim                          ! read dimension present or not
    integer , allocatable        :: idata2d(:,:)                     ! temporary data
    integer , allocatable        :: idata1d(:)                       ! temporary data
    real(r8), allocatable        :: rdata1d(:)                       ! temporary data
    character(len=32)            :: subname = 'surfrd_get_grid_conn' ! subroutine name

    !-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          call endrun( msg=' ERROR: filename is empty)'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Check if the dimensions are present

    call ncd_inqdid(ncid,'nCells',dimid, readdim)

    if ( .not.readdim ) then
       call endrun( msg=' ERROR: Dimension nCells missing in '//filename// &
            errMsg(__FILE__, __LINE__))
    end if
    ier = pio_inq_dimlen(ncid, dimid, nCells)

    call ncd_inqdid(ncid,'maxEdges',dimid,readdim)

    if ( .not.readdim ) then
       call endrun( msg=' ERROR: Dimension maxEdges missing in '//filename// &
            errMsg(__FILE__, __LINE__))
    end if
    ier = pio_inq_dimlen(ncid, dimid, maxEdges)

    call ncd_inqdid(ncid,'nEdges',dimid,readdim)
    if ( .not.readdim ) then
       call endrun( msg=' ERROR: Dimension nEdges missing in '//filename// &
            errMsg(__FILE__, __LINE__))
    end if
    ier = pio_inq_dimlen(ncid, dimid, nEdges)

    ! Determine the size of local array that needs to be saved.
    nCells_loc = nCells/npes
    remainder  = nCells - nCells_loc*npes
    if (iam < remainder) nCells_loc = nCells_loc + 1

    nEdges_loc = nEdges/npes
    remainder  = nEdges - nEdges_loc*npes
    if (iam < remainder) nEdges_loc = nEdges_loc + 1

    ! Determine the beginning and ending index of the data to
    ! be saved
    ibeg_c = 0
    iend_c = 0
    call MPI_Scan(nCells_loc, ibeg_c, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    call MPI_Scan(  nCells_loc, iend_c, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    ibeg_c = ibeg_c + 1 - nCells_loc

    ibeg_e = 0
    iend_e = 0
    call MPI_Scan(nEdges_loc, ibeg_e, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    call MPI_Scan(  nEdges_loc, iend_e, 1, MPI_INTEGER, MPI_SUM, mpicom, ier)
    ibeg_e = ibeg_e + 1 - nEdges_loc

    ! Allocate memory
    allocate(cellsOnCell   (maxEdges, nCells_loc))
    allocate(edgesOnCell   (maxEdges, nCells_loc))
    allocate(nEdgesOnCell  (nCells_loc          ))
    allocate(areaCell      (nCells_loc          ))
    allocate(dcEdge        (nEdges_loc          ))
    allocate(dvEdge        (nEdges_loc          ))

    ! Read the data independently (i.e. each MPI-proc reads in the entire
    ! dataset)

    allocate(idata2d(maxEdges, nCells))

    ! Read cellsOnCell
    call ncd_io(ncid=ncid, varname='cellsOnCell', data=idata2d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: cellsOnCell not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    cellsOnCell(:,:) = idata2d(:,ibeg_c:iend_c)

    ! Read edgesOnCell
    call ncd_io(ncid=ncid, varname='edgesOnCell', data=idata2d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: edgesOnCell not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    edgesOnCell(:,:) = idata2d(:,ibeg_c:iend_c)

    deallocate(idata2d)

    ! Read nEdgesOnCell
    allocate(idata1d(nCells))
    call ncd_io(ncid=ncid, varname='nEdgesOnCell', data=idata1d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: areaCell not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    nEdgesOnCell(:) = idata1d(ibeg_c:iend_c)
    deallocate(idata1d)

    ! Read areaCell
    allocate(rdata1d(nCells))
    call ncd_io(ncid=ncid, varname='areaCell', data=rdata1d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: areaCell not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    areaCell(:) = rdata1d(ibeg_c:iend_c)
    deallocate(rdata1d)

    ! Read dcEdge
    allocate(rdata1d(nEdges))
    call ncd_io(ncid=ncid, varname='dcEdge', data=rdata1d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: dcEdge not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    dcEdge(:) = rdata1d(ibeg_e:iend_e)

    ! Read dvEdge
    call ncd_io(ncid=ncid, varname='dvEdge', data=rdata1d, flag='read', readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: dvEdge not found in the file'//errMsg(__FILE__, __LINE__))
    end if
    dvEdge(:) = rdata1d(ibeg_e:iend_e)

    deallocate(rdata1d)

    ! Perform cleanup
    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_grid_conn
  
  !-----------------------------------------------------------------------------------------------------
  subroutine surfrd_topounit_data(begg, endg, lfsurdat)
    !
    ! !DESCRIPTION:
    ! Read topounit surface properties data
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t, var_desc_t, ncd_pio_openfile, ncd_pio_closefile
    use ncdio_pio       , only : ncd_io, check_var, ncd_inqfdims, check_dim, ncd_inqdid, ncd_inqdlen
    use elm_varctl      , only: fsurdat
    use fileutils       , only : getfil   
	use GridcellType    , only : grc_pp
    use elm_varsur      , only : wt_tunit, elv_tunit, slp_tunit, asp_tunit,num_tunit_per_grd
    use topounit_varcon ,  only : max_topounits, has_topounit
    
    !
    ! !ARGUMENTS:
    integer          , intent(in)    :: begg, endg 
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename

    !
    ! !LOCAL VARIABLES:
    type(var_desc_t)  :: vardesc
    integer  :: n,t                ! indices
    character(len=256):: locfn                ! local file name
  
    logical  :: readvar
    integer :: dimid
    type(file_desc_t)     :: ncid         ! netcdf id
   
    real(r8),pointer :: maxTopoElv(:)            ! Maximum topounit elevation
    integer ,pointer :: numTopoPerGrid(:)        ! Number of topounits per grid
    real(r8),pointer :: TopounitFracArea(:,:)    ! Topounit fractional area
    real(r8) ,pointer :: TopounitElv(:,:)         ! Topounit elevation
    real(r8),pointer :: TopounitSlope(:,:)       ! Topounit slope 
    integer ,pointer :: TopounitAspect(:,:)      ! Topounit aspect
    integer ,pointer :: num_topo_per_grid(:)      ! Topounit aspect
    real(r8),pointer :: GridElevation(:)      ! Topounit aspect
!    integer ,pointer :: TopounitIndices(:,:)     ! Topounit indices in each grid
	
    character(len=32) :: subname = 'surfrd_topounit_data'  ! subroutine name
!-----------------------------------------------------------------------  

    allocate(maxTopoElv(begg:endg))
    allocate(GridElevation(begg:endg))
    allocate(numTopoPerGrid(begg:endg))
    allocate(TopounitFracArea(begg:endg,max_topounits))
    allocate(TopounitElv(begg:endg,max_topounits))
    allocate(TopounitSlope(begg:endg,max_topounits))
    allocate(TopounitAspect(begg:endg,max_topounits))
    allocate(num_topo_per_grid(begg:endg))
!    allocate(TopounitIndices(begg:endg,max_topounits))
    
    ! Read surface data
    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
	
    !call check_dim(ncid, 'nlevsoi', nlevsoifl)
    call check_var(ncid=ncid, varname='MaxTopounitElv', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='MaxTopounitElv', flag='read', data=maxTopoElv, &
         dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='topoPerGrid', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='topoPerGrid', flag='read', data=numTopoPerGrid, &
         dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='TopounitFracArea', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='TopounitFracArea', flag='read', data=TopounitFracArea, &
         dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='TopounitAveElv', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='TopounitAveElv', flag='read', data=TopounitElv, &
         dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='TopounitSlope', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='TopounitSlope', flag='read', data=TopounitSlope, &
         dim1name=grlnd, readvar=readvar)
    endif

    call check_var(ncid=ncid, varname='TopounitAspect', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='TopounitAspect', flag='read', data=TopounitAspect, &
         dim1name=grlnd, readvar=readvar)
    endif
    
    call check_var(ncid=ncid, varname='topoPerGrid', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='topoPerGrid', flag='read', data=num_topo_per_grid, &
         dim1name=grlnd, readvar=readvar)
    endif
    
    call check_var(ncid=ncid, varname='TOPO2', vardesc=vardesc, readvar=readvar)
    if (readvar) then
       call ncd_io(ncid=ncid, varname='TOPO2', flag='read', data=GridElevation, &
         dim1name=grlnd, readvar=readvar)
    endif
    if (readvar) then
        do n = begg,endg          
           grc_pp%MaxElevation(n) = maxTopoElv(n) 	
           grc_pp%elevation(n) = GridElevation(n) 
           num_tunit_per_grd(n) = num_topo_per_grid(n)
           grc_pp%ntopounits(n) = numTopoPerGrid(n)
           do t = 1, max_topounits
              wt_tunit(n,t) = TopounitFracArea(n,t)
              elv_tunit(n,t) = TopounitElv(n,t)
        !      slp_tunit(n,t) = TopounitSlope(n,t)
        !      asp_tunit(n,t) = TopounitAspect(n,t)              
           end do
        end do		
     endif	
    deallocate(maxTopoElv,TopounitFracArea,TopounitElv,TopounitSlope,TopounitAspect,GridElevation)
    
    call ncd_pio_closefile(ncid)
    
  end subroutine surfrd_topounit_data


!-----------------------------------------------------------------------
  subroutine surfrd_get_topo_for_solar_rad(domain,filename)
! !DESCRIPTION:
! Read the topography parameters for TOP solar radiation parameterization:
! Assume domain has already been initialized and read

! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil

! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Dalei Hao
!
! !LOCAL VARIABLES:
!EOP
    type(file_desc_t)   :: ncid             ! netcdf file id
    integer             :: n                ! indices
    integer             :: ni,nj,ns         ! size of grid on file
    integer             :: dimid,varid      ! netCDF id's
    integer             :: ier              ! error status
    real(r8)            :: eps = 1.0e-12_r8 ! lat/lon error tolerance
    integer             :: beg,end          ! local beg,end indices
    logical             :: isgrid2d         ! true => file is 2d lat/lon
    real(r8),pointer    :: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn            ! local file name
    logical             :: readvar          ! is variable on file
    character(len=32)   :: subname = 'surfrd_get_topo_for_solar_rad'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       else
          write(iulog,*) 'Attempting to read topography parameters from fsurdat ',trim(filename)
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

    if (domain%ns /= ns) then
       write(iulog,*) trim(subname),' ERROR: fsurdat file mismatch ns',&
            domain%ns,ns
       call endrun()
    endif
    
    beg = domain%nbeg
    end = domain%nend

    allocate(latc(beg:end),lonc(beg:end))

    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on fsurdat file' )

    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LATIXY  NOT on fsurdat file' )

    do n = beg,end
       if (abs(latc(n)-domain%latc(n)) > eps .or. &
           abs(lonc(n)-domain%lonc(n)) > eps) then
          write(iulog,*) trim(subname),' ERROR: fsurdat file mismatch lat,lon',latc(n),&
               domain%latc(n),lonc(n),domain%lonc(n),eps
          call endrun()
       endif
    enddo

    call ncd_io(ncid=ncid, varname='STDEV_ELEV', flag='read', data=domain%stdev_elev, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
         write(iulog,*) trim(subname),' WARNING: STDEV_ELEV  NOT on fsurdat file. Try to use STD_ELEV instead.'
         call ncd_io(ncid=ncid, varname='STD_ELEV', flag='read', data=domain%stdev_elev, &
              dim1name=grlnd, readvar=readvar)
         if (.not. readvar) call endrun( trim(subname)//' ERROR: BOTH STD_ELEV and STDEV_ELEV NOT on fsurdat file' )
    endif
    call ncd_io(ncid=ncid, varname='SKY_VIEW', flag='read', data=domain%sky_view, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: SKY_VIEW  NOT on fsurdat file' )
    call ncd_io(ncid=ncid, varname='TERRAIN_CONFIG', flag='read', data=domain%terrain_config, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: TERRAIN_CONFIG  NOT on fsurdat file' )
    call ncd_io(ncid=ncid, varname='SINSL_COSAS', flag='read', data=domain%sinsl_cosas, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: SINSL_COSAS  NOT on fsurdat file' )
    call ncd_io(ncid=ncid, varname='SINSL_SINAS', flag='read', data=domain%sinsl_sinas, &
         dim1name=grlnd, readvar=readvar)
    If (.not. readvar) call endrun( trim(subname)//' ERROR: SINSL_SINAS  NOT on fsurdat file' )

    deallocate(latc,lonc)

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_topo_for_solar_rad


  subroutine surfrd_fates_nocropmod( ncid, begg, endg )

    !--------------------------------------------------------------------------
    !     This routine evaluates the natural and crop functional
    !     type fractions in the surface file and returns them to
    !     a single, concatenated vector.  These weights
    !     are only used for a satellite phenology run.
    !     Note that FATES will actually allocate a different number of patches
    !     and will use a mapping table to connect its own pft and cft
    !     definitions to those it finds in the surface file.
    !--------------------------------------------------------------------------
    
    ! !USES:
    use elm_varsur      , only : wt_nat_patch, wt_lunit
    use elm_varpar      , only : cft_size, surfpft_lb, surfpft_ub
    use landunit_varcon , only : istsoil, istcrop
    use topounit_varcon , only : max_topounits
    
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid         ! netcdf id
    integer          , intent(in)    :: begg, endg
    
    !
    ! !LOCAL VARIABLES:
    logical  :: readvar                              ! is variable on dataset
    real(r8),pointer :: array3d_pft(:,:,:)           ! local array
    real(r8),pointer :: array3d_cft(:,:,:)           ! local array
    integer :: g,p,t
    integer :: cft_dimlen,surfpft_dimlen,dimid
    
    character(len=32) :: subname = 'surfrd_fates_nocropmod'! subroutine name
    
    call ncd_inqdlen(ncid, dimid, cft_dimlen, 'cft')
    call ncd_inqdlen(ncid, dimid, surfpft_dimlen, 'natpft')

    ! double check that cft_dimlen+surfpft_dimlen = surfpft_size
    if((cft_dimlen+surfpft_dimlen).ne.(surfpft_ub-surfpft_lb+1))then
       call endrun( msg=' ERROR: PCT+CFT dimlen does not match array size for wt_nat_patch when fates is on'&
            //errMsg(__FILE__, __LINE__))
    end if
    
    allocate( array3d_cft(begg:endg,1:max_topounits,1:cft_dimlen) )
    allocate( array3d_pft(begg:endg,1:max_topounits,1:surfpft_dimlen) )
    
    call ncd_io(ncid=ncid, varname='PCT_CFT', flag='read', data=array3d_cft, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_CFT NOT on surfdata file'//errMsg(__FILE__, __LINE__))
    
    call ncd_io(ncid=ncid, varname='PCT_NAT_PFT', flag='read', data=array3d_pft, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( msg=' ERROR: PCT_NAT_PFT NOT on surfdata file'//errMsg(__FILE__, __LINE__))

    ! In fates, all the weights in both the cft and pfts go into this array
    ! It is only used by SP mode, and it can choose what PFTs to align with

    wt_nat_patch(begg:,:,0:surfpft_dimlen-1) = array3d_pft(begg:,:,:)
    wt_nat_patch(begg:,:,surfpft_dimlen:surfpft_dimlen+cft_dimlen-1) = array3d_cft(begg:,:,:)
    
    do g = begg, endg
       do t = 1, max_topounits
          if ( wt_lunit(g,t,istcrop) > 0.0_r8 )then
             ! Move CFT over to PFT and do weighted average of the crop and soil parts
             wt_nat_patch(g,t,0:surfpft_dimlen-1) = wt_nat_patch(g,t,0:surfpft_dimlen-1) * wt_lunit(g,t,istsoil)
             wt_nat_patch(g,t,surfpft_dimlen:surfpft_dimlen+cft_dimlen-1)       = &
                  wt_nat_patch(g,t,surfpft_dimlen:surfpft_dimlen+cft_dimlen-1) * wt_lunit(g,t,istcrop)
             wt_lunit(g,t,istsoil) = (wt_lunit(g,t,istsoil) + wt_lunit(g,t,istcrop)) ! Add crop landunit to soil landunit
             wt_nat_patch(g,t,:)   =  wt_nat_patch(g,t,:) / wt_lunit(g,t,istsoil)
             wt_lunit(g,t,istcrop) = 0.0_r8                ! Zero out crop CFT's
          else
             wt_nat_patch(g,t,surfpft_dimlen:surfpft_dimlen+cft_dimlen-1) = 0.0_r8    ! Make sure generic crops are zeroed out
          end if
       end do
    end do
    
    deallocate(array3d_cft,array3d_pft)
    
    
  end subroutine surfrd_fates_nocropmod
  
end module surfrdMod
