module surfrdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: surfrdMod
!
! !DESCRIPTION:
! Contains methods for reading in surface data file and determining
! two-dimensional subgrid weights as well as writing out new surface
! dataset. When reading in the surface dataset, determines array
! which sets the PFT for each of the [maxpatch] patches and
! array which sets the relative abundance of the PFT.
! Also fills in the PFTs for vegetated portion of each grid cell.
! Fractional areas for these points pertain to "vegetated"
! area not to total grid area. Need to adjust them for fraction of grid
! that is vegetated. Also fills in urban, lake, wetland, and glacier patches.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : nlevsoifl, numpft, &
                           maxpatch_pft, numcft, maxpatch, &
                           npatch_urban_tbd, npatch_urban_hd, npatch_urban_md, &
                           numurbl, npatch_lake, npatch_wet, &
                           npatch_glacier,maxpatch_urb, npatch_glacier_mec, &
                           maxpatch_glcmec
  use clm_varctl  , only : glc_topomax, iulog, scmlat, scmlon, single_column, &
                           create_glacier_mec_landunit
  use clm_varsur  , only : wtxy, vegxy, topoxy, pctspec
  use decompMod   , only : get_proc_bounds
  use clmtype
  use spmdMod                         
  use ncdio_pio
  use pio
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
!
! !PUBLIC DATA MEMBERS:
  logical, public :: crop_prog = .false. ! If prognostic crops is turned on
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by T Craig
! Updated by Mariana Vertenstein Jan 2011
!
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_wtxy_special
  private :: surfrd_wtxy_veg_all
  private :: surfrd_wtxy_veg_dgvm
!
! !PRIVATE DATA MEMBERS:
  ! default multiplication factor for epsilon for error checks
  real(r8), private, parameter :: eps_fact = 2._r8
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_globmask
!
! !INTERFACE:
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
    implicit none
    character(len=*), intent(in)    :: filename  ! grid filename
    integer         , pointer       :: mask(:)   ! grid mask 
    integer         , intent(out)   :: ni, nj    ! global grid sizes
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
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
          call endrun()
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
    if (.not. readvar) call endrun( trim(subname)//' ERROR: landmask not on fatmlndfrc file' )

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_globmask

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_grid
!
! !INTERFACE:
  subroutine surfrd_get_grid(ldomain, filename, glcfilename)
!
! !DESCRIPTION:
! THIS IS CALLED AFTER THE DOMAIN DECOMPOSITION HAS BEEN CREATED
! Read the surface dataset grid related information:
! o real latitude  of grid cell (degrees)
! o real longitude of grid cell (degrees)
!
! !USES:
    use clm_varcon, only : spval, re
    use domainMod , only : domain_type, domain_init, domain_clean, lon1d, lat1d
    use decompMod , only : get_proc_bounds
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: ldomain   ! domain to init
    character(len=*) ,intent(in)    :: filename  ! grid filename
    character(len=*) ,optional, intent(in) :: glcfilename ! glc mask filename
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
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
    character(len=16) :: vname              ! temporary
    character(len=256):: locfn              ! local file name
    integer :: n                            ! indices
    real(r8):: eps = 1.0e-12_r8             ! lat/lon error tolerance
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       if (filename == ' ') then
          write(iulog,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif
    end if

    call getfil( filename, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Determine dimensions

    call ncd_inqfdims(ncid, isgrid2d, ni, nj, ns)

    ! Determine isgrid2d flag for domain

    call get_proc_bounds(beg, end)
    call domain_init(ldomain, isgrid2d=isgrid2d, ni=ni, nj=nj, nbeg=beg, nend=end)
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
       if (.not. readvar) call endrun( trim(subname)//' ERROR: area NOT on file' )
       
       call ncd_io(ncid=ncid, varname= 'xc', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: xc NOT on file' )
       
       call ncd_io(ncid=ncid, varname= 'yc', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: yc NOT on file' )
    else
       call ncd_io(ncid=ncid, varname= 'AREA', flag='read', data=ldomain%area, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: AREA NOT on file' )
       
       call ncd_io(ncid=ncid, varname= 'LONGXY', flag='read', data=ldomain%lonc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY NOT on file' )
       
       call ncd_io(ncid=ncid, varname= 'LATIXY', flag='read', data=ldomain%latc, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: LATIXY NOT on file' )
    end if

    if (isgrid2d) then
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
    end if

    ! Check lat limited to -90,90

    if (minval(ldomain%latc) < -90.0_r8 .or. &
        maxval(ldomain%latc) >  90.0_r8) then
       write(iulog,*) trim(subname),' WARNING: lat/lon min/max is ', &
            minval(ldomain%latc),maxval(ldomain%latc)
       ! call endrun( trim(subname)//' ERROR: lat is outside [-90,90]' )
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
          call endrun( trim(subname)//' ERROR: LANDMASK NOT on fracdata file' )
       end if
    end if

    call ncd_io(ncid=ncid, varname='LANDFRAC', flag='read', data=ldomain%frac, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call ncd_io(ncid=ncid, varname='frac', flag='read', data=ldomain%frac, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) then
          call endrun( trim(subname)//' ERROR: LANDFRAC NOT on fracdata file' )
       end if
    end if

    call ncd_pio_closefile(ncid)

    if (present(glcfilename)) then
       if (masterproc) then
          if (glcfilename == ' ') then
             write(iulog,*) trim(subname),' ERROR: glc filename must be specified '
             call endrun()
          endif
       end if
       call getfil( glcfilename, locfn, 0 )
       call ncd_pio_openfile (ncidg, trim(locfn), 0)

       ldomain%glcmask(:) = 0
       call ncd_io(ncid=ncidg, varname='GLCMASK', flag='read', data=ldomain%glcmask, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: GLCMASK NOT in file' )

       ! Make sure the glc mask is a subset of the land mask
       do n = beg,end
          if (ldomain%glcmask(n)==1 .and. ldomain%mask(n)==0) then
             write(iulog,*)trim(subname),&
                  'initialize1: landmask/glcmask mismatch'
             write(iulog,*)trim(subname),&
                  'glc requires input where landmask = 0, gridcell index', n
             call endrun()
          endif
       enddo
       call ncd_pio_closefile(ncidg)
    endif   ! present(glcfilename)

  end subroutine surfrd_get_grid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_topo
!
! !INTERFACE:
  subroutine surfrd_get_topo(domain,filename)
!
! !DESCRIPTION:
! Read the topo dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
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
          call endrun()
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
       call endrun()
    endif
    
    beg = domain%nbeg
    end = domain%nend

    allocate(latc(beg:end),lonc(beg:end))

    call ncd_io(ncid=ncid, varname='LONGXY', flag='read', data=lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    call ncd_io(ncid=ncid, varname='LATIXY', flag='read', data=latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    do n = beg,end
       if (abs(latc(n)-domain%latc(n)) > eps .or. &
           abs(lonc(n)-domain%lonc(n)) > eps) then
          write(iulog,*) trim(subname),' ERROR: topo file mismatch lat,lon',latc(n),&
               domain%latc(n),lonc(n),domain%lonc(n),eps
          call endrun()
       endif
    enddo

    call ncd_io(ncid=ncid, varname='TOPO', flag='read', data=domain%topo, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: LONGXY  NOT on topodata file' )

    deallocate(latc,lonc)

    call ncd_pio_closefile(ncid)

  end subroutine surfrd_get_topo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_data
!
! !INTERFACE:
  subroutine surfrd_get_data (ldomain, lfsurdat)
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
    use clm_varctl  , only : allocate_all_vegpfts, create_crop_landunit
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
    use domainMod   , only : domain_type, domain_init, domain_clean
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in) :: ldomain     ! land domain associated with wtxy
    character(len=*), intent(in) :: lfsurdat    ! surface dataset filename
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    type(var_desc_t)  :: vardesc              ! pio variable descriptor
    type(domain_type) :: surfdata_domain      ! local domain associated with surface dataset
    character(len=256):: locfn                ! local file name
    integer           :: n                    ! loop indices
    integer           :: ni,nj,ns             ! domain sizes
    character(len=16) :: lon_var, lat_var     ! names of lat/lon on dataset
    logical           :: readvar              ! true => variable is on dataset
    real(r8)          :: rmaxlon,rmaxlat      ! local min/max vars
    type(file_desc_t) :: ncid                 ! netcdf id
    integer           :: begg,endg            ! beg,end gridcell indices
    logical           :: istype_domain        ! true => input file is of type domain
    logical           :: isgrid2d             ! true => intut grid is 2d 
    character(len=32) :: subname = 'surfrd_get_data'    ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          write(iulog,*)'lfsurdat must be specified'; call endrun()
       endif
    endif

    call get_proc_bounds(begg,endg)
    allocate(pctspec(begg:endg))

    vegxy(:,:) = noveg
    wtxy(:,:)  = 0._r8
    pctspec(:) = 0._r8
    if (allocated(topoxy)) topoxy(:,:) = 0._r8

    ! Read surface data

    call getfil( lfsurdat, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Read in pft mask - this variable is only on the surface dataset - but not
    ! on the domain dataset

    call ncd_io(ncid=ncid, varname= 'PFTDATA_MASK', flag='read', data=ldomain%pftm, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: pftm NOT on surface dataset' )

    ! Check if fsurdat grid is "close" to fatmlndfrc grid, exit if lats/lon > 0.001

    call check_var(ncid=ncid, varname='xc', vardesc=vardesc, readvar=readvar) 
    if (readvar) then
       istype_domain = .true.
    else
       call check_var(ncid=ncid, varname='LONGXY', vardesc=vardesc, readvar=readvar) 
       if (readvar) then
          istype_domain = .false.
       else
          call endrun( trim(subname)//' ERROR: unknown domain type')
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
    call domain_init(surfdata_domain, isgrid2d, ni, nj, begg, endg, clmlevel=grlnd)

    call ncd_io(ncid=ncid, varname=lon_var, flag='read', data=surfdata_domain%lonc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: lon var NOT on surface dataset' )

    call ncd_io(ncid=ncid, varname=lat_var, flag='read', data=surfdata_domain%latc, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: lat var NOT on surface dataset' )

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
       write(iulog,*) trim(subname)//': surfdata/fatmgrid lon/lat mismatch error',&
            rmaxlon,rmaxlat
       call endrun(trim(subname))
    end if
    call domain_clean(surfdata_domain)

    ! Obtain special landunit info

    call surfrd_wtxy_special(ncid, ldomain%ns)

    ! Obtain vegetated landunit info

#if (defined CNDV)
    if (create_crop_landunit) then ! CNDV means allocate_all_vegpfts = .true.
       call surfrd_wtxy_veg_all(ncid, ldomain%ns, ldomain%pftm)
    end if
    call surfrd_wtxy_veg_dgvm()
#else
    if (allocate_all_vegpfts) then
       call surfrd_wtxy_veg_all(ncid, ldomain%ns, ldomain%pftm)
    else
       call endrun (trim(subname) // 'only allocate_all_vegpfts is supported')
    end if
#endif

    call ncd_pio_closefile(ncid)

    if ( masterproc )then
       write(iulog,*) 'Successfully read surface boundary data'
       write(iulog,*)
    end if

  end subroutine surfrd_get_data

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_special 
!
! !INTERFACE:
  subroutine surfrd_wtxy_special(ncid, ns)
!
! !DESCRIPTION:
! Determine weight with respect to gridcell of all special "pfts" as well
! as soil color and percent sand and clay
!
! !USES:
    use pftvarcon     , only : noveg
    use UrbanInputMod , only : urbinp
    use clm_varpar    , only : maxpatch_glcmec, nlevurb
    use clm_varcon    , only : udens_base, udens_tbd, udens_hd, udens_md 
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    integer          , intent(in)    :: ns     ! domain size
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: n,nl,nurb,g                ! indices
    integer  :: begg,endg                  ! gcell beg/end
    integer  :: dimid,varid                ! netCDF id's
    real(r8) :: nlevsoidata(nlevsoifl)
    logical  :: found                      ! temporary for error check
    integer  :: nindx                      ! temporary for error check
    integer  :: ier                        ! error status
    integer  :: nlev                       ! level
    integer  :: npatch
    logical  :: readvar
    real(r8),pointer :: pctgla(:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb(:,:)      ! percent of grid cell is urbanized
    real(r8),pointer :: pctglc_mec(:,:)   ! percent of grid cell is glacier_mec (in each elev class)
    real(r8),pointer :: pctglc_mec_tot(:) ! percent of grid cell is glacier (sum over classes)
    real(r8),pointer :: pcturb_tot(:)  ! percent of grid cell is urban (sum over density classes)
    integer  :: dindx                      ! temporary for error check
    real(r8),pointer :: topoglc_mec(:,:)  ! surface elevation in each elev class
    character(len=32) :: subname = 'surfrd_wtxy_special'  ! subroutine name
    real(r8) closelat,closelon
!!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg,numurbl),pcturb_tot(begg:endg))
    if (create_glacier_mec_landunit) then
       allocate(pctglc_mec(begg:endg,maxpatch_glcmec))
       allocate(pctglc_mec_tot(begg:endg))
       allocate(topoglc_mec(begg:endg,maxpatch_glcmec))
       allocate(glc_topomax(0:maxpatch_glcmec))
    endif

    call check_dim(ncid, 'nlevsoi', nlevsoifl)

       ! Obtain non-grid surface properties of surface dataset other than percent pft

    call ncd_io(ncid=ncid, varname='PCT_WETLAND', flag='read', data=pctwet, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_WETLAND  NOT on surfdata file' )

    call ncd_io(ncid=ncid, varname='PCT_LAKE'   , flag='read', data=pctlak, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_LAKE NOT on surfdata file' )

    call ncd_io(ncid=ncid, varname='PCT_GLACIER', flag='read', data=pctgla, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLACIER NOT on surfdata file' )

    ! If PCT_URBAN is not multi-density then set pcturb and nlevurb to zero 
    if (nlevurb == 0) then
      pcturb = 0._r8
      write(iulog,*)'PCT_URBAN is not multi-density, pcturb set to 0'
    else
      call ncd_io(ncid=ncid, varname='PCT_URBAN'  , flag='read', data=pcturb, &
           dim1name=grlnd, readvar=readvar)
      if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_URBAN NOT on surfdata file' )
    end if
    if ( nlevurb == 0 )then
       if ( any(pcturb > 0.0_r8) ) call endrun( trim(subname)//' ERROR: PCT_URBAN MUST be zero when nlevurb=0' )
    end if

    pcturb_tot(:) = 0._r8
    do n = 1, numurbl
       do nl = begg,endg
          pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
       enddo
    enddo

    if (create_glacier_mec_landunit) then          ! call ncd_io_gs_int2d

       call check_dim(ncid, 'nglcec',   maxpatch_glcmec   )
       call check_dim(ncid, 'nglcecp1', maxpatch_glcmec+1 )

       call ncd_io(ncid=ncid, varname='GLC_MEC', flag='read', data=glc_topomax, &
            readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//'ERROR: GLC_MEC was NOT on the input surfdata file' )

       call ncd_io(ncid=ncid, varname='PCT_GLC_MEC', flag='read', data=pctglc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_GLC_MEC NOT on surfdata file' )

       call ncd_io(ncid=ncid, varname='TOPO_GLC_MEC',  flag='read', data=topoglc_mec, &
            dim1name=grlnd, readvar=readvar)
       if (.not. readvar) call endrun( trim(subname)//' ERROR: TOPO_GLC_MEC NOT on surfdata file' )

       pctglc_mec_tot(:) = 0._r8
       do n = 1, maxpatch_glcmec
          do nl = begg,endg
             pctglc_mec_tot(nl) = pctglc_mec_tot(nl) + pctglc_mec(nl,n)
          enddo
       enddo

       ! Make sure sum of pctglc_mec = pctgla (approximately), then correct pctglc_mec so
       ! that its sum more exactly equals pctgla, then zero out pctgla
       !
       ! WJS (9-28-12): The reason for the correction piece of this is: in the surface
       ! dataset, pctgla underwent various minor corrections that make it the trusted
       ! value (as opposed to sum(pctglc_mec). sum(pctglc_mec) can differ from pctgla by
       ! single precision roundoff. This difference can cause problems later (e.g., in the
       ! consistency check in pftdynMod), so we do this correction here. It might be
       ! better to do this correction in mksurfdata_map, but because of time constraints,
       ! which make me unable to recreate surface datasets, I have to do it here instead
       ! (and there are arguments for putting it here anyway).

       do nl = begg,endg
          ! We need to use a fairly high threshold for equality (2.0e-5) because pctgla
          ! and pctglc_mec are computed using single precision inputs. Note that this
          ! threshold agrees with the threshold in the error checks in mkglcmecMod:
          ! mkglcmec in mksurfdata_map. 
          if (abs(pctgla(nl) - pctglc_mec_tot(nl)) > 2.0e-5) then
             write(iulog,*) ' '
             write(iulog,*) 'surfrd error: pctgla not equal to sum of pctglc_mec for nl=', nl
             write(iulog,*) 'pctgla =', pctgla(nl)
             write(iulog,*) 'pctglc_mec_tot =', pctglc_mec_tot(nl)
             call endrun()
          endif

          ! Correct the remaining minor differences in pctglc_mec so sum more exactly
          ! equals pctglc (see notes above for justification)
          if (pctglc_mec_tot(nl) > 0.0_r8) then
             pctglc_mec(nl,:) = pctglc_mec(nl,:) * pctgla(nl)/pctglc_mec_tot(nl)
          end if

          ! Now recompute pctglc_mec_tot, and confirm that we are now much closer to pctgla
          pctglc_mec_tot(nl) = 0._r8
          do n = 1, maxpatch_glcmec
             pctglc_mec_tot(nl) = pctglc_mec_tot(nl) + pctglc_mec(nl,n)
          end do

          if (abs(pctgla(nl) - pctglc_mec_tot(nl)) > 1.0e-13) then
             write(iulog,*) ' '
             write(iulog,*) 'surfrd error: after correction, pctgla not equal to sum of pctglc_mec for nl=', nl
             write(iulog,*) 'pctgla =', pctgla(nl)
             write(iulog,*) 'pctglc_mec_tot =', pctglc_mec_tot(nl)
             call endrun()
          endif

          ! Finally, zero out pctgla
          pctgla(nl) = 0._r8
       enddo

       ! If pctglc_mec_tot is very close to 100%, round to 100%

       do nl = begg,endg
          ! The inequality here ( <= 2.0e-5 ) is designed to be the complement of the
          ! above check that makes sure pctglc_mec_tot is close to pctgla: so if pctglc=
          ! 100 (exactly), then exactly one of these conditionals will be triggered.
          ! Update 9-28-12: Now that there is a rescaling of pctglc_mec to bring it more
          ! in line with pctgla, we could probably decrease this tolerance again (the
          ! point about exactly one of these conditionals being triggered no longer holds)
          ! - or perhaps even get rid of this whole block of code. But I'm keeping this as
          ! is for now because that's how I tested it, and I don't think it will hurt
          ! anything to use this larger tolerance.
          if (abs(pctglc_mec_tot(nl) - 100._r8) <= 2.0e-5) then
             pctglc_mec(nl,:) = pctglc_mec(nl,:) * 100._r8 / pctglc_mec_tot(nl)
             pctglc_mec_tot(nl) = 100._r8
          endif
       enddo

       pctspec = pctwet + pctlak + pcturb_tot + pctglc_mec_tot

       if ( masterproc ) write(iulog,*) '   elevation limits = ', glc_topomax

    else

       pctspec = pctwet + pctlak + pcturb_tot + pctgla
 
    endif

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg,endg
       if (pctspec(nl) > 100._r8+1.e-04_r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: PFT cover>100 for nl=',nindx
       call endrun()
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg,endg

       vegxy(nl,npatch_lake)   = noveg
       wtxy(nl,npatch_lake)    = pctlak(nl)/100._r8

       vegxy(nl,npatch_wet)    = noveg
       wtxy(nl,npatch_wet)     = pctwet(nl)/100._r8

       vegxy(nl,npatch_glacier)= noveg
       wtxy(nl,npatch_glacier) = pctgla(nl)/100._r8

       ! Initialize urban tall building district weights
       n = udens_tbd - udens_base
       do nurb = npatch_urban_tbd, npatch_urban_hd-1 
          vegxy(nl,nurb) = noveg
          wtxy(nl,nurb)  = pcturb(nl,n) / 100._r8
       end do
       if ( pcturb(nl,n) > 0.0_r8 )then
          wtxy(nl,npatch_urban_tbd)   = wtxy(nl,npatch_urban_tbd)*urbinp%wtlunit_roof(nl,n)
          wtxy(nl,npatch_urban_tbd+1) = wtxy(nl,npatch_urban_tbd+1)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_tbd+2) = wtxy(nl,npatch_urban_tbd+2)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_tbd+3) = wtxy(nl,npatch_urban_tbd+3)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * (1.-urbinp%wtroad_perv(nl,n))
          wtxy(nl,npatch_urban_tbd+4) = wtxy(nl,npatch_urban_tbd+4)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * urbinp%wtroad_perv(nl,n)
       end if

       ! Initialize urban high density weights
       n = udens_hd - udens_base
       do nurb = npatch_urban_hd, npatch_urban_md-1 
          vegxy(nl,nurb) = noveg
          wtxy(nl,nurb)  = pcturb(nl,n) / 100._r8
       end do
       if ( pcturb(nl,n) > 0.0_r8 )then
          wtxy(nl,npatch_urban_hd)   = wtxy(nl,npatch_urban_hd)*urbinp%wtlunit_roof(nl,n)
          wtxy(nl,npatch_urban_hd+1) = wtxy(nl,npatch_urban_hd+1)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_hd+2) = wtxy(nl,npatch_urban_hd+2)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_hd+3) = wtxy(nl,npatch_urban_hd+3)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * (1.-urbinp%wtroad_perv(nl,n))
          wtxy(nl,npatch_urban_hd+4) = wtxy(nl,npatch_urban_hd+4)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * urbinp%wtroad_perv(nl,n)
       end if

       ! Initialize urban medium density weights
       n = udens_md - udens_base
       do nurb = npatch_urban_md, npatch_lake-1 
          vegxy(nl,nurb) = noveg
          wtxy(nl,nurb)  = pcturb(nl,n) / 100._r8
       end do
       if ( pcturb(nl,n) > 0.0_r8 )then
          wtxy(nl,npatch_urban_md)   = wtxy(nl,npatch_urban_md)*urbinp%wtlunit_roof(nl,n)
          wtxy(nl,npatch_urban_md+1) = wtxy(nl,npatch_urban_md+1)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_md+2) = wtxy(nl,npatch_urban_md+2)*(1 - urbinp%wtlunit_roof(nl,n))/3
          wtxy(nl,npatch_urban_md+3) = wtxy(nl,npatch_urban_md+3)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * (1.-urbinp%wtroad_perv(nl,n))
          wtxy(nl,npatch_urban_md+4) = wtxy(nl,npatch_urban_md+4)*(1 - urbinp%wtlunit_roof(nl,n))/3 &
                                        * urbinp%wtroad_perv(nl,n)
       end if

    end do

    ! Check to make sure we have valid urban data for each urban patch

    found = .false.
    do nl = begg,endg
      do n = 1, numurbl
        if ( pcturb(nl,n) > 0.0_r8 ) then
          if (urbinp%canyon_hwr(nl,n)            .le. 0._r8 .or. &
              urbinp%em_improad(nl,n)            .le. 0._r8 .or. &
              urbinp%em_perroad(nl,n)            .le. 0._r8 .or. &
              urbinp%em_roof(nl,n)               .le. 0._r8 .or. &
              urbinp%em_wall(nl,n)               .le. 0._r8 .or. &
              urbinp%ht_roof(nl,n)               .le. 0._r8 .or. &
              urbinp%thick_roof(nl,n)            .le. 0._r8 .or. &
              urbinp%thick_wall(nl,n)            .le. 0._r8 .or. &
              urbinp%t_building_max(nl,n)        .le. 0._r8 .or. &
              urbinp%t_building_min(nl,n)        .le. 0._r8 .or. &
              urbinp%wind_hgt_canyon(nl,n)       .le. 0._r8 .or. &
              urbinp%wtlunit_roof(nl,n)          .le. 0._r8 .or. &
              urbinp%wtroad_perv(nl,n)           .le. 0._r8 .or. &
              any(urbinp%alb_improad_dir(nl,n,:) .le. 0._r8) .or. &
              any(urbinp%alb_improad_dif(nl,n,:) .le. 0._r8) .or. &
              any(urbinp%alb_perroad_dir(nl,n,:) .le. 0._r8) .or. &
              any(urbinp%alb_perroad_dif(nl,n,:) .le. 0._r8) .or. &
              any(urbinp%alb_roof_dir(nl,n,:)    .le. 0._r8) .or. &
              any(urbinp%alb_roof_dif(nl,n,:)    .le. 0._r8) .or. &
              any(urbinp%alb_wall_dir(nl,n,:)    .le. 0._r8) .or. &
              any(urbinp%alb_wall_dif(nl,n,:)    .le. 0._r8) .or. &
              any(urbinp%tk_roof(nl,n,:)         .le. 0._r8) .or. &
              any(urbinp%tk_wall(nl,n,:)         .le. 0._r8) .or. &
              any(urbinp%cv_roof(nl,n,:)         .le. 0._r8) .or. &
              any(urbinp%cv_wall(nl,n,:)         .le. 0._r8)) then
            found = .true.
            nindx = nl
            dindx = n
            exit
          else
            if (urbinp%nlev_improad(nl,n) .gt. 0) then
               nlev = urbinp%nlev_improad(nl,n)
               if (any(urbinp%tk_improad(nl,n,1:nlev) .le. 0._r8) .or. &
                   any(urbinp%cv_improad(nl,n,1:nlev) .le. 0._r8)) then
                  found = .true.
                  nindx = nl
                  dindx = n
                  exit
               end if
            end if
          end if
          if (found) exit
        end if
      end do
    end do
    if ( found ) then
       write(iulog,*)'surfrd error: no valid urban data for nl=',nindx
       write(iulog,*)'density type:    ',dindx
       write(iulog,*)'canyon_hwr:      ',urbinp%canyon_hwr(nindx,dindx)
       write(iulog,*)'em_improad:      ',urbinp%em_improad(nindx,dindx)
       write(iulog,*)'em_perroad:      ',urbinp%em_perroad(nindx,dindx)
       write(iulog,*)'em_roof:         ',urbinp%em_roof(nindx,dindx)
       write(iulog,*)'em_wall:         ',urbinp%em_wall(nindx,dindx)
       write(iulog,*)'ht_roof:         ',urbinp%ht_roof(nindx,dindx)
       write(iulog,*)'thick_roof:      ',urbinp%thick_roof(nindx,dindx)
       write(iulog,*)'thick_wall:      ',urbinp%thick_wall(nindx,dindx)
       write(iulog,*)'t_building_max:  ',urbinp%t_building_max(nindx,dindx)
       write(iulog,*)'t_building_min:  ',urbinp%t_building_min(nindx,dindx)
       write(iulog,*)'wind_hgt_canyon: ',urbinp%wind_hgt_canyon(nindx,dindx)
       write(iulog,*)'wtlunit_roof:    ',urbinp%wtlunit_roof(nindx,dindx)
       write(iulog,*)'wtroad_perv:     ',urbinp%wtroad_perv(nindx,dindx)
       write(iulog,*)'alb_improad_dir: ',urbinp%alb_improad_dir(nindx,dindx,:)
       write(iulog,*)'alb_improad_dif: ',urbinp%alb_improad_dif(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dir: ',urbinp%alb_perroad_dir(nindx,dindx,:)
       write(iulog,*)'alb_perroad_dif: ',urbinp%alb_perroad_dif(nindx,dindx,:)
       write(iulog,*)'alb_roof_dir:    ',urbinp%alb_roof_dir(nindx,dindx,:)
       write(iulog,*)'alb_roof_dif:    ',urbinp%alb_roof_dif(nindx,dindx,:)
       write(iulog,*)'alb_wall_dir:    ',urbinp%alb_wall_dir(nindx,dindx,:)
       write(iulog,*)'alb_wall_dif:    ',urbinp%alb_wall_dif(nindx,dindx,:)
       write(iulog,*)'tk_roof:         ',urbinp%tk_roof(nindx,dindx,:)
       write(iulog,*)'tk_wall:         ',urbinp%tk_wall(nindx,dindx,:)
       write(iulog,*)'cv_roof:         ',urbinp%cv_roof(nindx,dindx,:)
       write(iulog,*)'cv_wall:         ',urbinp%cv_wall(nindx,dindx,:)
       if (urbinp%nlev_improad(nindx,dindx) .gt. 0) then
          nlev = urbinp%nlev_improad(nindx,dindx)
          write(iulog,*)'tk_improad: ',urbinp%tk_improad(nindx,dindx,1:nlev)
          write(iulog,*)'cv_improad: ',urbinp%cv_improad(nindx,dindx,1:nlev)
       end if
       call endrun()
    end if
    ! Initialize glacier_mec weights

    if (create_glacier_mec_landunit) then
       do n = 1, maxpatch_glcmec
          npatch = npatch_glacier_mec - maxpatch_glcmec + n
          do nl = begg, endg
             vegxy (nl,npatch) = noveg
             wtxy  (nl,npatch) = pctglc_mec(nl,n)/100._r8
             topoxy(nl,npatch) = topoglc_mec(nl,n)
          enddo   ! nl
       enddo      ! maxpatch_glcmec
       deallocate(pctglc_mec, pctglc_mec_tot, topoglc_mec)
    endif          ! create_glacier_mec_landunit

    deallocate(pctgla,pctlak,pctwet,pcturb,pcturb_tot)

  end subroutine surfrd_wtxy_special

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_all
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_all(ncid, ns, pftm)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use clm_varctl  , only : create_crop_landunit, fpftdyn, irrigate
    use pftvarcon   , only : nc3crop, nc3irrig, npcropmin, &
                             ncorn, ncornirrig, nsoybean, nsoybeanirrig, &
                             nscereal, nscerealirrig, nwcereal, nwcerealirrig
    use spmdMod     , only : mpicom, MPI_LOGICAL, MPI_LOR
!
! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! netcdf id
    integer          ,intent(in)    :: ns     ! domain size
    integer          ,pointer       :: pftm(:)
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: m,mp7,mp8,mp11,n,nl            ! indices
    integer  :: begg,endg                      ! beg/end gcell index
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: ier                            ! error status	
    logical  :: readvar                        ! is variable on dataset
    real(r8) :: sumpct                         ! sum of %pft over maxpatch_pft
    real(r8),allocatable :: pctpft(:,:)        ! percent of vegetated gridcell area for PFTs
    real(r8),pointer :: arrayl(:,:)            ! local array
    real(r8) :: numpftp1data(0:numpft)         
    logical  :: crop = .false.                 ! if crop data on this section of file
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    allocate(pctpft(begg:endg,0:numpft))

    call check_dim(ncid, 'lsmpft', numpft+1)

    allocate(arrayl(begg:endg,0:numpft))
    call ncd_io(ncid=ncid, varname='PCT_PFT', flag='read', data=arrayl, &
         dim1name=grlnd, readvar=readvar)
    if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_PFT NOT on surfdata file' )
    pctpft(begg:endg,0:numpft) = arrayl(begg:endg,0:numpft)
    deallocate(arrayl)

    do nl = begg,endg
       if (pftm(nl) >= 0) then

          ! Error check: make sure PFTs sum to 100% cover for vegetated landunit 
          ! (convert pctpft from percent with respect to gridcel to percent with 
          ! respect to vegetated landunit)

          ! THESE CHECKS NEEDS TO BE THE SAME AS IN pftdynMod.F90!
          if (pctspec(nl) < 100._r8 * (1._r8 - eps_fact*epsilon(1._r8))) then  ! pctspec not within eps_fact*epsilon of 100
             sumpct = 0._r8
             do m = 0,numpft
                sumpct = sumpct + pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
             end do
             if (abs(sumpct - 100._r8) > 0.1e-4_r8) then
                write(iulog,*) trim(subname)//' ERROR: sum(pct) over numpft+1 is not = 100.'
                write(iulog,*) sumpct, sumpct-100._r8, nl
                call endrun()
             end if
             if (sumpct < -0.000001_r8) then
                write(iulog,*) trim(subname)//' ERROR: sum(pct) over numpft+1 is < 0.'
                write(iulog,*) sumpct, nl
                call endrun()
             end if
             do m = npcropmin, numpft
                if ( pctpft(nl,m) > 0.0_r8 ) crop = .true.
             end do
          end if

          ! Set weight of each pft wrt gridcell (note that maxpatch_pft = numpft+1 here)

          do m = 1,numpft+1
             vegxy(nl,m)  = m - 1 ! 0 (bare ground) to numpft
             wtxy(nl,m) = pctpft(nl,m-1) / 100._r8
          end do

       end if
    end do

    call mpi_allreduce(crop,crop_prog,1,MPI_LOGICAL,MPI_LOR,mpicom,ier)
    if (ier /= 0) then
       write(iulog,*) trim(subname)//' mpi_allreduce error = ',ier
       call endrun( trim(subname)//' ERROR: error in reduce of crop_prog' )
    endif
    if (crop_prog .and. .not. create_crop_landunit) then
       call endrun( trim(subname)//' ERROR: prognostic crop '// &
                    'PFTs require create_crop_landunit=.true.' )
    end if
    if (crop_prog .and. fpftdyn /= ' ') &
       call endrun( trim(subname)//' ERROR: prognostic crop '// &
                    'is incompatible with transient landuse' )
    if (.not. crop_prog .and. irrigate) then
       call endrun( trim(subname)//' ERROR surfrdMod: irrigate = .true. requires CROP model active.' )
    end if

    if (masterproc .and. crop_prog .and. .not. irrigate) then
       write(iulog,*) trim(subname)//' crop=.T. and irrigate=.F., so merging irrigated pfts with rainfed' ! in the following do-loop
    end if

! repeat do-loop for error checking and for rainfed crop case

    do nl = begg,endg
       if (pftm(nl) >= 0) then
          if (pctspec(nl) < 100._r8 * (1._r8 - eps_fact*epsilon(1._r8))) then  ! pctspec not within eps_fact*epsilon of 100
             if (.not. crop_prog .and. wtxy(nl,nc3irrig+1) > 0._r8) then
                call endrun( trim(subname)//' ERROR surfrdMod: irrigated crop PFT requires CROP model active.' )
             end if
             if (crop_prog .and. .not. irrigate) then
                wtxy(nl,nc3crop+1)       = wtxy(nl,nc3crop+1)  + wtxy(nl,nc3irrig+1)
                wtxy(nl,nc3irrig+1)      = 0._r8
                wtxy(nl,ncorn+1)         = wtxy(nl,ncorn+1)    + wtxy(nl,ncornirrig+1)
                wtxy(nl,ncornirrig+1)    = 0._r8
                wtxy(nl,nscereal+1)      = wtxy(nl,nscereal+1) + wtxy(nl,nscerealirrig+1)
                wtxy(nl,nscerealirrig+1) = 0._r8
                wtxy(nl,nwcereal+1)      = wtxy(nl,nwcereal+1) + wtxy(nl,nwcerealirrig+1)
                wtxy(nl,nwcerealirrig+1) = 0._r8
                wtxy(nl,nsoybean+1)      = wtxy(nl,nsoybean+1) + wtxy(nl,nsoybeanirrig+1)
                wtxy(nl,nsoybeanirrig+1) = 0._r8
             end if
          end if
       end if
    end do


    deallocate(pctpft)

  end subroutine surfrd_wtxy_veg_all

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_dgvm
!
! !INTERFACE:
  subroutine surfrd_wtxy_veg_dgvm()
!
! !DESCRIPTION:
! Determine wtxy and vegxy for CNDV mode.
!
! !USES:
    use pftvarcon , only : noveg, crop
    use clm_varctl, only : create_crop_landunit
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine surfrd in this module
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/04
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: m,nl         ! indices
    integer  :: begg,endg   ! beg/end gcell index
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    do nl = begg,endg
       do m = 1, maxpatch_pft ! CNDV means allocate_all_vegpfts = .true.
          if (create_crop_landunit) then ! been through surfrd_wtxy_veg_all
             if (crop(m-1) == 0) then    ! so update natural vegetation only
                wtxy(nl,m) = 0._r8       ! crops should have values >= 0.
             end if
          else                   ! not been through surfrd_wtxy_veg_all
             wtxy(nl,m) = 0._r8  ! so update all vegetation
             vegxy(nl,m) = m - 1 ! 0 (bare ground) to maxpatch_pft-1 (= 16)
          end if
       end do
       ! bare ground weights
       wtxy(nl,noveg+1) = max(0._r8, 1._r8 - sum(wtxy(nl,:)))
       if (abs(sum(wtxy(nl,:)) - 1._r8) > 1.e-5_r8) then
          write(iulog,*) 'all wtxy =', wtxy(nl,:)
          call endrun()
       end if ! error check
    end do

  end subroutine surfrd_wtxy_veg_dgvm
   
end module surfrdMod
