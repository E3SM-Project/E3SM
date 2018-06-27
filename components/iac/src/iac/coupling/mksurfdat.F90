module mksurfdat

    use shr_kind_mod , only : r8 => shr_kind_r8, r4 => shr_kind_r4
    use shr_sys_mod  , only : shr_sys_getenv
    use shr_file_mod , only : shr_file_getunit, shr_file_freeunit
    use shr_timer_mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mksurfdat
!
! !DESCRIPTION:
! Creates land model surface dataset from original "raw" data files.
! Surface dataset contains model grid, pfts, inland water, glacier,
! soil texture, soil color, LAI and SAI, urban fraction, and urban
! parameters.
!
!EOP

! !USES:
    use mkfileutils  , only : getfil, opnfil, get_filename, getavu
    use mklaiMod     , only : mklai
    use mkpftMod     , only : pft_idx, pft_frc, mkpft, mkpftInit, mkpft_parse_oride, &
                              mkirrig
    use mksoilMod    , only : soil_sand, soil_clay, mksoiltex, mksoiltexInit, &
                              soil_color, mksoilcol, mksoilcolInit, mkorganic
    use mkvocefMod   , only : mkvocef
    use mklanwatMod  , only : mklanwat
    use mkfmaxMod    , only : mkfmax
    use mkglcmecMod  , only : nglcec, mkglcmec, mkglcmecInit, mkglacier
    use mkharvestMod , only : mkharvest, mkharvest_init, mkharvest_fieldname, &
                              mkharvest_numtypes, mkharvest_parse_oride
    use mkurbanparMod, only : mkurban, mkurbanpar, mkelev
    use mkfileMod    , only : mkfile
    use mkvarpar     , only : numpft, nlevsoi, elev_thresh
    use mkvarctl
    use mkncdio      , only : check_ret
    use mknanMod
    use mkdomainMod

    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Authors: Gordon Bonan, Sam Levis and Mariana Vertenstein
! Revised: Nan Rosenbloom to add fmax processing.
! 3/18/08: David Lawrence added organic matter processing
! 1/22/09: Keith Oleson added urban parameter processing
!
!

! !LOCAL VARIABLES:
    integer  :: nsoicol                     ! number of model color classes
    integer  :: k,m,n                       ! indices
    integer  :: ni,nj,ns_o                  ! indices
    integer  :: ier                         ! error status
    integer,save  :: ndiag                       ! unit numbers
    integer  :: ncid                        ! netCDF id
    integer  :: omode                       ! netCDF output mode
    integer  :: varid                       ! netCDF variable id
    integer  :: ndims                       ! netCDF number of dimensions
    integer  :: beg(4),len(4),dimids(4)     ! netCDF dimension sizes 
    integer  :: ret                         ! netCDF return status
    integer,save  :: ntim                        ! time sample for dynamic land use
    logical  :: all_veg                     ! if gridcell will be 100% vegetated land-cover
    logical,save :: first_call = .true.         ! first call logical
    real(r8) :: suma                        ! sum for error check
    character(len=256) :: fgrddat           ! grid data file
    character(len=256) :: fdyndat           ! dynamic landuse data file name
    character(len=256) :: fdyndat1          ! dynamic landuse data file name
    character(len=256) :: fdyndat2          ! dynamic landuse data file name
    character(len=256) :: fname             ! generic filename
    character(len=256) :: loc_fn            ! local file name
    character(len=256) :: fsurdat           ! output surface data file name
    character(len=256) :: fsurdat1          ! output surface data file name
    character(len=256) :: fsurdat2          ! output surface data file name
    character(len=256) :: fsurlog           ! output surface log file name
    integer  :: t1                          ! timer
    real(r8),parameter :: p5  = 0.5_r8      ! constant
    real(r8),parameter :: p25 = 0.25_r8     ! constant
    integer  :: nunit                       ! unit number

    real(r8), allocatable  :: landfrac_pft(:)    ! PFT data: % land per gridcell
    real(r8), allocatable  :: pctlnd_pft(:)      ! PFT data: % of gridcell for PFTs
    real(r8), allocatable  :: pctlnd_pft_dyn(:)  ! PFT data: % of gridcell for dyn landuse PFTs
    integer , allocatable  :: pftdata_mask(:)    ! mask indicating real or fake land type
    real(r8), pointer      :: pctpft(:,:)        ! PFT data: land fraction per gridcell
    real(r8), pointer      :: harvest(:,:)       ! harvest data: normalized harvesting
    real(r8), pointer      :: pctpft_i(:,:)      ! PFT data: % fraction on input grid
    real(r8), allocatable  :: pctgla(:)          ! percent of grid cell that is glacier  
    real(r8), allocatable  :: pctglcmec(:,:)     ! glacier_mec pct coverage in each gridcell and class
    real(r8), allocatable  :: topoglcmec(:,:)    ! glacier_mec sfc elevation in each gridcell and class
    real(r8), allocatable  :: thckglcmec(:,:)    ! glacier_mec ice sheet thcknss in each gridcell and class
    real(r8), allocatable  :: elevclass(:)       ! glacier_mec elevation classes
    real(r8), allocatable  :: pctlak(:)          ! percent of grid cell that is lake     
    real(r8), allocatable  :: pctwet(:)          ! percent of grid cell that is wetland  
    real(r8), allocatable  :: pctirr(:)          ! percent of grid cell that is irrigated  
    real(r8), allocatable  :: pcturb(:)          ! percent of grid cell that is urbanized
    real(r8), allocatable  :: elev(:)            ! glc elevation (m)
    real(r8), allocatable  :: topo(:)            ! land elevation (m)
    real(r8), allocatable  :: fmax(:)            ! fractional saturated area
    integer , allocatable  :: soicol(:)          ! soil color                            
    real(r8), allocatable  :: pctsand(:,:)       ! soil texture: percent sand            
    real(r8), allocatable  :: pctclay(:,:)       ! soil texture: percent clay            
    real(r8), allocatable  :: ef1_btr(:)         ! Isoprene emission factor for broadleaf
    real(r8), allocatable  :: ef1_fet(:)         ! Isoprene emission factor for fine/everg
    real(r8), allocatable  :: ef1_fdt(:)         ! Isoprene emission factor for fine/dec
    real(r8), allocatable  :: ef1_shr(:)         ! Isoprene emission factor for shrubs
    real(r8), allocatable  :: ef1_grs(:)         ! Isoprene emission factor for grasses
    real(r8), allocatable  :: ef1_crp(:)         ! Isoprene emission factor for crops
    real(r8), allocatable  :: organic(:,:)       ! organic matter density (kg/m3)            

    type(domain1_type) :: ldomain

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

    subroutine mksurfdat_run(year,plodata)

!
! !ARGUMENTS:
    implicit none
    integer,intent(in)  :: year                        ! year for dynamic land use
    real(r8),pointer,intent(in) :: plodata(:,:)

    logical :: exists
    character(len=32) :: subname = 'mksrfdat_run'  ! name

    namelist /mksurfnml/              &
	 mksrf_fgrid,              &	
	 mksrf_gridnm,             &	
	 mksrf_gridtype,           &	
         mksrf_fvegtyp,            &
	 mksrf_fsoitex,            &
         mksrf_forganic,           &
         mksrf_fsoicol,            &
         mksrf_fvocef,             &
         mksrf_flanwat,            &
         mksrf_fglacier,           &
         mksrf_fglctopo,           &
         mksrf_flndtopo,           &
         mksrf_ffrac,              &
         mksrf_fmax,               &
         mksrf_furban,             &
         mksrf_flai,               &
         mksrf_firrig,             &
         mksrf_fdynuse,            &
         nglcec,                   &
         soil_color,               &
         soil_sand,                &
         soil_clay,                &
         pft_idx,                  &
         pft_frc,                  &
         all_urban,                &
         map_fpft,                 &
         map_flanwat,              &
         map_fglacier,             &
         map_fsoitex,              &
         map_fsoicol,              &
         map_furban,               &
         map_fglctopo,             &
         map_flndtopo,             &
         map_fmax,                 &
         map_forganic,             &
         map_fvocef,               &
         map_fglcmec_t2g,          &
         map_fglcmec_g2g,          &
         map_flai,                 &
         map_fharvest,             &
         map_firrig,               &
         outnc_large_files,        &
         outnc_double,             &
         outnc_dims,               &
         fsurdat,                  &
         fdyndat,                  &   
         fsurlog     

!-----------------------------------------------------------------------

    ! ======================================================================
    ! Read input namelist
    ! ======================================
    ! Must specify settings for the output grid:
    ! ======================================
    !    mksrf_fgrid -- Grid dataset
    ! ======================================
    ! Must specify settings for input high resolution datafiles
    ! ======================================
    !    mksrf_ffrac ---- land fraction and land mask dataset
    !    mksrf_fglacier - Glacier dataset
    !    mksrf_flai ----- Leaf Area Index dataset
    !    mksrf_flanwat -- Land water dataset
    !    mksrf_forganic - Organic soil carbon dataset
    !    mksrf_fmax ----- Max fractional saturated area dataset
    !    mksrf_fsoicol -- Soil color dataset
    !    mksrf_fsoitex -- Soil texture dataset
    !    mksrf_fglctopo-- Topography dataset (for glacier multiple elevation classes)
    !                     (and for limiting urban areas)
    !    mksrf_furban --- Urban dataset
    !    mksrf_fvegtyp -- PFT vegetation type dataset
    !    mksrf_fvocef  -- Volatile Organic Compund Emission Factor dataset
    ! ======================================
    ! Optionally specify setting for:
    ! ======================================
    !    mksrf_firrig ------ Irrigation dataset
    !    mksrf_fdynuse ----- input filename (deprecated, now use data from interface, tc)
    !    mksrf_gridtype ---- Type of grid (default is 'global')
    !    mksrf_gridnm ------ Name of output grid resolution
    !    outnc_double ------ If output should be in double precision
    !    outnc_large_files - If output should be in NetCDF large file format
    !    nglcec ------------ If you want to change the number of Glacier elevation classes
    ! ======================================
    ! Optional settings to change values for entire area
    ! ======================================
    !    all_urban --------- If entire area is urban
    !    soil_color -------- If you want to change the soil_color to this value everywhere
    !    soil_clay --------- If you want to change the soil_clay % to this value everywhere
    !    soil_sand --------- If you want to change the soil_sand % to this value everywhere
    !    pft_idx ----------- If you want to change to 100% veg covered with given PFT indices
    !    pft_frc ----------- Fractions that correspond to the pft_idx above
    ! ==================
    ! ======================================================================

    write(6,*) 'Attempting to initialize control settings .....'

    mksrf_gridtype    = 'global'
    outnc_large_files = .false.
    outnc_double      = .true.
    all_urban         = .false.

    nunit = shr_file_getUnit()
    open(nunit,file="iac_in",status="old",action="read")
    read(nunit, mksurfnml, iostat=ier)
    if (ier /= 0) then
       write(6,*)'error: namelist input resulted in error code ',ier
       call abort()
    endif
    close(nunit)
    call shr_file_freeUnit(nunit)

    write (6,*) 'Attempting to create surface boundary data .....'
    write (6,'(72a1)') ("-",n=1,60)

    ! ----------------------------------------------------------------------
    ! Error check namelist input
    ! ----------------------------------------------------------------------
    
    if (mksrf_fgrid /= ' ')then
       fgrddat = mksrf_fgrid
       write(6,*)'mksrf_fgrid = ',mksrf_fgrid
    else
       write (6,*)'must specify mksrf_fgrid'
       call abort()
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    else
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
       write (6,*)'illegal mksrf_gridtype, must be global or regional '
       call abort()
    endif
    if ( outnc_large_files )then
       write(6,*)'Output file in NetCDF 64-bit large_files format'
    end if
    if ( outnc_double )then
       write(6,*)'Output ALL data in file as 64-bit'
    end if
    if ( all_urban )then
       write(6,*) 'Output ALL data in file as 100% urban'
    end if

    allocate ( elevclass(nglcec+1) )

    call mksoiltexInit( )
    call mksoilcolInit( )
    call mkpftInit( all_urban, all_veg )
    call mkglcmecInit (elevclass)

    if ( all_veg )then
       write(6,*) 'Output ALL data in file as 100% vegetated'
    end if

    ! ----------------------------------------------------------------------
    ! Determine land model grid, fractional land and land mask
    ! ----------------------------------------------------------------------
    
    write(6,*)'calling domain_read'
    call domain1_read_map(ldomain, fgrddat)
    write(6,*)'finished domain_read'
    
    ! Invalidate mask and frac for ldomain 

    ldomain%mask = bigint
    ldomain%frac = nan

    ! Determine if will have 1d output

    if (ldomain%ni /= -9999 .and. ldomain%nj /= -9999) then
       write(6,*)'fsurdat is 2d lat/lon grid'
       write(6,*)'nlon= ',ldomain%ni,' nlat= ',ldomain%nj
       if (outnc_dims == 1) then
          write(6,*)' writing output file in 1d gridcell format'
       end if
    else
       write(6,*)'fsurdat is 1d gridcell grid'
       outnc_dims = 1
    end if

    ! ----------------------------------------------------------------------
    ! Allocate and initialize dynamic memory
    ! ----------------------------------------------------------------------

    ns_o = ldomain%ns
    allocate ( landfrac_pft(ns_o)      , &
               pctlnd_pft(ns_o)        , & 
               pftdata_mask(ns_o)      , & 
               pctpft(ns_o,0:numpft)   , & 
               pctgla(ns_o)            , & 
               pctlak(ns_o)            , & 
               pctwet(ns_o)            , & 
               pcturb(ns_o)            , & 
               pctsand(ns_o,nlevsoi)   , & 
               pctclay(ns_o,nlevsoi)   , & 
               soicol(ns_o)              )
    landfrac_pft(:) = spval 
    pctlnd_pft(:)   = spval
    pftdata_mask(:) = -999
    pctpft(:,:)     = spval
    pctgla(:)       = spval
    pctlak(:)       = spval
    pctwet(:)       = spval
    pcturb(:)       = spval
    pctsand(:,:)    = spval
    pctclay(:,:)    = spval
    soicol(:)       = -999

    ! ----------------------------------------------------------------------
    ! Open diagnostic output log file
    ! ----------------------------------------------------------------------
    
    if (fsurlog == ' ') then
       write(6,*)' must specify fsurlog in namelist'
       stop
    else
       if (first_call) then
          ndiag = getavu(); call opnfil (fsurlog, ndiag, 'f')
       endif
    end if
    
    if (mksrf_fgrid /= ' ')then
       write (ndiag,*)'using fractional land data from file= ', &
            trim(mksrf_fgrid),' to create the surface dataset'
    endif

    if (trim(mksrf_gridtype) == 'global' .or. &
        trim(mksrf_gridtype) == 'regional') then
       write(6,*)'mksrf_gridtype = ',trim(mksrf_gridtype)
    endif

    write(ndiag,*) 'PFTs from:                  ',trim(mksrf_fvegtyp)
    write(ndiag,*) 'fmax from:                  ',trim(mksrf_fmax)
    write(ndiag,*) 'glaciers from:              ',trim(mksrf_fglacier)
    write(ndiag,*) 'glc topography from:        ',trim(mksrf_fglctopo)
    write(ndiag,*) '           with:            ', nglcec, ' glacier elevation classes'
    write(ndiag,*) 'land topography from:       ',trim(mksrf_flndtopo)
    write(ndiag,*) 'urban from:                 ',trim(mksrf_furban)
    write(ndiag,*) 'inland water from:          ',trim(mksrf_flanwat)
    write(ndiag,*) 'soil texture from:          ',trim(mksrf_fsoitex)
    write(ndiag,*) 'soil organic from:          ',trim(mksrf_forganic)
    write(ndiag,*) 'soil color from:            ',trim(mksrf_fsoicol)
    write(ndiag,*) 'VOC emission factors from:  ',trim(mksrf_fvocef)
    if (mksrf_firrig /= ' ') then
       write(ndiag,*)&
                   'irrigated area from:        ',trim(mksrf_firrig)
    endif
    write(ndiag,*)' mapping for pft             ',trim(map_fpft)
    write(ndiag,*)' mapping for lanwat          ',trim(map_flanwat)
    write(ndiag,*)' mapping for glacier         ',trim(map_fglacier)
    write(ndiag,*)' mapping for soil texture    ',trim(map_fsoitex)
    write(ndiag,*)' mapping for soil color      ',trim(map_fsoicol)
    write(ndiag,*)' mapping for soil organic    ',trim(map_forganic)
    write(ndiag,*)' mapping for urban           ',trim(map_furban)
    write(ndiag,*)' mapping for fmax            ',trim(map_fmax)
    write(ndiag,*)' mapping for VOC pct emis    ',trim(map_fvocef)
    write(ndiag,*)' mapping for harvest         ',trim(map_fharvest)
    write(ndiag,*)' mapping for irrigation      ',trim(map_firrig)
    write(ndiag,*)' mapping for lai/sai         ',trim(map_flai)
    write(ndiag,*)' mapping for glc topography  ',trim(map_fglctopo)
    write(ndiag,*)' mapping for land topography ',trim(map_flndtopo)
    write(ndiag,*)' mapping for topo to raw glacier          ',trim(map_fglcmec_t2g)
    write(ndiag,*)' mapping for raw glacier to model glacier ',trim(map_fglcmec_g2g)

    ! ----------------------------------------------------------------------
    ! Make surface dataset fields
    ! ----------------------------------------------------------------------

    ! Make irrigated area fraction [pctirr] from [firrig] dataset if requested in namelist

    if (mksrf_firrig /= ' ') then
       allocate ( pctirr(ns_o) )
       pctirr(:) = spval
       call mkirrig(ldomain, mapfname=map_firrig, datfname=mksrf_firrig,&
            ndiag=ndiag, irrig_o=pctirr)
    endif

    ! Make PFTs [pctpft] from dataset [fvegtyp] (1/2 degree PFT data)

    nullify(pctpft_i)
    call mkpft(ldomain, mapfname=map_fpft, fpft=mksrf_fvegtyp, &
         firrig=mksrf_firrig, ndiag=ndiag, pctlnd_o=pctlnd_pft, pctirr_o=pctirr, &
         pctpft_o=pctpft, pct_pft_i=pctpft_i)

    ! Make inland water [pctlak, pctwet] from Cogley's one degree data [flanwat]

    call mklanwat (ldomain, mapfname=map_flanwat, datfname=mksrf_flanwat, &
         ndiag=ndiag, zero_out=all_urban.or.all_veg, lake_o=pctlak, swmp_o=pctwet)

    ! Make glacier fraction [pctgla] from [fglacier] dataset

    call mkglacier (ldomain, mapfname=map_fglacier, datfname=mksrf_fglacier, &
         ndiag=ndiag, zero_out=all_urban.or.all_veg, glac_o=pctgla)

    ! Make soil texture [pctsand, pctclay] from IGBP 5 minute data [fsoitex]

    call mksoiltex (ldomain, mapfname=map_fsoitex, datfname=mksrf_fsoitex, &
         ndiag=ndiag, pctglac_o=pctgla, sand_o=pctsand, clay_o=pctclay)

    ! Make soil color classes [soicol] from BATS T42 data [fsoicol]

    call mksoilcol (ldomain, mapfname=map_fsoicol, datfname=mksrf_fsoicol, &
         ndiag=ndiag, pctglac_o=pctgla, soil_color_o=soicol, nsoicol=nsoicol)

    ! Make urban fraction [pcturb] from [furban] dataset

    call mkurban (ldomain, mapfname=map_furban, datfname=mksrf_furban, &
         ndiag=ndiag, zero_out=all_veg, urbn_o=pcturb)


    ! Make elevation [elev] from [ftopo, ffrac] dataset
    ! Used only to screen pcturb  
    ! Screen pcturb by elevation threshold from elev dataset

    allocate(elev(ns_o))
    elev(:) = spval
    call mkelev (ldomain, mapfname=map_fglctopo, datfname=mksrf_fglctopo, &
         varname='TOPO_ICE', ndiag=ndiag, elev_o=elev)

    if ( .not. all_urban )then
       where (elev .gt. elev_thresh)
         pcturb = 0._r8
       end where
    end if
    deallocate(elev)
    
    ! Determine topography

    allocate(topo(ns_o))
    call mkelev (ldomain, mapfname=map_flndtopo, datfname=mksrf_flndtopo, &
         varname='TOPO', ndiag=ndiag, elev_o=topo)

    ! Make fmax [fmax] from [fmax] dataset

    allocate(fmax(ns_o))
    fmax(:) = spval
    call mkfmax (ldomain, mapfname=map_fmax, datfname=mksrf_fmax, &
         ndiag=ndiag, fmax_o=fmax)

    ! Make organic matter density [organic] from Global Soil Data Task [forganic]
    allocate (organic(ns_o,nlevsoi))
    organic(:,:) = spval
    call mkorganic (ldomain, mapfname=map_forganic, datfname=mksrf_forganic, &
         ndiag=ndiag, organic_o=organic)

    ! Make VOC emission factors for isoprene &
    ! [ef1_btr,ef1_fet,ef1_fdt,ef1_shr,ef1_grs,ef1_crp] from 
    ! MEGAN data [fvocef] (heald, 04/06)

    allocate ( ef1_btr(ns_o) , & 
               ef1_fet(ns_o) , & 
               ef1_fdt(ns_o) , & 
               ef1_shr(ns_o) , & 
               ef1_grs(ns_o) , & 
               ef1_crp(ns_o) )
    ef1_btr(:) = 0._r8
    ef1_fet(:) = 0._r8
    ef1_fdt(:) = 0._r8
    ef1_shr(:) = 0._r8
    ef1_grs(:) = 0._r8
    ef1_crp(:) = 0._r8
    
    call mkvocef (ldomain, mapfname=map_fvocef, datfname=mksrf_fvocef, ndiag=ndiag, &
         ef_btr_o=ef1_btr, ef_fet_o=ef1_fet, ef_fdt_o=ef1_fdt,  &
         ef_shr_o=ef1_shr, ef_grs_o=ef1_grs, ef_crp_o=ef1_crp)

    ! Do landuse changes such as for the poles, Ross ice-shelf and etc.
    call change_landuse( ldomain, dynpft=.false. )

    do n = 1,ns_o

       ! Assume wetland and/or lake when dataset landmask implies ocean 
       ! (assume medium soil color (15) and loamy texture).
       ! Also set pftdata_mask here
       
       if (pctlnd_pft(n) < 1.e-6_r8) then
          pftdata_mask(n)  = 0
          soicol(n)        = 15
          pctwet(n)        = 100._r8 - pctlak(n)
          pcturb(n)        = 0._r8
          if (mksrf_firrig /= ' ') pctirr(n) = 0._r8
          pctgla(n)        = 0._r8
          pctpft(n,:)      = 0._r8
          pctsand(n,:)     = 43._r8
          pctclay(n,:)     = 18._r8
          organic(n,:)   = 0._r8
       else
          pftdata_mask(n) = 1
       end if

       ! Truncate all percentage fields on output grid. This is needed to
       ! insure that wt is not nonzero (i.e. a very small number such as
       ! 1e-16) where it really should be zero
       
       do k = 1,nlevsoi
          pctsand(n,k) = float(nint(pctsand(n,k)))
          pctclay(n,k) = float(nint(pctclay(n,k)))
       end do
       pctlak(n) = float(nint(pctlak(n)))
       pctwet(n) = float(nint(pctwet(n)))
       pctgla(n) = float(nint(pctgla(n)))
       if (mksrf_firrig /= ' ') pctirr(n) = float(nint(pctirr(n)))
       
       ! Make sure sum of land cover types does not exceed 100. If it does,
       ! subtract excess from most dominant land cover.
       
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma > 130._r4) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb and pctgla is greater than 130%'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n)
          call abort()
       else if (suma > 100._r4) then
          pctlak(n) = pctlak(n) * 100._r8/suma
          pctwet(n) = pctwet(n) * 100._r8/suma
          pcturb(n) = pcturb(n) * 100._r8/suma
          pctgla(n) = pctgla(n) * 100._r8/suma
       end if
       
    end do

    call normalizencheck_landuse(ldomain)

    ! Write out sum of PFT's

    do k = 0,numpft
       suma = 0._r8
       do n = 1,ns_o
          suma = suma + pctpft(n,k)
       enddo
       write(6,*) 'sum over domain of pft ',k,suma
    enddo
    write(6,*)

    ! Make glacier multiple elevation classes [pctglcmec,topoglcmec] from [fglacier,ftopo] dataset
    ! This call needs to occur after pctgla has been adjusted for the final time

    allocate (pctglcmec(ns_o,nglcec),&
              topoglcmec(ns_o,nglcec),&
              thckglcmec(ns_o,nglcec))
    pctglcmec(:,:)  = spval
    topoglcmec(:,:) = spval
    thckglcmec(:,:) = spval
    call mkglcmec (ldomain, mapfname_t2g=map_fglcmec_t2g, mapfname_g2g=map_fglcmec_g2g, &
                   datfname_fglctopo=mksrf_fglctopo, datfname_fglacier=mksrf_fglacier, &
                   ndiag=ndiag, pctglac_o=pctgla, pctglcmec_o=pctglcmec, &
                   topoglcmec_o=topoglcmec, thckglcmec_o=thckglcmec )

    ! Determine fractional land from pft dataset

    do n = 1,ns_o
       landfrac_pft(n) = pctlnd_pft(n)/100._r8
    end do

    ! ----------------------------------------------------------------------
    ! Create surface dataset
    ! ----------------------------------------------------------------------

    ! Create netCDF surface dataset.  

    if (fsurdat == ' ') then
       write(6,*)' must specify fsurdat in namelist'
       stop
    end if

    fsurdat1 = fsurdat(1:len_trim(fsurdat)-2)
    fsurdat2 = fsurdat(len_trim(fsurdat)-1:len_trim(fsurdat))
    if (trim(fsurdat2) == 'nc') then
       write(fsurdat,'(a,i4.4,a)') trim(fsurdat1),year,'.nc'
    endif
    write(6,*) 'fsurdat output file = ',trim(fsurdat),':',trim(fsurdat1),':',trim(fsurdat2)

    call mkfile(ldomain, trim(fsurdat), dynlanduse = .false.)

    call domain1_write(ldomain, fsurdat)

    call check_ret(nf_open(trim(fsurdat), nf_write, ncid), subname)
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Write fields OTHER THAN lai, sai, heights, and urban parameters to netcdf surface dataset

    call check_ret(nf_inq_varid(ncid, 'PFTDATA_MASK', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, pftdata_mask), subname)

    call check_ret(nf_inq_varid(ncid, 'LANDFRAC_PFT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, landfrac_pft), subname)

    call check_ret(nf_inq_varid(ncid, 'mxsoil_color', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, nsoicol), subname)

    call check_ret(nf_inq_varid(ncid, 'SOIL_COLOR', varid), subname)
    call check_ret(nf_put_var_int(ncid, varid, soicol), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_SAND', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctsand), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_CLAY', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctclay), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctwet), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctlak), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctgla), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_GLC_MEC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctglcmec), subname)

    call check_ret(nf_inq_varid(ncid, 'GLC_MEC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, elevclass), subname)

    call check_ret(nf_inq_varid(ncid, 'TOPO_GLC_MEC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, topoglcmec), subname)

    call check_ret(nf_inq_varid(ncid, 'TOPO', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, topo), subname)

    call check_ret(nf_inq_varid(ncid, 'THCK_GLC_MEC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, thckglcmec), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pcturb), subname)

    call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctpft), subname)

    call check_ret(nf_inq_varid(ncid, 'FMAX', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, fmax), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_BTR', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_btr), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_FET', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_fet), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_FDT', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_fdt), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_SHR', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_shr), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_GRS', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_grs), subname)

    call check_ret(nf_inq_varid(ncid, 'EF1_CRP', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, ef1_crp), subname)

    call check_ret(nf_inq_varid(ncid, 'ORGANIC', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, organic), subname)

    if (mksrf_firrig /= ' ') then
    call check_ret(nf_inq_varid(ncid, 'PCT_IRRIG', varid), subname)
    call check_ret(nf_put_var_double(ncid, varid, pctirr), subname)
    endif

    ! Deallocate arrays NOT needed for dynamic-pft section of code

    deallocate ( organic )
    deallocate ( ef1_btr, ef1_fet, ef1_fdt, ef1_shr, ef1_grs, ef1_crp )
    deallocate ( pctglcmec, topoglcmec, thckglcmec, elevclass )
    deallocate ( fmax )
    deallocate ( pctsand, pctclay )
    deallocate ( soicol )
    deallocate ( topo )

    ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

    call check_ret(nf_sync(ncid), subname)

    ! ----------------------------------------------------------------------
    ! Make Urban Parameters from 1/2 degree data and write to surface dataset 
    ! Write to netcdf file is done inside mkurbanpar routine
    ! Only call this routine if pcturb is greater than zero somewhere.  Raw urban
    ! datasets will have no associated parameter fields if there is no urban 
    ! (e.g., mksrf_urban.060929.nc).
    ! ----------------------------------------------------------------------

    write(6,*)'calling mkurbanpar'
    if (any(pcturb > 0._r8)) then
       call mkurbanpar(ldomain, mapfname=map_furban, datfname=mksrf_furban, &
            ndiag=ndiag, ncido=ncid)
    else
       write(6,*) 'PCT_URBAN is zero everywhere, no urban parameter fields will be created'
    end if

    ! ----------------------------------------------------------------------
    ! Make LAI and SAI from 1/2 degree data and write to surface dataset 
    ! Write to netcdf file is done inside mklai routine
    ! ----------------------------------------------------------------------

    write(6,*)'calling mklai'
    if ( .not. associated(pctpft_i) )then
       write(6,*)'error: pctpft_i is not allocated at this point'
       call abort()
    end if
    call mklai(ldomain, mapfname=map_flai, datfname=mksrf_flai, &
         firrig=mksrf_firrig, ndiag=ndiag, ncido=ncid, pctpft_i=pctpft_i)
    deallocate( pctpft_i )

    ! Close surface dataset

    call check_ret(nf_close(ncid), subname)

    write (6,'(72a1)') ("-",n=1,60)
!    write (6,'(a,f5.1,a4,f5.1,a5)') 'land model surface data set successfully created for ', &
!         'grid of size ',ns_o
    write (6,*) 'land model surface data set successfully created for ', &
         'grid of size ',ns_o

    ! ----------------------------------------------------------------------
    ! Create dynamic land use dataset if appropriate
    ! ----------------------------------------------------------------------

    if (mksrf_fdynuse /= ' ') then

       write(6,*)'creating dynamic land use dataset'

       allocate(pctlnd_pft_dyn(ns_o))
       call mkharvest_init( ns_o, spval, harvest, mksrf_fvegtyp )

       if (fdyndat == ' ') then
          write(6,*)' must specify fdyndat in namelist if mksrf_fdynuse is not blank'
          stop
       end if

       fdyndat1 = fdyndat(1:len_trim(fdyndat)-2)
       fdyndat2 = fdyndat(len_trim(fdyndat)-1:len_trim(fdyndat))
!tcxof       if (trim(fdyndat2) == 'nc') then
!          write(fdyndat,'(a,i4.4,a)') trim(fdyndat1),year,'.nc'
!       endif
       write(6,*) 'fdyndat output file = ',trim(fdyndat),':',trim(fdyndat1),':',trim(fdyndat2)
       inquire(file=trim(fdyndat),exist=exists)

       write(6,*) subname,'fdyndat exists value ',exists

       if (.not.exists) then

          ! Define dimensions and global attributes

          call mkfile(ldomain, trim(fdyndat), dynlanduse=.true.)

          ! Write fields other pft to dynamic land use dataset

          call domain1_write(ldomain, fdyndat)

          call check_ret(nf_open(trim(fdyndat), nf_write, ncid), subname)
          call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

          call check_ret(nf_inq_varid(ncid, 'PFTDATA_MASK', varid), subname)
          call check_ret(nf_put_var_int(ncid, varid, pftdata_mask), subname)

          call check_ret(nf_inq_varid(ncid, 'LANDFRAC_PFT', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, landfrac_pft), subname)

          call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, pctwet), subname)

          call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, pctlak), subname)

          call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, pctgla), subname)

          call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
          call check_ret(nf_put_var_double(ncid, varid, pcturb), subname)

          ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

          call check_ret(nf_sync(ncid), subname)

          ! Read in each dynamic pft landuse dataset

       else

          call check_ret(nf_open(trim(fdyndat), nf_write, ncid), subname)
          call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)
          call check_ret(nf_sync(ncid), subname)

       endif
!tcx no loop       do 
          ! Read input pft data

          !
          ! If pft fraction override is set, than intrepret mksrf_fdynuse as PFT and harvesting override values
          !
          if ( any(pft_frc > 0.0_r8 ) )then
             fname = ' '
             call mkpft_parse_oride(mksrf_fdynuse)
             call mkharvest_parse_oride(mksrf_fdynuse)
	     write(6,*)'PFT and harvesting values are ',trim(mksrf_fdynuse),' year is ',year
          !
          ! Otherwise intrepret mksrf_fdynuse as a filename with PFT and harvesting values in it
          !
          else
             fname = mksrf_fdynuse
	     write(6,*)'input pft dynamic dataset is  ',trim(fname),' year is ',year
          end if

          ! Create pctpft data at model resolution

          call mkpft(ldomain, mapfname=map_fpft, fpft=fname, firrig=mksrf_firrig, &
               ndiag=ndiag, pctlnd_o=pctlnd_pft_dyn, pctirr_o=pctirr, &
               pctpft_o=pctpft, pct_pft_i=pctpft_i, plodata=plodata)

          ! Create harvesting data at model resolution

          call mkharvest( ldomain, mapfname=map_fharvest, datfname=fname, &
               ndiag=ndiag, harv_o=harvest , plodata=plodata)

          ! Consistency check on input land fraction

          do n = 1,ns_o
             if (pctlnd_pft_dyn(n) /= pctlnd_pft(n)) then
                write(6,*) subname,' error: pctlnd_pft for dynamics data = ',&
                     pctlnd_pft_dyn(n), ' not equal to pctlnd_pft for surface data = ',&
                     pctlnd_pft(n),' at n= ',n
                if ( trim(fname) == ' ' )then
                   write(6,*) ' PFT mksrf_fdynuse = ', mksrf_fdynuse
                else
                   write(6,*) ' PFT file = ', fname
                end if
                call abort()
             end if
          end do
          call change_landuse(ldomain, dynpft=.true.)

          call normalizencheck_landuse(ldomain)

          ! Output pctpft data for current year

          call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
          call check_ret(nf_inq_varndims(ncid, varid, ndims), subname)
          call check_ret(nf_inq_vardimid(ncid, varid, dimids), subname)
          beg(1:ndims-1) = 1
          do n = 1,ndims
             call check_ret(nf_inq_dimlen(ncid, dimids(n), len(n)), subname)
          end do
          ntim = len(ndims) + 1
          len(ndims) = 1
          beg(ndims) = ntim

          call check_ret(nf_put_vara_double(ncid, varid, beg, len, pctpft), subname)

!  The YEAR and time represent the PCTPFT YEAR and time.
          call check_ret(nf_inq_varid(ncid, 'YEAR', varid), subname)
          call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

          call check_ret(nf_inq_varid(ncid, 'time', varid), subname)
          call check_ret(nf_put_vara_int(ncid, varid, ntim, 1, year), subname)

          call check_ret(nf_inq_varid(ncid, 'input_pftdata_filename', varid), subname)
          call check_ret(nf_put_vara_text(ncid, varid, (/ 1, ntim /), (/ len_trim(mksrf_fdynuse), 1 /), trim(mksrf_fdynuse) ), subname)

          do k = 1, mkharvest_numtypes()
             call check_ret(nf_inq_varid(ncid, trim(mkharvest_fieldname(k)), varid), subname)
             call check_ret(nf_inq_varndims(ncid, varid, ndims), subname)
             call check_ret(nf_inq_vardimid(ncid, varid, dimids), subname)
             beg(1:ndims-1) = 1
             do n = 1,ndims-1
                call check_ret(nf_inq_dimlen(ncid, dimids(n), len(n)), subname)
             end do
             len(ndims) = 1
!The harvest rate represents the rate from curr_year - 1 to curr_year.  This is the GLM convention. 
!jt             beg(ndims) = ntim
             beg(ndims) = ntim - 1
             call check_ret(nf_put_vara_double(ncid, varid, beg, len, harvest(:,k)), subname)
          end do

	  ! Synchronize the disk copy of a netCDF dataset with in-memory buffers

	  call check_ret(nf_sync(ncid), subname)

!tcx no loop       end do   ! end of read loop

       deallocate(pctlnd_pft_dyn)
       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-create dynamic landust dataset   

    deallocate (landfrac_pft     , &
               pctlnd_pft        , & 
               pftdata_mask      , & 
               pctpft            , & 
               pctgla            , & 
               pctlak            , & 
               pctwet            , & 
               pcturb            )
    if (mksrf_firrig /= ' ') deallocate(pctirr)

    ! ----------------------------------------------------------------------
    ! Close diagnostic dataset
    ! ----------------------------------------------------------------------

!tcx    close (ndiag)
    write (6,*)
    write (6,*) 'Surface data output file = ',trim(fsurdat)
    write (6,*) '   This file contains the land model surface data'
    write (6,*) 'Diagnostic log file      = ',trim(fsurlog)
    write (6,*) '   See this file for a summary of the dataset'
    write (6,*)

    first_call = .false.

    end subroutine mksurfdat_run

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: change_landuse
!
! !INTERFACE:
subroutine change_landuse( ldomain, dynpft )
!
! !DESCRIPTION:
!
! Do landuse changes such as for the poles, Ross ice-shelf and etc.
!
! !USES:
    use mkdomainMod
    implicit none
!
! !ARGUMENTS:
    type(domain1_type)   :: ldomain
    logical, intent(in)  :: dynpft   ! if part of the dynpft section of code

!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    logical  :: first_time = .true.         ! flag if this is the first pass through or not
    integer ,parameter :: bdtroptree = 6    ! Index for broadleaf decidious tropical tree
    integer ,parameter :: bdtemptree = 7    ! Index for broadleaf decidious temperate tree
    integer ,parameter :: bdtempshrub = 10  ! Index for broadleaf decidious temperate shrub
    real(r8),parameter :: troplat = 23.5_r8 ! Latitude to define as tropical
    real(r8),parameter :: rosslat = -79._r8 ! Latitude to define as Ross ice-shelf
    integer  :: n,ns_o                       ! indices
    character(len=32) :: subname = 'change_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    ns_o = ldomain%ns
    do n = 1,ns_o

       ! Set pfts 7 and 10 to 6 in the tropics to avoid lais > 1000
       ! Using P. Thornton's method found in surfrdMod.F90 in clm3.5
       
       if (abs(ldomain%latc(n))<troplat .and. pctpft(n,bdtemptree)>0._r8) then
          pctpft(n,bdtroptree) = pctpft(n,bdtroptree) + pctpft(n,bdtemptree)
          pctpft(n,bdtemptree) = 0._r8
          if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', &
               bdtemptree, ' now added to pft ', bdtroptree
       end if
       if (abs(ldomain%latc(n))<troplat .and. pctpft(n,bdtempshrub)>0._r8) then
          pctpft(n,bdtroptree) = pctpft(n,bdtroptree) + pctpft(n,bdtempshrub)
          pctpft(n,bdtempshrub) = 0._r8
          if ( first_time ) write (6,*) subname, ' Warning: all wgt of pft ', &
               bdtempshrub, ' now added to pft ', bdtroptree
       end if
       first_time = .false.
       
       ! Set land values on Ross ice shelf to glacier
       ! non-dynamic-PFT part of the code
       
       if (ldomain%latc(n) < rosslat) then
          pctpft(n,:) = 0._r8
          pctlak(n)   = 0._r8
          pctwet(n)   = 0._r8
          pcturb(n)   = 0._r8
          pctgla(n)   = 100._r8
          if ( .not. dynpft )then
             ef1_btr(n)     = 0._r8
             ef1_fet(n)     = 0._r8
             ef1_fdt(n)     = 0._r8
             ef1_shr(n)     = 0._r8
             ef1_grs(n)     = 0._r8
             ef1_crp(n)     = 0._r8
             organic(n,:)   = 0._r8
             soicol(n)      = 0
             if (mksrf_firrig /= ' ') pctirr(n) = 0._r8
             pctsand(n,:)   = 0._r8
             pctclay(n,:)   = 0._r8
          end if
       end if

       ! If have pole points on grid - set south pole to glacier
       ! north pole is assumed as non-land
       
       if (abs((ldomain%latc(n) - 90._r8)) < 1.e-6_r8) then
          pctpft(n,:) = 0._r8
          pctlak(n)   = 0._r8
          pctwet(n)   = 0._r8
          pcturb(n)   = 0._r8
          pctgla(n)   = 100._r8
          if ( .not. dynpft )then
             organic(n,:)   = 0._r8
             ef1_btr(n)     = 0._r8
             ef1_fet(n)     = 0._r8
             ef1_fdt(n)     = 0._r8
             ef1_shr(n)     = 0._r8
             ef1_grs(n)     = 0._r8
             ef1_crp(n)     = 0._r8
             soicol(n)      = 0
             if (mksrf_firrig /= ' ') pctirr(n) = 0._r8
             pctsand(n,:)   = 0._r8
             pctclay(n,:)   = 0._r8
          end if
       end if

    end do

end subroutine change_landuse

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: normalizencheck_landuse
!
! !INTERFACE:
subroutine normalizencheck_landuse(ldomain)
!
! !DESCRIPTION:
!
! Normalize land use and make sure things add up to 100% as well as
! checking that things are as they should be.
!
! !USES:
    use mkdomainMod
    implicit none
! !ARGUMENTS:
    type(domain1_type)   :: ldomain
!
! !REVISION HISTORY:
! 9/10/09: Erik Kluzek spin off subroutine from original embedded code
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m,k,n,ns_o                  ! indices
    real(r8) :: suma                        ! sum for error check
    real(r8) :: bare_urb_diff               ! difference between bare soil and urban %
    real(r8) :: pcturb_excess               ! excess urban % not accounted for by bare soil
    real(r8) :: sumpft                      ! sum of non-baresoil pfts
    real(r8) :: sum8, sum8a                 ! sum for error check
    real(r4) :: sum4a                       ! sum for error check
    character(len=32) :: subname = 'normalizencheck_landuse'  ! subroutine name
!-----------------------------------------------------------------------

    ns_o = ldomain%ns
    do n = 1,ns_o
       if (pcturb(n) .gt. 0._r8) then
          
          ! Replace bare soil preferentially with urban
          suma = pctlak(n)+pctwet(n)+pctgla(n)
          bare_urb_diff = 0.01_r8 * pctpft(n,0) * (100._r8 - suma) - pcturb(n)
          pctpft(n,0) = max(0._r8,bare_urb_diff)
          pcturb_excess = abs(min(0._r8,bare_urb_diff))
          
          ! Normalize pctpft to be the remainder of [100 - (special landunits)]
          ! including any urban not accounted for by bare soil above
          sumpft = sum(pctpft(n,1:numpft))
          if (sumpft > 0._r8) then
             suma = pctlak(n)+pctwet(n)+pctgla(n)
             do m = 1, numpft
                pctpft(n,m) = 0.01_r8 * pctpft(n,m) * (100._r8 - suma) - &
                     pcturb_excess*pctpft(n,m)/sumpft
             end do
          end if

       else
          
          ! Normalize pctpft to be the remainder of [100 - (special landunits)]
          suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
          do m = 0, numpft
             pctpft(n,m) = 0.01_r8 * pctpft(n,m) * (100._r8 - suma)
          end do

       end if
       
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       do m = 0,numpft
          suma = suma + pctpft(n,m)
       end do
       
       if (suma < 90._r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla and pctpft is less than 90'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctpft(n,:)
          call abort()
       else if (suma > 100._r8 + 1.e-4_r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla and pctpft is greater than 100'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctpft,sum= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctpft(n,:),suma
          call abort()
       else
          pctlak(n)   = pctlak(n)   * 100._r8/suma
          pctwet(n)   = pctwet(n)   * 100._r8/suma
          pcturb(n)   = pcturb(n)   * 100._r8/suma
          pctgla(n)   = pctgla(n)   * 100._r8/suma
          pctpft(n,:) = pctpft(n,:) * 100._r8/suma
       end if
       
       ! Roundoff error fix
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma < 100._r8 .and. suma > (100._r8 - 100._r8*epsilon(suma))) then
          write (6,*) 'Special land units near 100%, but not quite for n,suma =',n,suma
          write (6,*) 'Adjusting special land units to 100%'
          if (pctlak(n) >= 25._r8) then
             pctlak(n) = 100._r8 - (pctwet(n) + pcturb(n) + pctgla(n))
          else if (pctwet(n) >= 25._r8) then
             pctwet(n) = 100._r8 - (pctlak(n) + pcturb(n) + pctgla(n))
          else if (pcturb(n) >= 25._r8) then
             pcturb(n) = 100._r8 - (pctlak(n) + pctwet(n) + pctgla(n))
          else if (pctgla(n) >= 25._r8) then
             pctgla(n) = 100._r8 - (pctlak(n) + pctwet(n) + pcturb(n))
          else
             write (6,*) subname, 'Error: sum of special land units nearly 100% but none is >= 25% at ', &
                  'n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctpft(n,:),suma = ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),pctpft(n,:),suma
             call abort()
          end if
          pctpft(n,:) = 0._r8
       end if
       
       suma = pctlak(n) + pctwet(n) + pcturb(n) + pctgla(n)
       if (suma < 100._r8-epsilon(suma) .and. suma > (100._r8 - 4._r8*epsilon(suma))) then
          write (6,*) subname, 'n,pctlak,pctwet,pcturb,pctgla,pctpft= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
               pctpft(n,:)
          call abort()
       end if
       do m = 0,numpft
          suma = suma + pctpft(n,m)
       end do
       if ( abs(suma-100._r8) > 1.e-10_r8) then
          write (6,*) subname, ' error: sum of pctlak, pctwet,', &
               'pcturb, pctgla and pctpft is NOT equal to 100'
          write (6,*)'n,pctlak,pctwet,pcturb,pctgla,pctpft,sum= ', &
               n,pctlak(n),pctwet(n),pcturb(n),pctgla(n),&
               pctpft(n,:), sum8
          call abort()
       end if
       
    end do

    ! Check that when pctpft identically zero, sum of special landunits is identically 100%

    if ( .not. outnc_double )then
       do n = 1,ns_o
          sum8  =        real(pctlak(n),r4)
          sum8  = sum8 + real(pctwet(n),r4)
          sum8  = sum8 + real(pcturb(n),r4)
          sum8  = sum8 + real(pctgla(n),r4)
          sum4a = 0.0_r4
          do k = 0,numpft
             sum4a = sum4a + real(pctpft(n,k),r4)
          end do
          if ( sum4a==0.0_r4 .and. sum8 < 100._r4-2._r4*epsilon(sum4a) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctpft==0 sum = ', sum8
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n), pctpft(n,:)
             call abort()
          end if
       end do
    else
       do n = 1,ns_o
          sum8  =        pctlak(n)
          sum8  = sum8 + pctwet(n)
          sum8  = sum8 + pcturb(n)
          sum8  = sum8 + pctgla(n)
          sum8a = 0._r8
          do k = 0,numpft
             sum8a = sum8a + pctpft(n,k)
          end do
          if ( sum8a==0._r8 .and. sum8 < (100._r8-4._r8*epsilon(sum8)) )then
             write (6,*) subname, ' error: sum of pctlak, pctwet,', &
                  'pcturb, pctgla is < 100% when pctpft==0 sum = ', sum8
             write (6,*) 'Total error, error/epsilon = ',100._r8-sum8, ((100._r8-sum8)/epsilon(sum8))
             write (6,*)'n,pctlak,pctwet,pcturb,pctgla,epsilon= ', &
                  n,pctlak(n),pctwet(n),pcturb(n),pctgla(n), pctpft(n,:), epsilon(sum8)
             call abort()
          end if
       end do
    end if

end subroutine normalizencheck_landuse

end module mksurfdat
