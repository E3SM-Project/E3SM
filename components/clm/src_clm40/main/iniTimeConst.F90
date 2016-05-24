!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iniTimeConst
!
! !INTERFACE:
subroutine iniTimeConst
!
! !DESCRIPTION:
! Initialize time invariant clm variables
! 1) removed references to shallow lake - since it is not used
! 2) ***Make z, zi and dz allocatable depending on if you
!    have lake or soil
! 3) rootfr only initialized for soil points
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
  use decompMod   , only : get_proc_bounds, get_proc_global
  use decompMod   , only : gsMap_lnd_gdc2glo
  use clm_atmlnd  , only : clm_a2l
  use clm_varpar  , only : nlevsoi, nlevgrnd, nlevlak, numpft, numrad, nlevurb
  use clm_varcon  , only : istice, istdlak, istwet, isturb, istice_mec,  &
                           icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv, &
                           zlak, dzlak, zsoi, dzsoi, zisoi, spval, &
                           albsat, albdry, secspday
  use clm_varctl  , only : fsurdat,scmlon,scmlat,single_column, iulog, use_cn, use_cndv
  use pftvarcon   , only : noveg, ntree, roota_par, rootb_par,  &
                           smpso, smpsc, fnitr, nbrdlf_dcd_brl_shrub, &
                           z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                           qe25, mp, c3psn, slatop, dsladlai, leafcn, flnr, woody, &
                           lflitcn, frootcn, livewdcn, deadwdcn, froot_leaf, stem_leaf, croot_stem, &
                           flivewd, fcur, lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig, &
                           leaf_long, evergreen, stress_decid, season_decid, &
                           resist, pftpar20, pftpar28, pftpar29, pftpar30, pftpar31, &
                           allom1s, allom2s, &
                           allom1 , allom2 , allom3  , reinickerp, dwood
  use pftvarcon       , only : graincn
  use clm_time_manager, only : get_step_size
  use abortutils      , only : endrun
  use fileutils       , only : getfil
  use organicFileMod  , only : organicrd 
  use spmdMod         , only : mpicom, MPI_INTEGER, masterproc
  use clm_varctl      , only : fsnowoptics, fsnowaging
  use SNICARMod       , only : SnowAge_init, SnowOptics_init
  use shr_scam_mod    , only : shr_scam_getCloseLatLon
  use ncdio_pio       , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile

!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod.
!
! !REVISION HISTORY:
! Created by Gordon Bonan.
! Updated to clm2.1 data structrues by Mariana Vertenstein
! 4/26/05, Peter Thornton: Eliminated exponential decrease in saturated hydraulic
!   conductivity (hksat) with depth. 
! Updated: Colette L. Heald (05/06) reading in VOC emission factors
! 27 February 2008: Keith Oleson; Qing Liu (2004) saturated hydraulic conductivity 
! and matric potential
! 29 February 2008: David Lawrence; modified soil thermal and hydraulic properties to
! account for organic matter
! 18 March 2008: David Lawrence; nlevgrnd changes
! 03/28/08 Mark Flanner, read in netcdf files for SNICAR parameters
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: ivt(:)             !  vegetation type index
  integer , pointer :: pcolumn(:)         ! column index of corresponding pft
  integer , pointer :: pgridcell(:)       ! gridcell index of corresponding pft
  integer , pointer :: clandunit(:)       ! landunit index of column
  integer , pointer :: cgridcell(:)       ! gridcell index of column
  integer , pointer :: ctype(:)           ! column type index
  integer , pointer :: ltype(:)           ! landunit type index
  real(r8), pointer :: thick_wall(:)      ! total thickness of urban wall
  real(r8), pointer :: thick_roof(:)      ! total thickness of urban roof
  real(r8), pointer :: lat(:)             ! gridcell latitude (radians)
!
! local pointers to implicit out arguments
!
  real(r8), pointer :: z(:,:)             ! layer depth (m)
  real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
  real(r8), pointer :: dz(:,:)            ! layer thickness depth (m)
  real(r8), pointer :: rootfr(:,:)        ! fraction of roots in each soil layer
  real(r8), pointer :: rootfr_road_perv(:,:) ! fraction of roots in each soil layer for urban pervious road
  real(r8), pointer :: rresis(:,:)        !root resistance by layer (0-1)  (nlevgrnd)	
  real(r8), pointer :: dewmx(:)           ! maximum allowed dew [mm]
  real(r8), pointer :: bsw(:,:)           ! Clapp and Hornberger "b" (nlevgrnd)  
  real(r8), pointer :: bsw2(:,:)          ! Clapp and Hornberger "b" for CN code
  real(r8), pointer :: psisat(:,:)        ! soil water potential at saturation for CN code (MPa)
  real(r8), pointer :: vwcsat(:,:)        ! volumetric water content at saturation for CN code (m3/m3)
  real(r8), pointer :: watsat(:,:)        ! volumetric soil water at saturation (porosity) (nlevgrnd) 
  real(r8), pointer :: watfc(:,:)         ! volumetric soil water at field capacity (nlevsoi)
  real(r8), pointer :: watdry(:,:)        ! btran parameter for btran=0
  real(r8), pointer :: watopt(:,:)        ! btran parameter for btran = 1
  real(r8), pointer :: hksat(:,:)         ! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd) 
  real(r8), pointer :: sucsat(:,:)        ! minimum soil suction (mm) (nlevgrnd) 
  real(r8), pointer :: csol(:,:)          ! heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 
  real(r8), pointer :: tkmg(:,:)          ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
  real(r8), pointer :: tkdry(:,:)         ! thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
  real(r8), pointer :: tksatu(:,:)        ! thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
  real(r8), pointer :: wtfact(:)          ! maximum saturated fraction for a gridcell
  real(r8), pointer :: smpmin(:)          ! restriction for min of soil potential (mm) (new)
  real(r8), pointer :: hkdepth(:)         ! decay factor (m)
  integer , pointer :: isoicol(:)         ! soil color class
  real(r8), pointer :: gwc_thr(:)         ! threshold soil moisture based on clay content
  real(r8), pointer :: mss_frc_cly_vld(:) ! [frc] Mass fraction clay limited to 0.20
  real(r8), pointer :: efisop(:,:)        ! emission factors for isoprene (ug isoprene m-2 h-1)
  real(r8), pointer :: max_dayl(:)        ! maximum daylength (s)
  real(r8), pointer :: sandfrac(:)
  real(r8), pointer :: clayfrac(:)
!
!
! !OTHER LOCAL VARIABLES:
!EOP
  type(file_desc_t)  :: ncid   ! netcdf id
  integer  :: n,j,ib,lev,bottom! indices
  integer  :: g,l,c,p          ! indices
  integer  :: m                ! vegetation type index
  real(r8) :: bd               ! bulk density of dry soil material [kg/m^3]
  real(r8) :: tkm              ! mineral conductivity
  real(r8) :: xksat            ! maximum hydraulic conductivity of soil [mm/s]
  real(r8) :: scalez = 0.025_r8   ! Soil layer thickness discretization (m)
  real(r8) :: clay,sand        ! temporaries
  real(r8) :: slope,intercept        ! temporary, for rooting distribution
  real(r8) :: temp, max_decl   ! temporary, for calculation of max_dayl
  integer  :: begp, endp       ! per-proc beginning and ending pft indices
  integer  :: begc, endc       ! per-proc beginning and ending column indices
  integer  :: begl, endl       ! per-proc beginning and ending landunit indices
  integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
  integer  :: numg             ! total number of gridcells across all processors
  integer  :: numl             ! total number of landunits across all processors
  integer  :: numc             ! total number of columns across all processors
  integer  :: nump             ! total number of pfts across all processors

  real(r8),pointer :: temp_ef(:)        ! read in - temporary EFs
  real(r8),pointer :: efisop2d(:,:)     ! read in - isoprene emission factors

  real(r8),pointer :: arrayl(:)      ! generic global array
  integer ,pointer :: irrayg(:)      ! generic global array
  integer ,pointer :: soic2d(:)      ! read in - soil color
  real(r8),pointer :: sand3d(:,:)    ! read in - soil texture: percent sand
  real(r8),pointer :: clay3d(:,:)    ! read in - soil texture: percent clay
  real(r8),pointer :: organic3d(:,:) ! read in - organic matter: kg/m3
  real(r8),pointer :: gti(:)         ! read in - fmax
  real(r8) :: om_frac                ! organic matter fraction
  real(r8) :: om_watsat    = 0.9_r8  ! porosity of organic soil
  real(r8) :: om_hksat     = 0.1_r8  ! saturated hydraulic conductivity of organic soil [mm/s]
  real(r8) :: om_tkm       = 0.25_r8 ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
  real(r8) :: om_sucsat    = 10.3_r8 ! saturated suction for organic matter (Letts, 2000)
  real(r8) :: om_csol      = 2.5_r8  ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
  real(r8) :: om_tkd       = 0.05_r8 ! thermal conductivity of dry organic soil (Farouki, 1981)
  real(r8) :: om_b         = 2.7_r8  ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
  real(r8) :: organic_max  = 130._r8 ! organic matter (kg/m3) where soil is assumed to act like peat 
  real(r8) :: csol_bedrock = 2.0e6_r8 ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
  real(r8) :: pc           = 0.5_r8   ! percolation threshold
  real(r8) :: pcbeta       = 0.139_r8 ! percolation exponent
  real(r8) :: perc_frac               ! "percolating" fraction of organic soil
  real(r8) :: perc_norm               ! normalize to 1 when 100% organic soil
  real(r8) :: uncon_hksat             ! series conductivity of mineral/organic soil
  real(r8) :: uncon_frac              ! fraction of "unconnected" soil
  integer  :: varid                  ! netCDF id's
  integer  :: ret

  integer  :: ier                                ! error status
  character(len=256) :: locfn                    ! local filename
  character(len= 32) :: subname = 'iniTimeConst' ! subroutine name
  integer :: mxsoil_color                        ! maximum number of soil color classes
  real(r8), allocatable :: zurb_wall(:,:)        ! wall (layer node depth)
  real(r8), allocatable :: zurb_roof(:,:)        ! roof (layer node depth)
  real(r8), allocatable :: dzurb_wall(:,:)       ! wall (layer thickness)
  real(r8), allocatable :: dzurb_roof(:,:)       ! roof (layer thickness)
  real(r8), allocatable :: ziurb_wall(:,:)       ! wall (layer interface)
  real(r8), allocatable :: ziurb_roof(:,:)       ! roof (layer interface)
  logical :: readvar 
!------------------------------------------------------------------------

  integer :: closelatidx,closelonidx
  real(r8):: closelat,closelon
  integer :: iostat

!------------------------------------------------------------------------

  if (masterproc) write(iulog,*) 'Attempting to initialize time invariant variables'

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  allocate(soic2d(begg:endg), gti(begg:endg))
  allocate(sand3d(begg:endg,nlevsoi), clay3d(begg:endg,nlevsoi))
  allocate(organic3d(begg:endg,nlevsoi))

  allocate(temp_ef(begg:endg),efisop2d(6,begg:endg))

  efisop          => gve%efisop

  ! Assign local pointers to derived subtypes components (gridcell-level)
  lat             => grc%lat
     
  ! Assign local pointers to derived subtypes components (landunit-level)

  ltype               => lun%itype
  thick_wall          => lps%thick_wall
  thick_roof          => lps%thick_roof

  ! Assign local pointers to derived subtypes components (column-level)

  ctype           => col%itype
  clandunit       => col%landunit
  cgridcell       => col%gridcell
  z               => cps%z
  dz              => cps%dz
  zi              => cps%zi
  bsw             => cps%bsw
  bsw2            => cps%bsw2
  psisat          => cps%psisat
  vwcsat          => cps%vwcsat
  watsat          => cps%watsat
  watfc           => cps%watfc
  watdry          => cps%watdry  
  watopt          => cps%watopt  
  rootfr_road_perv => cps%rootfr_road_perv
  hksat           => cps%hksat
  sucsat          => cps%sucsat
  tkmg            => cps%tkmg
  tksatu          => cps%tksatu
  tkdry           => cps%tkdry
  csol            => cps%csol
  smpmin          => cps%smpmin
  hkdepth         => cps%hkdepth
  wtfact          => cps%wtfact
  isoicol         => cps%isoicol
  gwc_thr         => cps%gwc_thr
  mss_frc_cly_vld => cps%mss_frc_cly_vld
  max_dayl        => cps%max_dayl

  ! Assign local pointers to derived subtypes components (pft-level)

  ivt             => pft%itype
  pgridcell       => pft%gridcell
  pcolumn         => pft%column
  dewmx           => pps%dewmx
  rootfr          => pps%rootfr
  rresis          => pps%rresis
  sandfrac        => pps%sandfrac
  clayfrac        => pps%clayfrac

  allocate(zurb_wall(begl:endl,nlevurb),    &
           zurb_roof(begl:endl,nlevurb),    &
           dzurb_wall(begl:endl,nlevurb),   &
           dzurb_roof(begl:endl,nlevurb),   &
           ziurb_wall(begl:endl,0:nlevurb), &
           ziurb_roof(begl:endl,0:nlevurb), &
           stat=ier)
  if (ier /= 0) then
     call endrun( 'iniTimeConst: allocation error for zurb_wall,zurb_roof,dzurb_wall,dzurb_roof,ziurb_wall,ziurb_roof' )
  end if

  ! --------------------------------------------------------------------
  ! Read soil color, sand and clay from surface dataset 
  ! --------------------------------------------------------------------

  if (masterproc) then
     write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
  end if

  call getfil (fsurdat, locfn, 0)
  call ncd_pio_openfile (ncid, locfn, 0)

  ! Determine number of soil color classes - if number of soil color classes is not
  ! on input dataset set it to 8

  if (single_column) then
     call shr_scam_getCloseLatLon(locfn,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
  end if
  call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, &
              readvar=readvar)
  if ( .not. readvar ) mxsoil_color = 8  

  ! Read fmax

  call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun( trim(subname)//' ERROR: FMAX NOT on surfdata file') 

  ! Read in soil color, sand and clay fraction

  call ncd_io(ncid=ncid, varname='SOIL_COLOR', flag='read', data=soic2d, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun( trim(subname)//' ERROR: SOIL_COLOR NOT on surfdata file' ) 

  ! Read in emission factors

  call ncd_io(ncid=ncid, varname='EF1_BTR', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_BTR')
  efisop2d(1,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='EF1_FET', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_FET')
  efisop2d(2,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='EF1_FDT', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_FDT')
  efisop2d(3,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='EF1_SHR', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_SHR')
  efisop2d(4,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='EF1_GRS', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_GRS')
  efisop2d(5,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='EF1_CRP', flag='read', data=temp_ef, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun('iniTimeConst: errror reading EF1_CRP')
  efisop2d(6,:)=temp_ef(:)

  call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_SAND NOT on surfdata file' ) 

  call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd, readvar=readvar)
  if (.not. readvar) call endrun( trim(subname)//' ERROR: PCT_CLAY NOT on surfdata file' ) 

  call ncd_pio_closefile(ncid)

  if (masterproc) then
     write(iulog,*) 'Successfully read fmax, soil color, sand and clay boundary data'
     write(iulog,*)
  endif

  ! Determine saturated and dry soil albedos for n color classes and 
  ! numrad wavebands (1=vis, 2=nir)

  allocate(albsat(mxsoil_color,numrad), albdry(mxsoil_color,numrad), stat=ier)
  if (ier /= 0) then
     write(iulog,*)'iniTimeConst: allocation error for albsat, albdry'
     call endrun()
  end if

  if (mxsoil_color == 8) then
     albsat(1:8,1) = (/0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8/)
     albsat(1:8,2) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
     albdry(1:8,1) = (/0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8/)
     albdry(1:8,2) = (/0.48_r8,0.44_r8,0.40_r8,0.36_r8,0.32_r8,0.28_r8,0.24_r8,0.20_r8/)
  else if (mxsoil_color == 20) then
     albsat(1:20,1) = (/0.25_r8,0.23_r8,0.21_r8,0.20_r8,0.19_r8,0.18_r8,0.17_r8,0.16_r8,&
                        0.15_r8,0.14_r8,0.13_r8,0.12_r8,0.11_r8,0.10_r8,0.09_r8,0.08_r8,0.07_r8,0.06_r8,0.05_r8,0.04_r8/)
     albsat(1:20,2) = (/0.50_r8,0.46_r8,0.42_r8,0.40_r8,0.38_r8,0.36_r8,0.34_r8,0.32_r8,&
                        0.30_r8,0.28_r8,0.26_r8,0.24_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     albdry(1:20,1) = (/0.36_r8,0.34_r8,0.32_r8,0.31_r8,0.30_r8,0.29_r8,0.28_r8,0.27_r8,&
                        0.26_r8,0.25_r8,0.24_r8,0.23_r8,0.22_r8,0.20_r8,0.18_r8,0.16_r8,0.14_r8,0.12_r8,0.10_r8,0.08_r8/)
     albdry(1:20,2) = (/0.61_r8,0.57_r8,0.53_r8,0.51_r8,0.49_r8,0.48_r8,0.45_r8,0.43_r8,&
                        0.41_r8,0.39_r8,0.37_r8,0.35_r8,0.33_r8,0.31_r8,0.29_r8,0.27_r8,0.25_r8,0.23_r8,0.21_r8,0.16_r8/)
  else
     write(iulog,*)'maximum color class = ',mxsoil_color,' is not supported'
     call endrun
  end if
  
  do p = begp,endp
     g = pgridcell(p)
     sandfrac(p) = sand3d(g,1)/100.0_r8
     clayfrac(p) = clay3d(g,1)/100.0_r8
  end do

  ! --------------------------------------------------------------------
  ! If a organic matter dataset has been specified, read it
  ! --------------------------------------------------------------------

  call organicrd(organic3d)

  ! --------------------------------------------------------------------
  ! Initialize time constant arrays of ecophysiological constants and
  ! arrays of dgvm ecophysiological constants
  ! --------------------------------------------------------------------

   do m = 0,numpft
      if (m <= ntree) then
         pftcon%tree(m) = 1
      else
         pftcon%tree(m) = 0
      end if
      pftcon%z0mr(m) = z0mr(m)
      pftcon%displar(m) = displar(m)
      pftcon%dleaf(m) = dleaf(m)
      pftcon%xl(m) = xl(m)
      do ib = 1,numrad
         pftcon%rhol(m,ib) = rhol(m,ib)
         pftcon%rhos(m,ib) = rhos(m,ib)
         pftcon%taul(m,ib) = taul(m,ib)
         pftcon%taus(m,ib) = taus(m,ib)
      end do
      pftcon%qe25(m) = qe25(m)
      pftcon%mp(m) = mp(m)
      pftcon%c3psn(m) = c3psn(m)
      pftcon%slatop(m) = slatop(m)
      pftcon%dsladlai(m) = dsladlai(m)
      pftcon%leafcn(m) = leafcn(m)
      pftcon%flnr(m) = flnr(m)
      pftcon%smpso(m) = smpso(m)
      pftcon%smpsc(m) = smpsc(m)
      pftcon%fnitr(m) = fnitr(m)
      pftcon%woody(m) = woody(m)
      pftcon%lflitcn(m) = lflitcn(m)
      pftcon%frootcn(m) = frootcn(m)
      pftcon%livewdcn(m) = livewdcn(m)
      pftcon%deadwdcn(m) = deadwdcn(m)
      pftcon%graincn(m) = graincn(m)
      pftcon%froot_leaf(m) = froot_leaf(m)
      pftcon%stem_leaf(m) = stem_leaf(m)
      pftcon%croot_stem(m) = croot_stem(m)
      pftcon%flivewd(m) = flivewd(m)
      pftcon%fcur(m) = fcur(m)
      pftcon%lf_flab(m) = lf_flab(m)
      pftcon%lf_fcel(m) = lf_fcel(m)
      pftcon%lf_flig(m) = lf_flig(m)
      pftcon%fr_flab(m) = fr_flab(m)
      pftcon%fr_fcel(m) = fr_fcel(m)
      pftcon%fr_flig(m) = fr_flig(m)
      pftcon%leaf_long(m) = leaf_long(m)
      pftcon%evergreen(m) = evergreen(m)
      pftcon%stress_decid(m) = stress_decid(m)
      pftcon%season_decid(m) = season_decid(m)
      pftcon%resist(m) = resist(m)
      pftcon%dwood(m) = dwood
   end do

   if (use_cndv) then
      do m = 0,numpft
         dgv_pftcon%crownarea_max(m) = pftpar20(m)
         dgv_pftcon%tcmin(m) = pftpar28(m)
         dgv_pftcon%tcmax(m) = pftpar29(m)
         dgv_pftcon%gddmin(m) = pftpar30(m)
         dgv_pftcon%twmax(m) = pftpar31(m)
         dgv_pftcon%reinickerp(m) = reinickerp
         dgv_pftcon%allom1(m) = allom1
         dgv_pftcon%allom2(m) = allom2
         dgv_pftcon%allom3(m) = allom3
         ! modification for shrubs by X.D.Z
         if (m > ntree .and. m <= nbrdlf_dcd_brl_shrub ) then 
            dgv_pftcon%allom1(m) = allom1s
            dgv_pftcon%allom2(m) = allom2s
         end if
      end do
   end if

   ! --------------------------------------------------------------------
   ! Define layer structure for soil, lakes, urban walls and roof 
   ! Vertical profile of snow is not initialized here 
   ! --------------------------------------------------------------------

   ! Lake layers (assumed same for all lake patches)

   dzlak(1) = 0.1_r8
   dzlak(2) = 1._r8
   dzlak(3) = 2._r8
   dzlak(4) = 3._r8
   dzlak(5) = 4._r8
   dzlak(6) = 5._r8
   dzlak(7) = 7._r8
   dzlak(8) = 7._r8
   dzlak(9) = 10.45_r8
   dzlak(10)= 10.45_r8

   zlak(1) =  0.05_r8
   zlak(2) =  0.6_r8
   zlak(3) =  2.1_r8
   zlak(4) =  4.6_r8
   zlak(5) =  8.1_r8
   zlak(6) = 12.6_r8
   zlak(7) = 18.6_r8
   zlak(8) = 25.6_r8
   zlak(9) = 34.325_r8
   zlak(10)= 44.775_r8

   ! Soil layers and interfaces (assumed same for all non-lake patches)
   ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil

   do j = 1, nlevgrnd
      zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
   enddo

   dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
   do j = 2,nlevgrnd-1
      dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
   enddo
   dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)

   zisoi(0) = 0._r8
   do j = 1, nlevgrnd-1
      zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
   enddo
   zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8*dzsoi(nlevgrnd)

   ! Column level initialization for urban wall and roof layers and interfaces
   do l = begl, endl

   ! "0" refers to urban wall/roof surface and "nlevsoi" refers to urban wall/roof bottom
    if (ltype(l)==isturb) then
   
      do j = 1, nlevurb
        zurb_wall(l,j) = (j-0.5)*(thick_wall(l)/float(nlevurb))  !node depths
      end do
      do j = 1, nlevurb
        zurb_roof(l,j) = (j-0.5)*(thick_roof(l)/float(nlevurb))  !node depths
      end do

      dzurb_wall(l,1) = 0.5*(zurb_wall(l,1)+zurb_wall(l,2))    !thickness b/n two interfaces
      do j = 2,nlevurb-1
        dzurb_wall(l,j)= 0.5*(zurb_wall(l,j+1)-zurb_wall(l,j-1)) 
      enddo
      dzurb_wall(l,nlevurb) = zurb_wall(l,nlevurb)-zurb_wall(l,nlevurb-1)

      dzurb_roof(l,1) = 0.5*(zurb_roof(l,1)+zurb_roof(l,2))    !thickness b/n two interfaces
      do j = 2,nlevurb-1
        dzurb_roof(l,j)= 0.5*(zurb_roof(l,j+1)-zurb_roof(l,j-1)) 
      enddo
      dzurb_roof(l,nlevurb) = zurb_roof(l,nlevurb)-zurb_roof(l,nlevurb-1)

      ziurb_wall(l,0) = 0.
      do j = 1, nlevurb-1
        ziurb_wall(l,j) = 0.5*(zurb_wall(l,j)+zurb_wall(l,j+1))          !interface depths
      enddo
      ziurb_wall(l,nlevurb) = zurb_wall(l,nlevurb) + 0.5*dzurb_wall(l,nlevurb)

      ziurb_roof(l,0) = 0.
      do j = 1, nlevurb-1
        ziurb_roof(l,j) = 0.5*(zurb_roof(l,j)+zurb_roof(l,j+1))          !interface depths
      enddo
      ziurb_roof(l,nlevurb) = zurb_roof(l,nlevurb) + 0.5*dzurb_roof(l,nlevurb)
    end if
   end do

   ! Grid level initialization
   do g = begg, endg

      ! VOC emission factors
      ! Set gridcell and landunit indices
      efisop(:,g)=efisop2d(:,g)

   end do


   ! --------------------------------------------------------------------
   ! Initialize soil and lake levels
   ! Initialize soil color, thermal and hydraulic properties
   ! --------------------------------------------------------------------

   ! Column level initialization
   do c = begc, endc

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)
      
      ! initialize maximum daylength, based on latitude and maximum declination
      ! maximum declination hardwired for present-day orbital parameters, 
      ! +/- 23.4667 degrees = +/- 0.409571 radians, use negative value for S. Hem
      max_decl = 0.409571
      if (lat(g) .lt. 0._r8) max_decl = -max_decl
      temp = -(sin(lat(g))*sin(max_decl))/(cos(lat(g)) * cos(max_decl))
      temp = min(1._r8,max(-1._r8,temp))
      max_dayl(c) = 2.0_r8 * 13750.9871_r8 * acos(temp)

      ! Initialize restriction for min of soil potential (mm)
      smpmin(c) = -1.e8_r8

      ! Decay factor (m)
      hkdepth(c) = 1._r8/2.5_r8

      ! Maximum saturated fraction
      wtfact(c) = gti(g)

      ! Soil color
      isoicol(c) = soic2d(g)

      ! Soil hydraulic and thermal properties
        ! Note that urban roof, sunwall and shadewall thermal properties used to 
        ! derive thermal conductivity and heat capacity are set to special 
        ! value because thermal conductivity and heat capacity for urban 
        ! roof, sunwall and shadewall are prescribed in SoilThermProp.F90 in 
        ! SoilTemperatureMod.F90
      if (ltype(l)==istdlak .or. ltype(l)==istwet .or. ltype(l)==istice .or. ltype(l)==istice_mec) then
         do lev = 1,nlevgrnd
            bsw(c,lev)    = spval
            bsw2(c,lev)   = spval
            psisat(c,lev) = spval
            vwcsat(c,lev) = spval
            watsat(c,lev) = spval
            watfc(c,lev)  = spval
            hksat(c,lev)  = spval
            sucsat(c,lev) = spval
            tkmg(c,lev)   = spval
            tksatu(c,lev) = spval
            tkdry(c,lev)  = spval
            if (ltype(l)==istwet .and. lev > nlevsoi) then
               csol(c,lev) = csol_bedrock
            else
               csol(c,lev)= spval
            endif
            watdry(c,lev) = spval 
            watopt(c,lev) = spval 
         end do
      else if (ltype(l)==isturb .and. (ctype(c) /= icol_road_perv) .and. (ctype(c) /= icol_road_imperv) )then
         ! Urban Roof, sunwall, shadewall properties set to special value
         do lev = 1,nlevurb
            watsat(c,lev) = spval
            watfc(c,lev)  = spval
            bsw(c,lev)    = spval
            bsw2(c,lev)   = spval
            psisat(c,lev) = spval
            vwcsat(c,lev) = spval
            hksat(c,lev)  = spval
            sucsat(c,lev) = spval
            tkmg(c,lev)   = spval
            tksatu(c,lev) = spval
            tkdry(c,lev)  = spval
            csol(c,lev)   = spval
            watdry(c,lev) = spval 
            watopt(c,lev) = spval 
         end do
      else  ! soil columns of both urban and non-urban types
         do lev = 1,nlevgrnd
            ! duplicate clay and sand values from 10th soil layer
            if (lev .le. nlevsoi) then
               clay    = clay3d(g,lev)
               sand    = sand3d(g,lev)
               om_frac = (organic3d(g,lev)/organic_max)**2._r8
            else
               clay    = clay3d(g,nlevsoi)
               sand    = sand3d(g,nlevsoi)
               om_frac = 0._r8
            endif
            ! No organic matter for urban
            if (ltype(l)==isturb) then
              om_frac = 0._r8
            end if
            ! Note that the following properties are overwritten for urban impervious road 
            ! layers that are not soil in SoilThermProp.F90 within SoilTemperatureMod.F90
            watsat(c,lev) = 0.489_r8 - 0.00126_r8*sand
            bsw(c,lev)    = 2.91 + 0.159*clay
            sucsat(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
            bd            = (1._r8-watsat(c,lev))*2.7e3_r8 
            watsat(c,lev) = (1._r8 - om_frac)*watsat(c,lev) + om_watsat*om_frac
            tkm           = (1._r8-om_frac)*(8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac ! W/(m K)
            bsw(c,lev)    = (1._r8-om_frac)*(2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
            bsw2(c,lev)   = -(3.10_r8 + 0.157_r8*clay - 0.003_r8*sand)
            psisat(c,lev) = -(exp((1.54_r8 - 0.0095_r8*sand + 0.0063_r8*(100.0_r8-sand-clay))*log(10.0_r8))*9.8e-5_r8)
            vwcsat(c,lev) = (50.5_r8 - 0.142_r8*sand - 0.037_r8*clay)/100.0_r8
            sucsat(c,lev) = (1._r8-om_frac)*sucsat(c,lev) + om_sucsat*om_frac  
            xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

            ! perc_frac is zero unless perf_frac greater than percolation threshold
            if (om_frac > pc) then
               perc_norm=(1._r8 - pc)**(-pcbeta)
               perc_frac=perc_norm*(om_frac - pc)**pcbeta
            else
               perc_frac=0._r8
            endif
            ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
            uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac
            ! uncon_hksat is series addition of mineral/organic conductivites
            if (om_frac .lt. 1._r8) then
              uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                   +((1._r8-perc_frac)*om_frac)/om_hksat)
            else
              uncon_hksat = 0._r8
            end if
            hksat(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

            tkmg(c,lev)   = tkm ** (1._r8- watsat(c,lev))           
            tksatu(c,lev) = tkmg(c,lev)*0.57_r8**watsat(c,lev)
            tkdry(c,lev)  = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd))*(1._r8-om_frac) + &
                            om_tkd*om_frac  
            csol(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) +   &
                           om_csol*om_frac)*1.e6_r8  ! J/(m3 K)  
            if (lev .gt. nlevsoi) then
               csol(c,lev) = csol_bedrock
            endif
            watdry(c,lev) = watsat(c,lev) * (316230._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
            watopt(c,lev) = watsat(c,lev) * (158490._r8/sucsat(c,lev)) ** (-1._r8/bsw(c,lev)) 
            !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
            ! water content at field capacity, defined as hk = 0.1 mm/day
            ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
            watfc(c,lev) = watsat(c,lev) * (0.1_r8 / (hksat(c,lev)*secspday))**(1._r8/(2._r8*bsw(c,lev)+3._r8))
         end do
         !
         ! Urban pervious and impervious road
         !
         ! Impervious road layers -- same as above except set watdry and watopt as missing
         if (ctype(c) == icol_road_imperv) then
            do lev = 1,nlevgrnd
               watdry(c,lev) = spval 
               watopt(c,lev) = spval 
            end do
         ! pervious road layers -- same as above except also set rootfr_road_perv
         ! Currently, pervious road has same properties as soil
         else if (ctype(c) == icol_road_perv) then 
            do lev = 1, nlevgrnd
               rootfr_road_perv(c,lev) = 0._r8
            enddo
            do lev = 1,nlevsoi
               rootfr_road_perv(c,lev) = 0.1_r8  ! uniform profile
            end do
         end if
      endif

      ! Define lake or non-lake levels, layers and interfaces
      if (ltype(l) == istdlak) then
         z(c,1:nlevlak)  = zlak(1:nlevlak)
         dz(c,1:nlevlak) = dzlak(1:nlevlak)
      else if (ltype(l) == isturb) then
         if (ctype(c)==icol_sunwall .or. ctype(c)==icol_shadewall) then
            z(c,1:nlevurb)  = zurb_wall(l,1:nlevurb)
            zi(c,0:nlevurb) = ziurb_wall(l,0:nlevurb)
            dz(c,1:nlevurb) = dzurb_wall(l,1:nlevurb)
         else if (ctype(c)==icol_roof) then
            z(c,1:nlevurb)  = zurb_roof(l,1:nlevurb)
            zi(c,0:nlevurb) = ziurb_roof(l,0:nlevurb)
            dz(c,1:nlevurb) = dzurb_roof(l,1:nlevurb)
         else
            z(c,1:nlevurb)  = zsoi(1:nlevurb)
            zi(c,0:nlevurb) = zisoi(0:nlevurb)
            dz(c,1:nlevurb) = dzsoi(1:nlevurb)
         end if
      else
         z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
         zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
         dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
      end if

      ! Initialize terms needed for dust model
      clay = clay3d(g,1)
      gwc_thr(c) = 0.17_r8 + 0.14_r8*clay*0.01_r8
      mss_frc_cly_vld(c) = min(clay*0.01_r8, 0.20_r8)

   end do

   ! pft level initialization
   do p = begp, endp

      ! Initialize maximum allowed dew

      dewmx(p)  = 0.1_r8

      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).

      c = pcolumn(p)
      if (ivt(p) /= noveg) then
         do lev = 1, nlevgrnd
            rootfr(p,lev) = 0._r8
         enddo
         do lev = 1, nlevsoi-1
            rootfr(p,lev) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,lev-1))  &
                               + exp(-rootb_par(ivt(p)) * zi(c,lev-1))  &
                               - exp(-roota_par(ivt(p)) * zi(c,lev  ))  &
                               - exp(-rootb_par(ivt(p)) * zi(c,lev  )) )
         end do
         rootfr(p,nlevsoi) = .5_r8*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
                                + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )
         rootfr(p,nlevsoi+1:nlevgrnd) =  0.0_r8

!if (use_cn) then
!        ! replacing the exponential rooting distribution
!        ! with a linear decrease, going to zero at the bottom of the lowest
!        ! soil layer for woody pfts, but going to zero at the bottom of
!        ! layer 8 for non-woody pfts.  This corresponds to 3.43 m for woody
!        ! bottom, vs 1.38 m for non-woody bottom.
!        if (woody(ivt(p)) == 1) then
!           bottom = nlevsoi
!           slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._r8/zi(c,bottom)
!           do lev = 1, bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._r8
!              end do
!           end if
!        else
!           bottom = 8
!           slope = -2._r8/(zi(c,bottom)*zi(c,bottom))
!           intercept   = 2._r8/zi(c,bottom)
!           do lev=1,bottom
!              rootfr(p,lev) = dz(c,lev) * 0.5_r8 * ((intercept+slope*zi(c,lev-1)) + (intercept+slope*zi(c,lev)))
!           end do
!           if (bottom < nlevsoi) then
!              do lev=bottom+1,nlevgrnd
!                 rootfr(p,lev) = 0._r8
!              end do
!           end if
!        end if
! endif
      else
         rootfr(p,1:nlevsoi) = 0._r8
      endif
      
      ! initialize rresis, for use in ecosystemdyn
      do lev = 1,nlevgrnd
         rresis(p,lev) = 0._r8
      end do

   end do ! end pft level initialization
   
   if (use_cn) then
      ! initialize the CN variables for special landunits, including lake points
      call CNiniSpecial()
   end if

   deallocate(soic2d,sand3d,clay3d,gti,organic3d)
   deallocate(temp_ef,efisop2d)
   deallocate(zurb_wall, zurb_roof, dzurb_wall, dzurb_roof, ziurb_wall, ziurb_roof)


   ! Initialize SNICAR optical and aging parameters:
   call SnowOptics_init( )

   call SnowAge_init( )

   if (masterproc) write(iulog,*) 'Successfully initialized time invariant variables'

end subroutine iniTimeConst
