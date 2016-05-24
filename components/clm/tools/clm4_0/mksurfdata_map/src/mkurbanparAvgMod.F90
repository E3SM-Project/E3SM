module mkurbanparAvgMod
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mkurbanparAvg
!
! !DESCRIPTION:
! Make Urban Parameter data, using an average of input cells
!
! !REVISION HISTORY:
! Author: Keith Oleson, Mariana Vertenstein
! Feb 2012: Bill Sacks: pulled out some functionality to mkurbanparCommonMod
!
!-----------------------------------------------------------------------
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_sys_mod , only : shr_sys_flush
  use mkdomainMod , only : domain_checksame
  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:
  public  :: mkurban       ! Get the urban percentage
  public  :: mkurbanpar    ! Make the urban parameters

!EOP

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurban
!
! !INTERFACE:
subroutine mkurban(ldomain, mapfname, datfname, ndiag, zero_out, urbn_o)
!
! !DESCRIPTION:
! make percent urban
!
! !USES:
  use mkdomainMod , only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkurbanparCommonMod, only : mkurban_pct, mkurban_pct_diagnostics, MIN_DENS
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*) , intent(in) :: mapfname  ! input mapping file name
  character(len=*) , intent(in) :: datfname  ! input data file name
  integer          , intent(in) :: ndiag     ! unit number for diag out
  logical          , intent(in) :: zero_out  ! if should zero urban out
  real(r8)         , intent(out):: urbn_o(:) ! output grid: %urban
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type) :: tdomain               ! local domain
  type(gridmap_type) :: tgridmap              ! local gridmap
  real(r8), allocatable :: urbn_i(:)          ! input grid: percent urbn
  integer  :: no,ns                           ! indices
  integer  :: ncid,dimid,varid                ! input netCDF id's
  integer  :: ier                             ! error status
  character(len=32) :: subname = 'mkurban'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make %urban .....'

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  call domain_read(tdomain, datfname)

  ns = tdomain%ns
  allocate(urbn_i(ns), stat=ier)
  if (ier/=0) call abort()

  write (6,*) 'Open urban file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncid), subname)
  call check_ret(nf_inq_varid (ncid, 'PCT_URBAN', varid), subname)
  call check_ret(nf_get_var_double (ncid, varid, urbn_i), subname)
  call check_ret(nf_close(ncid), subname)

  ! Compute local fields _o

  if ( all_urban .or. zero_out )then

     do no = 1, ldomain%ns
        if (all_urban )then
           urbn_o(no) = 100.00_r8
        else 
           urbn_o(no) =   0.00_r8
        end if
     enddo

  else

     call gridmap_mapread(tgridmap, mapfname)

     call mkurban_pct(ldomain, tdomain, tgridmap, urbn_i, urbn_o)

     do no = 1, ldomain%ns
        if (urbn_o(no) < MIN_DENS) then
           urbn_o(no) = 0._r8
        end if
     end do

     call mkurban_pct_diagnostics(ldomain, tdomain, tgridmap, urbn_i, urbn_o, ndiag)

     write (6,*) 'Successfully made %urban'
     write (6,*)

  end if

  ! Deallocate dynamic memory

  !  call domain_clean(tdomain)
  if ( .not. all_urban .and. .not. zero_out )then
     call gridmap_clean(tgridmap)
  end if

  deallocate (urbn_i)

end subroutine mkurban

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkurbanpar
!
! !INTERFACE:
subroutine mkurbanpar(ldomain, mapfname, datfname, ndiag, ncido)
!
! !DESCRIPTION:
! Make Urban Parameter data
!
! !USES:
  use mkdomainMod  , only : domain_type, domain_clean, domain_read
  use mkgridmapMod
  use mkvarpar	
  use mkvarctl    
  use mkncdio
!
! !ARGUMENTS:
  implicit none
  type(domain_type), intent(in) :: ldomain
  character(len=*)  , intent(in) :: mapfname  ! input mapping file name
  character(len=*)  , intent(in) :: datfname  ! input data file name
  integer           , intent(in) :: ndiag     ! unit number for diag out
  integer           , intent(in) :: ncido     ! output netcdf file id
!
! !CALLED FROM:
! subroutine mksrfdat in module mksrfdatMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
  type(domain_type)    :: tdomain              ! local domain 
  type(domain_type)    :: tdomain_mask         ! local domain that contains "mask"
  type(gridmap_type)    :: tgridmap             ! local gridmap
  real(r8), allocatable :: canyon_hwr_o(:)      ! canyon height to width ratio out
  real(r8), allocatable :: canyon_hwr_i(:)      ! canyon_height to width ratio in
  real(r8), allocatable :: em_improad_o(:)      ! emissivity of impervious road out
  real(r8), allocatable :: em_improad_i(:)      ! emissivity of impervious road in
  real(r8), allocatable :: em_perroad_o(:)      ! emissivity of pervious road out
  real(r8), allocatable :: em_perroad_i(:)      ! emissivity of pervious road in
  real(r8), allocatable :: em_roof_o(:)         ! emissivity of roof out
  real(r8), allocatable :: em_roof_i(:)         ! emissivity of roof in
  real(r8), allocatable :: em_wall_o(:)         ! emissivity of wall out
  real(r8), allocatable :: em_wall_i(:)         ! emissivity of wall in
  real(r8), allocatable :: ht_roof_o(:)         ! height of roof out
  real(r8), allocatable :: ht_roof_i(:)         ! height of roof in
  real(r8), allocatable :: thick_roof_o(:)      ! thickness of roof out
  real(r8), allocatable :: thick_roof_i(:)      ! thickness of roof in
  real(r8), allocatable :: thick_wall_o(:)      ! thickness of wall out
  real(r8), allocatable :: thick_wall_i(:)      ! thickness of wall in
  real(r8), allocatable :: t_building_max_o(:)  ! maximum interior building temperature out
  real(r8), allocatable :: t_building_max_i(:)  ! maximum interior building temperature in
  real(r8), allocatable :: t_building_min_o(:)  ! minimum interior building temperature out
  real(r8), allocatable :: t_building_min_i(:)  ! minimum interior building temperature in
  real(r8), allocatable :: wind_hgt_canyon_o(:) ! height of wind in canyon out
  real(r8), allocatable :: wind_hgt_canyon_i(:) ! height of wind in canyon in
  real(r8), allocatable :: wtlunit_roof_o(:)    ! fraction of roof out
  real(r8), allocatable :: wtlunit_roof_i(:)    ! fraction of roof in
  real(r8), allocatable :: wtroad_perv_o(:)     ! fraction of pervious road out
  real(r8), allocatable :: wtroad_perv_i(:)     ! fraction of pervious road in
  real(r8), allocatable :: alb_improad_o(:,:,:) ! albedo of impervious road out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_improad_i(:,:,:) ! albedo of impervious road in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_o(:,:,:) ! albedo of pervious road out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_perroad_i(:,:,:) ! albedo of pervious road in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_o(:,:,:)    ! albedo of roof out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_roof_i(:,:,:)    ! albedo of roof in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_o(:,:,:)    ! albedo of wall out (lon,lat,numrad,numsolar)
  real(r8), allocatable :: alb_wall_i(:,:,:)    ! albedo of wall in (lon,lat,numrad,numsolar)
  real(r8), allocatable :: tk_roof_o(:,:)       ! thermal conductivity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_roof_i(:,:)       ! thermal conductivity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_o(:,:)       ! thermal conductivity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_wall_i(:,:)       ! thermal conductivity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_o(:,:)    ! thermal conductivity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: tk_improad_i(:,:)    ! thermal conductivity of impervious road in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_o(:,:)       ! volumetric heat capacity of roof out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_roof_i(:,:)       ! volumetric heat capacity of roof in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_o(:,:)       ! volumetric heat capacity of wall out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_wall_i(:,:)       ! volumetric heat capacity of wall in (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_o(:,:)    ! volumetric heat capacity of impervious road out (lon,lat,nlevurb)
  real(r8), allocatable :: cv_improad_i(:,:)    ! volumetric heat capacity of impervious road in (lon,lat,nlevurb)
  integer,  allocatable :: nlev_improad_o(:)    ! number of impervious road layers out
  real(r8), allocatable :: mask_i(:)            ! input grid: mask (0, 1)
  integer  :: nlevurb_i                         ! input  grid: number of urban vertical levels
  integer  :: numsolar_i                        ! input  grid: number of solar type (DIR/DIF)
  integer  :: numrad_i                          ! input  grid: number of solar bands (VIS/NIR)
  integer  :: ns_i,ns_o                         ! indices
  integer  :: k,l,n,m,ni,no                     ! indices
  integer  :: nsolar,nrad,nurb                  ! indices
  integer  :: ncidi,dimid,varid                 ! input netCDF id's
  integer  :: numlev                            ! number of valid impervious road layers
  integer  :: ier                               ! error status
  real(r8) :: relerr = 0.00001                  ! max error: sum overlap wts ne 1
  character(len=256) :: name                    ! name of attribute
  character(len=256) :: unit                    ! units of attribute
  character(len= 32) :: subname = 'mkurbanpar'
!-----------------------------------------------------------------------

  write (6,*) 'Attempting to make Urban Parameters .....'
  call shr_sys_flush(6)

  ! -----------------------------------------------------------------
  ! Read input file
  ! -----------------------------------------------------------------

  ! Obtain input grid info, read local fields

  !***NOTE - this determines the tdomain mask based on the "mask"
  ! file variable NOT the LANDMASK file variable - and the "mask"
  ! variable is a SUBSET of the LANDMASK variable - so we need
  ! to call the interplation with the mask_src optional argument 

  call domain_read(tdomain     , datfname)
  call domain_read(tdomain_mask, datfname, readmask=.true.)

  write (6,*) 'Open urban parameter file: ', trim(datfname)
  call check_ret(nf_open(datfname, 0, ncidi), subname)
  call check_ret(nf_inq_dimid(ncidi, 'nlevurb', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, nlevurb_i), subname)
  call check_ret(nf_inq_dimid(ncidi, 'numsolar', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numsolar_i), subname)
  call check_ret(nf_inq_dimid(ncidi, 'numrad', dimid), subname)
  call check_ret(nf_inq_dimlen(ncidi, dimid, numrad_i), subname)

  if (nlevurb_i /= nlevurb) then
     write(6,*)'MKURBANPAR: parameter nlevurb= ',nlevurb, &
          'does not equal input dataset nlevurb= ',nlevurb_i
     stop
  endif
  if (numsolar_i /= numsolar) then
     write(6,*)'MKURBANPAR: parameter numsolar= ',numsolar, &
          'does not equal input dataset numsolar= ',numsolar_i
     stop
  endif
  if (numrad_i /= numrad) then
     write(6,*)'MKURBANPAR: parameter numrad= ',numrad, &
          'does not equal input dataset numrad= ',numrad_i
     stop
  endif

  if (tdomain%ns /= tdomain_mask%ns) then
     write(6,*)'MKURBANPAR: inconsistent sizes for tdomain and tdomain_mask'
     write(6,*)' domain size tdomain     = ',tdomain%ns
     write(6,*)' domain size tdomain_mask= ',tdomain_mask%ns
     stop
  else
     do n = 1,tdomain%ns
        if (tdomain%mask(n) == 0 .and. tdomain_mask%mask(n) == 1) then 
           write(6,*)'tdomain_mask is not a submask tdomain at n= ',n
           stop
        end if
     end do
  end if

  ! Allocation

  ns_i = tdomain%ns
  allocate(canyon_hwr_i(ns_i), &
           em_improad_i(ns_i), &
           em_perroad_i(ns_i), &
           em_roof_i(ns_i), &
           em_wall_i(ns_i), &
           ht_roof_i(ns_i), &
           thick_roof_i(ns_i), &
           thick_wall_i(ns_i), &
           t_building_max_i(ns_i), &
           t_building_min_i(ns_i), &
           wind_hgt_canyon_i(ns_i), &
           wtlunit_roof_i(ns_i), &
           wtroad_perv_i(ns_i), &
           alb_improad_i(ns_i,numrad_i,numsolar_i), &
           alb_perroad_i(ns_i,numrad_i,numsolar_i), &
           alb_roof_i(ns_i,numrad_i,numsolar_i), &
           alb_wall_i(ns_i,numrad_i,numsolar_i), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ns_o = ldomain%ns
  allocate(canyon_hwr_o(ns_o), &
           em_improad_o(ns_o), &
           em_perroad_o(ns_o), &
           em_roof_o(ns_o), &
           em_wall_o(ns_o), &
           ht_roof_o(ns_o), &
           thick_roof_o(ns_o), &
           thick_wall_o(ns_o), &
           t_building_max_o(ns_o), &
           t_building_min_o(ns_o), &
           wind_hgt_canyon_o(ns_o), &
           wtlunit_roof_o(ns_o), &
           wtroad_perv_o(ns_o), &
           alb_improad_o(ns_o,numrad,numsolar), &
           alb_perroad_o(ns_o,numrad,numsolar), &
           alb_roof_o(ns_o,numrad,numsolar), &
           alb_wall_o(ns_o,numrad,numsolar), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ! Compute local fields _o
  ! Area average and then deallocate input data

  call gridmap_mapread(tgridmap, mapfname)

  ! Error checks for domain and map consistencies

  call domain_checksame( tdomain, ldomain, tgridmap )
     
  ! Determine urban variables on output grid

  ! IMPORTANT - first create a mask for the mapping that is based on 
  ! tdomain_mask%mask rather than tdomain%mask

  allocate(mask_i(ns_i))
  do n = 1,tgridmap%ns
     ni = tgridmap%src_indx(n)
     if (tdomain_mask%mask(ni) .eq. 0.) then
        mask_i(ni) = 0.
     else
        mask_i(ni) = 1.
     end if
  end do

  call check_ret(nf_inq_varid (ncidi, 'CANYON_HWR', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, canyon_hwr_i), subname)
  call gridmap_areaave(tgridmap,canyon_hwr_i,canyon_hwr_o,mask_src=mask_i)
  deallocate (canyon_hwr_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_improad_i), subname)
  call gridmap_areaave(tgridmap,em_improad_i,em_improad_o,mask_src=mask_i)
  deallocate (em_improad_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_PERROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_perroad_i), subname)
  call gridmap_areaave(tgridmap,em_perroad_i,em_perroad_o,mask_src=mask_i)
  deallocate (em_perroad_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_roof_i), subname)
  call gridmap_areaave(tgridmap,em_roof_i,em_roof_o,mask_src=mask_i)
  deallocate (em_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'EM_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, em_wall_i), subname)
  call gridmap_areaave(tgridmap,em_wall_i,em_wall_o,mask_src=mask_i)
  deallocate (em_wall_i)

  call check_ret(nf_inq_varid (ncidi, 'HT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, ht_roof_i), subname)
  call gridmap_areaave(tgridmap,ht_roof_i,ht_roof_o,mask_src=mask_i)
  deallocate (ht_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'THICK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_roof_i), subname)
  call gridmap_areaave(tgridmap,thick_roof_i,thick_roof_o,mask_src=mask_i)
  deallocate (thick_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'THICK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, thick_wall_i), subname)
  call gridmap_areaave(tgridmap,thick_wall_i,thick_wall_o,mask_src=mask_i)
  deallocate (thick_wall_i)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_max_i), subname)
  call gridmap_areaave(tgridmap,t_building_max_i,t_building_max_o,mask_src=mask_i)
  deallocate (t_building_max_i)

  call check_ret(nf_inq_varid (ncidi, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, t_building_min_i), subname)
  call gridmap_areaave(tgridmap,t_building_min_i,t_building_min_o,mask_src=mask_i)
  deallocate (t_building_min_i)

  call check_ret(nf_inq_varid (ncidi, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wind_hgt_canyon_i), subname)
  call gridmap_areaave(tgridmap,wind_hgt_canyon_i,wind_hgt_canyon_o,mask_src=mask_i)
  deallocate (wind_hgt_canyon_i)

  call check_ret(nf_inq_varid (ncidi, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtlunit_roof_i), subname)
  call gridmap_areaave(tgridmap,wtlunit_roof_i,wtlunit_roof_o,mask_src=mask_i)
  deallocate (wtlunit_roof_i)

  call check_ret(nf_inq_varid (ncidi, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, wtroad_perv_i), subname)
  call gridmap_areaave(tgridmap,wtroad_perv_i,wtroad_perv_o,mask_src=mask_i)
  deallocate (wtroad_perv_i)

  do n = 1,numrad_i
  do m = 1,numsolar_i     
     call check_ret(nf_inq_varid (ncidi, 'ALB_IMPROAD', varid), subname)
     call check_ret(nf_get_var_double (ncidi, varid, alb_improad_i), subname)
     call gridmap_areaave(tgridmap, alb_improad_i(:,n,m), alb_improad_o(:,n,m),&
          mask_src=mask_i)

     call check_ret(nf_inq_varid (ncidi, 'ALB_PERROAD', varid), subname)
     call check_ret(nf_get_var_double (ncidi, varid, alb_perroad_i), subname)
     call gridmap_areaave(tgridmap,alb_perroad_i(:,n,m),alb_perroad_o(:,n,m),&
          mask_src=mask_i)

     call check_ret(nf_inq_varid (ncidi, 'ALB_ROOF', varid), subname)
     call check_ret(nf_get_var_double (ncidi, varid, alb_roof_i), subname)
     call gridmap_areaave(tgridmap,alb_roof_i(:,n,m),alb_roof_o(:,n,m),&
          mask_src=mask_i)
     
     call check_ret(nf_inq_varid (ncidi, 'ALB_WALL', varid), subname)
     call check_ret(nf_get_var_double (ncidi, varid, alb_wall_i), subname)
     call gridmap_areaave(tgridmap,alb_wall_i(:,n,m),alb_wall_o(:,n,m),&
          mask_src=mask_i)
  end do
  end do
  deallocate (alb_improad_i)
  deallocate (alb_perroad_i)
  deallocate (alb_roof_i)
  deallocate (alb_wall_i)

  ! Now write output data to the file and then deallocate
  call check_ret(nf_inq_varid(ncido, 'CANYON_HWR', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, canyon_hwr_o),subname)
  deallocate (canyon_hwr_o)

  call check_ret(nf_inq_varid(ncido, 'EM_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_improad_o), subname)
  deallocate (em_improad_o)

  call check_ret(nf_inq_varid(ncido, 'EM_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_perroad_o), subname)
  deallocate (em_perroad_o)

  call check_ret(nf_inq_varid(ncido, 'EM_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_roof_o), subname)
  deallocate (em_roof_o)

  call check_ret(nf_inq_varid(ncido, 'EM_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, em_wall_o), subname)
  deallocate (em_wall_o)

  call check_ret(nf_inq_varid(ncido, 'HT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, ht_roof_o), subname)
  deallocate (ht_roof_o)

  call check_ret(nf_inq_varid(ncido, 'THICK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_roof_o), subname)
  deallocate (thick_roof_o)

  call check_ret(nf_inq_varid(ncido, 'THICK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, thick_wall_o), subname)
  deallocate (thick_wall_o)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MAX', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_max_o), subname)
  deallocate (t_building_max_o)

  call check_ret(nf_inq_varid(ncido, 'T_BUILDING_MIN', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, t_building_min_o), subname)
  deallocate (t_building_min_o)

  call check_ret(nf_inq_varid(ncido, 'WIND_HGT_CANYON', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wind_hgt_canyon_o), subname)
  deallocate (wind_hgt_canyon_o)

  call check_ret(nf_inq_varid(ncido, 'WTLUNIT_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtlunit_roof_o), subname)
  deallocate (wtlunit_roof_o)

  call check_ret(nf_inq_varid(ncido, 'WTROAD_PERV', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, wtroad_perv_o), subname)
  deallocate (wtroad_perv_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_improad_o), subname)
  deallocate (alb_improad_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_PERROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_perroad_o), subname)
  deallocate (alb_perroad_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_roof_o), subname)
  deallocate (alb_roof_o)

  call check_ret(nf_inq_varid(ncido, 'ALB_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, alb_wall_o), subname)
  deallocate (alb_wall_o)
  !
  ! 3D nlevurb fields
  !
  ! First allocate data
  allocate(cv_improad_i(ns_i,nlevurb), &
           tk_roof_i(ns_i,nlevurb), &
           tk_wall_i(ns_i,nlevurb), &
           tk_improad_i(ns_i,nlevurb), &
           cv_roof_i(ns_i,nlevurb), &
           cv_wall_i(ns_i,nlevurb), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if
           
  allocate(tk_roof_o(ns_o,nlevurb), &
           tk_wall_o(ns_o,nlevurb), &
           tk_improad_o(ns_o,nlevurb), &
           cv_roof_o(ns_o,nlevurb), &
           cv_wall_o(ns_o,nlevurb), &
           cv_improad_o(ns_o,nlevurb), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  ! Do the areaaveraging and then deallocate input data

  call check_ret(nf_inq_varid (ncidi, 'TK_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_roof_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'TK_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_wall_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_ROOF', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_roof_i), subname)
  call check_ret(nf_inq_varid (ncidi, 'CV_WALL', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_wall_i), subname)

  do n = 1,nlevurb
     call gridmap_areaave(tgridmap,tk_roof_i(:,n),tk_roof_o(:,n),&
          mask_src=mask_i)
     call gridmap_areaave(tgridmap,tk_wall_i(:,n),tk_wall_o(:,n),&
          mask_src=mask_i)
     call gridmap_areaave(tgridmap,cv_roof_i(:,n),cv_roof_o(:,n),&
          mask_src=mask_i)
     call gridmap_areaave(tgridmap,cv_wall_i(:,n),cv_wall_o(:,n),&
          mask_src=mask_i)
  end do

  deallocate (tk_roof_i)
  deallocate (tk_wall_i)
  deallocate (cv_roof_i)
  deallocate (cv_wall_i)
  deallocate (mask_i)

  ! Write output data then deallocate output data
  call check_ret(nf_inq_varid(ncido, 'TK_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_wall_o), subname)
  deallocate (tk_wall_o)

  call check_ret(nf_inq_varid(ncido, 'TK_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_roof_o), subname)
  deallocate (tk_roof_o)

  call check_ret(nf_inq_varid(ncido, 'CV_WALL', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_wall_o), subname)
  deallocate (cv_wall_o)

  call check_ret(nf_inq_varid(ncido, 'CV_ROOF', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_roof_o), subname)
  deallocate (cv_roof_o)

  ! Get fields from input file
  call check_ret(nf_inq_varid (ncidi, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, cv_improad_i), subname)

  call check_ret(nf_inq_varid (ncidi, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_get_var_double (ncidi, varid, tk_improad_i), subname)

  ! Impervious road thermal conductivity and heat capacity need to be
  ! handled differently because of varying levels of data.

  allocate(mask_i(ns_i))
  do nurb = 1,nlevurb
      ! Create mask for input data from missing values
      do n = 1,tgridmap%ns
         ni = tgridmap%src_indx(n)
         if (tk_improad_i(ni,nurb) .eq. -999.) then
            mask_i(ni) = 0.
         else
            mask_i(ni) = 1.
         end if
      end do
      call gridmap_areaave(tgridmap, tk_improad_i(:,nurb), tk_improad_o(:,nurb), &
           mask_src=mask_i)
      call gridmap_areaave(tgridmap, cv_improad_i(:,nurb), cv_improad_o(:,nurb), &
           mask_src=mask_i)
  end do
  deallocate(cv_improad_i)
  deallocate(tk_improad_i)
  deallocate(mask_i)

  call check_ret(nf_inq_varid(ncido, 'TK_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, tk_improad_o), subname)

  call check_ret(nf_inq_varid(ncido, 'CV_IMPROAD', varid), subname)
  call check_ret(nf_put_var_double(ncido, varid, cv_improad_o), subname)

  allocate(nlev_improad_o(ns_o), stat=ier)
  if (ier /= 0) then
     write(6,*)'mkurbanpar allocation error'; call abort()
  end if

  nlev_improad_o(:) = 0
  do no = 1,ns_o
     numlev = 0
     do nurb = 1,nlevurb
        if (tk_improad_o(no,nurb) > 0. .and. cv_improad_o(no,nurb) > 0.) then
           numlev = numlev+1 
        end if
     end do
     nlev_improad_o(no) = numlev
  end do

  call check_ret(nf_inq_varid(ncido, 'NLEV_IMPROAD', varid), subname)
  call check_ret(nf_put_var_int(ncido, varid, nlev_improad_o), subname)

  call check_ret(nf_sync(ncido), subname)
  call check_ret(nf_close(ncidi), subname)

  ! Deallocate dynamic memory

  deallocate (tk_improad_o)
  deallocate (cv_improad_o)
  deallocate (nlev_improad_o)
  call domain_clean(tdomain)
  call domain_clean(tdomain_mask)
  call gridmap_clean(tgridmap)

  write (6,*) 'Successfully made Urban Parameters'
  write (6,*)
  call shr_sys_flush(6)

end subroutine mkurbanpar

!-----------------------------------------------------------------------

end module mkurbanparAvgMod
