#define DEBUG
Module gcam2emisfile_mod

  !---------------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: gcam2emisfile_mod
  !
  !  
  !
  ! !DESCRIPTION:  The routine performs downscaling of the regional GCAM
  !                CO2 emissions as well as updating a dynamic CO2 flux
  !                file with the downscaled fluxes.  The CO2 boundary file
  !                is used as boundry data for CAM.  This routine is based
  !                on downscaling code provided by Yuyu Zhou at PNNL.
  !
  ! !USES:

  use iac_fields_mod
  use shr_file_mod, only: shr_file_getunit, shr_file_freeunit
  use shr_cal_mod
  use shr_sys_mod
  use netcdf
  use shr_kind_mod, only : r8 => shr_kind_r8,r4 => shr_kind_r4
  use shr_const_mod, only: pi=>shr_const_pi,re=>shr_const_rearth

  implicit none
  SAVE
  private                              ! By default make data private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: gcam2emisfile_init_mod               ! clm initialization
  public :: gcam2emisfile_run_mod                ! clm run phase
  public :: gcam2emisfile_final_mod              ! clm finalization/cleanup
  public :: handle_ncerr
  public :: interp2d
  public :: cremapbin
  ! !PUBLIC DATA MEMBERS: 


  real(r8), parameter ::  d2r  = pi/180._r8               ! radians to degrees
  real(r8), parameter ::  r2d  = 180._r8/pi               ! degrees to radians 
  real(r8), parameter ::  twopi=pi*2._r8
  real(r8), parameter ::  zero=0._r8
  real(r8), parameter ::  mssng=1.e36

  character(len=512)  :: &
       base2000yrco2file,&
       grid720x360,      &
       grid288x192,      &
       dynco2fluxfile,   &
       ship2000baseco2

  character*128 ::                   &
       emisname,gcam2emis_downscale, &
       gcam2emis_co2_allsteps,       &
       gcam2emis_co2base2000,        &
       gcam_lut_720x360_mapping

  character*15 ::                   &
       regionname(14)

  character*3         :: &
       source(10)

  integer             :: &
       cid,              &
       column,           &
       country,          &
       cum_day,          &
       date1(1),      &
       datesec1(1),      &
       dimid,            &
       dom(12),          &
       doy(12),          &
       dpm(12),          &
       latdimid,         &
       latind,           &
       latshipdimid,     &
       line,             &
       londimid,         &
       lonind,           &
       lonshipdimid,     &
       nc,               &
       ncidallsteps,     &
       ncidco2,          &
       ncidco2base,      &
       ncidds,           &
       ncidlatlon,       &
       ncidmap,          &
       ncidship,         &
       ncidtmp,          &
       ng,               &
       ns,               &
       ntime,            &
       numco2basetime,   &
       numco2lat,        &
       numco2lon,        &
       numcountry,       &
       numcountrygridpts,&
       numcountrymappts, &
       numgas,           &
       numlat,           &
       numlon,           &
       numreg,           &
       numregion,        &
       numsector,        &
       numshipco2lat,    &
       numshipco2lon,    &
       numshipco2time,   &
       numshiplat,       &
       numshiplon,       &
       numsource,        &
       numtime,          &
       numtmplat,        &
       numtmplon,        &
       nlatbaseremap,    &
       nlonbaseremap,    &
       numtmptime,       &
       position,         &
       region,           &
       regionind,        &
       rid,              &
       status,           &
       timedimid,        &
       timeshipdimid,    &
       tmpdimid,         &
       varid,            &
       varids(10),       &
       varidstest(10),   &
       year

  integer, allocatable, dimension(:)  :: &
       Grid05ID,        &
       Country05ID,     &
       CountryID,       &
       regionid,        &
       years

  logical             :: &
       masterproc = .true.

  real(r8)            :: &
       area,                        &
       gcamship2005,                &
       latbaseremap(360),           &
       lonbaseremap(720),           &
       lonp180(720),                &
       scalefactor(360,180),        &
       scalefactor09x125(288,192),  &
       scalefactorhalfdeg(720,360)

  real(r8), pointer ::              &
       area09x125(:,:),             &
       area720x360(:,:)

  real(r8), allocatable, dimension(:) :: &
       co2gridcarlat,                  &
       co2gridcarlon,                  &
       co2lat,                      &
       co2lon,                      &
       shiplat,                     &
       shiplon

  real(r8), allocatable, dimension(:,:,:) :: &
       co2ship2000base,             &
       co2ship2005base,             &
       co2ship2005basekgms,         &
       co2ship2000basehalfdeg,      &
       co2ship2000basehalfdegkggrdmth, &
       ship_future,                 &
       ship_season

  real(r8), allocatable, dimension(:,:,:) :: emiss_prev,emiss_new,emiss_newCO2,emiss_future
  real(r8), allocatable, dimension(:,:) :: co2_interp,emisgcam_fluxfilegrid,ship_interp

  real(r8) ::                       &
       emissionCheck,               &
       globalTotalAfter,            &
       globalTotalBefore,           &
       temp,                        &
       totalArea,                   &
       totalBaseEmissions

  real(r8), allocatable, dimension(:,:)  :: &
       Diff,                        &
       GDP_C_2year,                 &
       GDP_C_Final,                 &
       GDP_C_Preliminary,           &
       GDP_Increase,                &
       GDP_R,                       &
       GDP_R_GCAM,                  &
       GDP_R_GCAM_pred,             &
       GDP_R_IEA,                   &
       GHG_C2000,                   &
       GHG_R_GCAM_allsteps,         &
       tmp_GHG_R_GCAM_allsteps,     &
       POP_R_GCAM,                  &
       Pop_C,                       &
       Pop_R_IIASA,                 &
       annual_avg,                  &
       countryEmission,             &
       countryEmission2,            &
       emis_future_grid_rcp45,      &
       tmp1,      &
       emis_future_grid_rcp45p180  

  real(r8), allocatable, dimension(:,:,:) :: &
       co2base2000,                 &
       co2cam2000base,              &
       co2cam2000basehalfdeg,       &
       co2cam2000basehalfdegtg,     &
       emiss_season,                &
       share,                       &
       tmpco2

  real(r8), allocatable, dimension(:) :: &
       Area_sumGrid,                &
       Diff2,                       &
       GDPpc_Grow,                  &
       GHG_C_Final,                 &
       GHG_C_Preliminary,           &
       GHG_R_GCAM,                  &
       GHG_R_Pre,                   &
       GHG_Share_R,                 &
       GHGpg_Grow,                  &
       areaCountry,                 &
       areaGrid,                    &
       areaIntersect,               &
       emissionInCountryIn12year,   &
       totalBaseEmissionsCountry
  
  real(r8)  :: regiontot(14)

  integer, allocatable, dimension(:) :: &
       country2region

  ! !REVISION HISTORY:
  ! Author: J Truesdale


  ! !PRIVATE DATA MEMBERS:

  !EOP
  !===============================================================
contains
  !===============================================================

  !---------------------------------------------------------------------------
  !BOP

  ! !IROUTINE: gcam2emisfile_init_mod

  ! !INTERFACE:
  subroutine gcam2emisfile_init_mod( EClock, cdata, gcamoemis)

    ! !DESCRIPTION:
    ! Initialize interface for gcam2emisfile

    ! !USES:
    implicit none

    ! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real(r8), pointer :: gcamoemis(:,:)


    ! !LOCAL VARIABLES:
    logical :: restart_run,lexist
    logical :: initial_run
    integer :: iu,iun,tmpyears(2),tmp(1),ii,n,nn,maxlonlat(2)
    integer       ::i,j,idx2000,idx2005,idx2080,idx2090,idx2100,mth
    integer, dimension(2) :: start2,count2
    integer, dimension(3) :: start3,count3
    !jt    real(r4) :: annual_avg
    character(len=*),parameter :: subname='(gcam2emisfile_init_mod)'
    character(len=512) :: gcam2emis_co2base2000
    character(len=512) :: gcam2emis_halfdeg_mapping
    integer, allocatable, dimension(:) :: &
         datearr
    integer            ::skip,endidx
    ! !REVISION HISTORY:
    ! Author: JET          ! rewrite of matlat preprocessing script new_grids_matchforest4.m

    !EOP
    !-----------------------------------------------------------------------

    regionname=(/  &
         'USA            ', &
         'Canada         ', &
         'Western Europe ', &
         'Japan          ', &
         'Australia_NZ   ', &
         'USSR           ', &
         'China          ', &
         'Middle East    ', &
         'Africa         ', &
         'Latin America  ', &
         'South East Asia', &
         'Eastern Europe ', &
         'Korea          ', &
         'India          '/)

    iu  = cdata%i(iac_cdatai_logunit)
#ifdef DEBUG
    write(iu,*) subname,' starting subroutine '
#endif

    restart_run  = cdata%l(iac_cdatal_rest)
    !    gcamsize = cdata%i(iac_cdatai_gcam_naez)*cdata%i(iac_cdatai_gcam_nreg)
    initial_run = cdata%l(iac_cdatal_initrun)
    !    gcam2emis_co2base2000 = trim(cdata%c(iac_cdatac_gcam2emis_co2base2000))
    !    gcam2emis_halfdeg_mapping = trim(cdata%c(iac_cdatac_gcam2emis_lut_0.5deg_mapping))

    base2000yrco2file=trim(cdata%c(iac_cdatac_gcam2emisfile_co2base2000))
    grid720x360=trim(cdata%c(iac_cdatac_gcam2emisfile_grid720x360))
    grid288x192=trim(cdata%c(iac_cdatac_gcam2emisfile_grid288x192))
    ship2000baseco2=trim(cdata%c(iac_cdatac_gcam2emisfile_co2shipbase2000))
    gcam_lut_720x360_mapping=trim(cdata%c(iac_cdatac_gcam2emisfile_lut720x360map))
    gcam2emis_downscale=trim(cdata%c(iac_cdatac_gcam2emisfile_downscaleinfo))
    gcam2emis_co2_allsteps=trim(cdata%c(iac_cdatac_gcam2emisfile_rcp45allsteps))

    !
    ! setup 
    !
    source(:)= (/"agr","awb","dom","ene","ind","lcf","sav","slv","tra","wst"/)

    !coming from gcam plus air and ship source(:)= (/"agr","air","awb","dom","ene","ind","lcf","sav","ship","slv","tra","wst"/)
    !jt    numsource=10 for now jut using co2 aggregated.

    numsource=1
    numcountry=231

    !-------------------------------------------------------
    ! time setup
    dpm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) ! days per month
    
    dom = int(dpm(:) / 2.0)
    doy(1) = int(dpm(1) / 2.0)
    cum_day = dpm(1)
    do mth=2,12
       doy(mth) = cum_day + int(dpm(mth) / 2.0)
       cum_day = cum_day + dpm(mth)
    end do
    
    !-------------------------------------------------------
    ! create datesec
    
    datesec1(1) = 0
    
    allocate (area720x360(720,360))
    allocate (area09x125(288,192))
    allocate (co2base2000(720,360,numsource))

    call handle_ncerr(nf90_open(grid720x360,nf90_nowrite,ncidtmp),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidtmp, "cell_area", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidtmp,varid,area720x360),subname,__LINE__)
    call handle_ncerr( nf90_close(ncidtmp),subname,__LINE__)

    call handle_ncerr(nf90_open(grid288x192,nf90_nowrite,ncidtmp),subname,__LINE__)
    call handle_ncerr( nf90_inq_dimid( ncidtmp,  'lon', tmpdimid ),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidtmp, tmpdimid, len = numtmplon),subname,__LINE__)
    call handle_ncerr( nf90_inq_dimid( ncidtmp,  'lat', tmpdimid ),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidtmp, tmpdimid, len = numtmplat),subname,__LINE__)

    start2(1)=1
    count2(1)=numtmplon
    start2(2)=1
    count2(2)=numtmplat
    call handle_ncerr(nf90_inq_varid(ncidtmp, "cell_area", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidtmp,varid,area09x125,start2,count2),subname,__LINE__)
    call handle_ncerr( nf90_close(ncidtmp),subname,__LINE__)


    ! calculate seasonal average of co2 using last 10 years of data (1995-2004)
    !
    ! Open and read the co2base emissions file
    !                                                                                                                                              
    !    dynco2fluxfile='/project/amp/jet/inputdata/atm/cam/ggas/co2flux_fossil_1751-2006-monthly_0.9x1.25_c20100204.nc'
    dynco2fluxfile='./co2flux_iESM_dyn.nc'
    call handle_ncerr(nf90_open(dynco2fluxfile,nf90_nowrite,ncidco2),subname,__LINE__)
    call handle_ncerr( nf90_inq_dimid( ncidco2,  'lon', londimid ),subname,__LINE__)
    call handle_ncerr( nf90_inq_dimid( ncidco2,  'lat', latdimid ),subname,__LINE__)
    call handle_ncerr( nf90_inq_dimid( ncidco2,  'time', timedimid ),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidco2, latdimid, len = numco2lat),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidco2, londimid, len = numco2lon),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidco2, timedimid, len = numtime),subname,__LINE__)

    allocate (co2lon(numco2lon))
    allocate (co2lat(numco2lat))
    allocate (co2cam2000base(288,192,12))
    allocate (co2cam2000basehalfdeg(720,360,12))
    allocate (co2cam2000basehalfdegtg(720,360,12))
    allocate (emiss_season(numco2lon,numco2lat,12))
    allocate ( emiss_prev(numco2lon,numco2lat,12))
    allocate ( emiss_new(numco2lon,numco2lat,12))
    allocate ( emiss_newCO2(numco2lon,numco2lat,12))
    allocate ( emiss_future(numco2lon,numco2lat,12))
    allocate ( co2_interp(numco2lon,numco2lat))
    allocate ( ship_interp(numco2lon,numco2lat))
    allocate ( emisgcam_fluxfilegrid(numco2lon,numco2lat))



    allocate (tmpco2(numco2lon,numco2lat,12))
    allocate (tmp1(720,360))
    allocate (datearr(numtime))
    allocate (annual_avg(numco2lon,numco2lat))

    call handle_ncerr(nf90_inq_varid(ncidco2, "lon", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidco2,varid,co2lon),subname,__LINE__)
    call handle_ncerr(nf90_inq_varid(ncidco2, "lat", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidco2,varid,co2lat),subname,__LINE__)
    call handle_ncerr(nf90_inq_varid(ncidco2,'date',varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidco2,varid,datearr),subname,__LINE__)

    start3(1)=1
    count3(1)=numco2lon
    start3(2)=1
    count3(2)=numco2lat
    tmp=MAXLOC(datearr,mask=datearr.EQ.19960115)
    start3(3)=tmp(1)
    count3(3)=12
    emiss_season=0.
    scalefactor09x125=area09x125*365.*86400.*1.e-9
    do n=1,10
       tmpco2=1.e36
       start3(3)=tmp(1)+(n-1)*12
       write(6,*)'start3,count3,numtime,tmp(1),yrnum=',start3,count3,numtime,tmp(1),n
       call handle_ncerr(nf90_inq_varid(ncidco2, "CO2_flux", varid),subname,__LINE__)
       call handle_ncerr(nf90_get_var(ncidco2,varid,tmpco2,start3,count3),subname,__LINE__)
       annual_avg = sum(tmpco2,DIM=3)/size(tmpco2,DIM=3)
       do nn=1,12
          where(annual_avg.gt.0.)
             emiss_season(:,:,nn) = emiss_season(:,:,nn)+tmpco2(:,:,nn)/annual_avg
          end where
       end do

    end do
    ! average over 10 years
    emiss_season = emiss_season/10._r8
    write(6,*)'sum emiss_season=',sum(emiss_season)

    emiss_season=0.

    do n=1,10
       tmpco2=1.e36
       start3(3)=tmp(1)+(n-1)*12
       write(6,*)'start3,count3,numtime,tmp(1),yrnum=',start3,count3,numtime,tmp(1),n
       call handle_ncerr(nf90_inq_varid(ncidco2, "CO2_flux", varid),subname,__LINE__)
       call handle_ncerr(nf90_get_var(ncidco2,varid,tmpco2,start3,count3),subname,__LINE__)
       do nn=1,12
          scalefactor09x125=area09x125*dpm(nn)*86400
          tmpco2(:,:,nn)=tmpco2(:,:,nn)*scalefactor09x125
       end do 
       tmpco2=tmpco2/sum(tmpco2)
       emiss_season=emiss_season+tmpco2
    end do 
    emiss_season = emiss_season/sum(emiss_season)
    write(6,*)'sum emiss_season new =',sum(emiss_season)

    ! Read in 2000 co2 base year data from dynamic CO2 flux file.  Convert base units from CO2 to C

    start3(1)=1
    count3(1)=numco2lon
    start3(2)=1
    count3(2)=numco2lat
    tmp=MAXLOC(datearr,mask=datearr.EQ.20000115)
    start3(3)=tmp(1)
    count3(3)=12
    call handle_ncerr(nf90_inq_varid(ncidco2, "CO2_flux", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidco2,varid,co2cam2000base,start3,count3),subname,__LINE__)
    write(6,*)'co2cam2000base start3,count3 =',start3,count3
    do i=1,12
       write(6,*)'co2cam2000base mth max min sum =',i,maxval(co2cam2000base(:,:,i)),minval(co2cam2000base(:,:,i)),sum(co2cam2000base(:,:,i))
    end do
    co2cam2000base=co2cam2000base*12.0/(12.0+(2.*16.0))
    deallocate (datearr)
    deallocate (annual_avg)
    call handle_ncerr(nf90_close(ncidco2),subname,__LINE__)

    !  regrid co2cam2000 base to half degree 
    lonbaseremap(1)=-179.75
    do i=2,720
    lonbaseremap(i)=lonbaseremap(i-1)+.5
    end do
    write(6,*)'lonbaseremap=',lonbaseremap
    latbaseremap(1)=-89.75
    do i=2,360
    latbaseremap(i)=latbaseremap(i-1)+.5
    end do
    write(6,*)'latbaseremap=',latbaseremap
    nlonbaseremap=720
    nlatbaseremap=360

    do i =1,12
       call cremapbin(1  ,360   ,720   ,192    ,288 , &
            co2cam2000base(:,:,i) ,co2cam2000basehalfdeg(:,:,i) ,co2lat    ,co2lon    ,latbaseremap, &
            lonbaseremap  ,numco2lat    ,360   ,1.0d0   , &
            mssng                              )
       scalefactorhalfdeg=area720x360*float(dpm(i))*86400.*1.e-9
       scalefactor09x125=area09x125*float(dpm(i))*86400.*1.e-9
       co2cam2000basehalfdegtg(:,:,i)=co2cam2000basehalfdeg(:,:,i)*scalefactorhalfdeg
       write(6,*)'co2cam2000base global total co2 Tg/mth for month ',i,' = ',sum(co2cam2000base(:,:,i)*scalefactor09x125)
       write(6,*)'co2cam2000basehalfdeg global total co2 Tg/mth for month ',i,' = ',sum(co2cam2000basehalfdegtg(:,:,i))
    end do
    write(6,*)'co2cam2000basehalfdeg global annual total co2 Tg/Yr = ',sum(co2cam2000basehalfdegtg)
    write(6,*)'co2cam2000basehalfdeg global annual total co2 Tg/Yr = ',sum(sum(co2cam2000basehalfdegtg,dim=3))
    scalefactorhalfdeg=area720x360*365.0d0*86400.0d0*1.0e-9
    tmp1(:,:)=sum(co2cam2000basehalfdegtg,dim=3)/scalefactorhalfdeg
    deallocate (tmpco2)
    deallocate (tmp1)


    !  Read in 2000 co2 base year ship data.

    call handle_ncerr(nf90_open(ship2000baseco2,nf90_nowrite,ncidship),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid( ncidship,  'lon', lonshipdimid ),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid( ncidship,  'lat', latshipdimid ),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid( ncidship,  'time', timeshipdimid ),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidship, latshipdimid, len = numshipco2lat),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidship, lonshipdimid, len = numshipco2lon),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidship, timeshipdimid, len = numshipco2time),subname,__LINE__)

    allocate (shiplon(numshipco2lon))
    allocate (shiplat(numshipco2lat))

    call handle_ncerr(nf90_inq_varid(ncidship, "lon", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidship,varid,shiplon),subname,__LINE__)
    call handle_ncerr(nf90_inq_varid(ncidship, "lat", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidship,varid,shiplat),subname,__LINE__)

    allocate (co2ship2000base(numco2lon,numco2lat,12))
    allocate (co2ship2005base(numco2lon,numco2lat,12))
    allocate (co2ship2005basekgms(numco2lon,numco2lat,12))
    allocate (co2ship2000basehalfdeg(numshipco2lon,numshipco2lat,12))
    allocate (co2ship2000basehalfdegkggrdmth(numshipco2lon,numshipco2lat,12))
    allocate (ship_season(numco2lon,numco2lat,12))
    allocate (ship_future(numco2lon,numco2lat,12))

    start3(1)=1
    count3(1)=numshipco2lon
    start3(2)=1
    count3(2)=numshipco2lat
    start3(3)=1
    count3(3)=12

    call handle_ncerr(nf90_inq_varid(ncidship, "emiss_shp", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidship,varid,co2ship2000basehalfdeg,start3,count3),subname,__LINE__)
    call handle_ncerr(nf90_close(ncidship),subname,__LINE__)
    !regrid
    do nn=1,12
       scalefactorhalfdeg=area720x360*dpm(nn)*86400.
       co2ship2000basehalfdegkggrdmth(:,:,nn)=co2ship2000basehalfdeg(:,:,nn)*scalefactorhalfdeg

       call cremapbin(1,  numco2lat,  numco2lon,  numshipco2lat,  numshipco2lon, &
            co2ship2000basehalfdeg(:,:,nn)     ,co2ship2000base(:,:,nn)      ,shiplat    ,shiplon    ,co2lat, &
            co2lon  ,numshipco2lat    ,numco2lat   ,1.0d0   , &
            mssng                              )
       !convert output ship file on fluxfile grid from kg/m2/s to kg/g/month
       scalefactor09x125=area09x125*dpm(nn)*86400.
       co2ship2000base(:,:,nn)=scalefactor09x125*co2ship2000base(:,:,nn)
    end do

    deallocate (shiplat)
    deallocate (shiplon)

    !divide by kg/y = ship scaling
    ship_season(:,:,:)=co2ship2000base/sum(co2ship2000base)
    write(6,*)'global original co2ship2000basehalfdeg in Tg/yr ',sum(co2ship2000basehalfdegkggrdmth)*1.e-9
    write(6,*)'global co2ship2000base in Tg/yr ',sum(co2ship2000base)*1.e-9
    write(6,*)'sum ship_season is ',sum(ship_season)


    ! Handle steps 1 and 2 of downscaling code provided by yuyu
    ! 1) read in base year emissions and create base gcam country emissions
    ! 2) read in base year pop and GDP and calculate GPD per capita 

    !  open 0.5 deg gridded co2 file and read in emission flux in Tg/m2/sec
    !  file created from 1 deg yr 2000 gridcar files http://cdiac.ornl.gov/epubs/ndp/ndp058/ndp058_v2013.html
    !     cat > g360x180.des << EOF
    !     gridtype = lonlat
    !     xsize = 360
    !     ysize = 180
    !     xfirst = -179.5
    !     xinc = 1.
    !     yfirst = 89.5
    !     yinc = -1.
    !     EOF
    !     cdo -f nc input,g360x180.des gridcar_2000.nc < gridcar.2000
    !     cdo remapcon,r720x360 gridcar_2000.nc gridcar_2000.720x360.remapcon.nc
    !     ncks -vgridbox_area gridbox_area_720x360.nc gridcar_2000.720x360.remapcon.nc
    !     ncap2 -A -s "emission_flux=var1/gridbox_area/(86400.*365.)" gridcar_2000.720x360.remapcon.nc
    !     ncatted -a units,emission_flux,c,c,"Tg(species)/m^2/sec" gridcar_2000.720x360.remapcon.nc
    !  open LUT mapping file to map .5 deg gridded data to country.
    !  open Base population files (GCAM and IIASA_UN)
    !  open Base GDP file (GCAM and GDP_IEA_Scaled2GCAM) 
    !  open 1990-2100 GCAM all steps file

    !==========================================================================

    !                                                                                                                                              
    ! Open and read the country2grid mapping file                                                                                                  
    !                                                                                                                                              
    call handle_ncerr(nf90_open(gcam_lut_720x360_mapping,nf90_nowrite,ncidmap),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid(ncidmap, "CountryGridPts", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidmap, dimid, len = numcountrygridpts),subname,__LINE__)

    allocate (Grid05ID(numcountrygridpts))
    allocate (Country05id(numcountrygridpts))
    allocate (areaGrid(numcountrygridpts))
    allocate (areaCountry(numcountrygridpts))
    allocate (areaIntersect(numcountrygridpts))
    allocate (Area_sumGrid(numcountrygridpts))

    call handle_ncerr(nf90_inq_varid(ncidmap, "Grid05ID", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,grid05id),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidmap, "CountryID", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,country05id),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidmap, "areaGrid", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,areaGrid),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidmap, "areaCountry", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,areaCountry),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidmap, "areaIntersect", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,areaIntersect),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidmap, "Area_sumGrid", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidmap,varid,area_sumgrid),subname,__LINE__)

    call handle_ncerr(nf90_close(ncidmap),subname,__LINE__)



    allocate (totalBaseEmissionsCountry(numcountry))
    allocate (countryEmission(numcountry,numsource))
    allocate (countryEmission2(numcountry,numsource))

    countryEmission=0.
    totalArea=0

    !  step1_grid2country                                                                                                                          

    numgas=1
    do ng=1,numgas
       do ns=1,numsource

          countryEmission(:,ns)=0.

          !jt          emisname=source(ns)

          emisname='CO2_flux'

          !**********************************************************
          ! cam2000base from dyn cam file init - convert from monthly tg/grid/month to yearly kg/m2/s
          scalefactorhalfdeg=area720x360*365.*86400.*1.e-9
          co2base2000(:,:,ns)=sum(co2cam2000basehalfdegtg,dim=3)/scalefactorhalfdeg
          !**********************************************************

          ! scale base emissions from kg/m2/s to kg/m2/s *1e9
          !!!fix Check with Yuyu on this  - removed further down for countryEmission
          !!! but the scaled version of co2base2000 is used in step4.  There must be a compensating
          !!! scaling taking place in step 3.
          co2base2000(:,:,ns) = co2base2000(:,:,ns)*1.e9
          !                                                                                                                                              
          ! base emissions in kg/m^2/s  apportion to each country                                                                                        
          !                                                                                                                                              
          do i=1,numcountrygridpts
             latind=nlatbaseremap-(Grid05ID(i)-1)/nlonbaseremap
             lonind=mod(Grid05ID(i)-1,nlonbaseremap)+1
	     ! country Emission kg/country/yr
             countryEmission(country05id(i),ns) = countryEmission(country05id(i),ns)+  co2base2000(lonind,latind,ns)*areaGrid(i)*1000*1000*365*86400*areaIntersect(i)/Area_sumGrid(i)
             if (country05id(i).eq.102) then 
                write(6,*)"GridID=",Grid05ID(i)," emissionsFliped = ",co2base2000(lonind,latind,ns)," GridArea = ", areaGrid(i)," intersect=",areaIntersect(i),"sumgrig=",Area_sumGrid(i)," cemis= ",countryEmission(country05id(i),ns),"latind,lonind=",latind,lonind
             end if

             if(countryEmission(country05id(i),ns)==0) countryEmission(country05id(i),ns)=0.0000000001
             totalBaseEmissionsCountry(country05id(i)) = totalBaseEmissionsCountry(country05id(i)) + co2base2000(lonind,latind,ns)*areaIntersect(i)/Area_sumGrid(i)*86400*365*areaGrid(i)
          end do
          !                                                                                                                                              
          ! Get rid of 1e9 scaling and save country_emissions as kg/country/yr
          !                                                                                                                                              
          countryEmission(:,ns)=countryEmission(:,ns)/1.e9
       end do
    end do

    !                                                                                                                                              
    ! Step2 PopGDP                                                                                                                                 

    !                                                                                                                                              
    ! Open and read the pop gdp downscaling file                                                                                                   
    !                                                                                                                                              

    call handle_ncerr(nf90_open(gcam2emis_downscale,nf90_nowrite,ncidds),subname,__LINE__)

    call handle_ncerr(nf90_inq_dimid(ncidds, "region", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidds, dimid, len = numregion),subname,__LINE__)

    call handle_ncerr(nf90_inq_dimid(ncidds, "country", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidds, dimid, len = numcountry),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid(ncidds, "time", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidds, dimid, len = numtime),subname,__LINE__)

    allocate (regionid(numregion))
    allocate (countryid(numcountry))
    allocate (years(numtime))


    call handle_ncerr(nf90_inq_varid(ncidds, "region", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,regionid),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "country", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,countryid),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "year", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,years),subname,__LINE__)

    allocate (country2region(numcountry))
    allocate (GDPpc_Grow(numcountry))

    allocate (GDP_C_2year(numcountry,numtime))
    allocate (Pop_C(numcountry,numtime))
    allocate (GDP_C_Preliminary(numcountry,numtime))
    allocate (GDP_C_Final(numcountry,numtime))
    allocate (Pop_R_GCAM(numregion,numtime))
    allocate (Pop_R_IIASA(numregion,numtime))
    allocate (GDP_R_GCAM(numregion,numtime))
    allocate (GDP_R_IEA(numregion,numtime))
    allocate (GDP_R_GCAM_pred(numregion,numtime))
    allocate (Diff(numregion,numtime))
    allocate (GDP_Increase(numregion,numtime))

    Pop_C(:,:)=0.0
    GDP_C_Preliminary(:,:)=0.0
    GDP_C_Final(:,:)=0.0
    Pop_R_GCAM(:,:)=0.0
    GDP_R_IEA(:,:)=0.0
    Pop_R_IIASA(:,:)=0.0
    GDP_R_GCAM(:,:)=0.0
    GDP_R_GCAM_pred(:,:)=0.0
    Diff(:,:)=0.0
    GDP_Increase(:,:)=0.0
    country2region(:)=0

    call handle_ncerr(nf90_inq_varid(ncidds, "GDP_GCAM", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,GDP_R_GCAM),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "GDP_IEA", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,GDP_C_2year),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "POP_GCAM", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,Pop_R_GCAM),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "POP_IIASA_UN", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,Pop_C),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidds, "country2region", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidds,varid,country2region),subname,__LINE__)
    call handle_ncerr(nf90_close(ncidds),subname,__LINE__)

    ! Population downscaling                                                                                                                   
    Pop_R_IIASA(:,:)=0
    do i=1,numcountry
       regionind=country2region(countryid(i))
       do j=1,numtime
          Pop_R_IIASA(regionind,j) =  Pop_R_IIASA(regionind,j) + Pop_C(countryid(i),j)
       end do
    end do

    do i=1,numcountry
       regionind=country2region(countryid(i))
       do j=1,numtime
          Pop_C(countryid(i),j) = Pop_C(countryid(i),j) * Pop_R_GCAM(regionind,j)/Pop_R_IIASA(regionind,j)
       end do
    end do

    ! GDP downscaling                                                                                                                          

    idx2100=21
    idx2090=19
    idx2080=17
    idx2000=1
    idx2005=2

    !GCAM regional GDP                                                                                                                         

    !IEA country GDP scaled to GCAM already
    do i=1,numcountry
       regionind=country2region(countryid(i))
       GDP_R_IEA(regionind,:idx2005) =  GDP_R_IEA(regionind,:idx2005) + GDP_C_2year(countryid(i),:idx2005)
    end do

    !preliminary GPD per capita                                                                                                                
    do i=1,numcountry
       regionind=country2region(countryid(i))
       if(regionind.ge.1 .and. regionind.le.14) then

          temp =  ((GDP_R_GCAM(regionind,idx2100)*1000000000.0/Pop_R_GCAM(regionind,idx2100))/(GDP_R_GCAM(regionind,idx2090)*1000000000.0/Pop_R_GCAM(regionind,idx2090) ))**0.1-1      !growth rate from 2100 to 2090
          temp = GDP_R_GCAM(regionind,idx2100)*1000000000.0/Pop_R_GCAM(regionind,idx2100)*((1+temp)**(2150-2100)) !GDP per capita in 2150                  
          GDPpc_Grow(countryid(i)) = (temp/(GDP_C_2year(countryid(i),idx2005) *1000000000.0/Pop_C(countryid(i),idx2005)))**(1.00/(2150.0-2005.0))   ! use 2005 as base year 
          GDP_C_Preliminary(countryid(i),idx2000) = GDP_C_2year(countryid(i),idx2000)*1000000000.0 /Pop_C(countryid(i),idx2000) !2000                                     
          GDP_C_Preliminary(countryid(i),idx2005) = GDP_C_2year(countryid(i),idx2005)*1000000000.0 /Pop_C(countryid(i),idx2005) !2005 base year                           

          GDP_R_GCAM_pred(regionind,idx2000) = GDP_R_GCAM_pred(regionind,idx2000) +  GDP_C_Preliminary(countryid(i),idx2000)*Pop_C(countryid(i),idx2000)/1000000000.0
          GDP_R_GCAM_pred(regionind,idx2005) = GDP_R_GCAM_pred(regionind,idx2005) +  GDP_C_Preliminary(countryid(i),idx2005)*Pop_C(countryid(i),idx2005)/1000000000.0

          do j=3,numtime
	     ! The following line is hardcoded for 5 year time slices. 
             GDP_C_Preliminary(countryid(i),j) = GDP_C_2year(countryid(i),idx2005)*1000000000.0 / Pop_C(countryid(i),idx2005)*(GDPpc_Grow(countryid(i))**((j-3)*5+5)) !from 2005 
             GDP_R_GCAM_pred(regionind,j) = GDP_R_GCAM_pred(regionind,j) +  GDP_C_Preliminary(countryid(i),j)*Pop_C(countryid(i),j)/1000000000.0
             !regional GDP increase between interval                                                                                           
             GDP_Increase(regionind,j) = GDP_Increase(regionind,j) +  (GDP_C_Preliminary(countryid(i),j)*Pop_C(countryid(i),j)/1000000000.0 - GDP_C_Preliminary(countryid(i),j-1)*Pop_C(countryid(i),j-1)/1000000000.0)
          end do
       end if
    end do

    !Regional difference between preliminary GPD and GCAM estimation                                                                              
    do i=1,numregion
       do j=1,numtime
          if(GDP_R_GCAM(i,j).gt.0) then
             Diff(i,j) = GDP_R_GCAM(i,j)-GDP_R_GCAM_pred(i,j)
          end if
       end do
    end do

    !Final GDP in each country                                                                                                                 

    do i=1,numcountry
       regionind = country2region(countryid(i))
       if(regionind.ge.1 .and. regionind.le.14) then
          GDP_C_Final(countryid(i),idx2000) =GDP_C_Preliminary(countryid(i),idx2000)
          GDP_C_Final(countryid(i),idx2005) =GDP_C_Preliminary(countryid(i),idx2005)
          GDP_C_Final(countryid(i),idx2000) = GDP_C_Final(countryid(i),idx2000) * Pop_C(countryid(i),idx2000)/1000000000.0
          GDP_C_Final(countryid(i),idx2005) = GDP_C_Final(countryid(i),idx2005) * Pop_C(countryid(i),idx2005)/1000000000.0

          do j=3,numtime
             GDP_C_Final(countryid(i),j) =  GDP_C_Preliminary(countryid(i),j) +  Diff(regionind,j)*1000000000.0 * &
                  (GDP_C_Preliminary(countryid(i),j)*Pop_C(countryid(i),j)/1000000000.0 - GDP_C_Preliminary(countryid(i),j-1)*Pop_C(countryid(i),j-1)/1000000000.0)/GDP_Increase(regionind,j)/Pop_C(countryid(i),j)
             GDP_C_Final(countryid(i),j) = GDP_C_Final(countryid(i),j) * Pop_C(countryid(i),j)/1000000000.0
          end do
       end if
    end do


    deallocate (years)
    deallocate (regionid)
    deallocate (GDPpc_Grow)
    deallocate (GDP_C_2year)
    deallocate (GDP_C_Preliminary)
    deallocate (GDP_Increase)
    deallocate (GDP_R_GCAM_pred)
    deallocate (GDP_R_GCAM)
    deallocate (GDP_R_IEA)
    deallocate (Pop_C)
    deallocate (POP_R_GCAM)
    deallocate (Pop_R_IIASA)
    deallocate (Diff)
    deallocate (co2ship2000basehalfdegkggrdmth)

  end subroutine gcam2emisfile_init_mod

  !---------------------------------------------------------------------------
  !BOP
  
  ! !IROUTINE: gcam2emisfile_run_mod
  
  ! !INTERFACE:
  subroutine gcam2emisfile_run_mod( EClock, cdata, gcamoemis)

    ! !DESCRIPTION:
    ! Run interface for gcam2emisfile

    ! !USES:
    implicit none

    ! !ARGUMENTS:
    integer, pointer :: EClock(:)
    type(iac_cdata_type) :: cdata
    real(r8), pointer :: gcamoemis(:,:)
    ! !LOCAL VARIABLES:

    integer       ::i,k,allstepyridx,ii,iii,allstep2100idx,allstep2090idx,allstep2080idx,newdate(1)
    integer :: iu,tmp(1),n
    integer :: ymd, tod, dt,mth
    integer :: startyr,endyr,yr,ind,yeararr(12),montharr(12),start3(3),count3(3),datearr(12),ierr,curryr,exp_rate,start(1),count(1),year
    integer, allocatable, dimension(:) :: alldates
    real(r8)  :: startday,endday,calday(1),tot,startday2000
    ! !PARAMETERS:
    integer                            :: date(60),nr,numcat
    real(r8)                           :: time1(60),countrytest,gcamshipnew,gcam_ship_scalar,anntot
    character(len=*),parameter :: subname='(gcam2emisfile_run_mod)'
    !!!fix  set this to number of emissions 
    logical :: emismask(12)
    


    ! !REVISION HISTORY:
    ! Author: J Truesdale, Yuyu Zhou

    !EOP
    !-----------------------------------------------------------------------
    iu  = cdata%i(iac_cdatai_logunit)
    ymd = EClock(iac_EClock_ymd)
    tod = EClock(iac_EClock_tod)
    dt  = EClock(iac_EClock_dt)
    curryr=ymd/10000
    exp_rate=curryr-2000

#ifdef DEBUG
    write(iu,*) trim(subname),' date (same as gcam) = ',ymd,tod
#endif

    ! Use CO2 input from GCAM and downscale data from init to convert from region to country
    
    !
    ! Open and read the GHG_GCAM_allsteps file
    !
    
    call handle_ncerr(nf90_open(gcam2emis_co2_allsteps,nf90_nowrite,ncidallsteps),subname,__LINE__)
    call handle_ncerr(nf90_inq_dimid(ncidallsteps, "region", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidallsteps, dimid, len = numreg),subname,__LINE__)

    call handle_ncerr(nf90_inq_dimid(ncidallsteps, "time", dimid),subname,__LINE__)
    call handle_ncerr(nf90_inquire_dimension(ncidallsteps, dimid, len = ntime),subname,__LINE__)

    allocate (regionid(numreg))
    allocate (years(ntime))

    call handle_ncerr(nf90_inq_varid(ncidallsteps, "region", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidallsteps,varid,regionid),subname,__LINE__)

    call handle_ncerr(nf90_inq_varid(ncidallsteps, "time", varid),subname,__LINE__)
    call handle_ncerr(nf90_get_var(ncidallsteps,varid,years),subname,__LINE__)
    write(6,*)'numreg,numcountry,ntime,nlonbaseremap,nlatbaseremap=',numreg,numcountry,ntime,nlonbaseremap,nlatbaseremap

    allocate (Diff2(numreg))
    allocate (GDP_R(numreg,ntime))
    allocate (GHG_C2000(numcountry,numsource))
    allocate (GHG_C_Final(numcountry))
    allocate (GHG_C_Preliminary(numcountry))
    allocate (GHG_R_GCAM(numreg))
    allocate (GHG_R_GCAM_allsteps(numreg,ntime))
    allocate (GHG_R_Pre(numreg))
    allocate (GHG_Share_R(numreg))
    allocate (GHGpg_Grow(numcountry))
    allocate (emis_future_grid_rcp45(nlonbaseremap,nlatbaseremap))
    allocate (emis_future_grid_rcp45p180(nlonbaseremap,nlatbaseremap))
    allocate (emissionInCountryIn12year(numcountry))
    allocate (share(nlonbaseremap,nlatbaseremap,numcountry))
    allocate (tmp_GHG_R_GCAM_allsteps(numreg,ntime))

    tmp=MAXLOC(years,mask=years.EQ.curryr)
    allstepyridx=tmp(1)
    tmp=MAXLOC(years,mask=years.EQ.2100)
    allstep2100idx=tmp(1)
    tmp=MAXLOC(years,mask=years.EQ.2090)
    allstep2090idx=tmp(1)
    tmp=MAXLOC(years,mask=years.EQ.2080)
    allstep2080idx=tmp(1)

    GDP_R=0.
    GHG_C2000=0

    where(GDP_C_Final==0) GDP_C_Final=10
    do i=1,numcountry
       rid = country2region(countryid(i))
       GDP_R(rid,:ntime) = GDP_R(rid,:ntime)+GDP_C_Final(countryid(i),:ntime)
    end do

    GHG_C2000=countryEmission    !kg/country/yr
    where(GHG_C2000<1.e-10) GHG_C2000=10

    !step 3

    !jt       numsector=10
    numsector=1
    numgas=1
    do ng=1,numgas
       do ns=1,numsector
          ! convert 2000 base country emission to regions                                                                                                                                                       
          GHG_R_GCAM=0.
          GHG_C_Preliminary=0.
          GHG_C_Final=0.
          GHG_R_Pre=0.
          GHG_share_R=0.
          Diff2=0.

	  !!!fix - should just read allsteps once in init not everytime its run

          !jt             emisname='emiss_'//source(ns)
          emisname='CO2_flux'
          call handle_ncerr(nf90_inq_varid(ncidallsteps, emisname,varid),subname,__LINE__)
          call handle_ncerr(nf90_get_var(ncidallsteps,varid,tmp_GHG_R_GCAM_allsteps),subname,__LINE__)
          call handle_ncerr( nf90_close(ncidallsteps),subname,__LINE__)

	  !!!fix GHG_R_GCAM_allsteps this is ordered by region id = netcdf array should be using regionid as coordinate
	  !!!fix GHG_R_GCAM_allsteps region order needs to match GHG_R_GCAM (gcamoemis)
          do i=1,14
	     GHG_R_GCAM_allsteps(i,:)=tmp_GHG_R_GCAM_allsteps(regionid(i),:)
	  end do
          where(GHG_R_GCAM_allsteps==0.) GHG_R_GCAM_allsteps=0.0000000000000000000001
          GHG_R_GCAM_allsteps(:,1)=0.

          ! unpack gcamoemis into GHG_R_GCAM and gcamshipnew and gcamship2005

          gcamshipnew=0.
          numcat=12
          
          do nr=1,numreg
             gcamshipnew=gcamshipnew+gcamoemis(1,(nr-1)*numcat+9)
             write(6,*)'gcamshipnew from region ',nr,' is ',gcamoemis(1,(nr-1)*numcat+9)
          end do
          if (curryr.eq.2005) then
             gcamship2005=gcamshipnew
             write(6,*)'setting gcamship2005',gcamship2005
          endif
          write(6,*)'gcamshipnew total from all regions=',gcamshipnew

          ! elements 2 and 9 of gcamo emission categories correspond to air and shipping co2.  
          ! Use shipping for ship co2 scalar
          !!!fix need to have mapping for these categories in iac_share
          emismask=.false.
          emismask(4)=.true.  !DOM sector
          emismask(5)=.true.  !ENE sector
          emismask(6)=.true.  !IND sector
          emismask(11)=.true. !TRA sector
          
          if (curryr.eq.2000) then
             do i=1,numcountry
                rid = country2region(countryid(i))
                GHG_R_GCAM(rid) = GHG_R_GCAM(rid) + GHG_C2000(countryid(i),ns)/1000000000
             end do
             do nr=1,numreg
                write(6,*)curryr,' GHG_R_GCAM(nr)=',GHG_R_GCAM(nr)
             end do
          else
             write(6,*)'GHG_R_GCAM_allsteps(:,allstepyridx)=',GHG_R_GCAM_allsteps(:,allstepyridx)
             write(6,*)'gcamoemis=',gcamoemis(1,:)
             do nr=1,numreg
                write(6,*)'region ',nr
                GHG_R_GCAM(nr)=sum(gcamoemis(1,(nr-1)*numcat+1:nr*numcat),mask=emismask)
                write(6,*)curryr,'GHG_R_GCAM(nr)=',GHG_R_GCAM(nr)
             end do
          end if
          write(6,*)curryr,'GHG_R_GCAM Total for all Regions=',sum(GHG_R_GCAM(:))

          do i=1,numcountry
             cid=countryid(i)
             rid = country2region(cid)

             !deal with negtive emissions in 2100                                                                                                                                                 
             if(GHG_R_GCAM_allsteps(rid,allstep2100idx)/GHG_R_GCAM_allsteps(rid,allstep2090idx)<0) then
                temp = (((GHG_R_GCAM_allsteps(rid,allstep2090idx)*1000000000.0/GDP_R(rid,allstep2090idx))/(GHG_R_GCAM_allsteps(rid,allstep2080idx)*1000000000.0/GDP_R(rid,allstep2080idx)))**0.1) - 1
                temp = GHG_R_GCAM_allsteps(rid,allstep2090idx)*1000000000.0/GDP_R(rid,allstep2090idx)*((1+temp)**(2150-2090)) !GDP per capita in 2150 from 2090                                                           
             else
                temp = (((GHG_R_GCAM_allsteps(rid,allstep2100idx)*1000000000.0/GDP_R(rid,allstep2100idx))/(GHG_R_GCAM_allsteps(rid,allstep2090idx)*1000000000.0/GDP_R(rid,allstep2090idx)))**0.1) - 1
                temp = GHG_R_GCAM_allsteps(rid,allstep2100idx)*1000000000.0/GDP_R(rid,allstep2100idx)*((1+temp)**(2150-2100)) !GDP per capita in 2150                                                                     
             end if     !growth rate from 2090 to 2100                                                                                                                                            

             GHGpg_Grow(cid) = (temp/(GHG_C2000(cid,ns)/GDP_C_Final(cid,1)))**(1.00/(2150.0-2000.0))                !use 2000 as base year
             GHG_C_Preliminary(cid) = GHG_C2000(cid,ns)/GDP_C_Final(cid,1)*(GHGpg_Grow(cid)**(exp_rate)) !from yearly rate
             GHG_R_Pre(rid) = GHG_R_Pre(rid) +  GHG_C_Preliminary(cid)*GDP_C_Final(cid,allstepyridx) !emissions                                                                                           
             !regional GDP increase between interval                                                                                                                                           
             GHG_share_R(rid) =GHG_share_R(rid) +  GDP_C_Final(cid,allstepyridx)*GHG_C_Preliminary(cid)
             if(.not.(GHG_R_Pre(rid).gt.0 .and. GHG_R_Pre(rid).lt.9999999999999)) write(6,*)'stop cid,ghg_r_pre,rid=',cid,GHG_R_Pre(rid),rid
          end do


          !Regional difference between preliminary Emissions and GCAM estimation                                                                                                                  
          do i=1,numreg
             if(GHG_R_GCAM(i).gt.0) then
                Diff2(i) = GHG_R_GCAM(i)*1000000000-GHG_R_Pre(i)
             end if
          end do

          do i=1,numcountry
             cid=countryid(i)
             rid = country2region(cid)

             if(rid .gt. 0 .and. rid .lt. 15) then
                GHG_C_Final(cid) =GHG_C_Preliminary(cid)*GDP_C_Final(cid,allstepyridx)+  Diff2(rid) * (GHG_C_Preliminary(cid)*GDP_C_Final(cid,allstepyridx) )/GHG_share_R(rid)
                if(curryr.eq.2000 .and. GHG_C_Final(cid).eq.10) GHG_C_Final(cid)=0
                write(6,*)emisname,'GHG_C_Final(',cid,')=',GHG_C_Final(cid)
             end if
          end do

          !step 4 country to grid

          emissionInCountryIn12year=0
          share=0.
          emis_future_grid_rcp45=0.
          globalTotalBefore=0.
          globalTotalAfter=0.
          countryEmission2=0.
          totalBaseEmissionsCountry(:)=0.   
          do i=1,numcountrygridpts
             latind=nlatbaseremap-(Grid05ID(i)-1)/nlonbaseremap
             lonind=mod(Grid05ID(i)-1,nlonbaseremap)+1
             countryEmission2(country05id(i),ns) = countryEmission2(country05id(i),ns)+  co2base2000(lonind,latind,ns)*areaGrid(i)*1000*1000*365*86400*areaIntersect(i)/Area_sumGrid(i)
             if(countryEmission2(country05id(i),ns)==0) countryEmission2(country05id(i),ns)=0.0000000001
             totalBaseEmissionsCountry(country05id(i)) = totalBaseEmissionsCountry(country05id(i)) + co2base2000(lonind,latind,ns)*areaGrid(i)*1000*1000*365*86400*areaIntersect(i)/Area_sumGrid(i)
          end do
          write(6,*)'sum total difference between country emissions is ',sum(countryEmission2-countryEmission)
         
          do i=1,numcountrygridpts
             latind=nlatbaseremap-(Grid05ID(i)-1)/nlonbaseremap
             lonind=mod(Grid05ID(i)-1,nlonbaseremap)+1
             share(lonind,latind,country05id(i))=share(lonind,latind,country05id(i)) +(co2base2000(lonind,latind,ns)*areaGrid(i)*1000*1000*365*86400*areaIntersect(i)/Area_sumGrid(i))/countryEmission2(country05id(i),ns)
          end do

          emisname='emiss_'//source(ns)

          do i=1,numcountry
             cid=countryid(i)
             emissionInCountryIn12year(cid)= GHG_C_Final(cid)
             globalTotalBefore = globalTotalBefore+GHG_C_Final(cid)
          end do
          countrytest=0.
          do latind=nlatbaseremap,1,-1
             do lonind=1,nlonbaseremap
                temp=0
                do k=1,numcountry
                   temp = temp + emissionInCountryIn12year(k)*share(lonind,latind,k)
                   if (k.eq.4) countrytest=countrytest+emissionInCountryIn12year(k)*share(lonind,latind,k)
                end do
                globalTotalAfter = globalTotalAfter+temp
                emis_future_grid_rcp45(lonind,latind)=temp
             end do
          end do

          regiontot=0.
          do i=1,numcountrygridpts
             latind=nlatbaseremap-(Grid05ID(i)-1)/nlonbaseremap
             lonind=mod(Grid05ID(i)-1,nlonbaseremap)+1
             regionind=country2region(country05id(i))
             regiontot(regionind)=regiontot(regionind)+emis_future_grid_rcp45(lonind,latind)/1.e9
          end do

          do i=1,numregion
             write(6,*)curryr,' gcam downscaled CO2 Tg/yr for ',regionname(i),' is ',regiontot(i)
          end do
             
          write(6,*)'before=',globalTotalBefore,' after = ',globalTotalAfter,' chinatest=',countrytest
          write(6,*)'sum of emis_future_grid_rcp45 =',sum(emis_future_grid_rcp45)
       end do
    end do
    emiss_new=0.
    emiss_newCO2=0.
    emiss_future=0.
    gcam_ship_scalar=0.
    co2_interp=0.
    ship_interp=0.
    write(6,*)'emis_future_grid_rcp45 Tg/yr = ',sum(emis_future_grid_rcp45)/1.e9

    !
    !---------------------------------------------------------------------------
    !           ... Check units of emission flux. If necessary, convert to kg/m2/s
    !           ... Original units in kg/grid/yr
    !---------------------------------------------------------------------------

    write(6,*)'emis_rcp_720x360 kg/grid/y global sum=',sum(emis_future_grid_rcp45)
    scalefactorhalfdeg=area720x360*365.*86400.
    emis_future_grid_rcp45=emis_future_grid_rcp45/scalefactorhalfdeg 

    !---------------------------------------------------------------------------
    !       ... Set the transform from 0.5 deg lats to dynco2fluxfile lats
    !       ... swap longitudes of emis_future_grid to match co2flux_file
    !---------------------------------------------------------------------------
    lonp180(:360)=lonbaseremap(361:)
    lonp180(361:)=lonbaseremap(:360)
    where(lonp180.lt.0._r8) lonp180=lonp180+360._r8
    emis_future_grid_rcp45p180(:360,:)=emis_future_grid_rcp45(361:,:)
    emis_future_grid_rcp45p180(361:,:)=emis_future_grid_rcp45(:360,:)

    call cremapbin(1  ,numco2lat   ,numco2lon   ,nlatbaseremap    ,nlonbaseremap , &
         dble(emis_future_grid_rcp45p180)     ,emisgcam_fluxfilegrid  ,latbaseremap    ,lonp180    ,co2lat, &
         co2lon  ,nlatbaseremap    ,numco2lat   ,1.0d0   , &
         mssng                              )

    scalefactor09x125=area09x125*365.*86400.*1.e-9
    write(6,*)'emisgcam_fluxfilegrid after remapping to 0.9x1.25 Tg/Yr ',sum(emisgcam_fluxfilegrid*scalefactor09x125)

    !---------------------------------------------------------------------------
    ! ... Add in annual cycle  convert to kg/grd/yr before multiplying by emiss_season and then back to kg/m2/sec
    !---------------------------------------------------------------------------
    scalefactor09x125=area09x125*365.*86400.
    emisgcam_fluxfilegrid=emisgcam_fluxfilegrid*scalefactor09x125
    write(6,*)'emisgcam_fluxfilegrid after scaling to Kg/Yr ',sum(emisgcam_fluxfilegrid)
    do mth=1,12
       write(6,*)'sum emiss_season for mth ',mth,' is ',sum(emiss_season(:,:,mth))
       emiss_new(:,:,mth) = sum(emisgcam_fluxfilegrid(:,:)) * emiss_season(:,:,mth)
       emiss_newCO2(:,:,mth) = sum(emisgcam_fluxfilegrid(:,:)) * emiss_season(:,:,mth)*((12.0+2.*16.0)/12.)
       write(6,*)'emisgcam_fluxfilegrid for mth ',mth,' in Kg/mth is ',sum(emiss_new(:,:,mth))
       write(6,*)'CO2 for mth ',mth,' in Kg/mth is ',sum(emiss_newCO2(:,:,mth))
    end do
    write(6,*)'emisgcam_fluxfilegrid for the year in Kg/yr is ',sum(emiss_new)
    write(6,*)'CO2 for the year in Kg/yr is ',sum(emiss_newCO2)
    !---------------------------------------------------------------------------
    ! ... back to kg/m2/sec
    !---------------------------------------------------------------------------

    do mth=1,12
       scalefactor09x125=area09x125*dpm(mth)*86400._r8
       emiss_new(:,:,mth) = emiss_new(:,:,mth)/scalefactor09x125
       emiss_newCO2(:,:,mth) = emiss_newCO2(:,:,mth)/scalefactor09x125
    end do

    !---------------------------------------------------------------------------
    ! Update the CO2fluxfile with downscaled gcam emission data
    !---------------------------------------------------------------------------

    if (ymd/10000.ge.2010) then
       
       !---------------------------------------------------------------------------
       !... gcam only has shipping emissions from 2005 onward, need to extrapolate
       !... 2000 base ship data to 2005 base year
       !---------------------------------------------------------------------------
       write(6,*)'using gcamship2005 ',gcamship2005,ymd
       gcam_ship_scalar=gcamshipnew/gcamship2005

       if (curryr.eq.2010) then
          write(6,*)'setting co2ship2005base sum=',sum(co2ship2000base),ymd
          co2ship2005base=co2ship2000base*gcam_ship_scalar
       end if

       !  Scale 2005 co2 observed ship emissions by change in gcam ship emissions since 2005
       ship_future=sum(co2ship2005base)*gcam_ship_scalar*ship_season

       write(6,*)'co2ship2005base Tg/yr=',sum(co2ship2005base)*1.e-9,' gcam_ship_scalar=',gcam_ship_scalar,' curryr=',curryr
       write(6,*)'co2ship2005base*gcam_ship_scalar Tg/yr=',sum(co2ship2005base)*gcam_ship_scalar*1.e-9
       write(6,*)'co2ship2005base*gcam_ship_scalar blown out on grid Tg/yr=',sum(ship_future)*1.e-9
!  Convert Kg/g/mth to kg/m2/s for interpolation and adding to 
       do mth=1,12
          scalefactor09x125=area09x125*dpm(mth)*86400.
          ship_future(:,:,mth)=ship_future(:,:,mth)/scalefactor09x125
	  co2ship2005basekgms(:,:,mth)=co2ship2005base(:,:,mth)/scalefactor09x125
       end do

       !---------------------------------------------------------------------------
       !           ... make sure curr gcam date is 5 years advanced 
       !           ... from last dynco2fluxfile date
       !---------------------------------------------------------------------------
       call handle_ncerr(nf90_open(dynco2fluxfile,nf90_write,ncidco2),subname,__LINE__)
       call handle_ncerr( nf90_inq_dimid( ncidco2,  'lon', londimid ),subname,__LINE__)
       call handle_ncerr( nf90_inq_dimid( ncidco2,  'lat', latdimid ),subname,__LINE__)
       call handle_ncerr( nf90_inq_dimid( ncidco2,  'time', timedimid ),subname,__LINE__)
       call handle_ncerr(nf90_inquire_dimension(ncidco2, latdimid, len = numco2lat),subname,__LINE__)
       call handle_ncerr(nf90_inquire_dimension(ncidco2, londimid, len = numco2lon),subname,__LINE__)
       call handle_ncerr(nf90_inquire_dimension(ncidco2, timedimid, len = numtime),subname,__LINE__)
       allocate (alldates(numtime))
       call handle_ncerr(nf90_inq_varid(ncidco2,'date',varid),subname,__LINE__)
       call handle_ncerr(nf90_get_var(ncidco2,varid,alldates),subname,__LINE__)
       year=ymd/10000
       year=(year-5)*10000+115
       tmp=MAXLOC(alldates,mask=alldates.EQ.year)
       start=tmp(1)
       count=12
       write(6,*)'alldates=',alldates
       write(6,*)'start,count,year=',start,count,ymd,year
       call handle_ncerr(nf90_get_var(ncidco2,varid,datearr,start,count),subname,__LINE__)
       deallocate (alldates)
       yeararr=datearr/10000
       montharr=(datearr-yeararr*10000)/100

       !  Grab off 12 months for interpolation .. assume 5 years from current date.
       
       start3(1)=1
       count3(1)=numco2lon
       start3(2)=1
       count3(2)=numco2lat
       start3(3)=tmp(1)
       count3(3)=12
       call handle_ncerr(nf90_inq_varid(ncidco2, "CO2_flux", varid),subname,__LINE__)
       call handle_ncerr(nf90_get_var(ncidco2,varid,emiss_prev,start3,count3),subname,__LINE__)      
       scalefactor09x125=area09x125*365.*86400.*1.e-9
       anntot=0.
       do mth=1,12
          write(6,*)'CO2 (wo/ship) for mth ',mth,' in Tg/yr is ',sum(emiss_prev(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
	  anntot=anntot+sum(emiss_prev(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
       end do
       write(6,*)'CO2 annual total (wo/ship) in Tg/yr is ',anntot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Use this for raw downscaled GCAM emissions
!
       emiss_future=emiss_newCO2
!

!
!  First year add ship base to emis_prev  - this will be used for interpolation only not written back to dataset
!
       if (curryr.eq.2010) then
	  emiss_prev=emiss_prev+co2ship2005basekgms
       end if
!
!  Add ship future to emis_future for interpolation
!
       emiss_future=emiss_future+ship_future


       anntot=0.
       do mth=1,12
          write(6,*)'CO2 (w/ship) for mth ',mth,' in Tg/yr is ',sum(emiss_prev(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
	  anntot=anntot+sum(emiss_prev(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
       end do
       write(6,*)'CO2 prev annual (w/ship) in Tg/yr is ',anntot

       anntot=0.
       do mth=1,12
          write(6,*)'CO2 (w/ship) for mth ',mth,' in Tg/yr is ',sum(emiss_future(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
	  anntot=anntot+sum(emiss_future(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
       end do
       write(6,*)'CO2 future annual (w/ship) in Tg/yr is ',anntot


       startyr=yeararr(1)+1
       endyr=ymd/10000
       ind=start3(3)+11
       do yr = startyr, endyr
          anntot=0.
          do mth=1,12
             startday=yeararr(1)*365+doy(mth)
             startday2000=2000*365+doy(mth)
             endday=ymd/10000*365+doy(mth)
             calday(1) = yr*365.0  + doy(mth)
             date1(1) = yr*10000 + (mth)*100 + dom(mth)
             scalefactor09x125=area09x125*365.*86400.*1.e-9
             write(6,*)'interpolating co2 for yr,mth,startyr,endyr,startday,calday,endday=',yr,mth,startyr,endyr,startday,calday,endday
             write(6,*)'CO2 prev (w/ship) for mth ',mth,' in Tg/yr is ',sum(emiss_prev(:,:,mth)*scalefactor09x125*dpm(mth)/365.)
             write(6,*)'CO2 fut (w/ship) for mth ',mth,' in Tg/yr is ',sum(emiss_future(:,:,mth)*scalefactor09x125*dpm(mth)/365.)

             ! don't have to interpolate for endyr thats given.
	     if (yr.eq.endyr) then
                co2_interp=emiss_future(:,:,mth)
             else
                ! interpolate co2
                call interp2d( startday, endday, calday(1), &
                     emiss_prev(:,:,mth), &
                     emiss_future(:,:,mth), &
                     co2_interp(:,:) )
             end if
             write(6,*)'Interpolated CO2 for mth ',mth,' in Tg/mth is ',sum(co2_interp(:,:)*scalefactor09x125*dpm(mth)/365.)
             anntot=anntot+sum(co2_interp(:,:)*scalefactor09x125*dpm(mth)/365.)
             ind=ind+1
             start3(1)=1
             count3(1)=numco2lon
             start3(2)=1
             count3(2)=numco2lat
             start3(3)=ind
             count3(3)=1
             start=ind
             count=1
             write(6,*)'writing interpolated co2 month',mth,start3,count3
             call handle_ncerr(nf90_inq_varid(ncidco2, "CO2_flux", varid),subname,__LINE__)
             call handle_ncerr( nf90_put_var(ncidco2, varid,co2_interp,start3,count3),subname,__LINE__)
             call handle_ncerr(nf90_inq_varid(ncidco2, "date", varid),subname,__LINE__)
             call handle_ncerr( nf90_put_var(ncidco2, varid,date1,start,count),subname,__LINE__)
             call handle_ncerr(nf90_inq_varid(ncidco2, "datesec", varid),subname,__LINE__)
             call handle_ncerr( nf90_put_var(ncidco2, varid,datesec1,start,count),subname,__LINE__)
             call handle_ncerr(nf90_inq_varid(ncidco2, "time", varid),subname,__LINE__)
             call handle_ncerr( nf90_put_var(ncidco2, varid,calday,start,count),subname,__LINE__)
          end do
	  write(6,*)'annual total for year ',yr,' is ',anntot
       end do
       call handle_ncerr( nf90_close(ncidco2),subname,__LINE__)
    end if


    deallocate (regionid)
    deallocate (years)
    deallocate(Diff2)
    deallocate(GDP_R)
    deallocate(GHG_C2000)
    deallocate(GHG_C_Final)
    deallocate(GHG_C_Preliminary)
    deallocate(GHG_R_GCAM)
    deallocate(GHG_R_GCAM_allsteps)
    deallocate(GHG_R_Pre)
    deallocate(GHG_Share_R)
    deallocate(GHGpg_Grow)
    deallocate(emis_future_grid_rcp45)
    deallocate(emis_future_grid_rcp45p180)
    deallocate(emissionInCountryIn12year)
    deallocate(share)
    deallocate(tmp_GHG_R_GCAM_allsteps)

  end subroutine gcam2emisfile_run_mod


  !---------------------------------------------------------------------------
  !BOP

  ! !IROUTINE: gcam2emisfile_final_mod

  ! !INTERFACE:
  subroutine gcam2emisfile_final_mod( )

    ! !DESCRIPTION:
    ! Finalize gcam2emisfile model
    ! !USES:
    implicit none

    ! !ARGUMENTS:

    ! !LOCAL VARIABLES:
    integer :: iu
    character(len=*),parameter :: subname='(gcam2emisfile_final_mod)'

    ! !REVISION HISTORY:
    ! Author: J Truesdale, Yuyu Zhou

    !EOP

    !---------------------------------------------------------------------------

    !    iu  = nint(cdata(iac_cdata_logunit))
    !    write(iu,*) trim(subname)

    deallocate (Area_sumGrid)
    deallocate (Country05id)
    deallocate (Grid05ID)
    deallocate (area09x125)
    deallocate (area720x360)
    deallocate (areaCountry)
    deallocate (areaGrid)
    deallocate (areaIntersect)
    deallocate (co2_interp)
    deallocate (co2base2000)
    deallocate (co2cam2000base)
    deallocate (co2cam2000basehalfdeg)
    deallocate (co2cam2000basehalfdegtg)
    deallocate (co2lat)
    deallocate (co2lon)
    deallocate (co2ship2000base)
    deallocate (co2ship2005base)
    deallocate (co2ship2000basehalfdeg)
    deallocate (co2ship2005basekgms)
    deallocate (countryEmission)
    deallocate (countryEmission2)
    deallocate (countryid)
    deallocate (country2region)
    deallocate (emisgcam_fluxfilegrid)
    deallocate (emiss_future)
    deallocate (emiss_new)
    deallocate (emiss_newCO2)
    deallocate (emiss_prev)
    deallocate (ship_future)
    deallocate (ship_interp)
    deallocate (ship_season)
    deallocate (emiss_season)
    deallocate (totalBaseEmissionsCountry)
    deallocate (GDP_C_Final)
  end subroutine gcam2emisfile_final_mod


  subroutine handle_ncerr( ret, mes, line )

    !----------------------------------------------------------------------- 
    ! Purpose: 
    ! Check netCDF library function return code.  If error detected 
    ! issue error message then abort.
    !
    ! Author: B. Eaton
    !----------------------------------------------------------------------- 

    !-----------------------------------------------------------------------
    use netcdf
    !-----------------------------------------------------------------------

    integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
    character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
    integer, intent(in), optional :: line
    !-----------------------------------------------------------------------

    if ( ret .ne. NF90_NOERR ) then
       if(present(line)) then
          write(6,*) mes, line
       else   
          write(6,*) mes
       end if
       write(6,*) nf90_strerror( ret )
       !jt        call endrun ('HANDLE_NCERR')
       call abort
    endif

    return

  end subroutine handle_ncerr

  !====================================================================================

  subroutine interp2d( t1, t2, tint, f1, f2, fint )
    !-----------------------------------------------------------------------
    ! 	... Linearly interpolate between f1(t1) and f2(t2) to fint(tint).
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    ! 	... Dummy arguments
    !-----------------------------------------------------------------------
    real(r8), intent(in) :: &
         t1, &            ! time level of f1
         t2, &            ! time level of f2
         tint             ! interpolant time
    real(r8), dimension(:,:), intent(in) :: &
         f1, &            ! field at time t1
         f2               ! field at time t2

    real(r8), intent(out) :: &
         fint(:,:) ! field at time tint

    !-----------------------------------------------------------------------
    ! 	... Local variables
    !-----------------------------------------------------------------------
    integer  :: j, plat
    real(r8) :: factor

    plat = size(f1,2)
    factor = (tint - t1)/(t2 - t1)
    do j = 1,plat
       fint(:,j) = f1(:,j) + (f2(:,j) - f1(:,j))*factor
    end do

  end subroutine interp2d
function nearest_index (value, array)
    !=======================================================================
    !
    !     nearest_index = index of nearest data point within "array" corresponding to
    !            "value".
    !
    !     inputs:
    !
    !     value  = arbitrary data...same units as elements in "array"
    !     array  = array of data points  (must be monotonically increasing)
    !     ia     = dimension of "array"
    !
    !     output:
    !
    !     nearest_index =  index of nearest data point to "value"
    !             if "value" is outside the domain of "array" then nearest_index = 1
    !             or "ia" depending on whether array(1) or array(ia) is
    !             closest to "value"
    !
    !             note: if "array" is dimensioned array(0:ia) in the calling
    !                   program, then the returned index should be reduced
    !                   by one to account for the zero base.
    !
    !     example:
    !
    !     let model depths be defined by the following:
    !     parameter (km=5)
    !     dimension z(km)
    !     data z /5.0, 10.0, 50.0, 100.0, 250.0/
    !
    !     k1 = nearest_index (12.5, z, km)
    !     k2 = nearest_index (0.0, z, km)
    !
    !     k1 would be set to 2, and k2 would be set to 1 so that
    !     z(k1) would be the nearest data point to 12.5 and z(k2) would
    !     be the nearest data point to 0.0
    !
    !=======================================================================

    integer :: nearest_index, ia, i, ii
    real(r8) :: value
    real(r8), dimension(:) :: array
    logical keep_going
    character(len=*),parameter :: subname='(nearest_index)'

    ia = size(array)

    do i=2,ia
       if (array(i) < array(i-1)) then
          write (6,*) '=> Error: "nearest_index" array must be monotonically increasing when searching for nearest value to ',value
          write (6,*) '          array(i) < array(i-1) for i=',i 
          write (6,*) '          array(i) for i=1..ia follows:'
          do ii=1,ia
             write (6,*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          call shr_sys_abort(subname//' ERROR nearest_index array must be monotonically increasing.')
       endif
    enddo
    if (value < array(1) .or. value > array(ia)) then
       if (value < array(1))  nearest_index = 1
       if (value > array(ia)) nearest_index = ia
    else
       i=1
       keep_going = .true.
       do while (i <= ia .and. keep_going)
          i = i+1
          if (value <= array(i)) then
             nearest_index = i
             if (array(i)-value > value-array(i-1)) nearest_index = i-1
             keep_going = .false.
          endif
       enddo
    endif
  end function nearest_index

! NCLFORTSTART
      subroutine cremapbin(plev   ,plato   ,plono   ,plat    ,plon , &
                           xx     ,yy      ,clat    ,clon    ,clato, &
                           clono  ,nlat    ,nlato   ,bin_factor    , &
                           xxmsg                                   )
!
!--------1---------2---------3---------4---------5---------6---------7--
!
! Grid-Box Binning
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plev      ! vertical dimension of input/output field
      integer plato     ! latitude dimension of output field
      integer plono     ! longitude dimension of output field
      integer plat      ! latitude dimension of input field
      integer plon      ! longitude dimension of input field
!
      double precision xx(plon ,plat ,plev) ! input analysis field
      double precision yy(plono,plato,plev) ! horizontally interpolated
!                                           ! (output) field
      double precision clat (plat )         ! Input latitude in degrees
!                                           ! oriented S->N
      double precision clon (plon )         ! Input longitude in degrees
!                                           ! oriented W->E
      double precision clato(plato)         ! Output latitude in degrees
!                                           ! oriented S->N
      double precision clono(plono)         ! Output longitude in degrees
!                                           ! oriented W->E
      integer nlat                          ! Number of Global Gaussian latitudes (input)
      integer nlato                         ! Number of Global Gaussian latitudes (output)
      double precision bin_factor           ! bin-box area expansion/contraction factor
!                                           ! relative to output grid-box area.
      double precision xxmsg
!
!-----------------------------------------------------------------------
!
! NCLEND
!
!---------------------------Local workspace-----------------------------
!
!                                               ! Max # of box segments
      integer max_segs
      parameter (max_segs = 100000 )
!                                               
      integer i, j, ii, jj, k, jfirst, jfirsto ! Indices
      integer nx, ny, nx_max, ny_max          

      integer plon2, plonhalf
      integer i_in(max_segs),i_out(max_segs)
      integer j_in(max_segs),j_out(max_segs)
      integer grid_flag, grido_flag ! grid flags: 0=Regular, 1=Gaussian
      double precision xx_loc(plon*2, plat, plev)
      double precision pi, pio180, pio2, factor
      double precision flat  (plat)   ,flon  (plon*2  )
      double precision flato (plato  ),flono (plono   )
      double precision flati (plat+1), floni (plon*2+1)
      double precision flatoi(plato+1),flonoi(plono+1  )
      double precision tmps, tmpn, tmp(plono,plato)
      double precision edge_w (plon*2), edge_e (plon*2), edge_s (plat ),edge_n (plat )
      double precision edgeo_w(plono ), edgeo_e(plono ), edgeo_s(plato),edgeo_n(plato)
      double precision sin_s (plat ),sin_n (plat )
      double precision sino_s(plato),sino_n(plato)
      double precision dx(max_segs), dy(max_segs)
      double precision distmin, dist, zero, three
      double precision dlat, dlato, eps
!
! - the following are only relevant for grids that are Gaussian
!
      double precision, allocatable :: flat_glob (:)  ! Global Gaussian latitudes (based on input  grid resolution)
!                                                     ! (radians)
      double precision, allocatable :: flato_glob(:)  ! Global Gaussian latitudes (based on output grid resolution)
!                                                     ! (radians)
      double precision, allocatable :: gw_glob   (:)  ! Global Gaussian weights   (based on input  grid resolution)
      double precision, allocatable :: gwo_glob  (:)  ! Global Gaussian weights   (based on output grid resolution)
      integer ierror
!
!-----------------------------------------------------------------------
!
      zero     = 0.d0
      three    = 3.d0
      plon2    = plon*2
      plonhalf = plon/2
      pi       = 4.d0*atan(1.d0)
      pio180   = pi/180.d0
      pio2     = pi/2.d0
      eps      = 1.d-5
!
! Sanity checks
!
! djs if(bin_factor .lt. 0.05d0) then
      if(bin_factor .lt. 1.00d0) then
         write(6,*) 'ERROR ("CREMAPBIN"): binning factor out of range'
         write(6,*) 'bin_factor = ', bin_factor
         call abort
      end if
      if(clat(3) - clat(2) .lt. zero) then
         write(6,*) 'ERROR ("CREMAPBIN"): Input latitudes oriented'
         write(6,*) '                     N->S. Should be S->N'
         call abort
      end if
      if(clato(3) - clato(2) .lt. zero) then
         write(6,*) 'ERROR ("CREMAPBIN"):  Output latitudes oriented'
         write(6,*) '                      N->S. Should be S->N'
         call abort
      end if
      if(clon(3) - clon(2) .lt. zero) then
         write(6,*) 'ERROR ("CREMAPBIN"): Input longitudes oriented'
         write(6,*) '                     E->W. Should be W->E'
         call abort
      end if
      if(clono(3) - clono(2) .lt. zero) then
         write(6,*) 'ERROR ("CREMAPBIN"): Output longitudes oriented'
         write(6,*) '                     E->W. Should be W->E'
         call abort
      end if
!
! Determine if input/output grids are Regular or Gaussian
!
      dlat  = ( clat (plat ) - clat (1) ) /(plat -1)
      dlato = ( clato(plato) - clato(1) ) /(plato-1)
      grid_flag  = 0
      grido_flag = 0
      do j = 1,plat-1
         if( abs (clat (j+1) - clat (j) - dlat ) .gt. eps) grid_flag = 1
      end do
      do j = 1,plato-1
         if( abs (clato(j+1) - clato(j) - dlato) .gt. eps) grido_flag =1
      end do
!
! Get global lats/weights for those grids that are Gaussian
!
      allocate ( flat_glob(nlat) )
      allocate ( gw_glob  (nlat) )
      if(grid_flag .eq. 1) then
         if(nlat .lt. plat) then
            write(6,*) 'ERROR ("CREMAPBIN"): number of latitudes for '
            write(6,*) 'the input grid cannot be greater than the  '
            write(6,*) 'global number of latitudes for that grid   '
            write(6,*) 'resolution'
            write(6,*) 'nlat, plat = ', nlat, plat
            call abort
         end if
         call binning_get_global_lats_wgts(nlat, flat_glob, gw_glob)
      end if
!
      allocate ( flato_glob(nlato) )
      allocate ( gwo_glob  (nlato) )
      if(grido_flag .eq. 1) then
         if(nlato .lt. plato) then
            write(6,*) 'ERROR ("CREMAPBIN"): number of latitudes for '
            write(6,*) 'the output grid cannot be greater than the '
            write(6,*) 'global number of latitudes for that grid   '
            write(6,*) 'resolution'
            write(6,*) 'nlato, plato = ', nlato, plato
            call abort
         end if
         call binning_get_global_lats_wgts(nlato, flato_glob, gwo_glob)
      end if
!
! Copy input data to wrap-around array (wrap half-way around globe 
! at each end of x-direction)
!
      do k = 1,plev
         do j = 1,plat
            ii = plonhalf
            do i = 1,plon2
               ii = ii + 1
               if(ii .gt. plon) ii = 1
               xx_loc(i,j,k) = xx(ii,j,k)
            end do
         end do
      end do
!
! Convert input/output grid coordinates to radians (wrap half-way around
! globe at each end of x-direction of input grid)
!
      ii   = plonhalf
      do i = 1,plon2
         ii = ii + 1
         if(ii .gt. plon    )      ii = 1
         if(i  .le. plonhalf)      flon(i) = clon(ii)*pio180-4*pio2
         if(i  .gt. plonhalf+plon) flon(i) = clon(ii)*pio180+4*pio2
         if(i  .gt. plonhalf .and. i .le. plonhalf+plon)  flon(i) = clon(ii)*pio180
      end do

      do j = 1,plat
         flat (j) = clat (j)*pio180
      end do

      do i = 1,plono
         flono(i) = clono(i)*pio180
      end do

      do j = 1,plato
         flato(j) = clato(j)*pio180
      end do
!
! Map "regional" latitudes into global latitude arrays for input/output grids
!
      if(grid_flag .eq. 1) then
         call binning_map_lats(nlat , plat , flat , flat_glob , jfirst )
      end if
      if(grido_flag .eq. 1) then
         call binning_map_lats(nlato, plato, flato, flato_glob, jfirsto)
      end if
!
! Compute box edges for input and output grids
!
      call binning_map_edges(plat      , plon2 , nlat      , jfirst  , &
                             flon      , flat  , gw_glob   ,           &
                             grid_flag , floni , flati     )
      call binning_map_edges(plato     , plono , nlato     , jfirsto ,  &
                             flono     , flato , gwo_glob  ,            &
                             grido_flag, flonoi, flatoi    )
!
! Copy grid interfaces to "edge" arrays
!
      do i = 1,plon*2
        edge_w(i) = floni(i  )
        edge_e(i) = floni(i+1)
      end do

      do j = 1,plat
        edge_s(j) = flati(j  )
        edge_n(j) = flati(j+1)
        sin_s (j) = sin(edge_s(j))
        sin_n (j) = sin(edge_n(j))
      end do
!
! Expand/contract bin box area for each output grid box by "bin_factor"
!
      factor = sqrt(bin_factor)

      do i = 1,plono
        edgeo_w(i) = flono(i) - ( flono (i  ) - flonoi(i) )*factor
        edgeo_e(i) = flono(i) + ( flonoi(i+1) - flono (i) )*factor
      end do

      do j = 1,plato
        tmps       = flato(j) - ( flato (j  ) - flatoi(j) )*factor
        tmpn       = flato(j) + ( flatoi(j+1) - flato (j) )*factor
        edgeo_s(j) = max( tmps, -pio2) - max( ( tmpn - pio2), zero)
        edgeo_n(j) = min( tmpn,  pio2) + max( (-pio2 - tmps), zero)
        sino_s (j) = sin(edgeo_s(j))
        sino_n (j) = sin(edgeo_n(j))
      end do
!
! Make vector of box segments in x-direction
!
      nx = 0
      do i = 1,plono
         do ii = 1,plon*2
            if(edge_e (ii) .gt. edgeo_w( i) .and. edgeo_e( i) .gt. edge_w (ii) ) then
               nx = nx + 1
               if(nx .gt. max_segs) then
                  write(6,*) 'ERROR  ("CREMAPBIN"):  number of box'
                  write(6,*) 'segments greater than "max_segs"'
                  call abort
               end if
               i_in (nx) = ii
               i_out(nx) = i
               dx   (nx) = min(min(min(edge_e(ii)-edge_w(ii),   &
                                       edgeo_e(i)-edgeo_w(i) ), &
                                       edge_e(ii)-edgeo_w(i) ), &
                                       edgeo_e(i)-edge_w(ii) )
            end if
            if(edge_w (ii) .ge. edgeo_e( i)) exit
         end do
      end do
!
! Make vector of box segments in y-direction
!
      ny = 0
      do j = 1,plato
         do jj = 1,plat
            if(edge_n (jj) .gt. edgeo_s( j) .and. edgeo_n( j) .gt. edge_s (jj) ) then
               ny = ny + 1
               if(ny .gt. max_segs) then
                  write(6,*) 'ERROR  ("CREMAPBIN"):  number of box'
                  write(6,*) 'segments greater than "max_segs"'
                  call abort
               end if
               j_in (ny) = jj
               j_out(ny) = j
               distmin   = edge_n(jj)-edge_s(jj)
               dy(ny)    = sin_n (jj)-sin_s (jj)
               dist      = edgeo_n(j)-edgeo_s(j)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sino_n(j)-sino_s(j)
               end if
               dist      = edge_n(jj)-edgeo_s(j)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sin_n(jj)-sino_s(j)
               end if
               dist      = edgeo_n(j)-edge_s(jj)
               if(dist .lt. distmin) then
                  distmin = dist
                  dy(ny)  = sino_n(j)-sin_s(jj)
               end if
            end if
            if(edge_s (jj) .ge. edgeo_n( j)) exit
         end do
      end do

      nx_max = nx
      ny_max = ny
!
! Begin weighted binning
!
      do k = 1,plev
         do j = 1,plato
            do i = 1,plono
               yy(i,j,k) = 0.
            end do
         end do
      end do

      do k = 1,plev
         do ny = 1,ny_max
            j  = j_out(ny)
            jj = j_in (ny)
            do nx = 1,nx_max
               i  = i_out(nx)
               ii = i_in (nx)
               yy(i,j,k) = yy(i,j,k) + xx_loc(ii,jj,k)*dx(nx)*dy(ny)
            end do
         end do
      end do
!
! Normalize
!
      do j = 1,plato
         do i = 1,plono
            tmp(i,j) = (edgeo_e(i) - edgeo_w(i))*(sino_n(j) - sino_s(j))
         end do
      end do

      do k = 1,plev
         do j = 1,plato
            do i = 1,plono
               yy(i,j,k) = yy(i,j,k)/tmp(i,j)
            end do
         end do
      end do
!
      deallocate ( flat_glob  )
      deallocate ( gw_glob    )
      deallocate ( flato_glob )
      deallocate ( gwo_glob   )
!
! CRUDE .... 
! .   At any level where the input "xx" has a missing value
! .   set the corresponding "yy" level to missing.
!
      do k = 1,plev
         do j = 1,plat
            do i = 1,plon 
               if (xx(i,j,k).eq.xxmsg) then
                   do jj = 1,plato
                      do ii = 1,plono
                         yy(ii,jj,k) = xxmsg
                      end do
                   end do
               end if
               go to 100
            end do
         end do
  100    continue
      end do
    
      return
    end subroutine cremapbin
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine binning_get_global_lats_wgts(nlat, lat_glob, gw_glob)
!
!--------1---------2---------3---------4---------5---------6---------7--
!
! Compute Global Gaussian latitudes/weights based upon # of latitudes
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer nlat                      ! Number of Global Gaussian latitudes
      double precision lat_glob (nlat)  ! Global Gaussian latitudes (radians)
      double precision gw_glob  (nlat)  ! Global Gaussian weights
!
!---------------------------Local workspace-----------------------------
!
      integer ierror, lwork
      double precision pio2
      double precision, allocatable :: work(:) ! Work array
!
!-----------------------------------------------------------------------
!
      pio2 = 2.d0*atan(1.d0)
!
      if(nlat .le. 2) then
        write(6,*) 'Error in "cremapbin": Not enough Gaussian latitudes'
        write(6,*) 'nlat = ', nlat
        call abort
      end if
!
      lwork = 4*nlat*(nlat+1)+2
      allocate ( work(lwork) )
      call gaqdncl(nlat,lat_glob,gw_glob,work,lwork,ierror)
      deallocate ( work )

      if(ierror .ne. 0) then
         write(6,*)
         write(6,*) 'Error: in call to routine "gaqdncl", ierror = ',ierror
         if(ierror .eq. 1) then
            write(6,*) "Not enough work space declared for number of"
            write(6,*) "Gaussian latitudes" 
            write(6,*) 'lwork, nlat     = ', lwork,nlat
            write(6,*) 'lwork should be = ', 4*nlat*(nlat+1)+2
         end if
         call abort
      end if
!
      lat_glob(:) = lat_glob(:) - pio2
!
      return
    end subroutine binning_get_global_lats_wgts
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine binning_map_lats(nlat, plat, flat, flat_glob, jfirst)
!
!--------1---------2---------3---------4---------5---------6---------7--
!
! Map "regional" latitudes into global latitude arrays for input/output grids
! and check that the grid latitudes are an identical subset of the global array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer nlat                      ! Number of Global Gaussian latitudes
      integer plat                      ! Number of grid   Gaussian latitudes
      integer jfirst                    ! index of Global lat array that maps
!                                       ! into the first grid lat array
      double precision flat_glob(nlat)  ! Global Gaussian latitudes (radians)
      double precision flat     (plat)  ! grid Gaussian latitudes (radians)
!
!---------------------------Local workspace-----------------------------
!
      integer j, jj
      double precision eps
      logical found
!
!-----------------------------------------------------------------------
!
      eps = 1.d-5
!
! Find latitude in Global array that corresponds to the first latitude
! of the grid array.
!
      found  = .false.
      jfirst = 0
      do j = 1,nlat
         if( abs(flat_glob(j) - flat(1)) .lt. eps ) then
            found  = .true.
            jfirst = j
            exit
         end if
      end do
!
      if(.not. found) then
         write(6,*) 'Error in "cremapbin":'
         write(6,*) "Could not map global lat array into grid array"
         call abort
      end if
!
      if(plat+jfirst-1 .gt. nlat) then
         write(6,*) 'Error in "cremapbin":'
         write(6,*) "Stepping out of bounds of the global lat array"
         call abort
      end if
!
! Test that subsequent grid lats all match the global lat array
!
      do j = 2,plat
         if( abs(flat_glob(j+jfirst-1) - flat(j)) .gt. eps ) then
            write(6,*) 'Error in "cremapbin":'
            write(6,*) "Gaussian latitudes in grid array do not"
            write(6,*) "match those in the global array"
            call abort
         end if
      end do
!
      return
    end subroutine binning_map_lats
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine binning_map_edges(plat     , plon , nlat     , jfirst , &
                                   flon     , flat , gw_glob  ,          &
                                   grid_flag, floni, flati    )
!
!--------1---------2---------3---------4---------5---------6---------7--
!
! Based on input grid, compute grid-box edges for either Gaussian or
! Regular (evenly-spaced) grids.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer plat                       ! latitude  dimension of input field
      integer plon                       ! longitude dimension of input field
      integer nlat                       ! Number of Global Gaussian latitudes
      integer jfirst                     ! Index of global Gaussian lat array mapped into the first grid lat array
!
      double precision flon     (plon)   ! longitudes in radians oriented W->E
      double precision flat     (plat)   ! latitudes  in radians oriented S->N
      double precision gw_glob  (nlat)   ! Global Gaussian weights
      integer grid_flag                  ! grid flags: 0=Regular, 1=Gaussian
      double precision floni    (plon+1) ! longitudes of box edges in radians oriented W->E
      double precision flati    (plat+1) ! latitudes  of box edges in radians oriented S->N
!
!---------------------------Local workspace-----------------------------
!
      integer i, j, platp1               ! Indices
      double precision sum
      double precision pi, pio2, half, one, two, three
!
!-----------------------------------------------------------------------
!
      platp1   = plat + 1
      half     = 0.5d0
      one      = 1.d0
      two      = 2.d0
      three    = 3.d0
      pi       = 4.d0*atan(one)
      pio2     = pi/two
!
! Compute longitudes of box edges
!
      floni(     1) = ( three*flon(   1) - flon(     2) )*half
      floni(plon+1) = ( three*flon(plon) - flon(plon-1) )*half
      do i = 2,plon
         floni(i) = half*(flon(i-1) + flon(i))
      end do
!
! If Regular grid, use algebraic mean to determine latitudes of box edges (extrapolation for endpoints)
! Else, if Gaussian grid, use partial sums of Gaussian weights.
!
      if(grid_flag .eq. 0) then
         flati(1     ) = ( three*flat(1) - flat(2) )*half
         flati(1     ) = max( flati(1), -pio2)
         flati(platp1) = ( three*flat(plat) - flat(plat-1) )*half
         flati(platp1) = min( flati(platp1),  pio2)
         do j = 1,plat-1
            flati(j+1) = half*(flat(j) + flat(j+1))
         end do
      else
!
! Sum Gaussian weights up to first latitude of data grid to get first box edge
!
         sum = 0.d0
         if(jfirst .le. 1) then
            flati(1) = -pio2
         else
            do j = 1,jfirst-1
               sum = sum + gw_glob(j)
            end do
            flati(1) = asin( sum-one )
         end if
!
! Determine subsequent box edges
!
         do j = 1,plat
            sum = sum + gw_glob(jfirst+j-1)
            flati(j+1) = asin( min (one,(sum-one) ) )
         end do
      end if
!
      return
    end subroutine binning_map_edges

      SUBROUTINE GAQDNCL(NLAT,THETA,WTS,WORK,LWORK,IERROR)
!     SUBROUTINE GAQDNCL COMPUTES GAUSSIAN POINTS (IN RADIANS) AND WEIGHTS
!     ON THE SPHERE IN THE INTERVAL (0,PI).  (THESE CAN BE USED IN
!     GAUSSIAN QUADRATURE FOR ESTIMATING INTEGRALS ON THE SPHERE)
!
!
!
!
      implicit none
      integer lwork,ierror,i1,i2,i3
      integer nlat
      integer n
      DIMENSION WORK(LWORK),THETA(NLAT),WTS(NLAT)
      DOUBLE PRECISION WORK,THETA,WTS,X
      N = NLAT
      IERROR = 1
!     CHECK WORK SPACE LENGTH
      IF (LWORK.LT.4*N*(N+1)+2) RETURN
      IERROR = 2
      IF (N.LE.0) RETURN
      IERROR = 0
      IF (N.GT.2) THEN
!     PARTITION WORK SPACE FOR DOUBLE PRECISION EIGENVALUE(VECTOR COMPUTATION)
      I1 = 1
      I2 = I1+2*N
      I3 = I2+2*N
      CALL GAQDNCL1(N,THETA,WTS,WORK(I1),WORK(I2),WORK(I3),IERROR)
      IF (IERROR.NE.0) THEN
      IERROR = 3
      RETURN
      END IF
      RETURN
      ELSE IF (N.EQ.1) THEN
      WTS(1) = 2.0D0
      THETA(1) = DACOS(0.0D0)
      ELSE IF (N.EQ.2) THEN
!     COMPUTE WEIGHTS AND POINTS ANALYTICALLY WHEN N=2
      WTS(1) = 1.0D0
      WTS(2) = 1.0D0
      X = DSQRT(1.0D0/3.0D0)
      THETA(1) = DACOS(X)
      THETA(2) = DACOS(-X)
      RETURN
      END IF
    END SUBROUTINE GAQDNCL
      SUBROUTINE GAQDNCL1(N,THETA,WTS,W,E,WRK,IER)
      integer n,j,matz,indx,np1,n2,i,ier
      DIMENSION THETA(N),WTS(N),W(N),E(N),WRK(1)
      DOUBLE PRECISION THETA,WTS,TEMP,W,E,WRK
!     SET SYMMETRIC TRIDIAGNONAL MATRIX SUBDIAGONAL AND DIAGONAL
!     COEFFICIENTS FOR MATRIX COMING FROM COEFFICIENTS IN THE
!     RECURSION FORMULA FOR LEGENDRE POLYNOMIALS
!     A(N)*P(N-1)+B(N)*P(N)+C(N)*P(N+1) = 0.
      WRK(1)=0.D0
      WRK(N+1) = 0.D0
      W(1)=0.D0
      E(1) = 0.D0
      DO 100 J=2,N
      WRK(J)= (J-1.D0)/DSQRT((2.D0*J-1.D0)*(2.D0*J-3.D0))
      WRK(J+N)=0.D0
      E(J) = WRK(J)
      W(J) = 0.D0
  100 CONTINUE
!     COMPUTE EIGENVALUES  OF MATRIX
      MATZ = 1
      INDX = 2*N+1
      NP1=N+1
      CALL DRSTNCL(N,N,W,E,MATZ,WRK(INDX),IER)
      IF (IER.NE.0) RETURN
!     COMPUTE GAUSSIAN WEIGHTS AND POINTS
      DO 101 J=1,N
      THETA(J) = DACOS(W(J))
!     SET GAUSSIAN WEIGHTS AS 1ST COMPONENTS OF EIGENVECTORS SQUARED
      INDX = (J+1)*N+1
      WTS(J) = 2.0D0*WRK(INDX)**2
  101 CONTINUE
!     REVERSE ORDER OF GAUSSIAN POINTS TO BE
!     MONOTONIC INCREASING IN RADIANS
      N2 = N/2
      DO 102 I=1,N2
      TEMP = THETA(I)
      THETA(I) = THETA(N-I+1)
      THETA(N-I+1) = TEMP
  102 CONTINUE
      RETURN
      END subroutine GAQDNCL1


      SUBROUTINE DRSTNCL(NM,N,W,E,MATZ,Z,IERR)
!     DRSTNCL IS A DOUBLE PRECISION MODIFICATION OF RST OFF EISPACK
!     TO BE USED  TO COMPUTE GAUSSIAN POINTS AND WEIGHTS

!
      INTEGER I,J,N,NM,IERR,MATZ
      DOUBLE PRECISION W(N),E(N),Z(NM,N)

!
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
!
	 DO 30 J = 1, N
	    Z(J,I) = 0.0D0
   30    CONTINUE
!
	 Z(I,I) = 1.0D0
   40 CONTINUE
!
      CALL  DINTQLNCL(NM,N,W,E,Z,IERR)
      RETURN
   END subroutine DRSTNCL

      SUBROUTINE DINTQLNCL(NM,N,D,E,Z,IERR)
!     DINTQLNCL IS A DOUBLE PRECISION MODIFICATION OF INTQL2 OFF
!     EISPACK TO BE USED BY GAQDNCL IN SPHEREPACK FOR COMPUTING
!     GAUSSIAN WEIGHTS AND POINTS
!
      INTEGER I,J,K,L,M,N,II,NM,MML,IERR,nm1,mdo
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
!      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,DPYTHANCL
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      E(N) = 0.0D0
!
      DO 240 L = 1, N
	 J = 0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
!   
  105    NM1 = N-1     
         IF(L .GT. NM1) GO TO 111
         DO 110 MDO = L, NM1
            M = MDO
	    TST1 = DABS(D(M)) + DABS(D(M+1))
	    TST2 = TST1 + DABS(E(M))
	    IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
  111    M = N
!
  120    P = D(L)
	 IF (M .EQ. L) GO TO 240
	 IF (J .EQ. 30) GO TO 1000
	 J = J + 1
!     .......... FORM SHIFT ..........
	 G = (D(L+1) - P) / (2.0D0 * E(L))
!	 R = DPYTHANCL(G,1.0D0)
	 call DPYTHANCL(G,1.0D0,R)
	 G = D(M) - P + E(L) / (G + SIGN(R,G))
	 S = 1.0D0
	 C = 1.0D0
	 P = 0.0D0
	 MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
	 DO 200 II = 1, MML
	    I = M - II
	    F = S * E(I)
	    B = C * E(I)
!	    R = DPYTHANCL(F,G)
	    call DPYTHANCL(F,G,R)
	    E(I+1) = R
	    IF (R .EQ. 0.0D0) GO TO 210
	    S = F / R
	    C = G / R
	    G = D(I+1) - P
	    R = (D(I) - G) * S + 2.0D0 * C * B
	    P = S * R
	    D(I+1) = G + P
	    G = C * R - B
!     .......... FORM VECTOR ..........
	    DO 180 K = 1, N
	       F = Z(K,I+1)
	       Z(K,I+1) = S * Z(K,I) + C * F
	       Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
!
  200    CONTINUE
!
	 D(L) = D(L) - P
	 E(L) = G
	 E(M) = 0.0D0
	 GO TO 105
!     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
	 E(M) = 0.0D0
	 GO TO 105
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
	 I = II - 1
	 K = I
	 P = D(I)
!
	 DO 260 J = II, N
	    IF (D(J) .GE. P) GO TO 260
	    K = J
	    P = D(J)
  260    CONTINUE
!
	 IF (K .EQ. I) GO TO 300
	 D(K) = D(I)
	 D(I) = P
!
	 DO 280 J = 1, N
	    P = Z(J,I)
	    Z(J,I) = Z(J,K)
	    Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
   END subroutine  DINTQLNCL

!!$  !===============================================================
!!$      DOUBLE PRECISION FUNCTION DPYTHANCL(A,B)
!!$      DOUBLE PRECISION A,B
!!$!     DPYTHANCL IS A DOUBLE PRECISION MODIFICATION OF PYTHAG OFF EISPACK
!!$!     FOR USE BY DIMTQL
!!$
!!$!
!!$!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!!$!
!!$      DOUBLE PRECISION P,R,S,T,U
!!$      P = DABS(A)
!!$      IF (DABS(B).GE.DABS(A)) P = DABS(B)
!!$!     P = AMAX1(DABS(A),DABS(B))
!!$      IF (P .EQ. 0.0D0) GO TO 20
!!$      R = (DABS(A)/P)**2
!!$      IF (DABS(B).LT.DABS(A)) R = (DABS(B)/P)**2
!!$!     R = (AMIN1(DABS(A),DABS(B))/P)**2
!!$   10 CONTINUE
!!$	 T = 4.0D0 + R
!!$	 IF (T .EQ. 4.0D0) GO TO 20
!!$	 S = R/T
!!$	 U = 1.0D0 + 2.0D0*S
!!$	 P = U*P
!!$	 R = (S/U)**2 * R
!!$      GO TO 10
!!$   20 DPYTHANCL = P
!!$      RETURN
!!$      END function DPYTHANCL
!!$  !====================================================================================
  !===============================================================
      subroutine DPYTHANCL(A,B,DPYRET)
      DOUBLE PRECISION A,B,DPYRET
!     DPYTHANCL IS A DOUBLE PRECISION MODIFICATION OF PYTHAG OFF EISPACK
!     FOR USE BY DIMTQL

!
!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!
      DOUBLE PRECISION P,R,S,T,U
      P = DABS(A)
      IF (DABS(B).GE.DABS(A)) P = DABS(B)
!     P = AMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DABS(A)/P)**2
      IF (DABS(B).LT.DABS(A)) R = (DABS(B)/P)**2
!     R = (AMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
	 T = 4.0D0 + R
	 IF (T .EQ. 4.0D0) GO TO 20
	 S = R/T
	 U = 1.0D0 + 2.0D0*S
	 P = U*P
	 R = (S/U)**2 * R
      GO TO 10
   20 DPYRET = P
      RETURN
    END subroutine DPYTHANCL
  !====================================================================================


  subroutine create_test_file1(ncid,varids)

    use netcdf

    integer,intent(out)   ::   ncid
    integer,dimension(10),intent(out)   ::   varids
    character(len=*),parameter :: subname='(create_test_file1)'

    ! error status return
    integer stat
    ! netCDF ncid
    integer  ind,tmpvarid

    ! dimension lengths
    integer latitude_len
    parameter (latitude_len = 360)
    integer longitude_len
    parameter (longitude_len = 720)
    integer time_len
    parameter (time_len = 25)
    ! dimension ids
    integer latitude_dim
    integer longitude_dim
    integer time_dim

    ! variable ids
    integer latitude_id
    integer longitude_id
    integer time_id
    integer var1_id
    integer var2_id
    integer var3_id


    ! rank (number of dimensions) for each variable
    integer latitude_rank
    parameter (latitude_rank = 1)
    integer longitude_rank
    parameter (longitude_rank = 1)
    integer time_rank
    parameter (time_rank = 1)
    integer var1_rank
    parameter (var1_rank = 3)
    integer var2_rank
    parameter (var2_rank = 3)
    integer var3_rank
    parameter (var3_rank = 3)

    ! variable shapes
    integer latitude_dims(latitude_rank)
    integer longitude_dims(longitude_rank)
    integer time_dims(time_rank)
    integer var1_dims(var1_rank)
    integer var2_dims(var2_rank)
    integer var3_dims(var3_rank)

    ! variable declarations
    ! attribute vectors
    integer textval(1)
    real realval(1)


    ! enter define mode
    call handle_ncerr( nf90_create('testfile.nc', nf90_clobber, ncid),subname,__LINE__)
    ! define dimensions
    call handle_ncerr( nf90_def_dim(ncid, 'lat', latitude_len, latitude_dim),subname,__LINE__)
    call handle_ncerr( nf90_def_dim(ncid, 'lon', longitude_len, longitude_dim),subname,__LINE__)
    call handle_ncerr( nf90_def_dim(ncid, 'time', time_len, time_dim),subname,__LINE__)

    ! define variables

    latitude_dims(1) = latitude_dim
    call handle_ncerr( nf90_def_var(ncid, 'lat', nf90_float, latitude_dims, latitude_id),subname,__LINE__)

    longitude_dims(1) = longitude_dim
    call handle_ncerr( nf90_def_var(ncid, 'lon', nf90_float, longitude_dims, longitude_id),subname,__LINE__)

    time_dims(1) = time_dim
    call handle_ncerr( nf90_def_var(ncid, 'time', nf90_float, time_dims, time_id),subname,__LINE__)

    var1_dims(1) = longitude_dim
    var1_dims(2) = latitude_dim
    var1_dims(3) = time_dim
    call handle_ncerr( nf90_def_var(ncid, 'var1', nf90_float, var1_dims, var1_id),subname,__LINE__)
    call handle_ncerr( nf90_def_var(ncid, 'var2', nf90_float, var1_dims, var2_id),subname,__LINE__)
    call handle_ncerr( nf90_def_var(ncid, 'var3', nf90_float, var1_dims, var3_id),subname,__LINE__)

    ! assign per-variable attributes
    ! define units
    call handle_ncerr( nf90_put_att(ncid, latitude_id, 'units', 'degrees_north'),subname,__LINE__)
    ! define _CoordinateAxisType
    call handle_ncerr( nf90_put_att(ncid, latitude_id, '_CoordinateAxisType', 'Lat'),subname,__LINE__)
    ! define comments
    call handle_ncerr( nf90_put_att(ncid, latitude_id, 'comments',  'center of cell'),subname,__LINE__)
    ! define units
    call handle_ncerr( nf90_put_att(ncid, longitude_id, 'units', 'degrees_east'),subname,__LINE__)
    ! define _CoordinateAxisType
    call handle_ncerr( nf90_put_att(ncid, longitude_id, '_CoordinateAxisType', 'Lon'),subname,__LINE__)
    ! define comments
    call handle_ncerr( nf90_put_att(ncid, longitude_id, 'comments',  'center of cell'),subname,__LINE__)
    ! define units
    call handle_ncerr( nf90_put_att(ncid, time_id, 'units',  'days since 0000-01-01 00:00'),subname,__LINE__)
    ! define _CoordinateAxisType
    call handle_ncerr( nf90_put_att(ncid, time_id, '_CoordinateAxisType', 'Time'),subname,__LINE__)
    ! define calendar
    call handle_ncerr( nf90_put_att(ncid, time_id, 'calendar', 'gregorian'),subname,__LINE__)
    ! define units
    call handle_ncerr( nf90_put_att(ncid, var1_id, 'units',  'TgC/grid/year'),subname,__LINE__)
    ! define missing_value
    realval(1) = -999.
    call handle_ncerr( nf90_put_att(ncid, var1_id, 'missing_value',realval),subname,__LINE__)
    ! define reference
    call handle_ncerr( nf90_put_att(ncid, var1_id, 'reference',  'JGCRI'),subname,__LINE__)
    ! define long_name
    call handle_ncerr( nf90_put_att(ncid, var1_id, 'long_name',  'RCP 4.5 Stab BC Aviation'),subname,__LINE__)
    ! define comments
    call handle_ncerr( nf90_put_att(ncid, var1_id, 'comments',  'Created by Yuyu Zhou at JGCRI'),subname,__LINE__)

    ! leave define mode
    call handle_ncerr( nf90_enddef(ncid),subname,__LINE__)
    varids=(/var1_id,var2_id,var3_id,0,0,0,0,0,0,0/)
    call handle_ncerr(nf90_inq_varid(ncid, "lon", tmpvarid),subname,__LINE__)
    call handle_ncerr( nf90_put_var(ncid, tmpvarid,lonbaseremap),subname,__LINE__)
    call handle_ncerr(nf90_inq_varid(ncid, "lat", tmpvarid),subname,__LINE__)
    call handle_ncerr( nf90_put_var(ncid, tmpvarid,latbaseremap),subname,__LINE__)

  end subroutine create_test_file1



end module gcam2emisfile_mod

