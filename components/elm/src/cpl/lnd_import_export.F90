module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type
  use lnd2atmType  , only: lnd2atm_type
  use lnd2glcMod   , only: lnd2glc_type
  use atm2lndType  , only: atm2lnd_type
  use glc2lndMod   , only: glc2lnd_type
  use GridcellType , only: grc_pp          ! for access to gridcell topology
  use TopounitDataType , only: top_as, top_af  ! atmospheric state and flux variables  
  use elm_cpl_indices
  use mct_mod
  !
  implicit none
  !===============================================================================

contains

  !===============================================================================
  subroutine lnd_import( bounds, x2l, atm2lnd_vars, glc2lnd_vars, lnd2atm_vars)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
    ! !USES:
    use elm_varctl       , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit, &
                                 metdata_type, metdata_bypass, metdata_biases, co2_file, aero_file, use_atm_downscaling_to_topunit
    use elm_varctl       , only: const_climate_hist, add_temperature, add_co2, use_cn, use_fates
    use elm_varctl       , only: startdate_add_temperature, startdate_add_co2
    use elm_varcon       , only: rair, o2_molar_const, c13ratio
    use clm_time_manager , only: get_nstep, get_step_size, get_curr_calday, get_curr_date 
    use controlMod       , only: NLFilename
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use domainMod        , only: ldomain
    use shr_kind_mod     , only: r8 => shr_kind_r8, CL => shr_kind_CL
    use fileutils        , only: getavu, relavu
    use spmdmod          , only: masterproc, mpicom, iam, npes, MPI_REAL8, MPI_INTEGER, MPI_STATUS_SIZE
    use elm_nlUtilsMod   , only : find_nlgroup_name
    use FrictionVelocityMod, only: implicit_stress, atm_gustiness
    use lnd_disagg_forc
    use lnd_downscale_atm_forcing
    use netcdf
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars      ! clm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars      ! clm internal input data type
    type(lnd2atm_type) , intent(in)    :: lnd2atm_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,topo,i,m,thism,nstep,ier  ! indices, number of steps, and error code
    integer status(MPI_STATUS_SIZE)
    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: e, ea                ! vapor pressure (Pa)
    real(r8) :: qsat                 ! saturation specific humidity (kg/kg)
    real(r8) :: forc_t               ! atmospheric temperature (Kelvin)
    real(r8) :: forc_q               ! atmospheric specific humidity (kg/kg)
    real(r8) :: forc_pbot            ! atmospheric pressure (Pa)
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: co2_ppmv_diag        ! temporary
    real(r8) :: co2_ppmv_prog        ! temporary
    real(r8) :: co2_ppmv_val         ! temporary
    integer  :: co2_type_idx         ! integer flag for co2_type options
    real(r8) :: esatw                ! saturation vapor pressure over water (Pa)
    real(r8) :: esati                ! saturation vapor pressure over ice (Pa)
    real(r8) :: a0,a1,a2,a3,a4,a5,a6 ! coefficients for esat over water
    real(r8) :: b0,b1,b2,b3,b4,b5,b6 ! coefficients for esat over ice
    real(r8) :: tdc, t               ! Kelvins to Celcius function and its input
    real(r8) :: vp                   ! water vapor pressure (Pa)
    integer  :: thisng, np, num, nu_nml, nml_error                 
    integer  :: ng_all(100000)
    real(r8) :: swndf, swndr, swvdf, swvdr, ratio_rvrf, frac, q
    real(r8) :: thiscosz, avgcosz, szenith
    integer  :: swrad_period_len, swrad_period_start, thishr, thismin
    real(r8) :: timetemp(2)
    real(r8) :: latixy(500000), longxy(500000)
    integer ::  ierr, varid, dimid, yr, mon, day, tod, nindex(2), caldaym(13)
    integer ::  ncid, met_ncids(14), mask_ncid, thisncid, ng, tm
    integer ::  aindex(2), tindex(14,2), starti(3), counti(3)
    integer ::  grid_map(500000), zone_map(500000)
    integer ::  met_nvars, nyears_spinup, nyears_trans, starti_site, endi_site
    real(r8) :: smap05_lat(360), smap05_lon(720)
    real(r8) :: smapt62_lat(94), smapt62_lon(192)
    real(r8) :: smap2_lat(96), smap2_lon(144)
    real(r8) :: thisdist, mindist, thislon
    real(r8) :: tbot, tempndep(1,1,158), thiscalday, wt1(14), wt2(14), thisdoy
    real(r8) :: site_metdata(14,12)
    real(r8) :: var_month_mean(12)
    !real(r8) :: hdm1(720,360,1), hdm2(720,360,1) 
    !real(r8) :: lnfm1(192,94,2920)
    !real(r8) :: ndep1(144,96,1), ndep2(144,96,1)
    !real(r8) :: aerodata(14,144,96,14)
    integer  :: lnfmind(2)
    integer  :: var_month_count(12)
    integer*2 :: temp(1,500000)
    integer :: xtoget, ytoget, thisx, thisy, calday_start
    integer :: sdate_addt, sy_addt, sm_addt, sd_addt
    integer :: sdate_addco2, sy_addco2, sm_addco2, sd_addco2
    character(len=200) metsource_str, thisline
    character(len=*), parameter :: sub = 'lnd_import_mct'
    integer :: av, v, n, nummetdims, g3, gtoget, ztoget, line, mystart, tod_start, thistimelen  
    character(len=20) aerovars(14), metvars(14)
    character(len=3) zst
    integer :: stream_year_first_lightng, stream_year_last_lightng, model_year_align_lightng
    integer :: stream_year_first_popdens, stream_year_last_popdens, model_year_align_popdens
    integer :: stream_year_first_ndep,    stream_year_last_ndep,    model_year_align_ndep
    character(len=CL)  :: metdata_fname  
    character(len=CL)  :: lightngmapalgo = 'bilinear'! Mapping alogrithm
    character(len=CL)  :: popdensmapalgo = 'bilinear' 
    character(len=CL)  :: ndepmapalgo    = 'bilinear' 
    character(len=CL)  :: stream_fldFileName_lightng ! lightning stream filename to read
    character(len=CL)  :: stream_fldFileName_popdens ! poplulation density stream filename
    character(len=CL)  :: stream_fldFileName_ndep    ! nitrogen deposition stream filename
    logical :: use_sitedata, has_zonefile, use_daymet, use_livneh
    data caldaym / 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /    

    ! Constants to compute vapor pressure
    parameter (a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8)

    parameter (b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8)
    !
    ! function declarations
    !
    tdc(t) = min( 50._r8, max(-50._r8,(t-SHR_CONST_TKFRZ)) )
    esatw(t) = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    esati(t) = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    !---------------------------------------------------------------------------

    namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

    namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

    namelist /ndepdyn_nml/        &
        stream_year_first_ndep,  &
    stream_year_last_ndep,   &
        model_year_align_ndep,   &
        ndepmapalgo,             &
        stream_fldFileName_ndep

    stream_fldFileName_lightng = ' '
    stream_fldFileName_popdens = ' '
   
    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

    thisng = bounds%endg - bounds%begg + 1
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       
       ! Determine flooding input, sign convention is positive downward and
       ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
       ! change the sign to indicate addition of water to system.

       atm2lnd_vars%forc_flood_grc(g)   = -x2l(index_x2l_Flrr_flood,i)  

       atm2lnd_vars%volr_grc(g)   = x2l(index_x2l_Flrr_volr,i) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%volrmch_grc(g)= x2l(index_x2l_Flrr_volrmch,i) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%supply_grc(g) = x2l(index_x2l_Flrr_supply,i)
       atm2lnd_vars%deficit_grc(g) = x2l(index_x2l_Flrr_deficit,i)

       if (index_x2l_Sr_h2orof /= 0) then
         atm2lnd_vars%h2orof_grc(g)      = x2l(index_x2l_Sr_h2orof,i)
         atm2lnd_vars%frac_h2orof_grc(g) = x2l(index_x2l_Sr_frac_h2orof,i)
       endif

       ! Determine required receive fields

#ifdef CPL_BYPASS
        !read forcing data directly, bypass coupler
        atm2lnd_vars%forc_flood_grc(g)   = 0._r8
        atm2lnd_vars%volr_grc(g)   = 0._r8 

        !Get meteorological data, concatenated to include whole record
        !Note we only do this at the first timestep and keep the whole forcing dataset in the memory
       
  !-----------------------------------Meteorological forcing  -----------------------------------

        call get_curr_date( yr, mon, day, tod )
        thiscalday = get_curr_calday()
        nstep = get_nstep()

        !on first timestep, read all the met data for relevant gridcell(s) and store in array.
        !   Met data are held in short integer format to save memory.
        !   Each node must have enough memory to hold these data.
        met_nvars=7
        if (metdata_type(1:3) == 'cpl') met_nvars=14

        if (atm2lnd_vars%loaded_bypassdata == 0) then
          !meteorological forcing
          if (index(metdata_type, 'qian') .gt. 0) then 
            atm2lnd_vars%metsource = 0   
          else if (index(metdata_type,'cru') .gt. 0) then
            atm2lnd_vars%metsource = 1  
          else if (index(metdata_type,'site') .gt. 0) then 
            atm2lnd_vars%metsource = 2
          else if (index(metdata_type,'princeton') .gt. 0) then 
            atm2lnd_vars%metsource = 3
          else if (index(metdata_type,'gswp3') .gt. 0) then
            atm2lnd_vars%metsource = 4
          else if (index(metdata_type,'cpl') .gt. 0) then 
            atm2lnd_vars%metsource = 5
          else
            call endrun( sub//' ERROR: Invalid met data source for cpl_bypass' )
          end if

          use_livneh = .false.
          use_daymet = .false.
          if(index(metdata_type, 'livneh') .gt. 0) then 
              use_livneh = .true.
          else if (index(metdata_type, 'daymet') .gt. 0) then 
              use_daymet = .true.
          end if
 
          metvars(1) = 'TBOT'
          metvars(2) = 'PSRF'
          metvars(3) = 'QBOT'
          if (atm2lnd_vars%metsource .eq. 2) metvars(3) = 'RH'
          if (atm2lnd_vars%metsource .ne. 5) metvars(4) = 'FSDS'
          if (atm2lnd_vars%metsource .ne. 5) metvars(5) = 'PRECTmms'
          if (atm2lnd_vars%metsource .ne. 5) metvars(6) = 'WIND'
          metvars(4) = 'FSDS'
          metvars(5) = 'PRECTmms'
          metvars(6) = 'WIND'
          metvars(7) = 'FLDS'
          if (atm2lnd_vars%metsource .eq. 5) then 
              metvars(4) = 'SWNDF'
              metvars(5) = 'RAINC'
              metvars(6) = 'U'
              metvars(8) = 'SWNDR'
              metvars(9) = 'SWVDF'
              metvars(10) = 'SWVDR'
              metvars(11) = 'RAINL'
              metvars(12) = 'SNOWC'
              metvars(13) = 'SNOWL'
              metvars(14) = 'V'
          else
              metvars(4) = 'FSDS'
              metvars(5) = 'PRECTmms'
              metvars(6) = 'WIND'
          end if

          !set defaults
          atm2lnd_vars%startyear_met       = 1901
          atm2lnd_vars%endyear_met_spinup  = 1920
          if (atm2lnd_vars%metsource == 0) then 
            metsource_str = 'qian'
            atm2lnd_vars%startyear_met       = 1948
            atm2lnd_vars%endyear_met_spinup  = 1972
            atm2lnd_vars%endyear_met_trans   = 2004
          else if (atm2lnd_vars%metsource == 1) then 
            metsource_str = 'cruncep'
            atm2lnd_vars%endyear_met_trans  = 2016
          else if (atm2lnd_vars%metsource == 2) then
            metsource_str = 'site'
            !get year information from file
            ierr = nf90_open(trim(metdata_bypass) // '/all_hourly.nc', nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'start_year', varid) 
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%startyear_met)
            ierr = nf90_inq_varid(ncid, 'end_year', varid)
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%endyear_met_spinup)
            ierr = nf90_close(ncid)
            atm2lnd_vars%endyear_met_trans = atm2lnd_vars%endyear_met_spinup
          else if (atm2lnd_vars%metsource == 3) then 
            metsource_str = 'princeton'
            atm2lnd_vars%endyear_met_trans = 2012 
          else if (atm2lnd_vars%metsource == 4) then 
            atm2lnd_vars%endyear_met_trans  = 2014
          else if (atm2lnd_vars%metsource == 5) then
            atm2lnd_vars%startyear_met      = 566 !76
            atm2lnd_vars%endyear_met_spinup = 590 !100
            atm2lnd_vars%endyear_met_trans  = 590 !100
          end if

          if (use_livneh) then 
              atm2lnd_vars%startyear_met      = 1950
              atm2lnd_vars%endyear_met_spinup = 1969
          else if (use_daymet) then 
              atm2lnd_vars%startyear_met      = 1980
              atm2lnd_vars%endyear_met_spinup = atm2lnd_vars%endyear_met_trans
          end if

          nyears_spinup = atm2lnd_vars%endyear_met_spinup - &
                             atm2lnd_vars%startyear_met + 1
          nyears_trans  = atm2lnd_vars%endyear_met_trans - &
                             atm2lnd_vars%startyear_met  + 1

          !check for site data in run directory (monthly mean T, precip)
          inquire(file=trim(metdata_biases), exist=use_sitedata)

          !get grid lat/lon information, zone mappings
          inquire(file=trim(metdata_bypass) // '/zone_mappings.txt', exist=has_zonefile)
          if (has_zonefile) then
            open(unit=13, file=trim(metdata_bypass) // '/zone_mappings.txt')
          else if (atm2lnd_vars%metsource .ne. 2) then
            call endrun( sub//' ERROR: Zone mapping file does not exist for cpl_bypass' )
          end if

          if (atm2lnd_vars%metsource .ne. 2) then 
            ng = 0     !number of points
            do v=1,500000
              read(13,*, end=10), longxy(v), latixy(v), zone_map(v), grid_map(v)
              ng = ng + 1
            end do
10          continue
            close(unit=13)

            !Figure out the closest point and which zone file to open
            mindist=99999
            do g3 = 1,ng
              thisdist = 100*((latixy(g3) - ldomain%latc(g))**2 + &
                              (longxy(g3) - ldomain%lonc(g))**2)**0.5
              if (thisdist .lt. mindist) then 
                mindist = thisdist
                ztoget = zone_map(g3)
                gtoget = grid_map(g3)
              end if
            end do
          else
            gtoget = 1
          end if

          !get the site metdata for bias correction if they exist (lat/lons must match domain file)
          if (use_sitedata) then 
            open(unit=9, file=trim(metdata_biases),status='old')
            read(9,*) thisline
            site_metdata(:,:)=-999._r8
            do while ((site_metdata(1,1) .lt. ldomain%lonc(g) - 0.01 .or. &
                    site_metdata(1,1) .gt. ldomain%lonc(g) + 0.01) .and. &
                      (site_metdata(2,1) .lt. ldomain%latc(g) - 0.01 .or. &
                       site_metdata(2,1) .gt. ldomain%latc(g) + 0.01))
              read(9,*) site_metdata(1:7,1)
              if (site_metdata(1,1) .lt. 0) site_metdata(1,1) = site_metdata(1,1)+360._r8
            end do
            do line=2,12
              read(9,*) site_metdata(1:7,line)
            end do
            close(unit=9)
          end if

          do v=1,met_nvars
            write(zst, '(I3)') 100+ztoget
            if (atm2lnd_vars%metsource == 0) then 
                metdata_fname =  trim(metsource_str) // '_' // trim(metvars(v)) // '_z' // zst(2:3) // '.nc'
            else if (atm2lnd_vars%metsource == 1) then 
                metdata_fname = 'CRUNCEP.v5_' // trim(metvars(v)) // '_1901-2013_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'CRUNCEP5_Livneh_' // trim(metvars(v)) // '_1950-2013_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'CRUNCEP5_Daymet3_' // trim(metvars(v)) // '_1980-2013_z' // zst(2:3) // '.nc'
                end if
            else if (atm2lnd_vars%metsource == 2) then
                metdata_fname = 'all_hourly.nc'
            else if (atm2lnd_vars%metsource == 3) then 
               metdata_fname = 'Princeton_' // trim(metvars(v)) // '_1901-2012_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Livneh_' // trim(metvars(v)) // '_1950-2012_z' // zst(2:3) // '.nc'
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
                    metdata_fname = 'Princeton_Daymet3_' // trim(metvars(v)) // '_1980-2012_z' // zst(2:3) // '.nc'
                end if
            else if (atm2lnd_vars%metsource == 4) then 
                metdata_fname = 'GSWP3_' // trim(metvars(v)) // '_1901-2014_z' // zst(2:3) // '.nc'
                if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'GSWP3_Livneh_' // trim(metvars(v)) // '_1950-2010_z' // zst(2:3) // '.nc'                
                else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then 
                    metdata_fname = 'GSWP3_Daymet3_' // trim(metvars(v)) // '_1980-2010_z' // zst(2:3) // '.nc' 
                end if
            else if (atm2lnd_vars%metsource == 5) then 
                    !metdata_fname = 'WCYCL1850S.ne30_' // trim(metvars(v)) // '_0076-0100_z' // zst(2:3) // '.nc'
                    metdata_fname = 'CBGC1850S.ne30_' // trim(metvars(v)) // '_0566-0590_z' // zst(2:3) // '.nc'
            end if
  
            ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, met_ncids(v))
            if (ierr .ne. 0) call endrun(msg=' ERROR: Failed to open cpl_bypass input meteorology file' )
       
            !get timestep information
            ierr = nf90_inq_dimid(met_ncids(v), 'DTIME', dimid)
            ierr = nf90_Inquire_Dimension(met_ncids(v), dimid, len = atm2lnd_vars%timelen(v))

            starti(1) = 1
            counti(1) = 2
            ierr = nf90_inq_varid(met_ncids(v), 'DTIME', varid)
            ierr = nf90_get_var(met_ncids(v), varid, timetemp, starti(1:1), counti(1:1))   
            atm2lnd_vars%timeres(v)        = (timetemp(2)-timetemp(1))*24._r8
            atm2lnd_vars%npf(v)            = 86400d0*(timetemp(2)-timetemp(1))/get_step_size()  
            atm2lnd_vars%timelen_spinup(v) = nyears_spinup*(365*nint(24./atm2lnd_vars%timeres(v)))
    
            ierr = nf90_inq_varid(met_ncids(v), trim(metvars(v)), varid)
            !get the conversion factors
            ierr = nf90_get_att(met_ncids(v), varid, 'scale_factor', atm2lnd_vars%scale_factors(v))
            ierr = nf90_get_att(met_ncids(v), varid, 'add_offset', atm2lnd_vars%add_offsets(v))
            !get the met data         
            starti(1) = 1
            starti(2) = gtoget
            counti(1) = atm2lnd_vars%timelen_spinup(v)
            counti(2) = 1
            if (.not. const_climate_hist .and. (yr .ge. 1850 .or. use_sitedata)) counti(1) = atm2lnd_vars%timelen(v)

            if (i == 1 .and. v == 1)  then 
              allocate(atm2lnd_vars%atm_input       (met_nvars,bounds%begg:bounds%endg,1,1:counti(1)))
            end if 

            ierr = nf90_get_var(met_ncids(v), varid, atm2lnd_vars%atm_input(v,g:g,1,1:counti(1)), starti(1:2), counti(1:2))
            ierr = nf90_close(met_ncids(v))
    
            if (use_sitedata .and. v == 1) then 
                starti_site = max((nint(site_metdata(4,1))-atm2lnd_vars%startyear_met) * &
                                     365*nint(24./atm2lnd_vars%timeres(v))+1,1)
                endi_site   = (min(atm2lnd_vars%endyear_met_trans,nint(site_metdata(5,1))) - &
                                     atm2lnd_vars%startyear_met+1)*(365*nint(24./atm2lnd_vars%timeres(v)))
            end if
             
            atm2lnd_vars%var_offset(v,g,:) = 0._r8
            atm2lnd_vars%var_mult(v,g,:)   = 1._r8

            if (use_sitedata) then 
              !Compute monthly biases for site vs. reanalysis
              var_month_mean(:)  = 0._r8
              var_month_count(:) = 0
              do i=starti_site, endi_site
                thisdoy = mod(i,365*nint(24./atm2lnd_vars%timeres(v)))/(nint(24./atm2lnd_vars%timeres(v)))+1
                do m=1,12
                  if (thisdoy .ge. caldaym(m) .and. thisdoy .lt. caldaym(m+1)) thism = m
                enddo
                var_month_mean(thism) = var_month_mean(thism) + (atm2lnd_vars%atm_input(v,g,1,i)* &
                                          atm2lnd_vars%scale_factors(v) + atm2lnd_vars%add_offsets(v))
                var_month_count(thism) = var_month_count(thism)+1
              end do
     
              do m = 1,12
                var_month_mean(m) = var_month_mean(m)/var_month_count(m)
                !calculate offset and linear bias factors for temperature and precipitation
                if (v .eq. 1) atm2lnd_vars%var_offset(v,g,m) = (site_metdata(6,m)+SHR_CONST_TKFRZ) - var_month_mean(m)
                if (v .eq. 5 .and. var_month_mean(m) .gt. 0) &     
                      atm2lnd_vars%var_mult(v,g,m) = (site_metdata(7,m))/(caldaym(m+1)-caldaym(m))/24._r8/ &
                                                      3600._r8 / var_month_mean(m)
              end do
            end if
        
            !Align spinups and transient simulations
            !figure out which year to start with (assuming spinups always use integer multiple of met cycles)
            mystart = atm2lnd_vars%startyear_met
            do while (mystart > 1850)
              mystart = mystart - nyears_spinup
            end do
            if (atm2lnd_vars%metsource == 5) mystart=1850

            if (yr .lt. 1850) then 
              atm2lnd_vars%tindex(g,v,1) = (mod(yr-1,nyears_spinup) + (1850-mystart)) * 365 * nint(24./atm2lnd_vars%timeres(v))
            else if (yr .le. atm2lnd_vars%endyear_met_spinup) then
              atm2lnd_vars%tindex(g,v,1) = (mod(yr-1850,nyears_spinup) + (1850-mystart)) * 365 * nint(24./atm2lnd_vars%timeres(v))
            else
              atm2lnd_vars%tindex(g,v,1) = (yr - atm2lnd_vars%startyear_met) * 365 * nint(24./atm2lnd_vars%timeres(v))
            end if
            !adjust for starts not at beginning of year (but currently MUST begin at hour 0)
            atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1) + (caldaym(mon)+day-2)* &
                                         nint(24./atm2lnd_vars%timeres(v))
            
            atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,1) + 1
            if (atm2lnd_vars%tindex(g,v,1) == 0) then 
              atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen(v)
              if (yr .le. atm2lnd_vars%endyear_met_spinup) atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen_spinup(v)
             end if
          end do    !end variable loop        
        else
          do v=1,met_nvars
            if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then 
              if (v .eq. 4 .or. v .eq. 5 .or. (v .ge. 8 .and. v .le. 13)) then    !rad/Precipitation
                if (mod(tod/get_step_size(),nint(atm2lnd_vars%npf(v))) == 1 .and. nstep .gt. 3) then
                  atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+1
                  atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+1
                end if
              else  
                if (mod(tod/get_step_size()-1,nint(atm2lnd_vars%npf(v))) <= atm2lnd_vars%npf(v)/2._r8 .and. &
                    mod(tod/get_step_size(),nint(atm2lnd_vars%npf(v))) > atm2lnd_vars%npf(v)/2._r8) then 
                  atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+1
                  atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+1
                end if
              end if
            else
              atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%tindex(g,v,1)+nint(1/atm2lnd_vars%npf(v))
              atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%tindex(g,v,2)+nint(1/atm2lnd_vars%npf(v))  
            end if

            if (const_climate_hist .or. yr .le. atm2lnd_vars%startyear_met) then
              if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,1) = 1
              if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,2) = 1
            else if (yr .gt. atm2lnd_vars%endyear_met_trans) then
              if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen(v)) then
                 atm2lnd_vars%tindex(g,v,1) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
              end if
              if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen(v)) then
                 atm2lnd_vars%tindex(g,v,2) = atm2lnd_vars%timelen(v)-atm2lnd_vars%timelen_spinup(v)+1
              end if
            end if

            !if (yr .gt. atm2lnd_vars%startyear_met) then 
            !  if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen(v)) atm2lnd_vars%tindex(g,v,1) = 1
            !  if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen(v)) atm2lnd_vars%tindex(g,v,2) = 1
            !else
            !  if (atm2lnd_vars%tindex(g,v,1) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,1) = 1
            !  if (atm2lnd_vars%tindex(g,v,2) .gt. atm2lnd_vars%timelen_spinup(v)) atm2lnd_vars%tindex(g,v,2) = 1
            !end if
          end do
        end if

        tindex = atm2lnd_vars%tindex(g,:,:)

        !get weights for linear interpolation 
        do v=1,met_nvars
          if (atm2lnd_vars%npf(v) - 1._r8 .gt. 1e-3) then
               wt1(v) = 1._r8 - (mod((tod+86400)/get_step_size()-atm2lnd_vars%npf(v)/2._r8, &
                   atm2lnd_vars%npf(v))*1._r8)/atm2lnd_vars%npf(v)
               wt2(v) = 1._r8 - wt1(v)
          else
             wt1(v) = 0._r8    
             wt2(v) = 1._r8
          end if
        end do

        !Air temperature
        atm2lnd_vars%forc_t_not_downscaled_grc(g)  = min(((atm2lnd_vars%atm_input(1,g,1,tindex(1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(1,2))* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)             
        atm2lnd_vars%forc_th_not_downscaled_grc(g) = min(((atm2lnd_vars%atm_input(1,g,1,tindex(1,1))*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,g,1,tindex(1,2))* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)
       
        tbot = atm2lnd_vars%forc_t_not_downscaled_grc(g)

        !Air pressure
        atm2lnd_vars%forc_pbot_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(2,g,1,tindex(2,1))*atm2lnd_vars%scale_factors(2)+ &
                                                        atm2lnd_vars%add_offsets(2))*wt1(2) + (atm2lnd_vars%atm_input(2,g,1,tindex(2,2)) &
                                                        *atm2lnd_vars%scale_factors(2)+atm2lnd_vars%add_offsets(2))*wt2(2)) * &
                                                        atm2lnd_vars%var_mult(2,g,mon) + atm2lnd_vars%var_offset(2,g,mon), 4e4_r8)       
        !Specific humidity
        atm2lnd_vars%forc_q_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(3,g,1,tindex(3,1))*atm2lnd_vars%scale_factors(3)+ &
                                                     atm2lnd_vars%add_offsets(3))*wt1(3) + (atm2lnd_vars%atm_input(3,g,1,tindex(3,2)) &
                                                     *atm2lnd_vars%scale_factors(3)+atm2lnd_vars%add_offsets(3))*wt2(3)) * &
                                                     atm2lnd_vars%var_mult(3,g,mon) + atm2lnd_vars%var_offset(3,g,mon), 1e-9_r8)

        if (atm2lnd_vars%metsource == 2) then  !convert RH to qbot                             
          if (tbot > SHR_CONST_TKFRZ) then
            e = esatw(tdc(tbot))
          else
            e = esati(tdc(tbot))
          end if
          qsat           = 0.622_r8*e / (atm2lnd_vars%forc_pbot_not_downscaled_grc(g) - 0.378_r8*e)
          atm2lnd_vars%forc_q_not_downscaled_grc(g) = qsat * atm2lnd_vars%forc_q_not_downscaled_grc(g) / 100.0_r8
        end if

        !use longwave from file if provided
        atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ((atm2lnd_vars%atm_input(7,g,1,tindex(7,1))*atm2lnd_vars%scale_factors(7)+ &
                                                        atm2lnd_vars%add_offsets(7))*wt1(7) + (atm2lnd_vars%atm_input(7,g,1,tindex(7,2)) &
                                                        *atm2lnd_vars%scale_factors(7)+atm2lnd_vars%add_offsets(7))*wt2(7)) * &
                                                        atm2lnd_vars%var_mult(7,g,mon) + atm2lnd_vars%var_offset(7,g,mon)  
        if (atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .le. 50 .or. atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .ge. 600) then 
        !Longwave radiation (calculated from air temperature, humidity)
            e =  atm2lnd_vars%forc_pbot_not_downscaled_grc(g) * atm2lnd_vars%forc_q_not_downscaled_grc(g) / &
                 (0.622_R8 + 0.378_R8 * atm2lnd_vars%forc_q_not_downscaled_grc(g) )
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ea * SHR_CONST_STEBOL * tbot**4
        end if 

        !Shortwave radiation (cosine zenith angle interpolation)
        thishr = (tod-get_step_size()/2)/3600
        if (thishr < 0) thishr=thishr+24
        thismin = mod((tod-get_step_size()/2)/60, 60)
        thiscosz = max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr,thismin,0)* &
                        3.14159265358979/180.0d0), 0.001d0)
        avgcosz = 0d0
        if (atm2lnd_vars%npf(4) - 1._r8 .gt. 1e-3) then 
          swrad_period_len   = get_step_size()*nint(atm2lnd_vars%npf(4))
          swrad_period_start = ((tod-get_step_size()/2)/swrad_period_len) * swrad_period_len
          !set to last period if first model timestep of the day
          if (tod-get_step_size()/2 < 0) swrad_period_start = ((86400-get_step_size()/2)/swrad_period_len) * swrad_period_len   

          do tm=1,nint(atm2lnd_vars%npf(4))  
            !Get the average cosine zenith angle over the time resolution of the input data
            thishr  = (swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/3600
            if (thishr > 23) thishr=thishr-24  
            thismin = mod((swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/60, 60) 
            avgcosz  = avgcosz + max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr, thismin, 0) &
                       *3.14159265358979/180.0d0), 0.001d0)/atm2lnd_vars%npf(4)
          end do
        else
          avgcosz = thiscosz
        end if
        if (thiscosz > 0.001d0) then 
          wt2(4) = min(thiscosz/avgcosz, 10.0_r8)
        else
          wt2(4) = 0d0
        end if
        
        if (atm2lnd_vars%metsource == 5) then 
            wt2(4)=1.0   !cosz interp not working 
            wt2(8:10)=1.0
            swndf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                     atm2lnd_vars%add_offsets(4))*wt2(4)), 0.0_r8)
            swndr = max(((atm2lnd_vars%atm_input(8,g,1,tindex(8,2))*atm2lnd_vars%scale_factors(8)+ &
                                     atm2lnd_vars%add_offsets(8))*wt2(8)), 0.0_r8)
            swvdf = max(((atm2lnd_vars%atm_input(9,g,1,tindex(9,2))*atm2lnd_vars%scale_factors(9)+ &
                                     atm2lnd_vars%add_offsets(9))*wt2(9)), 0.0_r8)
            swvdr = max(((atm2lnd_vars%atm_input(10,g,1,tindex(10,2))*atm2lnd_vars%scale_factors(10)+ &
                                     atm2lnd_vars%add_offsets(10))*wt2(10)), 0.0_r8)
            atm2lnd_vars%forc_solad_grc(g,2) = swndr
            atm2lnd_vars%forc_solad_grc(g,1) = swvdr
            atm2lnd_vars%forc_solai_grc(g,2) = swndf
            atm2lnd_vars%forc_solai_grc(g,1) = swvdf
        else
            swndr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                     atm2lnd_vars%add_offsets(4))*wt2(4)) * 0.50_R8, 0.0_r8)
            swndf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdr = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            swvdf = max(((atm2lnd_vars%atm_input(4,g,1,tindex(4,2))*atm2lnd_vars%scale_factors(4)+ &
                                    atm2lnd_vars%add_offsets(4))*wt2(4))*0.50_R8, 0.0_r8)
            ratio_rvrf =   min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr &
                           -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
            atm2lnd_vars%forc_solad_grc(g,2) = ratio_rvrf*swndr
            atm2lnd_vars%forc_solai_grc(g,2) = (1._R8 - ratio_rvrf)*swndf
            ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                               -9.0039e-06_R8*swvdr**2 +8.1351e-09_R8*swvdr**3,0.01_R8))
            atm2lnd_vars%forc_solad_grc(g,1) = ratio_rvrf*swvdr
            atm2lnd_vars%forc_solai_grc(g,1) = (1._R8 - ratio_rvrf)*swvdf
        end if
        !Rain and snow
        if (atm2lnd_vars%metsource == 5) then 
          forc_rainc = max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = max((((atm2lnd_vars%atm_input(11,g,1,tindex(11,2))*atm2lnd_vars%scale_factors(11)+ &
                                        atm2lnd_vars%add_offsets(11)))*atm2lnd_vars%var_mult(11,g,mon) + &
                                        atm2lnd_vars%var_offset(11,g,mon)), 0.0_r8)
          forc_snowc = max((((atm2lnd_vars%atm_input(12,g,1,tindex(12,2))*atm2lnd_vars%scale_factors(12)+ &
                                        atm2lnd_vars%add_offsets(12)))*atm2lnd_vars%var_mult(12,g,mon) + &
                                        atm2lnd_vars%var_offset(12,g,mon)), 0.0_r8)
          forc_snowl = max((((atm2lnd_vars%atm_input(13,g,1,tindex(13,2))*atm2lnd_vars%scale_factors(13)+ &
                                        atm2lnd_vars%add_offsets(13)))*atm2lnd_vars%var_mult(13,g,mon) + &
                                          atm2lnd_vars%var_offset(13,g,mon)), 0.0_r8)
        else
          frac = (atm2lnd_vars%forc_t_not_downscaled_grc(g) - SHR_CONST_TKFRZ)*0.5_R8       ! ramp near freezing
          frac = min(1.0_R8,max(0.0_R8,frac))           ! bound in [0,1]
          !Don't interpolate rainfall data
          forc_rainc = 0.1_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = 0.9_R8 * frac * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                                         atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                         atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8) 
          forc_snowc = 0.1_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)  
          forc_snowl = 0.9_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,g,1,tindex(5,2))*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5))) * atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8) 
        end if
        !Wind
        atm2lnd_vars%forc_u_grc(g) = (atm2lnd_vars%atm_input(6,g,1,tindex(6,1))*atm2lnd_vars%scale_factors(6)+ &
                                     atm2lnd_vars%add_offsets(6))*wt1(6) + (atm2lnd_vars%atm_input(6,g,1,tindex(6,2))* &
                                     atm2lnd_vars%scale_factors(6)+atm2lnd_vars%add_offsets(6))*wt2(6)
        if (atm2lnd_vars%metsource == 5) then 
          atm2lnd_vars%forc_v_grc(g) = (atm2lnd_vars%atm_input(14,g,1,tindex(14,1))*atm2lnd_vars%scale_factors(14)+ &
                                     atm2lnd_vars%add_offsets(14))*wt1(14) + (atm2lnd_vars%atm_input(14,g,1,tindex(14,2))* &
                                     atm2lnd_vars%scale_factors(14)+atm2lnd_vars%add_offsets(14))*wt2(14)
        else
            atm2lnd_vars%forc_v_grc(g) = 0.0_R8 
        end if
        atm2lnd_vars%forc_hgt_grc(g) = 30.0_R8 !(atm2lnd_vars%atm_input(8,g,1,tindex(1))*wt1 + &
                                             !atm2lnd_vars%atm_input(8,g,1,tindex(2))*wt2)    ! zgcmxy  Atm state, default=30m

  !------------------------------------Fire data -------------------------------------------------------
 
        nindex(1) = yr-1848
        nindex(2) = nindex(1)+1
        if (yr .lt. 1850 .or. const_climate_hist) nindex(1:2) = 2
        if (yr .ge. 2010 .and. .not. const_climate_hist) nindex(1:2) = 161
      
        model_filter: if (use_cn .or. use_fates) then 
          if (atm2lnd_vars%loaded_bypassdata == 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then  
            if (masterproc .and. i .eq. 1) then 
              ! Read pop_dens streams namelist to get filename
              nu_nml = getavu()
              open(nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
              if (nml_error == 0) then
                  read(nu_nml, nml=popd_streams,iostat=nml_error)
                  if (nml_error /= 0) then
                      call endrun(msg='ERROR reading popdens namelist')
                  end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_popdens), NF90_NOWRITE, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lon)
              ierr = nf90_inq_varid(ncid, 'hdm', varid)
              starti(1:2) = 1 
              starti(3)   = nindex(1)
              counti(1) = 720
              counti(2) = 360
              counti(3) = 1       
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm1, starti, counti)
              starti(3) = nindex(2)
              if (nindex(1) .ne. nindex(2)) then 
                  ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm2, starti, counti)
              else
                  atm2lnd_vars%hdm2 = atm2lnd_vars%hdm1 
              end if
              ierr = nf90_close(ncid)
            end if

            if (i .eq. 1) then 
              call mpi_bcast (atm2lnd_vars%hdm1, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (atm2lnd_vars%hdm2, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lon, 720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lat, 360, MPI_REAL8, 0, mpicom, ier)
            end if
          end if

          !figure out which point to get
          if (atm2lnd_vars%loaded_bypassdata == 0) then 
            mindist=99999
            do thisx = 1,720
              do thisy = 1,360
                  if (ldomain%lonc(g) .lt. 0) then
                      if (smap05_lon(thisx) >= 180) smap05_lon(thisx) = smap05_lon(thisx)-360._r8
                  else if (ldomain%lonc(g) .ge. 180) then
                      if (smap05_lon(thisx) < 0) smap05_lon(thisx) = smap05_lon(thisx) + 360._r8
                  end if
                  thisdist = 100*((smap05_lat(thisy) - ldomain%latc(g))**2 + &
                          (smap05_lon(thisx) - ldomain%lonc(g))**2)**0.5
                  if (thisdist .lt. mindist) then
                      mindist = thisdist
                      atm2lnd_vars%hdmind(g,1) = thisx
                      atm2lnd_vars%hdmind(g,2) = thisy
                  end if
              end do
            end do
          end if
          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
          atm2lnd_vars%forc_hdm(g) = atm2lnd_vars%hdm1(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt1(1) + &
                                     atm2lnd_vars%hdm2(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt2(1)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. masterproc .and. i .eq. 1) then 
            ! Read light_streams namelist to get filename
            nu_nml = getavu()
            open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
            call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
            if (nml_error == 0) then
              read(nu_nml, nml=light_streams,iostat=nml_error)
              if (nml_error /= 0) then
                call endrun(msg='ERROR reading light_streams namelist')
              end if
            end if
            close(nu_nml)
            call relavu( nu_nml )

            !Get all of the data (master processor only)
            allocate(atm2lnd_vars%lnfm_all       (192,94,2920))
            ierr = nf90_open(trim(stream_fldFileName_lightng), NF90_NOWRITE, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)
            ierr = nf90_get_var(ncid, varid, smapt62_lon)
            ierr = nf90_inq_varid(ncid, 'lnfm', varid)
            ierr = nf90_get_var(ncid, varid, atm2lnd_vars%lnfm_all)
            ierr = nf90_close(ncid)
          end if
          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. i .eq. 1) then
            call mpi_bcast (smapt62_lon, 192, MPI_REAL8, 0, mpicom, ier)
            call mpi_bcast (smapt62_lat, 94, MPI_REAL8, 0, mpicom, ier)
          end if
          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,192
              do thisy = 1,94
                if (ldomain%lonc(g) .lt. 0) then 
                  if (smapt62_lon(thisx) >= 180) smapt62_lon(thisx) = smapt62_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then 
                  if (smapt62_lon(thisx) < 0) smapt62_lon(thisx) = smapt62_lon(thisx) + 360._r8
                end if
                thisdist = 100*((smapt62_lat(thisy) - ldomain%latc(g))**2 + &
                            (smapt62_lon(thisx) - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  lnfmind(1) = thisx
                  lnfmind(2) = thisy
                end if
              end do
            end do
            if (masterproc) then
              atm2lnd_vars%lnfm(g,:) = atm2lnd_vars%lnfm_all(lnfmind(1),lnfmind(2),:)
              do np = 1,npes-1
                if (i == 1) then 
                  call mpi_recv(thisng,  1, MPI_INTEGER, np, 100000+np, mpicom, status, ier)
                  ng_all(np) = thisng
                end if
                if (i <= ng_all(np)) then 
                  call mpi_recv(lnfmind, 2, MPI_INTEGER, np, 200000+np, mpicom, status, ier)
                  call mpi_send(atm2lnd_vars%lnfm_all(lnfmind(1),lnfmind(2),:), 2920, &
                            MPI_REAL8, np, 300000+np, mpicom, ier)
                end if
              end do
            else
              if (i == 1)  call mpi_send(thisng,  1, MPI_INTEGER, 0, 100000+iam, mpicom, ier)
              call mpi_send(lnfmind, 2, MPI_INTEGER, 0, 200000+iam, mpicom, ier) 
              call mpi_recv(atm2lnd_vars%lnfm(g,:), 2920, MPI_REAL8, 0, 300000+iam, mpicom, status, ier)
            end if
          end if

          !Lightning data is 3-hourly.  Does not currently interpolate.
          atm2lnd_vars%forc_lnfm(g) = atm2lnd_vars%lnfm(g, ((int(thiscalday)-1)*8+tod/(3600*3))+1)

   !------------------------------------Nitrogen deposition----------------------------------------------

          !DMR note - ndep will NOT be correct if more than 1850 years of model
          !spinup (model year > 1850)
          nindex(1) = min(max(yr-1848,2), 168)
          nindex(2) = min(nindex(1)+1, 168)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then 
            if (masterproc .and. i .eq. 1) then 
              nu_nml = getavu()
              open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'ndepdyn_nml', status=nml_error)
              if (nml_error == 0) then
                read(nu_nml, nml=ndepdyn_nml,iostat=nml_error)
                if (nml_error /= 0) then
                  call endrun(msg='ERROR reading ndep namelist')
                end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_ndep), nf90_nowrite, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)      
              ierr = nf90_get_var(ncid, varid, smap2_lon)
              ierr = nf90_inq_varid(ncid, 'NDEP_year', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1)   = 144
              counti(2)   = 96
              counti(3)   = 1
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep1, starti, counti)
              if (nindex(1) .ne. nindex(2)) then 
                starti(3) = nindex(2)
                ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep2, starti, counti)
              else
                atm2lnd_vars%ndep2 = atm2lnd_vars%ndep1
              end if
              ierr = nf90_close(ncid)
             end if
             if (i .eq. 1) then
               call mpi_bcast (atm2lnd_vars%ndep1, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (atm2lnd_vars%ndep2, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lon, 144, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lat, 96, MPI_REAL8, 0, mpicom, ier)
             end if
          end if

          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then 
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then 
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then 
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
          end if

          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
  
          atm2lnd_vars%forc_ndep_grc(g)    = (atm2lnd_vars%ndep1(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt1(1) + &
                                              atm2lnd_vars%ndep2(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt2(1)) / (365._r8 * 86400._r8)
        end if model_filter

   !------------------------------------Aerosol forcing--------------------------------------------------
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then 
          if (masterproc .and. i .eq. 1) then 
            aerovars(1) = 'BCDEPWET'
            aerovars(2) = 'BCPHODRY'
            aerovars(3) = 'BCPHIDRY'
            aerovars(4) = 'OCDEPWET'
            aerovars(5) = 'OCPHODRY'
            aerovars(6) = 'OCPHIDRY'
            aerovars(7) = 'DSTX01DD'
            aerovars(8) = 'DSTX02DD'
            aerovars(9) = 'DSTX03DD'
            aerovars(10) = 'DSTX04DD'
            aerovars(11) = 'DSTX01WD'
            aerovars(12) = 'DSTX02WD'
            aerovars(13) = 'DSTX03WD'
            aerovars(14) = 'DSTX04WD'
            ierr = nf90_open(trim(aero_file), nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)      
            ierr = nf90_get_var(ncid, varid, smap2_lon)
            starti(1:2) = 1
            starti(3)   = max((min(yr,2100)-1849)*12+1, 13)-1
            counti(1)   = 144
            counti(2)   = 96
            counti(3)   = 14
            do av=1,14
              ierr = nf90_inq_varid(ncid, trim(aerovars(av)), varid)
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%aerodata(av,:,:,:), starti, counti)
            end do
            ierr = nf90_close(ncid)
          end if
          if (i .eq. 1) then 
             call mpi_bcast (atm2lnd_vars%aerodata, 14*144*96*14, MPI_REAL8, 0, mpicom, ier)
          end if
        end if

        !Use ndep grid indices since they're on the same grid
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. (.not. (use_fates .or. use_cn) )   ) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
        end if

        !get weights for interpolation (note this method doesn't get the month boundaries quite right..)
        aindex(1) = mon+1
        if (thiscalday .le. (caldaym(mon+1)+caldaym(mon))/2._r8) then 
           wt1(1) = 0.5_r8 + (thiscalday-caldaym(mon))/(caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)-1
        else
           wt1(1) = 1.0_r8 - (thiscalday-(caldaym(mon+1)+caldaym(mon))/2._r8)/   &
                          (caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)+1
        end if
        wt2(1) = 1._r8 - wt1(1)

        do av = 1,14
          atm2lnd_vars%forc_aer_grc(g,av)  =  atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(1))*wt1(1)+atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(2))*wt2(1)
        end do    

       !Parse startdate for adding temperature
       if (startdate_add_temperature .ne. '') then 
         call get_curr_date( yr, mon, day, tod )
         read(startdate_add_temperature,*) sdate_addt
         sy_addt     = sdate_addt/10000
         sm_addt     = (sdate_addt-sy_addt*10000)/100
         sd_addt     = sdate_addt-sy_addt*10000-sm_addt*100
         read(startdate_add_co2,*) sdate_addco2
         sy_addco2     = sdate_addco2/10000
         sm_addco2     = (sdate_addco2-sy_addco2*10000)/100
         sd_addco2     = sdate_addco2-sy_addco2*10000-sm_addt*100
       end if 
       if (startdate_add_temperature .ne. '') then
         if ((yr == sy_addt .and. mon == sm_addt .and. day >= sd_addt) .or. &
             (yr == sy_addt .and. mon > sm_addt) .or. (yr > sy_addt)) then
           atm2lnd_vars%forc_t_not_downscaled_grc(g) = atm2lnd_vars%forc_t_not_downscaled_grc(g) + add_temperature
           atm2lnd_vars%forc_th_not_downscaled_grc(g) = atm2lnd_vars%forc_th_not_downscaled_grc(g) + add_temperature
         end if
       end if

       !set the topounit-level atmospheric state and flux forcings (bypass mode)
       do topo = grc_pp%topi(g), grc_pp%topf(g)
         ! first, all the state forcings
         top_as%tbot(topo)    = atm2lnd_vars%forc_t_not_downscaled_grc(g)      ! forc_txy  Atm state K
         top_as%thbot(topo)   = atm2lnd_vars%forc_th_not_downscaled_grc(g)     ! forc_thxy Atm state K
         top_as%pbot(topo)    = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)   ! ptcmxy    Atm state Pa
         top_as%qbot(topo)    = atm2lnd_vars%forc_q_not_downscaled_grc(g)      ! forc_qxy  Atm state kg/kg
         top_as%ubot(topo)    = atm2lnd_vars%forc_u_grc(g)                     ! forc_uxy  Atm state m/s
         top_as%vbot(topo)    = atm2lnd_vars%forc_v_grc(g)                     ! forc_vxy  Atm state m/s
         if (implicit_stress) then
            top_as%wsresp(topo)  = 0._r8                                       !           Atm state m/s/Pa
            top_as%tau_est(topo) = 0._r8                                       !           Atm state Pa
         end if
         top_as%ugust(topo) = 0._r8                                            !           Atm state m/s
         top_as%zbot(topo)    = atm2lnd_vars%forc_hgt_grc(g)                   ! zgcmxy    Atm state m
         ! assign the state forcing fields derived from other inputs
         ! Horizontal windspeed (m/s)
         top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
         if (atm_gustiness) then
            top_as%windbot(topo) = top_as%windbot(topo) + top_as%ugust(topo)
         end if
         ! Relative humidity (percent)
         if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
            e = esatw(tdc(top_as%tbot(topo)))
         else
            e = esati(tdc(top_as%tbot(topo)))
         end if
         qsat           = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
         top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
         ! partial pressure of oxygen (Pa)
         top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
         ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
         vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
         top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))

         ! second, all the flux forcings
         top_af%rain(topo)    = forc_rainc + forc_rainl            ! sum of convective and large-scale rain
         top_af%snow(topo)    = forc_snowc + forc_snowl            ! sum of convective and large-scale snow
         top_af%solad(topo,2) = atm2lnd_vars%forc_solad_grc(g,2)   ! forc_sollxy  Atm flux  W/m^2
         top_af%solad(topo,1) = atm2lnd_vars%forc_solad_grc(g,1)   ! forc_solsxy  Atm flux  W/m^2
         top_af%solai(topo,2) = atm2lnd_vars%forc_solai_grc(g,2)   ! forc_solldxy Atm flux  W/m^2
         top_af%solai(topo,1) = atm2lnd_vars%forc_solai_grc(g,1)   ! forc_solsdxy Atm flux  W/m^2
         top_af%lwrad(topo)   = atm2lnd_vars%forc_lwrad_not_downscaled_grc(g)     ! flwdsxy Atm flux  W/m^2
         ! derived flux forcings
         top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
       end do
     
  !-----------------------------------------------------------------------------------------------------
#else

       atm2lnd_vars%forc_hgt_grc(g)     = x2l(index_x2l_Sa_z,i)         ! zgcmxy  Atm state m
       atm2lnd_vars%forc_u_grc(g)       = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
       atm2lnd_vars%forc_v_grc(g)       = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
       atm2lnd_vars%forc_solad_grc(g,2) = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solad_grc(g,1) = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,2) = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
       atm2lnd_vars%forc_solai_grc(g,1) = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2

       atm2lnd_vars%forc_th_not_downscaled_grc(g)    = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
       atm2lnd_vars%forc_q_not_downscaled_grc(g)     = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
       atm2lnd_vars%forc_pbot_not_downscaled_grc(g)  = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy  Atm state Pa
       atm2lnd_vars%forc_t_not_downscaled_grc(g)     = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
       atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2

       forc_rainc                                    = x2l(index_x2l_Faxa_rainc,i)   ! mm/s
       forc_rainl                                    = x2l(index_x2l_Faxa_rainl,i)   ! mm/s
       forc_snowc                                    = x2l(index_x2l_Faxa_snowc,i)   ! mm/s
       forc_snowl                                    = x2l(index_x2l_Faxa_snowl,i)   ! mm/s

       ! atmosphere coupling, for prognostic/prescribed aerosols
       atm2lnd_vars%forc_aer_grc(g,1)  =  x2l(index_x2l_Faxa_bcphidry,i)
       atm2lnd_vars%forc_aer_grc(g,2)  =  x2l(index_x2l_Faxa_bcphodry,i)
       atm2lnd_vars%forc_aer_grc(g,3)  =  x2l(index_x2l_Faxa_bcphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,4)  =  x2l(index_x2l_Faxa_ocphidry,i)
       atm2lnd_vars%forc_aer_grc(g,5)  =  x2l(index_x2l_Faxa_ocphodry,i)
       atm2lnd_vars%forc_aer_grc(g,6)  =  x2l(index_x2l_Faxa_ocphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,7)  =  x2l(index_x2l_Faxa_dstwet1,i)
       atm2lnd_vars%forc_aer_grc(g,8)  =  x2l(index_x2l_Faxa_dstdry1,i)
       atm2lnd_vars%forc_aer_grc(g,9)  =  x2l(index_x2l_Faxa_dstwet2,i)
       atm2lnd_vars%forc_aer_grc(g,10) =  x2l(index_x2l_Faxa_dstdry2,i)
       atm2lnd_vars%forc_aer_grc(g,11) =  x2l(index_x2l_Faxa_dstwet3,i)
       atm2lnd_vars%forc_aer_grc(g,12) =  x2l(index_x2l_Faxa_dstdry3,i)
       atm2lnd_vars%forc_aer_grc(g,13) =  x2l(index_x2l_Faxa_dstwet4,i)
       atm2lnd_vars%forc_aer_grc(g,14) =  x2l(index_x2l_Faxa_dstdry4,i)
       
       !set the topounit-level atmospheric state and flux forcings
       if (use_atm_downscaling_to_topunit) then
          if(atm_gustiness) then
             call endrun("Error: atm_gustiness not yet supported with multiple topounits")
          end if
         do topo = grc_pp%topi(g) , grc_pp%topf(g)
            top_as%ugust(topo) = 0._r8
         end do

         call downscale_atm_forcing_to_topounit(g, i, x2l, lnd2atm_vars)
       else
         do topo = grc_pp%topi(g), grc_pp%topf(g)
           ! first, all the state forcings
           top_as%tbot(topo)    = x2l(index_x2l_Sa_tbot,i)      ! forc_txy  Atm state K
           top_as%thbot(topo)   = x2l(index_x2l_Sa_ptem,i)      ! forc_thxy Atm state K
           top_as%pbot(topo)    = x2l(index_x2l_Sa_pbot,i)      ! ptcmxy    Atm state Pa
           top_as%qbot(topo)    = x2l(index_x2l_Sa_shum,i)      ! forc_qxy  Atm state kg/kg
           top_as%ubot(topo)    = x2l(index_x2l_Sa_u,i)         ! forc_uxy  Atm state m/s
           top_as%vbot(topo)    = x2l(index_x2l_Sa_v,i)         ! forc_vxy  Atm state m/s
           top_as%zbot(topo)    = x2l(index_x2l_Sa_z,i)         ! zgcmxy    Atm state m
         if (implicit_stress) then
            top_as%wsresp(topo)  = x2l(index_x2l_Sa_wsresp,i) !           Atm state m/s/Pa
            top_as%tau_est(topo) = x2l(index_x2l_Sa_tau_est,i)!           Atm state Pa
         end if
         if (atm_gustiness) then
            top_as%ugust(topo)  = x2l(index_x2l_Sa_ugust,i)   !           Atm state m/s
         else
            top_as%ugust(topo) = 0._r8
         end if
           ! assign the state forcing fields derived from other inputs
           ! Horizontal windspeed (m/s)
           top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
         if (atm_gustiness) then
            top_as%windbot(topo) = top_as%windbot(topo) + top_as%ugust(topo)
         end if
           ! Relative humidity (percent)
           if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
            e = esatw(tdc(top_as%tbot(topo)))
           else
            e = esati(tdc(top_as%tbot(topo)))
           end if
           qsat = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
           top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
           ! partial pressure of oxygen (Pa)
           top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
           ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
           vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
           top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))
         
           ! second, all the flux forcings
           top_af%rain(topo)    = forc_rainc + forc_rainl       ! sum of convective and large-scale rain
           top_af%snow(topo)    = forc_snowc + forc_snowl       ! sum of convective and large-scale snow
           top_af%solad(topo,2) = x2l(index_x2l_Faxa_swndr,i)   ! forc_sollxy  Atm flux  W/m^2
           top_af%solad(topo,1) = x2l(index_x2l_Faxa_swvdr,i)   ! forc_solsxy  Atm flux  W/m^2
           top_af%solai(topo,2) = x2l(index_x2l_Faxa_swndf,i)   ! forc_solldxy Atm flux  W/m^2
           top_af%solai(topo,1) = x2l(index_x2l_Faxa_swvdf,i)   ! forc_solsdxy Atm flux  W/m^2
           top_af%lwrad(topo)   = x2l(index_x2l_Faxa_lwdn,i)    ! flwdsxy Atm flux  W/m^2
           ! derived flux forcings
           top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
         end do
       end if  
#endif

       ! Determine optional receive fields
       ! CO2 (and C13O2) concentration: constant, prognostic, or diagnostic
       if (co2_type_idx == 0) then                    ! CO2 constant, value from namelist
         co2_ppmv_val = co2_ppmv
       else if (co2_type_idx == 1) then               ! CO2 prognostic, value from coupler field
         co2_ppmv_val = x2l(index_x2l_Sa_co2prog,i)
       else if (co2_type_idx == 2) then               ! CO2 diagnostic, value from coupler field
         co2_ppmv_val = x2l(index_x2l_Sa_co2diag,i)
       else
         call endrun( sub//' ERROR: Invalid co2_type_idx, must be 0, 1, or 2 (constant, prognostic, or diagnostic)' )
       end if
       ! Assign to topounits, with conversion from ppmv to partial pressure (Pa)
       ! If using C13, then get the c13ratio from elm_varcon (constant value for pre-industrial atmosphere)

       do topo = grc_pp%topi(g), grc_pp%topf(g)
         top_as%pco2bot(topo) = co2_ppmv_val * 1.e-6_r8 * top_as%pbot(topo)
         if (use_c13) then
            top_as%pc13o2bot(topo) = top_as%pco2bot(topo) * c13ratio;
         end if
       end do
       ! CH4
       if (index_x2l_Sa_methane /= 0) then
          do topo = grc_pp%topi(g), grc_pp%topf(g)
            top_as%pch4bot(topo) = x2l(index_x2l_Sa_methane,i)
          end do
       endif

       if (index_x2l_Sa_co2prog /= 0) then
          co2_ppmv_prog = x2l(index_x2l_Sa_co2prog,i)   ! co2 atm state prognostic
       else
          co2_ppmv_prog = co2_ppmv
       end if

       if (index_x2l_Sa_co2diag /= 0) then
          co2_ppmv_diag = x2l(index_x2l_Sa_co2diag,i)   ! co2 atm state diagnostic
       else
          co2_ppmv_diag = co2_ppmv
       end if

       if (index_x2l_Sa_methane /= 0) then
          atm2lnd_vars%forc_pch4_grc(g) = x2l(index_x2l_Sa_methane,i)
       endif

       ! Determine derived quantities for required fields

       forc_t = atm2lnd_vars%forc_t_not_downscaled_grc(g)
       forc_q = atm2lnd_vars%forc_q_not_downscaled_grc(g)
       forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)
       
       atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of wind [m]
       atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of temperature [m]
       atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)    !observational height of humidity [m]
       atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot  / (0.622_r8 + 0.378_r8 * forc_q)
       atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
            (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
       atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
       atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
       atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                        atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)
       
       atm2lnd_vars%forc_rain_not_downscaled_grc(g)  = forc_rainc + forc_rainl
       atm2lnd_vars%forc_snow_not_downscaled_grc(g)  = forc_snowc + forc_snowl
       if (forc_t > SHR_CONST_TKFRZ) then
          e = esatw(tdc(forc_t))
       else
          e = esati(tdc(forc_t))
       end if
       qsat           = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
       atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)
       ! Make sure relative humidity is properly bounded
       ! atm2lnd_vars%forc_rh_grc(g) = min( 100.0_r8, atm2lnd_vars%forc_rh_grc(g) )
       ! atm2lnd_vars%forc_rh_grc(g) = max(   0.0_r8, atm2lnd_vars%forc_rh_grc(g) )

       ! Determine derived quantities for optional fields
       ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
       ! Note that forc_pbot is in Pa

#ifdef CPL_BYPASS
       co2_type_idx = 2
#endif

       if (co2_type_idx == 1) then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_type_idx == 2) then
#ifdef CPL_BYPASS
        !atmospheric CO2 (to be used for transient simulations only)
        if (atm2lnd_vars%loaded_bypassdata .eq. 0) then 
          ierr = nf90_open(trim(co2_file), nf90_nowrite, ncid)
          ierr = nf90_inq_dimid(ncid, 'time', dimid)
          ierr = nf90_Inquire_Dimension(ncid, dimid, len = thistimelen)
          ierr = nf90_inq_varid(ncid, 'CO2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%co2_input(:,:,1:thistimelen))
          ierr = nf90_inq_varid(ncid, 'C13O2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%c13o2_input(:,:,1:thistimelen))
          ierr = nf90_close(ncid)
        end if

        !get weights/indices for interpolation (assume values represent annual averages)
        nindex(1) = min(max(yr,1850),2100)-1764
        if (thiscalday .le. 182.5) then 
          nindex(2) = nindex(1)-1  
        else
          nindex(2) = nindex(1)+1
        end if
        wt1(1) = 1._r8 - abs((182.5 - (thiscalday -1._r8))/365._r8)
        wt2(1) = 1._r8 - wt1(1)

        co2_ppmv_val = atm2lnd_vars%co2_input(1,1,nindex(1))*wt1(1) + atm2lnd_vars%co2_input(1,1,nindex(2))*wt2(1)
        if (startdate_add_co2 .ne. '') then
          if ((yr == sy_addco2 .and. mon == sm_addco2 .and. day >= sd_addco2) .or. &
              (yr == sy_addco2 .and. mon > sm_addco2) .or. (yr > sy_addco2)) then
            co2_ppmv_val=co2_ppmv_val + add_co2
          end if
        end if

        if (use_c13) then 
          atm2lnd_vars%forc_pc13o2_grc(g) = (atm2lnd_vars%c13o2_input(1,1,nindex(1))*wt1(1) + &
               atm2lnd_vars%c13o2_input(1,1,nindex(2))*wt2(1)) * 1.e-6_r8 * forc_pbot
        end if
        co2_type_idx = 1
#else
          co2_ppmv_val = co2_ppmv_diag 
           if (use_c13) then
             atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
           end if
#endif
       else
          co2_ppmv_val = co2_ppmv
          if (use_c13) then
            atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
          end if
       end if
       atm2lnd_vars%forc_pco2_grc(g)   = co2_ppmv_val * 1.e-6_r8 * forc_pbot 

#ifdef CPL_BYPASS
       do topo = grc_pp%topi(g), grc_pp%topf(g)
         top_as%pco2bot(topo) = atm2lnd_vars%forc_pco2_grc(g)
         if (use_c13) then
            top_as%pc13o2bot(topo) = atm2lnd_vars%forc_pc13o2_grc(g)
         end if
       end do
#endif
      
       ! glc coupling 

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             glc2lnd_vars%frac_grc(g,num)  = x2l(index_x2l_Sg_frac(num),i)
             glc2lnd_vars%topo_grc(g,num)  = x2l(index_x2l_Sg_topo(num),i)
             glc2lnd_vars%hflx_grc(g,num)  = x2l(index_x2l_Flgg_hflx(num),i)
          end do
          glc2lnd_vars%icemask_grc(g)  = x2l(index_x2l_Sg_icemask,i)
          glc2lnd_vars%icemask_coupled_fluxes_grc(g)  = x2l(index_x2l_Sg_icemask_coupled_fluxes,i)
       end if

    end do     
#ifdef CPL_BYPASS
    atm2lnd_vars%loaded_bypassdata = 1
#endif

  end subroutine lnd_import

  !===============================================================================

  subroutine lnd_export( bounds, lnd2atm_vars, lnd2glc_vars, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the elm model to the coupler 
    ! 
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use elm_varctl         , only : iulog, create_glacier_mec_landunit
    use clm_time_manager   , only : get_nstep, get_step_size  
    use domainMod          , only : ldomain
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_vars ! clm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_vars ! clm land to atmosphere exchange data type
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i   ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step   
    integer  :: num   ! counter
    character(len=*), parameter :: sub = 'lnd_export_mct'
    !---------------------------------------------------------------------------

    dtime = get_step_size()

    ! cesm sign convention is that fluxes are positive downward

    l2x(:,:) = 0.0_r8

    do g = bounds%begg,bounds%endg
       i = 1 + (g-bounds%begg)
       l2x(index_l2x_Sl_t,i)        =  lnd2atm_vars%t_rad_grc(g)
       l2x(index_l2x_Sl_snowh,i)    =  lnd2atm_vars%h2osno_grc(g)
       l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_vars%albd_grc(g,1)
       l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_vars%albd_grc(g,2)
       l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_vars%albi_grc(g,1)
       l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_vars%albi_grc(g,2)
       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_vars%t_ref2m_grc(g)
       l2x(index_l2x_Sl_qref,i)     =  lnd2atm_vars%q_ref2m_grc(g)
       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_vars%u_ref10m_grc(g)
       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_vars%taux_grc(g)
       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_vars%tauy_grc(g)
       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_vars%eflx_lh_tot_grc(g)
       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_vars%eflx_sh_tot_grc(g)
       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_vars%eflx_lwrad_out_grc(g)
       l2x(index_l2x_Fall_evap,i)   = -lnd2atm_vars%qflx_evap_tot_grc(g)
       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_vars%fsa_grc(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x(index_l2x_Fall_fco2_lnd,i) = -lnd2atm_vars%nee_grc(g)  
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i)     =  lnd2atm_vars%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)       =  lnd2atm_vars%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)    =  lnd2atm_vars%h2osoi_vol_grc(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -lnd2atm_vars%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -lnd2atm_vars%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -lnd2atm_vars%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -lnd2atm_vars%flxdst_grc(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               lnd2atm_vars%ddvel_grc(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -lnd2atm_vars%flxvoc_grc(g,:shr_megan_mechcomps_n)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x(index_l2x_Fall_methane,i) = -lnd2atm_vars%flux_ch4_grc(g) 
       endif

       ! sign convention is positive downward with 
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x(index_l2x_Flrl_rofi,i) = lnd2atm_vars%qflx_rofice_grc(g)
       l2x(index_l2x_Flrl_rofsur,i) = lnd2atm_vars%qflx_rofliq_qsur_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsurp_grc(g)   !  surface ponding
       l2x(index_l2x_Flrl_rofsub,i) = lnd2atm_vars%qflx_rofliq_qsub_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsubp_grc(g)   !  perched drainiage
       l2x(index_l2x_Flrl_rofgwl,i) = lnd2atm_vars%qflx_rofliq_qgwl_grc(g)
  
       l2x(index_l2x_Flrl_demand,i) =  lnd2atm_vars%qflx_irr_demand_grc(g)   ! needs to be filled in
       if (l2x(index_l2x_Flrl_demand,i) > 0.0_r8) then
           write(iulog,*)'lnd2atm_vars%qflx_irr_demand_grc is',lnd2atm_vars%qflx_irr_demand_grc(g)
           write(iulog,*)'l2x(index_l2x_Flrl_demand,i) is',l2x(index_l2x_Flrl_demand,i)
           call endrun( sub//' ERROR: demand must be <= 0.')
       endif
       l2x(index_l2x_Flrl_Tqsur,i)  = lnd2atm_vars%Tqsur_grc(g)
       l2x(index_l2x_Flrl_Tqsub,i)  = lnd2atm_vars%Tqsub_grc(g)
       l2x(index_l2x_coszen_str,i)  = lnd2atm_vars%coszen_str(g)
       l2x(index_l2x_Flrl_rofmud,i) = lnd2atm_vars%qflx_rofmud_grc(g)
       l2x(index_l2x_Flrl_wslake,i) = lnd2atm_vars%wslake_grc(g)/dtime

       if (index_l2x_Flrl_inundinf /= 0) then
          l2x(index_l2x_Flrl_inundinf,i) = lnd2atm_vars%qflx_h2orof_drain_grc(g)
       endif
       
       ! glc coupling

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             l2x(index_l2x_Sl_tsrf(num),i)   = lnd2glc_vars%tsrf_grc(g,num)
             l2x(index_l2x_Sl_topo(num),i)   = lnd2glc_vars%topo_grc(g,num)
             l2x(index_l2x_Flgl_qice(num),i) = lnd2glc_vars%qice_grc(g,num)
          end do
       end if

    end do

  end subroutine lnd_export

end module lnd_import_export



double precision function szenith(xcoor, ycoor, ltm, jday, hr, min, offset)     
  !Function to calcualte solar zenith angle
  !Used in coupler bypass mode to compute inerpolation for incoming solar

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  implicit none
  !inputs
  real(r8) xcoor, ycoor, offset_min
  integer jday, hr, min, ltm, offset
  !working variables
  real(r8) d2r, r2d, lsn, latrad, decrad, decdeg, ha
  real(r8) hangle, harad, saltrad, saltdeg, sazirad, sazideg
  real(r8) szendeg,szenrad
  
  real pi
  parameter(pi = 3.14159265358979)
  offset_min = offset/60d0   !note assumes 1hr or smaller timestep
  min = min - offset_min  
   
  !adjust time for offsets
  if (min < 0) then
    hr = hr - 1
    min = min+60
  end if
  if (min >= 60) then 
    hr = hr+1
    min = min-60
  end if
  if (hr < 0) then  
    hr = hr+24
    jday = jday-1
  end if
  if (hr >= 24) then
    hr = hr-24
    jday = jday+1
  end if
    
  if (jday < 1) jday = 1
  if (xcoor > 180d0) xcoor = xcoor-360d0

  d2r     = pi/180d0
  r2d     = 1/d2r
  lsn     = 12.0d0+((ltm-xcoor)/15.0d0)
  latrad  = ycoor*d2r
  decrad  = 23.45*d2r*sin(d2r*360d0*(284d0+jday)/365d0)
  decdeg  = decrad*r2d
  ha      = hr+min/60.0d0 
  hangle  = (lsn-ha)*60.0d0               
  harad   = hangle*0.0043633d0       
  
  saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad) &
       *cos(harad)))
  saltdeg = saltrad * r2d
  sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
  sazideg = sazirad * r2d
  
  IF (saltdeg.LT.0.0d0 .OR. saltrad.GT.180.0d0) THEN  ! sun is below horizon
     saltdeg = 0.0d0
     saltrad = 0.0d0
     szendeg = 90.0d0
     szenrad = 90.0d0*d2r
     !mass    = 1229d0**.5d0             ! if solaralt=0 -> sin(0)=0
  ELSE
     szendeg = 90d0-saltdeg
     szenrad = szendeg*d2r
  ENDIF
  szenith = szendeg
  
end function szenith
