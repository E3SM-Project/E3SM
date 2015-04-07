module SoilStateInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Set hydraulic and thermal properties 
  !
  ! !USES
  use SoilStateType , only : soilstate_type
  use LandunitType  , only : lun                
  use ColumnType    , only : col                
  use PatchType     , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: SoilStateInitTimeConst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ReadNL
  !
  ! !PRIVATE DATA:
  ! Control variables (from namelist)
  logical, private :: organic_frac_squared ! If organic fraction should be squared (as in CLM4.5)
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine ReadNL( nlfilename )
    !
    ! !DESCRIPTION:
    ! Read namelist for SoilStateType
    !
    ! !USES:
    use shr_mpi_mod    , only : shr_mpi_bcast
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use fileutils      , only : getavu, relavu, opnfil
    use clm_nlUtilsMod , only : find_nlgroup_name
    use clm_varctl     , only : iulog
    use spmdMod        , only : mpicom, masterproc
    use abortUtils     , only : endrun    
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: nlfilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'SoilState_readnl'  ! subroutine name
    !-----------------------------------------------------------------------

    character(len=*), parameter :: nl_name  = 'clm_soilstate_inparm'  ! Namelist name
                                                                      ! MUST agree with name in namelist and read
    namelist / clm_soilstate_inparm / organic_frac_squared

    ! preset values

    organic_frac_squared = .false.

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in '//nl_name//' namelist'
       call opnfil (nlfilename, unitn, 'F')
       call find_nlgroup_name(unitn, nl_name, status=ierr)
       if (ierr == 0) then
          read(unit=unitn, nml=clm_soilstate_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading '//nl_name//' namelist"//errmsg(__FILE__, __LINE__))
          end if
       end if
       call relavu( unitn )

    end if

    call shr_mpi_bcast(organic_frac_squared, mpicom)

  end subroutine ReadNL

  !-----------------------------------------------------------------------
  subroutine SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) 
    !
    ! !USES:
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use shr_log_mod         , only : errMsg => shr_log_errMsg
    use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
    use decompMod           , only : bounds_type
    use abortutils          , only : endrun
    use spmdMod             , only : masterproc
    use ncdio_pio           , only : file_desc_t, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use ncdio_pio           , only : ncd_pio_openfile, ncd_pio_closefile, ncd_inqdlen
    use clm_varpar          , only : more_vertlayers, numpft, numrad 
    use clm_varpar          , only : nlevsoi, nlevgrnd, nlevlak, nlevsoifl, nlayer, nlayert, nlevurb, nlevsno
    use clm_varcon          , only : zsoi, dzsoi, zisoi, spval
    use clm_varcon          , only : secspday, pc, mu, denh2o, denice, grlnd
    use clm_varctl          , only : use_cn, use_lch4, use_ed
    use clm_varctl          , only : iulog, fsurdat, paramfile 
    use landunit_varcon     , only : istice, istdlak, istwet, istsoil, istcrop, istice_mec
    use column_varcon       , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv 
    use fileutils           , only : getfil
    use organicFileMod      , only : organicrd 
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use RootBiophysMod      , only : init_vegrootfr
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds  
    type(soilstate_type) , intent(inout) :: soilstate_inst
    character(len=*)     , intent(in)    :: nlfilename ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer            :: p, lev, c, l, g, j            ! indices
    real(r8)           :: om_frac                       ! organic matter fraction
    real(r8)           :: om_tkm         = 0.25_r8      ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8)           :: om_watsat_lake = 0.9_r8       ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     ! porosity of organic soil
    real(r8)           :: om_hksat                      ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_csol        = 2.5_r8       ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8)           :: om_tkd         = 0.05_r8      ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b                          ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8       ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: csol_bedrock   = 2.0e6_r8     ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: pc_lake        = 0.5_r8       ! percolation threshold
    real(r8)           :: perc_frac                     ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm                     ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat                   ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac                    ! fraction of "unconnected" soil
    real(r8)           :: bd                            ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm                           ! mineral conductivity
    real(r8)           :: xksat                         ! maximum hydraulic conductivity of soil [mm/s]
    real(r8)           :: clay,sand                     ! temporaries
    real(r8)           :: organic_max                   ! organic matter (kg/m3) where soil is assumed to act like peat
    integer            :: dimid                         ! dimension id
    logical            :: readvar 
    type(file_desc_t)  :: ncid                          ! netcdf id
    real(r8) ,pointer  :: zsoifl (:)                    ! Output: [real(r8) (:)]  original soil midpoint 
    real(r8) ,pointer  :: zisoifl (:)                   ! Output: [real(r8) (:)]  original soil interface depth 
    real(r8) ,pointer  :: dzsoifl (:)                   ! Output: [real(r8) (:)]  original soil thickness 
    real(r8) ,pointer  :: gti (:)                       ! read in - fmax 
    real(r8) ,pointer  :: sand3d (:,:)                  ! read in - soil texture: percent sand (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: clay3d (:,:)                  ! read in - soil texture: percent clay (needs to be a pointer for use in ncdio)
    real(r8) ,pointer  :: organic3d (:,:)               ! read in - organic matter: kg/m3 (needs to be a pointer for use in ncdio)
    character(len=256) :: locfn                         ! local filename
    integer            :: ipedof  
    integer            :: begp, endp
    integer            :: begc, endc
    integer            :: begg, endg
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    do c = begc,endc
       soilstate_inst%smpmin_col(c) = -1.e8_r8
    end do

    ! --------------------------------------------------------------------
    ! Read namelist
    ! --------------------------------------------------------------------

    call ReadNL( nlfilename )

    ! --------------------------------------------------------------------
    ! Initialize root fraction (computing from surface, d is depth in meter):
    ! --------------------------------------------------------------------

    ! Currently pervious road has same properties as soil
    do c = begc,endc
       l = col%landunit(c)

       if (lun%urbpoi(l) .and. col%itype(c) == icol_road_perv) then 
          do lev = 1, nlevgrnd
             soilstate_inst%rootfr_road_perv_col(c,lev) = 0._r8
          enddo
          do lev = 1,nlevsoi
             soilstate_inst%rootfr_road_perv_col(c,lev) = 0.1_r8  ! uniform profile
          end do
       end if
    end do

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          soilstate_inst%rootfr_col (c,nlevsoi+1:nlevgrnd) = 0._r8
       else
          ! Inactive CH4 columns 
          ! (Also includes (lun%itype(l)==istdlak .and.  allowlakeprod), which used to be
          ! in a separate branch of the conditional)
          soilstate_inst%rootfr_col (c,:) = spval
       end if
    end do

    ! Initialize root fraction 
    ! Note that ED has its own root fraction root fraction routine and should not
    ! use the following since it depends on patch%itype - which ED should not use

    if (.not. use_ed) then
        call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%rootfr_patch(begp:endp,1:nlevgrnd))
     end if

    ! --------------------------------------------------------------------
    ! dynamic memory allocation
    ! --------------------------------------------------------------------

    allocate(sand3d(begg:endg,nlevsoifl))
    allocate(clay3d(begg:endg,nlevsoifl))

    ! Determine organic_max from parameter file

    call getfil (paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call ncd_io(ncid=ncid, varname='organic_max', flag='read', data=organic_max, readvar=readvar)
    if ( .not. readvar ) call endrun(msg=' ERROR: organic_max not on param file'//errMsg(__FILE__, __LINE__))
    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! Read surface dataset
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    end if

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    call ncd_inqdlen(ncid,dimid,nlevsoifl,name='nlevsoi')
    if ( .not. more_vertlayers )then
       if ( nlevsoifl /= nlevsoi )then
          call endrun(msg=' ERROR: Number of soil layers on file does NOT match the number being used'//&
               errMsg(__FILE__, __LINE__))
       end if
    else
       ! read in layers, interpolate to high resolution grid later
    end if

    ! Read in organic matter dataset 

    allocate(organic3d(begg:endg,nlevsoifl))
    call organicrd(organic3d)

    ! Read in sand and clay data

    call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_SAND NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if

    call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: PCT_CLAY NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if

    do p = begp,endp
       g = patch%gridcell(p)
       if ( sand3d(g,1)+clay3d(g,1) == 0.0_r8 )then
          if ( any( sand3d(g,:)+clay3d(g,:) /= 0.0_r8 ) )then
             call endrun(msg='found depth points that do NOT sum to zero when surface does'//&
                  errMsg(__FILE__, __LINE__)) 
          end if
          sand3d(g,:) = 1.0_r8
          clay3d(g,:) = 1.0_r8
       end if
       if ( any( sand3d(g,:)+clay3d(g,:) == 0.0_r8 ) )then
          call endrun(msg='after setting, found points sum to zero'//errMsg(__FILE__, __LINE__)) 
       end if

       soilstate_inst%sandfrac_patch(p) = sand3d(g,1)/100.0_r8
       soilstate_inst%clayfrac_patch(p) = clay3d(g,1)/100.0_r8
    end do

    ! Read fmax

    allocate(gti(begg:endg))
    call ncd_io(ncid=ncid, varname='FMAX', flag='read', data=gti, dim1name=grlnd, readvar=readvar)
    if (.not. readvar) then
       call endrun(msg=' ERROR: FMAX NOT on surfdata file'//errMsg(__FILE__, __LINE__)) 
    end if
    do c = begc, endc
       g = col%gridcell(c)
       soilstate_inst%wtfact_col(c) = gti(g)
    end do
    deallocate(gti)

    ! Close file

    call ncd_pio_closefile(ncid)

    ! --------------------------------------------------------------------
    ! get original soil depths to be used in interpolation of sand and clay
    ! --------------------------------------------------------------------

    allocate(zsoifl(1:nlevsoifl), zisoifl(0:nlevsoifl), dzsoifl(1:nlevsoifl))
    do j = 1, nlevsoifl
       zsoifl(j) = 0.025*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
    enddo

    dzsoifl(1) = 0.5_r8*(zsoifl(1)+zsoifl(2))             !thickness b/n two interfaces
    do j = 2,nlevsoifl-1
       dzsoifl(j)= 0.5_r8*(zsoifl(j+1)-zsoifl(j-1))
    enddo
    dzsoifl(nlevsoifl) = zsoifl(nlevsoifl)-zsoifl(nlevsoifl-1)

    zisoifl(0) = 0._r8
    do j = 1, nlevsoifl-1
       zisoifl(j) = 0.5_r8*(zsoifl(j)+zsoifl(j+1))         !interface depths
    enddo
    zisoifl(nlevsoifl) = zsoifl(nlevsoifl) + 0.5_r8*dzsoifl(nlevsoifl)

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: non-lake
    ! --------------------------------------------------------------------

    !   urban roof, sunwall and shadewall thermal properties used to 
    !   derive thermal conductivity and heat capacity are set to special 
    !   value because thermal conductivity and heat capacity for urban 
    !   roof, sunwall and shadewall are prescribed in SoilThermProp.F90 
    !   in SoilPhysicsMod.F90


    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istwet .or. lun%itype(l)==istice .or. lun%itype(l)==istice_mec) then

          do lev = 1,nlevgrnd
             soilstate_inst%bsw_col(c,lev)    = spval
             soilstate_inst%watsat_col(c,lev) = spval
             soilstate_inst%watfc_col(c,lev)  = spval
             soilstate_inst%hksat_col(c,lev)  = spval
             soilstate_inst%sucsat_col(c,lev) = spval
             soilstate_inst%watdry_col(c,lev) = spval 
             soilstate_inst%watopt_col(c,lev) = spval 
             soilstate_inst%bd_col(c,lev)     = spval 
             if (lev <= nlevsoi) then
                soilstate_inst%cellsand_col(c,lev) = spval
                soilstate_inst%cellclay_col(c,lev) = spval
                soilstate_inst%cellorg_col(c,lev)  = spval
             end if
          end do

          do lev = 1,nlevgrnd
             soilstate_inst%tkmg_col(c,lev)   = spval
             soilstate_inst%tksatu_col(c,lev) = spval
             soilstate_inst%tkdry_col(c,lev)  = spval
             if (lun%itype(l)==istwet .and. lev > nlevsoi) then
                soilstate_inst%csol_col(c,lev) = csol_bedrock
             else
                soilstate_inst%csol_col(c,lev)= spval
             endif
          end do

       else if (lun%urbpoi(l) .and. (col%itype(c) /= icol_road_perv) .and. (col%itype(c) /= icol_road_imperv) )then

          ! Urban Roof, sunwall, shadewall properties set to special value
          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = spval
             soilstate_inst%watfc_col(c,lev)  = spval
             soilstate_inst%bsw_col(c,lev)    = spval
             soilstate_inst%hksat_col(c,lev)  = spval
             soilstate_inst%sucsat_col(c,lev) = spval
             soilstate_inst%watdry_col(c,lev) = spval 
             soilstate_inst%watopt_col(c,lev) = spval 
             soilstate_inst%bd_col(c,lev) = spval 
             if (lev <= nlevsoi) then
                soilstate_inst%cellsand_col(c,lev) = spval
                soilstate_inst%cellclay_col(c,lev) = spval
                soilstate_inst%cellorg_col(c,lev)  = spval
             end if
          end do

          do lev = 1,nlevgrnd
             soilstate_inst%tkmg_col(c,lev)   = spval
             soilstate_inst%tksatu_col(c,lev) = spval
             soilstate_inst%tkdry_col(c,lev)  = spval
             soilstate_inst%csol_col(c,lev)   = spval
          end do

       else

          do lev = 1,nlevgrnd

             if ( more_vertlayers )then ! duplicate clay and sand values from last soil layer

                if (lev .eq. 1) then
                   clay = clay3d(g,1)
                   sand = sand3d(g,1)
                   om_frac = organic3d(g,1)/organic_max 
                else if (lev <= nlevsoi) then
                   do j = 1,nlevsoifl-1
                      if (zisoi(lev) >= zisoifl(j) .AND. zisoi(lev) < zisoifl(j+1)) then
                         clay = clay3d(g,j+1)
                         sand = sand3d(g,j+1)
                         om_frac = organic3d(g,j+1)/organic_max    
                      endif
                   end do
                else
                   clay = clay3d(g,nlevsoifl)
                   sand = sand3d(g,nlevsoifl)
                   om_frac = 0._r8
                endif
             else
                if (lev <= nlevsoi) then ! duplicate clay and sand values from 10th soil layer
                   clay = clay3d(g,lev)
                   sand = sand3d(g,lev)
                   if ( organic_frac_squared )then
                      om_frac = (organic3d(g,lev)/organic_max)**2._r8
                   else
                      om_frac = organic3d(g,lev)/organic_max
                   end if
                else
                   clay = clay3d(g,nlevsoi)
                   sand = sand3d(g,nlevsoi)
                   om_frac = 0._r8
                endif
             end if

             if (lun%itype(l) == istdlak) then

                if (lev <= nlevsoi) then
                   soilstate_inst%cellsand_col(c,lev) = sand
                   soilstate_inst%cellclay_col(c,lev) = clay
                   soilstate_inst%cellorg_col(c,lev)  = om_frac*organic_max
                end if

             else if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types

                if (lun%urbpoi(l)) then
                   om_frac = 0._r8 ! No organic matter for urban
                end if

                if (lev <= nlevsoi) then
                   soilstate_inst%cellsand_col(c,lev) = sand
                   soilstate_inst%cellclay_col(c,lev) = clay
                   soilstate_inst%cellorg_col(c,lev)  = om_frac*organic_max
                end if

                ! Note that the following properties are overwritten for urban impervious road 
                ! layers that are not soil in SoilThermProp.F90 within SoilTemperatureMod.F90

                !determine the type of pedotransfer function to be used based on soil order
                !I will use the following implementation to further explore the ET problem, now
                !I set soil order to 0 for all soils. Jinyun Tang, Mar 20, 2014

                ipedof=get_ipedof(0)
                call pedotransf(ipedof, sand, clay, &
                     soilstate_inst%watsat_col(c,lev), soilstate_inst%bsw_col(c,lev), soilstate_inst%sucsat_col(c,lev), xksat)

                om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
                om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
                om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
                om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), 0.0001_r8)

                soilstate_inst%bd_col(c,lev)        = (1._r8 - soilstate_inst%watsat_col(c,lev))*2.7e3_r8 
                soilstate_inst%watsat_col(c,lev)    = (1._r8 - om_frac) * soilstate_inst%watsat_col(c,lev) + om_watsat*om_frac
                tkm                                 = (1._r8-om_frac) * (8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac ! W/(m K)
                soilstate_inst%bsw_col(c,lev)       = (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
                soilstate_inst%sucsat_col(c,lev)    = (1._r8-om_frac) * soilstate_inst%sucsat_col(c,lev) + om_sucsat*om_frac  
                soilstate_inst%hksat_min_col(c,lev) = xksat

                ! perc_frac is zero unless perf_frac greater than percolation threshold
                if (om_frac > pcalpha) then
                   perc_norm=(1._r8 - pcalpha)**(-pcbeta)
                   perc_frac=perc_norm*(om_frac - pcalpha)**pcbeta
                else
                   perc_frac=0._r8
                endif

                ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
                uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

                ! uncon_hksat is series addition of mineral/organic conductivites
                if (om_frac < 1._r8) then
                   uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                        +((1._r8-perc_frac)*om_frac)/om_hksat)
                else
                   uncon_hksat = 0._r8
                end if
                soilstate_inst%hksat_col(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

                soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8- soilstate_inst%watsat_col(c,lev))           

                soilstate_inst%tksatu_col(c,lev) = soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)

                soilstate_inst%tkdry_col(c,lev)  = ((0.135_r8*soilstate_inst%bd_col(c,lev) + 64.7_r8) / &
                     (2.7e3_r8 - 0.947_r8*soilstate_inst%bd_col(c,lev)))*(1._r8-om_frac) + om_tkd*om_frac  

                soilstate_inst%csol_col(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
                     om_csol*om_frac)*1.e6_r8  ! J/(m3 K)

                if (lev > nlevsoi) then
                   soilstate_inst%csol_col(c,lev) = csol_bedrock
                endif

                soilstate_inst%watdry_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (316230._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 
                soilstate_inst%watopt_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (158490._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

                !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
                ! water content at field capacity, defined as hk = 0.1 mm/day
                ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / secspday (day/sec)
                soilstate_inst%watfc_col(c,lev) = soilstate_inst%watsat_col(c,lev) * &
                     (0.1_r8 / (soilstate_inst%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))
             end if
          end do

          ! Urban pervious and impervious road
          if (col%itype(c) == icol_road_imperv) then
             ! Impervious road layers -- same as above except set watdry and watopt as missing
             do lev = 1,nlevgrnd
                soilstate_inst%watdry_col(c,lev) = spval 
                soilstate_inst%watopt_col(c,lev) = spval 
             end do
          else if (col%itype(c) == icol_road_perv) then 
             ! pervious road layers  - set in UrbanInitTimeConst
          end if

       end if
    end do

    ! --------------------------------------------------------------------
    ! Set soil hydraulic and thermal properties: lake
    ! --------------------------------------------------------------------

    do c = begc, endc
       g = col%gridcell(c)
       l = col%landunit(c)

       if (lun%itype(l)==istdlak) then

          do lev = 1,nlevgrnd
             if ( lev <= nlevsoi )then
                clay    =  soilstate_inst%cellclay_col(c,lev)
                sand    =  soilstate_inst%cellsand_col(c,lev)
                if ( organic_frac_squared )then
                   om_frac = (soilstate_inst%cellorg_col(c,lev)/organic_max)**2._r8
                else
                   om_frac = soilstate_inst%cellorg_col(c,lev)/organic_max
                end if
             else
                clay    = soilstate_inst%cellclay_col(c,nlevsoi)
                sand    = soilstate_inst%cellsand_col(c,nlevsoi)
                om_frac = 0.0_r8
             end if

             soilstate_inst%watsat_col(c,lev) = 0.489_r8 - 0.00126_r8*sand

             soilstate_inst%bsw_col(c,lev)    = 2.91 + 0.159*clay

             soilstate_inst%sucsat_col(c,lev) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )

             bd = (1._r8-soilstate_inst%watsat_col(c,lev))*2.7e3_r8

             soilstate_inst%watsat_col(c,lev) = (1._r8 - om_frac)*soilstate_inst%watsat_col(c,lev) + om_watsat_lake * om_frac

             tkm = (1._r8-om_frac)*(8.80_r8*sand+2.92_r8*clay)/(sand+clay) + om_tkm * om_frac ! W/(m K)

             soilstate_inst%bsw_col(c,lev)    = (1._r8-om_frac)*(2.91_r8 + 0.159_r8*clay) + om_frac * om_b_lake

             soilstate_inst%sucsat_col(c,lev) = (1._r8-om_frac)*soilstate_inst%sucsat_col(c,lev) + om_sucsat_lake * om_frac

             xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

             ! perc_frac is zero unless perf_frac greater than percolation threshold
             if (om_frac > pc_lake) then
                perc_norm = (1._r8 - pc_lake)**(-pcbeta)
                perc_frac = perc_norm*(om_frac - pc_lake)**pcbeta
             else
                perc_frac = 0._r8
             endif

             ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
             uncon_frac = (1._r8-om_frac) + (1._r8-perc_frac)*om_frac

             ! uncon_hksat is series addition of mineral/organic conductivites
             if (om_frac < 1._r8) then
                xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s
                uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat + ((1._r8-perc_frac)*om_frac)/om_hksat_lake)
             else
                uncon_hksat = 0._r8
             end if

             soilstate_inst%hksat_col(c,lev)  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat_lake
             soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8- soilstate_inst%watsat_col(c,lev))
             soilstate_inst%tksatu_col(c,lev) = soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)
             soilstate_inst%tkdry_col(c,lev)  = ((0.135_r8*bd + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd))*(1._r8-om_frac) + &
                                       om_tkd * om_frac
             soilstate_inst%csol_col(c,lev)   = ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) +   &
                                       om_csol * om_frac)*1.e6_r8  ! J/(m3 K)
             if (lev > nlevsoi) then
                soilstate_inst%csol_col(c,lev) = csol_bedrock
             endif

             soilstate_inst%watdry_col(c,lev) = soilstate_inst%watsat_col(c,lev) &
                  * (316230._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev))
             soilstate_inst%watopt_col(c,lev) = soilstate_inst%watsat_col(c,lev) &
                  * (158490._r8/soilstate_inst%sucsat_col(c,lev)) ** (-1._r8/soilstate_inst%bsw_col(c,lev))

             !! added by K.Sakaguchi for beta from Lee and Pielke, 1992
             ! water content at field capacity, defined as hk = 0.1 mm/day
             ! used eqn (7.70) in CLM3 technote with k = 0.1 (mm/day) / (# seconds/day)
             soilstate_inst%watfc_col(c,lev) = soilstate_inst%watsat_col(c,lev) * (0.1_r8 / &
                  (soilstate_inst%hksat_col(c,lev)*secspday))**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))
          end do
       endif

    end do

    ! --------------------------------------------------------------------
    ! Initialize threshold soil moisture and mass fracion of clay limited to 0.20
    ! --------------------------------------------------------------------

    do c = begc,endc
       g = col%gridcell(c)

       soilstate_inst%gwc_thr_col(c) = 0.17_r8 + 0.14_r8 * clay3d(g,1) * 0.01_r8
       soilstate_inst%mss_frc_cly_vld_col(c) = min(clay3d(g,1) * 0.01_r8, 0.20_r8)
    end do

    ! --------------------------------------------------------------------
    ! Deallocate memory
    ! --------------------------------------------------------------------

    deallocate(sand3d, clay3d, organic3d)
    deallocate(zisoifl, zsoifl, dzsoifl)

  end subroutine SoilStateInitTimeConst

end module SoilStateInitTimeConstMod
