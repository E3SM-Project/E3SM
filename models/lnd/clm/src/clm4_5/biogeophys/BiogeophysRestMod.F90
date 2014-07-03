module BiogeophysRestMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or biogeophysics restart/initial data
  !
  ! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : shr_infnan_isnan
  use spmdMod          , only : masterproc
  use decompMod        , only : bounds_type, get_proc_global
  use clm_time_manager , only : is_first_step
  use abortutils       , only : endrun
  use restUtilMod
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  ! save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BiogeophysRest
  !
  ! !PUBLIC DATA:
  logical, public  :: bound_h2osoi = .true.
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BiogeophysRest(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write biogeophysics information to/from restart file.
    !
    ! !USES:
    use clmtype
    use clm_varpar      , only : nlevgrnd, nlevsno, nlevlak, nlevurb, nlevsoi, nlevcan
    use clm_varcon      , only : denice, denh2o, istcrop, istdlak, istsoil, pondmx, watmin, spval  
    use clm_varcon      , only : icol_roof, icol_sunwall, icol_shadewall, zsoi
    use clm_varctl      , only : nsrest, nsrContinue, nsrStartup, nsrBranch, fpftdyn, iulog 
    use clm_varctl      , only : use_cndv, use_snicar_frc 
    use clm_atmlnd      , only : clm_a2l
    use SNICARMod       , only : snw_rds_min
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    real(r8)              :: maxwatsat    ! maximum porosity    
    real(r8)              :: excess       ! excess volumetric soil water
    real(r8)              :: totwat       ! total soil water (mm)
    integer               :: p,c,l,g,j,iv ! indices
    integer               :: numg_global  ! total number of grid cells, globally
    integer               :: numl_global  ! total number of landunits, globally
    integer               :: numc_global  ! total number of columns, globally
    integer               :: nump_global  ! total number of pfts, globally
    integer               :: nlevs        ! number of layers
    logical               :: readvar      ! determine if variable is on initial file
    logical               :: do_io        ! whether to do i/o for the given variable
    integer               :: dimlen       ! dimension length
    integer               :: err_code     ! error code
    real(r8), pointer     :: temp2d(:,:)  ! temporary for zisno
    character(len=7)      :: filetypes(0:3)
    character(len=*), parameter :: sub="BiogeophysRest"
    !-----------------------------------------------------------------------

    filetypes(:)           = "missing"
    filetypes(nsrStartup)  = "finidat"
    filetypes(nsrContinue) = "restart"
    filetypes(nsrBranch)   = "nrevsn"

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='zsoi', xtype=ncd_double, &
            dim1name='levgrnd', long_name='coordinate soil levels', units='m')
    else if (flag == 'write') then
       call ncd_io(ncid=ncid, varname='zsoi', data=zsoi, flag='write')
    end if

    ! Get expected total number of points, for later error checks
    call get_proc_global(numg_global, numl_global, numc_global, nump_global)

    ! Note - for the snow interfaces, are only examing the snow interfaces
    ! above zi=0 which is why zisno and zsno have the same level dimension below
    ! (Note - for zisno, zi(0) is set to 0 in routine iniTimeConst)

    ! pft energy flux - eflx_lwrad_out
    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=pef%eflx_lwrad_out)

    ! column water state variable - snow levels
    call restartvar(ncid=ncid, flag=flag, varname='SNLSNO', xtype=ncd_int,  & 
         dim1name='column', &
         long_name='number of snow layers', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=cps%snl)

    ! column water state variable - snow_depth
    call restartvar(ncid=ncid, flag=flag, varname='SNOW_DEPTH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='snow depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=cps%snow_depth) 

    ! column water state variable - int_snow
    call restartvar(ncid=ncid, flag=flag, varname='INT_SNOW', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='accuumulated snow', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=cws%int_snow)
    if (flag=='read' .and. .not. readvar) then
       cws%int_snow(:) = 0.0_r8
    end if
    
    ! perennial snow persistence - snow_persistence
    call restartvar(ncid=ncid, flag=flag, varname='SNOW_PERS', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='continuous snow cover time', units='sec', &
         interpinic_flag='interp', readvar=readvar, data=cps%snow_persistence)    
    if (flag=='read' .and. .not. readvar) then
         cps%snow_persistence(:) = 0.0_r8
    end if
       
    ! column water state variable - wa
    call restartvar(ncid=ncid, flag=flag, varname='WA', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water in the unconfined aquifer', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=cws%wa)

    ! column water state variable - zwt
    call restartvar(ncid=ncid, flag=flag, varname='ZWT', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=cws%zwt)

    ! column water state variable - frost_table
    call restartvar(ncid=ncid, flag=flag, varname='FROST_TABLE', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='frost table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=cws%frost_table)
    if (flag == 'read' .and. .not. readvar) then
       cws%frost_table(bounds%begc:bounds%endc) = cps%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

    ! column water state variable - zwt_perched
    call restartvar(ncid=ncid, flag=flag, varname='ZWT_PERCH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='perched water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=cws%zwt_perched)
    if (flag == 'read' .and. .not. readvar) then
       cws%zwt_perched(bounds%begc:bounds%endc) = cps%zi(bounds%begc:bounds%endc,nlevsoi)
    end if

    ! column type physical state variable - frac_sno_eff
    call restartvar(ncid=ncid, flag=flag, varname='frac_sno_eff', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='fraction of ground covered by snow (0 to 1)',units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=cps%frac_sno_eff)
    if (flag == 'read' .and. .not. readvar) then
       cps%frac_sno_eff(bounds%begc:bounds%endc) = 0.0_r8
    end if

    ! column type physical state variable - frac_sno
    call restartvar(ncid=ncid, flag=flag, varname='frac_sno', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='fraction of ground covered by snow (0 to 1)',units='unitless',&
         interpinic_flag='interp', readvar=readvar, data=cps%frac_sno)

    ! column type physical state variable - dzsno
    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = cps%dz(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='DZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer thickness', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       cps%dz(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    ! column type physical state variable - zsno
    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno+1:0))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) = cps%z(bounds%begc:bounds%endc,-nlevsno+1:0)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZSNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       cps%z(bounds%begc:bounds%endc,-nlevsno+1:0) = temp2d(bounds%begc:bounds%endc,-nlevsno+1:0) 
    end if
    deallocate(temp2d)

    ! column type physical state variable - zisno
    allocate(temp2d(bounds%begc:bounds%endc,-nlevsno:-1))
    if (flag == 'write') then
       temp2d(bounds%begc:bounds%endc,-nlevsno:-1) = cps%zi(bounds%begc:bounds%endc,-nlevsno:-1)
    end if
    call restartvar(ncid=ncid, flag=flag, varname='ZISNO', xtype=ncd_double,  & 
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno, upperb2=-1, &
         long_name='snow interface depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=temp2d)
    if (flag == 'read') then
       cps%zi(bounds%begc:bounds%endc,-nlevsno:-1) = temp2d(bounds%begc:bounds%endc,-nlevsno:-1) 
    end if
    deallocate(temp2d)

    ! column type physical state variable - coszen
    call restartvar(ncid=ncid, flag=flag, varname='coszen', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='cosine of solar zenith angle', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=cps%coszen)

    ! landunit type physical state variable - sabs_roof_dir
    call restartvar(ncid=ncid, flag=flag, varname='sabs_roof_dir', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='direct solar absorbed by roof per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_roof_dir)

    ! landunit type physical state variable - sabs_roof_dif
    call restartvar(ncid=ncid, flag=flag, varname='sabs_roof_dif', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='diffuse solar absorbed by roof per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_roof_dif)

    ! landunit type physical state variable - sabs_sunwall_dir
    call restartvar(ncid=ncid, flag=flag, varname='sabs_sunwall_dir', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='direct solar absorbed by sunwall per unit wall area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_sunwall_dir)

    ! landunit type physical state variable - sabs_sunwall_dif
    call restartvar(ncid=ncid, flag=flag, varname='sabs_sunwall_dif', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='diffuse solar absorbed by sunwall per unit wall area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_sunwall_dif)

    ! landunit type physical state variable - sabs_shadewall_dir
    call restartvar(ncid=ncid, flag=flag, varname='sabs_shadewall_dir', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='direct solar absorbed by shadewall per unit wall area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_shadewall_dir)

    ! landunit type physical state variable - sabs_shadewall_dif
    call restartvar(ncid=ncid, flag=flag, varname='sabs_shadewall_dif', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='diffuse solar absorbed by shadewall per unit wall area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_shadewall_dif)

    ! landunit type physical state variable - sabs_improad_dir
    call restartvar(ncid=ncid, flag=flag, varname='sabs_improad_dir', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='direct solar absorbed by impervious road per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_improad_dir)

    ! landunit type physical state variable - sabs_improad_dif
    call restartvar(ncid=ncid, flag=flag, varname='sabs_improad_dif', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='diffuse solar absorbed by impervious road per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_improad_dif)

    ! landunit type physical state variable - sabs_perroad_dir
    call restartvar(ncid=ncid, flag=flag, varname='sabs_perroad_dir', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='direct solar absorbed by pervious road per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_perroad_dir)

    ! landunit type physical state variable - sabs_perroad_dif
    call restartvar(ncid=ncid, flag=flag, varname='sabs_perroad_dif', xtype=ncd_double,  & 
         dim1name='landunit', dim2name='numrad', switchdim=.true., &
         long_name='diffuse solar absorbed by pervious road per unit ground area per unit incident flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=lps%sabs_perroad_dif)

    ! landunit type physical state variable - taf
    call restartvar(ncid=ncid, flag=flag, varname='taf', xtype=ncd_double,  & 
         dim1name='landunit', &
         long_name='urban canopy air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=lps%taf)

    ! landunit type physical state variable - qaf
    call restartvar(ncid=ncid, flag=flag, varname='qaf', xtype=ncd_double,  & 
         dim1name='landunit', &
         long_name='urban canopy specific humidity', units='kg/kg', &
         interpinic_flag='interp', readvar=readvar, data=lps%qaf)

    ! pft type physical state variable - albd
    call restartvar(ncid=ncid, flag=flag, varname='albd', xtype=ncd_double,  & 
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%albd)

    ! pft type physical state variable - albi
    call restartvar(ncid=ncid, flag=flag, varname='albi', xtype=ncd_double,  & 
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (diffuse) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%albi)

    ! column type physical state variable - albgrd
    call restartvar(ncid=ncid, flag=flag, varname='albgrd', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='ground albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%albgrd) 

    ! column type physical state variable - albgri
    call restartvar(ncid=ncid, flag=flag, varname='albgri', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='ground albedo (indirect) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%albgri)

    ! column type physical state variable - albsod
    call restartvar(ncid=ncid, flag=flag, varname='albsod', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='soil albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%albsod)

    ! column type physical state variable - albsoi
    call restartvar(ncid=ncid, flag=flag, varname='albsoi', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='soil albedo (indirect) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%albsoi)

    if (use_snicar_frc) then

       ! column type physical state variable - albgrd_bc
       call restartvar(ncid=ncid, flag=flag, varname='albgrd_bc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without BC (direct) (0 to 1)', units='', &
            interpinic_flag='interp',readvar=readvar, data=cps%albgrd_bc)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_bc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_bc to albgrd"
          cps%albgrd_bc(bounds%begc:bounds%endc,:) = cps%albgrd(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgri_bc
       call restartvar(ncid=ncid, flag=flag, varname='albgri_bc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without BC (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgri_bc)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_bc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_bc to albgri"
          cps%albgri_bc(bounds%begc:bounds%endc,:) = cps%albgri(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgrd_pur
       call restartvar(ncid=ncid, flag=flag, varname='albgrd_pur', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='pure snow ground albedo (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgrd_pur)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_pur in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_pur to albgrd"
          cps%albgrd_pur(bounds%begc:bounds%endc,:) = cps%albgrd(bounds%begc:bounds%endc,:)
       end if
       
       ! column type physical state variable - albgri_pur
       call restartvar(ncid=ncid, flag=flag, varname='albgri_pur', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='pure snow ground albedo (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgri_pur)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_pur in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_pur to albgri"
          cps%albgri_pur(bounds%begc:bounds%endc,:) = cps%albgri(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgrd_oc
       call restartvar(ncid=ncid, flag=flag, varname='albgrd_oc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without OC (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgrd_oc)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_oc in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_oc to albgrd"
          cps%albgrd_oc(bounds%begc:bounds%endc,:) = cps%albgrd(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgri_oc
       call restartvar(ncid=ncid, flag=flag, varname='albgri_oc', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without OC (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgri_oc)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_oc in restart (or initial) file..."
          if (masterproc) write(iulog,*) "Initialize albgri_oc to albgri"
          cps%albgri_oc(bounds%begc:bounds%endc,:) = cps%albgri(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgrd_dst
       call restartvar(ncid=ncid, flag=flag, varname='albgrd_dst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without dust (direct) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgrd_dst)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgrd_dst in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgrd_dst to albgrd"
          cps%albgrd_dst(bounds%begc:bounds%endc,:) = cps%albgrd(bounds%begc:bounds%endc,:)
       end if

       ! column type physical state variable - albgri_dst
       call restartvar(ncid=ncid, flag=flag, varname='albgri_dst', xtype=ncd_double,  &
            dim1name='column', dim2name='numrad', switchdim=.true., &
            long_name='ground albedo without dust (diffuse) (0 to 1)', units='', &
            interpinic_flag='interp', readvar=readvar, data=cps%albgri_dst)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "SNICAR: can't find albgri_dst in initial file..."
          if (masterproc) write(iulog,*) "Initialize albgri_dst to albgri"
          cps%albgri_dst(bounds%begc:bounds%endc,:) = cps%albgri(bounds%begc:bounds%endc,:)
       end if

    end if  ! end of if-use_snicar_frc 

    ! column water state variable - h2osfc
    call restartvar(ncid=ncid, flag=flag, varname='H2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=cws%h2osfc)
    if (flag=='read' .and. .not. readvar) then
       cws%h2osfc(bounds%begc:bounds%endc) = 0.0_r8
    end if

    ! column type physical state variable - frac_h2osfc
    call restartvar(ncid=ncid, flag=flag, varname='FH2OSFC', xtype=ncd_double,  &
         dim1name='column',&
         long_name='fraction of ground covered by h2osfc (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%frac_h2osfc)
    if (flag == 'read' .and. .not. readvar) then
       cps%frac_h2osfc(bounds%begc:bounds%endc) = 0.0_r8
    end if

    ! column energy state variable - t_h2osfc
    call restartvar(ncid=ncid, flag=flag, varname='TH2OSFC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='surface water temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_h2osfc)
    if (flag=='read' .and. .not. readvar) then
       ces%t_h2osfc(bounds%begc:bounds%endc) = 274.0_r8
    end if

    ! column water state variable - h2osno
    call restartvar(ncid=ncid, flag=flag, varname='H2OSNO', xtype=ncd_double,  &
         dim1name='column', &
         long_name='snow water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=cws%h2osno)

    ! column water state variable - h2osoi_liq
    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_LIQ', xtype=ncd_double,  &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='liquid water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=cws%h2osoi_liq)

    ! column water state variable - h2osoi_ice
    call restartvar(ncid=ncid, flag=flag, varname='H2OSOI_ICE', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='ice lens', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=cws%h2osoi_ice)
         
    ! column energy state variable - t_grnd
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND', xtype=ncd_double,  &
         dim1name='column', &
         long_name='ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_grnd)

    ! column urban energy state variable - eflx_urban_ac
    call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC', xtype=ncd_double,  &
         dim1name='column', &
         long_name='urban air conditioning flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=cef%eflx_urban_ac)

    ! column urban energy state variable - eflx_urban_heat
    call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT', xtype=ncd_double,  &
         dim1name='column', &
         long_name='urban heating flux', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=cef%eflx_urban_heat)

    ! pft energy state variable - t_ref2m_min
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min)

    ! pft energy state variable - t_ref2m_max
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max)

    ! pft energy state variable - t_ref2m_min_inst
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min_inst)

    ! pft energy state variable - t_ref2m_max_inst
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max_inst)

    ! pft energy state variable - t_ref2m_u
    call restartvar(ncid=ncid, flag=flag, varname="T_REF2M_U", xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Urban 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_u)

    ! column energy state variable - t_grnd_u
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_U', xtype=ncd_double,  &
         dim1name='column', &
         long_name='urban ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_grnd_u)

    ! pft energy state variable - t_ref2m_min_u
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_U', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='urban daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min_u)

    ! pft energy state variable - t_ref2m_max_u
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_U', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='urban daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max_u)

    ! pft energy state variable - t_ref2m_min_inst_u
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_U', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='urban instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min_inst_u)

    ! pft energy state variable - t_ref2m_max_inst_u
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_U', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='urban instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max_inst_u)

    ! pft energy state variable - t_ref2m_r
    call restartvar(ncid=ncid, flag=flag, varname="T_REF2M_R", xtype=ncd_double,  &
         dim1name='pft', &
         long_name='Rural 2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_r)

    ! column energy state variable - t_grnd_r
    call restartvar(ncid=ncid, flag=flag, varname='T_GRND_R', xtype=ncd_double,  &
         dim1name='column', &
         long_name='rural ground temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_grnd_r)
         
    ! pft energy state variable - t_ref2m_min_r
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily minimum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min_r)
         
    ! pft energy state variable - t_ref2m_max_r
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural daily maximum of average 2 m height surface air temperature (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max_r)
         
    ! pft energy state variable - t_ref2m_min_inst_r
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MIN_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily min of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_min_inst_r)

    ! pft energy state variable - t_ref2m_max_inst_r
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M_MAX_INST_R', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='rural instantaneous daily max of average 2 m height surface air temp (K)', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m_max_inst_r)

    ! column energy state variable - t_soisno
    call restartvar(ncid=ncid, flag=flag, varname='T_SOISNO', xtype=ncd_double,   &
         dim1name='column', dim2name='levtot', switchdim=.true., &
         long_name='soil-snow temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_soisno)

    ! column type energy state variable - t_lake
    call restartvar(ncid=ncid, flag=flag, varname='T_LAKE', xtype=ncd_double,  &
         dim1name='column', dim2name='levlak', switchdim=.true., &
         long_name='lake temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=ces%t_lake)

    ! pft physical state variable - frac_veg_nosno_alb
    call restartvar(ncid=ncid, flag=flag, varname='FRAC_VEG_NOSNO_ALB', xtype=ncd_int,  &
         dim1name='pft',&
         long_name='fraction of vegetation not covered by snow (0 or 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%frac_veg_nosno_alb)

    ! pft type physical state variable - fwet
    call restartvar(ncid=ncid, flag=flag, varname='FWET', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='fraction of canopy that is wet (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fwet)

    ! pft type physical state variable - tlai
    call restartvar(ncid=ncid, flag=flag, varname='tlai', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='one-sided leaf area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%tlai)

    ! pft type physical state variable - tsai
    call restartvar(ncid=ncid, flag=flag, varname='tsai', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='one-sided stem area index, no burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%tsai)

    ! pft type physical state variable - tlai_z
    call restartvar(ncid=ncid, flag=flag, varname='tlai_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='tlai increment for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%tlai_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find tlai_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize tlai_z to tlai/nlevcan" 
       do iv=1,nlevcan
          pps%tlai_z(bounds%begp:bounds%endp,iv) = pps%tlai(bounds%begp:bounds%endp)/nlevcan
       end do
    end if

    ! pft type physical state variable - tsai_z
    call restartvar(ncid=ncid, flag=flag, varname='tsai_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='tsai increment for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%tsai_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find tsai_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize tsai_z to tsai/nlevcan" 
       do iv=1,nlevcan
          pps%tsai_z(bounds%begp:bounds%endp,iv) = pps%tsai(bounds%begp:bounds%endp)/nlevcan
       end do
    end if

    ! pft type physical state variable - ncan
    call restartvar(ncid=ncid, flag=flag, varname='ncan', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='number of canopy layers', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%ncan)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find ncan in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize ncan to nlevcan" 
       pps%ncan(bounds%begp:bounds%endp) = nlevcan
    end if

    ! pft type physical state variable - nrad
    call restartvar(ncid=ncid, flag=flag, varname='nrad', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='number of canopy layers, above snow for radiative transfer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%nrad)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find nrad in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize nrad to nlevcan" 
       pps%nrad(bounds%begp:bounds%endp) = nlevcan
    end if
    
    ! pft type physica state variable - mlaidiff
    call restartvar(ncid=ncid, flag=flag, varname='mlaidiff', xtype=ncd_double,  &
         dim1name='pft',&
         long_name='difference between lai month one and month two', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%mlaidiff)

    ! pft type physical state variable - elai
    call restartvar(ncid=ncid, flag=flag, varname='elai', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='one-sided leaf area index, with burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%elai)
    
    ! pft type physical state variable - esai
    call restartvar(ncid=ncid, flag=flag, varname='esai', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='one-sided stem area index, with burying by snow', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%esai)
    
    ! pft type physical state variable - fsun
    call restartvar(ncid=ncid, flag=flag, varname='fsun', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='sunlit fraction of canopy', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fsun)
    if (flag=='read' )then
       do p = bounds%begp,bounds%endp
          if (shr_infnan_isnan( pps%fsun(p)) ) then
             pps%fsun(p) = spval
          end if
       end do
    end if
    
    ! pft type physical state variable - vcmaxcintsun
    call restartvar(ncid=ncid, flag=flag, varname='vcmaxcintsun', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='sunlit canopy scaling coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%vcmaxcintsun)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find vcmaxcintsun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize vcmaxcintsun to 1"
       pps%vcmaxcintsun(bounds%begp:bounds%endp) = 1._r8
    end if

    ! pft type physical state variable - vcmaxcintsha
    call restartvar(ncid=ncid, flag=flag, varname='vcmaxcintsha', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='shaded canopy scaling coefficient', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%vcmaxcintsha)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find vcmaxcintsha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize vcmaxcintsha to 1"
       pps%vcmaxcintsha(bounds%begp:bounds%endp) = 1._r8
    end if

    ! pft type physical state variable - htop
    call restartvar(ncid=ncid, flag=flag, varname='htop', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='canopy top', units='m', &
         interpinic_flag='interp', readvar=readvar, data=pps%htop)

    ! pft type physical state variable - hbot
    call restartvar(ncid=ncid, flag=flag, varname='hbot', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='canopy botton', units='m', &
         interpinic_flag='interp', readvar=readvar, data=pps%hbot)

    ! pft type physical state variable - fabd
    call restartvar(ncid=ncid, flag=flag, varname='fabd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabd)

    ! pft type physical state variable - fabi
    call restartvar(ncid=ncid, flag=flag, varname='fabi', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by veg per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabi)

    ! pft type physical state variable - fabd_sun
    call restartvar(ncid=ncid, flag=flag, varname='fabd_sun', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by sunlit leaf per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabd_sun)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sun to fabd/2"
       pps%fabd_sun(bounds%begp:bounds%endp,:) = pps%fabd(bounds%begp:bounds%endp,:)/2._r8
    end if

    ! pft type physical state variable - fabd_sha
    call restartvar(ncid=ncid, flag=flag, varname='fabd_sha', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by shaded leaf per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabd_sha)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sha to fabd/2"
       pps%fabd_sha(bounds%begp:bounds%endp,:) = pps%fabd(bounds%begp:bounds%endp,:)/2._r8
    end if

    ! pft type physical state variable - fabi_sun
    call restartvar(ncid=ncid, flag=flag, varname='fabi_sun', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by sunlit leaf per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabi_sun)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sun in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sun to fabi/2"
       pps%fabi_sun(bounds%begp:bounds%endp,:) = pps%fabi(bounds%begp:bounds%endp,:)/2._r8
    end if

    ! pft type physical state variable - fabi_sha
    call restartvar(ncid=ncid, flag=flag, varname='fabi_sha', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='flux absorbed by shaded leaf per unit diffuse flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabi_sha)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sha in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sha to fabi/2"
       pps%fabi_sha(bounds%begp:bounds%endp,:) = pps%fabi(bounds%begp:bounds%endp,:)/2._r8
    end if

    ! pft type physical state variable - fabd_sun_z
    call restartvar(ncid=ncid, flag=flag, varname='fabd_sun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed sunlit leaf direct PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabd_sun_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sun_z to (fabd/2)/nlevcan" 
       do iv=1,nlevcan
          pps%fabd_sun_z(bounds%begp:bounds%endp,iv) = (pps%fabd(bounds%begp:bounds%endp,1)/2._r8)/nlevcan
       end do
    end if

    ! pft type physical state variable - fabd_sha_z
    call restartvar(ncid=ncid, flag=flag, varname='fabd_sha_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed shaded leaf direct PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabd_sha_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabd_sha_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabd_sha_z to (fabd/2)/nlevcan" 
       do iv=1,nlevcan
          pps%fabd_sha_z(bounds%begp:bounds%endp,iv) = (pps%fabd(bounds%begp:bounds%endp,1)/2._r8)/nlevcan
       end do
    end if

    ! pft type physical state variable - fabi_sun_z
    call restartvar(ncid=ncid, flag=flag, varname='fabi_sun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed sunlit leaf diffuse PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabi_sun_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sun_z to (fabi/2)/nlevcan"
       do iv=1,nlevcan
          pps%fabi_sun_z(bounds%begp:bounds%endp,iv) = (pps%fabi(bounds%begp:bounds%endp,1)/2._r8)/nlevcan
       end do
    end if

    ! pft type physical state variable - fabi_sha_z
    call restartvar(ncid=ncid, flag=flag, varname='fabi_sha_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='absorbed shaded leaf diffuse PAR (per unit lai+sai) for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fabi_sha_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fabi_sha_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fabi_sha_z to (fabi/2)/nlevcan"
       do iv=1,nlevcan
          pps%fabi_sha_z(bounds%begp:bounds%endp,iv) = (pps%fabi(bounds%begp:bounds%endp,1)/2._r8)/nlevcan
       end do
    end if

    ! pft type physical state variable - fsun_z
    call restartvar(ncid=ncid, flag=flag, varname='fsun_z', xtype=ncd_double,  &
         dim1name='pft', dim2name='levcan', switchdim=.true., &
         long_name='sunlit fraction for canopy layer', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%fsun_z)
    if (flag=='read' .and. .not. readvar) then
       if (masterproc) write(iulog,*) "can't find fsun_z in restart (or initial) file..."
       if (masterproc) write(iulog,*) "Initialize fsun_z to 0"
       do iv=1,nlevcan
          pps%fsun_z(bounds%begp:bounds%endp,iv) = 0._r8
       end do
    end if

    ! pft type physical state variable - ftdd
    call restartvar(ncid=ncid, flag=flag, varname='ftdd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down direct flux below veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%ftdd)

    ! pft type physical state variable - ftid
    call restartvar(ncid=ncid, flag=flag, varname='ftid', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down diffuse flux below veg per unit direct flux', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%ftid)

    ! pft type physical state variable - ftii
    call restartvar(ncid=ncid, flag=flag, varname='ftii', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='down diffuse flux below veg per unit diffuse flux', units='', &      
         interpinic_flag='interp', readvar=readvar, data=pps%ftii)

    ! pft energy state variable - t_veg
    call restartvar(ncid=ncid, flag=flag, varname='T_VEG', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='vegetation temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_veg)

    ! pft energy state variable - t_ref2m
    call restartvar(ncid=ncid, flag=flag, varname='T_REF2M', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='2m height surface air temperature', units='K', &
         interpinic_flag='interp', readvar=readvar, data=pes%t_ref2m)
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! pft type water state variable - h2ocan
    call restartvar(ncid=ncid, flag=flag, varname='H2OCAN', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='canopy water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=pws%h2ocan)
    
    ! pft irrigation variable - n_irrig_steps_left
    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
          readvar = .false.
       end if
    else if (flag == 'define' .or. do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=pps%n_irrig_steps_left)
       if (flag=='read' .and. .not. readvar) then
          pps%n_irrig_steps_left = 0
       end if
    end if

    ! pft irrigation variable - irrig_rate
    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
          readvar = .false.
       end if
    else if (flag == 'define' .or. do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=pps%irrig_rate)
       if (flag=='read' .and. .not. readvar) then
          pps%irrig_rate = 0.0_r8
       end if
    end if

    ! ------------------------------------------------------------
    ! Determine volumetric soil water (for read only)
    ! ------------------------------------------------------------

    if (flag == 'read' ) then

       do c = bounds%begc, bounds%endc
          l = col%landunit(c)
          if ( col%itype(c) == icol_sunwall   .or. &
               col%itype(c) == icol_shadewall .or. &
               col%itype(c) == icol_roof )then
             nlevs = nlevurb
          else
             nlevs = nlevgrnd
          end if
          ! NOTE: THIS IS A MEMORY INEFFICIENT COPY
          if ( lun%itype(l) /= istdlak ) then ! This calculation is now done for lakes in initSLake.
             do j = 1,nlevs
                cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                    + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
             end do
          end if
       end do
    end if

    ! ------------------------------------------------------------
    ! If initial run -- ensure that water is properly bounded (read only)
    ! ------------------------------------------------------------

    if (flag == 'read' ) then
       if ( is_first_step() .and. bound_h2osoi) then
          do c = bounds%begc, bounds%endc
             l = col%landunit(c)
             if ( col%itype(c) == icol_sunwall .or. col%itype(c) == icol_shadewall .or. &
                  col%itype(c) == icol_roof )then
                nlevs = nlevurb
             else
                nlevs = nlevgrnd
             end if
             do j = 1,nlevs
                l = col%landunit(c)
                if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
                   cws%h2osoi_liq(c,j) = max(0._r8,cws%h2osoi_liq(c,j))
                   cws%h2osoi_ice(c,j) = max(0._r8,cws%h2osoi_ice(c,j))
                   cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                       + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
                   if (j == 1) then
                      maxwatsat = (cps%watsat(c,j)*cps%dz(c,j)*1000.0_r8 + pondmx) / &
                                  (cps%dz(c,j)*1000.0_r8)
                   else
                      maxwatsat = cps%watsat(c,j)
                   end if
                   if (cws%h2osoi_vol(c,j) > maxwatsat) then 
                      excess = (cws%h2osoi_vol(c,j) - maxwatsat)*cps%dz(c,j)*1000.0_r8
                      totwat = cws%h2osoi_liq(c,j) + cws%h2osoi_ice(c,j)
                      cws%h2osoi_liq(c,j) = cws%h2osoi_liq(c,j) - &
                                           (cws%h2osoi_liq(c,j)/totwat) * excess
                      cws%h2osoi_ice(c,j) = cws%h2osoi_ice(c,j) - &
                                           (cws%h2osoi_ice(c,j)/totwat) * excess
                   end if
                   cws%h2osoi_liq(c,j) = max(watmin,cws%h2osoi_liq(c,j))
                   cws%h2osoi_ice(c,j) = max(watmin,cws%h2osoi_ice(c,j))
                   cws%h2osoi_vol(c,j) = cws%h2osoi_liq(c,j)/(cps%dz(c,j)*denh2o) &
                                       + cws%h2osoi_ice(c,j)/(cps%dz(c,j)*denice)
                end if
             end do
          end do
       end if

    endif   ! end if if-read flag

    !--------------------------------
    ! variables needed for SNICAR
    !--------------------------------

    ! column type physical state variable - snw_rds
    call restartvar(ncid=ncid, flag=flag, varname='snw_rds', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer effective radius', units='um', &
         interpinic_flag='interp', readvar=readvar, data=cps%snw_rds)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize snw_rds
       if (masterproc) then
          write(iulog,*) "SNICAR: This is an initial run (not a restart), and grain size/aerosol " // &
               "mass data are not defined in initial condition file. Initialize snow " // &
               "effective radius to fresh snow value, and snow/aerosol masses to zero."
       endif
       do c= bounds%begc, bounds%endc
          if (cps%snl(c) < 0) then
             cps%snw_rds(c,cps%snl(c)+1:0) = snw_rds_min
             cps%snw_rds(c,-nlevsno+1:cps%snl(c)) = 0._r8
             cps%snw_rds_top(c) = snw_rds_min
             cps%sno_liq_top(c) = cws%h2osoi_liq(c,cps%snl(c)+1) / &
                                 (cws%h2osoi_liq(c,cps%snl(c)+1)+cws%h2osoi_ice(c,cps%snl(c)+1))
          elseif (cws%h2osno(c) > 0._r8) then
             cps%snw_rds(c,0) = snw_rds_min
             cps%snw_rds(c,-nlevsno+1:-1) = 0._r8
             cps%snw_rds_top(c) = spval
             cps%sno_liq_top(c) = spval
          else
             cps%snw_rds(c,:) = 0._r8
             cps%snw_rds_top(c) = spval
             cps%sno_liq_top(c) = spval
          endif
       enddo
    endif

    ! column type physical state variable - mss_bcpho
    call restartvar(ncid=ncid, flag=flag, varname='mss_bcpho', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophobic black carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_bcpho)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_bcpho to zero
       cps%mss_bcpho(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! column type physical state variable - mss_bcphi
    call restartvar(ncid=ncid, flag=flag, varname='mss_bcphi', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophilic black carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_bcphi)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_bcphi to zero
       cps%mss_bcphi(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! column type physical state variable - mss_ocpho
    call restartvar(ncid=ncid, flag=flag, varname='mss_ocpho', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophobic organic carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_ocpho)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_ocpho to zero
       cps%mss_ocpho(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! column type physical state variable - mss_ocphi
    call restartvar(ncid=ncid, flag=flag, varname='mss_ocphi', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer hydrophilic organic carbon mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_ocphi)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_ocphi to zero
       cps%mss_ocphi(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! column type physical state variable - mss_dst1
    call restartvar(ncid=ncid, flag=flag, varname='mss_dst1', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 1 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_dst1)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst1 to zero
       cps%mss_dst1(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if
    
    ! column type physical state variable - mss_dst2
    call restartvar(ncid=ncid, flag=flag, varname='mss_dst2', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 2 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_dst2)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst2 to zero
       cps%mss_dst2(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    ! column type physical state variable - mss_dst3
    call restartvar(ncid=ncid, flag=flag, varname='mss_dst3', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 3 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar,  data=cps%mss_dst3)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst3 to zero
       cps%mss_dst3(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    ! column type physical state variable - mss_dst4
    call restartvar(ncid=ncid, flag=flag, varname='mss_dst4', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer dust species 4 mass', units='kg m-2', &
         interpinic_flag='interp', readvar=readvar, data=cps%mss_dst4)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize mss_dst4 to zero
       cps%mss_dst4(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    end if

    ! column type physical state variable - flx_absdv
    call restartvar(ncid=ncid, flag=flag, varname='flx_absdv', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (direct, VIS)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=cps%flx_absdv)

    ! column type physical state variable - flx_absdn
    call restartvar(ncid=ncid, flag=flag, varname='flx_absdn', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (direct, NIR)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=cps%flx_absdn)

    ! column type physical state variable - flx_absiv
    call restartvar(ncid=ncid, flag=flag, varname='flx_absiv', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (diffuse, VIS)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=cps%flx_absiv)

    ! column type physical state variable - flx_absin
    call restartvar(ncid=ncid, flag=flag, varname='flx_absin', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno1', switchdim=.true., lowerb2=-nlevsno+1, upperb2=1, &
         long_name='snow layer flux absorption factors (diffuse, NIR)', units='fraction', &
         interpinic_flag='interp', readvar=readvar, data=cps%flx_absin)

    ! column type physical state variable - albsnd_hst
    call restartvar(ncid=ncid, flag=flag, varname='albsnd_hst', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='snow albedo (direct) (0 to 1)', units='proportion', &
         interpinic_flag='interp', readvar=readvar, data=cps%albsnd_hst)

    ! column type physical state variable - albsni_hst
    call restartvar(ncid=ncid, flag=flag, varname='albsni_hst', xtype=ncd_double,  &
         dim1name='column', dim2name='numrad', switchdim=.true., &
         long_name='snow albedo (diffuse) (0 to 1)', units='proportion', &
         interpinic_flag='interp', readvar=readvar, data=cps%albsni_hst)

    ! column type water flux variable - qflx_snofrz_lyr
    call restartvar(ncid=ncid, flag=flag, varname='qflx_snofrz_lyr', xtype=ncd_double,  &
         dim1name='column', dim2name='levsno', switchdim=.true., lowerb2=-nlevsno+1, upperb2=0, &
         long_name='snow layer ice freezing rate', units='kg m-2 s-1', &
         interpinic_flag='interp', readvar=readvar, data=cwf%qflx_snofrz_lyr)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snofrz_lyr to zero
       cwf%qflx_snofrz_lyr(bounds%begc:bounds%endc,-nlevsno+1:0) = 0._r8
    endif

    ! column type water flux variable - qflx_snow_melt
    call restartvar(ncid=ncid, flag=flag, varname='qflx_snow_melt', xtype=ncd_double,  &
         dim1name='column', &
         long_name='net snow melt', units='mm/s', &
         interpinic_flag='interp', readvar=readvar, data=cwf%qflx_snow_melt)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, not restart: initialize qflx_snow_melt to zero
       cwf%qflx_snow_melt = 0._r8
    endif

    ! gridcell type water flux variable - qflx_floodg 
    !TODO - what should this be in terms of interpolation? For now set it to skip
    call restartvar(ncid=ncid, flag=flag, varname='qflx_floodg', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='flood water flux', units='mm/s', &
         interpinic_flag='skip', readvar=readvar, data=clm_a2l%forc_flood)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, readvar=readvar, not restart: initialize flood to zero
       clm_a2l%forc_flood = 0._r8
    endif

    ! initialize other variables that are derived from those
    ! stored in the restart buffer. (there may be a more appropriate
    ! place to do this, but functionally this works)
    if (flag == 'read' ) then
       do j = -nlevsno+1,0
          do c = bounds%begc, bounds%endc
             ! mass concentrations of aerosols in snow
             if (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j) > 0._r8) then
                cps%mss_cnc_bcpho(c,j) = cps%mss_bcpho(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_bcphi(c,j) = cps%mss_bcphi(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_ocpho(c,j) = cps%mss_ocpho(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_ocphi(c,j) = cps%mss_ocphi(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))

                cps%mss_cnc_dst1(c,j) = cps%mss_dst1(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst2(c,j) = cps%mss_dst2(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst3(c,j) = cps%mss_dst3(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
                cps%mss_cnc_dst4(c,j) = cps%mss_dst4(c,j) / (cws%h2osoi_ice(c,j)+cws%h2osoi_liq(c,j))
             else
                cps%mss_cnc_bcpho(c,j) = 0._r8
                cps%mss_cnc_bcphi(c,j) = 0._r8
                cps%mss_cnc_ocpho(c,j) = 0._r8
                cps%mss_cnc_ocphi(c,j) = 0._r8

                cps%mss_cnc_dst1(c,j) = 0._r8
                cps%mss_cnc_dst2(c,j) = 0._r8
                cps%mss_cnc_dst3(c,j) = 0._r8
                cps%mss_cnc_dst4(c,j) = 0._r8
             endif
          enddo
       enddo
    endif
    !-- SNICAR variables

  end subroutine BiogeophysRest

end module BiogeophysRestMod
