module mkfileMod

contains

!-----------------------------------------------------------------------
  subroutine mkfile(domain, fname, dynlanduse)

    use shr_kind_mod , only : r8 => shr_kind_r8
    use shr_sys_mod  , only : shr_sys_getenv
    use fileutils    , only : get_filename
    use mkvarpar     , only : nlevsoi, nlevurb, numrad
    use mkvarctl
    use mkurbanparMod, only : numurbl
    use mkglcmecMod  , only : nglcec
    use mkpftMod     , only : mkpftAtt
    use mksoilMod    , only : mksoilAtt
    use mkharvestMod , only : mkharvest_fieldname, mkharvest_numtypes, mkharvest_longname
    use mkncdio      , only : check_ret, ncd_defvar
    use mkdomainMod  

    implicit none
    include 'netcdf.inc'
    type(domain_type) , intent(in) :: domain
    character(len=*)  , intent(in) :: fname
    logical           , intent(in) :: dynlanduse

    integer :: ncid
    integer :: j                    ! index
    integer :: dimid                ! temporary
    integer :: values(8)            ! temporary
    character(len=256) :: str       ! global attribute string
    character(len=256) :: name      ! name of attribute
    character(len=256) :: unit      ! units of attribute
    character(len= 18) :: datetime  ! temporary
    character(len=  8) :: date      ! temporary
    character(len= 10) :: time      ! temporary
    character(len=  5) :: zone      ! temporary
    integer            :: ier       ! error status
    integer            :: omode     ! netCDF output mode
    integer            :: xtype     ! external type
    character(len=32)  :: subname = 'mkfile'  ! subroutine name
!-----------------------------------------------------------------------

    call check_ret(nf_create(trim(fname), ior(nf_clobber,nf_64bit_offset), &
                                ncid), subname)

    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Define dimensions.

    if (outnc_1d) then
       call check_ret(nf_def_dim (ncid, 'gridcell', domain%ns, dimid), subname)
    else
       call check_ret(nf_def_dim (ncid, 'lsmlon'  , domain%ni, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lsmlat'  , domain%nj, dimid), subname)
    end if

    if (.not. dynlanduse) then
       if ( nglcec > 0 )then
          call check_ret(nf_def_dim (ncid, 'nglcec'  , nglcec      , dimid), subname)
          call check_ret(nf_def_dim (ncid, 'nglcecp1', nglcec+1    , dimid), subname)
       end if
    end if
    call check_ret(nf_def_dim (ncid, 'numurbl' , numurbl     , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'nlevurb' , nlevurb     , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'numrad'  , numrad      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'nchar'   , 256         , dimid), subname)

    ! Create global attributes.

    str = 'NCAR-CSM'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Conventions', len_trim(str), trim(str)), subname)

    call date_and_time (date, time, zone, values)
    datetime(1:8) =        date(5:6) // '-' // date(7:8) // '-' // date(3:4)
    datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '
    str = 'created on: ' // datetime
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'History_Log', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('LOGNAME', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Logname', len_trim(str), trim(str)), subname)

    call shr_sys_getenv ('HOST', str, ier)
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Host', len_trim(str), trim(str)), subname)

    str = 'Community Land Model: CLM4'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Source', len_trim(str), trim(str)), subname)

    str = &
'$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_5_1_r085/models/lnd/clm/tools/clm4_5/mksurfdata_map/src/mkfileMod.F90 $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id: mkfileMod.F90 47951 2013-06-12 11:13:58Z sacks $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Revision_Id', len_trim(str), trim(str)), subname)

#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif

    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Compiler_Optimized', len_trim(str), trim(str)), subname)

    if ( all_urban )then
       str = 'TRUE'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'all_urban', len_trim(str), trim(str)), subname)
    end if

    if ( no_inlandwet )then
       str = 'TRUE'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'no_inlandwet', len_trim(str), trim(str)), subname)
    end if

    call check_ret(nf_put_att_int(ncid, NF_GLOBAL, &
         'nglcec', nf_int, 1, nglcec), subname)

    ! Raw data file names

    str = get_filename(mksrf_fgrid)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Input_grid_dataset', len_trim(str), trim(str)), subname)

    str = trim(mksrf_gridtype)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Input_gridtype', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_fvocef)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'VOC_EF_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if

    str = get_filename(mksrf_flakwat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Inland_lake_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fwetlnd)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Inland_wetland_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fglacier)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Glacier_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_furbtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Urban_Topography_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_flndtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Land_Topography_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_furban)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Urban_raw_data_file_name', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_flai)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if

    str = get_filename(mksrf_fabm)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'agfirepkmon_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fgdp)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'gdp_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fpeat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'peatland_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_ftopostats)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'topography_stats_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fvic)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'vic_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fch4)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'ch4_params_raw_data_file_name', len_trim(str), trim(str)), subname)

    ! Mapping file names

    str = get_filename(map_fpft)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_pft_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(map_flakwat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_lakwat_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fwetlnd)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_wetlnd_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fglacier)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_glacier_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fsoitex)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_texture_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fsoicol)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_color_file', len_trim(str), trim(str)), subname)
  
     str = get_filename(map_fsoiord)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_order_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_forganic)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_organic_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_furban)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_urban_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fmax)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_fmax_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fvocef)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_VOC_EF_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fharvest)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_harvest_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_flai)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_lai_sai_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_furbtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_urban_topography_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_flndtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_land_topography_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fabm)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_agfirepkmon_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fgdp)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_gdp_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fpeat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_peatland_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_ftopostats)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_topography_stats_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fvic)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_vic_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fch4)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_ch4_params_file', len_trim(str), trim(str)), subname)

    ! ----------------------------------------------------------------------
    ! Define variables
    ! ----------------------------------------------------------------------
 
    if ( .not. outnc_double )then
       xtype = nf_float
    else
       xtype = nf_double
    end if

    call mksoilAtt( ncid, dynlanduse, xtype )

    call mkpftAtt(  ncid, dynlanduse, xtype )

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='AREA' , xtype=nf_double, &
            dim1name='gridcell',&
            long_name='area', units='km^2')
    else
       call ncd_defvar(ncid=ncid, varname='AREA' , xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='area', units='km^2')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LONGXY', xtype=nf_double, &
            dim1name='gridcell',&
            long_name='longitude', units='degrees east')
    else
       call ncd_defvar(ncid=ncid, varname='LONGXY', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='longitude', units='degrees east')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LATIXY', xtype=nf_double, &
            dim1name='gridcell',&
            long_name='latitude', units='degrees north')
    else
       call ncd_defvar(ncid=ncid, varname='LATIXY', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='latitude', units='degrees north')
    end if

    if (.not. dynlanduse) then
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_BTR', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF btr (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_BTR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF btr (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_FET', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF fet (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_FET', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF fet (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_FDT', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF fdt (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_FDT', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF fdt (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_SHR', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF shr (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_SHR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF shr (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_GRS', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF grs (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_GRS', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF grs (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EF1_CRP', xtype=xtype, &
               dim1name='gridcell',&
               long_name='EF crp (isoprene)', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EF1_CRP', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='EF crp (isoprene)', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CANYON_HWR', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='canyon height to width ratio', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='CANYON_HWR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='canyon height to width ratio', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='emissivity of impervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='emissivity of impervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_PERROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='emissivity of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_PERROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='emissivity of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='emissivity of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='emissivity of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='emissivity of wall', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='emissivity of wall', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='HT_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='height of roof', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='HT_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='height of roof', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='THICK_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='thickness of roof', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='THICK_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='thickness of roof', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='THICK_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='thickness of wall', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='THICK_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='thickness of wall', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MAX', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='maximum interior building temperature', units='K')
       else
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MAX', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='maximum interior building temperature', units='K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MIN', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='minimum interior building temperature', units='K')
       else
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MIN', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='minimum interior building temperature', units='K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WIND_HGT_CANYON', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='height of wind in canyon', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='WIND_HGT_CANYON', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='height of wind in canyon', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WTLUNIT_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='fraction of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='WTLUNIT_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='fraction of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WTROAD_PERV', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='fraction of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='WTROAD_PERV', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='fraction of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD_DIR', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='direct albedo of impervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD_DIR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='direct albedo of impervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD_DIF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='diffuse albedo of impervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD_DIF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='diffuse albedo of impervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD_DIR', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='direct albedo of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD_DIR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='direct albedo of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD_DIF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='diffuse albedo of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD_DIF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='diffuse albedo of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF_DIR', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='direct albedo of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF_DIR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad',&
               long_name='direct albedo of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF_DIF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='diffuse albedo of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF_DIF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad',&
               long_name='diffuse albedo of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_WALL_DIR', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='direct albedo of wall', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_WALL_DIR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='direct albedo of wall', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_WALL_DIF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='numrad', &
               long_name='diffuse albedo of wall', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_WALL_DIF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='numrad', &
               long_name='diffuse albedo of wall', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='thermal conductivity of roof', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='thermal conductivity of roof', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='thermal conductivity of wall', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='thermal conductivity of wall', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='thermal conductivity of impervious road', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='thermal conductivity of impervious road', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='volumetric heat capacity of roof', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='volumetric heat capacity of roof', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='volumetric heat capacity of wall', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='volumetric heat capacity of wall', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numurbl', dim3name='nlevurb', &
               long_name='volumetric heat capacity of impervious road', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', dim4name='nlevurb', &
               long_name='volumetric heat capacity of impervious road', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='NLEV_IMPROAD', xtype=nf_int, &
               dim1name='gridcell', dim2name='numurbl', &
               long_name='number of impervious road layers', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='NLEV_IMPROAD', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
               long_name='number of impervious road layers', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='peatf', xtype=xtype, &
               dim1name='gridcell',&
               long_name='peatland fraction', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='peatf', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='peatland fraction', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='abm', xtype=nf_int, &
               dim1name='gridcell',&
               long_name='agricultural fire peak month', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='abm', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='agricultural fire peak month', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='gdp', xtype=xtype, &
               dim1name='gridcell',&
               long_name='gdp', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='gdp', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='gdp', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='SLOPE', xtype=xtype, &
               dim1name='gridcell',&
               long_name='mean topographic slope', units='degrees')
       else
          call ncd_defvar(ncid=ncid, varname='SLOPE', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='mean topographic slope', units='degrees')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='STD_ELEV', xtype=xtype, &
               dim1name='gridcell',&
               long_name='standard deviation of elevation', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='STD_ELEV', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='standard deviation of elevation', units='m')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='binfl', xtype=xtype, &
               dim1name='gridcell',&
               long_name='VIC b parameter for the Variable Infiltration Capacity Curve', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='binfl', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='VIC b parameter for the Variable Infiltration Capacity Curve', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='Ws', xtype=xtype, &
               dim1name='gridcell',&
               long_name='VIC Ws parameter for the ARNO curve', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='Ws', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='VIC Ws parameter for the ARNO Curve', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='Dsmax', xtype=xtype, &
               dim1name='gridcell',&
               long_name='VIC Dsmax parameter for the ARNO curve', units='mm/day')
       else
          call ncd_defvar(ncid=ncid, varname='Dsmax', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='VIC Dsmax parameter for the ARNO curve', units='mm/day')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='Ds', xtype=xtype, &
               dim1name='gridcell',&
               long_name='VIC Ds parameter for the ARNO curve', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='Ds', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='VIC Ds parameter for the ARNO curve', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='LAKEDEPTH', xtype=xtype, &
               dim1name='gridcell',&
               long_name='lake depth', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='LAKEDEPTH', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='lake depth', units='m')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='F0', xtype=xtype, &
               dim1name='gridcell',&
               long_name='maximum gridcell fractional inundated area', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='F0', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='maximum gridcell fractional inundated area', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='P3', xtype=xtype, &
               dim1name='gridcell',&
               long_name='coefficient for qflx_surf_lag for finundated', units='s/mm')
       else
          call ncd_defvar(ncid=ncid, varname='P3', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='coefficient for qflx_surf_lag for finundated', units='s/mm')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ZWT0', xtype=xtype, &
               dim1name='gridcell',&
               long_name='decay factor for finundated', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='ZWT0', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='decay factor for finundated', units='m')
       end if

    endif

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PCT_WETLAND', xtype=xtype, &
            dim1name='gridcell',&
            long_name='percent wetland', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_WETLAND', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent wetland', units='unitless')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PCT_LAKE', xtype=xtype, &
            dim1name='gridcell',&
            long_name='percent lake', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_LAKE', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent lake', units='unitless')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PCT_GLACIER', xtype=xtype, &
            dim1name='gridcell',&
            long_name='percent glacier', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_GLACIER', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent glacier', units='unitless')
    end if

    if (.not. dynlanduse) then
       if ( nglcec > 0 )then
          call ncd_defvar(ncid=ncid, varname='GLC_MEC', xtype=xtype, &
               dim1name='nglcecp1', long_name='Glacier elevation class', units='m')
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC', xtype=xtype, &
                  dim1name='gridcell', dim2name='nglcec', &
                  long_name='percent glacier for each glacier elevation class (% of landunit)', units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
                  long_name='percent glacier for each glacier elevation class (% of landunit)', units='unitless')
          end if
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
                  dim1name='gridcell', dim2name='nglcec', &
                  long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC_GIC', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
                  long_name='percent smaller glaciers and ice caps for each glacier elevation class (% of landunit)', units='unitless')
          end if
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
                  dim1name='gridcell', dim2name='nglcec', &
                  long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC_ICESHEET', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
                  long_name='percent ice sheet for each glacier elevation class (% of landunit)', units='unitless')
          end if
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_GIC', xtype=xtype, &
                  dim1name='gridcell', &
                  long_name='percent ice caps/glaciers (% of landunit)', units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_GIC', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', &
                  long_name='percent ice caps/glaciers (% of landunit)', units='unitless')
          end if
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
                  dim1name='gridcell', &
                  long_name='percent ice sheet (% of landunit)', units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname='PCT_GLC_ICESHEET', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', &
                  long_name='percent ice sheet (% of landunit)', units='unitless')
          end if
          
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname='TOPO_GLC_MEC', xtype=xtype, &
                  dim1name='gridcell', dim2name='nglcec', &
                  long_name='mean elevation on glacier elevation classes', units='m')
          else
             call ncd_defvar(ncid=ncid, varname='TOPO_GLC_MEC', xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
                  long_name='mean elevation on glacier elevation classes', units='m')
          end if
          
       end if
    end if

    if (.not. dynlanduse) then
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TOPO', xtype=xtype, &
               dim1name='gridcell', &
               long_name='mean elevation on land', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='TOPO', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='mean elevation on land', units='m')
       end if
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=xtype, &
            dim1name='gridcell', dim2name='numurbl', &
            long_name='percent urban for each density type', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', dim3name='numurbl', &
            long_name='percent urban for each density type', units='unitless')
    end if

    if (.not. dynlanduse) then

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='URBAN_REGION_ID', xtype=nf_int, &
               dim1name='gridcell',&
               long_name='urban region ID', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='URBAN_REGION_ID', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='urban region ID', units='unitless')
       end if

    end if

    if (dynlanduse) then
       do j = 1, mkharvest_numtypes()
          if (outnc_1d) then
             call ncd_defvar(ncid=ncid, varname=mkharvest_fieldname(j), xtype=xtype, &
                  dim1name='gridcell', dim2name='time', &
                  long_name=mkharvest_longname(j), units='unitless')
          else
             call ncd_defvar(ncid=ncid, varname=mkharvest_fieldname(j), xtype=xtype, &
                  dim1name='lsmlon', dim2name='lsmlat', dim3name='time', &
                  long_name=mkharvest_longname(j), units='unitless')
          end if
       end do
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='APATITE_P', xtype=xtype, &
            dim1name='gridcell',&
            long_name='Apatite Phosphorus ', units='gP/m2')
    else
       call ncd_defvar(ncid=ncid, varname='APATITE_P', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='Apatite Phosphorus', units='gP/m2')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LABILE_P', xtype=xtype, &
            dim1name='gridcell',&
            long_name='Labile Inorganic Phosphorus', units='gP/m2')
    else
       call ncd_defvar(ncid=ncid, varname='LABILE_P', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='Labile Inorganic Phosphorus', units='gP/m2')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='OCCLUDED_P', xtype=xtype, &
            dim1name='gridcell',&
            long_name='Occluded Phosphorus', units='gP/m2')
    else
       call ncd_defvar(ncid=ncid, varname='OCCLUDED_P', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='Occluded Phosphorus', units='gP/m2')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='SECONDARY_P', xtype=xtype, &
            dim1name='gridcell',&
            long_name='Secondary Mineral Phosphorus', units='gP/m2')
    else
       call ncd_defvar(ncid=ncid, varname='SECONDARY_P', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='Secondary Mineral Phosphorus', units='gP/m2')
    end if

    ! End of define mode

    call check_ret(nf_enddef(ncid), subname)
    call check_ret(nf_close(ncid), subname)

  end subroutine mkfile

end module mkfileMod
