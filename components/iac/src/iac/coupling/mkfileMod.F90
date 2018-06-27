module mkfileMod

contains

!-----------------------------------------------------------------------
  subroutine mkfile(domain, fname, dynlanduse)

    use shr_kind_mod, only : r8 => shr_kind_r8
    use shr_sys_mod , only : shr_sys_getenv
    use mkfileutils , only : get_filename
    use mkvarpar    , only : nlevsoi, numpft, nlevurb, numsolar, numrad
    use mkvarctl
    use mkglcmecMod , only : nglcec
    use mkharvestMod, only : mkharvest_fieldname, mkharvest_numtypes, mkharvest_longname
    use mkncdio     , only : check_ret, ncd_defvar, ncd_convl2i
    use mkdomainMod  

    implicit none
    include 'netcdf.inc'
    type(domain1_type), intent(in) :: domain
    character(len=*)  , intent(in) :: fname
    logical           , intent(in) :: dynlanduse

    integer :: ncid
    integer :: j                    ! index
    integer :: pftsize              ! size of lsmpft dimension
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
    logical            :: outnc_1d
    character(len=32)  :: subname = 'mkfile'  ! subroutine name
!-----------------------------------------------------------------------

    outnc_1d = .false.
    if ((domain%ni == -9999 .and. domain%nj == -9999) .or. outnc_dims ==1) then
       outnc_1d = .true.
    end if

    if ( .not. outnc_large_files )then
       call check_ret(nf_create(trim(fname), nf_clobber, ncid), subname)
    else
       call check_ret(nf_create(trim(fname), ior(nf_clobber,nf_64bit_offset), &
                                ncid), subname)
    end if
    call check_ret(nf_set_fill (ncid, nf_nofill, omode), subname)

    ! Define dimensions.

    pftsize = numpft + 1

    if (outnc_1d) then
       call check_ret(nf_def_dim (ncid, 'gridcell', domain%ns, dimid), subname)
    else
       call check_ret(nf_def_dim (ncid, 'lsmlon'  , domain%ni, dimid), subname)
       call check_ret(nf_def_dim (ncid, 'lsmlat'  , domain%nj, dimid), subname)
    end if

    if (.not. dynlanduse) then
       call check_ret(nf_def_dim (ncid, 'nlevsoi' , nlevsoi     , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'nglcec'  , nglcec      , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'nglcecp1', nglcec+1    , dimid), subname)
    end if
    call check_ret(nf_def_dim (ncid, 'nlevurb' , nlevurb     , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'numsolar', numsolar    , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'numrad'  , numrad      , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'lsmpft'  , pftsize     , dimid), subname)
    call check_ret(nf_def_dim (ncid, 'time'    , nf_unlimited, dimid), subname)
    call check_ret(nf_def_dim (ncid, 'nchar'   , 128         , dimid), subname)

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

#ifdef OPT
    str = 'TRUE'
#else
    str = 'FALSE'
#endif
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Compiler_Optimized', len_trim(str), trim(str)), subname)

    str = &
'$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/clm4_0_22/models/lnd/clm/tools/mksurfdata/mkfileMod.F90 $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Version', len_trim(str), trim(str)), subname)

    str = '$Id: mkfileMod.F90 25173 2010-10-15 23:39:33Z erik $'
    call check_ret(nf_put_att_text (ncid, NF_GLOBAL, &
         'Revision_Id', len_trim(str), trim(str)), subname)

    !call check_ret(nf_put_att_int1(ncid, NF_GLOBAL, &
    !     'all_urban', ncd_convl2i(all_urban)), subname)

    str = get_filename(mksrf_fgrid)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Input_grid_dataset', len_trim(str), trim(str)), subname)

    str = trim(mksrf_gridtype)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Input_gridtype', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fvegtyp)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Vegetation_type_raw_data_filename', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_fsoitex)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_texture_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fsoicol)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Soil_color_raw_data_file_name', len_trim(str), trim(str)), subname)

       str = get_filename(mksrf_fvocef)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'VOC_EF_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if

    str = get_filename(mksrf_forganic)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Organic_matter_raw_data_file_name', len_trim(str), trim(str)), subname)
       
    str = get_filename(mksrf_flanwat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Inland_water_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fglacier)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Glacier_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fglctopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'GLC_Topography_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_flndtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Land_Topography_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_ffrac)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Fracdata_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_fmax)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Fmax_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_furban)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Urban_raw_data_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(mksrf_firrig)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'Irrig_raw_data_file_name', len_trim(str), trim(str)), subname)

    if (.not. dynlanduse) then
       str = get_filename(mksrf_flai)
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
            'Lai_raw_data_file_name', len_trim(str), trim(str)), subname)
    end if

    str = get_filename(map_fpft)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_pft_file_name', len_trim(str), trim(str)), subname)

    str = get_filename(map_flanwat)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_lanwat_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fglacier)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_glacier_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fsoitex)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_texture_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fsoicol)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_soil_color_file', len_trim(str), trim(str)), subname)

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

    str = get_filename(map_firrig)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_irrigation_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_flai)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_lai_sai_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fglctopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_glc_topography_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_flndtopo)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_land_topography_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fglcmec_t2g)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_topo_to_raw_glacier_file', len_trim(str), trim(str)), subname)

    str = get_filename(map_fglcmec_g2g)
    call check_ret(nf_put_att_text(ncid, NF_GLOBAL, &
         'map_raw_glacier_to_model_glacier_file', len_trim(str), trim(str)), subname)

    ! ----------------------------------------------------------------------
    ! Define variables
    ! ----------------------------------------------------------------------
 
    if ( .not. outnc_double )then
       xtype = nf_float
    else
       xtype = nf_double
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LATN' , xtype=nf_double, &
            dim1name='gridcell', &
            long_name='latitude of north edge', units='degrees north')
    else
       call ncd_defvar(ncid=ncid, varname='LATN' , xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='latitude of north edge', units='degrees north')
    end if
       
    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LONE' , xtype=nf_double, &
            dim1name='gridcell', &
            long_name='longitude of east edge', units='degrees east')
    else
       call ncd_defvar(ncid=ncid, varname='LONE' , xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='longitude of east edge', units='degrees east')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LATS' , xtype=nf_double, &
            dim1name='gridcell',&
            long_name='latitude of south edge', units='degrees north')
    else
       call ncd_defvar(ncid=ncid, varname='LATS' , xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='latitude of south edge', units='degrees north')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LONW' , xtype=nf_double, &
            dim1name='gridcell',&
            long_name='longitude of west edge', units='degrees east')
    else
       call ncd_defvar(ncid=ncid, varname='LONW' , xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='longitude of west edge', units='degrees east')
    end if

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

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
            dim1name='gridcell',&
            long_name='land fraction from pft dataset', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='LANDFRAC_PFT', xtype=nf_double, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='land fraction from pft dataset', units='unitless')
    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
            dim1name='gridcell',&
            long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PFTDATA_MASK', xtype=nf_int, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='land mask from pft dataset, indicative of real/fake points', units='unitless')
    end if

    if (.not. dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='mxsoil_color', xtype=nf_int, &
            long_name='maximum numbers of soil colors', units='unitless')

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
               dim1name='gridcell',&
               long_name='soil color', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='SOIL_COLOR', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='soil color', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevsoi', &
               long_name='percent sand', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_SAND', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
               long_name='percent sand', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevsoi', &
               long_name='percent clay', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_CLAY', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
               long_name='percent clay', units='unitless')
       end if

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
          call ncd_defvar(ncid=ncid, varname='ORGANIC', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevsoi', &
               long_name='organic matter density at soil levels', &
               units='kg/m3 (assumed carbon content 0.58 gC per gOM)')
       else
          call ncd_defvar(ncid=ncid, varname='ORGANIC', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevsoi', &
               long_name='organic matter density at soil levels', &
               units='kg/m3 (assumed carbon content 0.58 gC per gOM)')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CANYON_HWR', xtype=xtype, &
               dim1name='gridcell',&
               long_name='canyon height to width ratio', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='CANYON_HWR', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='canyon height to width ratio', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_IMPROAD', xtype=xtype, &
               dim1name='gridcell',&
               long_name='emissivity of impervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='emissivity of impervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_PERROAD', xtype=xtype, &
               dim1name='gridcell',&
               long_name='emissivity of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_PERROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='emissivity of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_ROOF', xtype=xtype, &
               dim1name='gridcell',&
               long_name='emissivity of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='emissivity of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='EM_WALL', xtype=xtype, &
               dim1name='gridcell',&
               long_name='emissivity of wall', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='EM_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='emissivity of wall', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='HT_ROOF', xtype=xtype, &
               dim1name='gridcell',&
               long_name='height of roof', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='HT_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='height of roof', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='THICK_ROOF', xtype=xtype, &
               dim1name='gridcell',&
               long_name='thickness of roof', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='THICK_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='thickness of roof', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='THICK_WALL', xtype=xtype, &
               dim1name='gridcell',&
               long_name='thickness of wall', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='THICK_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='thickness of wall', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MAX', xtype=xtype, &
               dim1name='gridcell',&
               long_name='maximum interior building temperature', units='K')
       else
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MAX', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='maximum interior building temperature', units='K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MIN', xtype=xtype, &
               dim1name='gridcell',&
               long_name='minimum interior building temperature', units='K')
       else
          call ncd_defvar(ncid=ncid, varname='T_BUILDING_MIN', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='minimum interior building temperature', units='K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WIND_HGT_CANYON', xtype=xtype, &
               dim1name='gridcell',&
               long_name='height of wind in canyon', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='WIND_HGT_CANYON', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='height of wind in canyon', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WTLUNIT_ROOF', xtype=xtype, &
               dim1name='gridcell',&
               long_name='fraction of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='WTLUNIT_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='fraction of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='WTROAD_PERV', xtype=xtype, &
               dim1name='gridcell',&
               long_name='fraction of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='WTROAD_PERV', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='fraction of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numrad', dim3name='numsolar', &
               long_name='albedo of impervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numrad', dim4name='numsolar', &
               long_name='albedo of impervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='numrad', dim3name='numsolar', &
               long_name='albedo of pervious road', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_PERROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numrad', dim4name='numsolar', &
               long_name='albedo of pervious road', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='numrad', dim3name='numsolar', &
               long_name='albedo of roof', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numrad',dim4name='numsolar', &
               long_name='albedo of roof', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='ALB_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='numrad', dim3name='numsolar', &
               long_name='albedo of wall', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='ALB_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='numrad', dim4name='numsolar', &
               long_name='albedo of wall', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='thermal conductivity of roof', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='thermal conductivity of roof', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='thermal conductivity of wall', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='thermal conductivity of wall', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TK_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='thermal conductivity of impervious road', units='W/m*K')
       else
          call ncd_defvar(ncid=ncid, varname='TK_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='thermal conductivity of impervious road', units='W/m*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_ROOF', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='volumetric heat capacity of roof', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_ROOF', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='volumetric heat capacity of roof', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_WALL', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='volumetric heat capacity of wall', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_WALL', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='volumetric heat capacity of wall', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='CV_IMPROAD', xtype=xtype, &
               dim1name='gridcell', dim2name='nlevurb', &
               long_name='volumetric heat capacity of impervious road', units='J/m^3*K')
       else
          call ncd_defvar(ncid=ncid, varname='CV_IMPROAD', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nlevurb', &
               long_name='volumetric heat capacity of impervious road', units='J/m^3*K')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='NLEV_IMPROAD', xtype=nf_int, &
               dim1name='gridcell',&
               long_name='number of impervious road layers', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='NLEV_IMPROAD', xtype=nf_int, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='number of impervious road layers', units='unitless')
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
       call ncd_defvar(ncid=ncid, varname='GLC_MEC', xtype=xtype, &
            dim1name='nglcecp1', long_name='Glacier elevation class', units='m')
   
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC', xtype=xtype, &
               dim1name='gridcell', dim2name='nglcec', &
               long_name='percent for each glacier elevation class', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_GLC_MEC', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
               long_name='percent for each glacier elevation class', units='unitless')
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

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='TOPO', xtype=xtype, &
               dim1name='gridcell', &
               long_name='mean elevation on land', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='TOPO', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='mean elevation on land', units='m')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='THCK_GLC_MEC', xtype=xtype, &
               dim1name='gridcell', dim2name='nglcec', &
               long_name='mean ice sheet thickness on glacier elevation classes', units='m')
       else
          call ncd_defvar(ncid=ncid, varname='THCK_GLC_MEC', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='nglcec', &
               long_name='mean ice sheet thickness on glacier elevation classes', units='m')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='FMAX', xtype=xtype, &
               dim1name='gridcell', &
               long_name='maximum fractional saturated area', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='FMAX', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='maximum fractional saturated area', units='unitless')
       end if

    end if

    if (outnc_1d) then
       call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=xtype, &
            dim1name='gridcell',&
            long_name='percent urban', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='PCT_URBAN', xtype=xtype, &
            dim1name='lsmlon', dim2name='lsmlat', &
            long_name='percent urban', units='unitless')
    end if

    if (mksrf_firrig /= ' ') then
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_IRRIG', xtype=xtype, &
               dim1name='gridcell',&
               long_name='percent irrigated area', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_IRRIG', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', &
               long_name='percent irrigated area', units='unitless')
       end if
    end if

    if (.not. dynlanduse) then
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
               dim1name='gridcell', dim2name='lsmpft', &
               long_name='percent plant functional type of gridcell', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', &
               long_name='percent plant functional type of gridcell', units='unitless')
       end if
    else
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
               dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
               long_name='percent plant functional type of gridcell', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='PCT_PFT', xtype=xtype, &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
               long_name='percent plant functional type of gridcell', units='unitless')
       end if

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

    if (.not. dynlanduse) then
       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=xtype,  &
               dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
               long_name='monthly leaf area index', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='MONTHLY_LAI', xtype=xtype,  &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
               long_name='monthly leaf area index', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=xtype,  &
               dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
               long_name='monthly stem area index', units='unitless')
       else
          call ncd_defvar(ncid=ncid, varname='MONTHLY_SAI', xtype=xtype,  &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
               long_name='monthly stem area index', units='unitless')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
               dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
               long_name='monthly height top', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_TOP', xtype=xtype,  &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
               long_name='monthly height top', units='meters')
       end if

       if (outnc_1d) then
          call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
               dim1name='gridcell', dim2name='lsmpft', dim3name='time', &
               long_name='monthly height bottom', units='meters')
       else
          call ncd_defvar(ncid=ncid, varname='MONTHLY_HEIGHT_BOT', xtype=xtype,  &
               dim1name='lsmlon', dim2name='lsmlat', dim3name='lsmpft', dim4name='time', &
               long_name='monthly height bottom', units='meters')
       end if

    end if
       
    if (dynlanduse) then
       call ncd_defvar(ncid=ncid, varname='YEAR', xtype=nf_int,  &
            dim1name='time', &
            long_name='Year of PFT data', units='unitless')
       call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='year', units='unitless')
       call ncd_defvar(ncid=ncid, varname='input_pftdata_filename', xtype=nf_char,  &
            dim1name='nchar', dim2name='time',  &
            long_name='Input filepath for PFT values for this year', units='unitless')
    else
       call ncd_defvar(ncid=ncid, varname='time', xtype=nf_int,  &
            dim1name='time', &
            long_name='Calendar month', units='month')
    end if

    ! End of define mode

    call check_ret(nf_enddef(ncid), subname)
    call check_ret(nf_close(ncid), subname)

  end subroutine mkfile

end module mkfileMod
