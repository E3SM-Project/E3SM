module initInterpMod

  !----------------------------------------------------------------------- 
  ! Interpolate initial conditions file from one resolution and/or landmask
  ! to another resolution and/or landmask
  !----------------------------------------------------------------------- 

  use shr_kind_mod   , only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_const_mod  , only: SHR_CONST_PI, SHR_CONST_REARTH
  use shr_sys_mod    , only: shr_sys_flush
  use shr_infnan_mod , only: shr_infnan_isnan
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use elm_varctl     , only: iulog
  use abortutils     , only: endrun
  use spmdMod        , only: masterproc
  use restUtilMod    , only: iflag_interp, iflag_copy, iflag_skip
  use restUtilMod    , only: iflag_noswitchdim, iflag_switchdim
  use elm_varcon     , only: spval, re
  use ncdio_pio
  use pio

  implicit none
  private
  save

  ! Public methods

  public :: initInterp

  ! Private methods

  private :: check_dim_subgrid
  private :: check_dim_level
  private :: findMinDist
  private :: set_subgrid_dist
  private :: set_subgrid_glob
  private :: set_mindist
  private :: is_sametype
  private :: is_baresoil
  private :: interp_0d_copy
  private :: interp_1d_double
  private :: interp_1d_int
  private :: interp_2d_double

  ! Private data
 
  integer          :: ipft_not_vegetated 
  integer          :: icol_vegetated_or_bare_soil
  integer          :: ilun_vegetated_or_bare_soil
  integer          :: ilun_landice_multiple_elevation_classes 
  character(len=8) :: created_glacier_mec_landunits
  logical          :: override_missing = .true. ! override missing types with closest bare-soil
  
  type, public :: subgrid_type
     character(len=16) :: name               ! pft, column, landunit
     integer , pointer :: ptype(:) => null() ! used for pft type 
     integer , pointer :: ctype(:) => null() ! used for pft or col type
     integer , pointer :: ltype(:) => null() ! used for pft, col or lun type
     real(r8), pointer :: topoglc(:) => null()
     real(r8), pointer :: lat(:)
     real(r8), pointer :: lon(:)
     real(r8), pointer :: coslat(:)
  end type subgrid_type

contains

  !=======================================================================

  subroutine initInterp (filei, fileo, bounds)

    !----------------------------------------------------------------------- 
    ! Read initial data from netCDF instantaneous initial data history file 
    !-----------------------------------------------------------------------

    use decompMod, only: bounds_type

    ! --------------------------------------------------------------------
    ! arguments
    character(len=*)  , intent(in) :: filei     !input  initial dataset
    character(len=*)  , intent(in) :: fileo    !output initial dataset
    type(bounds_type) , intent(in) :: bounds
    !
    ! local variables
    integer            :: i,j,k,l,m,n     ! loop indices    
    integer            :: begi, endi      ! beginning/ending indices 
    integer            :: bego, endo      ! beginning/ending indices 
    integer            :: begp_i, endp_i  ! input file pft bounds
    integer            :: begp_o, endp_o  ! output file pft bounds
    integer            :: begc_i, endc_i  ! input file column bounds
    integer            :: begc_o, endc_o  ! output file column bounds
    integer            :: begl_i, endl_i  ! input file landunit bounds
    integer            :: begl_o, endl_o  ! output file landunit bounds
    integer            :: begg_i, endg_i  ! input file gridcell bounds
    integer            :: begg_o, endg_o  ! output file gridcell bounds
    integer            :: nlevi,nlevo     ! input/output number of levels
    type(file_desc_t)  :: ncidi, ncido    ! input/output pio fileids 
    integer            :: dimleni,dimleno ! input/output dimension length       
    integer            :: nvars           ! number of variables
    character(len=256) :: varname         ! variable name
    character(len=16)  :: vec_dimname     ! subgrid dimension name
    character(len=16)  :: lev_dimname     ! level dimension name
    type(Var_desc_t)   :: vardesc         ! pio variable descriptor
    integer            :: levdimi,levdimo ! input/output level dimension (2d fields only)
    integer            :: vecdimi,vecdimo ! input/output non-level dimension (2d fields only)
    integer            :: xtypeo          ! netCDF variable type
    integer            :: varido          ! netCDF variable id
    integer            :: ndimso          ! netCDF number of dimensions
    integer            :: dimidso(3) = -1 ! netCDF dimension ids
    integer            :: dimidsi(3) = -1 ! netCDF dimension ids
    integer            :: status          ! return code
    integer            :: switchdim_flag  ! 1 => level dimension is first dimension
    integer            :: iflag_interpinic
    logical            :: switchdimi, switchdimo
    real(r8)           :: rvalue
    integer            :: ivalue 
    integer            :: spinup_state_i, spinup_state_o
    integer            :: decomp_cascade_state_i, decomp_cascade_state_o 
    integer            :: npftsi, ncolsi, nlunsi, ngrcsi 
    integer            :: npftso, ncolso, nlunso, ngrcso 
    integer , pointer  :: pftindx(:)     
    integer , pointer  :: colindx(:)     
    integer , pointer  :: lunindx(:)     
    integer , pointer  :: grcindx(:) 
    logical , pointer  :: pft_activei(:), pft_activeo(:) 
    logical , pointer  :: col_activei(:), col_activeo(:) 
    logical , pointer  :: lun_activei(:), lun_activeo(:) 
    logical , pointer  :: grc_activei(:), grc_activeo(:) 
    integer , pointer  :: sgridindex(:)
    logical , pointer  :: activei(:), activeo(:)
    !--------------------------------------------------------------------

    if (masterproc) then
       write (iulog,*) '**** Mapping clm initial data from input ',trim(filei),&
            '  to output ',trim(fileo),' ****'
    end if

    ! --------------------------------------------
    ! Open input and output initial conditions files (both just for reading now)
    ! --------------------------------------------

    call ncd_pio_openfile (ncidi, trim(filei) , 0)
    call ncd_pio_openfile (ncido, trim(fileo),  ncd_write)

    ! --------------------------------------------
    ! Determine dimensions and error checks on dimensions
    ! --------------------------------------------

    call check_dim_subgrid(ncidi, ncido, dimname ='pft'     , dimleni=npftsi, dimleno=npftso)
    call check_dim_subgrid(ncidi, ncido, dimname ='column'  , dimleni=ncolsi, dimleno=ncolso)
    call check_dim_subgrid(ncidi, ncido, dimname ='landunit', dimleni=nlunsi, dimleno=nlunso)
    call check_dim_subgrid(ncidi, ncido, dimname ='gridcell', dimleni=ngrcsi, dimleno=ngrcso)

    if (masterproc) then
       write (iulog,*) 'input gridcells = ',ngrcsi,' output gridcells = ',ngrcso
       write (iulog,*) 'input landuntis = ',nlunsi,' output landunits = ',nlunso
       write (iulog,*) 'input columns   = ',ncolsi,' output columns   = ',ncolso
       write (iulog,*) 'input pfts      = ',npftsi,' output pfts      = ',npftso
    end if

    call check_dim_level(ncidi, ncido, dimname='levsno' )
    call check_dim_level(ncidi, ncido, dimname='levsno1')
    call check_dim_level(ncidi, ncido, dimname='levcan' )
    call check_dim_level(ncidi, ncido, dimname='levlak' )
    call check_dim_level(ncidi, ncido, dimname='levtot' )
    call check_dim_level(ncidi, ncido, dimname='levgrnd')
    call check_dim_level(ncidi, ncido, dimname='numrad' )

    ! --------------------------------------------
    ! Determine input file global attributes that are needed 
    ! --------------------------------------------

    status = pio_get_att(ncidi, pio_global, &
         'ipft_not_vegetated',                      ipft_not_vegetated)
    status = pio_get_att(ncidi, pio_global, &
         'icol_vegetated_or_bare_soil',             icol_vegetated_or_bare_soil)
    status = pio_get_att(ncidi, pio_global, &
         'ilun_vegetated_or_bare_soil',             ilun_vegetated_or_bare_soil)
    status = pio_get_att(ncidi, pio_global, &
         'ilun_landice_multiple_elevation_classes', ilun_landice_multiple_elevation_classes)
    status = pio_get_att(ncidi, pio_global, &
         'created_glacier_mec_landunits',           created_glacier_mec_landunits)

    if (masterproc) then
       write(iulog,*)'ipft_not_vegetated                      = ' ,ipft_not_vegetated
       write(iulog,*)'icol_vegetated_or_bare_soil             = ' ,icol_vegetated_or_bare_soil
       write(iulog,*)'ilun_vegetated_or_bare_soil             = ' ,ilun_vegetated_or_bare_soil
       write(iulog,*)'ilun_landice_multiple_elevation_classes = ' ,ilun_landice_multiple_elevation_classes
       write(iulog,*)'create_glacier_mec_landunits            = ',trim(created_glacier_mec_landunits)
    end if

    ! --------------------------------------------
    ! Find closest values for pfts, cols, landunits
    ! --------------------------------------------

    begp_i = 1 ;  endp_i = npftsi
    begc_i = 1 ;  endc_i = ncolsi
    begl_i = 1 ;  endl_i = nlunsi
    begg_i = 1 ;  endg_i = ngrcsi

    begp_o = bounds%begp ;  endp_o = bounds%endp
    begc_o = bounds%begc ;  endc_o = bounds%endc
    begl_o = bounds%begl ;  endl_o = bounds%endl
    begg_o = bounds%begg ;  endg_o = bounds%endg

    allocate(pft_activei(begp_i:endp_i))
    allocate(col_activei(begc_i:endc_i))
    allocate(lun_activei(begl_i:endl_i))
    allocate(grc_activei(begg_i:endg_i))

    allocate(pft_activeo(begp_o:endp_o))
    allocate(col_activeo(begc_o:endc_o))
    allocate(lun_activeo(begl_o:endl_o))
    allocate(grc_activeo(begg_o:endg_o))

    allocate(pftindx(begp_o:endp_o))
    allocate(colindx(begc_o:endc_o))
    allocate(lunindx(begl_o:endl_o))
    allocate(grcindx(begg_o:endg_o))

    ! For each output pft, find the input pft, pftindx, that is closest

    if (masterproc) then
       write(iulog,*)'finding minimum distance for pfts'
    end if
    vec_dimname = 'pft'
    call findMinDist(vec_dimname, begp_i, endp_i, begp_o, endp_o, ncidi, ncido, &
         pft_activei, pft_activeo, pftindx )

    ! For each output column, find the input column, colindx, that is closest

    if (masterproc) then
       write(iulog,*)'finding minimum distance for columns'
    end if
    vec_dimname = 'column'
    call findMinDist(vec_dimname, begc_i, endc_i, begc_o, endc_o, ncidi, ncido, &
         col_activei, col_activeo, colindx )

    ! For each output landunit, find the input landunit, lunindx, that is closest

    if (masterproc) then
       write(iulog,*)'finding minimum distance for landunits'
    end if
    vec_dimname = 'landunit'
    call findMinDist(vec_dimname, begl_i, endl_i, begl_o, endl_o, ncidi, ncido, &
         lun_activei, lun_activeo, lunindx )

    ! For each output gridcell, find the input gridcell, grcindx, that is closest

    if (masterproc) then
       write(iulog,*)'finding minimum distance for gridcells'
    end if
    vec_dimname = 'gridcell'
    call findMinDist(vec_dimname, begg_i, endg_i, begg_o, endg_o, ncidi, ncido, &
         grc_activei, grc_activeo, grcindx )

    !------------------------------------------------------------------------          
    ! Read input initial data and write output initial data
    !------------------------------------------------------------------------          

    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    ! Read input initial data and write output initial data
    ! Only examing the snow interfaces above zi=0 => zisno and zsno have
    ! the same level dimension below

    if (masterproc) then
       write(iulog,*)'reading in initial dataset'
    end if
    ! Get number of output variables and loop over them 
    status = pio_inquire(ncido, nVariables=nvars)
    do varido = 1, nvars

       !---------------------------------------------------          
       ! Given varido, get out variable data 
       !---------------------------------------------------          

       status = pio_inquire_variable(ncido, varid=varido, name=varname, &
            xtype=xtypeo, ndims=ndimso, dimids=dimidso)

       !---------------------------------------------------          
       ! If variable is zsoi, SKIP this variable
       !---------------------------------------------------          

       if  ( trim(varname) == 'zsoi' ) then
          if (masterproc) then
             write(iulog,*) 'Skipping     : ',trim(varname)
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! If interpinic flag is set to skip on output file
       ! SKIP this variable
       !---------------------------------------------------          

       status = pio_inq_varid (ncido, trim(varname), vardesc)
       status = pio_get_att(ncido, vardesc, 'interpinic_flag', iflag_interpinic)
       if (iflag_interpinic == iflag_skip) then
          if (masterproc) then
             write (iulog,*) 'Skipping     : ', trim(varname)
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! If variable is on output file - but not on input file
       ! SKIP this variable
       !---------------------------------------------------          

       call pio_seterrorhandling(ncidi, PIO_BCAST_ERROR)
       status = pio_inq_varid(ncidi, name=varname, vardesc=vardesc)
       call pio_seterrorhandling(ncidi, PIO_INTERNAL_ERROR)
       if (status /= PIO_noerr) then
          if (masterproc) then
             write (iulog,*) 'Skipping     : ', trim(varname), ' variable is NOT on input file'
          end if
          CYCLE
       end if

       !---------------------------------------------------          
       ! For scalar outut variables 
       !---------------------------------------------------          

       if ( ndimso == 0 ) then

          if  ( trim(varname) .eq. 'spinup_state' ) then

             ! since we are copying soil variables, need to also copy spinup state 
             ! since otherwise if they are different then it will break the spinup procedure
             status = pio_inq_varid(ncidi, trim(varname), vardesc)
             status = pio_get_var(ncidi, vardesc, spinup_state_i)
             status = pio_inq_varid(ncido, trim(varname), vardesc)
             status = pio_get_var(ncido, vardesc, spinup_state_o)
             if ( spinup_state_i  /=  spinup_state_o ) then
                if (masterproc) then
                   write (iulog,*) 'Spinup states are different: Copying: ', trim(varname)
                end if
                status = pio_put_var(ncido, vardesc, spinup_state_i)
             else
                if (masterproc) then
                   write (iulog,*) 'Spinup states match: Skipping: ', trim(varname)
                end if
             endif
             
          else if  ( trim(varname) .eq. 'decomp_cascade_state' ) then

             ! ditto for the decomposition cascade
             status = pio_inq_varid(ncidi, trim(varname), vardesc)
             status = pio_get_var(ncidi, vardesc, decomp_cascade_state_i)
             status = pio_inq_varid(ncido, trim(varname), vardesc)
             status = pio_get_var(ncido, vardesc, decomp_cascade_state_o)
             if ( decomp_cascade_state_i  /=  decomp_cascade_state_o ) then
                call endrun(msg='ERROR: Decomposition cascade states are different'//errMsg(__FILE__, __LINE__))
             else
                if (masterproc) then
                   write (iulog,*) 'Decomposition cascade states match: Skipping: ', trim(varname)
                end if
             endif

          else if (iflag_interpinic == iflag_copy) then

             call interp_0d_copy(varname, xtypeo, ncidi, ncido)

          else if (iflag_interpinic == iflag_skip) then

             if (masterproc) then
                write(iulog,*) 'Skipping     : ',trim(varname)
             end if

          end if

       !---------------------------------------------------          
       ! For 1D output variables
       !---------------------------------------------------          

       else if ( ndimso == 1 ) then

          status = pio_inq_dimname(ncido, dimidso(1), vec_dimname)
          if ( vec_dimname == 'pft' )then
             begi = begp_i
             endi = endp_i
             bego = begp_o
             endo = endp_o
             activei => pft_activei
             activeo => pft_activeo
             sgridindex => pftindx
          else if ( vec_dimname  == 'column' )then
             begi = begc_i
             endi = endc_i
             bego = begc_o
             endo = endc_o
             activei => col_activei
             activeo => col_activeo
             sgridindex => colindx
          else if ( vec_dimname == 'landunit' )then
             begi = begl_i
             endi = endl_i
             bego = begl_o
             endo = endl_o
             activei => lun_activei
             activeo => lun_activeo
             sgridindex => lunindx
          else if ( vec_dimname == 'gridcell' )then
             begi = begg_i
             endi = endg_i
             bego = begg_o
             endo = endg_o
             activei => grc_activei
             activeo => grc_activeo
             sgridindex => grcindx
          else
             call endrun(msg='ERROR interpinic: 1D variable '//trim(varname)//&
                  'with unknown subgrid dimension: '//trim(vec_dimname)//&
                  errMsg(__FILE__, __LINE__))
          end if

          if ( xtypeo == pio_int )then
             call interp_1d_int ( varname, vec_dimname, begi, endi, bego, endo, &
                  ncidi, ncido, activei, activeo, sgridindex )
          else if ( xtypeo == pio_double )then                                 
             call interp_1d_double( varname, vec_dimname, begi, endi, bego, endo, &
                  ncidi, ncido, activei, activeo, sgridindex )
          else
             call endrun(msg='ERROR interpinic: 1D variable with unknown type: '//&
                  trim(varname)//errMsg(__FILE__, __LINE__))
          end if

       !---------------------------------------------------          
       ! For 2D output variables
       !---------------------------------------------------          

       else if ( ndimso == 2 )then

          if ( xtypeo /= pio_double )then
             call endrun(msg='ERROR interpinic: 2D variable with unknown type: '//&
                  trim(varname)//errMsg(__FILE__, __LINE__))
          end if
          ! Determine order of level and subgrid dimension in restart file
          status = pio_inq_varid (ncidi, trim(varname), vardesc)
          status = pio_inquire_variable(ncidi, vardesc=vardesc, dimids=dimidsi)
          status = pio_get_att(ncidi, vardesc, 'switchdim_flag', switchdim_flag)
          if (switchdim_flag > 0) then
             levdimi = dimidsi(1)  
             vecdimi = dimidsi(2)  
             switchdimi = .true.
          else
             levdimi = dimidsi(2) 
             vecdimi = dimidsi(1) 
             switchdimi = .false.
          end if
          status = pio_inq_dimlen(ncidi,levdimi,nlevi)

          status = pio_inq_varid (ncido, trim(varname), vardesc)
          status = pio_get_att(ncido, vardesc, 'switchdim_flag', switchdim_flag)
          if (switchdim_flag > 0) then
             levdimo = dimidso(1)  
             vecdimo = dimidso(2)  
             switchdimo = .true.
          else
             levdimo = dimidso(2) 
             vecdimo = dimidso(1) 
             switchdimo = .false.
          end if
          status = pio_inq_dimlen(ncido,levdimo,nlevo)
          status = pio_inq_dimname(ncido, vecdimo, vec_dimname)
          status = pio_inq_dimname(ncido, levdimo, lev_dimname)

          if ( vec_dimname == 'pft' )then
             begi = begp_i
             endi = endp_i
             bego = begp_o
             endo = endp_o
             activei => pft_activei
             activeo => pft_activeo
             sgridindex => pftindx
          else if ( vec_dimname  == 'column' )then
             begi = begc_i
             endi = endc_i
             bego = begc_o
             endo = endc_o
             activei => col_activei
             activeo => col_activeo
             sgridindex => colindx
          else if ( vec_dimname == 'landunit' )then
             begi = begl_i
             endi = endl_i
             bego = begl_o
             endo = endl_o
             activei => lun_activei
             activeo => lun_activeo
             sgridindex => lunindx
          else
             call endrun(msg='ERROR interpinic: 2D variable with unknown subgrid dimension: '//&
                  trim(varname)//errMsg(__FILE__, __LINE__))
          end if
          call interp_2d_double(varname, vec_dimname, lev_dimname, &
               begi, endi, bego, endo, nlevi, nlevo, switchdimi, switchdimo, &
               ncidi, ncido, activei, activeo, sgridindex)

       else

          call endrun(msg='ERROR interpinic: variable NOT scalar, 1D or 2D: '//&
               trim(varname)//errMsg(__FILE__, __LINE__))

       end if
       call shr_sys_flush(iulog)

    end do
    ! Close input and output files

    call pio_closefile(ncidi)
    call pio_closefile(ncido)

    if (masterproc) then
       write (iulog,*) ' Successfully created initial condition file mapped from input IC file'
    end if

  end subroutine initInterp

  !=======================================================================

  subroutine findMinDist( dimname, begi, endi, bego, endo, ncidi, ncido, &
       activei, activeo, minindx)

    ! --------------------------------------------------------------------
    !
    ! Find the PFT distances based on the column distances already calculated
    !
    ! arguments
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    logical           , intent(out)   :: activei(begi:endi)
    logical           , intent(out)   :: activeo(bego:endo)
    integer           , intent(out)   :: minindx(bego:endo)         
    !
    ! local variables 
    type(subgrid_type)   :: subgridi 
    type(subgrid_type)   :: subgrido 
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*)'calling set_subgrid_glob for ',trim(dimname)
    end if
    call set_subgrid_glob(begi, endi, dimname, ncidi, activei, subgridi)

    if (masterproc) then
       write(iulog,*)'calling set_subgrid_dist',trim(dimname)
    end if
    call set_subgrid_dist(bego, endo, dimname, ncido, activeo, subgrido)

    if (masterproc) then
       write(iulog,*)'calling set_mindist for ',trim(dimname)
    end if
    call set_mindist(begi, endi, bego, endo, activei, activeo, subgridi, subgrido, minindx)

    deallocate(subgridi%lat, subgridi%lon, subgridi%coslat)
    deallocate(subgrido%lat, subgrido%lon, subgrido%coslat)
    
  end subroutine findMinDist

 !=======================================================================

  subroutine set_subgrid_dist(beg, end, dimname, ncid, active, subgrid)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)    :: beg, end
    type(file_desc_t)  , intent(inout) :: ncid
    character(len=*)   , intent(inout) :: dimname
    logical            , intent(inout) :: active(beg:end)    
    type(subgrid_type) , intent(inout) :: subgrid
    !
    ! local variables
    integer              :: n
    integer, pointer     :: itemp(:) 
    real(r8), parameter  :: deg2rad  = SHR_CONST_PI/180._r8
    !-----------------------------------------------------------------------

    subgrid%name = dimname

    allocate(itemp(beg:end))
    allocate(subgrid%lat(beg:end), subgrid%lon(beg:end), subgrid%coslat(beg:end))
    if (dimname == 'pft') then
       allocate(subgrid%ptype(beg:end), subgrid%ctype(beg:end), subgrid%ltype(beg:end))
    else if (dimname == 'column') then
       allocate(subgrid%ctype(beg:end), subgrid%ltype(beg:end))
    else if (dimname == 'landunit') then
       allocate(subgrid%ltype(beg:end))
    end if

    ! determine if is_glcmec from global attributes
    if (trim(created_glacier_mec_landunits) == 'true') then
       if  (dimname == 'pft' .or. dimname == 'column') then
          allocate(subgrid%topoglc(beg:end))
       end if
    end if

    if (dimname == 'pft') then
       call ncd_io(ncid=ncid, varname='pfts1d_lon'    , flag='read', data=subgrid%lon  , dim1name='pft') 
       call ncd_io(ncid=ncid, varname='pfts1d_lat'    , flag='read', data=subgrid%lat  , dim1name='pft') 
       call ncd_io(ncid=ncid, varname='pfts1d_itypveg', flag='read', data=subgrid%ptype, dim1name='pft')
       call ncd_io(ncid=ncid, varname='pfts1d_itypcol', flag='read', data=subgrid%ctype, dim1name='pft')
       call ncd_io(ncid=ncid, varname='pfts1d_ityplun', flag='read', data=subgrid%ltype, dim1name='pft')
       call ncd_io(ncid=ncid, varname='pfts1d_active' , flag='read', data=itemp        , dim1name='pft')
       if (associated(subgrid%topoglc)) then
          call ncd_io(ncid=ncid, varname='pfts1d_topoglc', flag='read', data=subgrid%topoglc, dim1name='pft')
       end if
    else if (dimname == 'column') then
       call ncd_io(ncid=ncid, varname='cols1d_lon'    , flag='read', data=subgrid%lon  , dim1name='column') 
       call ncd_io(ncid=ncid, varname='cols1d_lat'    , flag='read', data=subgrid%lat  , dim1name='column')  
       call ncd_io(ncid=ncid, varname='cols1d_ityp'   , flag='read', data=subgrid%ctype, dim1name='column') 
       call ncd_io(ncid=ncid, varname='cols1d_ityplun', flag='read', data=subgrid%ltype, dim1name='column') 
       call ncd_io(ncid=ncid, varname='cols1d_active' , flag='read', data=itemp        , dim1name='column')
       if (associated(subgrid%topoglc)) then
          call ncd_io(ncid=ncid, varname='cols1d_topoglc', flag='read', data=subgrid%topoglc, dim1name='column') 
       end if
    else if (dimname == 'landunit') then
       call ncd_io(ncid=ncid, varname='land1d_lon'    , flag='read', data=subgrid%lon  , dim1name='landunit') 
       call ncd_io(ncid=ncid, varname='land1d_lat'    , flag='read', data=subgrid%lat  , dim1name='landunit') 
       call ncd_io(ncid=ncid, varname='land1d_ityplun', flag='read', data=subgrid%ltype, dim1name='landunit')
       call ncd_io(ncid=ncid, varname='land1d_active' , flag='read', data=itemp        , dim1name='landunit') 
    else if (dimname == 'gridcell') then
       call ncd_io(ncid=ncid, varname='grid1d_lon'    , flag='read', data=subgrid%lon  , dim1name='gridunit') 
       call ncd_io(ncid=ncid, varname='grid1d_lat'    , flag='read', data=subgrid%lat  , dim1name='gridunit') 
       itemp(:) = 1
    end if

    do n = beg,end
       if (itemp(n) > 0) then
          active(n) = .true.
       else
          active(n) = .false.
       end if
       subgrid%lat(n) = subgrid%lat(n)*deg2rad
       subgrid%lon(n) = subgrid%lon(n)*deg2rad
       subgrid%coslat(n) = cos(subgrid%lat(n))
    end do

    deallocate(itemp)

  end subroutine set_subgrid_dist

  !=======================================================================

  subroutine set_subgrid_glob(beg, end, dimname, ncid, active, subgrid)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)    :: beg, end
    type(file_desc_t)  , intent(inout) :: ncid
    character(len=*)   , intent(inout) :: dimname
    logical            , intent(out)   :: active(beg:end)    
    type(subgrid_type) , intent(inout) :: subgrid
    !
    ! local variables
    integer              :: n
    integer , pointer    :: itemp(:) 
    real(r8), parameter  :: deg2rad  = SHR_CONST_PI/180._r8
    !-----------------------------------------------------------------------

    subgrid%name = dimname

    allocate(itemp(beg:end))
    allocate(subgrid%lat(beg:end), subgrid%lon(beg:end), subgrid%coslat(beg:end))
    if (dimname == 'pft') then
       allocate(subgrid%ptype(beg:end), subgrid%ctype(beg:end), subgrid%ltype(beg:end))
    else if (dimname == 'column') then
       allocate(subgrid%ctype(beg:end), subgrid%ltype(beg:end))
    else if (dimname == 'landunit') then
       allocate(subgrid%ltype(beg:end))
    end if

    ! determine if is_glcmec from global attributes
    if (trim(created_glacier_mec_landunits) == 'true') then
       if  (dimname == 'pft' .or. dimname == 'column') then
          allocate(subgrid%topoglc(beg:end))
       end if
    end if

    if (dimname == 'pft') then
       call ncd_io(ncid=ncid, varname='pfts1d_lon'    , flag='read', data=subgrid%lon  ) 
       call ncd_io(ncid=ncid, varname='pfts1d_lat'    , flag='read', data=subgrid%lat  ) 
       call ncd_io(ncid=ncid, varname='pfts1d_itypveg', flag='read', data=subgrid%ptype)
       call ncd_io(ncid=ncid, varname='pfts1d_itypcol', flag='read', data=subgrid%ctype)
       call ncd_io(ncid=ncid, varname='pfts1d_ityplun', flag='read', data=subgrid%ltype)
       call ncd_io(ncid=ncid, varname='pfts1d_active' , flag='read', data=itemp)
       if (associated(subgrid%topoglc)) then
          call ncd_io(ncid=ncid, varname='pfts1d_topoglc', flag='read', data=subgrid%topoglc)
       end if
    else if (dimname == 'column') then
       call ncd_io(ncid=ncid, varname='cols1d_lon'    , flag='read', data=subgrid%lon) 
       call ncd_io(ncid=ncid, varname='cols1d_lat'    , flag='read', data=subgrid%lat) 
       call ncd_io(ncid=ncid, varname='cols1d_ityp'   , flag='read', data=subgrid%ctype)
       call ncd_io(ncid=ncid, varname='cols1d_ityplun', flag='read', data=subgrid%ltype)
       call ncd_io(ncid=ncid, varname='cols1d_active' , flag='read', data=itemp)
       if (associated(subgrid%topoglc)) then
          call ncd_io(ncid=ncid, varname='cols1d_topoglc', flag='read', data=subgrid%topoglc)
       end if
    else if (dimname == 'landunit') then
       call ncd_io(ncid=ncid, varname='land1d_lon'    , flag='read', data=subgrid%lon  ) 
       call ncd_io(ncid=ncid, varname='land1d_lat'    , flag='read', data=subgrid%lat  ) 
       call ncd_io(ncid=ncid, varname='land1d_ityplun', flag='read', data=subgrid%ltype)
       call ncd_io(ncid=ncid, varname='land1d_active' , flag='read', data=itemp) 
    else if (dimname == 'gridcell') then
       call ncd_io(ncid=ncid, varname='grid1d_lon'    , flag='read', data=subgrid%lon  )
       call ncd_io(ncid=ncid, varname='grid1d_lat'    , flag='read', data=subgrid%lat  )
    end if

    do n = beg,end
       if (itemp(n) > 0) then
          active(n) = .true.
       else
          active(n) = .false.
       end if
       subgrid%lat(n) = subgrid%lat(n)*deg2rad
       subgrid%lon(n) = subgrid%lon(n)*deg2rad
       subgrid%coslat(n) = cos(subgrid%lat(n))
    end do

    deallocate(itemp)

  end subroutine set_subgrid_glob

  !=======================================================================

  subroutine set_mindist(begi, endi, bego, endo, activei, activeo, subgridi, subgrido, mindist_index)

    ! --------------------------------------------------------------------
    ! arguments
    integer            , intent(in)  :: begi, endi 
    integer            , intent(in)  :: bego, endo 
    logical            , intent(in)  :: activei(begi:endi) 
    logical            , intent(in)  :: activeo(bego:endo) 
    type(subgrid_type) , intent(in)  :: subgridi
    type(subgrid_type) , intent(in)  :: subgrido
    integer            , intent(out) :: mindist_index(bego:endo) 
    !
    ! local variables
    real(r8) :: dx,dy
    real(r8) :: distmin,dist,hgtdiffmin,hgtdiff    
    integer  :: nsizei, nsizeo
    integer  :: ni,no,nmin,ier,n,noloc
    logical  :: closest
    ! --------------------------------------------------------------------

    mindist_index(bego:endo) = 0
    distmin = spval

!$OMP PARALLEL DO PRIVATE (ni,no,n,nmin,distmin,dx,dy,dist,closest,hgtdiffmin,hgtdiff)
    do no = bego,endo
       
       ! If output type is contained in input dataset ...
       if (activeo(no)) then 

          nmin    = 0
          distmin = spval
          hgtdiffmin = spval
          do ni = begi,endi
             if (activei(ni)) then
                if (is_sametype(ni, no, subgridi, subgrido)) then
                   dy = abs(subgrido%lat(no)-subgridi%lat(ni))*re
                   dx = abs(subgrido%lon(no)-subgridi%lon(ni))*re * &
                        0.5_r8*(subgrido%coslat(no)+subgridi%coslat(ni))
                   dist = dx*dx + dy*dy
                   if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                      hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                   end if
                   closest = .false.
                   if ( dist < distmin ) then
                      closest = .true.
                      distmin = dist
                      nmin = ni
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiffmin = hgtdiff
                      end if
                   end if
                   if (.not. closest) then
                      if (associated(subgridi%topoglc) .and. associated(subgrido%topoglc)) then
                         hgtdiff = abs(subgridi%topoglc(ni) - subgrido%topoglc(no))
                         if ((dist == distmin) .and. (hgtdiff < hgtdiffmin)) then
                            closest = .true.
                            hgtdiffmin = hgtdiff
                            distmin = dist
                            nmin = ni
                         end if
                      end if
                   end if
                end if
             end if
          end do
          
          ! If output type is not contained in input dataset, then use closest bare soil 
          if ( override_missing .and. distmin == spval) then
             do ni = begi, endi
                if (activei(ni)) then
                   if ( is_baresoil(ni, subgridi)) then
                      dy = abs(subgrido%lat(no)-subgridi%lat(ni))*re
                      dx = abs(subgrido%lon(no)-subgridi%lon(ni))*re * &
                           0.5_r8*(subgrido%coslat(no)+subgridi%coslat(ni))
                      dist = dx*dx + dy*dy
                      if ( dist < distmin )then
                         distmin = dist
                         nmin = ni
                      end if
                   end if
                end if
             end do
          end if

          ! Error conditions
          if ( distmin == spval )then
             write(iulog,*) 'ERROR interpinic set_mindist: Cannot find the closest output ni,no,type= ',&
                  ni, no,subgridi%name
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if

          mindist_index(no) = nmin

       end if ! end if activeo block
    end do
!$OMP END PARALLEL DO
    
  end subroutine set_mindist

  !=======================================================================

  logical function is_sametype (ni, no, subgridi, subgrido)

    ! --------------------------------------------------------------------
    ! arguments
    integer           , intent(in)  :: ni 
    integer           , intent(in)  :: no 
    type(subgrid_type), intent(in)  :: subgridi
    type(subgrid_type), intent(in)  :: subgrido
    ! --------------------------------------------------------------------

    is_sametype = .false.

    if (trim(subgridi%name) == 'pft' .and. trim(subgrido%name) == 'pft') then
       if ( subgridi%ltype(ni) == ilun_landice_multiple_elevation_classes .and. &
            subgrido%ltype(no) == ilun_landice_multiple_elevation_classes) then
          is_sametype = .true.
       else if (subgridi%ptype(ni) == subgrido%ptype(no) .and. &
                subgridi%ctype(ni) == subgrido%ctype(no) .and. &
                subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'column' .and. trim(subgrido%name) == 'column') then
       if ( subgridi%ltype(ni) == ilun_landice_multiple_elevation_classes  .and. &
            subgrido%ltype(no) == ilun_landice_multiple_elevation_classes ) then
          is_sametype = .true.
       else if (subgridi%ctype(ni) == subgrido%ctype(no) .and. &
                subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'landunit' .and. trim(subgrido%name) == 'landunit') then
       if (subgridi%ltype(ni) == subgrido%ltype(no)) then
          is_sametype = .true.
       end if
    else if (trim(subgridi%name) == 'gridcell' .and. trim(subgrido%name) == 'gridcell') then
       is_sametype = .true.
    else 
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: is_sametype check on input and output type not supported'
          write(iulog,*)'typei = ',trim(subgridi%name)
          write(iulog,*)'typeo = ',trim(subgrido%name)
       end if
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end function is_sametype

  !=======================================================================

  logical function is_baresoil (n, subgrid)

    ! --------------------------------------------------------------------
    ! arguments
    integer           , intent(in)  :: n 
    type(subgrid_type), intent(in)  :: subgrid
    ! --------------------------------------------------------------------

    is_baresoil = .false.

    if (subgrid%name == 'pft') then
       if (subgrid%ptype(n) == ipft_not_vegetated) then
          is_baresoil = .true.
       end if
    else if (subgrid%name == 'column') then
       if (subgrid%ctype(n) == icol_vegetated_or_bare_soil) then
          is_baresoil = .true.
       end if
    else if (subgrid%name == 'landunit') then
       if (subgrid%ltype(n) == ilun_vegetated_or_bare_soil) then
          is_baresoil = .true.
       end if
    else 
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: is_baresoil subgrid type ',subgrid%name,' not supported'
       end if
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end function is_baresoil

  !=======================================================================

  subroutine interp_0d_copy (varname, xtype, ncidi, ncido)

    ! --------------------------------------------------------------------
    ! arguments
    character(len=*)   , intent(inout) :: varname
    integer            , intent(in)    :: xtype
    type(file_desc_t)  , intent(inout) :: ncidi
    type(file_desc_t)  , intent(inout) :: ncido
    !
    ! local variables
    integer :: ivalue
    real(r8):: rvalue
    ! --------------------------------------------------------------------
 
    if (masterproc) then
       write(iulog,*) 'Copying      : ',trim(varname)
    end if

    if (xtype == pio_int) then
       call ncd_io(ncid=ncidi, varname=trim(varname), flag='read' , data=ivalue)
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=ivalue)
    else if (xtype == pio_double) then
       call ncd_io(ncid=ncidi, varname=trim(varname), flag='read' , data=rvalue)
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rvalue)
    else
       if (masterproc) then
          write(iulog,*)'ERROR interpinic: unhandled case for var ',trim(varname),' stopping' 
       end if
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine interp_0d_copy

  !=======================================================================

  subroutine interp_1d_double (varname, dimname, begi, endi, bego, endo, ncidi, ncido, &
       activei, activeo, sgridindex)

    ! ------------------------ arguments ---------------------------------
    character(len=*)  , intent(inout) :: varname 
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi
    type(file_desc_t) , intent(inout) :: ncido
    logical           , intent(in)    :: activei(begi:endi)
    logical           , intent(in)    :: activeo(bego:endo)
    integer           , intent(in)    :: sgridindex(bego:endo)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer           :: no,ni          ! indices
    real(r8), pointer :: rbufsli(:)     ! input array
    real(r8), pointer :: rbufslo(:)     ! output array
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Interpolating: ',trim(varname)
    end if

    allocate (rbufsli(begi:endi), rbufslo(bego:endo))
    call ncd_io(ncid=ncidi, varname=trim(varname), flag='read', data=rbufsli)
    call ncd_io(ncid=ncido, varname=trim(varname), flag='read', data=rbufslo, &
         dim1name=dimname)

    do no = bego,endo
       if (activeo(no)) then
          ni = sgridindex(no)
          if (ni > 0) then
             if ( shr_infnan_isnan(rbufsli(ni)) ) then
                rbufslo(no) = spval
             else
                rbufslo(no) = rbufsli(ni)
             end if
          else
             rbufslo(no) = spval
          end if
       else
          if ( shr_infnan_isnan(rbufslo(no)) ) then
             rbufslo(no) = spval
          end if
       end if
    end do

    call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rbufslo, &
         dim1name=dimname)

    deallocate(rbufsli, rbufslo)

  end subroutine interp_1d_double

  !=======================================================================

  subroutine interp_1d_int (varname, dimname, begi, endi, bego, endo, ncidi, ncido, &
       activei, activeo, sgridindex)

    ! ------------------------ arguments ---------------------------------
    character(len=*)  , intent(inout) :: varname 
    character(len=*)  , intent(inout) :: dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    type(file_desc_t) , intent(inout) :: ncidi
    type(file_desc_t) , intent(inout) :: ncido
    logical           , intent(in)    :: activei(begi:endi)
    logical           , intent(in)    :: activeo(bego:endo)
    integer           , intent(in)    :: sgridindex(bego:endo)
    ! --------------------------------------------------------------------

    ! ------------------------ local variables --------------------------
    integer           :: no,ni          !indices
    integer , pointer :: ibufsli(:)     !input array
    integer , pointer :: ibufslo(:)     !output array
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Interpolating: ',trim(varname)
    end if

    allocate (ibufsli(begi:endi), ibufslo(bego:endo))

    call ncd_io(ncid=ncidi, varname=trim(varname), flag='read', &
         data=ibufsli)
    call ncd_io(ncid=ncido, varname=trim(varname), flag='read', &
         data=ibufslo, dim1name=dimname)

    do no = bego,endo
       if (activeo(no)) then
          ni = sgridindex(no)
          if (ni > 0) then
             ibufslo(no) = ibufsli(ni)  
          end if
       end if
    end do

    call ncd_io(ncid=ncido, varname=trim(varname), flag='write', &
         data=ibufslo, dim1name=dimname)

    deallocate (ibufsli, ibufslo)

  end subroutine interp_1d_int

  !=======================================================================

  subroutine interp_2d_double (varname, vec_dimname, lev_dimname, &
       begi, endi, bego, endo, nlevi, nlevo, &
       switchdimi, switchdimo, ncidi, ncido, activei, activeo, sgridindex)

    ! --------------------------------------------------------------------
    ! arguments
    character(len=*)  , intent(inout) :: varname     
    character(len=*)  , intent(inout) :: vec_dimname
    character(len=*)  , intent(inout) :: lev_dimname
    integer           , intent(in)    :: begi, endi
    integer           , intent(in)    :: bego, endo
    integer           , intent(in)    :: nlevi, nlevo
    logical           , intent(inout) :: switchdimi
    logical           , intent(inout) :: switchdimo
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    logical           , intent(in)    :: activei(begi:endi)        
    logical           , intent(in)    :: activeo(bego:endo)        
    integer           , intent(in)    :: sgridindex(bego:endo)
    !
    ! local variables
    integer             :: ni,no               ! indices
    integer             :: ji, jj, index_lower ! indices
    integer             :: status              ! netCDF return code
    integer             :: lev                 ! temporary
    integer             :: start(2), count(2)
    logical             :: copylevels, doneloop
    type(Var_desc_t)    :: vardesc             ! pio variable descriptor
    real(r8), pointer   :: rbuf2do(:,:)        ! output array
    real(r8), pointer   :: rbuf2di(:,:)        ! input array
    real(r8), pointer   :: rbuf1di(:)          ! input array
    real(r8), pointer   :: zsoii(:), zsoio(:) 
    real(r8), parameter :: eps = 1e-6
    ! --------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Interpolating: ',trim(varname)
    end if

    allocate(rbuf2do(bego:endo, nlevo))
    call ncd_io(ncid=ncido, varname=trim(varname), flag='read', data=rbuf2do, &
         dim1name=trim(vec_dimname), switchdim=switchdimo)

    if (nlevi == nlevo) then

       ! Read in 1 level of input array at a time and do interpolation for just that level

       status = pio_inq_varid(ncidi, trim(varname), vardesc)
       allocate(rbuf1di(begi:endi))
       do lev = 1,nlevo
          if (switchdimi) then
             start(1) = lev
             count(1) = 1
             start(2) = 1
             count(2) = endi-begi+1
          else
             start(1) = 1
             count(1) = endi-begi+1
             start(2) = lev
             count(2) = 1
          end if
          status = pio_get_var(ncidi, vardesc, start, count, rbuf1di)
          do no = bego,endo
             if (activeo(no)) then
                ni = sgridindex(no)
                if (ni > 0) then
                   rbuf2do(no,lev) = rbuf1di(ni)
                end if
             else
                if ( shr_infnan_isnan(rbuf2do(no,lev)) ) then
                   rbuf2do(no,lev) = spval
                end if
             end if
          end do
       end do
       deallocate(rbuf1di)
       
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rbuf2do, &
            dim1name=trim(vec_dimname), switchdim=switchdimo)
       
    else if (nlevi < nlevo) then

       ! for the special case when we are regridding from the standard grid to the more
       ! vertlayers grid, calculate a linear interpolation from the one grid to the other

       if (switchdimi) then
          allocate(rbuf2di(nlevi, begi:endi))
       else
          allocate(rbuf2di(begi:endi, nlevi))
       end if
       call ncd_io(ncid=ncidi, varname=trim(varname), flag='read', data=rbuf2di)

       if ( ( nlevi .eq. 15) .and. (nlevo .eq. 30) ) then
          !!! this is the case for variables on the levgrnd grid
          allocate(zsoii(nlevi), zsoio(nlevo))
          status = pio_inq_varid(ncidi, 'zsoi', vardesc)
          status = pio_get_var(ncidi, vardesc, zsoii)
          status = pio_inq_varid(ncido, 'zsoi', vardesc)
          status = pio_get_var(ncido, vardesc, zsoio)  
       else if ( ( nlevi .eq. 20) .and. (nlevo .eq. 35) ) then
          !!! this is the case for variables on the levtot grid
          allocate(zsoii(nlevi), zsoio(nlevo))
          status = pio_inq_varid(ncidi, 'zsoi', vardesc)
          status = pio_get_var(ncidi, vardesc, zsoii(6:20))
          status = pio_inq_varid(ncido, 'zsoi', vardesc)
          status = pio_get_var(ncido, vardesc, zsoio(6:35))
          !! need to put in dummy coordinates for the five snow levels
          zsoii(1:5) = (/ -5._r8, -4._r8, -3._r8, -2._r8, 0._r8 /)
          zsoio(1:5) = (/ -5._r8, -4._r8, -3._r8, -2._r8, 0._r8 /)
       else
          call endrun(msg='ERROR: vertical grid must be either levgrnd or levtot'//&
            errMsg(__FILE__, __LINE__))
       endif
       !
       do ji = 1, nlevo
          doneloop = .false.
          do jj = 1, nlevi
             if ( .not. doneloop) then
                if ( (abs(zsoio(ji) - zsoii(jj))  <  eps ) .or. (jj .eq. nlevi) ) then
                   doneloop = .true.
                   copylevels = .true.
                   index_lower = jj
                else if ( (zsoio(ji)  >  zsoii(jj)) .and. (zsoio(ji)  <  zsoii(jj+1)) ) then
                   doneloop = .true.
                   copylevels = .false.
                   index_lower = jj
                end if
             end if
          end do
          do no = bego,endo
             if (activeo(no)) then
                ni = sgridindex(no)
                if (ni > 0) then
                   if (switchdimi) then
                      if ( copylevels) then
                         rbuf2do(no,ji) = rbuf2di(index_lower,ni)
                      else
                         rbuf2do(no,ji) = &
                              rbuf2di(index_lower+1,ni) &
                              * (zsoio(ji) - zsoii(index_lower)) &
                              / (zsoii(index_lower+1) - zsoii(index_lower)) + &
                              rbuf2di(index_lower,ni) &
                              * (zsoii(index_lower+1) - zsoio(ji) ) &
                              / (zsoii(index_lower+1) - zsoii(index_lower))
                      end if
                   else
                      if ( copylevels) then
                         rbuf2do(no,ji) = rbuf2di(ni,index_lower)
                      else
                         rbuf2do(no,ji) = &
                              rbuf2di(ni,index_lower+1) &
                              * (zsoio(ji) - zsoii(index_lower)) &
                              / (zsoii(index_lower+1) - zsoii(index_lower)) + &
                              rbuf2di(ni,index_lower) &
                              * (zsoii(index_lower+1) - zsoio(ji) ) &
                              / (zsoii(index_lower+1) - zsoii(index_lower))
                      end if
                   end if
                end if
             end if
          end do
       end do
       deallocate(rbuf2di, zsoii, zsoio)
       !
       call ncd_io(ncid=ncido, varname=trim(varname), flag='write', data=rbuf2do, &
            dim1name=trim(vec_dimname), switchdim=switchdimo)

    else

       call endrun(msg='ERROR: new grid must have more vertical levels than old grid'//&
            errMsg(__FILE__, __LINE__))

    end if

    deallocate(rbuf2do)

  end subroutine interp_2d_double

  !=======================================================================

  subroutine check_dim_subgrid(ncidi, ncido, dimname, dimleni, dimleno)

    ! --------------------------------------------------------------------
    ! arguments
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    character(len=*)  , intent(in)    :: dimname
    integer           , intent(out)   :: dimleni
    integer           , intent(out)   :: dimleno
    !
    ! local variables
    integer :: status
    integer :: dimid
    ! --------------------------------------------------------------------

    status = pio_inq_dimid (ncidi, dimname, dimid)
    status = pio_inq_dimlen(ncidi, dimid  , dimleni)
    status = pio_inq_dimid (ncido, dimname, dimid)
    status = pio_inq_dimlen(ncido, dimid  , dimleno)

  end subroutine check_dim_subgrid

  !=======================================================================

  subroutine check_dim_level(ncidi, ncido, dimname)

    ! --------------------------------------------------------------------
    ! arguments
    type(file_desc_t) , intent(inout) :: ncidi         
    type(file_desc_t) , intent(inout) :: ncido         
    character(len=*)  , intent(in)    :: dimname
    !
    ! local variables
    integer :: status
    integer :: dimid
    integer :: dimleni, dimleno
    ! --------------------------------------------------------------------

    status = pio_inq_dimid (ncidi, dimname, dimid)
    status = pio_inq_dimlen(ncidi, dimid  , dimleni)
    status = pio_inq_dimid (ncido, dimname, dimid)
    status = pio_inq_dimlen(ncido, dimid  , dimleno)

    if ( (trim(dimname) == 'levgrnd') .or. (trim(dimname) == 'levtot') ) then
       if (dimleni /= dimleno) then
          if (masterproc) then
             write (iulog,*) 'input and output ',trim(dimname),' values disagree'
             write (iulog,*) 'input nlevgrnd = ',dimleni,' output nlevgrnd = ',dimleno
          end if
          if (dimleni > dimleno) then
             if (masterproc) then
                write (iulog,*) 'ERROR interpinic: input > output for variable ',trim(dimname),'; not supported, stopping'
             end if
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
       end if
    else if (dimleni/=dimleno) then
       !!if (masterproc) then
          write (iulog,*) 'ERROR interpinic: input and output ',trim(dimname),' values disagree'
          write (iulog,*) 'input dimlen = ',dimleni,' output dimlen = ',dimleno
       !!end if
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine check_dim_level

end module initInterpMod
