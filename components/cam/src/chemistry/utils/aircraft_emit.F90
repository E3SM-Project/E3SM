module aircraft_emit
  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! Manages reading and interpolation of aircraft aerosols
  !
  ! Authors: Chih-Chieh (Jack) Chen and Cheryl Craig -- February 2010
  ! Authors: Balwinder Singh and Bryce Harrop - added native grid forcing
  !          files code - September 2019
  !-----------------------------------------------------------------------

  !IMPORTANT:
  ! ** figure out if we need "delta days in time_coord%initialize ***

  !BSINGH - future work:
  ! 1. It will nice to have seprate arrays for forcings which needs to be native grid vs. the others,
  !    as of now, we allocate both fields
  ! 2. Unify the code calling "infld" routine with co2_data_flux.F90 code as much of the code is
  !    repeated
  ! 3. Add "called by" for each subroutine
  ! 4. Code has been modified subtantially, so there may be unused variables.

  use shr_kind_mod,     only: r8 => shr_kind_r8, cx =>SHR_KIND_CX, cl =>SHR_KIND_CL, &
       cs =>SHR_KIND_CS, cxx =>SHR_KIND_CXX
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use tracer_data,      only: trfld, trfile
  use cam_logfile,      only: iulog
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use input_data_utils, only: time_coordinate

  implicit none
  private
  save

  public :: aircraft_emit_init
  public :: aircraft_emit_adv
  public :: aircraft_emit_register
  public :: aircraft_emit_readnl

  !------------------------------------------------------------------
  !DEFINITION:
  !"native grid forcing file": A forcing file which has to be on the
  !same grid horizontally as the model is running on. For example,
  !if the model is running on ne30 grid, forcing file has to be on
  !ne30 grid horizontally. The vertical resolution can be different
  !from the model's vertical resolution.
  !------------------------------------------------------------------

  type :: forc_air_native_grid
     !------------------------------------------------------------------
     !"forc_air_native_grid" is forcing from files which has to be on the
     !native grid (only in horizontal,vertical resolution may be different
     !from the model's grid resolution)
     !That is, forcing files has to be on same grid as the grid used for
     !the model run
     !------------------------------------------------------------------

     !Number of levels in the 3D forcing file
     integer                               :: lev_frc

     !Data structure to store two time samples from a file to do time interpolation in the next step
     !(pcols,lev_frc,begchunk:endchunk,2)
     real(r8), pointer, dimension(:,:,:,:) :: native_grid_flds_tslices

     !Data structure to store data after time interpolation from two time samples
     !(pcols,lev_frc,begchunk:endchunk)
     real(r8), pointer, dimension(:,:,:)   :: native_grid_flds

     !Data structure to keep track of time
     type(time_coordinate) :: time_coord

     !specie name
     character( len = cx)  :: spc_name_ngrd


     !Level bounds read from input file
     real(r8), pointer, dimension(:,:) :: lev_bnds

     !Forcing file name
     character( len = cx)  :: input_file

     !Units of forcing data
     character( len = cs)  :: units

     !logical to control first data read
     logical               :: initialized

     !pbuf index to store read in data in pbuf
     integer               :: pbuf_ndx = -1
  end type forc_air_native_grid
  type(forc_air_native_grid),allocatable :: native_grid_frc_air(:)

  type :: forcing_air
     real(r8)              :: mw
     character(len=cl)     :: filelist
     character(len=cl)     :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=11)     :: species
     character(len=8)      :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld),pointer       :: fields(:)
     type(trfile)              :: file
  end type forcing_air

  type(forcing_air), allocatable :: forcings_air(:)

  integer, parameter :: N_AERO = 10

  character(len=11), parameter :: aero_names(N_AERO) = (/'ac_HC      ','ac_NOX     ','ac_PMNV    ',&
       'ac_PMSO    ','ac_PMFO    ','ac_FUELBURN','ac_CO2     ','ac_H2O     ',&
       'ac_SOX     ','ac_CO      '/)

  real(r8), parameter :: molmass(N_AERO) = 1._r8

  logical, parameter  :: advective_tracer(N_AERO) = (/.false., .false., .false., .false., .false., &
       .false., .false., .false., .false.,.false./)

  !Logical to control whether to force that horizontal grid is on native grid or not
  logical, parameter  :: horz_native(N_AERO) = (/.false., .false., .false., .false., .false., &
       .false., .true., .false., .false.,.false./)

  !Logical to control whether to convert units to mmr (mass mixing ratios) or not
  !sometimes unit coversion is handled at other places in the code
  ! (e.g. unit conversion of ac_CO2 is handled in co2_cycle.F90)
  logical, parameter  :: convert_to_mmr(N_AERO) = (/.true., .true., .true., .true., .true., &
       .true., .false., .true., .true.,.true./)

  !mixtyoe is only used for advected tracers
  character(len=3) :: mixtype(N_AERO) = (/'wet','wet','wet','wet','wet','wet','wet','wet','wet','wet'/)

  real(r8) :: cptmp = 666.0_r8
  real(r8) :: qmin = 0.0_r8
  logical  :: readiv = .false.
  logical  :: has_fixed_ubs = .false.
  logical  :: cam_outfld = .false.

  integer            :: index_map(N_AERO)
  character(len=cl)  :: air_specifier(N_AERO)=''
  character(len=cl)  :: air_datapath=''
  character(len=24)  :: air_type = 'CYCLICAL_LIST' ! 'CYCLICAL_LIST'

  logical            :: rmv_file = .false.
  logical            :: dimnames_set = .false.

  integer :: number_flds

  integer :: aircraft_cnt = 0

  character(len=16)  :: spc_name_list(N_AERO)
  character(len=cl)  :: spc_flist(N_AERO),spc_fname(N_AERO)
  character(len=8)   :: dim1name, dim2name
  integer, parameter :: huge_int = huge(1)

contains

  subroutine aircraft_emit_register()

    !------------------------------------------------------------------
    ! **** Add the aircraft aerosol data to the physics buffer ****
    ! called by:
    !------------------------------------------------------------------
    use physics_buffer, only: pbuf_add_field, dtype_r8
    use tracer_data,    only: incr_filename
    use constituents,   only: cnst_add
    use ppgrid,         only: pver, pcols

    integer            :: i,idx, mm, ind, n
    character(len=16)  :: spc_name
    character(len=cl)  :: filelist, curr_filename
    character(len=128) :: long_name
    logical            :: has_fixed_ubc=.false.
    logical            :: read_iv=.false.

    !------------------------------------------------------------------
    ! Return if air_specifier is blank (no aircraft data to process)
    !------------------------------------------------------------------
    if (air_specifier(1) == "") return

    ! count aircraft emission species used in the simulation
    count_emis: do n=1,N_AERO

       if( len_trim(air_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(air_specifier(n),'->')
       spc_name = trim(adjustl(air_specifier(n)(:i-1)))
       filelist = trim(adjustl(air_specifier(n)(i+2:)))

       mm = get_aircraft_ndx(spc_name)
       if( mm < 1 ) then
          call endrun(trim(spc_name)//' is not in the aircraft emission dataset '//errmsg(__FILE__,__LINE__))
       endif

       aircraft_cnt = aircraft_cnt + 1
       call pbuf_add_field(aero_names(mm),'physpkg',dtype_r8,(/pcols,pver/),idx)

       spc_flist(aircraft_cnt) = filelist
       spc_name_list(aircraft_cnt) = spc_name
       index_map(aircraft_cnt) = mm

       !"incr_filename" reads the filename mentioned inside file "spc_flist(aircraft_cnt)" at path "air_datapath"
       curr_filename=''
       spc_fname(aircraft_cnt) = incr_filename( curr_filename, filenames_list=spc_flist(aircraft_cnt), datapath=air_datapath)

       if( advective_tracer(mm) ) then
          long_name = 'aircraft_'//trim(spc_name)
          call cnst_add(aero_names(mm),molmass(mm),cptmp,qmin,ind,longname=long_name,readiv=read_iv, &
               mixtype=mixtype(mm),cam_outfld=cam_outfld,fixed_ubc=has_fixed_ubc)
       endif

    enddo count_emis
    ! count aircraft emission species used in the simulation

  endsubroutine aircraft_emit_register

  subroutine aircraft_emit_init(state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Initialize the aircraft aerosol data handling ****
    ! called by:
    !-------------------------------------------------------------------
    use cam_history,      only: addfld, add_default
    use tracer_data,      only: trcdata_init
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk, pcols
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_index
    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use cam_grid_support, only: cam_grid_id, cam_grid_check
    use cam_grid_support, only: cam_grid_get_dim_names
    use dyn_grid,         only: get_horiz_grid_dim_d
    use dycore,           only: dycore_is
    use cam_pio_utils,    only: cam_pio_openfile
    use pio,              only: file_desc_t, pio_nowrite, pio_closefile, pio_inq_dimid, pio_bcast_error, &
         pio_seterrorhandling, pio_noerr, pio_inquire_dimension, pio_get_att, pio_inq_varid, pio_get_var


    implicit none

    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !local vars
    type(file_desc_t)   :: fh
    character(len=16)   :: spc_name
    character(len=cxx)  :: err_str

    integer :: ndx, istat, i, astat, m, n, mm, c
    integer :: grid_id
    integer :: dimlevid, var_id, errcode, dim1id, dim2id, dim1len, dim2len
    integer :: dimbndid, nbnd
    integer :: hdim1_d, hdim2_d    ! model grid size

    real(r8) :: dtime

    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------
    if (aircraft_cnt == 0 ) return


    !------------------------------------------------------------------
    ! For forcing files which has to be on the native grid,dimensions
    ! are set in the following if condition
    !------------------------------------------------------------------
    if( any( horz_native(:) ) ) then
       if (.not. dimnames_set) then
          grid_id = cam_grid_id('physgrid')
          if (.not. cam_grid_check(grid_id)) then
             call endrun('no "physgrid" grid:'//errmsg(__FILE__,__LINE__))
          endif
          !dim1name and dim2name are populated here with the grid dimension the model is running on (e.g. ne30, lat, lon etc.)
          !For SE grid, dim1name = dim2name = "ncol"
          !For FV grid, dim1name = lon, dim2name = lat
          call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
          dimnames_set = .true.
       end if

       !--------------------------------------------------------------------------------
       ! allocate forcings type array for native grid forcing files
       !--------------------------------------------------------------------------------
       allocate( native_grid_frc_air(aircraft_cnt), stat=astat )
       if( astat /= 0 ) then
          write(err_str,*) 'failed to allocate native_grid_frc_air array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       end if
    endif

    if (masterproc) write(iulog,*) ' '

    if( any( .not. horz_native(:) ) ) then
       !-----------------------------------------------------------------------
       !       allocate forcings type array
       !-----------------------------------------------------------------------
       allocate( forcings_air(aircraft_cnt), stat=astat )
       if( astat/= 0 ) then
          write(err_str,*) 'failed to allocate forcings_air array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       end if
    endif


    !-----------------------------------------------------------------------
    !       setup the forcings_air type array
    !-----------------------------------------------------------------------
    species_loop : do m = 1,aircraft_cnt

       spc_name = spc_name_list(m)

       if( horz_native(index_map(m))) then
          !-----------------------------------------------------------------------
          !       initialize variables for native grid forcing files
          !-----------------------------------------------------------------------

          native_grid_frc_air(m)%spc_name_ngrd  = spc_name
          native_grid_frc_air(m)%input_file     = trim(air_datapath)//'/'//trim(spc_fname(m))

          native_grid_frc_air(m)%initialized    = .false.
          dtime = 1.0_r8 - 200.0_r8 / 86400.0_r8
          call native_grid_frc_air(m)%time_coord%initialize(trim(adjustl(native_grid_frc_air(m)%input_file)), &
               force_time_interp=.true., delta_days=dtime)

          !-----------------------------------------------------------------------
          !       Open file
          !-----------------------------------------------------------------------
          call cam_pio_openfile(fh, trim(adjustl(native_grid_frc_air(m)%input_file)), PIO_NOWRITE)

          !ask PIO to return the control if it experiences an error so that we can
          !handle it explicitly in the code
          call pio_seterrorhandling(fh, pio_bcast_error)

          !-----------------------------------------------------------------------
          !       Sanity checks for the native grid
          !-----------------------------------------------------------------------

          !if forcing file is on a different grid than the model grid
          !(e.g. model is running on an FV grid and forcing netcdf file is on an SE grid), exit with an error
          if(pio_inq_dimid(fh, trim(adjustl(dim1name)), dim1id) /= pio_noerr) then
             !pio_inq_dimid function tries to find dim1name in file with id "fh"
             !if it can't find dim1name, it means there is a mismacth in model and netcdf
             !file grid
             call endrun('grid mismatch, failed to find '//dim1name//' dimension in file:'&
                  ' '//trim(adjustl(native_grid_frc_air(m)%input_file))//' '&
                  ' '//errmsg(__FILE__,__LINE__))
          endif

          !find if the model and netcdf file has same grid resolution
          call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
          if( dycore_is('SE') )  then
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr) then
                if(dim1len /= hdim1_d ) then !compare model grid length with file's
                   write(err_str,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                        hdim1_d,'), netcdf file is:'//trim(adjustl(native_grid_frc_air(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_frc_air(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          elseif( dycore_is('LR')) then
             if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim2id)) then !obtain lat dimension of model
                call endrun('failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
                     ' '//trim(adjustl(native_grid_frc_air(m)%input_file))//' '//errmsg(__FILE__,__LINE__))
             endif
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr .and. &
                  pio_inquire_dimension(fh, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
                if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
                   write(err_str,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                        hdim1_d,' x ',hdim2_d,'), netcdf file is:'//trim(adjustl(native_grid_frc_air(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_frc_air(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          else
             call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
          endif

          !Find the value of vertical levels in the forcing file
          if( pio_inq_dimid(fh, 'lev', dimlevid) ==  pio_noerr ) then
             if ( pio_inquire_dimension(fh, dimlevid, len =  native_grid_frc_air(m)%lev_frc) /=  pio_noerr ) then
                write(err_str,*)'failed to obtain value of "lev" dimension from file:',&
                     trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
             !obtain level bounds needed for vertical interpolation
             if( pio_inq_varid(fh, 'lev_bnds', var_id) ==  pio_noerr ) then
                !get dimension "bound"
                if( pio_inq_dimid(fh, 'bound', dimbndid) ==  pio_noerr ) then
                   if ( pio_inquire_dimension(fh, dimbndid, len = nbnd) ==  pio_noerr ) then
                      !"nbnd" has to be 2 (it is obvious but adding a check here doesn't hurt)
                      if(nbnd /= 2) then
                         write(err_str,*)'"bound" should be equal to 2, bound=',nbnd,' in file:', &
                              trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                         call endrun(err_str)
                      endif
                      allocate(native_grid_frc_air(m)%lev_bnds(nbnd,native_grid_frc_air(m)%lev_frc))
                      if (pio_get_var(fh, var_id,native_grid_frc_air(m)%lev_bnds) /=  pio_noerr ) then
                         write(err_str,*)'failed to read "lev_bnds" variable from file:',&
                              trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                         call endrun(err_str)
                      endif
                   else
                      write(err_str,*)'failed to obtain value of "bound" dimension from file:',&
                           trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                      call endrun(err_str)
                   endif
                else
                   write(err_str,*)'failed to inquire "bound" dimension from file:',&
                        trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                   call endrun(err_str)
                endif
             else
                write(err_str,*)'failed to obtain "lev_bnds" variable from file:',&
                     trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'Dimension "lev" is not found in:',&
                  trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif

          !get units of the data in the forcing file
          if(pio_inq_varid( fh, spc_name, var_id ) == pio_noerr ) then
             if(pio_get_att( fh, var_id, 'units', native_grid_frc_air(m)%units) .ne. pio_noerr ) then
                write(err_str,*)'failed to obtain units of variable ',trim(spc_name),' in &
                     &file:',trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'variable ',trim(spc_name),' not found in:',trim(adjustl(native_grid_frc_air(m)%input_file)), &
                  ',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif

          !close file
          call pio_closefile(fh)

          !allocate arrays to stroe data for interpolation in time
          allocate(native_grid_frc_air(m)%native_grid_flds_tslices(pcols, native_grid_frc_air(m)%lev_frc, &
               begchunk:endchunk,2), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_frc_air(',m,')%native_grid_flds_tslices array;&
                  error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif

          !allocate arrays to hold data before the vertical interpolation
          allocate(native_grid_frc_air(m)%native_grid_flds(pcols, native_grid_frc_air(m)%lev_frc,begchunk:endchunk), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_frc_air(',m,')%native_grid_flds array; error = ',&
                  astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif

          !get pbuf index to store the field in pbuf
          native_grid_frc_air(m)%pbuf_ndx = pbuf_get_index(spc_name,errcode)
          if(errcode < 0 ) then
             write(err_str,*)'failed to get pbuf index for specie:',spc_name,' errorcode is:',errcode,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif

       else

          allocate( forcings_air(m)%sectors(1), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'aircraft_emit_init: failed to allocate forcings_air%sectors &
                  &array; error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          end if

          allocate( forcings_air(m)%fields(1), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'aircraft_emit_init: failed to allocate forcings_air%fields &
                  &array; error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          end if

          !-----------------------------------------------------------------------
          !         default settings
          !-----------------------------------------------------------------------
          forcings_air(m)%file%stepTime    = .true.  ! Aircraft data is not to be interpolated in time
          forcings_air(m)%file%cyclical_list    = .true.  ! Aircraft data cycles over the filename list
          forcings_air(m)%file%weight_by_lat     = .true.  ! Aircraft data -  interpolated with latitude weighting
          forcings_air(m)%file%conserve_column = .true. ! Aircraft data - vertically interpolated to conserve the total column
          forcings_air(m)%species          = spc_name
          forcings_air(m)%sectors          = spc_name ! Only one species per file for aircraft data
          forcings_air(m)%nsectors         = 1
          forcings_air(m)%filelist         = spc_flist(m)
          !         forcings_air(m)%file%curr_filename    = spc_fname(m)
          forcings_air(m)%filename         = spc_fname(m)
       endif

       call addfld( trim(spc_name), (/ 'lev' /), 'A',  '1/s',     &
            'aircraft emission '//trim(spc_name) )
       call add_default( trim(spc_name), 1, ' ' )

    end do species_loop

    if (masterproc) then
       !-----------------------------------------------------------------------
       !            diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'aircraft_emit_init: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'aircraft_emit timing specs'
       write(iulog,*) 'type = ',air_type
       write(iulog,*) ' '
       write(iulog,*) 'there are ',aircraft_cnt,' species of aircraft emission'
       do m = 1,aircraft_cnt
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',spc_name_list(m)
          write(iulog,*) 'filelist= ',spc_flist(m)
       end do
       write(iulog,*) ' '
    endif


    !------------------------------------------------------------------
    !       Initialize the aircraft file processing
    !------------------------------------------------------------------
    do m=1,aircraft_cnt

       number_flds = 0
       if(horz_native(index_map(m))) then
          if (associated(native_grid_frc_air(m)%native_grid_flds_tslices)) &
               number_flds = 1
          !read the forcing file once to initialize variable including time cordinate
          call advance_native_grid_data( native_grid_frc_air(m) )
          native_grid_frc_air(m)%initialized = .true.
       else
          allocate (forcings_air(m)%file%in_pbuf(size(forcings_air(m)%sectors)))
          forcings_air(m)%file%in_pbuf(:) = .true.
          if (associated(forcings_air(m)%fields)) number_flds = size( forcings_air(m)%fields )

          call trcdata_init( forcings_air(m)%sectors, forcings_air(m)%filename, forcings_air(m)%filelist, air_datapath, &
               forcings_air(m)%fields, forcings_air(m)%file, rmv_file, 0, 0, 0, air_type)
       endif

       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(err_str,*) 'There are no aircraft aerosols ',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
       end if
    end do
  end subroutine aircraft_emit_init


  subroutine aircraft_emit_adv( state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Advance to the next aircraft data ****
    ! called by:
    !-------------------------------------------------------------------

    use perf_mod,     only: t_startf, t_stopf
    use tracer_data,  only: advance_trcdata
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use ppgrid,       only: pcols, pver
    use string_utils, only: to_lower, GLC
    use cam_history,  only: outfld
    use physconst,    only: mwdry       ! molecular weight dry air ~ kg/kmole
    use physconst,    only: boltz                ! J/K/molecule
    ! C.-C. Chen
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer  :: ind, c, ncol, i, caseid, m, pbuf_ndx
    real(r8) :: to_mmr(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)
    character(len = cs) :: units_spc

    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------
    if (aircraft_cnt == 0 ) return
    call t_startf('All_aircraft_emit_adv')

    !-------------------------------------------------------------------
    !    For each field, read more data if needed and interpolate it to the current model time
    !-------------------------------------------------------------------
    do m = 1, aircraft_cnt

       if(horz_native(index_map(m))) then ! if horizontal grid is native
          units_spc = native_grid_frc_air(m)%units
          pbuf_ndx  = native_grid_frc_air(m)%pbuf_ndx

          !read in next time slice (if needed) and interpolate in time
          !following call just reads in time slices in horizontal
          !vertical interpolation is done in the next call
          call advance_native_grid_data( native_grid_frc_air(m) )

          !do vertical interpolation

          !following call needs state to get ncol for each chunk and
          !pbuf for storing the interpolated (time and vertically)
          !field in pbuf
          call vert_interp( state, pbuf_ndx, native_grid_frc_air(m), pbuf2d)

       else
          units_spc = forcings_air(m)%fields(i)%units
          pbuf_ndx  = forcings_air(m)%fields(i)%pbuf_ndx
          call advance_trcdata( forcings_air(m)%fields, forcings_air(m)%file, state, pbuf2d)
       endif

       !-------------------------------------------------------------------
       !    set the tracer fields with the correct units
       !-------------------------------------------------------------------
       do i = 1,number_flds

          !initialize caseid so that it is not used inadvertantly
          caseid = huge_int
          ind    = index_map(i)

          !only assign valid integer if we need unit conversion
          if ( convert_to_mmr(ind) ) then
             ! C.-C. Chen, adding case 4  for kg/sec
             select case ( to_lower(trim(units_spc(:GLC(units_spc)))) )
             case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
                caseid = 1
             case ('kg/kg','mmr')
                caseid = 2
             case ('mol/mol','mole/mole','vmr','fraction')
                caseid = 3
             case ('kg/kg/sec')
                caseid = 4
             case default
                print*, 'aircraft_emit_adv: units = ',trim(units_spc) ,' are not recognized'
                call endrun('aircraft_emit_adv: units are not recognized '//errmsg(__FILE__, __LINE__))
             end select
          endif

          !$OMP PARALLEL DO PRIVATE (C, NCOL, TO_MMR, tmpptr, pbuf_chnk)
          do c = begchunk,endchunk
             ncol = state(c)%ncol

             !initialize to_mmr to 1.0 for cases where unit conversion is not required
             !(i.e., caseid is not assigned an integer value)
             to_mmr(:ncol,:) = 1.0_r8

             if (caseid == 1) then ! change it to select-case?
                to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
             elseif(caseid == 2) then
                to_mmr(:ncol,:) = 1.0_r8
             elseif(caseid == 3) then
                to_mmr(:ncol,:) = molmass(ind)/mwdry
             elseif(caseid == 4) then
                to_mmr(:ncol,:) = 1.0_r8
             endif

             pbuf_chnk => pbuf_get_chunk(pbuf2d, c)

             call pbuf_get_field(pbuf_chnk, pbuf_ndx , tmpptr )

             tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

             call outfld( trim(spc_name_list(m)), &
                  tmpptr, ncol, state(c)%lchnk )
          enddo
       enddo
    enddo

    call t_stopf('All_aircraft_emit_adv')
  end subroutine aircraft_emit_adv


  subroutine advance_native_grid_data( native_grid_strct )
    !-------------------------------------------------------------------
    !    This subroutine reads the data from the native grid and
    !    interpolates in time
    ! called by:
    !-------------------------------------------------------------------

    use ppgrid,         only: begchunk, endchunk, pcols
    use ncdio_atm,      only: infld
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, pio_nowrite, pio_closefile

    implicit none

    !args
    type(forc_air_native_grid), intent (inout) :: native_grid_strct

    !local vars
    type(file_desc_t) :: fh
    character(len=cs) :: spc_name

    logical  :: read_data
    integer  :: indx2_pre_adv
    logical  :: found

    integer :: nstep, i, j, k

    !obtain name of the specie
    spc_name = native_grid_strct%spc_name_ngrd

    !Decide whether to read new data or not (e.g. data may needs to be read on month boundaries )
    read_data = native_grid_strct%time_coord%read_more() .or. .not. native_grid_strct%initialized

    !Find time index to decide whether to read new data or recycle previously read data
    indx2_pre_adv = native_grid_strct%time_coord%indxs(2)

    !compute weights for time interpolation (time_coord%wghts) by advancing in time
    call native_grid_strct%time_coord%advance()

    if ( read_data ) then

       !open file
       call cam_pio_openfile(fh, trim(adjustl(native_grid_strct%input_file)), PIO_NOWRITE)

       ! read time-level 1
       if (native_grid_strct%initialized .and. native_grid_strct%time_coord%indxs(1) == indx2_pre_adv) then
          ! skip the read if the needed vals for time level 1 are present in time-level 2
          native_grid_strct%native_grid_flds_tslices(:,:,:,1) = native_grid_strct%native_grid_flds_tslices(:,:,:,2)
       else
          !NOTE: infld call doesn't do any interpolation in space, it just reads in the data
          call infld(trim(spc_name), fh, dim1name, dim2name, 'lev',&
               1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
               native_grid_strct%native_grid_flds_tslices(:,:,:,1), found, &
               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          if (.not. found) then
             call endrun(trim(spc_name) // ' not found '//errmsg(__FILE__,__LINE__))
          endif
       endif

       ! read time level 2
       call infld(trim(spc_name), fh, dim1name, dim2name, 'lev',&
            1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
            native_grid_strct%native_grid_flds_tslices(:,:,:,2), found, &
            gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(2))




       if (.not. found) then
          call endrun(trim(spc_name) // ' not found '//errmsg(__FILE__,__LINE__))
       endif

       !close file
       call pio_closefile(fh)
    endif

    ! interpolate between time-levels
    ! If time:bounds is in the dataset, and the dataset calendar is compatible with EAM's,
    ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
    ! generating fluxes that are piecewise constant in time.

    if (native_grid_strct%time_coord%wghts(2) == 0.0_r8) then
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1)
    else
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1) + &
            native_grid_strct%time_coord%wghts(2) * (native_grid_strct%native_grid_flds_tslices(:,:,:,2) - &
            native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    endif

  end subroutine advance_native_grid_data


  subroutine vert_interp( state, pbuf_ndx, native_grid_strct, pbuf2d )

    !-------------------------------------------------------------------
    !    This subroutine interpolates in vertical direction. Finally the
    !    interpolated data is stored  in pbuf
    ! called by:
    !-------------------------------------------------------------------

    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use ppgrid,         only: pver, pcols


    !args
    type(physics_state), intent(in)        :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer     :: pbuf2d(:,:)
    type(forc_air_native_grid), intent(in) :: native_grid_strct
    integer, intent(in)                    :: pbuf_ndx

    !local vars
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer  :: ic, ncol, icol
    real(r8) :: vrt_interp_field(pcols, pver,begchunk:endchunk)
    real(r8),pointer :: tmpptr_native_grid(:,:)

    real(r8) :: nzil, nzir, mzil, mzir                       ! level bounds
    real(r8) :: ovrl, ovrr, ovrf                             ! overlap bounds
    real(r8) :: ovrmat(pver, native_grid_strct%lev_frc) ! overlap fractions matrix
    integer  :: kinp, kmdl                              ! vertical indexes

    ! Vertical interpolation follows Joeckel (2006) ACP method for conservation
    do ic = begchunk,endchunk
       ncol = state(ic)%ncol

       do icol = 1, ncol
          do kinp = 1, native_grid_strct%lev_frc
             nzil = native_grid_strct%lev_bnds(1,kinp)
             nzir = native_grid_strct%lev_bnds(2,kinp)
             do kmdl = 1, pver
                mzil               = state(ic)%zi(icol, kmdl)
                mzir               = state(ic)%zi(icol, kmdl+1)
                ovrl               = MAX( MIN(nzil, nzir), MIN(mzil, mzir) )
                ovrr               = MIN( MAX(nzil, nzir), MAX(mzil, mzir) )
                ovrf               = INT( SIGN(0.5_r8, ovrr-ovrl) + 1._r8 )
                ovrmat(kmdl,kinp) = ABS(ovrr - ovrl) * ovrf / ABS(nzir - nzil)
             end do
          end do
          vrt_interp_field(icol,:,ic)  = MATMUL( ovrmat, native_grid_strct%native_grid_flds(icol,:,ic) )
       end do
    enddo

    !future work: these two (above abd below) do loops should be combined into one.

    !add field to pbuf
    !$OMP PARALLEL DO PRIVATE (IC, NCOL, tmpptr_native_grid, pbuf_chnk, vrt_interp_field)
    do ic = begchunk,endchunk
       ncol = state(ic)%ncol
       pbuf_chnk => pbuf_get_chunk(pbuf2d, ic)
       call pbuf_get_field(pbuf_chnk, pbuf_ndx , tmpptr_native_grid )
       tmpptr_native_grid(:ncol,:) = vrt_interp_field(:ncol,:,ic)
    enddo

  end subroutine vert_interp

  subroutine aircraft_emit_readnl(nlfile)
    !-------------------------------------------------------------------
    ! **** Read in the aircraft_emit namelist *****
    ! called by:
    !-------------------------------------------------------------------
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aircraft_emit_readnl'

    character(len=cl)  :: aircraft_specifier(N_AERO)
    character(len=cl)  :: aircraft_datapath
    character(len=24)  :: aircraft_type

    namelist /aircraft_emit_nl/  aircraft_specifier, aircraft_datapath, aircraft_type
    !-----------------------------------------------------------------------------

    ! Initialize namelist variables from local module variables.
    aircraft_specifier= air_specifier
    aircraft_datapath = air_datapath
    aircraft_type     = air_type

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aircraft_emit_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aircraft_emit_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist '//errmsg(__FILE__,__LINE__))
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(aircraft_specifier,len(aircraft_specifier(1))*N_AERO,     mpichar, 0, mpicom)
    call mpibcast(aircraft_datapath, len(aircraft_datapath),                mpichar, 0, mpicom)
    call mpibcast(aircraft_type,     len(aircraft_type),                    mpichar, 0, mpicom)
#endif

    ! Update module variables with user settings.
    air_specifier  = aircraft_specifier
    air_datapath   = aircraft_datapath
    air_type       = aircraft_type

  end subroutine aircraft_emit_readnl

  integer function get_aircraft_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_aircraft_ndx = 0
    do i = 1,N_AERO
       if ( trim(name) == trim(aero_names(i)) ) then
          get_aircraft_ndx = i
          return
       endif
    enddo

  end function get_aircraft_ndx

end module aircraft_emit
