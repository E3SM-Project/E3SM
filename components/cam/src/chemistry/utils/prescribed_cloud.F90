module prescribed_cloud

!--------------------------------------------------------------------------
! Purpose:
!
! Reads cloud-related fields, puts the them into the physics buffer for use
! by radiation
!
!--------------------------------------------------------------------------

!++BEH
!  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_kind_mod,     only : r8 =>shr_kind_r8, cx =>SHR_KIND_CX, cl =>SHR_KIND_CL, &
                               cs =>SHR_KIND_CS, cxx =>SHR_KIND_CXX
!--BEH
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog
!++BEH
  use input_data_utils, only : time_coordinate
  use shr_log_mod ,     only : errMsg => shr_log_errMsg
!--BEH

  implicit none
  private
  save

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_cloud_init
  public :: prescribed_cloud_adv
  public :: write_prescribed_cloud_restart
  public :: read_prescribed_cloud_restart
  public :: has_prescribed_cloud
  public :: prescribed_cloud_register
  public :: init_prescribed_cloud_restart
  public :: prescribed_cloud_readnl

  logical :: has_prescribed_cloud = .false.
!! JGOmod
!  integer          , parameter :: nflds             = 10
!  character(len=16), parameter :: cloud_name(nflds) = (/'DEI_in'   ,'MU_in'  ,'LAMBDAC_in' ,'ICIWP_in' ,'ICLWP_in' ,'DES_in' , &
!                                                        'ICSWP_in' ,'CLD_in' ,'CLDLIQ_in'  ,'CLDICE_in'  /)
!
!  character(len=16)  :: fld_name(nflds)             = (/'DEI_rad'  ,'MU_rad' ,'LAMBDAC_rad','ICIWP_rad','ICLWP_rad','DES_rad', &
!                                                        'ICSWP_rad','CLD_rad','CLDLIQ_rad' ,'CLDICE_rad' /)
! BPM: added 1 to nflds; added CLDFSNOW_* into *_name
  integer          , parameter :: nflds             = 10
! BEH no need for _in, just use _rad so that we don't have to rename variables to be read in
!  character(len=16), parameter :: cloud_name(nflds) = (/'DEI_in'   ,'MU_in'  ,'LAMBDAC_in' ,'ICIWP_in' ,'ICLWP_in' ,'DES_in' , &
!                                                        'ICSWP_in' ,'CLD_in', 'CLDFSNOW_in', 'CONCLD_in' /)
  character(len=16), parameter :: cloud_name(nflds) = (/'DEI_rad'   ,'MU_rad'  ,'LAMBDAC_rad' ,'ICIWP_rad' ,'ICLWP_rad' ,'DES_rad' , &
                                                        'ICSWP_rad' ,'CLD_rad', 'CLDFSNOW_rad', 'CONCLD_rad' /)

  character(len=16)  :: fld_name(nflds)             = (/'DEI_rad'  ,'MU_rad' ,'LAMBDAC_rad','ICIWP_rad','ICLWP_rad','DES_rad', &
                                                        'ICSWP_rad','CLD_rad', 'CLDFSNOW_rad', 'CONCLD_rad'/)
!! JGOmod
  character(len=256) :: filename                    = ' '
  character(len=256) :: filelist                    = ' '
  character(len=256) :: datapath                    = ' '
!BEH  character(len=32)  :: data_type                   = 'SERIAL'
  logical            :: rmv_file                    = .false.
  integer            :: cycle_yr                    = 0
  integer            :: fixed_ymd                   = 0
  integer            :: fixed_tod                   = 0
  character(len=32)  :: specifier(nflds)            = ''


!++BEH
  logical, parameter :: horz_native = .true. ! they will all have the same behavior if they
                                             ! all come from the same file
  logical            :: dimnames_set = .false. 
  integer            :: number_flds
  character(len=16)  :: spc_name_list(nflds)
  character(len=cl)  :: spc_flist(nflds),spc_fnames(nflds)
  character(len=8)   :: dim1name, dim2name
  character(len=32)  :: data_type = 'CYCLICAL' ! 'CYCLICAL'
  real(r8)           :: num_file_years = 0._r8
  character(len=32)  :: air_type  = 'CYCLICAL_LIST'
!  logical, parameter :: horz_native(nflds) = (/.true., .true., .true., .true., .true., &
!       .true., .true., .true., .true., .true./)
!  integer            :: index_map(nflds)
  !------------------------------------------------------------------
  !DEFINITION:
  !"native grid forcing file": A forcing file which has to be on the 
  !same grid horizontally as the model is running on. For example, 
  !if the model is running on ne30 grid, forcing file has to be on
  !ne30 grid horizontally. The vertical resolution can be different 
  !from the model vertical resolution.
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
     character( len = cx)  :: spc_cname_ngrd
     character( len = cx)  :: spc_fname_ngrd


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
  type(forc_air_native_grid),allocatable :: native_grid_cloud(:)

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
!--BEH
contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_cloud_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: i,idx

    if (has_prescribed_cloud) then
       do i = 1,nflds
          call pbuf_add_field(cloud_name(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_cloud_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_cloud_init_orig()

    use tracer_data, only : trcdata_init

    implicit none

    integer :: ndx, istat, i

    if ( has_prescribed_cloud ) then
       if ( masterproc ) then
          write(iulog,*) 'cloud is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    do i = 1,nflds
       specifier(i) = trim(cloud_name(i))//':'//trim(fld_name(i))
    end do


    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    file%stepTime   = .true.
    file%xyzint     = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

  end subroutine prescribed_cloud_init_orig

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_cloud_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_cloud_readnl'

   character(len=256) :: prescribed_cloud_file
   character(len=256) :: prescribed_cloud_filelist
   character(len=256) :: prescribed_cloud_datapath
   character(len=32)  :: prescribed_cloud_type
   integer            :: prescribed_cloud_cycle_yr
   integer            :: prescribed_cloud_fixed_ymd
   integer            :: prescribed_cloud_fixed_tod
   real(r8)           :: prescribed_cloud_num_file_years

   namelist /prescribed_cloud_nl/ &
      prescribed_cloud_file,      &
      prescribed_cloud_filelist,  &
      prescribed_cloud_datapath,  &
      prescribed_cloud_type,      &
      prescribed_cloud_cycle_yr,  &
      prescribed_cloud_fixed_ymd, &
      prescribed_cloud_fixed_tod, &
      prescribed_cloud_num_file_years
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_cloud_file     = filename
   prescribed_cloud_filelist = filelist
   prescribed_cloud_datapath = datapath
   prescribed_cloud_type     = data_type
   prescribed_cloud_cycle_yr = cycle_yr
   prescribed_cloud_fixed_ymd= fixed_ymd
   prescribed_cloud_fixed_tod= fixed_tod
   prescribed_cloud_num_file_years= num_file_years

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_cloud_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_cloud_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_cloud_file,     len(prescribed_cloud_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_filelist, len(prescribed_cloud_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_datapath, len(prescribed_cloud_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_type,     len(prescribed_cloud_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_cloud_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_cloud_fixed_tod,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_cloud_num_file_years, 1, mpir8, 0, mpicom)
#endif

   ! Update module variables with user settings.
   filename   = prescribed_cloud_file
   filelist   = prescribed_cloud_filelist
   datapath   = prescribed_cloud_datapath
   data_type  = prescribed_cloud_type
   cycle_yr   = prescribed_cloud_cycle_yr
   fixed_ymd  = prescribed_cloud_fixed_ymd
   fixed_tod  = prescribed_cloud_fixed_tod
   num_file_years = prescribed_cloud_num_file_years

   ! Turn on prescribed cloud if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_cloud = .true.

end subroutine prescribed_cloud_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_cloud_adv_orig( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole

    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk, pbuf_get_field, pbuf_set_field

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if( .not. has_prescribed_cloud ) return

    call advance_trcdata( fields, file, state, pbuf2d )

  end subroutine prescribed_cloud_adv_orig

!-------------------------------------------------------------------

  subroutine init_prescribed_cloud_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_cloud', piofile, file )

  end subroutine init_prescribed_cloud_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_cloud_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_cloud_restart

!-------------------------------------------------------------------
  subroutine read_prescribed_cloud_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_cloud', piofile, file )

  end subroutine read_prescribed_cloud_restart
!================================================================================================


  subroutine prescribed_cloud_init(state, pbuf2d)
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
    
!BEH: check to see if state and pbuf2d really need to be passed in.  I didn't use them
!     in the init routine for prescrbed surface fluxes.
    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    !local vars
    type(file_desc_t)   :: fh
    character(len=16)   :: spc_name
    character(len=16)   :: spc_cname
    character(len=16)   :: spc_fname
    character(len=cxx)  :: err_str
    
    integer :: ndx, istat, i, astat, m, n, mm, c
    integer :: grid_id
    integer :: dimlevid, var_id, errcode, dim1id, dim2id, dim1len, dim2len
    integer :: dimbndid, nbnd
    integer :: hdim1_d, hdim2_d    ! model grid size
    real(r8) :: dtime
    logical  :: fixed, cyclical
    
    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------
!++BEH
!    if (aircraft_cnt == 0 ) return
    if ( has_prescribed_cloud ) then
       if ( masterproc ) then
          write(iulog,*) 'now cloud is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    do i = 1,nflds
!       specifier(i) = trim(cloud_name(i))//':'//trim(fld_name(i))
       specifier(i) = trim(cloud_name(i))
       if ( masterproc ) then
          write(iulog,*) 'A specifier:'//specifier(i)
       endif
    end do
!--BEH
    
    
    !------------------------------------------------------------------
    ! For forcing files which has to be on the native grid,dimensions
    ! are set in the following if condition
    !------------------------------------------------------------------
!    if( any( horz_native(:) ) ) then ! reduce horz_native to single boolean
    if( horz_native ) then
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
!++BEH
!       allocate( native_grid_frc_air(aircraft_cnt), stat=astat )
       allocate( native_grid_cloud(nflds), stat=astat )
!--BEH
       if( astat /= 0 ) then 
          write(err_str,*) 'failed to allocate native_grid_cloud array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       end if
    endif
    
    if (masterproc) write(iulog,*) ' '
    
!    if( any( .not. horz_native(:) ) ) then
    if( .not. horz_native ) then
       !-----------------------------------------------------------------------
       !       allocate forcings type array
       !-----------------------------------------------------------------------
!++BEH
!       allocate( forcings_air(aircraft_cnt), stat=astat )
       allocate( forcings_air(nflds), stat=astat )
!--BEH
       if( astat/= 0 ) then
          write(err_str,*) 'failed to allocate forcings_air array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str) 
       end if
    endif
    
    
    !-----------------------------------------------------------------------
    !       setup the forcings_air type array
    !-----------------------------------------------------------------------
!++BEH
!    species_loop : do m = 1,aircraft_cnt
    species_loop : do m = 1,nflds
       
!       spc_name = spc_name_list(m)
       spc_name = specifier(m)

       spc_cname = trim(cloud_name(m))
!       spc_fname = trim(fld_name(m))
       spc_fname = trim(cloud_name(m))
       if ( masterproc ) then
          write(iulog,*) 'spc_name is: '//trim(spc_name)
          write(iulog,*) 'spc_cname is: '//trim(spc_cname)
          write(iulog,*) 'spc_fname is: '//trim(spc_fname)
       endif
!--BEH
       
!       if( horz_native(index_map(m))) then
       if( horz_native ) then
          !-----------------------------------------------------------------------
          !       initialize variables for native grid forcing files
          !-----------------------------------------------------------------------

!++BEH          
          native_grid_cloud(m)%spc_name_ngrd  = spc_name
          native_grid_cloud(m)%spc_cname_ngrd = spc_cname
          native_grid_cloud(m)%spc_fname_ngrd = spc_fname

          fixed    = .false.
          cyclical = .false.
          select case ( data_type )
          case( 'FIXED' )
             fixed = .true.
          case( 'CYCLICAL' )
             cyclical = .true.
             !file%cyc_yr = data_cycle_yr
             ! BEH - should I be setting num_file_years for this?
          case( 'SERIAL' )
             ! Do nothing
          case default 
             write(iulog,*) 'prescribed_cloud: invalid data type: '//trim(data_type)//' file: '//trim(filename)
             write(iulog,*) 'prescribed_cloud: valid data types: SERIAL | CYCLICAL | FIXED '
             call endrun('prescribed_cloud: invalid data type: '//trim(data_type)//' file: '//trim(filename))
          endselect
!--BEH
!++BEH
!          native_grid_cloud(m)%input_file     = trim(datapath)//'/'//trim(spc_fnames(m)) !BEH air_datapath -> datapath
          native_grid_cloud(m)%input_file     = trim(datapath)//'/'//trim(filename)
!--BEH

          native_grid_cloud(m)%initialized    = .false.
          !dtime = 1.0_r8 - 200.0_r8 / 86400.0_r8
          !dtime = -1.0_r8
          dtime  = 0.0_r8
          call native_grid_cloud(m)%time_coord%initialize(trim(adjustl(native_grid_cloud(m)%input_file)), &
               force_time_interp=.true., delta_days=dtime, fixed=fixed, &
               cyclical=cyclical, num_file_years=num_file_years)

          !-----------------------------------------------------------------------
          !       Open file
          !-----------------------------------------------------------------------
          call cam_pio_openfile(fh, trim(adjustl(native_grid_cloud(m)%input_file)), PIO_NOWRITE) 

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
                  ' '//trim(adjustl(native_grid_cloud(m)%input_file))//' '&
                  ' '//errmsg(__FILE__,__LINE__))
          endif
          
          !find if the model and netcdf file has same grid resolution
          call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
          if( dycore_is('SE') )  then
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr) then
                if(dim1len /= hdim1_d ) then !compare model grid length with file's
                   write(err_str,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                        hdim1_d,'), netcdf file is:'//trim(adjustl(native_grid_cloud(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_cloud(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          elseif( dycore_is('LR')) then
             if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim2id)) then !obtain lat dimension of model
                call endrun('failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
                     ' '//trim(adjustl(native_grid_cloud(m)%input_file))//' '//errmsg(__FILE__,__LINE__))
             endif
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr .and. &
                  pio_inquire_dimension(fh, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
                if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
                   write(err_str,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                        hdim1_d,' x ',hdim2_d,'), netcdf file is:'//trim(adjustl(native_grid_cloud(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_cloud(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          else
             call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
          endif
          
          !Find the value of vertical levels in the forcing file
          if( pio_inq_dimid(fh, 'lev', dimlevid) ==  pio_noerr ) then
             if ( pio_inquire_dimension(fh, dimlevid, len =  native_grid_cloud(m)%lev_frc) /=  pio_noerr ) then
                write(err_str,*)'failed to obtain value of "lev" dimension from file:',&
                     trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
             !obtain level bounds needed for vertical interpolation
!++BEH --- don't need lev_bnds since there is no vertical interpolation
!             if( pio_inq_varid(fh, 'lev_bnds', var_id) ==  pio_noerr ) then
!                !get dimension "bound"
!                if( pio_inq_dimid(fh, 'bound', dimbndid) ==  pio_noerr ) then
!                   if ( pio_inquire_dimension(fh, dimbndid, len = nbnd) ==  pio_noerr ) then
!                      !"nbnd" has to be 2 (it is obvious but adding a check here doesn't hurt)
!                      if(nbnd /= 2) then
!                         write(err_str,*)'"bound" should be equal to 2, bound=',nbnd,' in file:', &
!                              trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
!                         call endrun(err_str)
!                      endif
!                      allocate(native_grid_frc_air(m)%lev_bnds(nbnd,native_grid_frc_air(m)%lev_frc))
!                      if (pio_get_var(fh, var_id,native_grid_frc_air(m)%lev_bnds) /=  pio_noerr ) then
!                         write(err_str,*)'failed to read "lev_bnds" variable from file:',&
!                              trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
!                         call endrun(err_str)
!                      endif
!                   else
!                      write(err_str,*)'failed to obtain value of "bound" dimension from file:',&
!                           trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
!                      call endrun(err_str)
!                   endif
!                else
!                   write(err_str,*)'failed to inquire "bound" dimension from file:',&
!                        trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
!                   call endrun(err_str)
!                endif
!             else
!                write(err_str,*)'failed to obtain "lev_bnds" variable from file:',&
!                     trim(adjustl(native_grid_frc_air(m)%input_file)),',',errmsg(__FILE__, __LINE__)
!                call endrun(err_str)                
!             endif
!--BEH
          else
             write(err_str,*)'Dimension "lev" is not found in:',&
                  trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !get units of the data in the forcing file
!BEH === do I need units for these things?  Looks like all the vars have the units attribute
          if(pio_inq_varid( fh, spc_fname, var_id ) == pio_noerr ) then
             if(pio_get_att( fh, var_id, 'units', native_grid_cloud(m)%units) .ne. pio_noerr ) then
                write(err_str,*)'failed to obtain units of variable ',trim(spc_fname),' in &
                     &file:',trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'variable ',trim(spc_fname),' not found in:',trim(adjustl(native_grid_cloud(m)%input_file)), &
                  ',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !close file
          call pio_closefile(fh)
          
          !allocate arrays to store data for interpolation in time
          allocate(native_grid_cloud(m)%native_grid_flds_tslices(pcols, native_grid_cloud(m)%lev_frc, &
               begchunk:endchunk,2), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_cloud(',m,')%native_grid_flds_tslices array;&
                  error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !allocate arrays to hold data before the vertical interpolation
          allocate(native_grid_cloud(m)%native_grid_flds(pcols, native_grid_cloud(m)%lev_frc,begchunk:endchunk), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_cloud(',m,')%native_grid_flds array; error = ',&
                  astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !get pbuf index to store the field in pbuf
          native_grid_cloud(m)%pbuf_ndx = pbuf_get_index(spc_cname,errcode)
          if(errcode < 0 ) then
             write(err_str,*)'failed to get pbuf index for specie:',spc_cname,' errorcode is:',errcode,',',errmsg(__FILE__, __LINE__)
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
          !         forcings_air(m)%file%curr_filename    = spc_fnames(m)
          forcings_air(m)%filename         = spc_fnames(m)
       endif

!++BEH --- I don't need these, but may want to keep them for debugging purposes
!       call addfld( trim(spc_cname)//'_debug', (/ 'lev' /), 'A',  '1/s',     &
!            'aircraft emission '//trim(spc_cname) )
!       call add_default( trim(spc_cname)//'_debug', 1, ' ' )
!--BEH
       
    end do species_loop
    
!++BEH not sure what to do this, so for now commenting out
!    if (masterproc) then
!       !-----------------------------------------------------------------------
!       !            diagnostics
!       !-----------------------------------------------------------------------
!       write(iulog,*) ' '
!       write(iulog,*) 'aircraft_emit_init: diagnostics'
!       write(iulog,*) ' '
!       write(iulog,*) 'aircraft_emit timing specs'
!       write(iulog,*) 'type = ',air_type
!       write(iulog,*) ' '
!       write(iulog,*) 'there are ',aircraft_cnt,' species of aircraft emission'
!       do m = 1,aircraft_cnt
!          write(iulog,*) ' '          
!          write(iulog,*) 'forcing type ',m
!          write(iulog,*) 'species = ',spc_name_list(m)
!          write(iulog,*) 'filelist= ',spc_flist(m)
!       end do
!       write(iulog,*) ' '
!    endif
!--BEH
    
    
    !------------------------------------------------------------------
    !       Initialize the aircraft file processing
    !------------------------------------------------------------------
!++BEH
!    do m=1,aircraft_cnt
    do m=1,nflds
!--BEH

       number_flds = 0
!       if(horz_native(index_map(m))) then
       if(horz_native) then
          if (associated(native_grid_cloud(m)%native_grid_flds_tslices)) &
               number_flds = 1
          !read the forcing file once to initialize variable including time cordinate
          call advance_native_grid_data( native_grid_cloud(m) )
          native_grid_cloud(m)%initialized = .true.
       else
          allocate (forcings_air(m)%file%in_pbuf(size(forcings_air(m)%sectors)))
          forcings_air(m)%file%in_pbuf(:) = .true.
          if (associated(forcings_air(m)%fields)) number_flds = size( forcings_air(m)%fields )
          
!++BEH
!          call trcdata_init( forcings_air(m)%sectors, forcings_air(m)%filename, forcings_air(m)%filelist, air_datapath, &
!               forcings_air(m)%fields, forcings_air(m)%file, rmv_file, 0, 0, 0, air_type)
          call trcdata_init( forcings_air(m)%sectors, forcings_air(m)%filename, forcings_air(m)%filelist, datapath, &
               forcings_air(m)%fields, forcings_air(m)%file, rmv_file, 0, 0, 0, air_type)
!--BEH
       endif
       
       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(err_str,*) 'There are no aircraft aerosols ',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
       end if
    end do
  end subroutine prescribed_cloud_init

  subroutine prescribed_cloud_adv( state, pbuf2d)
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
    
    integer             :: ind, c, ncol, i, caseid, m, pbuf_ndx
!    real(r8)            :: to_mmr(pcols,pver)
    real(r8), pointer   :: tmpptr(:,:)
!    real(r8), pointer   :: tmpptr(:)
    character(len = cs) :: units_spc
    
    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------
!++BEH
!    if (aircraft_cnt == 0 ) return
    if( .not. has_prescribed_cloud ) return
!--BEH
    call t_startf('All_aircraft_emit_adv')
    
    !-------------------------------------------------------------------
    !    For each field, read more data if needed and interpolate it to the current model time
    !-------------------------------------------------------------------
!++BEH
!    do m = 1, aircraft_cnt
    do m = 1,nflds
!--BEH
       
!       if(horz_native(index_map(m))) then ! if horizontal grid is native
       if (horz_native) then
          units_spc = native_grid_cloud(m)%units
          pbuf_ndx  = native_grid_cloud(m)%pbuf_ndx
          
          !read in next time slice (if needed) and interpolate in time
          !following call just reads in time slices in horizontal
          !vertical interpolation is done in the next call
          call advance_native_grid_data( native_grid_cloud(m) )
 
          !do vertical interpolation
          
          !following call needs state to get ncol for each chunk and
          !pbuf for storing the interpolated (time and vertically)
          !field in pbuf
!++BEH --- don't do vertical interpolation anymore
!          call vert_interp( state, pbuf_ndx, native_grid_frc_air(m), pbuf2d)
!--BEH

!++BEH ---->  Testing putting data into pbuf
          !Need to store the native_grid_cloud data into the physics buffer.
          !$OMP PARALLEL DO PRIVATE (C, NCOL, TMPPTR, PBUF_CHNK)
          do c = begchunk, endchunk
             ncol = state(c)%ncol
             pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
             call pbuf_get_field(pbuf_chnk, pbuf_ndx, tmpptr)
             tmpptr(:ncol,:) = native_grid_cloud(m)%native_grid_flds(:,:,c)
          enddo
!--BEH
          
       else          
          units_spc = forcings_air(m)%fields(i)%units
          pbuf_ndx  = forcings_air(m)%fields(i)%pbuf_ndx
          call advance_trcdata( forcings_air(m)%fields, forcings_air(m)%file, state, pbuf2d)
       endif

       !-------------------------------------------------------------------
       !    set the tracer fields with the correct units
       !-------------------------------------------------------------------
!++BEH I do not need to adjust units for the cloud fields.  All this should go.
!       do i = 1,number_flds
!          
!          !initialize caseid so that it is not used inadvertantly
!          caseid = huge_int
!          ind    = index_map(i)
!          
!          !only assign valid integer if we need unit conversion
!          if ( convert_to_mmr(ind) ) then
!             ! C.-C. Chen, adding case 4  for kg/sec
!             select case ( to_lower(trim(units_spc(:GLC(units_spc)))) )
!             case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
!                caseid = 1
!             case ('kg/kg','mmr')
!                caseid = 2
!             case ('mol/mol','mole/mole','vmr','fraction')
!                caseid = 3
!             case ('kg/kg/sec')
!                caseid = 4
!             case default
!                print*, 'aircraft_emit_adv: units = ',trim(units_spc) ,' are not recognized'
!                call endrun('aircraft_emit_adv: units are not recognized '//errmsg(__FILE__, __LINE__))
!             end select
!          endif
!          
!          !$OMP PARALLEL DO PRIVATE (C, NCOL, TO_MMR, tmpptr, pbuf_chnk)
!          do c = begchunk,endchunk
!             ncol = state(c)%ncol
!             
!             !initialize to_mmr to 1.0 for cases where unit conversion is not required
!             !(i.e., caseid is not assigned an integer value)
!             to_mmr(:ncol,:) = 1.0_r8
!             
!             if (caseid == 1) then ! change it to select-case?
!                to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
!             elseif(caseid == 2) then
!                to_mmr(:ncol,:) = 1.0_r8
!             elseif(caseid == 3) then
!                to_mmr(:ncol,:) = molmass(ind)/mwdry
!             elseif(caseid == 4) then
!                to_mmr(:ncol,:) = 1.0_r8
!             endif
!             
!             pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
!             
!             call pbuf_get_field(pbuf_chnk, pbuf_ndx , tmpptr )
!             
!             tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)
!             
!             call outfld( trim(spc_name_list(m)), &
!                  tmpptr, ncol, state(c)%lchnk )
!          enddo
!       enddo
!--BEH
    enddo
    
    call t_stopf('All_aircraft_emit_adv')
  end subroutine prescribed_cloud_adv

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
    character(len=cs) :: spc_cname
    character(len=cs) :: spc_fname
    
    logical  :: read_data
    integer  :: indx2_pre_adv
    logical  :: found

    !obtain name of the specie
    spc_name  = native_grid_strct%spc_name_ngrd
    spc_cname = native_grid_strct%spc_cname_ngrd
    spc_fname = native_grid_strct%spc_fname_ngrd
    
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
          call infld(trim(spc_fname), fh, dim1name, dim2name, 'lev',&
               1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
               native_grid_strct%native_grid_flds_tslices(:,:,:,1), found, &
               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          if (.not. found) then
             call endrun(trim(spc_fname) // ' not found '//errmsg(__FILE__,__LINE__))
          endif
       endif
       
       ! read time level 2
       call infld(trim(spc_fname), fh, dim1name, dim2name, 'lev',&
            1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
            native_grid_strct%native_grid_flds_tslices(:,:,:,2), found, &
            gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(2))
       
       if (.not. found) then
          call endrun(trim(spc_fname) // ' not found '//errmsg(__FILE__,__LINE__))
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

end module prescribed_cloud
