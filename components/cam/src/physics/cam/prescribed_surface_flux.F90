
module prescribed_sfc_flux

!------------------------------------------------------------------------------------------------
! Purpose:
! Reads in surface fluxes from a file, puts them into the physics buffer.
!
! Author: Bryce Harrop
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8, cx => shr_kind_cx, cl => shr_kind_cl, &
                              cs => shr_kind_cs, cxx => shr_kind_cxx
  use ppgrid,           only: begchunk, endchunk, pcols
  use cam_abortutils,   only: endrun
  use dycore,           only: dycore_is
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use input_data_utils, only: time_coordinate
  use cam_pio_utils,    only: cam_pio_openfile
  use pio,              only: file_desc_t, pio_nowrite, pio_closefile, pio_inq_dimid, pio_bcast_error, &
                              pio_seterrorhandling, pio_noerr, pio_inquire_dimension
  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog
  
  implicit none

! public type

  public presc_sfc_flux_type

! public interface
  public presc_sfc_flux_readnl
  public presc_sfc_flux_register
  public presc_sfc_flux_init
  public presc_sfc_flux_adv
  public presc_sfc_flux_overwrite
  public has_presc_sfc_flux

  logical :: has_presc_sfc_flux = .false.

! private data

  private
 
!++BEH
   ! Need to have multiple fields read in:
   ! shf, cflx, lhf, lwup, asdir, aldir, asdif, aldif
   integer          , parameter :: nflds              = 8
!   character(len=16), parameter :: sflx_name(nflds)   = (/'shf_in', 'cflx_in', 'lhf_in', 'lwup_in', 'asdir_in', 'aldir_in', 'asdif_in', 'aldif_in' /)
   character(len=16), parameter :: sflx_name(nflds)   = (/ 'SHFLX', 'QFLX', 'LHFLX', 'FLUS', 'ASDIR', 'ALDIR', 'ASDIF', 'ALDIF' /)
   character(len=16), parameter :: sflx_pname(nflds)  = (/ 'p_SHFLX', 'p_QFLX', 'p_LHFLX', 'p_FLUS', 'p_ASDIR', 'p_ALDIR', 'p_ASDIF', 'p_ALDIF' /)
   logical, parameter           :: sflx_dtimes(nflds) = (/ .false., .false., .false., .true., .true., .true., .true., .true. /)
   character(len=256)           :: filename           = ' '
   character(len=256)           :: filelist           = ' '
   character(len=256)           :: datapath           = ' '
   character(len=32)            :: data_type          = 'CYCLICAL'
   real(r8)                     :: num_file_years     = 0._r8
   real(r8)                     :: input_dtime        = 0._r8
   logical                      :: rmv_file           = .false.
   integer                      :: cycle_yr           = 0
   integer                      :: fixed_ymd          = 0
   integer                      :: fixed_tod          = 0
 
   logical, parameter           :: horz_native        = .true. ! they will all have the same behavior if they
                                                               ! all come from the same file
   logical                      :: dimnames_set       = .false. 
   integer                      :: number_flds
   character(len=16)            :: spc_name_list(nflds)
   character(len=cl)            :: spc_flist(nflds), spc_fnames(nflds)
   character(len=8)             :: dim1name, dim2name
   character(len=24)            :: air_type           = 'CYCLICAL_LIST' ! 'CYCLICAL_LIST'
!--BEH

!--------------------------------------------------------------------------------------------------
type :: presc_sfc_flux_type          
!++BEH
   !(pcols,begchunk:endchunk,2)
!B3   real(r8), pointer, dimension(:,:,:) :: native_grid_flds_tslices
   real(r8), pointer, dimension(:,:,:,:) :: native_grid_flds_tslices
   integer                               :: lev_frc
   
   !Data structure to store data after time interpolation from two time samples
   !(pcols,begchunk:endchunk)
!B3   real(r8), pointer, dimension(:,:)   :: native_grid_flds
   real(r8), pointer, dimension(:,:,:)   :: native_grid_flds
!--BEH

   !Forcing file name
   character(len=cl) :: input_file  

   !Data structure to keep track of time
   type(time_coordinate) :: time_coord

   !specie name
   character(len=cl)     :: spc_name_ngrd
   character(len=cl)     :: spc_pname_ngrd
   
   !logical to control first data read
   logical               :: initialized

   !Units of forcing data
   character(len=cl)     :: units

   !pbuf index to store read in data in pbuf
   integer               :: pbuf_ndx = -1
 
end type presc_sfc_flux_type

type(presc_sfc_flux_type), dimension(nflds) :: natgrid_sfc_flux_in

!===============================================================================
contains
!===============================================================================

subroutine presc_sfc_flux_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'presc_sfc_flux_readnl'

   character(len=256) :: presc_sfc_flux_file
   character(len=256) :: presc_sfc_flux_filelist
   character(len=256) :: presc_sfc_flux_datapath
   character(len=32)  :: presc_sfc_flux_type
   logical            :: presc_sfc_flux_rmfile
   integer            :: presc_sfc_flux_cycle_yr
   integer            :: presc_sfc_flux_fixed_ymd
   integer            :: presc_sfc_flux_fixed_tod
   real(r8)           :: presc_sfc_flux_num_file_years
   real(r8)           :: presc_sfc_flux_input_dtime

   namelist /presc_sfc_flux_nl/ &
      presc_sfc_flux_file,      &
      presc_sfc_flux_filelist,  &
      presc_sfc_flux_datapath,  &
      presc_sfc_flux_type,      &
      presc_sfc_flux_rmfile,    &
      presc_sfc_flux_cycle_yr,  &
      presc_sfc_flux_fixed_ymd, &
      presc_sfc_flux_fixed_tod, &
      presc_sfc_flux_num_file_years, &
      presc_sfc_flux_input_dtime
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   presc_sfc_flux_file      = filename
   presc_sfc_flux_filelist  = filelist
   presc_sfc_flux_datapath  = datapath
   presc_sfc_flux_type      = data_type
   presc_sfc_flux_rmfile    = rmv_file
   presc_sfc_flux_cycle_yr  = cycle_yr
   presc_sfc_flux_fixed_ymd = fixed_ymd
   presc_sfc_flux_fixed_tod = fixed_tod
   presc_sfc_flux_num_file_years = num_file_years
   presc_sfc_flux_input_dtime    = input_dtime

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'presc_sfc_flux_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, presc_sfc_flux_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(presc_sfc_flux_file,     len(presc_sfc_flux_file),     mpichar, 0, mpicom)
   call mpibcast(presc_sfc_flux_filelist, len(presc_sfc_flux_filelist), mpichar, 0, mpicom)
   call mpibcast(presc_sfc_flux_datapath, len(presc_sfc_flux_datapath), mpichar, 0, mpicom)
   call mpibcast(presc_sfc_flux_type,     len(presc_sfc_flux_type),     mpichar, 0, mpicom)
   call mpibcast(presc_sfc_flux_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(presc_sfc_flux_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(presc_sfc_flux_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(presc_sfc_flux_fixed_tod,1, mpiint,  0, mpicom)
   call mpibcast(presc_sfc_flux_num_file_years, 1, mpir8, 0, mpicom)
   call mpibcast(presc_sfc_flux_input_dtime,    1, mpir8, 0, mpicom)
#endif

   ! Update module variables with user settings.
   filename   = presc_sfc_flux_file
   filelist   = presc_sfc_flux_filelist
   datapath   = presc_sfc_flux_datapath
   data_type  = presc_sfc_flux_type
   rmv_file   = presc_sfc_flux_rmfile
   cycle_yr   = presc_sfc_flux_cycle_yr
   fixed_ymd  = presc_sfc_flux_fixed_ymd
   fixed_tod  = presc_sfc_flux_fixed_tod
   num_file_years = presc_sfc_flux_num_file_years
   input_dtime    = presc_sfc_flux_input_dtime

   ! Turn on prescribed surface fluxes if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_presc_sfc_flux = .true.

end subroutine presc_sfc_flux_readnl



subroutine presc_sfc_flux_register()
!B3   use ppgrid,         only: pcols
   use ppgrid,         only: pcols, pver
   use physics_buffer, only: pbuf_add_field, dtype_r8

   integer :: i,idx

   if (has_presc_sfc_flux) then
      do i = 1,nflds
!         call pbuf_add_field(sflx_name(i), 'physpkg', dtype_r8, (/pcols/), idx)
!B3         call pbuf_add_field(sflx_pname(i), 'physpkg', dtype_r8, (/pcols/), idx)
         call pbuf_add_field(sflx_pname(i), 'physpkg', dtype_r8, (/pcols,pver/), idx)
      enddo
   endif

end subroutine presc_sfc_flux_register



subroutine presc_sfc_flux_init()

!-------------------------------------------------------------------------------
! Initialize presc_sfc_flux_type instance
!   including initial read of input and interpolation to the current timestep
!-------------------------------------------------------------------------------

   use cam_history,      only: addfld, add_default, horiz_only
   use ppgrid,           only: begchunk, endchunk, pcols
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_index
   use physics_buffer,   only: physics_buffer_desc
   use cam_grid_support, only: cam_grid_id, cam_grid_check
   use cam_grid_support, only: cam_grid_get_dim_names
   use dyn_grid,         only: get_horiz_grid_dim_d
   use dycore,           only: dycore_is
   use cam_pio_utils,    only: cam_pio_openfile
   use pio,              only: file_desc_t, pio_nowrite, pio_closefile, & 
                               pio_inq_dimid, pio_bcast_error, pio_seterrorhandling, & 
                               pio_noerr, pio_inquire_dimension, pio_get_att, & 
                               pio_inq_varid, pio_get_var

   implicit none

   ! Local variables
!++BEH
   type(file_desc_t)  :: fh
   character(len=16)  :: spc_name
   character(len=16)  :: spc_pname
   character(len=cxx) :: err_str
    
   integer            :: ndx, istat, i, astat, m, n, mm, c
   integer            :: dimbndid, nbnd, var_id, errcode
   logical            :: fixed, cyclical
!--BEH
   integer            :: grid_id, dim1len, dim2len, dim1id, dim2id ! netcdf file ids and sizes
   integer            :: hdim1_d, hdim2_d    ! model grid size
   real(r8)           :: dtime
!B3++B3
   integer            :: dimlevid
!B3--B3
   !----------------------------------------------------------------------------

!++BEH
   if ( has_presc_sfc_flux ) then
      if ( masterproc ) then
         write(iulog,*) 'Prescribed surface fluxes are in: '//trim(filename)
      endif
   else
      return
   endif

   do i = 1,nflds
      if ( masterproc ) then
         write(iulog,*) 'Overwriting flux: '//trim(sflx_name(i))
      endif
   end do
!--BEH

   if (.not. dimnames_set) then
      grid_id = cam_grid_id('physgrid')
      if (.not. cam_grid_check(grid_id)) then
         call endrun('ERROR: no "physgrid" grid:'//errmsg(__FILE__,__LINE__))
      endif
      !dim1name and dim2name are populated here with the grid dimension the model is running on (e.g. ne30, lat, lon etc.)
      !For SE grid, dim1name = dim2name = "ncol"
      !For FV grid, dim1name = lon, dim2name = lat
      call cam_grid_get_dim_names(grid_id, dim1name, dim2name) 
      dimnames_set = .true.
   end if

!++BEH -- don't need this anymore.  natgrid_sfc_flux_in will always be nflds in dimension, so specify above.
!   allocate( natgrid_sfc_flux_in(nflds), stat=astat )
!   if( astat /= 0 ) then 
!      write(err_str,*) 'failed to allocate natgrid_sfc_flux_in array; error = ',astat,',',errmsg(__FILE__, __LINE__)
!      call endrun(err_str)
!   end if
!--BEH

   flux_loop : do m = 1,nflds
      spc_name  = trim(sflx_name(m))
      spc_pname = trim(sflx_pname(m))
      if ( masterproc ) then
         write(iulog,*) 'spc_name is: '//trim(spc_name)
      endif

      ! Set species name
      natgrid_sfc_flux_in(m)%spc_name_ngrd  = spc_name
      natgrid_sfc_flux_in(m)%spc_pname_ngrd = spc_pname

      ! Set logicals for how to read flux files
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
         write(iulog,*) 'prescribed_surface_flux: invalid data type: '//trim(data_type)//' file: '//trim(filename)
         write(iulog,*) 'prescribed_surface_flux: valid data types: SERIAL | CYCLICAL | FIXED '
         call endrun('prescribed_surface_flux: invalid data type: '//trim(data_type)//' file: '//trim(filename))
      endselect

      ! Set input file
      natgrid_sfc_flux_in(m)%input_file     = trim(datapath)//'/'//trim(filename)

      ! Set initialized state to false, to force the initialization procedure
      natgrid_sfc_flux_in(m)%initialized    = .false.

      ! dtime is the offset time between the model and file.  Modify if not synchronized properly.
      if (sflx_dtimes(m)) then
!         dtime = (1.0_r8 / 24._r8) ! Read from the previous time step (hour) for dtime
         dtime = input_dtime / 86400. ! This is the timestep length for the input data converted from seconds to days
         !BEH -- should 86400 be converted to some named constant?
      else
         dtime = 0.0_r8
      end if

      ! Initialize the time coordinate of the native grid variable
      call natgrid_sfc_flux_in(m)%time_coord%initialize(trim(adjustl(natgrid_sfc_flux_in(m)%input_file)), & 
                                                        force_time_interp=.true., delta_days=dtime, &
                                                        fixed=fixed, cyclical=cyclical, &
                                                        num_file_years=num_file_years)

      ! ------------------------------------------------------------
      !                          Open file
      ! ------------------------------------------------------------

      call cam_pio_openfile(fh, trim(adjustl(natgrid_sfc_flux_in(m)%input_file)), PIO_NOWRITE)
      
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
                     ' '//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))//' '&
                     ' '//errmsg(__FILE__,__LINE__))
      endif
          
      !find if the model and netcdf file has same grid resolution
      call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
      if( dycore_is('SE') )  then
         if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr) then
            if(dim1len /= hdim1_d ) then !compare model grid length with file's
               write(err_str,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                        hdim1_d,'), netcdf file is:'//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))
               call endrun(err_str//errmsg(__FILE__,__LINE__))
            endif
         else
            call endrun('failed while inquiring dimensions of file:'//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
         endif
      elseif( dycore_is('LR')) then
         if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim2id)) then !obtain lat dimension of model
            call endrun('failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
                     ' '//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))//' '//errmsg(__FILE__,__LINE__))
         endif
         if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr .and. &
                  pio_inquire_dimension(fh, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
            if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
               write(err_str,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                        hdim1_d,' x ',hdim2_d,'), netcdf file is:'//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))
               call endrun(err_str//errmsg(__FILE__,__LINE__))
            endif
         else
            call endrun('failed while inquiring dimensions of file:'//trim(adjustl(natgrid_sfc_flux_in(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
         endif
      else
         call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
      endif

!B3++B3
          !Find the value of vertical levels in the forcing file
          if( pio_inq_dimid(fh, 'lev', dimlevid) ==  pio_noerr ) then
             if ( pio_inquire_dimension(fh, dimlevid, len =  natgrid_sfc_flux_in(m)%lev_frc) /=  pio_noerr ) then
                write(err_str,*)'failed to obtain value of "lev" dimension from file:',&
                     trim(adjustl(natgrid_sfc_flux_in(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'Dimension "lev" is not found in:',&
                  trim(adjustl(natgrid_sfc_flux_in(m)%input_file)),',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
!B3--B3
          
          !get units of the data in the forcing file
!BEH === do I need units for these things?  Looks like all the vars have the units attribute
      if(pio_inq_varid( fh, spc_name, var_id ) == pio_noerr ) then
         if(pio_get_att( fh, var_id, 'units', natgrid_sfc_flux_in(m)%units) .ne. pio_noerr ) then
            write(err_str,*)'failed to obtain units of variable ',trim(spc_name),' in &
                 &file:',trim(adjustl(natgrid_sfc_flux_in(m)%input_file)),',',errmsg(__FILE__, __LINE__)
            call endrun(err_str)
         endif
      else
         write(err_str,*)'variable ',trim(spc_name),' not found in:',trim(adjustl(natgrid_sfc_flux_in(m)%input_file)), &
                  ',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
      endif
          
      !close file
      call pio_closefile(fh)
          
      !allocate arrays to store data for interpolation in time
!B3      allocate(natgrid_sfc_flux_in(m)%native_grid_flds_tslices(pcols, begchunk:endchunk,2), stat=astat )
      allocate(natgrid_sfc_flux_in(m)%native_grid_flds_tslices(pcols, natgrid_sfc_flux_in(m)%lev_frc, begchunk:endchunk,2), stat=astat )
      if( astat/= 0 ) then
         write(err_str,*) 'failed to allocate natgrid_sfc_flux_in(',m,')%native_grid_flds_tslices array;&
                  error = ',astat,',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
      endif
          
      !allocate arrays to hold data before the vertical interpolation
!B3      allocate(natgrid_sfc_flux_in(m)%native_grid_flds(pcols, begchunk:endchunk), stat=astat )
      allocate(natgrid_sfc_flux_in(m)%native_grid_flds(pcols, natgrid_sfc_flux_in(m)%lev_frc, begchunk:endchunk), stat=astat )
      if( astat/= 0 ) then
         write(err_str,*) 'failed to allocate natgrid_sfc_flux_in(',m,')%native_grid_flds array; error = ',&
                  astat,',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
      endif
          
      !get pbuf index to store the field in pbuf
      natgrid_sfc_flux_in(m)%pbuf_ndx = pbuf_get_index(spc_pname,errcode)
      if(errcode < 0 ) then
         write(err_str,*)'failed to get pbuf index for specie:',spc_pname,' errorcode is:',errcode,',',errmsg(__FILE__, __LINE__)
         call endrun(err_str)
      endif

      call addfld( 'INFLX_'//trim(spc_name), horiz_only, 'A',  'W/m2',     &
            'Input flux for '//trim(spc_name) )

   end do flux_loop


    !------------------------------------------------------------------
    !       Initialize the aircraft file processing
    !------------------------------------------------------------------
!++BEH
!    do m=1,aircraft_cnt
    do m=1,nflds
!--BEH

       number_flds = 0
!       if(horz_native(index_map(m))) then
       if (associated(natgrid_sfc_flux_in(m)%native_grid_flds_tslices)) &
            number_flds = 1
       !read the forcing file once to initialize variable including time cordinate
       call advance_native_grid_data( natgrid_sfc_flux_in(m) )
       natgrid_sfc_flux_in(m)%initialized = .true.
       
       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(err_str,*) 'There are no surface fluxes or properties ',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
       end if
    end do

end subroutine presc_sfc_flux_init



subroutine presc_sfc_flux_adv(state, pbuf2d)

!-------------------------------------------------------------------------------
! Advance the contents of a presc_sfc_flux_type instance
!   including reading new data, if necessary
!-------------------------------------------------------------------------------

   use perf_mod,       only: t_startf, t_stopf
   use physics_types,  only: physics_state
   use ppgrid,         only: begchunk, endchunk
   use ppgrid,         only: pcols, pver
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
   use cam_history,    only: outfld

   implicit none

   type(physics_state), intent(in)    :: state(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   ! Local variables
   integer                            :: c, ncol, m, pbuf_ndx
   character(len=cs)                  :: units_spc
!B3   real(r8), pointer, dimension(:)    :: tmpptr_native_grid
   real(r8), pointer, dimension(:,:)  :: tmpptr_native_grid
   character(len=16)                  :: spc_name

   !----------------------------------------------------------------------------

   !------------------------------------------------------------------
   ! Return if no surface flux file
   !------------------------------------------------------------------
!++BEH
!    if (aircraft_cnt == 0 ) return
    if( .not. has_presc_sfc_flux ) return
!--BEH
    call t_startf('All_prescribed_surface_fluxes')

    !-------------------------------------------------------------------
    !    For each field, read more data if needed and interpolate it to the current model time
    !-------------------------------------------------------------------

    do m = 1,nflds
       units_spc = natgrid_sfc_flux_in(m)%units
       pbuf_ndx  = natgrid_sfc_flux_in(m)%pbuf_ndx
       spc_name  = natgrid_sfc_flux_in(m)%spc_name_ngrd
       
       !read in next time slice (if needed) and interpolate in time
       !following call just reads in time slices in horizontal
       !vertical interpolation is done in the next call
       call advance_native_grid_data( natgrid_sfc_flux_in(m) )

       !$OMP PARALLEL DO PRIVATE (C, NCOL, TMPPTR_NATIVE_GRID, PBUF_CHNK)
       do c = begchunk, endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_chnk, pbuf_ndx, tmpptr_native_grid)
!B3          tmpptr_native_grid(:ncol) = natgrid_sfc_flux_in(m)%native_grid_flds(:,c)
          tmpptr_native_grid(:ncol,:) = natgrid_sfc_flux_in(m)%native_grid_flds(:,:,c)

!B3          call outfld( 'INFLX_'//trim(spc_name), tmpptr_native_grid(:ncol), &
!B3                       ncol, state(c)%lchnk )
          call outfld( 'INFLX_'//trim(spc_name), tmpptr_native_grid(:ncol, natgrid_sfc_flux_in(m)%lev_frc), &
                       ncol, state(c)%lchnk )
       enddo

    enddo



    call t_stopf('All_prescribed_surface_fluxes')

end subroutine presc_sfc_flux_adv



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

    !arguments
    type(presc_sfc_flux_type), intent (inout) :: native_grid_strct 
    
    !local vars
    type(file_desc_t) :: fh
    character(len=cs) :: spc_name
    logical           :: read_data
    integer           :: indx2_pre_adv
    logical           :: found

    !obtain name of the specie
    spc_name  = native_grid_strct%spc_name_ngrd
    
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
!B3          native_grid_strct%native_grid_flds_tslices(:,:,1) = native_grid_strct%native_grid_flds_tslices(:,:,2)
          native_grid_strct%native_grid_flds_tslices(:,:,:,1) = native_grid_strct%native_grid_flds_tslices(:,:,:,2)
       else
          !NOTE: infld call doesn't do any interpolation in space, it just reads in the data
          if(masterproc) then
             write(iulog, *) 'First infld call in presc_sfc_flux'
             write(iulog, *) 'dim1name is ', dim1name, ' dim2name is ', dim2name
          endif
!B3          call infld(trim(spc_name), fh, dim1name, dim2name, &
!B3               1, pcols, begchunk, endchunk, &
!B3               native_grid_strct%native_grid_flds_tslices(:,:,1), found, &
!B3               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          call infld(trim(spc_name), fh, dim1name, dim2name, 'lev', &
               1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
               native_grid_strct%native_grid_flds_tslices(:,:,:,1), found, &
               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          if (.not. found) then
             call endrun(trim(spc_name) // ' not found '//errmsg(__FILE__,__LINE__))
          endif
       endif
       
       ! read time level 2
       if(masterproc) then
          write(iulog, *) 'Second infld call in presc_sfc_flux'
       endif
!B3       call infld(trim(spc_name), fh, dim1name, dim2name, &
!B3            1, pcols, begchunk, endchunk, &
!B3            native_grid_strct%native_grid_flds_tslices(:,:,2), found, &
!B3            gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(2))
       call infld(trim(spc_name), fh, dim1name, dim2name, 'lev', &
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
!B3       native_grid_strct%native_grid_flds(:,:) = native_grid_strct%native_grid_flds_tslices(:,:,1)
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1)
    else
!B3       native_grid_strct%native_grid_flds(:,:) = native_grid_strct%native_grid_flds_tslices(:,:,1) + &
!B3            native_grid_strct%time_coord%wghts(2) * (native_grid_strct%native_grid_flds_tslices(:,:,2) - &
!B3            native_grid_strct%native_grid_flds_tslices(:,:,1))
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1) + &
            native_grid_strct%time_coord%wghts(2) * (native_grid_strct%native_grid_flds_tslices(:,:,:,2) - &
            native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    endif

!++BEH a test for zeros, write out useful stuff
    if ( all(native_grid_strct%native_grid_flds(:,:,:) .eq. 0.0_r8) ) then
       call endrun(trim(spc_name) // ' is all zeros.  Oh no!'//errmsg(__FILE__,__LINE__))
    endif
!--BEH
    
end subroutine advance_native_grid_data

subroutine presc_sfc_flux_overwrite( cam_in, pbuf2d )
!--------------------------------------------
! Based on get_met_srf2
!--------------------------------------------
    use camsrfexch,          only: cam_in_t    
    use phys_grid,           only: get_ncols_p
    use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, &
                                   pbuf_get_chunk, pbuf_get_index

    implicit none

    type(cam_in_t), intent(inout), dimension(begchunk:endchunk) :: cam_in
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer   :: c, ncol, err
    integer   :: index_shf,   index_cflx,  index_lhf,   index_lwup
    integer   :: index_asdir, index_aldir, index_asdif, index_aldif
    real(r8), pointer, dimension(:) :: shf, cflx, lhf, lwup, asdir, aldir, asdif, aldif

! BEH Add in some nudging options
    logical, parameter  :: nudge_turbulent_fluxes = .true.
!    real(r8), parameter :: nts = 3600._r8
    real(r8), parameter :: nts = 2._r8
    real(r8)            :: factor
! BEH NEED TO INCLUDE MODEL TIMESTEP IN FACTOR

    if ( .not. has_presc_sfc_flux ) return

    do c = begchunk, endchunk
       pbuf_chnk  => pbuf_get_chunk(pbuf2d, c)

       index_shf   = pbuf_get_index('p_SHFLX')
       index_cflx  = pbuf_get_index('p_QFLX')
       index_lhf   = pbuf_get_index('p_LHFLX')
       index_lwup  = pbuf_get_index('p_FLUS')
       index_asdir = pbuf_get_index('p_ASDIR')
       index_aldir = pbuf_get_index('p_ALDIR')
       index_asdif = pbuf_get_index('p_ASDIF')
       index_aldif = pbuf_get_index('p_ALDIF')

       call pbuf_get_field(pbuf_chnk, index_shf,   shf)
       call pbuf_get_field(pbuf_chnk, index_cflx,  cflx)
       call pbuf_get_field(pbuf_chnk, index_lhf,   lhf)
       call pbuf_get_field(pbuf_chnk, index_lwup,  lwup)
       call pbuf_get_field(pbuf_chnk, index_asdir, asdir)
       call pbuf_get_field(pbuf_chnk, index_aldir, aldir)
       call pbuf_get_field(pbuf_chnk, index_asdif, asdif)
       call pbuf_get_field(pbuf_chnk, index_aldif, aldif)

!    do c=begchunk,endchunk
       ncol = get_ncols_p(c)

       if (nudge_turbulent_fluxes) then
          factor = 1._r8 / nts !BEH: THIS IS WRONG.  IT SHOULD BE MODEL TIMESTEP / NUDGING TIME SCALE
!          factor = 0.95_r8
!          factor = 0.999_r8
!          factor = 0.1_r8
          cam_in(c)%shf(:ncol)    = cam_in(c)%shf(:ncol)    + factor * (shf(:ncol)  - cam_in(c)%shf(:ncol))
          cam_in(c)%cflx(:ncol,1) = cam_in(c)%cflx(:ncol,1) + factor * (cflx(:ncol) - cam_in(c)%cflx(:ncol,1))
          cam_in(c)%lhf(:ncol)    = cam_in(c)%lhf(:ncol)    + factor * (lhf(:ncol)  - cam_in(c)%lhf(:ncol))
       else
          cam_in(c)%shf(:ncol)    = shf(:ncol)
          cam_in(c)%cflx(:ncol,1) = cflx(:ncol)
          cam_in(c)%lhf(:ncol)    = lhf(:ncol)
       end if
       cam_in(c)%lwup(:ncol)   = lwup(:ncol)
       cam_in(c)%asdir(:ncol)  = asdir(:ncol)
       cam_in(c)%aldir(:ncol)  = aldir(:ncol)
       cam_in(c)%asdif(:ncol)  = asdif(:ncol)
       cam_in(c)%aldif(:ncol)  = aldif(:ncol)
!!       cam_in(c)%shf(:ncol)    = natgrid_sfc_flux_in(1)%native_grid_flds(:ncol, c)
!!       cam_in(c)%cflx(:ncol,1) = natgrid_sfc_flux_in(2)%native_grid_flds(:ncol, c)
!!       cam_in(c)%lhf(:ncol)    = natgrid_sfc_flux_in(3)%native_grid_flds(:ncol, c)
!!       cam_in(c)%lwup(:ncol)   = natgrid_sfc_flux_in(4)%native_grid_flds(:ncol, c)
!!       cam_in(c)%asdir(:ncol)  = natgrid_sfc_flux_in(5)%native_grid_flds(:ncol, c)
!!       cam_in(c)%aldir(:ncol)  = natgrid_sfc_flux_in(6)%native_grid_flds(:ncol, c)
!!       cam_in(c)%asdif(:ncol)  = natgrid_sfc_flux_in(7)%native_grid_flds(:ncol, c)
!!       cam_in(c)%aldif(:ncol)  = natgrid_sfc_flux_in(8)%native_grid_flds(:ncol, c)

       if (masterproc) then
          write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf),   minval(shf): ', &
               maxval(cam_in(c)%shf),   minval(cam_in(c)%shf)
          write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf2),   minval(shf2): ', &
               maxval(shf(:ncol)),   &
               minval(shf(:ncol))
          write(iulog,*) 'PRESC_SFC_FLUX   maxval(lhf),   minval(lhf): ', &
               maxval(cam_in(c)%lhf),   minval(cam_in(c)%lhf)
          write(iulog,*) 'PRESC_SFC_FLUX   maxval(lwup),  minval(lwup): ', &
               maxval(cam_in(c)%lwup),  minval(cam_in(c)%lwup)
          write(iulog,*) 'PRESC_SFC_FLUX   maxval(asdir), minval(asdir): ', &
               maxval(cam_in(c)%asdir), minval(cam_in(c)%asdir)
       endif

       !++BEH a test for zeros, write out useful stuff
       if ( all(cam_in(c)%shf .eq. 0.0_r8) ) then
          call endrun('cam_in SHF is all zeros.  Oh no!'//errmsg(__FILE__,__LINE__))
       endif
       !--BEH

    end do

!    if (masterproc) then
!!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf),   minval(shf): ', &
!!                       maxval(cam_in(c)%shf),   minval(cam_in(c)%shf)
!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf),   minval(shf): ', &
!                       maxval(cam_in(c)%shf),   minval(cam_in(c)%shf)
!!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf2),   minval(shf2): ', &
!!                       maxval(shf(:ncol)),   &
!!                       minval(shf(:ncol))
!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(shf2),   minval(shf2): ', &
!                       maxval(shf(:)),   &
!                       minval(shf(:))
!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(lhf),   minval(lhf): ', &
!                       maxval(cam_in(c)%lhf),   minval(cam_in(c)%lhf)
!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(lwup),  minval(lwup): ', &
!                       maxval(cam_in(c)%lwup),  minval(cam_in(c)%lwup)
!       write(iulog,*) 'PRESC_SFC_FLUX   maxval(asdir), minval(asdir): ', &
!                       maxval(cam_in(c)%asdir), minval(cam_in(c)%asdir)
!    endif



    

end subroutine presc_sfc_flux_overwrite


end module prescribed_sfc_flux

