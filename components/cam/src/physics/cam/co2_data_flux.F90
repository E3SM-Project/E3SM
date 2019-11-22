!BALLI: take care of the units of CO2!!, do cleanup

module co2_data_flux

!------------------------------------------------------------------------------------------------
! for data reading and interpolation                                           
!------------------------------------------------------------------------------------------------

  use shr_kind_mod,     only : r8 => shr_kind_r8,cx => shr_kind_cx, cl => shr_kind_cl
  use ppgrid,           only : begchunk, endchunk, pcols
  use cam_abortutils,   only : endrun
  use dycore,           only: dycore_is
  use shr_log_mod ,     only: errMsg => shr_log_errMsg
  use input_data_utils, only: time_coordinate
  use cam_pio_utils,    only: cam_pio_openfile
  use pio,              only: file_desc_t, pio_nowrite, pio_closefile, pio_inq_dimid, pio_bcast_error, &
       pio_seterrorhandling, pio_noerr, pio_inquire_dimension
  use time_manager,   only: get_nstep

#ifdef CO2_BILIN_REGRID
  use tracer_data,      only : trfld, trfile, trcdata_init, advance_trcdata
#endif

  use spmd_utils,       only: masterproc

  implicit none

! public type

  public co2_data_flux_type

! public interface
  public co2_data_flux_init
  public co2_data_flux_advance

! private data

  private
 
!--------------------------------------------------------------------------------------------------
type :: co2_data_flux_type          
   
   !To store two time samples from a file to do time interpolation in the next step
   !(pcols,begchunk:endchunk,2)
   real(r8), pointer, dimension(:,:,:) :: co2bdy

   !To store data after time interpolation from two time samples
   !(pcols,begchunk:endchunk)
   real(r8), pointer, dimension(:,:)   :: co2flx

   !Forcing file name
   character(len=cl) :: filename   

   !Data structure to keep track of time
   type(time_coordinate) :: time_coord

   !specie name
   character(len=cl)     :: spec_name
   
   !logical to control first data read
   logical               :: initialized

#ifdef CO2_BILIN_REGRID
   type(trfld),  pointer, dimension(:)   :: fields
   type(trfile)                          :: file
#endif
 
end type co2_data_flux_type

! dimension names for physics grid (physgrid)
logical           :: dimnames_set = .false.
character(len=8)  :: dim1name, dim2name

!===============================================================================
contains
!===============================================================================

subroutine co2_data_flux_init (input_file, spec_name, xin)

!-------------------------------------------------------------------------------
! Initialize co2_data_flux_type instance
!   including initial read of input and interpolation to the current timestep
!-------------------------------------------------------------------------------

   use ioFileMod,        only: getfil
   use cam_grid_support, only: cam_grid_id, cam_grid_check
   use cam_grid_support, only: cam_grid_get_dim_names
   use dyn_grid,         only: get_horiz_grid_dim_d

   ! Arguments
   character(len=*),          intent(in)    :: input_file
   character(len=*),          intent(in)    :: spec_name
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   character(len = cx) :: msg
   integer  :: grid_id, ierr, dim1len, dim2len, dim1id, dim2id ! netcdf file ids and sizes
   integer  :: hdim1_d, hdim2_d    ! model grid size
   real(r8) :: dtime
   type(file_desc_t) :: fh_co2_data_flux
   !----------------------------------------------------------------------------

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

   !find if the "input_file" exists locally and update xin%filename with the input file path
   call getfil(input_file, xin%filename)

   !Do some sanity checks before proceeding futher
   call cam_pio_openfile(fh_co2_data_flux, trim(xin%filename), pio_nowrite)

   !Ask PIO to return the control if it experiences an error so that we can handle it explicitly here
   call pio_seterrorhandling(fh_co2_data_flux, pio_bcast_error)

   !if input file is on a different grid than the model grid
   !(e.g. model is running on an FV grid and input netcdf file is on an SE grid), exit with an error
   if(pio_inq_dimid(fh_co2_data_flux, trim(adjustl(dim1name)), dim1id) /= pio_noerr) then
      !pio_inq_dimid function tries to find dim1name in file with id "fh_co2_data_flux"
      !if it can't find dim1name, it means there is a mismacth in model and netcdf 
      !file grid
      call endrun('ERROR: grid mismatch, failed to find '//dim1name//' dimension in file:'//input_file//' '&
           ' '//errmsg(__FILE__,__LINE__))
   endif
   
   !find if the model and netcdf file has same grid resolution
   call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
   if( dycore_is('SE') )  then
      if(pio_inquire_dimension(fh_co2_data_flux, dim1id, len = dim1len) ==  pio_noerr) then
         if(dim1len /= hdim1_d ) then !compare model grid length with file's
            write(msg,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                 hdim1_d,'), netcdf file is:'//input_file
            call endrun(msg//errmsg(__FILE__,__LINE__))
         endif
      else
         call endrun('ERROR: failed while inquiring dimensions of file:'//input_file//' '//errmsg(__FILE__,__LINE__))
      endif
   elseif( dycore_is('LR')) then
      if(pio_inq_dimid(fh_co2_data_flux, trim(adjustl(dim2name)), dim2id)) then !obtain lat dimension of model
         call endrun('ERROR: failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
              ' '//input_file//' '//errmsg(__FILE__,__LINE__))
      endif
      if(pio_inquire_dimension(fh_co2_data_flux, dim1id, len = dim1len) ==  pio_noerr .and. &
         pio_inquire_dimension(fh_co2_data_flux, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
         if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
            write(msg,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                 hdim1_d,' x ',hdim2_d,'), netcdf file is:'//input_file
            call endrun(msg//errmsg(__FILE__,__LINE__))
         endif
      else
          call endrun('ERROR: failed while inquiring dimensions of file:'//input_file//' '//errmsg(__FILE__,__LINE__))
      endif
   else
      call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
   endif

   !Sanity checks end

   !close file
   call pio_closefile(fh_co2_data_flux)


   !Populate xin data structure
   xin%spec_name = spec_name
   xin%initialized = .false.

   !BALLI: Not sure why dtime is defined this way...try a run without sending this optional arg to see if the results are different
   !initialze time coord, this is used....BALLI??
   dtime = 1.0_r8 - 200.0_r8 / 86400.0_r8
   call xin%time_coord%initialize(input_file, force_time_interp=.true., delta_days=dtime)

   !xin%co2bdy will store values of co2 for two time levels 
   !xin%co2flx is the co2 interpolated in time based on model time
   allocate( xin%co2bdy(pcols,begchunk:endchunk,2), &
             xin%co2flx(pcols,begchunk:endchunk)    )

   !Read the file and populate xin%co2flx once
   call co2_data_flux_advance(xin)

   xin%initialized = .true.

end subroutine co2_data_flux_init

!============================================================================================================

subroutine co2_data_flux_advance (xin)

!-------------------------------------------------------------------------------
! Advance the contents of a co2_data_flux_type instance
!   including reading new data, if necessary
!-------------------------------------------------------------------------------

   use ncdio_atm,        only: infld

   ! Arguments
   type(co2_data_flux_type),  intent(inout) :: xin

   ! Local variables
   logical           :: read_data
   integer           :: indx2_pre_adv
   type(file_desc_t) :: fh_co2_data_flux
   logical           :: found
   integer :: nstep, i, j, k

   real(r8), dimension(pcols,begchunk:endchunk) :: co2bdy_dum
   !----------------------------------------------------------------------------

   read_data = xin%time_coord%read_more() .or. .not. xin%initialized

   indx2_pre_adv = xin%time_coord%indxs(2)
   !if(masterproc)write(102,*)'time_adv_surf'
   call xin%time_coord%advance()

   if ( read_data ) then

      call cam_pio_openfile(fh_co2_data_flux, trim(xin%filename), PIO_NOWRITE)

      ! read time-level 1
      ! skip the read if the needed vals are present in time-level 2
      if (xin%initialized .and. xin%time_coord%indxs(1) == indx2_pre_adv) then
         xin%co2bdy(:,:,1) = xin%co2bdy(:,:,2)
      else
         !NOTE: infld call doesn't do any interpolation in space, it just reads in the data
         call infld(trim(xin%spec_name), fh_co2_data_flux, dim1name, dim2name, &
              1, pcols, begchunk, endchunk, xin%co2bdy(:,:,1), found, &
              gridname='physgrid', timelevel=xin%time_coord%indxs(1))

         if (.not. found) then
            call endrun('ERROR: ' // trim(xin%spec_name) // ' not found'//errmsg(__FILE__,__LINE__))
         endif
      endif

      ! read time-level 2
      call infld(trim(xin%spec_name), fh_co2_data_flux, dim1name, dim2name, &
           1, pcols, begchunk, endchunk, xin%co2bdy(:,:,2), found, &
           gridname='physgrid', timelevel=xin%time_coord%indxs(2))

      do j = -10,12222
         call infld(trim(xin%spec_name), fh_co2_data_flux, dim1name, dim2name, &
              1, pcols, begchunk, endchunk, co2bdy_dum, found, &
              gridname='physgrid', timelevel=j)
         i = j+1
      enddo





      if (.not. found) then
         call endrun('ERROR: ' // trim(xin%spec_name) // ' not found'//errmsg(__FILE__,__LINE__))
      endif

      call pio_closefile(fh_co2_data_flux)
   endif

   ! interpolate between time-levels
   ! If time:bounds is in the dataset, and the dataset calendar is compatible with CAM's,
   ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
   ! generating fluxes that are piecewise constant in time.

   if(masterproc)then
      nstep = get_nstep()
      
      do j = begchunk,endchunk
         do i = 1, pcols
            if(xin%co2bdy(i,j,1) == xin%co2bdy(i,j,2)) then
               write(102,*)'SAME:',nstep,xin%time_coord%indxs(1),xin%time_coord%indxs(2)
            else
               write(102,*)'DIFF:',nstep,xin%time_coord%indxs(1),xin%time_coord%indxs(2)
            endif
            write(102,*)i,j,xin%co2bdy(i,j,1),xin%co2bdy(i,j,2)
         enddo
      enddo
   endif


   if (xin%time_coord%wghts(2) == 0.0_r8) then
      !if(masterproc)write(102,*)'No weights'
      xin%co2flx(:,:) = xin%co2bdy(:,:,1)!real(nstep,r8) * 10e-7_r8 !xin%co2bdy(:,:,1)
   else
      !if(masterproc)write(102,*)'weights are:', xin%time_coord%wghts(2)
      xin%co2flx(:,:) = xin%time_coord%wghts(2)! real(nstep,r8) * 10e-7_r8 !xin%co2bdy(:,:,1) + &
           !xin%time_coord%wghts(2) * (xin%co2bdy(:,:,2) - xin%co2bdy(:,:,1))
   endif
   !BRYCE-unit conversion?????
   ! atm_import_export.F90 wants surface fluxes in kgCO2/m2/s so no conversion necessary

end subroutine co2_data_flux_advance

end module co2_data_flux

