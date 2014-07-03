module ocn_comp
  !-----------------------------------------------------------------------
  !
  ! Method:  CAM Data Ocean Model
  !
  !-----------------------------------------------------------------------
  use shr_kind_mod,         only: r8 => shr_kind_r8
  use physconst,            only: tmelt
  use shr_sys_mod,          only: shr_sys_abort

  use cam_logfile,          only: iulog
  use cam_control_mod,      only: nsrest, aqua_planet
  use cam_history,          only: outfld
  use filenames,            only: caseid
  use ppgrid,               only: pcols, begchunk, endchunk
  use phys_grid,            only: read_chunk_from_field, write_field_from_chunk, &
                                  gather_chunk_to_field, get_ncols_p, scatter_field_to_chunk
  use units,                only: getunit, freeunit
  use ioFileMod,            only: opnfil, getfil


  use ocn_time_manager,     only: get_step_size, get_nstep, get_curr_calday, &
                                  is_end_curr_day, is_first_step, &
                                  timemgr_write_restart, timemgr_read_restart, &
                                  timemgr_restart, timemgr_init
  use ocn_types,            only: ocn_out_t, frac_t, ocn_types_alloc
  use ocn_spmd
  use sst_data,             only: sstini, sstint, sst, sstcyc
  use perf_mod
  use pio
  use cam_pio_utils,        only : cam_pio_openfile
  !
  implicit none
  private                  ! By default everything private to this module
#include <mpif.h>
  save
  !
  ! Public methods
  !
  public ocn_init          ! Initialization method
  public ocn_run           ! Run method
  public ocn_final         ! Finalization method
  public ocn_write_restart
  public ocn_read_restart
  !
  ! Private module data
  !
  integer :: nrf  = -1                    ! logical unit number for ocn restart dataset
  integer :: nrpf = -1                    ! logical unit number for ocn restart pointer file
  integer, parameter  :: nlen = 256       ! Length of character strings
  character(len=nlen) :: dom_branch_file = ' ' ! full pathname of restart file to branch from (nsrest=3)
  character(len=nlen) :: rest_pfile = './rpointer.ocn' ! restart pointer file contains name of most recently
  character(len=nlen) :: bndtvs           ! sst file
  character(len=nlen) :: focndomain       ! ocn domain file

  type(frac_t), allocatable, public :: frac(:) 
  type(file_desc_t) :: ncid_sst                     ! netcdf file handle for sst file


!======================================================================= 
contains
!======================================================================= 

  subroutine ocn_init( mpicom_ocn, ocn_out, &
       start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
       perpetual_run, perpetual_ymd, calendar)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Initialization method for CAM data ocean model.
    ! 
    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    integer, intent(in)      :: mpicom_ocn	
    type(ocn_out_t), pointer :: ocn_out(:)
    integer, intent(in)      :: start_ymd       ! Start date (YYYYMMDD)
    integer, intent(in)      :: start_tod       ! Start time of day (sec)
    integer, intent(in)      :: ref_ymd         ! Reference date (YYYYMMDD)
    integer, intent(in)      :: ref_tod         ! Reference time of day (sec)
    integer, intent(in)      :: stop_ymd        ! Stop date (YYYYMMDD)
    integer, intent(in)      :: stop_tod        ! Stop time of day (sec)
    logical, intent(in)      :: perpetual_run   ! If in perpetual mode or not
    integer, intent(in)      :: perpetual_ymd   ! Perpetual date (YYYYMMDD)
    character(len=256), intent(in) :: calendar  ! Calendar type
    !
    ! Local variables
    !
    integer  :: ncol                ! number of columns
    integer  :: c                   ! Chunk loop index
    integer  :: i                   ! Column loop index
    integer  :: ierr                ! error status 
    integer  :: unitn               ! namelist unit
    character(len=256) :: locfn     ! netcdf local filename to open
    !
    namelist /dom_inparm/ sstcyc, dom_branch_file, rest_pfile, bndtvs, focndomain
    !--------------------------------------------------------------------------------
    !
    ! Allocate dynamic memory
    ! 
    call ocn_types_alloc( ocn_out )
    call ocn_alloc( )
    !
    ! Initialize ocn MPI communicator 
    !
    call ocn_spmd_init( mpicom_ocn )
    !
    ! Read ocn namelist
    !
    sstcyc = .true.
    if (masterproc) then
       unitn = getunit()
       write(iulog,*) 'Read in camdom namelist from file= ocn_in'
       open( unitn, file='ocn_in', status='old' )
       ierr = 1
       do while ( ierr /= 0 )
          read(unitn, dom_inparm, iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( 'ocn_comp encountered end-of-file on namelist read' )
          endif
       end do
       call freeunit( unitn )
       if (sstcyc) then
          write(iulog,*)'SST dataset will be reused for each model year'
       else
          write(iulog,*)'SST dataset will not be cycled'
       end if
    end if
    call mpi_bcast(sstcyc, 1, MPI_LOGICAL,0,mpicom_ocn,ierr)
    call mpi_bcast(focndomain, nlen, MPI_CHARACTER, 0, mpicom, ierr)
    call mpi_bcast(dom_branch_file, nlen, MPI_CHARACTER, 0, mpicom, ierr)
    call mpi_bcast(rest_pfile, nlen, MPI_CHARACTER, 0, mpicom, ierr)    
    call mpi_bcast(bndtvs, nlen, MPI_CHARACTER, 0, mpicom, ierr)
   !
   ! Data ocean model
   !
   if (nsrest == 0) then
      call ocn_read_inidat( )
   else
      call ocn_read_restart( ocn_out, stop_ymd, stop_tod )
   end if
   !
   ! Obtain time-variant sst datatset
   !
   if (.not. aqua_planet) then
      call getfil(bndtvs, locfn)
      call cam_pio_openfile(ncid_sst, locfn, 0)
   endif
   !
   ! Initialize ocean surface datasets
   !
   call sstini(ncid_sst)
   if (is_first_step()) then
      call sstint(ncid_sst, prev_timestep=.true.)
   else
      call sstint(ncid_sst, prev_timestep=.false.)
   end if
   !
   ! Determine initial surface temperature
   !
   if (is_first_step()) then
      do c = begchunk,endchunk
         ncol = get_ncols_p(c)
         do i = 1, ncol
            ocn_out(c)%ts(i) = sst(i,c) + tmelt
         end do
      end do
   end if
   
 end subroutine ocn_init
 
 !
 !----------------------------------------------------------------------- 
 !
 
  subroutine ocn_run( ocn_out )

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Run method for CAM ocean surface fluxes.
    ! 
    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(ocn_out_t), intent(inout) :: ocn_out(begchunk:endchunk)
    !
    ! Local variables
    !
    integer :: ncol           ! number of columns in chunk
    integer :: i              ! Column index
    integer :: c              ! chunk index
    !-----------------------------------------------------------------------
    !
    ! Output initial data
    !
    call t_startf ('ocn_IC_output')
    call ocn_IC_output( ) 
    call t_stopf ('ocn_IC_output')
    !
    ! Get SST from dataset
    !
    call t_startf ('sstint')
    call sstint (ncid_sst)
    call t_stopf ('sstint')
    !
    ! convert ts units from degC to degK over open ocean
    !
    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       do i = 1,ncol
          ocn_out(c)%ts(i)   = sst(i,c) + tmelt
       end do
    end do

  end subroutine ocn_run

  !
  !-----------------------------------------------------------------------
  !

  subroutine ocn_final( ocn_out )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Finalization of cam data ocean
    !
    !-----------------------------------------------------------------------
    type(ocn_out_t), pointer :: ocn_out(:)

    deallocate (ocn_out)

  end subroutine ocn_final

  !
  !----------------------------------------------------------------------- 
  !

  subroutine ocn_alloc()
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! CAM Data Ocean module data allocatation.
    !
    !-----------------------------------------------------------------------
    integer :: c

    allocate (frac(begchunk:endchunk))
    do c = begchunk,endchunk
       frac(c)%land(:) = 0._r8
    end do

  end subroutine ocn_alloc

  !
  !----------------------------------------------------------------------- 
  !

  subroutine ocn_read_inidat()
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! CAM Data model read of initial data.
    !
    !-----------------------------------------------------------------------
    use scamMod,         only: scmlon,scmlat,single_column
    use shr_scam_mod,    only: shr_scam_GetCloseLatLon
    use ncdio_atm, only : infld
    !
    ! Local variables
    !
    logical :: readvar              ! inquiry:  true => variable exists on file
    type(file_desc_t) :: ncid_landfrac        ! netcdf file id	
    integer :: c, ncols             ! column indices
    type(file_desc_t) :: ncid_dom             ! netcdf file id
    real(r8), pointer :: arr2d(:,:) ! temporary 2D array
    character(len=16) :: fieldname  ! field name
    character(len=256) :: locfn     ! netcdf local filename to open
    integer :: closelatidx
    integer :: closelonidx
    integer :: fracid,rcode
    real(r8) :: closelat,closelon

    !-----------------------------------------------------------------------

    allocate ( arr2d(1:pcols,begchunk:endchunk) )

    if (aqua_planet) then
       arr2d(:,:) = 1._r8
    else
       call getfil(focndomain, locfn)
       call cam_pio_openfile(ncid_dom, locfn, 0)
       if (single_column) then
          call shr_scam_getCloseLatLon(ncid_dom%fh,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          rcode = pio_inq_varid(ncid_dom, 'frac', fracid)
          rcode = pio_get_var(ncid_dom,fracid,(/closelonidx,closelatidx/),(/1,1/),arr2d)
       else
          fieldname = 'frac'
          call infld(fieldname, ncid_dom, 'ni', 'nj', 1, pcols, begchunk,&
               endchunk, arr2d ,readvar, grid_map='PHYS')
!          call read_domain (fieldname, ncid_dom, 'ni', 'nj', 'xc','yc',&
!               1, pcols, begchunk, endchunk, arr2d, readvar)
          if(.not. readvar) call shr_sys_abort('dom: error in reading LANDFRAC')
       end if
       call pio_closefile(ncid_dom)
    end if
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       ! first convert from ocn fraction to land fraction
       frac(c)%land(:ncols) = 1._r8 - arr2d(:ncols,c)
    end do

    deallocate ( arr2d )

  end subroutine ocn_read_inidat

  !
  !-----------------------------------------------------------------------
  !

  subroutine ocn_write_restart( fname, ocn_out )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Write the CAM ocean surface restart file out.
    !
    !-----------------------------------------------------------------------
    !
    ! Input arguments
    !
    character(len=nlen), intent(in) :: fname            ! restart filename
    type(ocn_out_t)    , intent(in) :: ocn_out(begchunk:endchunk)
    !
    ! Local workspace
    !
    real(r8) :: tmpfield(pcols,begchunk:endchunk) ! temporary
    integer  :: i                                 ! loop index
    !-----------------------------------------------------------------------

    if (masterproc) then
       if ( nrf == -1 ) nrf = getunit()
       call opnfil(fname, nrf, 'u')
    endif

    do i=begchunk,endchunk
       tmpfield(:,i) = frac(i)%land(:)
    end do
    call write_field_from_chunk(nrf,1,1,1,tmpfield)

    ! write time manager restart
    
    if (masterproc) then
	call timemgr_write_restart(nrf)
     end if
  
    if (masterproc) then
       close (nrf)

       if ( nrpf == -1 ) nrpf = getunit()
       call opnfil(rest_pfile, nrpf, 'f')
       rewind nrpf
       write (nrpf,'(a)') trim(fname)
       close(nrpf)
       write(iulog,*)'(ocn_write_restart): successfully wrote local restart pointer file ',trim(rest_pfile)
    end if

  end subroutine ocn_write_restart

  !
  !-----------------------------------------------------------------------
  !

  subroutine ocn_read_restart(ocn_out, stop_ymd, stop_tod)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Read the CAM ocean surface restart file in.
    !
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ocn_out_t), intent(inout) :: ocn_out(begchunk:endchunk)
    integer        , intent(in)    :: stop_ymd      ! Stop date (YYYYMMDD)
    integer        , intent(in)    :: stop_tod      ! Stop time of day (sec)
    !
    ! Local workspace
    !
    real(r8) :: tmpfield(pcols,begchunk:endchunk)
    integer  :: i                          ! loop index
    integer  :: ncol                       ! number of vertical columns
    character(len=nlen) :: fname           ! restart filename
    character(len=nlen) :: pname           ! restart full pathname
    !-----------------------------------------------------------------------

    ! restart run =>nsrest=1) and branch run=>nsrest=3.  
    ! Only read the restart pointer file for a restart run.
    
    if (masterproc) then
       if (nsrest == 1) then
          nrpf = getunit()
          call opnfil (rest_pfile, nrpf, 'f', status="old")
          read (nrpf,'(a)') pname
          close(nrpf)
       else
          if ( trim(dom_branch_file) == '' )then
             call shr_sys_abort ('ocn_read_restart: dom_branch_file is empty')
          end if
          if ( dom_branch_file(1:1) /= '/' )then
             call shr_sys_abort ('ocn_read_restart: dom_branch_file is not an absolute pathname')
          end if
          if ( len_trim(dom_branch_file) > nlen )then
             call shr_sys_abort ('ocn_read_restart: dom_branch_file is too long :'//dom_branch_file)
          end if
          pname = trim(dom_branch_file)
       end if
       call getfil(pname, fname)
       nrf = getunit()
       call opnfil(fname, nrf, 'u')
    endif

    call read_chunk_from_field(nrf,1,1,1,tmpfield)
    do i = begchunk,endchunk
       ncol = get_ncols_p(i)
       frac(i)%land(:ncol) = tmpfield(:ncol,i)
    end do

    ! Restart the time manager.
    
    if (masterproc) then
       call timemgr_read_restart(nrf)
    end if
    call timemgr_restart( stop_ymd=stop_ymd, stop_tod=stop_tod )

    if (masterproc) then
       close(nrf)
    end if

  end subroutine ocn_read_restart

!
!----------------------------------------------------------------------- 
!

  subroutine read_domain(varname, ncid , dimlonnam, dimlatnam, lonnam, latnam, &
       dim1b, dim1e, dim2b, dim2e, field, readvar)
    use abortutils, only : endrun
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in)  :: varname     ! variable name
    type(file_desc_t)         , intent(inout)  :: ncid        ! input unit
    character(len=*), intent(in)  :: dimlonnam   ! name of longitude dimension of field on file
    character(len=*), intent(in)  :: dimlatnam   ! name of latitude  dimension of field on file
    character(len=*), intent(in)  :: lonnam      ! name of longitude variable  of field on file
    character(len=*), intent(in)  :: latnam      ! name of latitude  variable  of field on file
    integer         , intent(in)  :: dim1b       ! start of first  dimension of array to be returned
    integer         , intent(in)  :: dim1e       ! end   of first  dimension of array to be returned
    integer         , intent(in)  :: dim2b       ! start of second dimension of array to be returned
    integer         , intent(in)  :: dim2e       ! end   of second dimension of array to be returned
    real(r8)        , intent(out) :: field(dim1b:dim1e,dim2b:dim2e) ! array to be returned (decomposed or global)
    logical         , intent(out) :: readvar     ! true => variable is on initial dataset
    !
    ! local variables
    !
    integer :: i,j                      ! index
    integer :: ierr                      ! error status
    integer :: varid                    ! variable id
    integer :: dimlon, dimlat           ! lon, lat, lev dimension lengths
    integer :: tmptype
    integer :: ndims                    ! number of dimensions
    integer :: dims(PIO_MAX_VAR_DIMS)    ! variable shape
    integer :: londimid, latdimid       ! Dimension ID's
    integer :: strt(3)                  ! start lon, lat, time indices for netcdf 2-d
    integer :: cnt (3)                  ! lon, lat, time counts for netcdf 2-d
    data strt/3*1/                      ! 
    data cnt /1,1,1/                    ! 2-d arrs
    real(r8), pointer :: tmp(:,:)       ! input data
    logical :: readvar_tmp              ! if true, variable is on tape
    character(len=PIO_MAX_NAME) tmpname
    character(len=32) :: subname='read_domain' ! subroutine name
    
    !-----------------------------------------------------------------------
    !
!    call check_var(ncid, varname, varid, readvar_tmp)
!    if (readvar_tmp) then
!       ierr = pio_inq_dimid  (ncid, dimlonnam, londimid)
!       ierr = pio_inq_dimlen (ncid, londimid , dimlon)
!       ierr = pio_inq_dimid  (ncid, dimlatnam, latdimid)
!       ierr = pio_inq_dimlen (ncid, latdimid , dimlat)

       ! Check order of dimensions in variable
!       ierr = pio_inq_varndims(ncid, varid,ndims)
!       ierr = pio_inq_vardimid (ncid, varid, dims(1:ndims))
!       if (dims(1) /= londimid .or. dims(2) /= latdimid .or. ndims > 3) then
!          write(iulog,*) trim(subname), ' Error: Bad number of dims or ordering while reading field ', trim(varname)
!          call endrun()
!       end if

       ! Allocate memory and read variable
!       cnt(1) = dimlon
!       cnt(2) = dimlat
!       allocate ( tmp(dimlon,dimlat) )
!       ierr = pio_get_var (ncid, varid, strt, cnt, tmp)
  
       
     

!    end if  ! end of readvar_tmp

  end subroutine read_domain
  !
  !----------------------------------------------------------------------- 
  !
  subroutine check_var(ncid, varname, varid, readvar)
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(file_desc_t), intent(inout)          :: ncid
    character(len=*), intent(in) :: varname
    integer, intent(out)         :: varid
    logical, intent(out)         :: readvar 
    !
    ! Local Variables
    !
    integer :: ret     ! return value
    !-----------------------------------------------------------------------
    readvar = .true.
    call pio_seterrorhandling(ncid, pio_bcast_error)
    ret = pio_inq_varid (ncid, varname, varid)
    call pio_seterrorhandling(ncid, pio_internal_error)
    if (ret/=PIO_NOERR) then
       if (masterproc) then
          write(iulog,*)'CHECK_VAR Warning:  variable ',trim(varname),' is not on initial dataset'
       end if
       readvar = .false.
    end if
  end subroutine check_var
  !
  !----------------------------------------------------------------------- 
  !
  subroutine ocn_IC_output( )

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Output Ice initial condition data to history files.
    !
    !----------------------------------------------------------------------- 
    integer :: c       ! Chunk index
    !
    ! IC Outfield calls
    !
    do c = begchunk, endchunk
       call outfld('TSOCN&IC   ', sst(1,c), pcols, c) 
    end do
    
  end subroutine ocn_IC_output
  
end module ocn_comp
