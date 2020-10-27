
module metdata
!----------------------------------------------------------------------- 
!
! BOP
!
! !MODULE: metdata
!
! !DESCRIPTION
! Handles reading and interpolating offline meteorological data which
! is used to drive the dynamics.
!
! !USES
  use shr_kind_mod,       only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use shr_cal_mod,        only: shr_cal_gregorian
  use time_manager,       only: get_curr_date, get_step_size, timemgr_is_caltype
  use spmd_utils,         only: masterproc
  use ppgrid,             only: pcols, pver, begchunk, endchunk
  use time_manager,       only: get_curr_calday, get_curr_date, get_step_size
  use cam_abortutils,         only: endrun
  use dynamics_vars,      only: T_FVDYCORE_GRID

#if ( defined SPMD )
  use mpishorthand,       only: mpicom, mpir8, mpiint,mpichar
  use mod_comm,           only: mp_sendirr,mp_recvirr
#endif
  use perf_mod
  use cam_logfile,        only: iulog
  use pio, only: file_desc_t, pio_put_att, pio_global, pio_get_att, pio_inq_att, pio_inq_dimid, pio_inq_dimlen, &
       pio_closefile, pio_get_var, pio_inq_varid
  use cam_pio_utils,      only: cam_pio_openfile
    

  implicit none

  private  ! all unless made public
  save 

! !PUBLIC MEMBERS

  public :: metdata_dyn_init  ! subroutine to open files, allocate blocked arrays, etc
  public :: metdata_phys_init ! subroutine to allocate chunked arrays
  public :: advance_met       ! subroutine to read more data and interpolate
  public :: get_met_fields    ! interface to set meteorology fields
  public :: get_met_srf1
  public :: get_met_srf2
  public :: get_us_vs       
  public :: metdata_readnl
  public :: met_winds_on_walls
  public :: write_met_restart
  public :: read_met_restart
  public :: met_rlx
  public :: met_fix_mass
  public :: met_srf_feedback

  interface write_met_restart
     Module procedure write_met_restart_bin
     Module procedure write_met_restart_pio
  end interface

  interface read_met_restart
     Module procedure read_met_restart_bin
     Module procedure read_met_restart_pio
  end interface


  !------------------------------------------------------------------
  ! Interface to access the meteorology fields.  Possible invocations
  ! are as follows:
  !   call get_met_fields( physics_state, us, vs , tend, dt )
  !   call get_met_fields( u, v )
  !   call get_met_fields( cam_in_t )
  !------------------------------------------------------------------
  Interface get_met_fields                       ! overload accessors
     Module Procedure get_dyn_flds
     Module Procedure get_uv_centered
     Module Procedure get_ps    
     Module Procedure get_ocn_ice_frcs 
  End Interface
  
  real(r8), allocatable :: met_ps_next(:,:)   ! PS interpolated to next timestep
  real(r8), allocatable :: met_ps_curr(:,:)   ! PS interpolated to next timestep

  logical :: met_cell_wall_winds = .false.  ! true => met data winds are cell centered
  logical :: met_remove_file = .false.  ! delete metdata file when finished with it

  character(len=16) :: met_shflx_name = 'SHFLX'
  character(len=16) :: met_qflx_name = 'QFLX'
  real(r8) :: met_snowh_factor = 1._r8
  real(r8) :: met_shflx_factor = 1._r8
  real(r8) :: met_qflx_factor  = 1._r8
  logical  :: met_srf_feedback = .true.
  logical  :: met_srf_nudge_flux = .true. ! wsx, wsy, shf, and cflx nudged rather than forced.
                                          ! This is done primarily to prevent unrealistic
                                          ! surface temperatures.

! !REVISION HISTORY:
!   31 Oct 2003  Francis Vitt     Creation
!   05 Feb 2004  F Vitt   Removed reading/inperpolating PS for current timestep
!                         -- only met_ps_next is needed
!   10 Nov 2004  F Vitt   Implemented ability to read from series of files
!   16 Dec 2004  F Vitt   Added offline_met_defaultopts and offline_met_setopts
!   14 Jul 2005  W Sawyer Removed pmgrid, spmd_dyn dependencies
!   12 Apr 2006  W Sawyer Removed unneeded ghosting of met_us, met_vs
!   08 Apr 2010  J Edwards Replaced serial netcdf calls with pio interface  
!
! EOP
!----------------------------------------------------------------------- 
! $Id$
! $Author$
!----------------------------------------------------------------------- 

  type input2d
     real(r8), dimension(:,:), pointer :: data => null()
  endtype input2d

  type input3d
     real(r8), dimension(:,:,:), pointer :: data => null()
  endtype input3d

  real(r8), allocatable :: met_t(:,:,:)  ! interpolated temperature 
  real(r8), allocatable :: met_u(:,:,:)  ! interpolated zonal wind
  real(r8), allocatable :: met_v(:,:,:)  ! interpolated meridional wind
  real(r8), allocatable :: met_us(:,:,:) ! interpolated zonal wind -staggered
  real(r8), allocatable :: met_vs(:,:,:) ! interpolated meridional wind -staggered
  real(r8), allocatable :: met_q(:,:,:)  ! interpolated water vapor

  real(r8), allocatable :: met_shflx(:,:)! interpolated surface pressure
  real(r8), allocatable :: met_qflx(:,:) ! interpolated water vapor flux
  real(r8), allocatable :: met_taux(:,:) ! interpolated 
  real(r8), allocatable :: met_tauy(:,:) ! interpolated
  real(r8), allocatable :: met_snowh(:,:) ! interpolated snow height

  real(r8), allocatable :: met_ts(:,:) ! interpolated

  type(input3d) :: met_ti(2)
  type(input3d) :: met_ui(2)
  type(input3d) :: met_vi(2)
  type(input3d) :: met_usi(2)
  type(input3d) :: met_vsi(2)
  type(input3d) :: met_qi(2)

  type(input2d) :: met_psi_next(2)
  type(input2d) :: met_psi_curr(2)
  type(input2d) :: met_shflxi(2)
  type(input2d) :: met_qflxi(2)
  type(input2d) :: met_tauxi(2)
  type(input2d) :: met_tauyi(2)
  type(input2d) :: met_tsi(2)
  type(input2d) :: met_snowhi(2)

  integer :: dateid           ! var id of the date in the netCDF
  integer :: secid            ! var id of the sec data 
  real(r8) :: datatimem = -1.e36_r8     ! time of prv. values read in
  real(r8) :: datatimep = -1.e36_r8     ! time of nxt. values read in
  real(r8) :: datatimemn = -1.e36_r8    ! time of prv. values read in for next timestep
  real(r8) :: datatimepn  = -1.e36_r8   ! time of nxt. values read in for next timestep

  integer, parameter :: nm=1    ! array index for previous (minus) data
  integer, parameter :: np=2    ! array indes for next (plus) data

  real(r8) :: curr_mod_time ! model time - calendar day
  real(r8) :: next_mod_time ! model time - calendar day - next time step

  character(len=256) :: curr_filename, next_filename, met_data_file
  character(len=256) :: met_filenames_list = ''
  character(len=256) :: met_data_path = ''
  type(file_desc_t) :: curr_fileid, next_fileid     ! the id of the NetCDF file
  real(r8), pointer, dimension(:) :: curr_data_times => null()
  real(r8), pointer, dimension(:) :: next_data_times => null()

  real(r8) :: alpha = 1.0_r8 ! don't read in water vapor  
  !   real(r8), private :: alpha = 0.0  ! read in water vaper each time step

  real(r8), parameter ::  D0_0                    =     0.0_r8
  real(r8), parameter ::  D0_5                    =     0.5_r8
  real(r8), parameter ::  D0_75                   =     0.75_r8
  real(r8), parameter ::  D1_0                    =     1.0_r8  
  real(r8), parameter ::  days_per_month          =    30.6_r8
  real(r8), parameter ::  days_per_non_leapyear   =   365.0_r8
  real(r8), parameter ::  days_per_year           =   365.25_r8
  real(r8), parameter ::  seconds_per_day         = 86400.0_r8
  real(r8), parameter ::  fill_value              = -9999.0_r8

  logical :: online_test = .false.
  logical :: debug = .false.

  real(r8) :: met_rlx(pver)
  integer  :: met_levels
  integer  :: num_met_levels
  real(r8) :: met_rlx_top = 60._r8
  real(r8) :: met_rlx_bot = 50._r8
  real(r8) :: met_max_rlx = 1._r8

  logical  :: met_fix_mass = .true.
  logical  :: has_ts = .false.
 
contains

!-------------------------------------------------------------------------
! meteorology data options
!-------------------------------------------------------------------------
  subroutine metdata_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'metdata_readnl'

   namelist /metdata_nl/ &
        met_data_file, &
        met_data_path, & 
        met_remove_file, & 
        met_cell_wall_winds, &
        met_filenames_list, &
        met_rlx_top, &
        met_rlx_bot, &
        met_max_rlx, &
        met_fix_mass, &
        met_shflx_name, &
        met_shflx_factor, &
        met_snowh_factor, &
        met_qflx_name, &
        met_qflx_factor, &
        met_srf_feedback, &
        met_srf_nudge_flux

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'metdata_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, metdata_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#if ( defined SPMD )

   ! Broadcast namelist variables

   call mpibcast (met_data_file  ,len(met_data_file) ,mpichar,0,mpicom)
   call mpibcast (met_data_path  ,len(met_data_path) ,mpichar,0,mpicom)
   call mpibcast (met_remove_file    ,1 ,mpilog, 0, mpicom )
   call mpibcast (met_cell_wall_winds,1 ,mpilog, 0, mpicom )
   call mpibcast (met_filenames_list ,len(met_filenames_list),mpichar,0,mpicom)
   call mpibcast (met_rlx_top,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_rlx_bot,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_max_rlx,        1 ,mpir8,  0, mpicom )
   call mpibcast (met_fix_mass,       1 ,mpilog, 0, mpicom )
   call mpibcast (met_qflx_name      ,len(met_qflx_name),     mpichar,0,mpicom)
   call mpibcast (met_shflx_name     ,len(met_shflx_name),    mpichar,0,mpicom)
   call mpibcast (met_qflx_factor    ,1, mpir8,  0, mpicom )
   call mpibcast (met_shflx_factor   ,1, mpir8,  0, mpicom )
   call mpibcast (met_snowh_factor   ,1, mpir8,  0, mpicom )
   call mpibcast (met_srf_feedback   ,1 ,mpilog, 0, mpicom )
#endif

   if (masterproc) then
       write(iulog,*)'Time-variant meteorological dataset (met_data_file) is: ', trim(met_data_file)
       write(iulog,*)'Meteorological data file will be removed (met_remove_file): ', met_remove_file
       write(iulog,*)'Meteorological winds are on cell walls (met_cell_wall_winds): ', met_cell_wall_winds
       write(iulog,*)'Meteorological file names list file: ', trim(met_filenames_list) 
       write(iulog,*)'Meteorological relaxation top is (km): ', met_rlx_top
       write(iulog,*)'Meteorological relaxation bottom is (km): ', met_rlx_bot
       write(iulog,*)'Meteorological maximum relaxation is : ',met_max_rlx
       write(iulog,*)'Offline driver mass fixer is trurned on (met_fix_mass): ',met_fix_mass
       write(iulog,*)'Meteorological shflx field name : ', trim(met_shflx_name) 
       write(iulog,*)'Meteorological shflx multiplication factor : ', met_shflx_factor
       write(iulog,*)'Meteorological qflx field name : ', trim(met_qflx_name) 
       write(iulog,*)'Meteorological qflx multiplication factor : ', met_qflx_factor 
       write(iulog,*)'Meteorological snowh multiplication factor : ', met_snowh_factor
       write(iulog,*)'Meteorological allow srf models feedbacks : ', met_srf_feedback
    endif

 end subroutine metdata_readnl

!--------------------------------------------------------------------------
! Opens file, allocates arrays
!--------------------------------------------------------------------------
  subroutine metdata_dyn_init(grid)
    use infnan,          only : nan, assignment(=)
    use cam_control_mod, only : nsrest
    implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid


    integer            :: im, km, jfirst, jlast, kfirst, klast 
    integer            :: ng_d, ng_s

    im        = grid%im
    km        = grid%km
    jfirst    = grid%jfirst
    jlast     = grid%jlast
    kfirst    = grid%kfirst
    klast     = grid%klast
    ng_d      = grid%ng_d
    ng_s      = grid%ng_s


    if (nsrest/=1) then ! initial run or branch run
       curr_filename = met_data_file
       next_filename = ''
    else
       ! restart run
       ! curr_filename & next_filename already set by restart_dynamics
    endif

    call open_met_datafile( curr_filename, curr_fileid, curr_data_times, met_data_path, check_dims=.true., grid=grid)

    if ( len_trim(next_filename) > 0 ) &
         call open_met_datafile( next_filename, next_fileid, next_data_times, met_data_path )

!
! allocate space for data arrays ...
! 
    ! dynamics grid

    allocate( met_psi_next(nm)%data(im, jfirst:jlast) )
    allocate( met_psi_next(np)%data(im, jfirst:jlast) )
    allocate( met_psi_curr(nm)%data(im, jfirst:jlast) )
    allocate( met_psi_curr(np)%data(im, jfirst:jlast) )
    allocate( met_ps_next(im, jfirst:jlast) )
    allocate( met_ps_curr(im, jfirst:jlast) )

    allocate( met_us(im, jfirst-ng_d:jlast+ng_s, kfirst:klast) )
    allocate( met_vs(im, jfirst-ng_s:jlast+ng_d, kfirst:klast) )

    met_us = nan
    met_vs = nan

    if (met_cell_wall_winds) then
       allocate( met_usi(nm)%data(im, jfirst:jlast, kfirst:klast) )
       allocate( met_usi(np)%data(im, jfirst:jlast, kfirst:klast) )
       allocate( met_vsi(nm)%data(im, jfirst:jlast, kfirst:klast) )
       allocate( met_vsi(np)%data(im, jfirst:jlast, kfirst:klast) )
    endif

    if (.not. met_cell_wall_winds) then

       allocate( met_u(im, jfirst-ng_d:jlast+ng_d, kfirst:klast) )
       allocate( met_ui(nm)%data(im, jfirst:jlast, kfirst:klast) )
       allocate( met_ui(np)%data(im, jfirst:jlast, kfirst:klast) )

       allocate( met_v(im, jfirst-ng_s:jlast+ng_d, kfirst:klast) )
       allocate( met_vi(nm)%data(im, jfirst:jlast, kfirst:klast) )
       allocate( met_vi(np)%data(im, jfirst:jlast, kfirst:klast) )

    endif

 end subroutine metdata_dyn_init

!=================================================================================

 subroutine metdata_phys_init
   use infnan,      only : nan, assignment(=)
   use cam_history, only : addfld, horiz_only

   call addfld ('MET_RLX',(/ 'lev' /), 'A','     ','Meteorology relax function',gridname='fv_centers')
   call addfld ('MET_TAUX',horiz_only, 'A','        ','Meteorology taux', gridname='physgrid')
   call addfld ('MET_TAUY',horiz_only, 'A','        ','Meteorology tauy', gridname='physgrid')
   call addfld ('MET_SHFX',horiz_only, 'A','        ','Meteorology shflx', gridname='physgrid')
   call addfld ('MET_QFLX',horiz_only, 'A','        ','Meteorology qflx', gridname='physgrid')
   call addfld ('MET_PS',horiz_only, 'A','          ','Meteorology PS',gridname='fv_centers')
   call addfld ('MET_T',(/ 'lev' /), 'A','        ','Meteorology T', gridname='physgrid')
   call addfld ('MET_U',(/ 'lev' /), 'A','        ','Meteorology U',gridname='fv_centers')
   call addfld ('MET_V',(/ 'lev' /), 'A','        ','Meteorology V',gridname='fv_centers')
   call addfld ('MET_SNOWH',horiz_only,    'A','    ','Meteorology snow height', gridname='physgrid')

   call addfld ('MET_TS',horiz_only, 'A','K','Meteorology TS', gridname='physgrid')
   call addfld ('MET_OCNFRC',horiz_only, 'A','fraction','Ocean frac derived from met TS', gridname='physgrid')
   call addfld ('MET_ICEFRC',horiz_only, 'A','fraction','Sea ice frac derived from met TS', gridname='physgrid')

! allocate chunked arrays   

   allocate( met_ti(nm)%data(pcols,pver,begchunk:endchunk) )
   allocate( met_ti(np)%data(pcols,pver,begchunk:endchunk) )
   allocate( met_t(pcols,pver,begchunk:endchunk) )

   allocate( met_qi(nm)%data(pcols,pver,begchunk:endchunk) )
   allocate( met_qi(np)%data(pcols,pver,begchunk:endchunk) )
   allocate( met_q(pcols,pver,begchunk:endchunk) )

   allocate( met_shflxi(nm)%data(pcols,begchunk:endchunk) )
   allocate( met_shflxi(np)%data(pcols,begchunk:endchunk) )
   allocate( met_shflx(pcols,begchunk:endchunk) )

   allocate( met_qflxi(nm)%data(pcols,begchunk:endchunk) )
   allocate( met_qflxi(np)%data(pcols,begchunk:endchunk) )
   allocate( met_qflx(pcols,begchunk:endchunk) )

   allocate( met_tauxi(nm)%data(pcols,begchunk:endchunk) )
   allocate( met_tauxi(np)%data(pcols,begchunk:endchunk) )
   allocate( met_taux(pcols,begchunk:endchunk) )

   allocate( met_tauyi(nm)%data(pcols,begchunk:endchunk) )
   allocate( met_tauyi(np)%data(pcols,begchunk:endchunk) )
   allocate( met_tauy(pcols,begchunk:endchunk) )

   allocate( met_tsi(nm)%data(pcols,begchunk:endchunk) )
   allocate( met_tsi(np)%data(pcols,begchunk:endchunk) )
   allocate( met_ts(pcols,begchunk:endchunk) )
   met_ts(:,:) = nan
   
   if(.not.met_srf_feedback) then
      allocate( met_snowhi(nm)%data(pcols,begchunk:endchunk) )
      allocate( met_snowhi(np)%data(pcols,begchunk:endchunk) )
      allocate( met_snowh(pcols,begchunk:endchunk) )
      met_snowh(:,:) = nan
   endif

   call set_met_rlx()

   ! initialize phys surface fields...
   call get_model_time()
   call check_files()
   call read_phys_srf_flds()
   call interp_phys_srf_flds()
   datatimem = -1.e36_r8
   datatimep = -1.e36_r8
 end subroutine metdata_phys_init


!-----------------------------------------------------------------------
! Reads more data if needed and interpolates data to current model time 
!-----------------------------------------------------------------------
 subroutine advance_met(grid)
   use cam_history, only : outfld
   implicit none

   type (T_FVDYCORE_GRID), intent(in) :: grid

   real(r8) ::  met_rlx_2d(grid%ifirstxy:grid%ilastxy,grid%km)
   integer :: i,j,k, idim

   call t_startf('MET__advance')

!
!
    call get_model_time()

    if ( ( curr_mod_time > datatimep ) .or. &
         ( next_mod_time > datatimepn ) ) then
       call check_files()
    endif

    if ( curr_mod_time > datatimep ) then
       call read_next_metdata(grid)
    end if

    if ( next_mod_time > datatimepn ) then
       call read_next_ps(grid)
    end if

! need to inperpolate the data, regardless !
! each mpi tasks needs to interpolate
    call interpolate_metdata(grid)

    call t_stopf('MET__advance')

    idim = grid%ilastxy - grid%ifirstxy + 1
    do j = grid%jfirstxy, grid%jlastxy
       do k = 1, grid%km
          do i = grid%ifirstxy, grid%ilastxy
             met_rlx_2d(i,k) = met_rlx(k)
          enddo
       enddo
       call outfld('MET_RLX',met_rlx_2d, idim, j)
    enddo
  end subroutine advance_met

!-------------------------------------------------------------------
! Method to get some the meteorology data. 
! Sets the following cam_in_t member fields to the 
! meteorology data :
!   qflx
!   shflx
!   taux
!   tauy
!   snowh
!-------------------------------------------------------------------
  subroutine get_met_srf2( cam_in )
    use camsrfexch,          only: cam_in_t    
    use phys_grid,           only: get_ncols_p
    use cam_history,         only: outfld
    use shr_const_mod,       only: shr_const_stebol

    implicit none

    type(cam_in_t), intent(inout), dimension(begchunk:endchunk) :: cam_in

    integer :: c,ncol,i

    if (met_srf_nudge_flux) then
       do c=begchunk,endchunk
          ncol = get_ncols_p(c)
          cam_in(c)%wsx(:ncol)     = (1._r8-met_rlx(pver)) * cam_in(c)%wsx(:ncol)    + met_rlx(pver) * met_taux(:ncol,c)
          cam_in(c)%wsy(:ncol)     = (1._r8-met_rlx(pver)) * cam_in(c)%wsy(:ncol)    + met_rlx(pver) * met_tauy(:ncol,c)
          cam_in(c)%shf(:ncol)     = (1._r8-met_rlx(pver)) * cam_in(c)%shf(:ncol)    + &
               met_rlx(pver) * (met_shflx(:ncol,c) * met_shflx_factor)
          cam_in(c)%cflx(:ncol,1)  = (1._r8-met_rlx(pver)) * cam_in(c)%cflx(:ncol,1) + &
               met_rlx(pver) * (met_qflx(:ncol,c)  * met_qflx_factor)
       end do                    ! Chunk loop
    else
       do c=begchunk,endchunk
          ncol = get_ncols_p(c)
          cam_in(c)%wsx(:ncol)     = met_taux(:ncol,c)
          cam_in(c)%wsy(:ncol)     = met_tauy(:ncol,c)
          cam_in(c)%shf(:ncol)     = met_shflx(:ncol,c) * met_shflx_factor
          cam_in(c)%cflx(:ncol,1)  = met_qflx(:ncol,c)  * met_qflx_factor
       end do                    ! Chunk loop
    end if

    if (debug) then
       if (masterproc) then
          write(iulog,*)'METDATA   maxval(met_taux),minval(met_taux): ',maxval(met_taux),minval(met_taux)
          write(iulog,*)'METDATA   maxval(met_tauy),minval(met_tauy): ',maxval(met_tauy),minval(met_tauy)
          write(iulog,*)'METDATA maxval(met_shflx),minval(met_shflx): ',maxval(met_shflx),minval(met_shflx)
          write(iulog,*)'METDATA   maxval(met_qflx),minval(met_qflx): ',maxval(met_qflx),minval(met_qflx)
       endif
    endif
    
    do c = begchunk, endchunk
       call outfld('MET_TAUX',cam_in(c)%wsx , pcols   ,c   )
       call outfld('MET_TAUY',cam_in(c)%wsy , pcols   ,c   )
       call outfld('MET_SHFX',cam_in(c)%shf , pcols   ,c   )
       call outfld('MET_QFLX',cam_in(c)%cflx(:,1) , pcols   ,c   )
    enddo

  end subroutine get_met_srf2

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_met_srf1( cam_in )
    use camsrfexch,          only: cam_in_t    
    use phys_grid,           only: get_ncols_p
    use cam_history,         only: outfld
    use shr_const_mod,       only: shr_const_stebol

    implicit none

    type(cam_in_t), intent(inout), dimension(begchunk:endchunk) :: cam_in

    integer :: c,ncol,i
    
    if (met_srf_feedback) return
    if (.not.has_ts) then
       call endrun('The meteorolgy input must have TS to run with met_srf_feedback set to FALSE')
    endif

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       cam_in(c)%ts(:ncol)     = met_ts(:ncol,c)
       do i = 1,ncol
          cam_in(c)%snowhland(i) = met_snowh(i,c)*cam_in(c)%landfrac(i) * met_snowh_factor
       enddo
    end do ! Chunk loop

    if (debug) then
       if (masterproc) then
          write(iulog,*)'METDATA       maxval(met_ts),minval(met_ts): ',maxval(met_ts),minval(met_ts)
          write(iulog,*)'METDATA maxval(met_snowh),minval(met_snowh): ',maxval(met_snowh),minval(met_snowh)
       endif
    endif

    do c = begchunk, endchunk
       call outfld('MET_SNOWH',cam_in(c)%snowhland, pcols   ,c   )
    enddo

  end subroutine get_met_srf1

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_ocn_ice_frcs( lndfrc, ocnfrc, icefrc, lchnk, ncol )

    use shr_const_mod, only: SHR_CONST_TKFRZSW
    use shr_const_mod, only: SHR_CONST_TKFRZ
    use cam_history,   only: outfld

    ! args
    real(r8), intent( in) :: lndfrc (pcols)
    real(r8), intent(out) :: ocnfrc (pcols)
    real(r8), intent(out) :: icefrc (pcols)

    integer, intent(in)   :: lchnk 
    integer, intent(in)   :: ncol

    ! local vars
    integer :: i
    
    if (.not.has_ts) then
       if (masterproc) then
          write(iulog,*) 'get_ocn_ice_frcs: TS is not in the met dataset and cannot set ocnfrc and icefrc'
          write(iulog,*) ' try setting drydep_method to xactive_atm or table'
          call endrun('get_ocn_ice_frcs: TS is not in the met dataset')
       endif
    endif

    do i = 1,ncol

       if ( met_ts(i,lchnk) < SHR_CONST_TKFRZ-2._r8 ) then
          ocnfrc(i) = 0._r8
          icefrc(i) = 1._r8 - lndfrc(i)
       else
          icefrc(i) = 0._r8
          ocnfrc(i) = 1._r8 - lndfrc(i)
       endif

    enddo

    call outfld('MET_TS', met_ts(:ncol,lchnk) , ncol   ,lchnk   )
    call outfld('MET_OCNFRC', ocnfrc(:ncol) , ncol   ,lchnk   )
    call outfld('MET_ICEFRC', icefrc(:ncol) , ncol   ,lchnk   )

  endsubroutine get_ocn_ice_frcs

!-------------------------------------------------------------------
! allows access to physics state fields 
!   q  : water vapor
!   ps : surface pressure
!   t  : temperature
!-------------------------------------------------------------------
  subroutine get_dyn_flds( state, tend, dt )

    use physics_types,  only: physics_state, physics_tend, physics_dme_adjust
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use phys_grid,      only: get_ncols_p
    use cam_history,    only: outfld 

    implicit none

    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: tend
    real(r8),            intent(in   ) :: dt                  ! model physics timestep

    integer :: lats(pcols)           ! array of latitude indices
    integer :: lons(pcols)           ! array of longitude indices
    integer :: c, ncol, i,j,k
    real(r8):: qini(pcols,pver)   ! initial specific humidity

    real(r8) :: tmp(pcols,pver)

    call t_startf('MET__GET_DYN2')
    
    do c = begchunk, endchunk
       ncol = get_ncols_p(c)
       do k=1,pver
          do i=1,ncol
             state(c)%t(i,k) = (1._r8-met_rlx(k))*state(c)%t(i,k) + met_rlx(k)*met_t(i,k,c)

             qini(i,k) = state(c)%q(i,k,1)

             ! at this point tracer mixing ratios have already been
             ! converted from dry to moist
!!$             if (  moist_q_mmr .and. (.not. online_test)) then
                state(c)%q(i,k,1) = alpha*state(c)%q(i,k,1) + &
                     (1-alpha)*met_q(i,k,c)
!!$             else 
!!$                ! dry-to-moist conversion
!!$                state(c)%q(i,k,1) = alpha*state(c)%q(i,k,1) + &
!!$                     (1.-alpha)*met_q(i,c,k) &
!!$                     * state(c)%pdeldry(i,k)/state(c)%pdel(i,k)
!!$             endif

             if ((state(c)%q(i,k,1) < D0_0).and. (alpha .ne. D1_0 )) state(c)%q(i,k,1) = D0_0

          end do

       end do

       ! now adjust mass of each layer now that water vapor has changed
       if (( .not. online_test ) .and. (alpha .ne. D1_0 )) then
          call physics_dme_adjust(state(c), tend(c), qini, dt)
       endif

    end do

    if (debug) then
    if (masterproc) then
      write(iulog,*)'METDATA maxval(met_t),minval(met_t): ', maxval(met_t),minval(met_t)
      write(iulog,*)'METDATA maxval(met_ps_next),minval(met_ps_next): ',  maxval(met_ps_next),minval(met_ps_next)
    endif
    endif
    
    do c = begchunk, endchunk
       call outfld('MET_T  ',state(c)%t , pcols   ,c   )
    enddo
    call t_stopf('MET__GET_DYN2')

  end subroutine get_dyn_flds

!------------------------------------------------------------------------
! get the meteorological winds on the grid cell centers (A-grid)
!------------------------------------------------------------------------
  subroutine get_uv_centered( grid, u, v )

    use cam_history,    only: outfld

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(out) :: u(grid%im, grid%jfirst-grid%ng_d:grid%jlast+grid%ng_d, &
                               grid%kfirst:grid%klast)  ! u-wind on A-grid
    real(r8), intent(out) :: v(grid%im, grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d, &
                               grid%kfirst:grid%klast)  ! v-wind on A-grid

    integer :: i,j,k

    integer :: jm, jfirst, jlast, jfirstxy, jlastxy, kfirst, klast, ng_d, ng_s, ifirstxy, ilastxy

    real(r8) :: u3s_tmp(grid%ifirstxy:grid%ilastxy,grid%km)
    real(r8) :: v3s_tmp(grid%ifirstxy:grid%ilastxy,grid%km)

    jm      = grid%jm
    jfirstxy= grid%jfirstxy
    jlastxy = grid%jlastxy
    jfirst  = grid%jfirst
    jlast   = grid%jlast
    kfirst  = grid%kfirst
    klast   = grid%klast
    ifirstxy= grid%ifirstxy
    ilastxy = grid%ilastxy

    ng_d    = grid%ng_d
    ng_s    = grid%ng_s

    u(:,:,:) = D0_0
    v(:,:,:) = D0_0

    u(         :,  max(1,jfirst-ng_d):min(jm,jlast+ng_d), kfirst:klast ) = &
         met_u(:,  max(1,jfirst-ng_d):min(jm,jlast+ng_d), kfirst:klast )

    v(         :,  max(1,jfirst-ng_s):min(jm,jlast+ng_d), kfirst:klast ) = &
         met_v(:,  max(1,jfirst-ng_s):min(jm,jlast+ng_d), kfirst:klast )

    if (masterproc) then
    if (debug) write(iulog,*)'METDATA maxval(u),minval(u),maxval(v),minval(v) : ',&
         maxval(u(:,  max(1,jfirst-ng_d):min(jm,jlast+ng_d), kfirst:klast )),&
         minval(u(:,  max(1,jfirst-ng_d):min(jm,jlast+ng_d), kfirst:klast )),&
         maxval(v(:,  max(1,jfirst-ng_s):min(jm,jlast+ng_d), kfirst:klast )),&
         minval(v(:,  max(1,jfirst-ng_s):min(jm,jlast+ng_d), kfirst:klast ))
    endif

    if ( grid%twod_decomp .eq. 0 ) then
       do j = jfirst, jlast
          do k = kfirst, klast
             do i = 1, grid%im
                u3s_tmp(i,k) = u(i,j,k)
                v3s_tmp(i,k) = v(i,j,k)
             enddo
          enddo
          call outfld ('MET_U ', u3s_tmp, grid%im, j )
          call outfld ('MET_V ', v3s_tmp, grid%im, j )
       enddo
    endif

  end subroutine get_uv_centered

!------------------------------------------------------------------------
! get the meteorological surface pressure interp to dyn substep
!------------------------------------------------------------------------
  subroutine get_ps( grid, ps, nsubsteps, n )

    use cam_history,    only: outfld

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(out) :: ps(grid%im, grid%jfirst:grid%jlast)
    integer, intent(in) :: nsubsteps
    integer, intent(in) :: n

    real(r8) :: num1, num2
    integer  :: j

    num1 = n
    num2 = nsubsteps

    ps(:,:) = met_ps_curr(:,:) + num1*(met_ps_next(:,:)-met_ps_curr(:,:))/num2

    if ( grid%twod_decomp .eq. 0 ) then
       do j = grid%jfirst, grid%jlast
          call outfld('MET_PS',ps(:,j), grid%im   ,j   )
       enddo
    endif
  end subroutine get_ps

!------------------------------------------------------------------------
! get the meteorological winds on the grid cell walls (vorticity winds)
!   us : staggered zonal wind
!   vs : staggered meridional wind
!------------------------------------------------------------------------
  subroutine get_us_vs( grid, us, vs )

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    real(r8), intent(inout) :: us(grid%im, grid%jfirst-grid%ng_d:grid%jlast+grid%ng_s, &
                                  grid%kfirst:grid%klast)    ! u-wind on d-grid
    real(r8), intent(inout) :: vs(grid%im, grid%jfirst-grid%ng_s:grid%jlast+grid%ng_d, &
                                  grid%kfirst:grid%klast)    ! v-wind on d-grid

    integer :: i,j,k

    integer :: jm, jfirst, jlast, kfirst, klast, ng_d, ng_s

    jm      = grid%jm
    jfirst  = grid%jfirst
    jlast   = grid%jlast
    kfirst  = grid%kfirst
    klast   = grid%klast
    ng_d    = grid%ng_d
    ng_s    = grid%ng_s

    call t_startf('MET__get_us_vs')

    ! vertical relaxation (blending) occurs in dyn_run (dyn_comp.F90)

    us(:,:,:) = 1.e36_r8
    vs(:,:,:) = 1.e36_r8
    us(          :, max(2,jfirst):  min(jm,jlast), kfirst:klast) =      &
         met_us( :, max(2,jfirst):  min(jm,jlast), kfirst:klast)
    vs(          :, max(1,jfirst):  min(jm,jlast), kfirst:klast) =      &
         met_vs( :, max(1,jfirst):  min(jm,jlast), kfirst:klast)
    if (masterproc) then
    if (debug) write(iulog,*)grid%iam,': METDATA maxval(us),minval(us),maxval(vs),minval(vs) : ',&
         maxval(us(          :, max(2,jfirst):  min(jm,jlast),    kfirst:klast)),&
         minval(us(          :, max(2,jfirst):  min(jm,jlast),    kfirst:klast)),&
         maxval(vs(          :, max(1,jfirst):  min(jm,jlast),    kfirst:klast)),&
         minval(vs(          :, max(1,jfirst):  min(jm,jlast),    kfirst:klast))
    endif

!!$    if (debug) then
!!$        u3s_tmp = 1.e36
!!$       do j = jfirst, jlast
!!$          do k = kfirst, klast
!!$             do i = 1, im
!!$                if (j >= 2) u3s_tmp(i,k) = us(i,j,k)
!!$                v3s_tmp(i,k) = vs(i,j,k)
!!$             enddo
!!$          enddo
!!$          call outfld ('MET_US ', u3s_tmp, im, j )
!!$          call outfld ('MET_VS ', v3s_tmp, im, j )
!!$       enddo
!!$    endif
!!$
    call t_stopf('MET__get_us_vs')

  end subroutine get_us_vs

!-------------------------------------------------------------------------
! writes file names to restart file
!-------------------------------------------------------------------------

  subroutine write_met_restart_pio(File)
    type(file_desc_t), intent(inout) :: File
    integer :: ierr    
    ierr =  pio_put_att(File, PIO_GLOBAL, 'current_metdata_filename', curr_filename)
    ierr =  pio_put_att(File, PIO_GLOBAL, 'next_metdata_filename', next_filename)

  end subroutine write_met_restart_pio
  subroutine read_met_restart_pio(File)
    type(file_desc_t), intent(inout) :: File
    
    integer :: ierr, xtype, slen

    ierr = pio_inq_att(File, PIO_GLOBAL, 'current_metdata_filename',xtype, slen)
    ierr = pio_get_att(File, PIO_GLOBAL, 'current_metdata_filename', curr_filename)
    curr_filename(slen+1:256) = ''

    ierr = pio_inq_att(File, PIO_GLOBAL, 'next_metdata_filename',xtype, slen)
    ierr = pio_get_att(File, PIO_GLOBAL, 'next_metdata_filename', next_filename)
    next_filename(slen+1:256) = ''

  end subroutine read_met_restart_pio

  subroutine write_met_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    integer :: ioerr   ! error status

    if (masterproc) then
       write(nrg, iostat=ioerr) curr_filename
       if (ioerr /= 0 ) then
          write(iulog,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('WRITE_RESTART_DYNAMICS')
       end if
       write(nrg, iostat=ioerr) next_filename
       if (ioerr /= 0 ) then
          write(iulog,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('WRITE_RESTART_DYNAMICS')
       end if
    end if
  end subroutine write_met_restart_bin

!-------------------------------------------------------------------------
! reads file names from restart file
!-------------------------------------------------------------------------
  subroutine read_met_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    integer :: ioerr   ! error status

    if (masterproc) then
       read(nrg, iostat=ioerr) curr_filename
       if (ioerr /= 0 ) then
          write(iulog,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('READ_RESTART_DYNAMICS')
       end if
       read(nrg, iostat=ioerr) next_filename
       if (ioerr /= 0 ) then
          write(iulog,*) 'READ ioerror ',ioerr,' on i/o unit = ',nrg
          call endrun ('READ_RESTART_DYNAMICS')
       end if
    end if

#if ( defined SPMD )
    call mpibcast ( curr_filename ,len(curr_filename) ,mpichar,0,mpicom)    
    call mpibcast ( next_filename ,len(next_filename) ,mpichar,0,mpicom)    
#endif
  end subroutine read_met_restart_bin

!-------------------------------------------------------------------------
! returns true if the met winds are defined on cell walls
!-------------------------------------------------------------------------
  function met_winds_on_walls()
    logical :: met_winds_on_walls

    met_winds_on_walls = met_cell_wall_winds
  end function met_winds_on_walls

! internal methods :

!-------------------------------------------------------------------------
! transfers cell-centered winds to cell walls
!-------------------------------------------------------------------------
  subroutine transfer_windsToWalls(grid)

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    integer :: i,j,k
    integer :: im, jfirst, jlast, kfirst, klast

    im      = grid%im
    jfirst  = grid%jfirst
    jlast   = grid%jlast
    kfirst  = grid%kfirst
    klast   = grid%klast

    call t_startf('MET__transfer_windsToWalls')

!$omp parallel do private (i, j, k)
    do k = kfirst, klast

       do j = jfirst+1,jlast
          do i = 1,im
             met_us(i,j,k) = ( met_u(i,j,k) + met_u(i,j-1,k) )*D0_5
          end do
       end do

#if defined( SPMD )
       if ( jfirst .gt. 1 ) then
          do i = 1, im
             ! met_u is alread ghosted at this point
             met_us(i,jfirst,k) = ( met_u(i,jfirst,k) + met_u(i,jfirst-1,k) )*D0_5
          enddo
       endif
#endif

       do j = jfirst,jlast
          met_vs(1,j,k) = ( met_v(1,j,k) + met_v(im,j,k) )*D0_5
          do i = 2,im
             met_vs(i,j,k) = ( met_v(i,j,k) + met_v(i-1,j,k) )*D0_5
          end do
       end do
    end do

    call t_stopf('MET__transfer_windsToWalls')

  end subroutine transfer_windsToWalls

  subroutine get_model_time()
    implicit none
    integer yr, mon, day, ncsec  ! components of a date

    call t_startf('MET__get_model_time')

    call get_curr_date(yr, mon, day, ncsec)

    curr_mod_time = get_time_float( yr, mon, day, ncsec )
    next_mod_time = curr_mod_time + get_step_size()/seconds_per_day

    call t_stopf('MET__get_model_time')

  end subroutine get_model_time

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine check_files()

    use shr_sys_mod, only: shr_sys_system
    use ioFileMod, only: getfil

    implicit none

!-----------------------------------------------------------------------
! 	... local variables
!-----------------------------------------------------------------------
    character(len=256) :: ctmp
    character(len=256) :: loc_fname   
    integer            :: istat


    if (next_mod_time > curr_data_times(size(curr_data_times))) then
       if ( .not. associated(next_data_times) ) then
          ! open next file...
          next_filename = incr_filename( curr_filename )
          call open_met_datafile( next_filename, next_fileid, next_data_times, met_data_path )
       endif
    endif

    if ( associated(next_data_times) ) then
       if (curr_mod_time >= next_data_times(1)) then

          ! close current file ...
          call pio_closefile( curr_fileid )
          if (masterproc) then
             ! remove if requested
             if( met_remove_file ) then
                call getfil( curr_filename, loc_fname, 0 )
                write(iulog,*) 'check_files: removing file = ',trim(loc_fname) 
                ctmp = 'rm -f ' // trim(loc_fname) 
                write(iulog,*) 'check_files: fsystem issuing command - '
                write(iulog,*) trim(ctmp)
                call shr_sys_system( ctmp, istat )
             end if
          endif

          curr_filename = next_filename
          curr_fileid = next_fileid

          deallocate( curr_data_times )
          allocate( curr_data_times( size( next_data_times ) ) )
          curr_data_times(:) = next_data_times(:)

          next_filename = ''

          deallocate( next_data_times )
          nullify(  next_data_times )

       endif
    endif

  end subroutine check_files

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  function incr_filename( filename )

    !-----------------------------------------------------------------------
    ! 	... Increment or decrement a date string withing a filename
    !           the filename date section is assumed to be of the form
    !           yyyy-dd-mm
    !-----------------------------------------------------------------------

    use string_utils, only : incstr
    use shr_file_mod, only : shr_file_getunit, shr_file_freeunit

    implicit none


    character(len=*), intent(in) :: filename ! present dynamical dataset filename
    character(len=256) :: incr_filename      ! next filename in the sequence

    ! set new next_filename ...

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: pos, pos1, istat
    character(len=256) :: fn_new, line
    character(len=6)   :: seconds
    character(len=5)   :: num
    integer :: ios,unitnumber

    if ( len_trim(met_filenames_list) .eq. 0) then
       !-----------------------------------------------------------------------
       !	... ccm type filename
       !-----------------------------------------------------------------------
       pos = len_trim( filename )
       fn_new = filename(:pos)
       write(iulog,*) 'incr_flnm: old filename = ',trim(fn_new)
       if( fn_new(pos-2:) == '.nc' ) then
          pos = pos - 3
       end if
       istat = incstr( fn_new(:pos), 1 )
       if( istat /= 0 ) then
          write(iulog,*) 'incr_flnm: incstr returned ', istat
          write(iulog,*) '           while trying to decrement ',trim( fn_new )
          call endrun
       end if

    else

       ! open met_filenames_list
       write(iulog,*) 'incr_flnm: old filename = ',trim(filename)
       write(iulog,*) 'incr_flnm: open met_filenames_list : ',met_filenames_list 
       unitnumber = shr_file_getUnit()
       open( unit=unitnumber, file=met_filenames_list, iostat=ios, status="OLD")
       if (ios /= 0) then
          call endrun('not able to open met_filenames_list file: '//met_filenames_list)
       endif

       ! read file names
       read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          call endrun('not able to increment file name from met_filenames_list file: '//met_filenames_list)
       endif
       do while( trim(line) /= trim(filename) )
          read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
          if (ios /= 0) then
             call endrun('not able to increment file name from met_filenames_list file: '//met_filenames_list)
          endif
       enddo

       read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          call endrun('not able to increment file name from met_filenames_list file: '//met_filenames_list)
       endif
       fn_new = trim(line)

       close(unit=unitnumber)
       call shr_file_freeUnit(unitnumber)
    endif
    incr_filename = trim(fn_new)
    write(iulog,*) 'incr_flnm: new filename = ',incr_filename

  end function incr_filename

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine find_times( itms, fids, datatm, datatp, time )

    implicit none

    integer, intent(out) :: itms(2) ! record numbers that bracket time
    type(file_desc_t), intent(out) :: fids(2) ! ids of files that contains these recs
    real(r8), intent(in) :: time    ! time of interest
    real(r8), intent(out):: datatm, datatp

    integer np1        ! current forward time index of dataset
    integer n,i      ! 
    integer :: curr_tsize, next_tsize, all_tsize

    real(r8), allocatable, dimension(:):: all_data_times

    curr_tsize = size(curr_data_times)
    next_tsize = 0
    if ( associated(next_data_times)) next_tsize = size(next_data_times)

    all_tsize = curr_tsize + next_tsize

    allocate( all_data_times( all_tsize ) )

    all_data_times(:curr_tsize) = curr_data_times(:)
    if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = next_data_times(:)

    ! find bracketing times 
    do n=1, all_tsize-1
       np1 = n + 1
       datatm = all_data_times(n)
       datatp = all_data_times(np1)
       if ( (time .ge. datatm) .and. (time .le. datatp) ) then
          goto 20
       endif
    enddo

    write(iulog,*)'FIND_TIMES: Failed to find dates bracketing desired time =', time
    write(iulog,*)' datatm = ',datatm
    write(iulog,*)' datatp = ',datatp
    write(iulog,*)' all_data_times = ',all_data_times

    call endrun

20  continue

    deallocate( all_data_times )
  
    itms(1) = n
    itms(2) = np1
    fids(:) = curr_fileid
  
    do i=1,2
       if ( itms(i) > curr_tsize ) then 
          itms(i) = itms(i) - curr_tsize 
          fids(i) = next_fileid
       endif
    enddo

  end subroutine find_times

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine read_next_ps(grid)
    use ncdio_atm,          only: infld
    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid

    integer :: recnos(2)
    type(file_desc_t) :: fids(2)       
    character(len=8) :: varname
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy

    real(r8) :: wrk_xy(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy )

    logical :: readvar

    if(online_test) then
       varname='arch_PS'
    else
       varname='PS'
    end if

    jfirstxy= grid%jfirstxy
    jlastxy = grid%jlastxy
    ifirstxy= grid%ifirstxy
    ilastxy = grid%ilastxy

    call find_times( recnos, fids, datatimemn, datatimepn, next_mod_time )

    call infld(varname, fids(1), 'lon', 'lat',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
         wrk_xy, readvar, gridname='fv_centers', timelevel=recnos(1))

    ! transpose xy -> yz decomposition
    call transpose_xy2yz_2d( wrk_xy, met_psi_next(nm)%data, grid )

    call infld(varname, fids(2), 'lon', 'lat',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
         wrk_xy, readvar, gridname='fv_centers', timelevel=recnos(2))

    ! transpose xy -> yz decomposition
    call transpose_xy2yz_2d( wrk_xy, met_psi_next(np)%data, grid )

    if(masterproc) write(iulog,*)'READ_NEXT_PS: Read meteorological data  '

  end subroutine read_next_ps

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine transpose_xy2yz_2d( xy_2d, yz_2d, grid )

    implicit none
    type (T_FVDYCORE_GRID), intent(in)  :: grid
    real(r8),               intent(in)  :: xy_2d(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy)
    real(r8),               intent(out) :: yz_2d(1:grid%im, grid%jfirst:grid%jlast)

    real(r8) :: xy3(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, 1:grid%npr_z )
    integer :: i,j,k

    if (grid%iam .lt. grid%npes_xy) then
       if ( grid%twod_decomp .eq. 1 ) then

#if defined( SPMD )
!$omp parallel do private(i,j,k)
          do k=1,grid%npr_z
             do j=grid%jfirstxy,grid%jlastxy
                do i=grid%ifirstxy,grid%ilastxy
                   xy3(i,j,k) = xy_2d(i,j)
                enddo
             enddo
          enddo

          call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc, &
               grid%xy2d_to_yz2d%RecvDesc, xy3, yz_2d,             &
               modc=grid%modc_dynrun )
          call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc, &
               grid%xy2d_to_yz2d%RecvDesc, xy3, yz_2d,             &
               modc=grid%modc_dynrun )
#endif

       else
          yz_2d(:,:) = xy_2d(:,:)
       endif
    endif  !  (grid%iam .lt. grid%npes_xy)

  end subroutine transpose_xy2yz_2d

!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine transpose_xy2yz_3d( xy_3d, yz_3d, grid )

    implicit none
    type (T_FVDYCORE_GRID), intent(in)  :: grid
    real(r8),               intent(in)  :: xy_3d(grid%ifirstxy:grid%ilastxy, grid%jfirstxy:grid%jlastxy, 1:grid%km )
    real(r8),               intent(out) :: yz_3d(1:grid%im, grid%jfirst:grid%jlast, grid%kfirst:grid%klast)

    if (grid%iam .lt. grid%npes_xy) then
       if ( grid%twod_decomp .eq. 1 ) then
#if defined( SPMD )
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc, &
               grid%ijk_xy_to_yz%RecvDesc, xy_3d, yz_3d,            &
               modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc, &
               grid%ijk_xy_to_yz%RecvDesc, xy_3d, yz_3d,            &
               modc=grid%modc_dynrun )
#endif
       else
          yz_3d(:,:,:) = xy_3d(:,:,:)
       endif
    endif  !  (grid%iam .lt. grid%npes_xy)

  end subroutine transpose_xy2yz_3d



!------------------------------------------------------------------------
!------------------------------------------------------------------------
  subroutine read_next_metdata(grid)
    use ncdio_atm,          only: infld
    use cam_grid_support,   only: cam_grid_check, cam_grid_id
    use cam_grid_support,   only: cam_grid_get_dim_names

    implicit none

    type (T_FVDYCORE_GRID), intent(in) :: grid
    integer recnos(2),  i      ! 
    type(file_desc_t) :: fids(2)

    character(len=8) :: Uname, Vname, Tname, Qname, psname
    character(len=8) :: dim1name, dim2name
    integer :: im, jm, km
    logical :: readvar
    integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
    real(r8), allocatable :: wrk2_xy(:,:)
    real(r8), allocatable :: wrk3_xy(:,:,:)
    real(r8), allocatable :: tmp_data(:,:,:)
    integer :: elev1,blev1, elev2,blev2
    integer :: elev3,blev3, elev4,blev4
    integer :: grid_id  ! grid ID for data mapping

    call t_startf('MET__read_next_metdata')

    jfirstxy= grid%jfirstxy
    jlastxy = grid%jlastxy
    ifirstxy= grid%ifirstxy
    ilastxy = grid%ilastxy

    im      = grid%im
    jm      = grid%jm
    km      = grid%km

    call find_times( recnos, fids, datatimem, datatimep, curr_mod_time )
    !
    ! Set up hyperslab corners
    !

    if(online_test) then
       Tname='arch_T'
       Qname='arch_Q'
       PSname='arch_PS'
       if(met_cell_wall_winds) then
          Uname='arch_US'
          Vname='arch_VS'
       else
          Uname='arch_U'
          Vname='arch_V'
       end if
    else
       Tname='T'
       Qname='Q'
       PSname='PS'
       if(met_cell_wall_winds) then
          Uname='US'
          Vname='VS'
       else
          Uname='U'
          Vname='V'
       end if
       
    end if


    if ( num_met_levels>km ) then

       blev1 = 1
       elev1 = km

       blev2 = num_met_levels-km+1
       elev2 = num_met_levels

       blev3 = num_met_levels-km+1
       elev3 = num_met_levels

       blev4 = 1
       elev4 = num_met_levels

    else

       blev1 = km-num_met_levels+1
       elev1 = km

       blev2 = 1
       elev2 = num_met_levels

       blev3 = 1
       elev3 = km

       blev4 = km-num_met_levels+1
       elev4 = km

    endif

    allocate(tmp_data(pcols, 1:num_met_levels, begchunk:endchunk))
    allocate(wrk2_xy(ifirstxy:ilastxy, jfirstxy:jlastxy))
    allocate(wrk3_xy(ifirstxy:ilastxy, jfirstxy:jlastxy, 1:max(km,num_met_levels)))

    ! physgrid intput for FV is probably always lon/lat but let's be pedantic
    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
       call endrun('read_next_metdata: Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

    do i=1,2

       met_ti(i)%data = 0._r8

       call infld(Tname, fids(i), dim1name, 'lev', dim2name,  1, pcols, 1,num_met_levels , &
            begchunk, endchunk, tmp_data, readvar, gridname='physgrid',timelevel=recnos(i))

       met_ti(i)%data(:,blev1:elev1,:) = tmp_data(:, blev2:elev2, :)

       met_qi(i)%data = 0._r8

       call infld(Qname, fids(i), dim1name, 'lev', dim2name,  1, pcols, 1,num_met_levels, &
            begchunk, endchunk, tmp_data, readvar, gridname='physgrid',timelevel=recnos(i))

       met_qi(i)%data(:,blev1:elev1,:) = tmp_data(:, blev2:elev2, :)

       if (met_cell_wall_winds) then 

          wrk3_xy = 0._r8
          met_usi(i)%data(:,:,:) = 0._r8
          call infld(Uname, fids(i), 'lon', 'slat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1,num_met_levels, wrk3_xy(:,:,blev4:elev4), readvar, gridname='fv_u_stagger',timelevel=recnos(i))

          ! transpose xy -> yz decomposition
          call transpose_xy2yz_3d( wrk3_xy(:,:,blev3:elev3), met_usi(i)%data(:,:,:), grid )

          wrk3_xy = 0._r8
          met_vsi(i)%data(:,:,:) = 0._r8
          call infld(Vname, fids(i), 'slon', 'lat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1,num_met_levels, wrk3_xy(:,:,blev4:elev4), readvar, gridname='fv_v_stagger',timelevel=recnos(i))

          ! transpose xy -> yz decomposition
          call transpose_xy2yz_3d( wrk3_xy(:,:,blev3:elev3), met_vsi(i)%data(:,:,:), grid )

       else

          ! read into lower portion of the array...

          wrk3_xy = 0._r8
          met_ui(i)%data = 0._r8
          call infld(Uname, fids(i), 'lon', 'lat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1,num_met_levels, wrk3_xy(:,:,blev4:elev4), readvar, gridname='fv_centers',timelevel=recnos(i))

          ! transpose xy -> yz decomposition
          call transpose_xy2yz_3d( wrk3_xy(:,:,blev3:elev3), met_ui(i)%data(:,:,:), grid )

          wrk3_xy = 0._r8
          met_vi(i)%data = 0._r8
          call infld(Vname, fids(i), 'lon', 'lat', 'lev',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
               1,num_met_levels, wrk3_xy(:,:,blev4:elev4), readvar, gridname='fv_centers',timelevel=recnos(i))

          ! transpose xy -> yz decomposition
          call transpose_xy2yz_3d( wrk3_xy(:,:,blev3:elev3), met_vi(i)%data(:,:,:), grid )

       endif ! met_cell_wall_winds

       call infld(PSname, fids(i), 'lon', 'lat',  ifirstxy, ilastxy, jfirstxy, jlastxy, &
            wrk2_xy, readvar, gridname='fv_centers', timelevel=recnos(i))

       ! transpose xy -> yz decomposition
        call transpose_xy2yz_2d( wrk2_xy, met_psi_curr(i)%data, grid )

    enddo

    deallocate(tmp_data)
    deallocate(wrk3_xy)
    deallocate(wrk2_xy)

    ! 2-D feilds
    call read_phys_srf_flds( )

    if(masterproc) write(iulog,*)'READ_NEXT_METDATA: Read meteorological data '

    call t_stopf('MET__read_next_metdata')

  end subroutine read_next_metdata

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine read_phys_srf_flds( )
    use ncdio_atm,          only: infld

    integer :: i, recnos(2)
    type(file_desc_t) :: fids(2)
    logical :: readvar

    call find_times( recnos, fids, datatimem, datatimep, curr_mod_time )
    do i=1,2

       call infld(met_shflx_name, fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
            met_shflxi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       call infld(met_qflx_name, fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
            met_qflxi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       call infld('TAUX', fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
            met_tauxi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       call infld('TAUY', fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
            met_tauyi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       if ( .not.met_srf_feedback ) then
          call infld('SNOWH', fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
               met_snowhi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       endif
       if (has_ts) then
          call infld('TS', fids(i), 'lon', 'lat',  1, pcols, begchunk, endchunk, &
               met_tsi(i)%data, readvar, gridname='physgrid',timelevel=recnos(i))
       endif
    enddo
  end subroutine read_phys_srf_flds

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine interp_phys_srf_flds( )

    real(r4) :: fact1, fact2
    real(r8) :: deltat

    deltat = datatimep - datatimem
    fact1 = (datatimep - curr_mod_time)/deltat
    fact2 = D1_0-fact1

    met_shflx(:,:) = fact1*met_shflxi(nm)%data(:,:) + fact2*met_shflxi(np)%data(:,:)
    met_qflx(:,:) = fact1*met_qflxi(nm)%data(:,:) + fact2*met_qflxi(np)%data(:,:)
    met_taux(:,:) = fact1*met_tauxi(nm)%data(:,:) + fact2*met_tauxi(np)%data(:,:)
    met_tauy(:,:) = fact1*met_tauyi(nm)%data(:,:) + fact2*met_tauyi(np)%data(:,:)
    if ( .not.met_srf_feedback ) then
       met_snowh(:,:) = fact1*met_snowhi(nm)%data(:,:) + fact2*met_snowhi(np)%data(:,:)
    endif
    if (has_ts) then
       met_ts(:,:) = fact1*met_tsi(nm)%data(:,:) + fact2*met_tsi(np)%data(:,:)
    endif

  end subroutine interp_phys_srf_flds
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine interpolate_metdata(grid)

#if defined( SPMD )
    use mod_comm, only : mp_send4d_ns, mp_recv4d_ns
#endif

    implicit none
    type (T_FVDYCORE_GRID), intent(in) :: grid

    real(r4) fact1, fact2
    real(r4) nfact1, nfact2
    real(r8) deltat,deltatn
    integer :: i,j,k
    integer :: im, jm, km, jfirst, jlast, kfirst, klast, ng_d, ng_s

    im      = grid%im
    jm      = grid%jm
    km      = grid%km
    jfirst  = grid%jfirst
    jlast   = grid%jlast
    kfirst  = grid%kfirst
    klast   = grid%klast
    ng_d    = grid%ng_d
    ng_s    = grid%ng_s

    call t_startf('MET__interpolate_metdata')

    deltat = datatimep - datatimem
    deltatn = datatimepn - datatimemn

    fact1 = (datatimep - curr_mod_time)/deltat
!    fact2 = (curr_mod_time - datatimem)/deltat
    fact2 = D1_0-fact1

    nfact1 = (datatimepn - next_mod_time)/deltatn
!    nfact2 = (next_mod_time - datatimemn)/deltatn
    nfact2 = D1_0-nfact1

    met_t(:,:,:) = fact1*met_ti(nm)%data(:,:,:) + fact2*met_ti(np)%data(:,:,:)
    met_q(:,:,:) = fact1*met_qi(nm)%data(:,:,:) + fact2*met_qi(np)%data(:,:,:)

    if (.not. online_test) where (met_q .lt. D0_0) met_q = D0_0

    met_ps_next(:,:) =  nfact1*met_psi_next(nm)%data(:,:) + nfact2*met_psi_next(np)%data(:,:)
    met_ps_curr(:,:) =   fact1*met_psi_curr(nm)%data(:,:) +  fact2*met_psi_curr(np)%data(:,:)

    call interp_phys_srf_flds()

    if (has_ts) then
       met_ts(:,:) = fact1*met_tsi(nm)%data(:,:) + fact2*met_tsi(np)%data(:,:)
    endif

    if ( .not. met_cell_wall_winds ) then

       met_u(1:im,jfirst:jlast,kfirst:klast) = fact1*met_ui(nm)%data(1:im,jfirst:jlast,kfirst:klast) &
                                                 + fact2*met_ui(np)%data(1:im,jfirst:jlast,kfirst:klast)
       met_v(1:im,jfirst:jlast,kfirst:klast) = fact1*met_vi(nm)%data(1:im,jfirst:jlast,kfirst:klast) &
                                                 + fact2*met_vi(np)%data(1:im,jfirst:jlast,kfirst:klast)

       ! ghost u,v
#if defined( SPMD )
       call mp_send4d_ns( grid%commxy, im, jm, km, 1, jfirst, jlast,       &
            kfirst, klast, ng_d, ng_d, met_u )
       call mp_send4d_ns( grid%commxy, im, jm, km, 1, jfirst, jlast,       &
            kfirst, klast, ng_d, ng_s, met_v )
       call mp_recv4d_ns( grid%commxy, im, jm, km, 1, jfirst, jlast,       &
            kfirst, klast, ng_d, ng_d, met_u )
       call mp_recv4d_ns( grid%commxy, im, jm, km, 1, jfirst, jlast,       &
            kfirst, klast, ng_d, ng_s, met_v )
#endif

       ! average to cell walls (vorticity winds)
       call transfer_windsToWalls(grid)
    else
       met_us(:,jfirst:jlast,kfirst:klast) = fact1*met_usi(nm)%data(:,jfirst:jlast,kfirst:klast) + &
            fact2*met_usi(np)%data(:,jfirst:jlast,kfirst:klast)
       met_vs(:,jfirst:jlast,kfirst:klast) = fact1*met_vsi(nm)%data(:,jfirst:jlast,kfirst:klast) + &
            fact2*met_vsi(np)%data(:,jfirst:jlast,kfirst:klast)

    endif

    ! ghost staggered u,v
    ! WS 2006.04.11: not necessary here since it will be done in cd_core

!    write(iulog,*)'INTERPOLATE_METDATA: complete.'

    call t_stopf('MET__interpolate_metdata')

  end subroutine interpolate_metdata

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine get_dimension( fid, dname, dsize )
    implicit none
    type(file_desc_t), intent(in) :: fid
    character(*), intent(in) :: dname
    integer, intent(out) :: dsize

    integer :: dimid, ierr

    ierr = pio_inq_dimid( fid, dname, dimid )
    ierr = pio_inq_dimlen( fid, dimid, dsize )

  end subroutine get_dimension

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine open_met_datafile( fname, fileid, times, datapath, check_dims, grid )

    use ioFileMod, only: getfil
    use pio,       only: pio_seterrorhandling, PIO_INTERNAL_ERROR, PIO_BCAST_ERROR, PIO_NOERR

    implicit none

    character(*), intent(in) :: fname
    type(file_desc_t), intent(inout) :: fileid
    real(r8), pointer, intent(inout) :: times(:)
    character(*), intent(in) :: datapath
    logical, optional, intent(in) :: check_dims
    type (T_FVDYCORE_GRID), optional, intent(in) :: grid

    character(len=256) :: filepath   
    character(len=256) :: filen   
    integer :: year, month, day, dsize, i, timesize
    integer :: dateid,secid
    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: ierr
    integer :: im, jm, km
    integer :: varid

    !
    ! open file and get fileid
    !
    if (len_trim( datapath )>0) then
       filepath = trim(datapath)//'/'//trim(fname)
    else
       filepath = trim(fname)
    endif
    call getfil( filepath, filen, 0 )

    call cam_pio_openfile( fileid, filen, 0 )

    call pio_seterrorhandling(fileid, PIO_BCAST_ERROR)

    ierr = pio_inq_varid( fileid, 'TS', varid )
    has_ts = ierr==PIO_NOERR

    call pio_seterrorhandling(fileid, PIO_INTERNAL_ERROR)

    if (masterproc) write(iulog,*) 'open_met_datafile: ',trim(filen)

    call get_dimension( fileid, 'time', timesize )

    if ( associated(times) ) deallocate(times)
    allocate( times(timesize) )

    allocate( dates(timesize) )
    allocate( datesecs(timesize) )

    ierr = pio_inq_varid( fileid, 'date',    dateid  )
    ierr = pio_inq_varid( fileid, 'datesec', secid  )

    ierr = pio_get_var( fileid, dateid, dates )
    ierr = pio_get_var( fileid, secid,  datesecs  )

    do i=1,timesize
       year = dates(i) / 10000
       month = mod(dates(i),10000)/100
       day = mod(dates(i),100)
       times(i) = get_time_float( year, month, day, datesecs(i) )
    enddo

    deallocate( dates )
    deallocate( datesecs )       


!
! check that the data dim sizes match models dimensions
!
    if (present(check_dims) .and. present(grid)) then
       im         = grid%im
       jm         = grid%jm
       km         = grid%km

       if (check_dims) then

          call get_dimension( fileid, 'lon', dsize )
          if (dsize /= im) then
             write(iulog,*)'open_met_datafile: lonsiz=',dsize,' must = ',im
             call endrun
          endif
          call get_dimension( fileid, 'lat', dsize )
          if (dsize /= jm) then
             write(iulog,*)'open_met_datafile: latsiz=',dsize,' must = ',jm
             call endrun
          endif
          call get_dimension( fileid, 'lev', dsize )
          met_levels = min( dsize, km )
          num_met_levels = dsize
       endif
    endif

  end subroutine open_met_datafile

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  function get_time_float( year, month, day, sec )

! returns float representation of time -- number of days
! since 1 jan 0001 00:00:00.000

    implicit none

    integer, intent(in) :: year, month, day
    integer, intent(in) :: sec
    real(r8) :: get_time_float

! ref date is 1 jan 0001

    integer :: refyr, refmn, refdy
    real(r8) :: refsc, fltdy
    integer :: doy(12)

!              jan feb mar apr may jun jul aug sep oct nov dec
!              31  28  31  30  31  30  31  31  31  31  30  31
    data doy /  1, 32, 60, 91,121,152,182,213,244,274,305,335 /

    refyr = 1
    refmn = 1
    refdy = 1
    refsc = D0_0

    if ( timemgr_is_caltype(trim(shr_cal_gregorian))) then
       fltdy = greg2jday(year, month, day) - greg2jday(refyr,refmn,refdy)
    else ! assume no_leap (all years are 365 days)
       fltdy = (year - refyr)*days_per_non_leapyear + &
               (doy(month)-doy(refmn)) + &
               (day-refdy) 
    endif

    get_time_float = fltdy + ((sec-refsc)/seconds_per_day)

  endfunction get_time_float

!-----------------------------------------------------------------------
! 	... Return Julian day number given Gregorian date.
!
! Algorithm from Hatcher,D.A., Simple Formulae for Julian Day Numbers
! and Calendar Dates, Q.Jl.R.astr.Soc. (1984) v25, pp 53-55.
!-----------------------------------------------------------------------
  function greg2jday( year, month, day )

    implicit none

    integer, intent(in) :: year, month, day
    integer :: greg2jday

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer :: ap, mp
    integer :: y, d, n, g

    !-----------------------------------------------------------------------
    !     	... Modify year and month numbers
    !-----------------------------------------------------------------------
    ap = year - (12 - month)/10
    mp = MOD( month-3,12 )
    if( mp < 0 ) then
       mp = mp + 12
    end if

    !-----------------------------------------------------------------------
    !     	... Julian day
    !-----------------------------------------------------------------------
    y = INT( days_per_year*( ap + 4712 ) )
    d = INT( days_per_month*mp + D0_5 )
    n = y + d + day  + 59
    g = INT( D0_75*INT( ap/100 + 49 ) ) - 38
    greg2jday = n - g

  end function greg2jday

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine set_met_rlx( )

    use pmgrid
    use hycoef, only: hypm, ps0

    integer :: k, k_cnt, k_top
    real(r8), parameter :: h0 = 7._r8   ! scale height (km)
    real(r8) :: p_top, p_bot

996 format( 'set_met_rlx: ',a15, I10 )
997 format( 'set_met_rlx: ',a15, E10.2 )
998 format( 'set_met_rlx: ',a15, PLEV(E10.2))
999 format( 'set_met_rlx: ',a15, PLEV(F10.5))

    met_rlx(:) = 999._r8

    p_top = ps0 * exp( - met_rlx_top/h0 )
    p_bot = ps0 * exp( - met_rlx_bot/h0 )

    if (masterproc) then
       write(iulog,fmt=997) 'p_top = ',p_top
       write(iulog,fmt=997) 'p_bot = ',p_bot
    endif

    if ( p_bot < hypm( pver-met_levels+1 ) .and. ( met_levels < pver ) ) then
       call endrun(  'set_met_rlx: met_rlx_bot is too high ' )
    endif

    where( hypm < p_top )
       met_rlx = 0._r8
    endwhere

    where( hypm > p_bot )
       met_rlx = met_max_rlx
    endwhere

    if ( any( met_rlx(:) /= met_max_rlx) ) then
       k_top = max(plev - met_levels, 1)

       do while ( met_rlx(k_top) /= 999._r8 )
          k_top = k_top + 1
          if ( k_top == pver ) then
             call endrun ( 'set_met_rlx: cannot find ramped region ')
          endif
       enddo

       met_rlx(1:k_top) = 0._r8

       k_cnt = count( met_rlx == 999._r8 )

       if (masterproc) then
          write(iulog,fmt=996) 'k_cnt = ',k_cnt
          write(iulog,fmt=996) 'k_top = ',k_top 
       endif

       do k = k_top,k_top+k_cnt
          met_rlx(k) = met_max_rlx*real( k - k_top ) / real(k_cnt)
       enddo
    endif

    if (masterproc) then
       write(iulog,fmt=996) '  met_levels = ',met_levels
       write(iulog,fmt=996) 'non-zero terms:',count( met_rlx /= 0._r8 )
    endif

    if ( met_levels < count( met_rlx /= 0._r8 ) ) then
       call endrun('set_met_rlx: met_rlx_top is too high for the meteorology data')
    endif

    if (masterproc) then
       write(iulog,fmt=998) 'press levels = ',hypm
       write(iulog,fmt=999) '     met_rlx = ',met_rlx
    endif

    if ( any( (met_rlx > 1._r8) .or. (met_rlx < 0._r8) ) ) then
      call endrun('Offline meteorology relaxation function not set correctly.')
    endif

  end subroutine set_met_rlx

end module metdata
