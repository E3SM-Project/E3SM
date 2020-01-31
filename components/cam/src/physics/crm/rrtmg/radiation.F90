module radiation

!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to RRTMG
!
! Revision history:
! May  2004, D. B. Coleman,  Initial version of interface module.
! July 2004, B. Eaton,       Use interfaces from new shortwave, longwave, and ozone modules.
! Feb  2005, B. Eaton,       Add namelist variables and control of when calcs are done.
! May  2008, Mike Iacono     Initial version for RRTMG
! June, 2009, Minghuai Wang,   MMF cam
!              The MMF treatment is added to the subroutine of radiation_tend. 
!              These modifications are based on the spcam3.5, which was developled
!              by Marat Khairoutdinov. The spcam3.5 only have one radiation package 
!              (camrt). See comments in radiation_tend for the details. 
!July, 2009, Minghuai Wang: 
!             For the Morrison's two momenent microphysics in SAM, droplet and ice crystal effective radius
!             used in the radiation code are calcualted at each CRM column by calling m2005_effradius
!October, 2009, Minghuai Wang:
!             CRM-scale aerosol water is used to calculate aerosol optical depth
!
! Nov  2010, J. Kay          Add COSP simulator calls
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8
use spmd_utils,      only: masterproc, iam, npes
use ppgrid,          only: pcols, pver, pverp, begchunk, endchunk
use physics_types,   only: physics_state, physics_ptend
use physconst,       only: cpair, cappa
use time_manager,    only: get_nstep, is_first_restart_step
use cam_abortutils,  only: endrun
use error_messages,  only: handle_err
use cam_control_mod, only: lambm0, obliqr, mvelpp, eccen
use scamMod,         only: scm_crm_mode, single_column,have_cld,cldobs,&
                           have_clwp,clwpobs,have_tg,tground,swrad_off,&
                           lwrad_off
use perf_mod,        only: t_startf, t_stopf
use cam_logfile,     only: iulog

use rad_constituents, only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info
use radconstants,     only: rrtmg_sw_cloudsim_band, rrtmg_lw_cloudsim_band, nswbands, nlwbands

implicit none
private
save

public :: &
   radiation_register,    &! registers radiation physics buffer fields
   radiation_nextsw_cday, &! calendar day of next radiation calculation
   radiation_do,          &! query which radiation calcs are done this timestep
   radiation_init,        &! calls radini
   radiation_readnl,      &! read radiation namelist
   radiation_tend          ! moved from radctl.F90

integer,public, allocatable :: cosp_cnt(:)       ! counter for cosp
integer,public              :: cosp_cnt_init = 0 !initial value for cosp counter

integer, public, parameter   :: kiss_seed_num = 4
integer, public, allocatable :: rad_randn_seedrst(:,:,:), tot_chnk_till_this_prc(:) !total number of chunks till this processor

! Private module data
integer :: qrs_idx      = 0 
integer :: qrl_idx      = 0 
integer :: su_idx       = 0 
integer :: sd_idx       = 0 
integer :: lu_idx       = 0 
integer :: ld_idx       = 0 
integer :: cldfsnow_idx = 0 
integer :: cld_idx      = 0 
integer :: concld_idx   = 0
integer :: rel_idx      = 0
integer :: rei_idx      = 0
integer :: dei_idx      = 0

! Default values for namelist variables

! Frequency of shortwave and longwave calculations in time steps (positive) or
! in hours (negative).
integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).

integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                           ! or hours (negative) SW/LW radiation will be
                           ! run continuously from the start of an
                           ! initial or restart run
logical :: spectralflux  = .false. ! calculate fluxes (up and down) per band.

logical :: use_rad_dt_cosz  = .false. ! if true, uses the radiation dt for all cosz calculations !BSINGH - Added for solar insolation calc.
character(len=16) :: microp_scheme  ! microphysics scheme

character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

logical :: dohirs = .false. ! diagnostic  brightness temperatures at the top of the
                            ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
                            ! channels (1,2,3,4).
integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

integer, allocatable :: clm_rand_seed(:,:,:)

real(r8) :: dt_avg=0.0_r8  ! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true !BSINGH - Added for solar insolation calc.

logical :: pergro_mods = .false. ! for activating pergro mods
integer :: firstblock, lastblock      ! global block indices

!===============================================================================
contains
!===============================================================================

subroutine radiation_readnl(nlfile, dtime_in)
!-------------------------------------------------------------------------------
! Purpose: Read radiation_nl namelist group.
!-------------------------------------------------------------------------------

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_logical, &
                              mpi_character, masterproc
   use time_manager,    only: get_step_size

   ! File containing namelist input
   character(len=*), intent(in) :: nlfile
   integer, intent(in), optional :: dtime_in

   ! Local variables
   integer :: unitn, ierr
   integer :: dtime  ! timestep size
   character(len=*), parameter :: subroutine_name = 'radiation_readnl'

   ! Variables defined in namelist
   namelist /radiation_nl/ iradsw, iradlw, irad_always, &
                           use_rad_dt_cosz, spectralflux

   ! Read the namelist, only if called from master process
   ! TODO: better documentation and cleaner logic here?
   if (masterproc) then
      unitn = getunit()
      open(unitn, file=trim(nlfile), status='old')
      call find_group_name(unitn, 'radiation_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, radiation_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subroutine_name // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   call mpibcast(spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
#endif

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (present(dtime_in)) then
      dtime = dtime_in
   else
      dtime  = get_step_size()
   end if
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Print runtime options to log.
   if (masterproc) then
      call radiation_printopts()
   end if

end subroutine radiation_readnl


  subroutine radiation_register
!-----------------------------------------------------------------------
! 
! Register radiation fields in the physics buffer
!
!-----------------------------------------------------------------------

    use physics_buffer,  only: pbuf_add_field, dtype_r8

    call pbuf_add_field('QRS' , 'global',dtype_r8,(/pcols,pver/), qrs_idx) ! shortwave radiative heating rate 
    call pbuf_add_field('QRL' , 'global',dtype_r8,(/pcols,pver/), qrl_idx) ! longwave  radiative heating rate 

    ! If the namelist has been configured for preserving the spectral fluxes, then create
    ! physics buffer variables to store the results.
    if (spectralflux) then
      call pbuf_add_field('SU'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), su_idx) ! shortwave upward flux (per band)
      call pbuf_add_field('SD'  , 'global',dtype_r8,(/pcols,pverp,nswbands/), sd_idx) ! shortwave downward flux (per band)
      call pbuf_add_field('LU'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), lu_idx) ! longwave upward flux (per band)
      call pbuf_add_field('LD'  , 'global',dtype_r8,(/pcols,pverp,nlwbands/), ld_idx) ! longwave downward flux (per band)
    end if

  end subroutine radiation_register

!================================================================================================
subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, irad_always_out, spectralflux_out, use_rad_dt_cosz_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   logical, intent(out), optional :: spectralflux_out
   logical, intent(out), optional :: use_rad_dt_cosz_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = -999
   if ( present(irad_always_out) ) irad_always_out = irad_always
   if ( present(spectralflux_out) ) spectralflux_out = spectralflux
   if ( present(use_rad_dt_cosz_out) ) use_rad_dt_cosz_out = use_rad_dt_cosz
end subroutine radiation_defaultopts

!================================================================================================


subroutine radiation_setopts(dtime, nhtfrq, iradsw_in, iradlw_in, iradae_in, &
   irad_always_in, spectralflux_in, use_rad_dt_cosz_in)
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! *** NOTE *** This routine needs information about dtime (init by dycore) 
!              and nhtfrq (init by history) to do its checking.  Being called
!              from runtime_opts provides these values possibly before they
!              have been set in the modules responsible for them.
!-----------------------------------------------------------------------

   integer, intent(in)           :: dtime           ! timestep size (s)
   integer, intent(in)           :: nhtfrq          ! output frequency of primary history file
   integer, intent(in), optional :: iradsw_in
   integer, intent(in), optional :: iradlw_in
   integer, intent(in), optional :: iradae_in
   integer, intent(in), optional :: irad_always_in
   logical, intent(in), optional :: spectralflux_in
   logical, intent(in), optional :: use_rad_dt_cosz_in
   
   ! Local
   integer :: ntspdy   ! no. timesteps per day
   integer :: nhtfrq1  ! local copy of input arg nhtfrq
   integer :: iradae   ! not used by RRTMG
!-----------------------------------------------------------------------

   if ( present(iradsw_in) )      iradsw = iradsw_in
   if ( present(iradlw_in) )      iradlw = iradlw_in
   if ( present(iradae_in) )      iradae = iradae_in
   if ( present(irad_always_in) ) irad_always = irad_always_in
   if ( present(spectralflux_in) ) spectralflux = spectralflux_in
   if ( present(use_rad_dt_cosz_in) ) use_rad_dt_cosz = use_rad_dt_cosz_in

   ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
   if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
   if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
   if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

   ! Has user specified iradae?
   if (iradae /= -999) then
      call endrun('radiation_setopts: iradae not used by RRTMG.')
   end if

end subroutine radiation_setopts

!===============================================================================

subroutine radiation_get(iradsw_out, iradlw_out, iradae_out, irad_always_out, spectralflux_out)
!----------------------------------------------------------------------- 
! Purpose: Provide access to private module data.  (This should be eliminated.)
!-----------------------------------------------------------------------

   integer, intent(out), optional :: iradsw_out
   integer, intent(out), optional :: iradlw_out
   integer, intent(out), optional :: iradae_out
   integer, intent(out), optional :: irad_always_out
   logical, intent(out), optional :: spectralflux_out
   !-----------------------------------------------------------------------

   if ( present(iradsw_out) )      iradsw_out = iradsw
   if ( present(iradlw_out) )      iradlw_out = iradlw
   if ( present(iradae_out) )      iradae_out = -999
   if ( present(irad_always_out) ) irad_always_out = irad_always
   if ( present(spectralflux_out) ) spectralflux_out = spectralflux_out

end subroutine radiation_get

!================================================================================================

subroutine radiation_printopts
!----------------------------------------------------------------------- 
! Purpose: Print runtime options to log.
!-----------------------------------------------------------------------


   if(irad_always /= 0) write(iulog,10) irad_always
   write(iulog,20) iradsw,iradlw
10 format(' Execute SW/LW radiation continuously for the first ',i5, ' timestep(s) of this run')
20 format(' Frequency of Shortwave Radiation calc. (IRADSW)     ',i5/, &
          ' Frequency of Longwave Radiation calc. (IRADLW)      ',i5)

end subroutine radiation_printopts

!================================================================================================

function radiation_do(op, timestep)
!----------------------------------------------------------------------- 
! Purpose: Returns true if the specified operation is done this timestep.
!-----------------------------------------------------------------------

   character(len=*), intent(in) :: op             ! name of operation
   integer, intent(in), optional:: timestep
   logical                      :: radiation_do   ! return value

   ! Local variables
   integer :: nstep             ! current timestep number
   !-----------------------------------------------------------------------

   if (present(timestep)) then
      nstep = timestep
   else
      nstep = get_nstep()
   end if

   select case (op)

   case ('sw') ! do a shortwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradsw == 1                     &
                    .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('lw') ! do a longwave heating calc this timestep?
      radiation_do = nstep == 0  .or.  iradlw == 1                     &
                    .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                    .or. nstep <= irad_always

   case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
      ! for RRTMG there is no abs/ems restart file
      radiation_do = .false.
         
   case default
      call endrun('radiation_do: unknown operation:'//op)

   end select
end function radiation_do

!================================================================================================

real(r8) function radiation_nextsw_cday()
  
!----------------------------------------------------------------------- 
! Purpose: Returns calendar day of next sw radiation calculation
!-----------------------------------------------------------------------

   use time_manager, only: get_curr_calday, get_nstep, get_step_size

   ! Local variables
   integer :: nstep      ! timestep counter
   logical :: dosw       ! true => do shosrtwave calc   
   integer :: offset     ! offset for calendar day calculation
   integer :: dTime      ! integer timestep size
   real(r8):: calday     ! calendar day of 
   !-----------------------------------------------------------------------

   radiation_nextsw_cday = -1._r8
   dosw   = .false.
   nstep  = get_nstep()
   dtime  = get_step_size()
   offset = 0
   do while (.not. dosw)
      nstep = nstep + 1
      offset = offset + dtime
      if (radiation_do('sw', nstep)) then
         radiation_nextsw_cday = get_curr_calday(offset=offset) 
         dosw = .true.
      end if
   end do
   if(radiation_nextsw_cday == -1._r8) then
      call endrun('error in radiation_nextsw_cday')
   end if
        
end function radiation_nextsw_cday

!================================================================================================

  subroutine radiation_init(phys_state)
!-----------------------------------------------------------------------
!
! Initialize the radiation parameterization, add fields to the history buffer
! 
!-----------------------------------------------------------------------
    use physics_buffer, only: pbuf_get_index
    use phys_grid,      only: npchunks, get_ncols_p, chunks, knuhcs, ngcols, dyn_to_latlon_gcol_map
    use cam_history,    only: addfld, horiz_only, add_default
    use constituents,   only: cnst_get_ind
    use physconst,      only: gravit, stebol, &
                              pstd, mwdry, mwco2, mwo3
    use phys_control,   only: phys_getopts
    use cospsimulator_intr, only: docosp, cospsimulator_intr_init
    use radsw,          only: radsw_init
    use radlw,          only: radlw_init
    use hirsbt,         only: hirsbt_init
    use hirsbtpar,      only: hirsname, msuname

    use radiation_data, only: init_rad_data
    use modal_aer_opt, only: modal_aer_opt_init
    use rrtmg_state,   only: rrtmg_state_init
    use time_manager,   only: get_step_size
    use dyn_grid,       only: get_block_bounds_d
#ifdef SPMD
    use mpishorthand,   only: mpi_integer, mpicom, mpi_comm_world
#endif
    use crmdims,        only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
    use rad_constituents, only: oldcldoptics
    use cloud_rad_props, only: cloud_rad_props_init

    type(physics_state), intent(in) :: phys_state(begchunk:endchunk)

    integer :: icall, nmodes
    logical :: active_calls(0:N_DIAG)
    integer :: nstep                       ! current timestep number
    logical :: history_amwg                ! output the variables used by the AMWG diag package
    logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
    logical :: history_budget              ! output tendencies and state variables for CAM4
                                           ! temperature, water vapor, cloud ice and cloud
                                           ! liquid budgets.
    integer :: history_budget_histfile_num ! output history file number for budget fields
    integer :: err
    integer :: dtime                        ! time step

    logical :: use_SPCAM                      ! SPCAM flag
    character(len=16) :: SPCAM_microp_scheme  ! SPCAM microphysics scheme

    !variables for pergro_mods
    character (len=250) :: errstr
    integer, allocatable, dimension(:,:,:) :: clm_id_mstr
    integer, allocatable, dimension(:,:) :: clm_id
    integer :: id, lchnk, ncol, ilchnk, astat, iseed, ipes, ipes_tmp
    integer :: igcol, chunkid, icol, iown, tot_cols, ierr, max_chnks_in_blk 
    !-----------------------------------------------------------------------
    
    call rrtmg_state_init()

    call init_rad_data() ! initialize output fields for offline driver

    call phys_getopts( use_SPCAM_out           = use_SPCAM           )
    call phys_getopts( SPCAM_microp_scheme_out = SPCAM_microp_scheme )
    call phys_getopts( microp_scheme_out       = microp_scheme       )

    call radsw_init()
    call radlw_init()

    ! Initialize cloud optics
    call cloud_rad_props_init()

    ! Set the radiation timestep for cosz calculations if requested using the adjusted iradsw value from radiation
    if (use_rad_dt_cosz)  then
       dtime  = get_step_size()
       dt_avg = iradsw*dtime
    end if

    call phys_getopts(history_amwg_out   = history_amwg,    &
                      history_vdiag_out  = history_vdiag,   &
                      history_budget_out = history_budget,  &
                      history_budget_histfile_num_out = history_budget_histfile_num, &
                      pergro_mods_out    = pergro_mods)
   
    ! Determine whether modal aerosols are affecting the climate, and if so
    ! then initialize the modal aerosol optics module
    call rad_cnst_get_info(0, nmodes=nmodes)
    if (nmodes > 0) call modal_aer_opt_init()

    call hirsbt_init()

    ! "irad_always" is number of time steps to execute radiation continuously from start of
    ! initial OR restart run

    nstep = get_nstep()
    if ( irad_always > 0) then
       nstep       = get_nstep()
       irad_always = irad_always + nstep
    end if


    if (docosp) call cospsimulator_intr_init

    
    allocate(cosp_cnt(begchunk:endchunk))

    !Modification needed by pergro_mods for generating random numbers
    if (pergro_mods) then
       max_chnks_in_blk = maxval(npchunks(:))  !maximum of the number for chunks in each procs
       allocate(clm_rand_seed(pcols,kiss_seed_num,max_chnks_in_blk), stat=astat)
       if( astat /= 0 ) then
          write(iulog,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_rand_seed; error = ',astat
          call endrun
       end if

       allocate(tot_chnk_till_this_prc(0:npes-1), stat=astat )
       if( astat /= 0 ) then
          write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate tot_chnk_till_this_prc variable; error = ',astat
          call endrun (errstr)
       end if
       
       !BSINGH - Build lat lon relationship to chunk and column
       !Compute maximum number of chunks each processor have
       if(masterproc) then
          tot_chnk_till_this_prc(0:npes-1) = huge(1)
          do ipes = 0, npes - 1
             tot_chnk_till_this_prc(ipes) = 0
             do ipes_tmp = 0, ipes-1
                tot_chnk_till_this_prc(ipes) = tot_chnk_till_this_prc(ipes) + npchunks(ipes_tmp)
             enddo
          enddo
       endif
#ifdef SPMD
       !BSINGH - Ideally we should use mpi_scatter but we are using this variable
       !in "if(masterproc)" below in phys_run1, so broadcast is iused here
       call mpibcast(tot_chnk_till_this_prc,npes, mpi_integer, 0, mpicom)
#endif
       call get_block_bounds_d(firstblock,lastblock)
       
       allocate(clm_id(pcols,max_chnks_in_blk), stat=astat)
       if( astat /= 0 ) then
          write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_id; error = ',astat
          call endrun(errstr)
       end if
       
       allocate(clm_id_mstr(pcols,max_chnks_in_blk,npes), stat=astat)
       if( astat /= 0 ) then
          write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_id_mstr; error = ',astat
          call endrun(errstr)
       end if
       !compute all clm ids on masterproc and then scatter it ....
       if(masterproc) then
          do igcol = 1, ngcols
             if (dyn_to_latlon_gcol_map(igcol) .ne. -1) then                
                chunkid  = knuhcs(igcol)%chunkid
                icol = knuhcs(igcol)%col
                iown  = chunks(chunkid)%owner
                ilchnk = (chunks(chunkid)%lcid - lastblock) - tot_chnk_till_this_prc(iown)
                clm_id_mstr(icol,ilchnk,iown+1) = igcol
             endif
          enddo
       endif
       
#ifdef SPMD
       !Scatter
       tot_cols = pcols*max_chnks_in_blk
       call MPI_Scatter( clm_id_mstr, tot_cols,  mpi_integer, &
            clm_id,    tot_cols,  mpi_integer, 0,             &
            mpicom,ierr)
#else
       !BSINGH - Haven't tested it.....               
       call endrun('radiation.F90(rrtmg)-radiation_init: non-mpi compiles are not tested yet for pergro test...')
#endif       
    endif
       
    if (is_first_restart_step()) then
       cosp_cnt(begchunk:endchunk)=cosp_cnt_init
       if (pergro_mods) then
          !--------------------------------------
          !Read seeds from restart file
          !--------------------------------------
          !For restart runs, rad_randn_seedrst array  will already be allocated in the restart_physics.F90
          
          do ilchnk = 1, max_chnks_in_blk
             lchnk = begchunk + (ilchnk -1)
             ncol = phys_state(lchnk)%ncol
             do iseed = 1, kiss_seed_num
                do icol = 1, ncol                
                   clm_rand_seed(icol,iseed,ilchnk) = rad_randn_seedrst(icol,iseed,lchnk)
                enddo
             enddo
          enddo
       endif
    else
       cosp_cnt(begchunk:endchunk)=0           
       if (pergro_mods) then
          !---------------------------------------
          !create seeds based off of column ids
          !---------------------------------------
          !allocate array rad_randn_seedrst for initial run for  maintaining exact restarts
          !For restart runs, it will already be allocated in the restart_physics.F90
          allocate(rad_randn_seedrst(pcols,kiss_seed_num,begchunk:endchunk), stat=astat)
          if( astat /= 0 ) then
             write(iulog,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate rad_randn_seedrst; error = ',astat
             call endrun
          end if
          do ilchnk = 1, max_chnks_in_blk
             lchnk = begchunk + (ilchnk -1)
             ncol = phys_state(lchnk)%ncol
             do iseed = 1, kiss_seed_num
                do icol = 1, ncol
                   id = clm_id(icol,ilchnk)
                   clm_rand_seed(icol,iseed,ilchnk) = id + (iseed -1)
                enddo
             enddo
          enddo
       endif
    end if


    ! Shortwave radiation
    call addfld('TOT_CLD_VISTAU', (/ 'lev' /), 'A',   '1', 'Total gbx cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
    call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Total in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
    call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Liquid in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
    call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', 'Ice in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)



    call add_default('TOT_CLD_VISTAU',  1, ' ')
    call add_default('TOT_ICLD_VISTAU', 1, ' ')

    ! get list of active radiation calls
    call rad_cnst_get_call_list(active_calls)

    do icall = 0, N_DIAG

       if (active_calls(icall)) then
          call addfld('SOLIN'//diag(icall),  horiz_only,     'A',   'W/m2', 'Solar insolation', sampling_seq='rad_lwsw')
          call addfld('SOLL'//diag(icall),  horiz_only,     'A',    'W/m2', 'Solar downward near infrared direct  to surface',&
           sampling_seq='rad_lwsw')
          call addfld('SOLS'//diag(icall),  horiz_only,     'A',    'W/m2', 'Solar downward visible direct  to surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('SOLLD'//diag(icall),  horiz_only,     'A',   'W/m2', 'Solar downward near infrared diffuse to surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('SOLSD'//diag(icall),  horiz_only,     'A',   'W/m2', 'Solar downward visible diffuse to surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('QRS'//diag(icall),   (/ 'lev' /),  'A',     'K/s', 'Solar heating rate', sampling_seq='rad_lwsw')
          call addfld('QRSC'//diag(icall),   (/ 'lev' /),  'A',    'K/s', 'Clearsky solar heating rate', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNS'//diag(icall),  horiz_only,     'A',    'W/m2', 'Net solar flux at surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNT'//diag(icall),  horiz_only,     'A',    'W/m2', 'Net solar flux at top of model', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNTOA'//diag(icall),  horiz_only,     'A',  'W/m2', 'Net solar flux at top of atmosphere', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSUTOA'//diag(icall),  horiz_only,     'A',  'W/m2', 'Upwelling solar flux at top of atmosphere', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNTOAC'//diag(icall),  horiz_only,     'A', 'W/m2', 'Clearsky net solar flux at top of atmosphere', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSUTOAC'//diag(icall),  horiz_only,     'A',  'W/m2', 'Clearsky upwelling solar flux at top of atmosphere', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSN200'//diag(icall),  horiz_only,     'A',  'W/m2', 'Net shortwave flux at 200 mb', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSN200C'//diag(icall),  horiz_only,     'A', 'W/m2', 'Clearsky net shortwave flux at 200 mb', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNTC'//diag(icall),  horiz_only,     'A',   'W/m2', 'Clearsky net solar flux at top of model', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSNSC'//diag(icall),  horiz_only,     'A',   'W/m2', 'Clearsky net solar flux at surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSDSC'//diag(icall),  horiz_only,     'A',   'W/m2', 'Clearsky downwelling solar flux at surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FSDS'//diag(icall),  horiz_only,     'A',    'W/m2', 'Downwelling solar flux at surface', &
                                                                                 sampling_seq='rad_lwsw')
          call addfld('FUS'//diag(icall),  (/ 'ilev' /), 'I',     'W/m2', 'Shortwave upward flux')
          call addfld('FDS'//diag(icall),  (/ 'ilev' /), 'I',     'W/m2', 'Shortwave downward flux')
          call addfld('FUSC'//diag(icall),  (/ 'ilev' /), 'I',    'W/m2', 'Shortwave clear-sky upward flux')
          call addfld('FDSC'//diag(icall),  (/ 'ilev' /), 'I',    'W/m2', 'Shortwave clear-sky downward flux')
          call addfld('FSNIRTOA'//diag(icall),  horiz_only,     'A','W/m2',&
           'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
          call addfld('FSNRTOAC'//diag(icall),  horiz_only,     'A','W/m2', &
                      'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', sampling_seq='rad_lwsw')
          call addfld('FSNRTOAS'//diag(icall),  horiz_only,     'A','W/m2', &
          'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', sampling_seq='rad_lwsw')
          call addfld ('SWCF'//diag(icall),  horiz_only,     'A',   'W/m2', 'Shortwave cloud forcing', sampling_seq='rad_lwsw')

          if (history_amwg) then
             call add_default('SOLIN'//diag(icall),   1, ' ')
             call add_default('QRS'//diag(icall),     1, ' ')
             call add_default('FSNS'//diag(icall),    1, ' ')
             call add_default('FSNT'//diag(icall),    1, ' ')
             call add_default('FSNTOA'//diag(icall),  1, ' ')
             call add_default('FSUTOA'//diag(icall),  1, ' ')
             call add_default('FSNTOAC'//diag(icall), 1, ' ')
             call add_default('FSUTOAC'//diag(icall), 1, ' ')
             call add_default('FSNTC'//diag(icall),   1, ' ')
             call add_default('FSNSC'//diag(icall),   1, ' ')
             call add_default('FSDSC'//diag(icall),   1, ' ')
             call add_default('FSDS'//diag(icall),    1, ' ')
             call add_default('SWCF'//diag(icall),    1, ' ')
          endif

       end if
    end do


    if (single_column .and. scm_crm_mode) then
       call add_default ('FUS     ', 1, ' ')
       call add_default ('FUSC    ', 1, ' ')
       call add_default ('FDS     ', 1, ' ')
       call add_default ('FDSC    ', 1, ' ')
    endif


    ! Longwave radiation

    do icall = 0, N_DIAG

       if (active_calls(icall)) then
          call addfld('QRL'//diag(icall),  (/ 'lev' /), 'A',     'K/s', 'Longwave heating rate', sampling_seq='rad_lwsw')
          call addfld('QRLC'//diag(icall),  (/ 'lev' /), 'A',    'K/s', 'Clearsky longwave heating rate', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLDS'//diag(icall), horiz_only,    'A',    'W/m2', 'Downwelling longwave flux at surface', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLDSC'//diag(icall), horiz_only,    'A',   'W/m2', 'Clearsky Downwelling longwave flux at surface', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLNS'//diag(icall), horiz_only,    'A',    'W/m2', 'Net longwave flux at surface', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLNT'//diag(icall), horiz_only,    'A',    'W/m2', 'Net longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLUT'//diag(icall), horiz_only,    'A',    'W/m2', 'Upwelling longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLUTC'//diag(icall), horiz_only,    'A',   'W/m2', 'Clearsky upwelling longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLNTC'//diag(icall), horiz_only,    'A',   'W/m2', 'Clearsky net longwave flux at top of model', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('LWCF'//diag(icall), horiz_only,    'A',    'W/m2', 'Longwave cloud forcing', sampling_seq='rad_lwsw')
          call addfld('FLN200'//diag(icall), horiz_only,    'A',  'W/m2', 'Net longwave flux at 200 mb', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLN200C'//diag(icall), horiz_only,    'A', 'W/m2', 'Clearsky net longwave flux at 200 mb', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FLNSC'//diag(icall), horiz_only,    'A',   'W/m2', 'Clearsky net longwave flux at surface', &
                                                                           sampling_seq='rad_lwsw')
          call addfld('FUL'//diag(icall), (/ 'ilev' /),'I',     'W/m2', 'Longwave upward flux')
          call addfld('FDL'//diag(icall), (/ 'ilev' /),'I',     'W/m2', 'Longwave downward flux')
          call addfld('FULC'//diag(icall), (/ 'ilev' /),'I',    'W/m2', 'Longwave clear-sky upward flux')
          call addfld('FDLC'//diag(icall), (/ 'ilev' /),'I',    'W/m2', 'Longwave clear-sky downward flux')
 
         if (history_amwg) then
            call add_default('QRL' //diag(icall),   1, ' ')
            call add_default('FLNS'//diag(icall),  1, ' ')
            call add_default('FLDS'//diag(icall),  1, ' ')
            call add_default('FLNT'//diag(icall),  1, ' ')
            call add_default('FLUT'//diag(icall),  1, ' ')
            call add_default('FLUTC'//diag(icall), 1, ' ')
            call add_default('FLNTC'//diag(icall), 1, ' ')
            call add_default('FLNSC'//diag(icall), 1, ' ')
            call add_default('LWCF'//diag(icall),  1, ' ')
         endif

         ! Add cloud-scale radiative quantities
         if (use_SPCAM) then
            call addfld ('CRM_QRAD', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'A', 'K/s', 'Radiative heating tendency')
            call addfld ('CRM_QRS ', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', 'CRM Shortwave radiative heating rate')
            call addfld ('CRM_QRSC', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', 'CRM Clearsky shortwave radiative heating rate')
            call addfld ('CRM_QRL ', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', 'CRM Longwave radiative heating rate' )
            call addfld ('CRM_QRLC', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', 'CRM Longwave radiative heating rate' )
            call addfld ('CRM_CLD_RAD', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'fraction', 'CRM cloud fraction' )
         end if
       end if
    end do

    call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave emissivity')

    if (single_column.and.scm_crm_mode) then
       call add_default ('FUL     ', 1, ' ')
       call add_default ('FULC    ', 1, ' ')
       call add_default ('FDL     ', 1, ' ')
       call add_default ('FDLC    ', 1, ' ')
    endif

    ! HIRS/MSU diagnostic brightness temperatures
    if (dohirs) then
       call addfld (hirsname(1),horiz_only,'A','K','HIRS CH2 infra-red brightness temperature')
       call addfld (hirsname(2),horiz_only,'A','K','HIRS CH4 infra-red brightness temperature')
       call addfld (hirsname(3),horiz_only,'A','K','HIRS CH6 infra-red brightness temperature')
       call addfld (hirsname(4),horiz_only,'A','K','HIRS CH8 infra-red brightness temperature')
       call addfld (hirsname(5),horiz_only,'A','K','HIRS CH10 infra-red brightness temperature')
       call addfld (hirsname(6),horiz_only,'A','K','HIRS CH11 infra-red brightness temperature')
       call addfld (hirsname(7),horiz_only,'A','K','HIRS CH12 infra-red brightness temperature')
       call addfld (msuname(1),horiz_only,'A','K','MSU CH1 microwave brightness temperature')
       call addfld (msuname(2),horiz_only,'A','K','MSU CH2 microwave brightness temperature')
       call addfld (msuname(3),horiz_only,'A','K','MSU CH3 microwave brightness temperature')
       call addfld (msuname(4),horiz_only,'A','K','MSU CH4 microwave brightness temperature')

       call add_default (hirsname(1), 1, ' ')
       call add_default (hirsname(2), 1, ' ')
       call add_default (hirsname(3), 1, ' ')
       call add_default (hirsname(4), 1, ' ')
       call add_default (hirsname(5), 1, ' ')
       call add_default (hirsname(6), 1, ' ')
       call add_default (hirsname(7), 1, ' ')
       call add_default (msuname(1), 1, ' ')
       call add_default (msuname(2), 1, ' ')
       call add_default (msuname(3), 1, ' ')
       call add_default (msuname(4), 1, ' ')
    end if

    ! Heating rate needed for d(theta)/dt computation
    call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

    if ( history_budget .and. history_budget_histfile_num > 1 ) then
       call add_default ('QRL     ', history_budget_histfile_num, ' ')
       call add_default ('QRS     ', history_budget_histfile_num, ' ')
    end if

    if (history_vdiag) then
       call add_default('FLUT', 2, ' ')
       call add_default('FLUT', 3, ' ')
    end if


    ! (Almost) net radiative flux at surface, does not have lwup.
    ! call addfld ('SRFRAD  ','W/m2    ',1,    'A','Net radiative flux at surface',phys_decomp)
    ! call add_default ('SRFRAD  ', 1, ' ')

    ! call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)

    if ( history_budget .and. history_budget_histfile_num > 1 ) then
       call add_default ('QRL     ', history_budget_histfile_num, ' ')
       call add_default ('QRS     ', history_budget_histfile_num, ' ')
    end if

     if (history_vdiag) then
       call add_default('FLUT', 2, ' ')
       call add_default('FLUT', 3, ' ')
    end if 

    cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)
    cld_idx      = pbuf_get_index('CLD')
    concld_idx   = pbuf_get_index('CONCLD')
    rel_idx      = pbuf_get_index('REL')
    rei_idx      = pbuf_get_index('REI')
    dei_idx      = pbuf_get_index('DEI')

    if (use_SPCAM .and. SPCAM_microp_scheme .eq. 'sam1mom') then
      cldfsnow_idx = 0
    end if

    if (cldfsnow_idx > 0) then
       call addfld ('CLDFSNOW',(/ 'lev' /),'I','1','CLDFSNOW',flag_xyfill=.true.)
       call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', 'Snow in-cloud extinction visible sw optical depth', &
                                                       sampling_seq='rad_lwsw', flag_xyfill=.true.)
    endif

    call addfld('COSZRS', horiz_only, 'I', '1', &
                'Cosine of solar zenith angle', &
                sampling_seq='rad_lwsw', flag_xyfill=.true.)

     ! Sanity check on cloud optics specified
     if (use_SPCAM .and. SPCAM_microp_scheme == 'sam1mom') then
        if (.not. oldcldoptics) then
           call endrun('radiation_init: must use oldcldoptics with sam1mom')
        end if
     end if

  end subroutine radiation_init

!===============================================================================
  
  subroutine radiation_tend( state, ptend,pbuf, &
       cam_out, cam_in, &
       landfrac,landm,icefrac,snowh, &
       fsns,    fsnt, flns,    flnt,  &
       fsds, net_flx, is_cmip6_volc)

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Driver for radiation computation.
    ! 
    ! Method: 
    ! Radiation uses cgs units, so conversions must be done from
    ! model fields to radiation fields.
    !
    ! Revision history:
    ! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
    ! 2004-08-09  B. Eaton         Add pointer variables for constituents.
    ! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
    ! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
    ! 2007-11-05  M. Iacono        Install rrtmg_lw and sw as radiation model.
    ! 2007-12-27  M. Iacono        Modify to use CAM cloud optical properties with rrtmg.
    ! 2009-06     Minghuai Wang,   add treatments for MMF CAM
    !              These modifications are based on the spcam3.5, which was developled
    !              by Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu). The spcam3.5 
    !              only have one radiation package (camrt). 
    !              Short wave and long wave radiation codes are called for every column
    !              of the CRM domain. CRM fields are named as *_crm, and domain-averaged fields are
    !              named as *_m. The domain-averaged fields and CRM fields are outputed at the end
    !              of the CRM domain loop (last_column=.true.).
    !              Several variables in state are updated with those from CRM output
    !              (liquid water, qc_rad; ice water, qi_rad; water vapor, qv_rad; 
    !              and temperature, t_rad). Several variables in pbuf are also updated 
    !              with those in CRM domain (cld, cicewp, cliqwp, csnowp, cldfsnow).  
    !              At the end of the radiation calculation, state and those in pbuf are
    !              restored to the old values. 
    !              Finally, a new cloud simulator are called, which takes account of cloud fileds 
    !              from the CRM domain. 
    ! 2009-07-13, Minghuai Wang: MMF CAM
    !             When Morrison's two momenent microphysics is used in SAM, droplet and ice crystal effective radius
    !             used in this radiation code are calcualted at each CRM column by calling m2005_effradius
    ! 2009-10-21, Minghuai Wang: MMF CAM
    !             CRM-scale aerosol water is used to calculate aerosol optical depth
    !----------------------------------------------------------------------------------------------


    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx, pbuf_get_index, pbuf_set_field
    
    use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
    use physics_types,   only: physics_state, physics_ptend
    use cospsimulator_intr, only: docosp, cospsimulator_intr_run,cosp_nradsteps
    use time_manager,    only: get_curr_calday
    use camsrfexch,      only: cam_out_t, cam_in_t
    use cam_history,     only: outfld
    use cam_history_support, only: fillvalue
    use parrrtm,         only: nbndlw
    use parrrsw,         only: nbndsw
    use hirsbt,          only: hirsrtm
    use hirsbtpar,       only: pnb_hirs, pnf_msu, hirsname, msuname
    use radheat,         only: radheat_tend
    use ppgrid
    use pspect
    use physconst,        only: cpair, stebol
    use radconstants,     only: nlwbands,idx_sw_diag
    use radsw,            only: rad_rrtmg_sw
    use radlw,            only: rad_rrtmg_lw
    use rad_constituents, only: rad_cnst_get_gas, rad_cnst_out, oldcldoptics, &
                                liqcldoptics, icecldoptics
    use aer_rad_props,    only: aer_rad_props_sw, aer_rad_props_lw
    use interpolate_data, only: vertinterp
    use cloud_rad_props,  only: get_ice_optics_sw, get_liquid_optics_sw, liquid_cloud_get_rad_props_lw, &
               ice_cloud_get_rad_props_lw, cloud_rad_props_get_lw, snow_cloud_get_rad_props_lw, get_snow_optics_sw
    use slingo,           only: slingo_liq_get_rad_props_lw, slingo_liq_optics_sw
    use ebert_curry,      only: ec_ice_optics_sw, ec_ice_get_rad_props_lw
    use rad_solar_var,    only: get_variability
    use radiation_data,   only: output_rad_data
    use rrtmg_state,      only: rrtmg_state_create, rrtmg_state_update, rrtmg_state_destroy, rrtmg_state_t, num_rrtmg_levs
    use crmdims,              only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
    use physconst,            only: gravit
    use constituents,         only: cnst_get_ind
    use radconstants,         only: idx_sw_diag
    use crm_physics,          only: m2005_effradius
#ifdef MODAL_AERO
    use modal_aero_data,       only: ntot_amode
#endif
    use phys_control,         only: phys_getopts
    use orbit,                only: zenith
    use output_aerocom_aie,   only: do_aerocom_ind3
    use pkg_cldoptics,        only: cldefr  ! for sam1mom microphysics


    ! Arguments
    logical,  intent(in)    :: is_cmip6_volc    ! true if cmip6 style volcanic file is read otherwise false 
    real(r8), intent(in)    :: landfrac(pcols)  ! land fraction
    real(r8), intent(in)    :: landm(pcols)     ! land fraction ramp
    real(r8), intent(in)    :: icefrac(pcols)   ! land fraction
    real(r8), intent(in)    :: snowh(pcols)     ! Snow depth (liquid water equivalent)
#ifdef MODAL_AERO
#endif
    real(r8), intent(inout) :: fsns(pcols)      ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)      ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)      ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)      ! Surface solar down flux
    real(r8), intent(inout) :: net_flx(pcols)

    type(physics_state), intent(inout), target :: state

    type(physics_ptend), intent(out)        :: ptend
    
    type(physics_buffer_desc), pointer      :: pbuf(:)
    type(cam_out_t),     intent(inout)      :: cam_out
    type(cam_in_t),      intent(in)         :: cam_in

    ! Local variables

    type(physics_state), target :: statein_copy
    logical :: dosw, dolw
    integer nstep                       ! current timestep number
    real(r8) britemp(pcols,pnf_msu)     ! Microwave brightness temperature
    real(r8) tb_ir(pcols,pnb_hirs)      ! Infrared brightness temperature
    real(r8) ts(pcols)                  ! surface temperature
    real(r8) pintmb(pcols,pverp)        ! Model interface pressures (hPa)
    real(r8) oro(pcols)                 ! Land surface flag, sea=0, land=1

    real(r8),pointer :: nc_rad(:,:,:,:) ! rad cloud water droplet number (#/kg)
    real(r8),pointer :: ni_rad(:,:,:,:) ! rad cloud ice crystal nubmer (#/kg)
    real(r8),pointer :: qs_rad(:,:,:,:) ! rad cloud snow crystal mass (kg/kg)
    real(r8),pointer :: ns_rad(:,:,:,:) ! rad cloud snow crystal number (#/kg)

    real(r8),pointer :: t_rad (:,:,:,:) ! rad temperuture
    real(r8),pointer :: qv_rad(:,:,:,:) ! rad vapor
    real(r8),pointer :: qc_rad(:,:,:,:) ! rad cloud water
    real(r8),pointer :: qi_rad(:,:,:,:) ! rad cloud ice
    real(r8),pointer :: cld_rad(:,:,:,:) ! 3D cloud fraction averaged over CRM integration
    real(r8),pointer :: crm_qrad(:,:,:,:) ! rad heating

    real(r8),pointer :: qaerwat_crm(:,:,:,:,:) ! aerosol water
    real(r8),pointer :: dgnumwet_crm(:,:,:,:,:) ! wet mode dimaeter

    integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
    real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
                                               !    maximally overlapped region.
                                               !    0->pmxrgn(i,1) is range of pressure for
                                               !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
                                               !    2nd region, etc
    real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
    real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables

    ! combined cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: c_cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    real(r8) :: c_cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    real(r8) :: c_cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    real(r8) :: c_cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    real(r8) :: c_cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: cld_tau    (nbndsw,pcols,pver) ! cloud extinction optical depth
    real(r8) :: cld_tau_w  (nbndsw,pcols,pver) ! cloud single scattering albedo * tau
    real(r8) :: cld_tau_w_g(nbndsw,pcols,pver) ! cloud assymetry parameter * w * tau
    real(r8) :: cld_tau_w_f(nbndsw,pcols,pver) ! cloud forward scattered fraction * w * tau
    real(r8) :: cld_lw_abs (nbndlw,pcols,pver) ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: ice_tau    (nbndsw,pcols,pver) ! ice extinction optical depth
    real(r8) :: ice_tau_w  (nbndsw,pcols,pver) ! ice single scattering albedo * tau
    real(r8) :: ice_tau_w_g(nbndsw,pcols,pver) ! ice assymetry parameter * tau * w
    real(r8) :: ice_tau_w_f(nbndsw,pcols,pver) ! ice forward scattered fraction * tau * w
    real(r8) :: ice_lw_abs (nbndlw,pcols,pver)   ! ice absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: snow_tau    (nbndsw,pcols,pver) ! snow extinction optical depth
    real(r8) :: snow_tau_w  (nbndsw,pcols,pver) ! snow single scattering albedo * tau
    real(r8) :: snow_tau_w_g(nbndsw,pcols,pver) ! snow assymetry parameter * tau * w
    real(r8) :: snow_tau_w_f(nbndsw,pcols,pver) ! snow forward scattered fraction * tau * w
    real(r8) :: snow_lw_abs (nbndlw,pcols,pver)   ! snow absorption optics depth (LW)
    real(r8) :: gb_snow_tau        (pcols,pver) ! grid-box mean snow_tau for COSP only
    real(r8) :: gb_snow_lw         (pcols,pver) ! grid-box mean LW snow optical depth for COSP only

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: liq_tau    (nbndsw,pcols,pver) ! liquid extinction optical depth
    real(r8) :: liq_tau_w  (nbndsw,pcols,pver) ! liquid single scattering albedo * tau
    real(r8) :: liq_tau_w_g(nbndsw,pcols,pver) ! liquid assymetry parameter * tau * w
    real(r8) :: liq_tau_w_f(nbndsw,pcols,pver) ! liquid forward scattered fraction * tau * w
    real(r8) :: liq_lw_abs (nbndlw,pcols,pver) ! liquid absorption optics depth (LW)

    real(r8) :: tot_cld_vistau(pcols,pver)  ! tot gbx cloud visible sw optical depth for output on history files
    real(r8) :: tot_icld_vistau(pcols,pver) ! tot in-cloud visible sw optical depth for output on history files
    real(r8) :: liq_icld_vistau(pcols,pver) ! liq in-cloud visible sw optical depth for output on history files
    real(r8) :: ice_icld_vistau(pcols,pver) ! ice in-cloud visible sw optical depth for output on history files
    real(r8) :: snow_icld_vistau(pcols,pver) ! snow in-cloud visible sw optical depth for output on history files

    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: rel      ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei      ! ice effective drop size (microns)
    real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cldfsnow ! cloud fraction of just "snow clouds- whatever they are"
    real(r8) :: cldfprime(pcols,pver)             ! combined cloud fraction (snow plus regular)
    real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction
    real(r8), pointer, dimension(:,:) :: qrs      ! shortwave radiative heating rate 
    real(r8), pointer, dimension(:,:) :: qrl      ! longwave  radiative heating rate 
    real(r8) :: qrsc(pcols,pver)                  ! clearsky shortwave radiative heating rate 
    real(r8) :: qrlc(pcols,pver)                  ! clearsky longwave  radiative heating rate 

    integer lchnk, ncol, lw
    real(r8) :: calday                        ! current calendar day
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)
    real(r8) :: coszrs(pcols)                 ! Cosine solar zenith angle
    logical  :: conserve_energy = .true.      ! flag to carry (QRS,QRL)*dp across time steps

    ! Local variables from radctl
    integer :: i, k, iseed, ilchnk                  ! index
    integer :: istat
    integer :: clm_seed (pcols,kiss_seed_num)
    real(r8) solin(pcols)         ! Solar incident flux
    real(r8) fsntoa(pcols)        ! Net solar flux at TOA
    real(r8) fsutoa(pcols)        ! Upwelling solar flux at TOA
    real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
    real(r8) fsutoac(pcols)       ! Clear sky upwelling solar flux at TOA
    real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
    real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
    real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
    real(r8) flut(pcols)          ! Upward flux at top of model
    real(r8) lwcf(pcols)          ! longwave cloud forcing
    real(r8) swcf(pcols)          ! shortwave cloud forcing
    real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
    real(r8) flntc(pcols)         ! Clear sky lw flux at model top
    real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
    real(r8) fldsc(pcols)         ! Clear sky lw flux at srf (down)
    real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) fns(pcols,pverp)     ! net shortwave flux
    real(r8) fcns(pcols,pverp)    ! net clear-sky shortwave flux
    real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
    real(r8) fsn200c(pcols)       ! fcns interpolated to 200 mb
    real(r8) fnl(pcols,pverp)     ! net longwave flux
    real(r8) fcnl(pcols,pverp)    ! net clear-sky longwave flux
    real(r8) qtot
    real(r8) factor_xy
    real(r8) cld_save   (pcols,pver)
    real(r8) fice       (pcols,pver)
    real(r8) cliqwp_crm (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) cicewp_crm (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) rel_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) rei_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) cld_tau_crm(pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) emis_crm   (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) qrl_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) qrs_crm    (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) qrlc_crm   (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) qrsc_crm   (pcols, crm_nx_rad, crm_ny_rad, crm_nz)
    real(r8) crm_fsnt   (pcols, crm_nx_rad, crm_ny_rad)   ! net shortwave fluxes at TOA at CRM grids
    real(r8) crm_fsntc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky shortwave fluxes at TOA at CRM grids
    real(r8) crm_fsns   (pcols, crm_nx_rad, crm_ny_rad)   ! net shortwave fluxes at surface at CRM grids
    real(r8) crm_fsnsc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky shortwave fluxes at surface at CRM grids
    real(r8) crm_flnt   (pcols, crm_nx_rad, crm_ny_rad)   ! net longwave fluxes at TOA at CRM grids
    real(r8) crm_flntc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky longwave fluxes at TOA at CRM grids
    real(r8) crm_flns   (pcols, crm_nx_rad, crm_ny_rad)   ! net longwave fluxes at surface at CRM grids
    real(r8) crm_flnsc  (pcols, crm_nx_rad, crm_ny_rad)   ! net clear-sky longwave fluxes at surface at CRM grids
    real(r8) crm_aodvisz(pcols, crm_nx_rad, crm_ny_rad, crm_nz)   ! layer aerosol optical depth at 550nm at CRM grids
    real(r8) crm_aodvis (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 550nm at CRM grids
    real(r8) crm_aod400 (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 400nm at CRM grids
    real(r8) crm_aod700 (pcols, crm_nx_rad, crm_ny_rad)   ! AOD at 700nm at CRM grids
    real(r8) aod400     (pcols)   ! AOD at 400nm at CRM grids
    real(r8) aod700     (pcols)   ! AOD at 700nm at CRM grids

    integer :: nct_tot_icld_vistau(pcols,pver) ! the number of CRM columns that has in-cloud visible sw optical depth 
    integer :: nct_liq_icld_vistau(pcols,pver) ! the number of CRM column that has liq in-cloud visible sw optical depth 
    integer :: nct_ice_icld_vistau(pcols,pver) ! the number of CRM column that has ice in-cloud visible sw optical depth 
    integer :: nct_snow_icld_vistau(pcols,pver) ! the number of CRM column that has snow in-cloud visible sw optical depth 

    real(r8) solin_m(pcols, 0:N_DIAG)         ! Solar incident flux
    real(r8) fsntoa_m(pcols, 0:N_DIAG)        ! Net solar flux at TOA
    real(r8) fsutoa_m(pcols, 0:N_DIAG)        ! upwelling solar flux at TOA
    real(r8) fsntoac_m(pcols, 0:N_DIAG)       ! Clear sky net solar flux at TOA
    real(r8) fsnirt_m(pcols, 0:N_DIAG)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc_m(pcols, 0:N_DIAG)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq_m(pcols, 0:N_DIAG)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc_m(pcols, 0:N_DIAG)         ! Clear sky total column abs solar flux
    real(r8) fsnsc_m(pcols, 0:N_DIAG)         ! Clear sky surface abs solar flux
    real(r8) fsdsc_m(pcols, 0:N_DIAG)         ! Clear sky surface downwelling solar flux
    real(r8) flut_m(pcols, 0:N_DIAG)          ! Upward flux at top of model
    real(r8) flutc_m(pcols, 0:N_DIAG)         ! Upward Clear Sky flux at top of model
    real(r8) flntc_m(pcols, 0:N_DIAG)         ! Clear sky lw flux at model top
    real(r8) flnsc_m(pcols, 0:N_DIAG)         ! Clear sky lw flux at srf (up-down)
    real(r8) fldsc_m(pcols, 0:N_DIAG)         ! Clear sky lw flux at srf (down)
    real(r8) flwds_m(pcols, 0:N_DIAG)          ! Down longwave flux at surface
    real(r8) fsns_m(pcols, 0:N_DIAG)          ! Surface solar absorbed flux
    real(r8) fsnt_m(pcols, 0:N_DIAG)          ! Net column abs solar flux at model top
    real(r8) flns_m(pcols, 0:N_DIAG)          ! Srf longwave cooling (up-down) flux
    real(r8) flnt_m(pcols, 0:N_DIAG)          ! Net outgoing lw flux at model top
    real(r8) fsds_m(pcols, 0:N_DIAG)          ! Surface solar down flux
    real(r8) fln200_m(pcols, 0:N_DIAG)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c_m(pcols, 0:N_DIAG)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) fsn200_m(pcols, 0:N_DIAG)        ! fns interpolated to 200 mb
    real(r8) fsn200c_m(pcols, 0:N_DIAG)       ! fcns interpolated to 200 mb
    real(r8) sols_m(pcols, 0:N_DIAG)          ! Solar downward visible direct  to surface
    real(r8) soll_m(pcols, 0:N_DIAG)          ! Solar downward near infrared direct  to surface
    real(r8) solsd_m(pcols, 0:N_DIAG)         ! Solar downward visible diffuse to surface
    real(r8) solld_m(pcols, 0:N_DIAG)         ! Solar downward near infrared diffuse to surface
    real(r8) qrs_m(pcols,pver, 0:N_DIAG)
    real(r8) qrl_m(pcols,pver, 0:N_DIAG)
    real(r8) qrsc_m(pcols,pver, 0:N_DIAG)
    real(r8) qrlc_m(pcols,pver, 0:N_DIAG)
    logical :: first_column
    logical :: last_column
    integer :: ii,jj,m
#ifdef MAML
    !define local cam_out fluxes for each CRM column  
    real(r8) lwup_loc 
    real(r8) :: sols_loc(pcols)
    real(r8) :: soll_loc(pcols)
    real(r8) :: solsd_loc(pcols)
    real(r8) :: solld_loc(pcols)
    real(r8) :: fsns_loc
    real(r8) :: flwds_loc(pcols)
#endif

    integer :: ixcldliq, ixcldice
    integer :: i_iciwp, i_iclwp, i_icswp
    real(r8), pointer, dimension(:, :) :: cicewp
    real(r8), pointer, dimension(:, :) :: cliqwp
    real(r8), pointer, dimension(:, :) :: csnowp
    real(r8) :: cicewp_save(pcols, pver)
    real(r8) :: cliqwp_save(pcols, pver)
    real(r8) :: csnowp_save(pcols, pver)
    real(r8) :: cldfsnow_save(pcols, pver)
    real(r8) :: rel_save(pcols, pver)
    real(r8) :: rei_save(pcols, pver)
    real(r8) :: effl    ! droplet effective radius [micrometer]
    real(r8) :: effi    ! ice crystal effective radius [micrometer]
    real(r8) :: effl_fn  ! effl for fixed number concentration of nlic = 1.e8

    real(r8) :: deffi    ! ice effective diameter for optics (radiation)
    real(r8) :: lamc     ! slope of droplet distribution for optics (radiation)
    real(r8) :: pgam     ! gamma parameter for optics (radiation)
    real(r8) :: dest     ! snow crystal effective diameters for optics (radiation) (micro-meter)
    real(r8), pointer, dimension(:, :) :: dei     ! ice effective diameter for optics (radiation)
    real(r8), pointer, dimension(:, :) :: mu      ! gamma parameter for optics (radiation)
    real(r8), pointer, dimension(:, :) :: lambdac ! slope of droplet distribution for optics (radiation)
    real(r8), pointer, dimension(:, :) :: des     ! snow crystal diameter for optics (mirometer, radiation)
    real(r8),allocatable ::  dei_save(:,:)
    real(r8),allocatable ::  mu_save(:,:)
    real(r8),allocatable ::  lambdac_save(:,:)
    real(r8),allocatable ::  des_save(:,:)
    real(r8),allocatable ::  dei_crm(:,:,:,:)     ! cloud scale ice effective diameter for optics
    real(r8),allocatable ::  mu_crm(:,:,:,:)      ! cloud scale gamma parameter for optics
    real(r8),allocatable ::  lambdac_crm(:,:,:,:) ! cloud scale slope of droplet distribution for optics
    real(r8),allocatable ::  des_crm(:,:,:,:)     ! cloud scale snow crystal diameter (micro-meter)
#ifdef MODAL_AERO
    real(r8), pointer, dimension(:,:,:) :: dgnumwet ! number mode diameter
    real(r8), pointer, dimension(:,:,:) :: qaerwat ! aerosol water 
    real(r8)  dgnumwet_save(pcols, pver, ntot_amode) 
    real(r8)  qaerwat_save(pcols, pver, ntot_amode)
#endif
 
    real(r8) pbr(pcols,pver)      ! Model mid-level pressures (dynes/cm2)
    real(r8) pnm(pcols,pverp)     ! Model interface pressures (dynes/cm2)
    real(r8) eccf                 ! Earth/sun distance factor
    real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

    real(r8) dy                   ! Temporary layer pressure thickness
    real(r8) tint(pcols,pverp)    ! Model interface temperature
    real(r8) :: sfac(1:nswbands)  ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band

    real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
    real(r8), dimension(pcols) :: co2_col_mean  ! co2 column mean mmr
    real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity

    real(r8), pointer, dimension(:,:,:) :: su => NULL()  ! shortwave spectral flux up
    real(r8), pointer, dimension(:,:,:) :: sd => NULL()  ! shortwave spectral flux down
    real(r8), pointer, dimension(:,:,:) :: lu => NULL()  ! longwave  spectral flux up
    real(r8), pointer, dimension(:,:,:) :: ld => NULL()  ! longwave  spectral flux down
    
    real(r8), allocatable :: su_m(:,:,:,:)  ! shortwave spectral flux up
    real(r8), allocatable :: sd_m(:,:,:,:)  ! shortwave spectral flux down
    real(r8), allocatable :: lu_m(:,:,:,:)  ! longwave  spectral flux up
    real(r8), allocatable :: ld_m(:,:,:,:)  ! longwave  spectral flux down

    ! Aerosol radiative properties
    real(r8) :: aer_tau    (pcols,0:pver,nbndsw) ! aerosol extinction optical depth
    real(r8) :: aer_tau_w  (pcols,0:pver,nbndsw) ! aerosol single scattering albedo * tau
    real(r8) :: aer_tau_w_g(pcols,0:pver,nbndsw) ! aerosol assymetry parameter * w * tau
    real(r8) :: aer_tau_w_f(pcols,0:pver,nbndsw) ! aerosol forward scattered fraction * w * tau
    real(r8) :: aer_lw_abs (pcols,pver,nbndlw)   ! aerosol absorption optics depth (LW)

    ! Gathered indicies of day and night columns 
    !  chunk_column_index = IdxDay(daylight_column_index)
    integer :: Nday                      ! Number of daylight columns
    integer :: Nnite                     ! Number of night columns
    integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
    integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

    integer :: icall                     ! index through climate/diagnostic radiation calls
    integer :: crm_nc_rad_idx, crm_ni_rad_idx, crm_qs_rad_idx, crm_ns_rad_idx
    integer :: crm_t_rad_idx, crm_qi_rad_idx, crm_qc_rad_idx, crm_qv_rad_idx, crm_qrad_idx
    integer :: dgnumwet_crm_idx, qaerwat_crm_idx

    logical :: active_calls(0:N_DIAG)
    logical :: use_SPCAM

    type(rrtmg_state_t), pointer :: r_state ! contains the atm concentratiosn in layers needed for RRTMG

! AeroCOM IND3 output +++mhwang
!MRN: Already defined above -- gives errors with GNU
!    real(r8) ::  aod400(pcols)        ! AOD at 400 nm 
!    real(r8) ::  aod700(pcols)        ! AOD at 700 nm
    real(r8) ::  angstrm(pcols)       ! Angstrom coefficient
    real(r8) ::  aerindex(pcols)      ! Aerosol index
    integer aod400_idx, aod700_idx, cld_tau_idx

    character(*), parameter :: name = 'radiation_tend'
    character(len=16)       :: SPCAM_microp_scheme  ! SPCAM_microphysics scheme
!----------------------------------------------------------------------
  
    call phys_getopts( use_SPCAM_out           = use_SPCAM )
    call phys_getopts( SPCAM_microp_scheme_out = SPCAM_microp_scheme)
    first_column = .false.
    last_column  = .false.

    ! In order to populate data structures with CRM state variables, we modify
    ! the physics_state object in-place and then restore the input physics_state
    ! at the end of the routine. So here we create a copy that we can restore.
    if (use_SPCAM) then 
      statein_copy = state
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
    endif

    lchnk = state%lchnk
    ncol  = state%ncol
    
    if(pergro_mods) then
       ilchnk = (lchnk - lastblock) - tot_chnk_till_this_prc(iam)
       clm_seed(1:pcols,1:kiss_seed_num) = clm_rand_seed (1:pcols,1:kiss_seed_num,ilchnk)       
    else
       ! For default simulation, clm_seed should never be used, assign it a value which breaks the simulation if used.
       clm_seed(1:pcols,1:kiss_seed_num) = huge(1)
    endif

    itim = pbuf_old_tim_idx()

    if (cldfsnow_idx > 0) then
       call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    endif
    call pbuf_get_field(pbuf, cld_idx,      cld,      start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, concld_idx,   concld,   start=(/1,1,itim/), kount=(/pcols,pver,1/)  )
    call pbuf_get_field(pbuf, qrs_idx,      qrs)
    call pbuf_get_field(pbuf, qrl_idx,      qrl)
    call pbuf_get_field(pbuf, rel_idx,      rel)
    call pbuf_get_field(pbuf, rei_idx,      rei)
    call pbuf_get_field(pbuf, dei_idx,      dei)

    if (spectralflux) then
      call pbuf_get_field(pbuf, su_idx, su)
      call pbuf_get_field(pbuf, sd_idx, sd)
      call pbuf_get_field(pbuf, lu_idx, lu)
      call pbuf_get_field(pbuf, ld_idx, ld)
      if(use_SPCAM) then
         allocate(su_m(pcols,pverp,nswbands,0:N_DIAG))
         allocate(sd_m(pcols,pverp,nswbands,0:N_DIAG))
         allocate(lu_m(pcols,pverp,nswbands,0:N_DIAG))
         allocate(ld_m(pcols,pverp,nswbands,0:N_DIAG))
      end if ! use_SPCAM
    end if
    
    if (use_SPCAM) then 
#ifdef MODAL_AERO
       dgnumwet_crm_idx = pbuf_get_index('CRM_DGNUMWET')
       qaerwat_crm_idx  = pbuf_get_index('CRM_QAERWAT')
       call pbuf_get_field(pbuf, dgnumwet_crm_idx, dgnumwet_crm)
       call pbuf_get_field(pbuf, qaerwat_crm_idx,  qaerwat_crm)
#endif

       if (SPCAM_microp_scheme .eq. 'm2005') then
          ! ifld = pbuf_get_index('DEI')
          ! call pbuf_get_field(pbuf, ifld, dei)
          ifld = pbuf_get_index('MU')
          call pbuf_get_field(pbuf, ifld, mu)
          ifld = pbuf_get_index('LAMBDAC')
          call pbuf_get_field(pbuf, ifld, LAMBDAC)
          ifld = pbuf_get_index('DES')
          call pbuf_get_field(pbuf, ifld, des)
       endif
    end if ! use_SPCAM
 
    if (do_aerocom_ind3) then
      cld_tau_idx = pbuf_get_index('cld_tau')
    end if
   
    ! For CRM, make cloud equal to input observations:
    if (single_column.and.scm_crm_mode.and.have_cld) then
      do k = 1,pver
         cld(:ncol,k)= cldobs(k)
      enddo
    endif

    if (cldfsnow_idx > 0) then
      call outfld('CLDFSNOW',cldfsnow,pcols,lchnk)
    endif

    ! Cosine solar zenith angle for current time step
    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    calday = get_curr_calday()
    call zenith (calday, clat, clon, coszrs, ncol, dt_avg)
    
    ! We can bypass the shortwave calculation by setting the cosine of the solar
    ! zenith angle to zero for all columns, because the shortwave code collapses
    ! the inputs to daytime-only arrays. In case the swrad_off flag is set then,
    ! we force shortwave to not be calculated by setting coszrs = 0.
    if (swrad_off) then
      coszrs(:)=0._r8
    endif    

    ! Output cosine solar zenith angle
    call outfld('COSZRS', coszrs(1:ncol), ncol, lchnk)

    ! The output_rad_data routine is intended to output all of the data needed
    ! to run the radiation code offline. This functionality may or may not be
    ! supported, and definitely is NOT for SP simulations.
    call output_rad_data(  pbuf, state, cam_in, landm, coszrs )

    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1,ncol
      if ( coszrs(i) > 0.0_r8 ) then
        Nday = Nday + 1
        IdxDay(Nday) = i
      else
        Nnite = Nnite + 1
        IdxNite(Nnite) = i
      end if
    end do ! i

    ! Allocate "save" variables that will be used to restore fields that are
    ! modified in-place in pbuf to populate with each crm column one at a time
    if (use_SPCAM) then
      allocate(dei_save(pcols, pver))
      allocate(dei_crm(pcols, crm_nx_rad, crm_ny_rad, crm_nz))
      if (SPCAM_microp_scheme .eq. 'm2005') then 
        allocate(mu_save      (pcols, pver))
        allocate(lambdac_save (pcols, pver))
        allocate(des_save     (pcols, pver))
        allocate(mu_crm       (pcols, crm_nx_rad, crm_ny_Rad, crm_nz))
        allocate(lambdac_crm  (pcols, crm_nx_rad, crm_ny_Rad, crm_nz))
        allocate(des_crm      (pcols, crm_nx_rad, crm_ny_Rad, crm_nz))
      end if
    end if

    ! Figure out if we are doing radiation at this timestep. For SP-CAM, these
    ! should ALWAYS return true...this is handled in radiation_init() by setting
    ! iradsw = iradlw = 1
    dosw = radiation_do('sw')      ! do shortwave heating calc this timestep?
    dolw = radiation_do('lw')      ! do longwave heating calc this timestep?

    ! Initialize averages over CRM columns to zero. These are aggregated over
    ! the loop over CRM columns below.
    if (use_SPCAM) then 

       ! Get CRM radiative heating from the physics buffer; we need to do this regardless of whether
       ! or not we are going to do radiative calculations this timestep, because this is still
       ! accessed outside the dosw .or. dolw logical block.
       crm_qrad_idx = pbuf_get_index('CRM_QRAD')
       call pbuf_get_field(pbuf, crm_qrad_idx, crm_qrad)

       ! Only zero SP fields when we are going to update the longwave or
       ! shortwave in case we are NOT going to update the radiation each
       ! timestep.
       if (dosw .or. dolw) then
         solin_m    = 0.   ; fsntoa_m   = 0. 
         fsutoa_m   = 0.   ; fsntoac_m  = 0.
         fsnirt_m   = 0.   ; fsnrtc_m   = 0.
         fsnirtsq_m = 0.   ; fsntc_m    = 0. 
         fsnsc_m    = 0.   ; fsdsc_m    = 0.
         flut_m     = 0.   ; flutc_m    = 0.
         flntc_m    = 0.   ; flnsc_m    = 0.
         fldsc_m    = 0.   ; flwds_m    = 0.
         fsns_m     = 0.   ; fsnt_m     = 0.
         flns_m     = 0.   ; flnt_m     = 0.
         fsds_m     = 0.
         fln200_m   = 0.   ; fln200c_m  = 0.
         fsn200_m   = 0.   ; fsn200c_m  = 0.
         sols_m     = 0.   ; soll_m     = 0.
         solsd_m    = 0.   ; solld_m    = 0.
         qrs_m      = 0.   ; qrl_m      = 0.
         qrsc_m     = 0.   ; qrlc_m     = 0.
         qrs_crm    = 0.   ; qrl_crm    = 0.
         qrsc_crm    = 0.  ; qrlc_crm   = 0.
         emis_crm   = 0.   ; cld_tau_crm= 0.
         crm_aodvisz= 0.   ; crm_aodvis = 0.
         crm_aod400 = 0.   ; crm_aod700 = 0.
         aod400     = 0.   ; aod700     = 0.
         crm_fsnt   = 0.   ; crm_fsntc  = 0. 
         crm_fsns   = 0.   ; crm_fsnsc  = 0.
         crm_flnt   = 0.   ; crm_flntc  = 0.
         crm_flns   = 0.   ; crm_flnsc  = 0.
         tot_cld_vistau   = 0
         tot_icld_vistau  = 0. ; nct_tot_icld_vistau  = 0.
         liq_icld_vistau  = 0. ; nct_liq_icld_vistau  = 0.
         ice_icld_vistau  = 0. ; nct_ice_icld_vistau  = 0.
         snow_icld_vistau = 0. ; nct_snow_icld_vistau = 0.
         if (spectralflux) then
           su_m = 0. ; sd_m = 0.
           lu_m = 0. ; ld_m = 0.
         end if

         i_iciwp  = pbuf_get_index('ICIWP')
         i_iclwp  = pbuf_get_index('ICLWP')
         i_icswp  = pbuf_get_index('ICSWP')
         call pbuf_get_field(pbuf, i_iciwp, cicewp)
         call pbuf_get_field(pbuf, i_iclwp, cliqwp)
         ! call pbuf_get_field(pbuf, i_icswp, csnowp)
         cicewp_save = cicewp     ! save to restore later
         cliqwp_save = cliqwp     ! save to restore later
         ! csnowp_save = csnowp     ! save to restore later

         if (SPCAM_microp_scheme .eq. 'm2005') then 
           call pbuf_get_field(pbuf, i_icswp, csnowp)
           csnowp_save = csnowp     ! save to restore later
         end if

         crm_t_rad_idx   = pbuf_get_index('CRM_T_RAD')
         crm_qc_rad_idx  = pbuf_get_index('CRM_QC_RAD')
         crm_qi_rad_idx  = pbuf_get_index('CRM_QI_RAD')
         crm_qv_rad_idx  = pbuf_get_index('CRM_QV_RAD')
         call pbuf_get_field(pbuf, crm_t_rad_idx,  t_rad)
         call pbuf_get_field(pbuf, crm_qc_rad_idx, qc_rad)
         call pbuf_get_field(pbuf, crm_qi_rad_idx, qi_rad)
         call pbuf_get_field(pbuf, crm_qv_rad_idx, qv_rad)

         ! Zero out radiative heating
         crm_qrad=0.

         if (SPCAM_microp_scheme .eq. 'm2005') then 
           crm_nc_rad_idx  = pbuf_get_index('CRM_NC_RAD')
           call pbuf_get_field(pbuf, crm_nc_rad_idx, nc_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
           crm_ni_rad_idx  = pbuf_get_index('CRM_NI_RAD')
           call pbuf_get_field(pbuf, crm_ni_rad_idx, ni_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
           crm_qs_rad_idx  = pbuf_get_index('CRM_QS_RAD')
           call pbuf_get_field(pbuf, crm_qs_rad_idx, qs_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
           crm_ns_rad_idx  = pbuf_get_index('CRM_NS_RAD')
           call pbuf_get_field(pbuf, crm_ns_rad_idx, ns_rad, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
         endif

         ! Get cloud fraction averaged over the CRM time integration
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_CLD_RAD'), cld_rad)
         call outfld('CRM_CLD_RAD', cld_rad, state%ncol, state%lchnk)

         cicewp(1:ncol,1:pver) = 0.  
         cliqwp(1:ncol,1:pver) = 0.

         factor_xy = 1./real( crm_nx_rad*crm_ny_rad ,r8)

         cld_save = cld  ! save to restore later
         rel_save = rel  ! save to restroe later
         rei_save = rei  ! save to restore later
         dei_save = dei  ! save to restore later
         cld = 0.0_r8
         rel = 0.0_r8
         rei = 0.0_r8
         dei = 0.0_r8
         if (cldfsnow_idx > 0) then
           cldfsnow  = 0.0_r8
           cldfsnow_save = cldfsnow
         end if
         if (SPCAM_microp_scheme .eq. 'm2005') then 
           mu_save      = mu
           lambdac_save = lambdac
           des_save     = des
           mu      = 0.0_r8
           lambdac = 0.0_r8
           des     = 0.0_r8
         endif

#ifdef MODAL_AERO
         ifld = pbuf_get_index('DGNUMWET')
         call pbuf_get_field(pbuf, ifld, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/) )
         ifld  = pbuf_get_index( 'QAERWAT' )
         call pbuf_get_field(pbuf, ifld, qaerwat, start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/) )
         dgnumwet_save = dgnumwet
         qaerwat_save = qaerwat
#endif /*MODAL_AERO*/
      end if ! dosw .or. dolw
    endif ! SPCAM

    if (dosw .or. dolw) then

      ! construct an RRTMG state object
      r_state => rrtmg_state_create( state, cam_in )

      ! For CRM, make cloud liquid water path equal to input observations
      if(single_column.and.scm_crm_mode.and.have_clwp)then
        call endrun('cloud water path must be passed through radiation interface')
        !do k=1,pver
        !   cliqwp(:ncol,k) = clwpobs(k)
        !end do
      endif

      ! calculate effective radius - moved outside of ii,jj loops for 1-moment microphysics
      if (SPCAM_microp_scheme .eq. 'sam1mom') then 
        call cldefr( lchnk, ncol, landfrac, state%t, rel, rei, state%ps, state%pmid, landm, icefrac, snowh )
      end if

      ! Start loop over CRM columns; the strategy here is to loop over each CRM
      ! column and separately call the radiative transfer codes with optical
      ! properties calculated from CRM fields for each of those columns. Note
      ! that here we loop over "crm_nx_rad" rather than "crm_nx". This is to
      ! allow the flexibility for the radiation to be calculated on a reduced
      ! resolution relative to the CRM by grouping (averaging) adjacent columns
      ! together.
      do jj=1,crm_ny_rad
        do ii=1,crm_nx_rad

          if (use_SPCAM) then 
            first_column = ii.eq.1.and.jj.eq.1
            last_column = ii.eq.crm_nx_rad.and.jj.eq.crm_ny_rad

            do m=1,crm_nz
              k = pver-m+1
              do i=1,ncol

                ! Overwrite cloud fraction with CRM cloud fraction
                cld(i,k) = cld_rad(i,ii,jj,m)

                ! Calculate water paths and fraction of ice
                if (cld(i,k) > 0) then
                  qtot = qc_rad(i,ii,jj,m) + qi_rad(i,ii,jj,m)
                  fice(i,k) = qi_rad(i,ii,jj,m)/qtot
                  cicewp(i,k) = qi_rad(i,ii,jj,m)*state%pdel(i,k)/gravit    &
                           / max(0.01_r8,cld(i,k)) ! In-cloud ice water path.
                  cliqwp(i,k) = qc_rad(i,ii,jj,m)*state%pdel(i,k)/gravit     &
                           / max(0.01_r8,cld(i,k)) ! In-cloud liquid water path. 
                else
                  fice(i,k)= 0.
                  cicewp(i,k) = 0.           ! In-cloud ice water path.
                  cliqwp(i,k) = 0.           ! In-cloud liquid water path.
                end if

                ! snow water-related variables: 
                ! snow water is an important component in m2005 microphysics, and is therefore taken
                ! into account in the radiative calculation (snow water path is several times larger
                ! than ice water path in m2005 globally). 
                if (SPCAM_microp_scheme .eq. 'm2005') then 
                  if( qs_rad(i, ii, jj, m).gt.1.0e-7) then
                    cldfsnow(i,k) = 0.99_r8   
                    csnowp(i,k) = qs_rad(i,ii,jj,m)*state%pdel(i,k)/gravit    &
                           / max(0.001_r8,cldfsnow(i,k)) ! In-cloud ice water path.
                  else
                    cldfsnow(i,k) = 0.0  
                    csnowp(i,k) = 0.0
                  end if
                end if

                ! Update ice water, liquid water, water vapor, and temperature in state
                state%q(i,k,ixcldice) =  qi_rad(i,ii,jj,m)
                state%q(i,k,ixcldliq) =  qc_rad(i,ii,jj,m)
                state%q(i,k,1)        =  max(1.e-9_r8,qv_rad(i,ii,jj,m))
                state%t(i,k)          =  t_rad(i, ii, jj, m)
#ifdef MODAL_AERO
                ! Use CRM scale aerosol water to calculate aerosol optical depth. Here we assume no
                ! aerosol water uptake for cloudy sky on CRM grids. This is not really physically
                ! correct, but if we assume 100% of relative humidity for aerosol water uptake, this
                ! will bias 'AODVIS' to be large, since 'AODVIS' is used to compare with observed
                ! clear sky AOD. In the future, AODVIS should be calculated from clear sky CRM AOD
                ! only. But before this is done, we will assume no water uptake on CRM grids for 
                ! cloudy conditions (The radiative effects of this assumption will be small, since 
                ! aerosol effects are small relative to cloud effects for cloudy sky anyway. 
                ! -Minghuai Wang (minghuai.wang@pnl.gov)
                qaerwat (i, k, 1:ntot_amode) =  qaerwat_crm(i, ii, jj, m, 1:ntot_amode)
                dgnumwet(i, k, 1:ntot_amode) = dgnumwet_crm(i, ii, jj, m, 1:ntot_amode)
#endif 
              end do ! i
            end do ! m

            ! update effective radius
            if (SPCAM_microp_scheme .eq. 'm2005') then 
              do m=1,crm_nz
                k = pver-m+1
                do i=1,ncol
                  call m2005_effradius(qc_rad(i,ii,jj,m), nc_rad(i,ii,jj,m), qi_rad(i,ii,jj,m), &
                                       ni_rad(i,ii,jj,m), qs_rad(i,ii,jj,m), ns_rad(i,ii,jj,m),  &
                                       1.0_r8, state%pmid(i,k), state%t(i,k), effl, effi, effl_fn, deffi, lamc, pgam, dest)
                  rel(i,k)     = effl
                  rei(i,k)     = effi
                  dei(i,k)     = deffi 
                  mu(i,k)      = pgam 
                  lambdac(i,k) = lamc 
                  des(i,k)     = dest
                  dei_crm(i,ii,jj,m)     = dei(i,k)
                  mu_crm(i,ii,jj,m)      = mu(i,k)
                  lambdac_crm(i,ii,jj,m) = lambdac(i,k)
                  des_crm(i,ii,jj,m)     = des(i,k)
                
                  rel_crm(i,ii,jj,m) = rel(i,k)
                  rei_crm(i,ii,jj,m) = rei(i,k)
                end do ! i
              end do ! m
            else if (SPCAM_microp_scheme .eq. 'sam1mom') then 
              ! for sam1mom, rel and rei are calcualted above, and are the same for all CRM columns
              do m=1,crm_nz
                k = pver-m+1
                rel_crm(:ncol,ii,jj,m)=rel(:ncol,k)
                rei_crm(:ncol,ii,jj,m)=rei(:ncol,k)
                
                dei(:ncol,k) = rei(:ncol,k) * 2.0_r8
                ! whannah - calculation of dei below is taken from m2005_effradius()
                ! dei(:ncol,k) = rei(:ncol,k) * 500._r8/917._r8 * 2._r8
                dei_crm(:ncol,ii,jj,m) = dei(:ncol,k)
              end do ! m
            end if ! sam1mom
          endif ! use_SPCAM

          call t_startf('cldoptics')
          if (dosw) then
            if(oldcldoptics) then
              call ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.false.)
              call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.false.)
            else
              select case (icecldoptics)
              case ('ebertcurry')
                call  ec_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp=.true.)
              case ('mitchell')
                call get_ice_optics_sw(state, pbuf, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
              case default
                call endrun('iccldoptics must be one either ebertcurry or mitchell')
              end select
              select case (liqcldoptics)
              case ('slingo')
                call slingo_liq_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp=.true.)
              case ('gammadist')
                call get_liquid_optics_sw(state, pbuf, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
              case default
                call endrun('liqcldoptics must be either slingo or gammadist')
              end select
            endif
            cld_tau    (:,1:ncol,:) = liq_tau    (:,1:ncol,:) + ice_tau    (:,1:ncol,:)
            cld_tau_w  (:,1:ncol,:) = liq_tau_w  (:,1:ncol,:) + ice_tau_w  (:,1:ncol,:)
            cld_tau_w_g(:,1:ncol,:) = liq_tau_w_g(:,1:ncol,:) + ice_tau_w_g(:,1:ncol,:)
            cld_tau_w_f(:,1:ncol,:) = liq_tau_w_f(:,1:ncol,:) + ice_tau_w_f(:,1:ncol,:)
 
            if (cldfsnow_idx > 0) then
              ! add in snow
              call get_snow_optics_sw(state, pbuf, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)
              do i=1,ncol
                do k=1,pver
                  cldfprime(i,k)=max(cld(i,k),cldfsnow(i,k))
                  if(cldfprime(i,k) > 0.)then
                    c_cld_tau    (1:nbndsw,i,k) = (cldfsnow(i,k)*snow_tau    (1:nbndsw,i,k) &
                                                      + cld(i,k)*cld_tau     (1:nbndsw,i,k))/cldfprime(i,k)
                    c_cld_tau_w  (1:nbndsw,i,k) = (cldfsnow(i,k)*snow_tau_w  (1:nbndsw,i,k) &
                                                      + cld(i,k)*cld_tau_w   (1:nbndsw,i,k))/cldfprime(i,k)
                    c_cld_tau_w_g(1:nbndsw,i,k) = (cldfsnow(i,k)*snow_tau_w_g(1:nbndsw,i,k) &
                                                      + cld(i,k)*cld_tau_w_g (1:nbndsw,i,k))/cldfprime(i,k)
                    c_cld_tau_w_f(1:nbndsw,i,k) = (cldfsnow(i,k)*snow_tau_w_f(1:nbndsw,i,k) &
                                                      + cld(i,k)*cld_tau_w_f (1:nbndsw,i,k))/cldfprime(i,k)
                  else
                    c_cld_tau    (1:nbndsw,i,k) = 0._r8
                    c_cld_tau_w  (1:nbndsw,i,k) = 0._r8
                    c_cld_tau_w_g(1:nbndsw,i,k) = 0._r8
                    c_cld_tau_w_f(1:nbndsw,i,k) = 0._r8
                  endif
                enddo
              enddo
              if (use_SPCAM) then 
                do m=1,crm_nz
                   k = pver-m+1
                   do i=1,ncol
                      cld_tau_crm(i,ii,jj,m) =  cld_tau(rrtmg_sw_cloudsim_band,i,k)
                   end do ! i
                end do ! m
              endif 
            else  ! cldfsnow_idx > 0
              c_cld_tau    (1:nbndsw,1:ncol,:) = cld_tau    (:,1:ncol,:)
              c_cld_tau_w  (1:nbndsw,1:ncol,:) = cld_tau_w  (:,1:ncol,:)
              c_cld_tau_w_g(1:nbndsw,1:ncol,:) = cld_tau_w_g(:,1:ncol,:)
              c_cld_tau_w_f(1:nbndsw,1:ncol,:) = cld_tau_w_f(:,1:ncol,:)
            endif  ! cldfsnow_idx > 0

            if(do_aerocom_ind3) then
              call pbuf_set_field(pbuf,cld_tau_idx,cld_tau(rrtmg_sw_cloudsim_band, :, :))                   
            end if

          endif

          if (dolw) then
            if(oldcldoptics) then
              call cloud_rad_props_get_lw(state, pbuf, cld_lw_abs, oldcloud=.true.)
            else  ! oldcldoptics
              select case (icecldoptics)
              case ('ebertcurry')
                call    ec_ice_get_rad_props_lw(state, pbuf, ice_lw_abs, oldicewp=.true.)
              case ('mitchell')
                call ice_cloud_get_rad_props_lw(state, pbuf, ice_lw_abs)
              case default
                call endrun('iccldoptics must be one either ebertcurry or mitchell')
              end select
              select case (liqcldoptics)
              case ('slingo')
                call   slingo_liq_get_rad_props_lw(state, pbuf, liq_lw_abs, oldliqwp=.true.)
              case ('gammadist')
                call liquid_cloud_get_rad_props_lw(state, pbuf, liq_lw_abs)
              case default
                call endrun('liqcldoptics must be either slingo or gammadist')
              end select
              cld_lw_abs(:,1:ncol,:) = liq_lw_abs(:,1:ncol,:) + ice_lw_abs(:,1:ncol,:)
            endif  ! oldcldoptics

            if (cldfsnow_idx > 0) then
              ! add in snow
              call snow_cloud_get_rad_props_lw(state, pbuf, snow_lw_abs)
              do i=1,ncol
                do k=1,pver
                   cldfprime(i,k)=max(cld(i,k),cldfsnow(i,k))
                   if(cldfprime(i,k) > 0.)then
                      c_cld_lw_abs(1:nbndlw,i,k)= &
                           (cldfsnow(i,k)*snow_lw_abs(1:nbndlw,i,k) + cld(i,k)*cld_lw_abs(1:nbndlw,i,k))/cldfprime(i,k)
                   else
                      c_cld_lw_abs(1:nbndlw,i,k)= 0._r8
                   endif
                enddo
              enddo
              if (use_SPCAM) then
                  do m=1,crm_nz
                     k = pver-m+1
                     do i=1,ncol
                        emis_crm(i,ii,jj,m)=1._r8 - exp(-cld_lw_abs(rrtmg_lw_cloudsim_band,i,k))
                     end do ! i
                  end do ! m
              endif  ! use_SPCAM
            else  ! cldfsnow_idx > 0
               c_cld_lw_abs(1:nbndlw,1:ncol,:)=cld_lw_abs(:,1:ncol,:)
            endif  ! cldfsnow_idx > 0
          endif  ! dolw

          if (.not.(cldfsnow_idx > 0)) then
            cldfprime(1:ncol,:)=cld(1:ncol,:)
          endif

          call t_stopf('cldoptics')

          ! construct cgs unit reps of pmid and pint and get "eccf" - earthsundistancefactor
          call radinp(ncol, state%pmid, state%pint, pbr, pnm, eccf)

          ! Calculate interface temperatures (following method
          ! used in radtpl for the longwave), using surface upward flux and
          ! stebol constant in mks units
          do i = 1,ncol
            tint(i,1) = state%t(i,1)
#ifdef MAML
            tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i,ii)/stebol))
#else
            tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i)/stebol))
#endif
            do k = 2,pver
               dy = (state%lnpint(i,k) - state%lnpmid(i,k)) / (state%lnpmid(i,k-1) - state%lnpmid(i,k))
               tint(i,k) = state%t(i,k) - dy * (state%t(i,k) - state%t(i,k-1))
            end do
          end do

          ! Solar radiation computation
          if (dosw) then
            call t_startf ('rad_sw')

            ! Calculate solar variability factor
            call get_variability(sfac)

            ! Get the active climate/diagnostic shortwave calculations
            call rad_cnst_get_call_list(active_calls)

            ! Loop over diagnostic cases (each of which can contain different
            ! radiative constituents. The climate (icall==0) calculation must 
            ! occur last, so we loop from N_DIAG to 0.
            do icall = N_DIAG, 0, -1

              if (active_calls(icall)) then

                ! Update the concentrations in the RRTMG state object
                call rrtmg_state_update( state, pbuf, icall, r_state )

                ! Calculate the aerosol optical properties
                call aer_rad_props_sw( icall, state, pbuf, nnite, idxnite, is_cmip6_volc, &
                                       aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)

                ! Run the shortwave radiation driver
                call t_startf ('rad_rrtmg_sw')
                call rad_rrtmg_sw( &
                     lchnk,        ncol,         num_rrtmg_levs, r_state,                    &
                     state%pmid,   cldfprime,                                                &
                     aer_tau,      aer_tau_w,    aer_tau_w_g,  aer_tau_w_f,                  &
                     eccf,         coszrs,       solin,        sfac,                         &
#ifdef MAML
                     cam_in%asdir(:ncol,ii), cam_in%asdif(:ncol,ii), cam_in%aldir(:ncol,ii), cam_in%aldif(:ncol,ii), &
#else
                     cam_in%asdir, cam_in%asdif, cam_in%aldir, cam_in%aldif,                 &
#endif
                     qrs,          qrsc,         fsnt,         fsntc,        fsntoa, fsutoa, &
                     fsntoac,      fsnirt,       fsnrtc,       fsnirtsq,     fsns,           &
#ifdef MAML
                     fsnsc,        fsdsc,        fsds,         sols_loc, soll_loc,   &
                     solsd_loc, solld_loc ,fns,          fcns,                               &
#else
                     fsnsc,        fsdsc,        fsds,         cam_out%sols, cam_out%soll,   &
                     cam_out%solsd,cam_out%solld,fns,          fcns,                         &
#endif

                     Nday,         Nnite,        IdxDay,       IdxNite,      clm_seed,       &
                     su,           sd,                                                       &
                     E_cld_tau=c_cld_tau, E_cld_tau_w=c_cld_tau_w, E_cld_tau_w_g=c_cld_tau_w_g, E_cld_tau_w_f=c_cld_tau_w_f, &
                     old_convert = .false.)
                call t_stopf ('rad_rrtmg_sw')
                   
                ! Output net fluxes at 200 mb
                call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcns, fsn200c)
                call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fns, fsn200)

                ! Calculate diagnostic quantities
                do i=1,ncol
                  swcf(i)=fsntoa(i) - fsntoac(i)
                  fsutoac(i) = solin(i) - fsntoac(i)
                end do

                ! Aggregate grid-mean averages from cloud-scale fluxes and heating rates 
                if (use_SPCAM) then 
                  do i = 1,ncol
                        qrs_m     (i,:pver, icall) =  qrs_m(i,:pver, icall) +  qrs(i,:pver)*factor_xy
                        qrsc_m    (i,:pver, icall) = qrsc_m(i,:pver, icall) + qrsc(i,:pver)*factor_xy
                        solin_m   (i, icall) = solin_m   (i, icall)+solin   (i)*factor_xy
                        fsds_m    (i, icall) = fsds_m    (i, icall)+fsds    (i)*factor_xy
                        fsnirt_m  (i, icall) = fsnirt_m  (i, icall)+fsnirt  (i)*factor_xy
                        fsnrtc_m  (i, icall) = fsnrtc_m  (i, icall)+fsnrtc  (i)*factor_xy
                        fsnirtsq_m(i, icall) = fsnirtsq_m(i, icall)+fsnirtsq(i)*factor_xy
                        fsnt_m    (i, icall) = fsnt_m    (i, icall)+fsnt    (i)*factor_xy
                        fsns_m    (i, icall) = fsns_m    (i, icall)+fsns    (i)*factor_xy
                        fsntc_m   (i, icall) = fsntc_m   (i, icall)+fsntc   (i)*factor_xy
                        fsnsc_m   (i, icall) = fsnsc_m   (i, icall)+fsnsc   (i)*factor_xy
                        fsdsc_m   (i, icall) = fsdsc_m   (i, icall)+fsdsc   (i)*factor_xy
                        fsntoa_m  (i, icall) = fsntoa_m  (i, icall)+fsntoa  (i)*factor_xy
                        fsutoa_m  (i, icall) = fsutoa_m  (i, icall)+fsutoa  (i)*factor_xy
                        fsntoac_m (i, icall) = fsntoac_m (i, icall)+fsntoac (i)*factor_xy
#ifdef MAML
                        sols_m    (i, icall) = sols_m    (i, icall)+sols_loc(i)*factor_xy
                        soll_m    (i, icall) = soll_m    (i, icall)+soll_loc(i)*factor_xy
                        solsd_m   (i, icall) = solsd_m   (i, icall)+solsd_loc(i)*factor_xy
                        solld_m   (i, icall) = solld_m   (i, icall)+solld_loc(i)*factor_xy
#else
                        sols_m    (i, icall) = sols_m    (i, icall)+cam_out%sols (i)*factor_xy
                        soll_m    (i, icall) = soll_m    (i, icall)+cam_out%soll (i)*factor_xy
                        solsd_m   (i, icall) = solsd_m   (i, icall)+cam_out%solsd(i)*factor_xy
                        solld_m   (i, icall) = solld_m   (i, icall)+cam_out%solld(i)*factor_xy
#endif
                        fsn200_m  (i, icall) = fsn200_m  (i, icall)+fsn200 (i)*factor_xy
                        fsn200c_m (i, icall) = fsn200c_m (i, icall)+fsn200c(i)*factor_xy
                        if (spectralflux) then
                           su_m(i,:,:,icall) = su_m(i,:,:,icall) + su(i,:,:)*factor_xy
                           sd_m(i,:,:,icall) = sd_m(i,:,:,icall) + sd(i,:,:)*factor_xy
                        end if
                  end do  ! i = 1,ncol

                  if(icall.eq.0) then  ! for the climate call
                       do i=1, ncol
                         crm_fsnt  (i,ii,jj) = fsnt(i)
                         crm_fsntc (i,ii,jj) = fsntc(i)
                         crm_fsns  (i,ii,jj) = fsns(i)
                         crm_fsnsc (i,ii,jj) = fsnsc(i)
#ifdef MAML
                         cam_out%netsw(i,ii) = fsns(i)
                         cam_out%sols(i,ii) = sols_loc(i)
                         cam_out%soll(i,ii) = soll_loc(i)
                         cam_out%solsd(i,ii) = solsd_loc(i)
                         cam_out%solld(i,ii) = solld_loc(i)
#endif
                         crm_aodvis(i,ii,jj) = sum(aer_tau(i, :, idx_sw_diag))
                         crm_aod400(i,ii,jj) = sum(aer_tau(i, :, idx_sw_diag+1))
                         crm_aod700(i,ii,jj) = sum(aer_tau(i, :, idx_sw_diag-1))
                         aod400(i) = aod400(i)+crm_aod400(i,ii,jj) * factor_xy
                         aod700(i) = aod700(i)+crm_aod700(i,ii,jj) * factor_xy
                       end do
                       do m=1,crm_nz
                         k = pver-m+1
                         qrs_crm(:ncol,ii,jj,m) = qrs(:ncol,k) / cpair
                         qrsc_crm(:ncol,ii,jj,m) = qrsc(:ncol,k) / cpair
                         crm_aodvisz(:ncol, ii, jj, m) = aer_tau(:ncol,k,idx_sw_diag)
                       end do

                       do i=1, ncol
                       do k=1, pver
                          if(c_cld_tau(idx_sw_diag,i,k).gt.1.0e-10) then
                            tot_icld_vistau(i,k)  = tot_icld_vistau(i,k) + c_cld_tau(idx_sw_diag,i,k)  
                            nct_tot_icld_vistau(i,k) = nct_tot_icld_vistau(i,k) + 1
                          end if
                          if(liq_tau(idx_sw_diag,i,k).gt.1.0e-10) then
                             liq_icld_vistau(i,k)  = liq_icld_vistau(i,k) + liq_tau(idx_sw_diag,i,k)
                             nct_liq_icld_vistau(i,k)  =   nct_liq_icld_vistau(i,k) + 1
                          end if
                          if(ice_tau(idx_sw_diag,i,k).gt.1.0e-10) then
                             ice_icld_vistau(i,k)  = ice_icld_vistau(i,k) + ice_tau(idx_sw_diag,i,k)
                             nct_ice_icld_vistau(i,k) = nct_ice_icld_vistau(i,k) + 1
                          end if
                          if(snow_tau(idx_sw_diag,i,k).gt.1.0e-10) then
                             snow_icld_vistau(i,k) = snow_icld_vistau(i,k) + snow_tau(idx_sw_diag,i,k)
                             nct_snow_icld_vistau(i,k) = nct_snow_icld_vistau(i,k) + 1
                          end if
                       end do
                       end do
                  end if  ! for the climate call

                  if(last_column) then
                      do i=1, ncol
                        qrs(i,:pver) = qrs_m(i,:pver, icall) 
                        qrsc(i,:pver) = qrsc_m(i,:pver, icall) 
                        solin(i) = solin_m(i, icall)
                        fsds(i)  = fsds_m(i, icall)
                        fsnirt(i)= fsnirt_m(i, icall)
                        fsnrtc(i)= fsnrtc_m(i, icall)
                        fsnirtsq(i)= fsnirtsq_m(i, icall)
                        fsnt(i)  = fsnt_m(i, icall)
                        fsns(i)  = fsns_m(i, icall)
                        fsntc(i) = fsntc_m(i, icall)
                        fsnsc(i) = fsnsc_m(i, icall)
                        fsdsc(i) = fsdsc_m(i, icall)
                        fsntoa(i)=fsntoa_m(i, icall)
                        fsutoa(i)=fsutoa_m(i, icall)
                        fsntoac(i)=fsntoac_m(i, icall)
#ifdef MAML
                        sols_loc(i)   =sols_m(i, icall)
                        soll_loc(i)   =soll_m(i, icall)
                        solsd_loc(i)   =solsd_m(i, icall)
                        solld_loc(i)   =solld_m(i, icall)
#else
                        cam_out%sols(i)   =sols_m(i, icall)
                        cam_out%soll(i)   =soll_m(i, icall)
                        cam_out%solsd(i)   =solsd_m(i, icall)
                        cam_out%solld(i)   =solld_m(i, icall)
#endif
                        fsn200(i)  = fsn200_m(i, icall)
                        fsn200c(i) = fsn200c_m(i, icall)
                        if (spectralflux) then
                           su(i,:,:) = su_m(i,:,:,icall)
                           sd(i,:,:) = sd_m(i,:,:,icall)
                        end if
                        swcf(i)=fsntoa(i) - fsntoac(i)
                        fsutoac(i) = solin(i) - fsntoac(i)
                      end do
                  end if  ! last_column
                end if ! (use_SPCAM)

                ! Dump shortwave radiation information to history tape buffer (diagnostics)
                if ( (use_SPCAM .and. last_column) .or. .not. use_SPCAM) then
                  ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
                  call outfld('QRS'//diag(icall),ftem  ,pcols,lchnk)
                  ftem(:ncol,:pver) = qrsc(:ncol,:pver)/cpair
                  call outfld('QRSC'//diag(icall),ftem  ,pcols,lchnk)
                  call outfld('SOLIN'//diag(icall),solin ,pcols,lchnk)
                  call outfld('FSDS'//diag(icall),fsds  ,pcols,lchnk)
                  call outfld('FSNIRTOA'//diag(icall),fsnirt,pcols,lchnk)
                  call outfld('FSNRTOAC'//diag(icall),fsnrtc,pcols,lchnk)
                  call outfld('FSNRTOAS'//diag(icall),fsnirtsq,pcols,lchnk)
                  call outfld('FSNT'//diag(icall),fsnt  ,pcols,lchnk)
                  call outfld('FSNS'//diag(icall),fsns  ,pcols,lchnk)
                  call outfld('FSNTC'//diag(icall),fsntc ,pcols,lchnk)
                  call outfld('FSNSC'//diag(icall),fsnsc ,pcols,lchnk)
                  call outfld('FSDSC'//diag(icall),fsdsc ,pcols,lchnk)
                  call outfld('FSNTOA'//diag(icall),fsntoa,pcols,lchnk)
                  call outfld('FSUTOA'//diag(icall),fsutoa,pcols,lchnk)
                  call outfld('FSNTOAC'//diag(icall),fsntoac,pcols,lchnk)
#ifdef MAML
                   call outfld('SOLS'//diag(icall),sols_loc  ,pcols,lchnk)
                   call outfld('SOLL'//diag(icall),soll_loc  ,pcols,lchnk)
                   call outfld('SOLSD'//diag(icall),solsd_loc ,pcols,lchnk)
                   call outfld('SOLLD'//diag(icall),solld_loc ,pcols,lchnk)
#else
                   call outfld('SOLS'//diag(icall),cam_out%sols  ,pcols,lchnk)
                   call outfld('SOLL'//diag(icall),cam_out%soll  ,pcols,lchnk)
                   call outfld('SOLSD'//diag(icall),cam_out%solsd ,pcols,lchnk)
                   call outfld('SOLLD'//diag(icall),cam_out%solld ,pcols,lchnk)
#endif
                  call outfld('FSN200'//diag(icall),fsn200,pcols,lchnk)
                  call outfld('FSN200C'//diag(icall),fsn200c,pcols,lchnk)
                  call outfld('SWCF'//diag(icall),swcf  ,pcols,lchnk)
                end if  ! (use_SPCAM .and. last_column) .or .not. use_SPCAM

                if(do_aerocom_ind3) then
                  aerindex = 0.0
                  angstrm = 0.0
                  aod400 = 0.0
                  aod700 = 0.0
                  do i=1, ncol
                     aod400(i) = sum(aer_tau(i, :, idx_sw_diag+1))
                     aod700(i) = sum(aer_tau(i, :, idx_sw_diag-1))
                     if(aod400(i).lt.1.0e4 .and. aod700(i).lt.1.e4  .and. &
                        aod400(i).gt.1.0e-10 .and. aod700(i).gt.1.0e-10) then
                        angstrm(i) = (log (aod400(i))-log(aod700(i)))/(log(0.700)-log(0.400))                               
                     else
                        angstrm(i) = fillvalue
                     end if
                     if(angstrm(i).ne.fillvalue) then 
                        aerindex(i) = angstrm(i)*sum(aer_tau(i,:,idx_sw_diag))
                     else 
                        aerindex(i) = fillvalue
                     end if
                  end do 
                  do i = 1, nnite
                     angstrm(idxnite(i)) = fillvalue
                     aod400(idxnite(i)) = fillvalue
                     aod700(idxnite(i)) = fillvalue
                     aerindex(idxnite(i)) = fillvalue
                  end do
                  if(icall.eq.0) then ! only for climatology run
                     call outfld('angstrm', angstrm, pcols, lchnk)
                     call outfld('aod400', aod400, pcols, lchnk)
                     call outfld('aod700', aod700, pcols, lchnk)
                     call outfld('aerindex', aerindex, pcols, lchnk)
                  end if
                end if  ! do_aerocom_ind3

              end if ! (active_calls(icall))
            end do ! icall

            if(use_SPCAM .and. last_column) then
              do i = 1, nnite 
                crm_aodvis(idxnite(i), :, :) = fillvalue
                crm_aod400(idxnite(i), :, :) = fillvalue
                crm_aod700(idxnite(i), :, :) = fillvalue
                aod400(idxnite(i)) = fillvalue
                aod700(idxnite(i)) = fillvalue
                crm_aodvisz(idxnite(i), :, :, :) = fillvalue
              end do
              call outfld('CRM_FSNT', crm_fsnt, pcols, lchnk)
              call outfld('CRM_FSNTC', crm_fsntc, pcols, lchnk)
              call outfld('CRM_FSNS', crm_fsns, pcols, lchnk)
              call outfld('CRM_FSNSC', crm_fsnsc, pcols, lchnk)
              call outfld('CRM_AODVIS', crm_aodvis, pcols, lchnk)
              call outfld('CRM_AOD400', crm_aod400, pcols, lchnk)
              call outfld('CRM_AOD700', crm_aod700, pcols, lchnk)
              call outfld('AOD400', aod400, pcols, lchnk)
              call outfld('AOD700', aod700, pcols, lchnk)
              call outfld('CRM_AODVISZ', crm_aodvisz, pcols, lchnk)

              do i=1,ncol
                do k=1,pver
                  tot_cld_vistau(i,k) = tot_icld_vistau(i,k) *  factor_xy
                  if(nct_tot_icld_vistau(i,k).ge.1) then
                    tot_icld_vistau(i,k)  = tot_icld_vistau(i,k)/nct_tot_icld_vistau(i,k)
                  else 
                    tot_icld_vistau(i,k)  = 0.0_r8
                  end if
                  if(nct_liq_icld_vistau(i,k).ge.1) then
                    liq_icld_vistau(i,k)  = liq_icld_vistau(i,k)/nct_liq_icld_vistau(i,k)
                  else
                    liq_icld_vistau(i,k)  = 0.0_r8
                  end if
                  if(nct_ice_icld_vistau(i,k).ge.1) then
                    ice_icld_vistau(i,k)  = ice_icld_vistau(i,k)/nct_ice_icld_vistau(i,k)
                  else
                    ice_icld_vistau(i,k)  = 0.0_r8
                  end if
                  if(nct_snow_icld_vistau(i,k).ge.1) then
                    snow_icld_vistau(i,k) = snow_icld_vistau(i,k)/nct_snow_icld_vistau(i,k)
                  else
                    snow_icld_vistau(i,k) = 0.0_r8
                  end if
                end do ! k
              end do ! i
            else if (.not. use_SPCAM) then
              ! Output cloud optical depth fields for the visible band
              tot_icld_vistau(:ncol,:)  = c_cld_tau(idx_sw_diag,:ncol,:)
              liq_icld_vistau(:ncol,:)  = liq_tau(idx_sw_diag,:ncol,:)
              ice_icld_vistau(:ncol,:)  = ice_tau(idx_sw_diag,:ncol,:)
              if (cldfsnow_idx > 0) then
                snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
              endif
              ! multiply by total cloud fraction to get gridbox value
              tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)
            endif  ! use_SPCAM .and. last_column

            ! add fillvalue for night columns
            if ( (use_SPCAM .and. last_column) .or. .not. use_SPCAM) then
              do i = 1, Nnite
                tot_cld_vistau(IdxNite(i),:)   = fillvalue
                tot_icld_vistau(IdxNite(i),:)  = fillvalue
                liq_icld_vistau(IdxNite(i),:)  = fillvalue
                ice_icld_vistau(IdxNite(i),:)  = fillvalue
                if (cldfsnow_idx > 0) then
                  snow_icld_vistau(IdxNite(i),:) = fillvalue
                endif
              end do

              call outfld('TOT_CLD_VISTAU', tot_cld_vistau, pcols, lchnk)       
              call outfld('TOT_ICLD_VISTAU', tot_icld_vistau, pcols, lchnk)
              call outfld('LIQ_ICLD_VISTAU', liq_icld_vistau, pcols, lchnk)
              call outfld('ICE_ICLD_VISTAU', ice_icld_vistau, pcols, lchnk)
              if (cldfsnow_idx > 0) then
                call outfld('SNOW_ICLD_VISTAU', snow_icld_vistau, pcols, lchnk)
              endif
            end if

            call t_stopf ('rad_sw')

          end if   ! dosw

          if( (use_SPCAM .and. last_column) .or. .not. use_SPCAM)  then
              ! Output aerosol mmr
              call rad_cnst_out(0, state, pbuf)
          end if  

          ! Longwave radiation computation
          if (dolw) then

            call t_startf ('rad_lw')

            ! Convert upward longwave flux units to CGS
           do i=1,ncol
#ifdef MAML
              lwupcgs(i) = cam_in%lwup(i,ii)*1000._r8
#else
              lwupcgs(i) = cam_in%lwup(i)*1000._r8
#endif
              if(single_column.and.scm_crm_mode.and.have_tg) &
                lwupcgs(i) = 1000*stebol*tground(1)**4
            end do

            ! Get the active climate/diagnostic shortwave calculations
            call rad_cnst_get_call_list(active_calls)

            ! Loop over diagnostic cases (each of which can contain different
            ! radiative constituents. The climate (icall==0) calculation must 
            ! occur last, so we loop from N_DIAG to 0.
            do icall = N_DIAG, 0, -1
              if (active_calls(icall)) then

                ! Update the concentrations in the RRTMG state object
                call rrtmg_state_update( state, pbuf, icall, r_state)

                ! Calculate aerosol optical properties
                call aer_rad_props_lw(is_cmip6_volc, icall, state, pbuf,  aer_lw_abs)
                    
                call t_startf ('rad_rrtmg_lw')
                call rad_rrtmg_lw( &
                         lchnk,        ncol,         num_rrtmg_levs,  r_state,                     &
                         state%pmid,   aer_lw_abs,   cldfprime,       c_cld_lw_abs,                &
                         qrl,          qrlc,                                                       &
#ifdef MAML
                         flns,         flnt,         flnsc,           flntc,        flwds_loc,     &
#else
                         flns,         flnt,         flnsc,           flntc,        cam_out%flwds, &
#endif
                         flut,         flutc,        fnl,             fcnl,         fldsc,         &
                         clm_seed,     lu,           ld                                            )
                call t_stopf ('rad_rrtmg_lw')

                if (lwrad_off) then
                  qrl(:,:) = 0._r8
                  qrlc(:,:) = 0._r8
                  flns(:) = 0._r8
                  flnt(:) = 0._r8
                  flnsc(:) = 0._r8
                  flntc(:) = 0._r8
#ifdef MAML
                  cam_out%flwds(:,ii) = 0._r8
#else
                  cam_out%flwds(:) = 0._r8
#endif
                  flut(:) = 0._r8
                  flutc(:) = 0._r8
                  fnl(:,:) = 0._r8
                  fcnl(:,:) = 0._r8
                  fldsc(:) = 0._r8
                end if !lwrad_off
    
                do i=1,ncol
                  lwcf(i)=flutc(i) - flut(i)
                end do

                !  Output fluxes at 200 mb
                call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fnl,  fln200)
                call vertinterp(ncol, pcols, pverp, state%pint, 20000._r8, fcnl, fln200c)

                ! Aggregate grid-mean averages from cloud-scale fluxes and heating rates 
                if (use_SPCAM) then
                  do i=1, ncol
                    qrl_m (i,:pver, icall) = qrl_m (i,:pver, icall) + qrl (i,:pver)*factor_xy
                    qrlc_m(i,:pver, icall) = qrlc_m(i,:pver, icall) + qrlc(i,:pver)*factor_xy
                    flnt_m   (i, icall) = flnt_m   (i, icall)+flnt(i)          *factor_xy
                    flut_m   (i, icall) = flut_m   (i, icall)+flut(i)          *factor_xy
                    flutc_m  (i, icall) = flutc_m  (i, icall)+flutc(i)         *factor_xy
                    flntc_m  (i, icall) = flntc_m  (i, icall)+flntc(i)         *factor_xy
                    flns_m   (i, icall) = flns_m   (i, icall)+flns(i)          *factor_xy
                    flnsc_m  (i, icall) = flnsc_m  (i, icall)+flnsc(i)         *factor_xy
                    fldsc_m  (i, icall) = fldsc_m  (i, icall)+fldsc(i)         *factor_xy
#ifdef MAML
                    flwds_m  (i, icall) = flwds_m  (i, icall)+flwds_loc(i) *factor_xy
#else
                    flwds_m  (i, icall) = flwds_m  (i, icall)+cam_out%flwds(i) *factor_xy
#endif
                    fln200_m (i, icall) = fln200_m (i, icall)+fln200(i)        *factor_xy
                    fln200c_m(i, icall) = fln200c_m(i, icall)+fln200c(i)       *factor_xy
                    if (spectralflux) then
                       lu_m(i,:,:,icall) = lu_m(i,:,:,icall) + lu(i,:,:)*factor_xy
                       ld_m(i,:,:,icall) = ld_m(i,:,:,icall) + ld(i,:,:)*factor_xy
                    end if

                    ! Only save the CRM fluxes for the case that affects the
                    ! climate (icall == 0)
                    if(icall.eq.0) then
                      crm_flnt (i,ii,jj) = flnt(i)
                      crm_flntc(i,ii,jj) = flntc(i)
                      crm_flns (i,ii,jj) = flns(i)
                      crm_flnsc(i,ii,jj) = flnsc(i)
#ifdef MAML
                      cam_out%flwds(i,ii) = flwds_loc(i)
#endif
                      do m=1,crm_nz
                         k = pver-m+1
                         qrl_crm(:ncol,ii,jj,m) = qrl(:ncol,k) / cpair
                         qrlc_crm(:ncol,ii,jj,m) = qrlc(:ncol,k) / cpair
                      end do
                    end if  ! icall == 0

                  end do  ! i = 1,ncol

                  ! Set GCM fluxes to grid-means of the cloud-scale fluxes if this
                  ! is the last column (since the aggregated averages will
                  ! represent the full grid-mean by now)
                  if(last_column) then
                    do i = 1,ncol
                      qrl (i,:pver) = qrl_m(i,:pver, icall)
                      qrlc(i,:pver) = qrlc_m(i,:pver, icall)
                      flnt (i) = flnt_m (i, icall)
                      flut (i) = flut_m (i, icall)
                      flutc(i) = flutc_m(i, icall)
                      flntc(i) = flntc_m(i, icall)
                      flns (i) = flns_m (i, icall)
                      flnsc(i) = flnsc_m(i, icall)
                      fldsc(i) = fldsc_m(i, icall)
#ifdef MAML
                      flwds_loc(i) = flwds_m(i, icall)
#else
                      cam_out%flwds(i) = flwds_m(i, icall)
#endif
                      fln200 (i) = fln200_m (i, icall)
                      fln200c(i) = fln200c_m(i, icall)
                      if (spectralflux) then
                         lu(i,:,:) = lu_m(i,:,:,icall)
                         ld(i,:,:) = ld_m(i,:,:,icall)
                      end if
                      lwcf(i)=flutc(i) - flut(i) 
                    end do  ! i = 1,ncol
                  endif  ! last_column

                endif ! use_SPACM

                ! Dump longwave radiation information to history tape buffer (diagnostics)
                if ( (use_SPCAM .and. last_column ) .or. .not. use_SPCAM) then
                  call outfld('QRL'//diag(icall),qrl (:ncol,:)/cpair,ncol,lchnk)
                  call outfld('QRLC'//diag(icall),qrlc(:ncol,:)/cpair,ncol,lchnk)
                  call outfld('FLNT'//diag(icall),flnt  ,pcols,lchnk)
                  call outfld('FLUT'//diag(icall),flut  ,pcols,lchnk)
                  call outfld('FLUTC'//diag(icall),flutc ,pcols,lchnk)
                  call outfld('FLNTC'//diag(icall),flntc ,pcols,lchnk)
                  call outfld('FLNS'//diag(icall),flns  ,pcols,lchnk)
              
                  call outfld('FLDSC'//diag(icall),fldsc ,pcols,lchnk)
                  call outfld('FLNSC'//diag(icall),flnsc ,pcols,lchnk)
                  call outfld('LWCF'//diag(icall),lwcf  ,pcols,lchnk)
                  call outfld('FLN200'//diag(icall),fln200,pcols,lchnk)
                  call outfld('FLN200C'//diag(icall),fln200c,pcols,lchnk)
#ifdef MAML
                  call outfld('FLDS'//diag(icall),flwds_loc ,pcols,lchnk)
#else
                  call outfld('FLDS'//diag(icall),cam_out%flwds ,pcols,lchnk)
#endif
                end if
                if (use_SPCAM .and. last_column ) then
                  if(icall.eq.0) then  ! the climate call
                    call outfld('CRM_FLNT', crm_flnt, pcols, lchnk)
                    call outfld('CRM_FLNTC', crm_flntc, pcols, lchnk)
                    call outfld('CRM_FLNS', crm_flns, pcols, lchnk)
                    call outfld('CRM_FLNSC', crm_flnsc, pcols, lchnk)
                  end if   ! the climate call
                end if

              end if  ! active_calls(icall)
            end do ! icall

            call t_stopf ('rad_lw')

          end if  !dolw
        end do ! ii = 1,crm_nx_rad
      end do ! jj = 1,crm_nx_rad

      ! Restore pbuf and state to values as input to this routine before we
      ! modified them in-place to populate with CRM column values
      if (use_SPCAM) then 
        cld    = cld_save
        cicewp = cicewp_save
        cliqwp = cliqwp_save
        if (cldfsnow_idx > 0) then
          csnowp   = csnowp_save
          cldfsnow = cldfsnow_save   
        end if
        rel   = rel_save
        rei   = rei_save
        state = statein_copy
        dei   = dei_save
        deallocate(dei_save)
        if (SPCAM_microp_scheme .eq. 'm2005') then
           mu      = mu_save
           lambdac = lambdac_save
           des     = des_save
           deallocate (mu_save, lambdac_save, des_save)
        endif
#ifdef MODAL_AERO
        qaerwat  = qaerwat_save
        dgnumwet = dgnumwet_save
#endif
      endif ! use_SPCAM

      ! Calculate net CRM heating rate from shortwave and longwave heating rates
      if (use_SPCAM) then 
        do m = 1,crm_nz
          do i = 1,ncol
            crm_qrad(i,:,:,m) = (qrs_crm(i,:,:,m) + qrl_crm(i,:,:,m))
          end do
        end do
      endif ! use_SPCAM

      ! deconstruct the RRTMG state object
      call rrtmg_state_destroy(r_state)

      ! mji/hirsrtm - Add call to HIRSRTM package
      ! HIRS brightness temperature calculation in 7 infra-red channels and 4 microwave
      ! channels as a diagnostic to compare to TOV/MSU satellite data.
      ! Done if dohirs set to .true. at time step frequency ihirsfq
      nstep = get_nstep()
      if ( dohirs .and. (mod(nstep-1,ihirsfq) .eq. 0) ) then

          do i= 1, ncol
#ifdef MAML
             lwup_loc =0._r8
             do ii =1, crm_nx
                lwup_loc = lwup_loc + cam_in%lwup(i,ii)
             end do
             lwup_loc = lwup_loc/crm_nx
             ts(i) = sqrt(sqrt(lwup_loc/stebol))
#else
             ts(i) = sqrt(sqrt(cam_in%lwup(i)/stebol))
#endif
             ! Set oro (land/sea flag) for compatibility with landfrac/icefrac/ocnfrac
             ! oro=0 (sea or ice); oro=1 (land)
             if (landfrac(i).ge.0.001) then
                oro(i)=1.
             else
                oro(i)=0.
             endif
             ! Convert pressure from Pa to hPa
             do k = 1, pver
                pintmb(i,k) = state%pint(i,k)*1.e-2_r8        
             end do
             pintmb(i,pverp) = state%pint(i,pverp)*1.e-2_r8 
          end do
          
          ! Get constituent mixing ratios (specific humidity, ozone mass mixing
          ! ratio, CO2 mass mixing ratio
          call rad_cnst_get_gas(0,'H2O', state, pbuf, sp_hum)
          call rad_cnst_get_gas(0,'O3',  state, pbuf, o3)
          call rad_cnst_get_gas(0,'CO2', state, pbuf, co2)

          call calc_col_mean(state, co2, co2_col_mean)

          ! Call the hirsrtm driver
          call t_startf ('hirstrm')
          call hirsrtm( lchnk  ,ncol , &
                        pintmb ,state%t  ,sp_hum ,co2_col_mean, &
                        o3     ,ts       ,oro    ,tb_ir  ,britemp )
          call t_stopf ('hirstrm')

          ! Send outputs to history buffer
          do i = 1, pnb_hirs
             call outfld(hirsname(i),tb_ir(1,i),pcols,lchnk)
          end do
          do i = 1, pnf_msu
             call outfld(msuname(i),britemp(1,i),pcols,lchnk)
          end do

      end if

      ! Run the CFMIP Observation Simulator Package (COSP)
      ! For the time being, the MMF stuff is not coupled with the COSP
      ! simulator, so bypass this code if we are using SP/MMF (for now)
      if (.not. use_SPCAM) then 
          !! initialize and calculate emis
          emis(:,:) = 0._r8
          emis(:ncol,:) = 1._r8 - exp(-cld_lw_abs(rrtmg_lw_cloudsim_band,:ncol,:))
          call outfld('EMIS      ',emis    ,pcols   ,lchnk   )

          !! compute grid-box mean SW and LW snow optical depth for use by COSP
          gb_snow_tau(:,:) = 0._r8
          gb_snow_lw(:,:) = 0._r8
          if (cldfsnow_idx > 0) then
             do i=1,ncol
                do k=1,pver
                   if(cldfsnow(i,k) > 0.)then
                      gb_snow_tau(i,k) = snow_tau(rrtmg_sw_cloudsim_band,i,k)*cldfsnow(i,k)
                      gb_snow_lw(i,k) = snow_lw_abs(rrtmg_lw_cloudsim_band,i,k)*cldfsnow(i,k)
                   end if
                enddo
             enddo
          end if

          if (docosp) then
             !! cosp_cnt referenced for each chunk... cosp_cnt(lchnk)
             !! advance counter for this timestep
             cosp_cnt(lchnk) = cosp_cnt(lchnk) + 1

             !! if counter is the same as cosp_nradsteps, run cosp and reset counter
              if (cosp_nradsteps .eq. cosp_cnt(lchnk)) then
                 !call should be compatible with camrt radiation.F90 interface too, should be with (in),optional
                 ! N.B.: For snow optical properties, the GRID-BOX MEAN shortwave and longwave optical depths are passed.

                 call t_startf ('cosp_run')
                 call cospsimulator_intr_run(state,  pbuf, cam_in, emis, coszrs, &
                      cld_swtau_in=cld_tau(rrtmg_sw_cloudsim_band,:,:),&
                      snow_tau_in=gb_snow_tau,snow_emis_in=gb_snow_lw)
                 cosp_cnt(lchnk) = 0  !! reset counter
                 call t_stopf ('cosp_run')

              end if
          end if
      endif  ! use_SPCAM

      if (use_SPCAM .and. SPCAM_microp_scheme .eq. 'm2005') then
          call outfld('CRM_MU    ', mu_crm     , pcols, lchnk)
          call outfld('CRM_DES   ', des_crm    , pcols, lchnk)
          call outfld('CRM_LAMBDA', lambdac_crm, pcols, lchnk)
          call outfld('CRM_TAU   ', cld_tau_crm, pcols, lchnk)
          deallocate(des_crm, mu_crm, lambdac_crm)
      endif

      if (use_SPCAM) then 
          call outfld('CRM_DEI  ', dei_crm, pcols, lchnk)
          deallocate(dei_crm)
          call outfld('CRM_REL  ', rel_crm, pcols, lchnk)
          call outfld('CRM_REI  ', rei_crm, pcols, lchnk)
          call outfld('CRM_QRL  ', qrl_crm, pcols, lchnk)
          call outfld('CRM_QRS  ', qrs_crm, pcols, lchnk)
          call outfld('CRM_QRLC ', qrlc_crm, pcols, lchnk)
          call outfld('CRM_QRSC ', qrsc_crm, pcols, lchnk)
      endif

    else  !  if (dosw .or. dolw) then

       ! If conserve_energy is true, then heating rates are multiplied by dp at
       ! the end of this routine to carry Q*dp across timesteps to conserve
       ! energy. Thus, if this flag is set, we need to divide by dp here before
       ! working with Q below because it was multiplied by dp in a previous
       ! call.
       if (conserve_energy) then
          do k = 1,pver
             do i = 1, ncol
                qrs(i,k) = qrs(i,k)/state%pdel(i,k)
                qrl(i,k) = qrl(i,k)/state%pdel(i,k)
             end do
          end do
          if (use_SPCAM) then
             do m = 1,crm_nz
                k = pver - m + 1
                do i = 1,ncol
                  crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) / state%pdel(i,k)
                end do
             end do
          end if
       end if

    end if  !  if (dosw .or. dolw)

    ! Compute net radiative heating tendency
    call t_startf ('radheat_tend')
    call radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                      fsnt, flns, flnt, cam_in%asdir, net_flx)
    call t_stopf ('radheat_tend')

    ! Compute heating rate for dtheta/dt
    do k=1,pver
       do i=1,ncol
          ftem(i,k) = (qrs(i,k) + qrl(i,k)) / cpair * (1.e5_r8/state%pmid(i,k))**cappa
       end do
    end do
    call outfld('HR', ftem, pcols, lchnk)

    ! convert radiative heating rates to Q*dp for energy conservation
    if (conserve_energy) then
       do k =1 , pver
          do i = 1, ncol
             qrs(i,k) = qrs(i,k)*state%pdel(i,k)
             qrl(i,k) = qrl(i,k)*state%pdel(i,k)
          end do
       end do
       if (use_SPCAM) then
          do m = 1,crm_nz
             k = pver - m + 1
             do i = 1,ncol
               crm_qrad(i,:,:,m) = crm_qrad(i,:,:,m) * state%pdel(i,k)
             end do
          end do
       end if
    end if

    ! write kissvec seeds for random numbers
    if (pergro_mods) then
       do iseed = 1, kiss_seed_num    
          do i = 1, ncol          
             rad_randn_seedrst(i,iseed,lchnk) = clm_seed(i,iseed)
          enddo
       enddo
    endif

    ! Compute net surface radiative flux for use by surface temperature code.
    ! Note that units have already been converted to mks in RADCTL.  Since
    ! fsns and flwds are in the buffer, array values will be carried across
    ! timesteps when the radiation code is not invoked.
    !cam_out%srfrad(:ncol) = fsns(:ncol) + cam_out%flwds(:ncol)
    !call outfld('SRFRAD  ',cam_out%srfrad,pcols,lchnk)
#ifndef MAML
     cam_out%netsw(:ncol) = fsns(:ncol)
#endif
  
  end subroutine radiation_tend

!===============================================================================

subroutine radinp(ncol, pmid, pint, pmidrd, pintrd, eccf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
! Convert model pressures to cgs.
! 
! Author: CCM1, CMS Contact J. Kiehl
!-----------------------------------------------------------------------
   use shr_orb_mod
   use time_manager, only: get_curr_calday

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! Pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! Pressure at model interfaces (pascals)
!
! Output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! Pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! Pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! Earth-sun distance factor

!
!---------------------------Local variables-----------------------------
!
   integer i                ! Longitude loop index
   integer k                ! Vertical loop index

   real(r8) :: calday       ! current calendar day
   real(r8) :: delta        ! Solar declination angle
!-----------------------------------------------------------------------
!
   calday = get_curr_calday()
   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf)

!
! Convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0_r8
         pintrd(i,k) = pint(i,k)*10.0_r8
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0_r8
   end do

end subroutine radinp

!===============================================================================

subroutine calc_col_mean(state, mmr_pointer, mean_value)
!----------------------------------------------------------------------- 
! 
! Compute the column mean mass mixing ratio.  
!
!-----------------------------------------------------------------------

   use cam_logfile,  only: iulog

   type(physics_state),        intent(in)  :: state
   real(r8), dimension(:,:),   pointer     :: mmr_pointer  ! mass mixing ratio (lev)
   real(r8), dimension(pcols), intent(out) :: mean_value   ! column mean mmr

   integer  :: i, k, ncol
   real(r8) :: ptot(pcols)
   !-----------------------------------------------------------------------

   ncol         = state%ncol
   mean_value   = 0.0_r8
   ptot         = 0.0_r8

   do k=1,pver
      do i=1,ncol
         mean_value(i) = mean_value(i) + mmr_pointer(i,k)*state%pdeldry(i,k)
         ptot(i)         = ptot(i) + state%pdeldry(i,k)
      end do
   end do
   do i=1,ncol
      mean_value(i) = mean_value(i) / ptot(i)
   end do

end subroutine calc_col_mean

!===============================================================================

end module radiation

