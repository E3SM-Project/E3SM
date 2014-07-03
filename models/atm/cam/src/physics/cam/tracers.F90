!======================================================================
! This is an interface for 3 test tracers: tr1,tr2,tr3
!
! This uses the tracers_suite module to initialize 
!   mixing ratios, fluxes and calculate tendencies.
!   Details of calling tree below. All of the detailed information about the
!  tracers should be store in the suite file, including the number & names of tracers. 
!
! Author B. Eaton
! History  D. Bundy, June 2003 modified to the format of physics interface
!        
!
!---------------------------------------------------------------
!
!  ------------  calling tree --------------
!  Register the tracers as advected fields, pass names to model
!  initindx.F90:			call tracers_register()
!
!  Initialize the tracer mixing ratio field
!  inidat.F90:read_inidat
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  inti.F90
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!  		-> addfld/add default for surface flux (SF)
!
!  Timestepping:
!  advnce.F90
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  tphysac.F90
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!  		-> tracers_suite.F90:tend_tr
!
!======================================================================


module tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog

  implicit none
  private
  save

! Public interfaces
  public tracers_register                  ! register constituent
  public tracers_implements_cnst           ! true if named constituent is implemented by this package
  public tracers_init_cnst                 ! initialize constituent field
  public tracers_init                      ! initialize history fields, datasets
  public tracers_timestep_tend             ! calculate tendencies
  public tracers_timestep_init             ! interpolate dataset for constituent each timestep

! Data from namelist variables
  logical, public :: tracers_flag  = .false.     ! true => turn on test tracer code, namelist variable

! Private module data

  integer :: trac_ncnst                    ! total number of test tracers
  integer :: ixtrct=-999                   ! index of 1st constituent
  logical :: debug = .false.
  
contains
!======================================================================
subroutine tracers_register
!----------------------------------------------------------------------- 
!
! Purpose: register advected tracers. Called by initindx.F90
!  The registration lets the model know what the tracer names are
!  and returns the index number ixtrct for the constituent array
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

   use physconst,    only: mwdry, cpair
   use constituents, only: cnst_add, cnst_num_avail
   use tracers_suite, only: get_tracer_name
   
   implicit none
!---------------------------Local workspace-----------------------------
   integer :: mm,m                                 ! dummy
   character(len=8) :: name   ! constituent name
   real(r8) minc

!-----------------------------------------------------------------------
   if ( tracers_flag ) then 
      minc = 0        ! min mixing ratio (normal setting)
      minc = -1.e36_r8   ! min mixing ratio (disable qneg3)
      
      ! Set the number of test tracers equal to the number of slots available
      ! in the constituent array
      trac_ncnst = cnst_num_avail()

      do m = 1,trac_ncnst 
         name = get_tracer_name(m)  ! get name from suite file
         
         ! add constituent name to list of advected, save index number ixtrct
         call cnst_add(name, mwdry, cpair, minc, mm, &  
              readiv=.false.,mixtype='dry')
         if ( m .eq. 1 ) ixtrct = mm  ! save index number of first tracer
         
      end do
   end if

end subroutine tracers_register
!======================================================================

function tracers_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

  use tracers_suite, only: get_tracer_name
  
  implicit none
!-----------------------------Arguments---------------------------------
  
  character(len=*), intent(in) :: name   ! constituent name
  logical :: tracers_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
   integer :: m
!-----------------------------------------------------------------------

   tracers_implements_cnst = .false.
   if ( tracers_flag ) then 
      do m = 1, trac_ncnst
         if (name == get_tracer_name(m)) then
            tracers_implements_cnst = .true.
            return
         end if
      end do
   end if
end function tracers_implements_cnst

!===============================================================================
subroutine tracers_init_cnst(name, q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose: initialize test tracers mixing ratio fields 
!  This subroutine is called at the beginning of an initial run ONLY
!
!-----------------------------------------------------------------------

  use tracers_suite,   only: init_cnst_tr, get_tracer_name

  implicit none

  character(len=*), intent(in) :: name
  real(r8), intent(out), dimension(:,:) :: q    ! kg tracer/kg dry air (gcol,plev)
  integer,  intent(in)                  :: gcid(:)  ! global column id
! Local
  integer m
  if ( tracers_flag ) then 
     do m = 1, trac_ncnst
        if (name ==  get_tracer_name(m))  then
           call init_cnst_tr(m,q, gcid)
        endif
     end do
  end if

end subroutine tracers_init_cnst

!===============================================================================
subroutine tracers_init

!----------------------------------------------------------------------- 
!
! Purpose: declare history variables, initialize data sets
!  This subroutine is called at the beginning of an initial or restart run
!
!-----------------------------------------------------------------------

   use tracers_suite,   only: init_tr, get_tracer_name
   use cam_history,     only: addfld, add_default, phys_decomp
   use ppgrid,          only: pver
   use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam

   ! Local
   integer m, mm
   character(len=8) :: name   ! constituent name

   if ( tracers_flag ) then     
     
      do m = 1,trac_ncnst 
         name = get_tracer_name(m)
         call cnst_get_ind(name, mm)
         call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
         call addfld (sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)

         call add_default (cnst_name(mm), 1, ' ')
         call add_default (sflxnam(mm),   1, ' ')
      end do
     
      ! initialize datasets, etc, needed for constituents.
      call init_tr  

  endif
     
end subroutine tracers_init

!======================================================================

subroutine tracers_timestep_init( phys_state )
!----------------------------------------------------------------------- 
!
! Purpose: At the beginning of a timestep, there are some things to do
! that just the masterproc should do. This currently just interpolates
! the emissions boundary data set to the current time step.
!
!-----------------------------------------------------------------------

  use tracers_suite, only: timestep_init_tr

  ! phys_state argument is unused in this version
  use ppgrid,         only: begchunk, endchunk
  use physics_types,  only: physics_state
  type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    
!-----------------------------------------------------------------------

  if ( tracers_flag ) then 
     
     call timestep_init_tr
     
     if (debug) write(iulog,*)'tracers_timestep_init done'
  endif

end subroutine tracers_timestep_init

!======================================================================

subroutine tracers_timestep_tend(state, ptend, cflx, landfrac, deltat)

!----------------------------------------------------------------------- 
!
! Purpose: During the timestep, compute test tracer mixing ratio 
! tendencies and surface fluxes.
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,        only: pcols, pver
  use constituents,  only: pcnst, sflxnam, cnst_cam_outfld
  use tracers_suite, only: flux_tr, tend_tr
  use cam_history,   only: outfld

  implicit none

  ! Arguments
   type(physics_state), intent(in)  :: state          ! state variables
   type(physics_ptend), intent(out) :: ptend          ! package tendencies
   real(r8),            intent(in)  :: deltat         ! timestep
   real(r8),            intent(in)  :: landfrac(pcols) ! Land fraction
   real(r8), intent(inout) :: cflx(pcols,pcnst) ! Surface constituent flux (kg/m^2/s)

! Local variables
   integer  :: m               ! tracer number (internal)

   logical  :: lq(pcnst)

!-----------------------------------------------------------------------

  if (.not. tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
  else
     lq(:)      = .FALSE.
     lq(ixtrct:ixtrct+trac_ncnst-1) = .TRUE.
     call physics_ptend_init(ptend, state%psetcols, 'tracers', lq=lq)
 
     do  m = 1,trac_ncnst 
        if (debug) write(iulog,*)'tracers.F90 calling for tracer ',m
        
        !calculate flux
        call flux_tr(m,state%ncol,state%lchnk, landfrac, cflx(:,ixtrct+m-1))
        
        !calculate tendency
        call tend_tr(m,state%ncol, state%q(:,:,ixtrct+m-1), deltat, ptend%q(:,:,ixtrct+m-1))
        
        !outfld calls could go here
        if ( cnst_cam_outfld(ixtrct+m-1) ) then
           call outfld (sflxnam(ixtrct+m-1),cflx(:,ixtrct+m-1),pcols,state%lchnk)
        end if
        
        
     end do
     
     if ( debug ) then 
        do  m = 1,trac_ncnst 
           write(iulog,*)'tracers_timestep_tend ixtrct,m,ixtrct+m-1',ixtrct,m,ixtrct+m-1
           write(iulog,*)'tracers_timestep_tend min max flux',minval(cflx(:,ixtrct+m-1)),maxval(cflx(:,ixtrct+m-1))
           write(iulog,*)'tracers_timestep_tend min max tend',minval(ptend%q(:,:,ixtrct+m-1)),maxval(ptend%q(:,:,ixtrct+m-1))
        end do
        write(iulog,*)'tracers_timestep_tend end'
     endif
  endif


end subroutine tracers_timestep_tend

end module tracers
