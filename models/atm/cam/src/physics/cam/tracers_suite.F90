
module tracers_suite

!---------------------------------------------------------------
!
! Implements artificial suite of tracers
!    1) low tracer with unit mixing ratio at about 800 mb
!    2) med tracer with unit mixing ratio at about 500 mb
!    3) high tracer with unit mixing ratio at about 200 mb
!    4) reverse med tracer with unit mixing ratio everywhere except about 500 mb
!    5) unit tracer with unit mixing ratio everywhere
!
!  D Bundy June 2003
!  modified Feb 2004 to include TT_UN and smoothing
!
!  A Mirin and B Eaton, August 2007
!  Modified to create up to 1000 distinct copies of the 5 basic tracers
!  by appending up to a 3 digit number to the base tracer name.
!  RESTRICTION - trac_ncnstmx cannot exceed 5000 unless the algorithm for
!                constructing new tracer names is extended.
!
! This is a module that contains a suite of tracers that is interfaced
! by tracers.F90. The details of the suite are contained entirely
! in this file, so the public routines are all very generic. 
!
!  ------------  calling tree --------------
!  Initialize the tracer mixing ratio field
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!  			-> init_cnst_lw
!  			-> init_cnst_md
!  			-> init_cnst_hi
!  			-> init_cnst_un
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!
!  Timestepping:
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!  			-> flux_lw  
!  			-> flux_md
!  			-> flux_hi
!  			-> flux_un
!
!  		-> tracers_suite.F90:tend_tr
!  			-> tend_lw
!  			-> tend_md
!  			-> tend_hi
!  			-> tend_un
!
!---------------------------------------------------------------


  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog

  implicit none
  private
  save

! Public interfaces
  public get_tracer_name  ! store names of tracers
  public init_cnst_tr      ! initialize tracer fields
  public init_tr      ! initialize data sets need for tracer
  public tend_tr      ! tracer tendency
  public flux_tr      ! surface flux of tracer
  public timestep_init_tr  ! interpolate tracer emissions data set


! Private module data
  integer, parameter :: trac_ncnstmx=5000 ! Max no. of tracers based on current algorithm for
                                          ! constructing tracer names.  This could easily be extended.
  integer, parameter :: trac_names=5      ! No. of base tracers

  logical, parameter :: smooth = .false.
  
contains

!======================================================================
function get_tracer_name(n)

!------------------------------------------------------------------------
! Purpose:
!
! The tracer names are only defined in this module. This function is for
! outside programs to grab the name for each tracer number. 
!    If n > trac_ncst calls endrun.
!
!------------------------------------------------------------------------

! -----------------------------Arguments---------------------------------

  integer, intent(in) :: n
  character(len=8) :: get_tracer_name  ! return value

! ----------------------------- Local ---------------------------------
  character(len=5), dimension(trac_names), parameter :: & ! constituent names
       tracer_names  =  (/ 'TT_LW', 'TT_MD', 'TT_HI', 'TTRMD' , 'TT_UN'/)
  
!-----------------------------------------------------------------------

  integer :: nbase  ! Corresponding base tracer index
  integer :: ncopy  ! No. of copies of base tracers
  character(len=1) :: c1
  character(len=2) :: c2
  character(len=3) :: c3

  if ( n > trac_ncnstmx ) then
     write(iulog,*) 'tracers_suite:get_tracer_name()','requested tracer',n
     write(iulog,*) 'only ',trac_ncnstmx,' tracers available'
     call endrun
  else
     nbase = mod(n-1, trac_names) + 1
     ncopy = (n-1)/trac_names
     if ( ncopy == 0 ) then
        get_tracer_name = tracer_names(nbase)
     else if ( ncopy >= 1  .and.  ncopy <= 9 ) then
        write (c1,'(i1)') ncopy
        get_tracer_name = tracer_names(nbase) // c1
     else if ( ncopy >= 10  .and.  ncopy <= 99 ) then
        write (c2,'(i2)') ncopy
        get_tracer_name = tracer_names(nbase) // c2
     else if ( ncopy >= 100  .and.  ncopy <= 999 ) then
        write (c3,'(i3)') ncopy
        get_tracer_name = tracer_names(nbase) // c3
     end if
  endif

  return

end function get_tracer_name


!======================================================================
!======================================================================
subroutine init_cnst_tr(m,q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose:
! calls initialization routine for tracer m, returns mixing ratio in q
!
! This routine must be consistent with trac_ncnstmx.
!
!----------------------------------------------------------------------- 

   implicit none

   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air (gcol, plev)
   integer ,intent(in) :: m           ! index of tracer
   integer, intent(in) :: gcid(:)     ! global column id   

   integer nbase ! Corresponding base tracer index

   if ( m > trac_ncnstmx ) then
      write(iulog,*) 'tracers_suite:init_cnst_tr()'
      write(iulog,*) ' asked to initialize tracer number ',m
      write(iulog,*) ' but there are only trac_ncnstmx = ',trac_ncnstmx,' tracers'
      call endrun
   endif

   nbase = mod(m-1,trac_names)+1

   if ( nbase == 1 ) then
      call init_cnst_lw(q, gcid)
   else if ( nbase == 2 ) then
      call init_cnst_md(q, gcid)
   else if ( nbase == 3 ) then
      call init_cnst_hi(q, gcid)
   else if ( nbase == 4 ) then
      call init_cnst_md(q, gcid, rev_in=1)
   else if ( nbase == 5 ) then
      call init_cnst_un(q, gcid)
   else
      write(iulog,*) 'tracers_suite:init_cnst_tr()'
      write(iulog,*) 'no initialization routine specified for tracer',nbase
      call endrun
   endif
      
end subroutine init_cnst_tr



!======================================================================
subroutine init_cnst_lw(q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 1: low
! 
!-----------------------------------------------------------------------
   implicit none

!Arguments
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air (gcol,plev)
   integer,  intent(in)  :: gcid(:)   ! global column id
! Local
  integer indx

!-----------------------------------------------------------------------
!
! Initialize low tracer to zero except at 800 level
!

   indx = setpindxtr(800._r8)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.876_r8)
  else 
     q = 0.0_r8
     q(:,indx) = 1.0_r8
  endif

!   write(iulog,*)'suite 1 init_cnst_lw min/max q',minval(q),maxval(q)

end subroutine init_cnst_lw

!======================================================================
subroutine init_cnst_md(q,gcid,rev_in)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 2: med
! 
!-----------------------------------------------------------------------
   implicit none

!Arguments
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air
   integer, intent(in) :: gcid(:)     ! global column id
   integer,  intent(in), optional :: rev_in         ! reverse the mixing ratio

! Local
  integer indx
  integer rev

!-----------------------------------------------------------------------
!
! Initialize med tracer to zero except at 500 level
!

  rev = 0
  if (present(rev_in)) then
     if (rev_in == 1) then
        rev = 1
     endif
  endif

  indx = setpindxtr(500._r8)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.876_r8,rev_in=rev)
  else
     if (rev == 1 ) then
        q = 1.0_r8
        q(:,indx) = 0.0_r8
     else
        q = 0.0_r8
        q(:,indx) = 1.0_r8
     endif
  endif
   
!   write(iulog,*)'suite 1 init_cnst_md min/max q',minval(q),maxval(q)

end subroutine init_cnst_md

!======================================================================
subroutine init_cnst_hi(q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 3: high
! 
!-----------------------------------------------------------------------

   implicit none

!Arguments
   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air
   integer, intent(in) :: gcid(:)     ! global column id
! Local
  integer indx

!-----------------------------------------------------------------------
!
! Initialize high tracer to zero except at 200 level
!

  indx = setpindxtr(200._r8)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.3_r8)
  else
     q = 0.0_r8
     q(:,indx) = 1.0_r8
  endif
  !   write(iulog,*)'suite 1 init_cnst_hi min/max q',minval(q),maxval(q)


end subroutine init_cnst_hi

!======================================================================
subroutine init_cnst_un(q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test unit tracer.
!    2) conserved unit tracer
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! copied from tracers.F90:initesttr()  D.Bundy Oct 15 2002
!-----------------------------------------------------------------------
   implicit none

   real(r8), intent(out) :: q(:,:)    ! kg tracer/kg dry air
   integer, intent(in)   :: gcid(:)   ! global column id
!-----------------------------------------------------------------------
! Initialize conserved unit tracer.

   q = 1.0_r8

end subroutine init_cnst_un


!======================================================================

subroutine init_tr

  !----------------------------------------------------------------------- 
  ! Purpose:
  !  Initialize any datasets that the constituents need. 
  !  Currently does nothing for this suite
  ! D Bundy, May 30, 2003
  !-----------------------------------------------------------------------


end subroutine init_tr

!======================================================================

subroutine timestep_init_tr
!----------------------------------------------------------------------- 
! 
! Purpose: call tracer timestep init processes
!----------------------------------------------------------------------- 


end subroutine timestep_init_tr

!======================================================================

subroutine setsmoothtr(indx,q,width,rev_in)
  use ref_pres, only : pref_mid
  implicit none


  !Arguments
  integer, intent(in)     :: indx               ! k index of pressure level
  real(r8), intent(inout) :: q(:,:)  ! kg tracer/kg dry air
  real(r8), intent(in)    :: width              ! eta difference from unit level where q = 0.1
  integer,  intent(in), optional :: rev_in      ! reverse the mixing ratio


  !Local variables
  integer k
  real(r8) alpha ! guassian width, determined by width, T
  real(r8) pdist ! pressure distance (eta.e4) from k=indx
  real(r8) T     ! desired m.r. in level specified by pdiff from k=indx
  integer rev  ! = 1 then reverse (q = 1, q(k=indx) = 0 )




  rev = 0
  if (present(rev_in)) then
     if (rev_in == 1) then
        rev = 1
     endif
  endif


!  write(iulog,*)'TR SMOOTH indx pref_mid(indx)',indx,pref_mid(indx)
!  write(iulog,67)'TR SMOOTH ','k','pref_mid(k)','pdist','-a*(pd^2)','e^-a*(pd^2)'

  T = 0.1_r8
  alpha = -log(T)/(width*1.e4_r8)**2  ! s.t. in level width from indx, mr = T

!  alpha = 3.e-8  ! m.r. ~ 0.1 in adjacent levels, where change eta ~ 0.08

  do k=1,pver
     pdist = pref_mid(k) - pref_mid(indx)

     if ( rev == 1 ) then
        q(:,k) = 1.0_r8 - exp(-alpha*(pdist**2))
     else
        q(:,k) =       exp(-alpha*(pdist**2))
  endif
!     write(iulog,66)'TR SMOOTH ',k,pref_mid(k),pdist,-alpha*pdist**2,q(1,k)
  end do
  
66 format (a15,i4,3f15.2,g15.6)
67 format (a15,a4,4a15)
  

end subroutine setsmoothtr


!======================================================================

integer function setpindxtr(pmb)

  ! find the index of layer nearest pmb

  use ref_pres, only : pref_mid

  implicit none
  
  
  real(r8) pmb
  real(r8) pmin, pdist
  integer indx, k
  
  indx = 0
  pmin = 1.e36_r8
  pdist = 1.e36_r8
  do k=1,pver
     pdist = abs(pref_mid(k) - pmb*100._r8)
     if (pdist < pmin) then
        indx = k
        pmin = pdist
     end if
  end do
  !write(iulog,*) ' index near', pmb, ' is ', indx
  setpindxtr = indx

end function setpindxtr

!======================================================================
!======================================================================

subroutine tend_tr(m,ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of radon tracer
! 
! Method:
!-------------------------Code History----------------------------------
!
! tracers.F90:rndecay()
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
!
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: m                 ! tracer number
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! radon mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg tracer /(s kg moist air))
!--------------------------Local Variables------------------------------

  integer nbase

  if ( m > trac_ncnstmx ) then
      write(iulog,*) 'tracers_suite:tend_tr()'
      write(iulog,*) ' asked to calculate tendency for tracer number ',m
      write(iulog,*) ' but there are only trac_ncnstmx = ',trac_ncnstmx,' tracers'
      call endrun('tracers_suite.F90:tend_tr() L457')
   endif

   nbase = mod(m-1,trac_names) + 1

   if ( nbase == 1 ) then
      call tend_lw(ncol, q, deltat, snk)
   else if ( nbase == 2 ) then
      call tend_md(ncol, q, deltat, snk)
   else if ( nbase == 3 ) then
      call tend_hi(ncol, q, deltat, snk)
   else if ( nbase == 4 ) then
      call tend_md(ncol, q, deltat, snk)
   else if ( nbase == 5 ) then
      call tend_un(ncol, q, deltat, snk)


   else
      write(iulog,*) 'tracers_suite:tend_tr()'
      write(iulog,*) 'no tendency routine specified for tracer',nbase
      call endrun
   endif
      


end subroutine tend_tr

!======================================================================

subroutine tend_lw(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 1 (low) (null)
! 
! Method: null
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! low mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg q /(s kg moist air))
!--------------------------Local Variables------------------------------
!
   snk = 0._r8
!   write(iulog,*)'suite 1 tend_lw min/max snk',minval(snk),maxval(snk)

end subroutine tend_lw

!======================================================================

subroutine tend_md(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 2 (med) = null
! 
! Method: null
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! low mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg q /(s kg moist air))
!--------------------------Local Variables------------------------------
   snk = 0._r8;
!   write(iulog,*)'suite 1 tend_md min/max snk',minval(snk),maxval(snk)


end subroutine tend_md

!======================================================================

subroutine tend_hi(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 3 (high) =  null
! 
! Method:  null 
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    !  mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg tracer /(s kg moist air))
!--------------------------Local Variables------------------------------

   snk = 0._r8

!   write(iulog,*)'suite 1 tend_hi min/max snk',minval(snk),maxval(snk)


end subroutine tend_hi

!======================================================================

subroutine tend_un(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 2 (unit) = null
! 
! Method: null
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! radon mixing ratio (kg/(kg dry air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg rn /(s kg dry air))
!--------------------------Local Variables------------------------------
   snk = 0._r8;
   !	write(iulog,*)'suite 2 tend_un min/max snk',minval(snk),maxval(snk)


end subroutine tend_un

!======================================================================
!======================================================================

subroutine flux_tr(m,ncol, lchnk, landfrac, flux )
!----------------------------------------------------------------------- 
!

!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments-------------------------------------

   integer,  intent(in)  :: m               ! tracer number
   integer,  intent(in)  :: ncol            ! number of atmospheric columns in chunk
   integer , intent(in)  :: lchnk           ! current identifier
   real(r8), intent(in)  :: landfrac(pcols) ! landfraction
   real(r8), intent(out) :: flux(pcols)     ! specified radon flux in kg/m^2/s

!--------------------------Local Variables------------------------------


   if ( m > trac_ncnstmx ) then
      write(iulog,*) 'tracers_suite:flux_tr()'
      write(iulog,*) ' asked to calculate flux for tracer number ',m
      write(iulog,*) ' but there are only trac_ncnstmx = ',trac_ncnstmx,' tracers'
      call endrun
   endif
   
   ! flux is null for all tracers
   flux = 0._r8


end subroutine flux_tr


end module tracers_suite
