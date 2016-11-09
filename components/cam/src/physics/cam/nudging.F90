
module nudging
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses.
!
! Author: Patrick Callaghan
!
! Description:
!    This module assumes that the user has {U,V,T,Q,PS} analyses which
!    have been preprocessed onto the current model grid and are stored
!    in individual files which are indexed with respect to year, month,
!    day, and second of the day. When the model is inbetween the given
!    begining and ending times, forcing is added to nudge the model toward
!    the appropriate analyses values. After the model passes the ending
!    analyses time, the forcing discontinues.
!
! Revisions:
!    01/14/13 - Modified to manage 'GAPS' in analyses data. For now the
!               approach is to coast through the gaps...  If a given
!               analyses file is missing, nudging is turned off for
!               that interval of time. Once an analyses file is found,
!               the Nudging is switched back on.
!    02/22/13 - Modified to add functionality for FV and EUL dynamical
!               cores.
!    03/03/13 - For ne120 runs, the automatic arrays used for reading in
!               U,V,T,Q,PS values were putting too much of a burden on the
!               stack memory. Until Parallel I/O is implemented, the impact
!               on the stack was reduced by using only one automatic array
!               to read in and scatter the data.
!    04/01/13 - Added Heaviside window function for localized nudging
!    04/10/13 - Modified call to physics_ptend_init() to accomodate the
!               new interface (in CESM1_2_BETA05).
!    05/06/13 - 'WRAP_NF' was modified from a generic interface so that
!               now it can only read in 1D arrays from netCDF files.
!               To eliminate errors from future meddling of this sort, all
!               refenences to the 'wrap_nf' module were removed and replaced
!               with direct nf90 calls.
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!
!    The nudging of the model toward the analyses data is controlled by
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. The strength of the nudging is
!    specified as a fractional coeffcient between [0,1]. The spatial distribution 
!    is specified with a profile index:
!
!        (U,V,T,Q) Profiles:      0 == OFF      (No Nudging of this variable)
!        -------------------      1 == CONSTANT (Spatially Uniform Nudging)
!                                 2 == HEAVISIDE WINDOW FUNCTION
!
!        (PS) Profiles:           0 == OFF (Not Implemented)
!        -------------------      1 == N/A (Not Implemented)
!
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 14 parameters:
!        Nudge_Hwin_lat0:     Provide the horizontal center of the window in degrees. 
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value leads to a smoother transition.
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Nudge_Vwin_Hindex:   range from [0,(NCOL+1)]. The transition lengths are also 
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NCOL+1), and the transition 
!                             lengths should be set to 0.1
!
!        Nudge_Hwin_lo:       For a given set of spatial parameters, the raw window 
!        Nudge_Hwin_hi:       function may not span the range [0,1], so those values are 
!        Nudge_Vwin_lo:       mapped to the range of values specified in by the user. 
!        Nudge_Vwin_hi:       The 'hi' values are mapped to the maximum of the raw window 
!                             function and 'lo' values are mapped to its minimum. 
!                             Typically the 'hi' values will be set equal to 1, and the 
!                             'lo' values set equal 0 or the desired window minimum. 
!                             Specifying the 'lo' value as 1 and the 'hi' value as 0 acts 
!                             to invert the window function. For a properly specified
!                             window its maximum should be equal to 1: MAX('lo','hi')==1
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lo = 0.               Nudge_Vwin_lo = 0.
!                        Nudge_Hwin_hi = 1.               Nudge_Vwin_hi = 1.
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.1
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.1
!                        Nudge_Hwin_lonWidth = 999.
!                        Nudge_Hwin_lonDelta = 1.0
!
!                 If on the other hand one desired to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_lo = 1.
!                        Nudge_Hwin_hi = 0.
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!      Nudge_Path          - CHAR path to the analyses files.
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.
!      Model_Times_Per_Day - INT Number of times to update the model state (used for nudging) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U.
!                                    Utau=(Nudge_Ucoef/analyses_timestep)
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V.
!                                    Vtau=(Nudge_Vcoef/analyses_timestep)
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T.
!                                    Ttau=(Nudge_Tcoef/analyses_timestep)
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q.
!                                    Qtau=(Nudge_Qcoef/analyses_timestep)
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS.
!                                    PStau=(Nudge_PScoef/analyses_timestep)
!      Nudge_Beg_Year      - INT nudging begining year.
!      Nudge_Beg_Month     - INT nudging begining month.
!      Nudge_Beg_Day       - INT nudging begining day.
!      Nudge_End_Year      - INT nudging ending year.
!      Nudge_End_Month     - INT nudging ending month.
!      Nudge_End_Day       - INT nudging ending day.
!      Nudge_Hwin_lo       - REAL value mapped to RAW horizontal window minimum. [0]
!      Nudge_Hwin_hi       - REAL value mapped to RAW horizontal window maximum. [1]
!      Nudge_Vwin_lo       - REAL value mapped to RAW vertical window minimum.   [0]
!      Nudge_Vwin_hi       - REAL value mapped to RAW vertical window maximum.   [1]
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length
!      Nudge_Vwin_Hdelta   - REAL HI transition length
!    /
!
!================
!  DIAG NOTE:
!================
!   The interface for reading and using analyses data is not complete for the FV
!   dynamical core. Wind values stored in the available data set are the values
!   on the staggered grid US,VS rather than U,V. To test the implementation of
!   the nudging for this case, the US,VS values were read in a loaded as if they
!   were U,V. The implementation of this hack is tagged with 'DIAG' where code
!   changed are needed to undo and fix what I have done.
!================
!
! TO DO:
! -----------
!    ** Currently the surface pressure is read in, but there is no forcing
!       meachnism implemented.
!    ** Analyses data is read in and then distributed to processing elements
!       via 'scatted_field_to_chunk' calls. The SE's want this to be changed
!       to parallel I/O calls.
!    ** Possibly implement time variation to nudging coeffcients, so that
!       rather than just bashing the model with a sledge hammer, the user has the
!       option to ramp up the nudging coefs over a startup time frame via a
!       heavyside step function.
!
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,only:timemgr_time_ge,timemgr_time_inc,get_curr_date,dtime
  use phys_grid   ,only:scatter_field_to_chunk
  use cam_abortutils  ,only:endrun
  use spmd_utils  ,only:masterproc
  use cam_logfile ,only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Nudge_Model,Nudge_ON
  public:: nudging_readnl
  public:: nudging_init
  public:: nudging_timestep_init
  public:: nudging_timestep_tend
  private::nudging_update_analyses_se
  private::nudging_update_analyses_eul
  private::nudging_update_analyses_fv
  private::nudging_set_PSprofile
  private::nudging_set_profile

  ! Nudging Parameters
  !--------------------
  logical::         Nudge_Model       =.false.
  logical::         Nudge_ON          =.false.
  logical::         Nudge_File_Present=.false.
  logical::         Nudge_Initialized =.false.
  character(len=cl) Nudge_Path
  character(len=cs) Nudge_File,Nudge_File_Template
  integer           Nudge_Times_Per_Day
  integer           Model_Times_Per_Day
  real(r8)          Nudge_Ucoef,Nudge_Vcoef
  integer           Nudge_Uprof,Nudge_Vprof
  real(r8)          Nudge_Qcoef,Nudge_Tcoef
  integer           Nudge_Qprof,Nudge_Tprof
  real(r8)          Nudge_PScoef
  integer           Nudge_PSprof
  integer           Nudge_Beg_Year ,Nudge_Beg_Month
  integer           Nudge_Beg_Day  ,Nudge_Beg_Sec
  integer           Nudge_End_Year ,Nudge_End_Month
  integer           Nudge_End_Day  ,Nudge_End_Sec
  integer           Nudge_Curr_Year,Nudge_Curr_Month
  integer           Nudge_Curr_Day ,Nudge_Curr_Sec
  integer           Nudge_Next_Year,Nudge_Next_Month
  integer           Nudge_Next_Day ,Nudge_Next_Sec
  integer           Nudge_Step
  integer           Model_Curr_Year,Model_Curr_Month
  integer           Model_Curr_Day ,Model_Curr_Sec
  integer           Model_Next_Year,Model_Next_Month
  integer           Model_Next_Day ,Model_Next_Sec
  integer           Model_Step
  real(r8)          Nudge_Hwin_lo
  real(r8)          Nudge_Hwin_hi
  real(r8)          Nudge_Hwin_lat0
  real(r8)          Nudge_Hwin_latWidth
  real(r8)          Nudge_Hwin_latDelta
  real(r8)          Nudge_Hwin_lon0
  real(r8)          Nudge_Hwin_lonWidth
  real(r8)          Nudge_Hwin_lonDelta
  real(r8)          Nudge_Vwin_lo
  real(r8)          Nudge_Vwin_hi
  real(r8)          Nudge_Vwin_Hindex
  real(r8)          Nudge_Vwin_Hdelta
  real(r8)          Nudge_Vwin_Lindex
  real(r8)          Nudge_Vwin_Ldelta
  real(r8)          Nudge_Hwin_latWidthH
  real(r8)          Nudge_Hwin_lonWidthH
  real(r8)          Nudge_Hwin_max
  real(r8)          Nudge_Hwin_min

  ! Nudging State Arrays
  !-----------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
!DIAG
  integer Nudge_slat
!DIAG
  real(r8),allocatable::Target_U(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q(:,:,:)     !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS(:,:)      !(pcols,begchunk:endchunk)
  real(r8),allocatable::Model_U(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_V(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_T(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_Q(:,:,:)      !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Model_PS(:,:)       !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Utau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Ttau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qtau(:,:,:)   !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PStau(:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable::Nudge_Ustep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Vstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Tstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_Qstep(:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Nudge_PSstep(:,:)   !(pcols,begchunk:endchunk)

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
   !
   ! NUDGING_READNL: Initialize default values controlling the Nudging
   !                 process. Then read namelist values to override
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                         Nudge_File_Template,Nudge_Times_Per_Day,      &
                         Model_Times_Per_Day,                          &
                         Nudge_Ucoef ,Nudge_Uprof,                     &
                         Nudge_Vcoef ,Nudge_Vprof,                     &
                         Nudge_Qcoef ,Nudge_Qprof,                     &
                         Nudge_Tcoef ,Nudge_Tprof,                     &
                         Nudge_PScoef,Nudge_PSprof,                    &
                         Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                         Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                         Nudge_Hwin_lo,Nudge_Hwin_hi,                  &
                         Nudge_Vwin_lo,Nudge_Vwin_hi,                  &
                         Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                         Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                         Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                         Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                         Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_ON          =.false.
   Nudge_File_Present=.false.
   Nudge_Beg_Sec=0
   Nudge_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model        =.false.
   Nudge_Path         ='./Data/YOTC_ne30np4_001/'
   Nudge_File_Template='YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Nudge_Times_Per_Day=4
   Model_Times_Per_Day=4
   Nudge_Ucoef  =0._r8
   Nudge_Vcoef  =0._r8
   Nudge_Qcoef  =0._r8
   Nudge_Tcoef  =0._r8
   Nudge_PScoef =0._r8
   Nudge_Uprof  =0
   Nudge_Vprof  =0
   Nudge_Qprof  =0
   Nudge_Tprof  =0
   Nudge_PSprof =0
   Nudge_Beg_Year =2008
   Nudge_Beg_Month=5
   Nudge_Beg_Day  =1
   Nudge_End_Year =2008
   Nudge_End_Month=9
   Nudge_End_Day  =1
   Nudge_Hwin_lo      =0.0_r8
   Nudge_Hwin_hi      =1.0_r8
   Nudge_Hwin_lat0    =0._r8
   Nudge_Hwin_latWidth=9999._r8
   Nudge_Hwin_latDelta=1.0_r8
   Nudge_Hwin_lon0    =180._r8
   Nudge_Hwin_lonWidth=9999._r8
   Nudge_Hwin_lonDelta=1.0_r8
   Nudge_Vwin_lo      =0.0_r8
   Nudge_Vwin_hi      =1.0_r8
   Nudge_Vwin_Hindex  =float(pver+1)
   Nudge_Vwin_Hdelta  =0.1_r8
   Nudge_Vwin_Lindex  =0.0_r8
   Nudge_Vwin_Ldelta  =0.1_r8

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'nudging_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Check for valid namelist values
   !----------------------------------
   if((max(Nudge_Hwin_lo,Nudge_Hwin_hi).ne.1.0).or. &
      (max(Nudge_Vwin_lo,Nudge_Vwin_hi).ne.1.0)   ) then
     write(iulog,*) 'NUDGING: The window function must have a maximum value of 1'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lo=',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Hwin_hi=',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING:  Nudge_Vwin_lo=',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING:  Nudge_Vwin_hi=',Nudge_Vwin_hi
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lat0.lt.-90.).or.(Nudge_Hwin_lat0.gt.+90.)) then
     write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0.lt.0.).or.(Nudge_Hwin_lon0.ge.360.)) then
     write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                         .or. &
      (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0.).or. &
      (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0.)   ) then
     write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta.le.0.).or.(Nudge_Hwin_lonDelta.le.0.).or. &
      (Nudge_Vwin_Hdelta  .le.0.).or.(Nudge_Vwin_Ldelta  .le.0.)    ) then
     write(iulog,*) 'NUDGING: Window Deltas must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if((Nudge_Hwin_latWidth.le.0.).or.(Nudge_Hwin_lonWidth.le.0.)) then
     write(iulog,*) 'NUDGING: Window widths must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Nudge_Path         ,len(Nudge_Path)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template,len(Nudge_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_File_Present , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Ucoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Tcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Qcoef    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PScoef   , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Uprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Vprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Tprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Qprof    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_PSprof   , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Year , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec  , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init
   !
   ! NUDGING_INIT: Allocate space and initialize Nudging values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld
   use shr_const_mod ,only: SHR_CONST_PI

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_PS',pcols*((endchunk-begchunk)+1))

   ! Allocate Space for spatial dependence of
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   allocate(Nudge_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Ttau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ttau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Nudge_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Tstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Tstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PSstep',pcols*((endchunk-begchunk)+1))

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld('Nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
   call addfld('Nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
   call addfld('Nudge_T',(/ 'lev' /),'A','cp*K/s' ,'T Nudging Tendency')
   call addfld('Nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')

   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Nudge_Step.
     !--------------------------------------------------------
     Model_Step=86400/Model_Times_Per_Day
     Nudge_Step=86400/Nudge_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'NUDGING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Nudge_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be more than Nudge_Step'
       write(iulog,*) 'NUDGING:  Setting Model_Step=Nudge_Step, Nudge_Step=',Nudge_Step
       write(iulog,*) ' '
       Model_Step=Nudge_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Nudge_nlon=hdim1_d
     Nudge_nlat=hdim2_d
     Nudge_ncol=hdim1_d*hdim2_d
     Nudge_nlev=pver
!DIAG
     Nudge_slat=Nudge_nlat-1
!DIAG

     ! Check the time relative to the nudging window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
     call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Nudge_End_Sec,Before_End)

     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Nudge_Next_Year =Year
       Nudge_Next_Month=Month
       Nudge_Next_Day  =Day
       Nudge_Next_Sec  =(Sec/Nudge_Step)*Nudge_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Nudge_Beg_Year
       Model_Next_Month=Nudge_Beg_Month
       Model_Next_Day  =Nudge_Beg_Day
       Model_Next_Sec  =Nudge_Beg_Sec
       Nudge_Next_Year =Nudge_Beg_Year
       Nudge_Next_Month=Nudge_Beg_Month
       Nudge_Next_Day  =Nudge_Beg_Day
       Nudge_Next_Sec  =Nudge_Beg_Sec
     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       Nudge_Model=.false.
       Nudge_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
       write(iulog,*) 'NUDGING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function
     !----------------------------------------
     lonp= 180.
     lon0=   0.
     lonn=-180.
     latp=  90.-Nudge_Hwin_lat0
     lat0=   0.
     latn= -90.-Nudge_Hwin_lat0

     Nudge_Hwin_lonWidthH=Nudge_Hwin_lonWidth/2.
     Nudge_Hwin_latWidthH=Nudge_Hwin_latWidth/2.

     Val1_p=(1.+tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2.
     Val2_p=(1.+tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2.
     Val3_p=(1.+tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2.
     Val4_p=(1.+tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2.

     Val1_0=(1.+tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2.
     Val2_0=(1.+tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2.
     Val3_0=(1.+tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2.
     Val4_0=(1.+tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2.

     Val1_n=(1.+tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2.
     Val2_n=(1.+tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2.
     Val3_n=(1.+tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2.
     Val4_n=(1.+tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2.

     Nudge_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Nudge_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialization is done,
     !--------------------------
     Nudge_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('NUDGING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'NUDGING: Nudge_Model=',Nudge_Model
     write(iulog,*) 'NUDGING: Nudge_Path=',Nudge_Path
     write(iulog,*) 'NUDGING: Nudge_File_Template=',Nudge_File_Template
     write(iulog,*) 'NUDGING: Nudge_Times_Per_Day=',Nudge_Times_Per_Day
     write(iulog,*) 'NUDGING: Model_Times_Per_Day=',Model_Times_Per_Day
     write(iulog,*) 'NUDGING: Nudge_Step=',Nudge_Step
     write(iulog,*) 'NUDGING: Model_Step=',Model_Step
     write(iulog,*) 'NUDGING: Nudge_Ucoef  =',Nudge_Ucoef
     write(iulog,*) 'NUDGING: Nudge_Vcoef  =',Nudge_Vcoef
     write(iulog,*) 'NUDGING: Nudge_Qcoef  =',Nudge_Qcoef
     write(iulog,*) 'NUDGING: Nudge_Tcoef  =',Nudge_Tcoef
     write(iulog,*) 'NUDGING: Nudge_PScoef =',Nudge_PScoef
     write(iulog,*) 'NUDGING: Nudge_Uprof  =',Nudge_Uprof
     write(iulog,*) 'NUDGING: Nudge_Vprof  =',Nudge_Vprof
     write(iulog,*) 'NUDGING: Nudge_Qprof  =',Nudge_Qprof
     write(iulog,*) 'NUDGING: Nudge_Tprof  =',Nudge_Tprof
     write(iulog,*) 'NUDGING: Nudge_PSprof =',Nudge_PSprof
     write(iulog,*) 'NUDGING: Nudge_Beg_Year =',Nudge_Beg_Year
     write(iulog,*) 'NUDGING: Nudge_Beg_Month=',Nudge_Beg_Month
     write(iulog,*) 'NUDGING: Nudge_Beg_Day  =',Nudge_Beg_Day
     write(iulog,*) 'NUDGING: Nudge_End_Year =',Nudge_End_Year
     write(iulog,*) 'NUDGING: Nudge_End_Month=',Nudge_End_Month
     write(iulog,*) 'NUDGING: Nudge_End_Day  =',Nudge_End_Day
     write(iulog,*) 'NUDGING: Nudge_Hwin_lo       =',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING: Nudge_Hwin_hi       =',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING: Nudge_Hwin_lat0     =',Nudge_Hwin_lat0
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidth =',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_latDelta =',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_lon0     =',Nudge_Hwin_lon0
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidth =',Nudge_Hwin_lonWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonDelta =',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_lo       =',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING: Nudge_Vwin_hi       =',Nudge_Vwin_hi
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hindex   =',Nudge_Vwin_Hindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hdelta   =',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Lindex   =',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Ldelta   =',Nudge_Vwin_Ldelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidthH=',Nudge_Hwin_latWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidthH=',Nudge_Hwin_lonWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_max      =',Nudge_Hwin_max
     write(iulog,*) 'NUDGING: Nudge_Hwin_min      =',Nudge_Hwin_min
     write(iulog,*) 'NUDGING: Nudge_Initialized   =',Nudge_Initialized
     write(iulog,*) ' '
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Step          , 1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Model         , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON            , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized   , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ncol          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlev          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlon          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlat          , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_max      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_min      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidthH, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidthH, 1, mpir8 , 0, mpicom)
!DIAG
   call mpibcast(Nudge_slat       , 1, mpiint, 0, mpicom)
!DIAG
#endif

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call nudging_set_profile(rlat,rlon,Nudge_Uprof,Wprof,pver)
       Nudge_Utau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Vprof,Wprof,pver)
       Nudge_Vtau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Tprof,Wprof,pver)
       Nudge_Ttau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Qprof,Wprof,pver)
       Nudge_Qtau(icol,:,lchnk)=Wprof(:)

       Nudge_PStau(icol,lchnk)=nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
     end do
     Nudge_Utau(:ncol,:pver,lchnk) =                             &
     Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef/float(Nudge_Step)
     Nudge_Vtau(:ncol,:pver,lchnk) =                             &
     Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef/float(Nudge_Step)
     Nudge_Ttau(:ncol,:pver,lchnk) =                             &
     Nudge_Ttau(:ncol,:pver,lchnk) * Nudge_Tcoef/float(Nudge_Step)
     Nudge_Qtau(:ncol,:pver,lchnk) =                             &
     Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef/float(Nudge_Step)
     Nudge_PStau(:ncol,lchnk)=                             &
     Nudge_PStau(:ncol,lchnk)* Nudge_PScoef/float(Nudge_Step)

     Nudge_Ustep(:pcols,:pver,lchnk)=0._r8
     Nudge_Vstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Tstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Qstep(:pcols,:pver,lchnk)=0._r8
     Nudge_PSstep(:pcols,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8
     Target_T(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8
     Target_PS(:pcols,lchnk)=0._r8
   end do

   ! End Routine
   !------------
   return
  end subroutine ! nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   !
   ! NUDGING_TIMESTEP_INIT:
   !                 Check the current time and update Model/Nudging
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
   !===============================================================
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use physconst    ,only: cpair

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Nudge,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_Initialized) then
     call endrun('nudging_timestep_init:: Nudging NOT Initialized')
   endif

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'NUDGING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Model_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
       Model_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
       Model_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
       Model_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
       Model_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
     end do
   endif

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_Next_Year*10000) + (Nudge_Next_Month*100) + Nudge_Next_Day
   call timemgr_time_ge(YMD1,Nudge_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)

   if((Before_End).and.(Update_Nudge)) then
     ! Increment the Nudge times by the current interval
     !---------------------------------------------------
     Nudge_Curr_Year =Nudge_Next_Year
     Nudge_Curr_Month=Nudge_Next_Month
     Nudge_Curr_Day  =Nudge_Next_Day
     Nudge_Curr_Sec  =Nudge_Next_Sec
     YMD1=(Nudge_Curr_Year*10000) + (Nudge_Curr_Month*100) + Nudge_Curr_Day
     call timemgr_time_inc(YMD1,Nudge_Curr_Sec,              &
                           YMD2,Nudge_Next_Sec,Nudge_Step,0,0)
     Nudge_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Nudge_Next_Year*10000)
     Nudge_Next_Month=(YMD2/100)
     Nudge_Next_Day  = YMD2-(Nudge_Next_Month*100)

     ! Update the Nudge arrays with analysis
     ! data at the NEXT time
     !-----------------------------------------------
     Nudge_File=interpret_filename_spec(Nudge_File_Template      , &
                                         yr_spec=Nudge_Next_Year , &
                                        mon_spec=Nudge_Next_Month, &
                                        day_spec=Nudge_Next_Day  , &
                                        sec_spec=Nudge_Next_Sec    )
     if(masterproc) then
      write(iulog,*) 'NUDGING: Reading analyses:',trim(Nudge_Path)//trim(Nudge_File)
     endif

     ! How to manage MISSING values when there are 'Gaps' in the analyses data?
     !  Check for analyses file existence. If it is there, then read data.
     !  If it is not, then issue a warning and switch off nudging to 'coast'
     !  thru the gap.
     !------------------------------------------------------------------------
     if(dycore_is('UNSTRUCTURED')) then
       call nudging_update_analyses_se(trim(Nudge_Path)//trim(Nudge_File))
     elseif(dycore_is('EUL')) then
       call nudging_update_analyses_eul(trim(Nudge_Path)//trim(Nudge_File))
     else !if(dycore_is('LR')) then
       call nudging_update_analyses_fv(trim(Nudge_Path)//trim(Nudge_File))
     endif
   endif

   !-------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between
   ! beginning and ending times, and the analyses file exists.
   !-------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(Nudge_File_Present) then
       Nudge_ON=.true.
     else
       Nudge_ON=.false.
       if(masterproc) then
         write(iulog,*) 'NUDGING: WARNING - analyses file NOT FOUND. Switching '
         write(iulog,*) 'NUDGING:           nudging OFF to coast thru the gap. '
       endif
     endif
   else
     Nudge_ON=.false.
   endif

   !-------------------------------------------------------
   ! HERE Implement time dependence of Nudging Coefs HERE
   !-------------------------------------------------------



   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Nudge).or.(Update_Model))) then
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Nudge_Ustep(:ncol,:pver,lchnk)=(  Target_U(:ncol,:pver,lchnk)  &
                                         -Model_U(:ncol,:pver,lchnk)) &
                                      *Nudge_Utau(:ncol,:pver,lchnk)
       Nudge_Vstep(:ncol,:pver,lchnk)=(  Target_V(:ncol,:pver,lchnk)  &
                                         -Model_V(:ncol,:pver,lchnk)) &
                                      *Nudge_Vtau(:ncol,:pver,lchnk)
       Nudge_Tstep(:ncol,:pver,lchnk)=(  Target_T(:ncol,:pver,lchnk)  &
                                         -Model_T(:ncol,:pver,lchnk)) &
                                      *Nudge_Ttau(:ncol,:pver,lchnk)*cpair
       Nudge_Qstep(:ncol,:pver,lchnk)=(  Target_Q(:ncol,:pver,lchnk)  &
                                         -Model_Q(:ncol,:pver,lchnk)) &
                                      *Nudge_Qtau(:ncol,:pver,lchnk)
       Nudge_PSstep(:ncol,     lchnk)=(  Target_PS(:ncol,lchnk)  &
                                         -Model_PS(:ncol,lchnk)) &
                                      *Nudge_PStau(:ncol,lchnk)
     end do

     !******************
     ! DIAG
     !******************
!    if(masterproc) then
!     write(iulog,*) 'PFC: Target_T(1,:pver,begchunk)=',Target_T(1,:pver,begchunk)  
!     write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!     write(iulog,*) 'PFC: Nudge_Tstep(1,:pver,begchunk)=',Nudge_Tstep(1,:pver,begchunk)
!     write(iulog,*) 'PFC: Nudge_Xstep arrays updated:'
!    endif
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Nudge_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Nudge_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Nudge_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Nudge_Tstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Nudge_Qstep(:ncol,:pver,lchnk)

     call outfld('Nudge_U',phys_tend%u          ,pcols,lchnk)
     call outfld('Nudge_V',phys_tend%v          ,pcols,lchnk)
     call outfld('Nudge_T',phys_tend%s          ,pcols,lchnk)
     call outfld('Nudge_Q',phys_tend%q(1,1,indw),pcols,lchnk)
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_se(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_SE:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid,varid
   real(r8) Xanal(Nudge_ncol,Nudge_nlev)
   real(r8) PSanal(Nudge_ncol)
   real(r8) Lat_anal(Nudge_ncol)
   real(r8) Lon_anal(Nudge_ncol)

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'ncol',varid)
!     call wrap_inq_dimlen(ncid,varid,ncol)
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     if((Nudge_ncol.ne.ncol).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_se: ncol=',ncol,' Nudge_ncol=',Nudge_ncol
      write(iulog,*) 'ERROR: nudging_update_analyses_se: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_se: analyses dimension mismatch')
     endif

     ! Read in and scatter data arrays
     !----------------------------------
!     call wrap_inq_varid    (ncid,'U',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal ,Target_U)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'V',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal ,Target_V)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'T',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal ,Target_T)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'Q',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_ncol,Xanal ,Target_Q)

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !-----------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
!  call scatter_field_to_chunk(1,         1,1,Nudge_ncol,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_se
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_eul(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_EUL:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'lon',varid)
!     call wrap_inq_dimlen(ncid,varid,nlon)
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_dimid (ncid,'lat',varid)
!     call wrap_inq_dimlen(ncid,varid,nlat)
     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_eul: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_eul: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
!     call wrap_inq_varid    (ncid,'U',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_U)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'V',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_V)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'T',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_T)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'Q',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_Q)

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !-----------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
!  call scatter_field_to_chunk(1,         1,1,Nudge_nlon,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_eul
  !================================================================


  !================================================================
  subroutine nudging_update_analyses_fv(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES_FV:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
!   use wrap_nf
   use ppgrid ,only: pver
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
!DIAG
   real(r8) Uanal(Nudge_nlon,Nudge_slat,Nudge_nlev)
!DIAG

   ! Check the existence of the analyses file; broadcast the file status to
   ! all the other MPI nodes. If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
!     call wrap_inq_dimid (ncid,'lon',varid)
!     call wrap_inq_dimlen(ncid,varid,nlon)
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_dimid (ncid,'lat',varid)
!     call wrap_inq_dimlen(ncid,varid,nlat)
     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_dimid (ncid,'lev',varid)
!     call wrap_inq_dimlen(ncid,varid,plev)
     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_varid(ncid,'lon',varid)
!     call wrap_get_var_realx(ncid,varid,Lon_anal)
     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

!     call wrap_inq_varid(ncid,'lat',varid)
!     call wrap_get_var_realx(ncid,varid,Lat_anal)
     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
!DIAG:  Dont have U, so jam US into U so tests can proceed:
!DIAG     call wrap_inq_varid    (ncid,'U',varid)
!DIAG     call wrap_get_var_realx(ncid,varid,Xanal)
!DIAG     do ilat=1,nlat
!DIAG     do ilev=1,plev
!DIAG     do ilon=1,nlon
!DIAG       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!DIAG     end do
!DIAG     end do
!DIAG     end do
!     call wrap_inq_varid    (ncid,'US',varid)
!     call wrap_get_var_realx(ncid,varid,Uanal)
     istat=nf90_inq_varid(ncid,'US',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Uanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,(nlat-1)
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Uanal(ilon,ilat,ilev)
     end do
     end do
     end do
     Xtrans(:,:,ilat)=Xtrans(:,:,ilat-1)
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_U)

   if(masterproc) then
!DIAG:  Dont have V, so jam VS into V so tests can proceed:
!DIAG     call wrap_inq_varid    (ncid,'V',varid)
!DIAG     call wrap_get_var_realx(ncid,varid,Xanal)
!DIAG     do ilat=1,nlat
!DIAG     do ilev=1,plev
!DIAG     do ilon=1,nlon
!DIAG       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
!DIAG     end do
!DIAG     end do
!DIAG     end do
!     call wrap_inq_varid    (ncid,'VS',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'VS',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_V)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'T',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_T)

   if(masterproc) then
!     call wrap_inq_varid    (ncid,'Q',varid)
!     call wrap_get_var_realx(ncid,varid,Xanal)
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans ,Target_Q)

   if(masterproc) then
!!    call wrap_inq_varid    (ncid,'PS',varid)
!!    call wrap_get_var_realx(ncid,varid,PSanal)
!    istat=nf90_inq_varid(ncid,'PS',varid)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif
!    istat=nf90_get_var(ncid,varid,PSanal)
!    if(istat.ne.NF90_NOERR) then
!      write(iulog,*) nf90_strerror(istat)
!      call endrun ('UPDATE_ANALYSES_SE')
!    endif

     ! Close the analyses file
     !-----------------------
!     call wrap_close(ncid)
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
!  call scatter_field_to_chunk(1,         1,1,Nudge_nlon,PSanal,Target_PS)

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_fv
  !================================================================


  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   !
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge_Hwin_max.le.Nudge_Hwin_min) then
       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge_Hwin_lo,Nudge_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge_Hwin_lat0
       lonx=rlon-Nudge_Hwin_lon0
       if(lonx.gt. 180.) lonx=lonx-360.
       if(lonx.le.-180.) lonx=lonx+360.

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge_Hwin_lonWidthH+lonx)/Nudge_Hwin_lonDelta
       lon_hi=(Nudge_Hwin_lonWidthH-lonx)/Nudge_Hwin_lonDelta
       lat_lo=(Nudge_Hwin_latWidthH+latx)/Nudge_Hwin_latDelta
       lat_hi=(Nudge_Hwin_latWidthH-latx)/Nudge_Hwin_latDelta
       Hcoef=((1.+tanh(lon_lo))/2.)*((1.+tanh(lon_hi))/2.) &
            *((1.+tanh(lat_lo))/2.)*((1.+tanh(lat_hi))/2.)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
       Hcoef=(1.-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
       lev_hi=(Nudge_Vwin_Hindex-float(ilev))/Nudge_Vwin_Hdelta
       Wprof(ilev)=((1.+tanh(lev_lo))/2.)*((1.+tanh(lev_hi))/2.)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if(Vmax.le.Vmin) then
       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_Vwin_lo + Wprof(:)*(Nudge_Vwin_hi-Nudge_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_set_profile
  !================================================================

  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   !
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0
   elseif(Nudge_PSprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_set_PSprofile
  !================================================================

end module nudging
