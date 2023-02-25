module RtmMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! Mosart Routing Model
!
! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_flush
  use shr_const_mod   , only : SHR_CONST_PI, SHR_CONST_CDAY
  use rof_cpl_indices , only : nt_rtm, rtm_tracers, KW, DW
  use RtmSpmd         , only : masterproc, npes, iam, mpicom_rof, ROFID, mastertask, &
                               MPI_REAL8,MPI_INTEGER,MPI_CHARACTER,MPI_LOGICAL,MPI_MAX
  use RtmVar          , only : re, spval, rtmlon, rtmlat, iulog, ice_runoff, &
                               frivinp_rtm, finidat_rtm, nrevsn_rtm,rstraflag,ngeom,nlayers,rinittemp, &
                               nsrContinue, nsrBranch, nsrStartup, nsrest, &
                               inst_index, inst_suffix, inst_name, wrmflag, inundflag, &
                               smat_option, decomp_option, barrier_timers, heatflag, sediflag, &
                               isgrid2d, data_bgc_fluxes_to_ocean_flag, use_lnd_rof_two_way, use_ocn_rof_two_way
  use RtmFileUtils    , only : getfil, getavu, relavu
  use RtmTimeManager  , only : timemgr_init, get_nstep, get_curr_date, advance_timestep
  use RtmHistFlds     , only : RtmHistFldsInit, RtmHistFldsSet 
  use RtmHistFile     , only : RtmHistUpdateHbuf, RtmHistHtapesWrapup, RtmHistHtapesBuild, &
                               rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq,     &
                               rtmhist_avgflag_pertape, rtmhist_avgflag_pertape, & 
                               rtmhist_fincl1, rtmhist_fincl2, rtmhist_fincl3,   &
                               rtmhist_fexcl1, rtmhist_fexcl2, rtmhist_fexcl3,   &
                               max_tapes, max_namlen
  use RtmRestFile     , only : RtmRestTimeManager, RtmRestGetFile, RtmRestFileRead, &
                               RtmRestFileWrite, RtmRestFileName
  use RunoffMod       , only : RunoffInit, rtmCTL, Tctl, Tunit, TRunoff, Tpara, Theat, &
                               gsmap_r, &
                               SMatP_dnstrm, avsrc_dnstrm, avdst_dnstrm, &
                               SMatP_upstrm, avsrc_upstrm, avdst_upstrm, &
                               SMatP_direct, avsrc_direct, avdst_direct
  use MOSART_physics_mod, only : Euler
  use MOSART_physics_mod, only : updatestate_hillslope, updatestate_subnetwork, &
                                 updatestate_mainchannel
  use WRM_type_mod    , only : ctlSubwWRM, WRMUnit, StorWater
  use WRM_subw_IO_mod , only : WRM_init, WRM_computeRelease
  use MOSARTinund_PreProcs_MOD, only : calc_chnlMannCoe, preprocess_elevProf
  use MOSARTinund_Core_MOD    , only : MOSARTinund_simulate, ManningEq, ChnlFPexchg
  use RtmIO
  use mct_mod
  use perf_mod
  use pio
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini          ! Initialize MOSART grid
  public Rtmrun          ! River routing model
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: RtmFloodInit

! !PRIVATE TYPES:

! MOSART tracers
  character(len=256) :: rtm_trstr   ! tracer string

! MOSART namelists
  integer, save :: coupling_period   ! mosart coupling period
  integer, save :: delt_mosart       ! mosart internal timestep (->nsub)

! MOSART constants
  real(r8) :: cfl_scale = 1.0_r8    ! cfl scale factor, must be <= 1.0
  real(r8) :: river_depth_minimum = 1.e-4 ! gridcell average minimum river depth [m]

!global (glo)
  integer , pointer :: ID0_global(:)  ! local ID index
  integer , pointer :: dnID_global(:) ! downstream ID based on ID0
  real(r8), pointer :: area_global(:) ! area
  real(r8), pointer :: DIN_global(:)  !
  real(r8), pointer :: DIP_global(:)  !
  real(r8), pointer :: DON_global(:)  !
  real(r8), pointer :: DOP_global(:)  !
  real(r8), pointer :: DOC_global(:)  !
  real(r8), pointer :: PP_global(:)   !
  real(r8), pointer :: DSi_global(:)  !
  real(r8), pointer :: POC_global(:)  !
  real(r8), pointer :: PN_global(:)   !
  real(r8), pointer :: DIC_global(:)  !
  real(r8), pointer :: Fe_global(:)   !
  integer , pointer :: IDkey(:)       ! translation key from ID to gindex

!local (gdc)
  real(r8), save, pointer :: evel(:,:)       ! effective tracer velocity (m/s)
  real(r8), save, pointer :: flow(:,:)       ! mosart flow (m3/s)
  real(r8), save, pointer :: eroup_lagi(:,:) ! erout previous timestep (m3/s)
  real(r8), save, pointer :: eroup_lagf(:,:) ! erout current timestep (m3/s)
  real(r8), save, pointer :: erowm_regi(:,:) ! erout previous timestep (m3/s)
  real(r8), save, pointer :: erowm_regf(:,:) ! erout current timestep (m3/s)
  real(r8), save, pointer :: eroutup_avg(:,:)! eroutup average over coupling period (m3/s)
  real(r8), save, pointer :: erlat_avg(:,:)  ! erlateral average over coupling period (m3/s)

  real(r8), save :: vol_chnl2fp                 ! Total volume of flows from main channels to floodplains (for all local grid cells and all sub-steps of coupling period) (m^3).
  real(r8), save :: vol_fp2chnl                 ! Total volume of flows from floodplains to main channels (for all local grid cells and all sub-steps of coupling period) (m^3).

  real(r8), save, pointer :: fa_fp_cplPeriod(:) ! Mean area of inundated floodplain (not including channel area) for all sub-steps of coupling period (m^2).

  real(r8), save, pointer :: totalVelo_down(:)  ! Total downward flow velocity (for all local grid cells and all sub-steps of coupling period) (m/s).
  real(r8), save, pointer :: totalVelo_up(:)    ! Total upward flow velocity (for all local grid cells and all sub-steps of coupling period) (m/s).
  integer, save, pointer :: chnlNum_down(:)     ! Total number of channels with downward flow velocities (for all local grid cells and all sub-steps of coupling period) (dimensionless).
  integer, save, pointer :: chnlNum_up(:)       ! Total number of channels with upward flow velocities (for all local grid cells and all sub-steps of coupling period) (dimensionless).

! global MOSART grid
  real(r8),pointer :: rlatc(:)    ! latitude of 1d grid cell (deg)
  real(r8),pointer :: rlonc(:)    ! longitude of 1d grid cell (deg)
  real(r8),pointer :: rlats(:)    ! latitude of 1d south grid cell edge (deg)
  real(r8),pointer :: rlatn(:)    ! latitude of 1d north grid cell edge (deg)
  real(r8),pointer :: rlonw(:)    ! longitude of 1d west grid cell edge (deg)
  real(r8),pointer :: rlone(:)    ! longitude of 1d east grid cell edge (deg)

  logical :: do_rtmflood
  logical :: do_rtm

  real(r8), save :: delt_save             ! previous delt 

!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmini
!
! !INTERFACE:
  subroutine Rtmini(rtm_active,flood_active)
!
! !DESCRIPTION:
! Initialize MOSART grid, mask, decomp
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    logical, intent(out) :: rtm_active
    logical, intent(out) :: flood_active
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
! Update: T Craig, Dec 2006
!
!
! !LOCAL VARIABLES:
!EOP
    real(r8) :: effvel0 = 10.0_r8             ! default velocity (m/s)
    real(r8) :: effvel(nt_rtm)                ! downstream velocity (m/s)
    real(r8) :: edgen                         ! North edge of the direction file
    real(r8) :: edgee                         ! East edge of the direction file
    real(r8) :: edges                         ! South edge of the direction file
    real(r8) :: edgew                         ! West edge of the direction file
    integer  :: i,j,k,n,ng,g,n2,nt,nn         ! loop indices
    integer  :: i1,j1,i2,j2
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8) :: lrtmarea                      ! tmp local sum of area
    real(r8),allocatable :: tempr(:,:)        ! temporary buffer
    integer ,allocatable :: itempr(:,:)       ! temporary buffer
    integer ,allocatable :: idxocn(:)         ! downstream ocean outlet cell
    integer ,allocatable :: nupstrm(:)        ! number of upstream cells including own cell
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of gridcells on a pe
    integer ,allocatable :: nba(:)            ! number of basins on each pe
    integer ,allocatable :: nrs(:)            ! begr on each pe
    integer ,allocatable :: basin(:)          ! basin to mosart mapping
    integer  :: nmos,nmos_chk                 ! number of mosart points
    integer  :: nout,nout_chk                 ! number of basin with outlets
    integer  :: nbas,nbas_chk                 ! number of basin/ocean points
    integer  :: nrof,nrof_chk                 ! num of active mosart points
    integer  :: baspe                         ! pe with min number of mosart cells
    integer  :: maxrtm                        ! max num of rtms per pe for decomp
    integer  :: minbas,maxbas                 ! used for decomp search
    integer  :: nl,nloops                     ! used for decomp search
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: numr                          ! tot num of roff pts on all pes
    integer  :: RoutingMethod                 ! 1 = KW, 2 = DW
    integer  :: DLevelH2R                     !
    integer  :: DLevelR                       !
    integer  :: OPT_inund                     ! Options for inundation
    integer  :: OPT_trueDW                    ! Options for diffusion wave channel routing method:
    integer  :: OPT_calcNr                    ! Options to calculate channel Manning roughness coefficients : 
    real(r8) :: nr_max                        ! Max Manning coefficient for channels (when OPT_calcNr = 1, 2, 3) ( s*m^(-1/3) ).
    real(r8) :: nr_min                        ! Min Manning coefficient for channels (when OPT_calcNr = 1, 2, 3) ( s*m^(-1/3) ).
    real(r8) :: nr_uniform                    ! The uniform Manning coefficient for all channels (when OPT_calcNr = 4) ( s*m^(-1/3) ).
    real(r8) :: rdepth_max                    ! Max channel depth (used when OPT_calcNr = 1, 2) (m).
    real(r8) :: rdepth_min                    ! Min channel depth (used when OPT_calcNr = 1, 2) (m). 
    real(r8) :: rwidth_max                    ! Max channel width (used when OPT_calcNr = 3) (m).
    real(r8) :: rwidth_min                    ! Min channel width (used when OPT_calcNr = 3) (m). 
    real(r8) :: rslp_assume                   ! Use this assumed riverbed slope when the input riverbed slope <= zero (dimensionless).
    real(r8) :: minL_tribRouting              ! Min tributary channel length for using tributary routing (m).  
    integer  :: OPT_elevProf                  ! Options of elevation profile data: 1 -- Use real data; 2 -- Use hypothetical values.
    integer  :: npt_elevProf                  ! Number of dividing points in the elevation profile.
    real(r8) :: threshold_slpRatio            ! Threshold of the ratio of the lowest section's slope to the second lowest section's slope in 
                                              ! the elevation profile (used to alleviate the effect of DEM pits on elevation profiles).
    real(r8) :: dtover,dtovermax              ! ts calc temporaries
    type(file_desc_t) :: ncid                 ! netcdf file id
    integer  :: dimid                         ! netcdf dimension identifier
    integer  :: nroflnd                       ! local number of land runoff 
    integer  :: nrofocn                       ! local number of ocn runoff
    integer  :: pid,np,npmin,npmax,npint      ! log loop control
    integer  :: na,nb,ns                      ! mct sizes
    integer  :: ni,no,go                      ! tmps
    integer ,pointer  :: rgdc2glo(:)          ! temporary for initialization
    integer ,pointer  :: rglo2gdc(:)          ! temporary for initialization
    integer ,pointer  :: gmask(:)             ! global mask
    logical           :: found                ! flag
    character(len=256):: fnamer               ! name of netcdf restart file 
    character(len=256):: pnamer               ! full pathname of netcdf restart file
    character(len=256):: locfn                ! local file name
    character(len=16384) :: rList             ! list of fields for SM multiply
    integer           :: unitn                ! unit for namelist file
    integer,parameter :: dbug = 1             ! 0 = none, 1=normal, 2=much, 3=max
    logical :: lexist                         ! File exists
    character(len= 7) :: runtyp(4)            ! run type
    integer ,allocatable :: gindex(:)         ! global index
    integer           :: cnt, lsize, gsize    ! counter
    integer           :: igrow,igcol,iwgt     ! mct field indices
    character(len=256):: nlfilename_rof       ! namelist filename
    type(mct_avect)   :: avtmp, avtmpG        ! temporary avects
    type(mct_sMat)    :: sMat                 ! temporary sparse matrix, needed for sMatP

!global (glo), temporary
    integer , pointer :: ID0_global(:)  ! local ID index
    integer , pointer :: dnID_global(:) ! downstream ID based on ID0
    integer , pointer :: nUp_global(:)  ! number of upstream units
    integer , pointer :: nUp_dstrm_global(:)  ! number of units flowing into the downstream unit
    real(r8), pointer :: area_global(:) ! area

    character(len=*),parameter :: subname = '(Rtmini) '
    integer           :: rtmn                 ! total number of cells

    real(r8) :: wd_chnl                       ! Channel water depth (m).
    real(r8) :: hydrR                         ! Hydraulic radius (= wet A / wet P) (m).
    real(r8) :: v_chnl                        ! Channel flow velocity (m/s).#endif

!-----------------------------------------------------------------------

    !-------------------------------------------------------
    ! Read in mosart namelist
    !-------------------------------------------------------

    namelist /mosart_inparm / ice_runoff, do_rtm, do_rtmflood, &
         frivinp_rtm, finidat_rtm, nrevsn_rtm, coupling_period, &
         rtmhist_ndens, rtmhist_mfilt, rtmhist_nhtfrq, &
         rtmhist_fincl1,  rtmhist_fincl2, rtmhist_fincl3, &
         rtmhist_fexcl1,  rtmhist_fexcl2, rtmhist_fexcl3, &
         rtmhist_avgflag_pertape, decomp_option, wrmflag,rstraflag,ngeom,nlayers,rinittemp, &
         inundflag, smat_option, delt_mosart, barrier_timers,          &
         RoutingMethod, DLevelH2R, DLevelR, sediflag, heatflag, data_bgc_fluxes_to_ocean_flag

    namelist /inund_inparm / opt_inund, &
         opt_truedw, opt_calcnr, nr_max, nr_min, &
         nr_uniform, rdepth_max, rdepth_min, rwidth_max, rwidth_min, &
         rslp_assume, minl_tribrouting, opt_elevprof, npt_elevprof, &
         threshold_slpratio

    ! Preset values
    do_rtm      = .true.
    do_rtmflood = .false.
    ice_runoff  = .true.
    wrmflag     = .false.
    rstraflag   = .false.
    rinittemp   = 283.15_r8
    ngeom       = 50  
    nlayers     = 30                  
    inundflag   = .false.
    sediflag    = .false.
    heatflag    = .false.
    barrier_timers = .false.
    finidat_rtm = ' '
    nrevsn_rtm  = ' '
    coupling_period   = -1
    delt_mosart = 3600
    decomp_option = 'basin'
    smat_option = 'opt'
    RoutingMethod = KW
    DLevelH2R = 5
    DLevelR = 3
    data_bgc_fluxes_to_ocean_flag = .false.

    OPT_inund = 0
    OPT_trueDW = 2
    OPT_calcNr = 1                
    nr_max = 0.05_r8              
    nr_min = 0.03_r8              
    nr_uniform = 0.04_r8          
    rdepth_max = 50.0_r8          
    rdepth_min = 1.0_r8           
    rwidth_max = 5000.0_r8        
    rwidth_min = 20.0_r8          
    rslp_assume = 0.00001_r8      
    minL_tribRouting = 10.0_r8    
    OPT_elevProf = 2              
    npt_elevProf = 11             
    threshold_slpRatio = 10.0_r8  

    nlfilename_rof = "mosart_in" // trim(inst_suffix)
    inquire (file = trim(nlfilename_rof), exist = lexist)
    if ( .not. lexist ) then
       write(iulog,*) subname // ' ERROR: nlfilename_rof does NOT exist:'&
            //trim(nlfilename_rof)
       call shr_sys_abort(trim(subname)//' ERROR nlfilename_rof does not exist')
    end if
    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in mosart_inparm namelist from: ', trim(nlfilename_rof)
       open( unitn, file=trim(nlfilename_rof), status='old' )
       ier = 1
       do while ( ier /= 0 )
          read(unitn, mosart_inparm, iostat=ier)
          if (ier < 0) then
             call shr_sys_abort( subname//' encountered end-of-file on mosart_inparm read' )
          endif
       end do
       call relavu( unitn )
       if (inundflag) then
          unitn = getavu()
          write(iulog,*) 'Read in inund_inparm namelist from: ', trim(nlfilename_rof)
          open( unitn, file=trim(nlfilename_rof), status='old' )
          ier = 1
          do while ( ier /= 0 )
             read(unitn, inund_inparm, iostat=ier)
             if (ier < 0) then
                call shr_sys_abort( subname//' encountered end-of-file on inund_inparm read' )
             endif
          end do
          call relavu( unitn )
       end if
       if (.not. inundflag .and. use_lnd_rof_two_way) then
          call shr_sys_abort(trim(subname)//' inundation model must be turned on for land river two way coupling')
       end if
    end if

    call mpi_bcast (coupling_period,   1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (delt_mosart    ,   1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (RoutingMethod  ,   1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (DLevelH2R      ,   1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (DLevelR        ,   1, MPI_INTEGER, 0, mpicom_rof, ier)

    call mpi_bcast (finidat_rtm  , len(finidat_rtm)  , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (frivinp_rtm  , len(frivinp_rtm)  , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (nrevsn_rtm   , len(nrevsn_rtm)   , MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (decomp_option, len(decomp_option), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (smat_option  , len(smat_option)  , MPI_CHARACTER, 0, mpicom_rof, ier)

    call mpi_bcast (do_rtm,         1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (do_rtmflood,    1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (ice_runoff,     1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (wrmflag,        1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (sediflag,       1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (heatflag,       1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (rstraflag,      1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (rinittemp,      1, MPI_REAL8, 0, mpicom_rof, ier)
    call mpi_bcast (ngeom,          1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (nlayers,        1, MPI_INTEGER, 0, mpicom_rof, ier)
    call mpi_bcast (inundflag,      1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (use_lnd_rof_two_way, 1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (heatflag,       1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (use_ocn_rof_two_way, 1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (barrier_timers, 1, MPI_LOGICAL, 0, mpicom_rof, ier)
    call mpi_bcast (data_bgc_fluxes_to_ocean_flag, 1, MPI_LOGICAL, 0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_nhtfrq, size(rtmhist_nhtfrq), MPI_INTEGER,   0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_mfilt , size(rtmhist_mfilt) , MPI_INTEGER,   0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_ndens , size(rtmhist_ndens) , MPI_INTEGER,   0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_fexcl1, (max_namlen+2)*size(rtmhist_fexcl1), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fexcl2, (max_namlen+2)*size(rtmhist_fexcl2), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fexcl3, (max_namlen+2)*size(rtmhist_fexcl3), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl1, (max_namlen+2)*size(rtmhist_fincl1), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl2, (max_namlen+2)*size(rtmhist_fincl2), MPI_CHARACTER, 0, mpicom_rof, ier)
    call mpi_bcast (rtmhist_fincl3, (max_namlen+2)*size(rtmhist_fincl3), MPI_CHARACTER, 0, mpicom_rof, ier)

    call mpi_bcast (rtmhist_avgflag_pertape, size(rtmhist_avgflag_pertape), MPI_CHARACTER, 0, mpicom_rof, ier)

    if (inundflag) then
       call mpi_bcast (OPT_inund,          1, MPI_INTEGER, 0, mpicom_rof, ier)
       call mpi_bcast (OPT_trueDW,         1, MPI_INTEGER, 0, mpicom_rof, ier)
       call mpi_bcast (OPT_calcNr,         1, MPI_INTEGER, 0, mpicom_rof, ier)
       call mpi_bcast (nr_max,             1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (nr_min,             1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (nr_uniform,         1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (rdepth_max,         1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (rdepth_min,         1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (rwidth_max,         1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (rwidth_min,         1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (rslp_assume,        1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (minL_tribRouting,   1, MPI_REAL8  , 0, mpicom_rof, ier)
       call mpi_bcast (OPT_elevProf,       1, MPI_INTEGER, 0, mpicom_rof, ier)
       call mpi_bcast (npt_elevProf,       1, MPI_INTEGER, 0, mpicom_rof, ier)
       call mpi_bcast (threshold_slpRatio, 1, MPI_REAL8  , 0, mpicom_rof, ier)
    end if

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'
    runtyp(nsrBranch   + 1) = 'branch '

    if ( use_ocn_rof_two_way ) then
       RoutingMethod = DW
    end if 

    Tctl%RoutingMethod = RoutingMethod
    Tctl%DLevelH2R     = DLevelH2R
    Tctl%DLevelR       = DLevelR
    if(.not.(Tctl%RoutingMethod==KW .or. Tctl%RoutingMethod==DW)) then 
	   call shr_sys_abort('Error in routing method setup! There are only 2 options available: 1==KW, 2==DW')
    end if

    if (inundflag) then
       Tctl%OPT_inund = OPT_inund     !
       Tctl%OPT_trueDW = OPT_trueDW   ! diffusion wave method
       Tctl%OPT_calcNr = OPT_calcNr   ! method to calculate channel Manning
       Tctl%nr_max = nr_max           ! Max Manning coefficient
       Tctl%nr_min = nr_min           ! Min Manning coefficient
       Tctl%nr_uniform = nr_uniform   ! uniform Manning for all channels
       Tctl%rdepth_max = rdepth_max   ! Max channel depth
       Tctl%rdepth_min = rdepth_min   ! Min channel depth
       Tctl%rwidth_max = rwidth_max   ! Max channel width
       Tctl%rwidth_min = rwidth_min   ! Min channel width
       Tctl%rslp_assume = rslp_assume ! assumed riverbed slope if input slope<=0
       Tctl%minL_tribRouting = minL_tribRouting
       Tctl%OPT_elevProf = OPT_elevProf
       Tctl%npt_elevProf = npt_elevProf
       Tctl%threshold_slpRatio = threshold_slpRatio
    end if

    if (masterproc) then
       write(iulog,*) 'define run:'
       write(iulog,*) '   run type              = ',runtyp(nsrest+1)
      !write(iulog,*) '   case title            = ',trim(ctitle)
      !write(iulog,*) '   username              = ',trim(username)
      !write(iulog,*) '   hostname              = ',trim(hostname)
       write(iulog,*) '   coupling_period       = ',coupling_period
       write(iulog,*) '   delt_mosart           = ',delt_mosart
       write(iulog,*) '   decomp_option         = ',trim(decomp_option)
       write(iulog,*) '   smat_option           = ',trim(smat_option)
       write(iulog,*) '   wrmflag               = ',wrmflag
       write(iulog,*) '   inundflag             = ',inundflag
       write(iulog,*) '   use_lnd_rof_two_way   = ',use_lnd_rof_two_way
       write(iulog,*) '   heatflag              = ',heatflag
       write(iulog,*) '   barrier_timers        = ',barrier_timers
       write(iulog,*) '   RoutingMethod         = ',Tctl%RoutingMethod
       write(iulog,*) '   DLevelH2R             = ',Tctl%DLevelH2R
       write(iulog,*) '   DLevelR               = ',Tctl%DLevelR
       write(iulog,*) '   data_bgc_fluxes_to_ocean_flag = ',data_bgc_fluxes_to_ocean_flag
       if (nsrest == nsrStartup .and. finidat_rtm /= ' ') then
          write(iulog,*) '   MOSART initial data   = ',trim(finidat_rtm)
       end if
       if (inundflag) then
          write(iulog,*) ' '
          write(iulog,*) 'inundation settings:'
          write(iulog,*) '  OPT_inund          = ',Tctl%OPT_inund
          write(iulog,*) '  OPT_trueDW         = ',Tctl%OPT_trueDW
          write(iulog,*) '  OPT_calcNr         = ',Tctl%OPT_calcNr
          write(iulog,*) '  nr_max             = ',Tctl%nr_max
          write(iulog,*) '  nr_min             = ',Tctl%nr_min
          write(iulog,*) '  nr_uniform         = ',Tctl%nr_uniform
          write(iulog,*) '  rdepth_max         = ',Tctl%rdepth_max
          write(iulog,*) '  rdepth_min         = ',Tctl%rdepth_min
          write(iulog,*) '  rwidth_max         = ',Tctl%rwidth_max
          write(iulog,*) '  rwidth_min         = ',Tctl%rwidth_min
          write(iulog,*) '  rslp_assume        = ',Tctl%rslp_assume
          write(iulog,*) '  minL_tribRouting   = ',Tctl%minL_tribRouting
          write(iulog,*) '  OPT_elevProf       = ',Tctl%OPT_elevProf
          write(iulog,*) '  npt_elevProf       = ',Tctl%npt_elevProf
          write(iulog,*) '  threshold_slpRatio = ',Tctl%threshold_slpRatio
       endif
    endif

    rtm_active = do_rtm
    flood_active = do_rtmflood
    
    if (do_rtm) then
       if (frivinp_rtm == ' ') then
          call shr_sys_abort( subname//' ERROR: do_rtm TRUE, but frivinp_rtm NOT set' )
       else
          if (masterproc) then
             write(iulog,*) '   MOSART river data       = ',trim(frivinp_rtm)
          endif
       end if
    else
       if (masterproc) then
          write(iulog,*)'MOSART will not be active '
       endif
       RETURN
    end if

    if (wrmflag .and. inundflag) then
       write(iulog,*) subname,' MOSART wrmflag and inundflag both set to be true'
    endif

    if (coupling_period <= 0) then
       write(iulog,*) subname,' ERROR MOSART coupling_period invalid',coupling_period
       call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
    endif

    if (delt_mosart <= 0) then
       write(iulog,*) subname,' ERROR MOSART delt_mosart invalid',delt_mosart
       call shr_sys_abort( subname//' ERROR: delt_mosart invalid' )
    endif
       
    do i = 1, max_tapes
       if (rtmhist_nhtfrq(i) == 0) then
          rtmhist_mfilt(i) = 1
       else if (rtmhist_nhtfrq(i) < 0) then
          rtmhist_nhtfrq(i) = nint(-rtmhist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*coupling_period))
       endif
    end do

    !-------------------------------------------------------
    ! Initialize MOSART time manager 
    !-------------------------------------------------------

    ! Intiialize MOSART pio
    call ncd_pio_init()

    ! Obtain restart file if appropriate
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
        (nsrest == nsrBranch  )) then
       call RtmRestGetfile( file=fnamer, path=pnamer )
    endif       

    ! Initialize time manager
    if (nsrest == nsrStartup) then  
       call timemgr_init(dtime_in=coupling_period)
    else
       call RtmRestTimeManager(file=fnamer)
    end if

    !-------------------------------------------------------
    ! Initialize rtm_trstr
    !-------------------------------------------------------

    rtm_trstr = trim(rtm_tracers(1))
    do n = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(n))
    enddo
    if (masterproc) then
       write(iulog,*)'MOSART tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !-------------------------------------------------------
    ! Read input data (river direction file)
    !-------------------------------------------------------

    ! Useful constants and initial values
    deg2rad = SHR_CONST_PI / 180._r8

    call t_startf('mosarti_grid')

    call getfil(frivinp_rtm, locfn, 0 )
    if (masterproc) then
       write(iulog,*) 'Read in MOSART file name: ',trim(frivinp_rtm)
       call shr_sys_flush(iulog)
    endif

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqfdims(ncid, isgrid2d, rtmlon, rtmlat, rtmn)

    if (masterproc) then
       write(iulog,*) 'Values for rtmlon/rtmlat: ',rtmlon,rtmlat
       write(iulog,*) 'Successfully read MOSART dimensions'
       if (isgrid2d) then
        write(iulog,*) 'MOSART input is 2d'
       else
        write(iulog,*) 'MOSART input is 1d'
       endif
       call shr_sys_flush(iulog)
    endif

    ! Allocate variables
    if (isgrid2d) then
      allocate(rlonc(rtmlon), rlatc(rtmlat), &
               rlonw(rtmlon), rlone(rtmlon), &
               rlats(rtmlat), rlatn(rtmlat), &
               rtmCTL%rlon(rtmlon),          &
               rtmCTL%rlat(rtmlat),          &
               stat=ier)
    else
      allocate(rlonc(rtmlon), rlatc(rtmlon), &
               rtmCTL%rlon(rtmlon),          &
               rtmCTL%rlat(rtmlon),          &
               stat=ier)
    endif
    if (ier /= 0) then
       write(iulog,*) subname,' : Allocation ERROR for rlon'
       call shr_sys_abort(subname//' ERROR alloc for rlon')
    end if

    ! reading the routing parameters
    allocate ( &
              ID0_global(rtmlon*rtmlat), area_global(rtmlon*rtmlat), &
              dnID_global(rtmlon*rtmlat), nUp_global(rtmlon*rtmlat), nUp_dstrm_global(rtmlon*rtmlat),&
              stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname, ' : Allocation error for ID0_global'
       call shr_sys_abort(subname//' ERROR alloc for ID0')
    end if

    if (data_bgc_fluxes_to_ocean_flag) then
      allocate ( &
                DIN_global(rtmlon*rtmlat), DIP_global(rtmlon*rtmlat), &
                DON_global(rtmlon*rtmlat), DOP_global(rtmlon*rtmlat), &
                DOC_global(rtmlon*rtmlat), PP_global(rtmlon*rtmlat) , &
                DSi_global(rtmlon*rtmlat), POC_global(rtmlon*rtmlat), &
                PN_global(rtmlon*rtmlat) , DIC_global(rtmlon*rtmlat), &
                Fe_global(rtmlon*rtmlat) ,  &
                stat=ier)
      if (ier /= 0) then
         write(iulog,*) subname, ' : Allocation error for global BGC variables'
         call shr_sys_abort(subname//' ERROR alloc for BGC variables')
      end if
    end if

    allocate(tempr(rtmlon,rtmlat))  
    allocate(itempr(rtmlon,rtmlat))  

    call ncd_io(ncid=ncid, varname='longxy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART longitudes')
    if (masterproc) write(iulog,*) 'Read longxy ',minval(tempr),maxval(tempr)
    do i=1,rtmlon
       rtmCTL%rlon(i) = tempr(i,1)
       rlonc(i) = tempr(i,1)
    enddo
    if (masterproc) write(iulog,*) 'rlonc ',minval(rlonc),maxval(rlonc)

    call ncd_io(ncid=ncid, varname='latixy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART latitudes')
    if (masterproc) write(iulog,*) 'Read latixy ',minval(tempr),maxval(tempr)
    if (isgrid2d) then
      do j=1,rtmlat
         rtmCTL%rlat(j) = tempr(1,j)
         rlatc(j) = tempr(1,j)
      end do
    else
      do j=1,rtmlon
         rtmCTL%rlat(j) = tempr(j,1)
         rlatc(j) = tempr(j,1)
      end do
    endif
    if (masterproc) write(iulog,*) 'rlatc ',minval(rlatc),maxval(rlatc)

    call ncd_io(ncid=ncid, varname='area', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART area')
    if (masterproc) write(iulog,*) 'Read area ',minval(tempr),maxval(tempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       area_global(n) = tempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'area ',minval(tempr),maxval(tempr)

    call ncd_io(ncid=ncid, varname='ID', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART ID')
    if (masterproc) write(iulog,*) 'Read ID ',minval(itempr),maxval(itempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       ID0_global(n) = itempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'ID ',minval(itempr),maxval(itempr)

    call ncd_io(ncid=ncid, varname='dnID', flag='read', data=itempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART dnID')
    if (masterproc) write(iulog,*) 'Read dnID ',minval(itempr),maxval(itempr)
    do j=1,rtmlat
    do i=1,rtmlon
       n = (j-1)*rtmlon + i
       dnID_global(n) = itempr(i,j)
    end do
    end do
    if (masterproc) write(iulog,*) 'dnID ',minval(itempr),maxval(itempr)

    if (data_bgc_fluxes_to_ocean_flag) then

      call ncd_io(ncid=ncid, varname='Ld_DIN_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DIN')
      if (masterproc) write(iulog,*) 'Read DIN ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DIN_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DIN ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DIP_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DIP')
      if (masterproc) write(iulog,*) 'Read DIP ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DIP_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DIP ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DON_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DON')
      if (masterproc) write(iulog,*) 'Read DON ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DON_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DON ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DOP_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DOP')
      if (masterproc) write(iulog,*) 'Read DOP ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DOP_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DOP ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DOC_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DOC')
      if (masterproc) write(iulog,*) 'Read DOC ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DOC_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DOC ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_PP_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART PP')
      if (masterproc) write(iulog,*) 'Read PP ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         PP_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'PP ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DSi_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DSi')
      if (masterproc) write(iulog,*) 'Read DSi ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DSi_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DSi ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_POC_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART POC')
      if (masterproc) write(iulog,*) 'Read POC ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         POC_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'POC ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_PN_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART PN')
      if (masterproc) write(iulog,*) 'Read PN ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         PN_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'PN ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_DIC_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART DIC')
      if (masterproc) write(iulog,*) 'Read DIC ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         DIC_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'DIC ',minval(tempr),maxval(tempr)

      call ncd_io(ncid=ncid, varname='Ld_Fe_mosart_half', flag='read', data=tempr, readvar=found)
      if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read MOSART Fe')
      if (masterproc) write(iulog,*) 'Read Fe ',minval(tempr),maxval(tempr)
      do j=1,rtmlat
      do i=1,rtmlon
         n = (j-1)*rtmlon + i
         Fe_global(n) = tempr(i,j)
      end do
      end do
      if (masterproc) write(iulog,*) 'Fe ',minval(tempr),maxval(tempr)

    end if  !  data_bgc_fluxes_to_ocean_flag

    deallocate(tempr)
    deallocate(itempr)             

    call ncd_pio_closefile(ncid)

    !-------------------------------------------------------
    ! RESET dnID indices based on ID0
    ! rename the dnID values to be consistent with global grid indexing.
    ! where 1 = lower left of grid and rtmlon*rtmlat is upper right.
    ! ID0 is the "key", modify dnID based on that.  keep the IDkey around
    ! for as long as needed.  This is a key that translates the ID0 value
    ! to the gindex value.  compute the key, then apply the key to dnID_global.
    ! As part of this, check that each value of ID0 is unique and within
    ! the range of 1 to rtmlon*rtmlat.
    !-------------------------------------------------------

    allocate(IDkey(rtmlon*rtmlat))
    IDkey = 0
    do n=1,rtmlon*rtmlat
       if (ID0_global(n) < 0 .or. ID0_global(n) > rtmlon*rtmlat) then
          write(iulog,*) subname,' ERROR ID0 out of range',n,ID0_global(n)
          call shr_sys_abort(subname//' ERROR error ID0 out of range')
       endif
       if (IDkey(ID0_global(n)) /= 0) then
          write(iulog,*) subname,' ERROR ID0 value occurs twice',n,ID0_global(n)
          call shr_sys_abort(subname//' ERROR ID0 value occurs twice')
       endif
       IDkey(ID0_global(n)) = n
    enddo
    if (minval(IDkey) < 1) then
       write(iulog,*) subname,' ERROR IDkey incomplete'
       call shr_sys_abort(subname//' ERROR IDkey incomplete')
    endif
    do n=1,rtmlon*rtmlat
       if (dnID_global(n) > 0 .and. dnID_global(n) <= rtmlon*rtmlat) then
          if (IDkey(dnID_global(n)) > 0 .and. IDkey(dnID_global(n)) <= rtmlon*rtmlat) then
             dnID_global(n) = IDkey(dnID_global(n))
          else
             write(iulog,*) subname,' ERROR bad IDkey',n,dnID_global(n),IDkey(dnID_global(n))
             call shr_sys_abort(subname//' ERROR bad IDkey')
          endif
       endif
    enddo
    
    nUp_global = 0
    nUp_dstrm_global = 0
    do n=1,rtmlon*rtmlat
       if(dnID_global(n) > 0  .and. dnID_global(n) <= rtmlon*rtmlat) then
              nUp_global(dnID_global(n)) =  nUp_global(dnID_global(n)) + 1
       end if
    enddo
    do n=1,rtmlon*rtmlat
       if(dnID_global(n) > 0  .and. dnID_global(n) <= rtmlon*rtmlat) then
           nUp_dstrm_global(n) = nUp_global(dnID_global(n))
       end if
    enddo
    deallocate(ID0_global)

    !-------------------------------------------------------
    ! Derive gridbox edges
    !-------------------------------------------------------

    if (isgrid2d) then

      ! assuming equispaced grid, calculate edges from rtmlat/rtmlon
      ! w/o assuming a global grid
      edgen = maxval(rlatc) + 0.5*abs(rlatc(1) - rlatc(2))
      edges = minval(rlatc) - 0.5*abs(rlatc(1) - rlatc(2))
      edgee = maxval(rlonc) + 0.5*abs(rlonc(1) - rlonc(2))
      edgew = minval(rlonc) - 0.5*abs(rlonc(1) - rlonc(2))

      if ( edgen .ne.  90._r8 )then
         if (masterproc) write(iulog,*) 'Regional grid: edgen = ', edgen
      end if
      if ( edges .ne. -90._r8 )then
         if (masterproc) write(iulog,*) 'Regional grid: edges = ', edges
      end if
      if ( edgee .ne. 180._r8 )then
         if (masterproc) write(iulog,*) 'Regional grid: edgee = ', edgee
      end if
      if ( edgew .ne.-180._r8 )then
         if (masterproc) write(iulog,*) 'Regional grid: edgew = ', edgew
      end if

      ! Set edge latitudes (assumes latitudes are constant for a given longitude)
      rlats(:) = edges
      rlatn(:) = edgen
      do j = 2, rtmlat
         if (rlatc(2) > rlatc(1)) then ! South to North grid
            rlats(j)   = (rlatc(j-1) + rlatc(j)) / 2._r8
            rlatn(j-1) = rlats(j)
         else  ! North to South grid
            rlatn(j)   = (rlatc(j-1) + rlatc(j)) / 2._r8
            rlats(j-1) = rlatn(j)
         end if
      end do

      ! Set edge longitudes
      rlonw(:) = edgew
      rlone(:) = edgee
      dx = (edgee - edgew) / rtmlon
      do i = 2, rtmlon
         rlonw(i)   = rlonw(i) + (i-1)*dx
         rlone(i-1) = rlonw(i)
      end do

    endif

    call t_stopf ('mosarti_grid')

    !!call t_stopf('mosarti_grid')   (Repetitive, removed on 6-1-17. --Inund.)

    !-------------------------------------------------------
    ! Determine mosart ocn/land mask (global, all procs)
    !-------------------------------------------------------

    call t_startf('mosarti_decomp')

    allocate (gmask(rtmlon*rtmlat), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname, ' : Allocation ERROR for gmask'
       call shr_sys_abort(subname//' ERROR alloc for gmask')
    end if

    !  1=land, 
    !  2=ocean,
    !  3=ocean outlet from land

    gmask = 2    ! assume ocean point
    do n=1,rtmlon*rtmlat         ! mark all downstream points as outlet
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= rtmlon*rtmlat)) then
          gmask(nr) = 3          ! <- nr
       end if
    enddo
    do n=1,rtmlon*rtmlat         ! now mark all points with downstream points as land
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= rtmlon*rtmlat)) then
          gmask(n) = 1           ! <- n
       end if
    enddo

    !-------------------------------------------------------
    ! Compute total number of basins and runoff points
    !-------------------------------------------------------

    nbas = 0
    nrof = 0
    nout = 0
    nmos = 0
    ! over all MOSART grid cells
    do nr=1,rtmlon*rtmlat
       ! if ocean outlet from land
       if (gmask(nr) == 3) then
          nout = nout + 1
          nbas = nbas + 1
          nmos = nmos + 1
          nrof = nrof + 1
       ! if ocean
       elseif (gmask(nr) == 2) then
          nbas = nbas + 1
          nrof = nrof + 1
       ! if land
       elseif (gmask(nr) == 1) then
          nmos = nmos + 1
          nrof = nrof + 1
       endif
    enddo
    if (masterproc) then
       write(iulog,*) 'Number of outlet basins (if ocean) = ',nout
       write(iulog,*) 'Number of total  basins (either ocean or land-ocean outlet) = ',nbas
       write(iulog,*) 'Number of mosart points (either land or land-ocean outlet) = ',nmos
       write(iulog,*) 'Number of runoff points (all grids)= ',nrof
    endif

    !-------------------------------------------------------
    ! Compute river basins, actually compute ocean outlet gridcell
    !-------------------------------------------------------

    ! idxocn = final downstream cell, index is global 1d ocean gridcell
    ! nupstrm = number of source gridcells upstream including self

    allocate(idxocn(rtmlon*rtmlat),nupstrm(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' : Allocation ERROR for ',&
            'idxocn,nupstrm'
       call shr_sys_abort(subname//' ERROR alloc for idxocn nupstrm')
    end if

    call t_startf('mosarti_dec_basins')
    idxocn  = 0
    nupstrm = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             nupstrm(n) = nupstrm(n) + 1
             n = dnID_global(n)
             g = g + 1
          end do
          if (gmask(n) == 3) then           ! found ocean outlet 
             nupstrm(n) = nupstrm(n) + 1    ! one more land cell for n
             idxocn(nr) = n                 ! set ocean outlet or nr to n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) subname,' ERROR closed basin found', &
               g,nr,gmask(nr),dnID_global(nr), &
               n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR closed basin found')
          elseif (gmask(n) == 2) then
             write(iulog,*) subname,' ERROR found invalid ocean cell ',nr
             call shr_sys_abort(subname//' ERROR found invalid ocean cell')
          else 
             write(iulog,*) subname,' ERROR downstream cell is unknown', &
               g,nr,gmask(nr),dnID_global(nr), &
               n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR downstream cell is unknown')
          endif
       elseif (gmask(n) >= 2) then  ! ocean, give to self
          nupstrm(n) = nupstrm(n) + 1
          idxocn(nr) = n
       endif
    enddo
    call t_stopf('mosarti_dec_basins')

    ! check

    nbas_chk = 0
    nrof_chk = 0
    do nr=1,rtmlon*rtmlat
!      !if (masterproc) write(iulog,*) 'nupstrm check ',nr,gmask(nr),nupstrm(nr),idxocn(nr)
       if (gmask(nr) >= 2 .and. nupstrm(nr) > 0) then
          nbas_chk = nbas_chk + 1
          nrof_chk = nrof_chk + nupstrm(nr)
       endif
    enddo

    if (nbas_chk /= nbas .or. nrof_chk /= nrof) then
       write(iulog,*) subname,' ERROR nbas nrof check',nbas,nbas_chk,nrof,nrof_chk
       call shr_sys_abort(subname//' ERROR nbas nrof check')
    endif

    !-------------------------------------------------------
    !--- Now allocate those basins to pes
    !-------------------------------------------------------

    call t_startf('mosarti_dec_distr')

    !--- this is the heart of the decomp, need to set pocn and nop by the end of this
    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of mosart cells/pe 

    allocate(pocn(rtmlon*rtmlat),     &  !global mosart array
             nop(0:npes-1), &
             nba(0:npes-1))

    pocn = -99
    nop = 0
    nba = 0

    if (trim(decomp_option) == 'basin') then
       baspe = 0
       maxrtm = int(float(nrof)/float(npes)*0.445) + 1
       nloops = 3
       minbas = nrof
       do nl=1,nloops
          maxbas = minbas - 1
          minbas = maxval(nupstrm)/(2**nl)
          if (nl == nloops) minbas = min(minbas,1)
          do nr=1,rtmlon*rtmlat
             if (gmask(nr) >= 2 .and. nupstrm(nr) > 0 .and. nupstrm(nr) >= minbas .and. nupstrm(nr) <= maxbas) then
                ! Decomp options
                !   find min pe (implemented but scales poorly)
                !   use increasing thresholds (implemented, ok load balance for l2r or calc)
                !   distribute basins using above methods but work from max to min basin size
                !
                !--------------
                ! find min pe
                !             baspe = 0
                !             do n = 1,npes-1
                !                if (nop(n) < nop(baspe)) baspe = n
                !             enddo
                !--------------
                ! find next pe below maxrtm threshhold and increment
                do while (nop(baspe) > maxrtm)
                   baspe = baspe + 1
                   if (baspe > npes-1) then
                      baspe = 0
                      maxrtm = max(maxrtm*1.5, maxrtm+1.0)   ! 3 loop, .445 and 1.5 chosen carefully
                   endif
                enddo
                !--------------
                if (baspe > npes-1 .or. baspe < 0) then
                   write(iulog,*) 'ERROR in decomp for MOSART ',nr,npes,baspe
                   call shr_sys_abort('ERROR mosart decomp')
                endif
                nop(baspe) = nop(baspe) + nupstrm(nr)
                nba(baspe) = nba(baspe) + 1
                pocn(nr) = baspe
             endif
          enddo ! nr
       enddo ! nl

       ! set pocn for land cells, was set for ocean above
       do nr=1,rtmlon*rtmlat
          if (idxocn(nr) > 0) then
             pocn(nr) = pocn(idxocn(nr))
             if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
                write(iulog,*) subname,' ERROR pocn lnd setting ',&
                   nr,idxocn(nr),idxocn(idxocn(nr)),pocn(idxocn(nr)),pocn(nr),npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    elseif (trim(decomp_option) == '1d') then
       ! distribute active points in 1d fashion to pes
       ! baspe is the pe assignment
       ! maxrtm is the maximum number of points to assign to each pe
       baspe = 0
       maxrtm = (nrof-1)/npes + 1
       do nr=1,rtmlon*rtmlat
          if (gmask(nr) >= 1) then
             pocn(nr) = baspe
             nop(baspe) = nop(baspe) + 1
             if (nop(baspe) >= maxrtm) then
                baspe = (mod(baspe+1,npes))
                if (baspe < 0 .or. baspe > npes-1) then
                   write(iulog,*) subname,' ERROR basepe ',baspe,npes
                   call shr_sys_abort(subname//' ERROR pocn lnd')
                endif
             endif
          endif
       enddo

    elseif (trim(decomp_option) == 'roundrobin') then
       ! distribute active points in roundrobin fashion to pes
       ! baspe is the pe assignment
       ! maxrtm is the maximum number of points to assign to each pe
       baspe = 0
       do nr=1,rtmlon*rtmlat
          if (gmask(nr) >= 1) then
             pocn(nr) = baspe
             nop(baspe) = nop(baspe) + 1
             baspe = (mod(baspe+1,npes))
             if (baspe < 0 .or. baspe > npes-1) then
                write(iulog,*) subname,' ERROR basepe ',baspe,npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    else
       write(iulog,*) subname,' ERROR decomp option unknown ',trim(decomp_option)
       call shr_sys_abort(subname//' ERROR pocn lnd')
    endif  ! decomp_option

    if (masterproc) then
       write(iulog,*) 'MOSART cells and basins total  = ',nrof,nbas
       write(iulog,*) 'MOSART cells per basin avg/max = ',nrof/nbas,maxval(nupstrm)
       write(iulog,*) 'MOSART cells per pe    min/max = ',minval(nop),maxval(nop)
       write(iulog,*) 'MOSART basins per pe   min/max = ',minval(nba),maxval(nba)
    endif

    deallocate(nupstrm)

    !-------------------------------------------------------
    !--- Count and distribute cells to rglo2gdc
    !-------------------------------------------------------

    rtmCTL%numr   = 0
    rtmCTL%lnumr  = 0

    do n = 0,npes-1
       if (iam == n) then
          rtmCTL%begr  = rtmCTL%numr  + 1
       endif
       rtmCTL%numr  = rtmCTL%numr  + nop(n)
       if (iam == n) then
          rtmCTL%lnumr = rtmCTL%lnumr + nop(n)
          rtmCTL%endr  = rtmCTL%begr  + rtmCTL%lnumr  - 1
       endif
    enddo

    allocate(rglo2gdc(rtmlon*rtmlat), &  !global mosart array
             nrs(0:npes-1))
    nrs = 0
    rglo2gdc = 0

    ! nrs is begr on each pe
    nrs(0) = 1
    do n = 1,npes-1
       nrs(n) = nrs(n-1) + nop(n-1)
    enddo

    ! reuse nba for nop-like counter here
    ! pocn -99 is unused cell
    nba = 0
    do nr = 1,rtmlon*rtmlat
       if (pocn(nr) >= 0) then
          rglo2gdc(nr) = nrs(pocn(nr)) + nba(pocn(nr))
          nba(pocn(nr)) = nba(pocn(nr)) + 1          
       endif
    enddo
    do n = 0,npes-1
       if (nba(n) /= nop(n)) then
          write(iulog,*) subname,' ERROR mosart cell count ',n,nba(n),nop(n)
          call shr_sys_abort(subname//' ERROR mosart cell count')
       endif
    enddo

    deallocate(nop,nba,nrs)
    deallocate(pocn)
    call t_stopf('mosarti_dec_distr')

    !-------------------------------------------------------
    !--- adjust area estimation from DRT algorithm for those outlet grids
    !--- useful for grid-based representation only
    !--- need to compute areas where they are not defined in input file
    !-------------------------------------------------------

    do n=1,rtmlon*rtmlat
       if (area_global(n) <= 0._r8) then
          i = mod(n-1,rtmlon) + 1
          j = (n-1)/rtmlon + 1
          dx = (rlone(i) - rlonw(i)) * deg2rad
          dy = sin(rlatn(j)*deg2rad) - sin(rlats(j)*deg2rad)
          area_global(n) = abs(1.e6_r8 * dx*dy*re*re)
          if (masterproc .and. area_global(n) <= 0) then
             write(iulog,*) 'Warning! Zero area for unit ', n, area_global(n),dx,dy,re
          end if
       end if
    end do

    call t_stopf('mosarti_decomp')

    !-------------------------------------------------------
    !--- Write per-processor runoff bounds depending on dbug level
    !-------------------------------------------------------

    call t_startf('mosarti_print')

    call shr_sys_flush(iulog)
    if (masterproc) then
       write(iulog,*) 'total runoff cells numr  = ',rtmCTL%numr
    endif
    call shr_sys_flush(iulog)
    call mpi_barrier(mpicom_rof,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    elseif (dbug == 3) then
       npint = 1
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)
       if (iam == pid) then
          write(iulog,'(2a,i9,a,i9,a,i9,a,i9)') &
             'MOSART decomp info',' proc = ',iam, &
             ' begr = ',rtmCTL%begr,&
             ' endr = ',rtmCTL%endr, &
             ' numr = ',rtmCTL%lnumr
       endif
       call shr_sys_flush(iulog)
       call mpi_barrier(mpicom_rof,ier)
    enddo

    call t_stopf('mosarti_print')

    !-------------------------------------------------------
    ! Allocate local flux variables
    !-------------------------------------------------------

    call t_startf('mosarti_vars')

    allocate (evel    (rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              flow    (rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroup_lagi(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroup_lagf(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erowm_regi(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erowm_regf(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              eroutup_avg(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              erlat_avg(rtmCTL%begr:rtmCTL%endr,nt_rtm), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' Allocation ERROR for flow'
       call shr_sys_abort(subname//' Allocationt ERROR flow')
    end if
    flow(:,:)        = 0._r8
    eroup_lagi(:,:)  = 0._r8
    eroup_lagf(:,:)  = 0._r8
    erowm_regi(:,:)  = 0._r8
    erowm_regf(:,:)  = 0._r8
    eroutup_avg(:,:) = 0._r8
    erlat_avg(:,:)   = 0._r8

    if (inundflag) then
       ! If inundation scheme is turned on :
       !if ( Tctl%OPT_inund .eq. 1 ) then
       allocate ( fa_fp_cplPeriod(rtmCTL%begr : rtmCTL%endr), stat=ier )
       if (ier /= 0) then
          write(iulog, *) subname,' Allocation ERROR for "fa_fp_cplPeriod".'
          call shr_sys_abort(subname//' Allocationt ERROR for "fa_fp_cplPeriod".')
       end if
       fa_fp_cplPeriod(:) = 0._r8
       !end if

       allocate ( totalVelo_down(nt_rtm), &
            totalVelo_up(nt_rtm), &
            chnlNum_down(nt_rtm), &
            chnlNum_up(nt_rtm), &
            stat=ier )
       if (ier /= 0) then
          write(iulog, *) subname,' Allocation ERROR for velocity variables.'
          call shr_sys_abort(subname//' Allocationt ERROR for velocity variables.')
       end if
       totalVelo_down(:) = 0._r8
       totalVelo_up(:) = 0._r8
       chnlNum_down(:) = 0
       chnlNum_up(:) = 0
    end if

    !-------------------------------------------------------
    ! Allocate runoff datatype 
    !-------------------------------------------------------

    call RunoffInit(rtmCTL%begr, rtmCTL%endr, rtmCTL%numr)

    !-------------------------------------------------------
    ! Initialize mosart flood - rtmCTL%fthresh and evel
    !-------------------------------------------------------

    if (do_rtmflood) then
       write(iulog,*) subname,' Flood not validated in this version, abort'
       call shr_sys_abort(subname//' Flood feature unavailable')
       call RtmFloodInit (frivinp_rtm, rtmCTL%begr, rtmCTL%endr, rtmCTL%fthresh, evel)
    else
       effvel(:) = effvel0  ! downstream velocity (m/s)
       rtmCTL%fthresh(:) = abs(spval)
       do nt = 1,nt_rtm
          do nr = rtmCTL%begr,rtmCTL%endr
             evel(nr,nt) = effvel(nt)
          enddo
       enddo
    end if

    !-------------------------------------------------------
    ! Initialize runoff data type
    !-------------------------------------------------------

    allocate(rgdc2glo(rtmCTL%numr), stat=ier)
    if (ier /= 0) then
       write(iulog,*) subname,' ERROR allocation of rgdc2glo'
       call shr_sys_abort(subname//' ERROR allocate of rgdc2glo')
    end if

    ! Set map from local to global index space
    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = rglo2gdc(n)
       if (nr > 0) then
          numr = numr + 1
          rgdc2glo(nr) = n
       endif
    end do
    end do
    if (numr /= rtmCTL%numr) then
       write(iulog,*) subname,'ERROR numr and rtmCTL%numr are different ',numr,rtmCTL%numr
       call shr_sys_abort(subname//' ERROR numr')
    endif

    ! Determine runoff datatype variables
    lrtmarea = 0.0_r8
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       rtmCTL%gindex(nr) = rgdc2glo(nr)
       rtmCTL%mask(nr) = gmask(rgdc2glo(nr))
       n = rgdc2glo(nr)
       i = mod(n-1,rtmlon) + 1
       j = (n-1)/rtmlon + 1
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(iulog,*) subname,' ERROR gdc2glo, nr,ng= ',nr,n
          call shr_sys_abort(subname//' ERROR gdc2glo values')
       endif
       rtmCTL%lonc(nr) = rtmCTL%rlon(i)
       if (isgrid2d) then 
          rtmCTL%latc(nr) = rtmCTL%rlat(j)
       else 
          rtmCTL%latc(nr) = rtmCTL%rlat(i)
       endif

       rtmCTL%outletg(nr) = idxocn(n)
       rtmCTL%area(nr) = area_global(n)
       lrtmarea = lrtmarea + rtmCTL%area(nr)
       if (dnID_global(n) <= 0) then
          rtmCTL%dsig(nr) = 0
       else
          if (rglo2gdc(dnID_global(n)) == 0) then
             write(iulog,*) subname,' ERROR glo2gdc dnID_global ',&
                  nr,n,dnID_global(n),rglo2gdc(dnID_global(n))
             call shr_sys_abort(subname//' ERROT glo2gdc dnID_global')
          endif
          cnt = cnt + 1
          rtmCTL%dsig(nr) = dnID_global(n)
          rtmCTL%iDown(nr) = rglo2gdc(dnID_global(n))
       endif
       rtmCTL%nUp_dstrm(nr) = nUp_dstrm_global(n)
       if (data_bgc_fluxes_to_ocean_flag) then
         rtmCTL%concDIN(nr) = DIN_global(n)
         rtmCTL%concDIP(nr) = DIP_global(n)
         rtmCTL%concDON(nr) = DON_global(n)
         rtmCTL%concDOP(nr) = DOP_global(n)
         rtmCTL%concDOC(nr) = DOC_global(n)
         rtmCTL%concPP(nr)  = PP_global(n)
         rtmCTL%concDSi(nr) = DSi_global(n)
         rtmCTL%concPOC(nr) = POC_global(n)
         rtmCTL%concPN(nr)  = PN_global(n)
         rtmCTL%concDIC(nr) = DIC_global(n)
         rtmCTL%concFe(nr)  = Fe_global(n)
       end if
    enddo
    
    rtmCTL%iUp = 0
    rtmCTL%nUp = 0
    do nr = rtmCTL%begr,rtmCTL%endr
        do n = rtmCTL%begr,rtmCTL%endr
            if(rtmCTL%iDown(n) == nr) then
                rtmCTL%nUp(nr) = rtmCTL%nUp(nr) + 1  ! initial value of rtmCTL%nUp is 0
                rtmCTL%iUp(nr,rtmCTL%nUp(nr)) = n
            end if
        end do
    enddo
    
    deallocate(gmask)
    deallocate(rglo2gdc)
    deallocate(rgdc2glo)
    deallocate(dnID_global,area_global)
    deallocate(nUp_global,nUp_dstrm_global)
    deallocate(idxocn)
    if (data_bgc_fluxes_to_ocean_flag) then
      deallocate (DIN_global,DIP_global)
      deallocate (DON_global,DOP_global)
      deallocate (DOC_global,PP_global)
      deallocate (DSi_global,POC_global)
      deallocate (PN_global ,DIC_global)
      deallocate (Fe_global)
    end if
    call shr_mpi_sum(lrtmarea,rtmCTL%totarea,mpicom_rof,'mosart totarea',all=.true.)
    if (masterproc) write(iulog,*) subname,'  earth area ',4.0_r8*shr_const_pi*1.0e6_r8*re*re
    if (masterproc) write(iulog,*) subname,' MOSART area ',rtmCTL%totarea
    if (minval(rtmCTL%mask) < 1) then
       write(iulog,*) subname,'ERROR rtmCTL mask lt 1 ',minval(rtmCTL%mask),maxval(rtmCTL%mask)
       call shr_sys_abort(subname//' ERROR rtmCTL mask')
    endif

    !-------------------------------------------------------
    ! Compute Sparse Matrix for downstream advection
    !-------------------------------------------------------

    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat
    allocate(gindex(lsize))
    do nr = rtmCTL%begr,rtmCTL%endr
       gindex(nr-rtmCTL%begr+1) = rtmCTL%gindex(nr)
    enddo
    call mct_gsMap_init( gsMap_r, gindex, mpicom_rof, ROFID, lsize, gsize )
    deallocate(gindex)

    if (smat_option == 'opt') then
       ! distributed smat initialization
       ! mct_sMat_init must be given the number of rows and columns that
       ! would be in the full matrix.  Nrows= size of output vector=nb.
       ! Ncols = size of input vector = na.

       cnt = 0
       do nr=rtmCTL%begr,rtmCTL%endr
          if(rtmCTL%dsig(nr) > 0) cnt = cnt + 1
       enddo

       call mct_sMat_init(sMat, gsize, gsize, cnt)
       igrow = mct_sMat_indexIA(sMat,'grow')
       igcol = mct_sMat_indexIA(sMat,'gcol')
       iwgt  = mct_sMat_indexRA(sMat,'weight')
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%dsig(nr) > 0) then
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%dsig(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          endif
       enddo

       call mct_sMatP_Init(sMatP_upstrm, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%dsig(nr) > 0) then
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%dsig(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%gindex(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          endif
       enddo

       call mct_sMatP_Init(sMatP_dnstrm, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

    elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then

       ! root initialization

       call mct_aVect_init(avtmp,rList='f1:f2',lsize=lsize)
       call mct_aVect_zero(avtmp)
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avtmp%rAttr(1,cnt) = rtmCTL%gindex(nr)
          avtmp%rAttr(2,cnt) = rtmCTL%dsig(nr)
       enddo
       call mct_avect_gather(avtmp,avtmpG,gsmap_r,mastertask,mpicom_rof)
       call mct_avect_clean(avtmp)
       if (masterproc) then
          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
             endif
          enddo

          call mct_sMat_init(sMat, gsize, gsize, cnt)
          igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
          igrow = mct_sMat_indexIA(sMat,'grow')
          iwgt  = mct_sMat_indexRA(sMat,'weight')

          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(2,n)
             endif
          enddo
       else
          call mct_sMat_init(sMat,1,1,1)
       endif

       call mct_sMatP_Init(sMatP_upstrm, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

       if (masterproc) then
          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(2,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(1,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             endif
          enddo
          call mct_avect_clean(avtmpG)
       endif

       call mct_sMatP_Init(sMatP_dnstrm, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

    else

       write(iulog,*) trim(subname),' MOSART ERROR: invalid smat_option '//trim(smat_option)
       call shr_sys_abort(trim(subname)//' ERROR invald smat option')

    endif

    ! initialize the AVs to go with sMatP
    write(rList,'(a,i3.3)') 'tr',1
    do nt = 2,nt_rtm
       write(rList,'(a,i3.3)') trim(rList)//':tr',nt
    enddo
    write(rList,'(a,i3.3)') trim(rList)//':tr',nt_rtm+1
    if (masterproc) write(iulog,*) trim(subname),' MOSART initialize avect ',trim(rList)
    call mct_aVect_init(avsrc_upstrm,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avdst_upstrm,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avsrc_dnstrm,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avdst_dnstrm,rList=rList,lsize=rtmCTL%lnumr)

    lsize = mct_smat_gNumEl(sMatP_upstrm%Matrix,mpicom_rof)
    if (masterproc) write(iulog,*) subname," Done initializing SmatP_upstrm, nElements = ",lsize
    lsize = mct_smat_gNumEl(sMatP_dnstrm%Matrix,mpicom_rof)
    if (masterproc) write(iulog,*) subname," Done initializing SmatP_dnstrm, nElements = ",lsize

    ! keep only sMatP
    call mct_sMat_clean(sMat)

    !-------------------------------------------------------
    ! Compute Sparse Matrix for direct to outlet transfer
    ! reuse gsmap_r
    !-------------------------------------------------------

    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat

    if (smat_option == 'opt') then
       ! distributed smat initialization
       ! mct_sMat_init must be given the number of rows and columns that
       ! would be in the full matrix.  Nrows= size of output vector=nb.
       ! Ncols = size of input vector = na.

       cnt = rtmCTL%endr - rtmCTL%begr + 1

       call mct_sMat_init(sMat, gsize, gsize, cnt)
       igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
       igrow = mct_sMat_indexIA(sMat,'grow')
       iwgt  = mct_sMat_indexRA(sMat,'weight')
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%outletg(nr) > 0) then
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%outletg(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          else
             cnt = cnt + 1
             sMat%data%iAttr(igcol,cnt) = rtmCTL%gindex(nr)
             sMat%data%iAttr(igrow,cnt) = rtmCTL%gindex(nr)
             sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
          endif
       enddo
       if (cnt /= rtmCTL%endr - rtmCTL%begr + 1) then
          write(iulog,*) trim(subname),' MOSART ERROR: smat cnt1 ',cnt,rtmCTL%endr-rtmCTL%begr+1
          call shr_sys_abort(trim(subname)//' ERROR smat cnt1')
       endif

       call mct_sMatP_Init(sMatP_direct, sMat, gsMap_r, gsMap_r, 0, mpicom_rof, ROFID)

    elseif (smat_option == 'Xonly' .or. smat_option == 'Yonly') then

       ! root initialization

       call mct_aVect_init(avtmp,rList='f1:f2',lsize=lsize)
       call mct_aVect_zero(avtmp)
       cnt = 0
       do nr = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avtmp%rAttr(1,cnt) = rtmCTL%gindex(nr)
          avtmp%rAttr(2,cnt) = rtmCTL%outletg(nr)
       enddo
       call mct_avect_gather(avtmp,avtmpG,gsmap_r,mastertask,mpicom_rof)
       if (masterproc) then

          cnt = rtmlon*rtmlat

          call mct_sMat_init(sMat, gsize, gsize, cnt)
          igcol = mct_sMat_indexIA(sMat,'gcol')   ! src
          igrow = mct_sMat_indexIA(sMat,'grow')
          iwgt  = mct_sMat_indexRA(sMat,'weight')

          cnt = 0
          do n = 1,rtmlon*rtmlat
             if (avtmpG%rAttr(2,n) > 0) then
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(2,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             else
                cnt = cnt + 1
                sMat%data%iAttr(igcol,cnt) = avtmpG%rAttr(1,n)
                sMat%data%iAttr(igrow,cnt) = avtmpG%rAttr(1,n)
                sMat%data%rAttr(iwgt ,cnt) = 1.0_r8
             endif
          enddo
          if (cnt /= rtmlon*rtmlat) then
             write(iulog,*) trim(subname),' MOSART ERROR: smat cnt2 ',cnt,rtmlon*rtmlat
             call shr_sys_abort(trim(subname)//' ERROR smat cnt2')
          endif
          call mct_avect_clean(avtmpG)
       else
          call mct_sMat_init(sMat,1,1,1)
       endif
       call mct_avect_clean(avtmp)

       call mct_sMatP_Init(sMatP_direct, sMat, gsMap_r, gsMap_r, smat_option, 0, mpicom_rof, ROFID)

    else

       write(iulog,*) trim(subname),' MOSART ERROR: invalid smat_option '//trim(smat_option)
       call shr_sys_abort(trim(subname)//' ERROR invald smat option')

    endif

    ! initialize the AVs to go with sMatP
    write(rList,'(a,i3.3)') 'tr',1
    do nt = 2,nt_rtm
       write(rList,'(a,i3.3)') trim(rList)//':tr',nt
    enddo
    write(rList,'(a,i3.3)') trim(rList)//':tr',nt_rtm+1
    if (masterproc) write(iulog,*) trim(subname),' MOSART initialize avect ',trim(rList)
    call mct_aVect_init(avsrc_direct,rList=rList,lsize=rtmCTL%lnumr)
    call mct_aVect_init(avdst_direct,rList=rList,lsize=rtmCTL%lnumr)
!    write(iulog,*) subname,' avsrc_direct lsize = ',iam,mct_aVect_lsize(avsrc_direct)
!    write(iulog,*) subname,' avdst_direct lsize = ',iam,mct_aVect_lsize(avdst_direct)

    lsize = mct_smat_gNumEl(sMatP_direct%Matrix,mpicom_rof)
    if (masterproc) write(iulog,*) subname," Done initializing SmatP_direct, nElements = ",lsize

    ! keep only sMatP
    call mct_sMat_clean(sMat)

    !-------------------------------------------------------
    ! Compute timestep and subcycling number
    !-------------------------------------------------------

! tcraig, old code based on cfl
!    dtover = 0._r8
!    dtovermax = 0._r8
!!    write(iulog,*) "tcx ddist ",minval(ddist),maxval(ddist)
!!    write(iulog,*) "tcx evel  ",minval(evel),maxval(evel)
!    do nt=1,nt_rtm
!       do nr=rtmCTL%begr,rtmCTL%endr
!          if (ddist(nr) /= 0._r8) then
!             dtover = evel(nr,nt)/ddist(nr)
!          else
!             dtover = 0._r8
!          endif
!          dtovermax = max(dtovermax,dtover)
!       enddo
!    enddo
!    dtover = dtovermax
!    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom_rof,ier)
!
!    write(iulog,*) "tcx dtover ",dtover,dtovermax
!
!    if (dtovermax > 0._r8) then
!       delt_mosart = (1.0_r8/dtovermax)*cfl_scale
!    else
!       write(iulog,*) subname,' ERROR in delt_mosart ',delt_mosart,dtover
!       call shr_sys_abort(subname//' ERROR delt_mosart')
!    endif
!
!    if (masterproc) write(iulog,*) 'mosart max timestep = ',delt_mosart,' (sec) for cfl_scale = ',cfl_scale
!    if (masterproc) call shr_sys_flush(iulog)
!
!    delt_mosart = 600._r8  ! here set the time interval for routing as 10 mins, which is sufficient for 1/8th degree resolution and coarser.
!    if (masterproc) write(iulog,*) 'mosart max timestep hardwired to = ',delt_mosart
!    if (masterproc) call shr_sys_flush(iulog)

    call t_stopf('mosarti_vars')

    !-------------------------------------------------------
    ! Initialize mosart
    !-------------------------------------------------------

    call t_startf('mosarti_mosart_init')

    !=== initialize MOSART related variables
!    if (masterproc) write(iulog,*) ' call mosart_init'
!    if (masterproc) call shr_sys_flush(iulog)
    call MOSART_init()

    call t_stopf('mosarti_mosart_init')

    if (wrmflag) then
       call t_startf('mosarti_wrm_init')
       if (wrmflag) then
          call WRM_init()
       endif
       call t_startf('mosarti_wrm_init')
    end if
    if (wrmflag .and. heatflag .and. rstraflag) then
       call regeom                    
    end if  

    !-------------------------------------------------------
    ! Read restart/initial info
    !-------------------------------------------------------

    call t_startf('mosarti_restart')

!    if (masterproc) write(iulog,*) ' call RtmRestFileRead'
!    if (masterproc) call shr_sys_flush(iulog)

    ! The call below opens and closes the file
    if ((nsrest == nsrStartup .and. finidat_rtm /= ' ') .or. &
        (nsrest == nsrContinue) .or. & 
        (nsrest == nsrBranch  )) then
       call RtmRestFileRead( file=fnamer )
       !write(iulog,*) ' MOSART init file is read'
       TRunoff%wh   = rtmCTL%wh
       TRunoff%wt   = rtmCTL%wt
       TRunoff%wr   = rtmCTL%wr
       TRunoff%mr   = rtmCTL%mr
       TRunoff%yr   = rtmCTL%yr
       TRunoff%pr   = rtmCTL%pr
       TRunoff%rr   = rtmCTL%rr
       TRunoff%erout= rtmCTL%erout

       if (inundflag) then
          ! If inundation scheme is turned on :
          if ( Tctl%OPT_inund .eq. 1 ) then
             !TRunoff%wf_ini(:) = rtmCTL%wf(:, 1)
             ! Innudation floodplain water volume (m3)
             TRunoff%wf_ini(:) = rtmCTL%inundwf(:)
             ! Inundation floodplain water depth (m)
             TRunoff%hf_ini(:) = rtmCTL%inundhf(:)
             ! Inundation floodplain fraction
             TRunoff%ff_ini(:) = rtmCTL%inundff(:)
             ! Inundation flooded fraction, including river and floodplain
             TRunoff%ffunit_ini(:) = rtmCTL%inundffunit(:)
          end if
       end if
       if (wrmflag) then 
           call WRM_computeRelease()
        endif

    else
!       do nt = 1,nt_rtm
!       do nr = rtmCTL%begr,rtmCTL%endr
!          TRunoff%wh(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 1.e-10_r8
!          TRunoff%wt(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 1.e-8_r8
!          TRunoff%wr(nr,nt) = rtmCTL%area(nr) * river_depth_minimum * 10._r8
!       enddo
!       enddo

        if (wrmflag) then
           call WRM_computeRelease()
        endif

    endif

    !----------------------------------  
    ! Set inital condition (water volumes and streamflows) :
    !----------------------------------  

    if (inundflag .and. nsrest .eq. nsrStartup) then

      do nr = rtmCTL%begr, rtmCTL%endr
        if ( rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
          ! Water volume (water depth) over hillslopes (m) :
          rtmCTL%wh(nr, 1) = 0.001_r8   ! Assumed water depth.

          ! Water volume in the subnetwork (tributary channels) (m^3) :
          rtmCTL%wt(nr, 1) = TUnit%tlen( nr ) * TUnit%twidth( nr ) * 1._r8   ! Assumed water depth.

          ! Water depth in the main channel (m) :
          wd_chnl = TUnit%rdepth( nr ) * 0.9_r8   ! Assumed water depth (90% full).

          ! Water volume in the main channel (m^3) :
          rtmCTL%wr(nr, 1) = TUnit%rlen( nr ) * TUnit%rwidth( nr ) * wd_chnl

          ! Hydraulic radius (m) :
          hydrR = ( TUnit%rwidth( nr ) * wd_chnl ) / ( TUnit%rwidth( nr ) + 2._r8 * wd_chnl )

          ! Estimate channel flow velocity with Manning equation (m/s) :
          v_chnl = ManningEq ( TUnit%rslp( nr ), TUnit%nr( nr ), hydrR )

          ! Main channel outflow (outflow is negative) (m^3/s) :
          rtmCTL%erout( nr, 1 ) = - v_chnl * TUnit%rwidth( nr ) * wd_chnl

          if ( Tctl%OPT_inund .eq. 1 ) then
            rtmCTL%inundwf(nr) = 0._r8
            rtmCTL%inundhf(nr) = 0._r8
            rtmCTL%inundff(nr) = 0._r8
            rtmCTL%inundffunit(nr) = 0._r8                                                           
          end if

        else

          rtmCTL%wh(nr, 1) = 0._r8
          rtmCTL%wt(nr, 1) = 0._r8
          rtmCTL%wr(nr, 1) = 0._r8
          rtmCTL%erout(nr, 1) = 0._r8
          if ( Tctl%OPT_inund .eq. 1 ) then
            rtmCTL%inundwf(nr) = 0._r8
            rtmCTL%inundhf(nr) = 0._r8
            rtmCTL%inundff(nr) = 0._r8
            rtmCTL%inundffunit(nr) = 0._r8
          end if

        end if
      end do

      ! For tracer 2 :
      rtmCTL%wh(:, 2) = 0._r8
      rtmCTL%wt(:, 2) = 0._r8
      rtmCTL%wr(:, 2) = 0._r8
      rtmCTL%erout(:, 2) = 0._r8

      TRunoff%wh   = rtmCTL%wh
      TRunoff%wt   = rtmCTL%wt
      TRunoff%wr   = rtmCTL%wr
      TRunoff%erout= rtmCTL%erout

      if ( Tctl%OPT_inund .eq. 1 ) then
        TRunoff%wf_ini(:) = rtmCTL%inundwf(:)
        TRunoff%hf_ini(:) = rtmCTL%inundhf(:)
        TRunoff%ff_ini(:) = rtmCTL%inundff(:)
        TRunoff%ffunit_ini(:) = rtmCTL%inundffunit(:)    
      end if

    end if
!#endif


 do nt = 1,nt_rtm
      do nr = rtmCTL%begr,rtmCTL%endr
      
       call UpdateState_hillslope(nr,nt)
       call UpdateState_subnetwork(nr,nt)   
                                          
        rtmCTL%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + TRunoff%wh(nr,nt)*rtmCTL%area(nr)*TUnit%frac(nr))  ! times "TUnit%frac( nr )" or not ?
        if (inundflag .and. nt == 1) then  
           rtmCTL%volr(nr,nt) = rtmCTL%volr(nr,nt) + TRunoff%wf_ini( nr )
        else        
           call UpdateState_mainchannel(nr,nt)
        endif
        
      enddo
 enddo

    call t_stopf('mosarti_restart')

    !-------------------------------------------------------
    ! Initialize mosart history handler and fields
    !-------------------------------------------------------

    call t_startf('mosarti_histinit')

!    if (masterproc) write(iulog,*) ' call RtmHistFldsInit'
!    if (masterproc) call shr_sys_flush(iulog)

    call RtmHistFldsInit()
    if (nsrest==nsrStartup .or. nsrest==nsrBranch) then
       call RtmHistHtapesBuild()
    end if
    call RtmHistFldsSet()
    delt_save = 0.0

    if (masterproc) write(iulog,*) subname,' done'
    if (masterproc) call shr_sys_flush(iulog)

    call t_stopf('mosarti_histinit')

  end subroutine Rtmini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmrun
!
! !INTERFACE:
  subroutine Rtmrun(rstwr,nlend,rdate)
!
! !DESCRIPTION:
! River routing model
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    logical ,         intent(in) :: rstwr          ! true => write restart file this step)
    logical ,         intent(in) :: nlend          ! true => end of run on this step
    character(len=*), intent(in) :: rdate          ! restart file time stamp for name
!
! !CALLED FROM:
! subroutine RtmMap in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i, j, n, nr, ns, nt, n2, nf, idam ! indices
    integer, parameter :: budget_terms_total = 80
    logical  :: output_all_budget_terms = .false.   ! output flag
    real(r8) :: budget_terms (budget_terms_total,nt_rtm)    ! local budget sums
    real(r8) :: budget_global(budget_terms_total,nt_rtm)    ! global budget sums

    real(r8) :: budget_glb_inund(budget_terms_total, nt_rtm)! Global values.
    real(r8) :: budget_diff                 ! Percentage difference between water volume change and flux difference (%).
    character( len = 9 ) :: tracerID

    real(r8),save :: budget_accum(nt_rtm)   ! BUDGET accumulator over run
    integer ,save :: budget_accum_cnt       ! counter for budget_accum
    real(r8) :: budget_input, budget_output, budget_volume, budget_total, &
                budget_other
    logical  :: budget_check, budget_write  ! do global budget check

    ! BUDGET term ids
    ! budget computed in m3 over each coupling period
    ! use delt_coupling to convert from rates to volumes
    ! General equation is as follows
    !   Vf = Vi + input - output + other
    !     or
    !   Vf-Vi - input + output - other == 0 for conservation
    ! budget_volume = Vf-Vi
    ! budget_input = input
    ! budget_output = output
    ! budget_other = other terms
    ! budget_total = budget_volume - budget_input + budget_output - budget_other
    ! bv_ generally accumulates a volume (m3)
    ! br_ generally accumulates a rate (m3/s)

    ! Storage/Volume TERMS (volumes, m3)
    integer,parameter :: bv_volt_i = 1  ! initial total volume
    integer,parameter :: bv_volt_f = 2  ! final   total volume
    integer,parameter :: bv_wt_i   = 3  ! initial wt volume
    integer,parameter :: bv_wt_f   = 4  ! final   wt volume
    integer,parameter :: bv_wr_i   = 5  ! initial wr volume
    integer,parameter :: bv_wr_f   = 6  ! final   wr volume
    integer,parameter :: bv_wh_i   = 7  ! initial wh volume
    integer,parameter :: bv_wh_f   = 8  ! final   wh volume
    integer,parameter :: bv_dstor_i= 9  ! initial dam storage
    integer,parameter :: bv_dstor_f= 10 ! final   dam storage

    integer,parameter :: bv_fp_i   = 11 ! Initial water volume over floodplains.
    integer,parameter :: bv_fp_f   = 12 ! Final water volume over floodplains.

    integer,parameter :: bmask_lndErr   = 15 ! Tunit%mask() is land, but rtmCTL%mask() is not land.
    integer,parameter :: bmask_outErr   = 16 ! Tunit%mask() is basin outlet, but rtmCTL%mask() is not basin outlet.
    integer,parameter :: bmask_ocnErr   = 17 ! Tunit%mask() is ocean, but rtmCTL%mask() is not ocean.

    integer,parameter :: bmask_uLndrLnd = 31 ! Tunit%mask() is land, rtmCTL%mask() is land.
    integer,parameter :: bmask_uLndrOut = 32 ! Tunit%mask() is land, rtmCTL%mask() is outlet.
    integer,parameter :: bmask_uLndrOcn = 33 ! Tunit%mask() is land, rtmCTL%mask() is ocean.

    integer,parameter :: bmask_uOutrLnd = 34 ! Tunit%mask() is outlet, rtmCTL%mask() is land.
    integer,parameter :: bmask_uOutrOut = 35 ! Tunit%mask() is outlet, rtmCTL%mask() is outlet.
    integer,parameter :: bmask_uOutrOcn = 36 ! Tunit%mask() is outlet, rtmCTL%mask() is ocean.

    integer,parameter :: bmask_uOcnrLnd = 37 ! Tunit%mask() is ocean, rtmCTL%mask() is land.
    integer,parameter :: bmask_uOcnrOut = 38 ! Tunit%mask() is ocean, rtmCTL%mask() is outlet.
    integer,parameter :: bmask_uOcnrOcn = 39 ! Tunit%mask() is ocean, rtmCTL%mask() is ocean.

    ! Input TERMS (rates, m3/s)
    integer,parameter :: br_qsur   = 20 ! input qsur
    integer,parameter :: br_qsub   = 21 ! input qsub
    integer,parameter :: br_qgwl   = 22 ! input qgwl
    integer,parameter :: br_qdto   = 23 ! input qdto
    integer,parameter :: br_qdem   = 24 ! input qdem

    ! Output TERMS (rates m3/s or volumes m3)
    integer,parameter :: br_ocnout = 40 ! runoff output to ocean
    integer,parameter :: br_lndout = 41 ! runoff output on non ocean points
    integer,parameter :: br_flood  = 42 ! flood term back to land
    integer,parameter :: br_direct = 43 ! direct output term
    integer,parameter :: bv_dsupp_i= 44 ! initial dam supply
    integer,parameter :: bv_dsupp_f= 45 ! final   dam supply

    integer,parameter :: bv_chnl2fp       = 46 ! Volume of flows from main channels to floodplains (m^3).
    integer,parameter :: bv_fp2chnl       = 47 ! Volume of flows from floodplains to main channels (m^3).
    integer,parameter :: br_landOutflow   = 48 ! Total streamflow (flow rate) from land to oceans (m^3/s).

    integer,parameter :: bi_landArea      = 50 ! Total land area (m^2).
    integer,parameter :: bi_floodedArea   = 51 ! Total flooded area (including channel area) (m^2).
    integer,parameter :: bi_mainChnlArea  = 52 ! Total channel surface area (m^2).

    integer,parameter :: bVelo_downward   = 55 ! Sum of all downward flow velocities (m/s).
    integer,parameter :: bVelo_downChnlNo = 56 ! Total number of main channels with downward flow velocities (dimensionless).
    integer,parameter :: bVelo_upward     = 57 ! Sum of all upward flow velocities (is negative) (m/s).
    integer,parameter :: bVelo_upChnlNo   = 58 ! Total number of channels with upward flow velocities (dimensionless).

    ! Other Diagnostic TERMS (rates, m3/s)
    integer,parameter :: br_erolpo = 60 ! erout lag ocn previous
    integer,parameter :: br_erolco = 61 ! erout lag ocn current
    integer,parameter :: br_erorpo = 62 ! erout lag ocn previous   (for WRM module. --Inund.)
    integer,parameter :: br_erorco = 63 ! erout lag ocn current   (for WRM module. --Inund.)
    integer,parameter :: br_eroutup= 64 ! erout upstream average
    integer,parameter :: br_erolpn = 65 ! erout lag non-ocn previous
    integer,parameter :: br_erolcn = 66 ! erout lag non-ocn current
    integer,parameter :: br_erorpn = 67 ! erout lag non-ocn previous   (for WRM module. --Inund.)
    integer,parameter :: br_erorcn = 68 ! erout lag non-ocn current   (for WRM module. --Inund.)
    integer,parameter :: br_erlat  = 69 ! erlateral 

    ! Accumuluation TERMS
    integer,parameter :: bv_naccum = 80 ! accumulated net budget

    !   volume = 2 - 1 + bv_dstor_f - bv_dstor_i
    !   input  = br_qsur + br_qsub + br_qgwl + br_qdto + br_qdem
    !   output = br_ocnout + br_flood + br_direct + 42
    !   total  = volume - input + output
    !   erlag  = br_erolpn - br_erolcn

    real(r8) :: volr_init                   ! temporary storage to compute dvolrdt
    real(r8),parameter :: budget_tolerance = 1.0e-6   ! budget tolerance, m3/day
    logical  :: abort                       ! abort flag
    real(r8) :: sum1,sum2
    integer  :: yr, mon, day, ymd, tod      ! time information
    integer  :: nsub                        ! subcyling for cfl
    real(r8) :: delt                        ! delt associated with subcycling
    real(r8) :: delt_coupling               ! real value of coupling_period
    integer , save :: nsub_save             ! previous nsub
    logical , save :: first_call = .true.   ! first time flag (for backwards compatibility)
    character(len=256) :: filer             ! restart file name
    integer  :: cnt                         ! counter for gridcells
    integer  :: ier                         ! error code
    integer,parameter  :: dbug = 1          ! local debug flag
!scs
! parameters used in negative runoff partitioning algorithm
    real(r8) :: river_volume_minimum        ! gridcell area multiplied by average river_depth_minimum [m3]
    real(r8) :: qgwl_volume                 ! volume of runoff during time step [m3]
!scs
    character(len=*),parameter :: subname = '(Rtmrun) '
!-----------------------------------------------------------------------

    call t_startf('mosartr_tot')

    delt_coupling = coupling_period*1.0_r8
    budget_check = .true.   ! leave this on all the time for now
    if (first_call) then
       budget_accum = 0._r8
       budget_accum_cnt = 0
       if (masterproc) write(iulog,'(2a,g20.12)') trim(subname),' MOSART coupling period ',delt_coupling
    end if

    budget_terms = 0._r8

    flow = 0._r8
    eroup_lagi = 0._r8
    eroup_lagf = 0._r8
    erowm_regi = 0._r8
    erowm_regf = 0._r8
    eroutup_avg = 0._r8
    erlat_avg = 0._r8
    rtmCTL%runoff = 0._r8              ! coupler return mosart basin derived flow [m3/s]
    rtmCTL%direct = 0._r8              ! coupler return direct flow [m3/s]
    rtmCTL%flood = 0._r8               ! coupler return flood water sent back to clm [m3/s]
    rtmCTL%runofflnd = spval        ! runoff masked for land (m3 H2O/s) ( spval = 1.e36_r8 )
    rtmCTL%runoffocn = spval       ! runoff masked for ocn  (m3 H2O/s)
    rtmCTL%dvolrdt = 0._r8            ! RTM change in storage (mm/s)
    rtmCTL%dvolrdtlnd = spval      ! dvolrdt masked for land (mm/s)
    rtmCTL%dvolrdtocn = spval     ! dvolrdt masked for ocn  (mm/s)

    if (inundflag) then
       ! If inundation scheme is turned on :
       if ( Tctl%OPT_inund .eq. 1 ) then
          vol_chnl2fp = 0._r8
          vol_fp2chnl = 0._r8
          fa_fp_cplPeriod(:) = 0._r8
       end if

       totalVelo_down(:) = 0._r8
       totalVelo_up(:) = 0._r8
       chnlNum_down(:) = 0
       chnlNum_up(:) = 0
    end if

    if(heatflag) then
      rtmCTL%templand_Tqsur = spval
      rtmCTL%templand_Tqsub = spval
      rtmCTL%templand_Ttrib = spval
      rtmCTL%templand_Tchanr = spval
      do n = rtmCTL%begr,rtmCTL%endr
          if (rtmCTL%mask(n) .eq. 1 .or. rtmCTL%mask(n) .eq. 3) then
              rtmCTL%templand_Tqsur(n) = 0._r8
              rtmCTL%templand_Tqsub(n) = 0._r8
              rtmCTL%templand_Ttrib(n) = 0._r8
              rtmCTL%templand_Tchanr(n) = 0._r8
          end if
      end do
    end if
    
    if (budget_check) then
       call t_startf('mosartr_budget')
       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(bv_volt_i,nt) = budget_terms( bv_volt_i,nt) + rtmCTL%volr(nr,nt)
          budget_terms(bv_wt_i,nt) = budget_terms(bv_wt_i,nt) + TRunoff%wt(nr,nt)
          budget_terms(bv_wr_i,nt) = budget_terms(bv_wr_i,nt) + TRunoff%wr(nr,nt)
          budget_terms(bv_wh_i,nt) = budget_terms(bv_wh_i,nt) + TRunoff%wh(nr,nt)*rtmCTL%area(nr)
          budget_terms(br_qsur,nt) = budget_terms(br_qsur,nt) + rtmCTL%qsur(nr,nt)*delt_coupling     ! (rtmCTL%qsur 's unit is m^3/s. --Inund.)
          budget_terms(br_qsub,nt) = budget_terms(br_qsub,nt) + rtmCTL%qsub(nr,nt)*delt_coupling
          budget_terms(br_qgwl,nt) = budget_terms(br_qgwl,nt) + rtmCTL%qgwl(nr,nt)*delt_coupling
          budget_terms(br_qdto,nt) = budget_terms(br_qdto,nt) + rtmCTL%qdto(nr,nt)*delt_coupling
          budget_terms(br_qdem,nt) = budget_terms(br_qdem,nt) + rtmCTL%qdem(nr,nt)*delt_coupling
       enddo
       enddo

       ! If inundation scheme is turned on :
       if (inundflag .and. Tctl%OPT_inund .eq. 1 ) then
         do nr = rtmCTL%begr, rtmCTL%endr

           !if ( TUnit%mask( nr ) .gt. 0 ) then      ! 0--Ocean; 1--Land; 2--Basin outlet.
           if ( rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).
             budget_terms(bv_fp_i, 1) = budget_terms(bv_fp_i, 1) + TRunoff%wf_ini( nr )
             !budget_terms(bv_fp_i, 1) = budget_terms(bv_fp_i, 1) + rtmCTL%inundwf(nr)        ! 17-6-7
           endif

           ! land river two way coupling, update floodplain inundation volume with drainage from lnd
           if (use_lnd_rof_two_way) then
             TRunoff%wf_ini(nr) = TRunoff%wf_ini(nr) - rtmCTL%inundinf(nr) * coupling_period

             if ( TRunoff%wf_ini(nr) < 0._r8 ) then
               TRunoff%wr(nr, 1) = TRunoff%wr(nr, 1) + TRunoff%wf_ini(nr)
               TRunoff%wf_ini(nr) = 0._r8
               TRunoff%yr(nr, 1) = TRunoff%wr(nr, 1) / TUnit%rlen(nr) / TUnit%rwidth(nr)
             endif
           endif

         end do
       end if

       if (wrmflag) then
          StorWater%supply = 0._r8             !initial supply at the start
          nt = 1
          do nr = rtmCTL%begr,rtmCTL%endr
             budget_terms(bv_dsupp_i,nt) = budget_terms(bv_dsupp_i,nt) + StorWater%supply(nr)
          enddo
          do idam = 1,ctlSubwWRM%LocalNumDam
             budget_terms(bv_dstor_i,nt) = budget_terms(bv_dstor_i,nt) + StorWater%storage(idam)
          enddo
       endif
       call t_stopf('mosartr_budget')
    endif ! budget_check

    ! data for euler solver, in m3/s here
    do nr = rtmCTL%begr,rtmCTL%endr
    do nt = 1,nt_rtm
       TRunoff%qsur(nr,nt) = rtmCTL%qsur(nr,nt)
       TRunoff%qsub(nr,nt) = rtmCTL%qsub(nr,nt)
       TRunoff%qgwl(nr,nt) = rtmCTL%qgwl(nr,nt)
       TRunoff%qdem(nr,nt) = rtmCTL%qdem(nr,nt)
    enddo
    enddo
  
    !-----------------------------------
    ! Compute flood
    ! Remove water from mosart and send back to clm
    ! Just consider land points and only remove liquid water 
    ! rtmCTL%flood is m3/s here
    !-----------------------------------

    if (.not.inundflag) then
       call t_startf('mosartr_flood')
       nt = 1 
       rtmCTL%flood = 0._r8
       do nr = rtmCTL%begr,rtmCTL%endr
          ! initialize rtmCTL%flood to zero
          if (rtmCTL%mask(nr) == 1) then
             if (rtmCTL%volr(nr,nt) > rtmCTL%fthresh(nr)) then 
                ! determine flux that is sent back to the land
                ! this is in m3/s
                rtmCTL%flood(nr) = &
                     (rtmCTL%volr(nr,nt)-rtmCTL%fthresh(nr)) / (delt_coupling)

                ! rtmCTL%flood will be sent back to land - so must subtract this 
                ! from the input runoff from land
                ! tcraig, comment - this seems like an odd approach, you
                !   might create negative forcing.  why not take it out of
                !   the volr directly?  it's also odd to compute this
                !   at the initial time of the time loop.  why not do
                !   it at the end or even during the run loop as the
                !   new volume is computed.  fluxout depends on volr, so
                !   how this is implemented does impact the solution.
                TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) - rtmCTL%flood(nr)
             endif
          endif
       enddo
       call t_stopf('mosartr_flood')
    endif

    !-----------------------------------
    ! DIRECT sMAT transfer to outlet point using sMat
    ! Remember to subract water from TRunoff forcing
    !-----------------------------------

    if (barrier_timers) then
       call t_startf('mosartr_SMdirect_barrier')
       call mpi_barrier(mpicom_rof,ier)
       call t_stopf ('mosartr_SMdirect_barrier')
    endif

    call t_startf('mosartr_SMdirect')

    !--- copy direct transfer fields to AV
    !--- convert kg/m2s to m3/s
    call mct_avect_zero(avsrc_direct)
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       cnt = cnt + 1
       do nt = 1,nt_rtm

          !-------------------------------
          ! This water is routed directly to the outlet and passed out
          ! Turn on and off terms as desired via commenting them out
          !-------------------------------

          !---- all dto water, none was going to TRunoff ---
          !---- *** DO NOT TURN THIS ONE OFF, conservation will fail *** ---
          avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + rtmCTL%qdto(nr,nt)

          !---- negative gwl water less than channel volume, remove from TRunoff ---
          !---- scs
          qgwl_volume = TRunoff%qgwl(nr,nt) * delt_mosart
          river_volume_minimum = river_depth_minimum * rtmCTL%area(nr)
          !---- if qgwl is negative, and adding it to the main channel would bring 
          !---- main channel storage below a threshold, send qgwl directly to ocean
          if (((qgwl_volume + TRunoff%wr(nr,nt)) < river_volume_minimum) &
             .and. (TRunoff%qgwl(nr,nt) < 0._r8)) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
             TRunoff%qgwl(nr,nt) = 0._r8
          endif
          !---- scs

          !---- negative qgwl water, remove from TRunoff ---
          if (TRunoff%qgwl(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
             TRunoff%qgwl(nr,nt) = 0._r8
          endif

          !---- all gwl water, remove from TRunoff ---
          avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qgwl(nr,nt)
          TRunoff%qgwl(nr,nt) = 0._r8

          !---- negative qsub water, remove from TRunoff ---
          if (TRunoff%qsub(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qsub(nr,nt)
             TRunoff%qsub(nr,nt) = 0._r8
          endif

          !---- negative qsur water, remove from TRunoff ---
          if (TRunoff%qsur(nr,nt) < 0._r8) then
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + TRunoff%qsur(nr,nt)
             TRunoff%qsur(nr,nt) = 0._r8
          endif

          !---- water outside the basin ---
          !---- *** DO NOT TURN THIS ONE OFF, conservation will fail *** ---
          if (inundflag) then
             if ( rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

             else
                avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + &
                     TRunoff%qsub(nr,nt) + &
                     TRunoff%qsur(nr,nt) + &
                     TRunoff%qgwl(nr,nt)
                TRunoff%qsub(nr,nt) = 0._r8
                TRunoff%qsur(nr,nt) = 0._r8
                TRunoff%qgwl(nr,nt) = 0._r8
             end if

          else
             if (TUnit%mask(nr) > 0) then
                ! mosart euler
             else
                avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + &
                     TRunoff%qsub(nr,nt) + &
                     TRunoff%qsur(nr,nt) + &
                     TRunoff%qgwl(nr,nt)
                TRunoff%qsub(nr,nt) = 0._r8
                TRunoff%qsur(nr,nt) = 0._r8
                TRunoff%qgwl(nr,nt) = 0._r8
             end if
          end if

          !---- all nt=2 water ---
          !---- can turn off euler_calc for this tracer ----
          if (nt == 2) then
             TUnit%euler_calc(nt) = .false.
             avsrc_direct%rAttr(nt,cnt) = avsrc_direct%rAttr(nt,cnt) + &
                  TRunoff%qsub(nr,nt) + &
                  TRunoff%qsur(nr,nt) + &
                  TRunoff%qgwl(nr,nt)
             TRunoff%qsub(nr,nt) = 0._r8
             TRunoff%qsur(nr,nt) = 0._r8
             TRunoff%qgwl(nr,nt) = 0._r8
          endif

       enddo
    enddo
    call mct_avect_zero(avdst_direct)

    call mct_sMat_avMult(avsrc_direct, sMatP_direct, avdst_direct)

    !--- copy direct transfer water from AV to output field ---
    cnt = 0
    do nr = rtmCTL%begr,rtmCTL%endr
       cnt = cnt + 1
       do nt = 1,nt_rtm
          rtmCTL%direct(nr,nt) = avdst_direct%rAttr(nt,cnt)
       enddo
    enddo
    call t_stopf('mosartr_SMdirect')

    !-----------------------------------
    ! MOSART Subcycling
    !-----------------------------------

    call t_startf('mosartr_subcycling')

    if (first_call .and. masterproc) then
       do nt = 1,nt_rtm
          write(iulog,'(2a,i6,l4)') trim(subname),' euler_calc for nt = ',nt,TUnit%euler_calc(nt)
       enddo
    endif

    nsub = coupling_period/delt_mosart
    if (nsub*delt_mosart < coupling_period) then
       nsub = nsub + 1
    end if
    delt = delt_coupling/float(nsub)
  
    if (delt /= delt_save) then
       if (masterproc) then
          write(iulog,'(2a,2g20.12,2i12)') trim(subname),' MOSART delt update from/to',delt_save,delt,nsub_save,nsub
       end if
    endif
    nsub_save = nsub
    delt_save = delt
    Tctl%DeltaT = delt

    !-----------------------------------
    ! mosart euler solver
    ! --- convert TRunoff fields from m3/s to m/s before calling Euler
    !-----------------------------------

    do nt = 1,nt_rtm
    do nr = rtmCTL%begr,rtmCTL%endr
       TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) / rtmCTL%area(nr)
       TRunoff%qsub(nr,nt) = TRunoff%qsub(nr,nt) / rtmCTL%area(nr)
       TRunoff%qgwl(nr,nt) = TRunoff%qgwl(nr,nt) / rtmCTL%area(nr)
       TRunoff%qdem(nr,nt) = TRunoff%qdem(nr,nt) / rtmCTL%area(nr) !m3 to m
    enddo
    enddo

    do ns = 1,nsub

       ! this advances the model time one coupling period when ns=1
       ! we want the advance in the subcycling because of the new_month flag
       ! this should really be advancing delt per ns, not coupling_period per call, but
       ! need to modify clock initialization and restart for that to be supported
       ! and make sure a rational timestep doesn't lead to accumulating errors
       call advance_timestep(ns)

       call get_curr_date(yr, mon, day, tod)
       ymd = yr*10000 + mon*100 + day
       if (tod == 0 .and. masterproc) then
          write(iulog,*) ' '
          write(iulog,'(2a,i4,a,i10,i6)') trim(subname),' subcycling=',ns,': model date=',ymd,tod
       endif
     
       !if (inundflag .and. wrmflag .eq. 0) then !use Luo's scheme when inundation is on and WM is off (keep it for now - tz)
       !   call t_startf('mosartr_inund_sim')
       !   call MOSARTinund_simulate ( )
       !   call t_stopf('mosartr_inund_sim')
       !else ! other cases
       
       call t_startf('mosartr_euler')
#ifdef DEBUG
       write(iulog,*) 'clm-mosart subT: (call Euler) ns=', ns
#endif
       call Euler()
       call t_stopf('mosartr_euler')


! tcraig - NOT using this now, but leave it here in case it's useful in the future
!   for some runoff terms.
!       !-----------------------------------
!       ! downstream advection using sMat
!       !-----------------------------------
!
!       if (barrier_timers) then
!          call t_startf('mosartr_SMdnstrm_barrier')
!          call mpi_barrier(mpicom_rof,ier)
!          call t_stopf ('mosartr_SMdnstrm_barrier')
!       endif
!
!       call t_startf('mosartr_SMdnstrm')
!
!       !--- copy fluxout into avsrc_upstrm ---
!       cnt = 0
!       do n = rtmCTL%begr,rtmCTL%endr
!          cnt = cnt + 1
!          do nt = 1,nt_rtm
!             avsrc_upstrm%rAttr(nt,cnt) = fluxout(n,nt)
!          enddo
!       enddo
!       call mct_avect_zero(avdst_upstrm)
!
!       call mct_sMat_avMult(avsrc_upstrm, sMatP_upstrm, avdst_upstrm)
!
!       !--- add mapped fluxout to sfluxin ---
!       cnt = 0
!       sfluxin = 0._r8
!       do n = rtmCTL%begr,rtmCTL%endr
!          cnt = cnt + 1
!          do nt = 1,nt_rtm
!             sfluxin(n,nt) = sfluxin(n,nt) + avdst_upstrm%rAttr(nt,cnt)
!          enddo
!       enddo
!       call t_stopf('mosartr_SMdnstrm')

       !-----------------------------------
       ! accumulate local flow field
       !-----------------------------------

       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          flow(nr,nt) = flow(nr,nt) + TRunoff%flow(nr,nt)
          eroup_lagi(nr,nt) = eroup_lagi(nr,nt) + TRunoff%eroup_lagi(nr,nt)
          eroup_lagf(nr,nt) = eroup_lagf(nr,nt) + TRunoff%eroup_lagf(nr,nt)
          erowm_regi(nr,nt) = erowm_regi(nr,nt) + TRunoff%erowm_regi(nr,nt)
          erowm_regf(nr,nt) = erowm_regf(nr,nt) + TRunoff%erowm_regf(nr,nt)
          eroutup_avg(nr,nt) = eroutup_avg(nr,nt) + TRunoff%eroutup_avg(nr,nt)
          erlat_avg(nr,nt) = erlat_avg(nr,nt) + TRunoff%erlat_avg(nr,nt)
       enddo
       enddo


       if (inundflag) then
          ! If 'budget_check' is true & inundation scheme is turned on :
          if ( inundflag .and. budget_check .and. Tctl%OPT_inund .eq. 1 ) then

             do nt = 1, 1
                do nr = rtmCTL%begr, rtmCTL%endr

                   !if (Tunit%mask(nr) .eq. 1 .or. Tunit%mask(nr) .eq. 2) then     ! 1 -- Land; 2 -- Basin outlet.
                   if (rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3) then    ! 1 -- Land; 3 -- Basin outlet.

                      ! Accumulate flow from main channel to floodplain (for all local grid cells and all sub-steps of coupling period):
                      if ( TRunoff%se_rf(nr) .gt. 0._r8 ) then
                         vol_chnl2fp = vol_chnl2fp + TRunoff%se_rf(nr)

                         ! Accumulate flow from floodplain to main channel (for all local grid cells and all sub-steps of coupling period):
                      elseif ( TRunoff%se_rf(nr) .lt. 0._r8 ) then
                         vol_fp2chnl = vol_fp2chnl - TRunoff%se_rf(nr)
                      end if

                      ! Accumulate inundated floodplain area for all sub-steps of coupling period (for each land grid cell):
                      fa_fp_cplPeriod(nr) = fa_fp_cplPeriod(nr) + TRunoff%fa_fp(nr)
                   end if

                end do
             end do
          end if

          if ( budget_check ) then

             do nt = 1, nt_rtm
                do nr = rtmCTL%begr, rtmCTL%endr

                   if (rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3) then    ! 1 -- Land; 3 -- Basin outlet.
                      !if (Tunit%mask(nr) .eq. 1 .or. Tunit%mask(nr) .eq. 2) then     ! 1 -- Land; 2 -- Basin outlet.

                      ! Flow velocity is downward or zero:
                      if ( TRunoff%vr( nr, nt ) .ge. 0._r8 ) then
                         totalVelo_down(nt) = totalVelo_down(nt) + TRunoff%vr( nr, nt )
                         chnlNum_down(nt) = chnlNum_down(nt) + 1

                         ! Flow velocity is upward (is negative):
                      else
                         totalVelo_up(nt) = totalVelo_up(nt) + TRunoff%vr( nr, nt )
                         chnlNum_up(nt) = chnlNum_up(nt) + 1
                      end if

                   end if
                end do
             end do
          end if
       end if

       if(heatflag) then
         do n = rtmCTL%begr,rtmCTL%endr
             if (rtmCTL%mask(n) .eq. 1 .or. rtmCTL%mask(n) .eq. 3) then
                 rtmCTL%templand_Tqsur(n) = rtmCTL%templand_Tqsur(n) + THeat%Tqsur(n)
                 rtmCTL%templand_Tqsub(n) = rtmCTL%templand_Tqsub(n) + THeat%Tqsub(n)
                 rtmCTL%templand_Ttrib(n) = rtmCTL%templand_Ttrib(n) + THeat%Tt_avg(n)
                 rtmCTL%templand_Tchanr(n) = rtmCTL%templand_Tchanr(n) + THeat%Tr_avg(n)
             end if
         enddo
       end if
       
    enddo ! nsub

    !-----------------------------------
    ! average flow over subcycling
    !-----------------------------------

    flow        = flow        / float(nsub)
    eroup_lagi  = eroup_lagi  / float(nsub)
    eroup_lagf  = eroup_lagf  / float(nsub)
    erowm_regi  = erowm_regi  / float(nsub)
    erowm_regf  = erowm_regf  / float(nsub)
    eroutup_avg = eroutup_avg / float(nsub)
    erlat_avg   = erlat_avg   / float(nsub)

    if (inundflag) then
       ! Mean inundated floodplain area for all sub-steps of coupling period (for each land grid cell):
       fa_fp_cplPeriod = fa_fp_cplPeriod / float(nsub)
    endif

    !-----------------------------------
    ! update states when subsycling completed
    !-----------------------------------

    rtmCTL%wh      = TRunoff%wh
    rtmCTL%wt      = TRunoff%wt
    rtmCTL%wr      = TRunoff%wr
    rtmCTL%mr      = TRunoff%mr
    rtmCTL%pr      = TRunoff%pr
    rtmCTL%yr      = TRunoff%yr
    rtmCTL%rr      = TRunoff%rr
    rtmCTL%erout   = TRunoff%erout

    ! If inundation scheme is turned on :
    if (inundflag .and. Tctl%OPT_inund .eq. 1 ) then
      rtmCTL%inundwf(:) = TRunoff%wf_ini(:)
      rtmCTL%inundhf(:) = TRunoff%hf_ini(:)
      rtmCTL%inundff(:) = TRunoff%ff_ini(:)
      rtmCTL%inundffunit(:) = TRunoff%ffunit_ini(:)
    end if

    if (heatflag) then
      rtmCTL%Tqsur   = THeat%Tqsur
      rtmCTL%Tqsub   = THeat%Tqsub
      rtmCTL%Tt      = THeat%Tt
      rtmCTL%Tr      = THeat%Tr
      rtmCTL%Ha_rout   = THeat%Ha_rout
    
      do n = rtmCTL%begr,rtmCTL%endr
         if(rtmCTL%mask(n) .eq. 1 .or. rtmCTL%mask(n) .eq. 3) then
            rtmCTL%templand_Tqsur(n) = rtmCTL%templand_Tqsur(n) / float(nsub)
            rtmCTL%templand_Tqsub(n) = rtmCTL%templand_Tqsub(n) / float(nsub)
            rtmCTL%templand_Ttrib(n) = rtmCTL%templand_Ttrib(n) / float(nsub)
            rtmCTL%templand_Tchanr(n) = rtmCTL%templand_Tchanr(n) / float(nsub)
         else
            rtmCTL%templand_Tqsur(n) = spval
            rtmCTL%templand_Tqsub(n) = spval
            rtmCTL%templand_Ttrib(n) = spval
            rtmCTL%templand_Tchanr(n) = spval
         end if
      end do
    end if

    do nt = 1,nt_rtm
    do nr = rtmCTL%begr,rtmCTL%endr
       volr_init = rtmCTL%volr(nr,nt)
       rtmCTL%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + &
                             TRunoff%wh(nr,nt)*rtmCTL%area(nr)) * TUnit%frac(nr)
       if (inundflag .and. Tctl%OPT_inund == 1 .and. nt == 1) then
          rtmCTL%volr(nr,nt) = rtmCTL%volr(nr,nt) + TRunoff%wf_ini(nr)
       endif
       rtmCTL%dvolrdt(nr,nt) = (rtmCTL%volr(nr,nt) - volr_init) / delt_coupling
       rtmCTL%runoff(nr,nt) = flow(nr,nt)

       rtmCTL%runofftot(nr,nt) = rtmCTL%direct(nr,nt)      ! coupler return direct flow [m3/s]
       if (rtmCTL%mask(nr) == 1) then                      ! 1=Land (rtmCTL%mask are consistent with Tunit%mask ? --Inund.)

          rtmCTL%runofflnd(nr,nt) = rtmCTL%runoff(nr,nt)   ! (rtmCTL%runoff is main channel outflow (m^3/s). --Inund.)
          rtmCTL%dvolrdtlnd(nr,nt)= rtmCTL%dvolrdt(nr,nt)  ! RTM change in storage (m^3/s, --Inund.)

       elseif (rtmCTL%mask(nr) >= 2) then                  ! 2=ocean, 3=outlet

          rtmCTL%runoffocn(nr,nt) = rtmCTL%runoff(nr,nt)   ! (For ocean, rtmCTL%runoff should be zero ? --Inund.)
          rtmCTL%runofftot(nr,nt) = rtmCTL%runofftot(nr,nt) + rtmCTL%runoff(nr,nt)   ! total runoff masked for ocn (direct flow + channel outflow, --Inund.)
          rtmCTL%dvolrdtocn(nr,nt)= rtmCTL%dvolrdt(nr,nt)  ! RTM change in storage (For ocean, storage [rtmCTL%volr] and storage change [rtmCTL%dvolrdt] should be zeros. --Inund.)
       endif
    enddo
    enddo

    call t_stopf('mosartr_subcycling')

    !-----------------------------------
    ! BUDGET
    !-----------------------------------

    budget_write = .false.
    if (day == 1 .and. mon == 1) budget_write = .true.
    if (tod == 0) budget_write = .true.

    if (budget_check) then
       call t_startf('mosartr_budget')
       do nt = 1,nt_rtm
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(bv_volt_f,nt) = budget_terms(bv_volt_f,nt) + rtmCTL%volr(nr,nt)                     ! (Final total volume at end of coupling period. --Inund.)
          budget_terms(bv_wt_f,nt) = budget_terms(bv_wt_f,nt) + TRunoff%wt(nr,nt)
          budget_terms(bv_wr_f,nt) = budget_terms(bv_wr_f,nt) + TRunoff%wr(nr,nt)
          budget_terms(bv_wh_f,nt) = budget_terms(bv_wh_f,nt) + TRunoff%wh(nr,nt)*rtmCTL%area(nr)
          budget_terms(br_direct,nt) = budget_terms(br_direct,nt) + rtmCTL%direct(nr,nt)*delt_coupling     ! (Volume of direct 'flow' to ocean. --Inund.)

          if (rtmCTL%mask(nr) >= 2) then    ! (2 -- Ocean; 3 -- Outlet. --Inund.)

             budget_terms(br_ocnout,nt) = budget_terms(br_ocnout,nt) + rtmCTL%runoff(nr,nt)*delt_coupling  ! (Volume of outflows to oceans. Note that rtmCTL%runoff is averge value of sub-steps used by MOSART. --Inund.)
             budget_terms(br_erolpo,nt) = budget_terms(br_erolpo,nt) + eroup_lagi(nr,nt)*delt_coupling     ! (Volume of outflows to oceans in the previous MOSART sub-step. 
                                                                                                           ! Also eroup_lagi is averge value of several previous MOSART sub-steps. --Inund.)
             budget_terms(br_erolco,nt) = budget_terms(br_erolco,nt) + eroup_lagf(nr,nt)*delt_coupling     ! (Volume of outflows to oceans. Actually eroup_lagf is same as the above rtmCTL%runoff. --Inund.)
             budget_terms(br_erorpo,nt) = budget_terms(br_erorpo,nt) + erowm_regi(nr,nt)*delt_coupling
             budget_terms(br_erorco,nt) = budget_terms(br_erorco,nt) + erowm_regf(nr,nt)*delt_coupling
          else       ! (1--Land. --Inund.)
             budget_terms(br_lndout,nt) = budget_terms(br_lndout,nt) + rtmCTL%runoff(nr,nt)*delt_coupling  ! (Sum up channel outflow volumes of all land cells. --Inund.)
             budget_terms(br_erolpn,nt) = budget_terms(br_erolpn,nt) + eroup_lagi(nr,nt)*delt_coupling     ! (Sum up channel outflow volumes of all land cells in the previous MOSART sub-step. --Inund.)
             budget_terms(br_erolcn,nt) = budget_terms(br_erolcn,nt) + eroup_lagf(nr,nt)*delt_coupling     ! (Sum up channel outflow volumes of all land cells. Actually eroup_lagf is same as the above rtmCTL%runoff. --Inund.)
             budget_terms(br_erorpn,nt) = budget_terms(br_erorpn,nt) + erowm_regi(nr,nt)*delt_coupling
             budget_terms(br_erorcn,nt) = budget_terms(br_erorcn,nt) + erowm_regf(nr,nt)*delt_coupling
          endif
          budget_terms(br_eroutup,nt) = budget_terms(br_eroutup,nt) - eroutup_avg(nr,nt)*delt_coupling     ! (Sum up upstream inflow volumes of all land cells. --Inund.)
          budget_terms(br_erlat,nt) = budget_terms(br_erlat,nt) - erlat_avg(nr,nt)*delt_coupling           ! (Sum up lateral inflow volumes of all land cells. --Inund.)
       enddo
       enddo
       nt = 1
       do nr = rtmCTL%begr,rtmCTL%endr
          budget_terms(br_flood,nt) = budget_terms(br_flood,nt) + rtmCTL%flood(nr)*delt_coupling           ! (This flood volume is from former method. It equals zero when the inundation scheme is turned on. --Inund.)
       enddo
       if (wrmflag) then
          nt = 1
          do nr = rtmCTL%begr,rtmCTL%endr
             budget_terms(bv_dsupp_f,nt) = budget_terms(bv_dsupp_f,nt) + StorWater%supply(nr)
             ! convert supply from m3 per coupling delta (3hrs)  to mm/s (N. Sun)
             if (StorWater%supply(nr) > 0) then            
               StorWater%supply(nr) = StorWater%supply(nr)/delt_coupling               
             endif
          end do
          do idam = 1,ctlSubwWRM%LocalNumDam
             budget_terms(bv_dstor_f,nt) = budget_terms(bv_dstor_f,nt) + StorWater%storage(idam)
          enddo
       endif

       if (inundflag) then
          do nt = 1, nt_rtm
             do nr = rtmCTL%begr, rtmCTL%endr
                if (rtmCTL%mask(nr) .eq. 3) then        ! 3 -- Basin outlet (downstream is ocean).
                   !if ( Tunit%mask(nr) .eq. 2 ) then      ! 2 -- Basin outlet (downstream is ocean).
                   ! Total streamflow (flow rate) from land to oceans (m^3/s):
                   budget_terms(br_landOutflow, nt) = budget_terms(br_landOutflow, nt) + rtmCTL%runoff(nr, nt)
                end if
             end do
          end do

          ! If inundation scheme is turned on :
          if (inundflag .and. Tctl%OPT_inund .eq. 1 ) then

             ! Total volume of flows from main channels to floodplains (for all local grid cells and all sub-steps of coupling period) (m^3):
             budget_terms(bv_chnl2fp, 1) = vol_chnl2fp

             ! Total volume of flows from floodplains to main channels (for all local grid cells and all sub-steps of coupling period) (m^3):
             budget_terms(bv_fp2chnl, 1) = vol_fp2chnl

             do nr = rtmCTL%begr, rtmCTL%endr
                if (rtmCTL%mask(nr) .eq. 1 .or. rtmCTL%mask(nr) .eq. 3) then      ! 1 -- Land; 3 -- Basin outlet.
                   !if ( Tunit%mask(nr) .eq. 1 .or. Tunit%mask(nr) .eq. 2 ) then     ! 1--Land; 2--Basin outlet (downstream is ocean).

                   ! Water volume over floodplains at coupling-period end:
                   if ( Tctl%OPT_inund .eq. 1 ) then
                      budget_terms(bv_fp_f, 1) = budget_terms(bv_fp_f, 1) + TRunoff%wf_ini(nr)
                   endif

                   ! Sum up channel surface area:
                   budget_terms(bi_mainChnlArea, 1) = budget_terms(bi_mainChnlArea, 1) + TUnit%area(nr) * TUnit%frac(nr) * TUnit%a_chnl(nr)

                   ! Sum up flooded area (including channel area):
                   budget_terms(bi_floodedArea, 1) = budget_terms(bi_floodedArea, 1) + fa_fp_cplPeriod(nr) + TUnit%area(nr) * TUnit%frac(nr) * TUnit%a_chnl(nr)

                   ! Sum up land area:
                   budget_terms(bi_landArea, 1) = budget_terms(bi_landArea, 1) + TUnit%area(nr) * TUnit%frac(nr)
                end if
             end do
          end if

          do nt = 1, nt_rtm
             ! Total downward flow velocity (for all local grid cells and all sub-steps of coupling period) (m/s):
             budget_terms(bVelo_downward, nt) = totalVelo_down(nt)

             ! Total upward flow velocity (for all local grid cells and all sub-steps of coupling period) (m/s):
             budget_terms(bVelo_upward, nt) = totalVelo_up(nt)

             ! Total number of channels with downward flow velocities (for all local grid cells and all sub-steps of coupling period):
             budget_terms(bVelo_downChnlNo, nt) = chnlNum_down(nt)

             ! Total number of channels with upward flow velocities (for all local grid cells and all sub-steps of coupling period):
             budget_terms(bVelo_upChnlNo, nt) = chnlNum_up(nt)
          end do

          !--- Debug: compare rtmCTL%mask() and Tunit%mask() ---
          do nr = rtmCTL%begr, rtmCTL%endr
             ! Tunit%mask() is land, but rtmCTL%mask() is not land :
             if ( Tunit%mask(nr) .eq. 1 .and. rtmCTL%mask(nr) .ne. 1 ) then
                budget_terms(bmask_lndErr, 1) = budget_terms(bmask_lndErr, 1) + 1
             end if

             ! Tunit%mask() is basin outlet, but rtmCTL%mask() is not basin outlet :
             if ( Tunit%mask(nr) .eq. 2 .and. rtmCTL%mask(nr) .ne. 3 ) then
                budget_terms(bmask_outErr, 1) = budget_terms(bmask_outErr, 1) + 1
             end if

             ! Tunit%mask() is ocean, but rtmCTL%mask() is not ocean :
             if ( Tunit%mask(nr) .eq. 0 .and. rtmCTL%mask(nr) .ne. 2 ) then
                budget_terms(bmask_ocnErr, 1) = budget_terms(bmask_ocnErr, 1) + 1
             end if

             ! Tunit%mask() is land, rtmCTL%mask() is land :
             if ( Tunit%mask(nr) .eq. 1 .and. rtmCTL%mask(nr) .eq. 1 ) then
                budget_terms(bmask_uLndrLnd, 1) = budget_terms(bmask_uLndrLnd, 1) + 1
             end if

             ! Tunit%mask() is land, rtmCTL%mask() is outlet :
             if ( Tunit%mask(nr) .eq. 1 .and. rtmCTL%mask(nr) .eq. 3 ) then
                budget_terms(bmask_uLndrOut, 1) = budget_terms(bmask_uLndrOut, 1) + 1
             end if

             ! Tunit%mask() is land, rtmCTL%mask() is ocean :
             if ( Tunit%mask(nr) .eq. 1 .and. rtmCTL%mask(nr) .eq. 2 ) then
                budget_terms(bmask_uLndrOcn, 1) = budget_terms(bmask_uLndrOcn, 1) + 1
             end if

             ! Tunit%mask() is outlet, rtmCTL%mask() is land :
             if ( Tunit%mask(nr) .eq. 2 .and. rtmCTL%mask(nr) .eq. 1 ) then
                budget_terms(bmask_uOutrLnd, 1) = budget_terms(bmask_uOutrLnd, 1) + 1
             end if

             ! Tunit%mask() is outlet, rtmCTL%mask() is outlet :
             if ( Tunit%mask(nr) .eq. 2 .and. rtmCTL%mask(nr) .eq. 3 ) then
                budget_terms(bmask_uOutrOut, 1) = budget_terms(bmask_uOutrOut, 1) + 1
             end if

             ! Tunit%mask() is outlet, rtmCTL%mask() is ocean :
             if ( Tunit%mask(nr) .eq. 2 .and. rtmCTL%mask(nr) .eq. 2 ) then
                budget_terms(bmask_uOutrOcn, 1) = budget_terms(bmask_uOutrOcn, 1) + 1
             end if

             ! Tunit%mask() is ocean, rtmCTL%mask() is land :
             if ( Tunit%mask(nr) .eq. 0 .and. rtmCTL%mask(nr) .eq. 1 ) then
                budget_terms(bmask_uOcnrLnd, 1) = budget_terms(bmask_uOcnrLnd, 1) + 1
             end if

             ! Tunit%mask() is ocean, rtmCTL%mask() is outlet :
             if ( Tunit%mask(nr) .eq. 0 .and. rtmCTL%mask(nr) .eq. 3 ) then
                budget_terms(bmask_uOcnrOut, 1) = budget_terms(bmask_uOcnrOut, 1) + 1
             end if

             ! Tunit%mask() is ocean, rtmCTL%mask() is ocean :
             if ( Tunit%mask(nr) .eq. 0 .and. rtmCTL%mask(nr) .eq. 2 ) then
                budget_terms(bmask_uOcnrOcn, 1) = budget_terms(bmask_uOcnrOcn, 1) + 1
             end if

          end do
          !--- Debug ---
       endif

       ! accumulate the budget total over the run to make sure it's decreasing on avg
       budget_accum_cnt = budget_accum_cnt + 1
       do nt = 1,nt_rtm
          budget_volume =  budget_terms(bv_volt_f,nt) - budget_terms(bv_volt_i,nt) + &
                           budget_terms(bv_dstor_f,nt) - budget_terms(bv_dstor_i,nt)             ! (Volume change during a coupling period. --Inund.)
          budget_input  =  budget_terms(br_qsur,nt) + budget_terms(br_qsub,nt) + &
                           budget_terms(br_qgwl,nt) + budget_terms(br_qdto,nt)
          budget_output =  budget_terms(br_ocnout,nt) + budget_terms(br_flood,nt) + &
                           budget_terms(br_direct,nt) + &
                           budget_terms(bv_dsupp_f,nt) - budget_terms(bv_dsupp_i,nt)
          budget_accum(nt) = budget_accum(nt) + budget_volume - budget_input + budget_output     ! (Accumulate the water balance errors. 
                                                                                                 ! 'budget_volume - budget_input + budget_output' should not be zero. 
                                                                                                 ! Because water is not balanced for grid cells of one computation node. --Inund.)
          budget_terms(bv_naccum,nt) = budget_accum(nt)/budget_accum_cnt                         ! (Average water balance error of all coupling periods up to now. --Inund.)
       enddo
       call t_stopf('mosartr_budget')
    endif  ! budget_check

    if (budget_check .and. budget_write) then
       call t_startf('mosartr_budget')
       !--- check budget

       ! convert terms from m3 to million m3
       budget_terms(:,:) = budget_terms(:,:) * 1.0e-6_r8

       ! global sum
       call shr_mpi_sum(budget_terms,budget_global,mpicom_rof,'mosart global budget',all=.false.)

       ! write budget
       if (masterproc) then
          write(iulog,'(2a,i10,i6)') trim(subname),' MOSART BUDGET diagnostics (million m3) for ',ymd,tod
          do nt = 1,nt_rtm
            budget_volume = (budget_global(bv_volt_f,nt) - budget_global(bv_volt_i,nt) + &
                             budget_global(bv_dstor_f,nt) - budget_global(bv_dstor_i,nt))   !(Global volume change during a coupling period. --Inund.)
            budget_input  = (budget_global(br_qsur,nt) + budget_global(br_qsub,nt) + &
                             budget_global(br_qgwl,nt) + budget_global(br_qdto,nt))
            budget_output = (budget_global(br_ocnout,nt) + budget_global(br_flood,nt) + &
                             budget_global(br_direct,nt) + &
                             budget_global(bv_dsupp_f,nt) - budget_global(bv_dsupp_i,nt))
            ! erout lag, need to remove current term and add in previous term, current term used in next timestep
            budget_other  = budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt) + &   !('previous MOSART sub-step channel outflow volume'-'current MOSART sub-step channel outflow volume'. --Inund.)
                            budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)       !(When WRM module is on: 'previous MOSART sub-step channel outflow volume'-'current MOSART sub-step channel outflow volume'. --Inund.)
            budget_total  = budget_volume - budget_input + budget_output - budget_other     !('budget_total' is supposed to be zero if water balance is perfect. --Inund.)

            write(iulog,'(2a)') trim(subname),'-----------------------------------------------------------------'
            write(iulog,'(2a,i4)')        trim(subname),'  tracer = ',nt
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wh    = ',nt,budget_global(bv_wh_f,nt)-budget_global(bv_wh_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wt    = ',nt,budget_global(bv_wt_f,nt)-budget_global(bv_wt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wr    = ',nt,budget_global(bv_wr_f,nt)-budget_global(bv_wr_i,nt)

            ! If inundation scheme is turned on :
            if (inundflag .and. Tctl%OPT_inund .eq. 1 ) then
              if ( nt .eq. 1 ) then
                write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume wf    = ',nt,budget_global(bv_fp_f,nt)-budget_global(bv_fp_i,nt)
              end if
            end if

            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   dvolume dstor = ',nt,budget_global(bv_dstor_f,nt)-budget_global(bv_dstor_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * dvolume total = ',nt,budget_volume   !(Global volume change during a coupling period. --Inund.)
          if (output_all_budget_terms) then           
            if (inundflag .and. Tctl%OPT_inund .eq. 1 .and. nt .eq. 1) then                                                                 
                write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x dvolume check = ',nt,budget_volume - &
                                                                             (budget_global(bv_wh_f,nt)-budget_global(bv_wh_i,nt) + &
                                                                              budget_global(bv_wt_f,nt)-budget_global(bv_wt_i,nt) + &
                                                                              budget_global(bv_wr_f,nt)-budget_global(bv_wr_i,nt) + &
                                                                              budget_global(bv_fp_f,nt)-budget_global(bv_fp_i,nt) + &
                                                                              budget_global(bv_dstor_f,nt)-budget_global(bv_dstor_i,nt)),&
                                                                              ' (should be zero)'
            else
                write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x dvolume check = ',nt,budget_volume - &
                                                                             (budget_global(bv_wh_f,nt)-budget_global(bv_wh_i,nt) + &
                                                                              budget_global(bv_wt_f,nt)-budget_global(bv_wt_i,nt) + &
                                                                              budget_global(bv_wr_f,nt)-budget_global(bv_wr_i,nt) + &
                                                                              budget_global(bv_dstor_f,nt)-budget_global(bv_dstor_i,nt)),&
                                                                              ' (should be zero)'
            endif
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volume   init = ',nt,budget_global(bv_volt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volume  final = ',nt,budget_global(bv_volt_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumeh  init = ',nt,budget_global(bv_wh_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumeh final = ',nt,budget_global(bv_wh_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumet  init = ',nt,budget_global(bv_wt_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumet final = ',nt,budget_global(bv_wt_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumer  init = ',nt,budget_global(bv_wr_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x volumer final = ',nt,budget_global(bv_wr_f,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x storage  init = ',nt,budget_global(bv_dstor_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x storage final = ',nt,budget_global(bv_dstor_f,nt)
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input surface = ',nt,budget_global(br_qsur,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input subsurf = ',nt,budget_global(br_qsub,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input gwl     = ',nt,budget_global(br_qgwl,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   input dto     = ',nt,budget_global(br_qdto,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * input total   = ',nt,budget_input
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x input check   = ',nt,budget_input - &
                                                                             (budget_global(br_qsur,nt)+budget_global(br_qsub,nt)+ &
                                                                              budget_global(br_qgwl,nt)+budget_global(br_qdto,nt)), &
                     ' (should be zero)'
                                                                             
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output runoff = ',nt,budget_global(br_ocnout,nt)   !(Outflows to oceans. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output direct = ',nt,budget_global(br_direct,nt)   !(Direct flows to oceans. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output flood  = ',nt,budget_global(br_flood,nt)    !(Former flood to land. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   output supply = ',nt,budget_global(bv_dsupp_f,nt)-budget_global(bv_dsupp_i,nt)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * output total  = ',nt,budget_output
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x output check  = ',nt,budget_output - &
                                                                             (budget_global(br_ocnout,nt) + budget_global(br_direct,nt) + &
                                                                              budget_global(br_flood,nt) + &
                                                                              budget_global(bv_dsupp_f,nt)-budget_global(bv_dsupp_i,nt)), &
                                                                             ' (should be zero)'
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
            !(--Inund. 'previous MOSART sub-step channel outflow volume'-'current MOSART sub-step channel outflow volume': )
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   other dwn lag = ',nt,budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt)

            !(--Inund. When WRM module is on: 'previous MOSART sub-step channel outflow volume'-'current MOSART sub-step channel outflow volume': )
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   other reg lnd = ',nt,budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)

            write(iulog,'(2a,i4,f22.6  )') trim(subname),' * other total   = ',nt,budget_other
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x other check   = ',nt,budget_other - &
                                                                            (budget_global(br_erolpn,nt) - budget_global(br_erolcn,nt) + &
                                                                             budget_global(br_erorpn,nt) - budget_global(br_erorcn,nt)), &
                                                                            ' (should be zero)'
          endif
            write(iulog,'(2a)') trim(subname),'----------------'
          if (output_all_budget_terms) then
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum dvolume   = ',nt,budget_volume     !(Global volume change during a coupling period. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum input     = ',nt,budget_input
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum output    = ',nt,budget_output
            write(iulog,'(2a,i4,f22.6  )') trim(subname),'   sum other     = ',nt,budget_other      !(Channel outflow volume difference between previous and current MOSART sub-step. --Inund.) 
          endif
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' * sum budget ** = ',nt,budget_total,' (should be zero, dv-in+out-oth)'
          if (output_all_budget_terms) then
            ! accum budget is just dv-i+o and should show that over time, the other terms go to zero (lag yes, reg land no)
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6,a)') trim(subname),' x accum budget  = ',nt,budget_global(bv_naccum,nt),' (should tend to zero over run, dv-in+out)'          ! (Average water balance error of all coupling periods up to now. --Inund.)
            write(iulog,'(2a)') trim(subname),'----------------'
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x runoff     ocn= ',nt,budget_global(br_ocnout,nt)      !(Outflows to oceans. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagi ocn= ',nt,budget_global(br_erolpo,nt)      !(Outflows to oceans in the previous MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagf ocn= ',nt,budget_global(br_erolco,nt)      !(Outflows to oceans in the current MOSART sub-step. 
                                                                                                                   !Actually 'eroup_lagf' is same as the above 'runoff'. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regi ocn= ',nt,budget_global(br_erorpo,nt)      !(When WRM module is on: Outflows to oceans in the previous MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regf ocn= ',nt,budget_global(br_erorco,nt)      !(When WRM module is on: Outflows to oceans in the current MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x other reg ocn = ',nt,budget_global(br_erorpo,nt) - budget_global(br_erorco,nt) !(Difference between the above two amounts. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x runoff     lnd= ',nt,budget_global(br_lndout,nt)      !(Total volume of channel outflows of all land cells. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagi lnd= ',nt,budget_global(br_erolpn,nt)      !(Total volume of all channel outflows in the previous MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroup_lagf lnd= ',nt,budget_global(br_erolcn,nt)      !(Total volume of all channel outflows in the current MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regi lnd= ',nt,budget_global(br_erorpn,nt)      !(When WRM module is on: Total volume of all channel outflows in the previous MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erowm_regf lnd= ',nt,budget_global(br_erorcn,nt)      !(When WRM module is on: Total volume of all channel outflows in the current MOSART sub-step. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x eroutup_avg   = ',nt,budget_global(br_eroutup,nt)     !(Total volume of upstream inflow amounts of all land cells. --Inund.)
            write(iulog,'(2a,i4,f22.6  )') trim(subname),' x erlateral     = ',nt,budget_global(br_erlat,nt)       !(Total volume of lateral inflow amounts of all main channels. --Inund.)
            write(iulog,'(2a)') trim(subname),'----------------'
          endif

            if ((budget_total) > 1.0e-6) then
               write(iulog,'(2a,i4)') trim(subname),' ***** BUDGET WARNING error gt 1. m3 for nt = ',nt
            endif
          enddo   ! (do nt = 1,nt_rtm   --Inund.)
          write(iulog,'(a)') '----------------------------------- '

          if (inundflag) then
             write(iulog,'(a)') ' '
             write(iulog,'(a)') '=================================================== '
             write(iulog,'(a)') ' '

             write(iulog,'(a)') 'MOSART-Inundation simulation information (global total or average values):'

             ! Unit conversion (from million m^3 to km^3 (cubic kilometer) ):
             !budget_glb_km3 = budget_global / 1000._r8

             ! Unit conversion:

             ! From million m^3 to km^3 :
             budget_glb_inund(bv_volt_i, 1) = budget_global(bv_volt_i, 1) / 1000._r8   ! Initial total surface-water volume.
             budget_glb_inund(bv_volt_f, 1) = budget_global(bv_volt_f, 1) / 1000._r8   ! Final total surface-water volume.
             budget_glb_inund(bv_wt_i, 1) = budget_global(bv_wt_i, 1) / 1000._r8       ! Initial water volume in subnetworks.
             budget_glb_inund(bv_wt_f, 1) = budget_global(bv_wt_f, 1) / 1000._r8       ! Final water volume in subnetworks.
             budget_glb_inund(bv_wr_i, 1) = budget_global(bv_wr_i, 1) / 1000._r8       ! Initial water volume in main channels.
             budget_glb_inund(bv_wr_f, 1) = budget_global(bv_wr_f, 1) / 1000._r8       ! Final water volume in main channels.
             budget_glb_inund(bv_wh_i, 1) = budget_global(bv_wh_i, 1) / 1000._r8       ! Inital water volume over hillslopes.
             budget_glb_inund(bv_wh_f, 1) = budget_global(bv_wh_f, 1) / 1000._r8       ! Final water volume over hillslopes.
             budget_glb_inund(br_qsur, 1) = budget_global(br_qsur, 1) / 1000._r8       ! Input surface runoff.
             budget_glb_inund(br_qsub, 1) = budget_global(br_qsub, 1) / 1000._r8       ! Input sub-surface runoff.
             budget_glb_inund(br_qgwl, 1) = budget_global(br_qgwl, 1) / 1000._r8       ! Input from glacier, wetland or lake.
             budget_glb_inund(br_qdto, 1) = budget_global(br_qdto, 1) / 1000._r8       ! Input direct-to-ocean runoff.
             budget_glb_inund(br_ocnout, 1) = budget_global(br_ocnout, 1) / 1000._r8   ! Output flows to oceans.
             budget_glb_inund(br_direct, 1) = budget_global(br_direct, 1) / 1000._r8   ! Output direct to oceans.
             budget_glb_inund(br_erolpn, 1) = budget_global(br_erolpn, 1) / 1000._r8   ! Output from previous sub-step
             budget_glb_inund(br_erolcn, 1) = budget_global(br_erolcn, 1) / 1000._r8   ! Output from current sub-step

             if (inundflag .and. Tctl%OPT_inund .eq. 1 ) then
                budget_glb_inund(bv_fp_i, 1) = budget_global(bv_fp_i, 1) / 1000._r8         ! From million m^3 to km^3 (initial water volume over floodplains).
                budget_glb_inund(bv_fp_f, 1) = budget_global(bv_fp_f, 1) / 1000._r8         ! From million m^3 to km^3 (final water volume over floodplains).
                budget_glb_inund(bv_chnl2fp, 1) = budget_global(bv_chnl2fp, 1) / 1000._r8   ! From million m^3 to km^3 (volume of flows from main channels to floodplains).
                budget_glb_inund(bv_fp2chnl, 1) = budget_global(bv_fp2chnl, 1) / 1000._r8   ! From million m^3 to km^3 (volume of flows from floodplains to main channels).
                budget_glb_inund(bi_landArea, 1) = budget_global(bi_landArea, 1) * 1.0e6_r8 / 1.0e9_r8         ! First recovered to m^2, then converted to thousand km^2 (total land area).
                budget_glb_inund(bi_floodedArea, 1) = budget_global(bi_floodedArea, 1) * 1.0e6_r8 / 1.0e9_r8   ! First recovered to m^2, then converted to thousand km^2 (total flooded area [including channel area] ).
                budget_glb_inund(bi_mainChnlArea, 1) = budget_global(bi_mainChnlArea, 1) * 1.0e6_r8 / 1.0e9_r8 ! First recovered to m^2, then converted to thousand km^2 (total channel surface area).
             end if

             budget_glb_inund(br_landOutflow, :) = budget_global(br_landOutflow, :) * 1.0e6_r8 / 1000._r8 ! First recovered to m^3/s, then converted to thousand m^3/s (total streamflow [flow rate] from land to oceans).

             budget_glb_inund(bVelo_downward, :) = budget_global(bVelo_downward, :) * 1.0e6_r8 ! Recovered to m/s (sum of all downward flow velocities).
             budget_glb_inund(bVelo_downChnlNo, :) = budget_global(bVelo_downChnlNo, :) * 1.0e6_r8 ! Recovered to number (total number of main channels with downward flow velocities).
             budget_glb_inund(bVelo_upward, :) = budget_global(bVelo_upward, :) * 1.0e6_r8     ! Recovered to m/s (sum of all upward flow velocities).
             budget_glb_inund(bVelo_upChnlNo, :) = budget_global(bVelo_upChnlNo, :) * 1.0e6_r8 ! Recovered to number (total number of main channels with upward flow velocities).

             budget_glb_inund(bmask_lndErr, :) = budget_global(bmask_lndErr, :) * 1.0e6_r8     ! Recovered to number (number of land cells not matched).
             budget_glb_inund(bmask_outErr, :) = budget_global(bmask_outErr, :) * 1.0e6_r8     ! Recovered to number (number of outlet cells not matched).
             budget_glb_inund(bmask_ocnErr, :) = budget_global(bmask_ocnErr, :) * 1.0e6_r8     ! Recovered to number (number of ocean cells not matched).

             budget_glb_inund(bmask_uLndrLnd, :) = budget_global(bmask_uLndrLnd, :) * 1.0e6_r8 ! Recovered to number (Tunit is land, rtmCTL is land).
             budget_glb_inund(bmask_uLndrOut, :) = budget_global(bmask_uLndrOut, :) * 1.0e6_r8 ! Recovered to number (Tunit is land, rtmCTL is outlet).
             budget_glb_inund(bmask_uLndrOcn, :) = budget_global(bmask_uLndrOcn, :) * 1.0e6_r8 ! Recovered to number (Tunit is land, rtmCTL is ocean).

             budget_glb_inund(bmask_uOutrLnd, :) = budget_global(bmask_uOutrLnd, :) * 1.0e6_r8 ! Recovered to number (Tunit is outlet, rtmCTL is land).
             budget_glb_inund(bmask_uOutrOut, :) = budget_global(bmask_uOutrOut, :) * 1.0e6_r8 ! Recovered to number (Tunit is outlet, rtmCTL is outlet).
             budget_glb_inund(bmask_uOutrOcn, :) = budget_global(bmask_uOutrOcn, :) * 1.0e6_r8 ! Recovered to number (Tunit is outlet, rtmCTL is ocean).

             budget_glb_inund(bmask_uOcnrLnd, :) = budget_global(bmask_uOcnrLnd, :) * 1.0e6_r8 ! Recovered to number (Tunit is ocean, rtmCTL is land).
             budget_glb_inund(bmask_uOcnrOut, :) = budget_global(bmask_uOcnrOut, :) * 1.0e6_r8 ! Recovered to number (Tunit is ocean, rtmCTL is outlet).
             budget_glb_inund(bmask_uOcnrOcn, :) = budget_global(bmask_uOcnrOcn, :) * 1.0e6_r8 ! Recovered to number (Tunit is ocean, rtmCTL is ocean).

             !--- Debug ---
             write(iulog,'(a)') ' '
             write(iulog,'(a)') '---------------------------------'
             write(iulog,'(a)') '   Compare Tunit%mask() and rtmCTL%mask() :'
             write(iulog,'(a)') '---------------------------------'

             write(iulog,'(a, f22.6)') '   Tunit%mask() is land, but rtmCTL%mask() is not land, number     =', budget_glb_inund(bmask_lndErr, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is outlet, but rtmCTL%mask() is not outlet, number =', budget_glb_inund(bmask_outErr, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is ocean, but rtmCTL%mask() is not ocean, number   =', budget_glb_inund(bmask_ocnErr, 1)

             write(iulog,'(a, f22.6)') '   Tunit%mask() is land, rtmCTL%mask() is land, number     =', budget_glb_inund(bmask_uLndrLnd, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is land, rtmCTL%mask() is outlet, number   =', budget_glb_inund(bmask_uLndrOut, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is land, rtmCTL%mask() is ocean, number    =', budget_glb_inund(bmask_uLndrOcn, 1)

             write(iulog,'(a, f22.6)') '   Tunit%mask() is outlet, rtmCTL%mask() is land, number   =', budget_glb_inund(bmask_uOutrLnd, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is outlet, rtmCTL%mask() is outlet, number =', budget_glb_inund(bmask_uOutrOut, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is outlet, rtmCTL%mask() is ocean, number  =', budget_glb_inund(bmask_uOutrOcn, 1)

             write(iulog,'(a, f22.6)') '   Tunit%mask() is ocean, rtmCTL%mask() is land, number    =', budget_glb_inund(bmask_uOcnrLnd, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is ocean, rtmCTL%mask() is outlet, number  =', budget_glb_inund(bmask_uOcnrOut, 1)
             write(iulog,'(a, f22.6)') '   Tunit%mask() is ocean, rtmCTL%mask() is ocean, number   =', budget_glb_inund(bmask_uOcnrOcn, 1)
             !--- Debug ---

             do nt = 1, nt_rtm

                write(iulog,'(a)') ' '
                write(iulog,'(a, i2)') 'tracer =', nt

                write (tracerID, '(a7, i1, a1)') '(tracer', nt, ')'

                write(iulog,'(a)') ' '
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                write(iulog,'(a)') trim(tracerID)//'   Surface-water balance check:'
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                budget_input = budget_glb_inund(br_qsur, nt) + budget_glb_inund(br_qsub, nt) + budget_glb_inund(br_qgwl, nt) + budget_glb_inund(br_qdto, nt)
                budget_output = budget_glb_inund(br_ocnout, nt) + budget_glb_inund(br_direct, nt)
                budget_other = budget_glb_inund(br_erolpn,nt) - budget_glb_inund(br_erolcn,nt)   !('previous MOSART sub-step channel outflow volume'-'current MOSART sub-step channel outflow volume'. --Inund.)                                                                                                                                                                                                

                if ( abs( budget_input - budget_output ) .lt. 1e-9_r8 ) then      ! 1e-9 km^3 = 1 m^3
                   write(iulog,'(a)') trim(tracerID)//'   Total input equals total output to oceans.'
                else
                   budget_diff = (budget_glb_inund(bv_volt_f, nt) - budget_glb_inund(bv_volt_i, nt) - budget_input + budget_output - budget_other)/abs(budget_input - budget_output) * 100._r8
                   write(iulog,'(a, f22.6)') trim(tracerID)//' * Surface-water balance error(%)             =', budget_diff
                end if

                write(iulog,'(a, f22.6)') trim(tracerID)//' * Total volume change (km^3)                 =', budget_glb_inund(bv_volt_f, nt) - budget_glb_inund(bv_volt_i, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//' * Total input minus output to oceans (km^3)  =', budget_input - budget_output + budget_other
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Input surface runoff (km^3)                =', budget_glb_inund(br_qsur, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Input sub-surface runoff (km^3)            =', budget_glb_inund(br_qsub, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Input from glacier, wetland or lake (km^3) =', budget_glb_inund(br_qgwl, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Input direct-to-ocean runoff (km^3)        =', budget_glb_inund(br_qdto, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Output flows to oceans (km^3)              =', budget_glb_inund(br_ocnout, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Output direct to oceans (km^3)             =', budget_glb_inund(br_direct, nt)

                write(iulog,'(a)') ' '
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                write(iulog,'(a)') trim(tracerID)//'   Initial surface water volumes (before coupling period):'
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                write(iulog,'(a, f22.6)') trim(tracerID)//' * Total surface water volume (km^3)          =', budget_glb_inund(bv_volt_i, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume over hillslopes (km^3)              =', budget_glb_inund(bv_wh_i, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume in subnetworks (km^3)               =', budget_glb_inund(bv_wt_i, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume in main channels (km^3)             =', budget_glb_inund(bv_wr_i, nt)

                ! If inundation scheme is on & the 1st tracer (liquid water) :
                if ( inundflag .and. Tctl%OPT_inund .eq. 1 .and. nt .eq. 1 ) then
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume over fps before cpl period (km^3)   =', budget_glb_inund(bv_fp_i, nt)
                end if

                write(iulog,'(a)') ' '
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                write(iulog,'(a)') trim(tracerID)//'   Final surface water volumes (after coupling period):'
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                write(iulog,'(a, f22.6)') trim(tracerID)//' * Total surface water volume (km^3)          =', budget_glb_inund(bv_volt_f, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume over hillslopes (km^3)              =', budget_glb_inund(bv_wh_f, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume in subnetworks (km^3)               =', budget_glb_inund(bv_wt_f, nt)
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume in main channels (km^3)             =', budget_glb_inund(bv_wr_f, nt)

                ! If inundation scheme is on & the 1st tracer (liquid water) :
                if ( inundflag .and. Tctl%OPT_inund .eq. 1 .and. nt .eq. 1 ) then
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Volume over floodplains (km^3)             =', budget_glb_inund(bv_fp_f, nt)
                end if

                ! If inundation scheme is on & the 1st tracer (liquid water) :
                if ( inundflag .and. Tctl%OPT_inund .eq. 1 .and. nt .eq. 1 ) then

                   write(iulog,'(a)') ' '
                   write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                   write(iulog,'(a)') trim(tracerID)//'   Floodplain water balance check:'
                   write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                   if ( abs( budget_glb_inund(bv_chnl2fp, nt) - budget_glb_inund(bv_fp2chnl, nt) ) .lt. 1e-9_r8 ) then      ! 1e-9 km^3 = 1 m^3
                      write(iulog,'(a)') trim(tracerID)//'   Total channel-to-floodplain flow equals total floodplain-to-channel flow.'
                   else
                      budget_diff = (budget_glb_inund(bv_fp_f, nt) - budget_glb_inund(bv_fp_i, nt) - budget_glb_inund(bv_chnl2fp, nt) + budget_glb_inund(bv_fp2chnl, nt))&
                           /abs(budget_glb_inund(bv_chnl2fp, nt) - budget_glb_inund(bv_fp2chnl, nt)) * 100._r8
                      write(iulog,'(a, f22.6)') trim(tracerID)//' * Floodplain water balance error(%)          =', budget_diff
                   end if

                   write(iulog,'(a, f22.6)') trim(tracerID)//' * Total change of floodplain volume (km^3)   =', budget_glb_inund(bv_fp_f, nt) - budget_glb_inund(bv_fp_i, nt)
                   write(iulog,'(a, f22.6)') trim(tracerID)//' * Total input minus output (km^3)            =', budget_glb_inund(bv_chnl2fp, nt) - budget_glb_inund(bv_fp2chnl, nt)
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Channel-to-floodplain flow (km^3)          =', budget_glb_inund(bv_chnl2fp, nt)
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Floodplain-to-channel flow (km^3)          =', budget_glb_inund(bv_fp2chnl, nt)

                   write(iulog,'(a)') ' '
                   write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                   write(iulog,'(a)') trim(tracerID)//'   Inundation statistics:'
                   write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Percentage of area flooded (incl. channel) =', budget_glb_inund(bi_floodedArea, nt)/budget_glb_inund(bi_landArea, nt)*100._r8
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Percentage of main channel surface area    =', budget_glb_inund(bi_mainChnlArea, nt)/budget_glb_inund(bi_landArea, nt)*100._r8
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Total land area (k*km^2)                   =', budget_glb_inund(bi_landArea, nt)
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Total flooded area (incl. channel; k*km^2) =', budget_glb_inund(bi_floodedArea, nt)
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Total main channel surface area (k*km^2)   =', budget_glb_inund(bi_mainChnlArea, nt)

                end if

                write(iulog,'(a)') ' '
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'
                write(iulog,'(a)') trim(tracerID)//'   Main channel flow statistics:'
                write(iulog,'(a)') trim(tracerID)//'---------------------------------'

                ! Number of channels with downward velocities (this is mean value for sub-steps in a coupling period):
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Number of channels with downward velocities=', budget_glb_inund(bVelo_downChnlNo, nt)/float(nsub)
                if ( budget_glb_inund(bVelo_downChnlNo, nt) .gt. 0 ) then
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Average downward flow velocity (m/s)       =', budget_glb_inund(bVelo_downward, nt)/budget_glb_inund(bVelo_downChnlNo, nt)
                else
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Average downward flow velocity (m/s)       =', 0._r8
                end if

                ! Number of channels with upward velocities (this is mean value for sub-steps in a coupling period):
                write(iulog,'(a, f22.6)') trim(tracerID)//'   Number of channels with upward velocities  =', budget_glb_inund(bVelo_upChnlNo, nt)/float(nsub)
                if ( budget_glb_inund(bVelo_upChnlNo, nt) .gt. 0 ) then
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Average upward flow velocity (m/s)         =', budget_glb_inund(bVelo_upward, nt)/budget_glb_inund(bVelo_upChnlNo, nt)
                else
                   write(iulog,'(a, f22.6)') trim(tracerID)//'   Average upward flow velocity (m/s)         =', 0._r8
                end if

                write(iulog,'(a, f22.6)') trim(tracerID)//' * Total land-to-ocean streamflow (k*m^3/s)   =', budget_glb_inund(br_landOutflow, nt)

             end do   ! end of 'do nt = 1, nt_rtm'

             write(iulog,'(a)') ' '
             write(iulog,'(a)') '=================================================== '
             write(iulog,'(a)') ' '
          endif
       endif   ! (End of if (masterproc). --Inund.)

       call t_stopf('mosartr_budget')
    endif  ! budget_write   (end of if (budget_check .and. budget_write). --Inund.)

    !-----------------------------------
    ! Write out MOSART history file
    !-----------------------------------

    call t_startf('mosartr_hbuf')
    call RtmHistFldsSet()
    call RtmHistUpdateHbuf()
    call t_stopf('mosartr_hbuf')

    call t_startf('mosartr_htapes')
    call RtmHistHtapesWrapup( rstwr, nlend )
    call t_stopf('mosartr_htapes')

    !-----------------------------------
    ! Write out MOSART restart file
    !-----------------------------------

    if (rstwr) then
       call t_startf('mosartr_rest')
       filer = RtmRestFileName(rdate=rdate)
       call RtmRestFileWrite( filer, rdate=rdate )
       call t_stopf('mosartr_rest')
    end if

    !-----------------------------------
    ! Done
    !-----------------------------------

    first_call = .false.

    call shr_sys_flush(iulog)
    call t_stopf('mosartr_tot')

  end subroutine Rtmrun

!-----------------------------------------------------------------------

  subroutine RtmFloodInit(frivinp, begr, endr, fthresh, evel )

    !-----------------------------------------------------------------------
    ! Uses

    ! Input variables
    character(len=*), intent(in) :: frivinp
    integer , intent(in)  :: begr, endr
    real(r8), intent(out) :: fthresh(begr:endr)
    real(r8), intent(out) :: evel(begr:endr,nt_rtm) 

    ! Local variables
    real(r8) , pointer :: rslope(:)   
    real(r8) , pointer :: max_volr(:)
    integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
    integer :: nt,n,cnt              ! indices
    logical :: readvar               ! read variable in or not
    integer :: ier                   ! status variable
    integer :: dids(2)               ! variable dimension ids 
    type(file_desc_t)  :: ncid       ! pio file desc
    type(var_desc_t)   :: vardesc    ! pio variable desc 
    type(io_desc_t)    :: iodesc     ! pio io desc
    character(len=256) :: locfn      ! local file name

    !MOSART Flood variables for spatially varying celerity
    real(r8) :: effvel(nt_rtm) = 0.7_r8   ! downstream velocity (m/s)
    real(r8) :: min_ev(nt_rtm) = 0.35_r8  ! minimum downstream velocity (m/s)
    real(r8) :: fslope = 1.0_r8           ! maximum slope for which flooding can occur
    character(len=*),parameter :: subname = '(RtmFloodInit) '
    !-----------------------------------------------------------------------

    allocate(rslope(begr:endr), max_volr(begr:endr), stat=ier)
    if (ier /= 0) call shr_sys_abort(subname // ' allocation ERROR')

    ! Assume that if SLOPE is on river input dataset so is MAX_VOLR and that
    ! both have the same io descriptor

    call getfil(frivinp, locfn, 0 )
    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ier = pio_inq_varid(ncid, 'SLOPE', vardesc)
    if (ier /= PIO_noerr) then
       if (masterproc) write(iulog,*) subname//' variable SLOPE is not on dataset'
       readvar = .false.
    else
       readvar = .true.
    end if
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
    if (readvar) then
       ier = pio_inq_vardimid(ncid, vardesc, dids)
       allocate(compdof(rtmCTL%lnumr))
       cnt = 0
       do n = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          compDOF(cnt) = rtmCTL%gindex(n)
       enddo
       call pio_initdecomp(pio_subsystem, pio_double, dids, compDOF, iodesc)
       deallocate(compdof)
! tcraig, there ia bug here, shouldn't use same vardesc for two different variable
       call pio_read_darray(ncid, vardesc, iodesc, rslope, ier)
       call pio_read_darray(ncid, vardesc, iodesc, max_volr, ier)
       call pio_freedecomp(ncid, iodesc)
    else
       rslope(:)   = 1._r8
       max_volr(:) = spval
    end if
    call pio_closefile(ncid)

    do nt = 1,nt_rtm
       do n = rtmCTL%begr, rtmCTL%endr
          fthresh(n) = 0.95*max_volr(n)*max(1._r8,rslope(n))
          ! modify velocity based on gridcell average slope (Manning eqn)
          evel(n,nt) = max(min_ev(nt),effvel(nt_rtm)*sqrt(max(0._r8,rslope(n)))) 
       end do
    end do

    deallocate(rslope, max_volr)

  end subroutine RtmFloodInit 

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
  subroutine MOSART_init
!
! !REVISION HISTORY:
! Author: Hongyi Li

! !DESCRIPTION:
! initialize MOSART variables
! 
! !USES:
! !ARGUMENTS:
  implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
  type(file_desc_t)  :: ncid       ! pio file desc
  type(var_desc_t)   :: vardesc    ! pio variable desc 
  type(io_desc_t)    :: iodesc_dbl ! pio io desc
  type(io_desc_t)    :: iodesc_int ! pio io desc
  integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
  integer :: ndims                 ! number of dimensions in the input
  integer :: dids(2)               ! variable dimension ids 
  integer :: dsizes(2)             ! variable dimension lengths
  integer :: ier                   ! error code
  integer :: begr, endr, iunit, nn, n, cnt, nr, nt, i, j
  integer :: numDT_r, numDT_t
  integer :: lsize, gsize
  integer :: igrow, igcol, iwgt, idim
  type(mct_aVect) :: avtmp, avtmpG ! temporary avects
  type(mct_aVect) :: avsrc, avdst  ! temporary
  type(mct_sMat)  :: sMat          ! temporary sparse matrix, needed for sMatP
  real(r8):: areatot_prev, areatot_tmp, areatot_new
  real(r8):: hlen_max, rlen_min
  integer :: tcnt
  character(len=16384) :: rList             ! list of fields for SM multiply
  character(len=1000) :: fname
  character(len=*),parameter :: subname = '(MOSART_init)'
  character(len=*),parameter :: FORMI = '(2A,2i10)'
  character(len=*),parameter :: FORMR = '(2A,2g15.7)'
  !real(r8), pointer          :: e_eprof_in2(:,:)
 
  begr = rtmCTL%begr
  endr = rtmCTL%endr
  
  if(endr >= begr) then
     ! routing parameters
     call ncd_pio_openfile (ncid, trim(frivinp_rtm), 0)
     call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
     allocate(compdof(rtmCTL%lnumr))
     cnt = 0
     do n = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        compDOF(cnt) = rtmCTL%gindex(n)
     enddo

     ! setup iodesc based on frac dids
     ier = pio_inq_varid(ncid, 'frac', vardesc)
     if (isgrid2d) then
        ndims = 2
     else
        ndims = 1
     endif
     ier = pio_inq_vardimid(ncid, vardesc, dids(1:ndims))
     do idim = 1, ndims
        ier = pio_inq_dimlen(ncid, dids(idim),dsizes(idim))
     enddo
     call pio_initdecomp(pio_subsystem, pio_double, dsizes(1:ndims), compDOF, iodesc_dbl)
     call pio_initdecomp(pio_subsystem, pio_int   , dsizes(1:ndims), compDOF, iodesc_int)

     deallocate(compdof)
     call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

     allocate(TUnit%euler_calc(nt_rtm))
     Tunit%euler_calc = .true.

     allocate(TUnit%frac(begr:endr))
     ier = pio_inq_varid(ncid, 'frac', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%frac, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read frac ',minval(Tunit%frac),maxval(Tunit%frac)
     call shr_sys_flush(iulog)
     
     if (wrmflag) then
       allocate(TUnit%domainfrac(begr:endr))
       ier = pio_inq_varid(ncid, 'domainfrac', vardesc)
       call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%domainfrac, ier)
       if (masterproc) write(iulog,FORMR) trim(subname),' read domainfrac ',minval(Tunit%domainfrac),maxval(Tunit%domainfrac)
       call shr_sys_flush(iulog)
     endif
     
     ! read fdir, convert to mask
     ! fdir <0 ocean, 0=outlet, >0 land
     ! tunit mask is 0=ocean, 1=land, 2=outlet for mosart calcs

     allocate(TUnit%mask(begr:endr))  
     ier = pio_inq_varid(ncid, 'fdir', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%mask, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read fdir mask ',minval(Tunit%mask),maxval(Tunit%mask)
     call shr_sys_flush(iulog)

     do n = rtmCtl%begr, rtmCTL%endr
        if (Tunit%mask(n) < 0) then
           Tunit%mask(n) = 0
        elseif (Tunit%mask(n) == 0) then
           Tunit%mask(n) = 2
           if (abs(Tunit%frac(n)-1.0_r8)>1.0e-9) then
              write(iulog,*) subname,' ERROR frac ne 1.0',n,Tunit%frac(n)
              call shr_sys_abort(subname//' ERROR frac ne 1.0')
           endif
        elseif (Tunit%mask(n) > 0) then
           Tunit%mask(n) = 1
           if (abs(Tunit%frac(n)-1.0_r8)>1.0e-9) then
              write(iulog,*) subname,' ERROR frac ne 1.0',n,Tunit%frac(n)
              call shr_sys_abort(subname//' ERROR frac ne 1.0')
           endif
        else
           call shr_sys_abort(subname//' Tunit mask error')
        endif
       
       if (wrmflag) then       
        if (Tunit%domainfrac(n) == 0) then
          Tunit%domainfrac(n) = 1
        elseif (Tunit%domainfrac(n) < 0) then
          write(iulog,*) subname,' ERROR domain frac < 0',n,Tunit%domainfrac(n)
          call shr_sys_abort(subname//' Tunit domainfrac error')
        endif
       endif
     enddo

     allocate(TUnit%ID0(begr:endr))  
     ier = pio_inq_varid(ncid, 'ID', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%ID0, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read ID0 ',minval(Tunit%ID0),maxval(Tunit%ID0)
     call shr_sys_flush(iulog)

     allocate(TUnit%dnID(begr:endr))  
     ier = pio_inq_varid(ncid, 'dnID', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%dnID, ier)
     if (masterproc) write(iulog,FORMI) trim(subname),' read dnID ',minval(Tunit%dnID),maxval(Tunit%dnID)
     call shr_sys_flush(iulog)

     !-------------------------------------------------------
     ! RESET ID0 and dnID indices using the IDkey to be consistent
     ! with standard gindex order to leverage gsmap_r
     !-------------------------------------------------------
     do n=rtmCtl%begr, rtmCTL%endr
        TUnit%ID0(n)  = IDkey(TUnit%ID0(n))
        if (Tunit%dnID(n) > 0 .and. TUnit%dnID(n) <= rtmlon*rtmlat) then
           if (IDkey(TUnit%dnID(n)) > 0 .and. IDkey(TUnit%dnID(n)) <= rtmlon*rtmlat) then
              TUnit%dnID(n) = IDkey(TUnit%dnID(n))
           else
              write(iulog,*) subname,' ERROR bad IDkey for TUnit%dnID',n,TUnit%dnID(n),IDkey(TUnit%dnID(n))
              call shr_sys_abort(subname//' ERROR bad IDkey for TUnit%dnID')
           endif
        endif
     enddo

     allocate(TUnit%area(begr:endr))  
     ier = pio_inq_varid(ncid, 'area', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%area, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read area ',minval(Tunit%area),maxval(Tunit%area)
     call shr_sys_flush(iulog)

     do n=rtmCtl%begr, rtmCTL%endr
        if (TUnit%area(n) <= 0._r8) TUnit%area(n) = rtmCTL%area(n)
        if (TUnit%area(n) /= rtmCTL%area(n)) then
           write(iulog,*) subname,' ERROR area mismatch',TUnit%area(n),rtmCTL%area(n)
           call shr_sys_abort(subname//' ERROR area mismatch')
        endif
     enddo

     allocate(TUnit%areaTotal(begr:endr))  
     ier = pio_inq_varid(ncid, 'areaTotal', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%areaTotal, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read areaTotal ',minval(Tunit%areaTotal),maxval(Tunit%areaTotal)
     call shr_sys_flush(iulog)

     allocate(TUnit%rlenTotal(begr:endr))
     TUnit%rlenTotal = 0._r8

     allocate(TUnit%nh(begr:endr))  
     ier = pio_inq_varid(ncid, 'nh', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nh, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read nh ',minval(Tunit%nh),maxval(Tunit%nh)
     call shr_sys_flush(iulog)

     allocate(TUnit%hslp(begr:endr))  
     ier = pio_inq_varid(ncid, 'hslp', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%hslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read hslp ',minval(Tunit%hslp),maxval(Tunit%hslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%hslpsqrt(begr:endr))  
     TUnit%hslpsqrt = 0._r8

     allocate(TUnit%gxr(begr:endr))  
     ier = pio_inq_varid(ncid, 'gxr', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%gxr, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read gxr ',minval(Tunit%gxr),maxval(Tunit%gxr)
     call shr_sys_flush(iulog)

     allocate(TUnit%hlen(begr:endr))
     TUnit%hlen = 0._r8

     allocate(TUnit%tslp(begr:endr))  
     ier = pio_inq_varid(ncid, 'tslp', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%tslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read tslp ',minval(Tunit%tslp),maxval(Tunit%tslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%tslpsqrt(begr:endr))  
     TUnit%tslpsqrt = 0._r8

     allocate(TUnit%tlen(begr:endr))
     TUnit%tlen = 0._r8

     allocate(TUnit%twidth(begr:endr))  
     ier = pio_inq_varid(ncid, 'twid', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%twidth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read twidth ',minval(Tunit%twidth),maxval(Tunit%twidth)
     call shr_sys_flush(iulog)

     allocate(TUnit%nt(begr:endr))  
     ier = pio_inq_varid(ncid, 'nt', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nt, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read nt ',minval(Tunit%nt),maxval(Tunit%nt)
     call shr_sys_flush(iulog)

     allocate(TUnit%rlen(begr:endr))  
     ier = pio_inq_varid(ncid, 'rlen', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rlen, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rlen ',minval(Tunit%rlen),maxval(Tunit%rlen)
     call shr_sys_flush(iulog)

     allocate(TUnit%rslp(begr:endr))  
     ier = pio_inq_varid(ncid, 'rslp', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rslp, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rslp ',minval(Tunit%rslp),maxval(Tunit%rslp)
     call shr_sys_flush(iulog)

     allocate(TUnit%rslpsqrt(begr:endr))  
     TUnit%rslpsqrt = 0._r8

     allocate(TUnit%rwidth(begr:endr))  
     ier = pio_inq_varid(ncid, 'rwid', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rwidth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rwidth ',minval(Tunit%rwidth),maxval(Tunit%rwidth)
     call shr_sys_flush(iulog)

     if (inundflag) then
        do n = rtmCtl%begr, rtmCTL%endr
           if ( rtmCTL%mask(n) .eq. 1 .or. rtmCTL%mask(n) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

              ! If Channel area >= unit area * 0.7 (note: 0.7 may be changed) :
              if ( TUnit%rwidth(n) * TUnit%rlen(n) .ge. TUnit%area(n) * 0.7_r8 ) then
                 TUnit%rwidth(n) = TUnit%area(n) * 0.7_r8 / TUnit%rlen(n)
              endif

           end if
        end do
     end if

     allocate(TUnit%rwidth0(begr:endr))  
     ier = pio_inq_varid(ncid, 'rwid0', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rwidth0, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rwidth0 ',minval(Tunit%rwidth0),maxval(Tunit%rwidth0)
     call shr_sys_flush(iulog)

     allocate(TUnit%rdepth(begr:endr))  
     ier = pio_inq_varid(ncid, 'rdep', vardesc)
     call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%rdepth, ier)
     if (masterproc) write(iulog,FORMR) trim(subname),' read rdepth ',minval(Tunit%rdepth),maxval(Tunit%rdepth)
     call shr_sys_flush(iulog)

     ! define outlets and relevant parameters where ocn rof two-way coupling is on
     if ( use_ocn_rof_two_way ) then
        allocate(TUnit%ocn_rof_coupling_ID(begr:endr))
        ier = pio_inq_varid(ncid, 'ocn_rof_coupling_ID', vardesc)
        call pio_read_darray(ncid, vardesc, iodesc_int, TUnit%ocn_rof_coupling_ID, ier)
        if (masterproc) write(iulog,FORMR) trim(subname),' read ocn_rof_coupling_ID',minval(Tunit%ocn_rof_coupling_ID),maxval(Tunit%ocn_rof_coupling_ID)
        call shr_sys_flush(iulog)

        allocate(TUnit%vdatum_conversion(begr:endr))
        ier = pio_inq_varid(ncid, 'vdatum_conversion', vardesc)
        call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%vdatum_conversion, ier)
        if (masterproc) write(iulog,FORMR) trim(subname),' read vdatum_conversion',minval(Tunit%vdatum_conversion),maxval(Tunit%vdatum_conversion)
        call shr_sys_flush(iulog)
     end if

     allocate(TUnit%nr(begr:endr))
  
     if (inundflag) then
        ! Calculate channel Manning roughness coefficients :
        call calc_chnlMannCoe ( )
     else
        !!allocate(TUnit%nr(begr:endr))   !(Repetitive, removed on 6-1-17. --Inund.)
        ier = pio_inq_varid(ncid, 'nr', vardesc)
        call pio_read_darray(ncid, vardesc, iodesc_dbl, TUnit%nr, ier)
        if (masterproc) write(iulog,FORMR) trim(subname),' read nr ',minval(Tunit%nr),maxval(Tunit%nr)
        call shr_sys_flush(iulog)
     endif

     if (inundflag) then

        !----------------------------------   
        ! Check input parameters :
        !----------------------------------   

       do n = rtmCtl%begr, rtmCTL%endr
          !if ( Tunit%mask(n) .eq. 1 .or. Tunit%mask(n) .eq. 2 ) then    ! 1-- Land; 2-- Basin outlet.
          if ( rtmCTL%mask(n) .eq. 1 .or. rtmCTL%mask(n) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

            if ( TUnit%area(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%area(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%area(n) <= 0 ')
            end if

            if ( TUnit%areaTotal(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%areaTotal(n) <= 0 for n=', n

              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%areaTotal(n) <= 0 ')
            end if

            if ( TUnit%nh(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%nh(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%nh(n) <= 0 ')
            end if

            if ( TUnit%hslp(n) .LT. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%hslp(n) < 0 for n=', n

              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%hslp(n) < 0 ')
            end if

            if ( TUnit%gxr(n) .LT. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%gxr(n) < 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%gxr(n) < 0 ')
            end if

            if ( TUnit%tslp(n) .LT. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%tslp(n) < 0 for n=', n

              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%tslp(n) < 0 ')
            end if

            if ( TUnit%twidth(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%twidth(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%twidth(n) <= 0 ')
            end if

            if ( TUnit%nt(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%nt(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%nt(n) <= 0 ')
            end if

            if ( TUnit%rlen(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%rlen(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%rlen(n) <= 0 ')
            end if

            if ( TUnit%rslp(n) .LT. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%rslp(n) < 0 for n=', n

              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%rslp(n) < 0 ')
            end if

            if ( TUnit%rwidth(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%rwidth(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%rwidth(n) <= 0 ')
            end if

            if ( TUnit%rdepth(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%rdepth(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%rdepth(n) <= 0 ')
            end if

            if ( TUnit%nr(n) .le. 0._r8 ) then
              write( iulog, * ) trim( subname ) // ' ERROR: TUnit%nr(n) <= 0 for n=', n
              call shr_sys_abort( trim( subname ) // ' ERROR: TUnit%nr(n) <= 0 ')
            end if

          end if
        end do

     end if

     if (inundflag) then
       if ( Tctl%RoutingMethod == DW ) then       ! Use diffusion wave method in channel routing computation.
          allocate (TUnit%rlen_dstrm(begr:endr))
          allocate (TUnit%rslp_dstrm(begr:endr))

          ! --------------------------------- 
          ! retrieve downstream values (TUnit%rlen_dstrm(:) and TUnit%rslp_dstrm(:)) for DW routing.
          ! --------------------------------- 

          call mct_aVect_init(avsrc,rList='rlen:rslp',lsize=rtmCTL%lnumr)
          call mct_aVect_init(avdst,rList='rlen:rslp',lsize=rtmCTL%lnumr)
          call mct_aVect_zero(avsrc)
          call mct_aVect_zero(avdst)
          cnt = 0
          do nr = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             avsrc%rAttr(1,cnt) = TUnit%rlen(nr)
             avsrc%rAttr(2,cnt) = TUnit%rslp(nr)
          enddo

          call mct_sMat_avMult(avsrc, sMatP_dnstrm, avdst)

          cnt = 0
          do nr = rtmCTL%begr,rtmCTL%endr
             cnt = cnt + 1
             TUnit%rlen_dstrm(nr) = avdst%rAttr(1,cnt)
             TUnit%rslp_dstrm(nr) = avdst%rAttr(2,cnt)
          enddo

          call mct_aVect_clean(avsrc)
          call mct_aVect_clean(avdst)

          if (masterproc) write(iulog,FORMR) trim(subname),' set rlen_dstrm ',minval(Tunit%rlen_dstrm),maxval(Tunit%rlen_dstrm)
          call shr_sys_flush(iulog)
          if (masterproc) write(iulog,FORMR) trim(subname),' set rslp_dstrm ',minval(Tunit%rslp_dstrm),maxval(Tunit%rslp_dstrm)
          call shr_sys_flush(iulog)
       end if

       if (Tctl%OPT_inund == 1) then
          allocate (TUnit%wr_bf(begr:endr))
          TUnit%wr_bf = 0.0_r8   

          allocate( TUnit%e_eprof_in2( Tctl%npt_elevProf, begr:endr ) )    

          ! --------------------------------- 
          ! (assign elevation values to TUnit%e_eprof_in2( :, : ) ).

          if ( Tctl%OPT_elevProf .eq. 1 ) then 
            call read_elevation_profile(ncid, 'ele', TUnit%e_eprof_in2)
          endif
          ! --------------------------------- 

          allocate (TUnit%a_eprof(begr:endr,12))
          TUnit%a_eprof = 0.0_r8

          allocate (TUnit%e_eprof(begr:endr,12))
          TUnit%e_eprof = 0.0_r8

          allocate (TUnit%a_chnl(begr:endr))
          TUnit%a_chnl = 0.0_r8

          allocate (TUnit%e_chnl(begr:endr))
          TUnit%e_chnl = 0.0_r8
 
          allocate (TUnit%ipt_bl_bktp(begr:endr))
          TUnit%ipt_bl_bktp = 0

          allocate (TUnit%a_eprof3(begr:endr, 13))
          TUnit%a_eprof3 = 0.0_r8

          allocate (TUnit%e_eprof3(begr:endr, 13))
          TUnit%e_eprof3 = 0.0_r8

          allocate (TUnit%npt_eprof3(begr:endr))
          TUnit%npt_eprof3 = 0

          allocate (TUnit%s_eprof3(begr:endr, 13))
          TUnit%s_eprof3 = 0.0_r8

          allocate (TUnit%alfa3(begr:endr, 13))
          TUnit%alfa3 = 0.0_r8

          allocate (TUnit%p3(begr:endr, 13))
          TUnit%p3 = 0.0_r8

          allocate (TUnit%q3(begr:endr, 13))
          TUnit%q3 = 0.0_r8      

          ! Pre-process elevation-profile parameters :
          call preprocess_elevProf ( )
       endif

     end if  ! inundflag
  
     ! initialize water states and fluxes
     allocate (TRunoff%wh(begr:endr,nt_rtm))
     TRunoff%wh = 0._r8

     allocate (TRunoff%dwh(begr:endr,nt_rtm))
     TRunoff%dwh = 0._r8

     allocate (TRunoff%yh(begr:endr,nt_rtm))
     TRunoff%yh = 0._r8

     allocate (TRunoff%qsur(begr:endr,nt_rtm))
     TRunoff%qsur = 0._r8

     allocate (TRunoff%qsub(begr:endr,nt_rtm))
     TRunoff%qsub = 0._r8

     allocate (TRunoff%qgwl(begr:endr,nt_rtm))
     TRunoff%qgwl = 0._r8
      
     allocate (TRunoff%qdem(begr:endr,nt_rtm)) 
     TRunoff%qdem = 0._r8
  
     allocate (TRunoff%ehout(begr:endr,nt_rtm))
     TRunoff%ehout = 0._r8

     allocate (TRunoff%tarea(begr:endr,nt_rtm))
     TRunoff%tarea = 0._r8

     allocate (TRunoff%wt(begr:endr,nt_rtm))
     TRunoff%wt= 0._r8

     allocate (TRunoff%dwt(begr:endr,nt_rtm))
     TRunoff%dwt = 0._r8

     allocate (TRunoff%yt(begr:endr,nt_rtm))
     TRunoff%yt = 0._r8

     allocate (TRunoff%mt(begr:endr,nt_rtm))
     TRunoff%mt = 0._r8

     allocate (TRunoff%rt(begr:endr,nt_rtm))
     TRunoff%rt = 0._r8

     allocate (TRunoff%pt(begr:endr,nt_rtm))
     TRunoff%pt = 0._r8

     allocate (TRunoff%vt(begr:endr,nt_rtm))
     TRunoff%vt = 0._r8

     allocate (TRunoff%tt(begr:endr,nt_rtm))
     TRunoff%tt = 0._r8

     allocate (TRunoff%etin(begr:endr,nt_rtm))
     TRunoff%etin = 0._r8

     allocate (TRunoff%etout(begr:endr,nt_rtm))
     TRunoff%etout = 0._r8

     allocate (TRunoff%rarea(begr:endr,nt_rtm))
     TRunoff%rarea = 0._r8

     allocate (TRunoff%wr(begr:endr,nt_rtm))
     TRunoff%wr = 0._r8

     allocate (TRunoff%dwr(begr:endr,nt_rtm))
     TRunoff%dwr = 0._r8

     allocate (TRunoff%yr(begr:endr,nt_rtm))
     TRunoff%yr = 0._r8

     allocate (TRunoff%mr(begr:endr,nt_rtm))
     TRunoff%mr = 0._r8

     allocate (TRunoff%rr(begr:endr,nt_rtm))
     TRunoff%rr = 0._r8

     allocate (TRunoff%pr(begr:endr,nt_rtm))
     TRunoff%pr = 0._r8

     allocate (TRunoff%vr(begr:endr,nt_rtm))
     TRunoff%vr = 0._r8

     allocate (TRunoff%tr(begr:endr,nt_rtm))
     TRunoff%tr = 0._r8

     allocate (TRunoff%erlg(begr:endr,nt_rtm))
     TRunoff%erlg = 0._r8

     allocate (TRunoff%erlateral(begr:endr,nt_rtm))
     TRunoff%erlateral = 0._r8

     allocate (TRunoff%erin(begr:endr,nt_rtm))
     TRunoff%erin = 0._r8

     allocate (TRunoff%erout(begr:endr,nt_rtm))
     TRunoff%erout = 0._r8

     allocate (TRunoff%eroup_lagi(begr:endr,nt_rtm))
     TRunoff%eroup_lagi = 0._r8

     allocate (TRunoff%eroup_lagf(begr:endr,nt_rtm))
     TRunoff%eroup_lagf = 0._r8

     allocate (TRunoff%erowm_regi(begr:endr,nt_rtm))
     TRunoff%erowm_regi = 0._r8

     allocate (TRunoff%erowm_regf(begr:endr,nt_rtm))
     TRunoff%erowm_regf = 0._r8

     allocate (TRunoff%eroutUp(begr:endr,nt_rtm))
     TRunoff%eroutUp = 0._r8

     allocate (TRunoff%eroutUp_avg(begr:endr,nt_rtm))
     TRunoff%eroutUp_avg = 0._r8

     allocate (TRunoff%erlat_avg(begr:endr,nt_rtm))
     TRunoff%erlat_avg = 0._r8

     allocate (TRunoff%ergwl(begr:endr,nt_rtm))
     TRunoff%ergwl = 0._r8

     allocate (TRunoff%flow(begr:endr,nt_rtm))
     TRunoff%flow = 0._r8
    
     allocate (TPara%c_twid(begr:endr))
     TPara%c_twid = 1.0_r8

     if ( Tctl%RoutingMethod == DW ) then       ! Use diffusion wave method in channel routing computation.
        allocate (TRunoff%rslp_energy(begr:endr))
        TRunoff%rslp_energy = 0.0_r8
        
        allocate (TRunoff%wr_dstrm(begr:endr, nt_rtm))
        TRunoff%wr_dstrm = 0.0_r8

        allocate (TRunoff%yr_dstrm(begr:endr))
        TRunoff%yr_dstrm = 0.0_r8    

        allocate (TRunoff%erin_dstrm(begr:endr,nt_rtm))
        TRunoff%erin_dstrm = 0.0_r8 

        do nr = rtmCTL%begr,rtmCTL%endr
            do i=1, rtmCTL%nUp(nr)
                n = rtmCTL%iUp(nr,i)            
                if(rtmCTL%iDown(n).ne.nr) then
                   write(iulog,*) trim(subname),' MOSART ERROR: upstream-downstream relationships ',nr
                   call shr_sys_abort(trim(subname)//' ERROR upstream-downstream relationships incorrect')
                end if
            enddo
        enddo
     end if
   
     if (inundflag) then
        !allocate (TRunoff%wr_ini(begr:endr))
        !TRunoff%wr_ini = 0.0

        !allocate (TRunoff%yr_ini(begr:endr))
        !TRunoff%yr_ini = 0.0

        allocate (TRunoff%wr_exchg(begr:endr))
        TRunoff%wr_exchg = 0.0_r8

        allocate (TRunoff%yr_exchg(begr:endr))
        TRunoff%yr_exchg = 0.0_r8

        if ( Tctl%RoutingMethod == DW ) then       ! Use diffusion wave method in channel routing computation.
           allocate (TRunoff%wr_exchg_dstrm(begr:endr))
           TRunoff%wr_exchg_dstrm = 0.0_r8

           allocate (TRunoff%yr_exchg_dstrm(begr:endr))
           TRunoff%yr_exchg_dstrm = 0.0_r8    
        end if

        !allocate (TRunoff%delta_wr(begr:endr))
        !TRunoff%delta_wr = 0.0

        allocate (TRunoff%wr_rtg(begr:endr))
        TRunoff%wr_rtg = 0.0_r8    

        allocate (TRunoff%yr_rtg(begr:endr))
        TRunoff%yr_rtg = 0.0_r8

        if ( Tctl%OPT_inund == 1 ) then

           allocate (TRunoff%wf_ini(begr:endr))
           TRunoff%wf_ini = 0.0_r8

           allocate (TRunoff%hf_ini(begr:endr))
           TRunoff%hf_ini = 0.0_r8
           
           allocate (TRunoff%ff_ini(begr:endr))
           TRunoff%ff_ini = 0.0_r8

           allocate (TRunoff%ffunit_ini(begr:endr))
           TRunoff%ffunit_ini = 0.0_r8
           allocate (TRunoff%netchange(begr:endr))
           TRunoff%netchange = 0.0_r8
           allocate (TRunoff%se_rf(begr:endr))
           TRunoff%se_rf = 0.0_r8

           allocate (TRunoff%ff_fp(begr:endr))
           TRunoff%ff_fp = 0.0_r8

           allocate (TRunoff%fa_fp(begr:endr))
           TRunoff%fa_fp = 0.0_r8

           allocate (TRunoff%wf_exchg(begr:endr))
           TRunoff%wf_exchg = 0.0_r8    

           allocate (TRunoff%hf_exchg(begr:endr))
           TRunoff%hf_exchg = 0.0_r8

           allocate (TRunoff%ff_unit(begr:endr))
           TRunoff%ff_unit = 0.0_r8
        endif
     end if
     
     if(heatflag) then
        ! initialize heat states and fluxes
        allocate (THeat%forc_t(begr:endr))
        THeat%forc_t = 273.15_r8
        allocate (THeat%forc_pbot(begr:endr))
        THeat%forc_pbot = 0._r8
        allocate (THeat%forc_vp(begr:endr))
        THeat%forc_vp = 0._r8
        allocate (THeat%forc_wind(begr:endr))
        THeat%forc_wind = 0._r8
        allocate (THeat%forc_lwrad(begr:endr))
        THeat%forc_lwrad = 0._r8
        allocate (THeat%forc_solar(begr:endr))
        THeat%forc_solar = 0._r8

        allocate (THeat%Tqsur(begr:endr))
        THeat%Tqsur = 273.15_r8
        allocate (THeat%Tqsub(begr:endr))
        THeat%Tqsub = 273.15_r8

        allocate (THeat%Tt(begr:endr))
        THeat%Tt = 273.15_r8
        allocate (THeat%Ha_h2t(begr:endr))
        THeat%Ha_h2t = 0._r8
        allocate (THeat%Ha_t2r(begr:endr))
        THeat%Ha_t2r = 0._r8
        allocate (THeat%Ha_lateral(begr:endr))
        THeat%Ha_lateral = 0._r8
        allocate (THeat%Hs_t(begr:endr))
        THeat%Hs_t = 0._r8
        allocate (THeat%Hl_t(begr:endr))
        THeat%Hl_t = 0._r8
        allocate (THeat%He_t(begr:endr))
        THeat%He_t = 0._r8
        allocate (THeat%Hh_t(begr:endr))
        THeat%Hh_t = 0._r8
        allocate (THeat%Hc_t(begr:endr))
        THeat%Hc_t = 0._r8
        allocate (THeat%deltaH_t(begr:endr))
        THeat%deltaH_t = 0._r8
        allocate (THeat%deltaM_t(begr:endr))
        THeat%deltaM_t = 0._r8

        allocate (THeat%Tr(begr:endr))
        THeat%Tr = 273.15_r8
        allocate (THeat%Ha_rin(begr:endr))
        THeat%Ha_rin = 0._r8
        allocate (THeat%Ha_rout(begr:endr))
        THeat%Ha_rout = 0._r8
        allocate (THeat%Ha_eroutUp(begr:endr))
        THeat%Ha_eroutUp = 0._r8
        allocate (THeat%Ha_eroutUp_avg(begr:endr))
        THeat%Ha_eroutUp_avg = 0._r8
        allocate (THeat%Ha_erlat_avg(begr:endr))
        THeat%Ha_erlat_avg = 0._r8
        allocate (THeat%Hs_r(begr:endr))
        THeat%Hs_r = 0._r8
        allocate (THeat%Hl_r(begr:endr))
        THeat%Hl_r = 0._r8
        allocate (THeat%He_r(begr:endr))
        THeat%He_r = 0._r8
        allocate (THeat%Hh_r(begr:endr))
        THeat%Hh_r = 0._r8
        allocate (THeat%Hc_r(begr:endr))
        THeat%Hc_r = 0._r8
        allocate (THeat%deltaH_r(begr:endr))
        THeat%deltaH_r = 0._r8
        allocate (THeat%deltaM_r(begr:endr))
        THeat%deltaM_r = 0._r8

        allocate (THeat%Tt_avg(begr:endr))
        THeat%Tt_avg = 273.15_r8
        allocate (THeat%Tr_avg(begr:endr))
        THeat%Tr_avg = 273.15_r8
        
       ! read the parameters for mosart-heat
        if(endr >= begr) then
            allocate(TPara%t_alpha(begr:endr))    
           TPara%t_alpha = 27.19_r8
            allocate(TPara%t_beta(begr:endr))
           TPara%t_beta = 13.63_r8
            allocate(TPara%t_gamma(begr:endr))
           TPara%t_gamma = 0.1576_r8
             allocate(TPara%t_mu(begr:endr))
           TPara%t_mu = 0.5278_r8
        end if
        
        allocate (THeat%coszen(begr:endr))
        THeat%coszen = 0._r8
        
     end if
    
     call pio_freedecomp(ncid, iodesc_dbl)
     call pio_freedecomp(ncid, iodesc_int)
     call pio_closefile(ncid)

   ! control parameters and some other derived parameters
   ! estimate derived input variables

     if (inundflag .and. Tctl%OPT_inund == 1) then
        do iunit = rtmCTL%begr, rtmCTL%endr
          if ( rtmCTL%mask(iunit) .eq. 1 .or. rtmCTL%mask(iunit) .eq. 3 ) then   ! 1--Land; 3--Basin outlet (downstream is ocean).

            ! Main channel storage capacity :
            TUnit%wr_bf( iunit ) = TUnit%rwidth( iunit ) * TUnit%rlen( iunit ) * TUnit%rdepth( iunit )

          end if 
        end do
     end if

     do iunit=rtmCTL%begr,rtmCTL%endr
        if(TUnit%Gxr(iunit) > 0._r8) then
           TUnit%rlenTotal(iunit) = TUnit%area(iunit)*TUnit%Gxr(iunit)
        end if
     end do

     do iunit=rtmCTL%begr,rtmCTL%endr
        if(TUnit%rlen(iunit) > TUnit%rlenTotal(iunit)) then
           TUnit%rlenTotal(iunit) = TUnit%rlen(iunit)
        end if
     end do     

     do iunit=rtmCTL%begr,rtmCTL%endr
      
        if(TUnit%rlen(iunit) > 0._r8) then
           TUnit%hlen(iunit) = TUnit%area(iunit) / TUnit%rlenTotal(iunit) / 2._r8
           hlen_max = max(1000.0_r8, sqrt(TUnit%area(iunit))) ! constrain the hillslope length
           if(TUnit%hlen(iunit) > hlen_max) then
              TUnit%hlen(iunit) = hlen_max   ! allievate the outlier in drainage density estimation. TO DO
           end if
           rlen_min = sqrt(TUnit%area(iunit))
           if(TUnit%rlen(iunit) < rlen_min) then
              TUnit%tlen(iunit) = TUnit%area(iunit) / rlen_min / 2._r8 - TUnit%hlen(iunit)
           else
              TUnit%tlen(iunit) = TUnit%area(iunit) / TUnit%rlen(iunit) / 2._r8 - TUnit%hlen(iunit)
           end if
  
           if(TUnit%twidth(iunit) < 0._r8) then
              TUnit%twidth(iunit) = 0._r8
           end if
           if(TUnit%tlen(iunit) > 0._r8) then
              if ((TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit) > 1._r8) then
                  TUnit%twidth(iunit) = TPara%c_twid(iunit)*TUnit%twidth(iunit)*((TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit))
              end if
           end if
          
           if(TUnit%tlen(iunit) > 0._r8 .and. TUnit%twidth(iunit) <= 0._r8) then
              TUnit%twidth(iunit) = 0._r8
           end if
        else
           TUnit%rlen(iunit) = 0._r8
           TUnit%hlen(iunit) = 0._r8
           TUnit%tlen(iunit) = 0._r8
           TUnit%twidth(iunit) = 0._r8
        end if
        
        if(TUnit%rslp(iunit) <= 0._r8) then

        if (inundflag) then
           TUnit%rslp(iunit) = Tctl%rslp_assume
        else
           TUnit%rslp(iunit) = 0.0001_r8
        endif

        end if
        if(TUnit%tslp(iunit) <= 0._r8) then
           TUnit%tslp(iunit) = 0.0001_r8
        end if
        if(TUnit%hslp(iunit) <= 0._r8) then
           TUnit%hslp(iunit) = 0.005_r8
        end if
        TUnit%rslpsqrt(iunit) = sqrt(Tunit%rslp(iunit))
        TUnit%tslpsqrt(iunit) = sqrt(Tunit%tslp(iunit))
        TUnit%hslpsqrt(iunit) = sqrt(Tunit%hslp(iunit))
     end do
  end if  ! endr >= begr

  ! retrieve the downstream channel attributes after some post-processing above
  if (Tctl%RoutingMethod == DW ) then       ! Use diffusion wave method in channel routing computation.
     allocate (TUnit%rlen_dstrm(begr:endr))
     allocate (TUnit%rslp_dstrm(begr:endr))

     ! --------------------------------- 
     ! Need code to retrieve values of TUnit%rlen_dstrm(:) and TUnit%rslp_dstrm(:) .
     ! --------------------------------- 
       call mct_aVect_zero(avsrc_dnstrm)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avsrc_dnstrm%rAttr(1,cnt) = TUnit%rlen(iunit)
       enddo
       call mct_aVect_zero(avdst_dnstrm)
       call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          TUnit%rlen_dstrm(iunit) = avdst_dnstrm%rAttr(1,cnt)
       enddo

       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          avsrc_dnstrm%rAttr(1,cnt) = TUnit%rslp(iunit)
       enddo
       call mct_aVect_zero(avdst_dnstrm)
       call mct_sMat_avMult(avsrc_dnstrm, sMatP_dnstrm, avdst_dnstrm)
       cnt = 0
       do iunit = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          TUnit%rslp_dstrm(iunit) = avdst_dnstrm%rAttr(1,cnt)
       enddo

     if (masterproc) write(iulog,FORMR) trim(subname),' set rlen_dstrm ',minval(Tunit%rlen_dstrm),maxval(Tunit%rlen_dstrm)
     call shr_sys_flush(iulog)
     if (masterproc) write(iulog,FORMR) trim(subname),' set rslp_dstrm ',minval(Tunit%rslp_dstrm),maxval(Tunit%rslp_dstrm)
     call shr_sys_flush(iulog)
  end if

  !--- compute areatot from area using dnID ---
  !--- this basically advects upstream areas downstream and
  !--- adds them up as it goes until all upstream areas are accounted for

  allocate(Tunit%areatotal2(rtmCTL%begr:rtmCTL%endr))
  Tunit%areatotal2 = 0._r8

  ! initialize avdst to local area and add that to areatotal2
  cnt = 0
  call mct_avect_zero(avdst_upstrm)
  do nr = rtmCTL%begr,rtmCTL%endr
     cnt = cnt + 1
     avdst_upstrm%rAttr(1,cnt) = rtmCTL%area(nr)
     Tunit%areatotal2(nr) = avdst_upstrm%rAttr(1,cnt)
  enddo

  tcnt = 0
  areatot_prev = -99._r8
  areatot_new = -50._r8
  do while (areatot_new /= areatot_prev .and. tcnt < rtmlon*rtmlat)

     tcnt = tcnt + 1

     ! copy avdst to avsrc for next downstream step
     cnt = 0
     call mct_avect_zero(avsrc_upstrm)
     do nr = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        avsrc_upstrm%rAttr(1,cnt) = avdst_upstrm%rAttr(1,cnt)
     enddo

     call mct_avect_zero(avdst_upstrm)

     call mct_sMat_avMult(avsrc_upstrm, sMatP_upstrm, avdst_upstrm)

     ! add avdst to areatot and compute new global sum
     cnt = 0
     areatot_prev = areatot_new
     areatot_tmp = 0._r8
     do nr = rtmCTL%begr,rtmCTL%endr
        cnt = cnt + 1
        Tunit%areatotal2(nr) = Tunit%areatotal2(nr) + avdst_upstrm%rAttr(1,cnt)
        areatot_tmp = areatot_tmp + Tunit%areatotal2(nr)
     enddo
     call shr_mpi_sum(areatot_tmp, areatot_new, mpicom_rof, 'areatot_new', all=.true.)

  enddo

  if (areatot_new /= areatot_prev) then
     write(iulog,*) trim(subname),' MOSART ERROR: areatot incorrect ',areatot_new, areatot_prev
     call shr_sys_abort(trim(subname)//' ERROR areatot incorrect')
  endif

  call SubTimestep ! prepare for numerical computation

  call shr_mpi_max(maxval(Tunit%numDT_r),numDT_r,mpicom_rof,'numDT_r',all=.false.)
  call shr_mpi_max(maxval(Tunit%numDT_t),numDT_t,mpicom_rof,'numDT_t',all=.false.)
  if (masterproc) then
     write(iulog,*) subname,' DLevelH2R = ',Tctl%DlevelH2R
     write(iulog,*) subname,' numDT_r   = ',minval(Tunit%numDT_r),maxval(Tunit%numDT_r)
     write(iulog,*) subname,' numDT_r max  = ',numDT_r
     write(iulog,*) subname,' numDT_t   = ',minval(Tunit%numDT_t),maxval(Tunit%numDT_t)
     write(iulog,*) subname,' numDT_t max  = ',numDT_t
  endif 

  end subroutine MOSART_init

!----------------------------------------------------------------------------

  subroutine read_elevation_profile(ncid, varname, e_eprof_in2)

    implicit none
    type(file_desc_t)      :: ncid       ! pio file desc
    character(len=*)       :: varname    ! variable name
    real(r8)               :: e_eprof_in2(:,:)

    character(len=*),parameter :: subname = '(read_elevation_profile)'

    type(var_desc_t)   :: vardesc    ! pio variable desc 
    logical            :: readvar    ! If variable exists or not
    type(io_desc_t)    :: iodesc     ! pio io desc
    integer            :: ndims      ! ndims for var      
    integer            :: dsizes(3)  ! dim sizes
    integer            :: dimids(3)  ! dim ids
    integer, pointer   :: compdof(:) ! computational degrees of freedom for pio 
    integer            :: begr, endr, cnt, m, n
    integer            :: ier
    integer            :: elesize    ! number of points for elevation profile
    real(r8), pointer  :: ele(:,:)
    character(len=2)   :: str

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    call check_var(ncid, varname, vardesc, readvar)

    if (readvar) then
      ier = pio_inq_varid(ncid, varname, vardesc)
      ier = pio_inq_varndims(ncid, vardesc, ndims)
      ier = pio_inq_vardimid(ncid, vardesc, dimids)

      do n = 1,ndims
        ier = pio_inq_dimlen(ncid,dimids(n),dsizes(n))
      enddo

      elesize = dsizes(ndims)
      if (elesize /= Tctl%npt_elevProf) then
        write(iulog,*) trim(subname),' MOSART ERROR: number of points in elevation profile is ', elesize
        call shr_sys_abort(trim(subname)//' ERROR number of points in elevation profile is not euqal to 11')
      endif

      allocate(compdof(rtmCTL%lnumr*elesize)) ! dims(ndims): 
      cnt = 0

      do n = 1, elesize
        do m = rtmCTL%begr,rtmCTL%endr
          cnt = cnt + 1
          compDOF(cnt) = (n-1)* rtmCTL%numr + rtmCTL%gindex(m)
        enddo
      enddo

      call pio_initdecomp(pio_subsystem, pio_double, dsizes(1:ndims), compDOF, iodesc)
      deallocate(compdof)

      allocate(ele(begr:endr,1:elesize))
      call pio_read_darray(ncid, vardesc, iodesc, ele, ier)

      do n = 1, elesize
        e_eprof_in2(n,:) = ele(:,n)
      enddo

      deallocate(ele)
      call pio_freedecomp(ncid, iodesc)

    else

      do n = 1, Tctl%npt_elevProf

        if (Tctl%npt_elevProf-1<10) then 
          write(str,'(I1)') Tctl%npt_elevProf-1
        elseif (Tctl%npt_elevProf<100) then
          write(str,'(I2)') Tctl%npt_elevProf-1
        endif

        call check_var(ncid, 'ele'//str, vardesc, readvar)

        if (readvar) then
          ier = pio_inq_varid(ncid, 'ele'//str, vardesc)
        else
          call shr_sys_abort(trim(subname)//' ERROR missing elevation profile data')
        endif

        if (n==1) then

          allocate(compdof(rtmCTL%lnumr))
          cnt = 0
          do m = rtmCTL%begr,rtmCTL%endr
            cnt = cnt + 1
            compDOF(cnt) = rtmCTL%gindex(m)
          enddo

          ier = pio_inq_varndims(ncid, vardesc, ndims)
          ier = pio_inq_vardimid(ncid, vardesc, dimids)

          do m = 1,ndims
            ier = pio_inq_dimlen(ncid,dimids(m),dsizes(m))
          enddo

          call pio_initdecomp(pio_subsystem, pio_double, dsizes(1:ndims), compDOF, iodesc)

        endif

        call pio_read_darray(ncid, vardesc, iodesc, e_eprof_in2(n,:), ier)

      enddo

      if (masterproc) then
         write(iulog,*) subname,'read elevattion profile successfully'
      endif 

    endif

  end subroutine read_elevation_profile

!----------------------------------------------------------------------------

  subroutine regeom
    
! Calculate reservoir layer average area (km2) 
 use shr_sys_mod   , only : shr_sys_flush
 use WRM_type_mod  , only : WRMUnit
 use RunoffMod     , only : rtmCTL
 use RtmVar         , only : iulog, ngeom, nlayers
    
 implicit none
 real(r8) :: M_W,M_L,gm_j,d_res,dd_in,C_a,C_v
 real(r8),dimension(ngeom+1) :: d_zi0,v_zti0,a_di0
 real(r8) ::dd_zz(ngeom),a_dd(ngeom+1),a_zi(ngeom+1),C_aa(ngeom+1), ar_f = 1.0e6                 ! Factor to convert area to m^2
 real(r8) :: pi      = 3.141593_r8      ! pi
 real(r8) :: d_s    = 1.0_r8        !0.60_r8                                                   ! Surface layer depth (m)
 real(r8) :: dv,da,dz
 real(r8) :: d_v(nlayers)                           ! Reservoir volume change at layer (m^3)
 real(r8) :: rho_z(nlayers)                           ! Reservoir layer density (kg/m^3), taken constant for furture revision 
 real(r8) :: delta_z                                ! depth change to calculate corresponding area/volume(m)
 integer :: i,j,k,iunit,damID                            ! indices
 real(r8) :: s_tin,s_lum,a_sur,fac                 ! Initial total storage (10*6m^3), Lumped storage for grids with multiple reservoirs (m^3)
 real(r8) :: A_cf,V_cf                        !Area and volume correcting factor for error from geometry estimation
 character(len=*),parameter :: subname = '(regeom)'        
                            
!**************************************************************************************************************************************************
    
 do iunit = rtmCTL%begr,rtmCTL%endr    
        
     damID = WRMUnit%INVicell(iunit)
     if (damID>0.0_r8 .and. WRMUnit%d_ns(damID) >= 1) then ! .and. WRMUnit%Depth(damID) > 0.0_r8 .and. WRMUnit%Height(damID) > 0.0_r8
          if (WRMUnit%Depth(damID) <= 2.0_r8)WRMUnit%Depth(damID) = 2._r8
          if (WRMUnit%Height(damID) <= 2.0_r8)WRMUnit%Height(damID) = WRMUnit%Depth(damID)
          WRMUnit%d_resrv(damID) = 0.95_r8*WRMUnit%Height(damID)
          WRMUnit%h_resrv(damID) = 0.95_r8*WRMUnit%Height(damID)
      
          s_lum = WRMUnit%StorCap(damID)
          s_tin = WRMUnit%V_str(damID)
          a_sur = WRMUnit%SurfArea(damID)
                
          if (s_lum <= 0.0_r8 .or. s_tin <= 0.0_r8 .or. ((s_lum/1.0e6) - s_tin)<0_r8) then
               fac = 0.0_r8
          elseif  (((s_lum/1.0e6) - s_tin)>=0_r8) then
               fac = ((s_lum/1.0e6) - s_tin)/s_tin
          endif
                
          !    Calculate layer depth    
          if (WRMUnit%Width_r(damID) <= 0.0_r8)WRMUnit%Width_r(damID) = 1._r8
          if (WRMUnit%Length_r(damID) <= 0.0_r8)WRMUnit%Length_r(damID) = 1._r8
                        
          ! Area and volume correcting factors for relative error as compared to GRanD
          if (WRMUnit%A_dfs(damID)>=0_r8) then
               A_cf = 1._r8 + (abs(WRMUnit%A_errs(damID))/100._r8)
          elseif(WRMUnit%A_dfs(damID)<0_r8) then
               A_cf = 1._r8 - (abs(WRMUnit%A_errs(damID))/100._r8)
          end if
          if (WRMUnit%V_dfs(damID)>=0_r8) then
               V_cf = 1._r8 + (abs(WRMUnit%V_errs(damID))/100._r8) + fac
          elseif(WRMUnit%V_dfs(damID)<0_r8) then
               V_cf = 1._r8 - (abs(WRMUnit%V_errs(damID))/100._r8) + fac
          end if
          if (A_cf<=0._r8)A_cf = 0.1_r8
          if (V_cf<=0._r8)V_cf = 0.1_r8
                
          M_W   = WRMUnit%Width_r(damID)
          M_L   = WRMUnit%Length_r(damID)
          if (M_W<1._r8)M_W = 1.0_r8
          if (M_L<1._r8)M_L = 1.0_r8
          gm_j  = WRMUnit%geometry(damID)
          d_res = WRMUnit%d_resrv(damID)        
          C_a   = WRMUnit%C_as(damID)
          C_v   = WRMUnit%C_vs(damID)
                
          if (WRMUnit%A_str(damID)<= 0._r8) WRMUnit%A_str(damID)= 1._r8
                
          if (WRMUnit%A_str(damID)<= 2._r8) then
               C_a  = 2.0_r8*WRMUnit%C_as(damID)
          end if
                            
          ! Uniform subsurface layer depth for initialization    and limit maximum layer thickness 
          dd_in = d_res/ngeom !bottom layers evenly descritized
                            
          ! Calculate reservoir geometry to establish depth-area-volume relationship            
                    
          ! Calculate depth area    
          !********** Curved Lake Bottom 
                
          do j = 1, ngeom!    
               if (gm_j == 1.0_r8) then    
                     a_dd(j) = C_a*M_L*M_W*(ar_f)*(1-((dd_in*(j-1))/d_res)**2._r8)
               else if (gm_j == 2.0_r8) then    
                     a_dd(j) = C_a*M_L*M_W*(ar_f)*(1-((dd_in*(j-1))/d_res)**2._r8)*(1-((dd_in*(j-1))/d_res))
               else if (gm_j == 3.0_r8) then    
                     a_dd(j) = C_a*M_L*M_W*(ar_f)*(1-((dd_in*(j-1))/d_res)**2._r8)*((d_res-(dd_in*(j-1)))/d_res)**0.5_r8
               else if (gm_j == 4.0_r8) then    
                     a_dd(j) = C_a*(2._r8/3._r8)*M_L*M_W*(ar_f)*(1-((dd_in*(j-1))/d_res)**2._r8)*(1-((dd_in*(j-1)))/d_res)
               else if (gm_j == 5.0_r8) then    
                     a_dd(j) = C_a*pi*0.25_r8*M_L*M_W*(ar_f)*(1-((dd_in*(j-1))/d_res)**2._r8)*((d_res-(dd_in*(j-1)))/d_res)**0.5_r8
               end if                    
          end do
                    
          a_dd(ngeom+1) = a_dd(ngeom)*0.001_r8            !Bottom area given non-zero value
                    
          ! Reverse indexing so that bottom is 1 and top is d_n+1 and convert to m2
          do j = 1,ngeom+1
               k =ngeom+2-j  
               a_di0(k) = ((a_dd(j)))
               if (a_di0(k)==0._r8) a_di0(k)=1._r8
               WRMUnit%a_di(damID,k)  = (A_cf*a_di0(k))     !Area corrected for error for optimal geometry
               if (WRMUnit%a_di(damID,k)==0._r8) WRMUnit%a_di(damID,k)=1._r8
          end do    
                
          ! Calculate layer depth,area,and volume from bottom 
          d_zi0(1) = 0._r8
          WRMUnit%d_zi(damID,1)  = d_zi0(1)
          do j = 2, ngeom+1    
               d_zi0(j) = d_zi0(j-1) + dd_in
               WRMUnit%d_zi(damID,j)  = d_zi0(j) 
          end do
                        
          ! Calculate layer average area,and total volume from bottom
          v_zti0(1) = (0.001_r8*0.5_r8*(WRMUnit%a_di(damID,1)+WRMUnit%a_di(damID,2))*dd_in)!    
          WRMUnit%v_zti(damID,1) = (V_cf*v_zti0(1))     !lower Volume corrected for error
          if (WRMUnit%v_zti(damID,1)==0._r8) WRMUnit%v_zti(damID,1)=1._r8
          do j = 2, ngeom+1     
               a_zi(j) = 0.5_r8*(WRMUnit%a_di(damID,j)+WRMUnit%a_di(damID,j-1)) !Area converted to m^2
               v_zti0(j) = v_zti0(j-1) + ((C_v*a_zi(j)*dd_in))
               WRMUnit%v_zti(damID,j) = (V_cf*v_zti0(j))     !Volume corrected for error
               if (WRMUnit%v_zti(damID,j)==0._r8) WRMUnit%v_zti(damID,j)=1._r8
               if (WRMUnit%a_di(damID,j)<0._r8 .or. WRMUnit%v_zti(damID,j)<0._r8) write(iulog,*) 'geom negative',iunit,WRMUnit%grandid(damID)
          end do
                
          !Initial reservoir storage 
          do j = WRMUnit%d_ns(damID),1,-1
               if (j == WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID) == 1) then
                     WRMUnit%dd_z(damID,j) = WRMUnit%d_resrv(damID)
               elseif (j == WRMUnit%d_ns(damID) .and. WRMUnit%d_ns(damID) > 1) then !top layer depth kept constant
                     WRMUnit%dd_z(damID,j) = d_s        !0.6_r8
               elseif ((WRMUnit%d_ns(damID)>1 .and.j < WRMUnit%d_ns(damID)).and.(WRMUnit%d_resrv(damID) - WRMUnit%dd_z(damID,WRMUnit%d_ns(damID)))>0._r8) then
                     WRMUnit%dd_z(damID,j) = (WRMUnit%d_resrv(damID) - WRMUnit%dd_z(damID,WRMUnit%d_ns(damID))) / (WRMUnit%d_ns(damID) - 1) !bottom layers evenly descritized
               endif
          end do    
            
          if (WRMUnit%d_ns(damID)>1 .and. WRMUnit%dd_z(damID,WRMUnit%d_ns(damID)-1)<d_s)then !layer thickness too small
               WRMUnit%d_ns(damID)=int((WRMUnit%d_resrv(damID)/d_s)+1) 
               do j=1,nlayers!WRMUnit%d_ns(damID)
                     WRMUnit%dd_z(damID,j) = 0._r8
               end do
                    
               ! Reinitialize layer thickness
               do j = WRMUnit%d_ns(damID),1,-1
                     if (j == WRMUnit%d_ns(damID)) then
                         WRMUnit%dd_z(damID,j) = d_s      !top layer depth
                     else
                         WRMUnit%dd_z(damID,j) = (WRMUnit%d_resrv(damID) - WRMUnit%dd_z(damID,WRMUnit%d_ns(damID)))/(WRMUnit%d_ns(damID) - 1) !bottom layers evenly descritized
                     end if
               end do    
          end if
        
          if (WRMUnit%d_ns(damID)>1) then
               WRMUnit%ddz_local(damID) = WRMUnit%dd_z(damID,WRMUnit%d_ns(damID)-1)
          else
               WRMUnit%ddz_local(damID) = WRMUnit%dd_z(damID,WRMUnit%d_ns(damID))
          end if
                
          ! Calculate layer depth (minimum at bottom)    
          WRMUnit%d_z(damID,1)=0._r8
          do j = 2, nlayers    
               if (j<=WRMUnit%d_ns(damID)+1) then
                     WRMUnit%d_z(damID,j) = WRMUnit%d_z(damID,j-1) + WRMUnit%dd_z(damID,j-1)
               else 
                     WRMUnit%d_z(damID,j) = 0._r8
               end if
          end do    
          ! Assign layer area and volume based on depth-area-volume relationship
          WRMUnit%a_d(damID,1)=(WRMUnit%a_di(damID,1))
          WRMUnit%v_zt(damID,1)=(WRMUnit%v_zti(damID,1))
          do i=2,nlayers
               if (i<=WRMUnit%d_ns(damID)+1) then
                     do j=2,ngeom+1
                          if (WRMUnit%d_z(damID,i)>WRMUnit%d_zi(damID,j-1).and.WRMUnit%d_z(damID,i)<=WRMUnit%d_zi(damID,j))then
                               delta_z = (WRMUnit%d_z(damID,i)-WRMUnit%d_zi(damID,j-1))/(WRMUnit%d_zi(damID,j)-WRMUnit%d_zi(damID,j-1))
                               WRMUnit%a_d(damID,i) = (delta_z*(WRMUnit%a_di(damID,j)-WRMUnit%a_di(damID,j-1))+WRMUnit%a_di(damID,j-1))
                               WRMUnit%v_zt(damID,i) = (delta_z*(WRMUnit%v_zti(damID,j)-WRMUnit%v_zti(damID,j-1))+WRMUnit%v_zti(damID,j-1))
                          else if (WRMUnit%d_z(damID,i) > WRMUnit%d_zi(damID,ngeom+1) .and. i<= WRMUnit%d_ns(damID)+1)then
                               delta_z = (WRMUnit%d_z(damID,i)-WRMUnit%d_zi(damID,ngeom))/(WRMUnit%d_zi(damID,ngeom+1)-WRMUnit%d_zi(damID,ngeom))
                               WRMUnit%a_d(damID,i) = (delta_z*(WRMUnit%a_di(damID,ngeom+1)-WRMUnit%a_di(damID,ngeom))+WRMUnit%a_di(damID,ngeom))
                               WRMUnit%v_zt(damID,i) = (delta_z*(WRMUnit%v_zti(damID,ngeom+1)-WRMUnit%v_zti(damID,ngeom))+WRMUnit%v_zti(damID,ngeom))
                          end if
                     end do
                     if (WRMUnit%v_zt(damID,i)<WRMUnit%v_zt(damID,i-1))write(iulog,*) 'check reservoir geometry for',WRMUnit%grandid(damID)
               else
                     WRMUnit%a_d(damID,i) = 0._r8
                     WRMUnit%v_zt(damID,i) = 0._r8
               end if
          end do
                                
          ! Calculate layer volume(m^3)
          do j = 1, nlayers    
               if (j<=WRMUnit%d_ns(damID)) then
                     WRMUnit%d_v(damID,j) = ((WRMUnit%v_zt(damID,j+1) - WRMUnit%v_zt(damID,j)))    
                     if (WRMUnit%d_v(damID,j)==0._r8) WRMUnit%d_v(damID,j)=1._r8    
                     WRMUnit%dd_z(damID,j) = WRMUnit%d_z(damID,j+1) - WRMUnit%d_z(damID,j)
                     if (WRMUnit%d_v(damID,j)<0.0_r8) then 
                         write(iulog,*) subname,'Layer volume negative:Check geometry data for dam',WRMUnit%grandid(damID)
                     end if
               else
                     WRMUnit%d_v(damID,j) = 0._r8
                     WRMUnit%dd_z(damID,j) = 0._r8
               end if
          end do        
                
          !     Intitialize layer temperature and total storage                    
          do j=1,WRMUnit%d_ns(damID)                    
               rho_z(j) = 1000._r8*( 1._r8 - 1.9549e-05*(abs(WRMUnit%temp_resrv(damID,j)-277._r8))**1.68_r8) 
               WRMUnit%v_zn(damID,j) = (WRMUnit%d_v(damID,j))
               WRMUnit%m_zo(damID,j) = (WRMUnit%d_v(damID,j)*rho_z(j)) 
               WRMUnit%m_zn(damID,j) = (WRMUnit%d_v(damID,j)*rho_z(j))
          end do
     end if
 end do    
 end subroutine regeom                    
!---------------------------------------------------------------------------- 

  subroutine SubTimestep
  ! !DESCRIPTION: predescribe the sub-time-steps for channel routing
    implicit none    
    integer :: iunit   !local index
    character(len=*),parameter :: subname = '(SubTimestep)'

    allocate(TUnit%numDT_r(rtmCTL%begr:rtmCTL%endr),TUnit%numDT_t(rtmCTL%begr:rtmCTL%endr))
    TUnit%numDT_r = 1
    TUnit%numDT_t = 1
    allocate(TUnit%phi_r(rtmCTL%begr:rtmCTL%endr),TUnit%phi_t(rtmCTL%begr:rtmCTL%endr))
    TUnit%phi_r = 0._r8
    TUnit%phi_t = 0._r8

    do iunit=rtmCTL%begr,rtmCTL%endr
       if(TUnit%mask(iunit) > 0 .and. TUnit%rlen(iunit) > 0._r8) then
          TUnit%phi_r(iunit) = TUnit%areaTotal2(iunit)*sqrt(TUnit%rslp(iunit))/(TUnit%rlen(iunit)*TUnit%rwidth(iunit))
          if(TUnit%phi_r(iunit) >= 10._r8) then
             TUnit%numDT_r(iunit) = (TUnit%numDT_r(iunit)*log10(TUnit%phi_r(iunit))*Tctl%DLevelR) + 1
          else 
             TUnit%numDT_r(iunit) = TUnit%numDT_r(iunit)*1.0_r8*Tctl%DLevelR + 1
          end if
       end if
       if(TUnit%numDT_r(iunit) < 1) TUnit%numDT_r(iunit) = 1
      
       if(TUnit%tlen(iunit) > 0._r8) then
          TUnit%phi_t(iunit) =      TUnit%area(iunit)*sqrt(TUnit%tslp(iunit))/(TUnit%tlen(iunit)*TUnit%twidth(iunit))
          if(TUnit%phi_t(iunit) >= 10._r8) then 
             TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*log10(TUnit%phi_t(iunit))*Tctl%DLevelR) + 1
          else 
             TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*1.0*Tctl%DLevelR) + 1
          end if
       end if
       if(TUnit%numDT_t(iunit) < 1) TUnit%numDT_t(iunit) = 1
    
    end do
  end subroutine SubTimestep

!-----------------------------------------------------------------------

end module RtmMod

