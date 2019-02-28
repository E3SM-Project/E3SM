!------------------------------------------------------------------------------
module m_xnldata
!------------------------------------------------------------------------------
!  module for computing the quadruplet interaction
!  Created by Gerbrant van Vledder
!
!  version 1.01   16/02/1999  Initial version
!          2.01   01/10/2001  various extensions added
!          3.1.01 01/10/2001  Array's for k4 -locus added
!          3.2    12/05/2002  Triplet data added
!          4.00   08/08/2002  Upgrade to version 4.0
!          4.01   19/08/2002  Various modifications for consistency reasons
!          5.01    9/09/2002  Length of strings aqname and bqname modified
!                             q_dstep added, step for BQF files
!                 11/09/2002  Filtering variables added
!          5.02   12/04/2003  Switch for triplet variables corrected
!          5.03   26/05/2003  Switch for lumping along locus added
!                 04/06/2003  Switch for Gauss-Legendre integration added
!                 06/06/2003  Switch iq_xdia added and NXDIA removed
!                 12/06/2003  Loop indices ik1,ia1,ik3,ia1 added
!                 16/06/2003  Switch IQ_SYM introduced
!                 04/09/2003  Version string set in subroutine q_version
!                 09/09/2003  Parameter id_facmax introduced
!------------------------------------------------------------------------------------
implicit none
!
character(len=60) q_version    ! version string
!
character(len=20) sub_name     ! Name of active subroutine
character(len=20) qbase        ! base name for I/O files
character(len=20) qf_error     ! name of file with error messages
!
integer iufind  ! Specifies handling of unit numbers, see Z_FILEIO
integer iscreen ! identifier for screen, set in XNL_INIT
!
!  unit numbers for I/O
!
integer luq_bqf  ! binary file storing and retrieving precomputed loci
integer luq_cfg  ! user defined configuration
integer luq_err  ! file with error messages
integer luq_fil  ! test output for filtering
integer luq_grd  ! ASCII file storing and retrieving precomputed loci
integer luq_int  ! test file for test output of integration
integer luq_loc  ! statistics about computed loci
integer luq_log  ! logging
integer luq_prt  ! general print file for quadruplets
integer luq_trf  ! testing transformation of loci
integer luq_tst  ! test file for quadruplets
integer luq_txt  ! reading (error) text file
integer luq_t13  ! test of basis integration
!------------------------------------------------------------------------------
!  physical coefficients, to be obtained through interface XNL_INIT
!------------------------------------------------------------------------------
real q_grav     ! gravitational acceleration (Earth = 9.81 m/s^2)
real qf_tail    ! power of spectral tail of E(f), e.g. -4,, -4.5, -5
!               ! these values must be set in the interface routine
!------------------------------------------------------------------------------
!  filtering coefficients
!------------------------------------------------------------------------------
real qf_krat     ! maximum ratio of the interacting wave numbers k1 and k3
real qf_dmax     ! maximum directional difference between k1 and k3
real qf_frac     ! fraction of maximum action density to filter
!
!  program switches, optionally to be reset in routine Q_SETCONFIG
!
integer iq_compact ! switch to compact data
!                 == 0, do not compact
!                 == 1, compact data by elimiting zero contribution along locus
!
integer iq_cple   ! type of coupling coefficient
!                 == 1, deep water coefficient of Webb
!                 == 2, deep water coefficient of Zakharov
!                 == 3, finite depth coefficient of Hasselmann & Herterich
!                 == 4, finite depth coefficient of Zakharov
!                 == 5, finite depth coefficient of Lin & Perrie
!
integer iq_disp   ! type of dispersion relation, viz. depth dependency
!                 == 1, deep water, possibly with geometric scaling
!                 == 2, linear dispersion relation, w^2 = g.k.tanh(kd)
!                 == 3, nonlinear dispersion relation
!
integer iq_dscale !  switch to activate depth scaling according to
                  !  Herterich and Hasselmann
!                 !  == 0, No depth scaling
!                 !  == 1, depth scaling activated
!
integer iq_filt   !  switch to activate filtering in wave number space
!                 !  ==0, no filtering
!                 !  ==1, filtering activated
!
integer iq_gauleg ! switch for Gauss-Legendre interpolation
!                 !  == 0, No Gauss-Legendre, default
!                 !  > 0   Gauss-Legendre, iq_gauleg is number of points
!
integer iq_geom   !  type of scaling
!                 == 0, no geometric scaling, only directional scaling of loci
!                 == 1, geometric scaling using Resio/Tracy method
!                       only possible in the case IQ_DISP=1
!
integer iq_grid   !  type of spectral grid
!                 == 1, sector & symmetric around zero
!                 == 2, sector & symmetric around zero & non-symmetric
!                 == 3, full circle & non-symmetric
!
integer iq_integ  ! option to output integration results
!                 !  ==0 no output of integration
!                 !  ==1 only sum per locus
!                 !  ==2 also information per point on locus
!                 !  ==3 only basic line integrals
!
integer iq_interp ! type of interpolation to retrieve action density
!                 !  == 1, bi-linear interpolation in discrete spectrum (default)
!                 !  == 2, take nearest bins, on the basis of maximum weight
!
integer iq_locus  !  Option for computation of locus
!                 !  ==1, explicit polar method with fixed k-step
!                 !  ==2, explicit polar method with adpative k-stepping
!                 !  ==3, explicit polar method with geometric k-spacing
!
integer iq_log    !  switch to activate logging to file QBASE//.LOG
!                 !  == 0, No print output
!                 !  == 1, print output
!
integer iq_lump   !  switch to activate lumping on locus
!                 !  == 0, No lumping
!                 !  == 1, Lumping along locus
!
integer iq_make   !  option to make quadruplet grid
!                    == 1, make when needed (default)
!                    == 2, always make quadruplet grid
!                    == 3, only make grid file
!
integer iq_mod    !  option to redistribute points on locus
!                 !  == 0, Points will be used as computed by tracing algortihm
!                 !  == 1, Equi-distant spacing on points along locus (NLOC1)
!
integer iq_prt    !  switch to activate print output, to file QBASE//.PRT
!                 !  == 0, No print output
!                 !  == 1, print output
!
integer iq_search ! switch to determine search for a proper grid
!                 == 0, no search is carried out
!                 == 1, search nearest (relative) interaction grid
!
integer iq_screen ! option to send output to the screen
!                 !  == 0, no output is send to screen
!                 !  == 1, output is send to screen
!
integer iq_sym    ! switch to activate use of symmetry reduction
!                 !  == 0, no symmetries are used
!                 !  == 1, symmetry activated (default)
!
integer iq_test   !  test level, output is directed to unit luqtst
!                 !  == 0, no test output
!                 !  == 1, output of basic I/O
!                 !  == 2, extensive test output
!
integer iq_trace  !  trace option
!                 !  == 0, no trace of subroutine calls
!                 !  > 0,  maximum number of traces per subroutine
!                 !  < 0,  as for >0 but now output is send to the screen
!
integer iq_trf    !  option to print transformed loci to special output file
!                 !  == 0, no output to data file unit luqtrf
!                 !  == 1, test output from routine Q_GETLOCUS
!
integer iq_t13    ! option to output T13 integration
!                 !  ==0, no output
!                 !  ==1, test output of T13 per locus
!
integer iq_xdia    ! switch to activate output to extended DIA data file
!                    == 0, no output
!                    >  0, output to data file, but only when lumping is also
!                          activated
!---------------------------------------------------------------------------------------
!
!
! grid administration
!
character(len=13) aqname       ! name of ASCII grid file
character(len=13) bqname       ! name of binary quadruplet grid file
character(len=13) lastquadfile ! name of last retrieved BQF file
character(len=21) q_header     ! header of Binary Quadruplet File as intended in BQF-file
character(len=21) r_header     ! header of Binary Quadruplet File as exists in BQF-file
logical lq_grid                ! flag to make (new) interaction grid
!
integer nkq     ! number of wave numbers of quad-grid
integer naq     ! number of angles of quad-grad
integer ncirc   ! number of angles on a full circle
!
integer ia_k1,ik_k1 ! indices of main loop variables
integer ia_k3,ik_k3 ! indices of main loop variables
!
real fqmin      ! lowest frequency in Hz
real fqmax      ! highest frequency in Hz
real q_sector   ! half plane width in degrees (for iq_grid=1,2)
real q_dstep    ! step size for generating BQF files
!
integer, parameter :: mq_stack=10 ! maximum number of elements in stack
!
integer mlocus   ! maximum number of points on locus for defining arrays
integer nlocus0  ! preferred number of points on locus
integer nlocus1  ! number of points on locus as computed in Q_CMPLOCUS
integer klocus   ! number of points on locus as stored in quadruplet database
                 ! based on nlocus0, iq_gauleg and iq_lump (without compacting)
                 ! used in Q_ALLOCATE to define size of data arrays
integer nlocus   ! number of points on locus, equal to klocus
integer nlocusx  ! number of points on locus for use in computation (nlocusx <= nlocus)
!
real kqmin       ! lowest wave number
real kqmax       ! highest wave number
real wk_max      ! maximum weight for wave number interpolation, set in Q_INIT
!
real k0x,k0y,dk0 ! components of initial wave number of locus,
real krefx,krefy ! components of reference wave number for quad-grid
real k1x,k1y     ! components of k1 wave number
real k2x,k2y     ! components of k2 wave number
real k3x,k3y     ! components of k3 wave number
real k4x,k4y     ! components of k4 wave number
real px,py       ! components of difference k1-k3 wave number
real pmag        ! magnitude of P-vector
real pang        ! angle related of P-vector, Pang = atan2(py,px), (radians)
real sang        ! angle of symmytry axis of locus, SANG = PANG +/ pi° (radians)
real xang        ! angle of locus for the case that w1=w3, Xang=atan2(-px,py), (radians)
real q           ! difference of radian frequencies, used in Resio-Tracy method
real kmin_loc    ! minimum wave number of locus along symmetry axis
real kmax_loc    ! maximum wave number of locus along symmetry axis
real kmid        ! wave number at midpoint of locus along symmetry axis
real kmidx       ! x-component of wave number at midpoint of locus along symmetry axis
real kmidy       ! y-component of wave number at midpoint of locus along symmetry axis
real loc_crf     ! circumference of locus in (kx,ky)-space
real loc_area    ! area of locus, measured in (kx-ky)- space
real loc_xz      ! x-coordinate of center of gravity of locus in (kx,ky)-space
real loc_yz      ! y-coordinate of center of gravity of locus in (kx,ky)-space
!
!  data for extended input k-grid, necessary when input grid is smaller than
!  internal k-grid.
!
! real fackx                            ! geometric spacing factor of input grid
! integer nkx                           ! new number of k-rings of extended input grid
! real, allocatable :: kx(:)            ! extended k-grid
! real, allocatable :: nspecx(:,:)      ! extended action density spectrum
!
!  information about pre_computed locus, only half the angles need to be saved
!
!
integer, allocatable :: quad_nloc(:,:)   ! number of points on locus
integer, allocatable :: quad_ik2(:,:,:)  ! lower wave number index of k2
integer, allocatable :: quad_ia2(:,:,:)  ! lower direction index of k2
integer, allocatable :: quad_ik4(:,:,:)  ! lower wave number index of k4
integer, allocatable :: quad_ia4(:,:,:)  ! lower direction index of k4
real, allocatable :: quad_w1k2(:,:,:)    ! weight 1 of k2
real, allocatable :: quad_w2k2(:,:,:)    ! weight 2 of k2
real, allocatable :: quad_w3k2(:,:,:)    ! weight 3 of k2
real, allocatable :: quad_w4k2(:,:,:)    ! weight 4 of k2
real, allocatable :: quad_w1k4(:,:,:)    ! weight 1 of k4
real, allocatable :: quad_w2k4(:,:,:)    ! weight 2 of k4
real, allocatable :: quad_w3k4(:,:,:)    ! weight 3 of k4
real, allocatable :: quad_w4k4(:,:,:)    ! weight 4 of k4
real, allocatable :: quad_zz  (:,:,:)    ! compound product of cple*ds*sym/jac
!
!  characteristic of computed locus
!
real, allocatable :: x2_loc(:)   ! k2x coordinates around locus
real, allocatable :: y2_loc(:)   ! k2y coordinates around locus
real, allocatable :: z_loc(:)    ! data value around locus
real, allocatable :: s_loc(:)    ! coordinate along locus
real, allocatable :: x4_loc(:)   ! k4x coordinates around locus
real, allocatable :: y4_loc(:)   ! k4y coordinates around locus
real, allocatable :: ds_loc(:)   ! step size around locus
real, allocatable :: jac_loc(:)  ! jacobian term around locus
real, allocatable :: cple_loc(:) ! coupling coefficient around locus
real, allocatable :: sym_loc(:)  ! factor for symmetry between k3 and k4
!
real, allocatable :: k_pol(:)  ! wave numbers during polar generation of locus
real, allocatable :: c_pol(:)  ! cosines during polar generation of locus
real, allocatable :: a_pol(:)  ! angles of polar locus
!
!  characteristics of modified locus, result
!
real, allocatable :: x2_mod(:)   ! k2x coordinates along locus
real, allocatable :: y2_mod(:)   ! k2y coordinates along locus
real, allocatable :: x4_mod(:)   ! k4x coordinates along locus
real, allocatable :: y4_mod(:)   ! k4y coordinates along locus
real, allocatable :: z_mod(:)    ! data value around locus
real, allocatable :: s_mod(:)    ! coordinate along locus
real, allocatable :: ds_mod(:)   ! step size around locus
real, allocatable :: jac_mod(:)  ! jacobian term around locus
real, allocatable :: cple_mod(:) ! coupling coefficient around locus
real, allocatable :: sym_mod(:)  ! factor for symmetry between k3 and k4
!
real, allocatable :: k2m_mod(:) ! k2 magnitude around locus
real, allocatable :: k2a_mod(:) ! k2 angle around locus
real, allocatable :: k4m_mod(:) ! k4 magnitude around locus
real, allocatable :: k4a_mod(:) ! k4 angle around locus
!
!   result of subroutine Q_weight
!
real, allocatable :: wk_k2(:)   ! position of k2 and k4 wave number
real, allocatable :: wk_k4(:)   ! w.r.t. discrete k-grid
real, allocatable :: wa_k2(:)   ! position of k2 and k4 wave number
real, allocatable :: wa_k4(:)   ! w.r.t. discrete a-grid
real, allocatable :: wt_k2(:)   ! weight factor in tail,
real, allocatable :: wt_k4(:)   ! wt==1 for wave numbers inside k-grid
!
integer, allocatable :: t_ik2(:)   ! transformed weight for k2-magnitude
integer, allocatable :: t_ia2(:)   ! transformed direction for k2
integer, allocatable :: t_ik4(:)   ! transformed tail factor for k2
integer, allocatable :: t_ia4(:)   ! transformed weight for k4
real, allocatable :: t_w1k2(:)  ! transformed weight 1 for k2
real, allocatable :: t_w2k2(:)  ! transformed weight 2 for k2
real, allocatable :: t_w3k2(:)  ! transformed weight 3 for k2
real, allocatable :: t_w4k2(:)  ! transformed weight 4 for k2
real, allocatable :: t_w1k4(:)  ! transformed weight 1 for k4
real, allocatable :: t_w2k4(:)  ! transformed weight 2 for k4
real, allocatable :: t_w3k4(:)  ! transformed weight 3 for k4
real, allocatable :: t_w4k4(:)  ! transformed weight 4 for k4
real, allocatable :: t_zz(:)    ! product term
!
!  corresponding declarations
!
integer, allocatable :: r_ik2(:)
integer, allocatable :: r_ia2(:)
integer, allocatable :: r_ik4(:)
integer, allocatable :: r_ia4(:)
real, allocatable :: r_w1k2(:),r_w2k2(:),r_w3k2(:),r_w4k2(:)
real, allocatable :: r_w1k4(:),r_w2k4(:),r_w3k4(:),r_w4k4(:)
real, allocatable :: r_zz(:),r_jac(:),r_cple(:),r_sym(:),r_ws(:)
!
real, allocatable :: dt13(:)    ! increment along locus
!
real, allocatable :: q_xk(:)    ! extended wave number array starting at index 0
real, allocatable :: q_sk(:)    ! step size of extended wave number array
real sk_max                     ! maximum wave number in extended array
!
real, allocatable :: q_k(:)     !  wave number grid [1/m]
real, allocatable :: q_dk(:)    !  width of wave number bins [1/m]
real, allocatable :: q_kpow(:)  !  wave number to a certain power, used in filtering
real, allocatable :: q_f(:)     !  frequencies accociated to wave number/depth
real, allocatable :: q_df(:)    !  step size of frequency grid
real, allocatable :: q_sig(:)   !  radian frequencies associated to wave number/depth
real, allocatable :: q_dsig(:)  !  step size of radian frequency grid
real, allocatable :: q_cg(:)    !  group velocity (m/s)
real, allocatable :: q_a(:)     !  directions of quadruplet grid in radians
real, allocatable :: q_ad(:)    !  directions of quadruplet grid in degrees
real, allocatable :: a(:,:)     !  Action density on wave number grid A(sigma,theta)
real, allocatable :: nspec(:,:) !  Action density on wave number grid N(kx,ky)
real, allocatable :: nk1d(:)    !  Internal 1d action density spectrum N(k)
real, allocatable :: qnl(:,:)   !  Nonlinear energy transfer Snl(k,theta)
!
integer id_facmax                   ! Factor for determining range of depth search (Q_SEARCHGRID)
real q_dird1,q_dird2                ! first and last direction of host model (via XNL_INIT) degrees
real q_depth                        ! local water depth in m
real q_maxdepth                     ! maximum water depth, set in XNL_INIT, used in Q_CTRGRID
real q_mindepth                     ! minimum water depth, set in XNL_INIT, used in Q_CTRGRID
real q_lambda                       ! geometric scaling factor for 'deep' water loci
real q_scale                        ! additional scale factor resulting from SEARCH for neasrest grid
!
real eps_q                          ! absolute accuracy for check of Q
real eps_k                          ! absolute accuracy for equality check of k
real rel_k                          ! relative accuracy for equality check of k
!
integer iq_stack                    ! Sequence number of stack with subroutine calls
character(len=21) cstack(mq_stack)  ! Stack with module names
!
!  characteristics of locus
!
real crf1     ! estimated circumference of locus
!---------------------------------------------------------------------------------
!
!  information about type of grid
!
integer iaref       !  index of first angle of reference wave numbers
integer iamax       !  maximum difference in indices for sector grids
integer iaq1,iaq2   !  indices of do-loop for directions
integer iag1,iag2   !  range of directions for precomputed interaction grid
real q_ang1,q_ang2  !  lower and upper angle of grid in degrees
real q_delta        !  directional spacing of angular grid in radians
real q_deltad       !  directional spacing of angular grid in degrees
!
real q_ffac         !  geometric factor between subsequent frequencies
real q_kfac         !  geometric factor between subsequent wave numbers
                    !  (only valid for IQ_IDISP==1)
real qk_tail        ! power of spectral tail of N(k), computed from qf_tail
!
!-----------------------------------------------------------------------------
!
!!/R real wq2(4)          ! interpolation weights for k2
!!/R real wq4(4)          ! interpolation weights for k4
!!/R real wqw             ! overall weight of contribution
!!/R real wtriq(40)       ! triplet weights
!!/R integer ikq2(4)      ! wave number index for k2
!!/R integer idq2(4)      ! angle index for k2
!!/R integer ikq4(4)      ! wave number index for k4
!!/R integer idq4(4)      ! angle index for k4
!!/R integer iktriq(40,3) ! k-indices of triplets
!!/R integer idtriq(40,3) ! direction indices of triplets

!
!============== General settings =================
!
integer iq_type   !  method for computing the nonlinear interactions
!                    depending on the value of iq_type a number of settings
!                    for other processes or schematizations are set in Q_COMPU
!  iq_type==1: deep water, symmetric spectrum, Webb coupling coefficient
!           2: deep water computation with WAM depth scaling based on Herterich
!              and Hasselmann (1980)
!           3: finite depth transfer
!
integer iq_err    !  counts the number of errors
!                    if no error occurred, IQ_ERR = 0
!                    for each occuring error, iq_err is incremented
!                    errors are always terminating
!                    routine Q_ERROR handles the reporting on the error
!
integer iq_warn   !  counts the number of warnings
!
!  indices for test output of actual integration
!  these values are set and optionally modified in Q_SETCONFIG
!
contains
!----------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine xnl_init(sigma,dird,nsigma,ndir,pftail,x_grav,depth,ndepth, &
&   iquad,iqgrid,iproc,ierror)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 4 September 2003
!   +---+ |   |  Release: 5.03
!         +---+
!
use m_fileio
use m_constants
! do not use m_xnldata
!
implicit none
!
!  0. Update history
!
!     10/01/2001   Initial version
!     14/02/2001   Release 3
!     06/11/2001   Depth forced to 1000 when iquad ==1,2
!     08/08/2002   Upgrade to version 4
!     19/08/2002   Extra test output
!     22/08/2002   User defined directions stored in Quad-system
!     26/08/2002   Minimum water depth in variable q_mindepth
!     09/09/2002   Initialized LASTQUADFILE
!     10/09/2002   Initialisation of Q_DSTEP
!     11/09/2002   Called of Q_ALLOC moved to location after Q_SETCFG
!                  Output unit luq_fil added
!     16/09/2002   Parameter IPROC added to take are of MPI
!     25/09/2002   Check added for directions of sector grid
!     25/04/2003   name q_alloc changed to q_allocate
!     04/06/2003   variable IQ_INT renamed IQ_INTEG
!     11/06/2003   Call to Q_SETCFG changed into Q_SETCONFIG
!                  Call to Q_CHKCFG changed into Q_CHKCONFIG
!                  Call to subroutine Q_SUMMARY added
!                  Compute size of points on locus, stored in KLOCUS
!     13/06/2003   Test parameters moved to Q_SETCONFIG
!     04/09/2003   Routine Q_SETVERSION added
!
!  1. Purpose:
!
!     Initialize coefficients, integration space, file i/o for computation
!     nonlinear quadruplet wave-wave interaction
!
!  2. Method
!
!     Set version number
!     Set unit unit numbers
!     Open quad related files
!     Optionally reset configuration by a back door option
!     Compute integration spaces for given water depths
!
!  3. Parameter list:
!
!Type   I/O              Name          Description
!------------------------------------------------------------------------------
integer, intent(in)  ::  nsigma        ! Number of sigma values
integer, intent(in)  ::  ndir          ! Number of directions
integer, intent(in)  ::  ndepth        ! Number of water depths
real, intent(in)     ::  sigma(nsigma) ! Radian frequencies
real, intent(in)     ::  dird(ndir)    ! Directions (degrees)
real, intent(in)     ::  pftail        ! power of spectral tail, e.g. -4 or -5
real, intent(in)     ::  depth(ndepth) ! depths for which integration space must be computed
real, intent(in)     ::  x_grav        ! gravitational acceleration
integer, intent(in)  ::  iquad         ! Type of method for computing nonlinear interactions
integer, intent(in)  ::  iqgrid        ! Type of grid for computing nonlinear interactions
integer, intent(in)  ::  iproc         ! Processor number, controls output file for MPI
integer, intent(out) ::  ierror        ! Error indicator. If no errors are detected IERR=0
!
!  4. Error messages
!
!     An error message is produced within the QUAD system.
!     If no errors are detected IERROR=0
!     otherwise IERROR > 0
!
!     ierr    Description of error
!     -------------------------------
!     1       Invalid value of iquad
!     2       Invalid value of iq_grid
!     3       Incompatability between iq_grid and input directions
!     4       Error in deleting *.ERR file
!     5       Error generated by Q_SETCONFIG
!     6       Error generated by Q_CHKCFG
!     7       Error generated by Q_CTRGRID
!
!  5. Called by:
!
!     host program, e.g. SWANQUAD4
!
!  6. Subroutines used:
!
!     Q_SETVERSION
!     Q_SETCONFIG
!     Q_CHKCFG
!     Q_SUMMARY
!     Q_ALLOCATE
!     Q_CTRGRID
!     Q_INIT
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!---------------------------------------------------------------------------------------
!
! Local parameters
!
integer iuerr     ! error indicator
integer idepth    ! index over water depths
integer igrid     ! status of quadruplet grid file
integer ia,ik     ! counters
!
real depmin       ! minimum water depth
real dstep        ! directional step
real dgap         ! directional gap between first and last direction
!
call q_setversion ! set version number
!------------------------------------------------------------------------------
! user defined settings
!------------------------------------------------------------------------------
q_mindepth = 0.1               ! Set minimum water depth
q_maxdepth = 2000              ! Set maximum water depth
q_dstep    = 0.1               ! Set minimum step for coding depth files
iscreen    = 6                 ! Identifier for screen output (UNIX=0, WINDOWS=6)
iufind     = 1                 ! search for unit numbers automatically
!----------------------------------------------------------------------------
! Initialisations
!
ierror     = 0                 ! set error condition
iq_stack   = 0                 ! initialize stack for tracing subroutines
qbase      = 'xnl4v5'          ! Base name for quadruplet files
qf_error   = 'xnl5_errors.txt' ! Text file with error messages
lastquadfile = 'quad?????.bqf' ! Initialize name of last retrieved quad file
!
!  set values of physical quantities
!  and store them in quad data area
!
q_grav    = x_grav        ! gravitational acceleration
qf_tail   = pftail        ! Power of parametric spectral tail
iq_type   = iquad         ! Type of method to compute transfer
iq_interp = 1             ! apply bi-linear interpolation
iq_prt    = 0             ! Print output on, to file /qbase/.prt
iq_test   = 0             ! test level
iq_trace  = 0             ! level of subroutine trace
iq_log    = 0             ! Set logging of q_routines off
iq_grid   = iqgrid        ! Grid type for computation of nonlinear transfer
iq_screen = 0             ! enable output to screen
!
!------------------------------------------------------------------------------
! Check input
!------------------------------------------------------------------------------
if(iq_type<1 .or. iq_type>3) then
  ierror = 1
  goto 9999
end if
!
if(iq_grid<1 .or. iq_grid>3) then
  ierror = 2
  goto 9999
end if
!
!  Retrieve size of spectral grid from input
!
fqmin = sigma(1)/(2.*pi)       ! minimum frequency in Hz
fqmax = sigma(nsigma)/(2.*pi)  ! maximum frequency in Hz
nkq   = nsigma                 ! number of frequencies/wave numbers
naq   = ndir                   ! number of directions
!
!  check if directions are given on full circle or in a symmetric sector
!----------------------------------------------------------------------------
!  1: compute directional step
!  2: compute gap between first and last
!  3: compare gap with step
!
dstep = dird(2)-dird(1)                          ! directional step
dgap = 180.- abs(180.- abs(dird(1)-dird(ndir)))  ! directional gap
!
!----------------------------------------------------------------------
if(iq_grid==1 .or. iq_grid==2) then
!
!  check if gap equal to step in the case of full circle
!
  if(abs(dstep-dgap) < 0.001) then
    ierror = 31
    goto 9999
  end if
!
!  check if sector is symmetric around zero in the case of sector grid
!
  if(abs(dird(1)+dird(ndir)) > 0.01) then
    ierror = 32
    goto 9999
  end if
end if
!
q_dird1 = dird(1)
q_dird2 = dird(ndir)
!
! assign unit numbers for QUAD related files
!
! If IUFIND=0, fixed prespecified unit numbers must be given
!    IUFIND=1, the numbers are searched automatically
!
if(iufind==0) then
  luq_err = 103
  luq_tst = 104
  luq_int = 105
  luq_log = 106
  luq_prt = 107
  luq_cfg = 108
  luq_bqf = 109
  luq_grd = 110
  luq_txt = 111
  luq_loc = 112
  luq_trf = 113
  luq_t13 = 114
  luq_fil = 117
end if
!
! delete old Error file, if it exists
!
call z_fileio(trim(qbase)//'.err','DF',iufind,luq_err,iuerr)
if(iuerr/=0) then
  call q_error('e','FILEIO','Problem in deleting error file *.ERR')
  ierror = 4
  goto 9999
end if
!
! create new files, first create logging file
!
call z_fileio(trim(qbase)//'.log','UF',iufind,luq_log,iuerr) ! logging
call z_fileio(trim(qbase)//'.prt','UF',iufind,luq_prt,iuerr) ! general print file
call z_fileio(trim(qbase)//'.tst','UF',iufind,luq_tst,iuerr) ! test output
!
!
write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.log connected to :',luq_log
write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.prt connected to :',luq_prt
write(luq_log,'(2a,i4)') 'XNL_INIT: ',trim(qbase)//'.tst connected to :',luq_tst
!
!
write(luq_prt,'(a)') '---------------------------------------------------------------'
write(luq_prt,'(a)') trim(q_version)
write(luq_prt,'(a)') 'Solution of Boltzmann integral using Webb/Resio/Tracy method'
write(luq_prt,'(a)') '---------------------------------------------------------------'
write(luq_prt,*)
write(luq_prt,'(a)') 'Initialisation'
write(luq_prt,*)
!
if(iproc >=0) write(luq_prt,'(a,i5)') '(MPI) processor number:',iproc
!---------------------------------------------------------------------------------
!   Reset configuration from file, using a backdoor
!---------------------------------------------------------------------------------
call q_setconfig(iquad)
if (iq_err /=0) then
  ierror = 5
  goto 9999
end if
!---------------------------------------------------------------------------------
!  check settings for inconsistencies
!---------------------------------------------------------------------------------
call q_chkconfig
if (iq_err /=0) then
  ierror = 6
  goto 9999
end if
!---------------------------------------------------------------------------------
! determine minimum size of number of points on locus as stored in database
!---------------------------------------------------------------------------------
klocus = nlocus0
if(iq_gauleg > 0) klocus = min(iq_gauleg,klocus)
if(iq_lump   > 0) klocus = min(iq_lump,klocus)
!---------------------------------------------------------------------------------
!  write summary of program settings
!---------------------------------------------------------------------------------
call q_summary
!----------------------------------------------------------------------------------
!  allocate data arrays
!-----------------------------------------------------------------------------
call q_allocate
!------------------------------------------------------------------------------
!  Generate interaction grid and coefficients for each valid water depth
!  Q_CTRGRID  controls grid generation
!------------------------------------------------------------------------------
do idepth=1,ndepth
  q_depth = depth(idepth)
!
  if(iquad==1 .or. iquad==2) q_depth = q_maxdepth
!
  if(q_depth < q_mindepth) then
    call q_error('w','DEPTH','Invalid depth')
    write(luq_err,'(a,e12.5,f10.2)') 'Incorrect depth & minimum:',q_depth,q_mindepth
  else
    call q_init
    call q_ctrgrid(2,igrid)
    if(iq_err /= 0) then
      ierror = 7
      goto 9999
    end if
  end if
!
  if(iquad==1 .and. ndepth > 0) then
    write(luq_prt,'(a)') 'XNL_INIT: For deep water only one grid suffices'
    exit
  end if
end do
!
!
!  Create or open triplet output data file if iq_triq > 0
!
!
9999 continue
!
!! if (iq_log ==0) call z_fileio(trim(qbase)//'.log','DF',iufind,luq_log,iuerr)
!! if (iq_prt ==0) call z_fileio(trim(qbase)//'.prt','DF',iufind,luq_prt,iuerr)
!
!
return
end subroutine
!-----------------------------------------------------------------------------!
subroutine xnl_main(aspec,sigma,angle,nsig,ndir,depth,iquad,xnl,diag, &
&   iproc, ierror)
!-----------------------------------------------------------------------------!
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant Ph. van Vledder
!   |   +---+
!   |   | +---+  Last update: 27 Sept. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use serv_xnl4v5
!
implicit none
!
!  0. Update history
!
!     25/02/1999 Initial version
!     25/07/1999 Various restructurations
!     12/10/1999 Error handling improved
!     15/10/1999 Existing test file deleted
!     25/10/1999 Parameter iq_filt added
!     29/10/1999 iq_call renamed to i_qmain and save statement added
!     08/12/1999 Call to CHKLAW included
!     16/12/1999 Extra output and check for wave number range
!     27/12/1999 Expansion of k-grid now in new subroutine Q_EXPAND
!     06/01/2000 Deallocate of KX and NSPECX removed
!     08/01/2000 Recontructed, subroutine Q_WRTVV splitted off
!     12/01/2000 Diagonal term added to interface
!     26/06/2000 Name changed to XNL_MAIN
!     01/02/2001 Interface change, spectrum must be given as Action density(sigma,theta)
!     06/11/2001 Water depth forced to 1000 when IQUAD==1,2
!     13/08/2002 Upgrade to release 4.0
!     20/08/2002 Action density array copied to internal A-array
!     09/09/2002 Upgrade to release 5
!     16/09/2002 Parameter IPROC included to take care of MPI processors
!     27/09/2002 Description of input argument SIGMA corrected
!
!  1. Purpose:
!
!     Compute nonlinear transfer for a given action density spectrum
!     on a given sigma and direction grid
!
!  2. Method
!
!     Webb/Resio/Tracy/Van Vledder
!
!  3. Parameter list:
!
! Type    I/O          Name               Description
!---------------------------------------------------------------------------------------------
integer,intent(in)  :: nsig             ! number of frequencies (sigma)
integer,intent(in)  :: ndir             ! number of directions
integer,intent(in)  :: iquad            ! method of computing nonlinear quadruplet interactions
integer, intent(in) :: iproc            ! MPI processor number
!
real,   intent(in)  :: aspec(nsig,ndir) ! Action density spectrum as a function of (sigma,theta)
real,   intent(in)  :: sigma(nsig)      ! radian frequencies
real,   intent(in)  :: angle(ndir)      ! directions in radians (sector or full circle)
real,   intent(in)  :: depth            ! water depth
real,   intent(out) :: xnl(nsig,ndir)   ! nonlinear quadruplet interaction computed with
!                                         a certain exact method (k,theta)
real,   intent(out) :: diag(nsig,ndir)  ! diagonal term for semi-implicit integration
integer, intent(out) :: ierror          ! error indicator
!
!--------------------------------------------------------------------------------
!
!  4. Error messages
!
!  5. Called by:
!
!     host program
!
!  6. Subroutines used
!
!     Q_DSCALE        WAM depth scaling
!     Q_XNL4V4        Xnl using Webb/Resio/Tracy/VanVledder
!     Q_STACK         Stack administration
!     Q_CHKCONS       Check conservation of energy, action and momentum
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!--------------------------------------------------------------------------------
!     local variables
!
integer, save :: i_qmain    ! counter number of calls of XNL_MAIN
integer i_qlast             ! value of iquad in last call
!
integer isig         ! counter for sigma values
integer idir         ! counter of directions
real q_dfac          ! depth scale factor for nonlinear transfer
!
real sum_e           ! sum of energy
real sum_a           ! sum of action
real sum_mx          ! sum of momentum in x-direction
real sum_my          ! sum of momentum in y-direction
!
data i_qmain /0/     ! keep track of number of calls of XNL_MAIN
data i_qlast /0/     ! keep track of last call with IQUAD
!
!--------------------------------------------------------------------------------
!
iq_stack =0          ! initialize stack order every time qmain is called
!
call q_stack('+xnl_main')
!
i_qmain = i_qmain + 1
!
if(iq_prt>=1) then
  write(luq_prt,*)
  write(luq_prt,'(a,i4,f16.3,i4)') 'XNL_MAIN: Input arguments: iquad depth iproc:',&
&  iquad,depth,iproc
end if
!
! initialisations for error handling
!
iq_err  = 0             ! No errors detected at start
q_depth = depth         ! water depth to be used in computation
!
if(iquad==1 .or. iquad==2) q_depth=1000.
!     !
!  check water depth to be used in computation
!
if(q_depth < q_mindepth) then
  xnl = 0.
  call q_error('w','DEPTH','Zero transfer returned')
  goto 9999
end if
!
!  check if iquad has changed since last call, this is no more allowed
!
!!if (iquad /= i_qlast .and. i_qmain/=1) then
!!  call q_error('e','IQUAD','Value of IQUAD differs from initial value')
!!  ierror = 1
!!  goto 9999
!!end if
!-----------------------------------------------------------------------------+
!  main choice between various options                                        |
!-----------------------------------------------------------------------------+
!
if(iquad>=1 .and. iquad <=3) then
!
  a = aspec
  call q_xnl4v4(aspec,sigma,angle,nsig,ndir,depth,xnl,diag,ierror)
!
  if(ierror/=0) then
    call q_error('e','wrtvv','Problem in Q_XNL4V4')
    goto 9999
  end if
!------------------------------------------------------------------------------
! compute scale factor to include WAM depth scaling
!------------------------------------------------------------------------------
!
  if(iq_dscale ==1) then
    call q_dscale(aspec,sigma,angle,nsig,ndir,depth,q_grav,q_dfac)
!
    xnl = xnl*q_dfac
!
    if(iq_prt >=1) write(luq_prt,'(a,f7.4)') 'XNL_MAIN depth scale factor:',q_dfac
  end if
end if
!
!  check conservation laws
!
  call q_chkcons(xnl,nsig,ndir,sum_e,sum_a,sum_mx,sum_my)
!
  if(iq_prt >= 1) then
    write(luq_prt,'(a)')        'XNL_MAIN: Conservation checks'
    write(luq_prt,'(a,4e13.5)') 'XNL_MAIN: E/A/MOMX/MOMY:',sum_e,sum_a,sum_mx,sum_my
  end if
!
9999 continue
!
ierror = iq_err
!
if(iq_log >= 1) then
  write(luq_log,*)
  write(luq_log,'(a,i4)') 'XNL_MAIN: Number of warnings:',iq_warn
  write(luq_log,'(a,i4)') 'XNL_MAIN: Number of errors  :',iq_err
end if
!
!!i_qlast = iquad
!
call q_stack('-xnl_main')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_allocate
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 August 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
implicit none
!
!
!  0. Update history
!
!     05/10/1999  initial version
!     25/11/1999  logging output added
!     10/01/2001  Inconsistency fixed
!     01/10/2001  Array's for k4-locus added
!     08/08/2002  Upgrade to release 4.0
!     20/08/2002  Internal action density spectrum added A
!     11/09/2002  q_kpow added
!     25/04/2003  Name modified from q_alloc -> q_allocate
!                 Triplet array added
!     02/05/2003  Bug fixed in allocate of triplet arrays
!     11/06/2003  Array SYM_LOC added
!                 Parameter KLOCUS introduced for actual maximum size of locus
!     16/06/2003  Loci information included, moved from Q_XNL4V4
!     08/08/2003  Value of MLOCUS modified, now 1.3*NLOCUS0, was 1.2
!
!  1. Purpose:
!
!     Check configuration for non-linear transfer
!
!  2. Method
!
!     Allocate data arrays
!
!  3. Parameters used
!
!  4. Error messaged
!
!  5. Called by:
!
!     XNL_INIT
!
!  6. Subroutines used
!
!     Q_STACK
!
!  7. Remarks
!
!  8. Stucture
!
!  9. Switches
!
!     /S  Subroutine tracing
!
! 10. Source code
!-------------------------------------------------------------------------------
!
!  Local variables
!-------------------------------------------------------------------------------
integer maq   ! number of theta elements in grid matrix
integer mkq   ! number of k-elements in grid matrix
!-------------------------------------------------------------------------------
call q_stack('+q_allocate')
!
if(iq_geom==0) then
  mkq = nkq*(nkq+1)/2
else
  mkq = nkq
end if
!
maq    = naq/2+1
mlocus = 1.3*nlocus0
!
if(iq_log >=1) write(luq_log,'(a,4i4)') &
&  'Q_ALLOCATE:  mkq maq mlocus klocus:',mkq,maq,mlocus,klocus
!
if (allocated (q_xk)) deallocate (q_xk)  ; allocate(q_xk(0:nkq))
if (allocated (q_sk)) deallocate (q_sk)  ; allocate(q_sk(0:nkq))
!
if (allocated (quad_nloc)) deallocate (quad_nloc) ;  allocate (quad_nloc(mkq,maq))
if (allocated (quad_ik2))  deallocate (quad_ik2)  ;  allocate (quad_ik2(mkq,maq,klocus))
if (allocated (quad_ia2))  deallocate (quad_ia2)  ;  allocate (quad_ia2(mkq,maq,klocus))
if (allocated (quad_ik4))  deallocate (quad_ik4)  ;  allocate (quad_ik4(mkq,maq,klocus))
if (allocated (quad_ia4))  deallocate (quad_ia4)  ;  allocate (quad_ia4(mkq,maq,klocus))
if (allocated (quad_w1k2)) deallocate (quad_w1k2) ;  allocate (quad_w1k2(mkq,maq,klocus))
if (allocated (quad_w2k2)) deallocate (quad_w2k2) ;  allocate (quad_w2k2(mkq,maq,klocus))
if (allocated (quad_w3k2)) deallocate (quad_w3k2) ;  allocate (quad_w3k2(mkq,maq,klocus))
if (allocated (quad_w4k2)) deallocate (quad_w4k2) ;  allocate (quad_w4k2(mkq,maq,klocus))
if (allocated (quad_w1k4)) deallocate (quad_w1k4) ;  allocate (quad_w1k4(mkq,maq,klocus))
if (allocated (quad_w2k4)) deallocate (quad_w2k4) ;  allocate (quad_w2k4(mkq,maq,klocus))
if (allocated (quad_w3k4)) deallocate (quad_w3k4) ;  allocate (quad_w3k4(mkq,maq,klocus))
if (allocated (quad_w4k4)) deallocate (quad_w4k4) ;  allocate (quad_w4k4(mkq,maq,klocus))
if (allocated (quad_zz))   deallocate (quad_zz)   ;  allocate (quad_zz(mkq,maq,klocus))
!
if (allocated(x2_loc))   deallocate(x2_loc)     ;  allocate (x2_loc(mlocus))
if (allocated(y2_loc))   deallocate(y2_loc)     ;  allocate (y2_loc(mlocus))
if (allocated(x4_loc))   deallocate(x4_loc)     ;  allocate (x4_loc(mlocus))
if (allocated(y4_loc))   deallocate(y4_loc)     ;  allocate (y4_loc(mlocus))
if (allocated(z_loc))    deallocate(z_loc)      ;  allocate (z_loc(mlocus))
if (allocated(s_loc))    deallocate(s_loc)      ;  allocate (s_loc(mlocus))
if (allocated(ds_loc))   deallocate(ds_loc)     ;  allocate (ds_loc(mlocus))
if (allocated(jac_loc))  deallocate(jac_loc)    ;  allocate (jac_loc(mlocus))
if (allocated(cple_loc)) deallocate(cple_loc)   ;  allocate (cple_loc(mlocus))
if (allocated(a_pol))    deallocate(a_pol)      ;  allocate (a_pol(mlocus))
if (allocated(c_pol))    deallocate(c_pol)      ;  allocate (c_pol(mlocus))
if (allocated(k_pol))    deallocate(k_pol)      ;  allocate (k_pol(mlocus))
if (allocated(sym_loc))  deallocate (sym_loc)   ;  allocate (sym_loc(mlocus))
!
if (allocated(x2_mod))   deallocate (x2_mod)    ;  allocate (x2_mod(mlocus))
if (allocated(y2_mod))   deallocate (y2_mod)    ;  allocate (y2_mod(mlocus))
if (allocated(x4_mod))   deallocate (x4_mod)    ;  allocate (x4_mod(mlocus))
if (allocated(y4_mod))   deallocate (y4_mod)    ;  allocate (y4_mod(mlocus))
if (allocated(z_mod))    deallocate (z_mod)     ;  allocate (z_mod(mlocus))
if (allocated(s_mod))    deallocate (s_mod)     ;  allocate (s_mod(mlocus))
if (allocated(ds_mod))   deallocate (ds_mod)    ;  allocate (ds_mod(mlocus))
if (allocated(jac_mod))  deallocate (jac_mod)   ;  allocate (jac_mod(mlocus))
if (allocated(cple_mod)) deallocate (cple_mod)  ;  allocate (cple_mod(mlocus))
if (allocated(sym_mod))  deallocate (sym_mod)   ;  allocate (sym_mod(mlocus))
!
if (allocated(k2m_mod)) deallocate (k2m_mod)   ;  allocate (k2m_mod(mlocus))
if (allocated(k2a_mod)) deallocate (k2a_mod)   ;  allocate (k2a_mod(mlocus))
if (allocated(k4m_mod)) deallocate (k4m_mod)   ;  allocate (k4m_mod(mlocus))
if (allocated(k4a_mod)) deallocate (k4a_mod)   ;  allocate (k4a_mod(mlocus))
!
if (allocated(wk_k2)) deallocate(wk_k2)   ;     allocate (wk_k2(mlocus))
if (allocated(wa_k2)) deallocate(wa_k2)   ;     allocate (wa_k2(mlocus))
if (allocated(wt_k2)) deallocate(wt_k2)   ;     allocate (wt_k2(mlocus))
if (allocated(wk_k4)) deallocate(wk_k4)   ;     allocate (wk_k4(mlocus))
if (allocated(wa_k4)) deallocate(wa_k4)   ;     allocate (wa_k4(mlocus))
if (allocated(wt_k4)) deallocate(wt_k4)   ;     allocate (wt_k4(mlocus))
!
!if (allocated(t_wk2))  deallocate (t_wk2) ;     allocate (t_wk2(mlocus))
!if (allocated(t_wa2))  deallocate (t_wa2) ;     allocate (t_wa2(mlocus))
!if (allocated(t_wt2))  deallocate (t_wt2) ;     allocate (t_wt2(mlocus))
!if (allocated(t_wk4))  deallocate (t_wk4) ;     allocate (t_wk4(mlocus))
!if (allocated(t_wa4))  deallocate (t_wa4) ;     allocate (t_wa4(mlocus))
!if (allocated(t_wt4))  deallocate (t_wt4) ;     allocate (t_wt4(mlocus))
!if (allocated(t_sym))  deallocate (t_sym) ;     allocate (t_sym(mlocus))
!if (allocated(t_grad)) deallocate (t_grad);     allocate (t_grad(mlocus))
!if (allocated(t_cple)) deallocate (t_cple);     allocate (t_cple(mlocus))
!if (allocated(t_s))    deallocate (t_s)   ;     allocate (t_s(mlocus))
!if (allocated(t_ds))   deallocate (t_ds)  ;     allocate (t_ds(mlocus))
!
!------------------------------------------------------------------------------
!  allocate data arrays for transformation of locus information
!  and integration along locus
!-----------------------------------------------------------------------------
if (allocated(r_ik2))  deallocate (r_ik2)  ;   allocate (r_ik2(klocus))
if (allocated(r_ia2))  deallocate (r_ia2)  ;   allocate (r_ia2(klocus))
if (allocated(r_ik4))  deallocate (r_ik4)  ;   allocate (r_ik4(klocus))
if (allocated(r_ia4))  deallocate (r_ia4)  ;   allocate (r_ia4(klocus))
if (allocated(r_w1k2)) deallocate (r_w1k2) ;   allocate (r_w1k2(klocus))
if (allocated(r_w2k2)) deallocate (r_w2k2) ;   allocate (r_w2k2(klocus))
if (allocated(r_w3k2)) deallocate (r_w3k2) ;   allocate (r_w3k2(klocus))
if (allocated(r_w4k2)) deallocate (r_w4k2) ;   allocate (r_w4k2(klocus))
if (allocated(r_w1k4)) deallocate (r_w1k4) ;   allocate (r_w1k4(klocus))
if (allocated(r_w2k4)) deallocate (r_w2k4) ;   allocate (r_w2k4(klocus))
if (allocated(r_w3k4)) deallocate (r_w3k4) ;   allocate (r_w3k4(klocus))
if (allocated(r_w4k4)) deallocate (r_w4k4) ;   allocate (r_w4k4(klocus))
if (allocated(r_zz))   deallocate (r_zz)   ;   allocate (r_zz(klocus))
!-------------------------------------------------------------------------------
if (allocated(t_ik2))  deallocate (t_ik2)  ;   allocate (t_ik2(klocus))
if (allocated(t_ia2))  deallocate (t_ia2)  ;   allocate (t_ia2(klocus))
if (allocated(t_ik4))  deallocate (t_ik4)  ;   allocate (t_ik4(klocus))
if (allocated(t_ia4))  deallocate (t_ia4)  ;   allocate (t_ia4(klocus))
if (allocated(t_w1k2)) deallocate (t_w1k2) ;   allocate (t_w1k2(klocus))
if (allocated(t_w2k2)) deallocate (t_w2k2) ;   allocate (t_w2k2(klocus))
if (allocated(t_w3k2)) deallocate (t_w3k2) ;   allocate (t_w3k2(klocus))
if (allocated(t_w4k2)) deallocate (t_w4k2) ;   allocate (t_w4k2(klocus))
if (allocated(t_w1k4)) deallocate (t_w1k4) ;   allocate (t_w1k4(klocus))
if (allocated(t_w2k4)) deallocate (t_w2k4) ;   allocate (t_w2k4(klocus))
if (allocated(t_w3k4)) deallocate (t_w3k4) ;   allocate (t_w3k4(klocus))
if (allocated(t_w4k4)) deallocate (t_w4k4) ;   allocate (t_w4k4(klocus))
if (allocated(t_zz))   deallocate (t_zz)   ;   allocate (t_zz(klocus))
!-------------------------------------------------------------------------------
if (allocated(dt13)) deallocate (dt13) ; allocate(dt13(klocus))
!
!-------------------  spectral grid data ------------------------
!
if (allocated(q_k))    deallocate (q_k)   ;     allocate (q_k(nkq))
if (allocated(q_dk))   deallocate (q_dk)  ;     allocate (q_dk(nkq))
if (allocated(q_kpow)) deallocate (q_kpow);     allocate (q_kpow(nkq))
if (allocated(q_f))    deallocate (q_f)   ;     allocate (q_f(nkq))
if (allocated(q_df))   deallocate (q_df)  ;     allocate (q_df(nkq))
if (allocated(q_sig))  deallocate (q_sig) ;     allocate (q_sig(nkq))
if (allocated(q_dsig)) deallocate (q_dsig);     allocate (q_dsig(nkq))
if (allocated(q_cg))   deallocate (q_cg);       allocate (q_cg(nkq))
if (allocated(q_a))    deallocate (q_a)   ;     allocate (q_a(naq))
if (allocated(q_ad))   deallocate (q_ad)  ;     allocate (q_ad(naq))
!
!
if (allocated(nspec)) deallocate (nspec)  ;     allocate (nspec(nkq,naq))
if (allocated(a)) deallocate (a)          ;     allocate (a(nkq,naq))
if (allocated(nk1d))  deallocate (nk1d)   ;     allocate (nk1d(nkq))
!! if (allocated(qnl))   deallocate (qnl)    ;     allocate (qnl(nkq,naq))
!
if(iq_log >=1) then
  write(luq_log,'(a)')    'Q_ALLOCATE: size of arrays'
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: mkq  :',mkq
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: maq  :',maq
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: nkq  :',nkq
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: naq  :',naq
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: mlocus:',mlocus
  write(luq_log,'(a,i4)') 'Q_ALLOCATE: klocus:',klocus
end if
!
call q_stack('-q_allocate')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_chkconfig
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 12 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
!
!  0. Update history
!
!     28/12/1999  Initial version
!     05/10/1999  Test for iq_filt included
!     01/11/1999  Implicit none introduced
!     12/11/1999  Update of tests
!     22/11/1999  Update of tests
!     28/12/1999  Check of IQ_LOCUS added
!     02/01/2000  IQ_START removed
!     05/01/2000  IQ_INT added
!     08/08/2002  Upgrade to release 4
!     22/08/2002  Extra checks included
!     04/06/2003  parameter IQ_INT renamed IQ_INTEG
!                 Switches IQ_GAULEG, IQ_LUMP added
!     11/06/2003  Name changed from Q_CHKCFG to Q_CHKCONFIG
!                 Parameter iq_space removed
!     12/06/2003  Extra check on IMOD, KLOCUS and NLOCUS0
!     16/06/2003  Switch IQ_SYM added
!
!  1. Purpose:
!
!     Check configuration for computation of non-linear transfer
!
!  2. Method
!
!     Check of each parameter setting
!
!  3. Parameters used
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_INIT
!
!  6. Subroutines used:
!
!     Q_ERROR
!
!  7. Remarkds
!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!
! 10. Source code
!-------------------------------------------------------------------------------------------
! do not use m_xnldata
implicit none
!
call q_stack('+Q_CHKCONFIG')
!
if(qf_tail > -1.) call q_error('e','CONFIG','Incorrect power of spectral: qf_tail')
!
if(iq_cple < 1 .or. iq_cple > 3) &
& call q_error('e','CONFIG','Invalid option for coupling coefficient iq_cple')
!
if(iq_compact < 0 .or. iq_compact > 1) &
& call q_error('e','CONFIG','iq_compact /= 0,1')
!
if(iq_filt < 0 .or. iq_filt > 1) &
& call q_error('e','CONFIG','iq_filt /= 0,1')
!
if(iq_gauleg < 0) &
&  call q_error('e','CONFIG','iq_gauleg <0')
!
if(iq_geom < 0 .or. iq_geom > 1) &
&  call q_error('e','CONFIG','iq_geom /= 0,1')
!
if(iq_interp < 1 .or. iq_interp > 2) &
&  call q_error('e','CONFIG','iq_interp /= 1,2')
!
if(iq_disp > 1 .and. iq_geom ==1) then
  call q_error('e','CONFIG','Invalid combination of iq_disp & iq_geom')
  write(luq_err,'(1x,a,2i4)') 'iq_disp iq_geom:',iq_disp,iq_geom
end if
!
if(iq_lump>0 .and. iq_gauleg>0) then
   call q_error('e','CONFIG','Lumping and Gauss-Legendre interpolation not together')
   write(luq_err,'(1x,a,2i4)') 'iq_lump iq_gauleg:',iq_lump,iq_gauleg
end if
!
if(iq_dscale < 0 .or. iq_dscale > 1) &
&  call q_error('e','CONFIG','Incorrect value for IQ_DSCALE, (0,1)')
!
if(iq_disp < 1 .or. iq_disp >2 ) &
&  call q_error('e','CONFIG','Incorrect value for IQ_DISP [DISP],(1,2) ')
!
if(iq_grid <1 .or. iq_grid > 3) &
&  call q_error('e','CONFIG','Incorrect value for IQ_GRID, (1,2,3)')
!
if(iq_integ < 0 .or. iq_make > 3) then
  call q_error('e','CONFIG','Invalid value for iq_integ')
  write(luq_err,'(1x,a,2i4)') 'iq_integ:',iq_integ
end if
!
if(iq_log < 0) &
&  call q_error('e','CONFIG','Incorrect value for IQ_LOG, (>=0) ')
!
if(iq_locus < 0 .or. iq_locus > 3) &
&  call q_error('e','CONFIG','Incorrect specifier for locus method')
!
if(iq_lump<0) then
  call q_error('e','CONFIG','Invalid value for iq_lump')
  write(luq_err,'(1x,a,2i4)') 'iq_lump:',iq_lump
end if
!
if(iq_make < 1 .or. iq_make > 3) then
  call q_error('e','CONFIG','Invalid value for iq_make')
  write(luq_err,'(1x,a,2i4)') 'iq_make:',iq_make
end if
!
if(iq_mod < 0 .or. iq_mod > 1) &
&  call q_error('e','CONFIG','Incorrect value for IQ_MOD [MOD] (0,1)')
!
if(iq_mod==0 .and. klocus<nlocus0) then
  call q_error('e','CONFIG','klocus < nlocus0')
  write(luq_err,'(a)') 'Lumping or Gauss-Integration enabled when IMOD=0'
end if
!
if(iq_prt < 0) &
&  call q_error('e','CONFIG','Incorrect value for IQ_PRT, (>=0) ')
!
!!if(iq_search==1 .and. iq_type/=3) &
!!&  call q_error('w','CONFIG','search option only active when IQUAD=3')
!
if(iq_sym <0 .or. iq_sym > 1) &
  call q_error('e','CONFIG','Incorrect value of IQ_SYM /=[0,1]')
!
if(iq_test < 0) &
&  call q_error('e','CONFIG','Incorrect value for IQ_TEST, (>=0) ')
!
if(iq_trf < 0 .or. iq_trf > 3) &
& call q_error('e','CONFIG','Incorrect value for IQ_TRF ')
!
!  parameter settings ------------------------------------------------
!
if(fqmin < 0) &
&  call q_error('e','CONFIG','Incorrect value for FQMIN')
!
if(fqmax < 0) &
&  call q_error('e','CONFIG','Incorrect value for FQMAX')
!
if(fqmax <= fqmin) &
&  call q_error('e','CONFIG','fmax <= fmin')
!
if(nkq < 1) &
& call q_error('e','CONFIG','Number of wave numbers NKQ < 0')
!
if(naq < 1) &
& call q_error('e','CONFIG','Number of directions NKQ < 0')
!
if(nlocus0 < 6) &
& call q_error('e','CONFIG','Preferred number of points on locus NLOCUS0 < 6')
!
if(q_sector < 40. .or. q_sector > 180.) &
&  call q_error('e','CONFIG','Sector too small (<40) or too large (>180)')
!
call q_stack('-Q_CHKCONFIG')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_chkcons(xnl,nk,ndir,sum_e,sum_a,sum_mx,sum_my)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last Update:  13 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
implicit none
!
!  0. Update history
!
!     29/07/1999 Initial version
!     01/11/1999 Implicit none introduced
!     08/12/1999 Bug fixed in definition os momentum sum
!     13/08/2002 Upgrade to release 4.0
!
!  1. Purpose:
!
!     Check conservation laws of non-linear transfer
!
!  2. Method
!
!     The following conservation laws should be fulfilled:
!
!     Wave Energy        SUME=0
!     Wave Action        SUMA=0
!     Momentum vector    SUMMX,SUMMY=0
!
!
!  3. Parameter list:
!
!Type   I/O          Name     Description
!
integer, intent(in) :: nk         ! number of wave numbers
integer, intent(in) :: ndir       ! number of directions
real, intent(in)  :: xnl(nk,ndir) ! transfer rate
real, intent(out) :: sum_e        ! sum of wave energy
real, intent(out) :: sum_a        ! sum of wave action
real, intent(out) :: sum_mx       ! sum of momentum in x-direction
real, intent(out) :: sum_my       ! sum of momentum in y-direction
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_MAIN
!
!  6. Subroutines used
!
!     Q_STACK
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real aa    ! action density
real ee    ! energy density
real kk    ! wave number
real momx  ! momentum in x-direction
real momy  ! momentum in y-direction
real qq    ! bin size
!
integer ia ! counter over directions
integer ik ! counter over wave numbers
!
!------------------------------------------------------------------------------
!
call q_stack('+q_chklaw')
!
!  initialize summations
!
sum_a  = 0.
sum_e  = 0.
sum_mx = 0.
sum_my = 0.
!
do ik=1,nkq
  qq = q_delta*q_dk(ik)
  kk = q_k(ik)
!
  do ia = 1,naq
    aa = xnl(ik,ia)
    ee = aa*q_sig(ik)
    momx = aa*kk*cos(q_a(ia))
    momy = aa*kk*sin(q_a(ia))
!
    sum_a  = sum_a  + aa*qq
    sum_e  = sum_e  + ee*qq
    sum_mx = sum_mx + momx*qq
    sum_my = sum_my + momy*qq
  end do
end do
!
call q_stack('-q_chklaw')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_chkres(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,dep,sum_kx,sum_ky,sum_w)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
implicit none
!
!  0. Update history
!
!     16/02/2001  Initial version
!     01/11/2001  Implicit none introduced
!     09/08/2002  Upgrade to release 4.0
!
!  1. Purpose:
!
!     Check resonance conditions of 4 interacting wave numbers
!     for a given water depth and dispersion relation
!
!  2. Method
!
!     The sum of wave number vectors and associated radian frequencies
!     are computed:
!
!     k1 + k2 - (k3 + k4)
!     w1 + w2 - (w3 + w4)
!
!     in which w_i = g k_i tanh(k_i d)
!
!  3. Parameter list:
!
! Type    I/O          Name       Description
!----------------------------------------------------------------------
real, intent(in)  ::   k1x     !  x-component of wave number vector k1
real, intent(in)  ::   k1y     !  y-component of wave number vector k1
real, intent(in)  ::   k2x     !  x-component of wave number vector k2
real, intent(in)  ::   k2y     !  y-component of wave number vector k2
real, intent(in)  ::   k3x     !  x-component of wave number vector k3
real, intent(in)  ::   k3y     !  y-component of wave number vector k3
real, intent(in)  ::   k4x     !  x-component of wave number vector k4
real, intent(in)  ::   k4y     !  y-component of wave number vector k4
real, intent(in)  ::   dep     !  depth in m
real, intent(out) ::   sum_kx  !  sum of x-components of quadruplet
real, intent(out) ::   sum_ky  !  sum of y-components of quadruplet
real, intent(out) ::   sum_w   !  sum of radian frequencies
!
!  4. Error messages
!
!  5. Subroutines used
!
!     X_DISPER
!
!  6. Called by:
!
!     Q_MAKEGRID
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real w1,w2,w3,w4            ! radian frequecies of wave numbers
!!real x_disper               ! dispersion relation
sum_kx = (k1x + k2x) - (k3x + k4x)
sum_ky = (k1y + k2y) - (k3y + k4y)
!
!  compute radian frequency on the basis of current dispersion relation
!
w1 = x_disper(sqrt(k1x**2 + k1y**2),dep)
w2 = x_disper(sqrt(k2x**2 + k2y**2),dep)
w3 = x_disper(sqrt(k3x**2 + k3y**2),dep)
w4 = x_disper(sqrt(k4x**2 + k4y**2),dep)
!
sum_w = w1 + w2 - (w3 + w4)
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_cmplocus(ka,kb,km,kw,loclen)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 August 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!-------------------------------------------------------------------------------
!  0. Update history
!
!     Date        Description
!
!     18/11/1999  Initial version
!     08/12/1999  Tracing of locus updated and Q_TRACE included
!     22/12/1999  Option Q_POLAR included
!     05/01/2000  LOCLEN added in interface
!     09/08/2002  Upgrade to release 4.0
!     10/09/2002  g added to interface with X_CPLE
!                 test output modified
!     12/06/2003  Call to Z_POYAREA added to check POLAR2
!     08/08/2003  Check on areas only for loci with k3m/k1m < 100
!                 Otherwise machine accuracy plays a role
!
!  1. Purpose:
!
!     Compute locus function used for the determination of the
!     resonnance condition
!
!  2. Method
!
!     See ALKYON, 1999
!
!  3. Parameter list:
!
!Type   I/O          Name          Description
!----------!----------------------------------------------------------------------------
real, intent(out) :: ka,kb       ! lowest and highest wave number magnitude of k2-locus
real, intent(out) :: km          ! wave number magnitude at mid point
real, intent(out) :: kw          ! half width of locus
real, intent(out) :: loclen      ! estimated length of locus
!
!  4. Error messages
!
!  5. Called by:
!
!    Q_MAKEGRID
!
!  6. Subroutines used
!
!     z_zero1
!
!
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!-------------------------------------------------------------------------------
!     Local variables
!-------------------------------------------------------------------------------
real k1m          ! magnitude of wave number k1
real k3m          ! magnitude of wave number k3
real pcos,psin    ! cosine and sine of normalize angle of P
real klen         ! total length of line locus for case w1=w3
!
real kx_beg       ! x-component at start point
real ky_beg       ! y-component at start point
real kx_end       ! x-component at end point
real ky_end       ! y-component at end point
!
real dsp,dsm      ! distances in plus and minus direction
real sum          ! total length of locus
!
real w1,w3        ! radian frequencies of wave numbers k1 and k3
real eps          ! local accuracy for determination of roots
real area1        ! area of locus as computed
real area2        ! area of locus as based on LOCPOS and ellipse
real ratio        ! maximum ratio between area1 and area2
!
integer ierr      ! local error level
integer iloc,jloc ! counters along locus
integer itest     ! local test level for test output to unit luqtst
integer ip1       ! index +1
integer im1       ! index -1
integer jj        ! counter
!------------------------------------------------------------------------------
!  function declarations
!
!!real x_disper                 ! dispersion relation
!!real x_cple                   ! coupling coefficient
!!real x_jacobian               ! Jacobian term
!------------------------------------------------------------------------------
call q_stack('+q_cmplocus')
!
!  set initial values
!
eps     = 10.*epsilon(1.)     ! set accuracy 10 times machine accuracy
itest   = iq_test             ! assign test level from overall setting
!! itest   = 1                ! (re)set local test level
!
! compute characteristics of configuration
!
px   = k1x - k3x
py   = k1y - k3y
pmag = sqrt(px**2 + py**2)
xang = atan2(-px,py)
pang = atan2(py,px)
k1m  = sqrt(k1x**2 + k1y**2)
k3m  = sqrt(k3x**2 + k3y**2)
w1   = x_disper(k1m,q_depth)
w3   = x_disper(k3m,q_depth)
q    = w1-w3
!
!  compute cosine and sine of direction of P-vector
!  reverse direction for the case q<0
!
if(q < 0) then
  sang = pang+pi
  pcos = cos(pang+pi)
  psin = sin(pang+pi)
else
  sang = pang
  pcos = cos(pang)
  psin = sin(pang)
end if
!
!
!  first solution along locus: k2 = k3
!
!  check for special case if q = 0
!
if (abs(q) < eps_q) then
!
  call q_loc_w1w3(k1x,k1y,k3x,k3y,nlocus0,x2_loc,y2_loc,x4_loc,y4_loc,s_loc)
  nlocus1 = nlocus0
  ds_loc  = s_loc(2)-s_loc(1)
  klen    = s_loc(nlocus0)
  ka      = 0.
  kb      = 0.
  km      = 0.
  kw      = 0.
  sang    = xang
!
else
!------------------------------------------------------------------------------
!  compute characteristics of locus, such as its position in
!  wave number space
!------------------------------------------------------------------------------
!
  call q_locpos(ka,kb,km,kw,loclen)
  if(iq_err/=0) goto 9999
!
!  compute position of start and end point for tracing
!  the locus
!
  kx_beg = ka*pcos
  ky_beg = ka*psin
  kx_end = kb*pcos
  ky_end = kb*psin
!
!  compute position of locus using polar method
!  see Van Vledder (2000)
!
!%  call q_polar(ka,kb,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
  call q_polar2(ka,kb,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
!
! check area of locus by a simple test  (added 12 June 2003)
!
  call z_polyarea(x2_loc,y2_loc,nlocus1,area1)
  area2 = pi*(kb-ka)*0.5*kw
  ratio = max(area1/area2,area2/area1)
!
!
  if(ratio>1.5 .and. k3m/k1m < 100.) then
    call q_error('e','LOCUS','Severe problem in POLAR2')
    write(luq_err,'(a)') 'Q_CMPLOCUS: ratio > 1.5'
!
    goto 9999
  end if
!
!  01/10/2001
!  compute position of k4 locus by a simple translation
!
  do iloc=1,nlocus1
    x4_loc(iloc) = x2_loc(iloc) + px
    y4_loc(iloc) = y2_loc(iloc) + py
  end do
!
end if
!
if (iq_test >=2) write(luq_tst,'(1x,a,4f12.5,i4)')&
&   'Q_CMPLOCUS: k1x/y k3x/y nlocus:',k1x,k1y,k3x,k3y,nlocus1
!----------------------------------------------------------------------------------
! compute characteristics around locus
!----------------------------------------------------------------------------------
!
s_loc(1) = 0.
sum      = 0
!
do iloc=1,nlocus1
!
!  compute step sizes
!
  if (abs(q) < eps_q) then
!
!  for this case the sum of ds_loc is unequal to s_loc(nlocus1)
!
    sum = s_loc(nlocus1)
  else
!
!   compute indices of previous and next index on locus
!
    ip1 = iloc+1
    if (ip1 > nlocus1) ip1 = 1
    im1 = iloc-1
    if (im1 < 1) im1 = nlocus1
!
    dsp = sqrt((x2_loc(iloc)-x2_loc(ip1))**2 + (y2_loc(iloc)-y2_loc(ip1))**2)
    dsm = sqrt((x2_loc(iloc)-x2_loc(im1))**2 + (y2_loc(iloc)-y2_loc(im1))**2)
    if(iloc < nlocus1) s_loc(iloc + 1) = s_loc(iloc) + dsp
    ds_loc(iloc) = 0.5*(dsp + dsm)
    sum = sum+ds_loc(iloc)
  end if
!
!  compute gradient/Jacobian terms along locus
!
   jac_loc(iloc) = x_jacobian(x2_loc(iloc),y2_loc(iloc),x4_loc(iloc),y4_loc(iloc))
!
!  compute coupling coefficients along locus
!
  k2x = x2_loc(iloc)
  k2y = y2_loc(iloc)
  k4x = x4_loc(iloc)
  k4y = y4_loc(iloc)
!
  cple_loc(iloc) = x_cple(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,iq_cple,q_depth,q_grav)
!
end do
!
!
9999 continue
!
call q_stack('-q_cmplocus')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_ctrgrid(itask,igrid)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 13 Sept. 2003
!   +---+ |   |  Release: 5.03
!         +---+
!
! do not use m_xnldata
use m_fileio
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     Version  Date    Modification
!
!     29/07/1999  Initial version
!     11/10/1999  Error messages via module
!     12/10/1999  File I/O of interaction grid files added, and consistency check
!     25/10/1999  Contents of q_header extended with iq_grid & iq_start & nlocus0
!     27/10/1999  Close statement after reading BQF file
!     01/11/1999  Close statments after call of Q_GRID
!     03/11/1999  Parameter IQ_MAKE included
!     22/11/1999  Use of z_fileio modified, use parameter IUERR if an attempt
!                 to open a non-existing file was made
!     30/11/1999  Extra messages to logging file when file are closed
!     03/01/2000  IQ_START replaced by IQ_LOCUS
!     12/01/2001  Output parameter IGRID added
!                 igrid=0: a proper grid exists or has been made, and will be read
!                      =1: grid file does not exist
!                      =2: grid file exists, but it is incorrect
!      8/08/2002  Upgrade to release 4.0
!     15/08/2002  Bug fixed ininitialisation of igrid
!     19/08/2002  header extended with parameter IQ_INTERP
!     20/08/2002  wave number array replaced by sigma array in grid file
!     22/08/2002  Extra i/o check when reading BQF file
!     23/08/2002  retrieve number of point on locus from BQF file
!     09/09/2002  aqfile and bqfile 5 units for depth
!     10/09/2002  new algorithm for coding depth
!                 Test added to avoid rereading of last read/generated BQF file
!     08/10/2002  Output to test file made conditional
!     05/09/2003  Water depth for creating and testing BQF file DSTEP dependend
!     09/09/2003  Bug fixed in assigning IGRID=0 when BQF still in memory
!     13/09/2003  When BFQ incorrupt, it is deleted and a new one is created
!                 Bug fixed in setting of s_depth  when iq_disp==1
!
!  1. Purpose:
!
!     Control of interaction grid administration
!
!  2. Method
!
!  3. Parameters used
!
integer, intent(in)  :: itask  !  task to perform by Q_CTRGRID
!                                 ==1: read and check header block
!                                 ==2: read and write grid file, according to
!                                      setting of IQ_MAKE
integer, intent(out) :: igrid  !  status of grid checking
!                                 ==0: a proper grid exists
!                                 ==1: grid file does not exist
!                                 ==2: grid file exists, but it is incorrect
!                                 ==3: read error in accessing grid information
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_INIT
!     Q_SEARCHGRID
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_ERROR
!     Q_MAKEGRID
!
!  7. Remarks
!
!     The generation of the database file depend on the control varaible of IQ_MAKE
!     if IQ_MAKE==1, make a grid when needed
!                 2, always make grid
!                 3, make a grid and stop, useful for test purposes
!
!     The maximum number of points on the locus, as stored in the BQF file
!     is read from the header and stored in the variable NLOCUS
!
!  8. Structure
!
!     Make header of grid file
!     Construct name of grid file
!     Check existence of grid file
!     if grid file exists
!       read header
!       check header
!       read first data block
!       check first data block
!       - directions
!       - wave numbers
!       - depth
!       set status of grid file
!     else
!       set status of grid file
!     end if
!
!     set status of generating/reading grid file
!
!     if make new grid
!       compute grid parameters
!       write grid parameters to file
!     else
!       read grid parameters from file
!     end if
!
!
!  9. Switches
!
! 10. Source code
!-------------------------------------------------------------------------------
!     Local variables
!
integer iaz,ikz,jaz,jkz                ! counters for checking header of BQF file
integer iz_geom,iz_disp,iz_cple        ! values as read in BQF file
integer naz                            ! number of directions in BQF file
integer nkz                            ! number of wave numbers in BQF file
integer idep,jdep                      ! coding for depth in BQF file
!
logical lbqf                           ! flag for existence of BQF file
real s_depth                           ! (stepped) depth
real q_depth_saved                     ! save input water depth, needed after call to Q_MAKEGRID
real z_depth                           ! water depth in BQF file
real, allocatable :: z_ad(:),z_sig(:)  ! directions and radian frequencies of grid in BQF file
integer ierr,iuerr                     ! error variables
!------------------------------------------------------------------------------
!
call q_stack('+q_ctrgrid')
!
!  echo input arguments
!
!
q_depth_saved = q_depth
!
!  generate header of BQF file
!
q_header = '000000000000000000000'
!           123456789012345678901
!                    1         2
write(q_header,'(3i3.3,6i2.2)') naq,nkq,nlocus0,&
&  iq_grid,iq_geom,iq_disp,iq_cple,iq_locus,iq_interp
!
if(iq_prt>=2) then
  write(luq_prt,'(2a)')    'Q_CTRGRID: header info:',trim(q_header)
  write(luq_prt,'(a,3i5)') 'Q_CTRGRID: naq nkq nlocus0:',naq,nkq,nlocus0
  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_grid,iq_geom,iq_disp:',iq_grid,iq_geom,iq_disp
  write(luq_prt,'(a,3i3)') 'Q_CTRGRID: iq_cple,iq_locus,iq_interp:',iq_cple,iq_locus,iq_interp
end if
!
!------------------------------------------------------------------------------
! construct name of grid file
!
if(iq_disp==1) then
  bqname = 'quad99999.bqf'
  s_depth = q_maxdepth
!
elseif(iq_disp==2) then
!
!---------------------------------------------------------------------------------------------
! generate code for actual depth
!---------------------------------------------------------------------------------------------
  idep = int(q_depth/q_dstep+0.5)
  jdep = idep*int(10.*q_dstep)
  jdep = max(1,jdep)
  jdep = min(99999,jdep)
!
  s_depth = real(idep)*q_dstep
!
!
  bqname = 'quad00000.bqf'
  write(bqname(5:9),'(i5.5)') min(int(q_maxdepth*10),jdep)
!
else
  call q_error('e','DISPER','Incorrect value for IQ_DISP')
  write(luq_err,'(a,i4)') 'IQ_DISP=',iq_disp
  goto 9999
end if
!
!
!-------------------------------------------------------------------------------------------
!  Compare LASTQUADFILE with bqname
!  if equal skip reading of BQF file, including checks of header
!-------------------------------------------------------------------------------------------
!
if(lastquadfile==bqname) then
  if(iq_screen>0) write(iscreen,'(2a)')   'Q_CTRGRID: Rereading of bqfile skipped: ',lastquadfile
  igrid = 0
  goto 9999
end if
!-------------------------------------------------------------------------------------------
if(iq_prt >= 2) then
  write(luq_prt,'(2a)') 'Q_CTRGRID: Header line of grid file:',trim(q_header)
  write(luq_prt,'(2a)') 'Q_CTRGRID: Name of BINARY grid file:',trim(bqname)
end if
!------------------------------------------------------------------------------
!
!  check if binary data file exists
!
call z_fileio(bqname,'OU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
!
if(iq_prt >= 2) write(luq_prt,'(2a,2x,2i4)') &
& 'Q_CTRGRID:  bqname:',trim(bqname),luq_bqf,iuerr
!
if(itask==2 .and. iq_make==2) luq_bqf=-1
!
!  if the file exists,
!     read header information
!     check header of file
!  end
!
!  If header is incorrect, set flag IQ_GRID to TRUE for generating new grid
!
if(luq_bqf > 0 .and. iuerr ==0) then
  if(iq_prt >= 2) then
    write(luq_prt,'(2a)')   'Q_CTRGRID: Binary grid file detected: ',trim(bqname)
    write(luq_prt,'(a,i4)') 'Q_CTRGRID: Connected to unit:',luq_bqf
  end if
!
!
! grid file exists, unless proven otherwise
!---------------------------------------------------------------------------------
!
  lq_grid = .false.
  igrid   = 0
  read(luq_bqf,iostat=ierr) r_header
  if(ierr/=0) then
    call q_error('w','READBQF','Read error for header in BQF file')
    write(luq_err,'(a)') 'BQF file deleted'
    call z_fileio(bqname,'DU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
    igrid = 3
    lq_grid = .true.
  else
    read(r_header,'(6x,i3)') nlocus
!
  end if
!-----------------------------------------------------------------------------
!   check header of grid file
!-----------------------------------------------------------------------------
!
  if(trim(r_header)/=trim(q_header).and. .not.lq_grid) then
    lq_grid = .true.
    igrid   = 2
    if(iq_prt >=2) then
      write(luq_prt,'(a,1x,a)') &
&     'Q_CTRGRID: Header in binary quad file         :',trim(r_header)
      write(luq_prt,'(a,1x,a)') &
&     'Q_CTRGRID: Expected header of binary quad file:',trim(q_header)
      write(luq_prt,'(a)') 'Q_CTRGRID: The file headers disagree'
      write(luq_prt,'(a)') 'Q_CTRGRID: A new grid will be generated'
    end if
  end if
!------------------------------------------------------------------------------
!  check other parts of binary grid file
!
  if(.not.lq_grid) then
    read(luq_bqf) naz,nkz
    allocate (z_sig(nkz),z_ad(naz))
    read(luq_bqf) z_sig
    read(luq_bqf) z_ad
    read(luq_bqf) iz_geom,iz_disp,iz_cple
    read(luq_bqf) z_depth
!
    if(iq_prt >=2) then
      write(luq_prt,'(a)')    'Q_CTRGRID: Contents of BQF file'
      write(luq_prt,'(2a)')   'Q_CTRGRID: Header:',trim(r_header)
      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NK:',nkz
      write(luq_prt,'(a,i4)') 'Q_CTRGRID: NA:',naz
    end if
  end if
!---------------------------------------------------------------------------------------
! check spectral interaction grid and depth for consistency
!---------------------------------------------------------------------------------------
  if(.not. lq_grid) then
az: do iaz = 1,naz
      if(abs(q_ad(iaz)-z_ad(iaz)) > 0.01) then
        write(luq_prt,'(a)') 'Q_CTRGRID: Directions do not agree'
        do jaz=1,naz
          write(luq_prt,'(1x,a,i4,2f10.3)') 'iaz q_ad z_ad:',jaz,q_ad(jaz),z_ad(jaz)
        end do
        lq_grid = .true.
        igrid = 2
        exit az
      end if
    end do az
  end if
!
  if(.not. lq_grid) then
ak: do ikz = 1,nkz
      if(abs(q_sig(ikz)-z_sig(ikz)) > 0.01) then
        write(luq_prt,'(a)') 'Q_CTRGRID: Wave numbers do not agree'
        do jkz=1,nkz
          write(luq_prt,'(1x,a,i4,2f10.3)') 'ikz q_k z_sig:',jkz,q_sig(jkz),z_sig(jkz)
        end do
        lq_grid = .true.
        igrid = 2
        exit ak
      end if
    end do ak
  end if
!
!  compare water depths
!
  if(abs(z_depth-s_depth) > 0.09 .and. iq_disp > 1 .and. .not. lq_grid) then
    write(luq_prt,'(a)') 'Q_CTRGRID: Water depths do not agree'
    write(luq_prt,'(a,2f16.2)') 'Q_CTRGRID: q_depth z_depth:',q_depth,z_depth
    lq_grid = .true.
    igrid = 2
  end if
!
  if(lq_grid) then
    close(luq_bqf)
    if(iq_log >= 1) write(luq_log,'(a)') 'Q_CTRGRID: Existing BQF-file invalid, it will be closed'
  end if
!
else
  lq_grid = .true.
  igrid = 1
end if
!------------------------------------------------------------------------------
if(itask==1) then
  if(luq_bqf>0) call z_fclose(luq_bqf)
  goto 9999
end if
!-----------------------------------------------------------------------------
!  if lq_grid==true  a new grid has to be generated
!                    if not, read the grid information into memory
!  or iq_make==2     always make an interaction grid
!  or iq_make==3     as 2, plus stop after making grid
!
if(lq_grid .or. iq_make==2 .or. iq_make==3) then
!
  if(luq_bqf>0) call z_fclose(luq_bqf)
  call z_fileio(bqname,'UU',iufind,luq_bqf,iuerr)   ! binary quadruplet file
!
  if(iq_log >= 1) then
    write(luq_log,*)
    write(luq_log,'(a)') 'Q_CTRGRID: New grid will be generated'
    write(luq_log,'(a,a)') 'Q_CTRGRID: Name of BQF file:',trim(bqname)
    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' connected to :',luq_bqf
  end if
!
  if(iq_screen >= 1) write(iscreen,'(2a)') &
& 'Q_CTRGRID: Generating wave number grid for quadruplet interactions: ',trim(bqname)
!
  q_depth = s_depth
  call q_makegrid
  q_depth = q_depth_saved
!
  if(iq_err /=0) then
     lastquadfile = 'quad_err_.bqf'
     goto 9999
  end if
!
  igrid = 0
!
  close(luq_bqf)
!
  if(iq_log >=1) then
    write(luq_log,'(a,i4)') 'Q_CTRGRID: '//trim(bqname)//' disconnected from:',luq_bqf
  end if
!
  if(iq_screen >=1) write(iscreen,'(a)') 'Q_CTRGRID: Grid generation completed succesfully'
!----------------------------------------------------------------------------------------
!
!  check of header and spectral grid succesfull
!  such that data can be read from BQF file
!----------------------------------------------------------------------------------------
!
else
  if(iq_screen >= 1) write(iscreen,'(2a)') 'Q_CTRGRID: Reading existing grid: ',trim(bqname)
  if(iq_prt >= 1)    write(luq_prt,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
  if(iq_log >= 1)    write(luq_log,'(2a)')  'Q_CTRGRID: Existing grid will be read:',trim(bqname)
!
  read(luq_bqf) quad_nloc
  read(luq_bqf) quad_ik2
  read(luq_bqf) quad_ia2
  read(luq_bqf) quad_ik4
  read(luq_bqf) quad_ia4
  read(luq_bqf) quad_w1k2
  read(luq_bqf) quad_w2k2
  read(luq_bqf) quad_w3k2
  read(luq_bqf) quad_w4k2
  read(luq_bqf) quad_w1k4
  read(luq_bqf) quad_w2k4
  read(luq_bqf) quad_w3k4
  read(luq_bqf) quad_w4k4
  read(luq_bqf) quad_zz
!
  lastquadfile = bqname
!
  close(luq_bqf)
!

end if
!
9999 continue
!
if (allocated(z_ad)) deallocate(z_ad,z_sig)
!
!
call q_stack('-q_ctrgrid')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_dscale(n,sigma,angle,nsig,nang,depth,grav,q_dfac)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 23 Aug. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
use serv_xnl4v5
implicit none
!
!  0. Update history
!
!      Date        Modification
!
!      25/02/1999  Initial version
!       2/12/1999  Result modified if total energy <= 0
!                  Cosmetic changes
!      13/08/2002  Upgrade to release 4.0
!      23/09/2002  Mean wave number multiplied by 0.75
!
!  1. Purpose:
!
!     Compute scaling factor for nonlinear transfer in finite depth
!
!  2. Method
!
!     Compute mean wave number km
!
!     Compute scale factor based on parameterized function of (km*d)
!     according to Herterich and Hasselmann
!     and parameterisation from WAM model
!
!
!  3. Interface parameter list:
!
! Type          I/O     Name           Description
!-------------------------------------------------------------------------
integer, intent (in) :: nsig          ! Number of sigma-values
integer, intent (in) :: nang          ! Number of directions
real, intent(in)     :: n(nsig,nang)  ! N(nsig,nang) Action density
real, intent(in)     :: sigma(nsig)   ! sigma values
real, intent(in)     :: angle(nang)   ! directions in (radians)
real, intent(in)     :: depth         ! Depth (m)
real, intent(in)     :: grav          ! Gravitational acceleration
real, intent(out)    :: q_dfac        ! scale factor
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_MAIN
!
!  6. Subroutines used
!
!     x_wnumb
!     z_steps
!     q_stack
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     local variables
!
real w          ! radian frequency
real kk         ! local wave number
real sqkk       ! square root of local wave number
real dnn        ! summation quantity
real kms        ! mean wave number
real kd         ! depth*mean wave number product
real sum0       ! summation variable for total energy
real sumk       ! summation variable for wave number
real delta      ! directional step, in radians
!
integer isig    ! counter over sigma loop
integer iang    ! counter over direction loop
!
! functions
!!!real z_wnumb    ! function to compute wave number
!
! temporary data
!
real dsigma(nsig) ! step size of sigma array, used for integration
!------------------------------------------------------------------------------
!
call q_stack('+q_dscale')
!
call z_steps(sigma,dsigma,nsig)   !  compute step size of sigma's
delta = angle(2)-angle(1)         !  compute directional step (radians)
!
sum0 = 0.
sumk = 0.
!
!  compute sums for total energy andwave number
!
do isig = 1,nsig
  w    = sigma(isig)
  kk   = z_wnumb(w,depth,grav)   ! compute wave number for given sigma,depth
  sqkk = sqrt(kk)
  do iang=1,nang
    dnn  = n(isig,iang)*dsigma(isig)*delta
    sum0 = sum0 + dnn
    sumk = sumk + 1./sqkk*dnn
  end do
end do
!
!  compute mean wave number and scale factor based
!  on the WAM approximation
!
if(sum0 > 0) then
  kms = (sum0/sumk)**2
  kd = max(0.5,0.75*kms*depth)
  q_dfac = 1+5.5/kd*(1.-5./6.*kd)*exp(-5./4.*kd)
!  pause
else
  kms = 0.
  kd  = 0.
  q_dfac = 1.
end if
!
call q_stack('-q_dscale')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_error(err_type,err_name,err_msg)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update 8 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
use m_fileio
implicit none
!
!  0. Update history
!
!     0.01  22/06/1999   Initial version
!     0.02  20/07/1999   Error message included
!     0.03  24/09/1999   Full message read from file Q_ERROR.TXT
!                        input argument ERR_NAME added
!     0.04  13/10/1999   Reading of multiple lines in Q_ERROR.TXT file
!     0.05  26/11/1999   Layout modified
!     0.06  30/11/1999   Extra output added, also to screen
!     4.01  08/08/2002   Upgrade to release 4
!
!  1. Purpose:
!
!     Error handling routine, produces a warning to an error
!     that has occured prints the error message and print
!     module stack to trace the origin of the error/
!
!  3. Parameter list:
!
!Type         I/O            Name           Description
!------------------------------------------------------------------------------
character(len=1), intent(in)   :: err_type   ! type of error
!                                              w or W: Warning or non-terminating error
!                                              e or E: terminating error
character(len=*), intent(in) :: err_name     !  reference to error message
character(len=*), intent(in) :: err_msg      !  Optional additional error message
!
!  4. Error messages
!
!  5. Called by:
!
!     All q_** subroutines
!
!  6. Subroutines used
!
!  7. Remarks
!
!     The reference to an error message is stored in the
!     string ERR_NAME. For each error number an associated
!     error is given.
!
!     No call is made to subroutine q_trace to avoid
!     infinite recursion
!
character(len=80) qline ! Input line from file with error messges
integer ntext           ! number of text line
integer iend            ! indicator for end of line
integer iutxt           ! unit number for text file
integer iuerr           ! indicator for error
integer j_stack         ! counter in printing stack
integer ispace          ! indicates that first character of line is space
!
call z_fileio(trim(qbase)//'.err','UF',iufind,luq_err,iuerr)   ! error messages
!
! logging of unit number
!
if(iq_log >= 1) write(luq_log,'(a,i4)') &
& 'Q_ERROR: '//trim(qbase)//'.ERR connected to unit:',luq_err
!
!  write general information, when the first error or
!
if(iq_warn ==0 .and. iq_err==0)  then
  write(luq_err,'(a)') q_version
  write(luq_err,'(a)')'--------------------------------------------------'
end if
!
!  check type of error
!
if(index('wW',err_type) > 0) then
  iq_warn = iq_warn + 1
  write(luq_err,'(a,i4)') 'Warning or non-terminating error:',iq_warn
  write(luq_err,'(a,a)')  'Name of error:',trim(err_name)
!
elseif(index('eE',err_type) > 0) then
  iq_err = iq_err + 1
  write(luq_err,'(a,i4)') 'Terminating error:',iq_err
  write(luq_err,'(a,a)')  'Name of error:',trim(err_name)
  write(*,'(1x,a,i4)') 'Terminating error:',iq_err
  write(*,'(1x,a,a)')  'Name of error:',trim(err_name)
end if
!
!  search explanation of error message in the file
!  QF_ERROR, set in XNL_INIT
!
ntext = len_trim(err_name)
!
if(ntext > 0) then
  call z_fileio(qf_error,'OF',iufind,luq_txt,iutxt)
!
  if(iutxt < 0) then
    if(iq_log > 0) write(luq_log,'(3a)') &
&   'Q_ERROR: File ',trim(qf_error),' does not exist in current directory'
!
  else
    if(iq_log >= 1) write(luq_log,'(a,i4)') &
&   'Q_ERROR: File Q_ERROR.TXT connected to unit:',luq_txt
    iend=0
!
!   scan all lines in the text file with error messages
!
    do while (iend==0)
      read(luq_txt,'(a)',iostat=iend) qline
      if(iend==0) then
!
!  check code word exists in text file
!
        if(qline(1:ntext) == err_name(1:ntext)) then
          write(luq_err,*)
          write(luq_err,'(a)') 'Explanation of error, and recommended action'
          write(luq_err,'(a)') '--------------------------------------------'
          write(luq_err,'(a)') trim(qline)
!
!  read following lines until end of file or a non-space in column 1
!
          ispace = 1
          do while (ispace ==1)
            read(luq_txt,'(a)',iostat=iend) qline
!
!  check conditions
!
            if(iend==0) then
              if(qline(1:1) == ' ') then
                write(luq_err,'(a)') trim(qline)
              else
                ispace = 0
              end if
            else
              ispace = 0
            end if
          end do
        end if
      end if
    end do
!
!  close text file with error messages
!
    close(luq_txt)
    if(iq_log >= 1) write(luq_log,'(3a,i4)') &
&   'Q_ERROR: File ',trim(qf_error),'  disconnected from unit:',luq_txt
  end if
end if
!
if(len_trim(err_msg) > 0) then
  write(luq_err,*)
  write(luq_err,'(a)') 'Additional message from point of occurrence:'
  write(luq_err,'(a)') trim(err_msg)
  write(luq_err,*)
end if
!
!  print stack of subroutines to trace the location where the
!  error occurred
!
write(luq_err,'(a)') 'Trace of error'
write(luq_err,'(a)') '--------------'
do j_stack=1,iq_stack
  write(luq_err,'(1x,i4,2x,a)') j_stack,trim(cstack(j_stack))
end do
!
write(luq_err,*)
!
if(iq_warn > 10) stop 'Too many warnings'
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_getlocus(ik1,ia1,ik3,ia3,ifnd)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 27 August 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
!----------------------------------------------------------------------------------
implicit none
!
!  0. Update history
!
!    25/02/1999  Initial version
!    15/04/1999  Extra parameter IFND to indicate if a
!                reference locus exists in the data base
!    19/07/1999  Restructured and some bugs removed
!    20/07/1999  Option added to compute locus directly, no database
!                or scaling involved. IFND < 0
!    15/10/1999  Information to transformation file updated
!    16/10/1999  Equation for computing address for storage updated
!    25/10/1999  Transformations updated
!    28/10/1999  Local variables ia1 and ia3 may not be changed, temp. variables
!                it1 and it3 included
!    29/10/1999  Use of IQ_TRF modified
!                EPSK introduced to check equality of loci
!    28/12/1999  A_CMPLOC renamed to Q_CMPLOC
!    03/01/2000  IQ_START replaced by IQ_LOCUS
!    05/01/2000  Interface with Q_CMPLOC modified
!    08/08/2002  Upgrade to release 4.0
!    13/08/2002  Indexing in terms of integers and real weights
!                upgrade to release 4.0
!    19/08/2002  Bug fixed in transforming CASE 6
!                Interpolation option added
!    12/06/2003  Parameter t_ws set equal to r_ws
!    13/06/2003  Parameter t_cple, t_jac and t_sym assigned
!                Bug fixed in nearest bin approach, symmetry regained
!    27/08/2003  Short-cut when number of points on locus is ZERO
!
!  1. Purpose:
!
!     Retrieve locus from basic locus as stored in the database
!
!  2. Method
!
!     In the case of geometric scaling, k-scaling is used using scale laws
!     described by Tracy
!
!     Directional transformation using linear transformations, shifting and mirror
!     imaging.
!
!
!  3. Parameter list:
!
!Type      I/O           name       Description
!------------------------------------------------------------------------------
integer, intent(in)  ::  ik1    !  k-index of wave number k1
integer, intent(in)  ::  ia1    !  theta-index of wave number k1
integer, intent(in)  ::  ik3    !  k-index of wave number k3
integer, intent(in)  ::  ia3    !  theta-index of wave number k3
integer, intent(out) :: ifnd    !  indicator if reference locus exists in database
!
!  4. Error messages
!
!  5. Called by
!
!     Q_T13V4
!
!  6. Subroutines used
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer it1,it3           ! work indices for directions, copy of ia1 and ia3
integer idir              ! switch to indicate if locus should be inverted
integer itrans            ! type of transformation
integer iloc,jloc         ! counters along locus
integer iadif,ikdif       ! difference in angular and k-index
integer ja1r,ja3r
integer imirror           ! extra step when locus is mirrorred
!
integer ibeta,kdif
integer nloc             ! number of points on locus
!
integer ierr
integer kmem                 ! index for storing 2-d matrix in 1-d array
integer amem                 ! index for storing direction of reference wave number k3
!
real lambda                  ! geometric scale factor
real j_lambda                ! scale factor for Jacobian term
real c_lambda                ! scale factor for coupling coefficient
real zz_lambda               ! combined scale factor
!
real xt2(nlocus),yt2(nlocus) ! xy-components of test k2-locus
real xt4(nlocus),yt4(nlocus) ! xy-components of test k4-locus
real wk,wa,vk,va
!! \A
!! real x_kfunc                 ! real function to compute wave number
!! \Z
!
integer ikmin,ja1,ja3,jk1,jk3,itmin
integer ibdif,nhalf
integer iaq,ikq              ! counters for loop over direction and wave numbers
!------------------------------------------------------------------------------
!
!! data i_getloc /0/  ! Initialise counter
!-------------------------------------------------------------------------------------
call q_stack('+q_getlocus')
!
!------------------------------------------------------------------------------
! initialisations
!------------------------------------------------------------------------------
!
it1 = ia1
it3 = ia3
!
imirror = 1
!
ikmin = min(ik1,ik3)    ! compute minimum of wave number index
ikdif = abs(ik1-ik3)    ! compute difference between wave number indices
!
if (iq_geom ==0) then
  jk1 = min(ik1,ik3)
  jk3 = max(ik1,ik3)
else
  jk1 = 1
  jk3 = ikdif + 1       ! compute k-index of wave number k3 relative to reference wave number
end if
!
itmin = min(it1,it3)    ! compute minimum angle of k1 and k3
iadif = abs(it1-it3)    ! difference index
ja1   = 1               ! index of direction of reference wave number k1
ja3   = iadif+iaref     ! compute theta-index of direction of wave number k3
!
!------------------------------------------------------------------------------
!  circle grid, modify ranges and transformation variables
!------------------------------------------------------------------------------
!
if (iq_grid==3) then
  nhalf = naq/2
  if (iadif > nhalf) then
    if(it1 > nhalf) it1 = it1 - naq
    if(it3 > nhalf) it3 = it3 - naq
  end if
  itmin = min(it1,it3)
  ibdif = (naq - abs(naq-2*abs(it1-it3)))/2   ! compute shortest difference in indices
                                              ! while taking care of periodicity
  ja3   = ibdif + 1
  iadif = ibdif
end if
!
ja1r = 1
ja3r = iadif + 1
amem = iadif + 1   ! compute index of reference wave number k3 in interaction grid
!
!------------------------------------------------------------------------------
! obtain k-index of reference wave number
!------------------------------------------------------------------------------
!
if(iq_geom==0) then
  kmem = (jk3-jk1+1) - (jk1-2*nkq-2)*(jk1-1)/2
else
  kmem = ikdif+1
end if
!
!
!------------------------------------------------------------------------------
!  check memory indexing
!------------------------------------------------------------------------------
!
if (amem > iamax) then
  ifnd = 0
  call q_error('e','MEMORY','Incorrect addres')
  goto 9999
end if
!
!-----------------------------------------------------------------------------
! retrieve info from reference locus in
! get actual number of valid points along locus (NLOCUSZ)
! depending on value of switch IQ_COMPACT
!------------------------------------------------------------------------------
!
nloc    = quad_nloc(kmem,amem)
nlocusx = nloc
!
!  short-cut when number of NON-ZERO points on locus is ZERO [27/8/2003]
!
if(nlocusx==0) goto 9999
!
r_ik2(1:nloc)  = quad_ik2(kmem,amem,1:nloc)
r_ia2(1:nloc)  = quad_ia2(kmem,amem,1:nloc)
r_ik4(1:nloc)  = quad_ik4(kmem,amem,1:nloc)
r_ia4(1:nloc)  = quad_ia4(kmem,amem,1:nloc)
!
r_w1k2(1:nloc) = quad_w1k2(kmem,amem,1:nloc)
r_w2k2(1:nloc) = quad_w2k2(kmem,amem,1:nloc)
r_w3k2(1:nloc) = quad_w3k2(kmem,amem,1:nloc)
r_w4k2(1:nloc) = quad_w4k2(kmem,amem,1:nloc)
!
r_w1k4(1:nloc) = quad_w1k4(kmem,amem,1:nloc)
r_w2k4(1:nloc) = quad_w2k4(kmem,amem,1:nloc)
r_w3k4(1:nloc) = quad_w3k4(kmem,amem,1:nloc)
r_w4k4(1:nloc) = quad_w4k4(kmem,amem,1:nloc)
!
r_zz(1:nloc)   = quad_zz(kmem,amem,1:nloc)
!
!------------------------------------------------------------------------------
kdif = ikmin - 1
if(iq_geom==0) then
  lambda = 1.
  kdif  = 0.
else
  lambda = q_kfac**(ikmin-1.)
end if
!
j_lambda = 1./sqrt(lambda)
c_lambda = lambda**6
!
!  compute combined scale factor
!
zz_lambda = lambda*c_lambda/j_lambda
!
!------------------------------------------------------------------------------
! select case to transform reference locus
!
! Transform of weigths reduces to an addition or subtraction
! because of log-spacing of wave numbers in the case of deep water
!
if(ik3 > ik1 .and. it3 >= it1) then      ! Case 1
  itrans = 1
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = r_ia2(1:nloc) + ibeta
  t_ia4(1:nloc)  = r_ia4(1:nloc) + ibeta
  idir   = 1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik3 > ik1 .and. it3 < it1)  then  ! Case 2
  itrans = 2
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = int(q_ad(ia1)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
  idir   = -1   ! according to theory
!  idir   = 1    ! as it should be to get symmetry
!
elseif(ik1 > ik3 .and. it3 >= it1) then      ! Case 3
  itrans = 3
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = int(q_ad(ia3)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
  idir   = 1
!
elseif(ik1 > ik3 .and. it1 > it3) then   ! Case 4
  itrans = 4
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = ibeta + r_ia4(1:nloc)
  t_ia4(1:nloc)  = ibeta + r_ia2(1:nloc)
  idir   = -1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik1==ik3 .and. it3 > it1) then  ! Case 5
  itrans = 5
  t_ik2(1:nloc)  = kdif + r_ik2(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik4(1:nloc)
  ibeta          = itmin-iaref
  t_ia2(1:nloc)  = r_ia2(1:nloc) + ibeta
  t_ia4(1:nloc)  = r_ia4(1:nloc) + ibeta
  idir   = 1
  t_w1k2(1:nloc)  = r_w1k2(1:nloc)
  t_w2k2(1:nloc)  = r_w2k2(1:nloc)
  t_w3k2(1:nloc)  = r_w3k2(1:nloc)
  t_w4k2(1:nloc)  = r_w4k2(1:nloc)
  t_w1k4(1:nloc)  = r_w1k4(1:nloc)
  t_w2k4(1:nloc)  = r_w2k4(1:nloc)
  t_w3k4(1:nloc)  = r_w3k4(1:nloc)
  t_w4k4(1:nloc)  = r_w4k4(1:nloc)
!
elseif(ik1==ik3 .and. it1 > it3) then  ! Case 6
  itrans = 6
  t_ik2(1:nloc)  = kdif + r_ik4(1:nloc)
  t_ik4(1:nloc)  = kdif + r_ik2(1:nloc)
  ibeta          = int(q_ad(ia1)/q_deltad+0.01)
  t_ia2(1:nloc)  = ibeta + 2.*iaref - r_ia2(1:nloc) -imirror
  t_ia4(1:nloc)  = ibeta + 2.*iaref - r_ia4(1:nloc) -imirror
!!  ibeta          = itmin-iaref
!!  t_ia2(1:nloc)  = r_ia4(1:nloc) + ibeta
!!  t_ia4(1:nloc)  = r_ia2(1:nloc) + ibeta
  idir   = -1
  t_w1k2(1:nloc)  = r_w3k2(1:nloc)
  t_w2k2(1:nloc)  = r_w4k2(1:nloc)
  t_w3k2(1:nloc)  = r_w1k2(1:nloc)
  t_w4k2(1:nloc)  = r_w2k2(1:nloc)
  t_w1k4(1:nloc)  = r_w3k4(1:nloc)
  t_w2k4(1:nloc)  = r_w4k4(1:nloc)
  t_w3k4(1:nloc)  = r_w1k4(1:nloc)
  t_w4k4(1:nloc)  = r_w2k4(1:nloc)
end if
!
t_zz(1:nloc)   = lambda*c_lambda/j_lambda * r_zz(1:nloc)
!
ifnd = 1
!
!------------------------------------------------------------------------------
!
9999 continue
!
call q_stack('-q_getlocus')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_init
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 25 Sep. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_fileio
use m_constants
use serv_xnl4v5
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!     25/02/1999 Initial version
!     13/10/1999 Error handling improved
!     18/10/1999 Test output to MATCHK.GRD added
!     21/10/1999 Extra output to MATCHK.GRD, iaref=1 for circle grids
!     01/11/1999 Allocatable arrays Q_XK and Q_SK added
!     14/02/2001 Version ready for WAVEWATCH III
!      8/08/2002 Release 4.
!     16/08/2002 Group velocity computed
!     22/08/2002 First and last used defined direction accounted for
!     11/09/2002 Call of Q_ALOC moved to higher level, viz. XNL_INIT
!                q_kpow initialized
!     25/09/2002 User defined directions used in the case of a sector grid
!
!  1. Purpose:
!
!     Initializing module for quadruplets
!     and setting default settings
!
!  2. Method
!
!     Conversion of power of spectral tail from E(f) to N(k) using the following
!     relations:
!
!       E(f) ~ f^qf_tail
!
!       N(k) ~ k^qk_tail
!
!       qk_tail = qf_tail/2 -1
!
!     See also Note 13 of G.Ph. van Vledder
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!
!  4. error meassage
!
!  5. Called by
!
!     XNL_INIT
!
!  6. Subroutines used
!
!     Q_STACK
!     Z_CMPCG
!     Z_STEPS
!     Z_WNUMB
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
!     /S Enable subroutine tracing
!
! 10. Source code
!------------------------------------------------------------------------------------------
!     Local variables
!
integer iaq,ikq    ! counters for loops over directions and wave numbers
real ff            ! frequency
!!!real z_wnumb       ! service function to compute wave number
!
integer iuerr      ! error indicator for i/o
!------------------------------------------------------------------------------
!
call q_stack('+q_init')
!
! set general settings
!
! convert power of E(f) f^qf_tail to power of N(k) k^qk_tail
! See Note 13 of G.Ph. van Vledder
!
qk_tail = (qf_tail-2.)/2. ! power of spectral tail, of N(k)
!
if(iq_prt >=2) then
  write(luq_prt,*)
  write(luq_prt,'(a,f6.1)') 'Q_INIT:  E(f)_tail: ',qf_tail
  write(luq_prt,'(a,f6.1)') 'Q_INIT:  N(k)_tail: ',qk_tail
end if
!
! set absolute and relative accuracies
!
eps_q   = 0.001           ! absolute accuracy for check of q==0
eps_k   = 1.e-5           ! absolute accuracy for equality check of k
rel_k   = 0.001           ! relative accuracy for equality check of k
!
sk_max = 50.              ! set maximum waver number
wk_max = real(nkq+0.9999) ! set maximum wave number index
!
! compute frequency and wave number grid
! assume that frequencies are always geometrically spaced,
! in the case of deep water this also holds for the wave numbers
!
q_ffac = (fqmax/fqmin)**real(1./(nkq-1.))         ! geometric spacing factor of frequencies
!
ff = fqmin                                        ! set minimum frequency
!
if(iq_prt>=2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Basic wave numbers, frequencies'
end if
!
do ikq=1,nkq                                       ! Generate wave number dependent variables
  q_f(ikq)    = ff                                 ! Frequency
  q_sig(ikq)  = ff*2.*pi                           ! Radian frequency
  q_k(ikq)    = z_wnumb(q_sig(ikq),q_depth,q_grav) ! compute wave number
  q_kpow(ikq) = (q_k(1)/q_k(ikq))**7.5             ! used in filtering
  ff          = ff*q_ffac                          ! Increase frequency
!
  call z_cmpcg(q_sig(ikq),q_depth,q_grav,q_cg(ikq))
  if(iq_prt >= 2) then
    write(luq_prt,'(a,i4,3f10.5,e12.4)') 'Q_INIT: ikq f sigma k k^p:', &
&    ikq,q_f(ikq),q_sig(ikq),q_k(ikq),q_kpow(ikq)
  end if
end do
!
! compute characteristics of extended k-array
!
if(iq_prt>=2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Extended wave numbers and spacing'
end if
!
do ikq=0,nkq
  if(ikq==0) then
    q_xk(ikq) = 0.
    q_sk(ikq) = q_k(1)
  elseif(ikq==nkq) then
    q_xk(ikq) = q_k(ikq)
    q_sk(ikq) = sk_max
  else
    q_xk(ikq) = q_k(ikq)
    q_sk(ikq) = q_k(ikq+1) - q_k(ikq)
  end if
!
end do
!
!
kqmin = q_k(1)
kqmax = q_k(nkq)
q_kfac = (kqmax/kqmin)**real(1./(nkq-1))  ! this value makes only sense in the
                                          ! case of deep water, IQ_DISP==1
!
! compute step size of frequency grids and wave number grid
!
call z_steps(q_f,  q_df,  nkq)           ! step size of frequencies
call z_steps(q_sig,q_dsig,nkq)           ! step size of radian frequencies
call z_steps(q_k,  q_dk,  nkq)           ! step size of wave numbers
!
if(iq_prt >= 2) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Q_INIT: Additional information'
  write(luq_prt,'(a,f8.1)')  'Q_depth (m):',q_depth
  write(luq_prt,'(a,i3)')    'Number of frequencies:',nkq
  write(luq_prt,'(a,f8.4)')  'Geometric f-spacing factor:',q_ffac
  write(luq_prt,'(a,f8.4)')  'Geometric k-spacing factor:',q_kfac
  write(luq_prt,'(a,2f8.3)') 'fmin fmax (Hz):',fqmin,fqmax
  write(luq_prt,'(a,2f8.3)') 'kmin kmax (Hz):',kqmin,kqmax
  write(luq_prt,*)
!
  write(luq_prt,*) '     i      f         df       sig      dsig       k         dk         cg'
!
  do ikq=1,nkq
    write(luq_prt,'(1x,i4,7f10.4)') &
 &  ikq,q_f(ikq),q_df(ikq),q_sig(ikq),q_dsig(ikq),q_k(ikq),q_dk(ikq),q_cg(ikq)
  end do
end if
!
! =============== D I R E C T I O N S ===============================================
!
! the directions in the array ANGLE are running from 1 to NAQ
! for a sector definition the middle direction has index IAREF
!
!  compute index IAREF of middle wave direction for sector grids
!
if(iq_grid ==1 .or. iq_grid==2) then
  iaref = (naq/2)+1
elseif(iq_grid==3) then
  iaref = 1
end if
!
if(iq_prt >= 2) write(luq_prt,'(a,i4)') &
&  'Q_INIT: Index of first direction for reference:',iaref
!
!  set loops indices
!
if(iq_grid==1) then    ! symmetric sector
  iaq1 = iaref
  iaq2 = naq
!
! non-symmetric sector and full circle
!
elseif(iq_grid==2 .or. iq_grid==3) then
  iaq1 = 1
  iaq2 = naq
end if
!
if(iq_prt >= 2) write(luq_prt,'(a,2i4)') &
&  'Q_INIT: Range of indices for loop over directions:',iaq1,iaq2
!
!  generate directions, given in degrees
!
q_sector = 0.5*(abs(q_dird1) + abs(q_dird2))
!
if(iq_grid==1 .or. iq_grid==2) then    ! define symmetric sector
  q_deltad = 2.*q_sector/real(naq-1.)  ! delta in degrees
  q_ang1 = -q_sector                   ! degrees
  q_ang2 =  q_sector                   ! degrees
  if(iq_prt>0) write(luq_prt,'(a)') 'Q_INIT: take care of q_dird1 and check if sector is OK'
!
elseif(iq_grid==3) then                ! full sector
  q_deltad = 360./real(naq)            ! degrees
  q_ang1 = 0                           ! degrees
  q_ang2 = 360.-q_delta                ! degrees
end if
!
q_delta = q_deltad*dera                ! directional step in radians
ncirc   = 2.00001*pi/q_delta           ! number of directions on circle
!
if(iq_prt >= 2) then
  write(luq_prt,'(a,3f10.3)')     'Q_INIT: d(1),d(n),dsector:',q_dird1,q_dird2,q_sector
  write(luq_prt,'(a,f6.2,a)')     'Q_INIT: Angular step     :',q_deltad,' degrees'
  write(luq_prt,'(a,2f8.2,i4,a)') 'Q_INIT: ang1 ang2 nang   :',q_ang1,q_ang2,naq,' degrees'
  write(luq_prt,'(a,i4)')         'Q_INIT: #Angles on circle:',ncirc
  write(luq_prt,*)
end if
!
!  generate directions arrays, given in degrees and radians
!
do iaq=1,naq
  q_ad(iaq) = q_ang1 + q_deltad*(iaq-1.)
  q_a(iaq)  = q_ad(iaq)*dera
  if(iq_prt >= 2) then
    write(luq_prt,'(a,i4,f10.4,f10.2)') 'Q_INIT: iaq q_a q_ad:',iaq,q_a(iaq),q_ad(iaq)
    if(iaq==naq) write(luq_prt,*)
  end if
end do
!
!  set loop indices for generation of grid
!  for sector grids and circle grids
!
if(iq_grid==1 .or. iq_grid==2) then
  iag1  = iaref
  iag2  = naq
!
!  circle grid
!
elseif(iq_grid==3) then
  iag1 = 1
  iag2 = naq/2+1
end if
!
iamax = iag2-iag1+1
!-------------------------------------------------------------------------
!
!
call q_stack('-q_init')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_locpos(ka,kb,km,kw,loclen)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 14 Oct. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5, only: z_root2
!
implicit none
!
!  0. Update history
!
!     Version Date        Description
!
!     03/12/1999  Initial version
!     09/08/2002  Upgrade to release 4.0
!     29/08/2002  Error handling z_root2 relaxed and some write statements modified
!     07/10/2002  Initialisation of QSQ replaced
!
!  1. Purpose:
!
!     Compute characteristics of locus used to optimize its acutal computation
!
!  2. Method
!
!  3. Parameter list:
!
!Type    I/O          Name     Description
!-----------------------------------------------------------------
real, intent (out) :: ka     ! minimum k along symmetry axis
real, intent (out) :: kb     ! maximum k along symmetry axis
real, intent (out) :: km     ! wave number at midpoint
real, intent (out) :: kw     ! half width of locus at midpoint
real, intent (out) :: loclen ! estimated length of locus
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!     z_zero2   Root finding method
!     x_locus1  Function of locus geometry, along symmetry axis
!     x_locus2  Function of locus geometry, perpendicular to symmetry axis
!     x_flocus  Locus function
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real kp           ! wave number at peak
real kpx,kpy      ! wave number at peak maximum
real zp           ! value of locus function at maximum
real za,zb        ! (test) value of locus function at kmin & kmax
real zz1,zz2      ! intermediate function values in interation process
real kk1,kk2      ! start values for finding root of locus equation
real kk1x,kk1y    ! wave number components at one side of root
real kk2x,kk2y    ! wave number components at other side of root
real beta1,beta2  ! parameters specifying cross component
real betaw        ! parameter specifying iterated cross component
real kwx,kwy      ! wave number at side of locus
real zw           ! function value at (kwx,kwy)
real a1,a2,b1,b2  ! constants in polynomial approximation of elliptic function
real aa,bb,mm,mm1 ! semi-major exis of ellips and derived parameters
!
real eps         ! local machine accuracy for determination of roots
real bacc        ! accuracy for determination of beta
real kacc        ! accuracy for determination of wave number roots
real qs          ! (w1-w3)/sqrt(g)
real qsq         ! gs^2
!
! Function declaration
!
!real z_root2     ! root finding using Ridders method
!
integer ierr     ! local error indicator, used in function Z-ZERO1
integer itest    ! local test level for test output
integer lutest   ! unit for test output in service routines
integer iter     ! local iteration number
integer maxiter  ! maximum number of iteration for determining starting points
!
!  function declarations
!!real, external :: x_locus2    ! locus function perpendicular to symmetry axis
!!real x_flocus                 ! 2-d locus function
!---------------------------------------------------------------------------------
!  assign test options
!
itest  = iq_test              ! assign test level
lutest = 0                    ! assign default, no test output in service routines
!
itest  = 0                    ! reset local test level
if(itest > 0) lutest=luq_tst   ! assign unit for test output
!
call q_stack('+q_locpos')
!
!  set initial values
!
eps     = epsilon(1.)         ! determine machine accurcy
maxiter = 20                  ! maximum number of iterations
!
! compute location of maximum, located at k_2 = P
!
kpx  = -px
kpy  = -py
kp   = sqrt(kpx**2 + kpy**2)
zp   = x_locus1(kp)
!
! find location of points A and B on locus
! for deep water, explicit relations are available
!
if(iq_disp==1) then
  qs = q/sqrtg
  qsq  = qs*qs
  if(qs < 0) then
    ka = 0.5*(-qs+sqrt(2.0*pmag-qsq))
    ka = ka**2
    kb = (pmag+qsq)/(2.*qs)
    kb = kb**2
    za = x_locus1(ka)
    zb = x_locus1(kb)
  else
    ka = 0.5*(-qs+sqrt(2.0*pmag-qsq))
    ka = -ka**2
    kb = (pmag-qsq)/(2.*qs)
    kb = kb**2
    za = x_locus1(ka)
    zb = x_locus1(kb)
  end if
!
!
!  find location of points A and B on locus
!  for water of finite depth, an iteration process is applied to
!  determine the zero-crossings of the locus function
!
else
!
  if(q<0) then
!
!   set two start points to locate position of wave number ka
!
    kk1 = 0.
    kk2 = kp
!
!   search root by Ridder's method
!
    kacc = 10.*max(kk1,kk2)*eps
    ka = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!
!
!   determine start points to locate position of wave number kb
!
    kk1 = kp
    kk2 = kp
    zz1 = zp
    zz2 = zp
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      call q_error('e','Start kb','Too many iterations needed')
      goto 9999
    end if
!
!   search root by Ridders method
!
    kacc = 10.*max(kk1,kk2)*eps
    kb = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!==================================================================
!   find positions for ka and kb for the case q > 0
!
  else
!
!   set two start points to locate position of wave number ka
!
    kk1  = 0.
    kk2  = -kp
    zz1  = x_locus1(kk1)
    zz2  = x_locus1(kk2)
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      call q_error('e','Start ka','Too many iterations needed')
      goto 9999
    end if
!
!   search root by Ridder's method
!
    kacc = 10.*max(abs(kk1),abs(kk2))*eps
    ka = z_root2(x_locus1,kk1,kk2,kacc,lutest,ierr)
!
!   determine start points to locate position of wave number kb
!
    kk1  = 0
    kk2  = kp
    zz1  = x_locus1(kk1)
    zz2  = x_locus1(kk2)
    iter = 0
!
!  ensure that two points are found on either side of zero-crossing
!
    do while (zz1*zz2 >= 0 .and. iter < maxiter)
      iter = iter + 1
      kk2 = kk2*2
      zz2 = x_locus1(kk2)
    end do
!
    if(iter>=maxiter) then
      call q_error('e','Start kb','Too many iterations needed')
      goto 9999
    end if
!
!   search root by Ridders method
!
    kacc = 10.*max(kk1,kk2)*eps
    kb = z_root2(x_locus1,kk1,kk2,kacc,luq_tst,ierr)
!
!   find positions for ka and kb for the case q > 0
!
  end if
!
  za = x_locus1(ka)
  zb = x_locus1(kb)
!
end if
!
! compute position of mid point
!
kmid = 0.5*(ka+kb)
km   = kmid
!
if(q < 0) then
  kmidx = kmid*cos(pang+pi)
  kmidy = kmid*sin(pang+pi)
else
  kmidx = kmid*cos(pang)
  kmidy = kmid*sin(pang)
end if
!
!
! compute width of locus near mid point of locus
!
! set starting values for determination of crossing point
!
beta1 = 0.
kk1x  = kmidx
kk1y  = kmidy
beta2 = 0.5
kk2x  = kmidx - beta2*py
kk2y  = kmidy + beta2*px
zz1   = x_flocus(kk1x,kk1y)
zz2   = x_flocus(kk2x,kk2y)
!
!
iter = 0
do while (zz1*zz2 > 0 .and. iter < maxiter)
  iter = iter + 1
  kk2x = kmidx - beta2*py
  kk2y = kmidy + beta2*px
  zz1  = x_flocus(kk1x,kk1y)
  zz2  = x_flocus(kk2x,kk2y)
  beta2 = beta2*2
end do
!
! call Ridders method to locate position of zero-crossing
!
!
bacc = 10.*max(beta1,beta2)*eps
betaw = z_root2(x_locus2,beta1,beta2,bacc,lutest,ierr)
!
!
kwx = kmidx - betaw*py
kwy = kmidy + betaw*px
zw  = x_flocus(kwx,kwy)
kw  = betaw*pmag
!
!
! estimate circumference of locus, assuming it to be an ellips
! estimate axis, this seems to be a rather good estimate
!
aa = 0.5*abs(ka-kb)
bb = kw
!
if (aa > bb) then
  mm = 1-(bb/aa)**2
else
  mm = 1-(aa/bb)**2
end if
!
mm1 = 1.-mm
a1 = 0.4630151;  a2 = 0.1077812;
b1 = 0.2452727;  b2 = 0.0412496;
!
if (mm1==0) then
  loclen = 4.*max(aa,bb)
else
  loclen = 4.*max(aa,bb)*((1. + a1*mm1 + a2*mm1**2) + (b1*mm1 + b2*mm1**2)*log(1/mm1))
end if
!
!
9999 continue
!
call q_stack('-q_locpos')
!
return
end subroutine
!
!------------------------------------------------------------------------------
subroutine q_makegrid
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 10 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
!
!  0. Update history
!
!     25/02/1999  Initial version
!     11/10/1999  Error handling improved; Bugs fixed when w1=w3
!     12/10/1999  Storage modified and non-geometric option included
!     16/10/1999  Equation for computing address of 2d array simplified
!     21/10/1999  Range of precomputed grid added to data file
!     22/10/1999  Renaming of some indices
!     25/10/1999  Header with grid info extended
!     12/11/1999  Output format modified of data to GRD file, adapted
!                 for use on UNIX systems at WES
!     08/12/1999  Interface with A_CMPLOC extended
!     28/12/1999  Routine A_CMPLOC renamed to Q_CMPLOC
!     03/01/2000  IQ_START replaced by IQ_LOCUS
!     05/01/2000  Interface with Q_CMPLOC modified
!     08/02/2000  Output to LUQLOC made conditional
!     09/08/2002  Name changed from Q_GRIDV1 to Q_MAKEGRID
!                 Upgrade to release 4.0
!     15/08/2002  Bug fixed in indexing bins below lowest wave number
!     20/08/2002  Sigma written to QUAD file, instead of wave numbers
!     22/08/2002  Data along locus compacted, elimate zero's
!     10/09/2002  Upgrade to release 5
!                 Value of LASTQUADFILE set
!     10/06/2003  Output to GRD file always without compacting
!
!  1. Purpose:
!
!     Set-up grid for computation of loci
!
!     Generate data file with basic loci for computation of
!     nonlinear quadruplet interactions
!
!  2. Method
!
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CTRGRID
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_CPMLOCUS
!     Q_MODIFY
!     Q_WEIGHT
!     Q_CHKRES
!     Q_NEAREST
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc,jloc            ! counters
integer iaq,ikq              ! counters
integer iaq3,ikq1,ikq3,nkq1  ! counters
integer jaq1,jaq3            ! counters
integer amem,kmem            ! index of angle and wave number in grid
real aa1,aa3,kk1,kk3         ! temporary wave number variables
!
integer nzloc                ! counter for non-zero contributions along locus
integer nztot1,nztot2        ! total number of zero and non-zero points on locus
integer ik2,ia2              ! index of wave number k2
integer ik4,ia4              ! index of wave number k4
!
real wk,wa                   ! weights
real w1k2,w2k2,w3k2,w4k2     ! interpolation weights
real w1k4,w2k4,w3k4,w4k4     ! interpolation weights
!
real ka,kb       ! lower and higher wave number magnitude
real km          ! wave number at mid point
real kw          ! half width of locus
!
real tfac        ! combined tail factor
!
logical lwrite   ! indicator if binary interaction grid has been written successfully
real smax        ! maximum s-value
!
real, allocatable :: xloc(:),yloc(:)
real qq
!-------------------------------------------------------------------------------
call q_stack('+q_makegrid')
!
! initializations
!
lwrite  = .false.
nztot1  = 0
nztot2  = 0
!%
quad_nloc = -1    ! number of points on all loci
!%
if(allocated(xloc)) deallocate(xloc) ; allocate (xloc(mlocus))
if(allocated(yloc)) deallocate(yloc) ; allocate (yloc(mlocus))
!
!  write header to grid file
!
!
!  set range of do loops for computing interaction grid
!
if(iq_geom==0 .or. iq_disp/=1) then
  nkq1 = nkq           ! loop over all k1 wave numbers, since no geometric scaling can be used
else
  nkq1 = 1             ! use only first wave number for k1, since geometric scaling can be used
end if
!
jaq1 = 1               ! index of direction of k1 in grid matrix
!-------------------------------------------------------------------------------------
!  compute components of reference wave number,
!  for setting up interaction grid
!-------------------------------------------------------------------------------------
k1: do ikq1=1,nkq1
!
  if(iq_screen==2) write(iscreen,*) 'k1-ring:',ikq1
!
  aa1   = q_ad(iaref)
  kk1   = q_k(ikq1)
  krefx = kk1*cos(q_ad(iaref)*dera)
  krefy = kk1*sin(q_ad(iaref)*dera)
!
  k1x  = krefx
  k1y  = krefy
!

k3: do ikq3 = ikq1,nkq   !
   if(iq_screen==2) write(iscreen,*) 'k1-k3 indices:',ikq1,ikq3
!
    kk3 = q_k(ikq3)
!
!
a3: do iaq3 = iag1,iag2
!
      if(iaq3 == iag1 .and. ikq3 == ikq1) cycle
!
      aa3 = q_ad(iaq3)
      k3x = kk3*cos(aa3*dera)
      k3y = kk3*sin(aa3*dera)
!------------------------------------------------------------------------------
!   compute locus for a specified combination of k1 and k3
!
!-----------------------------------------------------------------------------
      ia_k1 = iaq1; ik_k1 = ikq1
      ia_k3 = iaq3; ik_k3 = ikq3
      call q_cmplocus(ka,kb,km,kw,crf1)
!
      if(iq_err/=0) goto 9999
!------------------------------------------------------------------------------
!     redistibute or filter data points along locus
!
      call q_modify
      if(iq_err > 0) goto 9999
!------------------------------------------------------------------------------
!     compute weights for interpolation in computational grid
!
      call q_weight
      if(iq_err > 0) goto 9999
!------------------------------------------------------------------------------
!    special storing mechanism for interactions per combination of k1 and k3
!
      kmem  = (ikq3-ikq1+1) - (ikq1-2*nkq-2)*(ikq1-1)/2;
      jaq3  = iaq3-iaref+1        ! ensure that data stored in matrix start at index (1,1)
      amem  = jaq3                ! index of direction
!
!
!-------------------------------------------------------------------------------
!     Convert real indices to integer indexing and real weights
!
!    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
!    |    .      |              w2 = wk*(1-wa)
!    |. . + . . .| wa2   A      w3 = (1-wk)*wa
!    |    .      |       |      w4 = wk*wa
!    |    .      |       wa
!    |    .      |       |
!    1-----------2 ja2   V
!   jk2  wk2  jk2p
!
!    <-wk->
!
!-------------------------------------------------------------------------------
      nzloc = 0
!
loc:  do iloc = 1,nlocus
!
        ik2  = floor(wk_k2(iloc))
        ia2  = floor(wa_k2(iloc))
        wk   = wk_k2(iloc)-real(ik2)
        wa   = wa_k2(iloc)-real(ia2)
        w1k2 = (1.-wk)*(1.-wa)
        w2k2 = wk*(1.-wa)
        w3k2 = (1.-wk)*wa
        w4k2 = wk*wa
!
        ik4  = floor(wk_k4(iloc))
        ia4  = floor(wa_k4(iloc))
        wk   = wk_k4(iloc)-real(ik4)
        wa   = wa_k4(iloc)-real(ia4)
        w1k4 = (1.-wk)*(1.-wa)
        w2k4 = wk*(1.-wa)
        w3k4 = (1.-wk)*wa
        w4k4 = wk*wa
!
!  Take care of points that lie below lowest wave number
!  when no geometric scaling is applied, then modify weights
!  such that directional position is retained
!
        if(iq_geom==0) then
          if(ik2 ==0) then
            ik2  = 1
            w1k2 = w1k2 + w2k2
            w2k2 = 0.
            w3k2 = w3k2 + w4k2
            w4k2 = 0.
          end if
          if(ik4 ==0) then
            ik4  = 1
            w1k4 = w1k4 + w2k4
            w2k4 = 0.
            w3k4 = w3k4 + w4k4
            w4k4 = 0.
          end if
        end if
!
!  compute combined tail factor and product of coupling coefficient, step size,
!  symmetry factor, and tail factor divided by jacobian
!
        tfac = wt_k2(iloc)*wt_k4(iloc)
        quad_zz(kmem,amem,iloc)   = cple_mod(iloc)*ds_mod(iloc)*sym_mod(iloc)/jac_mod(iloc)*tfac
!
!----------------------------------------------------------------------------------------
!  compact data by elimating zero-contribution on locus
!----------------------------------------------------------------------------------------
!
        if(iq_compact==1 .and. abs(quad_zz(kmem,amem,iloc)) > 1.e-15) then
          nzloc = nzloc + 1
          jloc  = nzloc
          nztot1 = nztot1 + 1
        else
          jloc = iloc
        end if
        nztot2 = nztot2 + 1
!
!  shift data
!
        quad_zz(kmem,amem,jloc)  = quad_zz(kmem,amem,iloc)
!
        quad_ik2(kmem,amem,jloc) = ik2           ! lower wave number index of k2
        quad_ia2(kmem,amem,jloc) = ia2           ! lower direction index of k2
        quad_ik4(kmem,amem,jloc) = ik4           ! lower wave number index of k4
        quad_ia4(kmem,amem,jloc) = ia4           ! lower direction index of k4
!
        quad_w1k2(kmem,amem,jloc) = w1k2         ! weight 1 of k2
        quad_w2k2(kmem,amem,jloc) = w2k2         ! weight 2 of k2
        quad_w3k2(kmem,amem,jloc) = w3k2         ! weight 3 of k2
        quad_w4k2(kmem,amem,jloc) = w4k2         ! weight 4 of k2
!
        quad_w1k4(kmem,amem,jloc) = w1k4         ! weight 1 of k4
        quad_w2k4(kmem,amem,jloc) = w2k4         ! weight 2 of k4
        quad_w3k4(kmem,amem,jloc) = w3k4         ! weight 3 of k4
        quad_w4k4(kmem,amem,jloc) = w4k4         ! weight 4 of k4
!
!
      end do loc
!
      if(iq_compact==1) then
        quad_nloc(kmem,amem) = nzloc                ! store compacted number of points on locus
      else
        quad_nloc(kmem,amem) = nlocus               ! store number of points on locus
        nzloc = nlocus
      end if
!
!     write(luq_prt,'(a,4i5)') 'Q_MAKEGRID kmem amem nlocus:',kmem,amem,nlocus,nzloc
!
    end do a3
  end do k3
end do k1
!------------------------------------------------------------------------------
!  Write locus information to binary file
!------------------------------------------------------------------------------
!
write(luq_bqf) q_header
!
!------------------------------------------------------------------------------
! spectral interaction grid
!------------------------------------------------------------------------------
!
write(luq_bqf) naq,nkq
write(luq_bqf) q_sig
write(luq_bqf) q_ad
write(luq_bqf) iq_geom,iq_disp,iq_geom
write(luq_bqf) q_depth
!
!------------------------------------------------------------------------------
! interaction grid
!------------------------------------------------------------------------------
!
write(luq_bqf) quad_nloc
write(luq_bqf) quad_ik2
write(luq_bqf) quad_ia2
write(luq_bqf) quad_ik4
write(luq_bqf) quad_ia4
write(luq_bqf) quad_w1k2
write(luq_bqf) quad_w2k2
write(luq_bqf) quad_w3k2
write(luq_bqf) quad_w4k2
write(luq_bqf) quad_w1k4
write(luq_bqf) quad_w2k4
write(luq_bqf) quad_w3k4
write(luq_bqf) quad_w4k4
write(luq_bqf) quad_zz
!
!
lwrite = .true.
lastquadfile = bqname
!
if(iq_screen >= 1 .and. iq_test>=1) write(iscreen,'(2a)') 'q_makegrid: LASTQUADFILE: ',lastquadfile
!
9999 continue
!
if(allocated(xloc)) deallocate(xloc,yloc)
!
! check if BQF file has been written succesfully
! if not, deleted both the AQFILE and BQFILE
!
if(.not. lwrite) then
  close(luq_bqf,status='delete')
  if(iq_log > 0) then
    write(luq_log,*)
    write(luq_log,*) 'Q_MAKEGRID: Grid files ',trim(aqname),' and ',trim(bqname),' deleted'
    write(luq_log,*) 'Q_MAKEGRID: Since an error occurred during the generation'
    write(luq_log,*) 'Q_MAKEGRID: of the interaction grid'
  end if
end if
!-------------------------------------------------------------------------------
!  write statistics of compacting to print file
!
if(iq_prt >=1) then
  if(iq_compact==0) nztot1 = nztot2
  write(luq_prt,'(a,i10)') 'Total number of points on loci        :',nztot2
  write(luq_prt,'(a,i10)') 'Total number of stored points on locus:',nztot1
  write(luq_prt,'(a,i10)') 'Total number of zero points on locus  :',nztot2-nztot1
  write(luq_prt,'(a,f8.2)') 'Reduction factor (%):',real(nztot2-nztot1)/real(nztot2)*100.
end if
!
call q_stack('-q_makegrid')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_modify
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 11 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!       9/04/1999  Initial version
!      13/04/1999  New intermediate variables *_mod introduced
!      11/10/1999  Check on error messages in interpolation added
!      18/10/1999  Bug fixed in assigning new ds values to array DS_MOD
!      27/10/1999  Checked added on allocated of SOLD
!       8/12/1999  Test output added
!      29/12/1999  Bug fixed in assigning DS_MOD for first and last point on locus
!       1/10/2001  Components of k4-locus added
!                  No interpolation and modification if q==0
!       9/08/2002  Upgrade to version 4.0
!      15/08/2002  Step sizing improved
!       4/06/2003  Bug fixed in computing slen (length of locus)
!                  Locus closed to enable interpolation to finer resolution
!       6/06/2003  Activate output to XDIA configuration file
!      10/06/2003  Conversion to new indexing and lumping debugged
!      11/06/2003  Call to subroutine Q_SYMMETRY added
!
!  1. Purpose:
!
!     Modify points along the locus, such that they are evenly distributed
!     Only when intented, i.e. when IQ_LOCUS==2
!
!  2. Method
!
!     Compute new spacing along locus
!     Redistribute points and coefficient at new spacing using linear interpolation
!     Output DIA configuration when also lumping active
!
!     If no redistribution is needed, then copy relevant data
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_SYMMETRY
!     Z_INTP1
!
!  7. Remarks
!
!  8. structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local parameters
!
integer ierr,jerr    ! error indicators
integer nold,nnew    ! old and new number of points on locus
integer iold,inew    ! counter for loop along points
integer iloc         ! counter for loop along locus
integer jloc         ! counter for loop over lumped locus
integer itest        ! local test level, by default equal to IQ_TEST
!
real k2a,k2m         ! angle (deg) and wave number magnitude of wave number k2
real k4a,k4m         ! angle (deg) and wave number magnitude of wave number k4
real w2,w4           ! radian frequencies of wave numbers
!
!
real dk13,dk14       ! difference wave number
real dsnew,slen      ! new step size and length of locus
real zero            ! 0
real q_eps           ! accuracy to distinguish special case, with q=0
real diold           ! 'real' old number of indices between succeeding lumped bins
real dinew           ! 'real' new number of indices between succeeding lumped bins
!
!!real x_disper        ! evaluate dispersion relation
real, allocatable :: sold(:)     ! old coordinate along locus
real, allocatable :: snew(:)     ! new coordinate along locus
!--------------------------------------------------------------------------
call q_stack('+q_modify')
!
!  initialisations
!
zero  = 0.
q_eps = 1.e-5
itest = iq_test
!
! itest = 1   ! set local test level for test purposes
!
if(itest>=1) then
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_mod   :',iq_mod
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_xdia  :',iq_xdia
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_lump  :',iq_lump
  write(luq_tst,'(a,i4)') 'Q_MODIFY: iq_gauleg:',iq_gauleg
end if
!------------------------------------------------------------------------------
!  do not modify data when IQ_MOD==0
!------------------------------------------------------------------------------
!
if(iq_mod==0) then
  nlocus   = nlocus1
  x2_mod   = x2_loc
  y2_mod   = y2_loc
  x4_mod   = x4_loc
  y4_mod   = y4_loc
  s_mod    = s_loc
  ds_mod   = ds_loc
  jac_mod  = jac_loc
  cple_mod = cple_loc
  call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nlocus)
else
!------------------------------------------------------------------------------
! Modify spacing along locus
!------------------------------------------------------------------------------
  nold = nlocus1
!
! close locus by adding one point, equal to first point
! only for normal locus
!
  if(abs(q)>q_eps) nold  = nold+1
!
!------------------------------------------------------------------------------
! Determine new number of points along locus
!------------------------------------------------------------------------------
!
  if(iq_gauleg > 0) then
    nnew = iq_gauleg
  elseif(iq_lump > 0) then
    nnew = iq_lump
  else
    nnew  = nlocus0
  end if
!
!
  allocate (sold(nold),snew(nnew))
!------------------------------------------------------------------------------
!  Compute circumference of locus, distinguish 2 case, open or closed
!------------------------------------------------------------------------------
!
  if(abs(q)<q_eps) then
    slen = s_loc(nold)
    sold = s_loc
  else
    slen = 0
    do iold=1,nold-1               ! loop length minus one, since locus is closed
      sold(iold) = s_loc(iold)
      slen = slen + ds_loc(iold)
    end do
!
!------------------------------------------------------------------------------
!  close locus by copying first value in last value
!------------------------------------------------------------------------------
!
    sold(nold)     = slen
    x2_loc(nold)   = x2_loc(1)
    y2_loc(nold)   = y2_loc(1)
    x4_loc(nold)   = x4_loc(1)
    y4_loc(nold)   = y4_loc(1)
    jac_loc(nold)  = jac_loc(1)
    cple_loc(nold) = cple_loc(1)
  end if
!
!------------------------------------------------------------------------------
! compute new spacing along loci and coordinates along locus
! Gauss-Legendre integration
!------------------------------------------------------------------------------
!
  if(iq_gauleg > 0) then
    if(iq_gauleg > nnew) stop 'Q_MODIFY: iq_gauleg > nlocus0'
    nnew = iq_gauleg
    call y_gauleg(zero,slen,snew,ds_mod,nnew)
!
  else
    if(abs(q)>q_eps) then
      dsnew  = slen/real(nnew)
      do inew=1,nnew
        snew(inew) = (inew-1.)*dsnew
      end do
    else
      dsnew  = slen/real(nnew-1.)
      do inew=1,nnew
        snew(inew) = (inew-1)*dsnew
      end do
    end if
    ds_mod = dsnew
  end if
!
!
  jerr = 0
!------------------------------------------------------------------------------
! Compute characteristics of locus for special case q=0
!------------------------------------------------------------------------------
!
  if(abs(q)<1.e-5) then
    call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
 !
    call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
!
! ---  lumping along locus --------------------------------------------
!
    if(iq_lump>0) then
      diold  = slen/real(nold)
      dinew  = slen/real(nnew)
      ds_mod = 0.
      call q_symmetry(k1x,k1y,k3x,k3y,x4_loc,y4_loc,sym_loc,nold)
!
      do iloc=1,nlocus1
        jloc = floor((iloc-1.)*diold/dinew)+1
        ds_mod(jloc)   = ds_mod(jloc) + cple_loc(iloc)*ds_loc(iloc)/jac_loc(iloc)*sym_loc(iloc)
        jac_mod(jloc)  = 1.
        cple_mod(jloc) = 1.
      end do
!
      sym_mod = 1          ! symmetry already taken account in lumping proces
!
! --- No lumping -------------------------------------------------------------
!
    else
      call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
      jerr = jerr + ierr
!
      call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
      jerr = jerr + ierr
    end if
!------------------------------------------------------------------------------------------------
!  compute characteristics for closed locus
!------------------------------------------------------------------------------------------------
  else
    call z_intp1(sold,x2_loc,snew,x2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y2_loc,snew,y2_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,x4_loc,snew,x4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 x_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,y4_loc,snew,y4_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 y_loc, ierr=',ierr
    jerr = jerr + ierr
!
    call z_intp1(sold,s_loc,snew,s_mod,nold,nnew,ierr)
    if(ierr > 0) write(luq_err,*) 'Z_INTP1 s_loc, ierr=',ierr
    jerr = jerr + ierr
!
!
    call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym_mod,nnew)
!
!  ----- Lumping along locus -----------------------------------
!
    if(iq_lump>0) then
      diold  = slen/real(nold-1)
      dinew  = slen/real(nnew)
      ds_mod = 0.
      call q_symmetry(k1x,k1y,k3x,k3y,x4_loc,y4_loc,sym_loc,nold)
!
      do iloc=1,nold-1
        jloc = floor((iloc-1.)*diold/dinew + 1.49999)
        jloc = mod(jloc-1+nnew,nnew)+1
        ds_mod(jloc)   = ds_mod(jloc) + cple_loc(iloc)*ds_loc(iloc)/jac_loc(iloc)*sym_loc(iloc)
        jac_mod(jloc)  = 1.
        cple_mod(jloc) = 1.
      end do
!
      sym_mod = 1          ! symmetry already taken account in lumping proces
!
!------------  No lumping along locus  --------------------------------
!
    else
      call z_intp1(sold,jac_loc,snew,jac_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 jac_loc, ierr=',ierr
      jerr = jerr + ierr
!
      call z_intp1(sold,cple_loc,snew,cple_mod,nold,nnew,ierr)
      if(ierr > 0) write(luq_err,*) 'Z_INTP1 cp_loc, ierr=',ierr
      jerr = jerr + ierr
    end if
!
    if(jerr > 0) then
      iq_err = iq_err + 1
      call q_error('e','INTER','Problem in interpolation process')
      goto 9999
    end if
  end if
!
  nlocus = nnew
!
end if
!
!------------------------------------------------------------------------------
!
!
!------------------------------------------------------------------------------
!
!!  compute symmetry factor for reducing computational load
!!
!!call q_symmetry(k1x,k1y,k3x,k3y,x4_mod,y4_mod,sym,nnew)
!!
do iloc=1,nlocus
  k2x = x2_mod(iloc)
  k2y = y2_mod(iloc)
  k4x = x4_mod(iloc)
  k4y = y4_mod(iloc)
!
  k2m = sqrt(k2x**2 + k2y**2)
  k4m = sqrt(k4x**2 + k4y**2)
  k2a = atan2(k2y,k2x)*rade
  k4a = atan2(k4y,k4x)*rade
!
  k2m_mod(iloc) = k2m
  k4m_mod(iloc) = k4m
  k2a_mod(iloc) = k2a
  k4a_mod(iloc) = k4a
!
!
!
end do
!
!
9999 continue
!
if(allocated(sold)) deallocate(sold,snew)
!
call q_stack('-q_modify')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_polar2(kmin,kmax,kx_beg,ky_beg,kx_end,ky_end,loclen,ierr)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 8 Aug. 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5, only: z_wnumb
!
implicit none
!
!  0. Update history
!
!     Date        Description
!
!     03/12/1999  Initial version
!     09/08/2002  Geometric spacing of k added
!                 Upgrade to release 4.0
!     13/08/2002  reorganisation of loops generating points on locus
!     08/08/2003  Check included for maximum number of IPOL by using MPOL
!                    MPOL=MLOCUS/2+1-1  (-1 added regarding IPOL=IPOL+1 in Q_MODIFY)
!                 Check included on ARG=0 for IQ_LOCUS=2 and parameter dke added
!
!  1. Purpose:
!
!     Compute position of locus for given k1-k3 vector
!
!  2. Method
!
!     Explicit polar method, see Van Vledder 2000, Monterey paper
!     Optionally using a fixed k-step, geometric k-step or adaptive stepping
!
!  3. Parameters used:
!
!Type    I/O        Name             Description
!------------------------------------------------------------------------------
real, intent(in) :: kmin           ! minimum wave number on locus
real, intent(in) :: kmax           ! maximum wave number on locus
real, intent(in) :: kx_beg         ! x-coordinate of begin point
real, intent(in) :: ky_beg         ! y-coordinate of begin point
real, intent(in) :: kx_end         ! x-coordinate of end point
real, intent(in) :: ky_end         ! y-coordinate of end point
real, intent(in) :: loclen         ! estimated length of locus
integer, intent (out)  :: ierr     ! error condition
!
!     Parameters with module
!
!     nlocus0   Preferred number of points on locus
!     q         w1-w3, difference of radian frequencies
!     pmag      |k1-k3| (vector form)
!     pdir      direction of difference vector k1-k3
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CPMLOCUS
!
!  6. Subroutines used:
!
!     X_COSK
!
!  7. Remarks
!
!     The type of locus computation is controlled by the parameter IQ_LOCUS
!     Set in Q_SETCFG
!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variabels
!
integer ipol      ! counter
integer jpol      ! counter
integer iend      ! indicates end of locus computation
integer ipass     ! counter for passes
integer npol      ! number of points on locus
integer npass     ! number of passes in iteration process
integer mpol      ! maximum number of points on locus, related to MLOCUS
!
real kold         ! temporary wave number
real knew         ! temporary wave number
real cosold       ! 'old' cosine of angle
real cosnew       ! 'new' cosine of angle
real dkpol        ! step in wave number
real dkold        ! 'old' step in wave number
real ang1         ! 'old' angle
real ang2         ! 'new' angle
real kk1          ! 'old' wave number
real kk2          ! 'new' wave number
real kratio       ! ratio between succesive k-values when IQ_LOCUS=3
real arg          ! argument
real dk           ! step in wave number
real dke          ! estimate of new dk
real dsnew        ! new step size along locus
real dsz          ! estimated step size along locus
!
integer itest    ! local test level
integer lutest   ! unit number for test output in service routine
!
!  function declarations
!!!real    z_wnumb  ! compute wave number, via module SERV_XNL4V4
!!real    x_disper ! dispersion relation
!
!------------------------------------------------------------------------------
! initialisations
!------------------------------------------------------------------------------
call q_stack('+q_polar2')
!
ierr = 0                      ! set error code to zero
npol = (nlocus0+1)/2+1        ! first estimate of number k-values along symmetry axis
mpol = mlocus/2               ! set maximum number of points along locus axis
!
!-------------------------------------------------------------------------------
!
select case(iq_locus)
!------------------------------------------------------------------------------
! CASE = 1: Linear spacing of wave numbers along symmetry axis
!------------------------------------------------------------------------------
  case(1)
!
  dk = (kmax-kmin)/real(npol-1)
  do ipol=1,npol
    k_pol(ipol) = kmin + (ipol-1)*dk
    c_pol(ipol) = x_cosk(k_pol(ipol))
  end do
!------------------------------------------------------------------------------
!  Case = 2: Variable k-stepping along symmetry axis,
!            such that step along locus is more or less constant
!------------------------------------------------------------------------------
  case(2)
!
! set first point on locus
!
  ipol        = 1
  k_pol(ipol) = kmin
  c_pol(ipol) = -1.
  kold        = kmin
  cosold      = -1.
!
! compute initial step size of polar wave number
!
  dk0   = (kmax - kmin)/real(npol)      ! estimate of step size of equidistant radii
  dsz   = loclen/real(nlocus0)          ! estimate of step size along locus
  npass = 3                             ! set number of passes in iteration
  dk0   = dk0/2                         ! reduce initial step
  dk    = dk0
  iend  = 0
!
!
  do while (k_pol(ipol) < kmax .and. iend==0 .and. ipol < mpol)
    do ipass=1,npass
      knew  = min(kmax,k_pol(ipol)+dk)
      dkold = knew - k_pol(ipol)
      cosnew = x_cosk(knew)
      ang1  = pang + acos(cosold)
      ang2  = pang + acos(cosnew)
      kk1   = kold
      kk2   = knew
      arg   = kk1**2 + kk2**2 -2.*kk1*kk2*cos(ang1-ang2)
      dsnew = sqrt(abs(arg))
      if(dsnew>0) dke   = dk*dsz/dsnew
      dk    = dke
    end do
!----------------------------------------------------------------------------------------------
!  assign new estimate and check value of IPOL
!----------------------------------------------------------------------------------------------
    ipol        = ipol + 1
    k_pol(ipol) = k_pol(ipol-1) + dkold
    c_pol(ipol) = cosnew
    kold        = knew
    cosold      = cosnew
    if (abs(dkold) < 0.0005*(kmax-kmin)) iend=1
  end do
!
! fill last bin with coordinates of end point
!
  if(k_pol(ipol) < kmax .and. ipol <  mpol) then
    ipol = ipol + 1
    c_pol(ipol) = -1.
    k_pol(ipol) = kmax
  end if
!
!  update the number of k-points on symmetry axis
!
  npol = ipol
!
!-------------------------------------------------------------------------------
!  Case 3: Geometric spacing of wave numbers along symmetry axis
!-------------------------------------------------------------------------------
  case(3)
  kratio = (kmax/kmin)**(1./(npol-1.))
  do ipol=1,npol
    k_pol(ipol) = kmin*kratio**(ipol-1.)
    c_pol(ipol) = x_cosk(k_pol(ipol))
  end do
!
end select
!
!------------------------------------------------------------------------------
!
!  compute actual number of points on locus
!  this will always be an even number
!  mirror image the second half of the locus
!
nlocus1 = 2*npol-2
!
a_pol(1) = pang + acos(c_pol(1))
c_pol(1) = cos(a_pol(1))
!
do ipol=2,npol
  jpol = 2*npol-ipol
  a_pol(ipol) = pang + acos(c_pol(ipol))
  a_pol(jpol) = pang - acos(c_pol(ipol))
  c_pol(jpol) = cos(a_pol(jpol))
  k_pol(jpol) = k_pol(ipol)
end do
!
! compute x- and y-position along locus
!
do ipol=1,nlocus1
  x2_loc(ipol) = k_pol(ipol)*cos(a_pol(ipol))
  y2_loc(ipol) = k_pol(ipol)*sin(a_pol(ipol))
end do
!
!
9999 continue
!
call q_stack('-q_polar2')
!
return
end subroutine
!-----------------------------------------------------------------------------------
subroutine q_setconfig(iquad)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 16 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_fileio
use serv_xnl4v5
!--------------------------------------------------------------------------------
!
implicit none

!
!  0. Update history
!
!     20/07/1999  Initial version
!     11/10/1999  Option iq_geom added, consistency checks added
!     15/10/1999  Option iq_trf added, keyword TRANSF
!     02/11/1999  Close(LUQCFG) added, implicit none added
!     08/12/1999  Option IQ_MOD added
!     24/12/1999  Extra output when *.cfg does not exist, and IQ_PRT=1
!     08/02/2000  Error message included for IQUAD
!     12/05/2002  Triplet settings added
!     08/08/2002  Upgrade to release 4
!     19/08/2002  Inclusion of various test option and interpolation
!     22/08/2002  Switch to compact data included
!     09/09/2002  Parameter q_dstep added
!     11/09/2002  Parameter qf_frac added
!     26/05/2003  Parameter iq_lump added
!     04/06/2003  Parameter IQ_INT renamed IQ_INTEG
!                 Switch IQ_GAULEG added
!     11/06/2003  name changed from Q_SETCFG to Q_SETCONFIG
!                 Parameter IQ_SPACE removed
!     13/06/2003  Set test output, from XNL_INIT
!     16/06/2003  Switch IQ_SYM added
!     09/09/2003  Variable ID_FACMAX added
!
!  1. Purpose:
!
!     Set settings for computing the nonlinear interactions
!     set optimal basic settings
!     Set some settings based on the value of IQUAD
!
!  2. Method
!
!     Based on the value of IQUAD a number of settings are preset
!     In the case the file [qbase].CFG exists, this file
!     is analyzed and possibly some settings are reset
!
!  3. Parameter list:
!
!Type, I/O               Name      Description
!--------------------------------------------------------------------------
integer, intent(in) ::   iquad   ! Indicator for a specific choice of
!                                  settings for computing the nonlinear
!                                  interactions
!  4. Error messages
!
!  5. Called by:
!
!     XNL_INIT
!
!  6. Subroutines used
!
!  7. Remarks
!
!     IF no valid value for iquad is given, a default choice is
!     specified
!
!     The various options of the setting are specified in the general quads module
!
!  8. Structure
!
!  9. Switches
!
!     /S Enable subroutine tracing
!
! 10. Source code
!--------------------------------------------------------------------------------
! Local variables
!
integer iend               ! indicator for end of file
integer iuerr              ! error status of file io
character(len=10) cpar     ! character parameter
real rpar                  ! real parameter
!--------------------------------------------------------------------------------
!
call q_stack('+q_setconfig')
!--------------------------------------------------------------------------------
! default settings, which always work
!--------------------------------------------------------------------------------
nlocus0    = 30           ! Preferred number of points along locus
id_facmax  = 2            ! Factor for depth search in Q_SEARCHGRID
iq_filt    = 1            ! switch filtering on
iq_gauleg  = 0            ! No Gauss-Legendre interpolation
iq_locus   = 2            ! polar method, constant step with adaptive stepping
iq_lump    = 0            ! lumping disabled
iq_make    = 1            ! make grid at each new run
iq_mod     = 1            ! Modify spacing to equidistant spacing of points along locus
iq_compact = 1            ! Do not (yet) compact data along locus
iq_interp  = 1            ! bi-linear interpolation
iq_lump    = 0            ! no lumping of coefficient along locus
iq_search  = 0            ! No search is carried out for nearest quad grid
iq_sym     = 1            ! Activate symmetry reduction
!--------------------------------------------------------------------------------
!  set settings for special or test purposes
!--------------------------------------------------------------------------------
iq_trf     = 0            ! No output of transformed loci
iq_t13     = 0            ! No test output of T13 integration
!-------------------------------------------------------------------------------
! set filtering values for retricting integration space
!-------------------------------------------------------------------------------
qf_krat = 2.5             ! maximum ratio between wave numbers k1 and k3
qf_dmax = 75.0            ! difference in degrees between k1 and k3
qf_frac = 0.1             ! fraction of maximum energy density
!
q_sector = 120.           ! set size of half-plane direction sector (120)
!------------------------------------------------------------------------------
!
! Set specific parameter depending on IQUAD
!
!------------------------------------------------------------------------------
! deep water test version
!
if(iquad==1) then
  iq_geom   = 0            ! apply geometric scaling (Geometric scaling is disabled)
  iq_dscale = 0            ! no depth scaling
  iq_disp   = 1            ! deep water
  iq_cple   = 1            ! Webb's coupling coefficient
!
! 'deep' water computation and HH/WAM depth scaling
!
elseif(iquad==2) then
  iq_geom   = 0            ! apply geometric scaling
  iq_dscale = 1            ! put depth scaling on
  iq_disp   = 1            ! deep water
  iq_cple   = 1            ! Webb's coupling coefficient
!
!  full finite depth computation of interactions
!
elseif(iquad==3) then
  iq_dscale = 0            ! no depth scaling
  iq_disp   = 2            ! finite depth dispersion relation
  iq_geom   = 0            ! no geometric scaling
  iq_cple   = 2            ! finite depth coupling coefficient of H&H
else
  if(iq_screen>0) write(iscreen,'(a,i4)') 'Q_SETCONFIG: iquad=',iquad
  call q_error('e','IQUAD','No valid value of iquad has been given, default settings')
  write(luq_err,'(a,i4)') 'Q_SETCONFIG: Value of IQUAD:',iquad
  goto 9999
end if
!-------------------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------
! check if the configuration exists,
! and if so, override the settings
!----------------------------------------------------------------------------------
!
call z_fileio(trim(qbase)//'.cfg','OF',iufind,luq_cfg,iuerr)
if(luq_cfg > 0) then
  if(iq_log >= 1) then
    write(luq_log,*)
    write(luq_log,'(a)') 'Q_SETCONFIG: Configuration file '//trim(qbase)//'.cfg has been found'
    write(luq_log,'(a,i4)') 'Q_SETCONFIG: '//trim(qbase)//'.cfg connected to :',luq_cfg
  end if
!
  iend = 0
!
  do while (iend==0)
    read(luq_cfg,*,iostat=iend) cpar,rpar
!
    call z_upper(cpar) ! Convert string to upper case
!
    if(iend==0) then   ! process the command
!
      if(trim(cpar)=='DEPTH')    q_depth  = rpar
      if(trim(cpar)=='DSTEP')    q_dstep    = rpar
      if(trim(cpar)=='F_DMAX')   qf_dmax  = rpar
      if(trim(cpar)=='F_KRAT')   qf_krat  = rpar
      if(trim(cpar)=='F_FRAC')   qf_frac  = rpar
      if(trim(cpar)=='FMIN')     fqmin    = rpar
      if(trim(cpar)=='FMAX')     fqmax    = rpar
      if(trim(cpar)=='NLOCUS')   nlocus0    = int(rpar)
      if(trim(cpar)=='SECTOR')   q_sector   = rpar
!
      if(trim(cpar)=='GEOM') then
        iq_geom  = int(rpar)
        if(iq_geom==1) then
          iq_geom=0
          if(iq_screen>0) write(iscreen,'(a)') 'Q_SETCONFIG: geometric scaling disabled'
          if(iq_prt>=1)   write(luq_prt,'(a)') 'Q_SETCONFIG: geometric scaling disabled'
        end if
      end if
      if(trim(cpar)=='COMPACT')  iq_compact = int(rpar)
      if(trim(cpar)=='COUPLING') iq_cple    = int(rpar)
      if(trim(cpar)=='DISPER')   iq_disp    = int(rpar)
      if(trim(cpar)=='FILT')     iq_filt    = int(rpar)
      if(trim(cpar)=='GAULEG')   iq_gauleg  = int(rpar)
      if(trim(cpar)=='GRID')     iq_grid    = int(rpar)
      if(trim(cpar)=='INTEG')    iq_integ   = int(rpar)
      if(trim(cpar)=='INTERP')   iq_interp  = int(rpar)
      if(trim(cpar)=='LOCUS')    iq_locus   = int(rpar)
      if(trim(cpar)=='LOGGING')  iq_log     = int(rpar)
      if(trim(cpar)=='LUMPING')  iq_lump    = int(rpar)
      if(trim(cpar)=='MAKE')     iq_make    = int(rpar)
      if(trim(cpar)=='MODIFY')   iq_mod     = int(rpar)
      if(trim(cpar)=='PRINT')    iq_prt     = int(rpar)
      if(trim(cpar)=='PRINT')    iq_prt     = int(rpar)
      if(trim(cpar)=='SCREEN')   iq_screen  = int(rpar)
      if(trim(cpar)=='SEARCH')   iq_search  = int(rpar)
      if(trim(cpar)=='SYM')      iq_sym     = int(rpar)
      if(trim(cpar)=='T13')      iq_t13     = int(rpar)
      if(trim(cpar)=='TEST')     iq_test    = int(rpar)
      if(trim(cpar)=='TRACE')    iq_trace   = int(rpar)
      if(trim(cpar)=='TRANSF')   iq_trf     = int(rpar)
      if(trim(cpar)=='XDIA')     iq_xdia    = int(rpar)
    end if
  end do
!
  close(luq_cfg)
!
  if(iq_log >= 1) write(luq_log,'(a,i4)') &
& 'Q_SETCONFIG: '//trim(qbase)//'.cfg disconnected from :',luq_cfg
!
else
!  iq_prt = 1
  if(iq_log >= 1) then
    write(luq_log,*)
    write(luq_log,'(a)') 'Q_SETCONFIG: Configuration file '//trim(qbase)//'.CFG has not been found'
  end if
end if
!
9999 continue
!
call q_stack('-q_setconfig')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_searchgrid(depth,igrid)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 September 2003
!   +---+ |   |  Release: 5.03
!         +---+
!
! do not use m_xnldata
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     Version  Date    Modification
!
!     20/08/2002  Initial version
!     29/08/2002  Write statements made conditionsl
!      5/09/2003  Search algorithm improved
!     09/09/2003  factor ID_FACMAX introduced and extra test output created
!                 Input water depth saved for output
!
!  1. Purpose:
!
!     Search nearest valid grid, read grid file and scale factor
!
!  2. Method
!
!     Using the actual water depth
!     all possible interaction grids are checked
!     in upward and downward direction
!
!  3. Parameters used
!
real, intent(in)     :: depth  !  depth for which grid file must be found
integer, intent(out) :: igrid  !  status of grid checking
!                                 ==0: a proper grid exists
!                                 ==1: grid file does not exist
!                                 ==2: grid file exists, but it is incorrect
!                                 ==3: read error in accessing grid information
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_XNL4V4
!
!  6. Subroutines used
!
!     Q_CTRGRID
!     Q_STACK
!
!  7. Remarks
!
!
!  8. Structure
!
!
!  9. Switches
!
! 10. Source code
!---------------------------------------------------------------------------
!     Local variables
!
integer id        ! counter
integer idepth    ! integer depth
integer id_upper  ! upper limit of search
integer id_lower  ! lower limit of depth search
!
real d_lower      ! lower valid depth
real d_upper      ! upper valid depth
real r_lower      ! ratio with lower valid depth
real r_upper      ! ratio with upper valid depth
real s_depth      ! target depth in m, saved in this variable
real dfac1,dfac2  ! depth scale factors
real eps          ! accuracy
!------------------------------------------------------------------------------
!
call q_stack('+q_searchgrid')
!
eps = 0.0001
!
!------------------------------------------------------------------------------
!  check if a depth exists for current grid
!------------------------------------------------------------------------------
!
!
q_depth = depth + eps

call q_ctrgrid(1,igrid)
!
!
if(igrid==0) then
  if(iq_screen>=1) write(iscreen,'(a)') 'Q_SEARCHGRID: grid accepted, read whole database'
!
  call q_ctrgrid(2,igrid)
  goto 9999
end if
!
! save depth for which nearest grid file is to be found
!
s_depth  = depth
idepth   = int(s_depth*10+eps)
id_lower = int(q_mindepth*10+eps)
id_upper = int(q_maxdepth*10+eps)
!
id_upper = min(id_facmax*idepth,id_upper)
!
!  set 'not found' condition
!
d_lower = -1.
d_upper = -1.
!
!------------------------------------------------------------------------------
! search downwards until a valid grid is found
!------------------------------------------------------------------------------
!
do id = idepth-1,id_lower,-1
  q_depth = real(id)/10.+eps

  call q_ctrgrid(1,igrid)


  if(igrid==0) then
    d_lower = q_depth
    exit
  end if
end do
!
!------------------------------------------------------------------------------
!  seach upwards until a valid grid is found
!------------------------------------------------------------------------------
!
do id = idepth+1,id_upper
  q_depth = real(id)/10.+eps


  call q_ctrgrid(1,igrid)


  if(igrid==0) then
    d_upper = q_depth
    exit
  end if
end do
if(iq_prt>=1) write(luq_prt,*)
!------------------------------------------------------------------------------
!
!  determine nearest grid
!------------------------------------------------------------------------------
!
if(d_lower > 0) then
  r_lower = s_depth/d_lower
else
  r_lower = -1.
end if
!
if(d_upper > 0) then
  r_upper = d_upper/s_depth
else
  r_upper = -1.
end if
!
if(iq_prt>=1) then
  write(luq_prt,'(a,3f8.2)') 'Q_SEARCHGRID: d_lower d_target d_upper      :',d_lower,s_depth,d_upper
  write(luq_prt,'(a,2f8.2)') 'Q_SEARCHGRID: r_lower r_upper               :',r_lower,r_upper
end if
!------------------------------------------------------------------------------
!  select nearest valid grid
!------------------------------------------------------------------------------
if(r_lower>0 .and. r_upper>0) then
  if(r_lower < r_upper) then
    q_depth = d_lower
  else
    q_depth = d_upper
  end if
!
elseif(r_lower > 0 .and. r_upper <0 ) then
  q_depth = d_lower
elseif(r_lower < 0 .and. r_upper > 0) then
  q_depth = d_upper
else
  call q_error('e','SEARCHGRID','No valid nearest grid could be found')
  goto 9999
end if
!
!-----------------------------------------------------------------------------------------------
! compute depth scaling factors
!------------------------------------------------------------------------------
!
call q_dscale(a,q_sig,q_a,nkq,naq,s_depth,q_grav,dfac1)
call q_dscale(a,q_sig,q_a,nkq,naq,q_depth,q_grav,dfac2)
!
q_scale = dfac1/dfac2
!
if(iq_prt>=1) then
  write(luq_prt,'(a,2f8.4)') 'Q_SEARCHGRID: target and nearest scale factors:',dfac1,dfac2
  write(luq_prt,'(a,f8.4)')  'Q_SEARCHGRID: compound scale factor           :',q_scale
end if
!
!  Read BQF for nearest valid water depth
!
call q_ctrgrid(2,igrid)
if(iq_prt>=2) then
  write(luq_prt,'(a,f12.2)') 'Q_SEARCHGRID: Q_CTRGRID called with depth:',q_depth
  write(luq_prt,'(a,i4)') 'Q_SEARCHGRID: igrid of nearest grid operation:',igrid
end if
!
9999 continue
!
!  restore water depth
!
q_depth = s_depth
!
call q_stack('-q_searchgrid')
!
return
end subroutine
!-----------------------------------------------------------------
subroutine q_setversion
!-----------------------------------------------------------------
! do not use m_xnldata
!-----------------------------------------------------------------
! This subroutine has automatically been written by MODULE5
! Author: Gerbrant van Vledder
!
q_version ='GurboQuad  Version: 5.03 Build: 59 Date: 2003/09/15 [S]'
!
! Source code options:S
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_stack(mod_name)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 11 June 2003
!   +---+ |   |  Release: 5
!         +---+
!
! do not use m_xnldata
use m_fileio
implicit none
!
!
!  0. Update history
!
!     20/07/1999  Initial version
!     13/10/1999  Error handling improved
!     08/08/2002  Upgrade to release 4
!     11/06/2003  Extra check on output to print or test file
!
!  1. Purpose:
!
!     Add or remove mod_name name from module stack
!
!  2. Method
!
!     mod_name must be preceeded by a '+' , '-'
!     The module name is pushed to the stack when preceeded by '+'
!     and removed if mname starts with '-'.
!     In case an error is active,the module name is not removed
!     from the stack if mname starts with a '-'.The module is
!     always removed from the stack if mname starts with '!'.
!
!
!  3. Parameter list:
!
!Type           I/O          name          description
!-------------------------------------------------------
character(len=*), intent(in) :: mod_name    ! module name
!
!  4. Error messages
!
!  5. Called by
!
!     All q_** routines
!
!  6. Subroutines used
!
!     q_error
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!-------------------------------------------------------------------------------------
character(len=1) mod_task       ! task to do
integer mod_len                 ! length of mod_name
!
!!\A
if(iq_trace > 0) then
  if(iq_prt>0)  write(luq_prt,'(2a)') 'TRACE -> ',trim(mod_name)
  if(iq_test>0) write(luq_tst,'(2a)') 'TRACE -> ',trim(mod_name)
  if(iq_screen >= 2) write(iscreen,'(2a)') 'TRACE -> ',trim(mod_name)
end if
!
!  split MOD_NAME in two parts
!
!        MOD_TASK '+','-'
!
mod_len  = len_trim(mod_name)
mod_task = mod_name(1:1)
sub_name = mod_name(2:mod_len)
!
if(mod_task(1:1) == '+') then
  iq_stack = iq_stack + 1
!
  if(iq_stack > mq_stack) then
    call q_error('e','STACKMAX',' ')
    goto 9999
  else
    cstack(iq_stack) = mod_name(2:mod_len)
  end if
!------------------------------------------------------------------------
!  remove name from stack
!------------------------------------------------------------------------
elseif(mod_task(1:1) == '-') then
!
  if(mod_name(2:mod_len) == cstack(iq_stack)) then
    iq_stack = iq_stack - 1
  else
    write(luq_err,'(a)') 'Module name:',mod_name
    call q_error('e','STACKNAME',' ')
    goto 9999
  end if
else
  call q_error('e','STACKCALL',' ')
  goto 9999
end if
!
!!\Z
!
9999 continue
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_summary
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 16 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_fileio
use serv_xnl4v5
!--------------------------------------------------------------------------------
!
implicit none

!
!  0. Update history
!
!     11/06/2003  Initial version
!                 Parameter iq_space removed
!     16/06/2003  Switch IQ_SYM added
!
!  1. Purpose:
!
!     Write summary of GurboQuad settings to print file
!
!  2. Method
!
!     Based on the value of IQUAD a number of settings are preset
!     In the case the file [qbase].CFG exists, this file
!     is analyzed and possibly some settings are reset
!
!  3. Parameter list:
!
!Type, I/O               Name      Description
!--------------------------------------------------------------------------
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_INIT
!
!  6. Subroutines used
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
!     /S Enable subroutine tracing
!
! 10. Source code
!--------------------------------------------------------------------------------
! Local variables
!
!--------------------------------------------------------------------------------
!
call q_stack('+q_summary')
!--------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!  write summary of settings for computation of quadruplets
!  to print file
!
if (iq_prt > 0) then
  write(luq_prt,*)
  write(luq_prt,'(a)') 'Summary of settings for QUAD computation'
  write(luq_prt,'(a)') '------------------------------------------------'
  write(luq_prt,'(a,i4)')     'Number of wave numbers          :',nkq
  write(luq_prt,'(a,i4)')     'Number of directions            :',naq
  write(luq_prt,'(a,f10.5)')  'Minimum frequency (Hz)          :',fqmin
  write(luq_prt,'(a,f10.5)')  'Maximum frequency (Hz)          :',fqmax
  write(luq_prt,'(a,f10.2)')  'Water depth (m)                 :',q_depth
  write(luq_prt,'(a,i4)')     'Preferred number of locus points:',nlocus0
!
  write(luq_prt,*)
  write(luq_prt,'(a,f10.3)') 'Gravitational acceleration:',q_grav
!  write(luq_prt,'(a,f10.3)') 'Density of water          :',q_rhow
!  write(luq_prt,'(a,f10.2)') 'Power spectral tail E(f)  :',qf_tail
!  write(luq_prt,'(a,f10.2)') 'Power spectral tail N(k)  :',qk_tail
!
  write(luq_prt,*)
  if(iq_type==1) write(luq_prt,'(a)') 'IQUAD = 1: Deep water'
  if(iq_type==2) write(luq_prt,'(a)') 'IQUAD = 2: Deep water & WAM depth scaling'
  if(iq_type==3) write(luq_prt,'(a)') 'IQUAD = 3: Direct finite depth calculation'
  write(luq_prt,*)
!
  write(luq_prt,'(a,f5.2)') 'Step size in m of BQF coding:',q_dstep
  write(luq_prt,*)
!
  if(iq_grid==1) write(luq_prt,'(a)') 'Symmetric sector grid'
  if(iq_grid==2) write(luq_prt,'(a)') 'Non-symmetric sector grid'
  if(iq_grid==3) write(luq_prt,'(a)') 'Non-symmetric full circle grid'
!
  write(luq_prt,*)
  if(iq_compact==0) write(luq_prt,'(a)') 'No compacting of data along locus'
  if(iq_compact==1) write(luq_prt,'(a)') 'Compact data along locus by eliminating zero contributions'
!
  write(luq_prt,*)
  if(iq_dscale==0) write(luq_prt,'(a)') 'No WAM depth scaling'
  if(iq_dscale==1) write(luq_prt,'(a)') 'WAM depth scaling of transfer'
!
  write(luq_prt,*)
  if(iq_screen==0) write(luq_prt,'(a)') 'No output to screen'
  if(iq_screen>=1) write(luq_prt,'(a)') 'Intermediate output to screen'
  if(iq_screen>=2) write(luq_prt,'(a)') 'Intermediate output to screen + subroutine tracing'
  write(luq_prt,*)
!
  write(luq_prt,*)
  if(iq_search==0) write(luq_prt,'(a)') 'No search is carried out for nearest QUAD grid'
  if(iq_search==1) write(luq_prt,'(a)') 'A search is carried out for nearest QUAD grid'
!
  write(luq_prt,*)
  if(iq_gauleg==0) write(luq_prt,'(a)')    'Rectangular integration'
  if(iq_gauleg>0)  write(luq_prt,'(a,i4)') 'Gauss-Legendre integration with N=',iq_gauleg
!
  write(luq_prt,*)
  if(iq_cple==1) write(luq_prt,'(a)') 'Deep water coupling coefficient of Webb'
  if(iq_cple==2) write(luq_prt,'(a)') 'Finite depth coupling coefficient of H&H'
  if(iq_cple==3) write(luq_prt,'(a)') 'Finite depth coupling coefficient of Gorman'
  if(iq_cple==4) write(luq_prt,'(a)') 'Deep water coefficient of Zakharov'
  if(iq_cple==5) write(luq_prt,'(a)') 'Finite depth coefficient of Zakharov'
!
  write(luq_prt,*)
  if(iq_disp==1) write(luq_prt,'(a)') 'Deep water dispersion relation'
  if(iq_disp==2) write(luq_prt,'(a)') 'Finite depth linear dispersion relation'
  if(iq_disp==3) write(luq_prt,'(a)') 'Non linear finite depth dispersion'
!
  write(luq_prt,*)
  if(iq_filt==0) write(luq_prt,'(a)') 'Filtering of quadruplets off'
  if(iq_filt==1) then
     write(luq_prt,'(a)') 'Filtering of quadruplets on'
     write(luq_prt,*)
     write(luq_prt,'(a,f8.2)')  'Maximum ratio of k1 and k3        :',qf_krat
     write(luq_prt,'(a,f8.2)')  'Maximum directional difference    :',qf_dmax
     write(luq_prt,'(a,e12.3)') 'Fraction of maximum energy density:',qf_frac
  end if
!
!  write(luq_prt,*)
!  if(iq_geom==0) write(luq_prt,'(a)') 'Only directional scaling of loci'
!  if(iq_geom==1) write(luq_prt,'(a)') 'Geometric scaling of loci using R-T method'
!
  write(luq_prt,*)
  if(iq_locus==1) write(luq_prt,'(a)') 'Compute locus with polar method with fixed k-step'
  if(iq_locus==2) write(luq_prt,'(a)') 'Compute locus with polar method using adaptive k-step'
  if(iq_locus==3) write(luq_prt,'(a)') 'Compute locus with polar method using geometric k-step'
!
  write(luq_prt,*)
  if(iq_sym==0)  write(luq_prt,'(a)') 'Handling of symmetries disabled'
  if(iq_sym==1)  write(luq_prt,'(a)') 'Handling of symmetries enabled'
!
  write(luq_prt,*)
  if(iq_make==1) write(luq_prt,'(a)') 'Make quadruplet grid when necessary'
  if(iq_make==2) write(luq_prt,'(a)') 'Always make quadruplet grid'
  if(iq_make==3) write(luq_prt,'(a)') 'Stop after generation of quadruplet grid'
!
  write(luq_prt,*)
  if(iq_interp==1) write(luq_prt,'(a)') 'Apply bi-linear interpotion to retrieve action density'
  if(iq_interp==2) write(luq_prt,'(a)') 'Take nearest bin to retrieve action density'
!
  write(luq_prt,*)
  if(iq_lump==0) write(luq_prt,'(a)') 'Lumping of coefficients along locus disabled'
  if(iq_lump>0)  write(luq_prt,'(a)') 'Lumping of coefficients along locus enabled'
!
  write(luq_prt,*)
  if(iq_mod==0) write(luq_prt,'(a)') '?X? Spacing of point along locus as initially computed'
  if(iq_mod==1) write(luq_prt,'(a)') 'Equidistant spacing of points along locus'
!
  write(luq_prt,*)
  if(iq_trace==0) write(luq_prt,'(a)') 'Subroutine tracing disabled'
  if(iq_trace>0)  write(luq_prt,'(a)') 'Subroutine tracing enabled'
!
!
!
!
  write(luq_prt,*)
!
!  if(iq_disp==1 .and. iq_start==2) then
!    write(luqprt,'(a)') 'Start point for locus according to Resio&Tracy'
!  else
!    write(luqprt,'(a)') 'Start point for locus equal to k3'
!  end if
  write(luq_prt,*)
  write(luq_prt,'(a,i4)') 'Level of printed output        :',iq_prt
  write(luq_prt,'(a,i4)') 'Level of logging output        :',iq_log
  write(luq_prt,'(a,i4)') 'Level of test output           :',iq_test
  write(luq_prt,'(a,i4)') 'Level of trace output          :',iq_trace
  write(luq_prt,'(a,i4)') 'Level of transformation output :',iq_trf
  write(luq_prt,'(a)')   '----------------------------------------------'
end if
!
9999 continue
!
call q_stack('-q_summary')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_symmetry(k1x,k1y,k3x,k3y,k4x,k4y,symfac,nloc)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 16 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!     10/06/2003  Initial version
!     16/06/2003  Switch iq_sym added
!
!  1. Purpose:
!
!     Compute symmetry factor to reduce integration
!
!  2. Method
!
!     Compute distance between k1 and k3, and between k4 and k1
!
!  3. Parameter list:
!
! Type   i/o             Name           Description
!----------------------------------------------------------------------------------
integer, intent(in)   :: nloc         ! number of points in array with wave number
real, intent(in)      :: k1x          ! x-component  of wave number k1
real, intent(in)      :: k1y          ! y-component  of wave number k1
real, intent(in)      :: k3x          ! x-component  of wave number k3
real, intent(in)      :: k3y          ! y-component  of wave number k3
real, intent(in)      :: k4x(nloc)    ! x-components of wave number k4
real, intent(in)      :: k4y(nloc)    ! y-components of wave number k4
real, intent(out)     :: symfac(nloc) ! symmetry factor
!----------------------------------------------------------------------------------
!  4. Error messages
!
!  5. Called by:
!
!     Q_MODIFY
!
!  6. Subroutines used
!
!     Q_STACK
!
!  7. Remarks
!
!  8. structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
integer iloc      ! counter
real dk13         ! distance between k1 and k3
real dk14         ! distance between k1 and k4
!------------------------------------------------------------------------------
!
call q_stack('+q_symmetry')
!
!
! evaluate criterion |k3-k1| < |k4-k1|
! if true then symfac=1
!
symfac = 1.
if(iq_sym==1) then
  dk13 = (k1x-k3x)**2 + (k1y-k3y)**2
  do iloc=1,nloc
    dk14 = (k1x-k4x(iloc))**2 + (k1y-k4y(iloc))**2
    if (dk13 >= dk14) symfac(iloc) = 0.
  end do
end if
!
call q_stack('-q_symmetry')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_t13v4(ik1,ia1,ik3,ia3,t13,diagk1,diagk3)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 5 September 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
implicit none
!
!  0. Update history
!
!     25/02/1999  Initial version
!     14/04/1999  Extra check in GET_LOC if locus exists in internal database
!     12/10/1999  Error handling improved
!     15/01/2001  Interface extended with diagonal term
!     06/05/2002  Criterion f34_mod added to computational procedure
!     14/08/2002  Integration simplified
!     22/08/2002  Integration modified depending on actual number of non-zero points
!     26/09/2002  Boundary check for sector grid activated
!     15/04/2003  Bug fixed in handling of periodicity
!                 Nearest bin integration enabled, including diagonal term
!     25/04/2003  Output to triplet arrays for nearest bin
!     03/05/2003  Output of triplets for bi-linear interpolation enabled
!     04/06/2003  Parameter IQ_INT renamed IQ_INTEG
!     13/06/2003  Test of integration for case of nearest bin interpolation
!     25/06/2003  Bug fixed in computation of partial derivatives for contribution to
!                 diagonal term
!     27/08/2003  Short-cut when number of non-zero points on locus is ZERO
!     05/09/2003  Switches for test output in nearest bin approach modified
!
!  1. Purpose:
!
!     Compute the function T13, defined as a line integral around a locus
!
!  2. Method
!
!     See Tracy and Resio (1982) and Van Vledder (1999)
!
!  3. Parameter list:
!
! Type    I/O          Name             Description
!------------------------------------------------------------------------------
integer, intent(in) ::  ik1     !    Index of k-component of wave number k1
integer, intent(in) ::  ia1     !    Index of a-component of wave number k1
integer, intent(in) ::  ik3     !    Index of k-component of wave number k3
integer, intent(in) ::  ia3     !    Index of a-component of wave number k3
real, intent(out)   ::  t13     !    Value of line integral over a specific locus
real, intent(out)   ::  diagk1  !    Contribution to diagonal term of k1
real, intent(out)   ::  diagk3  !    Contribution to diagonal term of k3
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_XNL4V4
!
!  6. Subroutines used
!
!     Q_GETLOCUS
!     Q_PUT_BTRIPLETS
!     Q_PUT_NTRIPLETS
!
!  7. Remarks
!
!     The action density product term is given by:
!     P = n1.n2.(n3+n4)-(n1+n2).n3.n4
!
!     This term is rewritten as:
!
!     P = n1.n2.n3 + n1.n2.n4 - n1.n3.n4 - n2.n3.n4
!       = n1.n3.(n2-n4) + n2.n4.(n1-n3)

!
!  8. Structure
!
!  9. Switches
!
!     /S  enable subroutine tracing
!     /T  enable test output
!     /N  enable interpolation using nearest point
!
! 10. Source code:
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc          ! counter along locus
integer ifnd           ! indicator if correct locus is found
integer ja2,ja2p      ! direction indices for interpolation of k2
integer jk2,jk2p      ! wave number indices for interpolation of k2
integer ja4,ja4p      ! direction indices for interpolation of k4
integer jk4,jk4p      ! wave number indices for interpolation of k4
integer ikq,iaq       ! counters
!
real sumt13           ! sum along locus
real qn1,qn2,qn3,qn4  ! action densities at wave numbers k1, k2, k3 and k4
real nprod            ! wave number product
real t2,t4            ! tail factors for k2 and k4
real qd1,qd3          ! contribution to diagonal term
real rterm            ! product term along locus
!
real qn13p            ! product of N1 and N3
real qn13d            ! difference of N1 and N3
!
!
!------------------------------------------------------------------------------
call q_stack('+q_t13v4')
!
t13    = 0.
diagk1 = 0.
diagk3 = 0.
!
!
if(ik1==ik3 .and. ia1==ia3) goto 9999  ! skip routine if k1=k3
!
!  obtain information requested locus based on a information
!  about a precomputed locus, as stored in the database file
!
call q_getlocus(ik1,ia1,ik3,ia3,ifnd)
!
if(ifnd==0 .or. nlocusx==0) then
  t13 = 0.
  goto 9999
end if
!---------------------------------------------------------------------------------------
qn1 = nspec(ik1,ia1)
qn3 = nspec(ik3,ia3)
!
qn13p = qn1*qn3      ! compute product
qn13d = qn3-qn1      ! compute difference
!
sumt13 = 0
!
!    3-----------4 ja2p         w1 = (1-wk)*(1-wa)
!    |    .      |              w2 = wk*(1-wa)
!    |. . + . . .| wa2   A      w3 = (1-wk)*wa
!    |    .      |       |      w4 = wk*wa
!    |    .      |       wa
!    |    .      |       |
!    1-----------2 ja2   V
!   jk2  wk2  jk2p
!
!    <-wk->
!
!
t2 = 1.
t4 = 1.
!
!-----------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!  Main loop over the locus
!
do iloc=1,nlocusx
!
  jk2  = t_ik2(iloc)
  jk2p = min(jk2+1,nkq)
  ja2  = mod(t_ia2(iloc)-1+naq,naq)+1
  ja2p = mod(t_ia2(iloc)+naq,naq)+1
!
! compute tail parameters
!
!!  if(iq_geom==1) then
!!    jk2 = max(1,jk2)
!!    jk4 = max(1,jk4)
!!    t2 = max(1.,q_kfac**real(t_ik2(iloc)-nkq))
!!    t2 = t2**qk_tail
!!    t4 = max(1.,q_kfac**real(t_ik4(iloc)-nkq))
!!    t4 = t4**qk_tail
!!  end if
!---------------------------------------------------------------------------------------
!  check boundaries of sector grid
!
  if(iq_grid < 3) then
    ja2  = max(ja2,1)
    ja2  = min(ja2,naq)
    ja2p = max(ja2p,1)
    ja2p = min(ja2p,naq)
  end if
!
  qn2  = (t_w1k2(iloc)*nspec(jk2,ja2)  + t_w2k2(iloc)*nspec(jk2p,ja2) + &
&         t_w3k2(iloc)*nspec(jk2,ja2p) + t_w4k2(iloc)*nspec(jk2p,ja2p))*t2
!
  jk4  = t_ik4(iloc)
  jk4p = min(jk4+1,nkq)
  ja4  = mod(t_ia4(iloc)-1+naq,naq)+1
  ja4p = mod(t_ia4(iloc)+naq,naq)+1
!
!  special treatment for sector grids
!  limit range of indices
!  QQQ: in fact energy density should be set to ZERO
!
  if(iq_grid < 3) then
    ja4  = max(ja4,1)
    ja4  = min(ja4,naq)
    ja4p = max(ja4p,1)
    ja4p = min(ja4p,naq)
  end if
!
  qn4  = (t_w1k4(iloc)*nspec(jk4,ja4)  + t_w2k4(iloc)*nspec(jk4p,ja4) + &
&         t_w3k4(iloc)*nspec(jk4,ja4p) + t_w4k4(iloc)*nspec(jk4p,ja4p))*t4
!
!-------------------------------------------------------------------------------
!
  nprod     = qn13p*(qn4-qn2) + qn2*qn4*qn13d
  rterm     = t_zz(iloc)
  t13        = t13 + rterm*nprod
!
! output to triplets
!
!
!
!  add diagonal terms
!
!!  qd1    = qn3*(qn4-qn2) + qn2*qn4*qn3
!!  qd3    = qn1*(qn4-qn2) - qn2*qn4*qn1
!
  qd1    = qn3*(qn4-qn2) - qn2*qn4
  qd3    = qn1*(qn4-qn2) + qn2*qn4
  diagk1 = diagk1 + qd1*rterm
  diagk3 = diagk3 + qd3*rterm
!-----------------------------------------------------------------------------------
end do
!
!!/T if(iq_test>=4) then
!!/T   write(luq_tst,'(a)') 'Q_T13V4: NSPEC'
!!/T  do ikq=1,nkq
!!/T    write(luq_tst,'(100e12.4)') (nspec(ikq,iaq),iaq=1,naq)
!!T  end do
!!T end if
!!if(iq_integ==3) write(luq_int,'(4i3,i5,1000e13.5)') ik1,ia1,ik3,ia3,nloc, &
!!& t_s(nloc),t13,(dt13(iloc),iloc=1,nloc)
!
9999 continue
!
call q_stack('-q_t13v4')
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_weight
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 20 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
implicit none
!
!  0. Update history
!
!     13/04/1999  Initial version
!     27/10/1999  Weight computed in the case that k2m < k(1)
!     01/11/1999  Use of Q_XK and Q_SK added to compute weights if k > kmax
!     26/11/1999  Bug fixed when checking conversion
!      8/12/1999  Use of SK_MAX introduced to handle very large loci
!     09/08/2002  Modification of weights
!     13/08/2002  storage of log-spacing replace by linear spacing
!     20/08/2002  Bug fixed when geometric scaling is assumed
!
!  1. Purpose:
!
!     Compute interpolation weights of locus
!
!  2. Method
!
!     Compute position of wave number in wave number grid
!     Usable for linear interpolation
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_MAKEGRID
!
!  6. Subroutines used
!
!  7. Remarks
!
!     The tail factors wt_k2 and wt_k4 are valid for the decay of the action density spectrum
!     N(kx,ky). With p (qk_tail) the power p of the tail of the N(k) spectrum, and q the power
!     of the tail of the N(kx,ky) spectrum, we have q=p-1
!
!     Since N(k) = k N(kx,ky) with k the Jacobian
!     it follows that the tail functions are given by
!
!     k^p = k k^q   =>   k^p = k^(q+1) => p=q+1  => q=p-1
!
!  8. Structure
!
!     Initialisations
!     do for all points on locus
!       compute directional index for k2 and k4
!       if geometric scaling then
!         compute wave number index directly
!         convert log-scaling to linear scaling
!       else
!         search position of wave number in k-array
!         if k < kmin then
!           k-index = 1 and factor is 0
!         elsif k < kmax then
!           compute k-index for k2 and k4
!         else
!           compute tail factor
!         end if
!       end if
!     end do
!
!
!  9. Switches
!
!     /T  enable test output
!
! 10. Source code:
!------------------------------------------------------------------------------
!     Local variables
!
integer iloc      ! counter along locus loop
integer jpos      ! index for interpolation and tracking of position in wave numebr array
integer itest     ! local test level
real k2a,k2m      ! angle (radians) and magnitude of wave number k2
real k4a,k4m      ! angle (radians) and magnitude of wave number k2
real dk           ! difference between two wave numbers
real xtest        ! test value for checking computation of weights, by inversion test
real ff,gg        ! variables in transformation of log-spacing to linear spacing
!
! functions used
!
!!real x_kfunc      ! function to invert computation of wieghts
!------------------------------------------------------------------------------
call q_stack('+q_weight')
!
! initialisations
!
itest = iq_test        ! set local test level
itest = itest
!------------------------------------------------------------------------------
do iloc=1,nlocus
  k2m = k2m_mod(iloc)
  k2a = k2a_mod(iloc)
  k4m = k4m_mod(iloc)
  k4a = k4a_mod(iloc)
!
  wt_k2(iloc) = 0.
  wt_k4(iloc) = 0.
!
! compute directional weights
!
  wa_k2(iloc) = (k2a-q_ang1)/q_deltad+1
  wa_k4(iloc) = (k4a-q_ang1)/q_deltad+1
!------------------------------------------------------------------------------
!  compute position of k2 in wave number grid
!  and compute weight function
!-----------------------------------------------------------------------------
  if(iq_disp==1.and. iq_geom==1) then    ! deep water is assumed and loci have geometric scaling
!
    wk_k2(iloc) = 1.+alog(k2m/kqmin)/alog(q_kfac)
    wt_k2(iloc) = 1.
    wk_k4(iloc) = 1.+alog(k4m/kqmin)/alog(q_kfac)
    wt_k4(iloc) = 1.
!
!  Replace log-spacing by linear spacing
!
    ff = wk_k2(iloc)
    gg = floor(ff)
    wk_k2(iloc) = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
!
!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin2:', &
!!/T &  ff,gg,wk_k2(iloc),abs(wk_k2(iloc)-ff)/abs(ff)*100.
!
    ff = wk_k4(iloc)
    gg = floor(ff)
    wk_k4(iloc) = gg+(q_kfac**(ff-gg)-1.)/(q_kfac-1.)
!
!!/T    if(iq_test>=3) write(luq_tst,'(a,4f10.5)') 'Q_WEIGHT: wlog gg wlin4:', &
!!/T    ff,gg,wk_k4(iloc),abs(wk_k4(iloc)-ff)/abs(ff)*100.
!
!  for finite depth a search is carried out to compute
!  the position of the interacting wave number in the
!  non-geometric k-grid
!
  else
    jpos = 1
    do while (k2m > q_k(jpos))
      jpos = jpos + 1
      if(jpos > nkq) exit
    end do
!
    if(k2m <= q_k(1)) then
      wk_k2(iloc) = k2m/q_k(1)
      wt_k2(iloc) = 0.
    elseif(k2m < q_k(nkq) .and. k2m > q_k(1)) then
      dk          = q_k(jpos)-q_k(jpos-1)
      wk_k2(iloc) = real(jpos-1) + (k2m-q_k(jpos-1))/dk
      wt_k2(iloc) = 1.
    elseif(k2m >= q_k(nkq)) then
      wk_k2(iloc) = min(wk_max,real(nkq) + (k2m-q_k(nkq))/q_sk(nkq))
      wt_k2(iloc) = (k2m/q_k(nkq))**(qk_tail-1.)
!
! minus 1 to account for Jacobian from kx,ky to polar k-grid
!
    end if
!
!  compute position of k4 in wave number grid
!  and compute weight function
!
    jpos = 1
    do while (k4m > q_k(jpos))
      jpos = jpos + 1
      if(jpos > nkq) exit
    end do
!
    if(k4m <= q_k(1)) then
      wk_k4(iloc) = k4m/q_k(1)
      wt_k4(iloc) = 0.
    elseif(k4m < q_k(nkq) .and. k4m > q_k(1)) then
      dk   = q_k(jpos)-q_k(jpos-1)
      wk_k4(iloc) = real(jpos-1) + (k4m-q_k(jpos-1))/dk
      wt_k4(iloc) = 1.
    elseif(k4m >= q_k(nkq)) then
      wk_k4(iloc) = min(wk_max,real(nkq) + (k4m-q_k(nkq))/q_sk(nkq))
      wt_k4(iloc) = (k4m/q_k(nkq))**(qk_tail-1.)
    end if
!
  end if
!
!
end do
!
9999 continue
!
call q_stack('-q_weight')
!
return
end subroutine
!-----------------------------------------------------------------
subroutine q_loc_w1w3(k1x,k1y,k3x,k3y,npts,k2x,k2y,k4x,k4y,s)
!-----------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 11 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
!
implicit none
!
!  0. Update history
!
!     15/04/2002  Initial version
!     20/08/2002  Direction of k1 may be non-zero
!     27/08/2002  Singular solution crosses origin
!     11/06/2003  Length of locus fixed to 3
!
!  1. Purpose:
!
!     Compute locus for the special case w1=w3
!
!  2. Method
!
!     For this case, the k2-locus consists of a straight line
!
!  3. Parameter used:
!
integer, intent(in) :: npts      ! Number of points
real, intent(in)    :: k1x       ! x-component of wave number k1
real, intent(in)    :: k1y       ! y-component of wave number k1
real, intent(in)    :: k3x       ! x-component of wave number k3
real, intent(in)    :: k3y       ! y-component of wave number k3
!
real, intent(out)   :: k2x(npts) ! x-component of wave number k2
real, intent(out)   :: k2y(npts) ! y-component of wave number k2
real, intent(out)   :: k4x(npts) ! x-component of wave number k4
real, intent(out)   :: k4y(npts) ! y-component of wave number k4
real, intent(out)   :: s(npts)   ! distance along locus
!
!  4. Error messages
!
!  5. Caled by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!  7. Remarks
!
!     Routine based on modified version of routine SHLOCX of Resio and Tracy
!     On 15/4/2002 a bug fixed in computation of THR when angle of k3 is larger than 90°
!
!     In addition, the assumption that k1y=0 and thus dir1=0 is removed
!     In bug fix of 20/8/2002 this restriction is removed.
!
!  8. Structure
!
!     Compute angle of symmetry axis
!     Compute distance between 2 lines of solution
!     compute wave numbers along locus
!     rotate angles
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer ipt   ! counter of points along locus
!
real dirs     ! angle of symmetry axis
real dir1     ! direction of wave number k1
real dir3     ! direction of wave number k3
real dk0      ! step size along locus
real xk0      ! x-component
real yk0      ! y-component
real w2       ! radian frequency
real xx2,yy2  ! values along k2-locus
real xx4,yy4  ! values along k4-locus
real k1m      ! magnitude of wave number k1
!------------------------------------------------------------------------------
!
!    dirs is the angle of rotation from the x-axis to the "bisecting" angle
!
dir1 = atan2(k1y,k1x)
dir3 = atan2(k3y,k3x)
dirs = 0.5*(180-abs(180-abs(dir3-dir1)))
k1m  = sqrt(k1x**2 + k1y**2)
!
!     k1x is the total length of the wavenumber vector
!     xk0 is the length of this vector in the rotated coordinate system
!
xk0 = k1m * cos(dirs)
yk0 = k1m * sin(dirs)
!
! Specify step size for solution of singular case
!
!! dk0 = 0.11    ! Removed on 11/6/2003, this value is used in original WRT code
!! dk0 = kqmax/real(npts-1.)             this value depends on actual grid
dk0 = 3./real(npts-1.)                 !  this is test value
!
!  modify rotation angle
!
dirs = dirs + dir1
!
!  generate sequence of parallel lines
!  rotate lines over modified angle DIRS
!
do ipt=1,npts
!  w2       = real(ipt-1.)*dk0         ! removed on Aug. 27 2002
!
  w2       = 2.*real(ipt-npts/2)*dk0   ! create line on both sides of origin
  xx2      = w2*xk0
  yy2      = yk0
  k2x(ipt) = xx2*cos(dirs) - yy2*sin(dirs)
  k2y(ipt) = yy2*cos(dirs) + xx2*sin(dirs)
  xx4      = xx2
  yy4      = -yy2
  k4x(ipt) = xx4*cos(dirs) - yy4*sin(dirs)
  k4y(ipt) = yy4*cos(dirs) + xx4*sin(dirs)
  s(ipt)   = real(ipt-1)*dk0*xk0
end do
!
return
end subroutine
!------------------------------------------------------------------------------
subroutine q_xnl4v4(aspec,sigma,angle,nsig,nang,depth,xnl,diag,ierror)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 25 June 2003
!   +---+ |   |  Release: 5.0
!         +---+
!
! do not use m_xnldata
use m_constants
use serv_xnl4v5
implicit none
!------------------------------------------------------------------------------
!  0. Update history
!
!     08/01/2000 Initial version
!     12/01/2001 Updated interface
!     13/01/2001 Inclusion of diagonal term
!     14/02/2002 Upgrade to release 4.0, depth added to input
!     20/08/2002 quad depth adapted in the case of WAM-depth scaling
!                then deep water is assumed for conversion of A(sig,theta) -> N(kx,ky)
!                Search option for nearest grid included
!     23/08/2002 Allocation of work arrays set to fixed size
!     11/09/2002 Filtering of energy densities introduced and restructure
!     14/04/2003 Format of test write statement corrected
!     03/05/2003 Computation and output of triplets enabled
!     12/06/2003 Export spectral grid in case of Q_INTEG>1
!     16/06/2003 Switch IQ_SYM included
!                Allocation of dynamic data array's moved to Q_ALLOCATE
!     24/06/2003 Range of loop for IK3 made dependent on value of IQ_SYM
!     25/06/2003 Bug fixed in assigment of contribution of diagonal term
!
!  1. Purpose:
!
!     Compute nonlinear transfer for a given action density spectrum
!     on a given wave number and direction grid
!
!  2. Method
!
!     Compute nonlinear transfer in a surface gravity wave spectrum
!     due to resonant four wave-wave interactions
!
!     Methods: Webb/Resio/Tracy/VanVledder
!
!
!  3. Parameter list:
!
! Type    I/O          Name               Description
!------------------------------------------------------------------------------
integer,intent(in)  :: nsig             ! number of radian frequencies
integer,intent(in)  :: nang             ! number of directions
real,   intent(in)  :: aspec(nsig,nang) ! Action density spectrum as a function of (sigma,theta)
real,   intent(in)  :: sigma(nsig)      ! radian frequencies
real,   intent(in)  :: angle(nang)      ! directions in radians (sector or full circle)
real,   intent(in)  :: depth            ! water depth in m
real,   intent(out) :: xnl(nsig,nang)   ! nonlinear quadruplet interaction computed with
!                                         a certain exact method (k,theta)
real,   intent(out) :: diag(nsig,nang)  ! Diagonal term for WAM based implicit integration scheme
integer, intent(out) :: ierror          ! error indicator
!
!  4. Error messages
!
!  5. Called by:
!
!     XNL_MAIN
!
!  6. Subroutines used
!
!     Q_STACK
!     Q_INIT
!     Q_CTRGRID
!     Q_T13V4
!     Q_SEARCHGRID
!
!  7. Remarks
!
!     The external action density spectrum is given as N(sigma,dir)
!     The internal action density spectrum is given as N(kx,ky)
!
!     These 2 spectra are conected via the Jacobian transformation
!
!                cg
!     N(kx,ky) = -- N(sig,theta)
!                 k
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------------------
integer iaq      ! counter for directions
integer jaq      ! counter for directions
integer ikq      ! counter for wave numbers
integer iang     ! counter for directions
integer ia       ! counter for directions
integer ik       ! counter for wave numbers
integer idir1    ! direction in degrees of k1 (for integration test)
integer idir3    ! direction in degrees of k3 (for integration test)
real period      ! periodicity for direction, used in conversion of 2-spectra
real diagk1      ! diagonal term for k1
real diagk3      ! diagonal term for k3
!
real qn_max      ! maximum action density
real qn_min      ! minimum action density
!
real cg(nsig)         ! group velocity for conversion of spectrum and transfer
!
integer ia1,ia3,ja3   ! limits for directional loops
integer jk3           ! start of k3 loop
integer ik1,ik3       ! counters for wave number loop
integer nloc          ! number of points on locus
!
integer igrid         ! status of grid file
real t13              ! value of sub-integral
real k_rat            ! local ratio of wave numbers
real a_dif            ! directional difference
real jacobian         ! Jacobian
real qn1,qn3          ! action densities in k1 and k3
!
!  testing of diagonal term on a low level
!
real diagk1_0         ! saved value of diagk1
real diagk3_0         ! saved value of diagk3
real dq1              ! small change in action density of n1
real dq3              ! small change in action density of n3
real t13_0            ! Original estimated of diagonal term
real t13_1,t13_3      ! perturbed estimated of diagonal term
!
integer ifil_dir      ! indicator for filtering of directional criterion
integer ifil_krat     ! indicator for filtering of wave number ratio criterion
integer ifil_dens     ! indicator for filtering of action density criterion
integer ifil_tot      ! indicator for filtering due to any criterion
integer nfil_dir      ! counter to indicate filtering of directional criterion
integer nfil_krat     ! counter to indicate filtering of wave number criterion
integer nfil_dens     ! counter to indicate filtering of action density criterion
!
integer ntot_conf     ! total number of configurations
integer ntot_filt     ! total number of filtered configurations
!
!
!------------------------------------------------------------------------------
call q_stack('+q_xnl4v4')
!
! initialisations
!------------------------------------------------------------------------------
ierror = 0              ! error status
diag = 0                ! initialize output diagonal term
!
!
if(iq_type==3) then
  q_depth = depth         ! water depth to be used in computation
else
  q_depth = q_maxdepth
end if
!--------------------------------------------------------------------------
!  generate basic grid of loci and store loci in memory and to datafile
!--------------------------------------------------------------------------
if(iq_screen >= 1) write(iscreen,'(a)') 'Q_XNL4V4: Checking interaction grid '
!
if(iq_search==0 .or. iq_type/=3) then
  call q_init
  call q_ctrgrid(2,igrid)
  if(iq_err /= 0) goto 9999
!
  if(igrid/=0) then
    call q_error('e','NOGRID','No proper grid exists')
    goto 9999
  end if
!
  if(iq_make ==3) then
    call q_error('e','MAKEGRID','Only computation of grid')
    goto 9999
  end if
!------------------------------------------------------------------------------
!  set overall scale factor resulting from optional SEARCH for nearest grid
!------------------------------------------------------------------------------
!
  q_scale = 1.
!------------------------------------------------------------------------
else
!
!  search nearest valid grid and compute additional WAM scale factor
!  only active when IQ_SEARCH==1 .AND. IQ_TYPE==3
!
  call q_searchgrid(depth,igrid)


  if(igrid/=0) then
    call q_error('e','NOGRID','No proper grid exists')
    goto 9999
  end if
!
  if(iq_err /=0) goto 9999
end if
!
!------------------------------------------------------------------------------
!  convert input action density spectrum from A(sigma,theta) -> N(kx,ky)
!
do ikq=1,nkq
  call z_cmpcg(sigma(ikq),q_depth,q_grav,cg(ikq))
  do iaq=1,naq
    nspec(ikq,iaq) = aspec(ikq,iaq)/q_k(ikq)*cg(ikq)
  end do
end do
!
!------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------
!  integration over all possible configurations
!--------------------------------------------------------------------------------------
xnl = 0.
qn_max = maxval(nspec)
!
!--------------------------------------------------------------------------------------
do ik1 = 1,nkq
  if(iq_screen >= 1) write(iscreen,'(a,2i4,e12.3)') 'Q_XNL4V4: k1 nk d:',ik1,nkq,q_depth
  jk3 = ik1
  if(iq_sym==0) jk3 = 1
!
  do ia1 = iaq1,iaq2          ! loop over selected part of grid, set in q_init
!
    qn1 = nspec(ik1,ia1)
!
    do ik3 = jk3,nkq           ! compute only half-plane
      do ia3 = 1,naq           ! loop over all possible wave directions
        qn3 = nspec(ik3,ia3)
!
        if(iq_screen>=3) write(iscreen,'(a,4i4)') 'Q_XNL4V4: ik1 ia1 ik3 ia3:',ik1,ia1,ik3,ia3
!
!  computes distances in wave number space
!
        a_dif = 180. - abs(180. - abs(q_ad(ia1) - q_ad(ia3)))
        k_rat = max(q_k(ik1)/q_k(ik3), q_k(ik3)/q_k(ik1))
        qn_min = qf_frac*qn_max/(q_k(ik3)/q_k(1))**7.5
        qn_min = qf_frac*qn_max*q_kpow(ik3)
!
        ifil_dir  = 0
        ifil_krat = 0
        ifil_dens = 0
        ifil_tot  = 0
!
!  perform filtering
!
!  directional difference
!
        if(a_dif > qf_dmax) then
          ifil_dir = 1
        end if
!
!  wave number ratio
!
        if(k_rat > qf_krat) then
          ifil_krat = 1
        end if
!
!  energy density filtering
!
        if(qn1 < qn_min .and. qn3 < qn_min) then
          ifil_dens = 1
        end if
!
!
        if(ifil_dir==0 .and. ifil_krat==0 .and. ifil_dens==0 .or. iq_filt==0) then
!?        if(a_dif < qf_dmax .and. k_rat < qf_krat .or. iq_filt==0) then
!
!  perform integration along locus
!
          call q_t13v4(ik1,ia1,ik3,ia3,t13,diagk1,diagk3)
!
!
          if(iq_err /= 0) goto 9999
!
!  check contribution T13 with the computed with triplet method
!
!!/R           qt13 = 0.
!!/R           do iqtr = 1,ktriplets
!!/R             qt13 = qt13 + w_qtr(iqtr)*nspec(i_qtr(iqtr,1),i_qtr(iqtr,2))* &
!!/R &                   nspec(i_qtr(iqtr,3),i_qtr(iqtr,4))*nspec(i_qtr(iqtr,5),i_qtr(iqtr,6))
!!/R           end do
!!/R           write(iscreen,*) 'CHECK T13 QT13:',t13,qt13
!
!
!  take care of additional scale factor aring from search of nearest grid
!
          t13    = t13*q_scale
          diagk1 = diagk1*q_scale
          diagk3 = diagk3*q_scale
!
!  take care of symmetric storing of interactions
!  and factor 2 due to symmetry (if activated)
!
          if(iq_sym==1) then
            t13 = 2.*t13
            diagk1 = 2.*diagk1
            diagk3 = 2.*diagk3
          end if
!
          ja3 = ia3
          if(iq_grid==1 .and. ia3 < iaref) ja3 = naq-ia1+1
          xnl(ik1,ia1)  = xnl(ik1,ia1) + t13*q_k(ik3)*q_delta*q_dk(ik3)
          if(iq_sym==1) xnl(ik3,ja3)  = xnl(ik3,ja3) - t13*q_k(ik1)*q_delta*q_dk(ik1)
!
!  add diagonal term
!
          diag(ik1,ia1) = diag(ik1,ia1) + diagk1*q_k(ik3)*q_delta*q_dk(ik3)
          if(iq_sym==1) diag(ik3,ia3) = diag(ik3,ia3) - diagk3*q_k(ik1)*q_delta*q_dk(ik1)
!
        end if
!
!!/F         write(luq_fil,'(a,4i3,3e11.3,2f7.2,4i2)') &
!!/F &      'ik1 ia1 ik3 ia3 n1 n3 t13 adif krat fil1/2/3:', &
!!/F &       ik1,ia1,ik3,ia3,qn1,qn3,t13,a_dif,k_rat,&
!!/F &       ifil_dir,ifil_krat,ifil_dens,ifil_tot
      end do
    end do
  end do
end do
!
!
!
!
!
! write number of triplets that have been written
!
!
!------------------------------------------------------------------------------
! in the case of a symmetric sector, copy results to other part
!
! Examples: naq=5, iaref=3: 1,2,3,4,5 ->   Q(1)=Q(5)
!                                          Q(2)=Q(4)
!                                          Q(3)=Q(3)
!                                                     iaq+jaq=naq+1
!           naq=6, iaref=4: 1,2,3,4,5,6 -> Q(1)=Q(6)
!                                          Q(2)=Q(5)
!                                          Q(3)=Q(4)
!
if(iq_grid==1) then
  do ikq = 1,nkq
    do iaq=iaref,naq
      jaq = naq+1-iaq
      xnl(ikq,jaq) = xnl(ikq,iaq)
    end do
  end do
end if
!
!------------------------------------------------------------------------------
if(iq_screen>=2) write(iscreen,'(a)') 'Q_XNL4V4: Main computation ended'
!
!  Convert transfer from (kx,ky) grid to (sigma,theta) grid
!
do ikq=1,nkq
  jacobian = q_k(ikq)/cg(ikq)
  do iaq=1,naq
    xnl(ikq,iaq) = xnl(ikq,iaq)*jacobian
  end do
end do
!                                                   !
9999 continue
!
call q_stack('-q_xnl4v4')
!
return
end subroutine
!------------------------------------------------------------------------------
real function x_cosk(k)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 13 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
use serv_xnl4v5, only: z_wnumb
!
implicit none
!--------------------------------------------------------------------------------
!  0. Update history
!
!     Date        Description
!
!     13/08/2002  Initial version
!
!  1. Purpose:
!
!     Compute cosine of points on locus for given wave number k
!
!  2. Method
!
!     Explicit polar method, see Van Vledder 2000, Monterey paper
!     Optionally using a fixed k-step, geometric k-step or adaptive stepping
!
!  3. Parameters used:
!
real, intent(in) :: k  ! wave number along symmetry axis of locus
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_POLAR
!
!  6. Subroutines used:
!
!     Z_WNUMB   computation of wave number
!
!  7. Remarks
!
!     The variables q, pmag and q_depth are accessed from module m_xnldata
!     The variable q_grav is accessed from module m_constants
!
!  8. Structure
!
!  9. Switches
!
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real qq    ! constant in direct polar method qq=q/sqrt(g)
real wk    ! intemediate radian frequency
real kz    ! intermediate wave number
!------------------------------------------------------------------------------
select case(iq_disp)
!
case(1)   ! deep water
!
  qq = q/sqrt(q_grav)
  x_cosk = ((qq+sqrt(k))**4 - k**2 - pmag**2)/(2.*k*pmag)
!
case(2)   !  finite depth
!
  wk = q + x_disper(k,q_depth)
  kz = z_wnumb(wk,q_depth,q_grav)
  x_cosk = (kz**2-k**2 - pmag**2)/(2.*k*pmag)
!
end select
!
x_cosk = max(-1.,x_cosk)
x_cosk = min( 1.,x_cosk)
!
end function x_cosk
!------------------------------------------------------------------------------
real function x_cple(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,iq_cple,depth,grav)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 10 Sept. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
implicit none
!
!  0. Update history
!
!     25/02/1999  Initial version
!     25/10/1999  Names of some variables modified
!                 type and depth via interface
!     09/08/2002  Upgrade to release 4.0
!     10/09/2002  g included in interface
!
!  1. Purpose:
!
!     Compute coupling coefficient between a quadruplet of
!     interacting wave numbers
!
!  2. Method
!
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
!Type    I/O              Name      Description
!-----------------------------------------------------------------------------
real, intent(in) ::       k1x     !  x-component of wave number k1
real, intent(in) ::       k1y     !  y-component of wave number k1
real, intent(in) ::       k2x     !  x-component of wave number k2
real, intent(in) ::       k2y     !  y-component of wave number k2
real, intent(in) ::       k3x     !  x-component of wave number k3
real, intent(in) ::       k3y     !  y-component of wave number k3
real, intent(in) ::       k4x     !  x-component of wave number k4
real, intent(in) ::       k4y     !  y-component of wave number k4
integer, intent(in) ::    iq_cple !  Type of coupling coefficient
real, intent(in) ::       depth   !  Water depth in meters
real, intent(in) ::       grav    !  Gravitational acceleration
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used
!
!     X_WEBB
!     X_HH
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!-------------------------------------------------------------------------------
!     Local variables
!                        ! real functions to compute coupling coefficient
!!real xc_webb             ! Webb, deep water
!!real xc_hh               ! Herterich and Hasselmann, finite depth
!------------------------------------------------------------------------------
if (iq_cple < 1 .or. iq_cple > 4) then
  x_cple = 0.
  goto 9999
end if
!
select case(iq_cple)
!
!  1) Deep water coupling coefficient of Webb
!
case(1)
  x_cple = xc_webb(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,grav)
!
!  2) finite depth coupling coefficient of Herterich and Hasselmann
!     as implemented in the Resio-Tracy program SB5
!
case(2)
  x_cple = xc_hh(k4x,k4y,k3x,k3y,k2x,k2y,k1x,k1y,depth)
!
! x_cple = xc_hh2(k1x,k1y,k2x,k2y,k3x,k3y,depth,grav)
!
end select
!
9999 continue
!
return
end function
!------------------------------------------------------------------------------
real function x_flocus(kxx,kyy)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
use m_constants
implicit none
!
!  0. Update history
!
!     25/02/1999 Initial version
!     20/07/1999 Bug fixed when IDISP=2
!     09/10/1999 Values of w2 and w4 in double precision
!                to improve accuracy of computation of z
!     09/08/2002 Upgrade to release 4.0
!
!  1. Purpose:
!
!     Compute locus function used for the determination of the
!     resonance condition
!
!  2. Method
!
!     Explicit function evaluation
!
!  3. Parameter list:
!
!Type    I/O         Name       Description
!-----------------------------------------------------------
real, intent(in) ::  kxx      !  x-component of wave number
real, intent(in) ::  kyy      !  y-component of wave number
!
!  4. Error messages
!
!
!  5. Called by:
!
!     Q_LOCPOS
!
!  6. Subroutines used
!
!     X_DISPER
!
!  7. Remarks
!
!     if iq_disp not valid, then q_disper = -1
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!-------------------------------------------------------------------------------
!     Local variables
!
real z               ! diferrence
real k2m,k4m         ! wave number magnitudes
real (kind=8) w2,w4  ! radian frequencies
!!real x_disper
!------------------------------------------------------------------------------
!call q_stack('+x_flocus')
!
select case(iq_disp)
  case (1)
    w2 = sqrtg * (kxx**2 + kyy**2)**(0.25)
    w4 = sqrtg * ((kxx+px)**2 + (kyy+py)**2)**(0.25)
    z = q + w2 - w4
!
  case (2)
    k2m = sqrt(kxx**2+kyy**2)
    k4m = sqrt((kxx+px)**2 + (kyy+py)**2)
    w2 = x_disper(k2m,q_depth)
    w4 = x_disper(k4m,q_depth)
    z  = q + w2 - w4
!
  case default
    z = -1
  end select
!
x_flocus = z
!
!call q_stack('-x_flocus')
!
return
end function
!------------------------------------------------------------------------------
real function x_jacobian(x2,y2,x4,y4)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
!% use serv_xnl4v5, only: z_cmpcg
!
implicit none
!
!  0. Update history
!
!     25/02/1999 Initial version
!     12/10/1999 Overflow avoided by checking for argument k*d
!     29/10/1999 Bug fixed in computing finite depth gradient
!     27/12/1999 Factor SQRT(Grav) added
!     01/10/2001 Components of k4 wave number explicitly input
!                New version of function X_GRAD
!     09/08/2002 Computation of Jacobian replace by |cg2-cg4|
!                based on old routine X_GRAD2
!                Upgrade to release 4.0
!
!  1. Purpose:
!
!     Compute gradient/Jacobian term for a given point on the locus
!
!  2. Method
!
!     Explicit expressions for gradient term
!     Using expression of Rasmussen (1998)
!     J = |cg2-cg4|
!
!  3. Parameter list:
!
! Type    I/O        Name    Description
!--------------------------------------------------------------------
real, intent(in) ::  x2   !  x-component of wave number k2
real, intent(in) ::  y2   !  y-component of wave number k2
real, intent(in) ::  x4   !  x-component of wave number k4
real, intent(in) ::  y4   !  y-component of wave number k4
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_CMPLOCUS
!
!  6. Subroutines used:
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!------------------------------------------------------------------------------
! local variables
!
real k2m,k4m    ! wave number magnitudes
real k2md,k4md  ! k*d values
real ang2,ang4  ! directions
real cg2,cg4    ! group velocities
real sig2,sig4  ! radian frequencies
!------------------------------------------------------------------------------
k2m = sqrt(x2**2 + y2**2)
k4m = sqrt(x4**2 + y4**2)
!
ang2 = atan2(x2,y2)
ang4 = atan2(x4,y4)
!
sig2 = sqrt(q_grav*k2m*tanh(k2m*q_depth))
sig4 = sqrt(q_grav*k4m*tanh(k4m*q_depth))
!
k2md = k2m*q_depth
k4md = k4m*q_depth
!
if(k2md > 20) then
  cg2 = 0.5*q_grav/sig2
else
  cg2 = sig2/k2m*(0.5+k2md/sinh(2*k2md))
end if
!
if(k4md > 20) then
  cg4 = 0.5*q_grav/sig4
else
  cg4 = sig4/k4m*(0.5+k4md/sinh(2*k4md))
end if
!
x_jacobian = sqrt(cg2**2+cg4**2-2*cg2*cg4*cos(ang2-ang4))
!
return
end function
!------------------------------------------------------------------------------
real function x_disper(k,d)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
implicit none
!
!  0. Update history
!
!     16/02/1999  Initial version
!     25/02/1999  Short cut if kd > 10
!     09/08/2002  Upgrade to release 4.0
!
!  1. Purpose:
!
!     Compute radian frequency for a given wave number and water depth
!
!  2. Method
!
!     Depending on the value of the parameter iq_disp the radian
!     wave number is computed as:
!     1) deep water
!     2) finite depth linear dispersion relation
!     3) finited depth non-linear dispersion relation (NOT YET implemented)
!
!  3. Parameter list:
!
! Type    I/O          Name       Description
!----------------------------------------------------------------------
real, intent(in)   ::   k   !   wave number
real, intent(in)   ::   d   !   water depth in m
!
!  4. Error messages
!
!     if iq_type not valid, then q_disper = -1
!
!  5. Called by:
!
!     Q_CHKRES
!
!  6. Subroutines used
!
!
!  7. Remarks
!
!     Type of dispersion relation is determined by the parameter IQ_DISP:
!
!     IQ_DISP==1  deep water linear disperion relation is used
!              2  finite depth linear dispersion relation is used
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
integer id   ! copy of iq_type
real kd      ! k*d
!
kd = k * d
id = iq_disp
!
if (kd > 20.) id = 1
!
select case(id)
  case (1)                           ! deep water w^2=g k
    x_disper = sqrt(q_grav*k)
  case (2)                           ! finite depth w^2 = g k tanh(k d)
    x_disper = sqrt(q_grav*k*tanh(k*d))
  case default
    x_disper = -1.
end select
!
return
end function
!------------------------------------------------------------------------------
real function x_locus1(k2)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 AUg. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
use m_constants
implicit none
!
!  0. Update history
!
!     Date        Description
!
!     23/11/1999  Initial version
!      9/08/2002  Upgrade to release 4.0
!
!  1. Purpose:
!
!     Compute locus function along symmetry axis
!
!  2. Method
!
!     See ALKYON, 1999
!
!  3. Parameter list:
!
!Type   I/O         name    Description
!-------------------------------------------------------
real, intent(in) :: k2   !  Magnitude of wave number k2
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_LOCPOS
!
!  6. Subroutines used
!
!     x_disper
!
!  7. Remarks
!
!     The routine assumes that w1 < w3 or q<0
!     implying that the directions of k2 and P are opposite
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code
!------------------------------------------------------------------------------
!     Local variables
!
real k4       ! wave number magnitudes of k4
real w2,w4    ! radian frequencies of wave numbers k2 and k4
real z        ! function value
!
!!real x_disper
!
select case(iq_disp)
  case (1)
    w2 = sqrtg * sqrt(k2)
    w4 = sqrtg * sqrt(abs(-pmag+k2))
    z  = q + w2 - w4
!
  case (2)
    k4 = abs(-pmag+k2)
    w2 = x_disper(k2,q_depth)
    w4 = x_disper(k4,q_depth)
    z  = q + w2 - w4
!
  case default
    z = -1
  end select
!
x_locus1 = z
!
return
end function
!------------------------------------------------------------------------------
real function x_locus2(lambda)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 9 Aug. 2002
!   +---+ |   |  Release: 4.0
!         +---+
!
! do not use m_xnldata
use m_constants
implicit none
!
!  0. Update history
!
!     Date        Description
!
!     23/11/1999  Initial version
!     09/08/2002  Upgrade to release 4.0
!
!  1. Purpose:
!
!     Compute locus function perpendicluar to symmetry axis
!
!  2. Method
!
!     See ALKYON, 1999
!
!  3. Parameter list:
!
!     Name    I/O  Type  Description
!
real, intent(in) ::  lambda
!
!  4. Error messages
!
!  5. Called by:
!
!     Q_LOCPOS
!
!  6. Subroutines used
!
!     x_disper
!
!  7. Remarks
!
!     The routine assumes that w1 < w3 or q<0
!     implying that the directions of k2 and P are opposite
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!------------------------------------------------------------------------------
!     local variables
!
real kk2x,kk2y,kk2m  ! wave number components and magnitude for k2
real kk4x,kk4y,kk4m  ! wave number components and magnitude for k4
real w2,w4           ! radian frequencies of wave numbers k2 and k4
real z               ! function value
!!real x_disper
kk2x = kmidx - lambda*py
kk2y = kmidy + lambda*px
kk2m = sqrt(kk2x**2 + kk2y**2)
!
kk4x = kk2x + px
kk4y = kk2y + py
kk4m = sqrt(kk4x**2 + kk4y**2)
!
select case(iq_disp)
  case (1)
    w2 = sqrtg * sqrt(kk2m)
    w4 = sqrtg * sqrt(kk4m)
    z  = q + w2 - w4
!
  case (2)
!
    w2 = x_disper(kk2m,q_depth)
    w4 = x_disper(kk4m,q_depth)
    z  = q + w2 - w4
!
  case default
    z = -1
  end select
!
x_locus2 = z
!
return
end function
!------------------------------------------------------------------------------
real function xc_hh(w1x0,w1y0,w2x0,w2y0,w3x0,w3y0,z4x,z4y,h)
!------------------------------------------------------------------------------
!
!  factor EPS included
!
implicit none
!
real z4x,z4y  ! dummy arguments
!
real w1x0,w1y0,w2x0,w2y0,w3x0,w3y0,h,dsq
real om1,om2,om3,om4,scpl1,scpl2,scpl3,stot
real t1,t2,t3,t4,t5,tot1,tot2,tot3,tot4,tot5
real som1,som2,som3
real s1,s2,s3,z1,z2,z3,z4,z5
real p1,p2,p3,p4,di,tnz1,tnz2,tnz3,tnz23
real csz1,csz2,csz3,csz23
real e,g,gsq,omsq23,pi4
real dot123,dot23
!!real cosz,tanz
!
real eps
!

!     calculates coupling coefficient in shallow water given k1,k2,k3
      real k1,k2,k3,k1x,k2x,k3x,k1y,k2y,k3y,k23x,k23y,k23,k1x0,k1y0,   &
     &     k2x0,k2y0,k3x0,k3y0,k1zx,k1zy
      data pi4/0.785398163/
!
eps = 1.e-20
g = 9.81
!
z4x = z4x
z4y = z4y
!
!      print *,'entering cplesh depth = ',h
      k1x0=w1x0
      k1y0=w1y0
      k2x0=w2x0
      k2y0=w2y0
      k3x0=w3x0
      k3y0=w3y0

      tot1=0.
      tot2=0.
      tot3=0.
      tot4=0.
      tot5=0.
      z1=0.
      z2=0.
      z3=0.
      z4=0.
      z5=0.
      g=9.81
      gsq=g*g

      s1=1.
      s2=1.
      s3=-1.


      k1x=s1*k1x0
      k1y=s1*k1y0
      k2x=s2*k2x0
      k2y=s2*k2y0
      k3x=s3*k3x0
      k3y=s3*k3y0

      k1=sqrt(k1x**2+k1y**2)
      k2=sqrt(k2x**2+k2y**2)
      k3=sqrt(k3x**2+k3y**2)

       tnz1=tanz(k1*h)
       tnz2=tanz(k2*h)
       tnz3=tanz(k3*h)
       csz1=cosz(k1*h)
       csz2=cosz(k2*h)
       csz3=cosz(k3*h)
      om1=sqrt(g*k1*tnz1)
      om2=sqrt(g*k2*tnz2)
      om3=sqrt(g*k3*tnz3)

      som1=s1*om1
      som2=s2*om2
      som3=s3*om3
      dot23=k2x*k3x+k2y*k3y

      k23x=k2x+k3x
      k23y=k2y+k3y
      k23=sqrt(k23x**2+k23y**2)
      tnz23=tanz(k23*h)
      csz23=cosz(k23*h)

      omsq23=g*k23*tnz23
      dot123=k1x*k23x+k1y*k23y

!     note the "i**2" factor is included in this term
      DI=-(som2+som3)*(k2*k3*tnz2*tnz3-dot23)  &
     & +0.5*(som2*k3**2/(csz3)**2+som3*k2**2/(csz2)**2)

      E=0.5/g *(dot23-som2*som3/gsq*(om2**2+om3**2+som2*som3))

      p1=2.*(som1+som2+som3)*(om1**2.*omsq23/gsq-dot123)

      p2=-som1*(k23)**2/(csz23)**2

      p3=-(som2+som3)*k1**2/(csz1)**2

      z1=z1+di
      z2=z2+omsq23-(som2+som3)**2
      z3=z3+p1
      z4=z4+p2
      z5=z5+p3
      T1=DI/(omsq23-(som2+som3)**2 + eps ) * (p1+p2+p3)

      T2=-DI*som1/gsq *(om1**2+omsq23)

      p4=g*k1**2/(csz1)**2

      T3=E*(som1**3*(som2+som3)/g - g*dot123 - p4)

      T4=0.5*som1/gsq*dot23*((som1+som2+som3)*(om2**2+om3**2)  &
     &                           +som2*som3*(som2+som3))

      T5=-0.5*som1*om2**2*k3**2/gsq*(som1+som2+2.*som3)  &
     &   -0.5*som1*om3**2*k2**2/gsq*(som1+2.*som2+som3)

      scpl1=T1+T2+T3+T4+T5
      tot1=tot1+t1
      tot2=tot2+t2
      tot3=tot3+t3
      tot4=tot4+t4
      tot5=tot5+t5

      s1=1.
      s2=-1.
      s3=1.
      k1zx=k1x0
      k1zy=k1y0
      k1x0=k2x0
      k1y0=k2y0
      k2x0=k3x0
      k2y0=k3y0
      k3x0=k1zx
      k3y0=k1zy


      k1x=s1*k1x0
      k1y=s1*k1y0
      k2x=s2*k2x0
      k2y=s2*k2y0
      k3x=s3*k3x0
      k3y=s3*k3y0

      k1=sqrt(k1x**2+k1y**2)
      k2=sqrt(k2x**2+k2y**2)
      k3=sqrt(k3x**2+k3y**2)
      tnz1=tanz(k1*h)
      tnz2=tanz(k2*h)
      tnz3=tanz(k3*h)
       csz1=cosz(k1*h)
       csz2=cosz(k2*h)
       csz3=cosz(k3*h)
      om1=sqrt(g*k1*tnz1)
      om2=sqrt(g*k2*tnz2)
      om3=sqrt(g*k3*tnz3)
      som1=s1*om1
      som2=s2*om2
      som3=s3*om3
      dot23=k2x*k3x+k2y*k3y
      k23x=k2x+k3x
      k23y=k2y+k3y
      k23=sqrt(k23x**2+k23y**2)
      tnz23=tanz(k23*h)
      csz23=cosz(k23*h)
      omsq23=g*k23*tnz23
      dot123=k1x*k23x+k1y*k23y

!     note the "i**2" factor is included in this term
      DI=-(som2+som3)*(k2*k3*tnz2*tnz3-dot23)  &
     & +0.5*(som2*k3**2/(csz3)**2+som3*k2**2/(csz2)**2)

      E=0.5/g *(dot23-som2*som3/gsq *(om2**2+om3**2+som2*som3))

      p1=2.*(som1+som2+som3)*(om1**2.*omsq23/gsq-dot123)

      p2=-som1*(k23)**2/(csz23)**2

      p3=-(som2+som3)*k1**2/(csz1)**2
      z1=z1+di
      z2=z2+omsq23-(som2+som3)**2
      z3=z3+p1
      z4=z4+p2
      z5=z5+p3

      T1=DI/(omsq23-(som2+som3)**2) * (p1+p2+p3)

      T2=-DI*som1/gsq *(om1**2+omsq23)

      p4=g*k1**2/(csz1)**2

      T3=E*(som1**3*(som2+som3)/g - g*dot123 - p4)

      T4=0.5*som1/gsq*dot23*((som1+som2+som3)*(om2**2+om3**2) &
     &                           +som2*som3*(som2+som3))

      T5=-0.5*som1*om2**2*k3**2/gsq*(som1+som2+2.*som3)  &
     &   -0.5*som1*om3**2*k2**2/gsq*(som1+2.*som2+som3)

      scpl2=T1+T2+T3+T4+T5
      tot1=tot1+t1
      tot2=tot2+t2
      tot3=tot3+t3
      tot4=tot4+t4
      tot5=tot5+t5

      s1=-1.
      s2=1.
      s3=1.
      k1zx=k1x0
      k1zy=k1y0
      k1x0=k2x0
      k1y0=k2y0
      k2x0=k3x0
      k2y0=k3y0
      k3x0=k1zx
      k3y0=k1zy


      k1x=s1*k1x0
      k1y=s1*k1y0
      k2x=s2*k2x0
      k2y=s2*k2y0
      k3x=s3*k3x0
      k3y=s3*k3y0

      k1=sqrt(k1x**2+k1y**2)
      k2=sqrt(k2x**2+k2y**2)
      k3=sqrt(k3x**2+k3y**2)
      tnz1=tanz(k1*h)
      tnz2=tanz(k2*h)
      tnz3=tanz(k3*h)
      csz1=cosz(k1*h)
      csz2=cosz(k2*h)
      csz3=cosz(k3*h)
      om1=sqrt(g*k1*tnz1)
      om2=sqrt(g*k2*tnz2)
      om3=sqrt(g*k3*tnz3)
      som1=s1*om1
      som2=s2*om2
      som3=s3*om3
      dot23=k2x*k3x+k2y*k3y
      k23x=k2x+k3x
      k23y=k2y+k3y
      k23=sqrt(k23x**2+k23y**2)
      tnz23=tanz(k23*h)
       csz23=cosz(k23*h)
      omsq23=g*k23*tnz23
      dot123=k1x*k23x+k1y*k23y

!     note the "i**2" factor is included in this term
      DI=-(som2+som3)*(k2*k3*tnz2*tnz3-dot23) &
     & +0.5*(som2*k3**2/(csz3)**2+som3*k2**2/(csz2)**2)

      E=0.5/g *(dot23-som2*som3/gsq *(om2**2+om3**2+som2*som3))

      p1=2.*(som1+som2+som3)*(om1**2.*omsq23/gsq-dot123)

      p2=-som1*(k23)**2/(csz23)**2

      p3=-(som2+som3)*k1**2/(csz1)**2
      z1=z1+di
      z2=z2+omsq23-(som2+som3)**2
      z3=z3+p1
      z4=z4+p2
      z5=z5+p3

      T1=DI/(omsq23-(som2+som3)**2) * (p1+p2+p3)

      T2=-DI*som1/gsq*(om1**2+omsq23)

      p4=g*k1**2/(cosz(k1*h))**2

      T3=E*(som1**3*(som2+som3)/g - g*dot123 - p4)

      T4=0.5*som1/gsq*dot23*((som1+som2+som3)*(om2**2+om3**2) &
     &                           +som2*som3*(som2+som3))

      T5=-0.5*som1*om2**2*k3**2/gsq*(som1+som2+2.*som3) &
     &   -0.5*som1*om3**2*k2**2/gsq*(som1+2.*som2+som3)

      scpl3=T1+T2+T3+T4+T5
      tot1=tot1+t1
      tot2=tot2+t2
      tot3=tot3+t3
      tot4=tot4+t4
      tot5=tot5+t5

      stot=(scpl1+scpl2+scpl3)
      om4=om2+om3-om1
      dsq=stot*stot*pi4*gsq/(om1*om2*om3*om4+eps) ! eps by GVV
      xc_hh = dsq
!
!  possible bug fixed
!
      xc_hh = xc_hh*gsq
!
      RETURN
      end  function

      real function tanz(x)
      real x
!      print *,'inside tanz '
      if (x.gt.20.) x=25.
      tanz=tanh(x)
!      print *,'after def of tanz'
      return
      end  function

      real function cosz(x)
      real x
      if (x.gt.20.) x=25.
      cosz=cosh(x)
      return
      end  function


!------------------------------------------------------------------------------
real function xc_webb(k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,grav)
!------------------------------------------------------------------------------
!
!   +-------+    ALKYON Hydraulic Consultancy & Research
!   |       |    Gerbrant van Vledder
!   |   +---+
!   |   | +---+  Last update: 10  Sep. 2002
!   +---+ |   |  Release: 5.0
!         +---+
!
implicit none
!
!  0. Update history
!
!     25/02/1999 Initial version
!     10/09/2002 Upgrade of documention and interface
!
!  1. Purpose:
!
!     Compute deep water coupling coefficient for
!     non-linear quadruplet interactions
!
!
!  2. Method
!
!     Webb (1978) and modified and corrected by Dungey and Hui (1979)
!
!  3. Parameter list:
!
! Type    I/O        Name    Description
!--------------------------------------------------------------------
real, intent(in) ::  k1x   !  x-component of wave number k1
real, intent(in) ::  k1y   !  y-component of wave number k1
real, intent(in) ::  k2x   !  x-component of wave number k2
real, intent(in) ::  k2y   !  y-component of wave number k2
real, intent(in) ::  k3x   !  x-component of wave number k3
real, intent(in) ::  k3y   !  y-component of wave number k3
real, intent(in) ::  k4x   !  x-component of wave number k4
real, intent(in) ::  k4y   !  y-component of wave number k4
real, intent(in) ::  grav  !  gravitational acceleration m/s^2
!
!  4. Error messages
!
!  5. Called by:
!
!     X_CPLE
!
!  6. Subroutines used:
!
!  7. Remarks
!
!  8. Structure
!
!  9. Switches
!
! 10. Source code:
!------------------------------------------------------------------------------
! local variables
!
!
double precision wsqp12         ! derived variable
double precision wsqm13         ! derived variable
double precision wsq13          ! derived variable
double precision wsqm14         ! derived variable
double precision wsq14          ! derived variable
double precision wsq12          ! derived variable
real z,z12,z13,z14              ! derived variables
real dwebb                      ! final coefficient
real p1,p2,p3,p4,p5,p6,p7,p8,p9 ! partial summations
real w1,w2,w3,w4                ! radian frequencies
real k1,k2,k3,k4                ! wave number magnitudes
real dot12                      ! k1*k2
real dot13                      ! k1*k3
real dot14                      ! k1*k4
real dot23                      ! k2*k3
real dot24                      ! k2*k4
real dot34                      ! k3*k4
real pi                         ! pi
real pi4                        ! pi/4
real eps                        ! internal accuracy
!---------------------------------------------------------------------
! initialisations
!---------------------------------------------------------------------
pi  = 4.*atan(1.)
pi4 = 0.25*pi
!
eps = 1.0e-30
!
k1 = sqrt(k1x*k1x + k1y*k1y)
k2 = sqrt(k2x*k2x + k2y*k2y)
k3 = sqrt(k3x*k3x + k3y*k3y)
k4 = sqrt(k4x*k4x + k4y*k4y)
!
w1 = sqrt(k1)
w2 = sqrt(k2)
w3 = sqrt(k3)
w4 = sqrt(k4)
!
dot12 = k1x*k2x + k1y*k2y
dot13 = k1x*k3x + k1y*k3y
dot14 = k1x*k4x + k1y*k4y
dot23 = k2x*k3x + k2y*k3y
dot24 = k2x*k4x + k2y*k4y
dot34 = k3x*k4x + k3y*k4y
!
wsqp12= sqrt((k1x+k2x)*(k1x+k2x)+(k1y+k2y)*(k1y+k2y))
wsq12 = (w1+w2)*(w1+w2)
wsqm13= sqrt((k1x-k3x)*(k1x-k3x)+(k1y-k3y)*(k1y-k3y))
wsq13 = (w1-w3)*(w1-w3)
wsqm14= sqrt((k1x-k4x)*(k1x-k4x)+(k1y-k4y)*(k1y-k4y))
wsq14 = (w1-w4)*(w1-w4)
z12   = wsqp12-wsq12
z13   = wsqm13-wsq13
z14   = wsqm14-wsq14
z     = 2.*wsq12*(k1*k2-dot12)*(k3*k4-dot34)
p1    = z/(z12+eps)
z     = 2.*wsq13*(k1*k3+dot13)*(k2*k4+dot24)
p2    = z/(z13+eps)
z     = 2.*wsq14*(k1*k4+dot14)*(k2*k3+dot23)
p3    = z/(z14+eps)
p4    =  0.5 *(dot12*dot34 + dot13*dot24 + dot14*dot23)
p5    =  0.25*(dot13+dot24) * wsq13 * wsq13
p6    = -0.25*(dot12+dot34) * wsq12 * wsq12
p7    =  0.25*(dot14+dot23) * wsq14 * wsq14
p8    =  2.5*k1*k2*k3*k4
p9    = wsq12*wsq13*wsq14* (k1 + k2 + k3 + k4)
!
dwebb  = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
xc_webb = grav**2*pi4*dwebb*dwebb/(w1*w2*w3*w4+eps)
!
return
end  function
!
end module
