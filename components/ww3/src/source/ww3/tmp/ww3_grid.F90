#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      PROGRAM W3GRID
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |           J. H. Alves             |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2014 |
!/                  +-----------------------------------+
!/
!/    14-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    27-Jan-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/                  Add UNFORMATTED bath file option.
!/                  Read options with namelists.
!/    14-Feb-2000 : Adding exact Snl                    ( version 2.01 )
!/    04-May-2000 : Non central source term int.        ( version 2.03 )
!/    24-Jan-2001 : Flat grid option.                   ( version 2.06 )
!/    02-Feb-2001 : Xnl version 3.0                     ( version 2.07 )
!/    09-Feb-2001 : Third propagation scheme added.     ( version 2.08 )
!/    27-Feb-2001 : O0 output switch added.             ( version 2.08 )
!/    16-Mar-2001 : Fourth propagation scheme added.    ( version 2.09 )
!/    29-Mar-2001 : Sub-grid island treatment.          ( version 2.10 )
!/    20-Jul-2001 : Clean up.                           ( version 2.11 )
!/    12-Sep-2001 : Clean up.                           ( version 2.13 )
!/    09-Nov-2001 : Clean up.                           ( version 2.14 )
!/    11-Jan-2002 : Sub-grid ice treatment.             ( version 2.15 )
!/    17-Jan-2002 : DSII bug fix.                       ( version 2.16 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    26-Nov-2002 : Adding first version of NL-3/4.     ( version 3.01 )
!/                  Removed before distribution in 3.12.
!/    26-Dec-2002 : Relaxing CFL time step.             ( version 3.02 )
!/    01-Aug-2003 : Modify GSE correction for moving gr.( version 3.03 )
!/                  Add offset option for first direction.
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    04-May-2005 : Allow active points at edge.        ( version 3.07 )
!/    07-Jul-2005 : Add MAPST2 and map processing.      ( version 3.07 )
!/    09-Nov-2005 : Remove soft boundary options.       ( version 3.08 )
!/    23-Jun-2006 : Adding alternative source terms.    ( version 3.09 )
!/                  Module W3SLN1MD, dummy for others.
!/    28-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    28-Oct-2006 : Spectral partitioning.              ( version 3.09 )
!/    09-Jan-2007 : Correct edges of read mask.         ( version 3.10 )
!/    26-Mar-2007 : Add to spectral partitioning.       ( version 3.11 )
!/    14-Apr-2007 : Add Miche style limiter.            ( version 3.11 )
!/                  ( J. H. Alves )
!/    25-Apr-2007 : Battjes-Janssen Sdb added.          ( version 3.11 )
!/                  ( J. H. Alves )
!/    18-Sep-2007 : Adding WAM4 physics option.         ( version 3.13 )
!/                  ( F. Ardhuin )
!/    09-Oct-2007 : Adding bottom scattering SBS1.      ( version 3.13 )
!/                  ( F. Ardhuin )
!/    22-Feb-2008 : Initialize TRNX-Y properly.         ( version 3.13 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    23-Jul-2009 : Modification of ST3 namelist  .     ( version 3.14-SHOM )
!/    31-Mar-2010 : Addition of shoreline reflection    ( version 3.14-IFREMER )
!/    29-Jun-2010 : Adding Stokes drift profile output  ( version 3.14-IFREMER )
!/    30-Aug-2010 : Adding ST4 option                   ( version 3.14-IFREMER )
 
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    29-Oct-2010 : Clean up of unstructured grids      ( version 3.14.4 )
!/                  (A. Roland and F. Ardhuin)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid. Change GLOBAL
!/                  input in ww3_grid.inp to CSTRG.     ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    25-Jun-2011 : Adding movable bed friction         ( version 4.01 )
!/    16-Sep-2011 : Clean up.                           ( version 4.05 )
!/    01-Dec-2011 : New namelist for reflection         ( version 4.05 )
!/    01-Mar-2012 : Bug correction for NLPROP in ST2    ( version 4.05 )
!/    12-Jun-2012 : Add /RTD option or rotated grid.    ( version 4.06 )
!/                  (Jian-Guo Li, UK Met Office)
!/    13-Jul-2012 : Move data structures GMD (SNL3) and nonlinear
!/                  filter (SNLS) from 3.15 (HLT).      ( version 4.07 )
!/    02-Sep-2012 : Clean up of reflection and UG grids ( version 4.08 )
!/    12-Dec-2012 : Adding SMC grid.  JG_Li             ( version 4.08 )
!/    19-Dec-2012 : Add NOSWLL as namelist variable.    ( version 4.OF )
!/    05-Mar-2013 : Adjusted default roughness for rocks( version 4.09 )
!/    01-Jun-2013 : Adding namelist for spectral output ( version 4.10 )
!/    12-Sep-2013 : Adding Arctic part for SMC grid.    ( version 4.11 )
!/    01-Nov-2013 : Changed UG list name to UNST        ( version 4.12 )
!/    11-Nov-2013 : SMC and rotated grid incorporated in the main
!/                  trunk                               ( version 4.13 )
!/    13-Nov-2013 : Moved out reflection to W3UPDTMD    ( version 4.12 )
!/    27-Jul-2013 : Adding free infragravity waves      ( version 4.15 )
!/    02-Dec-2013 : Update of ST4                       ( version 4.16 )
!/    16-Feb-2014 : Adds wind bias correction: WCOR     ( version 5.00 )
!/    10-Mar-2014 : Adding namelist for IC2             ( version 5.01 )
!/    29-May-2014 : Adding namelist for IC3             ( version 5.01 )
!/    15 Oct-2015 : Change SMC grid input files. JGLi   ( version 5.09 )
!/
!/    Copyright 2009-2013 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     "Grid" preprocessing program, which writes a model definition
!     file containing the model parameter settigs and grid data.
!
!  2. Method :
!
!     Information is read from the file ww3_grid.inp (NDSI), or
!     preset in this program. A model definition file mod_def.ww3 is
!     then produced by W3IOGR. Note that the name of the model
!     definition file is set in W3IOGR.
!
!  3. Parameters :
!
!     Local parameters.
!     ----------------------------------------------------------------
!       NDSI    Int.  Input unit number ("ww3_grid.inp").
!       NDSS    Int.  Scratch file.
!       NDSG    Int.  Grid unit ( may be NDSI )
!       NDSTR   Int.  Sub-grid unit ( may be NDSI or NDSG )
!       VSC     Real  Scale factor.
!       VOF     Real  Add offset.
!       ZLIM    Real  Limiting bottom depth, used to define land.
!       IDLA    Int.  Layout indicator used by INA2R.
!       IDFM    Int.  Id. FORMAT indicator.
!       RFORM   C*16  Id. FORMAT.
!       FNAME   C*60  File name with bottom level data.
!       FROM    C*4   Test string for open, 'UNIT' or 'FILE'
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. W3GDATMD Set number of model.
!      W3SETG    Subr.   Id.    Point to selected model.
!      W3DIMS    Subr.   Id.    Set array dims for a spectral grid.
!      W3DIMX    Subr.   Id.    Set array dims for a spatial grid.
!      W3GRMP    Subr. W3GSRUMD Compute bilinear interpolation for point
!      W3NOUT    Subr. W3ODATMD Set number of model for output.
!      W3SETO    Subr.   Id.    Point to selected model for output.
!      W3DMO5    Subr.   Id.    Set array dims for output type 5.
!      ITRACE    Subr. W3SERVMD Subroutine tracing initialization.
!      STRACE    Subr.   Id.    Subroutine tracing.
!      NEXTLN    Subr.   Id.    Get next line from input filw
!      EXTCDE    Subr.   Id.    Abort program as graceful as possible.
!      DISTAB    Subr. W3DISPMD Make tables for solution of the
!                               dispersion relation.
!      READNL    Subr. Internal Read namelist.
!      INAR2R    Subr. W3ARRYMD Read in an REAL array.
!      PRTBLK    Subr.   Id.    Print plot of array.
!      W3IOGR    Subr. W3IOGRMD Reading/writing model definition file.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     None, stand-alone program.
!
!  6. Error messages :
!
!  7. Remarks :
!
!      Physical grid :
!     -----------------
!
!     The physical grid is defined by a grid counter IX defining the
!     discrete longitude and IY defining the discrete latitude as shown
!     below. For mathemathical convenience, these grid axes will
!     generally be denoted as the X and Y axes. Two-dimensional arrays
!     describing parameters on this grid are given as A(IY,IX).
!
!           IY=NY
!             ^  |      |      |      |      |      |            ^ N
!             |  |------|------|------|------|------|----        |
!             |  |  ::  |  25  |  26  |  27  |  28  |          --|--
!                |------|------|------|------|------|----        |
!           IY=3 |  ::  |  ::  |  9   |  10  |  11  |            |
!                |------|------|------|------|------|----
!           IY=2 |  ::  |   1  |   2  |  ::  |   3  |
!                |------|------|------|------|------|----
!           IY=1 |  ::  |  ::  |  ::  |  ::  |  ::  |
!                +------+------+------+------+------+----
!                  IX=1   IX=2   IX=3   IX=4   IX=5   ---> IX=NX
!
!                                        :: is a land point.
!
!     To reduce memory usage of the model, spectra are stored for sea
!     points only, in a one-dimensional grid with the length NSEA. This
!     grid is called the storage grid. The definition of the counter
!     in the storage grid is graphically depicted above. To transfer
!     data between the two grids, the maps MAPFS and MAPSF are
!     determined. MAPFS gives the counter of the storage grid ISEA
!     for every physical grid point (IY,IX), such that
!
!             MAPFS(IY,IX) = ISEA
!
!     ISEA = 0 corresponds to land points. The map MAPSF gives the grid
!     counters (IY,IX) for a given storage point ISEA.
!
!             MAPSF(ISEA,1) = IX
!             MAPSF(ISEA,2) = IY
!             MAPSF(ISEA,3) = IY+(IX-1)*NY  ( filled during reading )
!
!     Finally, a status maps MAPSTA and MAPST2 are determined, where
!     the status indicator ISTAT = MAPSTA(IY,IX) determines the type
!     of the grid point.
!
!         ISTAT  Means
!       ---------------------------------------------------
!           0    Point excluded from grid.
!        (-)1    Sea point
!        (-)2    "Active" boundary point (data prescribed)
!
!     For ISTAT=0, the secondary status counter ISTA2 is defined as
!
!         ISTA2  Means
!       ---------------------------------------------------
!           0    Land point.
!           1    Point excluded from grid.
!
!     Negative values of ISTAT identify points that are temporarily
!     taken out of the computation. For these points ISTA2 are
!     defined per bit
!
!         BIT    Means
!       ---------------------------------------------------
!          1     Ice flag (1 = ice coverage)
!          2     Dry flag (1 = dry point with depth 0)
!          3     Inferred land in multi-grid model.
!          4     Masking in multi-grid model.
!          5     land point flag for relocatable grid.
!
!      Thus ISTA2=0 for ISTAT<0 is in error, ISTA2=1 means ice cover,
!      ISTA2=3 means ice on dry point, etc.
!
!      Spectral grid :
!     -----------------
!
!     In the spectral grid (and in physical space in general),
!     the cartesian convention for directions is used, i.e., the
!     direction 0 corresponds to waves propagating in the positive
!     X-direction and 90 degr. corresponds to waves propagating in
!     the positive Y-direction. Similar definitions are used for the
!     internal description of winds and currents. Output can obviously
!     be transformed according to any preferred convention.
!
!          ITH=NTH
!             ^  |      |      |      |      |
!             |  |------|------|------|------|----
!             |  |      |      |      |      |      TH(3) = DTH*2.
!                |------|------|------|------|----
!          ITH=2 |      |      |      |      |      TH(2) = DTH
!                |------|------|------|------|----
!          ITH=1 |      |      |      |      |      TH(1) = 0.
!                +------+------+------+------+----
!                  IK=1   IK=2   IK=3   IK=4   ---> IK=NK
!
!     The spectral grid consists of NK wavenumbers. The first
!     wavenumber IK=1 corresponds to the longest wave. The wavenumber
!     grid varies in space, as given by an invariant relative freq.
!     grid and the local depth. The spectral grid furthermore contains
!     NTH directions, equally spaced over a full circle. the first
!     direction corresponds to the direction 0, etc.
!
! (Begin SMC description)
!
!      Spherical Multiple-Cell (SMC) grid
!     -----------------------------------
!
!     SMC grid is a multi-resolution grid using cells of multiple times
!     of each other.  It is similar to the lat-lon grid using rectangular
!     cells but only cells at sea points are retained.  All land points
!     have been removed from the model.  At high latitudes, cells are
!     merged longitudinally to relax the CFL resctiction on time steps.
!     Near coastlines, cells are divided into quarters in a few steps so
!     that high resolution is achieved to refine coastlines and resolve
!     small islands.  At present, three tiers of quarter cells are used.
!     For locating purpose, a usual x-y counter is setup by the smallest
!     cell size and starting from the south-west corner of the usual
!     rectuangular domain.  Each sea cell is then given a pair of x-y
!     index, plus a pair of increments.  These four index are stored in
!     the cell array IJKCel(NCel, 5), each row holds i, j, di, dj, ndps
!     where ndps is an integer depth in metre.  If precision higher than
!     a metre is required, it may use other unit (cm for instance) with a
!     conversion factor.
!
!     For transport calculation, two face arrays, IJKUFc(NUFc, 7) and
!     IJKVFc(NVFc,8), are also created to store the neighbouring cell
!     sequential numbers and the face location and size.  The 3 arrays
!     are calculated outside the wave model and input from text files.
!
!     Boundary condition is added for SMC grid so that it can be used for
!     regional model as well.  Most of the original boundary settings
!     are reclaimed as long as the boundary condition file is provided
!     by a lat-lon grid WW3 model, which will set the interpolation
!     parameters in the boundary condition file.  The NBI number is
!     reset with an input value because the NX-Y double loop overcount
!     the boundary cells for merged cells in the SMC grid.  ISBPI
!     boundary cell mapping array is fine as MAPFS uses duplicated cell
!     number in any merged cell.  From there, all original NBI loops are
!     reusable.
!
!     The whole Arctic can be included in the SMC grid if another option
!     ARC is activated along with the SMC option.  ARC option appends
!     the polar Arctic part above 86N to the existing SMC grid and uses
!     a map-east reference direction for this extra polar region.
!     Because the map-east direction changes with latitude and longitude
!     the wave spectra defined to the map-east direction could not mixed
!     up with the conventional spectra defined to the local east
!     direction.  A rotation sub is provided for convertion from one to
!     another.  Propagation part will be calculated together, including
!     the boundary cells.  The boundary cells are then updated by
!     assigning the corresponding inner cells to them after conversion.
!     Boundary cells are duplicated northmost 4 rows of the global part
!     and they can be excluded for source term and output if required.
!     For convenience, Arctic cellls are all base level cells and are
!     appended to the end of the global cells.  If refined cells were
!     used in the Arctic part, it would not be kept all together, making
!     the sub-loops much more complicated.
!
!     For more information about the SMC grid, please refer to
!     Li, J.G. (2012) Propagation of Ocean Surface Waves on a Spherical
!     Multiple-Cell Grid.  J. Comput. Phys., 231, 8262-8277.  online at
!     http://dx.doi.org/10.1016/j.jcp.2012.08.007
!
! (End SMC description)
!
!     ICEWIND is the scale factor for reduction of wind input by ice
!     concentration. Value specified corresponds to the fractional
!     input for 100% ice concentration. Default is 1.0, meaning that
!     100% ice concentration result in zero wind input.
!     Sin_in_ice=Sin_in_open_water * (1-ICE*ICEWIND)
 
!     -----------------------------------------------------------------*
!  8. Structure :
!
!     ----------------------------------------------------------------
!        1.   Set up grid storage structure.
!                               ( W3NMOD , W3NOUT , W3SETG , W3SETO )
!        2.a  I-O setup.
!          b  Print heading(s).
!        3.   Prepare int. table for dispersion relation   ( DISTAB )
!        4.   Read and process input file up to spectrum.
!          a  Get comment character
!          b  Name of grid
!          c  Define spectrum                              ( W3DIMS )
!        5.   Set-up discrete spectrum.
!          a  Directions.
!          b  Frequency for spectrum.
!        6.   Read and process input file up to numerical parameters
!          a  Set model flags and time steps
!          b  Set / select source term package
!          c  Pre-process namelists.
!          d  Wind input source term.
!          e  Nonlinear interactions.
!          f  Whitecapping term.
!          g  Bottom friction source term.
!          h  Depth indiced breaking source term.
!          i  Triad interaction source term.
!          j  Bottom scattering source term.
!          k  Undefined source term.
!          l  Set / select propagaton scheme
!          m  Parameters for propagation scheme.
!          n  Set misc. parameters (ice, seeding, ...)
!          o  End of namelist processing
!          p  Set various other variables
!        7.   Read and prepare grid.
!          a  Layout of grid
!          b  Storage of grid of grid
!          c  Read bottom depths
!          d  Set up temp map
!          e  Subgrid information
!            1 Info from input file
!            2 Open file and check if necessary
!            3 Read the data
!            4 Limit
!        8    Finalize status maps
!          a  Determine where to get the data
!             Get data in parts from input file
!             ----------------------------------------------------
!          b  Read and update TMPSTA with bound. and excl. points.
!          c  Finalize excluded points
!             ----------------------------------------------------
!             Read data from file
!             ----------------------------------------------------
!          d  Read data from file
!             ----------------------------------------------------
!          e  Get NSEA and other counters
!          f  Set up all maps                              ( W3DIMX )
!        9.   Prepare output boundary points.
!          a  Read
!          b  Update
!       10.   Write model definition file.                 ( W3IOGR )
!     ----------------------------------------------------------------
!
!  9. Switches :
!
!     !/FLX1  Stresses according to Wu (1980).
!     !/FLX2  Stresses according to T&C (1996).
!     !/FLX3  Stresses according to T&C (1996) with cap on Cd.
!     !/FLX4  Stresses according to Hwang (2011).
!
!     !/LN0   No linear input source term.
!     !/SEED  'Seeding' of lowest frequency for sufficiently strong
!             winds. Proxi for linear input.
!     !/LN1   Cavaleri and Melanotte-Rizzoli with Tolman filter.
!     !/LNX   Open slot.
!
!     !/ST0   No source terms included (input/dissipation)
!     !/ST1   WAM-3 physics package.
!     !/ST2   Tolman and Chalikov (1996) physics package.
!     !/ST3   WAM 4+ source terms from P.A.E.M. Janssen and J-R. Bidlot
!     !/ST4   Ardhuin et al. (2009,2010) input and dissipation
!     !/ST6   BYDRZ source term package featuring Donelan et al.
!             (2006) input and Babanin et al. (2001,2010) dissipation.
!     !/STX   Open slot.
!
!     !/NL0   No nonlinear interactions.
!     !/NL1   Discrete interaction approximation (DIA).
!     !/NL2   Exact interactions (WRT).
!     !/NL3   Generalized Multiple DIA (GMD).
!     !/NL4   Two Scale Approximation
!     !/NLX   Open slot.
!     !/NLS   Snl based HF filter.
!
!     !/BT0   No bottom friction included.
!     !/BT1   JONSWAP bottom friction package.
!     !/BT4   SHOWEX bottom friction using movable bed roughness
!                  (Tolman 1994, Ardhuin & al. 2003)
!     !/BTX   Open slot.
!
!     !/IC1   Sink term for interaction with ice (uniform k_i)
!     !/IC2   Sink term for under-ice boundary layer friction
!                  (Liu et al.    1991: JGR 96 (C3), 4605-4621)
!                  (Liu and Mollo 1988: JPO 18       1720-1712)
!     !/IC3   Sink term for interaction with ice (Wang and Shen method)
!                  (Wang and Shen JGR 2010)
!     !/IC4   Sink term for empirical, frequency-dependent attenuation
!                   in ice (Wadhams et al. 1988: JGR 93 (C6) 6799-6818)
!
!     !/DB0   No depth-induced breaking included.
!     !/DB1   Battjes-Janssen depth-limited breaking.
!     !/DBX   Open slot.
!     !/MLIM  Mich-style limiter.
!
!     !/TR0   No triad interactions included.
!     !/TRX   Open slot.
!
!     !/BS0   No bottom scattering included.
!     !/BS1   Routines from F. Ardhuin.
!     !/BSX   Open slot.
!
!     !/XX0   No unclasified source term included.
!     !/XXX   Open slot.
!
!     !/PR1   First order propagation scheme.
!     !/PR2   QUICKEST scheme with ULTIMATE limite and diffusion
!             correction for swell dispersion.
!     !/PR3   Averaging ULTIMATE QUICKEST scheme.
!
!     !/RTD   Rotated regular lat-lon grid.
!     !/SMC   UNO2 scheme on Spherical Multiple-Cell grid.
!     !/ARC   Append the Arctic part to the SMC grid.
!
!     !/MGG   GSE correction for moving grid.
!
!     !/S     Enable subroutine tracing.
!     !/T     Enable test output.
!     !/T0    Enable test output tables for boundary output.
!
!     !/O0    Print equivalent namelist setting to std out.
!     !/O1    Print tables with boundary points as part of output.
!     !/O2    Print MAPSTA as part of output.
!     !/O2a   Print land-sea mask in mask.ww3.
!     !/O2b   Print obstruction data.
!     !/O2c   Print extended status map.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!/
      USE W3TRIAMD
      USE W3GSRUMD, ONLY: W3GRMP
      USE W3ODATMD, ONLY: W3NOUT, W3SETO, W3DMO5
      USE W3IOGRMD, ONLY: W3IOGR
      USE W3SERVMD, ONLY: ITRACE, NEXTLN, EXTCDE
      USE W3ARRYMD, ONLY: INA2R, INA2I
      USE W3DISPMD, ONLY: DISTAB
!/
      USE W3GDATMD
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO
      USE W3ODATMD, ONLY: NBI, NBI2, NFBPO, NBO, NBO2, FLBPI, FLBPO,  &
                          IPBPO, ISBPO, XBPO, YBPO, RDBPO, FNMPRE,    &
                          IHMAX, HSPMIN, WSMULT, WSCUT, FLCOMB,       &
                          NOSWLL
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER      :: NFL = 6
      INTEGER                 :: NDSI, NDSI2, NDSS, NDSM, NDSG, NDSTR,&
                                 IERR, NDSTRC, NTRACE, ITH, IK, ITH0, &
                                 ISP, IYN(NFL), NRLIN, NRSRCE, NRNL,  &
                                 NRBT, NRDB, NRTR, NRBS, NRXX, NRPROP,&
                                 IDLA, IDFM, IX0, IXN, IX, IY, ISEA,  &
                                 IDX, IXO, IDY, IYO, IBA, NBA, ILOOP, &
                                 IFL, NBOTOT, NPO, IP, IX1, IX2, IY1, &
                                 IY2, J, JJ, IXR(4), IYR(4), ISEAI(4),&
                                 IST, NKI, NTHI, NRIC, NRIS,          &
                                 IDFT, NSTAT, NBT, NLAND, NOSW
      INTEGER             :: IBI, IP0, IPN, IPH, IPI
      INTEGER                 :: NCOL =  78
!     INTEGER                 :: NCOL = 130
!     INTEGER                 :: NCOL = 180
!
      INTEGER                 :: NMAP, IMAP
      INTEGER                 :: NMAPB, IMAPB
      INTEGER, ALLOCATABLE    :: TMPSTA(:,:), TMPMAP(:,:), READMP(:,:)
      REAL                    :: RXFR, RFR1, SIGMA, SXFR, FACHF,      &
                                 VSC, VOF,                            &
                                 ZLIM, X, Y, XP,  XO0, YO0, DXO, DYO, &
                                 XO, YO, RD(4), RDTOT,                &
                                 FACTOR, RTH0, FMICHE, RWNDC,         &
                                 WCOR1, WCOR2
      CHARACTER(LEN=4)        :: GSTRG, CSTRG
!
! Variables used to allow spectral output on full grid
!
      INTEGER                 :: P2SF,I1P2SF,I2P2SF
      INTEGER                 :: E3D,I1E3D,I2E3D
      INTEGER                 :: US3D,I1US3D,I2US3D,                  &
                                 TH1MF, I1TH1M, I2TH1M,               &
                                 STH1MF, I1STH1M, I2STH1M,            &
                                 TH2MF, I1TH2M, I2TH2M,               &
                                 STH2MF, I1STH2M, I2STH2M
!
      REAL                    :: CLIN, RFPM, RFHF
      REAL                    :: NLPROP
      REAL                    :: GAMMA
!
      REAL, ALLOCATABLE       :: XGRDIN(:,:), YGRDIN(:,:)
      REAL, ALLOCATABLE       :: ZBIN(:,:), OBSX(:,:), OBSY(:,:)
      REAL, ALLOCATABLE       :: REFD(:,:), REFD2(:,:), REFS(:,:)
      LOGICAL                 :: FLLIN, FLINDS, FLNL, FLBT, FLDB,     &
                                 FLTR, FLBS, FLXX, FLPROP, FLREF,     &
                                 FIRST, CONNCT, FLNEW, INGRID,FLIC,   &
                                 FLIS
      LOGICAL                 :: FLTC96 = .FALSE.
      LOGICAL                 :: FLNMLO = .FALSE.
      LOGICAL                 :: FLSTB2 = .FALSE.
      LOGICAL                 :: FLST4  = .FALSE.
      LOGICAL                 :: FLST6  = .FALSE.
 
!!Li  Add a logical variable to shelter regular grid lines from SMC grid.
      LOGICAL                 :: RGLGRD = .TRUE.
!!Li
      REAL                    :: FACBERG, REFSLOPE
!
 
 
 
      CHARACTER               :: COMSTR*1, PNAME*30, RFORM*16,        &
                                 FROM*4, FNAME*60, TNAME*60, LINE*80, &
                                 STATUS*20,FNAME2*60, PNAME2*40
      CHARACTER(LEN=6)        :: YESXNO(2)
!/ ------------------------------------------------------------------- /
!/ Namelists
!/
      INTEGER                 :: FLAGTR, IHM
      REAL                    :: CFLTM, CICE0, CICEN, PMOVE, XFILT,    &
                                 LICE, XSEED, XR, HSPM, WSM, WSC, STDX,&
                                 STDY, STDT, ICEHMIN, ICEHINIT, ICEWIND
      REAL(8)                 :: GSHIFT ! see notes in WMGHGH
      LOGICAL                 :: FLC, ICEDISP
!
      INTEGER                 :: SWELLFPAR,SDSISO,SDSBRFDF
      REAL                    :: ZWND, ALPHA0, Z0MAX, BETAMAX, SINTHP,&
                                 ZALP, Z0RAT, TAUWSHELTER, SWELLF,    &
                                 SWELLF2,SWELLF3,SWELLF4, SWELLF5,    &
                                 SWELLF6, SWELLF7, FXPM3, FXFM3,      &
                                 WNMEANPTAIL, WNMEANP, STXFTF,        &
                                 STXFTWN, SINBR, FXFMAGE, FXINCUT,    &
                                 FXDSCUT
      REAL                    :: STXFTFTAIL, SDSC1, SDSC2, SDSCUM,    &
                                 SDSC4, SDSC5, SDSC6, WHITECAPWIDTH,  &
                                 SDSSTRAIN, SDSBR, SDSP,              &
                                 SDSCOS, SDSDTH, SDSBCK, SDSABK,      &
                                 SDSPBK, SDSBINT, SDSHCK,             &
                                 SDSBR2, SDSBRF1,                     &
                                 SDSBM0, SDSBM1, SDSBM2, SDSBM3,      &
                                 SDSBM4, SDSLFGEN, SDSHFGEN
!
      REAL                    :: LAMBDA, KDCONV, KDMIN,               &
                                 SNLCS1, SNLCS2, SNLCS3
      REAL                    :: BJALFA, BJGAM
      LOGICAL                 :: BJFLAG
!
      REAL                    :: WDTHCG, WDTHTH
           REAL                    :: UGOBCDEPTH
           LOGICAL                 :: UGOBCAUTO, EXPFSN, EXPFSPSI,     &
                                      EXPFSFCT,IMPFSN
           INTEGER                 :: UNSTSCHEMES(4), UNSTSCHEME
!
      NAMELIST /SLN1/ CLIN, RFPM, RFHF
      NAMELIST /SIN4/ ZWND, ALPHA0, Z0MAX, BETAMAX, SINTHP, ZALP, &
                      TAUWSHELTER, SWELLFPAR, SWELLF,                 &
                      SWELLF2, SWELLF3, SWELLF4, SWELLF5, SWELLF6,    &
                      SWELLF7, Z0RAT, SINBR
      NAMELIST /SNL1/ LAMBDA, NLPROP, KDCONV, KDMIN,                  &
                      SNLCS1, SNLCS2, SNLCS3
      NAMELIST /SDS4/ SDSC1, WNMEANP, WNMEANPTAIL, FXPM3, FXFM3,      &
                      FXFMAGE, SDSC2, SDSCUM, SDSSTRAIN, SDSC4,       &
                      SDSC5, SDSC6, SDSBR, SDSBR2, SDSP, SDSISO,      &
                      SDSBCK, SDSABK, SDSPBK, SDSBINT, SDSHCK,        &
                      SDSDTH, SDSCOS, SDSBRF1, SDSBRFDF,              &
                      SDSBM0, SDSBM1, SDSBM2, SDSBM3, SDSBM4,         &
                      SDSHFGEN, SDSLFGEN, WHITECAPWIDTH, FXINCUT, FXDSCUT
      NAMELIST /SBT1/ GAMMA
      NAMELIST /SDB1/ BJALFA, BJGAM, BJFLAG
!
      NAMELIST /PRO3/ CFLTM, WDTHCG, WDTHTH
           NAMELIST /UNST/ UGOBCAUTO, UGOBCDEPTH,                     &
                           EXPFSN, EXPFSPSI, EXPFSFCT,                &
                           IMPFSN
           NAMELIST /MISC/ CICE0, CICEN, LICE, XSEED, FLAGTR, XP, XR, &
                      XFILT, PMOVE, IHM, HSPM, WSM, WSC, FLC, FMICHE, &
                      RWNDC, FACBERG, NOSW, GSHIFT, WCOR1, WCOR2,     &
                      STDX, STDY, STDT, ICEHMIN, ICEHINIT, ICEDISP,   &
                      ICEWIND
           NAMELIST /OUTS/ P2SF, I1P2SF, I2P2SF,                      &
                             US3D, I1US3D, I2US3D,                    &
                             E3D, I1E3D, I2E3D,                       &
                             TH1MF, I1TH1M, I2TH1M,                   &
                             STH1MF, I1STH1M, I2STH1M,                &
                             TH2MF, I1TH2M, I2TH2M,                   &
                             STH2MF, I1STH2M, I2STH2M
!/
!/ ------------------------------------------------------------------- /
!/
      DATA YESXNO / 'YES/--' , '---/NO' /
      FLNMLO = .TRUE.
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 1.  Set up grid storage structure
!
      CALL W3NMOD ( 1, 6, 6 )
      CALL W3SETG ( 1, 6, 6 )
      CALL W3NOUT (    6, 6 )
      CALL W3SETO ( 1, 6, 6 )
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 2.  IO set-up.
!
      NDSI   = 10
      NDSS   = 90
      NDSM   = 20
!
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSI,FILE=FNMPRE(:J)//'ww3_grid.inp',STATUS='OLD',        &
            ERR=2000,IOSTAT=IERR)
!
      NDSTRC =  6
      NTRACE =  10
      CALL ITRACE ( NDSTRC, NTRACE )
!
      WRITE (NDSO,900)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.a Interpolation table for dispersion relation.
!
      CALL DISTAB
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 3.b Table for friction factors
!
      CALL TABU_FW
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 4   Read and process input file up to spectrum
! 4.a Get comment character
!
      REWIND (NDSI)
      READ (NDSI,'(A)',END=2001,ERR=2002) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
! 4.b Name of grid :
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) GNAME
      WRITE (NDSO,902) GNAME
!
! 4.c Define spectrum
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) RXFR, RFR1, NKI, NTHI, RTH0
      NK     = NKI
      NK2    = NKI + 2
      NTH    = NTHI
      NSPEC  = NK * NTH
      XFR    = MAX ( RXFR , 1.00001 )
      FR1    = MAX ( RFR1 , 1.E-6 )
      DTH    = TPI / REAL(NTH)
      RTH0   = MAX ( -0.5 , MIN ( 0.5 , RTH0 ) )
      WRITE (NDSO,903) NTH, DTH*RADE
      WRITE (NDSO,904) 360./REAL(NTH)*RTH0
      WRITE (NDSO,905) NK, FR1, FR1*XFR**(NK-1), XFR
!
      CALL W3DIMS ( 1, NK, NTH, NDSE, NDST )
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 5.  Initialize spectral parameters.
! 5.a Directions :
!
      DO ITH=1, NTH
        TH  (ITH) = DTH * ( RTH0 + REAL(ITH-1) )
        ESIN(ITH) = SIN ( TH(ITH) )
        ECOS(ITH) = COS ( TH(ITH) )
        IF ( ABS(ESIN(ITH)) .LT. 1.E-5 ) THEN
          ESIN(ITH) = 0.
          IF ( ECOS(ITH) .GT. 0.5 ) THEN
            ECOS(ITH) =  1.
          ELSE
            ECOS(ITH) = -1.
            END IF
          END IF
        IF ( ABS(ECOS(ITH)) .LT. 1.E-5 ) THEN
          ECOS(ITH) = 0.
          IF ( ESIN(ITH) .GT. 0.5 ) THEN
            ESIN(ITH) =  1.
          ELSE
            ESIN(ITH) = -1.
            END IF
          END IF
        ES2 (ITH) = ESIN(ITH)**2
        EC2 (ITH) = ECOS(ITH)**2
        ESC (ITH) = ESIN(ITH)*ECOS(ITH)
        END DO
!
      DO IK=2, NK+1
        ITH0   = (IK-1)*NTH
        DO ITH=1, NTH
          ESIN(ITH0+ITH) = ESIN(ITH)
          ECOS(ITH0+ITH) = ECOS(ITH)
          ES2 (ITH0+ITH) = ES2 (ITH)
          EC2 (ITH0+ITH) = EC2 (ITH)
          ESC (ITH0+ITH) = ESC (ITH)
          END DO
        END DO
!
!   b Frequencies :
!
      SIGMA   = FR1 * TPI / XFR**2
      SXFR    = 0.5 * (XFR-1./XFR)
!
      DO IK=0, NK+1
        SIGMA    = SIGMA * XFR
        SIG (IK) = SIGMA
        DSIP(IK) = SIGMA * SXFR
        END DO
!
      DSII( 1) = 0.5 * SIG( 1) * (XFR-1.)
      DO IK=2, NK-1
        DSII(IK) = DSIP(IK)
        END DO
      DSII(NK) = 0.5 * SIG(NK) * (XFR-1.) / XFR
!
      DO IK=1, NK
        DDEN(IK) = DTH * DSII(IK) * SIG(IK)
        END DO
!
      DO ISP=1, NSPEC
        IK         = 1 + (ISP-1)/NTH
        SIG2 (ISP) = SIG (IK)
        DDEN2(ISP) = DDEN(IK)
        END DO
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 6   Read and process input file up to numerical parameters
! 6.a Set model flags and time steps
!
      WRITE (NDSO,910)
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002)                                 &
        FLDRY, FLCX, FLCY, FLCTH, FLCK, FLSOU
!
      IYN = 2
      IF ( FLDRY ) IYN(1) = 1
      IF ( FLCX  ) IYN(2) = 1
      IF ( FLCY  ) IYN(3) = 1
      IF ( FLCTH ) IYN(4) = 1
      IF ( FLCK  ) IYN(5) = 1
      IF ( FLSOU ) IYN(6) = 1
!
      WRITE (NDSO,911) (YESXNO(IYN(IFL)),IFL=1,NFL)
!
      IF ( .NOT. (FLDRY.OR.FLCX.OR.FLCY.OR.FLCK.OR.FLCTH.OR.FLSOU) ) THEN
          WRITE (NDSE,1010)
          CALL EXTCDE ( 2 )
        END IF
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) DTMAX, DTCFL, DTCFLI, DTMIN
 
      DTMAX  = MAX ( 1. , DTMAX )
!
! Commented to allow very high resolution zooms
!
!      DTCFL  = MAX ( 1. , DTCFL  )
!      DTCFLI = MIN ( DTMAX , MAX ( 1. , DTCFLI ) )
      DTMIN  = MIN ( DTMAX , MAX ( 0. , DTMIN  ) )
      WRITE (NDSO,912) DTMAX, DTCFL, DTCFLI, DTMIN
!
! 6.b Set / select source term package
!
      NRLIN  = 0
      NRSRCE = 0
      NRNL   = 0
      NRBT   = 0
      NRIC   = 0
      NRIS   = 0
      NRDB   = 0
      NRTR   = 0
      NRBS   = 0
      NRXX   = 0
!
      FLLIN  = .TRUE.
      FLINDS = .TRUE.
      FLNL   = .TRUE.
      FLBT   = .TRUE.
      FLIC   = .FALSE.
      FLIS   = .FALSE.
      FLDB   = .TRUE.
      FLTR   = .TRUE.
      FLBS   = .TRUE.
      FLREF  = .FALSE.
      FLXX   = .TRUE.
!
      NRLIN  = NRLIN + 1
!
      NRSRCE = NRSRCE + 1
      FLST4  = .TRUE.
!
      NRNL   = NRNL + 1
!
      NRBT   = NRBT + 1
!
      NRDB   = NRDB + 1
!
      NRTR   = NRTR + 1
      FLTR   = .FALSE.
!
      NRBS   = NRBS + 1
      FLBS   = .FALSE.
!
      NRXX   = NRXX + 1
      FLXX   = .FALSE.
!
      IF ( .NOT.FLLIN .AND.  .NOT.FLINDS .AND.  .NOT.FLNL .AND.        &
           .NOT.FLBT  .AND.  .NOT.FLIC   .AND.  .NOT.FLIS .AND.        &
           .NOT.FLDB  .AND.  .NOT.FLTR   .AND.  .NOT.FLBS .AND.        &
           .NOT.FLXX  .AND.  .NOT.FLREF  .AND.  FLSOU ) THEN
          WRITE (NDSE,1020)
          CALL EXTCDE ( 10 )
        END IF
!
      IF ( ( FLLIN .OR. FLINDS .OR. FLNL .OR. FLBT .OR. FLDB .OR.     &
             FLTR .OR. FLBS .OR. FLREF .OR. FLXX .OR. FLIC )          &
             .AND. .NOT.FLSOU ) THEN
          WRITE (NDSE,1021)
        END IF
!
      IF ( NRLIN .NE. 1 ) THEN
          WRITE (NDSE,1022) NRLIN
          CALL EXTCDE ( 11 )
        END IF
!
      IF ( NRSRCE .NE. 1 ) THEN
          WRITE (NDSE,1023) NRSRCE
          CALL EXTCDE ( 12 )
        END IF
!
      IF ( NRNL .NE. 1 ) THEN
          WRITE (NDSE,1024) NRNL
          CALL EXTCDE ( 13 )
        END IF
!
      IF ( NRBT .NE. 1 ) THEN
          WRITE (NDSE,1025) NRBT
          CALL EXTCDE ( 14 )
        END IF
!
      IF ( NRDB .NE. 1 ) THEN
          WRITE (NDSE,1026) NRDB
          CALL EXTCDE ( 15 )
        END IF
!
      IF ( NRTR .NE. 1 ) THEN
          WRITE (NDSE,1027) NRTR
          CALL EXTCDE ( 16 )
        END IF
!
      IF ( NRBS .NE. 1 ) THEN
          WRITE (NDSE,1028) NRBS
          CALL EXTCDE ( 17 )
        END IF
!
      IF ( NRXX .NE. 1 ) THEN
          WRITE (NDSE,1029) NRXX
          CALL EXTCDE ( 18 )
        END IF
!
      IF ( NRIC .GT. 1 ) THEN
          WRITE (NDSE,1034) NRIC
          CALL EXTCDE ( 19 )
        END IF
!
      IF ( NRIS .GT. 1 ) THEN
          WRITE (NDSE,1036) NRIS
          CALL EXTCDE ( 26 )
        END IF
!
! 6.c Pre-process namelists into scratch file
!
      WRITE (NDSO,915)
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSS,FILE=FNMPRE(:J)//'ww3_grid.scratch',FORM='FORMATTED')
!
      DO
        CALL NEXTLN ( COMSTR , NDSI , NDSE )
        READ (NDSI,'(A)',END=2001,ERR=2002) LINE
        IF ( LINE(1:16) .EQ. 'END OF NAMELISTS' ) THEN
            EXIT
          ELSE
            WRITE (NDSS,'(A)') LINE
          ENDIF
        END DO
!
      WRITE (NDSO,916)
!
! 6.d Define Sin.
! 6.d.1 Stresses
!
! 6.d.2 Linear input
!
      CLIN   = 80.
      RFPM   =  1.
      RFHF   =  0.5
!
      CALL READNL ( NDSS, 'SLN1', STATUS )
      WRITE (NDSO,820) STATUS
      CLIN   = MAX (0.,CLIN)
      RFPM   = MAX (0.,RFPM)
      RFHF   = MAX(0.,MIN (1.,RFHF))
      WRITE (NDSO,821) CLIN, RFPM, RFHF
      SLNC1  = CLIN * (DAIR/DWAT)**2 / GRAV**2
      FSPM   = RFPM
      FSHF   = RFHF
!
! 6.d.3 Exponential input
!
      ZWND   =   10.
      ALPHA0 = 0.0095
      Z0MAX = 0.0
      Z0RAT = 0.04
      BETAMAX   = 1.43
      SINTHP    = 2.
      SWELLF    = 0.66
      SWELLFPAR = 1
      SWELLF2 = -0.018
      SWELLF3 = 0.022
      SWELLF4 = 1.5E5
      SWELLF5 = 1.2
      SWELLF6 = 0.
      SWELLF7 = 360000.
      TAUWSHELTER = 0.3
      ZALP   = 0.006
      SINBR   = 0.
!
      CALL READNL ( NDSS, 'SIN4', STATUS )
      WRITE (NDSO,920) STATUS
      WRITE (NDSO,921) ALPHA0, BETAMAX, SINTHP, Z0MAX, ZALP, ZWND, TAUWSHELTER, &
           SWELLFPAR, SWELLF, SWELLF2, SWELLF3, SWELLF4, SWELLF5, &
           SWELLF6, SWELLF7, Z0RAT
      ZZWND  = ZWND
      AALPHA = ALPHA0
      BBETA  = BETAMAX
      SSINBR  = SINBR
      SSINTHP  = SINTHP
      ZZ0MAX  = Z0MAX
      ZZ0RAT  = Z0RAT
      ZZALP  = ZALP
      TTAUWSHELTER = TAUWSHELTER
      SSWELLF(1) = SWELLF
      SSWELLF(2) = SWELLF2
      SSWELLF(3) = SWELLF3
      SSWELLF(4) = SWELLF4
      SSWELLF(5) = SWELLF5
      SSWELLF(6) = SWELLF6
      SSWELLF(7) = SWELLF7
      SSWELLFPAR = SWELLFPAR
!
! 6.e Define Snl.
!
      LAMBDA =  0.25
      IF ( FLTC96 ) THEN
          NLPROP =  1.00E7
        ELSE IF ( FLST4 ) THEN
          NLPROP =  2.50E7
        ELSE IF ( FLST6 ) THEN
          NLPROP =  3.00E7
        ELSE
          NLPROP =  2.78E7
        END IF
!
      KDCONV =  0.75
      KDMIN  =  0.50
      SNLCS1 =  5.5
      SNLCS2 =  0.833
      SNLCS3 = -1.25
!
      CALL READNL ( NDSS, 'SNL1', STATUS )
      WRITE (NDSO,922) STATUS
      WRITE (NDSO,923) LAMBDA, NLPROP, KDCONV, KDMIN,            &
                       SNLCS1, SNLCS2, SNLCS3
      SNLC1  = NLPROP / GRAV**4
      LAM    = LAMBDA
      KDCON  = KDCONV
      KDMN   = KDMIN
      SNLS1  = SNLCS1
      SNLS2  = SNLCS2
      SNLS3  = SNLCS3
!
      FACHF  = 5.
!
!!/NL3      MSC    = MAX ( 0. , MIN ( 8. , MSC ) )  ! Disabled HLT ca. 2009
!
! 6.f Define Sds.
!
      SDSC1  = 0.0     ! not used in ST4, should be cleaned up
      WNMEANP = 0.5    ! taken from Bidlot et al. 2005
      FXFM3 = 2.5
      FXFMAGE = 0.
      FXINCUT = 0.
      FXDSCUT = 0.
      FXPM3 = 4.
      WNMEANPTAIL = -0.5
      SDSC2     = -2.2E-5
      SDSCUM     = -0.40344
      SDSC4     = 1.
      SDSC5     = 0.
      SDSC6     = 0.3
      SDSBR     = 0.90E-3
      SDSBRFDF  = 0
      SDSBRF1   = 0.5
      SDSP      = 2.   ! this is now fixed in w3sds4, should be cleaned up
      SDSDTH    = 80.
      SDSCOS    = 2.
      SDSISO    = 2
      SDSBM0    = 1.
      SDSBM1    = 0.
      SDSBM2    = 0.
      SDSBM3    = 0.
      SDSBM4    = 0.
      SDSBR2    = 0.8
      SDSBCK    = 0.
      SDSABK    = 1.5
      SDSPBK    = 4.
      SDSBINT   = 0.3
      SDSHCK    = 1.5
      WHITECAPWIDTH = 0.3
      SDSSTRAIN = 0.
      SDSHFGEN  = 0.
      SDSLFGEN  = 0.
!
      CALL READNL ( NDSS, 'SDS4', STATUS )
      WRITE (NDSO,924) STATUS
      WRITE (NDSO,925) SDSC2, SDSBCK, SDSCUM, WNMEANP
      SSDSC(1)   = SDSLFGEN
      SSDSC(2)   = SDSC2
      SSDSC(3)   = SDSCUM
      SSDSC(4)   = SDSC4
      SSDSC(5)   = SDSC5
      SSDSC(6)   = SDSC6
      SSDSC(7)   = WHITECAPWIDTH
      SSDSC(8)   = SDSSTRAIN ! Straining constant ...
      SSDSC(9)   = SDSHFGEN
      SSDSBR   = SDSBR
      SSDSBRF1 = SDSBRF1
      SSDSBRFDF= SDSBRFDF
      SSDSBM(0)   = SDSBM0
      SSDSBM(1)   = SDSBM1
      SSDSBM(2)   = SDSBM2
      SSDSBM(3)   = SDSBM3
      SSDSBM(4)   = SDSBM4
      SSDSBR2  = SDSBR2
      SSDSISO  = SDSISO
      SSDSCOS  = SDSCOS
      SSDSP    = SDSP
      SSDSDTH  = SDSDTH
      WWNMEANP   = WNMEANP
      FFXFM = FXFM3 * TPI
      FFXFA = FXFMAGE * TPI
      FFXFI = FXINCUT * TPI
      FFXFD = FXDSCUT * TPI
      FFXPM = FXPM3 * GRAV / 28.
      WWNMEANPTAIL   = WNMEANPTAIL
      SSDSBCK   = SDSBCK
      SSDSABK   = SDSABK
      SSDSPBK   = SDSPBK
      SSDSBINT  = SDSBINT
      SSDSHCK   = SDSHCK
!
! 6.g Define Sbt.
!
      GAMMA  = -0.067
!
      CALL READNL ( NDSS, 'SBT1', STATUS )
      WRITE (NDSO,926) STATUS
      WRITE (NDSO,927) GAMMA
      SBTC1  = 2. * GAMMA / GRAV
!
! 6.h Define Sdb.
!
      BJALFA = 1.
      BJGAM  = 0.73
      BJFLAG = .TRUE.
!
      CALL READNL ( NDSS, 'SDB1', STATUS )
      WRITE (NDSO,928) STATUS
      BJALFA = MAX ( 0. , BJALFA )
      BJGAM  = MAX ( 0. , BJGAM )
      WRITE (NDSO,929) BJALFA, BJGAM
      IF ( BJFLAG ) THEN
          WRITE (NDSO,*) '      Using Hmax/d ratio only.'
        ELSE
          WRITE (NDSO,*)                                       &
             '      Using Hmax/d in Miche style formulation.'
        END IF
      WRITE (NDSO,*)
      SDBC1  = 0.25 * BJALFA
      SDBC2  = BJGAM
      FDONLY = BJFLAG
!
! 6.i Define Str.
!
      WRITE (NDSO,930)
!
! 6.j Define Sbs.
!
      WRITE (NDSO,932)
!
! 6.k Define Sxx and Sic.
!
! !/XX0      WRITE (NDSO,934)
!
! 6.l Read unstructured data
! initialisation of logical related to unstructured grid
       UGOBCAUTO = .TRUE.
       UGOBCDEPTH= -10.
       EXPFSN    = .TRUE.
       EXPFSPSI  = .FALSE.
       EXPFSFCT  = .FALSE.
       IMPFSN    = .FALSE.
! read data from the unstructured devoted namelist
       CALL READNL ( NDSS, 'UNST', STATUS )
!
! 6.m Select propagation scheme
!
      WRITE (NDSO,950)
!
      NRPROP = 0
      FLPROP = .TRUE.
      PNAME  = '                              '
      PNAME  = '3rd order UQ'
       J = LEN_TRIM(PNAME)
      PNAME  = PNAME(1:J)//' + GSE averaging '
      NRPROP = NRPROP + 1
!
      IF ( (FLCX.OR.FLCY.OR.FLCTH.OR.FLCK) .AND. .NOT. FLPROP ) THEN
          WRITE (NDSE,1030)
          CALL EXTCDE ( 20 )
        END IF
!
      IF ( .NOT.(FLCX.OR.FLCY.OR.FLCTH.OR.FLCK) .AND. FLPROP ) THEN
          WRITE (NDSE,1031)
        END IF
!
      IF ( NRPROP.EQ.0 ) THEN
          WRITE (NDSE,1032)
          CALL EXTCDE ( 21 )
        END IF
!
      IF ( NRPROP .GT. 1 ) THEN
          WRITE (NDSE,1033) NRPROP
          CALL EXTCDE ( 22 )
        END IF
!
! 6.m Parameters for propagation scheme
!
      WRITE (NDSO,951) PNAME
!
      CFLTM  =  0.7
!
      WDTHCG = 1.5
      WDTHTH = WDTHCG
!
      CALL READNL ( NDSS, 'PRO3', STATUS )
      IF ( STATUS(18:18) .EQ. ':' ) STATUS(18:18) = ' '
       IF (GTYPE.NE.UNGTYPE) THEN
          WRITE (NDSO,952) STATUS(1:18)
      CFLTM  = MAX ( 0. , CFLTM )
          WRITE (NDSO,953) CFLTM, WDTHCG
      IF ( WDTHCG*(XFR-1.) .GT. 1. ) WRITE (NDSO,955) 1./(XFR-1.)
          WRITE (NDSO,954) WDTHTH
      IF ( WDTHTH*DTH .GT. 1. ) WRITE (NDSO,955) 1./DTH
          WRITE (NDSO,*)
       ENDIF
      WDCG   = WDTHCG
      WDTH   = WDTHTH
!
      CTMAX  = CFLTM
!
! 6.n Set miscellaneous parameters (ice, seeding, numerics ... )
!
      CICE0  = 0.5
      CICEN  = 0.5
      LICE   = 0.
      ICEHMIN= 0.2  ! the 0.2 value is arbitrary and needs to be tuned.
      ICEHINIT=0.5
      ICEWIND=1.0
      GSHIFT = 0.0D0
      PMOVE  = 0.5
      XSEED  = 1.
      FLAGTR = 0
      XP     = 0.15
      XR     = 0.10
      XFILT  = 0.05
      IHM    = 100
      HSPM   = 0.05
      WSM    = 1.7
      WSC    = 0.333
      FLC    = .TRUE.
      NOSW   = 5
      FMICHE = 1.6
      RWNDC  = 1.
      WCOR1  = 99.
      WCOR2  = 0.
! Variables for Space-Time Extremes
!  Default negative values make w3iogomd use local grid size
      STDX = -1.
      STDY = -1.
      STDT = -1.
      ICEDISP = .FALSE.
! Variables for 3D array output
      E3D=0
      I1E3D=1
      I2E3D=NK
      P2SF   = 0
      I1P2SF = 1
      I2P2SF = 15
      US3D=0
      I1US3D=1
      I2US3D=NK
      TH1MF=0
      I1TH1M=1
      I2TH1M=NK
      STH1MF=0
      I1STH1M=1
      I2STH1M=NK
      TH2MF=0
      I1TH2M=1
      I2TH2M=NK
      STH2MF=0
      I1STH2M=1
      I2STH2M=NK
!
      FACBERG=1.
      WRITE (NDSO,944)
!
!fixme: if USECGICE = .TRUE., don't allow use of IC3MAXTHK<100.0
 
 
!
      CALL READNL ( NDSS, 'OUTS', STATUS )
      WRITE (NDSO,4970) STATUS
!
! output of frequency spectra, th1m ...
!
      E3DF(1,1) = E3D
      E3DF(2,1) = MIN(MAX(1,I1E3D),NK)
      E3DF(3,1) = MIN(MAX(1,I2E3D),NK)
      E3DF(1,2) = TH1MF
      E3DF(2,2) = MIN(MAX(1,I1TH1M),NK)
      E3DF(3,2) = MIN(MAX(1,I2TH1M),NK)
      E3DF(1,3) = STH1MF
      E3DF(2,3) = MIN(MAX(1,I1STH1M),NK)
      E3DF(3,3) = MIN(MAX(1,I2STH1M),NK)
      E3DF(1,4) = TH2MF
      E3DF(2,4) = MIN(MAX(1,I1TH2M),NK)
      E3DF(3,4) = MIN(MAX(1,I2TH2M),NK)
      E3DF(1,5) = STH2MF
      E3DF(2,5) = MIN(MAX(1,I1STH2M),NK)
      E3DF(3,5) = MIN(MAX(1,I2STH2M),NK)
!
! output of microseismic source spectra
!
      P2MSF(1) = P2SF
      P2MSF(2) = MIN(MAX(1,I1P2SF),NK)
      P2MSF(3) = MIN(MAX(1,I2P2SF),NK)
!
! output of Stokes drift profile
!
      US3DF(1) = US3D
      US3DF(2) = I1US3D
      US3DF(3) = I2US3D
      US3DF(2) = MAX( 1 , MIN(NK,US3DF(2)) )
      US3DF(3) = MAX( 1 , MIN(NK,US3DF(3)) )
      WRITE (NDSO,4971) P2MSF(1:3)
      WRITE (NDSO,4972) US3DF(1:3)
      WRITE (NDSO,4973) E3DF(1:3,1)
!
      CALL READNL ( NDSS, 'MISC', STATUS )
      WRITE (NDSO,960) STATUS
!
      IF ( FLAGTR.LT.0 .OR. FLAGTR.GT.6 ) FLAGTR = 0
      CICEN  = MIN ( 1. , MAX ( 0. , CICEN ) )
      ICEWIND  = MIN ( 1. , MAX ( 0. , ICEWIND ) )
      FICEN  = CICEN
      GRIDSHIFT=GSHIFT
      ICE100WIND=ICEWIND
      CICE0  = MIN ( CICEN , MAX ( 0. , CICE0 ) )
      FICEL  = LICE
      IICEHMIN  = ICEHMIN
      IICEHINIT  = ICEHINIT
      IICEDISP= ICEDISP
      PMOVE  = MAX ( 0. , PMOVE )
      PFMOVE = PMOVE
!
! Notes: Presently, if we select CICE0.ne.CICEN requires an obstruction
!     grid, that is initialized with zeros as default.
      IF ( FLAGTR .LT. 3 ) THEN
        IF (CICE0.NE.CICEN) THEN
          CICE0 = CICEN
          IF (STATUS=='(user def. values) :')  WRITE (NDSO,2961)
          END IF
        END IF
      IF ( CICE0.EQ.CICEN .AND. FLAGTR.GE.3 ) FLAGTR = FLAGTR - 2
      WRITE (NDSO,961) CICE0, CICEN
      WRITE (NDSO,8972) ICEWIND
      FICE0  = CICE0
! Variables for Space-Time Extremes
      STEXU = STDX
      IF ( STDY .LE. 0. ) THEN
        STDY = STDX
      END IF
      STEYU = STDY
      STEDU = STDT
      IF ( STDX .GT. 0 ) THEN
         WRITE (NDSO,1040) STDX
         WRITE (NDSO,1041) STDY
      ELSE
         WRITE (NDSO,1042)
      END IF
      IF ( STDT .GT. 0 ) THEN
         WRITE (NDSO,1043) STDT
      ELSE
         WRITE (NDSO,1044)
      END IF
!
      FACSD  = XSEED
!
      XP     = MAX ( 1.E-6 , XP )
      XR     = MAX ( 1.E-6 , XR )
      XREL   = XR
      XFILT  = MAX ( 0. , XFILT )
      XFLT   = XFILT
      WRITE (NDSO,965) XP, XR, XFILT
      FACP   = XP / PI * 0.62E-3 * TPI**4 / GRAV**2
!
      IHMAX  = MAX ( 50, IHM )
      HSPMIN = MAX ( 0.0001 , HSPM )
      WSMULT = MAX ( 1. , WSM )
      WSCUT  = MIN ( 1.0001 , MAX ( 0. , WSC ) )
      FLCOMB = FLC
      NOSWLL = MAX ( 1 , NOSW )
        IF ( FLCOMB ) THEN
          J      = 1
        ELSE
          J      = 2
        END IF
      WRITE (NDSO,966) IHMAX, HSPMIN, WSMULT, WSCUT, YESXNO(J), NOSWLL
!!    WRITE (NDSO,966) IHMAX, HSPMIN, WSMULT, WSCUT, YESXNO(J)
!
      FHMAX  = MAX ( 0.01 , FMICHE )
      J      = 2
      J      = 1
      WRITE (NDSO,967) FHMAX, FHMAX/SQRT(2.), YESXNO(J)
      IF ( FHMAX.LT.0.50 .AND. J.EQ.1 ) WRITE (NDST,968)
      WRITE (NDSO,*)
!
! 6.x Read values for FLD stress calculation
!
! 6.o End of namelist processing
!
      CLOSE (NDSS,STATUS='DELETE')
!
      IF ( FLNMLO ) THEN
        WRITE (NDSO,917)
          WRITE (NDSO,2820) CLIN, RFPM, RFHF
        IF ( .NOT. FLSTB2 ) THEN
        ELSE
          END IF
!
          WRITE (NDSO,2920) ZWND, ALPHA0, Z0MAX, BETAMAX, SINTHP, ZALP,   &
            TAUWSHELTER, SWELLFPAR, SWELLF, SWELLF2, SWELLF3, SWELLF4, &
            SWELLF5, SWELLF6, SWELLF7, Z0RAT, SINBR
          WRITE (NDSO,2922) LAMBDA, NLPROP, KDCONV, KDMIN,       &
                            SNLCS1, SNLCS2, SNLCS3
          WRITE (NDSO,2924) SDSC1, SDSC2, SDSCUM, SDSSTRAIN, SDSC4, SDSC5, SDSC6, &
                    WNMEANP, FXPM3, FXFM3, FXFMAGE, FXINCUT, FXDSCUT, &
                    SDSBINT, SDSBCK, SDSABK, SDSPBK, SDSHCK,          &
                    SDSBR, SDSSTRAIN, SDSBR2, SDSP, SDSISO, SDSCOS,   &
                    SDSDTH, SDSBRF1, SDSBRFDF,                        &
                    SDSBM0, SDSBM1, SDSBM2, SDSBM3, SDSBM4,           &
                    WHITECAPWIDTH, SDSLFGEN, SDSHFGEN
          WRITE (NDSO,2926) GAMMA
          IF ( BJFLAG ) THEN
              WRITE (NDSO,2928) BJALFA, BJGAM, '.TRUE.'
            ELSE
              WRITE (NDSO,2928) BJALFA, BJGAM, '.FALSE.'
            END IF
          WRITE (NDSO,2953) CFLTM, WDTHCG, WDTHTH
!
        WRITE (NDSO,2956) UGOBCAUTO, UGOBCDEPTH,                      &
                          EXPFSN, EXPFSPSI, EXPFSFCT, IMPFSN
!
        WRITE (NDSO,2976)    P2SF, I1P2SF, I2P2SF, US3D, I1US3D, I2US3D,  &
                             E3D, I1E3D, I2E3D,                       &
                             TH1MF, I1TH1M, I2TH1M,                   &
                             STH1MF, I1STH1M, I2STH1M,                &
                             TH2MF, I1TH2M, I2TH2M,                   &
                             STH2MF, I1STH2M, I2STH2M
!
 
 
!
        IF ( FLCOMB ) THEN
          WRITE (NDSO,2966) CICE0, CICEN, LICE, PMOVE, XSEED, FLAGTR, &
                                XP, XR, XFILT, IHMAX, HSPMIN, WSMULT, &
                                WSCUT, '.TRUE.', NOSWLL, FHMAX,       &
                                RWNDC, WCOR1, WCOR2, FACBERG, GSHIFT, &
                                STDX, STDY, STDT, ICEHMIN, ICEHINIT,  &
                                ICEDISP, ICEWIND
        ELSE
          WRITE (NDSO,2966) CICE0, CICEN, LICE, PMOVE, XSEED, FLAGTR, &
                                XP, XR, XFILT, IHMAX, HSPMIN, WSMULT, &
                                WSCUT, '.FALSE.', NOSWLL, FHMAX,      &
                                RWNDC, WCOR1, WCOR2, FACBERG, GSHIFT, &
                                STDX, STDY, STDT, ICEHMIN, ICEHINIT,  &
                                ICEDISP, ICEWIND
 
          END IF
!
        WRITE (NDSO,918)
        END IF
!
! 6.p Set various other values ...
! ... Tail in integration       --> scale factor for A to E conv
!
      FTE    = 0.25 * SIG(NK)      * DTH * SIG(NK)
      FTF    = 0.20                * DTH * SIG(NK)
      FTWN   = 0.20 * SQRT(GRAV)   * DTH * SIG(NK)
      FTTR   = FTF
      FTWL   = GRAV / 6. / SIG(NK) * DTH * SIG(NK)
!
      STXFTF      = 1/(FACHF-1.-WNMEANP*2)                       &
                           * SIG(NK)**(2+WNMEANP*2) * DTH
      STXFTFTAIL  = 1/(FACHF-1.-WNMEANPTAIL*2)                   &
                           * SIG(NK)**(2+WNMEANPTAIL*2) * DTH
      STXFTWN = 1/(FACHF-1.-WNMEANP*2) * SIG(NK)**(2)            &
                 * (SIG(NK)/SQRT(GRAV))**(WNMEANP*2)   * DTH
      SSTXFTF     = STXFTF
      SSTXFTFTAIL = STXFTFTAIL
      SSTXFTWN    = STXFTWN
!
! ... High frequency cut-off
!
      FXFM   = 2.5
      FXPM   = 4.0
      FXPM   = FXPM * GRAV / 28.
      FXFM   = FXFM * TPI
      XFC    = 3.0
!
      FACTI1 = 1. / LOG(XFR)
      FACTI2 = 1. - LOG(TPI*FR1) * FACTI1
!
! Setting of FACHF moved to before !/NL2 set-up for consistency
!
      FACHFA = XFR**(-FACHF-2)
      FACHFE = XFR**(-FACHF)
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 7.  Read and prepare the grid.
! 7.a Type of grid
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) GSTRG, FLAGLL, CSTRG
      SELECT CASE (TRIM(GSTRG))
        CASE ('RECT')
          GTYPE = RLGTYPE
          WRITE (NDSO,3000) 'rectilinear'
        CASE ('CURV')
          GTYPE = CLGTYPE
          WRITE (NDSO,3000) 'curvilinear'
        CASE ('UNST')
          GTYPE = UNGTYPE
          WRITE (NDSO,3000) 'unstructured'
        CASE DEFAULT
          WRITE (NDSE,1007) TRIM(GSTRG)
          CALL EXTCDE ( 25 )
        END SELECT
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
!
      IF ( FLAGLL ) THEN
          FACTOR = 1.
          WRITE (NDSO,3001) 'spherical'
        ELSE
          FACTOR = 1.E-3
          WRITE (NDSO,3001) 'Cartesian'
        END IF
!
!     Only process grid closure string for logically rectangular grids.
!     Closure setting for unstructured grids is NONE.
      ICLOSE = ICLOSE_NONE
      IF ( GTYPE.NE.UNGTYPE ) THEN
          SELECT CASE (TRIM(CSTRG))
            CASE ('NONE')
              ICLOSE = ICLOSE_NONE
              WRITE (NDSO,3002) 'none'
            CASE ('SMPL')
              ICLOSE = ICLOSE_SMPL
              WRITE (NDSO,3002) 'simple'
            CASE ('TRPL')
              WRITE (NDSE,'(/2A)') ' *** WARNING WW3_GRID: TRIPOLE ',  &
              'GRID CLOSURE IMPLEMENTATION IS INCOMPLETE ***'
              ICLOSE = ICLOSE_TRPL
              WRITE (NDSO,3002) 'tripole'
              IF ( GTYPE.EQ.RLGTYPE ) THEN
                  WRITE (NDSE,1009)
                  CALL EXTCDE ( 25 )
                END IF
            CASE DEFAULT
              ! Check for old style GLOBAL input
              SELECT CASE (TRIM(CSTRG))
                CASE ('T','t','.TRU','.tru')
                  ICLOSE = ICLOSE_SMPL
                  WRITE (NDSO,3002) 'simple'
                  WRITE (NDSE,1013)
                CASE ('F','f','.FAL','.fal')
                  ICLOSE = ICLOSE_NONE
                  WRITE (NDSO,3002) 'none'
                  WRITE (NDSE,1013)
                CASE DEFAULT
                  WRITE (NDSE,1012) TRIM(CSTRG)
                  CALL EXTCDE ( 25 )
                END SELECT
            END SELECT
          IF ( ICLOSE.NE.ICLOSE_NONE .AND. .NOT.FLAGLL ) THEN
              WRITE (NDSE,1008)
              CALL EXTCDE ( 25 )
            END IF
        END IF !GTYPE.NE.UNGTYPE
!
! 7.b Size of grid
!
      IF ( GTYPE.NE.UNGTYPE) THEN
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) NX, NY
          NX     = MAX ( 3 , NX )
          NY     = MAX ( 3 , NY )
          WRITE (NDSO,3003) NX, NY
        ELSE
          NY =1
        END IF
!
! Propagation specific to unstructured grids
!
      IF ( GTYPE.EQ.UNGTYPE) THEN
        UNSTSCHEMES(:)=0
        IF (EXPFSN)   UNSTSCHEMES(1)=1
        IF (EXPFSPSI) UNSTSCHEMES(2)=1
        IF (EXPFSFCT) UNSTSCHEMES(3)=1
        IF (IMPFSN)   UNSTSCHEMES(4)=1
 
        DO IX=1,4
          IF (UNSTSCHEMES(IX).EQ.1) THEN
            UNSTSCHEME=IX
            EXIT
            END IF
          END DO
 
        SELECT CASE (UNSTSCHEME)
        CASE (1)
          FSN = EXPFSN
          PNAME2 = 'N Explicit (Fluctuation Splitting) '
        CASE (2)
          FSPSI = EXPFSPSI
          PNAME2 = 'PSI Explicit (Fluctuation Splitting)  '
        CASE (3)
          FSFCT = EXPFSFCT
          PNAME2 = ' Flux Corrected Transport Explicit'
        CASE (4)
          FSNIMP = IMPFSN
          PNAME2 = 'N Implicit (Fluctuation Splitting) '
          END SELECT
!
        IF (SUM(UNSTSCHEMES).GT.1) WRITE(NDSO,1035)
        WRITE (NDSO,2951) PNAME2
        END IF
!
! 7.c Grid coordinates (branch here based on grid type)
!
      IF ( GTYPE.NE.UNGTYPE) ALLOCATE ( XGRDIN(NX,NY), YGRDIN(NX,NY) )
      SELECT CASE ( GTYPE )
!
! 7.c.1 Rectilinear grid
!
        CASE ( RLGTYPE )
!
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) SX, SY, VSC
          VSC    = MAX ( 1.E-7 , VSC )
          SX     = SX / VSC
          SY     = SY / VSC
          SX     = MAX ( 1.E-7 , SX )
          SY     = MAX ( 1.E-7 , SY )
          IF ( ICLOSE.EQ.ICLOSE_SMPL ) SX = 360. / REAL(NX)
!
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) X0, Y0, VSC
          VSC    = MAX ( 1.E-7 , VSC )
          X0     = X0 / VSC
          Y0     = Y0 / VSC
!
          IF ( FLAGLL ) THEN
              WRITE (NDSO,3004) FACTOR*SX, FACTOR*SY,         &
                     FACTOR*X0, FACTOR*(X0+REAL(NX-1)*SX),    &
                     FACTOR*Y0, FACTOR*(Y0+REAL(NY-1)*SY)
            ELSE
              WRITE (NDSO,3005) FACTOR*SX, FACTOR*SY,         &
                     FACTOR*X0, FACTOR*(X0+REAL(NX-1)*SX),    &
                     FACTOR*Y0, FACTOR*(Y0+REAL(NY-1)*SY)
            END IF
!
          DO IY=1, NY
            DO IX=1, NX
              XGRDIN(IX,IY) = X0 + REAL(IX-1)*SX
              YGRDIN(IX,IY) = Y0 + REAL(IY-1)*SY
              END DO
            END DO
!
! 7.c.2 Curvilinear grid
!
        CASE ( CLGTYPE )
!
! 7.c.2.a Process x-coordinates
!
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) NDSG, VSC, VOF, &
                                          IDLA, IDFM, RFORM, FROM, FNAME
!
          IF (IDLA.LT.1 .OR. IDLA.GT.4) IDLA   = 1
          IF (IDFM.LT.1 .OR. IDFM.GT.3) IDFM   = 1
!
          WRITE (NDSO,3006) NDSG, VSC, VOF, IDLA, IDFM
          IF (IDFM.EQ.2) WRITE (NDSO,3008) TRIM(RFORM)
          IF (FROM.EQ.'NAME' .AND. NDSG.NE.NDSI) &
              WRITE (NDSO,3009) TRIM(FNAME)
!
          IF ( NDSG .EQ. NDSI ) THEN
              IF ( IDFM .EQ. 3 ) THEN
                  WRITE (NDSE,1004) NDSG
                  CALL EXTCDE (23)
                ELSE
                  CALL NEXTLN ( COMSTR , NDSI , NDSE )
                END IF
            ELSE
              IF ( IDFM .EQ. 3 ) THEN
                  IF (FROM.EQ.'NAME') THEN
                      OPEN (NDSG,FILE=TRIM(FNMPRE)//TRIM(FNAME),&
                            FORM='UNFORMATTED',                 &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG,                               &
                            FORM='UNFORMATTED',                 &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                ELSE
                  IF (FROM.EQ.'NAME') THEN
                      OPEN (NDSG,FILE=TRIM(FNMPRE)//TRIM(FNAME),&
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG,                               &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                END IF !IDFM
            END IF !NDSG
!
          CALL INA2R ( XGRDIN, NX, NY, 1, NX, 1, NY, NDSG, NDST, NDSE, &
                       IDFM, RFORM, IDLA, VSC, VOF)
!
! 7.c.2.b Process y-coordinates
!
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) NDSG, VSC, VOF, &
                                          IDLA, IDFM, RFORM, FROM, FNAME
!
          IF (IDLA.LT.1 .OR. IDLA.GT.4) IDLA   = 1
          IF (IDFM.LT.1 .OR. IDFM.GT.3) IDFM   = 1
!
          WRITE (NDSO,3007) NDSG, VSC, VOF, IDLA, IDFM
          IF (IDFM.EQ.2) WRITE (NDSO,3008) TRIM(RFORM)
          IF (FROM.EQ.'NAME' .AND. NDSG.NE.NDSI) &
              WRITE (NDSO,3009) TRIM(FNAME)
!
          IF ( NDSG .EQ. NDSI ) THEN
              IF ( IDFM .EQ. 3 ) THEN
                  WRITE (NDSE,1004) NDSG
                  CALL EXTCDE (23)
                ELSE
                  CALL NEXTLN ( COMSTR , NDSI , NDSE )
                END IF
            ELSE
              IF ( IDFM .EQ. 3 ) THEN
                  IF (FROM.EQ.'NAME') THEN
                      OPEN (NDSG,FILE=TRIM(FNMPRE)//TRIM(FNAME),&
                            FORM='UNFORMATTED',                 &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG,                               &
                            FORM='UNFORMATTED',                 &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                ELSE
                  IF (FROM.EQ.'NAME') THEN
                      OPEN (NDSG,FILE=TRIM(FNMPRE)//TRIM(FNAME),&
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG,                               &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                END IF !IDFM
            END IF !NDSG
!
          CALL INA2R ( YGRDIN, NX, NY, 1, NX, 1, NY, NDSG, NDST, NDSE, &
                       IDFM, RFORM, IDLA, VSC, VOF)
!
! 7.c.2.c Check for obvious errors in grid definition or input
!
! ....... Check for inverted grid (can result from wrong IDLA)
          IF ( (XGRDIN(2,1)-XGRDIN(1,1))*(YGRDIN(1,2)-YGRDIN(1,1)) .LT. &
               (YGRDIN(2,1)-YGRDIN(1,1))*(XGRDIN(1,2)-XGRDIN(1,1)) ) THEN
             WRITE (NDSE,1011) IDLA
!.........Notes: here, we are checking to make sure that the j axis is ~90 degrees
!................counter-clockwise from the i axis (the standard cartesian setup).
!................So, it is a check on the handedness of the grid.
!................We have confirmed for one case that a left-handed grid produces
!................errors in SCRIP. We have not confirmed that left-handed grids necessarily
!................produce errors in single-grid simulations, or that they necessarily
!................produce errors in all multi-grid simulations.
!................Note that transposing or flipping a grid will generally change the handedness.
             CALL EXTCDE (25)
          END IF
!
! 7.c.3 Unstructured grid
!
        CASE ( UNGTYPE )
!
             MAXX = 0.
             MAXY = 0.
             DXYMAX = 0.
             WRITE (NDSO,1150)
 
        END SELECT !GTYPE
!
! 7.d Depth information for grid
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) ZLIM, DMIN, NDSG, VSC, IDLA,    &
                                      IDFM, RFORM, FROM, FNAME
!
      DMIN    = MAX ( 1.E-3 , DMIN )
      IF (   ABS(VSC) .LT. 1.E-7  ) VSC    = 1.
      IF (IDLA.LT.1 .OR. IDLA.GT.4) IDLA   = 1
      IF (IDFM.LT.1 .OR. IDFM.GT.3) IDFM   = 1
!
      WRITE (NDSO,972) NDSG, ZLIM, DMIN, VSC, IDLA, IDFM
      IF (IDFM.EQ.2) WRITE (NDSO,973) TRIM(RFORM)
      IF (FROM.EQ.'NAME' .AND. NDSG.NE.NDSI) &
          WRITE (NDSO,974) TRIM(FNAME)
!
! 7.e Read bottom depths
!
      IF ( GTYPE.NE.UNGTYPE ) THEN
!
! Reading depths on structured grid
!
!Li Suspended for SMC grid, which uses depth stored in its cell array.
!Li               JGLi15Oct2014
        IF( RGLGRD ) THEN
!Li
          IF ( NDSG .EQ. NDSI ) THEN
              IF ( IDFM .EQ. 3 ) THEN
                  WRITE (NDSE,1004) NDSG
                  CALL EXTCDE (23)
                ELSE
                  CALL NEXTLN ( COMSTR , NDSI , NDSE )
                END IF
            ELSE  ! NDSG.NE.NDSI
              IF ( IDFM .EQ. 3 ) THEN
                  IF (FROM.EQ.'NAME') THEN
                      J = LEN_TRIM(FNMPRE)
                      OPEN (NDSG,FILE=TRIM(FNMPRE(:J))//TRIM(FNAME), &
                            FORM='UNFORMATTED',&
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG, FORM='UNFORMATTED',                &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                ELSE
                  IF (FROM.EQ.'NAME') THEN
                      J      = LEN_TRIM(FNMPRE)
                      OPEN (NDSG,FILE=TRIM(FNMPRE(:J))//TRIM(FNAME),  &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSG,                                     &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                END IF
            END IF  !( NDSG .EQ. NDSI )
!
!Li     End of RGLGRD block
        ENDIF
!Li
!
          ALLOCATE ( ZBIN(NX,NY), OBSX(NX,NY), OBSY(NX,NY) )
!
!         Initialize subgrid obstructions with zeros.
          ZBIN(:,:)=0.
          OBSX(:,:)=0.
          OBSY(:,:)=0.
 
!Li   Suspend read depth file.     JGLi15Oct2014
        IF( RGLGRD ) THEN
!Li
          CALL INA2R ( ZBIN, NX, NY, 1, NX, 1, NY, NDSG, NDST, NDSE,      &
                       IDFM, RFORM, IDLA, VSC, 0.0)
!Li     End of RGLGRD block
        ENDIF
!Li
!
        ELSE
!
! Reading depths on unstructured grid (this also sets number of mesh points, NX)
!
          CALL READMSH(NDSS,FNAME)
          ALLOCATE(ZBIN(NX, NY),OBSX(NX,NY),OBSY(NX,NY))
          ZBIN(:,1) = VSC*XYB(:,3)
!
! subgrid obstructions are not yet handled in unstructured grids
!
          OBSX(:,:)=0.
          OBSY(:,:)=0.
 
        END IF
!
! 7.f Set up temporary map
!
      ALLOCATE ( TMPSTA(NY,NX), TMPMAP(NY,NX) )
      TMPSTA = 0
!
      DO IY=1, NY
        DO IX=1, NX
          IF ( ZBIN(IX,IY) .LE. ZLIM ) TMPSTA(IY,IX) = 1
          END DO
        END DO
!
!Li   Suspended for SMC grid.  JGLi15Oct2014
      IF( RGLGRD ) THEN
!Li
!
! 7.g Subgrid information
!
      TRFLAG = FLAGTR
      IF ( TRFLAG.GT.6 .OR. TRFLAG.LT.0 ) TRFLAG = 0
!
      IF ( TRFLAG .EQ. 0 ) THEN
          WRITE (NDSO,976) 'Not available.'
          WRITE (NDSO,*)
        ELSE IF ( TRFLAG.EQ.1 .OR. TRFLAG.EQ.3 .OR. TRFLAG.EQ.5 ) THEN
          WRITE (NDSO,976) 'In between grid points.'
        ELSE
          WRITE (NDSO,976) 'At grid points.'
        END IF
!
      IF ( TRFLAG .NE. 0 ) THEN
!
! 7.g.1 Info from input file
!
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) NDSTR, VSC, IDLA, IDFT, RFORM, &
                                          FROM, TNAME
!
          IF (   ABS(VSC) .LT. 1.E-7  ) VSC    = 1.
          IF (IDLA.LT.1 .OR. IDLA.GT.4) IDLA   = 1
          IF (IDFT.LT.1 .OR. IDFT.GT.3) IDFT   = 1
!
          WRITE (NDSO,977) NDSTR, VSC, IDLA, IDFT
          IF (IDFT.EQ.2) WRITE (NDSO,973) RFORM
          IF (FROM.EQ.'NAME' .AND. NDSG.NE.NDSTR) WRITE (NDSO,974) TNAME
!
! 7.g.2 Open file and check if necessary
!
          IF ( NDSTR .EQ. NDSI ) THEN
              IF ( IDFT .EQ. 3 ) THEN
                  WRITE (NDSE,1004) NDSTR
                  CALL EXTCDE (23)
                ELSE
                  CALL NEXTLN ( COMSTR , NDSI , NDSE )
                END IF
            ELSE IF ( NDSTR .EQ. NDSG ) THEN
              IF ( ( IDFM.EQ.3 .AND. IDFT.NE.3 ) .OR.                 &
                   ( IDFM.NE.3 .AND. IDFT.EQ.3 ) ) THEN
                  WRITE (NDSE,1005) IDFM, IDFT
                  CALL EXTCDE (24)
                END IF
            ELSE
              IF ( IDFT .EQ. 3 ) THEN
                  IF (FROM.EQ.'NAME') THEN
                      J      = LEN_TRIM(FNMPRE)
                      OPEN (NDSTR,FILE=FNMPRE(:J)//TNAME,             &
                            FORM='UNFORMATTED',STATUS='OLD',ERR=2000, &
                            IOSTAT=IERR)
                    ELSE
                      OPEN (NDSTR,           FORM='UNFORMATTED',      &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                ELSE
                  IF (FROM.EQ.'NAME') THEN
                      J      = LEN_TRIM(FNMPRE)
                      OPEN (NDSTR,FILE=FNMPRE(:J)//TNAME,             &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSTR,                                    &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                END IF
            END IF
!
! 7.g.3 Read the data
!
          CALL INA2R ( OBSX, NX, NY, 1, NX, 1, NY, NDSTR, NDST, NDSE, &
                       IDFM, RFORM, IDLA, VSC, 0.0)
!
          IF ( NDSTR .EQ. NDSI ) CALL NEXTLN ( COMSTR , NDSI , NDSE )
!
          CALL INA2R ( OBSY, NX, NY, 1, NX, 1, NY, NDSTR, NDST, NDSE, &
                       IDFM, RFORM, IDLA, VSC, 0.0)
!
! 7.g.4 Limit
!
          DO IX=1, NX
            DO IY=1, NY
              OBSX(IX,IY) = MAX( 0. , MIN(1.,OBSX(IX,IY)) )
              OBSY(IX,IY) = MAX( 0. , MIN(1.,OBSY(IX,IY)) )
              END DO
            END DO
!
          WRITE (NDSO,*)
!
        END IF
!
!Li     End of RGLGRD block
        ENDIF
!Li
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 8.  Finalize status maps
! 8.a Automatic OBC detection for ug grids
!
      IF (GTYPE.EQ.UNGTYPE.AND.UGOBCAUTO)  &
        CALL UG_GETOPENBOUNDARY(TMPSTA,ZBIN,UGOBCDEPTH)
!
! 8.b Determine where to get the data
!
      CALL NEXTLN ( COMSTR , NDSI , NDSE )
      READ (NDSI,*,END=2001,ERR=2002) NDSTR, IDLA, IDFT, RFORM,     &
                                      FROM, TNAME
!
! ... Data to be read in parts
!
      IF ( FROM .EQ. 'PART' ) THEN
!
! 8.b Update TMPSTA with input boundary data (ILOOP=1)
!                        and excluded points (ILOOP=2)
!
        IF ( ICLOSE .EQ. ICLOSE_TRPL ) THEN
          WRITE(NDSE,*)'PROGRAM W3GRID STATUS MAP CALCULATION IS '//   &
          'NOT TESTED FOR TRIPOLE GRIDS FOR CASE WHERE USER OPTS '//   &
          'TO READ DATA IN PARTS. STOPPING NOW (107).'
          CALL EXTCDE ( 107 )
        END IF
 
        DO ILOOP=1, 2
!
          IF ( ILOOP .EQ. 1 ) THEN
              WRITE (NDSO,979) 'boundary points'
              NSTAT  = 2
            ELSE
              WRITE (NDSO,979) 'excluded points'
              NSTAT  = -1
            END IF
          FIRST  = .TRUE.
!
          DO
            CALL NEXTLN ( COMSTR , NDSI , NDSE )
            READ (NDSI,*,END=2001,ERR=2002) IX, IY, CONNCT
!
! ... Check if last point reached.
!
            IF (IX.EQ.0 .AND. IY.EQ.0) EXIT
!
! ... Check if point in grid.
!
            IF (GTYPE.EQ.UNGTYPE.AND.UGOBCAUTO) CYCLE
            IF (IX.LT.1 .OR. IX.GT.NX .OR.  IY.LT.1 .OR. IY.GT.NY) THEN
                WRITE (NDSO,981)
                WRITE (NDSO,*) '       ', IX, IY
                CYCLE
              END IF
!
! ... Check if intermediate points are to be added.
!
            IF ( CONNCT .AND. .NOT.FIRST ) THEN
                IDX    = IX - IXO
                IDY    = IY - IYO
                IF ( IDX.EQ.0 .OR. IDY.EQ.0 .OR.                      &
                    ABS(IDX).EQ.ABS(IDY) ) THEN
                    NBA    = MAX ( MAX(ABS(IDX),ABS(IDY))-1 , 0 )
                    IF (IDX.NE.0) IDX = SIGN(1,IDX)
                    IF (IDY.NE.0) IDY = SIGN(1,IDY)
                    IX     = IXO
                    IY     = IYO
                    DO IBA=1, NBA
                      IX     = IX + IDX
                      IY     = IY + IDY
                      IF ( TMPSTA(IY,IX).EQ.1 .OR. J.EQ.2 ) THEN
                          TMPSTA(IY,IX) = NSTAT
                      ELSE
                         WRITE(NDSO,*) 'WARNING: POINT (',IX,',',IY,  &
                                  ') CANNOT BE GIVEN THE STATUS ',NSTAT
                        END IF
                      END DO
                    IX     = IX + IDX
                    IY     = IY + IDY
                  ELSE
                    WRITE (NDSO,982)
                    WRITE (NDSO,*) '       ', IX , IY
                    WRITE (NDSO,*) '       ', IXO, IYO
                  END IF
              END IF
!
! ... Check if point itself is to be added
!
            IF ( TMPSTA(IY,IX).EQ.1 .OR. J.EQ.2 ) THEN
                TMPSTA(IY,IX) = NSTAT
              END IF
!
! ... Save data of previous point
!
            IXO    = IX
            IYO    = IY
            FIRST  = .FALSE.
!
! ... Branch back to read.
!
            END DO
!
! 8.c Final processing excluded points
!
          IF ( ILOOP .EQ. 2 ) THEN
!
              DO
                CALL NEXTLN ( COMSTR , NDSI , NDSE )
                READ (NDSI,*,END=2001,ERR=2002) IX, IY
!
! ... Check if last point reached.
!
                IF (IX.EQ.0 .AND. IY.EQ.0) EXIT
!
! ... Check if point in grid.
!
                IF (IX.LT.1 .OR. IX.GT.NX .OR. IY.LT.1 .OR. IY.GT.NY) THEN
                    WRITE (NDSO,981)
                    WRITE (NDSO,*) '       ', IX, IY
                    CYCLE
                  END IF
!
! ... Check if point already excluded
!
                IF ( TMPSTA(IY,IX) .EQ. NSTAT ) THEN
                    WRITE (NDSO,1981)
                    WRITE (NDSO,*) '       ', IX, IY
                    CYCLE
                  END IF
!
! ... Search for points to exclude
!
                TMPMAP = TMPSTA
                J      = 1
                IX1    = IX
                IY1    = IY
!
                JJ     = TMPSTA(IY,IX)
                TMPSTA(IY,IX) = NSTAT
!
                DO
                  NBT    = 0
!
                  DO IX=MAX(1,IX1-J), MIN(IX1+J,NX)
                    DO IY=MAX(1,IY1-J), MIN(IY1+J,NY)
                      IF ( TMPSTA(IY,IX) .EQ. JJ ) THEN
                          IF (IX.GT.1) THEN
                              IF (TMPSTA(IY  ,IX-1).EQ.NSTAT           &
                                      .AND. TMPMAP(IY  ,IX-1).EQ.JJ )  &
                              TMPSTA(IY,IX) = NSTAT
                            END IF
                          IF (IX.LT.NX) THEN
                              IF (TMPSTA(IY  ,IX+1).EQ.NSTAT           &
                                      .AND. TMPMAP(IY  ,IX+1).EQ.JJ )  &
                              TMPSTA(IY,IX) = NSTAT
                            END IF
                          IF (IY.LT.NY) THEN
                              IF (TMPSTA(IY+1,IX  ).EQ.NSTAT           &
                                      .AND. TMPMAP(IY+1,IX  ).EQ.JJ )  &
                              TMPSTA(IY,IX) = NSTAT
                            END IF
                          IF (IY.GT.1) THEN
                              IF (TMPSTA(IY-1,IX  ).EQ.NSTAT           &
                                      .AND. TMPMAP(IY-1,IX  ).EQ.JJ )  &
                               TMPSTA(IY,IX) = NSTAT
                            END IF
                          IF (TMPSTA(IY,IX).EQ.NSTAT) NBT = NBT + 1
                        END IF
                      END DO
                    END DO
!
                IF ( NBT .NE. 0 ) THEN
                    J      = J + 1
                  ELSE
                    EXIT
                  END IF
                END DO
!
              END DO
!
! ... Outer boundary excluded points
!
            IF ( GTYPE.NE.UNGTYPE ) THEN
 
                DO IX=1, NX
                    IF ( TMPSTA( 1,IX) .EQ. 1 ) TMPSTA( 1,IX) = NSTAT
                    IF ( TMPSTA(NY,IX) .EQ. 1 ) TMPSTA(NY,IX) = NSTAT
                  END DO
!
                IF ( ICLOSE.EQ.ICLOSE_NONE ) THEN
                    DO IY=2, NY-1
                        IF ( TMPSTA(IY, 1) .EQ. 1 ) TMPSTA(IY, 1) = NSTAT
                        IF ( TMPSTA(IY,NX) .EQ. 1 ) TMPSTA(IY,NX) = NSTAT
                      END DO
                  END IF
 
              END IF ! GTYPE
!
            END IF ! ILOOP .EQ. 2
!
! ... Branch back input / excluded points ( ILOOP in 8.b )
!
          END DO
!
        ELSE ! FROM .EQ. PART
!
! 8.d Read the map from file instead
!
          NSTAT  = -1
          IF (IDLA.LT.1 .OR. IDLA.GT.4) IDLA   = 1
          IF (IDFT.LT.1 .OR. IDFT.GT.3) IDFT   = 1
 
!!Li  Suspended for SMC grid though the file input line in  ww3_grid.inp
!!Li  is kept to divert the program into this block.  JGLi15Oct2014
!!Li
          IF( RGLGRD ) THEN
!!Li
!
          WRITE (NDSO,978) NDSTR, IDLA, IDFT
          IF (IDFT.EQ.2) WRITE (NDSO,973) RFORM
          IF (FROM.EQ.'NAME') WRITE (NDSO,974) TNAME
!
          IF ( NDSTR .NE. NDSI ) THEN
              IF ( IDFT .EQ. 3 ) THEN
                  WRITE (NDSE,1004) NDSTR
                  CALL EXTCDE (23)
                ELSE
                  CALL NEXTLN ( COMSTR , NDSI , NDSE )
                END IF
              IF ( IDFT .EQ. 3 ) THEN
                  IF (FROM.EQ.'NAME') THEN
                      J      = LEN_TRIM(FNMPRE)
                      OPEN (NDSTR,FILE=FNMPRE(:J)//TNAME,             &
                            FORM='UNFORMATTED',STATUS='OLD',ERR=2000, &
                            IOSTAT=IERR)
                    ELSE
                      OPEN (NDSTR,           FORM='UNFORMATTED',      &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                ELSE
                  IF (FROM.EQ.'NAME') THEN
                      J      = LEN_TRIM(FNMPRE)
                      OPEN (NDSTR,FILE=FNMPRE(:J)//TNAME,             &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    ELSE
                      OPEN (NDSTR,                                    &
                            STATUS='OLD',ERR=2000,IOSTAT=IERR)
                    END IF
                END IF
            END IF
!
          IF ( NDSTR .EQ. NDSI ) CALL NEXTLN ( COMSTR , NDSI , NDSE )
          ALLOCATE ( READMP(NX,NY) )
          CALL INA2I ( READMP, NX, NY, 1, NX, 1, NY, NDSTR, NDST,    &
                       NDSE, IDFM, RFORM, IDLA, 1, 0 )
!
          IF ( ICLOSE.EQ.ICLOSE_NONE ) THEN
              DO IY=2, NY-1
                IF ( READMP( 1,IY) .EQ. 1 ) READMP( 1,IY) = 3
                IF ( READMP(NX,IY) .EQ. 1 ) READMP(NX,IY) = 3
                END DO
            END IF
!
          DO IX=1, NX
            IF ( READMP(IX, 1) .EQ. 1 ) READMP(IX, 1) = 3
            IF ( READMP(IX,NY) .EQ. 1 .AND. ICLOSE .NE. ICLOSE_TRPL)   &
               READMP(IX,NY) = 3
            END DO
!
          DO IY=1, NY
            DO IX=1, NX
              IF ( READMP(IX,IY) .EQ. 3 ) THEN
                  TMPSTA(IY,IX) = NSTAT
                ELSE
                  TMPSTA(IY,IX) = READMP(IX,IY)
                END IF
              END DO
            END DO
          DEALLOCATE ( READMP )
!!Li
          ENDIF   !! RGLGRD
!!Li
!
        END IF !FROM .NE. 'PART'
!
! 8.e Get NSEA and other counters
!
      NSEA   = 0
      NLAND  = 0
      NBI    = 0
      NBT    = 0
!
      DO IX=1, NX
        DO IY=1, NY
          IF ( TMPSTA(IY,IX) .GT. 0 ) NSEA   = NSEA + 1
          IF ( TMPSTA(IY,IX) .EQ. 0 ) NLAND  = NLAND + 1
          IF ( TMPSTA(IY,IX) .LT. 0 ) NBT    = NBT + 1
          IF ( TMPSTA(IY,IX) .EQ. 2 ) NBI    = NBI + 1
          END DO
        END DO
!
      WRITE (NDSO,980)
      FLBPI  = NBI .GT. 0
      IF ( .NOT. FLBPI ) THEN
          WRITE (NDSO,985)
        ELSE
          WRITE (NDSO,986) NBI
          IF ( FLAGLL ) THEN
              WRITE (NDSO, 987)
            ELSE
              WRITE (NDSO,1987)
            END IF
          IBI    = 1
          DO IY=1, NY
            DO IX=1, NX
              IF (GTYPE.NE.UNGTYPE) THEN
                X = FACTOR * ( XGRDIN(IX,IY) )
                Y = FACTOR * ( YGRDIN(IX,IY) )
              ELSE
                X = FACTOR * XYB(IX,1)
                Y = FACTOR * XYB(IX,2)
                END IF
            IF ( TMPSTA(IY,IX).EQ.2 ) THEN
                  IF ( FLAGLL ) THEN
                      WRITE (NDSO, 988) IBI, IX, IY, X, Y
                    ELSE
                      WRITE (NDSO,1988) IBI, IX, IY, X, Y
                    END IF
                  IBI    = IBI + 1
                END IF
              END DO
            END DO
        END IF
!
      WRITE (NDSO,1980)
      IF ( NBT .EQ. 0 ) THEN
          WRITE (NDSO,1985)
        ELSE
          WRITE (NDSO,1986) NBT
        END IF
!
! 8.f Set up all maps
!
!!Li  CALL W3DIMX ( 1, NX, NY, NSEA, NDSE, NDST )
      CALL W3DIMX ( 1, NX, NY, NSEA, NDSE, NDST  &
                  )
!
! 8.g Activation of reflections and scattering
      FFACBERG=FACBERG
 
 
      IF (GTYPE.NE.UNGTYPE) THEN
        DO IY=1, NY
          DO IX=1, NX
            XGRD(IY,IX) = XGRDIN(IX,IY)
            YGRD(IY,IX) = YGRDIN(IX,IY)
            END DO
          END DO
          DEALLOCATE ( XGRDIN, YGRDIN )
          CALL W3GNTX ( 1, 6, 6 )
      ELSE
!
!FA:  This distinction  between structured and unstructured
! should be removed when XYB is replaced by XGRD and YGRD
!
        DO IX=1, NX
          XGRD(:,IX) = XYB(IX,1)
          YGRD(:,IX) = XYB(IX,2)
          END DO
        END IF   ! GTYPE
!
!!Li  MAPSTA = TMPSTA
!!Li  Shelter MAPSTA LLG definition for SMC by RGLGRD.
      IF( RGLGRD ) MAPSTA = TMPSTA
      MAPFS  = 0
!
      TRNX   = 0.
      TRNY   = 0.
!
!Li  Shelter MAPSTA etc LLG definitions for SMC by logical RGLGRD
      IF( RGLGRD ) THEN
      ISEA   = 0
      DO IY=1, NY
        DO IX=1, NX
          IF ( TMPSTA(IY,IX) .EQ. NSTAT ) THEN
              MAPSTA(IY,IX) = 0
              MAPST2(IY,IX) = 1
              TMPSTA(IY,IX) = 3
            ELSE
              MAPSTA(IY,IX) = TMPSTA(IY,IX)
              MAPST2(IY,IX) = 0
            END IF
          IF ( MAPSTA(IY,IX) .NE. 0 ) THEN
              ISEA           = ISEA + 1
              MAPFS (IY,IX)  = ISEA
              ZB(ISEA)       = ZBIN(IX,IY)
              MAPSF(ISEA,1)  = IX
              MAPSF(ISEA,2)  = IY
              IF ( FLAGLL ) THEN
                  Y              = YGRD(IY,IX)
                  CLATS(ISEA)    = COS(Y*DERA)
                  CLATIS(ISEA)   = 1. / CLATS(ISEA)
                  CTHG0S(ISEA)   = - TAN(DERA*Y) / RADIUS
                ELSE
                  CLATS(ISEA)    = 1.
                  CLATIS(ISEA)   = 1.
                  CTHG0S(ISEA)   = 0.
                END IF
            END IF
 
!/ ------------------------------------------------------------------- /
 
! notes: Oct 22 2012: I moved the following "if-then" statement from
! inside the  "IF ( MAPSTA(IY,IX) .NE. 0 )" statement to outside that
! statement. This is needed since later on, ATRNX is computed from
! TRNX(ix-1) , TRNX(ix) etc. which causes boundary effects if the
! MAPSTA=0 values are set to TRNX=0
 
              IF ( TRFLAG .NE. 0 ) THEN
                  TRNX(IY,IX) = 1. - OBSX(IX,IY)
                  TRNY(IY,IX) = 1. - OBSY(IX,IY)
                END IF
 
          END DO
        END DO
      ENDIF
!!Li End of RGLGRD IF block
!
      DO ISP=1, NSPEC+NTH
        MAPWN(ISP) = 1 + (ISP-1)/NTH
        MAPTH(ISP) = 1 + MOD(ISP-1,NTH)
        END DO
!
      NMAP   = 1 + (NX-1)/NCOL
      WRITE (NDSO,1100) NMAP
      DO IMAP=1, NMAP
        IX0    = 1 + (IMAP-1)*NCOL
        IXN    = MIN ( NX , IMAP*NCOL )
        DO IY=NY,1,-1
          WRITE (NDSO,1101) (TMPSTA(IY,IX),IX=IX0,IXN)
          END DO
        WRITE (NDSO,*) ' '
        END DO
      WRITE (NDSO,1102)
 
!
 
!
! 9.d Estimates shoreline direction for reflection
!     and shoreline treatment in general for UNST grids.
! NB: this is updated with moving water levels in W3ULEV
!
      IF (GTYPE.EQ.UNGTYPE) THEN
        CALL SETUGIOBP
        END IF
!
      DEALLOCATE ( ZBIN, TMPSTA, TMPMAP )
!
! 9.e Reads bottom information from file
!
 
 
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! 10.  Prepare output boundary points.
!     ILOOP = 1 to count NFBPO and NBO
!     ILOOP = 2 to fill data arrays
!
      WRITE (NDSO,990)
      J      = LEN_TRIM(FNMPRE)
      OPEN (NDSS,FILE=FNMPRE(:J)//'ww3_grid.scratch',FORM='FORMATTED')
!
      DO ILOOP = 1, 2
!
        IF ( ILOOP.EQ.2 ) CALL W3DMO5 ( 1, NDST, NDSE, 2 )
!
        NBOTOT = 0
        NFBPO  = 0
        NBO(0) = 0
        NBO2(0)= 0
        FIRST  = .TRUE.
        REWIND (NDSS)
        IF ( ILOOP .EQ. 1 ) THEN
            NDSI2 = NDSI
          ELSE
            NDSI2 = NDSS
          END IF
!
        DO
          CALL NEXTLN ( COMSTR , NDSI2 , NDSE )
          READ (NDSI2,*,END=2001,ERR=2002) XO0, YO0, DXO, DYO, NPO
!
          IF ( ILOOP .EQ. 1 ) THEN
              BACKSPACE (NDSI)
              READ (NDSI,'(A)') LINE
              WRITE (NDSS,'(A)') LINE
            END IF
!
! ... Check if new file to be used
!
          FIRST  = FIRST .OR. NPO.LE.0
          NPO    = ABS(NPO)
!
! ... Preparations for new output file including end check
!     and output for last output file
!
          IF ( FIRST ) THEN
!
              FIRST  = .FALSE.
!
              IF ( NFBPO.GE.1 .AND. ILOOP.EQ.2 ) THEN
                  WRITE (NDSO,991)  NFBPO, NBO(NFBPO) - NBO(NFBPO-1), &
                                          NBO2(NFBPO) - NBO2(NFBPO-1)
                  IF ( NBO(NFBPO) - NBO(NFBPO-1) .EQ. 1 ) THEN
                      IF ( FLAGLL ) THEN
                          WRITE (NDSO,992)
                        ELSE
                          WRITE (NDSO,2992)
                        END IF
                    ELSE
                      IF ( FLAGLL ) THEN
                          WRITE (NDSO,1992)
                        ELSE
                          WRITE (NDSO,3992)
                        END IF
                    END IF
                  IP0    = NBO(NFBPO-1)+1
                  IPN    = NBO(NFBPO)
                  IPH    = IP0 + (IPN-IP0-1)/2
                  IPI    = IPH -IP0 + 1 + MOD(IPN-IP0+1,2)
                  DO IP=IP0, IPH
                    IF ( FLAGLL ) THEN
                        WRITE (NDSO,1993) IP-NBO(NFBPO-1),     &
                                          FACTOR*XBPO(IP),     &
                                          FACTOR*YBPO(IP),     &
                                          IP+IPI-NBO(NFBPO-1), &
                                          FACTOR*XBPO(IP+IPI), &
                                          FACTOR*YBPO(IP+IPI)
                      ELSE
                        WRITE (NDSO,3993) IP-NBO(NFBPO-1),     &
                                          FACTOR*XBPO(IP),     &
                                          FACTOR*YBPO(IP),     &
                                          IP+IPI-NBO(NFBPO-1), &
                                          FACTOR*XBPO(IP+IPI), &
                                          FACTOR*YBPO(IP+IPI)
                      END IF
                    END DO
                  IF ( MOD(IPN-IP0+1,2) .EQ. 1 ) THEN
                      IF ( FLAGLL ) THEN
                          WRITE (NDSO, 993) IPH+1-NBO(NFBPO-1), &
                                            FACTOR*XBPO(IPH+1), &
                                            FACTOR*YBPO(IPH+1)
                        ELSE
                          WRITE (NDSO,2993) IPH+1-NBO(NFBPO-1), &
                                            FACTOR*XBPO(IPH+1), &
                                            FACTOR*YBPO(IPH+1)
                        END IF
                    END IF
                  WRITE (NDSO,*)
                END IF
!
              IF ( NPO .EQ. 0 ) EXIT
!
              NFBPO  = NFBPO + 1
              IF ( NFBPO .GT. 9 ) THEN
                  WRITE (NDSE,1006)
                  CALL EXTCDE ( 50 )
                END IF
              NBO2(NFBPO) = NBO2(NFBPO-1)
              NBO(NFBPO) = NBOTOT
!
            END IF
!
! ... Loop over line segment - - - - - - - - - - - - - - - - - - - - -
!
          DO IP=1, NPO
!
            XO     = XO0 + REAL(IP-1)*DXO
            YO     = YO0 + REAL(IP-1)*DYO
!
! ... Compute bilinear remapping weights
!
            INGRID = W3GRMP( GSU, XO, YO, IXR, IYR, RD )
!
!           Change cell-corners from counter-clockwise to column-major order
            IX     = IXR(3);  IY     = IYR(3);  X     = RD(3);
            IXR(3) = IXR(4);  IYR(3) = IYR(4);  RD(3) = RD(4);
            IXR(4) = IX    ;  IYR(4) = IY    ;  RD(4) = X    ;
!
! ... Check if point in grid
!
            IF ( INGRID ) THEN
!
! ... Check if point not on land
!
                IF ( ( MAPSTA(IYR(1),IXR(1)).GT.0 .AND.               &
                                          RD(1).GT.0.05 ) .OR.        &
                     ( MAPSTA(IYR(2),IXR(2)).GT.0 .AND.               &
                                          RD(2).GT.0.05 ) .OR.        &
                     ( MAPSTA(IYR(3),IXR(3)).GT.0 .AND.               &
                                          RD(3).GT.0.05 ) .OR.        &
                     ( MAPSTA(IYR(4),IXR(4)).GT.0 .AND.               &
                                          RD(4).GT.0.05 ) ) THEN
!
! ... Check storage and store coordinates
!
                    NBOTOT = NBOTOT + 1
                    IF ( ILOOP .EQ. 1 ) CYCLE
!
                    XBPO(NBOTOT) = XO
                    YBPO(NBOTOT) = YO
!
! ... Interpolation factors
!
                RDTOT = 0.
                DO J=1, 4
                  IF ( MAPSTA(IYR(J),IXR(J)).GT.0 .AND.               &
                                            RD(J).GT.0.05 ) THEN
                      RDBPO(NBOTOT,J) = RD(J)
                    ELSE
                      RDBPO(NBOTOT,J) = 0.
                    END IF
                    RDTOT = RDTOT + RDBPO(NBOTOT,J)
                  END DO
!
              DO J=1, 4
                RDBPO(NBOTOT,J) = RDBPO(NBOTOT,J) / RDTOT
                END DO
!
! ... Determine sea and interpolation point counters
!
                DO J=1, 4
                  ISEAI(J) = MAPFS(IYR(J),IXR(J))
                  END DO
!
                DO J=1, 4
                  IF ( ISEAI(J).EQ.0 .OR. RDBPO(NBOTOT,J).EQ. 0. ) THEN
                      IPBPO(NBOTOT,J) = 0
                    ELSE
                      FLNEW   = .TRUE.
                      DO IST=NBO2(NFBPO-1)+1, NBO2(NFBPO)
                        IF ( ISEAI(J) .EQ. ISBPO(IST) ) THEN
                            FLNEW  = .FALSE.
                            IPBPO(NBOTOT,J) = IST - NBO2(NFBPO-1)
                          END IF
                        END DO
                      IF ( FLNEW ) THEN
                          NBO2(NFBPO)        = NBO2(NFBPO) + 1
                          IPBPO(NBOTOT,J)    = NBO2(NFBPO) - NBO2(NFBPO-1)
                          ISBPO(NBO2(NFBPO)) = ISEAI(J)
                        END IF
                    END IF
                  END DO
!
! ... Error output
!
                  ELSE
                    WRITE (NDSE,995) FACTOR*XO, FACTOR*YO
                  END IF
              ELSE
                WRITE (NDSE,994) FACTOR*XO, FACTOR*YO
              END IF
!
            END DO
!
          NBO(NFBPO) = NBOTOT
!
! ... Branch back to read.
!
          END DO
!
! ... End of ILOOP loop
!
        END DO
!
      CLOSE ( NDSS, STATUS='DELETE' )
!
      FLBPO  = NBOTOT .GT. 0
      IF ( .NOT. FLBPO ) THEN
          WRITE (NDSO,996)
        ELSE
          WRITE (NDSO,997) NBOTOT, NBO2(NFBPO)
        END IF
!
!--- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!10.  Write model definition file.
!
      WRITE (NDSO,999)
      CALL W3IOGR ( 'WRITE', NDSM )
!
      CLOSE (NDSM)
!
      GOTO 2222
!
! Escape locations read errors :
!
 2000 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 60 )
!
 2001 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 61 )
!
 2002 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 62 )
!
 2222 CONTINUE
      IF ( GTYPE .NE. UNGTYPE) THEN
          IF ( NX*NY .NE. NSEA ) THEN
              WRITE (NDSO,9997) NX, NY, NX*NY, NSEA,                       &
                                100.*REAL(NSEA)/REAL(NX*NY), NBI, NLAND, NBT
            ELSE
              WRITE (NDSO,9998) NX, NY, NX*NY, NSEA, NBI, NLAND, NBT
            END IF
        ELSE IF ( GTYPE .EQ. UNGTYPE ) THEN
          IF ( NX*NY .NE. NSEA ) THEN
              WRITE (NDSO,9997)  0,  0, NX*NY, NSEA,                       &
                                100.*REAL(NSEA)/REAL(NX*NY), NBI, NLAND, NBT
            ELSE
              WRITE (NDSO,9998)  0,  0, NX*NY, NSEA, NBI, NLAND, NBT
            END IF
        ENDIF ! GTYPE .EQ. UNGTYPE
 
      WRITE (NDSO,9999)
 
!
! Formats
!
  900 FORMAT (/15X,'    *** WAVEWATCH III Grid preprocessor ***    '/ &
               15X,'==============================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
  902 FORMAT ( '  Grid name : ',A/)
  903 FORMAT (/'  Spectral discretization : '/                        &
               ' --------------------------------------------------'/ &
               '       Number of directions        :',I4/             &
               '       Directional increment (deg.):',F6.1)
  904 FORMAT ( '       First direction       (deg.):',F6.1)
  905 FORMAT ( '       Number of frequencies       :',I4/             &
               '       Frequency range        (Hz) :',F9.4,'-',F6.4/  &
               '       Increment factor            :',F8.3/)
!
  910 FORMAT (/'  Model definition :'/                                &
               ' --------------------------------------------------')
  911 FORMAT ( '       Dry run (no calculations)   :  ',A/            &
               '       Propagation in X-direction  :  ',A/            &
               '       Propagation in Y-direction  :  ',A/            &
               '       Refraction                  :  ',A/            &
               '       Current-induced k-shift     :  ',A/            &
               '       Source term calc. and int.  :  ',A/)
  912 FORMAT (/'  Time steps : '/                                     &
               ' --------------------------------------------------'/ &
               '       Maximum global time step      (s) :',F8.2/     &
               '       Maximum CFL time step X-Y     (s) :',F8.2/     &
               '       Maximum CFL time step k-theta (s) :',F8.2/     &
               '       Minimum source term time step (s) :',F8.2/)
  913 FORMAT (/ '  WARNING, TIME STEP LESS THAN 1 s, NITER:',I8 /)
  915 FORMAT ( '  Preprocessing namelists ...')
  916 FORMAT ( '  Preprocessing namelists finished.'/)
  917 FORMAT (/'  Equivalent namelists ...'/)
  918 FORMAT (/'  Equivalent namelists finished.'/)
!
  820 FORMAT (/'  Linear input (C&M-R 82) ',A/                   &
        ' --------------------------------------------------')
  821 FORMAT ( '       CLIN                        :',f8.2/      &
               '       Factor for fPM in filter    :',F8.2/      &
               '       Factor for fh in filter     :',F8.2/)
 2820 FORMAT ( '  &SLN1 CLIN =',F6.1,', RFPM =',F6.2,            &
               ', RFHF =',F6.2,' /')
!
  920 FORMAT (/'  Wind input (WAM 4+) ',A/                            &
        ' --------------------------------------------------')
  921 FORMAT ( '       minimum Charnock coeff.     :',F10.4/          &
               '       betamax                     :',F9.3/           &
               '       power of cos. in wind input :',F9.3/           &
               '       z0max                       :',F9.3/           &
               '       zalp                        :',F9.3/           &
               '       Height of input wind (m)    :',F8.2/           &
               '       wind stress sheltering      :',F9.3/           &
               '       swell attenuation param.    :',I5/             &
               '       swell attenuation factor    :',F9.3/           &
               '       swell attenuation factor2   :',F9.3/           &
               '       swell attenuation factor3   :',F9.3/           &
               '       critical Reynolds number    :',F9.1/           &
               '       swell attenuation factor5   :',F9.3/           &
               '       swell attenuation factor6   :',F9.3/           &
               '       swell attenuation factor7   :',F14.3/          &
               '       ratio of z0 for orb. & mean :',F9.3/)
 2920 FORMAT ( '  &SIN4 ZWND =',F5.1,', ALPHA0 =',F8.5,', Z0MAX =',F8.5,', BETAMAX =', &
                  F8.5,','/                                           &
              '        SINTHP =',F8.5,', ZALP =',F8.5,', TAUWSHELTER =',F8.5,           &
              ', SWELLFPAR =',I2,','/                                 &
              '        SWELLF =',F8.5,', SWELLF2 =',F8.5,             &
              ', SWELLF3 =',F8.5,', SWELLF4 =',F9.1,','/              &
              '        SWELLF5 =',F8.5,', SWELLF6 =',F8.5,            &
              ', SWELLF7 =',F12.2,', Z0RAT =',F8.5,', SINBR =',F8.5,'  /')
!
  922 FORMAT (/'  Nonlinear interactions (DIA) ',A/                   &
               ' --------------------------------------------------')
  923 FORMAT ( '       Lambda                      :',F8.2/      &
               '       Prop. constant              :',E10.3/     &
               '       kd conversion factor        :',F8.2/      &
               '       minimum kd                  :',F8.2/      &
               '       shallow water constants     :',F8.2,2F6.2/)
 2922 FORMAT ( '  &SNL1 LAMBDA =',F7.3,', NLPROP =',E10.3,       &
               ', KDCONV =',F7.3,', KDMIN =',F7.3,','/           &
               '        SNLCS1 =',F7.3,', SNLCS2 =',F7.3,        &
               ', SNLCS3 = ',F7.3,' /')
!
  924 FORMAT (/' Dissipation (Ardhuin et al. 2010) ',A/          &
        ' --------------------------------------------------')
  925 FORMAT ( '       SDSC2, SDSBCK, SDSCUM       :',3E11.3/    &
               '       Power of k in mean k        :',F8.2/)
 2924 FORMAT ( '  &SDS4 SDSC1 =',E12.4,', SDSC2 =',E12.4,        &
               ', SDSCUM =',F6.2,', SDSSTRAIN =',F4.1','/        &
               '        SDSC4 =',F6.2,', SDSC5 =',E12.4,         &
               ', SDSC6 =',E12.4,','/                            &
               '        WNMEANP =',F4.2,', FXPM3 =', F4.2,       &
               ', FXFM3 =',F4.2,', FXFMAGE =',F6.3,              &
               ', FXINCUT =',F6.3,', FXDSCUT =',F6.3,', '/       &
               '        SDSBINT =',E12.4,', SDSBCK =',E12.4,     &
               ', SDSABK =',F6.3,', SDSPBK =',F6.3,', '/         &
               '        SDSHCK =',F5.2,', SDSBR = ',E12.4,       &
               ', SDSSTRAIN =',F6.3,', '/                        &
               '        SDSBR2 =',F5.2,', SDSP =',F5.2,          &
               ', SDSISO =',I2, &
               ', SDSCOS =',F3.1,', SDSDTH =',F5.1,', '/         &
               '        SDSBRF1 = ',F5.2,', SDSBRFDF =',I2,', '/ &
               '        SDSBM0 = ',F5.2, ', SDSBM1 =',F5.2,      &
               ', SDSBM2 =',F5.2,', SDSBM3 =',F5.2,', SDSBM4 =', &
               F5.2,', '/,                                       &
               ',       WHITECAPWIDTH =',F5.2,', SDSLFGEN = ',   &
               F5.2,', SDSHFGEN = ',F5.2,' /')
!
  926 FORMAT (/'  Bottom friction (JONSWAP) ',A/                 &
               ' --------------------------------------------------')
  927 FORMAT ( '       gamma                       :',F8.4/)
 2926 FORMAT ( '  &SBT1 GAMMA =',E12.4,' /')
!
  928 FORMAT (/'  Surf breaking (B&J 1978) ',A/                  &
               ' --------------------------------------------------')
  929 FORMAT ( '       alpha                       :',F8.3/      &
               '       gamma                       :',F8.3)
 2928 FORMAT ( '  &SDB1 BJALFA =',F7.3,', BJGAM =',F7.3,         &
               ', BJFLAG = ',A,' /')
!
  930 FORMAT (/'  Triad interactions not defined.'/)
!
  932 FORMAT (/'  Bottom scattering not defined.'/)
!
  934 FORMAT (/'  Alternative source term slot not used.'/)
!
 944 FORMAT  (/'  Ice scattering not defined.'/)
!
  950 FORMAT (/'  Propagation scheme : '/                             &
               ' --------------------------------------------------')
  951 FORMAT ( '       Type of scheme (structured) :',1X,A)
 2951 FORMAT ( '       Type of scheme(unstructured):',1X,A)
  952 FORMAT ( '                                    ',1X,A)
!
  953 FORMAT ( '       CFLmax depth refraction     :',F9.3/      &
               '       Averaging area factor Cg    :',F8.2)
  954 FORMAT ( '       Averaging area factor theta :',F8.2)
  955 FORMAT ( '            **** Internal maximum .GE.',F6.2,' ****')
 2953 FORMAT ( '  &PRO3 CFLTM =',F5.2,                           &
                      ', WDTHCG = ',F4.2,', WDTHTH = ',F4.2,' /')
!
 2956 FORMAT ( '  &UNST UGOBCAUTO =',L3, ', UGOBCDEPTH =', F8.3/      &
               ',       EXPFSN =',L3,',EXPFSPSI =',L3, ',EXPFSFCT =', &
               L3,',IMPFSN =',L3         / )
!
  960 FORMAT (/'  Miscellaneous ',A/                                   &
               ' --------------------------------------------------')
 2961 FORMAT ( ' *** WAVEWATCH-III WARNING IN W3GRID :'/               &
               '     CICE0.NE.CICEN requires FLAGTR>2'/                &
               '     Parameters corrected: CICE0 = CICEN'/)
 2962 FORMAT (/' *** WAVEWATCH-III WARNING IN W3GRID : User requests', &
         'CICE0=CICEN corresponding to discontinuous treatment of ',   &
         'ice, so we will change FLAGTR')
 2963 FORMAT (/' *** WAVEWATCH-III WARNING IN W3GRID :'/               &
               '     Ice physics used, so we will change FLAGTR.')
  961 FORMAT ( '       Ice concentration cut-offs  :',F8.2,F6.2)
  965 FORMAT (/'    Dynamic source term integration scheme :'/        &
               '       Xp                      (-) :',F9.3/           &
               '       Xr                      (-) :',F9.3/           &
               '       Xfilt                   (-) :',F9.3)
  966 FORMAT (/'    Wave field partitioning :'/                       &
               '       Levels                  (-) :',I5/             &
               '       Minimum wave height     (m) :',F9.3/           &
               '       Wind area multiplier    (-) :',F9.3/           &
               '       Cut-off wind sea fract. (-) :',F9.3/           &
               '       Combine wind seas           :  ',A/            &
               '       Number of swells in fld out :',I5)
  967 FORMAT (/'    Miche-style limiting wave height :'/              &
               '       Hs,max/d factor         (-) :',F9.3/           &
               '       Hrms,max/d factor       (-) :',F9.3/           &
               '       Limiter activated           :  ',A)
  968 FORMAT ( '          *** FACTOR DANGEROUSLY LOW ***')
!
 8972 FORMAT ( '       Wind input reduction factor in presence of ', &
               /'         ice :',F6.2, &
               /'         (0.0==> no reduction and 1.0==> no wind', &
               /'         input with 100% ice cover)')
!
 4970 FORMAT (/'  Spectral output on full grid ',A/                   &
               ' --------------------------------------------------')
 4971 FORMAT ( '       Second order pressure at K=0:',3I4)
 4972 FORMAT ( '       Spectrum of Uss             :',3I4)
 4973 FORMAT ( '       Frequency spectrum          :',3I4)
!
 4980 FORMAT (/'  Coastal / iceberg reflection  ',A/                   &
               ' --------------------------------------------------')
 4981 FORMAT ( '       Coefficient for shorelines  :',F6.4)
 4989 FORMAT ( '          *** CURVLINEAR GRID: REFLECTION NOT IMPLEMENTED YET ***')
 2977 FORMAT ( '  &SIG1  IGMETHOD =',I2,', IGADDOUTP =',I2,', IGSOURCE =',I2, &
               ', IGSTERMS = ',I2,', IGBCOVERWRITE =', L3,','/        &
               '        IGSWELLMAX =', L3,', IGMAXFREQ =',F6.4,       &
               ', IGSOURCEATBP = ',I2,', IGKDMIN = ',F6.4,','/        &
               '        IGFIXEDDEPTH = ',F6.2', IGEMPIRICAL = ',F8.6,' /')
!
 2978 FORMAT ( '  &SIC2  IC2DISPER =',L3,', IC2TURB =',F8.2,          &
               ', IC2ROUGH  =',F10.6,','/                             &
               '        IC2REYNOLDS = ',F10.1,', IC2SMOOTH = ',F10.1, &
               ', IC2VISC =',F10.3,', IC2TURBS =',F8.2,' /')
!
 2979 FORMAT ( '  &SIC3 IC3MAXTHK =',F6.2, ', IC3MAXCNC =',F6.2,','/  &
               '        IC2TURB =',F8.2,                              &
               ', IC2ROUGH  =',F7.3,','/                              &
               '        IC2REYNOLDS = ',F10.1,', IC2SMOOTH = ',F10.1, &
               ', IC2VISC =',F10.3,','/                               &
               '        IC2TURBS =',F8.2,', IC3CHENG =',L3,           &
               ', USECGICE =',L3,', IC3HILIM = ',F6.2,','/            &
               '        IC3KILIM = ',E9.2,', IC3HICE = ',E9.2,        &
               ', IC3VISC = ',E9.2,','/                               &
               '        IC3DENS = ',E9.2,', IC3ELAS = ',E9.2,' /')
!
 2966 FORMAT ( '  &MISC CICE0 =',F6.3,', CICEN =',F6.3,               &
                     ', LICE = ',F8.1,', PMOVE =',F6.3,','/           &
               '        XSEED =',F6.3,', FLAGTR = ', I1,              &
                     ', XP =',F6.3,', XR =',F6.3,', XFILT =', F6.3 /  &
               '        IHM =',I5,', HSPM =',F6.3,', WSM =',F6.3,     &
                     ', WSC =',F6.3,', FLC = ',A/                     &
               '        NOSW =',I3,', FMICHE =',F6.3,', RWNDC =' ,    &
                        F6.3,', WCOR1 =',F6.2,', WCOR2 =',F6.2,','/   &
               '        FACBERG =',F4.1,', GSHIFT = ',E11.3,','/      &
               '        STDX =',F7.2,', STDY =',F7.2,', STDT =', F8.2,&
                     ', ICEHMIN = ',F5.2,', ICEHINIT = ',F5.2,','/   &
               '        ICEDISP = ',L3,', ICEWIND = ',F6.2,' /')
!
 2976 FORMAT ( '  &OUTS P2SF  =',I2,', I1P2SF =',I2,', I2P2SF =',I3,','/&
               '        US3D  =',I2,', I1US3D =',I3,', I2US3D =',I3,','/&
               '        E3D   =',I2,', I1E3D  =',I3,', I2E3D  =',I3,','/&
               '        TH1MF =',I2,', I1TH1M =',I3,', I2TH1M =',I3,','/&
               '        STH1MF=',I2,', I1STH1M=',I3,', I2STH1M=',I3,','/&
               '        TH2MF =',I2,', I1TH2M =',I3,', I2TH2M =',I3,','/&
               '        STH2MF=',I2,', I1STH2M=',I3,', I2STH2M=',I3,' /')
!
 2986 FORMAT ( '  &REF1 REFCOAST =',F5.2,', REFFREQ =',F5.2,', REFSLOPE =',F5.3, &
               ', REFMAP =',F4.1, ', REFMAPD =',F4.1, ', REFSUBGRID =',F5.2,','/ &
               '        REFRMAX=',F5.2,', REFFREQPOW =',F5.2,                    &
               ', REFICEBERG =',F5.2,', REFCOSP_STRAIGHT =',F4.1,' /')
!
 2987 FORMAT ( '  &FLD TAIL_ID =',I1,' TAIL_LEV =',F5.4,' TAILT1 =',F5.3,&
               ' TAILT2 =',F5.3,' /')
 
 3000 FORMAT (/'  The spatial grid: '/                                &
               ' --------------------------------------------------'/ &
              /'       Grid type                   : ',A)
 3001 FORMAT ( '       Coordinate system           : ',A)
 3002 FORMAT ( '       Index closure type          : ',A)
 3003 FORMAT ( '       Dimensions                  : ',I6,I8)
 3004 FORMAT (/'       Increments           (deg.) :',2F10.4/         &
               '       Longitude range      (deg.) :',2F10.4/         &
               '       Latitude range       (deg.) :',2F10.4)
 3005 FORMAT ( '       Increments             (km) :',2F8.2/          &
               '       X range                (km) :',2F8.2/          &
               '       Y range                (km) :',2F8.2)
 3006 FORMAT (/'       X-coordinate unit           :',I6/             &
               '       Scale factor                :',F10.4/           &
               '       Add offset                  :',E12.4/          &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
 3007 FORMAT (/'       Y-coordinate unit           :',I6/             &
               '       Scale factor                :',F10.4/           &
               '       Add offset                  :',E12.4/          &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
 3008 FORMAT ( '       Format                      : ',A)
 3009 FORMAT ( '       File name                   : ',A)
  972 FORMAT (/'       Bottom level unit           :',I6/             &
               '       Limiting depth          (m) :',F8.2/           &
               '       Minimum depth           (m) :',F8.2/           &
               '       Scale factor                :',F8.2/           &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
  973 FORMAT ( '       Format                      : ',A)
  974 FORMAT ( '       File name                   : ',A)
  976 FORMAT ( '       Sub-grid information        : ',A)
  977 FORMAT ( '       Obstructions unit           :',I6/             &
               '       Scale factor                :',F10.4/          &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
  978 FORMAT ( '       Mask information            : From file.'/     &
               '       Mask unit                   :',I6/             &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
 1977 FORMAT ( '       Shoreline slope             :',I6/             &
               '       Scale factor                :',F10.4/          &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
 1978 FORMAT ( '       Grain sizes                 :',I6/             &
               '       Scale factor                :',F10.4/          &
               '       Layout indicator            :',I6/             &
               '       Format indicator            :',I6)
!
  979 FORMAT ( '  Processing ',A)
  980 FORMAT (/'  Input boundary points : '/                          &
               ' --------------------------------------------------')
 1980 FORMAT (/'  Excluded points : '/                                &
               ' --------------------------------------------------')
  981 FORMAT ( '   *** POINT OUTSIDE GRID (SKIPPED), IX, IY =')
 1981 FORMAT ( '   *** POINT ALREADY EXCLUDED (SKIPPED), IX, IY =')
  982 FORMAT ( '   *** CANNOT CONNECT POINTS, IX, IY =')
  985 FORMAT ( '       No boundary points.'/)
  986 FORMAT ( '       Number of boundary points   :',I6/)
 1985 FORMAT ( '       No excluded points.'/)
 1986 FORMAT ( '       Number of excluded points   :',I6/)
  987 FORMAT ( '         Nr.|  IX |  IY |  Long.  |   Lat.  '/        &
               '       -----|-----|-----|---------|---------')
 1987 FORMAT ( '         Nr.|  IX |  IY |     X     |     Y     '/    &
               '       -----|-----|-----|-----------|-----------')
  988 FORMAT ( '       ',I4,2(' |',I4),2(' |',F8.2))
 1988 FORMAT ( '       ',I4,2(' |',I4),2(' |',F8.1,'E3'))
  989 FORMAT ( ' ')
!
  990 FORMAT (/'  Output boundary points : '/                         &
               ' --------------------------------------------------')
  991 FORMAT ( '       File nest',I1,'.ww3  Number of points  :',I6/  &
               '                       Number of spectra :',I6)
  992 FORMAT (/'         Nr.|  Long.  |   Lat.  '/               &
               '       -----|---------|---------')
 1992 FORMAT (/'         Nr.|  Long.  |   Lat.  ',               &
               '         Nr.|  Long.  |   Lat.  '/               &
               '       -----|---------|---------',               &
               '       -----|---------|---------')
  993 FORMAT ( '       ',I4,2(' |',F8.2))
 1993 FORMAT ( '       ',I4,2(' |',F8.2),                        &
              '        ',I4,2(' |',F8.2))
  994 FORMAT ( '   *** POINT OUTSIDE GRID (SKIPPED) : X,Y =',2F7.2)
  995 FORMAT ( '   *** POINT ON LAND      (SKIPPED) : X,Y =',2F7.2)
 2992 FORMAT (/'         Nr.|     X     |     Y     '/           &
               '       -----|-----------|-----------')
 3992 FORMAT (/'         Nr.|     X     |     Y     ',           &
                 '       Nr.|     X     |     Y     '/           &
               '       -----|-----------|-----------',           &
                 '     -----|-----------|-----------')
 2993 FORMAT ( '       ',I4,2(' |',F8.1,'E3'))
 3993 FORMAT ( '       ',I4,2(' |',F8.1,'E3'),                   &
                '      ',I4,2(' |',F8.1,'E3'))
 2994 FORMAT ( '   *** POINT OUTSIDE GRID (SKIPPED) : X,Y =',2(F8.1,'E3'))
 2995 FORMAT ( '   *** POINT ON LAND      (SKIPPED) : X,Y =',2(F8.1,'E3'))
  996 FORMAT ( '       No boundary points.'/)
  997 FORMAT ( '       Number of boundary points   :',I6/             &
               '       Number of spectra           :',I6/)
!
  999 FORMAT (/'  Writing model definition file ...'/)
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
         '     DISCREPANCY IN DECLARED AND READ SPECTRAL DIMENSIONS'/ &
               '     DECLARED NK, NTH : ',2I8/                         &
               '     READ NK, NTH     : ',2I8/)
!
 1004 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
               '     CANNOT READ UNFORMATTED (IDFM = 3) FROM UNIT',   &
               I4,' (ww3_grid.inp)'/)
!
 1005 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID : '/               &
               '     BOTTOM AND OBSTRUCTION DATA FROM SAME FILE '/    &
               '     BUT WITH INCOMPATIBLE FORMATS (',I1,',',I1,')'/)
!
 1006 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     TOO MANY NESTING OUTPUT FILES '/)
!
 1007 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'/                &
               '     ILLEGAL GRID TYPE:',A4)
!
 1008 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'/                &
               '     A CARTESIAN WITH CLOSURE IS NOT ALLOWED')
!
 1009 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'/                &
               '     A RECTILINEAR TRIPOLE GRID IS NOT ALLOWED')
!
 1010 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'//               &
               '     NO PROPAGATION + NO SOURCE TERMS = NO WAVE MODEL'// &
               '     ( USE DRY RUN FLAG TO TEMPORARILY SWITCH OFF ',  &
               'CALCULATIONS )'/)
!
 1011 FORMAT (/' *** WAVEWATCH-III WARNING IN W3GRID :'/              &
               '     LEFT-HANDED GRID -- POSSIBLE CAUSE IS WRONG '/   &
               '     IDLA:',I4,' . THIS MAY PRODUCE ERRORS '/         &
               '     (COMMENT THIS EXTCDE AT YOUR OWN RISK).')
!
 1012 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'/                &
               '     ILLEGAL GRID CLOSURE TYPE:',A4)
!
 1013 FORMAT (/' *** WAVEWATCH-III WARNING IN W3GRID :'/              &
               '     THE GLOBAL (LOGICAL) INPUT FLAG IS DEPRECATED'/  &
               '     AND REPLACED WITH A STRING INDICATING THE TYPE'/ &
               '     OF GRID INDEX CLOSURE (NONE, SMPL or TRPL).'/    &
               ' *** PLEASE UPDATE YOUR GRID INPUT FILE ACCORDINGLY ***'/)
!
 1020 FORMAT (/' *** WAVEWATCH-III ERROR IN W3GRID :'/                &
               '     SOURCE TERMS REQUESTED BUT NOT SELECTED'/)
 1021 FORMAT (/' *** WAVEWATCH III WARNING IN W3GRID :'/              &
               '     SOURCE TERMS SELECTED BUT NOT REQUESTED'/)
 1022 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/LNn OR SEED SWITCHES :',I3)
 1023 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/STn SWITCHES :',I3)
 1024 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/NLn SWITCHES :',I3)
 1025 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/BTn SWITCHES :',I3)
 1026 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/DBn SWITCHES :',I3)
 1027 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/TRn SWITCHES :',I3)
 1028 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/BSn SWITCHES :',I3)
 1029 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/XXn SWITCHES :',I3)
!
 1030 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     PROPAGATION REQUESTED BUT NO SCHEME SELECTED '/)
 1031 FORMAT (/' *** WAVEWATCH III WARNING IN W3GRID :'/              &
               '     NO PROPAGATION REQUESTED BUT SCHEME SELECTED '/)
 1032 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     NO PROPAGATION SCHEME SELECTED ( use !/PR0 ) '/)
 1033 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     MULTIPLE PROPAGATION SCHEMES SELECTED :',I3/     &
               '     CHECK !/PRn SWITCHES'/)
 1034 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/ICn SWITCHES :',I3)
 1035 FORMAT (/' *** WAVEWATCH III WARNING IN W3GRID :'/              &
               '     ONLY FIRST PROPAGATION SCHEME WILL BE USED: ')
 1036 FORMAT (/' *** WAVEWATCH III ERROR IN W3GRID :'/                &
               '     ILLEGAL NUMBER OF !/ISn SWITCHES :',I3)
!
 1040 FORMAT ( '       Space-time extremes DX      :',F10.2)
 1041 FORMAT ( '       Space-time extremes DX      :',F10.2)
 1042 FORMAT ( '       Space-time extremes DX-Y set to default 1000 m')
 1043 FORMAT ( '       Space-time extremes Dt      :',F8.2)
 1044 FORMAT ( '       Space-time extremes Dt set to default 1200 s')
!
 1100 FORMAT (/'  Status map, printed in',I3,' part(s) '/             &
               ' -----------------------------------'/)
 1101 FORMAT (2X,180I2)
 1102 FORMAT ( '  Legend : '/                                         &
               ' -----------------------------'/                      &
               '    0 : Land point            '/                      &
               '    1 : Sea point             '/                      &
               '    2 : Active boundary point '/                      &
               '    3 : Excluded point        '/)
 1103 FORMAT (/'  Obstruction map ',A1,', printed in',I3,' part(s) '/ &
               ' ---------------------------------------------'/)
 1104 FORMAT ( '  Legend : '/                                         &
               ' --------------------------------'/                   &
               '    fraction of obstruction * 10 '/)
 
 1105 FORMAT (/'  Shoreline slope, printed in',I3,' part(s) '/ &
               ' ---------------------------------------------'/)
 1106 FORMAT ( '  Legend : '/                                         &
               ' --------------------------------'/                   &
               '   Slope * 100'/)
 
 
 1150 FORMAT (/'  Reading unstructured grid definition files ...'/)
!
 9997 FORMAT (/'  Summary grid statistics : '/                        &
               ' --------------------------------------------------'/ &
               '       Number of longitudes      :',I10/              &
               '       Number of latitides       :',I10/              &
               '       Number of grid points     :',I10/              &
               '       Number of sea points      :',I10,' (',F4.1,'%)'/&
               '       Number of input b. points :',I10/              &
               '       Number of land points     :',I10/              &
               '       Number of excluded points :',I10/)
 9998 FORMAT (/'  Summary grid statistics : '/                        &
               ' --------------------------------------------------'/ &
               '       Number of longitudes      :',I10/              &
               '       Number of latitides       :',I10/              &
               '       Number of grid points     :',I10/              &
               '       Number of sea points      :',I10,' (100%)'/    &
               '       Number of input b. points :',I10/              &
               '       Number of land points     :',I10/              &
               '       Number of excluded points :',I10/)
 9999 FORMAT (/'  End of program '/                                   &
               ' ========================================'/           &
               '         WAVEWATCH III Grid preprocessor '/)
!
!/
!/ Internal function READNL ------------------------------------------ /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE READNL ( NDS, NAME, STATUS )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Jun-2013 |
!/                  +-----------------------------------+
!/
!  1. Purpose :
!
!     Read namelist info from file if namelist is found in file.
!
!  2. Method :
!
!     Look for namelist with name NAME in unit NDS and read if found.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   Data set number used for search.
!       NAME    C*4    I   Name of namelist.
!       STATUS  C*20   O   Status at end of routine,
!                            '(default values)  ' if no namelist found.
!                            '(user def. values)' if namelist read.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      EXTCDE    Subr. W3SERVMD Abort program as graceful as possible.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Program in which it is contained.
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDS
      CHARACTER, INTENT(IN)   :: NAME*4
      CHARACTER, INTENT(OUT)  :: STATUS*20
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IERR, I, J
      CHARACTER               :: LINE*80
!/
!/ ------------------------------------------------------------------- /
!/
!
      REWIND (NDS)
      STATUS  = '(default values) :  '
!
      DO
        READ (NDS,'(A)',END=800,ERR=800,IOSTAT=IERR) LINE
        DO I=1, 70
          IF ( LINE(I:I) .NE. ' ' ) THEN
              IF ( LINE(I:I) .EQ. '&' ) THEN
                  IF ( LINE(I+1:I+4) .EQ. NAME ) THEN
                      BACKSPACE (NDS)
                      SELECT CASE(NAME)
                        CASE('SLN1')
                          READ (NDS,NML=SLN1,END=801,ERR=802,IOSTAT=J)
                        CASE('SIN4')
                          READ (NDS,NML=SIN4,END=801,ERR=802,IOSTAT=J)
                        CASE('SNL1')
                          READ (NDS,NML=SNL1,END=801,ERR=802,IOSTAT=J)
                        CASE('SDS4')
                          READ (NDS,NML=SDS4,END=801,ERR=802,IOSTAT=J)
                        CASE('SBT1')
                          READ (NDS,NML=SBT1,END=801,ERR=802,IOSTAT=J)
                        CASE('SDB1')
                          READ (NDS,NML=SDB1,END=801,ERR=802,IOSTAT=J)
                        CASE('PRO3')
                          READ (NDS,NML=PRO3,END=801,ERR=802,IOSTAT=J)
                        CASE('UNST')
                          READ (NDS,NML=UNST,END=801,ERR=802,IOSTAT=J)
                        CASE('OUTS')
                          READ (NDS,NML=OUTS,END=801,ERR=802,IOSTAT=J)
                        CASE('MISC')
                          READ (NDS,NML=MISC,END=801,ERR=802,IOSTAT=J)
                        CASE DEFAULT
                          GOTO 803
                        END SELECT
                      STATUS  = '(user def. values) :'
                      RETURN
                    END IF
                ELSE
                  EXIT
                END IF
            ENDIF
          END DO
        END DO
!
  800 CONTINUE
      RETURN
!
  801 CONTINUE
      WRITE (NDSE,1001) NAME
      CALL EXTCDE(1)
      RETURN
!
  802 CONTINUE
      WRITE (NDSE,1002) NAME, J
      CALL EXTCDE(2)
      RETURN
!
  803 CONTINUE
      WRITE (NDSE,1003) NAME
      CALL EXTCDE(3)
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN READNL : '/          &
               '     PREMATURE END OF FILE IN READING ',A/)
 1002 FORMAT (/' *** WAVEWATCH III ERROR IN READNL : '/          &
               '     ERROR IN READING ',A,'  IOSTAT =',I8/)
 1003 FORMAT (/' *** WAVEWATCH III ERROR IN READNL : '/          &
               '     NAMELIST NAME ',A,' NOT RECOGNIZED'/)
!/
!/ End of READNL ----------------------------------------------------- /
!/
      END SUBROUTINE
!/
!/ End of W3GRID ----------------------------------------------------- /
!/
      END PROGRAM W3GRID
