#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!/
      MODULE W3SNL4MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Generic shallow-water Boltzmann integral (FBI or TSA)            !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
!! 1. Purpose :
!!
!!    Interface module for TSA type nonlinear interactions.
!!    Based on Resio and Perrie (2008) and Perrie and Resio (2009)
!!
!! 2. Variables and types :
!!
!!     Name      Type  Scope    Description
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!! 3. Subroutines and functions :
!!
!!     Name      Type  Scope    Description
!!    ------------------------------------------------------------------
!!     INSNL4    Subr. W3SNL4MD Corresponding initialization routine.
!!     ------
!!     W3SNL4    Subr. W3SNL4MD Main interface for TSA subroutines.
!!     ------                   Replaces main program "sboltz" in
!!                              "sbtsa-0-norm-Dec15-08.f" with
!!                              initialization done in subr. INSNL4
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!     shloxr    Subr. W3SNL4MD General locus solution
!!     shlocr    Subr. W3SNL4MD Locus solving routine - must converges
!!     cplshr    Subr. W3SNL4MD Computes Boltzmann coupling coefficient
!!     ------
!!op2
!!     Bash; Sections starting & ending with !!op2 are related to subr. optsa2
!!     optsa2    Subr. W3SNL4MD Converts the 2D Energy Density (f,theta)
!!     ------                   to Polar Action Density (k,theta) Norm. (in k)
!!                              then splits it into large and small scale
!!     --- - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!op2
!!     snlr_tsa  Subr. W3SNL4MD Computes dN(k,theta)/dt for TSA
!!     --------                 due to wave-wave inter. (set itsa = 1)
!!     snlr_fbi  Subr. W3SNL4MD Computes dN(k,theta)/dt for FBI
!!     --------                 due to wave-wave inter. (set itsa = 0)
!!
!!     wkfnc     fnc.  W3SNL4MD Compute wave number "k" for given
!!                              freq "f" (Hz) and water depth "d" (m)
!!                              or can use subr. WAVNU2
!!     cgfnc     fnc.  W3SNL4MD Compute group velocity "cg" for given
!!                              freq "f" (Hz), water depth "d" (m)
!!                              and phase speed "cvel" (m/s)
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines and functions used :
!!
!!     See subroutine documentation.
!!
!! 5. Remarks :
!!
!! 6. Switches :
!!
!!      !/S      Enable subroutine tracing.
!!      !/T(n)   Test output, see subroutines.
!!
!! 7. Source code :
!/
!!    --------------------------------------------------------------- &
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      PUBLIC
!!
!!
!!    ------------------------------------------------------------------
!!
!!
!!-0  Set these important run parameters here and declare them as PUBLIC
!! AC   itsa, ialt are set in mod_def and read here
!!      integer, parameter  :: itsa  =  1  !* = 1 for "snlr_tsa" or TSA
!!                           ****        !* = 0 for "snlr_fbi" or FBI
!!      integer, parameter  :: ialt  =  2  !* = 2 do    alternate in snlr's
!!                           ****        !* = 1 don't alternate in snlr's
      integer, parameter  :: ismo  =  1  !* = 1 do    smooth in interp2
!!                           ****        !* = 0 don't smooth in interp2
!!                                       !* interp2 is called only if ialt=2
      integer, parameter  :: npts  = 30  !* # of points on the locus
!!                           ****        !* can reduce npts for speed
      integer, parameter  :: ndep  = 37  !* # of depths in look-up tables
!!                           ****        !* can reduce ndep for speed
!!    ------------------------------------------------------------------
!!
!!
!!-1  Declare freq. related arrays & variables dim (nrng) and
!!            angle related arrays & variables dim (nang) as PUBLIC
      integer                            :: nrng, nzz, kzone, nb2fp
      integer                            :: nang, na2p1
      integer                            :: np2p1
      real                               :: dfrq, f0
      real                               :: ainc, twopi
      real,    allocatable, dimension(:) :: frqa, oma
      real,    allocatable, dimension(:) :: angl, sinan, cosan
      real,    allocatable, dimension(:) :: dep_tbl
!!    ------------------------------------------------------------------
!!
!!
!!-2  Declare gridsetr 11 look-up tables arrays dim (npts,nang,nzz,ndep)
!!    plus    pha_tbl array dim=(nrng,ndep) as PUBLIC
      integer, allocatable, dimension(:,:,:,:) :: kref2_tbl, kref4_tbl
      integer, allocatable, dimension(:,:,:,:) :: jref2_tbl, jref4_tbl
      real,    allocatable, dimension(:,:,:,:) :: wtk2_tbl,  wtk4_tbl
      real,    allocatable, dimension(:,:,:,:) :: wta2_tbl,  wta4_tbl
      real,    allocatable, dimension(:,:,:,:) :: tfac2_tbl, tfac4_tbl
      real,    allocatable, dimension(:,:,:,:) :: grad_tbl
      real,    allocatable, dimension(:,:)     :: pha_tbl
!!    ------------------------------------------------------------------
!!
!!
!!-3  Declare gridsetr 11 returned arrays dim (npts,nang,nzz) as PUBLIC
      integer, allocatable, dimension(:,:,:) :: kref2, kref4
      integer, allocatable, dimension(:,:,:) :: jref2, jref4
      real,    allocatable, dimension(:,:,:) :: wtk2,  wtk4
      real,    allocatable, dimension(:,:,:) :: wta2,  wta4
      real,    allocatable, dimension(:,:,:) :: tfac2, tfac4
      real,    allocatable, dimension(:,:,:) :: grad
!!    ------------------------------------------------------------------
!!
!!
!!-4  Declare shloxr/shlocr 5 returned arrays dim (npts) as PUBLIC
      real,    allocatable, dimension(:)     :: wk2x, wk2y
      real,    allocatable, dimension(:)     :: wk4x, wk4y, ds
!!    ------------------------------------------------------------------
!!
!!
!!-5  Declare w3snl4/optsa2 2 shared arrays dim (nrng,nang) & (nrng) as PUBLIC
!!    - ef2(nrng,nang) 'ww3' 2D Energy
!!    - ef1(nrng)      'ww3' 1D Energy from ef2()
      real,    allocatable, dimension(:,:)   :: ef2
      real,    allocatable, dimension(:)     :: ef1
!!    ------------------------------------------------------------------
!!
!!
!!-6  Declare optsa2 2 returned arrays dim (nrng,nang) as PUBLIC
!!    - dens1(nrng,nang)  'ww3' 2D Broad Scale Action
!!    - dens2(nrng,nang)  'ww3' 2D Small Scale Action
      real,    allocatable, dimension(:,:)   :: dens1, dens2
!!    ------------------------------------------------------------------
!!
!!
!!-7  Declare snlr's 4 returned arrays dim (nrng,nang) as PUBLIC
!!    tsa, diag  used for -tsa
!!    fbi, diag2 used for -fbi
      real,    allocatable, dimension(:,:)   :: tsa, diag
      real,    allocatable, dimension(:,:)   :: fbi, diag2
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
      CONTAINS
!!
!!
!!==============================================================================
!!
!!    ------------------------------------------------------------------
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE INSNL4
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!
!!    it returns: 11 look-up tables arrays dim=(npts,nang,nzz,ndep)
!!                kref2_tbl, kref4_tbl, jref2_tbl, jref4_tbl,
!!                wtk2_tbl,  wtk4_tbl,  wta2_tbl,  wta4_tbl,
!!                tfac2_tbl, tfac4_tbl & grad_tbl
!!                plus       pha_tbl  dim=(nrng,ndep)
!!                and        dep_tbl  dim=(ndep)
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!! 1. Purpose :
!!    It reads look-up tables (generated by gridsetr) if file exists
!!    otherwise it must generate the look-up tables file
!!
!! 2. Method :
!!    See subr gridsetr and subr W3SNL4 (or subr. W3IOGR)
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles
!!    npts      int.  Public    I   # of points on the locus
!!    ndep      int.  Public    I   # of depths in look-up tables
!!    dfrq      Real  Public    I   frequency multiplier for log freq. spacing
!!    dep_tbl   R.A.  Public    O   depthes in Look-up tables arrays dim=(ndep)
!!    grdfname  chr.  Local     -   Look-up tables filename (C*80)
!!    ------------------------------------------------------------------
!!
!!    *** The 11 look-up tables for grid integration geometry arrays
!!    *** at all selected 'ndep' depths defined in dep_tbl(ndep)' array
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2_tbl I.A.  Public    O   Index of reference wavenumber for k2
!!    kref4_tbl I.A.  Public    O   Idem for k4
!!    jref2_tbl I.A.  Public    O   Index of reference angle      for k2
!!    jref4_tbl I.A.  Public    O   Idem for k4
!!    wtk2_tbl  R.A.  Public    O   k2 Interpolation weigth along wavenumbers
!!    wtk4_tbl  R.A.  Public    O   Idem for k4
!!    wta2_tbl  R.A.  Public    O   k2 Interpolation weigth along angles
!!    wta4_tbl  R.A.  Public    O   Idem for k4
!!    tfac2_tbl R.A.  Public    O   Norm. for interp Action Density at k2
!!    tfac4_tbl R.A.  Public    O   Idem for k4
!!    grad_tbl  R.A.  Public    O   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!!    *** The 11 grid integration geometry arrays at one given depth
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2     I.A.  Public    O   Index of reference wavenumber for k2
!!    kref4     I.A.  Public    O   Idem for k4
!!    jref2     I.A.  Public    O   Index of reference angle      for k2
!!    jref4     I.A.  Public    O   Idem for k4
!!    wtk2      R.A.  Public    O   k2 Interpolation weigth along wavenumbers
!!    wtk4      R.A.  Public    O   Idem for k4
!!    wta2      R.A.  Public    O   k2 Interpolation weigth along angles
!!    wta4      R.A.  Public    O   Idem for k4
!!    tfac2     R.A.  Public    O   Norm. for interp Action Density at k2
!!    tfac4     R.A.  Public    O   Idem for k4
!!    grad      R.A.  Public    O   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     gridsetr  Subr. W3SERVMD Calc. the 11 grid geometry arrays for one depth
!!    ------------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     STRACE    Subr. W3SERVMD Subroutine tracing.
!!     W3SNL4    Subr. W3SNL4MD Interface for TSA  nonlinear interactions
!!     - or -    the option below was not used
!!     W3IOGR    Subr. W3INITMD Initialization (called by W3SHEL or WMINIT)
!!    ------------------------------------------------------------------
!!
!! 6. Error messages :
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
!/S      USE W3SERVMD, ONLY: STRACE
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO
!/MPI      USE WMMDATMD, ONLY: NMPSCR, IMPROC
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      IMPLICIT NONE
!!
!!    ==================================================================
!!
!!    Local variables & Parameters
!!    ----------------------------
!/S      INTEGER, SAVE           :: IENT = 0
!!
      integer              :: irng               !* dummy integer
      integer              :: nd                 !* dummy integer
!!
      logical              :: unavail = .TRUE.
      logical              :: file_exists
      character            :: grdfname*80
      integer              :: io_unit
!!
!!     
!!    Dimension wv# array and
!!    declare local var. dep2, cvel & cgnrng at nrng & depth 'nd'
      real                 :: wka2(nrng)  !* wv# array at depth nd (local)
      real                 :: pha2(nrng)  !* pha array at depth nd (local)
      real                 :: dep2        !* = dep_tbl(nd) depth at nd
      real                 :: cvel        !* Phase Velocity at (nrng,nd)
      real                 :: cgnrng      !* Group Velocity at (nrng,nd)
      real                 :: dwka        !* dummy storage for dk
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!/S      CALL STRACE (IENT, 'W3SNL4')
!!
!!    ==================================================================
!!
!!
!!-1  Make-up the filename from the main parameters
!!    ------------------------------------------------------------------
!b    example filename; grdfname = 'grd_1.1025_35_36_30_37.dat'
!!    grdfname = 'grd_dfrq_nr_na_np_nd.dat'
!!       where    fm = freq. mult. (dfrq) ex. dfrq = 1.1025 (F6.4)
!!                nr = # of rings  (nrng) ex. nrng = 35     (I2.2)
!!                na = # of angles (nang) ex. nang = 36     (I2.2)
!!                np = # of points (npts) ex. npts = 30     (I2.2)
!!                nd = # of depths (ndep) ex. nd   = 37     (I2.2)
      write(grdfname,'(A,F6.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A)')        &
      'grd_', dfrq,'_', nrng,'_', nang,'_',  npts,'_', ndep, '.dat'
!!    ==================================================================
!!
!!
!!-2  Check if the propre gridsetr Look-up tables file is available.
!!    if available read it and if not must generate it (by calling gridsetr)
!!    ------------------------------------------------------------------
      INQUIRE ( FILE=grdfname, EXIST=file_exists )
!!    Assign an unused UNIT number to io_unit.
!!    Note; It's important to look for an available unused number
      io_unit   = 60
      do while (unavail)
        io_unit = io_unit + 1
        INQUIRE ( io_unit, opened=unavail )
      enddo
!prt  print *, 'io_unit = ', io_unit
!!    ==================================================================
!!
!!
!!
!!
      IF ( file_exists ) THEN
!!
!!-3    File exists open it and read it
!!
!/MPI        if ( improc .eq. nmpscr ) then
        write ( ndso, 900 ) grdfname
!/MPI        end if
!!
        open (UNIT=io_unit, FILE=grdfname, STATUS='old',              &
              ACCESS='sequential', ACTION='read', FORM='unformatted')
        read (io_unit)  kref2_tbl, kref4_tbl, jref2_tbl, jref4_tbl,   &
                        wtk2_tbl,  wtk4_tbl,  wta2_tbl,  wta4_tbl,    &
                        tfac2_tbl, tfac4_tbl, grad_tbl,               &
                        pha_tbl,   dep_tbl
        close (io_unit)
!!      ----------------------------------------------------------------
!!
      ELSE      !* ELSE IF ( file_exists )
!!
!!
!!-4    File does not exist, create it here
!!
!/MPI        if ( improc .eq. nmpscr ) then
        write ( ndso, 901 ) grdfname
!/MPI        end if
!!      ----------------------------------------------------------------
!!
!!-4a   Define Look-up tables depth array 'dep_tbl(ndep)' for ndep=37
!!      with depths are +ve values
!!      ----------------------------------------------------------------
        dep_tbl(1:ndep) =                                             &
               (/  2.,  4.,  6.,  8., 10., 12., 14., 16., 18., 20.,   &
                  25., 30., 35., 40., 45., 50., 55., 60., 65., 70.,   &
                  80., 90.,100.,110.,120.,130.,140.,150.,160.,170.,   &
                 220.,270.,320.,370.,420.,470.,520.  /)
!prt    print *, ' ndep = ', ndep
!prt    print *, ' dep_tbl(1:ndep) = ', dep_tbl
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
        do 29 nd = 1,ndep
!!
!!
!!-4b     For given new depth dep2 = dep_tbl(nd) calculate
!!        a new array wka2(:) & new cgnrng corresp. to this depth
          dep2 = dep_tbl(nd)
          do irng=1,nrng
            wka2(irng) = wkfnc(frqa(irng),dep2)
          end do
          cvel    = oma(nrng)/wka2(nrng)         !* Phase Vel. at (nrng,nd)
          cgnrng  = cgfnc(frqa(nrng),dep2,cvel)  !* Group Vel. at (nrng,nd)
!!        --------------------------------------------------------------
!!        ==============================================================
!!
!!
!!-4c     Call gridsetr for this depth at nd
!!        --------------------------------------------------------------
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          call gridsetr ( dep2, wka2, cgnrng )
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!        it returns: 11 gridsetr arrays which are declared PUBLIC
!!                    kref2,kref4, jref2,jref4, wtk2,wtk4, wta2,wta4,
!!                    tfac2,tfac4  and   grad   all dim=(npts,nang,nzz)
!!        -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!        --------------------------------------------------------------
!!        ==============================================================
!!
!!-4d     Store in Look-up tables arrays at depth bin # 'nd'
          kref2_tbl(:,:,:,nd) = kref2(:,:,:)
          kref4_tbl(:,:,:,nd) = kref4(:,:,:)
          jref2_tbl(:,:,:,nd) = jref2(:,:,:)
          jref4_tbl(:,:,:,nd) = jref4(:,:,:)
          wtk2_tbl(:,:,:,nd)  = wtk2(:,:,:)
          wtk4_tbl(:,:,:,nd)  = wtk4(:,:,:)
          wta2_tbl(:,:,:,nd)  = wta2(:,:,:)
          wta4_tbl(:,:,:,nd)  = wta4(:,:,:)
          tfac2_tbl(:,:,:,nd) = tfac2(:,:,:)
          tfac4_tbl(:,:,:,nd) = tfac4(:,:,:)
          grad_tbl(:,:,:,nd)  = grad(:,:,:)
!!        --------------------------------------------------------------
!!        ==============================================================
!!
!!
!!-4e     Calculate pha2(:) at nd depth and store it in pha_tbl(:,:)
!!        pha2()=k*dk*dtheta, the base area at a grid intersection
!!        for use in integration of 2-D Density functions.
!!        --------------------------------------------------------------
!!        Below: variable dwka = dk centered at ring 1 (between 0 & 2)
!!        and computed    pha2(1) = k*dk*dtheta at ring 1
!!        with wkfnc(frqa(1)/dfrq,dep2) is like wka2(0)
!!        --assuming frqa(1)/dfrq is like frqa(0)
          dwka    = ( wka2(2) - wkfnc(frqa(1)/dfrq,dep2) ) / 2.
          pha2(1) = wka2(1)*dwka*ainc
!!
          do 23 irng=2,nrng-1
!!          Below: variable dwka = dk centered at irng (between irng-1 & irng+1)
!!          and computed    pha2(irng) = k*dk*dtheta at irng
            dwka       = ( wka2(irng+1) - wka2(irng-1) ) / 2.
            pha2(irng) = wka2(irng)*dwka*ainc
  23      continue
!!
!!        Below: variable dwka = dk centered at nrng (between nrng-1 & nrng+1)
!!        and computed    pha2(nrng) = k*dk*dtheta at nrng
!!        with wkfnc(dfrq*frqa(nrng),dep2) is like wka2(nrng+1)
!!        --assuming dfrq*frqa(nrng) is like frqa(nrng+1)
          dwka       = ( wkfnc(dfrq*frqa(nrng),dep2) - wka2(nrng-1) ) / 2.
          pha2(nrng) = wka2(nrng)*dwka*ainc
!!        --------------------------------------------------------------
!!        ==============================================================
!!
!!
!!-4f     Store pha2(:) at nd in  pha_tbl(:,nd) to be added to Look-up tables
          pha_tbl(1:nrng, nd) = pha2(1:nrng)
!!        --------------------------------------------------------------
!!        ==============================================================
!!
!!
  29    continue      !* End do 29 nd = 1,ndep
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-5    Ounce the Look-up tables arrays are full write it out to 'io_unit'
!!
!/MPI        if ( improc .eq. nmpscr ) then
        write( ndso,902 )
        open (UNIT=io_unit, FILE=grdfname, STATUS='new',              &
             ACCESS='sequential', ACTION='write', FORM='unformatted')
        write (io_unit) kref2_tbl, kref4_tbl, jref2_tbl, jref4_tbl,   &
                        wtk2_tbl,  wtk4_tbl,  wta2_tbl,  wta4_tbl,    &
                        tfac2_tbl, tfac4_tbl, grad_tbl,               &
                        pha_tbl,   dep_tbl
        close (io_unit)
        write( ndso,903 ) grdfname
!/MPI        end if
!!      ----------------------------------------------------------------
!!      ================================================================
!!
      ENDIF     !* End  IF ( file_exists )
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
      RETURN
!!
 900  format ( '  grdfname does exist = ',A/                          &
               '  open, read & close file ' )
!!
 901  format ( '  grdfname does not exist = ',A/                      &
               '  Generate look-up table arrays ' )
!!
 902  format ( '  Done generating look-up table arrays ----------- ' )
!!
 903  format ( '  Done writing & closing grdfname ', A )
!!
      END SUBROUTINE INSNL4
!!
!!==============================================================================
!!
!!    ------------------------------------------------------------------
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE W3SNL4 ( A, CG, WN, DEPTH,  S, D )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: S & D  dim = (NTH,NK) = (nang,nrng)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!! 1. Purpose :
!!
!!    Interface module for TSA type nonlinear interactions.
!!    Based on Resio and Perrie (2008) and Perrie and Resio (2009)
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    A         R.A.            I   2D Action Density A(NTH,NK) as function of
!!                                  direction (rad) and wavenumber (theta,k)
!!    CG        R.A.            I   Group velocities             dim=NK
!!    WN        R.A.            I   Wavenumbers                  dim=NK
!!    DEPTH     Real            I   Water depth (m)
!!    S         R.A.            O   Source term.                 dim=(NTH,NK)
!!    D         R.A.            O   Diagonal term of derivative. dim=(NTH,NK)
!!    ------------------------------------------------------------------
!!
!!    nrng      int.  Public    O   # of freq. or rings
!!    nang      int.  Public    O   # of angles
!!    npts      int.  Public    I   # of points on the locus
!!    ndep      int.  Public    I   # of depths in look-up tables
!!    nzz       int.  Public    O   linear irngxkrng = (NK*(NK+1))/2  
!!    kzone     int.  Public    O   zone of influence = INT(alog(4.0)/alog(dfrq))
!!    nb2fp     int.  Public    O   # of bins over fp => dfrq**nb2fp)*fp ~ 2.*fp
!!                                                    = INT(alog(2.0)/alog(dfrq))
!!    na2p1     int.  Public    O   = nang/2 + 1
!!    np2p1     int.  Public    O   = npts/2 + 1
!!    dfrq      Real  Public    O   frequency multiplier for log freq. spacing
!!    f0        Real  Public    O   = frqa(1); first freq. (Hz)
!!    ainc      Real  Public    O   = DTH; WW3 angle increment (radians)
!!    twopi     Real  Public    O   = TPI; WW3i 2*pi = 8.*atan(1.) (radians)
!!    oma       R.A.  Public    O   = SIG(1:NK)  WW3 waveumber array   dim=(nrng)
!!    frqa      R.A.  Public    O   = oma(:)/twopi WW3 frequency array dim=(nrng)
!!    angl      R.A.  Public    O   = TH(1:NTH);   WW3 angles array       dim=(nang)
!!    sinan     R.A.  Public    O   = ESIN(1:NTH); WW3 sin(angl(:)) array dim=(nang)
!!    cosan     R.A.  Public    O   = ECOS(1:NTH); WW3 cos(angl(:)) array dim=(nang)
!!    dep_tbl   R.A.  Public    I   depthes in Look-up tables arrays dim=(ndep)
!!    ------------------------------------------------------------------
!!
!!    *** The 11 look-up tables for grid integration geometry arrays
!!    *** at all selected 'ndep' depths defined in dep_tbl(ndep)' array
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2_tbl I.A.  Public    I   Index of reference wavenumber for k2
!!    kref4_tbl I.A.  Public    I   Idem for k4
!!    jref2_tbl I.A.  Public    I   Index of reference angle      for k2
!!    jref4_tbl I.A.  Public    I   Idem for k4
!!    wtk2_tbl  R.A.  Public    I   k2 Interpolation weigth along wavenumbers
!!    wtk4_tbl  R.A.  Public    I   Idem for k4
!!    wta2_tbl  R.A.  Public    I   k2 Interpolation weigth along angles
!!    wta4_tbl  R.A.  Public    I   Idem for k4
!!    tfac2_tbl R.A.  Public    I   Norm. for interp Action Density at k2
!!    tfac4_tbl R.A.  Public    I   Idem for k4
!!    grad_tbl  R.A.  Public    I   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!    
!!    *** The 11 grid integration geometry arrays at one given depth
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2     I.A.  Public    O   Index of reference wavenumber for k2
!!    kref4     I.A.  Public    O   Idem for k4
!!    jref2     I.A.  Public    O   Index of reference angle      for k2
!!    jref4     I.A.  Public    O   Idem for k4
!!    wtk2      R.A.  Public    O   k2 Interpolation weigth along wavenumbers
!!    wtk4      R.A.  Public    O   Idem for k4
!!    wta2      R.A.  Public    O   k2 Interpolation weigth along angles
!!    wta4      R.A.  Public    O   Idem for k4
!!    tfac2     R.A.  Public    O   Norm. for interp Action Density at k2
!!    tfac4     R.A.  Public    O   Idem for k4
!!    grad      R.A.  Public    O   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!!    ef2       R.A.  Public    O   2D Energy Density spectrum ef2(theta,f)
!!                                  = A(theta,k) * 2*pi*oma(f)/cga(f)
!!                                                     dim=(nrng,nang)
!!    ef1       R.A.  Public    O   1D Energy Density spectrum ef1(f)
!!                                                     dim=(nrng)
!!
!!    dens1     R.A.  Public    O   large-scale Action Density (k,theta)
!!                                                     dim=(nrng,nang)
!!    dens2     R.A.  Public    O   Small-scale Action Density (k,theta)
!!                                                     dim=(nrng,nang)
!!    ------------------------------------------------------------------
!!
!!    for -tsa; The 2 returned arrays tsa & diag dim=(nrng,nang)
!!    tsa       R.A.  Public    O   Snl-tsa = sumint + sumintsa
!!    diag      R.A.  Public    O   Snl-tsa diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!    
!!    for -fbi; The 2 returned arrays fbi & diag2 dim=(nrng,nang)
!!    fbi       R.A.  Public    O   Snl-fbi = sumint + sumintp  + sumintx
!!    diag2     R.A.  Public    O   Snl-fbi diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     STRACE    Subr. W3SERVMD Subroutine tracing.
!!     optsa2    Subr. W3SERVMD Converts the 2D Energy Density (f,theta)
!!     ------                   to Polar Action Density (k,theta) Norm. (in k)
!!                              then splits it into large and small scale
!!     snlr_tsa  Subr. W3SERVMD Computes dN(k,theta)/dt for TSA
!!     --------                 due to wave-wave inter. (set itsa = 1)
!!     snlr_fbi  Subr. W3SERVMD Computes dN(k,theta)/dt for FBI
!!     --------                 due to wave-wave inter. (set itsa = 0)
!!    ------------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     W3SRCE    Subr. w3srcemd Source term integration.
!!     W3EXPO    Subr. ww3_outp Point output post-processor.
!!     W3EXNC    Subr. ww3_ounp NetCDF Point output post-processor.
!!     GXEXPO    Subr. gx_outp  GrADS  Point output post-processor.
!!    ------------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S   Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      USE CONSTANTS, ONLY: TPI
      USE W3GDATMD,  ONLY: NK,  NTH,  XFR, DTH,  SIG, TH, ECOS, ESIN, &
                           ITSA, IALT
!!    dimension: SIG(0:NK+1),TH(NTH), ECOS(NSPEC+NTH), ESIN(NSPEC+NTH)
!!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3ODATMD, ONLY: NDSE, NDST, NDSO
!/S      USE W3SERVMD, ONLY: STRACE
!!    ==================================================================
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      REAL,    INTENT(IN)  :: A(NTH,NK), CG(NK), WN(NK), DEPTH
      REAL,    INTENT(OUT) :: S(NTH,NK), D(NTH,NK)
!!
      LOGICAL, SAVE        :: FIRST_TSA = .TRUE.
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!    Local Parameters & variables
!!    -----------------------------
!/S      INTEGER, SAVE           :: IENT = 0
      integer              :: irng, iang
      integer              :: nd3         !* bin # corresp. to ww3 dep  
      real                 :: dep         !* depth (m), get it from WW3 DEPTH
      real                 :: wka(NK)     !* from WW3 WN(1:NK) corresp. to "DEPTH"
      real                 :: cga(NK)     !* from WW3 CG(1:NK) corresp. to "DEPTH"
      real                 :: pha(NK)     !* k*dk*dtheta array corresp. to "DEPTH"
      real                 :: fac         !* twopi*oma()/cga()
      real                 :: sum1        !* dummy variable
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      integer              :: npk         !* bin# at peak frequency fpk
      integer              :: npk2        !* bin# of second peak frequency
      integer              :: npk0        !* dummy int. used in the shuffle of npk's
      integer              :: nsep        !* min # of bins that separates npk & npk2
                                          !* set nsep = 2
      integer              :: npeaks      !* # of peaks (=0, 1, or 2)
      integer              :: nfs         !* bin# of freq. separation
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      integer              :: nbins       !* actual # of bins > npk  (incl. nfs) or
!!                                        !* actual # of bins > npk2 (incl. nrng)
!!                                        !* to guarantee a min 1 bin in equi. range
      real                 :: fpk         !* peak frequency (Hz)
      real                 :: fpk2        !* second peak frequency (Hz)
      real                 :: e1max       !* 1D energy at 1st peak 'fpk'
      real                 :: e1max2      !* 1D energy at 2nd peak 'fpk2'
      real                 :: sumd1       !* sum dens1+dens2 at nfs
      real                 :: sumd2       !* sum dens1+dens2 at nfs+1
      real                 :: densat1     !* averaged dens1  at nfs
      real                 :: densat2     !* averaged dens1  at nfs+1
!!    --------------------------------------------------------------- &
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!/S      CALL STRACE (IENT, 'W3SNL4')
!!
!!ini
!!    Initialization of the output arrays
!!    before calling TSA subroutines.
!!    -----------------------------------
      S(:,:) = 0.0
      D(:,:) = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!    ==================================================================
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!
      IF ( FIRST_TSA ) THEN
!!
!!
!!-0    Set parameters & constants
!!      ---------------------------
        nrng  = NK                     !* nrng = NK  must be odd   <---
        nzz   = (NK * (NK+1)) / 2      !* linear irng, krng
        nang  = NTH                    !* nang = NTH must be even  <---
        na2p1 = nang/2 + 1             !* mid-angle or angle opposite to 1
        np2p1 = npts/2 + 1             !* mid-index of locus array
        twopi = TPI                    !* twopi = 8.*atan(1.)
                                       !* get it from WW3 TPI
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-1    Allocate freq & angle related array declared as PUBLIC
        if ( allocated (frqa) )    deallocate (frqa)
        if ( allocated (oma) )     deallocate (oma)
        if ( allocated (angl) )    deallocate (angl)
        if ( allocated (sinan) )   deallocate (sinan)
        if ( allocated (cosan) )   deallocate (cosan)
        if ( allocated (dep_tbl) ) deallocate (dep_tbl)
        allocate(frqa(nrng))
        allocate(oma(nrng))
        allocate(angl(nang))
        allocate(sinan(nang))
        allocate(cosan(nang))
        allocate(dep_tbl(ndep))
!!      ----------------------------------------------------------------
!!
!!-1a   Initialize frequency arrays and related parameters
!!      --------------------------------------------------
        oma(:)   = SIG(1:NK)           !* get it from WW3 SIG(1:NK)
        frqa(:)  = oma(:) / twopi
        f0       = frqa(1)
        dfrq     = XFR                 !* WW3 freq mult. for log freq
                                       !* get it from WW3 XFR
!!      ----------------------------------------------------------------
!!
!!-1b   Initialize direction arrays and related parameters
!!      --------------------------------------------------
        angl(:)  = TH(1:NTH)           !* get it from WW3 TH(1:NTH)
        cosan(:) = ECOS(1:NTH)         !* get it from WW3 ECOS(1:NTH)
        sinan(:) = ESIN(1:NTH)         !* get it from WW3 ESIN(1:NTH)
        ainc     = DTH                 !* WW3 angle increment (radians)
                                       !* get it from WW3 DTH
!!      ----------------------------------------------------------------
!!
!!-1c   Define kzone & nb2fp
!!kz
!!      kzone = zone of freq influence, function of dfrq
!!      for different values of x = 2,3,4 & 5
!!      So,    kzone(x) = INT( alog(x)/alog(dfrq) )
!!      +--------+----------+----------+----------+----------+
!!      | dfrq   | kzone(2) | kzone(3) | kzone(4) | kzone(5) |
!!      +--------+----------+----------+----------+----------+
!!      | 1.05   |    14    |    22    |    28    |    33    |
!!      +--------+----------+----------+----------+----------+
!!      | 1.07   |    10    |    16    |    20    |    24    |
!!      +--------+----------+----------+----------+----------+
!!      | 1.10   |     7    |    11    |    14    |    17    |
!!      +--------+----------+----------+----------+----------+
        kzone = INT( alog(2.0)/alog(dfrq) )  !* Bash; faster without loss of accuracy
!kz     kzone = INT( alog(3.0)/alog(dfrq) )  !* as in gridsetr & snlr_'s
!kz     kzone = INT( alog(4.0)/alog(dfrq) )  !* as in gridsetr & snlr_'s
!kz     kzone = INT( alog(5.0)/alog(dfrq) )  !* as in gridsetr & snlr_'s
!!kz---
!!
!!op2
!!      nb2fp = # of bins over fp (not incl. fp) - this depends on dfrq
!!              so that (dfrq**nb2fp)*fp ~ 2.*fp  (like kzone(2))
!!              used in 1 bin equi. range
        nb2fp = INT( alog(2.0)/alog(dfrq) )  !* for equi. range near 2*fp
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!op2
!!      ================================================================
!!
!!
!!-2    Allocate gridsetr 11 look-up tables arrays
!!      plus     pha_tbl array dim=(nrng,ndep) declared as PUBLIC
        if ( allocated (kref2_tbl) ) deallocate (kref2_tbl)
        if ( allocated (kref4_tbl) ) deallocate (kref4_tbl)
        allocate(kref2_tbl(npts,nang,nzz,ndep))
        allocate(kref4_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (jref2_tbl) ) deallocate (jref2_tbl)
        if ( allocated (jref4_tbl) ) deallocate (jref4_tbl)
        allocate(jref2_tbl(npts,nang,nzz,ndep))
        allocate(jref4_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (wtk2_tbl) ) deallocate (wtk2_tbl)
        if ( allocated (wtk4_tbl) ) deallocate (wtk4_tbl)
        allocate(wtk2_tbl(npts,nang,nzz,ndep))
        allocate(wtk4_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (wta2_tbl) ) deallocate (wta2_tbl)
        if ( allocated (wta4_tbl) ) deallocate (wta4_tbl)
        allocate(wta2_tbl(npts,nang,nzz,ndep))
        allocate(wta4_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (tfac2_tbl) ) deallocate (tfac2_tbl)
        if ( allocated (tfac4_tbl) ) deallocate (tfac4_tbl)
        allocate(tfac2_tbl(npts,nang,nzz,ndep))
        allocate(tfac4_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (grad_tbl) ) deallocate (grad_tbl)
        allocate(grad_tbl(npts,nang,nzz,ndep))
!!
        if ( allocated (pha_tbl) ) deallocate (pha_tbl)
        allocate(pha_tbl(nrng,ndep))
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-3    Allocate gridsetr 11 returned arrays declared as PUBLIC
        if ( allocated (kref2) ) deallocate (kref2)
        if ( allocated (kref4) ) deallocate (kref4)
        allocate(kref2(npts,nang,nzz))
        allocate(kref4(npts,nang,nzz))
!!
        if ( allocated (jref2) ) deallocate (jref2)
        if ( allocated (jref4) ) deallocate (jref4)
        allocate(jref2(npts,nang,nzz))
        allocate(jref4(npts,nang,nzz))
!!
        if ( allocated (wtk2) ) deallocate (wtk2)
        if ( allocated (wtk4) ) deallocate (wtk4)
        allocate(wtk2(npts,nang,nzz))
        allocate(wtk4(npts,nang,nzz))
!!
        if ( allocated (wta2) ) deallocate (wta2)
        if ( allocated (wta4) ) deallocate (wta4)
        allocate(wta2(npts,nang,nzz))
        allocate(wta4(npts,nang,nzz))
!!
        if ( allocated (tfac2) ) deallocate (tfac2)
        if ( allocated (tfac4) ) deallocate (tfac4)
        allocate(tfac2(npts,nang,nzz))
        allocate(tfac4(npts,nang,nzz))
!!
        if ( allocated (grad) ) deallocate (grad)
        allocate(grad(npts,nang,nzz))
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-4    Allocate shloxr/shlocr 5 returned arrays declared as PUBLIC
        if ( allocated (wk2x) ) deallocate (wk2x)
        if ( allocated (wk2y) ) deallocate (wk2y)
        allocate(wk2x(npts))
        allocate(wk2y(npts))
!!
        if ( allocated (wk4x) ) deallocate (wk4x)
        if ( allocated (wk4y) ) deallocate (wk4y)
        allocate(wk4x(npts))
        allocate(wk4y(npts))
!!
        if ( allocated (ds) ) deallocate (ds)
        allocate(ds(npts))
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-5    Allocate w3snlx/optsa2 2 shared arrays declared as PUBLIC
        if ( allocated (ef2) ) deallocate (ef2)
        if ( allocated (ef1) ) deallocate (ef1)
        allocate(ef2(nrng,nang))
        allocate(ef1(nrng))
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-6    Allocate optsa2 2 returned arrays declared as PUBLIC
        if ( allocated (dens1) ) deallocate (dens1)
        if ( allocated (dens2) ) deallocate (dens2)
        allocate(dens1(nrng,nang))
        allocate(dens2(nrng,nang))
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-7    Allocate snlr_??? 2 returned arrays declared as PUBLIC
        if ( itsa .eq. 1) then
!!        allocate tsa, diag  used for -tsa
          if ( allocated (tsa) ) deallocate (tsa)
          if ( allocated (diag) ) deallocate (diag)
          allocate(tsa(nrng,nang))
          allocate(diag(nrng,nang))
        elseif ( itsa .eq. 0) then
!!        allocate fbi, diag2 used for -fbi
          if ( allocated (fbi) ) deallocate (fbi)
          if ( allocated (diag2) ) deallocate (diag2)
          allocate(fbi(nrng,nang))
          allocate(diag2(nrng,nang))
        else
          write ( ndse,1000 ) itsa
          CALL EXTCDE ( 115 )
        endif
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
!!-8    Get the 11 look-up table arrays by calling INSNL4
!!      ----------------------------------------------------------------
!!
!!      ----------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call INSNL4
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!      it returns: 11 look-up tables arrays dim=(npts,nang,nzz,ndep)
!!                  kref2_tbl, kref4_tbl, jref2_tbl, jref4_tbl,
!!                  wtk2_tbl,  wtk4_tbl,  wta2_tbl,  wta4_tbl,
!!                  tfac2_tbl, tfac4_tbl & grad_tbl
!!                  plus       pha_tbl  dim=(nrng,ndep)
!!                  and        dep_tbl  dim=(ndep)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
        FIRST_TSA  = .FALSE.
!!
!!
      ENDIF    !! IF ( FIRST_TSA ) THEN
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    ##################################################################
!!
!!
!!
!!*i1 Map input ww3 "DEPTH" to "dep" and find corresp. depth bin # "nd3"
!!    ------------------------------------------------------------------
      dep = DEPTH                !* ww3 depth at a given time & loc.
      nd3 = MINLOC( abs(dep - dep_tbl(:)), dim=1 )
!prt  print *, 'DEPTH, corresp depth bin # (nd3)  = ', DEPTH, nd3
!!    ------------------------------------------------------------------
!!
!!
!!*i2 Map from Look-up tables the 11 gridsetr arrays corresp. to "nd3"
!!    kref2(:,:,:) -> grad(:,:,:) are used in subrs. "snlr_*"
!!    ------------------------------------------------------------------
      kref2(:,:,:) = kref2_tbl(:,:,:,nd3)
      kref4(:,:,:) = kref4_tbl(:,:,:,nd3)
      jref2(:,:,:) = jref2_tbl(:,:,:,nd3)
      jref4(:,:,:) = jref4_tbl(:,:,:,nd3)
      wtk2(:,:,:)  = wtk2_tbl(:,:,:,nd3)
      wtk4(:,:,:)  = wtk4_tbl(:,:,:,nd3)
      wta2(:,:,:)  = wta2_tbl(:,:,:,nd3)
      wta4(:,:,:)  = wta4_tbl(:,:,:,nd3)
      tfac2(:,:,:) = tfac2_tbl(:,:,:,nd3)
      tfac4(:,:,:) = tfac4_tbl(:,:,:,nd3)
      grad(:,:,:)  = grad_tbl(:,:,:,nd3)
      pha(:)       = pha_tbl(:,nd3)
!!    ------------------------------------------------------------------
!!
!!
!!*i3 Map input ww3 arrays "WN(:)" & "CG(:)" to "wka(:)" & "cga(:)"
!!    Note; Arrays wka(:) & cga(:) corresp to ww3 "DEPTH" & to be used in "optsa2"
!!    ------------------------------------------------------------------
      wka(:)  = WN(1:NK)               !* Wavenumber array at ww3 "DEPTH"
      cga(:)  = CG(1:NK)               !* Group velocity array at ww3 "DEPTH"
!!    ------------------------------------------------------------------
!!
!!
!!*i4 Convert input WW3 2D Action Density spectrum "A(theta,k)"
!!    to 2D Energy Density spectrum "ef2(theta,f)" & reverse indices
!!    ==>  ef2(f,theta) = A(theta,k) * 2*pi*oma(f)/cga(f)
!!    ------------------------------------------------------------------
      do 32 irng=1,nrng
        fac = twopi*oma(irng)/cga(irng)
        do 31 iang=1,nang
          ef2(irng,iang) = A(iang,irng) * fac
  31    continue
  32  continue
!!    ------------------------------------------------------------------
!!
!!
!!*i5 Calculte the 1D Energy Density "ef1(f)"
!!    ------------------------------------------------------------------
      do 42 irng=1,nrng
        sum1 = 0.0
        do 41 iang=1,nang
          sum1 = sum1 + ef2(irng,iang)
  41    continue
        ef1(irng) = sum1 * ainc
  42  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!op2
!!*   Bash;
!!*   Find 1 or 2 peaks that satisfy TSA min condition (below) ------- *
!!*   before calling TSA subrs. otherwise bailout (return) ----------- *
!!*   Bailout & return with init. values of S & D = 0.0 -------------- *
!!*   nsep   = min # of bins that separates between npk & npk2 (set=2) *
!!*   nbins  = actual # of bins > npk  (incl. nfs)  -- or --           *
!!*            actual # of bins > npk2 (incl. nrng)                    *
!!*            to guarantee a min 1 bin in equi. range                 *
!!*                                                                    *
!!*   ===>  In case of just 1 peak the TSA min condition ------------- *  *****
!!*   ===>  is relative to nrng and is satisfied when ---------------- *  <<<<<
!!*   ===>  npk.le.nrng-1, to guarantee min 1 bin (incl nrng) > npk -- *  <<<<<
!!*   ===>  we only need 1 bin in optsa2 to be in the equi. range ---- *  <<<<<
!!*   ===>  skip if condition is not met ie if  npk.gt.nrng-1  ------- *  <<<<<
!!*   ---------------------------------------------------------------- *
!!*                                                                    *
!!*   ===>  In case of 2 peaks the TSA min condition is applied twice: *  *****
!!*   ===>   *1) at low freq peak (npk), *2) at high freq peak (npk2)  *  *****
!!*   ===> *1) TSA min condition for the low freq peak (npk) --------- *  *****
!!*   ===>  is relative to nfs and is satisfied when ----------------- *  <<<<<
!!*   ===>  npk.le.nfs-1, to guarantee min 1 bin (incl nfs) > npk2 --- *  <<<<<
!!*   ===>  we only need 1 bin in optsa2 to be in the equi. range ---- *  <<<<<
!!*   ===>  skip if condition is not met ie if  npk.gt.nfs-1  -------- *  <<<<<
!!*                                                                    *
!!*   ===> *2) TSA min condition for the high freq peak (npk2) ------- *  *****
!!*   ===>  is relative to nrng and is satisfied when ---------------- *  <<<<<
!!*   ===>  npk2.le.nrng-1 to guarantee min 1 bin (incl nrng) > npk2   *  <<<<<
!!*   ===>  we only need 1 bin in optsa2 to be in the equi. range ---- *  <<<<<
!!*   ===>  skip if condition is not met ie if  npk2.gt.nrng-1 ------- *  <<<<<
!!*   ---------------------------------------------------------------- *
!!    ------------------------------------------------------------ !!op2
!!    ==================================================================
!!
!!op2 ctd
!!    First find the overall peak in ef1(:) with e1max must be > 0.000001
!!    Starting from low freq. find the Energy max "e1max" and
!!    corresp. peak freq. "fpk" and its freq. number "npk".
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      npk    = 0
      fpk    = 0.0
      e1max  = 0.0
      npeaks = 0
!!    Look in the freq range that works for TSA call (see condition below)
      do 43 irng=2,nrng-1             !* last peak loc. is at nrng-1      <<<<<
!!      Pick the 1st local abs. max in [2,nrng-1] using (ef1(irng).gt.e1max)
!!      so that if 2 equal adj. peaks are found it will pick the 1st e1max
!!      encountered (i.e. the lower freq. one)
!!      --------------------------------------------------------------!*  <<<<<
        if ( ef1(irng).gt.ef1(irng-1) .and. ef1(irng).gt.ef1(irng+1)  &
                                      .and. ef1(irng).gt.e1max ) then
!!      --------------------------------------------------------------!*  <<<<<
          npk    = irng               !* update npk
          fpk    = frqa(npk)          !* update fpk
          e1max  = ef1(npk)           !* update e1max
          npeaks = 1
        endif
  43  continue
!!    ------------------------------------------------------------------
!!
!!B   if a 1st peak is not found (npeaks=0 & e1max=0.0 < eps) or
!!B   if a 1st peak is found with a tiny peak energy (e1max < eps) or
!!B   if TSA min condition is not met rel. to nrng (npk.gt.nrng-1)   <<<<<
!!B   this spectrum is Not suitable for tsa, so don't call tsa
!!B   just return (don't stop) with init. values of S(:,:) and D(:,:)=0.0
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( e1max.lt.0.000001 ) return
!!    ------------------------------------------------------------ !!op2
!!    ==================================================================
!!
!!
!!op2 ctd
!!    Bash; if we are here (i.e. we did not return) then we must
!!    have found the 1st good peak (= overall peak) with e1max > eps
!!
!!    Now we look for a new 2nd peak that is at least 'nsep' bins away from
!!    the 1st peak (nsep=2) (i.e. iabs(irng-npk).gt.nsep) before 
!!    calling new "optsa2" (the 2nd peak will have e1max2 < e1max)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nsep   = 2
      npk2   = 0
      fpk2   = 0.0
      e1max2 = 0.0
!!    Again look in the freq range that is in line with TSA min condition
!!    and find the 2nd highest peak with  eps < e1max2 < e1max
      do 45 irng=2,nrng-1           !* last peak loc. is at nrng-1      <<<<<
!!      Pick the 2nd local abs. max in [2,nrng-1] that is at least 'nsep'
!!      bins away from the 1st peak using (ef1(irng).ge.e1max2) so that
!!      if 2 equal adj. peaks are found it will pick the 2nd e1max2
!!      encountered (i.e. the higher freq. one)
!!      --------------------------------------------------------------!*  <<<<<
        if ( ef1(irng).gt.ef1(irng-1) .and. ef1(irng).gt.ef1(irng+1)  &
         .and. ef1(irng).ge.e1max2 .and. iabs(irng-npk).gt.nsep )  then
!!      --------------------------------------------------------------!*  <<<<<
          npk2   = irng             !* update npk2
          fpk2   = frqa(npk2)       !* update fpk2
          e1max2 = ef1(npk2)        !* update e1max2
          npeaks = 2
        endif
  45  continue
!!    ------------------------------------------------------------------
!!
!!B   if a 2nd peak is not found (npeaks=1 & e1max2=0.0 < eps)
!!B   if a 2nd peak is found with a tiny peak energy (e1max2 < eps) or
!!B   if TSA min condition is not met rel. to nrng (npk2.gt.nrng-1)   <<<<<
!!B   This 2nd peak is not suitable for tsa, drop it and stay with just 1st peak.
      if ( e1max2.lt.0.000001 ) then
        npeaks = 1
        goto 200           !* skip the remaings tests goto 200
      endif
!!    ------------------------------------------------------------ !!op2
!!    ==================================================================
!!
!!
!!
!!
!!op2 ctd
      if ( npeaks.eq.2 ) then
!!-1    Shuffle the 2 peaks (if necessary) to keep npk to be always < npk2
!!      This says nothing about which peak is the dominant peak
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if ( npk2.lt.npk ) then
          npk0   = npk2
          npk2   = npk
          npk    = npk0                 !*  this way  npk < npk2  always
          fpk    = frqa(npk)
          fpk2   = frqa(npk2)
        endif
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!-2    here we have 2 peaks (npeaks=2) with  npk < npk2
!!      find the freq. separation "nfs" (that divide the freq. regime into 2)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nfs = INT ( (npk+npk2)   / 2.0 )  !* take the lower  bin # to be nfs
!b      nfs = INT ( (npk+npk2+1) / 2.0 )  !* take the higher bin # to be nfs
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endif   !! if ( npeaks.eq.2 )
!!
 200  continue
!!    ------------------------------------------------------------ !!op2
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
!!    Bash; With the new "optsa2" you are allowed one call (if 1 peak)
!!          or 2 calls (if 2 peaks) o account for spectra with double peaks.
!!    Note; when nrmn=1 & nrmx=nrng ==> optsa2 = the old optsa
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      if ( npeaks.eq.1 ) then
!!-1    one call to optsa2 for the whole freq. regime ( 1 --> nrng )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nbins = nrng - npk    !* # of bins in (npk, nrng] not incl. npk
        if ( nbins.gt.nb2fp ) nbins=nb2fp !* limit equi. range to ~2.0*fp
!!      ----------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call optsa2 ( 1,nrng,       npk, fpk,  nbins, wka, cga )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      It returns variables dens1(nrng,nang) and dens2(nrng,nang)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
      endif   !! if ( npeaks.eq.1 )
!!    ==================================================================
!!
!!
      if ( npeaks.eq.2 ) then
!!
!!-2    Now make two calls to new "optsa2" one for each freq regime.
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        nbins = nfs - npk     !* # of bins in (npk, nfs] not incl. npk
        if ( nbins.gt.nb2fp ) nbins=nb2fp  !* limit equi. range to ~2.0*fp
!!      ---------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call optsa2 ( 1,nfs,        npk, fpk,  nbins, wka, cga )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      It returns variables dens1(nrng,nang) and dens2(nrng,nang)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
!!
        nbins = nrng - npk2   !* # of bins in (npk2, nrng] not incl. npk2
        if ( nbins.gt.nb2fp ) nbins=nb2fp  !* limit equi. range to ~2.0*fp
!!      ----------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call optsa2 ( nfs+1,nrng,   npk2,fpk2, nbins, wka, cga )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      It returns variables dens1(nrng,nang) and dens2(nrng,nang)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!-3    Remove the step like jump (if exists) in dens1() between nfs & nfs+1
        do 440 iang=1,nang
          sumd1 = dens1(nfs,iang)   + dens2(nfs,iang)   !* sum at nfs
          sumd2 = dens1(nfs+1,iang) + dens2(nfs+1,iang) !* sum at nfs+1
!!
!!        do 3 bin average for dens1() at nfs   and store in densat1
          densat1 = ( dens1(nfs-1,iang) + dens1(nfs,iang) +           &
                                          dens1(nfs+1,iang) ) / 3.
!!        do 3 bin average for dens1() at nfs+1 and store in densat2
          densat2 = ( dens1(nfs,iang)   + dens1(nfs+1,iang) +         &
                                          dens1(nfs+2,iang) ) / 3.
!!
!!        subtitute back into dens1(nfs,iang) & dens1(nfs+1,iang)
          dens1(nfs,iang)   = densat1             ! dens1 at nfs
          dens1(nfs+1,iang) = densat2             ! dens1 at nfs+1
!!
!!        recalculate dens2(nfs,iang) & dens2(nfs+1,iang)
          dens2(nfs,iang)   = sumd1 - densat1     ! dens2 at nfs
          dens2(nfs+1,iang) = sumd2 - densat2     ! dens2 at nfs+1
 440    continue
!!
      endif   !! if ( npeaks.eq.2 )
!!    ==================================================================
!!    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Get Snl source term and its diagonal term  from "snlr"           !
!!    for -tsa  only       use "snlr_tsa"       itsa = 1               !
!!    for -fbi  only       use "snlr_fbi"       itsa = 0               !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
      if ( itsa .eq. 1) then
!!
!!      ----------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call snlr_tsa ( pha, ialt )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      It returns tsa(nrng,nang) & diag(nrng,nang)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
!!
!!      Pack results in proper format ---------------------------------- *
!!      S() & D() arrays are to be returned to WW3 in (k,theta) space
        do 52 irng=1,nrng
        do 51 iang=1,nang
!!        Convert the Norm. (in k) Polar tsa(k,theta) to Polar S(theta,k)
!!        and  reverse indices back to (iang,irng) as in WW3
          S(iang,irng) = tsa(irng,iang) * wka(irng)   !* <=============
          D(iang,irng) = diag(irng,iang)
!!        ---------------------------
  51    continue
  52    continue
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!
      elseif ( itsa .eq. 0) then
!!
!!
!!      ----------------------------------------------------------------
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call snlr_fbi ( pha, ialt )
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      It returns fbi(nrng,nang) & diag2(nrng,nang)
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!      ----------------------------------------------------------------
!!
!!      Pack results in proper format ---------------------------------- *
!!      S() & D() arrays are to be returned to WW3 in (k,theta) space
        do 54 irng=1,nrng
        do 53 iang=1,nang
!!        Convert the Norm. (in k) Polar fbi(k,theta) to Polar S(theta,k)
!!        and  reverse indices back to (iang,irng) as in WW3
          S(iang,irng) = fbi(irng,iang) * wka(irng)   !* <=============
          D(iang,irng) = diag2(irng,iang)
!!        --------------------------------
  53    continue
  54    continue
!!      -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      else
!!
         write( ndse,1000 ) itsa
         CALL EXTCDE ( 130 )
!!       --- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      endif
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
 1000 format ( ' W3SNL4 Error : Bad itsa value ',i4)
!!
      END SUBROUTINE W3SNL4
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE gridsetr ( dep, wka1, cgnrng1 )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: kref2,kref4, jref2,jref4, wtk2,wtk4, wta2,wta4,
!!                tfac2,tfac4  and    grad       all dim=(npts,nang,nzz)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -------------------------------------------------------------------------#
!!                                                                             !
!!    This routine sets up the geometric part of the Boltzmann integral        !
!!    based on a grid of wave frequencies and directions, with wave-           !
!!    numbers related to frequency and depth by linear dispersion.  It         !
!!    is adapted from Don's original code with changes to modify the           !
!!    indexing so there are fewer unused elements, and a number of algo-       !
!!    rithmic changes that are mathematically equivalent to Don's but          !
!!    take advantage of intrinsic functions to form smooth results with        !
!!    less reliance on if statements.                                          !
!!                                                                             !
!!    It calls locus-solving routines shloxr and shlocr and coupling           !
!!    coefficient routine cplshr.  If shlocr does not converge, ierr_gr        !
!!    will be something other than 0 and the routine will terminate,           !
!!    returning ierr_gr to the calling program (see shlocr).                   !
!!                                                                             !
!!    It returns array grad(,,), which is an estimate of the product           !
!!    C(k1,k2,k3,k4)*H(k1,k3,k4)*ds/|dW/dn| (where n and the k's are all       !
!!    vectors) as given, for example, by Eq.(7) of 'Nonlinear energy           !
!!    fluxes and the finite depth equilibrium range in wave spectra,'          !
!!    by Resio, Pihl, Tracy and Vincent (2001, JGR, 106(C4), p. 6985),         !
!!    as well as arrays for indexing, interpolating and weighting locus-       !
!!    based wavenumber vectors within the discrete solution grid.              !
!!    -------------------------------------------------------------------------#
!!
!!
!! 2. Method :
!!
!! 3. Parameters : 
!!
!!    Parameter list 
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles
!!    npts      int.  Public    I   # of points on the locus
!!    nzz       int.  Public    I   linear irngxkrng = (NK*(NK+1))/2
!!    kzone     int.  Public    I   zone of influence = INT(alog(4.0)/alog(dfrq))
!!    na2p1     int.  Public    I   = nang/2 + 1
!!    np2p1     int.  Public    I   = npts/2 + 1
!!    ------------------------------------------------------------------
!!
!!    dfrq      Real  Public    I   frequency multiplier for log freq. spacing
!!    f0        Real  Public    I   = frqa(1); first freq. (Hz)
!!    twopi     Real  Public    I   = TPI; WW3i 2*pi = 8.*atan(1.) (radians)
!!    ainc      Real  Public    I   = DTH; WW3 angle increment (radians)
!!    dep       Real  local     I   = depth (m)
!!    frqa      R.A.  Public    I   = oma(:)/twopi WW3 frequency array dim=(nrng)
!!    angl      R.A.  Public    I   = TH(1:NTH);   WW3 angles array       dim=(nang)
!!    sinan     R.A.  Public    I   = ESIN(1:NTH); WW3 sin(angl(:)) array dim=(nang)
!!    cosan     R.A.  Public    I   = ECOS(1:NTH); WW3 cos(angl(:)) array dim=(nang)
!!    wka1      R.A.  local     I   = wavenumber array at one depth dim=(nrng)
!!    cgnrng1   Real  local     I   = Group Vel. at nrng at one depth
!!    ------------------------------------------------------------------
!!
!!    *** The 11 grid integration geometry arrays at one given depth
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2     I.A.  Public    O   Index of reference wavenumber for k2
!!    kref4     I.A.  Public    O   Idem for k4
!!    jref2     I.A.  Public    O   Index of reference angle      for k2
!!    jref4     I.A.  Public    O   Idem for k4
!!    wtk2      R.A.  Public    O   k2 Interpolation weigth along wavenumbers
!!    wtk4      R.A.  Public    O   Idem for k4
!!    wta2      R.A.  Public    O   k2 Interpolation weigth along angles
!!    wta4      R.A.  Public    O   Idem for k4
!!    tfac2     R.A.  Public    O   Norm. for interp Action Density at k2
!!    tfac4     R.A.  Public    O   Idem for k4
!!    grad      R.A.  Public    O   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     shloxr    Subr. W3SERVMD General locus solution for input vectors
!!                              k1 & k3 when |k1| .eq. |k3|
!!     shlocr    Subr. W3SERVMD General locus solution for input vectors
!!                              k1 & k3 when |k1| .ne. |k3|
!!     cplshr    Subr. W3SERVMD Calculates Boltzmann coupling coefficient
!!                              in shallow water
!!    ------------------------------------------------------------------
!! 
!! 5. Called by : 
!! 
!!     Name      Type  Module   Description 
!!    ------------------------------------------------------------------
!!     insnl4    Subr. W3SNL4MD initialize the grid geometry
!!    ------------------------------------------------------------------
!! 
!! 6. Error messages : 
!! 
!!      None. 
!! 
!! 7. Remarks : 
!! 
!! 8. Structure : 
!! 
!!    See source code. 
!! 
!! 9. Switches : 
!! 
!!    !/S  Enable subroutine tracing. 
!! 
!!10. Source code : 
!! 
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: dep
      real,    intent(in)  :: wka1(nrng), cgnrng1  !* Use new names locally
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    Local Parameters & variables
!!    -----------------------------
      integer              :: irng,krng, iang,kang, ipt
      integer              :: iizz, izz, ir, i
      integer              :: kmax               !* = min(irng+kzone, nrng)
!!
      real                 :: g, gsq
      real                 :: alf0,aldfrq, wk1x,wk1y, wk3x,wk3y
!!
      real                 :: wn2,th2, wn2d,tnh2, om2,f2,cg2, tt2,w2
      real                 :: wn4,th4, wn4d,tnh4, om4,f4,cg4, tt4,w4
      real                 :: dWdnsq,dWdn, dif13,dif14, er
!!
!!hv  Bash; with !hv ON, move Heaviside section up & don't use var Heaviside
!hv   real                 :: Heaviside
!!hv---
!!
      real                 :: csq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    initial constants
!!    ------------------
      g      = 9.806                              !* set = GRAV as in CONSTANTS
      gsq    = 96.157636                          !* set = GRAV**2
!!
      alf0   = alog(frqa(1))                      !* ln(f0) for ir calc. below
      aldfrq = alog(dfrq)                         !* ln(dfrq)       "
!!
!!ini
!!    initialize array grad
!!    ----------------------
      grad(:,:,:)  = 0.0
      kref2(:,:,:) = 0
      kref4(:,:,:) = 0
      jref2(:,:,:) = 0
      jref4(:,:,:) = 0
      wtk2(:,:,:)  = 0.0
      wtk4(:,:,:)  = 0.0
      wta2(:,:,:)  = 0.0
      wta4(:,:,:)  = 0.0
      tfac2(:,:,:) = 0.0
      tfac4(:,:,:) = 0.0
!!ini---
!!------------------------------------------------------------------------------
!!
!!
!!    irng and iang are k1 parameters; krng and kang are k3 parameters
      iang = 1                               !* set = 1 and will remain = 1
!!
!!20
      do 20 irng=1,nrng
!!kz
        kmax = min(irng+kzone, nrng)   !* Bash; Sometimes a locus pt is outside nrng
!kz     kmax = min(irng+kzone, nrng-1) !* Bash; Taking 1 out will not affect kzone, try it
!!kz---
!!kz---
!!
        wk1x   = wka1(irng)
        wk1y   = 0.0                         !* set = 0.0 and will remain = 0.0
        iizz = (nrng-1)*(irng-1)-((irng-2)*(irng-1))/2
!!30
!!kz
        do 30 krng=irng,kmax
!!kz---
!kz     do 30 krng=irng,nrng
!!
!!        Bash; check1 - change this ratio from > 4 to > 3   and
!!              make it consistent with similar test done in subr. snlr_'s
!kz       if ( frqa(krng)/frqa(irng) .gt. 2. ) go to 30  !* Bash; use .gt. 2 for speed
!kz       if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 30  !* original snlr_'s
!kz       if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 30  !* original gridsetr
!!kz---
          izz = krng+iizz
!!40
          do 40 kang=1,nang
!!
            wk3x = wka1(krng)*cosan(kang)
            wk3y = wka1(krng)*sinan(kang)
!!
            if ( krng.eq.irng ) then              !* wn3 = wn1
!!
!!ba1         Bash; skip k1 but keep the opposite angle to k1 - orig setting
!!ba1               remember here iang = 1
              if ( kang .eq. 1 ) go to 40         !* th3 = th1
!!ba1---
!!            ----------------------------------------------------------
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              call shloxr ( dep,       wk1x,wk1y,wk3x,wk3y )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            it returns: wk2x, wk2y, wk4x, wk4y & ds all dim=(npts)
!!                        and all are PUBLIC
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
            else                                  !* wn3 > wn1
!!
!!            ----------------------------------------------------------
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              call shlocr ( dep,       wk1x,wk1y,wk3x,wk3y )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            it returns: wk2x, wk2y, wk4x, wk4y & ds all dim=(npts)
!!                        and all are PUBLIC
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
            end if   !!  if ( krng.eq.irng )
!!
!!          set the Heaviside coefficient
!b          dif13 = (wk1x-wk3x)**2   +  (wk1y-wk3y)**2        !* wk1y = 0.0
!b          dif13 = (wk1x-wk3x)**2   +       (wk3y)**2
            dif13 = (wk1x-wk3x)*(wk1x-wk3x) + wk3y*wk3y
!!50
            do 50 ipt=1,npts
!!
!!xlc1        Bash; skip k1 but keep the opposite angle to k1 - original setting
              if ( kang.eq.1 ) then                        !* th3=+th1, iang=1
                if (ipt.eq.1 .or. ipt.eq.np2p1) go to 50   !* skip x-axis loci
              end if
!!xlc1---
!!            ----------------------------------------------------------
!!
!!hv          Bash; with !hv ON, move Heaviside section from below to here
!!            Bash moved this section here. *** Check first compute after ***
!!            Skip first then compute only if Heaviside=1, without using it
!!            i.e. compute only if dif13.le.dif14 with Heaviside=1 omitted.
!!            Note; with !hv option is ON, you don't need to turn options
!!            ----  !k19p1 nor !cp4 ON.aYou only need one of the three.
!!            ----------------------------------------------------------
!!            set the Heaviside coefficient
!b            dif14 = (wk1x-wk4x(ipt))**2 + (wk1y-wk4y(ipt))**2 !* wk1y=0.0
!b            dif14 = (wk1x-wk4x(ipt))**2 +      (wk4y(ipt))**2
              dif14 = (wk1x-wk4x(ipt))*(wk1x-wk4x(ipt)) +             &
                                             wk4y(ipt)*wk4y(ipt)
!!
              if ( dif13 .gt. dif14 ) go to 50    !* skip, don't compute
!!
!b            if ( dif13 .gt. dif14 ) then
!b               Heaviside = 0.                   !* Eq(12) of RPTV
!b               go to 50
!b            else
!b               Heaviside = 1.                   !* Eq(11) of RPTV
!b            end if
!!hv---
!!            ----------------------------------------------------------
!!
!!            Set the coupling coefficient for ipt'th locus position
!!            ----------------------------------------------------------
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              call cplshr ( wk4x(ipt),wk4y(ipt), wk3x,wk3y,           &
                            wk2x(ipt),wk2y(ipt), dep, csq,            &
                            irng,krng,kang,ipt )
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            it returns: the coupling coefficient  csq
!!            -- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!            ----------------------------------------------------------
!!
!!wn2         Set parameters related to ipt'th locus wavenumber vector k2
!!            ----------------------------------------------------------
!b            wn2  = sqrt(wk2x(ipt)**2 + wk2y(ipt)**2)               !* k2
              wn2  = sqrt(wk2x(ipt)*wk2x(ipt) + wk2y(ipt)*wk2y(ipt)) !* k2
              th2  = atan2(wk2y(ipt),wk2x(ipt))                !* k2 direction
              if ( th2 .lt. 0. ) th2 = th2 + twopi             !* +ve in radians
              wn2d = wn2*dep                                   !* k2*depth
              tnh2 = tanh(wn2d)                                !* tanh(k2*depth)
              om2  = sqrt(g*wn2*tnh2)                          !* omega2 (rad)
!b            cg2  = 0.5*(om2/wn2)*(1.+wn2d*(1.-tnh2**2)/tnh2) !* group velocity
              cg2  = 0.5*(om2/wn2)*(1.+wn2d*((1./tnh2)-tnh2))  !* group velocity
              f2   = om2/twopi                                 !* f2 (Hz)
!!            ----------------------------------------------------------
!!
!!wn4         Set parameters related to ipt'th locus wavenumber vector k4
!!            ----------------------------------------------------------
!b            wn4  = sqrt(wk4x(ipt)**2 + wk4y(ipt)**2)
              wn4  = sqrt(wk4x(ipt)*wk4x(ipt) + wk4y(ipt)*wk4y(ipt))
              th4  = atan2(wk4y(ipt),wk4x(ipt))
              if ( th4 .lt. 0. ) th4 = th4 + twopi
              wn4d = wn4*dep
              tnh4 = tanh(wn4d)
              om4  = sqrt(g*wn4*tnh4)
!b            cg4  = 0.5*(om4/wn4)*(1.+wn4d*(1.-tnh4**2)/tnh4)
              cg4  = 0.5*(om4/wn4)*(1.+wn4d*((1./tnh4)-tnh4))
              f4   = om4/twopi
!!            ----------------------------------------------------------
!!
!!
!!hv          Bash; with !hv ON, move Heaviside section up
!!            Bash moved this section up. Check first compute after.
!!            ----------------------------------------------------------
!!            set the Heaviside coefficient
!b            dif14 = (wk1x-wk4x(ipt))**2 + (wk1y-wk4y(ipt))**2 !* wk1y=0.0
!b            dif14 = (wk1x-wk4x(ipt))**2 +      (wk4y(ipt))**2
!hv           dif14 = (wk1x-wk4x(ipt))*(wk1x-wk4x(ipt)) +             &
!hv                                          wk4y(ipt)*wk4y(ipt)
!hv           if ( dif13 .gt. dif14 ) then
!hv              Heaviside = 0.                   !* Eq(12) of RPTV
!hv           else
!hv              Heaviside = 1.                   !* Eq(11) of RPTV
!hv           end if
!!hv---
!!            ----------------------------------------------------------
!!
!!
!!            dWdn is the same as sqrt(zzsum) in Don's code, here reduced to a
!!            simpler but mathematically equivalent form that should vary
!!            smoothly between deep and intermediate water owing to identities
!!            using the computer's tanh() function
!!            ----------------------------------------------------------
!!
!!            set grad(,,);
!!            looks like the g^2 goes with csq (Webb'1978, eq. A2)
!!            ----------------------------------------------------------
!!
!b            dWdnsq = cg2**2  - 2.*cg2*cg4 * cos(th2-th4) + cg4**2
              dWdnsq = cg2*cg2 - 2.*cg2*cg4 * cos(th2-th4) + cg4*cg4
!!            ----------------------------------------------------------
!!
              dWdn               = sqrt(dWdnsq)
!!            ----------------------------------------------------------
!!
!!hv          Bash; with !hv ON, don't use var Heaviside (by here it's = 1.0)
              grad(ipt,kang,izz) =           ds(ipt)*csq*gsq/dWdn
!!hv---
!hv           grad(ipt,kang,izz) = Heaviside*ds(ipt)*csq*gsq/dWdn
!!hv---
!!            ----------------------------------------------------------
!!            ==========================================================
!!
!!
!!            Set interpolation, indexing and weight parameters for
!!            computations along wavenumber radials
!!            ----------------------------------------------------------
!!
!!f2          --------------------
              if ( f2 .lt. f0 ) then
                 wtk2(ipt,kang,izz)  = 1.
                 tfac2(ipt,kang,izz) = 0.
                 kref2(ipt,kang,izz) = 1
              else
                 ir = 1+int((alog(f2)-alf0)/aldfrq)
                 if ( ir+1 .gt. nrng ) then
                   wtk2(ipt,kang,izz)  = 0.
                   er = (wka1(nrng)/wn2)**(2.5)
                   tt2= er*(cg2/cgnrng1)*(frqa(nrng)/f2)*(wka1(nrng)/wn2)
                   tfac2(ipt,kang,izz) = tt2
                   kref2(ipt,kang,izz) = nrng - 1
                 else
                   w2 = (f2-frqa(ir))/(frqa(ir+1)-frqa(ir))
                   wtk2(ipt,kang,izz)  = 1. - w2
                   tfac2(ipt,kang,izz) = 1.
                   kref2(ipt,kang,izz) = ir
                 end if
              end if
!!            ----------------------------------------------------------
!!
!!f4          --------------------
              if ( f4 .lt. f0 ) then
                 wtk4(ipt,kang,izz)  = 1.
                 tfac4(ipt,kang,izz) = 0.
                 kref4(ipt,kang,izz) = 1
              else
                 ir = 1+int((alog(f4)-alf0)/aldfrq)
                 if ( ir+1 .gt. nrng ) then
                   wtk4(ipt,kang,izz)  = 0.
                   er = (wka1(nrng)/wn4)**2.5
                   tt4= er*(cg4/cgnrng1)*(frqa(nrng)/f4)*(wka1(nrng)/wn4)
                   tfac4(ipt,kang,izz) = tt4
                   kref4(ipt,kang,izz) = nrng - 1
                 else
                   w2 = (f4-frqa(ir))/(frqa(ir+1)-frqa(ir))
                   wtk4(ipt,kang,izz)  = 1. - w2
                   tfac4(ipt,kang,izz) = 1.
                   kref4(ipt,kang,izz) = ir
                 end if
              end if
!!            ----------------------------------------------------------
!!
!!
!!            Set indexing and weight parameters for computations around
!!            azimuths; it appears that jref2 and jref4 should be bounded
!!            between 0 and nang-1 so that when iang (=1,nang) is added in
!!            the integration section, the proper bin index will arise;
!!            the weights wta2 and wta4 seem to be the fractional bin
!!            widths between th2 or th4 and the next increasing
!!            directional bin boundary
!!            ----------------------------------------------------------
!!
              i = int(th2/ainc)
              wta2(ipt,kang,izz)  = 1. - abs(th2-i*ainc)/ainc
              if ( i .ge. nang )  i = i - nang
              jref2(ipt,kang,izz) = i
!mpc          jref2(ipt,kang,izz) = MOD(i,nang) !* is this better that the above two lines?
!!
              i = int(th4/ainc)
              wta4(ipt,kang,izz)  = 1. - abs(th4-i*ainc)/ainc
              if ( i .ge. nang )  i = i - nang
              jref4(ipt,kang,izz) = i
!mpc          jref4(ipt,kang,izz) = MOD(i,nang) !* is this better that the above two lines?
!!
  50        continue                              !* end of ipt loop
!!
  40      continue                                !* end of kang loop
!!
  30    continue                                  !* end of krng loop
!!
  20  continue                                    !* end of irng loop
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END SUBROUTINE gridsetr
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE shloxr ( dep, wk1x,wk1y, wk3x,wk3y )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: wk2x, wk2y, wk4x, wk4y & ds all dim=(npts)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -------------------------------------------------------------------------#
!!                                                                             !
!!    General locus solution for input vectors (wk1x,wk1y) and                 !
!!    (wk3x,wk3y) of the same magnitude but NOT in the same direction          !
!!    (or singularness will occur), output vectors (wk2x,wk2y) and             !
!!    (wk4x,wk4y), and element length ds along locus curve:                    !
!!                                                                             !
!!    With wavenumber vector n identified by (wknx,wkny), its magnitude        !
!!    given by wkn = sqrt(wknx**2+wkny**2) and its associated radian           !
!!    frequency given by sign = sqrt[g*wkn*tanh(wkn*dep)], where g is          !
!!    gravitational acceleration and dep is water depth, the four-wave         !
!!    resonance condition is satisfied along a locus of pts defined by         !
!!                                                                             !
!!    [1]  (wk1x,wk1y) + (wk2x,wk2y) - (wk3x,wk3y) - (wk4x,wk4y) = 0           !
!!                                                                             !
!!    [2]  sig1 + sig2 - sig3 - sig4 = 0                                       !
!!                                                                             !
!!    In the case where k1 [= sqrt(wk1x**2+wk1y**2)] is equal to k3            !
!!    [= sqrt(wk3x**2+wk3y**2)], we have by dispersion,                        !
!!                                                                             !
!!    [3]  sig1 = sqrt[g*k1*tanh(k1*h)] = sqrt[g*k3*tanh(k3*h)] = sig3         !
!!                                                                             !
!!    so sig1 - sig3 = 0 and [2] becomes sig2 = sig4, where, again by          !
!!    dispersion,                                                              !
!!                                                                             !
!!    [4]  sig2 = sqrt(g*k2*tanh(k2*h)] = sig4 = sqrt(g*k4*tanh(k4*h)]         !
!!                                                                             !
!!    and consequently k2 = k4.  This simplifies the locus solution            !
!!    considerably, and it can be shown that the (wk2x,wk2y) locus is          !
!!    along the perpendicular bisector of the (px,py) vector given by          !
!!                                                                             !
!!    [5]  (px,py) = (wk3x-wk1x,wk3y-wk1y)                                     !
!!                                                                             !
!!    and thereby from [1]                                                     !
!!    [6]  (wk4x,wk4y) = (wk2x,wk2y) - (px,py)                                 !
!!                                                                             !
!!    We note that these loci are independent of depth, although depth         !
!!    is used to set the length of the locus line by requiring that its        !
!!    range on either side of the p vector correspond to a wave with a         !
!!wkx freq four times that of the k1 vector (the locus line is made            !
!!    up of npts segments of length ds; the outer edges of the terminal        !
!!    segments satisfy the length constraint; vectors k2 and k4 extend         !
!!    to segment centers and will sufficiently approximate the length          !
!!    constraint).  As compared to srshlocr.f, we can do all                   !
!!    calculations here in dimensional space.                                  !
!!    -------------------------------------------------------------------------#
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    npts      int.  Public    I   # of points on the locus
!!    np2p1     int.  Public    I   = npts/2 + 1
!!    ----------------------------------------------------------------
!!
!!    *** arrays wk2x,wk2y, wk4x,wk4y & ds are related to locus solutioni
!!        for given vectors K1 & k3        all have dim=(npts)
!!    wk2x      R.A.  Public    O   x_cmp of vector k2 solution on the locus
!!    wk2y      R.A.  Public    O   y_cmp of vector k2 solution on the locus
!!    wk4x      R.A.  Public    O   x_cmp of vector k4 solution on the locus
!!    wk4y      R.A.  Public    O   y_cmp of vector k4 solution on the locus
!!    ds        R.A.  Public    O   element length along the locus
!!                                  (npts*ds circles the whole locus)
!!    ----------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!    ----------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ----------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
!!wvn
!wvn  USE W3DISPMD, ONLY: WAVNU2
!!wvn---
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: dep
      real,    intent(in)  :: wk1x,wk1y, wk3x,wk3y
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    Local Parameters & variables
!!    -----------------------------
      integer              :: m, n
      real                 :: g
      real                 :: wk1, f1, fx
      real                 :: wkx, db, px, py, p, thp, halfp
      real                 :: a, b, dth, a_halfp
!a    real                 :: wkfnc      !* real function or use WAVNU2
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    initial constants
!!    ------------------
      g     = 9.806                            !* set = GRAV as in CONSTANTS
!!
!!ini
!!    initial all returned arrays before they are computed
!!    -----------------------------------------------------
      wk2x(:) = 0.0
      wk2y(:) = 0.0
      wk4x(:) = 0.0
      wk4y(:) = 0.0
      ds(:)   = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!
!!wkx Bash; Try use wkx = wka(nrng) instead of  wkx = wkfnc(4.*f1,dep)
!b    wk1   = sqrt(wk1x**2+wk1y**2)           !* k1=wk1x since wk1y=0.0
      wk1   = wk1x
      f1    = sqrt(g*wk1*tanh(wk1*dep))/twopi !* f1=f(k1,dep)
      fx    = 4. * f1                         !* fx = 4*f1
!!
!!wvn Bash; Try use subr WAVNU2 instead of function wkfnc
!wvn  call WAVNU2(twopi*fx,dep,wkx,cgx)       !* => wkx=k(4*f1,dep) & cgx=Cg(4*f1,dep)
      wkx   = wkfnc(fx,dep)                   !* +> wkx=k(4*f1,dep)  (cgx not needed)
!!wvn---
!!    ------------------------------------------------------------------
!!
      db    = wkx/float(np2p1-1)              !* locus length increment
!!
!!
      px    = wk3x - wk1x
      py    = wk3y - wk1y                     !* wk1y=0.0 can be omitted
!b    p     = sqrt(px**2 + py**2)             !* argument never = 0.0
      p     = sqrt(px*px + py*py)
      thp   = atan2(py,px)
      halfp = 0.5*p
!!
!!
      do 10 n=np2p1,npts                      !* for npts = 30
!!                                            !* n = 16 --> 30
!!
!b       b   = 0.5 * db  + float(n-np2p1) * db
         b   = db * (0.5 + float(n-np2p1))
!b       a   = sqrt(1. + (2.*b/p)**2)
         a   = sqrt(1. + (2.*b/p)*(2.*b/p))
         dth = acos(1./a)
         a_halfp = a*halfp
!!
!b       wk2x(n) = a*halfp * cos(thp+dth)
!b       wk2y(n) = a*halfp * sin(thp+dth)
         wk2x(n) = a_halfp * cos(thp+dth)
         wk2y(n) = a_halfp * sin(thp+dth)
         wk4x(n) = wk2x(n) - px
         wk4y(n) = wk2y(n) - py
         ds(n)   = db
!!
         m       = npts - n + 1                !* m = 15 --> 1
!b       wk2x(m) = a*halfp * cos(thp-dth)
!b       wk2y(m) = a*halfp * sin(thp-dth)
         wk2x(m) = a_halfp * cos(thp-dth)
         wk2y(m) = a_halfp * sin(thp-dth)
         wk4x(m) = wk2x(m) - px
         wk4y(m) = wk2y(m) - py
         ds(m)   = db
!!
  10  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END SUBROUTINE shloxr
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE shlocr ( dep, wk1x,wk1y, wk3x,wk3y )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: wk2x, wk2y, wk4x, wk4y & ds all dim=(npts)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    With wavenumber vector n identified by (wknx,wkny), its magnitude!
!!    given by wkn = sqrt(wknx**2+wkny**2) and its associated radian   !
!!    frequency given by sign = sqrt[g*wkn*tanh(wkn*dep)], where g is  !
!!    gravitational acceleration and dep is water depth, the four-wave !
!!    resonance condition is satisfied along a locus of pts defined by !
!!                                                                     !
!!    [1]  (wk1x,wk1y) + (wk2x,wk2y) - (wk3x,wk3y) - (wk4x,wk4y) = 0   !
!!                                                                     !
!!    [2]  sig1 + sig2 - sig3 - sig4 = 0                               !
!!                                                                     !
!!    Because of the influence of depth, it is convenient to define new!
!!    vectors (wnx,wny) = (wknx*dep,wkny*dep) with magnitudes wn =     !
!!    sqrt(wnx**2+wny**2) = wkn*dep such that a dimensionless frequency!
!!    is sign*sqrt(dep/g) = sqrt[wkn*dep*tanh(wkn*dep)]                !
!!                        = sqrt[wn*tanh(wn)].                         !
!!    With these definitions and vectors (wk1x,wk1y) and (wk3x,wk3y)   !
!!    given as input, we can write (with some rearrangement            !
!!    of [1] and [2])                                                  !
!!                                                                     !
!!    [3]  w3x - w1x = px = w2x - w4x                                  !
!!    [4]  w3y - w1y = py = w2y - w4y                                  !
!!    [5]  sqrt[w3*tanh(w3)] - sqrt[w1*tanh(w1)] = q                   !
!!         = sqrt[w2*tanh(w2)] - sqrt[w4*tanh(w4)]                     !
!!                                                                     !
!!    With dimensionless vector (px,py) = (w3x-w1x,w3y-w1y) [magnitude !
!!    p = sqrt(px**2+py**2), direction atan2(py,px)] and dimensionless !
!!    frequency difference q = sqrt(w3*tanh(w3)] - sqrt(w1*tanh(w1)]   !
!!    defined by input parameters, we see from [3] and [4] that        !
!!    (w4x,w4y) = (w2x-px,w2y-py) [magnitude w4 = sqrt((w2x-px)**2 +   !
!!    (w2y-py)**2)] and thus from [5] we must basically find elements  !
!!    w2x and w2y that satisfy                                         !
!!                                                                     !
!!    [6] sqrt[sqrt(w2x**2+w2y**2)*tanh(sqrt(w2x**2+w2y**2))] -        !
!!        sqrt[sqrt((w2x-px)**2+(w2y-py)**2) *                         !
!!             tanh(sqrt((w2x-px)**2+(w2y-py)**2] = q                  !
!!                                                                     !
!!    The locus curve defined by the set of pts (w2x,w2y) crosses the  !
!!    p-vector axis at two points; one with magnitude w2=rmin*p with   !
!!    0.5 < rmin < 1.0 and one with magnitude w2=rmax*p with rmax > 1. !
!!    We first isolate rmin, rmax using various iterative algorithms,  !
!!    and then find locus pts that are not on the p-vector axis with   !
!!    another iterative scheme.  At the end, we un-normalize the w2    !
!!    and w4 vectors to find the wk2 and wk4 vectors.                  !
!!    -----------------------------------------------------------------#
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    npts      int.  Public    I   # of points on the locus
!!    np2p1     int.  Public    I   = npts/2 + 1
!!    ----------------------------------------------------------------
!!
!!    *** arrays wk2x,wk2y, wk4x,wk4y & ds are related to locus solutioni
!!        for given vectors K1 & k3        all have dim=(npts)
!!    wk2x      R.A.  Public    O   x_cmp of vector k2 solution on the locus
!!    wk2y      R.A.  Public    O   y_cmp of vector k2 solution on the locus
!!    wk4x      R.A.  Public    O   x_cmp of vector k4 solution on the locus
!!    wk4y      R.A.  Public    O   y_cmp of vector k4 solution on the locus
!!    ds        R.A.  Public    O   element length along the locus 
!!                                  (npts*ds circles the whole locus)
!!    ----------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!    ----------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ----------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3ODATMD, ONLY: NDSE
!!
!!
      IMPLICIT NONE
!!    
!!    Parameter list
!!    --------------
      real,    intent(in)  :: dep
      real,    intent(in)  :: wk1x,wk1y, wk3x,wk3y
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!    
!!    Local Parameters & variables
!!    -----------------------------
      integer              :: n, np, nnp, nplace
      integer              :: ierr_gr
!!
      real                 :: p,  px,  py,  q,    qrtp,  qsqp,        &
                              dr, dth, thp, dphi, cphi,               &
                              w1, w1x, w1y, wk1,  w3, w3x, w3y, wk3,  &
                              rold, rold1, rold2, rnew, rnew1, rnew2, &
                              pxod, pyod,  zpod,                      &
                              t, t1, t2, t3, tm, tp, ds1, ds2,        &
                              rmin, rmax, rcenter, rradius
!!
      double precision     :: dbt3,dbt4,dbt5,dbt6, dbz, dbp, dbqrtp,  &
                              cdthold, cdthnew, wate1, wate2
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!ini
!!    initial all returned arrays before they are computed
!!    -----------------------------------------------------
      wk2x(:) = 0.0
      wk2y(:) = 0.0
      wk4x(:) = 0.0
      wk4y(:) = 0.0
      ds(:)   = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!
!!
!b    wk1   = sqrt(wk1x**2 + wk1y**2)         !* k1=wk1x since wk1y=0.0
      wk1   = wk1x
!b    wk3   = sqrt(wk3x**2 + wk3y**2)
      wk3   = sqrt(wk3x*wk3x + wk3y*wk3y)
!!
      w1    = wk1  * dep
      w1x   = wk1x * dep
!b    w1y   = wk1y * dep
      w1y   = wk1y                            !* wk1y=0.0
!!
      w3    = wk3  * dep
      w3x   = wk3x * dep
      w3y   = wk3y * dep
!!
      px    = w3x - w1x
      py    = w3y - w1y
!b    p     = sqrt(px**2 + py**2)             !* argument never = 0.0
      p     = sqrt(px*px + py*py)             !* argument never = 0.0

      thp   = atan2(py,px)
      q     = sqrt(w3*tanh(w3)) - sqrt(w1*tanh(w1))
      qrtp  = q / sqrt(p)
      qsqp  = qrtp*qrtp
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    for (w2x,w2y) = rmin*(px,py) (locus crossing the p-vector axis   !
!!    nearest the origin), we have (w4x,w4y) = (w2x,w2y) - (px,py)     !
!!    = rmin*(px,py) - (px,py) = (rmin - 1)*(px,py); note that because !
!!    rmin < 1, the length of (w4x,w4y) is w4 = (1 - rmin)*p; then [6] !
!!    takes the simpler form                                           !
!!                                                                     !
!!    [7] sqrt(rmin*p*tanh(rmin*p)) - sqrt[(1-rmin)*p*tanh((1-rmin)*p)]!
!!        = q                                                          !
!!                                                                     !
!!    assuming the tanh() functions are slowly varying and can be      !
!!    treated as separate entities, [7] can be written as a quadratic  !
!!    in sqrt(rmin), i.e.,                                             !
!!                                                                     !
!!               2*qrtp*sqrt(rmin*p)                                   !
!!    [8] rmin - ------------------- sqrt(rmin) +                      !
!!                        T                                            !
!!                                                                     !
!!                                   qsqp-p*tanh((1-rmin)*p)           !
!!                                   -----------------------  =  0,    !
!!                                              T                      !
!!                                                                     !
!!    where   T    = tanh(rmin*p) + tanh((1-rmin)*p),                  !
!!            qrtp = q/sqrt(p)  and qsqp = qrtp**2                     !
!!                                                                     !
!!    the square of the most positive root of [8] can (with some       !
!!    algebra) be written                                              !
!!                                                                     !
!!    [9] rmin =                                                       !
!!   (1/T**2)*{qsqp*[tanh(rmin*p)-tanh((1-rmin)*p)]+T*tanh((1-rmin)*p) !
!!           +2*qrtp*sqrt[tanh(rmin*p)*tanh((1-rmin)*p)]*sqrt(T-qsqp)} !
!!                                                                     !
!!    setting rnew=rmin on the LHS and rold=rmin on the RHS in all     !
!!    instances allows the creation of an iterative algorithm for rmin;!
!!    convergence can be slow in general, so a coarse search for the   !
!!    crossing of [9] with the rnew=rold line is conducted first, then !
!!    a weighted iterative replacement loop is executed until the      !
!!    desired accuracy is achieved;                                    !
!!    Note that if p is sufficiently large, all tanh() -> 1            !
!!    and [9] becomes the analytic expression                          !
!!            rmin = 0.5*[1 + qrtp*sqrt(2-qsqp)]                       !
!!                                                                     !
!!    following is the coarse search using rold1 = 0.5,0.1,1.0,        !
!!                                         rold2 = rold1 + 0.1:        !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
      ierr_gr = 0
!!
      rold1 = 0.5
      tp    = tanh(rold1 * p)
      tm    = tanh((1.-rold1) * p)
      t     = tp + tm
      t1    = qsqp * (tp-tm)
      t2    = t  * tm
      t3    = 2. * qrtp * sqrt(tp*tm) * sqrt(t-qsqp)
      rnew1 = (t1 + t2 + t3) / (t*t)
!!
!!
      do 10 n=1,4
        rold2 = rold1 + 0.1
        tp    = tanh(rold2 * p)
        tm    = tanh((1.-rold2) * p)
        t     = tp + tm
        t1    = qsqp * (tp-tm)
        t2    = t  * tm
        t3    = 2. * qrtp * sqrt(tp*tm) * sqrt(t-qsqp)
        rnew2 = (t1 + t2 + t3) / (t*t)
        if ( rnew2 .lt. rold2 ) then
          rold = (rold2*rnew1-rold1*rnew2)/(rold2-rold1-rnew2+rnew1)
          go to 11
        end if
        rold1 = rold2
        rnew1 = rnew2
  10  continue
      rold = 0.9                       !* default if not otherwise found
  11  continue
!!    ------------------------------------------------------------------
!!
!!
!!    iterative replacement search for rmin
      do 20 n=1,50
        tp   = tanh(rold * p)
        tm   = tanh((1.-rold) * p)
        t    = tp + tm
        t1   = qsqp * (tp-tm)
        t2   = t  * tm
        t3   = 2. * qrtp*sqrt(tp*tm)*sqrt(t-qsqp)
        rnew = (t1 + t2 + t3) / (t*t)
        if ( abs(rnew-rold) .lt. 0.00001 ) then
          rmin = rnew
          go to 21
        end if
        rold = 0.5 * (rold + rnew)
  20  continue
      ierr_gr = ierr_gr + 1  !* set 1's flag in ierr_gr if no convergence
      rmin = rnew
  21  continue
!!    ------------------------------------------------------------------
!!
!!    set (dimensional) wavenumber components for this point on locus
      wk2x(1) =  rmin * px / dep
      wk2y(1) =  rmin * py / dep
      wk4x(1) = (rmin-1.) * px / dep
      wk4y(1) = (rmin-1.) * py / dep
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    for (w2x,w2y) = rmax*(px,py) (locus crossing the p-vector axis   !
!!    farthest from the origin), we have (w4x,w4y)=(w2x,w2y) - (px,py) !
!!    = rmax*(px,py) - (px,py) = (rmax - 1)*(px,py);                   !
!!    here, because rmax > 1, the length of (w4x,w4y) is               !
!!    w4 = (rmax - 1)*p; then [6] takes the form                       !
!!                                                                     !
!!    [10] sqrt(rmax*p*tanh(rmax*p))-sqrt[(rmax-1)*p*tanh((rmax-1)*p)] !
!!         = q                                                         !
!!                                                                     !
!!    rearranging terms, squaring both sides and again rearranging     !
!!    terms yields                                                     !
!!                                                                     !
!!    [11] rmax*p*[tanh(rmax*p) - tanh((rmax-1)*p)]                    !
!!         = 2*q*sqrt(tanh(rmax*p))*sqrt(rmax*p) +                     !
!!           q**2 + p*tanh((rmax-1)*p)                                 !
!!                                                                     !
!!    because the difference of the two tanh()'s on the LHS tend to    !
!!    make the whole term small, we solve for rmax from the rapidly    !
!!    varying part of the first term on the RHS, i.e.,                 !
!!                                                                     !
!!                         [tanh((rmax-1)*p)+qsqp + rmax*T]**2         !
!!    [12]          rmax = ----------------------------------- ,       !
!!                                 4*qsqp*tanh(rmax*p)                 !
!!                                                                     !
!!    where, in this algorithm, T = tanh(rmax*p) - tanh((rmax-1)*p);   !
!!    as for rmin in [9], setting rnew=rmax on the LHS and rold=rmax   !
!!    in all instances on the RHS allows the formation of an iterative !
!!    algorithm; initially, we only know rmax > 1 so we do a coarse    !
!!    search in the 10's place out to some reasonably big number to    !
!!    try to find the place where [12] crosses the rnew = rold line    !
!!    (if this fails, we set an error flag); in refining the estimate, !
!!    it appears that [12] can get a little squirrely, so we do a      ! 
!!    brute force successive decimation search to nplace decimal places!
!!    to home in on the answer; note that if p is big enough for the   !
!!    tanh()'s to reach unity, [12] becomes exact and                  !
!!    rmax = [(1 + qsqp)**2]/(4*qsqp)                                  !
!!    following is the coarse search with rold = 1,10,2001:            !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
      rold = 1.0
      do 30 n=1,200
        rold = rold + 10.
        tp   = tanh(rold * p)
        tm   = tanh((rold-1.) * p)
        t    = tp - tm
        t1   = tm + qsqp
        t2   = 4. * tp * qsqp
        rnew = ((t1+rold*t)**2) / t2
        if ( rnew .lt. rold ) then
          rold = rold - 10.
          go to 31
        end if
  30  continue
      ierr_gr = ierr_gr + 10    !* set 10's place in ierr_gr if no sol'n
  31  continue
!!    ------------------------------------------------------------------
!!
!!
!!    successive decimation search to refine rmax
      dr = 10.
      do 40 nplace=1,6
        dr = dr/10.
        do 50 n=1,10
          rold = rold + dr
          tp   = tanh(rold * p)
          tm   = tanh((rold-1.) * p)
          t    = tp - tm
          t1   = tm + qsqp
          t2   = 4. * tp * qsqp
          rnew = ((t1+rold*t)**2) / t2
          if ( rnew .lt. rold ) then
            rold = rold - dr
            go to 51
          end if
  50    continue
  51    continue
  40  continue
!!
      rmax = rold
!!
!!    set (dimensional) wavenumber components for this locus point
!!
      wk2x(np2p1) =  rmax * px / dep
      wk2y(np2p1) =  rmax * py / dep
      wk4x(np2p1) = (rmax-1.) * px / dep
      wk4y(np2p1) = (rmax-1.) * py / dep
!!
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    search for cos(dth) for off-p-vector solutions; use a circle     !
!!    centered on the p-vector axis at a distance                      ! 
!!    rcenter = 0.5*(rmax+rmin)  from the  origin with a               !
!!    radius  = 0.5*(rmax-rmin); radii from the center of the circle   !
!!    at successive angle increments np*dphi intersect the circle at   !
!!    distances r*p from the origin of the p vector such that          !
!!                                                                     !
!!    [13]  r**2 = rradius**2 + rcenter**2 -                           !
!!                 2*rcenter*rradius*cos(np*dphi)                      !
!!                                                                     !
!!    and makes an angle dth with the p vector satisfying              !
!!                                                                     !
!!    [14]   cdth = cos(dth) = (rcenter/r) - (rradius/r)*cos(np*dphi)  !
!!                                                                     !
!!    we then rotate this vector, holding its length=r*p constant and  !
!!    successively estimating cdth (using the above equation as an     !
!!    initial guess) until it intersects the locus curve; some algebra !
!!    yields the estimation equation as                                !
!!                                                                     !
!!                  r**2 + 1      [sqrt(r*tanh(rp)) - q/sqrt(p)]**4    !
!!    [15] cdthnew= ------- - ---------------------------------------- !
!!                    2*r     2*r*[tanh(p*sqrt(r**2-2*r*cdthold+1))]**2!
!!                                                                     !
!!    we use a weighted new estimate of cdthold with the weights based !
!!    on the argument of the tanh() function in the denominator        !
!!    (if the argument is big, tanh() -> 1 and cdthnew is found in one !
!!    pass; for small arguments, convergence is faster with equal      !
!!    weighting of old and new estimates; all this is empirical to try !
!!    to increase speed); double precision is used to gain enough      !
!!    accuracy when the arccos is taken; note that if p is big enough  !
!!    for all tanh()'s -> 1, [15] is exact and                         !
!!    cdthnew = cdth = [r**2 + 1 - (sqrt(r) - qrtp)**4]/(2*r)          !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!!
      rcenter = 0.5 * (rmax + rmin)
      rradius = 0.5 * (rmax - rmin)
      t1      = rradius**2 + rcenter**2
      t2      = 2. * rradius * rcenter
      dphi    = 6.283185308 / float(npts)
      pxod    = px / dep
      pyod    = py / dep
!!
      dbp     = dble(p)
      dbqrtp  = dble(qrtp)
!!
!!
      do 60 np=2,npts/2                        !* np = 2 --> 15
!!
        cphi    = cos(float(np-1)*dphi)
        dbz     = dsqrt(dble(t1-t2*cphi))
        cdthold = dble(rcenter-rradius*cphi) / dbz
        dbt3    = (dbz*dbz) + 1.d0
        dbt4    = dbt3 / (2.d0*dbz)
        dbt5    = ((dsqrt(dbz*dtanh(dbz*dbp))-dbqrtp)**4)/(2.d0*dbz)
        dbt6    = dbp * dsqrt(dbt3-2.d0*dbz*cdthold)
!!
        if ( dbt6 .gt. 0.55d0 ) then
          wate1 = dtanh(dbt6)
          wate2 = 1.d0 - wate1
        else
          wate1 = 0.5d0
          wate2 = 0.5d0
        end if
!!
        do 70 n=1,25
          cdthnew = dbt4 - dbt5 / ((dtanh(dbt6))**2)
          if ( dabs(cdthnew-cdthold) .lt. 0.0000001d0 ) go to 71
          cdthold = wate1 * cdthnew + wate2 * cdthold
          dbt6    = dbp * dsqrt(dbt3-2.d0*dbz*cdthold)
  70    continue
        ierr_gr = ierr_gr + 100   !* add to 100's place for every failure
  71    continue
!!
        dth  = sngl(dacos(cdthnew))
        zpod = sngl(dbz) * p / dep
!!
        wk2x(np) = zpod * cos(thp+dth)
        wk2y(np) = zpod * sin(thp+dth)
        wk4x(np) = wk2x(np) - pxod
        wk4y(np) = wk2y(np) - pyod
!!
        nnp = npts-np+2                        !* for npts = 30
!!                                             !* nnp = 30 --> 17
!!
        wk2x(nnp) = zpod * cos(thp-dth)
        wk2y(nnp) = zpod * sin(thp-dth)
        wk4x(nnp) = wk2x(nnp) - pxod
        wk4y(nnp) = wk2y(nnp) - pyod
!!
  60  continue
!!    ------------------------------------------------------------------
!!
!!
      if ( ierr_gr .ne. 0 ) then
        write ( ndse,1000 ) ierr_gr
        CALL EXTCDE ( 60 ) 
      endif
!!    ------------------------------------------------------------------
!!
!!
!!    set arc length ds as the sum of half the segment lengths on either
!!    side of a given point
!!
      ds1   = sqrt((wk2x(2)-wk2x(1))**2+(wk2y(2)-wk2y(1))**2)
      ds(1) = ds1
      do 80 np=3,npts/2+1
        ds2   = sqrt((wk2x(np)-wk2x(np-1))**2+(wk2y(np)-wk2y(np-1))**2)
        ds(np-1)      = 0.5*(ds1+ds2)
        ds(npts-np+3) = ds(np-1)
        ds1           = ds2
  80  continue
      ds(npts/2+1)    = ds2
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
 1000 format ( ' W3SNL4 Error : In shlocr. Error from gridset ',i10)
!!
      END SUBROUTINE shlocr
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE cplshr ( w1x0,w1y0, w2x0,w2y0, w3x0,w3y0,            &
                               h, csq, irng,krng, kang,ipt )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!
!!    it returns: the coupling coefficient  csq
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Calculates four-wave Boltzmann coupling coefficient in shallow   !
!!    water given k1,k2,k3 and following at least Hasselmann (1962)    !
!!    and probably Herterich and Hasselmann (1982).  Dimensional       !
!!    wavenumbers are (wnx0,wny0), n = 1,3, h = depth, csq = coupling  !
!!    coefficient.  This is the same as Don's cplesh, except within    !
!!    the algorithm, wavenumbers are made dimensionless with h and     !
!!    frequencies with sqrt(h/g), g = gravitational acceleration (the  !
!!    idea is to simplify and speed up the calculations while keeping  !
!!    a reasonable machine resolution of the result).  At the end,     !
!!    dimensionless csqhat is redimensioned as csq = csqhat/(h**6)     !
!!    so it is returned as a dimensional entity.                       !
!!                                                                     !
!!    This calculation can be a touchy bird, so we use double precision!
!!    for internal calculations, using single precision for input and  !
!!    output.                                                          !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ------------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      integer, intent(in)  :: irng,krng, kang,ipt
      real,    intent(in)  :: w1x0,w1y0, w2x0,w2y0, w3x0,w3y0
      real,    intent(in)  :: h          !* depth 'dep'
      real,    intent(out) :: csq
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    Local Parameters & variables
!!    -----------------------------
      integer              :: ipass
      double precision     :: hh
      double precision     :: s1,    s2,    s3
      double precision     :: k1x,   k2x,   k3x
      double precision     :: k1y,   k2y,   k3y
      double precision     :: k1,    k2,    k3
      double precision     :: om1,   om2,   om3
      double precision     :: som1,  som2,  som3
      double precision     :: om1sq, om2sq, om3sq
      double precision     :: k23,   k23x,  k23y
      double precision     :: dot23, dot123
      double precision     :: omsq23
!!mpc
      double precision     :: k1sq, k2sq, k3sq, k23sq
      double precision     :: tanh_k1, tanh_k2, tanh_k3, tanh_k23
!!mpc---
      double precision     :: k1x0,  k2x0,  k3x0,  k1zx
      double precision     :: k1y0,  k2y0,  k3y0,  k1zy
      double precision     :: di,    e
      double precision     :: p1,    p2,    p3,    p4
      double precision     :: t1,    t2,    t3,    t4,    t5
!!
      double precision     :: csqhatd, csqd
      double precision     :: scple
      double precision     :: pi4
!!
!!eps Bash; added +eps to avoid dividing by 0.0. Dividing by 0.0 causes NaN
!eps  double precision     :: eps
!!    
!!    Bash; Added domsq23 = denominator of  t1      in cplshr
!!          and   sumom   = denominator of  csqhatd in cplshr
      double precision     :: domsq23, sumom
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    initial constants
!!    ------------------
      hh    = dble(h)                     !* single to dbl precision
      pi4   = 0.785398175d0               !* Set = PI/4 as in CONSTANTS
!eps  eps   = 1.d-12                      !* set eps to a very small number
      scple = 0.d0                        !* initialize accumulator
!!ini
!!    initialize returned variable 'csq'
!!    ----------------------------------
      csq   = 0.d0
!!ini---
!!    ------------------------------------------------------------------
!!
      do 10 ipass=1,3
!p1
        if (ipass .eq. 1) then            !* initial pass (+1,+1,-1)
          s1   =  1.d0
          s2   =  1.d0
          s3   = -1.d0
          k1x0 = dble(w1x0) * hh          !* norm. k elements with h
          k1y0 = dble(w1y0) * hh
          k2x0 = dble(w2x0) * hh
          k2y0 = dble(w2y0) * hh
          k3x0 = dble(w3x0) * hh
          k3y0 = dble(w3y0) * hh
!p1
!p2
        else if (ipass .eq. 2) then       !* 1st permutation (+1,-1,+1)
          s1   =  1.d0
          s2   = -1.d0
          s3   =  1.d0
          k1zx = k1x0
          k1zy = k1y0
          k1x0 = k2x0
          k1y0 = k2y0
          k2x0 = k3x0
          k2y0 = k3y0
          k3x0 = k1zx
          k3y0 = k1zy
!p2
!p3
        else                              !* 2nd permutation (-1,+1,+1)
          s1   = -1.d0
          s2   =  1.d0
          s3   =  1.d0
          k1zx = k1x0
          k1zy = k1y0
          k1x0 = k2x0
          k1y0 = k2y0
          k2x0 = k3x0
          k2y0 = k3y0
          k3x0 = k1zx
          k3y0 = k1zy
!p3
        end if
!!k19p1
!!k19p1 Note: na2p1=nang/2+1   !* this is the angle index opposite to iang=1
!k19p1  if (krng.ne.irng .and. kang.eq.na2p1 .and. ipt.eq.1 .and.     &
!k19p1                                        ipass.eq.1) go to 10
!!k19p1---
!!
        k1x = s1 * k1x0                    !* sign the norm'ed k parts
        k1y = s1 * k1y0
        k2x = s2 * k2x0
        k2y = s2 * k2y0
        k3x = s3 * k3x0
        k3y = s3 * k3y0
!!mpc
!mpc    k1    = dsqrt(k1x**2 + k1y**2)     !* normalized |k|
!mpc    k2    = dsqrt(k2x**2 + k2y**2)
!mpc    k3    = dsqrt(k3x**2 + k3y**2)
!!mpc---
        k1sq  = (k1x*k1x + k1y*k1y)        !* normalized |k| **2
        k2sq  = (k2x*k2x + k2y*k2y)
        k3sq  = (k3x*k3x + k3y*k3y)
        k1    = dsqrt(k1sq)                !* normalized |k|
        k2    = dsqrt(k2sq)
        k3    = dsqrt(k3sq)
!!mpc---
!!
!!mpc
!mpc    om1   = dsqrt(k1*dtanh(k1))        !* norm. omega (by sqrt(h/g))
!mpc    om2   = dsqrt(k2*dtanh(k2))
!mpc    om3   = dsqrt(k3*dtanh(k3))
!mpc    om1sq = om1**2
!mpc    om2sq = om2**2
!mpc    om3sq = om3**2
!!mpc---
        tanh_k1 = dtanh(k1)
        tanh_k2 = dtanh(k2)
        tanh_k3 = dtanh(k3)
        om1sq   = k1*tanh_k1
        om2sq   = k2*tanh_k2
        om3sq   = k3*tanh_k3
        om1     = dsqrt(om1sq)             !* norm. omega (by sqrt(h/g))
        om2     = dsqrt(om2sq)
        om3     = dsqrt(om3sq)
!!mpc---
!!
        som1 = s1 * om1                    !* sign the norm'ed omega's
        som2 = s2 * om2
        som3 = s3 * om3
!!      ----------------------------------------------------------------
!!      ================================================================
!!
!!
        dot23  = k2x*k3x + k2y*k3y        !*  vector k2 dot vector k3
!!
        k23x   = k2x + k3x                !* (vector k2  +  vector k3)_x
        k23y   = k2y + k3y                !* (vector k2  +  vector k3)_y
!!
!!mpc
!mpc    k23    = dsqrt(k23x**2+k23y**2)   !* |vector k2  +  vector k3|
!!mpc---
        k23sq  = (k23x*k23x + k23y*k23y)
        k23    = dsqrt(k23sq)             !* |vector k2  +  vector k3|
!!mpc---
!!
!!mpc
!mpc    omsq23 = k23 * dtanh(k23)         !* norm sq frq of v.k2+v.k3
!!mpc---
        tanh_k23 = dtanh(k23)
        omsq23   = k23 * tanh_k23         !* norm sq frq of v.k2+v.k3
!!mpc---
!!
        dot123 = k1x*k23x + k1y*k23y      !* v.k1 dot (v.k2 + v.k3)
!!      ----------------------------------------------------------------
!!
!!      note: the "i**2" factor from some reference is included in this term
!!
!!mpc
!mpc    di = -(som2+som3)*(om2sq*om3sq-dot23)+0.5d0 *                 &
!mpc          (som2*(k3**2-om3sq**2)+som3*(k2**2-om2sq**2))
!!mpc---
        di = -(som2+som3)*(om2sq*om3sq-dot23)+0.5d0 *                 &
              (som2*(k3sq-om3sq*om3sq)+som3*(k2sq-om2sq*om2sq))
!!mpc---
!!
        e  = 0.5d0*(dot23-som2*som3*(om2sq+om3sq+som2*som3))
!!
        p1 = 2.d0 * (som1+som2+som3) * (om1sq*omsq23 - dot123)
!!
!!mpc
!mpc    p2 = -som1 * (k23**2 - omsq23**2)
!mpc    p3 = -(som2+som3) * (k1**2 - om1sq**2)
!mpc    p4 = k1**2 - om1sq**2
!!mpc---
!!      equation  p2 rewritten to preserve numerical precision
!!      equations p3, p4 rearranged to avoid recomputations.
        p2 = -som1 * (k23sq*(1 - tanh_k23*tanh_k23))
        p4 = (k1sq*(1-tanh_k1*tanh_k1))
        p3 = -(som2+som3) * p4
!!mpc---
!!      ----------------------------------------------------------------
!!
!!      Bash; added & used  variable domsq23 = denominator of t1
        domsq23 = omsq23 - ((som2+som3)**2)  !* Bash; needed for test below
!!      ----------------------------------------------------------------
!!
!!cp4   Bash; with !cp4 ON, test if ( domsq23 .eq. 0.d0 )
!cp4    if ( domsq23 .eq. 0.d0 ) then     !* Bash; this test was needed
!!                                        !* when !k19p1 & !hv were OFF
!!        domsq23=0.0 Dividing by 0.0 causes NaN; here we avoid it
!cp4      t1 = 0.d0
!eps      t1 = di * (p1+p2+p3) / (domsq23+eps)    !* Add eps to denominator
!!                                                !* and may be to numerator
!cp4    endif
!!cp4---
!!      Bash; with !cp4 OFF, don't test if ( domsq23 .eq. 0.d0 )
!!      domsq23 is not = 0.0 (when !k19p1 & !hv were OFF)
!b      t1 = di * (p1+p2+p3) / (omsq23 - ((som2+som3)**2))
        t1 = di * (p1+p2+p3) / (domsq23)
!!cp4---
!!      ----------------------------------------------------------------
!!
        t2 = -di * som1 * (om1sq+omsq23)
        t3 = e * ((som1**3) * (som2+som3) - dot123 - p4)
        t4 = 0.5d0 * som1 * dot23 *                                   &
          ((som1+som2+som3) * (om2sq+om3sq) + som2*som3*(som2+som3))
!!
!!mpc
!mpc    t5 = -0.5d0 * som1 *                                          &
!mpc        (om2sq * (k3**2) * (som1+som2 + 2.d0 * som3)  +           &
!mpc         om3sq * (k2**2) * (som1+som3 + 2.d0 * som2))
!!mpc---
        t5 = -0.5d0 * som1 *                                          &
            (om2sq * (k3sq) * (som1+som2 + 2.d0 * som3)  +            &
             om3sq * (k2sq) * (som1+som3 + 2.d0 * som2))
!!mpc---
!!
        scple = scple + t1 + t2 + t3 + t4 + t5
!!
  10  continue  !! end do 10 ipass=1,3
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!as  HH did division by 3 after adding 3 terms
!as   scple = scple/3.d0
!!as---
!!
!!    Bash; Added sumom   = denominator of  csqhatd in cplshr
      sumom = om1*om2*om3*(om2+om3-om1)
!b    csqhatd = scple*scple*pi4/(om1*om2*om3*(om2+om3-om1))  !* Bash; ok
      csqhatd = scple*scple*pi4/(sumom)
!!    ------------------------------------------------------------------
!!
      csqd    = csqhatd / (hh**6)
      csq     = sngl(csqd)                 !* from dbl to single precision
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END SUBROUTINE cplshr
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE optsa2 ( nrmn,nrmx,    npk,fpk, nbins, wka, cga )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!    ------------------------------------------------------------------
!!
!!    It returns variables dens1(nrng,nang) and dens2(nrng,nang)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    Splits the Action Density into two parts:
!!    (1) large-scale part dens1(nrng,nang)  and
!!    (2) small-scale part dens2(nrng,nang)
!!    dens1 & dens2 in Polar Action Density (k,theta) space Norm. (in k)
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!op2
!!    nrmn      int.  Local     I   number of first freq. bin in [1,nrng-1]
!!    nrmx      int.  Local     I   number of last  freq. bin in [2,nrng]
!!    npk       int.  Local     I   number of peak frequency  in [2,nrng-1]
!!    nbins     int.  Local     I   actual # of bins > npk  (incl. nfs)  or
!!                                  actual # of bins > npk2 (incl. nrng)
!!                                  to guarantee a min 1 bin in equi. range
!!                                  (see subr. W3SNL4)
!!    ------------------------------------------------------------ !!op2
!!
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles
!!
!!    dfrq      Real  Public    I   freq mult. for log freq spacing
!!    fpk       Real  Public    I   peak freq. [Hz] of initial freq spectrum
!!    oma       R.A.  Public    I   rel. freq. array (rad*Hz) ----- dim=(nrng)
!!    frqa      R.A.  Public    I   radian frequencies (Hz) ------- dim=(nrng)
!!
!!    ainc      Real  Public    I   angle increment (radians)
!!    angl      R.A.  Public    I   dir. array (rad) (full circle); dim=(nrng)
!!    cosan     R.A.  Public    I   cosine angles array ----------- dim=(nang)
!!    sinan     R.A.  Public    I   sine   angles array ----------- dim=(nang)
!!    ------------------------------------------------------------------
!!
!!    wka       R.A.  Local     I   wavenumbers array [1/m] ------- dim=(nrng)
!!    cga       R.A.  Local     I   group velocities array [m/s] -- dim=(nrng)
!!                                  wka & cga arrays are corrsp. to depth 'dep'
!!    ------------------------------------------------------------------
!!
!!    ef2       R.A.  Public    I   2D Energy Density spectrum ef2(theta,f)
!!                                  = A(theta,k)*2*pi*oma(f)/cga(f) dim=(nrng,nang)
!!    ef1       R.A.  Public    I   1D Energy Density spectrum ef1(f) dim=(nrng)
!!    ------------------------------------------------------------------
!!
!!    dens1     R.A.  Public    O   large-scale Action Density (k,theta)
!!                                                     dim=(nrng,nang)
!!    dens2     R.A.  Public    O   Small-scale Action Density (k,theta)
!!                                                     dim=(nrng,nang)
!!    ------------------------------------------------------------------
!! 
!! 4. Subroutines used :
!!     
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!! 
!! 5. Called by :
!!     
!!     Name      Type  Module   Description
!!    ------------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ------------------------------------------------------------------
!! 
!! 6. Error messages :
!!      
!!      None.
!! 
!! 7. Remarks :
!! 
!! 8. Structure :
!!    
!!    See source code.
!! 
!! 9. Switches :
!!    
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!    
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
!!
      IMPLICIT NONE
!!
!!
!!
!!    Parameter list
!!    --------------
!!op2 Bash; new for optsa2
      integer, intent(in)  :: nrmn, nrmx, nbins
!!    ------------------------------------------------------------ !!op2
!!
      integer, intent(in)  :: npk
      real,    intent(in)  :: fpk
      real,    intent(in)  :: wka(nrng),   cga(nrng)
!!    ------------------------------------------------------------------
!!
!!
!!    Local Parameters & variables
!!    -----------------------------
      integer              :: irng, iang
!!
!!p2
!!    Bash; Uses of original "psi2(:)", it was very bad (see below)
!p2   integer              :: n1,  n2,  m,   mm
!p2   integer              :: nn1, nn2, ii, idif
!p2   real                 :: q(16)
!p2   real                 :: emax
!p2   real                 :: y, qmin, adif
!!p2---
!!
!!p3
!!    Bash; This is an attempt to replace the original psi2(:)
!!          with a distr. based on sin()**mm  with 'newmaxang'
!!          - not good enough (see below)
!!          !!p4 is an override of !!p3 with mm=4
!p3   integer              :: n1,  n2,  m,   mm
!p3   real                 :: q(16)
!p3   real                 :: y, qmin, adif
!!    The var. below are needed to find 'newmaxang' used in !p3 & !p4
      integer              :: maxang,    newmaxang
      integer              :: maxangshift
      integer              :: halfangl,  halfangu
      real                 :: ef2maxrow(nang)
      real                 :: ef2shift(nang)
      real                 :: halfmax
!!p3---
!!
!!p4
!!    Bash; !!p4 is an override of !!p3 with mm=4
      integer              :: n1,  n2
      real                 :: q4
!!p4---
!!
!!eq
!!    Bash; Use variable equi. range suitable to TSA min condition
!!          simplifed to one point equi. range nearest to 2.*fp
!eq   integer              :: neq
!eq   real                 :: fovfp
!!eq---
!!
      integer              :: igam
      real                 :: sum1, fac
      real                 :: beta, gam
      real                 :: fdenp, fr, ratio, z, ddd
      real                 :: sigz               !* sigz  = 0.109
      real                 :: fk(nrng),    fknrm(nrng)
      real                 :: bscl1(nrng), fkscl1(nrng)
      real                 :: psi2(nang)
      real                 :: act2d(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!ini
!!    Bash; initialize psi2() array here
      psi2(:)    = 0.0
!!
!!    Initialize all the 1d & 2d arrays that are being used
!!    and especially those that are being returned
      fk(:)      = 0.0
      fknrm(:)   = 0.0
      fkscl1(:)  = 0.0
      bscl1(:)   = 0.0
      act2d(:,:) = 0.0
      dens1(:,:) = 0.0
      dens2(:,:) = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!
!!
!!*   Convert 2D Energy Density ef2(f,theta)
!!         to 2D Polar Action Density act2d(k,theta) Norm. (in k)
      do 25 irng=nrmn,nrmx
        fac   = cga(irng)/(twopi*oma(irng)*wka(irng))
        do 24 iang=1,nang
          act2d(irng,iang) = ef2(irng,iang) * fac
  24    continue
!!
!!*     Convert ef1(f) to fk(k); both are 1d Energy Density
        fk(irng)    = cga(irng)*ef1(irng)/twopi  !* fk(k) energy
!!
!!*     Normalize the 1d wavenumber Energy Density fk(k) to give fknrm(k)
        fknrm(irng) = fk(irng)*wka(irng)**2.5    !* fknrm(k) = Norm. fk(k)
  25  continue
!!    ------------------------------------------------------------------
!!
!!
!!    Fit parameters to spectrum
!!    --------------------------
!!eq
!eq   sum1 = 0.
!eq   neq  = 0
!eq   do 26 irng=nrmn,nrmx
!eq     fovfp = frqa(irng)/fpk
!!      Bash; check2 test equilibrium range
!b      if ( fovfp.ge.1.55.and.fovfp.le.2.45 ) then !* orig   equi range
!b      if ( fovfp.ge.1.20.and.fovfp.le.2.20 ) then !* wide   equi range
!b      if ( fovfp.ge.1.90.and.fovfp.le.2.20 ) then !* narrow equi range
!!      --------------------------------------------------------------!*  <<<<<
!!      Bash; select variable equi. range suitable to TSA min condition
!eq     if ( fovfp.ge.(dfrq**(nbins))-0.005 .and.                     &
!eq          fovfp.le.(dfrq**(nbins))+0.005 ) then  !* narrow equi range  <<<<<
!!      --------------------------------------------------------------!*  <<<<<
!eq       sum1 = sum1 + fknrm(irng)
!eq       neq  = neq + 1
!eq     endif
! 26  continue
!eq   beta = sum1 / neq
!eq   gam  = fknrm(npk) / beta
!!eq---
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    Simplify
      beta = fknrm(npk+nbins)
      gam  = fknrm(npk) / beta
!!eq---
!!    ------------------------------------------------------------------
!!
      do 226 irng=nrmn,nrmx
         fknrm(irng) = fknrm(irng) / beta
 226  continue
!!    ==================================================================
!!
!!
!!p2
!!    Construct Directional Distribution "psi2(:)" - original option
!!    ------------------------------------------------------------------
!!
!!    Solve for Normalizing Coefficient for Integral [1.0/(cos**m)]
!!    Note: n1, n2 spans half circle (from -pi/2 to +pi/2 going through 0.)
!p2   n1 = -nang/4 + 1
!p2   n2 =  nang/4 + 1
!p2   do 16 m=1,16
!p2     sum1 = 0.
!p2     do 15 iang=n1,n2
!p2       ii = iang
!p2       if ( iang .lt. 1 ) ii = iang + nang
!p2       sum1 = sum1 + cosan(ii)**m
! 15    continue
!p2     q(m) = 1./(sum1*ainc)
! 16  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Find peak direction "maxang" in ef2() at "npk" the peak in ef1()
!!    needed to define the energy spreading factor y=ef2(npk,maxang)/ef1(npk)
!!    Bash; Note; This original "psi2(:)" was simply very bad because the drift
!!          in "maxang" location causing the 2D Snl to lose symmetry
!p2   emax   = 0.
!p2   maxang = 0
!p2   do 27 iang=1,nang
!p2     if ( ef2(npk,iang).gt.emax ) then
!p2       emax   = ef2(npk,iang)
!p2       maxang = iang                          !* in [1,nang]
!p2     endif
! 27  continue
!p2   y  = ef2(npk,maxang)/ef1(npk)              !* Bash; Energy Spread
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Compare value of peak with q-array for closest fit to cos()**m at peak
!p2   mm   = 1
!p2   qmin = abs(q(1)-y)
!p2   do 28 m=2,16
!p2     adif = abs(q(m)-y)
!p2     if ( adif.lt.qmin ) then
!p2       qmin = adif
!p2       mm   = m
!p2     endif
! 28  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!p2   nn1 = maxang - nang/4           !* nn1 in [-8, 27], -ve/+ve (incl. 0)
!p2   nn2 = maxang + nang/4           !* nn2 in [10, 45], all +ve  (no  0)
!p2   do 29 iang=nn1,nn2              !* Bash; nn1 -> nn2 covers half circle
!p2     ii = iang                                !* ii always in range [1,nang]
!p2     if ( ii .lt.    1 ) ii = ii + nang       !* ""
!p2     if ( ii .gt. nang ) ii = ii - nang       !* ""
!p2     idif = iabs(maxang-iang) + 1             !* =10,9,..,2,1,2,..,9,10
!p2     psi2(ii) = q(mm) * cos(angl(idif))**mm   !* Normalized psi2 distr.
!  29  continue
!!p2---
!!    ==================================================================
!!
!!
!!p3
!!    Construct New Directional Distribution "psi2(:)"
!!    In an attempt to replace the original psi2(:) with
!!    a distribution based on sin()**mm with 'newmaxang' - not good enough
!!    ------------------------------------------------------------------
!!
!!    Solve for Normalizing Coefficient for Integral [1.0/(sin()**m)]
!!    Note: n1, n2 spans half circle (from 0 to +pi)
!p3   n1 =  1
!p3   n2 =  nang/2 + 1
!p3   do 36 m=1,16
!p3     sum1 = 0.
!p3     do 35 iang=n1,n2
!p3       sum1 = sum1 + sinan(iang)**m
! 35    continue
!p3     q(m) = 1./(sum1*ainc)
! 36  continue
!!p3---
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!p3
      ef2maxrow(:) = ef2(npk,:)
      maxang       = MAXLOC(ef2maxrow,1)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Shift the row so that the max is at location of 90 degress 
!!          Negative shift is to the right, Postive to the left
!!          halfangl - lower angular limit of the half maximum
!!          halfangu - upper angular limit of the half maximum
      ef2shift(:) = CSHIFT( ef2maxrow(:), (maxang-1-nang/4) )
      halfangu    = nang/4+2
      halfangl    = nang/4
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      halfmax = 0.5 * ef2(npk,maxang)
      do while((ef2shift(halfangu).gt.halfmax).and.(halfangu.lt.nang/2))
         halfangu = halfangu + 1
      enddo
      do while((ef2shift(halfangl).gt.halfmax).and.(halfangl.gt.1))
         halfangl = halfangl - 1
      enddo
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Convert angles indices with respect to peak
!!       e.g. halfangl should go to halfangl - (nang/4+1)
!!            halfangu should go to halfangu - (nang/4+1)
      halfangl = halfangl - (nang/4+1)
      halfangu = halfangu - (nang/4+1)
!!
!!    Now average the positions, round to nearest integer.
!!    -ve result means the centre is one greater than it should be.
      maxangshift = NINT( 0.5 * (halfangl + halfangu) )
      newmaxang   = maxang + maxangshift
      if (newmaxang .lt.    1) newmaxang = newmaxang + nang
      if (newmaxang .gt. nang) newmaxang = newmaxang - nang
!!p3---
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!p3
!!    Bash; need this section if you want to try sin()**mm with 'newmaxang'
!p3   y = ef2(npk,newmaxang) / ef1(npk)          !* New Energy Spread
!!
!!    Compare value of peak with q-array for closest fit to sin()**m at peak
!!    This !p3 section is needed for use with  sin()**mm
!!    Bash; Note; This new "psi2(:)" although better than original "psi2(:)"
!!          it was still not good enough:  the 2D Energy was OK but
!!          the 2D Snl, now with better symmetry, didn't always have the side lobes.
!p3   mm   = 1
!p3   qmin = abs(q(1)-y)
!p3   do 38 m=2,16
!p3     adif = abs(q(m)-y)
!p3     if ( adif.lt.qmin ) then
!p3       qmin = adif
!p3       mm   = m
!p3     endif
! 38  continue
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Final step, use 'mm' for sin()**mm
!p3   psi2(n1:n2) = (sinan(n1:n2))**mm           !* Un-norm. psi2 distr.
!p3   psi2(n1:n2) = q(mm) * psi2(n1:n2)          !* Norm. psi2 distr.
!!    Rotate peak to correct angle
!p3   psi2(:) = CSHIFT( psi2(:), newmaxang-1+nang/4 )
!!p3---
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!p4
!!    !!p4 is an override of !!p3 with mm=4, so go straight to Final step
!!    Note; all you need from !!p3 is the  "newmaxang"
!!    So it's a sin()**4 distr. shifted to "newmaxang" - worked very well
!!    ------------------------------------------------------------------
!!
!!    Solve for Normalizing Coefficient for Integral [1.0/(sin()**4)]
!!    Note: n1, n2 spans half circle (from 0 to +pi)
      n1 =  1
      n2 =  nang/2 + 1
!p4   sum1 = 0.
!p4   do 39 iang=n1,n2
!p4     sum1 = sum1 + sinan(iang)**4
! 39  continue
!p4   q4 = 1.0/(sum1*ainc)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    Change the angles that aren't zero  (0 deg to +180 deg)
      psi2(n1:n2) = (sinan(n1:n2))**4                 !* Un-norm. psi2 distr.
      q4          = 1.0/(SUM(psi2(n1:n2))*ainc)
      psi2(n1:n2) = q4 * psi2(n1:n2)                  !* Norm. psi2 distr.
!!
!!    Rotate peak to correct angle
      psi2(:) = CSHIFT( psi2(:), newmaxang-1+nang/4 )
!!p4---
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
!!    Estimate parametric spectrum and deviation from parametric spectrum
!!    ------------------------------------------------------------------
      igam  = (gam-0.4)*10 + 0.5
      sigz  = 0.109
      gam   = igam/10. + 0.4
      fdenp = gam * beta / wka(npk)**2.5
!!
!!
      do 40 irng=nrmn,nrmx
        fr = frqa(irng) / fpk
        if ( fr.le.1.0001 ) then
          if ( fr.ge.0.85 ) then
            ratio = 1.-(1.-fr)*0.7/0.15
          else
            ratio = 0.3*exp(-17.3*(0.85-fr))
          endif
          fkscl1(irng) = fdenp*ratio
          bscl1(irng)  = fkscl1(irng)/oma(irng)
        else
          z = 0.5*((fr-1.)/sigz)**1.2
          if ( z.gt.6. ) z = 6.
          ratio = 1.+exp(-z)*(gam-1.)
          fkscl1(irng) = beta*ratio/wka(irng)**2.5
          bscl1(irng)  = fkscl1(irng)/oma(irng)
        endif
!!
        do 41 iang=1,nang
          ddd = bscl1(irng) * psi2(iang) / wka(irng)  !* large-scale
          dens1(irng,iang) = ddd                      !* large-scale
          dens2(irng,iang) = act2d(irng,iang) - ddd   !* small-scale
  41    continue
  40  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END SUBROUTINE optsa2
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE snlr_fbi ( pha, ialt )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: fbi & diag2 all dim=(nrng,nang)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    For a given Action Density array dens1(k,theta), computes the    !
!!    rate-of-change array sumint(k,theta) = dN(k,theta)/dt owing to   !
!!    wave-wave interaction, as well as some ancillary arrays          !
!!    relating to positive and negative fluxes and their integrals.    !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Compute:                                                         !
!!    --------                                                         !
!!    for both -tsa and -fbi                                           !
!!    + sumint    contains scale 1 contibution for Snl -tsa & Snl -fbi !
!!                                                                     !
!!    for -tsa                                                         !
!!    + sumintsa  contains tsa approximation to Snl -tsa               !
!!                                                                     !
!!    for -fbi                                                         !
!!    + sumintp   contains scale 2 contribution to Snl -fbi            !
!!    + sumintx   contains cross interactions between scales 1 and 2   !
!!    -----------------------------------------------------------------#
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles
!!    npts      int.  Public    I   # of points on the locus
!!    nzz       int.  Public    I   linear irng x krng = (NK*(NK+1))/2
!!    ialt      int.  Public    I   integer switch ialt=2; do alternate
!!                                                 ialt=1; do not alternate
!!    kzone     int.  Public    I   zone of influence = INT(alog(4.0)/alog(dfrq))
!!    na2p1     int.  Public    I   = nang/2 + 1
!!    np2p1     int.  Public    I   = npts/2 + 1
!!    dfrq      real  Public    I   frequency multiplier for log freq. spacing
!!    frqa      R.A.  Public    I   radian frequencies (Hz);   dim=(nrng)
!!    pha       R.A.  local     I   pha = k*dk*dtheta      ;   dim=(nrng)
!!    ------------------------------------------------------------------
!!
!!    *** The 11 grid integration geometry arrays at one given depth
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2     I.A.  Public    I   Index of reference wavenumber for k2
!!    kref4     I.A.  Public    I   Idem for k4
!!    jref2     I.A.  Public    I   Index of reference angle      for k2
!!    jref4     I.A.  Public    I   Idem for k4
!!    wtk2      R.A.  Public    I   k2 Interpolation weigth along wavenumbers
!!    wtk4      R.A.  Public    I   Idem for k4
!!    wta2      R.A.  Public    I   k2 Interpolation weigth along angles
!!    wta4      R.A.  Public    I   Idem for k4
!!    tfac2     R.A.  Public    I   Norm. for interp Action Density at k2
!!    tfac4     R.A.  Public    I   Idem for k4
!!    grad      R.A.  Public    I   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!!    *** large & small scale Action Density from optsa dim=(nrng,nang)
!!    dens1     R.A.  Public    I   lrg-scl Action Density (k,theta);
!!    dens2     R.A.  Public    I   Sml-scl Action Density (k,theta);
!!    ------------------------------------------------------------------
!!
!!    for both -tsa and -fbi
!!    sumint    R.A.  local     O   contains scale 1 contribution to Snl
!!                                                       dim=(nrng,nang)
!!    for -tsa
!!    sumintsa  R.A.  local     O   contains tsa approximation to Snl -tsa
!!                                                       dim=(nrng,nang)
!!    for -fbi
!!    sumintp   R.A.  local     O   contains scale 2 contribution to Snl -fbi
!!                                                       dim=(nrng,nang)
!!    sumintx   R.A.  local     O   contains cross interactions " "   "  -fbi
!!                                                       dim=(nrng,nang)
!!    ------------------------------------------------------------------
!!
!!    for -tsa; The 2 returned arrays tsa & diag dim=(nrng,nang)
!!    tsa       R.A.  Public    O   Snl-tsa = sumint + sumintsa
!!    diag      R.A.  Public    O   Snl-tsa diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!
!!    for -fbi; The 2 returned arrays fbi & diag2 dim=(nrng,nang)
!!    fbi       R.A.  Public    O   Snl-fbi = sumint + sumintp  + sumintx
!!    diag2     R.A.  Public    O   Snl-fbi diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!    ----------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ----------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: pha(nrng)
      integer, intent(in)  :: ialt
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!    
!!    Local Parameters & variables
!!    ----------------------------- 
!!    for both -tsa and -fbi
      integer              :: irng,krng, iang,kang
      integer              :: ipt, iizz, izz
      integer              :: kmax
      integer              :: ia2, ia2p, k2, k2p
      integer              :: ia4, ia4p, k4, k4p
      integer              :: nref
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -tsa
!tsa  integer              :: nklimit, nalimit
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
      real                 :: d1,   d3,   d2,   d4
      real                 :: dp1,  dp3
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
!!    but for -tsa they are being calc. inside  if/endif if test is successful
!!    and for -fbi they are being calc. outside if/endif always
      real                 :: dz4,  dz5
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
      real                 :: dx13, ds13, dxp13, dsp13
      real                 :: dgm,  t31,  tr31
      real                 :: w2,   w2p,  wa2,  wa2p,  d2a,  d2b,  tt2
      real                 :: w4,   w4p,  wa4,  wa4p,  d4a,  d4b,  tt4
      real                 :: sumint(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -tsa
!tsa  real                 :: dz2a, dz3a,  ttsa,   trtsa
!tsa  real                 :: ddn1, ddn3,  diagk1, diagk3
!tsa  real                 :: sumintsa(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -fbi
      real                 :: dp2,   dp4
      real                 :: d2pa,  d4pa
      real                 :: d2pb,  d4pb
      real                 :: dz1,   dz2,   dz3,   dz6,   dz7,   dz8
      real                 :: dgmp,  tp31,  trp31, dzsum, txp31, trx31
!!
!!    for -fbi; Bash added 4 new terms for a full expression of diag2 term
      real                 :: ddp1,  ddp2,  ddp3,  ddp4   !* ddpi=di+dpi for i=1,4
      real                 :: dd2n1, dd2n3, diag2k1, diag2k3
      real                 :: sumintp(nrng,nang)
      real                 :: sumintx(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    for -tsa
!!    Bash; hardwire these two parameters
!tsa  nklimit = 6
!tsa  nalimit = 6
!!
!!
!!ini
!!    Bash; move initialization of all returned arrays from below to here
!!    ------------------------------------------------------------------
!!    for both -tsa and -fbi
!!    sumint is now initialized here instead of below!
      sumint(:,:)   = 0.0
!!
!!    for -tsa
!!    sumintsa are now initialized here instead of below!
!tsa  sumintsa(:,:) = 0.0
!tsa  tsa(:,:)      = 0.0
!tsa  diag(:,:)     = 0.0
!!
!!    for -fbi
!!    sumintp and sumintx are now initialized here instead of below!
      sumintp(:,:)  = 0.0
      sumintx(:,:)  = 0.0
      fbi(:,:)      = 0.0
      diag2(:,:)    = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!    ##################################################################
!!
!!
!!    for -tsa
!tsa  ddn1    = 0.0                           !* for -tsa diag  [dN/dn1]
!tsa  ddn3    = 0.0                           !* for -tsa diag  [dN/dn3]
!!
!!    for -fbi
      dd2n1   = 0.0                           !* for -fbi diag2 [dN/dn1]
      dd2n3   = 0.0                           !* for -fbi diag2 [dN/dn3]
!!
!!
!!
!!50
      do 50 irng=1,nrng,ialt
!!kz
        kmax = min(irng+kzone, nrng)   !* Bash; Sometimes a locus pt is outside nrng
!kz     kmax = min(irng+kzone, nrng-1) !* Bash; Taking 1 out will not affect kzone, try it
!!kz---
!!
        iizz = (nrng-1)*(irng-1) - ((irng-2)*(irng-1))/2
!!      ----------------------------------------------------------------
!!
!!
!!60
        do 60 iang=1,nang,ialt
!!
!!        for both -tsa and -fbi
          d1   = dens1(irng,iang)
          dp1  = dens2(irng,iang)
!!
!!        for -fbi
          ddp1 = d1+dp1             !! for full expression of diag2 term
!!
!!70
!!kz
!kz       do 70 krng=irng,nrng
          do 70 krng=irng,kmax,ialt
!!
!!          for both -tsa and -fbi
!!          Bash; check5  be consistent with gridsetr
!!          moved here from below (was after do 80 kang=1,nang)
!!          and changed go to 80 into go to 70 (i.e. go to next krng)
!kz         if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 70  !* original gridsetr
!kz         if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 70  !* original snlr_'s
!kz         if ( frqa(krng)/frqa(irng) .gt. 2. ) go to 70  !* Bash; use .gt. 2
!!kz---
!!
            izz = krng + iizz
!!          ------------------------------------------------------------
!!
!!80
            do 80 kang=1,nang,ialt
!!
!!            for both -tsa and -fbi
!!ba1         Bash; Remove self interaction
!!                  skip k1 but keep the opposite angle to k1 - original setting
              if ( krng.eq.irng ) then              !* wn3 = wn1
                if ( kang.eq.iang ) go to 80        !* th3 = th1
              endif
!!ba1---
!!            ----------------------------------------------------------
!!
!!            for both -tsa and -fbi
              d3   = dens1(krng,kang)
              dp3  = dens2(krng,kang)
!!
!!            for -fbi
              ddp3 = d3+dp3         !! for full expression of diag2 term
!!
!!
!!            for both -tsa and -fbi
              nref = kang - iang + 1
              if ( nref .lt. 1 ) nref = nref + nang
!!
!!
!!            for both -tsa and -fbi
!!            Bash; check5  be consistent with gridsetr
!!                  and move this test above right after do 70 krng=irng,nrng
!x            if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 80  !* gridsetr
!b            if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 80  !* original
!!
!!
!!            for both -tsa and -fbi
              t31     = 0.0             !* must be reset to 0.0
!!
!!            for -tsa
!tsa          ttsa    = 0.0             !* must be reset to 0.0
!tsa          diagk1  = 0.0             !* must be reset to 0.0
!tsa          diagk3  = 0.0             !* must be reset to 0.0
!!
!!            for -fbi
              tp31    = 0.0             !* must be reset to 0.0
              txp31   = 0.0             !* must be reset to 0.0
              diag2k1 = 0.0             !* must be reset to 0.0
              diag2k3 = 0.0             !* must be reset to 0.0
!!
!!            for both -tsa and -fbi
              dx13  = d1*d3
              ds13  = d3-d1
              dxp13 = dp1*dp3
              dsp13 = dp3-dp1
!!            ----------------------------------------------------------
!!
!!90
              do 90 ipt=1,npts
!!
!!              for both -tsa and -fbi
!!              save time by skipping insignificant contributions
!!e-30
!e-30           if ( grad(ipt,nref,izz) .lt. 1.e-30 ) go to 90
!!e-30---
                if ( grad(ipt,nref,izz) .lt. 1.e-15 ) go to 90
!!e-30---
!!              --------------------------------------------------------
!!
!!xlc1          Bash; skip k1 but keep the opposite angle to k1 - original setting
!xlc1           if ( kang.eq.iang ) then                     !* th3=+th1
!xlc1             if (ipt.eq.1 .or. ipt.eq.np2p1) go to 90   !* skip x-axis loci
!xlc1           end if
!!xlc1---
!!              --------------------------------------------------------
!!
!!
!!2             Estimation of Density for wave #2
!!
!!              for both -tsa and -fbi
                k2  = kref2(ipt,nref,izz)
                k2p = k2 + 1
                w2  = wtk2(ipt,nref,izz)
                w2p = 1. - w2
!!
!!              for both -tsa and -fbi
                ia2 = iang + jref2(ipt,nref,izz)
                if ( ia2 .gt. nang )  ia2  = ia2  - nang
!!
!!              for both -tsa and -fbi
                ia2p = ia2 + 1
                if ( ia2p .gt. nang ) ia2p = ia2p - nang
!!
!!              for both -tsa and -fbi
                wa2  = wta2(ipt,nref,izz)
                wa2p = 1. - wa2
                d2a  = w2 * dens1(k2,ia2)  + w2p * dens1(k2p,ia2)
                d2b  = w2 * dens1(k2,ia2p) + w2p * dens1(k2p,ia2p)
                tt2  = tfac2(ipt,nref,izz)
                d2   = (wa2*d2a  + wa2p*d2b)  * tt2
!!
!!              for -fbi
                d2pa = w2 * dens2(k2,ia2)  + w2p * dens2(k2p,ia2)
                d2pb = w2 * dens2(k2,ia2p) + w2p * dens2(k2p,ia2p)
!!
!!              for -fbi
                dp2  = (wa2*d2pa + wa2p*d2pb) * tt2   !* for -fbi
                ddp2 = d2+dp2       !! for full expression of diag2 term
!!              ========================================================
!!
!!
!!4             Estimation of Density for wave #4
!!
!!              for both -tsa and -fbi
                k4  = kref4(ipt,nref,izz)
                k4p = k4 + 1
                w4  = wtk4(ipt,nref,izz)
                w4p = 1. - w4
!!
!!              for both -tsa and -fbi
                ia4 = iang + jref4(ipt,nref,izz)
                if ( ia4 .gt. nang )  ia4  = ia4  - nang
!!
!!              for both -tsa and -fbi
                ia4p= ia4 + 1
                if ( ia4p .gt. nang ) ia4p = ia4p - nang
!!
!!              for both -tsa and -fbi
                wa4  = wta4(ipt,nref,izz)
                wa4p = 1. - wa4
                d4a  = w4*dens1(k4,ia4)  + w4p*dens1(k4p,ia4)
                d4b  = w4*dens1(k4,ia4p) + w4p*dens1(k4p,ia4p)
                tt4  = tfac4(ipt,nref,izz)
                d4   = (wa4*d4a  + wa4p*d4b)  * tt4
!!
!!              for -fbi
                d4pa = w4*dens2(k4,ia4)  + w4p*dens2(k4p,ia4)
                d4pb = w4*dens2(k4,ia4p) + w4p*dens2(k4p,ia4p)
!!
!!              for -fbi
                dp4  = (wa4*d4pa + wa4p*d4pb) * tt4   !* for -fbi
                ddp4 = d4+dp4       !! for full expression of diag2 term
!!              ========================================================
!!
!!
!!              for both -tsa and -fbi
                dgm  = dx13*(d4-d2) + ds13*d4*d2        !* dgm=B of R&P'08 eqn(8)
!!                                         !* represents Broad Scale interactions
                t31  = t31  + dgm  * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!
!!              for -fbi
                dgmp = dxp13*(dp4-dp2) + dsp13*dp4*dp2  !* dgmp=L of R&P'08 eqn(8)
!!                                          !* represents Local Scale interactions
                tp31 = tp31 + dgmp * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
!!              for -tsa : -diag
!!              use this expression for the diagonal term
!!              whose derivation neglect "dp2" & "dp4"
!tsa            ddn1   = (d3+dp3)*(d4-d2) - d4*d2              !* dN/dn1
!tsa            ddn3   = (d1+dp1)*(d4-d2) + d4*d2              !* dN/dn3
!tsa            diagk1 = diagk1 + ddn1 * grad(ipt,nref,izz)
!tsa            diagk3 = diagk3 + ddn3 * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!
!!              for -fbi : -diag2
!!              use the full expression for the diagonal terms
!!              whose derivation keeps all large + small scale
                dd2n1   = ddp3*(ddp4-ddp2) - ddp4*ddp2         !* dN/dn1
                dd2n3   = ddp1*(ddp4-ddp2) + ddp4*ddp2         !* dN/dn3
                diag2k1 = diag2k1 + dd2n1 * grad(ipt,nref,izz)
                diag2k3 = diag2k3 + dd2n3 * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
!!              for -fbi
                dz1  = dx13    * (dp4-dp2)
                dz2  = d1*dp3  * ((d4-d2)+(dp4-dp2))
                dz3  = d3*dp1  * ((d4-d2)+(dp4-dp2))
!!
!!              for -fbi (calc. dz4 & dz5 here)
                dz4  = dxp13   * (d4-d2)
                dz5  = d2*d4   *  dsp13
!!
!!              for -tsa
!!              Cross-interactions between parametric and perturbation
!!              that occur only when k3 is close enough to k1
!!              Bash; added an extra check on (nang-nalimit)
!b              if ( iabs(irng-krng).lt.nklimit .and.                 &
!b                   iabs(iang-kang).lt.nalimit )    then        !* original
!!
!tsa            if (     (krng-irng).lt.nklimit .and.                 &
!tsa               ( iabs(kang-iang).lt.nalimit .or.                  &
!tsa                 iabs(kang-iang).gt.(nang-nalimit) ) )  then !* Bash
!!
!!                for -tsa (calc. dz4 & dz5 here)
!tsa              dz4  = dxp13   * (d4-d2)
!tsa              dz5  = d2*d4   *  dsp13
!tsa              dz2a = d1*dp3 * (d4-d2)
!tsa              dz3a = d3*dp1 * (d4-d2)
!!
!tsa              ttsa = ttsa + (dz4+dz5+dz2a+dz3a)*grad(ipt,nref,izz)
!!
!tsa            endif
!!              --------------------------------------------------------
!!
!!              for -fbi
                dz6   = d2*dp4  * (ds13+dsp13)
                dz7   = d4*dp2  * (ds13+dsp13)
                dz8   = dp2*dp4 *  ds13
                dzsum = dz1 + dz2 + dz3 + dz4 + dz5 + dz6 + dz7 + dz8
                txp31 = txp31 + dzsum * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
  90          continue                        !* end of ipt (locus) loop
!!            ----------------------------------------------------------
!!
!!
!!            multiply the following components by factor 2. in here
!!
!!            for both -tsa and -fbi
              tr31  = 2. * t31
!!
!!            for -tsa
!tsa          trtsa = 2. * ttsa
!!
!!            for -fbi
              trp31 = 2. * tp31
              trx31 = 2. * txp31
!!
!!            for -tsa : -diag
!tsa          diagk1  = 2. * diagk1
!tsa          diagk3  = 2. * diagk3
!!
!!            for -fbi : -diag2
              diag2k1 = 2. * diag2k1
              diag2k3 = 2. * diag2k3
!!            ----------------------------------------------------------
!!
!!            for both -tsa and -fbi
              sumint(irng,iang)  = sumint(irng,iang)  + tr31*pha(krng)
              sumint(krng,kang)  = sumint(krng,kang)  - tr31*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -tsa
!tsa          sumintsa(irng,iang)= sumintsa(irng,iang)+ trtsa*pha(krng)
!tsa          sumintsa(krng,kang)= sumintsa(krng,kang)- trtsa*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -fbi
              sumintp(irng,iang) = sumintp(irng,iang) + trp31*pha(krng)
              sumintp(krng,kang) = sumintp(krng,kang) - trp31*pha(irng)
!!
!!            for -fbi
              sumintx(irng,iang) = sumintx(irng,iang) + trx31*pha(krng)
              sumintx(krng,kang) = sumintx(krng,kang) - trx31*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -tsa : -diag
!tsa          diag(irng,iang) = diag(irng,iang)  + diagk1*pha(krng)
!tsa          diag(krng,kang) = diag(krng,kang)  - diagk3*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -fbi : -diag2
              diag2(irng,iang) = diag2(irng,iang) + diag2k1*pha(krng)
              diag2(krng,kang) = diag2(krng,kang) - diag2k3*pha(irng)
!!            ----------------------------------------------------------
!!
  80        continue                              !* end of kang loop
!!
  70      continue                                !* end of krng loop
!!
  60    continue                                  !* end of iang loop
!!
  50  continue                                    !* end of irng loop
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    Final sum-up to get Snl and diag. term to be returned
!!
!!    for -tsa
!tsa  tsa(:,:)   = sumint(:,:) + sumintsa(:,:)
!b    diag(:,:)  = diag(:,:)   !* is Ok, already summed up
!!
!!    for -fbi
      fbi(:,:)   = sumint(:,:) + sumintp(:,:) + sumintx(:,:)
!b    diag2(:,:) = diag2(:,:)  !* is Ok, already summed up
!!    --------------------------------------------------------------------------
!!    ==========================================================================
!!
!!
!!alt Call interp2 only if ialt=2,
!!    Interpolate bi-linearly to fill in tsa/fbi & diag/diag2 arrays
!!    after alternating the irng, iang, krng & kang loops above
!!    ------------------------------------------------------------------
      if ( ialt.eq.2 ) then
!!      for -tsa
!tsa    call interp2 ( tsa )
!tsa    call interp2 ( diag )
!!
!!      for -fbi
        call interp2 ( fbi )
        call interp2 ( diag2 )
      endif
!!alt---
!!    --------------------------------------------------------------------------
!!    ==========================================================================
!!
      RETURN
!!
      END SUBROUTINE snlr_fbi
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE snlr_tsa ( pha, ialt )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    ------------------------------------------------------------------
!!
!!    it returns: tsa & diag  all dim=(nrng,nang)
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    For a given Action Density array dens1(k,theta), computes the    !
!!    rate-of-change array sumint(k,theta) = dN(k,theta)/dt owing to   !
!!    wave-wave interaction, as well as some ancillary arrays          !
!!    relating to positive and negative fluxes and their integrals.    !
!!                                                                     !
!!    -----------------------------------------------------------------#
!!                                                                     !
!!    Compute:                                                         !
!!    --------                                                         !
!!    for both -tsa and -fbi                                           !
!!    + sumint    contains scale 1 contibution for Snl -tsa & Snl -fbi !
!!                                                                     !
!!    for -tsa                                                         !
!!    + sumintsa  contains tsa approximation to Snl -tsa               !
!!                                                                     !
!!    for -fbi                                                         !
!!    + sumintp   contains scale 2 contribution to Snl -fbi            !
!!    + sumintx   contains cross interactions between scales 1 and 2   !
!!    -----------------------------------------------------------------#
!!
!! 2. Method :
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles
!!    npts      int.  Public    I   # of points on the locus
!!    nzz       int.  Public    I   linear irng x krng = (NK*(NK+1))/2
!!    ialt      int.  Public    I   integer switch ialt=2; do alternate
!!                                                 ialt=1; do not alternate
!!    kzone     int.  Public    I   zone of influence = INT(alog(4.0)/alog(dfrq))
!!    na2p1     int.  Public    I   = nang/2 + 1
!!    np2p1     int.  Public    I   = npts/2 + 1
!!    dfrq      real  Public    I   frequency multiplier for log freq. spacing
!!    frqa      R.A.  Public    I   radian frequencies (Hz);   dim=(nrng)
!!    pha       R.A.  local     I   pha = k*dk*dtheta      ;   dim=(nrng)
!!    ------------------------------------------------------------------
!!
!!    *** The 11 grid integration geometry arrays at one given depth
!!    *** from gridsetr.            dim=(npts,nang,nzz,ndep)
!!    kref2     I.A.  Public    I   Index of reference wavenumber for k2
!!    kref4     I.A.  Public    I   Idem for k4
!!    jref2     I.A.  Public    I   Index of reference angle      for k2
!!    jref4     I.A.  Public    I   Idem for k4
!!    wtk2      R.A.  Public    I   k2 Interpolation weigth along wavenumbers
!!    wtk4      R.A.  Public    I   Idem for k4
!!    wta2      R.A.  Public    I   k2 Interpolation weigth along angles
!!    wta4      R.A.  Public    I   Idem for k4
!!    tfac2     R.A.  Public    I   Norm. for interp Action Density at k2
!!    tfac4     R.A.  Public    I   Idem for k4
!!    grad      R.A.  Public    I   Coupling and gradient term in integral
!!                                  grad = C * H * g**2 * ds / |dW/dn|
!!    ------------------------------------------------------------------
!!
!!    *** large & small scale Action Density from optsa dim=(nrng,nang)
!!    dens1     R.A.  Public    I   lrg-scl Action Density (k,theta);
!!    dens2     R.A.  Public    I   Sml-scl Action Density (k,theta);
!!    ------------------------------------------------------------------
!!
!!    for both -tsa and -fbi
!!    sumint    R.A.  local     O   contains scale 1 contribution to Snl
!!                                                       dim=(nrng,nang)
!!    for -tsa
!!    sumintsa  R.A.  local     O   contains tsa approximation to Snl -tsa
!!                                                       dim=(nrng,nang)
!!    for -fbi
!!    sumintp   R.A.  local     O   contains scale 2 contribution to Snl -fbi
!!                                                       dim=(nrng,nang)
!!    sumintx   R.A.  local     O   contains cross interactions " "   "  -fbi
!!                                                       dim=(nrng,nang)
!!    ------------------------------------------------------------------
!!
!!    for -tsa; The 2 returned arrays tsa & diag dim=(nrng,nang)
!!    tsa       R.A.  Public    O   Snl-tsa = sumint + sumintsa
!!    diag      R.A.  Public    O   Snl-tsa diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!
!!    for -fbi; The 2 returned arrays fbi & diag2 dim=(nrng,nang)
!!    fbi       R.A.  Public    O   Snl-fbi = sumint + sumintp  + sumintx
!!    diag2     R.A.  Public    O   Snl-fbi diagonal term = [dN/dn1]
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!    ----------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!     gridsetr  Subr. W3SNL4MD Setup geometric integration grid
!!    ----------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: pha(nrng)
      integer, intent(in)  :: ialt
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!    
!!    Local Parameters & variables
!!    ----------------------------- 
!!    for both -tsa and -fbi
      integer              :: irng,krng, iang,kang
      integer              :: ipt, iizz, izz
      integer              :: kmax
      integer              :: ia2, ia2p, k2, k2p
      integer              :: ia4, ia4p, k4, k4p
      integer              :: nref
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -tsa
      integer              :: nklimit, nalimit
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
      real                 :: d1,   d3,   d2,   d4
      real                 :: dp1,  dp3
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
!!    but for -tsa they are being calc. inside  if/endif if test is successful
!!    and for -fbi they are being calc. outside if/endif always
      real                 :: dz4,  dz5
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for both -tsa and -fbi
      real                 :: dx13, ds13, dxp13, dsp13
      real                 :: dgm,  t31,  tr31
      real                 :: w2,   w2p,  wa2,  wa2p,  d2a,  d2b,  tt2
      real                 :: w4,   w4p,  wa4,  wa4p,  d4a,  d4b,  tt4
      real                 :: sumint(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -tsa
      real                 :: dz2a, dz3a,  ttsa,   trtsa
      real                 :: ddn1, ddn3,  diagk1, diagk3
      real                 :: sumintsa(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
!!    for -fbi
!fbi  real                 :: dp2,   dp4
!fbi  real                 :: d2pa,  d4pa
!fbi  real                 :: d2pb,  d4pb
!fbi  real                 :: dz1,   dz2,   dz3,   dz6,   dz7,   dz8
!fbi  real                 :: dgmp,  tp31,  trp31, dzsum, txp31, trx31
!!
!!    for -fbi; Bash added 4 new terms for a full expression of diag2 term
!fbi  real                 :: ddp1,  ddp2,  ddp3,  ddp4   !* ddpi=di+dpi for i=1,4
!fbi  real                 :: dd2n1, dd2n3, diag2k1, diag2k3
!fbi  real                 :: sumintp(nrng,nang)
!fbi  real                 :: sumintx(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    for -tsa
!!    Bash; hardwire these two parameters
      nklimit = 6
      nalimit = 6
!!
!!
!!ini
!!    Bash; move initialization of all returned arrays from below to here
!!    ------------------------------------------------------------------
!!    for both -tsa and -fbi
!!    sumint is now initialized here instead of below!
      sumint(:,:)   = 0.0
!!
!!    for -tsa
!!    sumintsa are now initialized here instead of below!
      sumintsa(:,:) = 0.0
      tsa(:,:)      = 0.0
      diag(:,:)     = 0.0
!!
!!    for -fbi
!!    sumintp and sumintx are now initialized here instead of below!
!fbi  sumintp(:,:)  = 0.0
!fbi  sumintx(:,:)  = 0.0
!fbi  fbi(:,:)      = 0.0
!fbi  diag2(:,:)    = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!    ------------------------------------------------------------------
!!    ##################################################################
!!
!!
!!    for -tsa
      ddn1    = 0.0                           !* for -tsa diag  [dN/dn1]
      ddn3    = 0.0                           !* for -tsa diag  [dN/dn3]
!!
!!    for -fbi
!fbi  dd2n1   = 0.0                           !* for -fbi diag2 [dN/dn1]
!fbi  dd2n3   = 0.0                           !* for -fbi diag2 [dN/dn3]
!!
!!
!!
!!50
      do 50 irng=1,nrng,ialt
!!kz
        kmax = min(irng+kzone, nrng)   !* Bash; Sometimes a locus pt is outside nrng
!kz     kmax = min(irng+kzone, nrng-1) !* Bash; Taking 1 out will not affect kzone, try it
!!kz---
!!
        iizz = (nrng-1)*(irng-1) - ((irng-2)*(irng-1))/2
!!      ----------------------------------------------------------------
!!
!!
!!60
        do 60 iang=1,nang,ialt
!!
!!        for both -tsa and -fbi
          d1   = dens1(irng,iang)
          dp1  = dens2(irng,iang)
!!
!!        for -fbi
!fbi      ddp1 = d1+dp1             !! for full expression of diag2 term
!!
!!70
!!kz
!kz       do 70 krng=irng,nrng
          do 70 krng=irng,kmax,ialt
!!
!!          for both -tsa and -fbi
!!          Bash; check5  be consistent with gridsetr
!!          moved here from below (was after do 80 kang=1,nang)
!!          and changed go to 80 into go to 70 (i.e. go to next krng)
!kz         if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 70  !* original gridsetr
!kz         if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 70  !* original snlr_'s
!kz         if ( frqa(krng)/frqa(irng) .gt. 2. ) go to 70  !* Bash; use .gt. 2
!!kz---
!!
            izz = krng + iizz
!!          ------------------------------------------------------------
!!
!!80
            do 80 kang=1,nang,ialt
!!
!!            for both -tsa and -fbi
!!ba1         Bash; Remove self interaction
!!                  skip k1 but keep the opposite angle to k1 - original setting
              if ( krng.eq.irng ) then              !* wn3 = wn1
                if ( kang.eq.iang ) go to 80        !* th3 = th1
              endif
!!ba1---
!!            ----------------------------------------------------------
!!
!!            for both -tsa and -fbi
              d3   = dens1(krng,kang)
              dp3  = dens2(krng,kang)
!!
!!            for -fbi
!fbi          ddp3 = d3+dp3         !! for full expression of diag2 term
!!
!!
!!            for both -tsa and -fbi
              nref = kang - iang + 1
              if ( nref .lt. 1 ) nref = nref + nang
!!
!!
!!            for both -tsa and -fbi
!!            Bash; check5  be consistent with gridsetr
!!                  and move this test above right after do 70 krng=irng,nrng
!x            if ( frqa(krng)/frqa(irng) .gt. 4. ) go to 80  !* gridsetr
!b            if ( frqa(krng)/frqa(irng) .gt. 3. ) go to 80  !* original
!!
!!
!!            for both -tsa and -fbi
              t31     = 0.0             !* must be reset to 0.0
!!
!!            for -tsa
              ttsa    = 0.0             !* must be reset to 0.0
              diagk1  = 0.0             !* must be reset to 0.0
              diagk3  = 0.0             !* must be reset to 0.0
!!
!!            for -fbi
!fbi          tp31    = 0.0             !* must be reset to 0.0
!fbi          txp31   = 0.0             !* must be reset to 0.0
!fbi          diag2k1 = 0.0             !* must be reset to 0.0
!fbi          diag2k3 = 0.0             !* must be reset to 0.0
!!
!!            for both -tsa and -fbi
              dx13  = d1*d3
              ds13  = d3-d1
              dxp13 = dp1*dp3
              dsp13 = dp3-dp1
!!            ----------------------------------------------------------
!!
!!90
              do 90 ipt=1,npts
!!
!!              for both -tsa and -fbi
!!              save time by skipping insignificant contributions
!!e-30
!e-30           if ( grad(ipt,nref,izz) .lt. 1.e-30 ) go to 90
!!e-30---
                if ( grad(ipt,nref,izz) .lt. 1.e-15 ) go to 90
!!e-30---
!!              --------------------------------------------------------
!!
!!xlc1          Bash; skip k1 but keep the opposite angle to k1 - original setting
!xlc1           if ( kang.eq.iang ) then                     !* th3=+th1
!xlc1             if (ipt.eq.1 .or. ipt.eq.np2p1) go to 90   !* skip x-axis loci
!xlc1           end if
!!xlc1---
!!              --------------------------------------------------------
!!
!!
!!2             Estimation of Density for wave #2
!!
!!              for both -tsa and -fbi
                k2  = kref2(ipt,nref,izz)
                k2p = k2 + 1
                w2  = wtk2(ipt,nref,izz)
                w2p = 1. - w2
!!
!!              for both -tsa and -fbi
                ia2 = iang + jref2(ipt,nref,izz)
                if ( ia2 .gt. nang )  ia2  = ia2  - nang
!!
!!              for both -tsa and -fbi
                ia2p = ia2 + 1
                if ( ia2p .gt. nang ) ia2p = ia2p - nang
!!
!!              for both -tsa and -fbi
                wa2  = wta2(ipt,nref,izz)
                wa2p = 1. - wa2
                d2a  = w2 * dens1(k2,ia2)  + w2p * dens1(k2p,ia2)
                d2b  = w2 * dens1(k2,ia2p) + w2p * dens1(k2p,ia2p)
                tt2  = tfac2(ipt,nref,izz)
                d2   = (wa2*d2a  + wa2p*d2b)  * tt2
!!
!!              for -fbi
!fbi            d2pa = w2 * dens2(k2,ia2)  + w2p * dens2(k2p,ia2)
!fbi            d2pb = w2 * dens2(k2,ia2p) + w2p * dens2(k2p,ia2p)
!!
!!              for -fbi
!fbi            dp2  = (wa2*d2pa + wa2p*d2pb) * tt2   !* for -fbi
!fbi            ddp2 = d2+dp2       !! for full expression of diag2 term
!!              ========================================================
!!
!!
!!4             Estimation of Density for wave #4
!!
!!              for both -tsa and -fbi
                k4  = kref4(ipt,nref,izz)
                k4p = k4 + 1
                w4  = wtk4(ipt,nref,izz)
                w4p = 1. - w4
!!
!!              for both -tsa and -fbi
                ia4 = iang + jref4(ipt,nref,izz)
                if ( ia4 .gt. nang )  ia4  = ia4  - nang
!!
!!              for both -tsa and -fbi
                ia4p= ia4 + 1
                if ( ia4p .gt. nang ) ia4p = ia4p - nang
!!
!!              for both -tsa and -fbi
                wa4  = wta4(ipt,nref,izz)
                wa4p = 1. - wa4
                d4a  = w4*dens1(k4,ia4)  + w4p*dens1(k4p,ia4)
                d4b  = w4*dens1(k4,ia4p) + w4p*dens1(k4p,ia4p)
                tt4  = tfac4(ipt,nref,izz)
                d4   = (wa4*d4a  + wa4p*d4b)  * tt4
!!
!!              for -fbi
!fbi            d4pa = w4*dens2(k4,ia4)  + w4p*dens2(k4p,ia4)
!fbi            d4pb = w4*dens2(k4,ia4p) + w4p*dens2(k4p,ia4p)
!!
!!              for -fbi
!fbi            dp4  = (wa4*d4pa + wa4p*d4pb) * tt4   !* for -fbi
!fbi            ddp4 = d4+dp4       !! for full expression of diag2 term
!!              ========================================================
!!
!!
!!              for both -tsa and -fbi
                dgm  = dx13*(d4-d2) + ds13*d4*d2        !* dgm=B of R&P'08 eqn(8)
!!                                         !* represents Broad Scale interactions
                t31  = t31  + dgm  * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!
!!              for -fbi
!fbi            dgmp = dxp13*(dp4-dp2) + dsp13*dp4*dp2  !* dgmp=L of R&P'08 eqn(8)
!!                                          !* represents Local Scale interactions
!fbi            tp31 = tp31 + dgmp * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
!!              for -tsa : -diag
!!              use this expression for the diagonal term
!!              whose derivation neglect "dp2" & "dp4"
                ddn1   = (d3+dp3)*(d4-d2) - d4*d2              !* dN/dn1
                ddn3   = (d1+dp1)*(d4-d2) + d4*d2              !* dN/dn3
                diagk1 = diagk1 + ddn1 * grad(ipt,nref,izz)
                diagk3 = diagk3 + ddn3 * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!
!!              for -fbi : -diag2
!!              use the full expression for the diagonal terms
!!              whose derivation keeps all large + small scale
!fbi            dd2n1   = ddp3*(ddp4-ddp2) - ddp4*ddp2         !* dN/dn1
!fbi            dd2n3   = ddp1*(ddp4-ddp2) + ddp4*ddp2         !* dN/dn3
!fbi            diag2k1 = diag2k1 + dd2n1 * grad(ipt,nref,izz)
!fbi            diag2k3 = diag2k3 + dd2n3 * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
!!              for -fbi
!fbi            dz1  = dx13    * (dp4-dp2)
!fbi            dz2  = d1*dp3  * ((d4-d2)+(dp4-dp2))
!fbi            dz3  = d3*dp1  * ((d4-d2)+(dp4-dp2))
!!
!!              for -fbi (calc. dz4 & dz5 here)
!fbi            dz4  = dxp13   * (d4-d2)
!fbi            dz5  = d2*d4   *  dsp13
!!
!!              for -tsa
!!              Cross-interactions between parametric and perturbation
!!              that occur only when k3 is close enough to k1
!!              Bash; added an extra check on (nang-nalimit)
!b              if ( iabs(irng-krng).lt.nklimit .and.                 &
!b                   iabs(iang-kang).lt.nalimit )    then        !* original
!!
                if (     (krng-irng).lt.nklimit .and.                 &
                   ( iabs(kang-iang).lt.nalimit .or.                  &
                     iabs(kang-iang).gt.(nang-nalimit) ) )  then !* Bash
!!
!!                for -tsa (calc. dz4 & dz5 here)
                  dz4  = dxp13   * (d4-d2)
                  dz5  = d2*d4   *  dsp13
                  dz2a = d1*dp3 * (d4-d2)
                  dz3a = d3*dp1 * (d4-d2)
!!
                  ttsa = ttsa + (dz4+dz5+dz2a+dz3a)*grad(ipt,nref,izz)
!!
                endif
!!              --------------------------------------------------------
!!
!!              for -fbi
!fbi            dz6   = d2*dp4  * (ds13+dsp13)
!fbi            dz7   = d4*dp2  * (ds13+dsp13)
!fbi            dz8   = dp2*dp4 *  ds13
!fbi            dzsum = dz1 + dz2 + dz3 + dz4 + dz5 + dz6 + dz7 + dz8
!fbi            txp31 = txp31 + dzsum * grad(ipt,nref,izz)
!!              --------------------------------------------------------
!!              ========================================================
!!
!!
  90          continue                        !* end of ipt (locus) loop
!!            ----------------------------------------------------------
!!
!!
!!            multiply the following components by factor 2. in here
!!
!!            for both -tsa and -fbi
              tr31  = 2. * t31
!!
!!            for -tsa
              trtsa = 2. * ttsa
!!
!!            for -fbi
!fbi          trp31 = 2. * tp31
!fbi          trx31 = 2. * txp31
!!
!!            for -tsa : -diag
              diagk1  = 2. * diagk1
              diagk3  = 2. * diagk3
!!
!!            for -fbi : -diag2
!fbi          diag2k1 = 2. * diag2k1
!fbi          diag2k3 = 2. * diag2k3
!!            ----------------------------------------------------------
!!
!!            for both -tsa and -fbi
              sumint(irng,iang)  = sumint(irng,iang)  + tr31*pha(krng)
              sumint(krng,kang)  = sumint(krng,kang)  - tr31*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -tsa
              sumintsa(irng,iang)= sumintsa(irng,iang)+ trtsa*pha(krng)
              sumintsa(krng,kang)= sumintsa(krng,kang)- trtsa*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -fbi
!fbi          sumintp(irng,iang) = sumintp(irng,iang) + trp31*pha(krng)
!fbi          sumintp(krng,kang) = sumintp(krng,kang) - trp31*pha(irng)
!!
!!            for -fbi
!fbi          sumintx(irng,iang) = sumintx(irng,iang) + trx31*pha(krng)
!fbi          sumintx(krng,kang) = sumintx(krng,kang) - trx31*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -tsa : -diag
              diag(irng,iang) = diag(irng,iang)  + diagk1*pha(krng)
              diag(krng,kang) = diag(krng,kang)  - diagk3*pha(irng)
!!            ----------------------------------------------------------
!!
!!            for -fbi : -diag2
!fbi          diag2(irng,iang) = diag2(irng,iang) + diag2k1*pha(krng)
!fbi          diag2(krng,kang) = diag2(krng,kang) - diag2k3*pha(irng)
!!            ----------------------------------------------------------
!!
  80        continue                              !* end of kang loop
!!
  70      continue                                !* end of krng loop
!!
  60    continue                                  !* end of iang loop
!!
  50  continue                                    !* end of irng loop
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
!!    Final sum-up to get Snl and diag. term to be returned
!!
!!    for -tsa
      tsa(:,:)   = sumint(:,:) + sumintsa(:,:)
!b    diag(:,:)  = diag(:,:)   !* is Ok, already summed up
!!
!!    for -fbi
!fbi  fbi(:,:)   = sumint(:,:) + sumintp(:,:) + sumintx(:,:)
!b    diag2(:,:) = diag2(:,:)  !* is Ok, already summed up
!!    --------------------------------------------------------------------------
!!    ==========================================================================
!!
!!
!!alt Call interp2 only if ialt=2,
!!    Interpolate bi-linearly to fill in tsa/fbi & diag/diag2 arrays
!!    after alternating the irng, iang, krng & kang loops above
!!    ------------------------------------------------------------------
      if ( ialt.eq.2 ) then
!!      for -tsa
        call interp2 ( tsa )
        call interp2 ( diag )
!!
!!      for -fbi
!fbi    call interp2 ( fbi )
!fbi    call interp2 ( diag2 )
      endif
!!alt---
!!    --------------------------------------------------------------------------
!!    ==========================================================================
!!
      RETURN
!!
      END SUBROUTINE snlr_tsa
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE interp2 ( X )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!
!! 1. Purpose :
!!
!!    Interpolate bi-linearly to fill in tsa/fbi & diag/diag2 arrays
!!    and then (optional) smooth the interior and the corners
!!    after alternating the irng, iang, krng & kang loops in snlr's
!!
!! 2. Method :
!! 
!! 3. Parameters :
!!    
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    nrng      int.  Public    I   # of freq. or rings
!!    nang      int.  Public    I   # of angles 
!!    ismo      int.  Local     I   switch;  ismo=0 skip smoothing
!!                                           ismo.ne.0 do smoothing
!!    X         R.A.  Local    I/O  Array to be ineterp. & smoothing
!!                                  its returned to snlr_tsa as tsa or diag
!!                                           and to snlr_fbi as fbi or diag2
!!                                           dim=(nrng,nang)
!!    ------------------------------------------------------------------
!!
!! 4. Subroutines used :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!    ----------------------------------------------------------------
!!
!! 5. Called by :
!!
!!     Name      Type  Module   Description
!!    ----------------------------------------------------------------
!!     snlr_tsa  Subr. W3SNL4MD Computes dN(k,theta)/dt for TSA
!!     --------                 due to wave-wave inter. (set itsa = 1)
!!     snlr_fbi  Subr. W3SNL4MD Computes dN(k,theta)/dt for FBI
!!     --------                 due to wave-wave inter. (set itsa = 0)
!!    ----------------------------------------------------------------
!!
!! 6. Error messages :
!!
!!      None.
!!
!! 7. Remarks :
!!
!! 8. Structure :
!!
!!    See source code.
!!
!! 9. Switches :
!!
!!    !/S  Enable subroutine tracing.
!!
!!10. Source code :
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!
!!    Parameter list
!!    --------------
      REAL,    INTENT(INOUT)  :: X(nrng,nang)
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!
!!
!!    Local Parameters
!!    ----------------
      integer              :: irng, iang
      real                 :: Y(nrng,nang)  !* dummy array used in smoothing
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!     
!!     
!!-0  Initial Y(:,:) array before it's computed
!!ini
      Y(:,:) = 0.0
!!ini---
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-1  Interpolate using simple 2 point averaging to fill in X
!!    Remeber: nang must be an even number ==> nang-2 is an even number
!!    and      nrng must be an odd  number ==> nrng-1 is an even number
!!    Example numbers used here are for nrng=35 & nang=36
!!    ------------------------------------------------------------------
!!
!!-1a For every calculated iang (1,3,5,..,nang-1=35)
!!    fill in missing irng's    (2,4,6,..,nrng-1=34)
      do 12 iang=1,nang-1,2      !* = 1,3,5,...,nang-1=35
      do 11 irng=2,nrng-1,2      !* = 2,4,6,...,nrng-1=34
        X(irng,iang) = 0.5 * ( X(irng-1,iang) + X(irng+1,iang) )
  11  continue
  12  continue
!!    ------------------------------------------------------------------
!!
!!-1b Now, for every irng (1,2,3,..,nrng  =35)
!!    fill missing iang's (2,4,6,..,nang-2=34)
      do 14 irng=1,nrng          !* 1,2,3,..,nrng  =35
      do 13 iang=2,nang-2,2      !* 2,4,6,..,nang-2=34
        X(irng,iang) = 0.5 * ( X(irng,iang-1) + X(irng,iang+1) )
  13  continue
  14  continue
!!    ------------------------------------------------------------------
!!
!!-1c for iang = nang (special case since nang is an even number)
      do 15 irng=1,nrng
        X(irng,nang) = 0.5 * ( X(irng,nang-1) + X(irng,1) )
  15  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!
!!    Skip smoothing only if ismo = 0
!!
!!
      if ( ismo.eq.0 ) goto 99
!!
!!
!!
!!-2  Smoothing the 2D array X into array Y
!!
!!-2a Smoothing the interior [2;nrng-1] x [2:nang-1]
!!-   Using 9 points averaged with equal weights.
!!-   Here use the dummy array so we don't spoil the original array.
      do 22 irng=2,nrng-1
      do 21 iang=2,nang-1
        Y(irng,iang)=(X(irng-1,iang-1)+X(irng-1,iang)+X(irng-1,iang+1) + &
                      X(irng,  iang-1)+X(irng,  iang)+X(irng,  iang+1) + &
                      X(irng+1,iang-1)+X(irng+1,iang)+X(irng+1,iang+1))/9.
  21  continue
  22  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-3  Smooth first & last line at iang=1 & iang=nang  (special cases)
!!
!!-3a Smooth line at iang = 1     (special case)
!!-   Using 9 points averaged with equal weights.
      do 31 irng=2,nrng-1
        Y(irng, 1) = (X(irng-1,nang) + X(irng-1, 1) + X(irng-1, 2) +  &
                      X(irng,  nang) + X(irng,   1) + X(irng,   2) +  &
                      X(irng+1,nang) + X(irng+1, 1) + X(irng+1, 2) )/9.
  31  continue
!!    ------------------------------------------------------------------
!!
!!-3b Smooth line at iang = nang  (special case)
!!-   Using 9 points averaged with equal weights.
      do 32 irng=2,nrng-1
        Y(irng,nang)=(X(irng-1,nang-1) +X(irng-1,nang) +X(irng-1,1) + &
                      X(irng,  nang-1) +X(irng,  nang) +X(irng,  1) + &
                      X(irng+1,nang-1) +X(irng+1,nang) +X(irng+1,1))/9.
  32  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-4  Smooth first & last col. at irng=1 & irng=nrng  (special cases)
!!
!!-4a Smooth col. at irng = 1     (low frq. can be skipped)
!!-   Using 6 points averaged with equal weights.
      do 33 iang=2,nang-1
        Y(1,iang)   = (X(1,iang-1) + X(1,iang) + X(1,iang+1) +        &
                       X(2,iang-1) + X(2,iang) + X(2,iang+1) )/6.
  33  continue
!!    ------------------------------------------------------------------
!!
!!-4b Smooth col. at irng = nrng  (high frq. can be skipped)
!!-   Using 6 points averaged with equal weights.
      do 34 iang=2,nang-1
        Y(nrng,iang)=(X(nrng-1,iang-1)+X(nrng-1,iang)+X(nrng-1,iang+1)+ &
                      X(nrng, iang-1)+X(nrng, iang)+X(nrng, iang+1) )/6.
  34  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-5  Smooth the 4 corners (optional):  <== Skip no sig. effect
!!-   Using 6 points averaged with equal weights
!!
!!-5a Corner (1, 1)
      Y(1, 1)      =( X(1,nang) + X(1, 1) + X(1, 2) +                 &
                      X(2,nang) + X(2, 1) + X(2, 2) )/6.0
!!    ------------------------------------------------------------------
!!
!!-5b Corner (nrng,1)
      Y(nrng,1)    =( X(nrng-1,nang) + X(nrng-1,1) + X(nrng-1,2) +    &
                      X(nrng,  nang) + X(nrng,  1) + X(nrng,  2) )/6.0
!!    ------------------------------------------------------------------
!!
!!-5c Corner (1,nang)
      Y(1,nang)    =( X(1,nang-1) + X(1,nang) + X(1, 1) +             &
                      X(2,nang-1) + X(2,nang) + X(2, 1) ) / 6.
!!    ------------------------------------------------------------------
!!
!!-5d Corner (nrng,nang)
      Y(nrng,nang) =( X(nrng-1,nang-1) +X(nrng-1,nang) +X(nrng-1,1) + &
                      X(nrng,  nang-1) +X(nrng,  nang) +X(nrng,  1) )/6.
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!!-6  Final, dump smoothed array Y(:,:) into X(:,:) to be returned
!!
!!-6a Done with X(:,:) re-initial before it's replaced by Y(:,:)
!!ini
      X(:,:) = 0.0
!!ini---
!!
!!-6b Dump smoothed array Y(:,:) into X(:,:) to be returned
      do 52 iang=1,nang
      do 51 irng=1,nrng
        X(irng,iang) = Y(irng,iang)
  51  continue
  52  continue
!!    Bash; can simplify in one line
!b    X(1:nrng, 1:nang) = Y(1:nrng, 1:nang)
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
  99  continue
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
      RETURN
!!
      END SUBROUTINE interp2
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL FUNCTION wkfnc ( f, dep )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!
!!    it returns: wavenumber 'k' in wkfnc
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!! 1. Purpose :
!!
!!    Calculate the Wavenumber k (rad/m) as function of
!!    frequency 'f' (Hz) and depth 'dep' (m).
!!
!! 2. Method :
!!
!!    Using what looks like a "Pade approximation" of an inversion
!!    of the linear wave dispersion relation.
!!        sigma^2 = gk*tanh(kd),  sigma = 2*pi*f
!!        Wavenumber k (rad/m) is returned in "wkfnc"
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    twopi     Real  Public    I   = TPI; WW3 2*pi=8.*atan(1.) (radians)
!!    ------------------------------------------------------------------
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: f, dep
!!
!!    Local variables
!!    ---------------
      real(KIND=8)         :: g, y, x
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
      g     = 9.806                   !* set = GRAV as in CONSTANTS
!!
      y     = ( (twopi*f)**2 ) * dep / g   !* sigma^2 d/g
!!
!!    --------------------------------------------------------------- &
      x     = y * ( y +                                               &
            1./(1.00000+y*(0.66667+y*(0.35550+y*(0.16084+y*(0.06320   &
            +y*(0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.00011)   &
            )))))))))
!!    --------------------------------------------------------------- &
!!
      x     = sqrt(x)                 !* kd
!!
      wkfnc = x / dep                 !* k
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END FUNCTION wkfnc
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL FUNCTION cgfnc ( f, dep, cvel )
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                 BIO |
!/                  |           Bash Toulany            |
!/                  |           Michael Casey           |
!/                  |           William Perrie          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         12-Apr-2016 |
!/                  +-----------------------------------+
!/
!/    01-Mar-2016 : Origination.                        ( version 5.13 )
!/
!!    it returns:  group velocity (m/s) in cgfnc
!!
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
!!
!! 1. Purpose :
!!
!!    Calculate the Group velocity Cg (m/s) as function of
!!    frequency 'f' (Hz), depth 'dep' (m) and phase speed 'cvel' (m/s)
!!
!! 2. Method :
!!
!!    This routine uses the identity
!!         sinh(2x) = 2*tanh(x)/(1-tanh(x)**2)
!!    to avoid extreme sinh(2x) for large x.
!!    thus,    2kd/sinh(2kd) = kd(1-tanh(kd)**2)/tanh(kd)
!!    Group velocity Cg (m/s) is returned in "cgfnc"
!!
!! 3. Parameters :
!!
!!    Parameter list
!!    ------------------------------------------------------------------
!!    Name     Type   Scope    I/O  Description
!!    ------------------------------------------------------------------
!!    twopi     Real  Public    I   = TPI; WW3 2*pi=8.*atan(1.) (radians)
!!    ------------------------------------------------------------------
!!
!!    --------------------------------------------------------------- &
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ----------------------------------------------------------------72
!!    ==================================================================
!!
!!
      IMPLICIT NONE
!!
!!    Parameter list
!!    --------------
      real,    intent(in)  :: f, dep, cvel
!!
!!    Local variables
!!    ---------------
      real                 :: wkd, tkd
!!    -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!    ---------------------::-----------------------------------------72
!!    ##################################################################
!!------------------------------------------------------------------------------
!!==============================================================================
!!
!!
      wkd   = twopi * f*dep/cvel           !* kd
      tkd   = tanh(wkd)                    !* tanh(kd)
      cgfnc = 0.5*cvel*(1.+wkd*(1.-tkd**2)/tkd)
!!    ------------------------------------------------------------------
!!    ==================================================================
!!
      RETURN
!!
      END FUNCTION cgfnc
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!
!!
      END MODULE W3SNL4MD
!!
!!==============================================================================
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
