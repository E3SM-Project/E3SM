! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_harmonic_analysis
!
!> \brief MPAS ocean analysis mode member: harmonic_analysis
!> \author Steven Brus
!> \date   July 2020
!> \details
!>  MPAS ocean analysis mode member: harmonic_analysis
!>
!>  This module contains subroutines to compute the harmonic decomposition
!>  of the sea-surface-height timeseries into an amplitude and phase for a 
!>  given number of tidal constituents.
!>
!>  The subroutines that compute the harmonic decomposition are taken from the
!>  ADCIRC model developed by R.A. Luettich, Jr. (University of North Carolina at Chapel Hill)
!>  and J.J. Westerink (University of Notre Dame). 
!>
!-----------------------------------------------------------------------

module ocn_harmonic_analysis

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_constants
   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_vel_tidal_potential, only: char_array,tidal_constituent_factors

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_harmonic_analysis, &
             ocn_compute_harmonic_analysis, &
             ocn_restart_harmonic_analysis, &
             ocn_finalize_harmonic_analysis

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type(char_array), dimension(37) :: constituentList
   integer :: leastSquaresSolution

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_harmonic_analysis
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_harmonic_analysis(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool, harmonicAnalysisAMPool 

      integer, pointer :: nAnalysisConstituents, nCellsSolve
      real (kind=RKIND), pointer :: harmonicAnalysisStart
      real (kind=RKIND), pointer :: harmonicAnalysisEnd
      real (kind=RKIND), dimension(:),   pointer :: ssh, sshInit
      real (kind=RKIND), dimension(:,:), pointer :: leastSquaresLHSMatrix
      real (kind=RKIND), dimension(:,:), pointer :: leastSquaresRHSVector
      real (kind=RKIND), dimension(:,:), pointer :: decomposedConstituentAmplitude 
      real (kind=RKIND), dimension(:,:), pointer :: decomposedConstituentPhase
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentFrequency 
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentNodalAmplitude
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentNodalPhase
      real (kind=RKIND), dimension(:), allocatable :: tidalConstituentAmplitude
      real (kind=RKIND), dimension(:), allocatable :: tidalConstituentLoveNumbers 
      real (kind=RKIND), dimension(:), allocatable :: tidalConstituentAstronomical 
      integer, dimension(:), allocatable :: tidalConstituentType

      type (MPAS_Time_Type) :: refTime
      integer :: iCon, iCell

      err = 0

      leastSquaresSolution = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'harmonicAnalysisAM', harmonicAnalysisAMPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

         call mpas_pool_get_array(harmonicAnalysisAMPool, 'sshInit', sshInit)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'nAnalysisConstituents', nAnalysisConstituents)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'leastSquaresLHSMatrix', leastSquaresLHSMatrix)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'leastSquaresRHSVector', leastSquaresRHSVector)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentFrequency', analysisConstituentFrequency)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentNodalAmplitude', analysisConstituentNodalAmplitude)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentNodalPhase', analysisConstituentNodalPhase)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'decomposedConstituentAmplitude', decomposedConstituentAmplitude)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'decomposedConstituentPhase', decomposedConstituentPhase)

         call mpas_set_time(refTime, dateTimeString=config_tidal_potential_reference_time)

         nAnalysisConstituents = 0
         if (config_AM_harmonicAnalysis_use_M2) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'M2'
         end if 

         if (config_AM_harmonicAnalysis_use_S2) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'S2'
         end if 
 
         if (config_AM_harmonicAnalysis_use_N2) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'N2'
         end if 

         if (config_AM_harmonicAnalysis_use_K2) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'K2'
         end if 

         if (config_AM_harmonicAnalysis_use_K1) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'K1'
         end if 

         if (config_AM_harmonicAnalysis_use_O1) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'O1'
         end if 

         if (config_AM_harmonicAnalysis_use_Q1) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'Q1'
         end if 

         if (config_AM_harmonicAnalysis_use_P1) then
           nAnalysisConstituents = nAnalysisConstituents + 1
           constituentList(nAnalysisConstituents)%constituent = 'P1'
         end if 

         allocate(tidalConstituentAmplitude(nAnalysisConstituents))
         allocate(tidalConstituentLoveNumbers(nAnalysisConstituents))
         allocate(tidalConstituentAstronomical(nAnalysisConstituents))
         allocate(tidalConstituentType(nAnalysisConstituents))

         call tidal_constituent_factors(constituentList,nAnalysisConstituents,refTime, &
                                        analysisConstituentFrequency, &
                                        tidalConstituentAmplitude, &
                                        tidalConstituentLoveNumbers, &
                                        analysisConstituentNodalAmplitude, &
                                        tidalConstituentAstronomical, &
                                        analysisConstituentNodalPhase, & 
                                        tidalConstituentType, &
                                        err)

         if (.not. config_do_restart) then
           leastSquaresRHSVector = 0.0_RKIND
           leastSquaresLHSMatrix = 0.0_RKIND
         end if
           
         do iCell = 1,nCellSSolve
           sshInit(iCell) = ssh(iCell)
         end do

         do iCon = 1,nAnalysisConstituents  
           call mpas_log_write('Constituent '//constituentList(iCon)%constituent)
           call mpas_log_write('  Frequency $r', realArgs=(/ analysisConstituentFrequency(iCon) /))
           call mpas_log_write('  Amplitude $r', realArgs=(/ tidalConstituentAmplitude(iCon) /))
           call mpas_log_write('  LoveNumbers $r', realArgs=(/ tidalConstituentLoveNumbers(iCon) /))
           call mpas_log_write('  NodalAmplitude $r', realArgs=(/ analysisConstituentNodalAmplitude(iCon) /))
           call mpas_log_write('  Astronomical argument $r', realArgs=(/ tidalConstituentAstronomical(iCon) /))
           call mpas_log_write('  NodalPhase $r', realArgs=(/ analysisConstituentNodalPhase(iCon) /))
           call mpas_log_write('  Type $i', intArgs=(/ tidalConstituentType(iCon) /))
           call mpas_log_write(' ')
         end do

         block => block % next
      end do

   end subroutine ocn_init_harmonic_analysis!}}}

!***********************************************************************
!
!  routine ocn_compute_harmonic_analysis
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_harmonic_analysis(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: harmonicAnalysisAMPool 

      integer, pointer :: nCellsSolve      
      integer, pointer :: nAnalysisConstituents
      real (kind=RKIND), dimension(:),   pointer :: ssh, sshInit
      real (kind=RKIND), dimension(:), allocatable :: sshDiff
      real (kind=RKIND), dimension(:,:), pointer :: leastSquaresLHSMatrix
      real (kind=RKIND), dimension(:,:), pointer :: leastSquaresRHSVector
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentFrequency 
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentNodalAmplitude
      real (kind=RKIND), dimension(:),   pointer :: analysisConstituentNodalPhase
      real (kind=RKIND), dimension(:,:), pointer :: decomposedConstituentAmplitude
      real (kind=RKIND), dimension(:,:), pointer :: decomposedConstituentPhase
      real (kind=RKIND), dimension(:),   pointer :: M2Amplitude
      real (kind=RKIND), dimension(:),   pointer :: M2Phase
      real (kind=RKIND), dimension(:),   pointer :: S2Amplitude
      real (kind=RKIND), dimension(:),   pointer :: S2Phase
      real (kind=RKIND), dimension(:),   pointer :: N2Amplitude
      real (kind=RKIND), dimension(:),   pointer :: N2Phase
      real (kind=RKIND), dimension(:),   pointer :: K2Amplitude
      real (kind=RKIND), dimension(:),   pointer :: K2Phase
      real (kind=RKIND), dimension(:),   pointer :: K1Amplitude
      real (kind=RKIND), dimension(:),   pointer :: K1Phase
      real (kind=RKIND), dimension(:),   pointer :: O1Amplitude
      real (kind=RKIND), dimension(:),   pointer :: O1Phase
      real (kind=RKIND), dimension(:),   pointer :: Q1Amplitude
      real (kind=RKIND), dimension(:),   pointer :: Q1Phase
      real (kind=RKIND), dimension(:),   pointer :: P1Amplitude
      real (kind=RKIND), dimension(:),   pointer :: P1Phase

      integer :: iCell
      integer :: iCon
      integer :: err_tmp
      real (kind=RKIND) :: time
      real (kind=RKIND) :: harmonicAnalysisEnd

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'harmonicAnalysisAM', harmonicAnalysisAMPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'sshInit', sshInit)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'nAnalysisConstituents', nAnalysisConstituents)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'leastSquaresLHSMatrix', leastSquaresLHSMatrix)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'leastSquaresRHSVector', leastSquaresRHSVector)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentFrequency', analysisConstituentFrequency)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentNodalAmplitude', analysisConstituentNodalAmplitude)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'analysisConstituentNodalPhase', analysisConstituentNodalPhase)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'decomposedConstituentAmplitude', decomposedConstituentAmplitude)
         call mpas_pool_get_array(harmonicAnalysisAMPool, 'decomposedConstituentPhase', decomposedConstituentPhase)

         ! get relevant time information
         time = daysSinceStartOfSim*86400_RKIND
         harmonicAnalysisEnd = config_AM_harmonicAnalysis_start_delay + config_AM_harmonicAnalysis_duration
         
         ! exit if harmonic analysis period has not begun yet
         if (daysSinceStartOfSim .lt. config_AM_harmonicAnalysis_start_delay) then
           call mpas_log_write('harmonicAnalysisAM exit: before HA period')
           return 
         end if

         ! update if within harmonic analysis period 
         if (daysSinceStartOfSim .le. harmonicAnalysisEnd) then
           allocate(sshDiff(nCellsSolve))

           do iCell = 1,nCellsSolve
             sshDiff(iCell) = ssh(iCell) - sshInit(iCell)
           end do

           CALL update_least_squares_LHS_matrix(nAnalysisConstituents, &
                                                time, &
                                                analysisConstituentFrequency, &
                                                leastSquaresLHSMatrix)

           CALL update_least_squares_RHS_vector(nAnalysisConstituents, &
                                                time, &
                                                nCellsSolve, &
                                                analysisConstituentFrequency, &
                                                sshDiff, &
                                                leastSquaresRHSVector)
           deallocate(sshDiff)
           call mpas_log_write('harmonicAnalysisAM update')
         end if

         ! solve harmonic analysis least squares system if harmonic analysis period is over
         if ((daysSinceStartOfSim .ge. harmonicAnalysisEnd) .and. (leastSquaresSolution == 0)) then
           call harmonic_analysis_solve(nCellsSolve, &
                                        nAnalysisConstituents, &
                                        leastSquaresLHSMatrix, &
                                        leastSquaresRHSVector, &
                                        analysisConstituentNodalAmplitude, &
                                        analysisConstituentNodalPhase, &
                                        decomposedConstituentAmplitude, &
                                        decomposedConstituentPhase)     

           ! copy amplitude and phase solutions into corresponding arrays
           do iCon = 1,nAnalysisConstituents
             if (constituentList(iCon)%constituent == 'M2') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'M2Amplitude', M2Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'M2Phase', M2Phase)
               M2Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               M2Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'S2') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'S2Amplitude', S2Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'S2Phase', S2Phase)
               S2Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               S2Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'N2') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'N2Amplitude', N2Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'N2Phase', N2Phase)
               N2Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               N2Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'K2') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'K2Amplitude', K2Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'K2Phase', K2Phase)
               K2Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               K2Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'K1') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'K1Amplitude', K1Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'K1Phase', K1Phase)
               K1Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               K1Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'O1') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'O1Amplitude', O1Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'O1Phase', O1Phase)
               O1Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               O1Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'Q1') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'Q1Amplitude', Q1Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'Q1Phase', Q1Phase)
               Q1Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               Q1Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
             if (constituentList(iCon)%constituent == 'P1') then
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'P1Amplitude', P1Amplitude)
               call mpas_pool_get_array(harmonicAnalysisAMPool, 'P1Phase', P1Phase)
               P1Amplitude(:) = decomposedConstituentAmplitude(iCon,:)
               P1Phase(:) = decomposedConstituentPhase(iCon,:)
             endif
           enddo

           call mpas_stream_mgr_write(domain % streamManager, streamID=config_AM_harmonicAnalysis_output_stream, forceWriteNow=.true., ierr=err_tmp)

           ! increment to indicate the system has been solved
           leastSquaresSolution = leastSquaresSolution + 1
           call mpas_log_write('harmonicAnalysisAM solve')

         else if ((daysSinceStartOfSim .ge. harmonicAnalysisEnd) .and. (leastSquaresSolution > 0)) then
           call mpas_log_write('harmonicAnalysisAM exit: past HA period')
         end if

         block => block % next
      end do

   end subroutine ocn_compute_harmonic_analysis!}}}

!***********************************************************************
!
!  routine ocn_restart_harmonic_analysis
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_harmonic_analysis(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_harmonic_analysis!}}}

!***********************************************************************
!
!  routine ocn_finalize_harmonic_analysis
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_harmonic_analysis(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      err = 0

   end subroutine ocn_finalize_harmonic_analysis!}}}

!***********************************************************************
!
!  routine update_least_squares_LHS_matrix
!
!> \brief   Update Left Hand Side Matrix 
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>   Adapted from the ADCIRC code by R.A. Luettich and J.J. Westerink 
!
!-----------------------------------------------------------------------

      SUBROUTINE update_least_squares_LHS_matrix(nfreq,TIMELOC,hafreq,ha) !{{{

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nfreq                           ! number of analysis constituents
      REAL(kind=RKIND), INTENT(IN) ::  TIMELOC               ! model time (sec)
      REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: hafreq   ! analysis contituent frequencies
      REAL(kind=RKIND), DIMENSION(:,:), INTENT(INOUT) :: ha  ! LHS matrix

      INTEGER :: I,J,I1,I2,J1,J2
      REAL(kind=RKIND) :: TF1,TF2

!     Update the Left Hand Side Matrix
!         Note: this is a symmetric matrix and therefore only store the
!         upper triangular part.  The lower part will be filled out in
!         subroutine least_squares_decompose() prior to decomposition 

      do i=1,nfreq
         do j=i,nfreq
            i1=2*i-1
            i2=i1+1
            j1=2*j-1
            j2=j1+1
            tf1=hafreq(i)*TIMELOC
            tf2=hafreq(j)*TIMELOC
            ha(i1,j1) = ha(i1,j1) + cos(tf1)*cos(tf2)
            ha(i1,j2) = ha(i1,j2) + cos(tf1)*sin(tf2)
            ha(i2,j2) = ha(i2,j2) + sin(tf1)*sin(tf2)
            if(i2.le.j1) ha(i2,j1) = ha(i2,j1) + sin(tf1)*cos(tf2)
         end do
      end do

      return
      end subroutine update_least_squares_LHS_matrix!}}}

!***********************************************************************
!
!  routine update_least_squares_RHS_vector
!
!> \brief   Update Right Hand Side Vector
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>   Adapted from the ADCIRC code by R.A. Luettich and J.J. Westerink 
!
!-----------------------------------------------------------------------

      SUBROUTINE update_least_squares_RHS_vector(nfreq,TIMEUD,NP,hafreq,GLOE,GLOELV) !{{{

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nfreq                                ! number of analysis consituents
      REAL(kind=RKIND), INTENT(IN) :: TIMEUD                      ! model time 
      INTEGER, INTENT(IN) :: NP                                   ! number of mesh cells
      REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: hafreq        ! analysis constituent frequencies
      REAL(kind=RKIND), DIMENSION(:), INTENT(IN) :: GLOE          ! sea surface height
      REAL(kind=RKIND), DIMENSION(:,:), INTENT(INOUT) :: GLOELV   ! RHS vector

      INTEGER I,N,I1,I2
      REAL(kind=RKIND) TF1,CTF1,STF1

!     Update the Right Hand Side Load Vectors
      do i=1,nfreq
         i1=2*i-1
         i2=i1+1
         tf1=hafreq(i)*TIMEUD
         ctf1 = cos(tf1)
         stf1 = sin(tf1)
         do n=1,np
            GLOELV(I1,N)=GLOELV(I1,N)+GLOE(N)*CTF1
            GLOELV(I2,N)=GLOELV(I2,N)+GLOE(N)*STF1
         end do
      end do

      return
      end subroutine update_least_squares_RHS_vector!}}}

!***********************************************************************
!
!  routine harmonic_analysis_solve 
!
!> \brief   Solve the least squares system 
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>   Adapted from the ADCIRC code by R.A. Luettich and J.J. Westerink 
!
!-----------------------------------------------------------------------

      SUBROUTINE harmonic_analysis_solve(MNP,nfreq,hmat,GLOELV,haff,haface,emagt,phaseden) !{{{

      IMPLICIT NONE

      integer, intent(in) :: MNP                                ! number of mesh cells
      integer, intent(in) :: nfreq                              ! number of analysis constituents
      real(kind=RKIND), dimension(:,:), intent(in) :: hmat      ! LHS matrix
      real(kind=RKIND), dimension(:,:), intent(in) :: GLOELV    ! RHS vector
      real(kind=RKIND), dimension(:), intent(in) :: haff        ! amplitude nodal factors
      real(kind=RKIND), dimension(:), intent(in) :: haface      ! phase nodal factors
      real(kind=RKIND), dimension(:,:), intent(out) :: emagt    ! constituent amplitudes
      real(kind=RKIND), dimension(:,:), intent(out) :: phaseden ! constituent phases

      integer J,N,K,I,I1,I2,IT,IFR
      integer mm
      real(kind=RKIND) :: convrd
      REAL(kind=RKIND),ALLOCATABLE :: PHASEE(:),EMAG(:)
      REAL(kind=RKIND),ALLOCATABLE :: hap(:),hax(:) 
      REAL(kind=RKIND),ALLOCATABLE :: ha(:,:)

      convrd=180.0_RKIND/pi

      mm = 2*nfreq
      ALLOCATE ( PHASEE(nfreq),EMAG(nfreq) )
      ALLOCATE ( hap(mm), hax(mm) )
      ALLOCATE ( ha(mm,mm) )

!     Copy LHS matrix and decompose
      do i = 1,mm
        do j = 1,mm
          ha(j,i) = hmat(j,i)
        end do
      end do
      call least_squares_decompose(nfreq,ha)

      DO N=1,MNP

!        At each node transfer the RHS vector and solve
         do k=1,mm
            hap(k) = GLOELV(k,n)
         end do
         call least_squares_solve(nfreq,ha,hap,hax)

!        Compute amplitude and phase for each constituent making sure that the
!        phase is between 0 and 360 deg.  
         do i=1,nfreq
               i1=2*i-1
               i2=i1+1
               emag(i)=sqrt(hax(i1)*hax(i1)+hax(i2)*hax(i2))
               emagt(i,n)=emag(i)/haff(i)
               if((hax(i1).eq.0.0_RKIND).and.(hax(i2).eq.0.0_RKIND)) then
                  phasee(i)=0.0_RKIND
               else
                  phasee(i) = atan2(hax(i2),hax(i1))
               endif
            phaseden(i,n)=phasee(i)+haface(i)
            phaseden(i,n)=convrd*phaseden(i,n)
            if(phaseden(i,n).lt.0.0_RKIND) phaseden(i,n)=phaseden(i,n)+360.0_RKIND
            if(phaseden(i,n).ge.360.0_RKIND) phaseden(i,n)=phaseden(i,n)-360.0_RKIND
         end do

      end do

      return
      end subroutine harmonic_analysis_solve !}}}


!***********************************************************************
!
!  routine least_squares_decompose 
!
!> \brief   Fill out symmetric matrix and decompose 
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>   Adapted from the ADCIRC code by R.A. Luettich and J.J. Westerink 
!
!-----------------------------------------------------------------------

      subroutine least_squares_decompose(nfreq,ha) !{{{

      implicit none

      integer, intent(in) :: nfreq                          ! number of analysis constituents
      real(kind=RKIND), dimension(:,:), intent(inout) :: ha ! LHS matrix

      integer i,j,ir,ire,k,jr,mm

      mm = 2*nfreq

!     Set up the lower triangular part of the LHS matrix

      do j=1,mm
         do i=j,mm
            ha(i,j)=ha(j,i)
         end do
      end do

!     Decomposition of matrix 

      do ir=1,mm
         ire=ir+1
         do j=ire,mm
            ha(ir,j)=ha(ir,j)/ha(ir,ir)
         end do
         if (ire.gt.mm) exit
         do j=ire,mm
           do k=ire,mm
              ha(k,j)=ha(k,j)-ha(k,ir)*ha(ir,j)
           end do
         end do
         do j=ire,mm
           ha(j,ir)=0.0_RKIND
         end do
      end do

      end subroutine least_squares_decompose !}}}

!***********************************************************************
!
!  routine least_squares_solve
!
!> \brief   Solves system a*x=b by l*d*l(tr) decomp in full storage mode 
!> \author  Steven Brus
!> \date    July 2020
!> \details
!>   Adapted from the ADCIRC code by R.A. Luettich and J.J. Westerink 
!
!-----------------------------------------------------------------------

      subroutine least_squares_solve(nfreq,ha,hap,hax) !{{{

      implicit none

      integer, intent(in) :: nfreq                        ! number of analysis constituents
      real(kind=RKIND), dimension(:,:), intent(in) :: ha  ! LHS matrix
      real(kind=RKIND), dimension(:), intent(in) :: hap   ! RHS vector
      real(kind=RKIND), dimension(:), intent(out) :: hax  ! solution vector

      integer idecom,i,j,ir,ire,k,jr,mm
      real(kind=RKIND),allocatable ::  c(:),y(:)

      mm = 2*nfreq

      allocate ( c(mm),y(mm) )

!     solve for y by forward substitution for l*y=p

      do ir=1,mm
         y(ir)=hap(ir)
         do jr=1,ir-1
            y(ir)=y(ir)-ha(jr,ir)*y(jr)
         end do
      end do

!     calculate c=d**(-1)*y

      do ir=1,mm
         c(ir)=y(ir)/ha(ir,ir)
      end do

!     solve for x by back-substituting for l(tr)*x=c

      do ir = mm,1,-1
         hax(ir)=c(ir)
         do jr=ir+1,mm
            hax(ir)=hax(ir)-ha(ir,jr)*hax(jr)
         end do
      end do

      end subroutine least_squares_solve!}}}

end module ocn_harmonic_analysis

! vim: foldmethod=marker
