! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added definitions of ISCCP axes
! Oct 2008 - H. Chepfer       - Added PARASOL_NREFL
! Jun 2010 - R. Marchand      - Modified to support quickbeam V3, added ifdef for hydrometeor definitions
! 
!
! 

#include "cosp_defs.h"
MODULE MOD_COSP_CONSTANTS
    IMPLICIT NONE

    character(len=32) :: COSP_VERSION='COSP v1.4'

    ! Indices to address arrays of LS and CONV hydrometeors
    integer,parameter :: I_LSCLIQ = 1
    integer,parameter :: I_LSCICE = 2
    integer,parameter :: I_LSRAIN = 3
    integer,parameter :: I_LSSNOW = 4
    integer,parameter :: I_CVCLIQ = 5
    integer,parameter :: I_CVCICE = 6
    integer,parameter :: I_CVRAIN = 7
    integer,parameter :: I_CVSNOW = 8
    integer,parameter :: I_LSGRPL = 9

    ! Missing value
    real,parameter :: R_UNDEF = -1.0E30

    ! Number of possible output variables
    integer,parameter :: N_OUT_LIST = 65 !+JEK 1.3 CESM mods changed from 45 to 47, I'll change this from 63 to 65
    integer,parameter :: N3D = 8
    integer,parameter :: N2D = 14
    integer,parameter :: N1D = 40

    ! Value for forward model result from a level that is under the ground
    real,parameter :: R_GROUND = -1.0E20

    ! Stratiform and convective clouds in frac_out
    integer, parameter :: I_LSC = 1, & ! Large-scale clouds
                          I_CVC = 2    ! Convective clouds

    ! Timing of different simulators, including statistics module
    integer, parameter :: N_SIMULATORS = 7
    integer,parameter :: I_RADAR = 1
    integer,parameter :: I_LIDAR = 2
    integer,parameter :: I_ISCCP = 3
    integer,parameter :: I_MISR  = 4
    integer,parameter :: I_MODIS = 5
    integer,parameter :: I_RTTOV = 6
    integer,parameter :: I_STATS = 7
    character*32, dimension(N_SIMULATORS) :: SIM_NAME = (/'Radar','Lidar','ISCCP','MISR ','MODIS','RTTOV','Stats'/)
!   use of timer causes pgi_acc to fail 
!   integer,dimension(N_SIMULATORS) :: tsim
!   data tsim/N_SIMULATORS*0/

    !--- Radar constants
    ! CFAD constants
    integer,parameter :: DBZE_BINS     =   15   ! Number of dBZe bins in histogram (cfad)
    real,parameter    :: DBZE_MIN      = -100.0 ! Minimum value for radar reflectivity
    real,parameter    :: DBZE_MAX      =   80.0 ! Maximum value for radar reflectivity
    real,parameter    :: CFAD_ZE_MIN   =  -50.0 ! Lower value of the first CFAD Ze bin
    real,parameter    :: CFAD_ZE_WIDTH =    5.0 ! Bin width (dBZe)


    !--- Lidar constants
    ! CFAD constants
    integer,parameter :: SR_BINS       =   15
    integer,parameter :: DPOL_BINS     =   6
    real,parameter    :: LIDAR_UNDEF   =   999.999

    ! Other constants
    integer,parameter :: LIDAR_NCAT    =   4
    integer,parameter :: PARASOL_NREFL =   5 ! parasol
    real,parameter,dimension(PARASOL_NREFL) :: PARASOL_SZA = (/0.0, 20.0, 40.0, 60.0, 80.0/)
    real,parameter    :: DEFAULT_LIDAR_REFF = 30.0e-6 ! Default lidar effective radius

    integer,parameter :: LIDAR_NTEMP = 40
    real,parameter,dimension(LIDAR_NTEMP) :: LIDAR_PHASE_TEMP=(/-91.5,-88.5,-85.5,-82.5,-79.5,-76.5,-73.5,-70.5,-67.5,-64.5, &
                   -61.5,-58.5,-55.5,-52.5,-49.5,-46.5,-43.5,-40.5,-37.5,-34.5, &
                   -31.5,-28.5,-25.5,-22.5,-19.5,-16.5,-13.5,-10.5, -7.5, -4.5, &
                    -1.5,  1.5,  4.5,  7.5, 10.5, 13.5, 16.5, 19.5, 22.5, 25.5/)
    real,parameter,dimension(2,LIDAR_NTEMP) :: LIDAR_PHASE_TEMP_BNDS=reshape(source=&
              (/-273.15,-90.,-90.,-87.,-87.,-84.,-84.,-81.,-81.,-78., &
                   -78.,-75.,-75.,-72.,-72.,-69.,-69.,-66.,-66.,-63., &
                   -63.,-60.,-60.,-57.,-57.,-54.,-54.,-51.,-51.,-48., &
                   -48.,-45.,-45.,-42.,-42.,-39.,-39.,-36.,-36.,-33., &
                   -33.,-30.,-30.,-27.,-27.,-24.,-24.,-21.,-21.,-18., &
                   -18.,-15.,-15.,-12.,-12., -9., -9., -6., -6., -3., &
                    -3.,  0.,  0.,  3.,  3.,  6.,  6.,  9.,  9., 12., &
                    12., 15., 15., 18., 18., 21., 21., 24., 24.,100./),shape=(/2,40/))

    !--- MISR constants
    integer,parameter :: MISR_N_CTH = 16

    !--- RTTOV constants
    integer,parameter :: RTTOV_MAX_CHANNELS = 20

    ! ISCCP tau-Pc axes
    real,parameter,dimension(7) :: ISCCP_TAU = (/0.15, 0.80, 2.45, 6.5, 16.2, 41.5, 100.0/)
    real,parameter,dimension(2,7) :: ISCCP_TAU_BNDS = reshape(source=(/0.0,0.3,0.3,1.30,1.30,3.6,3.6,9.4, &
                                                      9.4,23.0,23.0,60.0,60.0,100000.0/), shape=(/2,7/))

    real,parameter,dimension(7) :: ISCCP_PC = (/90000., 74000., 62000., 50000., 37500., 24500., 9000./)
    real,parameter,dimension(2,7) :: ISCCP_PC_BNDS = reshape(source=(/100000.0,80000.0,80000.0,68000.0,68000.0,56000.0 &
                               ,56000.0,44000.0,44000.0,31000.0,31000.0,18000.0,18000.0,0.0/), shape=(/2,7/))

    real,parameter,dimension(MISR_N_CTH) :: MISR_CTH = 1000.0*(/ 0., 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, &
                                            4.5, 6., 8., 10., 12., 14.5, 16., 18./)
    real,parameter,dimension(2,MISR_N_CTH) :: MISR_CTH_BNDS = 1000.0*reshape(source=(/ &
                                            -99.0,  0.0,       0.0,  0.5,       0.5,  1.0,      1.0,  1.5, &
                                              1.5,  2.0,       2.0,  2.5,       2.5,  3.0,      3.0,  4.0, &
                                              4.0,  5.0,       5.0,  7.0,       7.0,  9.0,      9.0, 11.0, &
                                             11.0, 13.0,      13.0, 15.0,      15.0, 17.0,     17.0, 99.0/), &
                                             shape=(/2,MISR_N_CTH/))


    !
    ! The following code was modifed by Roj with implementation of quickbeam V3
    !   (1) use ifdef to support more than one microphyscis scheme 
    !   (2) added constants  microphysic_scheme_name, LOAD_scale_LUTs, and SAVE_scale_LUTs 
    !

    ! directory where LUTs will be stored
    character*120 :: RADAR_SIM_LUT_DIRECTORY = './'

#ifdef MMF_V3_SINGLE_MOMENT

    !        
    !  Table hclass for quickbeam to support one-moment (bulk) microphysics scheme used by MMF V3.0 & V3.5
    !

    !
    ! NOTE:  if ANY value in this section of code is changed, the existing LUT 
    !        (i.e., the associated *.dat file) MUST be deleted so that a NEW
    !        LUT will be created !!!
    !
    character*120 :: RADAR_SIM_MICROPHYSICS_SCHEME_NAME = 'MMF_v3_single_moment'

    logical :: RADAR_SIM_LOAD_scale_LUTs_flag   = .false.
    logical :: RADAR_SIM_UPDATE_scale_LUTs_flag = .false.
    integer,parameter :: N_HYDRO = 9

    integer :: HCLASS_TYPE(N_HYDRO),HCLASS_PHASE(N_HYDRO)

    real :: HCLASS_DMIN(N_HYDRO),HCLASS_DMAX(N_HYDRO), &
            HCLASS_APM(N_HYDRO),HCLASS_BPM(N_HYDRO),HCLASS_RHO(N_HYDRO), &
            HCLASS_P1(N_HYDRO),HCLASS_P2(N_HYDRO),HCLASS_P3(N_HYDRO)

    ! HCLASS_CP is not used in the version of Quickbeam included in COSP
    !                   LSL    LSI      LSR     LSS   CVL    CVI   CVR     CVS   LSG
    data HCLASS_TYPE/    5,      1,      2,      2,     5,     1,   2,      2,    2/
    data HCLASS_PHASE/   0,      1,      0,      1,     0,     1,   0,      1,    1/
    data HCLASS_DMIN/   -1,     -1,     -1,     -1,    -1,    -1,   -1,    -1,   -1/
    data HCLASS_DMAX/   -1,     -1,     -1,     -1,    -1,    -1,   -1,    -1,   -1/
    data HCLASS_APM/   524,  110.8,    524,     -1,   524, 110.8,  524,    -1,   -1/
    data HCLASS_BPM/     3,   2.91,      3,     -1,     3,  2.91,    3,    -1,   -1/
    data HCLASS_RHO/    -1,     -1,     -1,    100,    -1,    -1,   -1,   100,  400/
    data HCLASS_P1/     -1,     -1,   8.e6,   3.e6,    -1,    -1, 8.e6,  3.e6, 4.e6/
    data HCLASS_P2/      6,     40,     -1,      -1,    6,    40,   -1,    -1,   -1/
    data HCLASS_P3/    0.3,      2,     -1,      -1,  0.3,     2,   -1,    -1,   -1/

    ! NOTES on HCLASS variables
    !
    ! TYPE - Set to
    ! 1 for modified gamma distribution,
    ! 2 for exponential distribution,
    ! 3 for power law distribution,
    ! 4 for monodisperse distribution,
    ! 5 for lognormal distribution.

    ! PHASE - Set to 0 for liquid, 1 for ice.

    ! DMIN - The minimum drop size for this class (micron), ignored for monodisperse.
    ! DMAX - The maximum drop size for this class (micron), ignored for monodisperse.
    ! Important note: The settings for DMIN and DMAX are
    ! ignored in the current version for all distributions except for power
    ! law. Except when the power law distribution is used, particle size
    ! is fixed to vary from zero to infinity, a restriction that is expected
    ! to be lifted in future versions. A placeholder must still be specified
    ! for each.

    ! Density of particles is given by apm*D^bpm or a fixed value rho. ONLY specify ONE of these two!!
    ! APM - The alpha_m coefficient in equation (1) (kg m**-beta_m )
    ! BPM - The beta_m coefficient in equation (1), see section 4.1.

    ! RHO - Hydrometeor density (kg m-3 ).

    ! P1, P2, P3 - are default distribution parameters that depend on the type
    ! of distribution (see quickmbeam documentation for more information)
    !
    ! Modified Gamma (must set P3 and one of P1 or P2)
    ! P1 - Set to the total particle number concentration Nt /rho_a (kg-1 ), where
    ! rho_a is the density of air in the radar volume.
    ! P2 - Set to the particle mean diameter D (micron).
    ! P3 - Set to the distribution width nu.
    !
    ! Exponetial (set one of)
    ! P1 - Set to a constant intercept parameter N0 (m-4).
    ! P2 - Set to a constant lambda (micron-1).
    !
    ! Power Law
    ! P1 - Set this to the value of a constant power law parameter br
    !
    ! Monodisperse
    ! P1 - Set to a constant diameter D0 (micron) = Re.
    !
    ! Log-normal (must set P3 and one of P1 or P2)
    ! P1 - Set to the total particle number concentration Nt /rho_a (kg-1 )
    ! P2 - Set to the geometric mean particle radius rg (micron).
    ! P3 - Set to the natural logarithm of the geometric standard deviation.
    !


    real,dimension(N_HYDRO) :: N_ax,N_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma_1,gamma_2,gamma_3,gamma_4

    ! Microphysical settings for the precipitation flux to mixing ratio conversion
    !                LSL    LSI       LSR       LSS   CVL    CVI       CVR       CVS      LSG
    data N_ax/       -1.,   -1.,     8.e6,     3.e6,  -1.,   -1.,     8.e6,     3.e6,     4.e6/
    data N_bx/       -1.,   -1.,      0.0,      0.0,  -1.,   -1.,      0.0,      0.0,      0.0/
    data alpha_x/    -1.,   -1.,      0.0,      0.0,  -1.,   -1.,      0.0,      0.0,      0.0/
    data c_x/        -1.,   -1.,    842.0,     4.84,  -1.,   -1.,    842.0,     4.84,     94.5/
    data d_x/        -1.,   -1.,      0.8,     0.25,  -1.,   -1.,      0.8,     0.25,      0.5/
    data g_x/        -1.,   -1.,      0.5,      0.5,  -1.,   -1.,      0.5,      0.5,      0.5/
    data a_x/        -1.,   -1.,    524.0,    52.36,  -1.,   -1.,    524.0,    52.36,   209.44/
    data b_x/        -1.,   -1.,      3.0,      3.0,  -1.,   -1.,      3.0,      3.0,      3.0/
    data gamma_1/    -1.,   -1., 17.83725, 8.284701,  -1.,   -1., 17.83725, 8.284701, 11.63230/
    data gamma_2/    -1.,   -1.,      6.0,      6.0,  -1.,   -1.,      6.0,      6.0,      6.0/
    data gamma_3/    -1.,   -1.,      2.0,      2.0,  -1.,   -1.,      2.0,      2.0,      2.0/
    data gamma_4/    -1.,   -1.,      6.0,      6.0,  -1.,   -1.,      6.0,      6.0,      6.0/



#endif


#ifdef MMF_V3p5_TWO_MOMENT

    !
    !  Table hclass for quickbeam to support two-moment "morrison" microphysics scheme used by V3.5 (SAM 6.8)
    !
    !  This Number concentriation Np in [1/kg] MUST be input to COSP/radar simulator
    !
    !  NOTE:  Be sure to check that the ice-density (rho) set it this tables matches what you used
    !

    !
    ! NOTE:  if ANY value in this section of code is changed, the existing LUT 
    !        (i.e., the associated *.dat file) MUST be deleted so that a NEW
    !        LUT will be created !!!
    !
    character*120 :: RADAR_SIM_MICROPHYSICS_SCHEME_NAME = 'MMF_v3.5_two_moment'

    logical :: RADAR_SIM_LOAD_scale_LUTs_flag   = .false.
    logical :: RADAR_SIM_UPDATE_scale_LUTs_flag = .false.

    integer,parameter :: N_HYDRO = 9

    integer :: HCLASS_TYPE(N_HYDRO),HCLASS_PHASE(N_HYDRO) 

    real :: HCLASS_DMIN(N_HYDRO),HCLASS_DMAX(N_HYDRO), &           
            HCLASS_APM(N_HYDRO),HCLASS_BPM(N_HYDRO),HCLASS_RHO(N_HYDRO), &
            HCLASS_P1(N_HYDRO),HCLASS_P2(N_HYDRO),HCLASS_P3(N_HYDRO)

    ! HCLASS_CP is not used in the version of Quickbeam included in COSP
    !                   LSL    LSI      LSR     LSS   CVL    CVI   CVR     CVS   LSG
    data HCLASS_TYPE/    1,      1,      1,      1,     1,     1,    1,      1,    1/
    data HCLASS_PHASE/   0,      1,      0,      1,     0,     1,    0,      1,    1/
    data HCLASS_DMIN/   -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/
    data HCLASS_DMAX/   -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/
    data HCLASS_APM/   524,     -1,    524,     -1,   524,    -1,  524,     -1,   -1/
    data HCLASS_BPM/     3,     -1,      3,     -1,     3,    -1,    3,     -1,   -1/
    data HCLASS_RHO/    -1,    500,     -1,    100,    -1,   500,   -1,    100,  900/
    data HCLASS_P1/     -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/
    data HCLASS_P2/     -1,     -1,     -1,     -1,    -1,    -1,   -1,     -1,   -1/
    data HCLASS_P3/     -2,      1,      1,      1,    -2,     1,    1,      1,    1/
    ! Note: value of "-2" for HCLASS_P3 uses martin 1994 parameteriztion of gamma function width with Number concentration
#endif

END MODULE MOD_COSP_CONSTANTS
