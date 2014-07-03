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
!
! 

MODULE MOD_COSP_CONSTANTS
    IMPLICIT NONE
    
    character(len=32) :: COSP_VERSION='COSP v1.3'
    
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
    integer,parameter :: N_OUT_LIST = 47   !modified by ZYY
    ! Value for forward model result from a level that is under the ground
    real,parameter :: R_GROUND = -1.0E20
    
    ! Stratiform and convective clouds in frac_out
    integer, parameter :: I_LSC = 1, & ! Large-scale clouds
                          I_CVC = 2    ! Convective clouds
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
    
    !--- MISR constants
    integer,parameter :: MISR_N_CTH = 16

    !--- RTTOV constants
    integer,parameter :: RTTOV_MAX_CHANNELS = 20
    
    ! ISCCP tau-Pc axes
    real,parameter,dimension(7) :: ISCCP_TAU = (/0.15, 0.80, 2.45, 6.5, 16.2, 41.5, 100.0/)
    real,parameter,dimension(2,7) :: ISCCP_TAU_BNDS = reshape(source=(/0.0,0.3,0.3,1.30,1.30,3.6,3.6,9.4, &
                                                      9.4,23.0,23.0,60.0,60.0,100000.0/), shape=(/2,7/))
!     real,parameter,dimension(7) :: ISCCP_PC = (/9000., 24500., 37500., 50000., 62000., 74000., 90000./)
!     real,parameter,dimension(2,7) :: ISCCP_PC_BNDS = reshape(source=(/0.0,18000.0,18000.0,31000.0,31000.0, &
!                                44000.0,44000.0,56000.0,56000.0,68000.0,68000.0,80000.0,80000.0,100000.0/), shape=(/2,7/))
   
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
            
    !  Table hclass for quickbeam
    integer,parameter :: N_HYDRO = 9
    real :: HCLASS_TYPE(N_HYDRO),HCLASS_COL(N_HYDRO),HCLASS_PHASE(N_HYDRO), &
            HCLASS_CP(N_HYDRO),HCLASS_DMIN(N_HYDRO),HCLASS_DMAX(N_HYDRO)
    real :: HCLASS_APM(N_HYDRO),HCLASS_BPM(N_HYDRO),HCLASS_RHO(N_HYDRO), &
            HCLASS_P1(N_HYDRO),HCLASS_P2(N_HYDRO),HCLASS_P3(N_HYDRO)
    real,dimension(N_HYDRO) :: N_ax,N_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma_1,gamma_2
    
     ! HCLASS_CP is not used in the version of Quickbeam included in COSP
!                      LSL    LSI      LSR     LSS   CVL    CVI   CVR     CVS   LSG
    data HCLASS_TYPE/    5,      1,      2,      2,     5,     1,   2,      2,    2/
    data HCLASS_COL/     1,      2,      3,      4,     5,     6,   7,      8,    9/
    data HCLASS_PHASE/   0,      1,      0,      1,     0,     1,   0,      1,    1/
    data HCLASS_CP/      0,      0,      0,      0,     0,     0,   0,      0,    0/
    data HCLASS_DMIN/   -1,     -1,     -1,     -1,    -1,    -1,   -1,    -1,   -1/
    data HCLASS_DMAX/   -1,     -1,     -1,     -1,    -1,    -1,   -1,    -1,   -1/
    data HCLASS_APM/   524,  110.8,    524,     -1,   524, 110.8,  524,    -1,   -1/
    data HCLASS_BPM/     3,   2.91,      3,     -1,     3,  2.91,    3,    -1,   -1/
    data HCLASS_RHO/    -1,     -1,     -1,    100,    -1,    -1,   -1,   100,  400/
    data HCLASS_P1/     -1,     -1,   8.e6,   3.e6,    -1,    -1, 8.e6,  3.e6, 4.e6/
    data HCLASS_P2/      6,     40,     -1,      -1,    6,    40,   -1,    -1,   -1/
    data HCLASS_P3/    0.3,      2,     -1,      -1,  0.3,     2,   -1,    -1,   -1/
    
    ! Microphysical settings for the precipitation flux to mixing ratio conversion
!                    LSL    LSI       LSR       LSS   CVL    CVI       CVR       CVS      LSG
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

! !                     LSL    LSI   LSR     LSS   CVL    CVI   CVR     CVS     LSG
!     data HCLASS_TYPE/   1,     1,    1,     -1,    1,     1,    1,      1,     -1/
!     data HCLASS_COL/    1,     2,    3,      4,    5,     6,    7,      8,      9/
!     data HCLASS_PHASE/  0,     1,    0,      1,    0,     1,    0,      1,      1/
!     data HCLASS_CP/     0,     0,    0,      0,    0,     0,    0,      0,      0/
!     data HCLASS_DMIN/  -1,    -1,   -1,     -1,   -1,    -1,   -1,     -1,     -1/
!     data HCLASS_DMAX/  -1,    -1,   -1,     -1,   -1,    -1,   -1,     -1,     -1/
!     data HCLASS_APM/   -1, 0.587,   -1, 0.0444,   -1, 0.587,   -1, 0.0444,  261.8/
!     data HCLASS_BPM/   -1,  2.45,   -1,    2.1,   -1,  2.45,   -1,    2.1,      3/
!     data HCLASS_RHO/ 1000,    -1, 1000,     -1, 1000,    -1, 1000,     -1,     -1/
!     data HCLASS_P1/    -1,    -1,   -1,     -1,   -1,    -1,   -1,     -1,     -1/
!     data HCLASS_P2/    10,    40, 1000,    120,   10,    40, 1000,    120,   1000/
!     data HCLASS_P3/     3,     1,    1,      1,    3,     1,    1,      1,    3.5/
! 
!     ! Microphysical settings for the precipitation flux to mixing ratio conversion
! !                     LSL    LSI      LSR    LSS   CVL    CVI       CVR       CVS      LSG
!     data N_ax/       -1.,   -1.,     26.2,   -1.,  -1.,   -1.,     26.2,     4.e6,     -1./
!     data N_bx/       -1.,   -1.,     1.57,   -1.,  -1.,   -1.,     1.57,      0.0,     -1./
!     data alpha_x/    -1.,   -1.,      0.0,   -1.,  -1.,   -1.,      0.0,      0.0,     -1./
!     data c_x/        -1.,   -1.,    386.8,   -1.,  -1.,   -1.,    386.8,     14.3,     -1./
!     data d_x/        -1.,   -1.,     0.67,   -1.,  -1.,   -1.,     0.67,    0.416,     -1./
!     data g_x/        -1.,   -1.,      0.4,   -1.,  -1.,   -1.,      0.4,      0.4,     -1./
!     data a_x/        -1.,   -1.,    523.6,   -1.,  -1.,   -1.,    523.6,   0.0444,     -1./
!     data b_x/        -1.,   -1.,      3.0,   -1.,  -1.,   -1.,      3.0,      2.1,     -1./
!     data gamma_1/    -1.,   -1., 14.78119,   -1.,  -1.,   -1., 14.78119, 3.382827,     -1./
!     data gamma_2/    -1.,   -1.,      6.0,   -1.,  -1.,   -1.,      6.0, 2.197659,     -1./

END MODULE MOD_COSP_CONSTANTS
