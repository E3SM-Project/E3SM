! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Original version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_CLOUDSAT_INTERFACE
  USE MOD_COSP_CONFIG, ONLY: DBZE_BINS,CFAD_ZE_MIN,CFAD_ZE_WIDTH,SR_BINS,DBZE_MAX,       &
                             DBZE_MIN
  USE COSP_KINDS,      ONLY: wp
  USE quickbeam,       ONLY: quickbeam_init,radar_cfg,Re_MAX_BIN,Re_BIN_LENGTH
  IMPLICIT NONE

  ! Directory where LUTs will be stored
  character(len=120) :: RADAR_SIM_LUT_DIRECTORY = './'
  logical :: RADAR_SIM_LOAD_scale_LUTs_flag   = .false.
  logical :: RADAR_SIM_UPDATE_scale_LUTs_flag = .false.

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cloudsat_IN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cloudsat_IN
     integer,pointer ::            &
          Npoints,         & ! Number of horizontal grid-points
          Nlevels,         & ! Number of vertical levels
          Ncolumns           ! Number of subcolumns
     real(wp),pointer ::   &
          hgt_matrix(:,:),   & ! Height of hydrometeors (km)
          z_vol(:,:,:),      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol(:,:,:),     & ! Attenuation coefficient hydro (dB/km)
          g_vol(:,:,:),      & ! Attenuation coefficient gases (dB/km)
          g_to_vol_in(:,:)     ! Gaseous atteunation, radar to vol (dB)
     type(radar_cfg),pointer :: rcfg   ! Radar simulator configuration
  end type cloudsat_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              SUBROUTINE cosp_cloudsat_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CLOUDSAT_INIT(radar_freq,k2,use_gas_abs,do_ray,undef,nhydro,   &
                                surface_radar,rcfg,cloudsat_micro_scheme,load_LUT)
    ! INPUTS
    real(wp),intent(in) :: &
         radar_freq,  & ! Radar frequency (GHz)
         k2,          & ! |K|^2, the dielectric constant
         undef          ! Undefined
    integer,intent(in) :: &
         use_gas_abs, & ! 1 = do gaseous abs calcs, 0=no gasesous absorbtion calculated,
                        ! 2 = calculate absorption for first profile on all profiles
         do_ray,      & !
         nhydro,      & !
         surface_radar
    logical,intent(in),optional :: &
         load_LUT
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme

    ! OUTPUTS
    type(radar_cfg) :: &
         rcfg           !

    ! LOCAL VARIABLES
    character(len=240) :: LUT_file_name
    logical       :: local_load_LUT
    integer       :: j

    if (present(load_LUT)) then
       local_load_LUT = load_LUT
    else
       local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
    endif
    
    ! LUT file name
    LUT_file_name = trim(RADAR_SIM_LUT_DIRECTORY) // &
         trim(cloudsat_micro_scheme)

    ! Initialize for NEW radar-configurarion derived type (radar_cfg)
    rcfg%freq                = radar_freq
    rcfg%k2                  = k2
    rcfg%use_gas_abs         = use_gas_abs
    rcfg%do_ray              = do_ray
    rcfg%nhclass             = nhydro
    rcfg%load_scale_LUTs     = local_load_LUT
    rcfg%update_scale_LUTs   = .false.
    rcfg%scale_LUT_file_name = LUT_file_name
    rcfg%N_scale_flag        = .false.
    rcfg%fc                  = undef
    rcfg%rho_eff             = undef
    rcfg%Z_scale_flag        = .false.
    rcfg%Ze_scaled           = 0._wp
    rcfg%Zr_scaled           = 0._wp
    rcfg%kr_scaled           = 0._wp

    ! Set up Re bin "structure" for z_scaling
    rcfg%base_list(1)=0
    do j=1,Re_MAX_BIN
       rcfg%step_list(j)=0.1_wp+0.1_wp*((j-1)**1.5)
       if(rcfg%step_list(j)>Re_BIN_LENGTH) then
          rcfg%step_list(j)=Re_BIN_LENGTH
       endif
       if(j>1) then
          rcfg%base_list(j)=rcfg%base_list(j-1)+floor(Re_BIN_LENGTH/rcfg%step_list(j-1))
       endif
    enddo

    ! Set flag denoting position of radar
    if (surface_radar == 1) then
       rcfg%radar_at_layer_one = .false.
    else
       rcfg%radar_at_layer_one = .true.
    endif

  END SUBROUTINE COSP_CLOUDSAT_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              	  END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CLOUDSAT_INTERFACE
